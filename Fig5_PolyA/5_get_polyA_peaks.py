import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gc  # Garbage collector
import re
import statsmodels.stats.multitest as multi
from scipy.stats import ks_2samp, gaussian_kde, mannwhitneyu, ttest_ind
from scipy.stats import pearsonr, spearmanr
from transcriptID2gene_names_biomart import *


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def kde_numpy(data, bandwidth=0.3, grid_size=200, grid_range=None):
    data = np.asarray(data)
    n = len(data)

    # Set the evaluation grid
    if grid_range is None:
        grid_min = data.min() - 3 * bandwidth
        grid_max = data.max() + 3 * bandwidth
    else:
        grid_min, grid_max = grid_range

    grid = np.linspace(grid_min, grid_max, grid_size)

    # Gaussian kernel for each data point
    kernel_values = np.zeros_like(grid)

    for xi in data:
        kernel_values += np.exp(-0.5 * ((grid - xi) / bandwidth) ** 2)

    # Normalize the result
    density = kernel_values / (n * bandwidth * np.sqrt(2 * np.pi))

    return grid, density


def get_max_peak_poly_tail_len(overwrite2=False):

    for file in glob.glob("DATA/*/polya_results_cdna.all.tsv"):
        print(file)
        file_name  = file.split("\\")[-2]
        file_out_temp = r"Max_PolyA_Peak_data/" + file_name + "_polyA_peak.tsv"
        print(file_out_temp)
        if not os.path.exists(file_out_temp) or overwrite2:
            print(file)
            condition = file.split("\\")[1].split("_")[0]
            print(condition)
            data_ = pd.read_table(file)

            data_ = data_[data_.qc_tag == "PASS"]
            coverage = []
            t_, pa_ = [], []
            t_idx = 0
            num_transcripts = len(data_["contig"].unique())

            length = 50
            for transcript, tbl in data_.groupby("contig"):
                t_idx += 1

                percent = t_idx / num_transcripts

                bar = '█' * int(length * percent) + '-' * (length - int(length * percent))
                print(f'\rLoading |{bar}| {percent:.0%}', end='', flush=True)
                x_vals, y_vals = kde_numpy(tbl["polya_length"],
                                           bandwidth=20,
                                           grid_size=200)
                peak_idx = np.argmax(y_vals)
                peak_x = x_vals[peak_idx]
                # print(peak_x)
                pa_.append(int(peak_x))
                t_.append(transcript)
                coverage.append(len(tbl["polya_length"].dropna()))

            out = pd.DataFrame({"contig": t_,
                                "polya_length": pa_,
                                "coverage": coverage})

            file = file.split("\\")[1]
            out["run"] = file[:5]
            out["condition"] = file[:2]
            print("\tDONE AND SAVE")
            print(out.head(5))
            out.to_csv(file_out_temp, sep="\t", index=False)


def check_protein_polya_diff_corr():
    polya = pd.read_table("polyA_gene_level_sumstats.tsv")
    print(polya.head(5))
    polya = pd.merge(polya, master, left_on="genes", right_on="gene_name", how="left")
    print(polya)

    plt.style.use('ggplot')
    plt.figure(figsize=(4, 4))
    num_dmr = 1
    polya = polya[
        (polya["dmr_count"] >= num_dmr)
        &
        ((polya["protein_qval"] <= 0.05)
        |
        (polya["ont_gene_qval"] <= 0.05))
    ]

    polya = polya.dropna(subset=["diff", "ont_gene_log2FC_C4vsC3", "protein_log2FC_C4vsC3"])
    r_rna, p_rna = pearsonr(polya["diff"], polya["ont_gene_log2FC_C4vsC3"])
    r_protein, p_protein = pearsonr(polya["diff"], polya["protein_log2FC_C4vsC3"])
    print("rna protein")
    print(r_rna, p_rna,
          r_protein, p_protein)
    sns.regplot(data=polya, x="diff", y="ont_gene_log2FC_C4vsC3",
                label=f"Gene level\nr={round(r_rna, 2)}, p={p_rna:.1e}",
                color="#909090",
                scatter_kws={"alpha": 0.8}, line_kws={"color": "#404040"})
    sns.regplot(data=polya, x="diff", y="protein_log2FC_C4vsC3",
                label=f"Protein level\nr={round(r_protein, 2)}, p={p_protein:.1e}",
                color="#303030",
                scatter_kws={"alpha": 0.8}, line_kws={"color": "k"})

    plt.legend(loc="lower left", fontsize=9)
    if num_dmr == 0:
        plt.title(f"Genes without DMRs")
    else:
        plt.title(f"Genes with ≥{num_dmr} DMRs")
    plt.xlabel("∆polyA length")
    plt.ylabel(f"Log$_2$(fold change)")
    plt.ylim(-2.1, 2.1)
    plt.yticks(np.arange(-2, 3, 1), np.arange(-2, 3, 1))
    plt.tight_layout()
    plt.savefig(f"polyA_diff_C4_C3_vs_logFC_C4_C3_num_dmr_{num_dmr}.svg")
    plt.show()
    return

def add_rowwise_zscores(
        df: pd.DataFrame,
        columns: list[str],
        suffix: str = "_zscore"
) -> pd.DataFrame:
    """
    Z-score selected columns across each row.

    For each transcript:
        z = (value - row mean) / row standard deviation

    The standard deviation is calculated with ddof=1.
    Rows with zero variance receive NaN z-scores.
    """
    result = df.copy()

    values = result[columns].apply(
        pd.to_numeric,
        errors="coerce"
    )

    row_mean = values.mean(axis=1, skipna=True)
    row_sd = values.std(axis=1, ddof=1, skipna=True)

    # Prevent division by zero for rows where all values are identical
    row_sd = row_sd.replace(0, np.nan)

    zscores = values.sub(row_mean, axis=0).div(row_sd, axis=0)

    zscore_columns = [
        f"{column}{suffix}"
        for column in columns
    ]

    zscores.columns = zscore_columns

    return pd.concat([result, zscores], axis=1)


if __name__ == '__main__':

    get_max_peak_poly_tail_len()

    overwrite = True
    coverage = 10
    filter_nans = True
    allow_one_nan_per_condition = False
    pvals, diffs = [], []

    master = pd.read_table("../0_TABLES/DMR_MASTER_TABLE_Gene_Level_C3_vs_C4_padj_DMR_025_PAPER.tsv")
    master = master.sort_values(by="protein_log2FC_C4vsC3")

    # combine data sets
    file_out = "polyA_max_peak_per_transctript_matrix.tsv"
    if not os.path.exists(file_out) or overwrite:
        print("overwrite")
        polya_data_ = pd.concat(pd.read_table(file) for file in glob.glob("Max_PolyA_Peak_data/*polyA_peak.tsv"))
        # polya_data_["contig"] = polya_data_["contig"].apply(lambda x: x.split(".")[0])
        polya_data_.condition = polya_data_.condition.apply(lambda x: x.replace("s", "C"))
        polya_data_.run = polya_data_.run.apply(lambda x: x.replace("s", "C").replace("r", "R"))
        temp = []
        cols = []
        for repl, tbl in polya_data_.groupby("run"):
            tbl = tbl[["contig", "polya_length", "coverage"]].set_index("contig")
            tbl = tbl[tbl.coverage > coverage]
            temp.append(tbl)
            cols.extend(["polya_length_" + repl, "coverage_" + repl])
        polya_data_ = pd.concat(temp, axis=1)
        polya_data_ = pd.DataFrame(polya_data_.to_numpy(), columns=cols, index=polya_data_.index)
        polya_data_ = polya_data_[sorted(polya_data_.columns)]
        print(polya_data_.head(5))
        polya_data_.to_csv(file_out, sep="\t", index=True)
    else:
        polya_data_ = pd.read_table(file_out, index_col=0)


    # filter out missing values > only one nan per condition and transctript is allowed
    if filter_nans:
        if allow_one_nan_per_condition:
            for cond in ["C3", "C4"]:
                polya_data_ = polya_data_[polya_data_[
                                  [f"polya_length_{cond}_R1",
                                   f"polya_length_{cond}_R2",
                                   f"polya_length_{cond}_R3"]].isna().sum(axis=1) <= 1]
        else:
            for cond in [# "C1", "C2",
                         "C3", "C4"]:
                polya_data_ = polya_data_[polya_data_[
                                  [f"polya_length_{cond}_R1",
                                   f"polya_length_{cond}_R2",
                                   f"polya_length_{cond}_R3"]].isna().sum(axis=1) == 0]


    print(polya_data_.head(5))

    polya_data_ = polya_data_.reset_index()
    polya_data_["contig"] = polya_data_["contig"].apply(lambda x: x.split(".")[0])

    # filter for xpore_data transcripts
    xpore_file = "../Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_GENENAMES.table"
    xpore_data = pd.read_table(xpore_file, low_memory=False)
    xpore_data["id"] = xpore_data["id"].apply(lambda x: x.split(".")[0])
    polya_data_ = polya_data_[polya_data_.contig.isin(xpore_data["id"])]

    # annotate gene names and group transcript by genes
    polya_data_ = get_gene_name_from_biomart(polya_data_, file_out=file_out.replace(".tsv", "_GENENAMES.tsv"),
                                             id_col="contig")
    polya_data_ = polya_data_.rename(columns={"Gene name": "gene_name"})
    print(polya_data_[polya_data_["gene_name"] == "NAMPT"])
    # polya_data_["gene_name"] = polya_data_["contig"].map(dict(zip(xpore_data["id"], xpore_data["Gene name"])))

    cols = [f"polya_length_C3_R1",
            f"polya_length_C3_R2",
            f"polya_length_C3_R3",
            f"polya_length_C4_R1",
            f"polya_length_C4_R2",
            f"polya_length_C4_R3"]

    polya_data_ = add_rowwise_zscores(
        df=polya_data_,
        columns=cols
    )



    print(len(polya_data_))
    polya_data_ = polya_data_.dropna(subset=["gene_name"])
    print(len(polya_data_))

    gene_level = pd.read_table(r"..\0_TABLES\DMR_MASTER_TABLE_Gene_Level_C3_vs_C4_padj_DMR_025_PAPER.tsv")
    gene_level["mean_gene_level"] = gene_level[["gene_level_ont_C3", "gene_level_ont_C4"]].mean(axis=1)
    print(gene_level["mean_gene_level"].describe())
    # Filter low-abundant transcripts > coverage > log2 >= 0 in ONT dataset
    polya_data_ = polya_data_[polya_data_["gene_name"].isin(gene_level.loc[gene_level["mean_gene_level"] > 2,
    "gene_name"].values)]

    print(polya_data_[polya_data_["gene_name"] == "NAMPT"])
    print(len(polya_data_))


    def weighted_mean(values, weights):
        """
        Weighted mean excluding rows with missing values,
        missing weights, or non-positive weights.
        """
        values = pd.to_numeric(values, errors="coerce")
        weights = pd.to_numeric(weights, errors="coerce")

        valid = (
                values.notna()
                & weights.notna()
                & (weights > 0)
        )

        if not valid.any():
            return np.nan

        return np.average(
            values[valid],
            weights=weights[valid]
        )

    search_genes = []
    for gene, data__per_gene in polya_data_.groupby("gene_name"):
        c3 = ["polya_length_C3_R1_zscore",  "polya_length_C3_R2_zscore",  "polya_length_C3_R3_zscore"]
        c4 = ["polya_length_C4_R1_zscore",  "polya_length_C4_R2_zscore",  "polya_length_C4_R3_zscore"]
        c3_weights = ["coverage_C3_R1",  "coverage_C3_R2",  "coverage_C3_R3"]
        c4_weights = ["coverage_C4_R1",  "coverage_C4_R2",  "coverage_C4_R3"]


        c3_weights =  [x for x in data__per_gene[c3_weights].values.ravel() if not np.isnan(x)]
        c3_weights = [x / max(c3_weights) for x in c3_weights]
        c4_weights = [x for x in data__per_gene[c4_weights].values.ravel() if not np.isnan(x)]
        c4_weights = [x / max(c4_weights) for x in c4_weights]
        c3 = [x for x in data__per_gene[c3].values.ravel() if not np.isnan(x)]
        c4 = [x for x in data__per_gene[c4].values.ravel() if not np.isnan(x)]
        # c4 = data__per_gene[c4].values.ravel()

        c3_weighted_mean = np.average(c3, weights=c3_weights)
        c4_weighted_mean = np.average(c4, weights=c4_weights)

        # delta = np.nanmean(c4) - np.nanmean(c3)
        delta = c4_weighted_mean - c3_weighted_mean

        if delta == 0:
            pass
        stat, p = ttest_ind([x for x in c3 if not np.isnan(x)],
                            [x for x in c4 if not np.isnan(x)],
                            equal_var=False
                            )  # alternative="greater"
        pvals.append(p)
        diffs.append(delta)
        search_genes.append(gene)

    df = pd.DataFrame(data={"genes": search_genes, "pval": pvals, "diff": diffs})
    print(df[df["genes"] == "NAMPT"])
    df = df.dropna()
    _, adj_val = multi.fdrcorrection(df["pval"])
    df["qval"] = adj_val
    df = df.sort_values(by="pval", ascending=True)
    print(df)
    print(df[df["pval"] <= 0.05])
    print(df[df["genes"] == "NAMPT"])
    # df = df[abs(df["diff"]) < 5]
    df.to_csv("polyA_gene_level_sumstats_NEW.tsv", sep="\t", index=False)