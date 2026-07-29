import pandas as pd
import glob
import statsmodels.stats.multitest as multi
import numpy as np


def filter_data_majority_direction_new(filter_mod_rate=True,
                                       filter_2_of_3=False,
                                       filter_kmers=True,
                                       filter_coverage=False,
                                       filter_p_adj=True):
    """
    Info:
    We can consider only one modification type per k-mer by finding the majority
    mod_assignment of each k-mer. For example, the majority of the modification means
    of GGACT (mu_mod) is lower than the non- modification counterpart (mu_unmod). We
    can filter out those positions whose mod_assigment values are not in line with those of the
    majority in order to restrict ourselves with one modification type per kmer in the analysis
    """

    print(file)
    groups = file.split("\\")[-2].split("_")[-1].split("-")
    print(groups)
    g1 = groups[0]
    g2 = groups[1]
    # take s3 as "control"
    data = pd.read_table(file)
    data["kmer"] = data["kmer"].apply(lambda x: x.replace("T", "U"))
    data["middle_base"] = data["kmer"].apply(lambda x: x[2] if len(x) == 5 else x)

    print("before filtering:", len(data))
    # filter nan >= 2/3
    if filter_2_of_3:
            data.iloc[:, 6:9] = data.iloc[:, 6:9].fillna(0)
            data.iloc[:, 9:12] = data.iloc[:, 9:12].fillna(0)
            data = data[data.iloc[:, 9:12].apply(lambda x: np.count_nonzero(x), axis=1) >= 2]
            data = data[data.iloc[:, 6:9].apply(lambda x: np.count_nonzero(x), axis=1) >= 2]
    else:
        # alternative filter:
        data = data.dropna()
    print("dropna:", len(data))

    if filter_coverage:
        # filter darasets coverage > 10  # 30 / 100
        coverage_cols = []
        for col in data.columns:
            if "coverage" in col:
                # tbl = tbl.dropna(subset=col)
                coverage_cols.append(col)
        # tbl["coverage"] = tbl[coverage_cols].min(axis=1)
        data = data[data[coverage_cols].min(axis=1) >= coverage]  # default 30

    print("coverage:", len(data))
    # considered dataset after qc

    # filter p_adj
    pcol = f"pval_{g1}_vs_{g2}"
    # Perform Benjamini-Hochberg correction / fdr correction
    # see https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
    _, adj_val = multi.fdrcorrection(data[pcol])
    padj_col = f"adj_pval_{g1}_vs_{g2}"
    data[padj_col] = adj_val
    # if file.endswith(".txt"):
    #     data.to_csv(file.replace(".txt", f"_before_padj_filtered_{coverage}.tsv"), sep="\t", index=False)
    # elif file.endswith(".table"):
    #     data.to_csv(file.replace(".table", f"_before_padj_filtered_{coverage}.tsv"), sep="\t", index=False)

    if filter_p_adj:
        data = data[data[padj_col] <= 0.05]
    else:
        data = data[data[pcol] <= 0.05]

    print("padj_col:", len(data))
    # selecting the middle base in the kmer column
    # filter by modification rate difference
    data.to_csv(folder + r"\majority_direction_kmer_diffmod_before_padj.tsv", sep="\t", index=False)

    if filter_mod_rate:
        data = data[abs(data[f"diff_mod_rate_{g1}_vs_{g2}"]) >= 0.25]

    print("diff_mod_rate:", len(data))
    # if file.endswith(".tsv"):
    #     data.to_csv(file.replace(".tsv", f"_filtered_pre_proc_{coverage}.tsv"), sep="\t", index=False)
    # if file.endswith(".txt"):
    #     data.to_csv(file.replace(".txt", f"_filtered_pre_proc_{coverage}.tsv"),  sep="\t", index=False)
    # if file.endswith(".table"):
    #     data.to_csv(file.replace(".table", f"_filtered_pre_proc_{coverage}.tsv"), sep="\t", index=False)


if __name__ == '__main__':
    coverage = 15
    wd = r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore"
    folder = wd + r"\direct_RNA_seq\Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4"
    file = folder + r"\majority_direction_kmer_diffmod.table"

    filter_data_majority_direction_new(filter_mod_rate=True,
                                       filter_2_of_3=False,
                                       filter_p_adj=True,
                                       filter_coverage=False)

