import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy import stats
import numpy as np
from matplotlib import cm
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import mutual_info_regression
import scipy.cluster.hierarchy as sch
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr, pearsonr, kendalltau


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def get_median_poly_tail_len():
    file_out = "direct_RNA_seq/PolyA/DATA/polyA_medians_per_cond.tsv"
    if not os.path.exists(file_out):
        out2 = []
        for file in glob.glob("direct_RNA_seq/PolyA/DATA/*/polya_results_cdna.all.tsv"):
            print(file)
            condition = file.split("\\")[1].split("_")[0]
            print(condition)
            data_ = pd.read_table(file)
            # print(tbl.head(5))
            data_ = data_[data_.qc_tag == "PASS"]
            out = data_.groupby("contig")["polya_length"].median().reset_index()
            # out = data_
            file = file.split("\\")[1]
            out["run"] = file[:5]
            out["condition"] = file[:2]
            out2.append(out)
        data_ = pd.concat(out2)
        data_.to_csv("direct_RNA_seq/PolyA/DATA/polyA_medians_per_cond_ungrouped.tsv", sep="\t", index=False)
        data_ = data_[["contig", "polya_length", "condition"]].groupby(["contig",
                                                                        "condition"
                                                                        ]).mean(numeric_only=True).reset_index()
        data_.to_csv(file_out, sep="\t", index=False)
    else:
        data_ = pd.read_table(file_out)
    data_["contig"] = data_["contig"].apply(lambda x: x.split(".")[0])
    data_.condition = data_.condition.apply(lambda x: x.replace("s", "C"))
    return data_


def create_xpore_master_table():

    # load tables
    xpore = pd.read_table(xpore_file)
    mrna_level = pd.read_table(mrna_level_file)
    mrna_level_ont = pd.read_table(mrna_level_ont_file)
    polya_tail_len = pd.read_table(polya_file)

    ##########################################################
    # filter and process tables
    ##########################################################
    xpore["id"] = xpore["id"].apply(lambda x: x.split(".")[0])
    xpore["kmer"] = xpore.kmer.apply(lambda x: x.replace("T", "U"))
    xpore["middle_base"] = xpore.kmer.apply(lambda x: x[2])
    xpore_counts = xpore.groupby("id").agg(  # gene_name=("gene_name", "first"),
        dmr_count=("kmer", "count"),
        cDNA_end_fiveUTR=("cDNA_end_fiveUTR", "first"),
        cDNA_start_fiveUTR=("cDNA_start_fiveUTR", "first"),
        cDNA_end_threeUTR=("cDNA_end_threeUTR", "first"),
        cDNA_start_threeUTR=("cDNA_start_threeUTR", "first"),
        fiveUTR_len=("five_prime_len", "first"),
        threeUTR_len=("three_prime_len", "first"),
        dmr_A_count=("middle_base", lambda x: (x == "A").sum()),
        dmr_C_count=("middle_base", lambda x: (x == "C").sum()),
        dmr_G_count=("middle_base", lambda x: (x == "G").sum()),
        dmr_U_count=("middle_base", lambda x: (x == "U").sum()),
    ).sort_values(by="dmr_count", ascending=False)

    xpore["modrate_C3"] = xpore[["mod_rate_s3-r1", "mod_rate_s3-r2", "mod_rate_s3-r3"]].mean(axis=1)
    xpore["modrate_C4"] = xpore[["mod_rate_s4-r1", "mod_rate_s4-r2", "mod_rate_s4-r3"]].mean(axis=1)

    modrates_lower = xpore[xpore["mod_assignment"] == "lower"].groupby("id").agg(
        modrate_C3_lower=("modrate_C3", "median"),
        modrate_C4_lower=("modrate_C4", "median"),
 )
    modrates_higher = xpore[xpore["mod_assignment"] == "higher"].groupby("id").agg(
        modrate_C3_higher=("modrate_C3", "median"),
        modrate_C4_higher=("modrate_C4", "median"),
    )

    modrates = xpore.groupby("id").agg({"modrate_C3": "median",
                                        "modrate_C4": "median"})
    moderates_3UTR = xpore[(xpore["relative_cdna_pos"] >= 1.8)].groupby("id").agg({"modrate_C3": "median",
                                                                                   "modrate_C4": "median"
                                                                                   })

    moderates_3UTR_higher = xpore[(xpore["mod_assignment"] == "higher")
                                  & (xpore["relative_cdna_pos"] >= 2)].groupby("id").agg(
        modrate_C3_3utr_higher=("modrate_C3", "median"),
        modrate_C4_3utr_higher=("modrate_C4", "median"),
    )

    moderates_3UTR_lower = xpore[(xpore["mod_assignment"] == "lower")
                                  & (xpore["relative_cdna_pos"] >= 2)].groupby("id").agg(
        modrate_C3_3utr_lower=("modrate_C3", "median"),
        modrate_C4_3utr_lower=("modrate_C4", "median"),
    )

    modrates["xpore_dmr_C4vsC3"] = modrates["modrate_C4"] - modrates["modrate_C3"]
    modrates_lower["xpore_dmr_C4vsC3_lower"] = modrates_lower["modrate_C4_lower"] - modrates_lower["modrate_C3_lower"]
    modrates_higher["xpore_dmr_C4vsC3_higher"] = modrates_higher["modrate_C4_higher"] - modrates_higher["modrate_C3_higher"]
    moderates_3UTR_higher["xpore_dmr_C4vsC3_3utr_higher"] = moderates_3UTR_higher["modrate_C4_3utr_higher"] - moderates_3UTR_higher["modrate_C3_3utr_higher"]
    moderates_3UTR_lower["xpore_dmr_C4vsC3_3utr_lower"] = moderates_3UTR_lower["modrate_C4_3utr_lower"] - moderates_3UTR_lower["modrate_C3_3utr_lower"]

    mrna_level["C3"] = mrna_level[[x for x in mrna_level.columns if "C3" in x]].mean(axis=1)
    mrna_level["C4"] = mrna_level[[x for x in mrna_level.columns if "C4" in x]].mean(axis=1)
    mrna_level["transcript"] = mrna_level["transcript"].apply(lambda x: x.split(".")[0])
    mrna_level = mrna_level[["transcript", "C3", "C4"]]
    mrna_level = mrna_level.rename(columns={"C3": "mrna_illumina_C3", "C4": "mrna_illumina_C4",
                                            "pval": "mrna_illumina_pval",
                                            "qval": "mrna_illumina_qval",
                                            })
    mrna_level["log2_mRNA_C4_vs_C3"] = mrna_level["mrna_illumina_C4"] - mrna_level["mrna_illumina_C3"]

    mrna_level_ont["mrna_ont_C3"] = np.log2(mrna_level_ont[["C3_R1", "C3_R2", "C3_R3"]].mean(axis=1) + 1)
    mrna_level_ont["mrna_ont_C4"] = np.log2(mrna_level_ont[["C4_R1", "C4_R2", "C4_R3"]].mean(axis=1) + 1)
    mrna_level_ont = mrna_level_ont.rename(columns={"HGNC symbol": "gene_name"})

    print(mrna_level_ont[mrna_level_ont["gene_name"] == "ITGB1"])


    mrna_level_ont = mrna_level_ont[["transcript_id", "gene_id", "mrna_ont_C3", "mrna_ont_C4", "gene_name", "transcript_length"]]
    mrna_level = mrna_level.drop_duplicates(subset="transcript", keep="first")

    polya_tail_len["contig"] = polya_tail_len["contig"].apply(lambda x: x.split(".")[0])
    polya_tail_len["polya_length_C3"] = polya_tail_len[["polya_length_s3_r1", "polya_length_s3_r2", "polya_length_s3_r3"]].mean(axis=1)
    polya_tail_len["polya_length_C4"] = polya_tail_len[["polya_length_s4_r1", "polya_length_s4_r2", "polya_length_s4_r3"]].mean(axis=1)
    polya_tail_len = polya_tail_len[["contig", "polya_length_C3", "polya_length_C4"]]
    #############################################
    # set transcrip_ID to index
    #############################################
    mrna_level = mrna_level.set_index("transcript")
    polya_tail_len = polya_tail_len.set_index("contig")
    mrna_level_ont = mrna_level_ont.set_index("transcript_id")
    # modrates = modrates.set_index("id") already set

    so = [xpore_counts,
          modrates,
          modrates_lower,
          modrates_higher,
          moderates_3UTR_higher,
          moderates_3UTR_lower,
          mrna_level,
          mrna_level_ont,
          polya_tail_len
          ]

    for tbl in so:
        print("\n")
        print(tbl.head(2))

    data2 = pd.concat(so, axis=1)
    data2 = data2.dropna(subset=["mrna_ont_C3", "mrna_ont_C4"])
    data2 = data2[data2["dmr_count"] > 0]

    print("SAVE DATA")
    data2 = data2.rename(columns={"pval": "illumina_pval",
                                "qval": "illumina_qval",
                                "log2_mRNA_C4_vs_C3": "mrna_illumina_log2FC_C4vsC3",
                                "modrate_C3": "xpore_median_modrate_C3",
                                "modrate_C4": "xpore_median_modrate_C4",
                                "kmer_count": "dmr_counts"
                                  })
    print(data2.info())
    for base in ["A", "C", "G", "U"]:
        data2[f"dmr_{base}_count"] = data2[f"dmr_{base}_count"].fillna(0).astype(int)
    data2["dmr_count"] = data2["dmr_count"].fillna(0).astype(int)
    data2["transcript_id"] = data2.index
    data2 = data2.dropna(subset=["transcript_id"])

    print(data2.columns)

    old = ['dmr_count', 'cDNA_end_fiveUTR', 'cDNA_start_fiveUTR', 'cDNA_end_threeUTR', 'cDNA_start_threeUTR',
           'dmr_A_count', 'dmr_C_count', 'dmr_G_count', 'dmr_U_count', 'fiveUTR_len', 'threeUTR_len',
           'xpore_median_modrate_C3', 'xpore_median_modrate_C4', 'xpore_dmr_C4vsC3', 'mrna_illumina_C3',
           'mrna_illumina_C4', 'mrna_illumina_log2FC_C4vsC3', 'gene_id', 'mrna_ont_C3',
           'mrna_ont_C4', 'gene_name', 'transcript_length', 'polya_length_C3', 'polya_length_C4', 'transcript_id']
    new = ['transcript_id', 'gene_id', 'gene_name',
           'dmr_count', 'dmr_A_count', 'dmr_C_count', 'dmr_G_count', 'dmr_U_count',
           'xpore_median_modrate_C3', 'xpore_median_modrate_C4', 'xpore_dmr_C4vsC3',
           "xpore_dmr_C4vsC3_lower",
           "xpore_dmr_C4vsC3_higher",
           "xpore_dmr_C4vsC3_3utr_lower",
           "xpore_dmr_C4vsC3_3utr_higher",
           'mrna_illumina_C3','mrna_illumina_C4', 'mrna_illumina_log2FC_C4vsC3', 'mrna_ont_C3',
           'mrna_ont_C4',  'transcript_length', 'fiveUTR_len', 'threeUTR_len',
           'polya_length_C3', 'polya_length_C4']
    data2 = data2.reset_index()
    data2 = data2.drop("index", axis=1)
    data2 = data2[new]
    print(data2)

    data2.to_csv(out_file, sep="\t", index=False)
    data2.to_excel(out_file.replace(".tsv", ".xlsx"),  index=False)
    print("Done")
    # print(data2.xpore_kmer_per_mrna.sum())
    return data2


if __name__ == '__main__':

    # define file_names
    # xpore_file = r"\majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR.tsv"
    xpore_file = r"\majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS.tsv"
    xpore_file = r"direct_RNA_seq\Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4" + xpore_file
    mrna_level_file = r"..\Illumina\hSMC_TGFB_PDGF\DATA\results_C3_C4_2025\tables\logcount-matrix\model_condition.logcount-matrix_postprocessed.tsv"
    mrna_level_ont_file = r"direct_RNA_seq/Nanocount/DESeq2_results_Tobias/ONT_transcript_level_C3_C4_deseq2_norm_tpm.tsv"
    gene_level_ont_file = r"direct_RNA_seq/Nanocount/DESeq2_results_Tobias/ONT_gene_level_C3_C4_deseq2_norm_tpm.tsv"
    degs_ont_deseq2_file = r"direct_RNA_seq/Nanocount/DESeq2_results_Tobias/ONT_gene_level_C3_C4_DEseq2_sum_stats.tsv"
    polya_file = r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore\direct_RNA_seq\PolyA\polya_C4_vs_C3_sumstats.tsv"

    # take DEseq2 normalized counts for Master table
    # degs_file_illumina = r"..\Illumina\hSMC_TGFB_PDGF\DATA\results_C3_C4_2025\tables\diffexp\model_condition.genes-representative.diffexp_postprocessed.tsv"
    overwrite = True
    tag = ""

    name_tag = "Transcript_Level_C3_vs_C4_padj_DMR_025_PAPER_higher_lower.tsv"
    out_file = f"direct_RNA_seq/0_TABLES/DMR_MASTER_TABLE_{tag}{name_tag}.tsv"

    if os.path.exists(out_file) and not overwrite:
        data = pd.read_table(out_file, index_col=0)
    else:
        data = create_xpore_master_table()
    quit()
