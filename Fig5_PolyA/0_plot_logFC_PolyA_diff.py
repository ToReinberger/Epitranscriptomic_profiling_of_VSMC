
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def check_protein_polya_diff_corr():
    # print(polya.head(5))
    polya = pd.merge(polya_in, master, left_on="genes", right_on="gene_name", how="right")
    print(polya)
    print(len(polya[polya["dmr_count"] >= num_dmr]))

    plt.style.use('ggplot')
    # plt.figure(figsize=(4, 4))

    if num_dmr == 0:
        polya = polya[
            (polya["dmr_count"] == 0)
            &
            (polya["protein_qval"] <= 0.05)
             # |
             # (polya["ont_gene_qval"] <= 0.05))
            ]
    else:
        polya = polya[
            (polya["dmr_count"] >= num_dmr)
            &
            ((polya["protein_qval"] <= 0.05)
            |
            (polya["ont_gene_qval"] <= 0.05))
        ]

    print(num_genes := len(polya))
    polya = polya.dropna(subset=["diff", "ont_gene_log2FC_C4vsC3", "protein_log2FC_C4vsC3"])
    r_rna, p_rna = pearsonr(polya["diff"], polya["ont_gene_log2FC_C4vsC3"])
    r_protein, p_protein = pearsonr(polya["diff"], polya["protein_log2FC_C4vsC3"])
    print("rna protein")
    print(r_rna, p_rna,
          r_protein, p_protein)

    plt.figure(figsize=(8, 5))
    plt.subplot(1, 2, 1)
    sns.regplot(data=polya, x="diff", y="ont_gene_log2FC_C4vsC3",
                label=f"Gene level\nr={round(r_rna, 2)}, p={p_rna:.1e}",
                color="#909090",
                scatter_kws={"alpha": 0.8,
                             },
                marker="D",
                line_kws={"color": "#808080", "alpha": 0.6},
                fit_reg=True)

    plt.xlabel("Δ z-score oftranscript-weighted\nmean poly(A) length\n(C4 − C3)")
    plt.ylabel(f"Log$_2$(fold change)\n(C4 vs. C3")
    plt.title("Gene level")
    plt.ylim(-2.6, 1.2)
    # plt.xlim(-2.0, 2.0)
    # plt.xlim(-70, 30)
    # plt.yticks(np.arange(-2, 3, 1), np.arange(-2, 3, 1))
    if show_names:
        for gene_name, x, y in zip(polya["gene_name"], polya["diff"], polya["ont_gene_log2FC_C4vsC3"]):
            plt.annotate(gene_name, xy=(x, y - 0.12),
                         xytext=(x, y - 0.12),
                         # textcoords='offset points',
                         ha='center', va='center',
                         fontsize=7, color="k", zorder=10000)

    plt.legend(loc=(0.05, 0.05), fontsize=9)
    if num_dmr == 0:
        plt.suptitle(f"{num_genes} genes without DMRs")
    else:
        plt.suptitle(f"{num_genes} genes with ≥{num_dmr} DMRs")
    plt.subplot(1, 2, 2)

    sns.regplot(data=polya, x="diff", y="protein_log2FC_C4vsC3",
                label=f"Protein level\nr={round(r_protein, 2)}, p={p_protein:.1e}",
                color="#909090",
                scatter_kws={"alpha": 0.8},
                line_kws={"color": "#808080", "alpha": 0.6},
                fit_reg=True)


    plt.ylim(-2.6, 1.2)
    # plt.xlim(-2.0, 2.0)

    plt.xlabel("Δ z-score of transcript-weighted\nmean poly(A) length\n(C4 − C3)")
    plt.ylabel(f"Log$_2$(fold change)\n(C4 vs. C3")
    plt.title("Protein level")
    if show_names:
        for gene_name, x, y in zip(polya["gene_name"], polya["diff"], polya["protein_log2FC_C4vsC3"]):
            plt.annotate(gene_name, xy=(x, y - 0.12),
                         xytext=(x, y - 0.12),
                         # textcoords='offset points',
                         ha='center', va='center',
                         fontsize=7, color="k", zorder=10000)

    # for y_p, y_g, x in zip(polya["protein_log2FC_C4vsC3"], polya["ont_gene_log2FC_C4vsC3"], polya["diff"]):
    #     dy = y_p - y_g
    #     plt.arrow(x, y_g - 0.1, 0,dy, color="k", width=0.01, head_width=0.1)

    plt.legend(loc=(0.05, 0.05), fontsize=9)


    plt.tight_layout()

    plt.savefig(f"polyA_z_Score_diff_C4_C3_vs_logFC_C4_C3_num_dmr_{num_dmr}_all_polya_data.svg")
    plt.show()
    return


if __name__ == '__main__':
    num_dmr = 1
    show_names = False
    polya_in = pd.read_table("polyA_gene_level_sumstats_NEW.tsv")
    # polya_in = polya_in[abs(polya_in["diff"]) >= 0.5]
    master = pd.read_table("../0_TABLES/DMR_MASTER_TABLE_Gene_Level_C3_vs_C4_padj_DMR_025_PAPER.tsv")
    master = master.sort_values(by="protein_log2FC_C4vsC3")
    sign_genes = polya_in.loc[polya_in["pval"] <= 0.05, "genes"].values
    # master = master[master["gene_name"].isin(sign_genes)]
    check_protein_polya_diff_corr()