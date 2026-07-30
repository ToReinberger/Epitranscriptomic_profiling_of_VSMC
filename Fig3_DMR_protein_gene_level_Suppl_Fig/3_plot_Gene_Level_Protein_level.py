import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy import stats
import numpy as np
from matplotlib import cm


pd.set_option('display.max_rows', 150)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def get_colormap_for_values(values, cmap_name):
    """
    :param values:
    :return colormap:
    """
    colors = []
    rgb = list(cm.get_cmap(name=cmap_name)(np.arange(0, 256)))
    max_val = max(values)
    min_val = min(values)
    if abs(min_val) > abs(max_val):
        max_val = abs(min_val)
    if max_val == 0:
        max_val = 0.0001
    for val in values:
        temp_ = val / max_val
        index = 127 + int(128 * temp_)
        if index == -1:
            index = 0
        colors.append(rgb[index])
    return colors, max_val, rgb


def plot_jointplot(data_,
                   size_col=None,
                   color_col=None,
                   group_col=None,
                   x_label="mRNA length",
                   y_label="Poly(A) length",
                   x_col='cdna_len',
                   y_col='polya_length',
                   annotate=False,
                   plot_line=False,
                   set_ax_lim=False,
                   fig_name="mRNA_len_polyA_tail.svg",
                   size=8
                   ):

    print(data_)
    show_genes = data_.loc[(data_["dmr_count"] >= 5)
                           & ((data["protein_qval"] <= 0.05) |
                              (data["ont_gene_qval"] <= 0.05))
    , "gene_name"].values


    if group_col is None:
        data_["group"] = False
        group_col = "group"
        r, p = stats.pearsonr(data[x_col],
                              data[y_col])
        print("is_DMR", r, p)
    else:
        pass
        r, p = stats.pearsonr(data.loc[data[group_col], x_col],
                              data.loc[data[group_col], y_col])
        print("is_DMR", r, p)

        r, p = stats.pearsonr(data.loc[~data[group_col], x_col],
                              data.loc[~data[group_col], y_col])
        print("non_DMR", r, p)

    # define scatter size
    if size_col is None:
        data_["ms"] = 20
    else:
        data_[size_col + "_temp"] = (data_[size_col] - data_[size_col].min()) / (
                data_[size_col].max() - data_[size_col].min())
        data_[size_col + "_temp"] = (data_[size_col + "_temp"] + 1) * 1.7
        # data_[size_col + "_temp"] = data_[size_col] / data_[size_col].max() * 20
        # data[size_col + "_temp"] = data[size_col] / 10
        data_["ms"] = [x ** 4.2 for x in data_[size_col + "_temp"].to_list()]
        data_.loc[data_.is_DMR == False, "ms"] = 20

    # define color gradient
    if color_col is None:
        data_.loc[data_.is_DMR == False, "colors"] = "#808080"
        data_.loc[data_.is_DMR, "colors"] = "firebrick"

        # data_.loc[data_.is_DMR == "pos", "colors"] = "firebrick"
        # data_.loc[data_.is_DMR == "neg", "colors"] = "steelblue"
        print(data_[data_["gene_name"].isin(["FN1"])])
        # quit()
    else:
        data_[color_col] = (data_[color_col] - data_[color_col].min()) / (
                    data_[color_col].max() - data_[color_col].min())
        color_vals = data_[color_col].to_list()
        # max_val = max(color_vals)
        # color_vals = [x / max_val for x in color_vals]
        colors, _, cmap = get_colormap_for_values(color_vals, "Reds")
        # mask = np.array([np.array_equal(c, cmap[0]) for c in colors])
        # gray = [x/255 for x in [80, 80, 80, 255]]
        #colors[mask] = gray
        data_["colors"] = colors
        data_.loc[~data_.is_DMR.isin(["neg", "pos"]), "colors"] = "gray"

    data_.rename(columns={x_col: x_label, y_col: y_label}, inplace=True)
    plt.style.use("ggplot")

    if set_ax_lim:
        g = sns.JointGrid(data=data_, x=x_label, y=y_label, height=size,
                          xlim=(-3.1, 3.1), ylim=(-3.1, 3.1)
                          )
    else:
        g = sns.JointGrid(data=data_, x=x_label, y=y_label, height=size,
                          # xlim=(-3.3, 3.3), ylim=(-3.3, 3.3)
                          )
    # Plot for each group
    for group, sub_data in data_.groupby(group_col):
        # Joint scatter and regression
        print(group, len(sub_data))
        s = sub_data.loc[sub_data[group_col] == group, "ms"].to_list()
        c = sub_data.loc[sub_data[group_col] == group, "colors"].to_list()

        # print(len(s), len(c), len(sub_data))
        if group:
            color = "firebrick"
        else:
            color = "gray"

        if group:
            rasterized = False
        else:
            rasterized = True

        sns.regplot(
            data=sub_data,
            truncate=False,
            x=x_label,
            y=y_label,
            ax=g.ax_joint,

            scatter_kws={'label': group,
                         'color': c,
                         "s": s,
                         "edgecolors": "white",
                         "linewidths": 0.4,
                         "rasterized": rasterized
                         },

            line_kws={'color': color},
            fit_reg=True

        )

        if plot_line:
            g.ax_joint.axvline(x=1, linestyle="dotted", color="#555555", alpha=0.6, zorder=5)
            g.ax_joint.axvline(x=-1, linestyle="dotted", color="#555555", alpha=0.6, zorder=5)

            g.ax_joint.axhline(y=-np.log2(1.5), linestyle="dotted", color="#555555")
            g.ax_joint.axhline(y=np.log2(1.5), linestyle="dotted", color="#555555")

        # Marginal histograms
        sns.kdeplot(sub_data[x_label], ax=g.ax_marg_x, color=color, alpha=0.4, multiple="stack", cut=0)
        sns.kdeplot(sub_data[y_label], ax=g.ax_marg_y, color=color, alpha=0.4, multiple="stack", vertical=True, cut=0)

    if annotate:
        # print(data_)
        data_["gene_name"].fillna(" ")
        print(data["dmr_count"].max())
        genes = []
        for p, m, k, gene in zip(data_[y_label],
                              data_[x_label],
                              data_["dmr_count"],
                              data_["gene_name"],
                            ):
            # print(p, m, k, gene)
            if type(gene) is not str:
                continue
            if gene not in show_genes:
                continue
            if "COL3A1" in gene:
                print("COL3A1")
            # if (((abs(p) > np.log2(1.5) or abs(m) > 1.0) and k > 0)
            #         or ((abs(p) > 1.5 and abs(m) > 3) and k > 1)
            #         or ((abs(p) > 2 or abs(m) > 3) and k > 2)
            #         # or ("COL3A1" in gene)
            #
            # ):
            if gene in genes or gene == np.nan:
                pass
            genes.append(gene)
            if m < 0:
                ha = "right"
            else:
                ha = "left"
            if abs(p) > 1.5 and abs(m) > 3 and k > 1:
                color = "gray"
            else:
                color = "k"
            gene = gene.split("-")[0] # + " | " + kmer
            g.ax_joint.text(x=m + np.sign(m) * 0.03, y=p + .03, s=gene,
                            fontsize=7, fontweight="bold", ha=ha, color="k",
                            bbox=dict(
                                facecolor="white",
                                alpha=0.4,
                                edgecolor="none",
                                pad=1.5
                            )
                            )
            # Increase fontsize of xticks and yticks
            g.ax_joint.tick_params(axis='x', labelsize=12)
            g.ax_joint.tick_params(axis='y', labelsize=12)
    # Add legend
    # g.ax_joint.legend(title='Group')

    plt.tight_layout()
    plt.savefig(fig_name)
    # plt.close()
    plt.show()


def plot_mrna_with_and_without_protein():
    rna_protein = data[(data[f"gene_level_ont_{cond}"] > 0) & (data[f"protein_level_{cond}"] > 0)]
    count_rna_protein_dmr = len(rna_protein[rna_protein["is_DMR"]])
    count_rna_protein = len(rna_protein[~rna_protein["is_DMR"]])
    rna_protein_ratio = count_rna_protein_dmr / (count_rna_protein_dmr + count_rna_protein) * 100

    no_protein = data[(data[f"gene_level_ont_{cond}"] > 0) & (data[f"protein_level_{cond}"] == 0)]
    print(no_protein)
    count_no_protein_dmr = len(no_protein[no_protein["is_DMR"]])
    count_no_protein = len(no_protein[~no_protein["is_DMR"]])
    no_protein_ratio = count_no_protein_dmr / (count_no_protein_dmr + count_no_protein) * 100
    print(count_no_protein_dmr, count_no_protein)
    no_rna = data[(data[f"gene_level_ont_{cond}"] == 0) & (data[f"protein_level_{cond}"] > 0)]
    count_no_rna_dmr = len(no_rna[no_rna["is_DMR"]])
    count_no_rna = len(no_rna[~no_rna["is_DMR"]])
    no_rna_ratio = count_no_rna_dmr / (count_no_rna_dmr + count_no_rna) * 100
    plt.style.use('ggplot')
    df = pd.DataFrame(data={
        "mRNA\nwith protein": [count_rna_protein_dmr, count_rna_protein],
        "mRNA\nwithout protein": [count_no_protein_dmr, count_no_protein],
        "Protein\nwithout mRNA": [count_no_rna_dmr, count_no_rna]
                            }, index=["DMR+", "DMR-"])

    df.T.plot(kind="bar", stacked=True, figsize=(6, 3), color=["firebrick", "gray"], alpha=0.7)
    plt.grid(False)
    plt.xticks(rotation=0)
    plt.ylabel("Counts")
    plt.text(x=0, y=count_rna_protein_dmr + 100, s=f"{round(rna_protein_ratio, 1)}%", ha="center", va="bottom", color="#555555")
    plt.text(x=1, y=count_no_protein_dmr + 100, s=f"{round(no_protein_ratio, 1)}%", ha="center", va="bottom", color="#555555")
    plt.text(x=2, y=count_no_rna_dmr + 100, s=f"{round(no_rna_ratio, 1)}%", ha="center", va="bottom", color="#555555")

    print(count_no_rna_dmr, count_no_rna)
    plt.tight_layout()
    plt.savefig(f"ONT_Gene_level_vs_Protein_level_mRNA_without_protein_{cond}.svg")
    plt.show()
    quit()


if __name__ == '__main__':

    cond = "C4"
    data = pd.read_table('direct_RNA_seq/0_TABLES/DMR_MASTER_TABLE_Gene_Level_C3_vs_C4_padj_DMR_025_PAPER.tsv')
    dmr_count_thresh = 3
    data["is_DMR"] = False
    data.loc[data["dmr_count"] >= dmr_count_thresh, "is_DMR"] = True

    print(data)



    # xpore = pd.read_table(majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS_GENENAMES.tsv
    trancr_data = pd.read_table("direct_RNA_seq/0_TABLES/DMR_MASTER_TABLE_Transcript_Level_C3_vs_C4_padj_DMR_025_PAPER_higher_lower.tsv.tsv")

    log_count_thresh = 2
    xpore_mode = ""

    trancr_data[f"xpore_dmr_C4vsC3{xpore_mode}"] = trancr_data[f"xpore_dmr_C4vsC3{xpore_mode}"].fillna(0)
    trancr_data_neg = trancr_data[(trancr_data[f"xpore_dmr_C4vsC3{xpore_mode}"] <= -0.25)
                                  & ((trancr_data["mrna_ont_C3"] > log_count_thresh)
                                     | (trancr_data["mrna_ont_C4"] > log_count_thresh))
                              & (trancr_data["dmr_count"] >= dmr_count_thresh) & (trancr_data["dmr_U_count"] >= 0)
                              ]

    trancr_data_pos = trancr_data[(trancr_data[f"xpore_dmr_C4vsC3{xpore_mode}"] >= 0.25)
                                  & ((trancr_data["mrna_ont_C3"] > log_count_thresh)
                                     | (trancr_data["mrna_ont_C4"] > log_count_thresh))
                              & (trancr_data["dmr_count"] >= dmr_count_thresh) & (trancr_data["dmr_U_count"] >= 0)
                              ]

    # trancr_data = trancr_data[(trancr_data["xpore_dmr_C4vsC3"] < -0.2) & ((trancr_data["mrna_ont_C3"] > 3)
    #                                                                    | (trancr_data["mrna_ont_C4"] > 3))
    #                           & (trancr_data["dmr_count"] >= dmr_count_thresh) & (trancr_data["dmr_U_count"] >= 1)
    # ]


    print(trancr_data[trancr_data["gene_name"].isin(["FN1", "THBS2"])])
    #
    # data.loc[data["gene_name"].isin(trancr_data_pos["gene_name"]), "is_DMR"] = "pos"
    # data.loc[data["gene_name"].isin(trancr_data_neg["gene_name"]), "is_DMR"] = "neg"

    print(len(data[data["is_DMR"] == "pos"]))
    print(len(data[data["is_DMR"] == "neg"]))

    #  data["DMR_direction"] = "none"
    # data.loc[data["gene_name"].isin(trancr_data_neg["gene_name"]), "DMR_direction"] = "neg"
    # data_temp_neg = data[data.is_DMR_neg]
    # data["is_DMR_pos"] = False
    # data.loc[(data["gene_name"].isin(trancr_data_pos["gene_name"]))
    #                                 & (~data["gene_name"].isin(trancr_data_neg["gene_name"])),
    #                                 "DMR_direction"] = "pos"
    # data_temp_pos = data[data.is_DMR_pos]

    print(data[data["gene_name"].isin(["FN1", "THBS2"])])

    #
    # sns.kdeplot(data=data[data["DMR_direction"] != "none"],
    #             hue="DMR_direction",
    #             y="protein_log2FC_C4vsC3",
    #             cut=0)
    # plt.show()
    # quit()
    # print(len(data))
    #
    # print(len(data[(data[f"ont_{cond}_gene_tpm"] > 0)
    #                & (data[f"protein_level_{cond}"] == 0)
    #                ]))
    #
    # print(len(data[(data[f"ont_{cond}_gene_tpm"] == 0)
    #                & (data[f"protein_level_{cond}"] > 0)
    #                ]))
    #
    #
    # # print(data[(data[f"ont_{cond}_gene_tpm"] == 0)
    # #                & (data[f"protein_level_{cond}"] > 0)
    # #                ])
    #
    # print(len(data[(data[f"ont_{cond}_gene_tpm"] > 0)
    #                & (data[f"protein_level_{cond}"] == 0) & data.is_DMR
    #                ]))
    #
    #
    # # print(data[(data[f"ont_{cond}_gene_tpm"] > 0)
    # #                & (data[f"protein_level_{cond}"] == 0) & data.is_DMR
    # #                ])
    #
    # print(len(data[(data[f"ont_{cond}_gene_tpm"] == 0)
    #                & (data[f"protein_level_{cond}"] > 0) & data.is_DMR
    #                ]))
    #
    # print(len(data[(data[f"ont_{cond}_gene_tpm"] > 0)
    #                & (data[f"protein_level_{cond}"] > 0) & data.is_DMR
    #                ]))
    #
    # print(len(data[(data[f"ont_{cond}_gene_tpm"] > 0)
    #                & (data[f"protein_level_{cond}"] > 0) & ~data.is_DMR
    #                ]))
    data = data[
        (data[f"gene_level_ont_{cond}"] > 0)
        & (data[f"protein_level_{cond}"] > 0)
        ]

    data["ont_gene_log2FC_C4vsC3"] = data["ont_gene_log2FC_C4vsC3"].fillna(0)
    data["protein_log2FC_C4vsC3"] = data["protein_log2FC_C4vsC3"].fillna(0)
    print(data.head(30))
    # for g in data.loc[data.xpore_kmer_per_mrna >= 5, "gene_name"]:
    #     print(g)

    data = data[
        (data["protein_qval"] <= 0.05)
        # |
        # (data["ont_gene_qval"] <= 0.05)
    ]

    # z-score
    # col = f'gene_level_ont_{cond}'
    # data[col] = (data[col] - data[col].mean()) / data[col].std()
    # col = f'protein_level_{cond}'
    # data[col] = (data[col] - data[col].mean()) / data[col].std()

    print(data[(data.protein_log2FC_C4vsC3 > 0) & data.is_DMR])
    print("\n")
    print(data[(data.protein_log2FC_C4vsC3 < 0) & data.is_DMR])


    # plot_mrna_with_and_without_protein()

    # plot_jointplot(data,
    #                color_col=None,
    #                size_col=None,
    #                group_col="is_DMR",
    #                x_label=f"{cond} - gene level| log2(TPM +1)",
    #                y_label=f"{cond} - protein level | log2(maxLFQ +1)",
    #                x_col=f'gene_level_ont_{cond}',
    #                y_col=f'protein_level_{cond}',
    #                annotate=True,
    #                set_ax_lim=False,
    #                size=4,
    #                fig_name=f"ONT_Gene_level_vs_Protein_level_{cond}.svg")


    # plot_jointplot(data,
    #                color_col="dmr_count",
    #                size_col="dmr_count",
    #                group_col="is_DMR",
    #                x_label=f"{cond} - gene level| log2(TPM +1)",
    #                y_label=f"{cond} - protein level | log2(maxLFQ +1)",
    #                x_col=f'gene_level_ont_{cond}',
    #                y_col=f'protein_level_{cond}',
    #                annotate=True,
    #                set_ax_lim=False,
    #                size=6,
    #                fig_name=f"ONT_Gene_level_vs_Protein_level_{cond}.svg")


    tag = "sign_protein"

    plot_jointplot(data,
                   color_col=None,
                   size_col="dmr_count",
                   group_col="is_DMR",
                   x_label="Gene level\nlog2(fold change)",
                   y_label="Protein level\nlog2(fold change)",
                   x_col=f'ont_gene_log2FC_C4vsC3',
                   y_col=f'protein_log2FC_C4vsC3',
                   annotate=True,
                   plot_line=True,
                   set_ax_lim=True,
                   fig_name=f"ONT_Gene_level_vs_Protein_level_C4_vs_C3_log2FC_{tag}_dmr_count{dmr_count_thresh}.pdf")

    # quit()
    print(data.head(5))


