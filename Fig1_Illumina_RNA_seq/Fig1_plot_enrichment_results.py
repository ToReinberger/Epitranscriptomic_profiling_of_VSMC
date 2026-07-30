import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import glob
import os

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 80)
pd.set_option('display.width', 1200)
pd.set_option('max_colwidth', 200)


def clean_terms(term):
    if "UV Response" in term:
        term = term.replace("Dn", "down")
    if " (GO" in term:
        term = term[:term.find(" (GO")]
    elif " R-HSA-" in term:
        term = term[:term.find(" R-HSA")]
    terms = term.split(" ")
    terms_temp = [terms[0]]
    for term in terms[1:]:
        if "NF-" in term or "A/1" in term or "I-" in term:
            terms_temp.append(term)
        else:
            terms_temp.append(term.lower())
    return " ".join(terms_temp)


def load_data(terms, color_map, overwrite=True):
    out = []
    out_temp = []
    if os.path.exists(out_file) and not overwrite:
        df = pd.read_excel(out_file, index_col=0)
        # df = df[df.Exclude != "x"]
        df = df.set_index(np.arange(0, len(df)))
        print(df)
    else:

        for file in glob.glob(folder):
            print(file)
            df = pd.read_csv(file, sep="\t")
            df["cluster"] = file[-6:-4]
            df = df[(df.num_genes >= 8)  # & (df["Odds ratio"] > 3)
            ]

            for pos_reg in df.loc[df["Term name"].str.contains("Positive regulation"), "Term name"]:
                df = df[df["Term name"] != pos_reg.replace("Positive r", "R")]
            for neg_reg in df.loc[df["Term name"].str.contains("Negative regulation"), "Term name"]:
                df = df[df["Term name"] != neg_reg.replace("Negative r", "R")]

            df["color"] = df["cluster"].map(color_map)
            df = df[df.library.isin(terms)]
            df = df.drop_duplicates(subset="Overlapping genes")

            out_temp.append(df)
            df = df.nsmallest(4, "Rank")

            out.append(df)
            """ for term, tbl in df.groupby("library"):
                if term not in terms:
                    continue
                # tbl = tbl[tbl["Adjusted p-value"] < 0.01]
                tbl = tbl.drop_duplicates(subset="Overlapping proteins")
                # tbl = tbl.sort_values(by="Odds ratio", ascending=False)

                # tbl = tbl.head(3)
                out.append(tbl.nsmallest(4, "Rank"))"""
        df = pd.concat(out)
        print(df)
        df_temp = pd.concat(out_temp)
        df = df.set_index(np.arange(0, len(df)))
        df_temp = df_temp.set_index(np.arange(0, len(df_temp)))
        df["Term name"] = df["Term name"].apply(clean_terms)
        df = df.drop_duplicates(subset="Term name")

        df_temp["Term name"] = df_temp["Term name"].apply(clean_terms)
        df_temp = df_temp.drop_duplicates(subset="Term name")
        df.to_excel(out_file)
        df_temp.to_excel(out_file)
        print("SELECT TERMS and RESTART")
        return df


def plot_enrichment_degs():
    plt.style.use("ggplot")
    cmap = sns.color_palette("colorblind").as_hex()[::-1]
    color_map = dict(zip(["C1",
                          "C2", "C3", "C4", "C5",
                          "C6",

                          ], cmap))
    fig2 = plt.figure(2)
    ms = 50
    legend_elements = []
    for k, c in color_map.items():
        legend_elements.append(plt.scatter([0], [0], marker='s', color=c, label=k[1], s=ms))
    plt.close(fig2)

    terms = ["MSigDB Hallmark 2020",
             "KEGG 2021 Human",
             "Reactome 2022",
             "GO Biological Process 2023",
             # "DSigDB",
             # "GWAS_Catalog_2023",
             # "Human_Gene_Atlas",
             # "GTEx_Tissues_V8_2023",
             ]
    libs = dict(zip([x.replace(" ", "_") for x in terms], ["MSig", "K", "R", "GO"]))
    use_all = True
    if use_all:
        terms = [x.replace(" ", "_") for x in terms]
        df = load_data(terms=terms, color_map=color_map, overwrite=True)
        cluster_name = ["C1", "C2", "C3", "C4", "C5", "C6",  # "C7"
                        ]
        df = df[~df["Term name"].str.contains("Positive regulation")]
        df["cluster_order"] = df["cluster"].map(dict(zip(cluster_name, np.arange(0, len(cluster_name))[::-1])))
        # df = df.sort_values(by=["cluster", "Odds ratio"], ascending=False)
        df = df.sort_values(by=["cluster_order", "Odds ratio"], ascending=True)
        df = df.set_index(np.arange(0, len(df)))
        print(df)

        plt.title("Enrichment")
        terms_temp = df["Term name"].to_list()

        add_break = False

        if add_break:
            terms = []
            for term in terms_temp:
                term = term.split(" ")
                center = int(len(term) / 2)
                if len(term) > 3:
                    terms.append(" ".join(term[:center]) + "\n" + " ".join(term[center:]))
                else:
                    terms.append(" ".join(term))
        else:
            terms = terms_temp
        # plt.yticks(df.index, df["Term name"],)
        # for x, y, s in zip(df["Odds ratio"], df.index, df.num_genes):
        #     plt.text(x - 0.4, y, s, fontsize=7.5, ha="left", va="center")
        # plt.gca().invert_yaxis()
        plt.xlabel("Odds ratio", fontsize=11)
        plt.ylim(-1, len(df))
        # plt.xlim(34, 0)
        plt.grid(False)
        # plt.gca().yaxis.tick_right()
        plt.barh(df.index, df["Odds ratio"], color=df.color, zorder=1, alpha=1)
        # plt.gca().tick_params(axis="y", direction="in", pad=-10, zorder=2)
        plt.yticks(np.arange(0, len(df)), terms, fontsize=10, ha="left", zorder=2)
        plt.yticks([])
        for y, name, num_genes, lib in zip(np.arange(0, len(df)), terms, df.num_genes, df.library):
            plt.text(0.3, y, name + f" ({num_genes}, {libs[lib]})", fontsize=9, ha="left", va="center")

        plt.legend(handles=legend_elements,
                   fontsize=9,
                   title='Gene\ncluster',
                   alignment="left",
                   title_fontsize=10,
                   frameon=False,
                   facecolor='white',
                   framealpha=1,
                   labelspacing=0.25,
                   handletextpad=0.15,
                   # bbox_to_anchor=(-2.8, 1.4),  # (left, bottom, width, height)
                   loc='best'
                   )
        plt.show()

    else:
        col_idx = 0
        for term in terms:
            term = term.replace(" ", "_")
            print(color_map)
            df = load_data(terms=[term], color_map=color_map, overwrite=True)
            overwrite = True
                # quit()
            cluster_name = ["C1", "C2", "C3", "C4", "C5",  "C6", # "C7"
                        ]
            df = df[~df["Term name"].str.contains("Positive regulation")]
            df["cluster_order"] = df["cluster"].map(dict(zip(cluster_name, np.arange(0, len(cluster_name))[::-1])))
            # df = df.sort_values(by=["cluster", "Odds ratio"], ascending=False)
            df = df.sort_values(by=["cluster_order", "Odds ratio"], ascending=True)
            df = df.set_index(np.arange(0, len(df)))
            print(df)
            # quit()

            """sns.barplot(x=df["Odds ratio"], y=df["Term name"], # tbl=df,
                            # hue="type"
                            )"""
            # df.loc[df.cluster == "C6", "color"] = '#ca9161'

            color_temp = ["red",  # C1
                          "lightgrey",  # C5
                          "green",  # C3
                          "steelblue",  # C4
                          "orange",  # C2
                          ]

            f = plt.figure(1, figsize=(14, 11))

            plt.subplots_adjust(left=0.1, bottom=0.15, top=0.95, right=0.42, wspace=0.1)
            plt.subplot2grid((1, 4), (0, col_idx))
            plot_expr = False
            if plot_expr:
                log_counts = pd.read_table(
                    r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Illumina\hSMC_TGFB_PDGF\DATA\results_conditions_NEW_2025\model_condition.logcount-matrix_temp.tsv",
                    index_col=0)
                conditions = ["C1", "C2", "C3", "C4", "C5"]
                for cond in conditions:
                    log_counts[cond] = log_counts[[f"{cond}_R1-1",
                                                   f"{cond}_R2-1",
                                                   f"{cond}_R3-1",
                                                   f"{cond}_R4-1",
                                                   f"{cond}_R5-1",  # carotids
                                                   f"{cond}_R6-1",
                                                   f"{cond}_R7-1",  # carotids
                                                   f"{cond}_R8-1"
                                                   ]].apply(lambda x: abs(np.mean(x)), axis=1)

                log_counts = log_counts[conditions]

                plt.subplot2grid((1, 5), (0, 0))
                # ax1 = f.add_subplot(121)
                plt.title("Expression")
                y_pos_ = 0
                conditions2 = [0, 4, 2, 3, 1]
                for cluster, gene_list in zip(df["cluster"], df["Overlapping genes"]):
                    expr_vals = []
                    gene_list = gene_list[2:-2].split("', '")
                    means = log_counts.loc[gene_list, :].median()
                    for cond in conditions2:
                        expr_vals.append(means[cond])
                    expr_vals = [2**x for x in expr_vals]
                    expr_vals = expr_vals / max(expr_vals)
                    plt.scatter(x=np.arange(0, 5), y=[y_pos_] * 5, s=expr_vals * 100, color=color_map[cluster])
                    y_pos_ += 1
                plt.ylim(-1, len(df))
                plt.xlim(-1, 5)
                plt.yticks([])
                for idx, label in enumerate(["+++", "---", "+--", "-+-", "++-"]):
                    plt.text(x=idx, y=-2.6, s=label, rotation=90, ha="center", va="bottom",
                           fontdict={"weight": "bold", "size": 12, "font": "Courier New"},
                             bbox=dict(facecolor=color_temp[idx], edgecolor='k', boxstyle='round,pad=0.2',
                                       alpha=0.7)
                             )

                """plt.xticks(np.arange(0, 5), ["+++", "---", "+--", "-+-", "++-"], rotation=90, ha="center",
                           fontdict={"weight": "bold", "size": 18, "font": "Courier New"}, color=color_temp)"""
                plt.xticks([])
                col_idx += 1
                # label = ax1.xaxis.get_ticklabels()
                # label.set_bbox(dict(facecolor=color_temp, edgecolor='none'))
            else:
                plt.subplot2grid((1, 4), (0, col_idx))

                # ax = f.add_subplot(122)
                plt.title(term)
                plt.barh(df.index, df["Odds ratio"], color=df.color)

                terms_temp = df["Term name"].to_list()
                terms = []
                for term in terms_temp:
                    term = term.split(" ")
                    center = int(len(term) / 2)
                    if len(term) > 3:
                        terms.append(" ".join(term[:center]) + "\n" + " ".join(term[center:]))
                    else:
                        terms.append(" ".join(term))

                plt.yticks(np.arange(0, len(df)), terms, fontsize=10)
                # plt.yticks(df.index, df["Term name"],)
                for x, y, s in zip(df["Odds ratio"], df.index, df.num_genes):
                    plt.text(x - 0.4, y, s, fontsize=7.5, ha="left", va="center")
                # plt.gca().invert_yaxis()
                plt.xlabel("Odds ratio", fontsize=11)
                plt.ylim(-1, len(df))
                plt.xlim(34, 0)
                plt.grid(False)
                plt.gca().yaxis.tick_right()

                plt.legend(handles=legend_elements,
                                            fontsize=9,
                                            title='Gene\ncluster',
                                            alignment="left",
                                            title_fontsize=10,
                                            frameon=False,
                                            facecolor='white',
                                            framealpha=1,
                                            labelspacing=0.25,
                                            handletextpad=0.15,
                                            # bbox_to_anchor=(-2.8, 1.4),  # (left, bottom, width, height)
                                            loc='best'
                                            )
                # ax.legend(custom_legend2)
                col_idx += 1
    plt.show()


def plot_enrichment_degs_up_down_between_groups():
    terms = ["MSigDB Hallmark 2020",
             "Reactome 2022",
             "KEGG 2021 Human",
             "GO Biological Process 2023",
             # "DSigDB",
            # "GWAS_Catalog_2023",
             # "Human_Gene_Atlas",
             # "GTEx_Tissues_V8_2023",
             # "GO_Cellular_Component_2023",
             "GO_Molecular_Function_2023"
             ]
    terms = [x.replace(" ", "_") for x in terms]

    out = []
    folder = rf"PLOTS/report_vsmc_covariates_2/*_thresh0.5.tsv"
    # folder = rf"PLOTS/report_vsmc_covariates_2/GeneFlow__C5_ignore_C3_ignore_C2_ignore_C1_up_thresh0.5_C1_C2_up.tsv"
    out_file = fr"PLOTS/report_vsmc_covariates_2/Enrichment_filtered_V2.xlsx"
    overwrite = True
    if os.path.exists(out_file) and not overwrite:
        df = pd.read_excel(out_file, index_col=0)
        if "Exclude" in df.columns:
            df = df[df.Exclude != "x"]
        df = df.set_index(np.arange(0, len(df)))
    else:
        for file in glob.glob(folder):
            dataset = file.split("\\")[-1][:-4].replace("GeneFlow_", "")
            df = pd.read_csv(file, sep="\t")
            print(df)
            df["cluster"] = dataset
            mode = dataset
            df = df[(df.num_genes >= 7)]
            # df["color"] = df["cluster"].apply(lambda x: "#93000C" if x.split("_")[-1] == "up" else "#1E309F")
            df = df[df.library.isin(terms)]
            df = df.sort_values(by=["Combined score", "Odds ratio"], ascending=False)
            df = df.sort_values(by=["Adjusted p-value"], ascending=True)
            out.append(df.head(6))
            """for term, tbl in df.groupby("library"):
                if term not in terms:
                    continue
                # tbl = tbl[tbl["Adjusted p-value"] < 0.01]
                print(tbl)
                tbl = tbl.sort_values(by=["Combined score", "Odds ratio"], ascending=False)
                tbl = tbl.drop_duplicates(subset="Overlapping proteins")
                tbl = tbl.head(3)
                out.append(tbl)"""
        df = pd.concat(out)
        print(df)
        df["Term name"] = df["Term name"].apply(clean_terms)
        # df = df.drop_duplicates(subset="Term name")
        df = df.set_index(np.arange(0, len(df)))
        df.to_excel(out_file)
        print("SELECT TERMS and RESTART")
        print(df)

    f = plt.figure(figsize=(10, 11))
    df.loc[df.cluster.str.contains("down"), "Odds ratio"] = df.loc[df.cluster.str.contains("down"), "Odds ratio"] * -1
    plt.subplots_adjust(left=0.1, bottom=0.1, top=0.95, right=0.4, wspace=0.1)
    row = 1
    conds = sorted(df.cluster.unique(), reverse=True)
    for c in conds:
        if "up" in c:
            color = "#93000C"
        else:
            color = "#1E309F"
        tbl = df[df.cluster == c]
        tbl = tbl.set_index((np.arange(len(tbl))))
        num_genes = df.loc[df.cluster == c, "num_genes"]
        ax = f.add_subplot(6, 1, row)
        plt.barh(tbl["Term name"], tbl["Odds ratio"], color=color, )
        for x, y, s in zip(tbl["Odds ratio"], tbl.index, tbl.num_genes):
            if "up" in c:
                plt.text(x - 0.6, y, s, fontsize=5.5, ha="right", va="center", color="white",
                         fontdict={"weight": "bold", "font": "DejaVu Sans"}
                         )
            else:
                plt.text(x + 0.6, y, s, fontsize=5.5, ha="left", va="center", color="white",
                         fontdict={"weight": "bold", "font": "DejaVu Sans"}
                         )

        if row == 1:
            plt.text(x=-10, y=len(tbl) + 0.1, s="Down", va="center", ha="center", color="#1E309F",
                     fontdict={"weight": "bold"})
            plt.text(x=10, y=len(tbl) + 0.1, s="Up", va="center", ha="center",  color="#93000C",
                     fontdict={"weight": "bold"})
        if row == 6:
            plt.xticks([-10, -20, 0, 10 , 20], [10, 20, 0, 10 , 20])
            plt.xlabel("Odds ratio")
        else:
            plt.xticks([-10, -20, 0, 10, 20], [" ", " ", " ", " ", " "])
            # plt.xticks([])
        # plt.yticks(np.arange(len(tbl["Term name"])), tbl["Term name"])
        row = row + 1
        plt.xlim(-25, 25)
        ax.yaxis.tick_right()
    plt.savefig("GeneFlow_Enrichr_v3.svg")
    plt.show()

    for c, tbl in df.groupby("cluster"):

        genes = []
        for gene_list in tbl["Overlapping proteins"]:
            gene_list = gene_list[2:-2].split("', '")
            for gene in gene_list:
                if gene not in genes:
                    genes.append(gene)
        genes = sorted(genes)
        temps = []
        if len(genes) > 18:
            steps = 3
        else:
            steps = 2

        for g in range(0, len(genes), steps):
            temp = ", ".join(genes[g:g + steps])
            temps.append(temp)
        show_genes = "\n".join(temps)

        print(tbl)
        f = plt.figure(1, figsize=(14, 3))
        plt.subplots_adjust(left=0.5, bottom=0.35, top=0.8, right=0.8, wspace=0.1)
        ax = f.add_subplot(111)
        plt.title(c, fontsize=5)
        tbl = tbl.sort_values(by="Odds ratio", ascending=True)
        plt.barh(tbl["Term name"], tbl["Odds ratio"], color="lightgray")
        tbl = tbl.set_index(np.arange(0, len(tbl)))

        for x, y, s in zip(tbl["Odds ratio"], tbl.index, tbl["Overlapping proteins"]):
            temps = []
            gene_list = sorted(s[2:-2].split("', '"))
            num_genes = len(gene_list)
            # gene_list = gene_list[:10]
            if x > 0.7 * tbl["Odds ratio"].max() and len(gene_list) < 10:
                show_genes = ", ".join(gene_list)
            elif x > 0.7 * tbl["Odds ratio"].max() and len(gene_list) >= 10:
                steps = int(len(gene_list) / 2) + 1
                for g in range(0, len(gene_list), steps):
                    temp = ", ".join(gene_list[g:g + steps])
                    temps.append(temp)
                show_genes = "\n".join(temps)
            elif len(gene_list) < 10:
                steps = int(len(gene_list) / 2) + 1
                for g in range(0, len(gene_list), steps):
                    temp = ", ".join(gene_list[g:g + steps])
                    temps.append(temp)
                show_genes = "\n".join(temps)
            else:
                gene_list2 = gene_list[:9]
                num_genes = num_genes - 8
                for g in range(0, len(gene_list2), 3):
                    temp = ", ".join(gene_list2[g:g + 3])
                    temps.append(temp)
                show_genes = "\n".join(temps)
                if num_genes != 0:
                    show_genes = show_genes + f" + {num_genes}"
            plt.text(0.02 * tbl["Odds ratio"].max(), y, show_genes, fontsize=6.5, ha="left", va="center")

        for x, y, s in zip(tbl["Odds ratio"], tbl.index, tbl.num_genes):
            if s > 100:
                plt.text(x - 0.6, y, s, fontsize=7.5, ha="right", va="center")

        # plt.text(x=tbl["Odds ratio"].max() * 1.1, y=(len(tbl) / 2) - 0.5, s=show_genes, fontsize=8, va="center", ha="left")
        # plt.gca().invert_yaxis()
        plt.xlabel("Odds ratio", fontsize=11)
        # plt.ylim(-1, len(tbl))
        # plt.xlim(38, 0)
        # plt.xlim(0 , tbl["Odds ratio"].max() * 1.3)
        # ax.yaxis.tick_right()
        plt.savefig("0" + "_Enrichr2.svg")
        plt.show()


if __name__ == '__main__':

    dataset = "results_conditions_NEW_2025_1000_var_logcount_mean_logcount_5_filter_degs"
    folder = rf"PLOTS/results_conditions_NEW_2025/6_cluster/{dataset}*.tsv"
    out_file = fr"PLOTS/results_conditions_NEW_2025/6_cluster/Enrichment.xlsx"
    out_file_unfiltered = fr"PLOTS/results_conditions_NEW_2025/6_cluster/Enrichment_unfiltered.xlsx"
    plot_enrichment_degs()
