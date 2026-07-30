import os
import pandas as pd

# from GENES.Gene_set_enrichment_analysis import enrichr_api
# Code source: Gaël Varoquaux
# License: BSD 3 clause
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from sklearn import decomposition
# from sklearn import datasets
import pandas as pd
import matplotlib
from scipy import stats
from scipy import cluster
import seaborn as sns
from collections import defaultdict
import statsmodels.formula.api as smf
from scipy.stats import ttest_ind, mannwhitneyu

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def pca_2d(X, y, show_annotation=False):
    fig = plt.figure(2, figsize=(8, 6))
    plt.subplots_adjust(right=0.8, bottom=0.2, hspace=0.4, wspace=0.4)
    num_pc: int = 5  # max 5
    pca = decomposition.PCA(n_components=num_pc)
    pca.fit(X)
    X = pca.transform(X)

    print(pca.explained_variance_ratio_)
    print(pca.singular_values_)
    pc_explained = [round(x * 100, 1) for x in pca.explained_variance_ratio_]
    # Reorder the labels to have colors_dict matching the cluster results
    y = np.choose(y, np.arange(0, num_groups + 1)).astype(float)
    # print(y)
    # print(X)
    pca_row, pca_col = 0, 0
    for pc in range(0, num_pc -1):
        plt.subplot2grid((2, 2), (pca_row, pca_col))
        plt.scatter(X[:, pc], X[:, pc + 1], c=color,  # c=y,
                    # cmap=plt.cm.Set1,
                    s=50, edgecolor="k", label=[0, 1], zorder=5)

        pca_col += 1
        if pca_col == 2:
            pca_row += 1
            pca_col = 0

        if show_annotation:
            for i in range(0, len(X[:, 0])):
                plt.text(x=X[i, pc], y=X[i, pc + 1], s=replicates[i], ha="center", va="center",
                         fontdict={"size": 8, "weight": "bold", "color": "white"})

        plt.axhline(xmin=0, xmax=1, y=0, c="gray", linestyle="--")
        plt.axvline(ymin=0, ymax=1, x=0, c="gray", linestyle="--")

        plt.xlabel(f"PC{pc + 1}: {pc_explained[pc]}% variance")
        plt.ylabel(f"PC{pc + 2}: {pc_explained[pc + 1]}% variance")

    legend_temp = []
    for i in range(num_groups):
        legend_temp.append(matplotlib.lines.Line2D([0],[0], linestyle="none", c=color_temp[i], marker='o'))
    plt.legend(legend_temp, labels, numpoints=1, loc=(1.1, 0.95))
    # ax.w_xaxis.set_ticklabels([])
    # ax.w_yaxis.set_ticklabels([])
    # ax.w_zaxis.set_ticklabels([])
    plt.savefig(f"PLOTS/{data_set}/{data_set}_{num_genes}_genes_{set_filter}_PCA.svg")


def calc_significance(x):
    # print(x)
    # print(x[:4], x[4:])
    t_, p_ = stats.mannwhitneyu(x[:4].values, x[4:].values)
    # t_, p_ = stats.kendalltau(x[:4].values, x[4:].values)
    # print(x[:4], x[3:], p_)
    # quit()
    # quit()
    return p_


def calc_significance_between_groups(x):
    # print(x)
    # print(x[:4], x[4:])
    t_, p_ = stats.mannwhitneyu(x[:num_replicates].values, x[num_replicates:].values)
    # t_, p_ = stats.kendalltau(x[:4].values, x[4:].values)
    # print(x[:4], x[3:], p_)
    # quit()
    # quit()
    return p_


def filter_data(data_temp, filter_options):

    # mean_logcount
    if filter_options["mean_logcount"]:
        data_temp["mean_logcount"] = data_temp.iloc[:, :].mean(axis=1)
        data_temp = data_temp[data_temp["mean_logcount"] > mean_logcount_tresh]  # default: >= 3
        data_temp = data_temp[data_temp.columns[:-1]]

    if filter_options["median_logcount"]:
        data_temp["median_logcount"] = (data_temp.iloc[:, :].median(axis=1))
        data_temp = data_temp[data_temp["median_logcount"] >= mean_logcount_tresh]  # default: >= 3
        data_temp = data_temp[data_temp.columns[:-1]]

    if filter_options["var_logcount"]:
        if filter_options["filter_degs"]:
            degs = pd.read_table(f"DATA/{data_set}/genes_representative.tsv")
            degs = degs.dropna(subset="ext_gene").dropna(subset="target_id")
            degs = degs.sort_values(by="qval", ascending=True)
            degs = degs.drop_duplicates(subset="ext_gene").drop_duplicates(subset="target_id")
            degs = degs[degs["ext_gene"].isin(data_temp.index)]
            degs = degs[degs.qval < 0.001]
            data_temp = data_temp[(data_temp.index.isin(degs["ext_gene"].values))]
        conditions = ["C1", "C2", "C3", "C4", "C5"]
        # var_logcount
        for cond in conditions:
            data_temp[cond] = data_temp[[f"{cond}_R1-1",
                                 f"{cond}_R2-1",
                                 f"{cond}_R3-1",
                                 f"{cond}_R4-1",
                                 f"{cond}_R5-1",  # carotids
                                 f"{cond}_R6-1",
                                 f"{cond}_R7-1",  # carotids
                                 f"{cond}_R8-1"
                                 ]].apply(lambda x: abs(np.nanmean(x)), axis=1)
        # polya_data_ = polya_data_[conditions]
        data_temp["var_logcount"] = data_temp.loc[:, conditions].std(axis=1)
        data_temp = data_temp.nlargest(num_genes, "var_logcount").drop(conditions, axis=1)
        data_temp = data_temp[data_temp.columns[:-1]]

    if filter_options["filter_degs"] and not filter_options["var_logcount"]:
        degs = pd.read_table(f"DATA/{data_set}/genes_representative.tsv")
        degs = degs.dropna(subset="ext_gene").dropna(subset="target_id")
        degs = degs.sort_values(by="qval", ascending=True)
        degs = degs.drop_duplicates(subset="ext_gene").drop_duplicates(subset="target_id")
        degs = degs[degs["ext_gene"].isin(data_temp.index)]
        degs = degs[degs.qval < 0.05]
        # degs = degs.nsmallest(num_genes, "qval")  # to excl. batch effects
        data_temp = data_temp[(data_temp.index.isin(degs["ext_gene"].values))]
        # print((set(degs["ext_gene"].values) - set(data_temp.index)))

    if filter_options["compare_groups"]:
        degs = pd.read_table(f"DATA/report_C3_C4/genes_representative.tsv")
        degs = degs[degs.pval < 0.05]
        degs = degs["ext_gene"].values
        data_temp = data_temp[data_temp.index.isin(degs)]

    return data_temp


def get_cluster_classes(den, label='ivl'):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))
    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l
    return cluster_classes


def perform_gene_clustering(g, matrix, do_gsea=True):
    # extract protein cluster
    plt.figure(5, figsize=(8, 6))
    plt.subplots_adjust(bottom=0.3)
    den = cluster.hierarchy.dendrogram(g.dendrogram_row.linkage,
                                       labels=matrix.index,
                                       color_threshold=cluster_thresh[num_genes])  # default: 3.9

    from collections import defaultdict

    clusters = get_cluster_classes(den)
    if do_gsea:
        for k, v in clusters.items():
            print(k)
            """gsea = enrichr_api.fetch_enrichr(gene_list=v, term="all")
            print(gsea)
            if gsea is not None:
                gsea.to_csv(f"PLOTS/{data_set}/{data_set}_{num_genes}_{set_filter}_genes_gene_cluster_gsea_{k}.tsv",
                            "\t", index=True)"""
            """for protein in v:
                print(protein)"""
            print("#######\n")

    plt.savefig(f"PLOTS/{data_set}/{data_set}_{num_genes}_{set_filter}_genes_gene_clusters.svg")
    return clusters


def plot_cluster_map(illumina_cluster):
    print(illumina_cluster)
    sns.set(font="DejaVu Sans")

    condition_colors = color_temp

    # sex_colors = pd.DataFrame(network_labels, index=networks.columns).map(network_lut)

    m = "lightsteelblue"
    f = "lightcoral"
    sex_colors = [m, m, f, m,  f,
                  m, f,
                  f] * len(condition_colors)
    # if filter_male:
    #     sex_colors = [m, m, m,
    #                   m,] * len(condition_colors)
    # else:
    #     sex_colors = [f, f, f,
    #                   f] * len(condition_colors)
    print(sex_colors)
    print(illumina_cluster.columns)

    col_colors = pd.DataFrame(data={"Condition": color, "Sex": sex_colors
                                    }, index=illumina_cluster.columns)
    # print(col_colors)
    # print(len(condition_colors), len(sex_colors))
    print("##########   Number of proteins in clustermap: ", len(illumina_cluster))

    if perform_clustering:
        fig_temp = plt.figure(10, figsize=(6, 6))
        clustergrid2 = sns.clustermap(illumina_cluster,
                                      method="ward",  # default average
                                      metric="euclidean",  # correlation or euclidean
                                      col_cluster=False,
                                      row_cluster=True,
                                      cmap="cividis",  # "cividis",
                                      # cbar_pos=(0.82, .15, .01, .3),  # (left, bottom, width, height)
                                      # col_colors=col_colors,
                                      # standard_scale=0,
                                      )

        plt.close(fig_temp)
        clusters_ = perform_gene_clustering(clustergrid2, illumina_cluster, do_gsea=True)
        den_tree_color = []
        den_gene_index = []
        a = clustergrid2.dendrogram_row.reordered_ind

        cmap = sns.color_palette("colorblind").as_hex()[::-1]

        idx = 0
        for c, genes in clusters_.items():
            print(c, genes)
            den_tree_color.extend([cmap[idx]] * (len(genes) - 2))
            for gene in genes:
                den_gene_index.append(gene)
            idx += 1
        row_colors = pd.DataFrame(data={" ": den_tree_color[:num_genes]}, index=illumina_cluster.iloc[a].index[:num_genes])
    else:
        row_colors = None
    # print(col_colors)
    # print(row_colors)
    clustergrid = sns.clustermap(illumina_cluster,
                                 figsize=(8, 13),
                                 method="ward",  # default average
                                 metric="euclidean",  # correlation or euclidean
                                 col_cluster=True,
                                 row_cluster=True,
                                 cmap="cividis",  # "cividis",
                                 # cbar_pos=(0.82, .15, .01, .3),  # (left, bottom, width, height)
                                 col_colors=col_colors,
                                 row_colors=row_colors,
                                 dendrogram_ratio=.08,
                                 colors_ratio=0.04,
                                 yticklabels=True,
                                 xticklabels=True,
                                 # z_score=0,
                                 standard_scale=0,
                                 # cbar=False,
                                 cbar_kws={"shrink": .2},
                                 # tree_kws={"colors_dict": den_tree_color}
                                 )
    # plt.show()

    clustergrid.fig.subplots_adjust(right=0.4)


    # title_font = font_manager.FontProperties(alignment="center")

    custom_legend = plt.legend(handles=legend_elements,
                               fontsize=9,
                               title='Condition (T+P/I+P)',
                               # title_fontproperties=title_font,
                               alignment="left",
                               title_fontsize=10,
                               frameon=False,
                               facecolor='white',
                               framealpha=1,
                               labelspacing=0.25,
                               handletextpad=0.15,
                               bbox_to_anchor=(-2.8, 1.9),  # (left, bottom, width, height)
                               loc='center left'
                               )

    custom_legend2 = plt.legend(handles=legend_elements_sex,
                                fontsize=9,
                                title='Sex',
                                alignment="left",
                                title_fontsize=10,
                                frameon=False,
                                facecolor='white',
                                framealpha=1,
                                labelspacing=0.25,
                                handletextpad=0.15,
                                bbox_to_anchor=(-2.8, 1.4),  # (left, bottom, width, height)
                                loc='center left'
                                )
    y_pos_cluster = 0
    """ num_cluster = len(row_colors[" "].unique())
    for col_idx, c in enumerate(cmap):
        if col_idx >= num_cluster:
            break
        len_cluster = len(row_colors[row_colors[" "] == c])
        clustergrid.ax_heatmap.text(x=-0.8, y=y_pos_cluster + int(len_cluster/ 2), s=col_idx + 1,
                                    ha="center", va="center", fontdict={"size": 8})
        y_pos_cluster = y_pos_cluster + len_cluster"""
    plt.gca().add_artist(custom_legend)
    plt.gca().add_artist(custom_legend2)

    # set cbar
    clustergrid.ax_cbar.set_position((0.82, .15, .01, .3))

    # clustergrid.ax_heatmap.set_xticks([])
    if num_genes > 200:
        pass
        # clustergrid.ax_heatmap.set_yticks([])
    clustergrid.ax_heatmap.tick_params(axis='y', labelsize=9)
    clustergrid.ax_heatmap.tick_params(axis='x', labelsize=6)
    clustergrid.ax_cbar.set_title(f"log$_2$(norm_counts)\nstandard scale"
                                  , fontdict={'fontsize': 10, "ha": "center"}, pad=10)
    clustergrid.ax_cbar.tick_params(axis='y', labelsize=6)
    # plt.xticks([])

    tag = "both_standard_scale_only_mRNA_C3_C4_mean_log2_counts_5"
    plt.savefig(f"PLOTS/{data_set}/{data_set}_{num_genes}_{set_filter}_genes_clustermap_{tag}.svg")


def read_data():

    if os.path.exists(file.replace(".tsv", "_temp.tsv")):
        rnaseq_data = pd.read_table(file.replace(".tsv", "_temp.tsv"), index_col=0)
        # print(rnaseq_data[rnaseq_data.index == "MYH11"])
        return rnaseq_data

    rnaseq_data = pd.read_table(file, # index_col="protein"
                                )  # ext_gene

    col_name = "protein"  # protein
    alt_name = "transcript"  # transcript
    rnaseq_data = rnaseq_data.fillna(".")
    print(rnaseq_data[rnaseq_data[col_name] == "."])
    print(rnaseq_data[rnaseq_data.transcript == "ENST00000300036.6"])

    rnaseq_data.loc[rnaseq_data[col_name] == ".", col_name] = rnaseq_data.loc[rnaseq_data[col_name] == ".", alt_name]

    reorder_conds = []
    for c in range(1, num_groups + 1):
        for r in range(1, num_replicates + 1):
            reorder_conds.append(f"C{c}_R{r}-1")
    # print(reorder_conds)
    # if readout != "logcount_matrix":
    #     for cond in reorder_conds:
    #         rnaseq_data[cond] = np.log2(rnaseq_data[cond] + 1)

    # rnaseq_data[reorder_conds] = rnaseq_data[reorder_conds].apply(lambda x: x + abs(min(x)) if min(x) < 0 else x, axis=1)
    rnaseq_data["mean_count"] = rnaseq_data[reorder_conds].apply(lambda x: abs(x).mean(), axis=1)
    rnaseq_data = rnaseq_data[rnaseq_data["mean_count"] >= 0]
    rnaseq_data = rnaseq_data.dropna()
    gene_repr = pd.read_table(fr"DATA/{data_set}/genes_representative.tsv", usecols=[3, 7])
    canonical_transcr = gene_repr.loc[gene_repr["canonical"], "target_id"].values
    rnaseq_data = rnaseq_data[rnaseq_data.transcript.isin(canonical_transcr)]

    rnaseq_data = rnaseq_data.sort_values(by="mean_count", ascending=False).drop("mean_count", axis=1)
    rnaseq_data = rnaseq_data.drop_duplicates(subset=col_name).dropna()
    rnaseq_data = rnaseq_data.set_index(rnaseq_data[col_name]).drop(col_name, axis=1)
    rnaseq_data = rnaseq_data[reorder_conds]
    rnaseq_data.to_csv(file.replace(".tsv", "_temp.tsv"), "\t", index=True)
    print(rnaseq_data)
    return rnaseq_data


def clean_data(df):
    X = df.copy()

    # 1) Ensure numeric and expose hidden NaNs
    X = X.apply(pd.to_numeric, errors='coerce')

    # 2) If you log-transform, do it safely (avoid -inf)
    # X = np.log2(X + 1)

    # 3) Replace explicit ±inf with NaN
    X = X.replace([np.inf, -np.inf], np.nan)

    # 4) Drop rows/cols that are all-NaN
    X = X.dropna(axis=0, how='all').dropna(axis=1, how='all')

    # 5) Drop zero-variance rows/cols (these break z-score / standard_scale)
    zero_var_rows = X.var(axis=1, ddof=0) == 0
    zero_var_cols = X.var(axis=0, ddof=0) == 0
    X = X.loc[~zero_var_rows, ~zero_var_cols]

    # 6) Option A: drop any remaining NaNs
    X_clean = X.dropna(axis=0).dropna(axis=1)
    return X_clean


def annotate_genes_manually():
    rnaseq_data_temp = rnaseq_data.copy()
    # rnaseq_data_temp = rnaseq_data_temp.reset_index()
    # print(rnaseq_data_temp)
    rnaseq_data_temp["gene_names"] = rnaseq_data_temp.index
    rnaseq_data_temp["RNA_mod"] = np.nan
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names.str.startswith("ALKBH"), "RNA_mod"] = "m6A eraser"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names == "FTO", "RNA_mod"] = "m6A eraser"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names == "FMR1", "RNA_mod"] = "m6A reader"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names.str.startswith("HNRNP"), "RNA_mod"] = "m6A reader"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names.str.startswith("YTHDC"), "RNA_mod"] = "m6A reader"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names.str.startswith("YTHDF"), "RNA_mod"] = "m6A reader"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names.str.startswith("PUS"), "RNA_mod"] = "Pseudo-U writer"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names == "NOP2", "RNA_mod"] = "m5C writer"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names == "IGF2BP1", "RNA_mod"] = "m6A reader"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names.str.startswith("METTL"), "RNA_mod"] = "m6A writer"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names == "WTAP", "RNA_mod"] = "m6A writer"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names == "ALYREF", "RNA_mod"] = "tRNA m7G writer"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names == "DNMT2", "RNA_mod"] = "m5C writer"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names.str.startswith("NSUN"), "RNA_mod"] = "m5C writer"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names == "FBL", "RNA_mod"] = "2′Ome writer"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names.isin(["TET1", "TET2", "TET3"]), "RNA_mod"] = "m5C erasers"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names.isin(["YBX1", "YBX2", "YBX3"]), "RNA_mod"] = "m5C reader"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names == "METTL1", "RNA_mod"] = "m7G writer"
    rnaseq_data_temp.loc[rnaseq_data_temp.gene_names.isin(["CMTR1", "CMTR2", "RNGTT", "RNMT"]), "RNA_mod"] = "m7GpppN writer"

    rnaseq_data_temp = rnaseq_data_temp[rnaseq_data_temp.gene_names != "NSUN5P2"].dropna()
    rnaseq_data_temp = rnaseq_data_temp[rnaseq_data_temp.gene_names != "NSUN5P1"]
    # rnaseq_data = rnaseq_data[rnaseq_data.index.isin(rnaseq_data_temp.dropna()["protein"])]
    rnaseq_data_temp = rnaseq_data_temp.sort_values(by=["RNA_mod", "gene_names"])
    rnaseq_data_temp = rnaseq_data_temp[~rnaseq_data_temp["gene_names"].str.contains("-")]

    print(rnaseq_data_temp[["RNA_mod", "gene_names"]])
    return rnaseq_data_temp


def plot_violine():
    tbl_collector = []
    for mod, tbl in rnaseq_data2[rnaseq_data2.columns].groupby("RNA_mod"):
        for c in ["C1", "C2", "C5"]:
            for col in tbl.columns:
                if c in col:
                    tbl = tbl.drop(col, axis=1)
        cols = tbl.columns[:16]
       # tbl[cols] = tbl[cols].div(tbl[cols].max(axis=1), axis=0)
        print(tbl)
        tbl = pd.melt(tbl.iloc[:, :16])
        for c in ["C1", "C2", "C3", "C4", "C5"]:
            tbl.loc[tbl["variable"].str.startswith(c), "variable"] = c
        tbl["RNAmod"] = mod
        print(tbl)
        tbl_collector.append(tbl)

    rna_mod_genes = pd.concat(tbl_collector)
    for mod, tbl in rna_mod_genes.groupby("RNAmod"):
        ax = sns.violinplot(data=rna_mod_genes, x="variable", y="value", palette=color_temp, width=0.9, cut=0)
        plt.title(f"{mod}")
        plt.show()


def plot_violine_only_degs():
    c3_c4 = pd.read_table(r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Illumina\hSMC_TGFB_PDGF\DATA\results_C3_C4_2025\tables\diffexp\model_condition.genes-representative.diffexp.tsv")
    c3_c4_degs = c3_c4.loc[c3_c4.pval <= 0.05, "ext_gene"]
    print(c3_c4_degs)

    tbl_collector = []
    for mod, tbl in rnaseq_data2.groupby("RNA_mod"):
        # print(mod)
        print(tbl)
        tbl = tbl.drop("RNA_mod", axis=1)
        # tbl = tbl[tbl.index.isin(c3_c4_degs)]
        if tbl.empty:
            continue
        # tbl = tbl.drop(["RNA_mod", "Link", "role in", "gene_RNAmod"], axis=1)
        tbl = tbl.T.reset_index()
        for c in ["C1", "C2", "C3", "C4", "C5"]:
            tbl.loc[tbl["index"].str.startswith(c), "index"] = c

        print(tbl)

        tbl["index"] = tbl["index"].apply(lambda x: x.split("_")[0])

        # plt.title(mod)
        num_genes = int(len(tbl.columns[1:]))
        if num_genes == 1:
            num_cols = 1
            num_rows = 1
            w, h = 3, 3
        elif 5 > num_genes > 1:
            num_cols = num_genes
            num_rows = 1
            w, h = int((3.5 * num_cols)), 3
        elif 11 > num_genes >= 5:
            num_cols = int(np.ceil(num_genes / 2))
            num_rows = int(np.ceil(num_genes / num_cols))
            w, h = int((4 * num_cols)), 6
        else:
            num_cols = int(np.ceil(num_genes / 3))
            num_rows = int(np.ceil(num_genes / num_cols))
            w, h = int((4 * num_cols)), 10

        plt.figure(1, figsize=(w, h))
        # plt.figure(1, figsize=(5, 8))
        plt.subplots_adjust(wspace=0.6, hspace=0.5, left=0.2, bottom=0.22, top=0.8)
        plt.suptitle(mod, fontsize=16, fontweight='bold', ha='center')
        n_col, n_row = 0, 0
        # print(num_rows, num_cols)
        vals, conds = [], []

        print(tbl)
        for col in sorted(tbl.columns[1:]):

            plt.subplot2grid((num_rows, num_cols), (n_row, n_col))
            ax = sns.violinplot(data=tbl, x="index", y=tbl[col], palette=color_temp, width=0.9, cut=0)
            n_col += 1
            if n_col == num_cols:
                n_col = 0
                n_row += 1
            ax.set_title(f"{col}")
            plt.ylim(min(tbl[col]) - 0.3, max(tbl[col]) + 0.3)
            plt.ylabel("log$_2$(counts + 1)")
            plt.xlabel("")
            tbl[col] = tbl[col] / tbl[col].max()
            vals.extend(tbl[col].values)
            conds.extend(list(tbl["index"].values))
        out = pd.DataFrame({"expr": vals, "condition": conds})
        out["RNAmod"] = mod
        tbl_collector.append(out)

        mod = mod.replace("/", "_")

        # plt.savefig(f"PLOTS/RNAmod_Genes/0_C1_C2_C3_C4_C5_{mod}_C3_C4_DEGS.png")
        # plt.savefig(f"PLOTS/RNAmod_Genes/0_C1_C2_C3_C4_C5_{mod}_C3_C4_DEGS.svg")
        plt.savefig(f"PLOTS/RNAmod_Genes/0_C1_C2_C3_C4_C5_{mod}.png")
        plt.savefig(f"PLOTS/RNAmod_Genes/0_C1_C2_C3_C4_C5_{mod}.svg")
       # plt.show()
        plt.close()
        # print(tbl.head(5))
        # print(tbl.T)
        print("\n")


if __name__ == '__main__':

    # only used for clustering > here not needed
    cluster_thresh = {2000: 12,
                      1000: 10,  # default 10
                      500: 7.6,
                      200: 5.5,
                      100: 4,
                      12707: 14,
                      17732: 14,
                      10096: 10,
                      10063: 10}

    num_genes = 1000
    mean_logcount_tresh = median_log_count = 5   # log2 matrix

    num_groups = 2
    num_replicates = 8
    data_set = "results_conditions_NEW_2025"
    readout = "model_condition.logcount-matrix"
    file = fr"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Illumina\hSMC_TGFB_PDGF\DATA\results_conditions_NEW_2025\model_condition.logcount-matrix.tsv"

    if not os.path.isdir(f"PLOTS/{data_set}"):
        os.mkdir(f"PLOTS/{data_set}")

    # rnaseq_data = pd.read_table(file)  # ext_gene

    rnaseq_data = read_data()
    filter_options_ = {"var_logcount": False,
                      "mean_logcount": False,
                      'median_logcount': True,
                      "filter_degs": False,
                      "compare_groups": False}

    set_filter_flags = []
    if filter_options_["median_logcount"]:
        filter_options_["mean_logcount"] = False

    for k, v in filter_options_.items():
        if v:
            if k == "mean_logcount":
                k = "mean_logcount" + "_" + str(mean_logcount_tresh)
            if k == "median_logcount":
                k = "median_logcount" + "_" + str(mean_logcount_tresh)
            set_filter_flags.append(k)

    set_filter = "_".join(set_filter_flags)
    rnaseq_data = filter_data(rnaseq_data,
                               filter_options_)

    data = []
    for col_idx in range(num_groups * num_replicates):
        data.append(rnaseq_data.iloc[:, col_idx])
    labels = []
    for i in range(1, num_groups + 1):
        labels.append(f"C{i}")
    replicates = []
    for i in range(1, num_replicates + 1):
        replicates.append(f"R{i}")
    replicates = replicates * num_groups

    groups = list()
    for i in range(num_groups):
        for j in range(num_replicates):
            groups.append(i)
    # print(groups)
    # groups = [0, 0, 0, 1, 1, 1]
    if num_groups > 2:
        color_temp = ["#BC80BD",  # C1
                            "#F7B6D2",  # C2
                            "#8DD3C7",  # C3
                            "#FDB462",  # C4
                            "#B3B3B3"  # C5
                            ]
    else:
        color_temp = [# "#BC80BD",  # C1
                      # "#F7B6D2",  # C2
                      "#8DD3C7",  # C3
                      "#FDB462",  # C4
                      # "#B3B3B3"  # C5
                      ]

    color = np.choose(groups, color_temp)
    # print(color)
    # groups = [0, 0, 1, 1]
    # pca_2d(tbl, groups, show_annotation=True)
    num_genes = len(rnaseq_data)
    fig3 = plt.figure(100)
    # add legend
    ms = 50
    legend_elements = [plt.scatter([0], [0], marker='s', color='#BC80BD', label='C1: +++', s=ms),
                       plt.scatter([0], [0], marker='s', color='#F7B6D2', label='C2: ++–', s=ms),
                       plt.scatter([0], [0], marker='s', color='#8DD3C7', label='C3: +––', s=ms),
                       plt.scatter([0], [0], marker='s', color='#FDB462', label='C4: –+–', s=ms),
                       plt.scatter([0], [0], marker='s', color='#B3B3B3', label='C5: –––', s=ms)
                       ]

    legend_elements_sex = [plt.scatter([0], [0], marker='s', color='lightcoral', label='female', s=ms),
                           plt.scatter([0], [0], marker='s', color='lightsteelblue', label='male', s=ms),
                           ]
    plt.close(fig3)

    mod_genes = pd.read_excel(r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\RNAmod_genes\human_RNA_mod_genes_expanded_PMID.xlsx",
                              skiprows=0)
    mod_genes = mod_genes[mod_genes["Main_RNA_substrate"].str.contains("mRNA")]
    mod_genes = mod_genes[~mod_genes["Modification"].str.contains("Cap")]
    mod_genes = mod_genes[~mod_genes["Modification"].str.contains("cap")]
    genes_temp = mod_genes["Gene_symbol"].unique()
    genes_temp = [x.strip() for x in genes_temp]

    plt.style.use('ggplot')

    rnaseq_data = rnaseq_data.reset_index()

    print(mod_genes)
    mod_genes = mod_genes.rename(columns={"Gene_symbol": "gene"})
    mod_genes["RNA_mod"] = mod_genes["Modification"] + " " + mod_genes["Role"]
    rnaseq_data = pd.merge(rnaseq_data, mod_genes[["gene", "Role", "Modification", "RNA_mod"]], on="gene", how="inner")

    # load RNAmod map
    rnaseq_data = rnaseq_data.dropna(subset="RNA_mod")

    # keep only mRNA related RNAmod genes
    rnaseq_data = rnaseq_data[~rnaseq_data["RNA_mod"].str.startswith("hnRNP")]
    rnaseq_data = rnaseq_data[~rnaseq_data["RNA_mod"].str.startswith("tRNA")]

    # use this as index
    rnaseq_data["gene_RNAmod"] = rnaseq_data["gene"] + " | " + rnaseq_data["RNA_mod"]

    perform_clustering = False
    rnaseq_data = rnaseq_data.fillna(0)

    # filter C3 and C4 if required
    new_cols = [x for x in rnaseq_data.columns if "C3" in x or "C4" in x]

    # analyze genesets:
    expression_df = rnaseq_data[new_cols].T.reset_index()
    print(expression_df)

    # rnaseq_data3 = rnaseq_data2.copy().T
    rnaseq_data2 = rnaseq_data.copy()
    print(rnaseq_data2.columns)
    rnaseq_data2 = rnaseq_data2.set_index("gene_RNAmod")
    rnaseq_data2 = rnaseq_data2.drop(["Role", "Modification", "gene", "RNA_mod"], axis=1)
    print(rnaseq_data2)

    # plot data
    plot_cluster_map(rnaseq_data2[new_cols])

    plot_violine()

    plot_violine_only_degs()





