import os
from sklearn.decomposition import PCA
import numpy as np
import enrichr_api
# Code source: Gaël Varoquaux
# License: BSD 3 clause
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import decomposition
from sklearn import datasets
import pandas as pd
import matplotlib
from scipy import stats
from scipy import cluster
import seaborn as sns
from collections import defaultdict


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def pca_3d(X, y):

    fig = plt.figure(1, figsize=(6, 4))
    plt.clf()
    ax = Axes3D(fig, # rect=[0, 0, 0.95, 1],
                elev=18, azim=125)

    plt.cla()
    pca = decomposition.PCA(n_components=3)
    pca.fit(X)
    X = pca.transform(X)

    print(pca.explained_variance_ratio_)
    print(pca.singular_values_)

    # Reorder the labels to have colors_dict matching the cluster results
    y = np.choose(y, np.arange(0, num_groups + 1)).astype(float)
    print(y)
    print(X)

    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=color, # c=y,
               # cmap=plt.cm.Set1,
               s=200,
               edgecolor="k", label=[0, 1])

    pc1, pc2, pc3 = [round(x * 100, 1) for x in pca.explained_variance_ratio_]
    ax.set_xlabel(f"PC1 ({pc1}%)")
    ax.set_ylabel(f"PC2 ({pc2}%)")
    ax.set_zlabel(f"PC3 ({pc3}%)")

    legend_temp = []
    for i in range(num_groups):
        legend_temp.append(matplotlib.lines.Line2D([0],[0], linestyle="none", c=color_temp[i], marker='o'))
    ax.legend(legend_temp, labels, numpoints=1)
    # ax.w_xaxis.set_ticklabels([])
    # ax.w_yaxis.set_ticklabels([])
    # ax.w_zaxis.set_ticklabels([])


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
        data_temp = data_temp[data_temp["mean_logcount"] >= mean_logcount_tresh]  # default: >= 3
        data_temp = data_temp[data_temp.columns[:-1]]

    if filter_options["median_logcount"]:
        data_temp["median_logcount"] = (data_temp.iloc[:, :].median(axis=1))
        data_temp = data_temp[data_temp["median_logcount"] >= mean_logcount_tresh]  # default: >= 3
        data_temp = data_temp[data_temp.columns[:-1]]

    if filter_options["var_logcount"]:
        if filter_options["filter_degs"]:
            print("filter_degs")
            degs = pd.read_table(f"DATA/{data_set}/genes_representative.tsv")
            degs = degs.dropna(subset="ext_gene").dropna(subset="target_id")
            degs = degs.sort_values(by="qval", ascending=True)
            degs = degs.drop_duplicates(subset="ext_gene").drop_duplicates(subset="target_id")
            degs = degs[degs["ext_gene"].isin(data_temp.index)]
            degs = degs[degs.qval < 0.001]
            data_temp = data_temp[(data_temp.index.isin(degs["ext_gene"].values))]
            print(len(data_temp))

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
        # data_ = data_[conditions]
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
        degs = degs.nsmallest(num_genes, "qval")  # to excl. batch effects
        data_temp = data_temp[(data_temp.index.isin(degs["ext_gene"].values))]
        print((set(degs["ext_gene"].values) - set(data_temp.index)))

    if filter_options["compare_groups"]:
        degs = pd.read_table(f"DATA/report_C1_C3/genes_representative.tsv")
        degs = degs[degs.qval < 0.001]
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
            gsea = enrichr_api.fetch_enrichr(gene_list=v, term="all")
            print(gsea)
            if gsea is not None:
                gsea.to_csv(f"PLOTS/{data_set}/{data_set}_{num_genes}_{set_filter}_genes_gene_cluster_gsea_{k}.tsv",
                            "\t", index=True)
            """for protein in v:
                print(protein)"""
            print("#######\n")

    plt.savefig(f"PLOTS/{data_set}/{data_set}_{num_genes}_{set_filter}_genes_gene_clusters.svg")
    return clusters


def plot_cluster_map(illumina_cluster):

    sns.set(font="DejaVu Sans")
    condition_colors = ["red",  # C1
                        "orange",  # C2
                        "green",  # C3
                        "steelblue",  # C4
                        "lightgrey"  # C5
                        ]

    # sex_colors = pd.DataFrame(network_labels, index=networks.columns).map(network_lut)

    m = "lightsteelblue"
    f = "lightcoral"
    temp = []
    sex_colors = [m, m, f, m, f, m, f, f] * len(condition_colors)

    print(sex_colors)
    print(illumina_cluster.columns)

    col_colors = pd.DataFrame(data={"Condition": color, "Sex": sex_colors}, index=illumina_cluster.columns)
    # print(col_colors)
    # print(len(condition_colors), len(sex_colors))
    print("##########   Number of proteins in clustermap: ", len(illumina_cluster))
    fig_temp = plt.figure(10)
    clustergrid2 = sns.clustermap(illumina_cluster,
                                 method="ward",  # default average
                                 metric="euclidean",  # correlation or euclidean
                                 col_cluster=True,
                                 row_cluster=True,
                                 cmap="cividis",  # "cividis",
                                 # cbar_pos=(0.82, .15, .01, .3),  # (left, bottom, width, height)
                                 # col_colors=col_colors,
                                 standard_scale=0,
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
        den_tree_color.extend([cmap[idx]] * (len(genes) - 1))
        for gene in genes:
            den_gene_index.append(gene)
        idx += 1

    row_colors = pd.DataFrame(data={" ": den_tree_color[:num_genes]}, index=illumina_cluster.iloc[a].index[:num_genes])
    # print(col_colors)
    # print(row_colors)
    clustergrid = sns.clustermap(illumina_cluster,
                                 figsize=(8, 6),
                                 method="ward",  # default average
                                 metric="euclidean",  # correlation or euclidean
                                 col_cluster=True,
                                 row_cluster=True,
                                 cmap="cividis",  # "cividis",
                                 # cbar_pos=(0.82, .15, .01, .3),  # (left, bottom, width, height)
                                 col_colors=col_colors,
                                 row_colors=row_colors,
                                 dendrogram_ratio=.1,
                                 yticklabels=True,
                                 xticklabels=True,
                                 standard_scale=0,
                                 # cbar=False,
                                 cbar_kws={"shrink": .2},
                                 # tree_kws={"colors_dict": den_tree_color}
                                 )

    clustergrid.fig.subplots_adjust(right=0.62)


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
    num_cluster = len(row_colors[" "].unique())
    for idx, c in enumerate(cmap):
        if idx >= num_cluster:
            break
        len_cluster = len(row_colors[row_colors[" "] == c])
        clustergrid.ax_heatmap.text(x=-0.8, y=y_pos_cluster + int(len_cluster/ 2), s=idx + 1,
                                    ha="center", va="center", fontdict={"size": 8})
        y_pos_cluster = y_pos_cluster + len_cluster
    plt.gca().add_artist(custom_legend)
    plt.gca().add_artist(custom_legend2)

    # set cbar
    clustergrid.ax_cbar.set_position((0.82, .15, .01, .3))

    clustergrid.ax_heatmap.set_xticks([])
    if num_genes > 100:
        clustergrid.ax_heatmap.set_yticks([])
    clustergrid.ax_heatmap.tick_params(axis='y', labelsize=6)
    clustergrid.ax_heatmap.tick_params(axis='x', labelsize=6)
    clustergrid.ax_cbar.set_title("Norm log2(counts)", fontdict={'fontsize': 10, "ha": "center"}, pad=10)
    clustergrid.ax_cbar.tick_params(axis='y', labelsize=6)
    # plt.xticks([])

    plt.savefig(f"PLOTS/{data_set}/{data_set}_{num_genes}_{set_filter}_genes_clustermap.svg")


def read_data():

    if os.path.exists(file.replace(".tsv", "_temp.tsv")):
        rnaseq_data = pd.read_table(file.replace(".tsv", "_temp.tsv"), index_col=0)
        print(rnaseq_data[rnaseq_data.index == "MYH11"])
        return rnaseq_data

    rnaseq_data = pd.read_table(file, # index_col="protein"
                                )  # ext_gene

    col_name = "gene"  # protein or gene
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
    if readout != "model_condition.logcount-matrix":
        for cond in reorder_conds:
            rnaseq_data[cond] = np.log2(rnaseq_data[cond] + 1)

    # rnaseq_data[reorder_conds] = rnaseq_data[reorder_conds].apply(lambda x: x + abs(min(x)) if min(x) < 0 else x, axis=1)
    rnaseq_data["mean_count"] = rnaseq_data[reorder_conds].apply(lambda x: abs(x).mean(), axis=1)
    rnaseq_data = rnaseq_data[rnaseq_data["mean_count"] >= 0]

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


if __name__ == '__main__':

    cluster_thresh = {2000: 12,
                      1000: 11.0,  # default 10
                      500: 7.6,
                      200: 5.5,
                      100: 4,
                      12707: 14}

    num_genes = 1000
    mean_logcount_tresh = median_log_count = 5   # log2 matrix

    num_groups = 5
    num_replicates = 8
    data_set = "results_conditions_NEW_2025"
    readout = "model_condition.logcount-matrix"
    file = fr"DATA/{data_set}/{readout}.tsv"

    if not os.path.isdir(f"PLOTS/{data_set}"):
        os.mkdir(f"PLOTS/{data_set}")

    # rnaseq_data = pd.read_table(file)  # ext_gene

    rnaseq_data = read_data()

    filter_options_ = {"var_logcount": True,
                      "mean_logcount": True,
                      'median_logcount': False,
                      "filter_degs": True,
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
    print(rnaseq_data)
    # rnaseq_data = rnaseq_data["mean_logcount"] = rnaseq_data.iloc[:, :].mean(axis=1)
    # rnaseq_data = rnaseq_data[rnaseq_data["mean_logcount"] >= 5].dropna()
    # rnaseq_data = rnaseq_data[rnaseq_data.columns[:-1]]

    # gene_list = gene_list[gene_list.Species == "Homo sapiens"]
    # gene_list.dropna(inplace=True)
    # gene_list = gene_list[gene_list.p_value < 0.05]

    data = []

    for col in range(num_groups * num_replicates):
        data.append(rnaseq_data.iloc[:, col])

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
    color_temp = ["#BC80BD",  # C1
                 "#F7B6D2",  # C2
                 "#8DD3C7",  # C3
                 "#FDB462",  # C4
                 "#B3B3B3"  # C5
                 ]
    color = np.choose(groups, color_temp)
    # print(color)
    # groups = [0, 0, 1, 1]
    pca_2d(data, groups)
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

    plot_cluster_map(rnaseq_data)
    plt.show()
