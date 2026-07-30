
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)



def load_data(file, filter_mod_rate=False):
    xpore_data = pd.read_table(file, low_memory=False)
    xpore_data_mrna = xpore_data[
        (xpore_data.relative_cdna_pos >= 0) &
        (xpore_data.relative_cdna_pos <= 3)]
    if filter_mod_rate:
        xpore_data_mrna = xpore_data_mrna[(xpore_data_mrna[["mod_rate_s3-r1",
                                                            "mod_rate_s3-r2",
                                                            "mod_rate_s3-r3"]].mean(axis=1) > 0.5)
                                          |
                                          (xpore_data_mrna[["mod_rate_s4-r1",
                                                            "mod_rate_s4-r2",
                                                            "mod_rate_s4-r3"]].mean(axis=1) > 0.5)
                                          ]

    ratio_dict = dict()
    for b, c in zip(bases, colors):
        ratio_dict[b] = []
        temp = xpore_data[(xpore_data.middle_base == b) & xpore_data.relative_cdna_pos.between(0, 3)]
        for i in range(0, 3):
            ratio_dict[b].append(round(len(temp[temp.relative_cdna_pos.between(i, i + 1)]) / len(temp) * 100, 1))
        # print(sum(ratio_dict[b]))
    print(ratio_dict)
    return xpore_data_mrna


def plot_data(xpore_data_mrna, col_idx, y_lim=32, set_ylim=True):
    plt.style.use('ggplot')
    color_dict = dict(zip(bases, colors))

    plt.figure(1, figsize=(6, 6))
    plt.subplots_adjust(hspace=0.1, wspace=0.7, left=0.18, right=0.85, top=0.85, bottom=0.15)
    num_rows_ = 8
    row_idx = 1

    # plot stacked bar plot
    plt.subplot2grid((num_rows_, 2), (0, col_idx))
    data_stacked = []
    for b in bases[::-1]:
        data_stacked.append(xpore_data_mrna.loc[xpore_data_mrna["middle_base"] == b, "relative_cdna_pos"].values)
    plt.hist(data_stacked, bins=bin_size, color=colors[::-1], zorder=1, stacked=True,
             rwidth=1)
    yticks = plt.gca().yaxis.get_major_ticks()
    yticks[0].set_visible(False)
    plt.grid(False)

    # plot line
    ax = sns.histplot(xpore_data_mrna.relative_cdna_pos.values, bins=bin_size, kde=True, color="#202020", alpha=0)
    plt.axvline(x=1, color='k', linestyle='dotted', alpha=0.6)
    plt.axvline(x=2, color='k', linestyle='dotted', alpha=0.6)
    plt.xticks([0, 1, 2, 3], [0, 1, 2, 3])
    plt.xticks([])
    plt.xlim(-0.05, 3.05)
    y_pad = -.16
    plt.grid(False)
    plt.ylabel("")

    # bin_dictionary = dict()
    for b, c in zip(bases, colors):

        plt.subplot2grid((num_rows_, 2), (row_idx, col_idx))
        # sns.histplot(xpore_data.loc[xpore_data["middle_base"] == b, "relative_cdna_pos"].values,
        #              bins=bin_size, kde=True, color="k", alpha=0)

        n_cdna, bins_cdna, _ = plt.hist(xpore_data_mrna.loc[
                                            (xpore_data_mrna["middle_base"] == b),
                                            "relative_cdna_pos"].values,
                                        bins=bin_size, color=c, zorder=1, stacked=False)
        # print(n_cdna, bins_cdna)
        # bin_dictionary[b] = [int(x) for x in n_cdna]

        plt.text(x=0.08, y=.7, s=b, ha="center", va="center", color=c,
                 fontsize=10, fontweight="bold",
                 transform=plt.gca().transAxes)

        plt.axvline(x=1, color='k', linestyle='dotted', alpha=0.6)
        plt.axvline(x=2, color='k', linestyle='dotted', alpha=0.6)
        plt.xticks([])
        plt.xlim(-0.05, 3.05)
        yticks = plt.gca().yaxis.get_major_ticks()
        yticks[0].set_visible(False)
        if set_ylim:
            plt.ylim(0, y_lim)

        plt.grid(False)
        if row_idx == 2:
            plt.ylabel(ylabel_name, fontsize=fs)

        if row_idx == 4:
            y_pad = -.22
            plt.text(x=1 / 6, y=y_pad, s="5'-UTR", ha="center", va="center", color="#555555", fontsize=10,
                     transform=plt.gca().transAxes)
            plt.text(x=1 / 2, y=y_pad, s="CDS", ha="center", va="center", color="#555555", fontsize=10,
                     transform=plt.gca().transAxes)
            plt.text(x=5 / 6, y=y_pad, s="3'-UTR", ha="center", va="center", color="#555555", fontsize=10,
                     transform=plt.gca().transAxes)
            plt.xlabel("mRNA position", labelpad=18, fontsize=fs)
        row_idx += 1


if __name__ == '__main__':
    save_fig = False
    bin_size = 100
    fs = 11
    bases, colors, zorder = ["A", "C", "G", "U"], ["#109648", "#255C99", "#F7B32B", "#D62839"], [1, 1, 1, 1]

    ylabel_name = "Counts"  # Modified bases  or "DMR counts\n(C3 vs. C4)"
    output_file = "direct_RNA_seq/0_PLOTS/UTR_analyses/C3_modified_bases_mod_rate_ALL_V2_FINAL.svg"
    folder = r"direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/"

    # xpore-filtered position
    file_in_all_sites = folder + "majority_direction_kmer_diffmod_GENOME_cDNA_UTR_REL_POS.tsv"

    # confident DMR position
    file_in_dmrs = folder + "majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS.tsv"

    data = load_data(file_in_all_sites)
    plot_data(data, 0, y_lim=20_000)
    data = load_data(file_in_dmrs)
    plot_data(data, 1)

    if save_fig:
        plt.savefig(output_file)

    plt.show()

