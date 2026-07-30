import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import matplotlib.ticker as ticker
import re
import os


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


if __name__ == '__main__':

    mode_dict = {
        # "FIVE_UTR": "5'-UTR",
        "CDS": "CDS",
        "THREE_UTR": "3'-UTR",
        }

    motifs = ["",
              "GUUU[UC]",
              "GUU[UC][UC]",
              "G[UC]U[UC][UC]",
              "GUUUU",
              "A[UC]U[GA][GU]",
              "ACUGG",
              "AAUCC"]

    selected_motif = motifs[0]
    overlapping_counts = True
    take_random_pos = False
    if take_random_pos:
        tag = "_RANDOM_CDS_3UTR"
    else:
        tag = "run2"

    name_flag = selected_motif
    name_flag2 = "_run2"  # _bases_combindied

    for window in [150]:
        folder_out = f"RNAfolds/{window}bp_window_run3_corrected"
        # folder_out = f"RNAfolds/{window}bp_window"

        file_out = f"RNAfold_x_fold_DMRs_{window}bp_window_{name_flag}{name_flag2}_overlap_counts_{overlapping_counts}" + tag
        print(file_out)
        folder = rf'{folder_out}/0_sumstats/'
        file = folder + f'sum_stats_kmer_as_forgi_kmer_in_structure_overlap_counts_{overlapping_counts}_{name_flag}{name_flag2}.csv'
        print(file)
        file2 = folder + rf"sum_stats_kmer_as_forgi_kmer_in_structure_overlap_counts_{overlapping_counts}_{name_flag}{name_flag2}_motif_ratio.csv"
        # "C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore\direct_RNA_seq\miRNA_binding\RNAfolds\170bp_window\0_sumstats"
        data = pd.read_csv(file)
        motif_counts_df = pd.read_csv(file2,index_col=0)
        # print(data)

        motif_counts_df = motif_counts_df[motif_counts_df.mean_ratio > 1]
        motif_counts = motif_counts_df.index
        forgi = sorted(data.loc[(data.q_val_fisher <= 0.05) |
                                (data.motif.isin( motif_counts_df[motif_counts_df.mean_ratio > 3].index)),
        "motif"].unique())

        print(forgi)

        if take_random_pos:
            forgi = ['shhhh', 'siiii', 'iiiii', 'hhhhh', 'sssii', 'iisss', 'mmmmm', 'ssiii', 'sssss', 'sisss', 'fffff']
        motif_counts_df = motif_counts_df[motif_counts_df.index.isin(forgi)]
        motif_counts_df = motif_counts_df.drop("mean_ratio", axis=1).reset_index()

        motif_counts_df = pd.melt(motif_counts_df, id_vars=["index"])

        motif_counts_df["mRNA_region"] = motif_counts_df["variable"].apply(lambda x: "_".join(x.split("_")[:-1]))
        motif_counts_df["Base"] = motif_counts_df["variable"].apply(lambda x: x.split("_")[-1])
        motif_counts_df = motif_counts_df.rename(columns={"index": "motif", "value": "ratio (%)"})

        data = data[data.motif.isin(forgi)]

        temp = []
        for base in ["A", "C", "G", "U"]:
            for k, v in mode_dict.items():
                for motif in forgi:
                    if len(data[(data.Base == base)
                                & (data.mRNA_region == k)
                                & (data.motif == motif)]) == 0:
                        temp.append([motif, 1, 1, 1, 1, 1, 0, 74082, 15, 638, 1, base, k])
        df_temp = pd.DataFrame(temp, columns=data.columns)
        data = pd.concat([data, df_temp], axis=0, ignore_index=True)
        data = data.sort_values(by="OR_fisher", ascending=False)
        # data.motif = data.motif.apply(lambda x: "5×" + x[0])

        data = pd.merge(data, motif_counts_df[["motif", "mRNA_region", "Base", "ratio (%)"]],
                        on=["motif", "mRNA_region", "Base"])
        print(data)
        print(forgi)
        col_idx = 0
        row_idx = 0
        plt.style.use('ggplot')
        plt.figure(figsize=(5, 7))
        plt.suptitle("Over-representation of DMRs\nin RNA secondary structure elements",
                     ha="center")
        plt.subplots_adjust(hspace=0.7, wspace=0.45, left=0.15, bottom=0.1, top=0.87)
        bases = ["A", "U", "C", "G", "All"]
        colors = dict(zip(bases, ["#109648", "#D62839", "#255C99", "#F7B32B", "#808080"]))
        for k, v in mode_dict.items():
            print(col_idx, row_idx)
            for base in [# "All",
                         "A", "C", "G", "U"]:
                if base == "G":
                    continue
                ax1 = plt.subplot2grid((25, 2), (row_idx, col_idx), rowspan=2)
                # ax2 = plt.twinx()
                data_temp = data[(data.Base == base) & (data.mRNA_region == k)]
                sns.barplot(data_temp, y="ratio (%)", x="motif", color=colors[base])

                if row_idx == 0:
                    plt.title(v, fontsize=11)

                plt.grid(False)
                plt.xticks([])
                plt.xlabel("")
                plt.ylim(0, 22)
                plt.yticks([2, 20], [2, 20], fontsize=9)
                # yticks = plt.gca().yaxis.get_major_ticks()
                # yticks[0].set_visible(False)
                if col_idx < 100:

                    plt.ylabel("%", fontsize=10)
                else:
                    plt.ylabel("")

                # ax2.set_ylim(0, 15)
                row_idx += 2
                ax2 = plt.subplot2grid((25, 2), (row_idx, col_idx), rowspan=4)

                sns.barplot(data_temp, y="OR_fisher", x="motif",
                            color=colors[base])

                # yticks = plt.gca().yaxis.get_major_ticks()
                # yticks[0].set_visible(False)
                plt.yticks([1, 2, 3, 4], [1, 2, 3, 4], fontsize=9)
                row_idx += 6

                x_pos = 0
                for m, q, x_fold in zip(data_temp.motif, data_temp.q_val_fisher, data_temp.OR_fisher):
                    if q <= 0.05:
                        plt.text(x=x_pos, y=x_fold + 0.1,
                                 s="*", va="center",
                                 ha="center", fontsize=11,
                                 color="#303030")
                    x_pos += 1

                # plt.ylim(0, 4.95)
                plt.ylim(0, 2.2)

                plt.xticks(rotation=90, font="Courier New", fontsize=10)
                if col_idx < 100:
                    plt.ylabel("OR", fontsize=10)
                else:
                    plt.ylabel("")
                # if col_idx == 1:
                #     # plt.yticks([])
                #     plt.ylabel("")
                # else:
                #     plt.ylabel(base, rotation=90, weight="bold")
                plt.grid(False)
                plt.axhline(y=1, linestyle="dotted", color="#303030", alpha=0.7)
                plt.text(x=0.9, y=.8, s=base, ha="center", va="center", color=colors[base],
                         fontsize=11.5, fontweight="bold",
                         transform=plt.gca().transAxes)

                plt.xlabel("")
                row_idx += 1
            col_idx += 1
            row_idx = 0
            #
            # col_idx = 0

        plt.savefig(fr"RNAfolds\{file_out}.svg")
        plt.savefig(fr"RNAfolds\{file_out}.png")
        os.startfile(fr"RNAfolds\{file_out}.png")

