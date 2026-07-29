import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import os
from scipy.stats import hypergeom
import re
import statsmodels.stats.multitest as multi
from scipy import stats
from scipy.stats import fisher_exact


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def get_kmer_forgi_elem(args):
    numbers = [str(x) for x in range(0, 10)]
    forgi_seq = []
    for elem in args["forgi_structure"]:
        if elem not in numbers:
            forgi_seq.append(elem)
    forgi_seq = "".join(forgi_seq)
    # forgi_seq = args["forgi_structure"][::2]
    p = args["pos"] - args["start"]
    return forgi_seq[p - 3: p + 2]


if __name__ == '__main__':

    motifs = ["",
              "G[UC]U[UC][UC]",
              "GUUUU",
              "GUUU[UC]",
              "GUU[UC][UC]",
              "A[UC]U[GA][GU]",
              "ACUGG",
              "AAUCC"]

    selected_motif = motifs[0]

    col_name_dict = {1: "kmer_in_structure",
                     2: "3mer_in_structure",
                     3: "middle_base_forgi_elem"}

    col_name = col_name_dict[1]  # default kmer_in_structure = 5-mer

    # window = 180
    overlapping_counts = True  # default True
    mirrow_motifs = False
    name_flag = selected_motif
    name_flag2 = "_C3_modRate_07"  # "C4_vs_C3_pos_zscore"
    overwrite = True
    take_random_pos = False

    for window in [150]:
        if take_random_pos:
            folder_out = f"RNAfolds/{window}bp_window_RANDOM_CDS_3UTR"
            tag = "_org"
        else:
            folder_out = f"RNAfolds/{window}bp_window_C3_modBASE_07"
            tag = ""

        if not os.path.exists(f"{folder_out}/0_sumstats"):
            os.mkdir(f"{folder_out}/0_sumstats")

        processed_file = f'{folder_out}/0_sumstats/seq_params_best_10.csv'
        print(processed_file)
        if not os.path.exists(processed_file) or overwrite:
            out = []

            # cdnas = pd.read_table(r"C:\Users\tobia\PycharmProjects\GENES\Genebank_Ensembl_GFF3\mRNA_cDNA_release_111.tsv",
                                  # low_memory=False)

            xpore_data = pd.read_table("../Nanopore/direct_RNA_seq/"
                                       "Xpore/DATA/"
                                       "xpore_cdna.all_group_compare_2025/"
                                       "diffmod_s3-s4/"
                                       "majority_direction_kmer_diffmod_C3_modRate_05_cDNA_UTR_REL_cDNA_POS.tsv")
            print(len(xpore_data))
            # xpore_data["id"] = xpore_data["id"].apply(lambda x: x.split(".")[0])
            """xpore_data = pd.merge(xpore_data, cdnas[["strand", "transcript_id"]], left_on="id", right_on="transcript_id",
                                  how="left")"""
            # xpore_data["kmer_rel_pos_in_transcript"] = xpore_data["position"] / xpore_data["len"]
            xpore_data = xpore_data[["id", "kmer", "position", "relative_cdna_pos", "z_score_s3_vs_s4", "z_score_s3_vs_s4"]]
            xpore_data = xpore_data.dropna(subset=["relative_cdna_pos"])

            print(len(xpore_data))
            # transcripts = xpore_data.id.unique()

            for file in glob.glob(f'{folder_out}/*/*/seq_params_best_10.csv'):
                temp = pd.read_csv(file)
                # print(file)
                transcript = file.split("\\")[-3]
                temp["transcript_ID"] = transcript
                # temp = temp.head(10)
                # motif_counts = temp[["kmer_in_structure", "kmer"]].groupby("kmer_in_structure").count().sort_values(by="kmer",
                #                                                                                                      ascending=False)
                # motif_ = motif_counts.index[0]
                # out.append(temp[temp.kmer_in_structure == motif_].head(1))
                out.append(temp.head(1))  # take first line with the lowest energy
            data = pd.concat(out)
            # tbl = tbl[tbl["mfe"] < - 20]
            print(data)
            # data = data.drop_duplicates(subset=["transcript_ID", "kmer", "pos", "kmer_in_structure"], keep='first')
            data = pd.merge(data, xpore_data,
                            left_on=["transcript_ID",
                                     "kmer" + tag, "pos" + tag
                                     ],
                            right_on=["id",
                                      "kmer", "position"
                                      ],
                            how="left")

            data = data.rename(columns={"kmer_x": "kmer"})
            # print(len(tbl))
            # print(tbl)
            # tbl["3mer_in_structure"] = tbl["kmer_in_structure"].apply(lambda x: x[1:-1])
            # tbl["kmer_in_structure"] = tbl["kmer_in_structure"].apply(lambda x: x[2])
            data["rel_pos_in_window"] = (data["pos"] - data["start"]) / (data["end"] - data["start"])
            data["kmer_in_structure"] = data[["pos", "start", "forgi_structure"]].apply(get_kmer_forgi_elem, axis=1)
            data.to_csv(processed_file, index=False)
        else:
            data = pd.read_csv(processed_file)

        data["3mer_in_structure"] = data["kmer_in_structure"].apply(lambda x: x[1:-1])
        data["middle_base_forgi_elem"] = data["kmer_in_structure"].apply(lambda x: x[2])
        data["middle_base"] = data["kmer"].apply(lambda x: x[2])

        data = data.drop_duplicates(subset=["transcript_ID", "position", "kmer_in_structure"])

        # data = data[data.z_score_s3_vs_s4 < 0]

        # print(data)
        # data_temp2 = data[data.kmer.str.contains("G[UC]U[UC][UC]",
        #                                                   flags=re.IGNORECASE, regex=True)]

        # data_temp2 = data.copy()
        # # print(data_temp2.sort_values(by="kmer"))
        # print(data[data.z_score_s3_vs_s4 < 0].groupby(["transcript_ID"]).count().sort_values(by="seq_window", ascending=False))

        #
        # print(a := data[(data.kmer_in_structure.str.contains("hhh") |
        #             data.kmer_in_structure.str.contains("iii"))].groupby(["transcript_ID"]).count().sort_values(by="seq_window", ascending=False))
        # for t in a[a["seq_window"] > 2].index:
        #     print(t)
        #
        # print(data[(data["id"] == "ENST00000225964")
        #            # & (data["kmer"] == "GUUUU")
        #            # & (data.kmer_in_structure.str.contains("hhh") |
        #            #    data.kmer_in_structure.str.contains("iii"))
        #            # & data["relative_cdna_pos"].between(2, 3)
        #            ])

        # print(data[(data["middle_base"] == "U")
        #            & (data["kmer"] == "GUUUU")
        #            & (data["kmer_in_structure"] == "iiiii")
        #            & data["relative_cdna_pos"].between(2, 3)
        #            ])  #  .groupby(["3mer_in_structure", "kmer"]).count()
        # print(a[a["seq_window"] > 5])

        sumstats_total_out = []
        motif_ratio_out = []
        motif_ratio_name = []
        RNA_element = {0: "FIVE_UTR", 1: "CDS", 2: "THREE_UTR", 3: "mRNA"}
        print(len(name_flag))
        for i in range(1, 3):
            for base in [# "All",
                         "A", "C", "G", "U"
                         ]:
                print(RNA_element[i], base)
                # tbl = tbl.drop_duplicates(subset=["transcript_ID", "kmer", "pos", "kmer_in_structure"])
                if len(name_flag) > 0:
                    if base != "U":
                        continue
                if i == 3:
                    data_temp = data.copy()
                else:
                    data_temp = data[data["relative_cdna_pos"].between(i, i + 1)]
                if base != "All":
                    data_temp = data_temp[data_temp.middle_base == base]
                if base == "U":
                    if len(name_flag) > 0:
                        print("filter for ", name_flag)
                        data_temp = data_temp[data_temp.kmer.str.contains(selected_motif,
                                                                          flags=re.IGNORECASE, regex=True)]
                    # data_temp = data_temp[data_temp.kmer == "GUUUU"]
                # print(data_temp)
                kmer_RNA_elem_counts = data_temp.groupby(col_name)["forgi_structure"].count().reset_index().sort_values(by="forgi_structure",
                                                                                                                   ascending=False).reset_index(drop=True)

                # print(kmer_RNA_elem_counts)
                if mirrow_motifs:
                    drop_rows = []
                    for i_idx in range(0, len(kmer_RNA_elem_counts) - 1):
                        for j in range(i_idx + 1, len(kmer_RNA_elem_counts)):
                            if kmer_RNA_elem_counts.iloc[i_idx, 0] == kmer_RNA_elem_counts.iloc[j, 0][::-1]:
                                # print(kmer_RNA_elem_counts.iloc[i, 0],  kmer_RNA_elem_counts.iloc[j, 0])
                                kmer_RNA_elem_counts.iloc[i_idx, 1] = kmer_RNA_elem_counts.iloc[i_idx, 1] + kmer_RNA_elem_counts.iloc[j, 1]
                                drop_rows.append(j)
                    kmer_RNA_elem_counts = kmer_RNA_elem_counts.drop(drop_rows).sort_values(by="forgi_structure", ascending=False)

                # kmer_RNA_elem_counts["forgi_structure"] = kmer_RNA_elem_counts["forgi_structure"] / kmer_RNA_elem_counts["forgi_structure"].sum() * 100

                numbers = [str(x) for x in range(0, 10)]
                rna_elements = kmer_RNA_elem_counts[col_name].to_list()

                data_temp["forgi_structure_str_only"] = data_temp["forgi_structure"].apply(lambda x: "".join([elem_
                                                                                                    for elem_
                                                                                                    in x if elem_ not in numbers]))
                #data_temp["len_forgi_structure1"] = data_temp["forgi_structure_str_only"].apply(lambda x: len(x))

                data_temp["forgi_structure_str_only"] = data_temp["forgi_structure_str_only"].apply(lambda x: x[10:-10])

                # data_temp["len_forgi_structure2"] = data_temp["forgi_structure_str_only"].apply(lambda x: len(x))
                # print(data_temp)

                """for elem in rna_elements:
                    tbl = tbl[tbl["forgi_structure_str_only"].str.contains(elem)]"""

                if overlapping_counts:
                    data_temp[rna_elements] = data_temp.apply(lambda x: [len([m.start() for m in re.finditer('(?={0})'.format(re.escape(e)),
                                                                                                   x["forgi_structure_str_only"])
                                                                    ])
                                                               for e in rna_elements],
                                                    axis=1,
                                                    result_type='expand')
                    # data_temp[rna_elements] = data_temp.apply(lambda x: x[rna_elements] / x[rna_elements].sum() * 100,
                    #                                 axis=1,
                    #                                 result_type='expand')

                else:
                    # count returns non-overlapping counts !
                    data_temp[rna_elements] = data_temp.apply(lambda x: [x["forgi_structure_str_only"].count(e) for e in rna_elements],
                                                    axis=1,
                                                    result_type='expand')

                    # data_temp[rna_elements] = data_temp.apply(lambda x: x[rna_elements] / x[rna_elements].sum() * 100,
                    #                                 axis=1,
                    #                                 result_type='expand')

                ###################################################################################################
                # Perform enrichment analysis
                ###################################################################################################
                """
                foreground_total = Anzahl aller RNA element (5m 5i, mmmssm) für middle base
                background_total = 
                """

                motif_counts_total = data_temp[rna_elements].sum().sort_values()
                motif_counts_total = motif_counts_total / motif_counts_total.sum() * 100
                motif_ratio_out.append(motif_counts_total)
                motif_ratio_name.append(RNA_element[i] + "_" + base)

                foreground_total = kmer_RNA_elem_counts["forgi_structure"].sum()
                background_total = data_temp[rna_elements].sum().sum()
                p_vals, p_adj_vals, x_folds = [], [], []
                p_vals_f, or_f = [], []
                stat_vals = []
                # print(data_temp)
                print("PERFORM ENRICHMENT")
                for motif, motif_count in zip(kmer_RNA_elem_counts[col_name], kmer_RNA_elem_counts["forgi_structure"]):
                    # print(motif, motif_count)
                    # ( data_temp[motif])

                    kmer_background = data_temp[motif].sum()
                    kmer_foreground = motif_count

                    M = background_total  # counts of all motifs in seq in all transcripts with >= 1 DMR
                    n = kmer_background  # counts of analyzed motif in seq in all transcripts with >= 1 DMR
                    N = foreground_total  # counts of motifs containing DMRs
                    x = kmer_foreground  # counts of analyzed motifs containing DMRs
                    stat_vals.append([M, n, N, x])
                    # Calculate p-value (survival function = P(X >= x))
                    p_value = hypergeom.sf(x - 1, M, n, N)

                    table = [
                        [x, N - x],
                        [n - x, M - N - n + x]
                    ]

                    odds_ratio, p_value_f = fisher_exact(table)
                    p_vals_f.append(p_value_f)
                    or_f.append(odds_ratio)
                    # Calculate X-Fold over-representation
                    expected_count = (n / M) * N
                    x_fold = (x / expected_count) if expected_count > 0 else float('inf')

                    # t_statistic, p_value = stats.ttest_1samp(data_temp[motif], motif_count)
                    # print("\n")
                    # print(motif)
                    # global_mean = data_temp[motif].mean()
                    # mean_motif_rations.append(global_mean)
                    # print("global ratio vs dmr ratio : ", global_mean, "vs.", motif_count)
                    # print(t_statistic, p_value)
                    # p_adj = p_value * len(data_temp) * window * len(kmer_RNA_elem_counts)
                    # print("p_adj: ", p_adj)
                    # print("x_folds: ",  motif_count / global_mean)
                    p_vals.append(p_value)
                    x_folds.append(x_fold)
                stat_vals = np.asarray(stat_vals)
                sum_stats = pd.DataFrame(data={"motif": kmer_RNA_elem_counts[col_name].values,
                                               # "mean_ratio_dmr_motif": kmer_RNA_elem_counts["forgi_structure"].values,
                                               # "mean_ratio_motif": mean_motif_rations,
                                               "p_val": p_vals,
                                               "x_fold": x_folds,
                                               "p_vals_fisher": p_vals_f,
                                               "OR_fisher": or_f,
                                               "total_num_motifs_M": stat_vals[:, 0],
                                               "total_num_analyzed_motif_n": stat_vals[:, 1],
                                               "num_all_motifs_with_DMR_N": stat_vals[:, 2],
                                               "num_analyzed_motif_with_DMR_x": stat_vals[:, 3],
                                               })

                _, adj_val = multi.fdrcorrection(sum_stats['p_val'])
                sum_stats["q_val"] = adj_val

                _, adj_val = multi.fdrcorrection(sum_stats['p_vals_fisher'])
                sum_stats["q_val_fisher"] = adj_val


                sum_stats = sum_stats[["motif", "p_val", "q_val", "x_fold", "p_vals_fisher",
                                       "q_val_fisher", "OR_fisher",
                                       "total_num_motifs_M",
                                       "total_num_analyzed_motif_n",
                                       "num_all_motifs_with_DMR_N",
                                       "num_analyzed_motif_with_DMR_x"
                                       ]]

                sum_stats["Base"] = base
                sum_stats["mRNA_region"] = RNA_element[i]

                sumstats_total_out.append(sum_stats)
                # sum_stats = sum_stats.sort_values(by="x_fold", ascending=False)
                # print(sum_stats[sum_stats.q_val <= 0.05])
                # print("\n")
                # file_out = f'RNAfolds/{window}bp_window/0_sumstats/sum_stats_kmer3_as_forgi_{base}_{RNA_element[i]}.csv'
                # sum_stats.to_csv(file_out, index=False)
                # print(sum_stats[(sum_stats["mean_ratio_dmr_motif"] > 3) & (sum_stats.p_adj_val < 0.05)].sort_values(by="odds_ratio", ascending=False))

        motif_counts_df = pd.concat(motif_ratio_out, axis=1, names=motif_ratio_name)
        motif_counts_df = motif_counts_df.rename(columns=dict(zip(range(12), motif_ratio_name)))
        motif_counts_df = motif_counts_df.fillna(0)
        motif_counts_df["mean_ratio"] = motif_counts_df.mean(axis=1)
        motif_counts_df = motif_counts_df.sort_values(by="mean_ratio", ascending=False)

        # print(motif_counts_df)
        # print(motif_counts_df.sum(axis=0))

        sumstats_total_out = pd.concat(sumstats_total_out)
        _, adj_val = multi.fdrcorrection(sumstats_total_out['p_val'])
        sumstats_total_out["q_val"] = adj_val
        sumstats_total_out = sumstats_total_out.sort_values(by=["mRNA_region", "Base", "x_fold"], ascending=False)
        file_out = (f'{folder_out}/0_sumstats/'
                    f'sum_stats_kmer_as_forgi_{col_name}_overlap_counts_{overlapping_counts}_{name_flag}{name_flag2}.csv')

        sumstats_total_out = sumstats_total_out.sort_values(by="p_vals_fisher", ascending=True)
        sumstats_total_out.to_csv(file_out, index=False)

        motif_counts_df.to_csv(file_out.replace(".csv", "_motif_ratio.csv"), index=True)

        print(sumstats_total_out[(sumstats_total_out.q_val_fisher <= 0.05)
                                 & (sumstats_total_out.motif.isin(motif_counts_df[motif_counts_df.mean_ratio > 0].index))])
        # print(sum_stats.sort_values(by="odds_ratio", ascending=False))

    """
            -------------------------------------------------
            |                  M                            |
            |          ______________                       |
            |          |             |                      |
            |          |     N    ---|-----------           |
            |          |         | x |          |           |
            |           _________|___|   n      |           |
            |                    |              |           |
            |                    |              |           |
            |                    ----------------           |
            -------------------------------------------------
            
            N = total n of positions / total number of background proteins
            M = n of all pos for RNA motif / number of protein in pathway
            n = n of all DMRs / DEGs
            x = n of DMRs with RNA motif
    """

