import pandas as pd
import numpy as np
import json
from scipy.stats import hypergeom
from scipy.stats import fisher_exact
import statsmodels.stats.multitest as multi
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
import re
import os


def count_all_possible_kmers():

    if os.path.exists(f"kmer_count_dict{tag}.json"):
        return
    filter_name = ("direct_RNA_seq/Xpore/DATA/"
                   "xpore_cdna.all_group_compare_2025/diffmod_s3-s4/"
                   "majority_direction_kmer_diffmod_padj_modrate_025_GENOME_cDNA_UTR.tsv")
    xpore_data = pd.read_table(filter_name)
    if tag == "_C4_vs_C3_pos_higher":
        xpore_data = xpore_data[(xpore_data.z_score_s3_vs_s4 < 0) & (xpore_data.mod_assignment == "higher")]
    if tag == "_C4_vs_C3_pos_lower":
        xpore_data = xpore_data[(xpore_data.z_score_s3_vs_s4 < 0) & (xpore_data.mod_assignment == "lower")]
    if tag == "_C4_vs_C3_neg_higher":
        xpore_data = xpore_data[(xpore_data.z_score_s3_vs_s4 > 0) & (xpore_data.mod_assignment == "higher")]
    if tag == "_C4_vs_C3_neg_lower":
        xpore_data = xpore_data[(xpore_data.z_score_s3_vs_s4 > 0) & (xpore_data.mod_assignment == "lower")]

    xpore_data.kmer = xpore_data.kmer.apply(lambda x: x.replace("T", "U"))
    xpore_data = xpore_data.dropna(subset=['cdna_seq'])
    # Create all combinations
    combinations = itertools.product('ACUG', repeat=5)
    kmers = np.array([''.join(p) for p in combinations])
    print(len(kmers))  # should print 1024
    # quit()
    kmer_dict = dict()
    for idx, kmer in enumerate(kmers):
        found = 0
        for transcript, tbl in xpore_data.groupby("id"):
            seq = tbl.cdna_seq.values[0].replace("T", "U")
            # print(seq)
            found = found + len([m.start() for m in re.finditer('(?={0})'.format(re.escape(kmer)), seq)])
        dmr_kmer = len(xpore_data[xpore_data.kmer == kmer])
        print(idx, kmer, dmr_kmer, found)
        kmer_dict[kmer] = (dmr_kmer, found)

    with open(f"kmer_count_dict{tag}.json", "w") as f:
        json.dump(kmer_dict, f)


def run_enrichment():

    if os.path.exists(f'kmer_dmr_enrichment{tag}_v2.tsv'):
        return
    # filter_name = "../Xpore/DATA/diffmod_s1_s2_s4_vs_s3__NEW_TR_filtered_coverage15_na_filter_pre_proc_FINAL.tsv"
    # xpore_data = xpore_data[xpore_data["col_idx"] == "s4"]
    filter_name = r"direct_RNA_seq/Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4\majority_direction_kmer_diffmod.table"
    xpore_data = pd.read_table(filter_name, usecols=[2, 24])
    xpore_data.kmer = xpore_data.kmer.apply(lambda x: x.replace("T", "U"))
    xpore_data = xpore_data.drop_duplicates()
    print(xpore_data.head(5))

    # 1024 kmer combinations > all

    # Analyze only kmers that were measured by xPore
    xpore_kmers = xpore_data.kmer.unique()
    print(len(xpore_kmers))
    # xpore_kmers = [x.replace("T", "U") for x in xpore_kmers]
    # GUUUU im ungefiltert Xpore
    # GUUUM DMRs in gefiltert xpore_data

    # Alle positionsn in  ungefiltert xpore_data
    # Alle DMRs  xpore_data

    # Load the dictionary
    with open(f"kmer_count_dict{tag}.json", 'r') as f:
        kmer_dict = json.load(f)
    print(len(kmer_dict.keys()))

    keys_to_remove = [k for k in kmer_dict if k not in xpore_kmers]
    for k in keys_to_remove:
        kmer_dict.pop(k)

    print(len(kmer_dict.keys()))
    counts = sorted([v[1] for v in kmer_dict.values()])
    print(dict(sorted(kmer_dict.items(), key=lambda item: item[1][1])))

    # plt.bar(np.arange(len(counts)), counts)
    # plt.show()

    # Parameters
    kmers, p, x_fold_ = [], [], []
    # Calculate totals
    foreground_total = sum([v[0] for v in kmer_dict.values()])
    background_total = sum([v[1] for v in kmer_dict.values()])
    kmer_foreground_ = []
    kmer_background_ = []
    for kmer in kmer_dict.keys():
        # kmer = 'GUUUU'

        # Extract counts
        kmer_foreground, kmer_background = kmer_dict[kmer]
        kmer_foreground_.append(kmer_foreground)
        kmer_background_.append(kmer_background)

        # Set hypergeometric parameters
        M = background_total        # total population size
        n = kmer_background          # total "success" states in the population
        N = foreground_total         # total number of draws
        x = kmer_foreground          # observed number of successes

        # Calculate p-value (survival function = P(X >= x))
        p_value = hypergeom.sf(x-1, M, n, N)

        neither = M - N - n + x
        # Equivalent one-sided Fisher exact test
        table = [
            [x, n - x],
            [N - x, neither]
        ]

        odds_ratio, p_value_fisher = fisher_exact(
            table,
            alternative="greater"
        )
        # Calculate X-Fold over-representation
        expected_overlap = n * N / M
        fold_enrichment = (x / n) / (N / M)

        # Calculate X-Fold over-representation
        # expected_count = (n / M) * N
        # x_fold = (x / expected_count) if expected_count > 0 else float('inf')
        kmers.append(kmer)
        p.append(p_value_fisher)
        x_fold_.append(fold_enrichment)

        # print(f"p-value for enrichment of {kmer}: {p_value:.4e}")
        # print(f"X-Fold over-representation of {kmer}: {x_fold:.2f}")

    df = pd.DataFrame({'kmer': kmers, 'pval': p, 'x_fold': x_fold_,
                       "background_total": background_total,
                       "foreground_total": foreground_total,
                       "kmer_background": kmer_background_,
                       "kmer_foreground": kmer_foreground_
                       })

    _, adj_val = multi.fdrcorrection(df['pval'])
    df["qval"] = adj_val
    df = df[["kmer", "pval", "qval", "x_fold",
             "background_total", "foreground_total",
             "kmer_background", "kmer_foreground"
             ]]
    df = pd.merge(df, xpore_data, how='left', on="kmer")
    df.to_csv(f'kmer_dmr_enrichment{tag}_v2.tsv', sep="\t", index=False)

    df = df[df.qval < 0.05]
    print(df.sort_values('x_fold', ascending=False))


def plot_enrichment():
    df = pd.read_table(f'kmer_dmr_enrichment{tag}_v2.tsv')
    df = df[df.pval < 0.05]
    # df = df[df.x_fold > 5]
    # df = df[df.kmer_foregr    ound > 5]
    df["color"] = "#909090"
    df.loc[df.qval < 0.05, "color"] = "#555555"
    df = df.nlargest(10, 'kmer_foreground')
    df = df.sort_values('x_fold', ascending=False)
    df.kmer = df.kmer.apply(lambda x: x[:2] + "-" + x[2] + "-" + x[3:])
    # df = df.sort_values(by="mod_assignment")
    plt.style.use('ggplot')
    plt.figure(figsize=(3, 4))
    plt.subplot(1, 2, 1)
    plt.subplots_adjust(left=0.45, bottom=0.15, top=0.95, right=0.86, wspace=0.3)
    sns.barplot(x="x_fold", y="kmer", # hue="mod_assignment",
                data=df, palette=df["color"].to_list()
                # palette={"lower": "#4586AC", "higher": "#CC5A49"},
                )

    xticks = plt.gca().xaxis.get_major_ticks()
    xticks[0].set_visible(False)
    plt.grid(False)

    # plt.legend(
    #     title="xPore mode",
    #     title_fontsize='small',
    #     fontsize='small',
    #     frameon=False,
    #     facecolor='white',
    #     framealpha=1,
    #     labelspacing=0.35,
    #     handletextpad=0.2,
    #     # bbox_to_anchor=(1.05, 1),
    #     # loc=2, borderaxespad=0.,
    # )

    plt.legend().set_visible(False)
    plt.ylabel("")
    plt.xlabel("OR", fontsize=13)
    # plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)

    plt.subplot(1, 2, 2)
    # df["kmer_background"] = df["kmer_background"] / df["background_total"] * 100
    df["kmer_foreground"] = df["kmer_foreground"] / num_sig_kmers * 100

    sns.barplot(x="kmer_foreground", y="kmer",
                data=df,
                palette=df["color"].to_list()
                )
    plt.axvline(x=100/1024, color='black', alpha=0.75, linestyle='dotted', linewidth=1.5)
    plt.ylabel("")
    plt.xlabel("%", fontsize=13)
    xticks = plt.gca().xaxis.get_major_ticks()
    xticks[0].set_visible(False)

    # plt.xticks(fontsize=13)
    plt.yticks([])
    plt.grid(False)

    plt.savefig(f'kmer_dmr_enrichment{tag}_v2.svg')
    plt.show()


if __name__ == '__main__':
    pass
    filter_name = ("direct_RNA_seq/Xpore/DATA/"
                   "xpore_cdna.all_group_compare_2025/diffmod_s3-s4/"
                   "majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS.tsv")
    xpore_data = pd.read_table(filter_name)
    print(len(xpore_data))
    num_sig_kmers = len(xpore_data)

    for mode in ["higher", "lower"]:
        for z in ["neg", "pos"]:
            tag = f"_C4_vs_C3_{z}_{mode}"
            # count_all_possible_kmers()  # C4 vs. C3
            # run_enrichment()
            plot_enrichment()
