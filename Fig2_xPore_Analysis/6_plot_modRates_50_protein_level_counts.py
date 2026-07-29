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
from adjustText import adjust_text
from scipy.signal import find_peaks
from mpl_toolkits.mplot3d import Axes3D
from sklearn.linear_model import LinearRegression

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def plot_modbase_ratio():
    # plot distributions
    plt.style.use('ggplot')
    plt.figure(figsize=(8, 4))
    print(data.columns)
    row_idx = 0
    col_idx = 0
    for base in ["A", "C", "G", "U"]:
        data["dmr_count_ratio_" + base] = data[f"dmr_{base}_count"] / data["transcript_length"] * 100
        if row_idx == 2:
            row_idx = 0
            col_idx += 1
        plt.subplot2grid((3, 2), (row_idx, col_idx))
        # idx += 1
        row_idx += 1

        sns.kdeplot(data.loc[data.is_DMR & (data[f"dmr_{base}_count"] > 0), "dmr_count_ratio_" + base],
                    color=colors[base], cut=0, fill=True,
                    alpha=0.75, label=base,
                    # label=f"{num_modbase_mrnas} dominant isoforms\nwith highest number of ∆-bases"
                    )
        plt.text(x=0.8, y=0.95, s=base,
                 transform=plt.gca().transAxes, color=colors[base], fontsize=10, va="top",
                 ha="right",
                 fontdict={"weight": "bold"})
        plt.xlim(0, 3)
        plt.xlabel("Ratio (%)")

    plt.subplot2grid((3, 2), (row_idx, 0), colspan=2)

    sns.kdeplot(data.loc[data.is_DMR, "dmr_count_ratio"],
                color="gray", cut=0, fill=True,
                alpha=0.75, legend=False,
                # label=f"{num_modbase_mrnas} dominant isoforms\nwith highest number of ∆-bases"
                )
    plt.text(x=0.8, y=0.95, s="All modified bases",
             transform=plt.gca().transAxes,
             color="gray",
             fontsize=10, va="top",
             ha="right",
             fontdict={"weight": "bold"})

    plt.xlabel("Ratio (%)")
    plt.xlim(0, 6)
    # plt.legend()
    # plt.tight_layout()
    plt.show()


def plot_corr_matrix():
    corr_features = ["dmr_count",
                     "dmr_count_ratio",
                     f"mrna_ont_{cond}",
                     "transcript_length",
                     "fiveUTR_len",
                     "threeUTR_len",
                     f"polya_length_{cond}",
                     # f"protein_level_{cond}",
                     # "protein_per_mRNA"
                     ]

    corr_matrix = data[corr_features].corr(method="spearman")  # kendall or spearman or pearson

    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    plt.figure(2, figsize=(7, 5))
    print(corr_matrix)
    sns.heatmap(corr_matrix, cmap="coolwarm", vmax=0.6, vmin=-0.6, center=0, mask=mask,
                cbar_kws={"shrink": 0.6, "location": "top"})

    xticks = plt.gca().get_xticks()
    xticks_label = features  # plt.gca().get_xticklabels()
    # xticks_label = ["polya_length", "DMR count", "mRNA level", "mRNA length", "5'-UTR length", "3'-UTR length"]
    plt.xticks(xticks[:-1], xticks_label[:-1], rotation=45, ha="right")
    yticks = plt.gca().get_yticks()
    yticks_label = features  # plt.gca().get_yticklabels()
    # yticks_label = ["mRNA level", "mRNA length", "5'-UTR length", "3'-UTR length", "Poly(A) length", "protein_level"]
    plt.yticks(yticks[1:], yticks_label[1:],  # rotation=-45, va="bottom"
               )
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    bases = ["A", "U", "C", "G"]
    colors = dict(zip(bases, ["#109648", "#D62839", "#255C99", "#F7B32B"]))

    cond = "C3"
    data = pd.read_table("direct_RNA_seq/0_TABLES/DMR_MASTER_TABLE_Transcript_Level_C3_modRate_05_PAPER.tsv.tsv")
    data["dmr_count_ratio"] = data["dmr_count"] / data["transcript_length"] * 100

    features = ["Number of ∆-bases",
                "Ratio of ∆-bases (%)",
                "mRNA level (log$_2$(TPM + 1))",
                "mRNA length",
                "5'-UTR length",
                "3'-UTR length",
                "Poly(A) length",
                # "Protein level",
                # "Protein/mRNA"
                ]

    old_names = ["dmr_count",
                 "dmr_count_ratio",
                 f"mrna_ont_{cond}",
                 "transcript_length",
                 "fiveUTR_len",
                 "threeUTR_len",
                 f"polya_length_{cond}",
                 # f"protein_level_{cond}",
                 # "protein_per_mRNA"
                 ]
    name_dict = dict(zip(old_names, features))

    data["mrna_ont_C3"] = data["mrna_ont_C3"].fillna(0)
    data["mrna_ont_C4"] = data["mrna_ont_C4"].fillna(0)
    dmr_count_thresh = 1
    data["is_DMR"] = False
    data.loc[data["dmr_count"] >= dmr_count_thresh, "is_DMR"] = True

    # plot_modbase_ratio()
    plot_corr_matrix()

    quit()
