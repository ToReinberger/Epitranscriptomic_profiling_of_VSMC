import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text

pd.set_option('display.max_rows', 200)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)

bases, colors, zorder = ["A", "C", "G", "U"], ["#109648", "#255C99", "#F7B32B", "#D62839"], [1, 1, 1, 1]

data = pd.read_table(
    r"direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS_GENENAMES.tsv")
data["id"] = data["id"].apply(lambda x: x.split(".")[0])
print(data)

data["diff_mod_rate_s4_vs_s3"] = data["diff_mod_rate_s3_vs_s4"] * -1
data["logP"] = np.log10(data["adj_pval_s3_vs_s4"]) * -1
# data["diff_mod_rate_s4_vs_s3"] = data["diff_mod_rate_s3_vs_s4"] * -1

plt.style.use('ggplot')
# data = data[data["middle_base"] == "U"]
# data = data[data["kmer"] == "GUUUU"]




row_idx = 0
col_idx = 0
for base, color in zip(["A", "C", "G", "U"], colors):
    y_thresh = 10
    if base == "G":
        continue
    if base == "U":
        col_span = 2
        y_thresh = 15
    else:
        col_span = 1
    print(row_idx, col_idx)
    if col_idx == 2:
        col_idx = 0
        row_idx += 1

    # plt.subplot2grid((2, 2),  (row_idx, col_idx), colspan=col_span)
    col_idx += 1

    data_temp = data[data["middle_base"] == base]
    max_y = data_temp.logP.max()
    plt.figure(figsize=(14, 8))
    sns.scatterplot(data=data_temp,
                    x="diff_mod_rate_s4_vs_s3",
                    y="logP",
                    color="gray")
    sns.scatterplot(data=data_temp[(abs(data_temp["diff_mod_rate_s4_vs_s3"]) >= 0.25)
                              & (data_temp["adj_pval_s3_vs_s4"] <= 0.05)],
                    x="diff_mod_rate_s4_vs_s3", y="logP",
                    color=color)
    plt.xlim(-0.9, 0.9)
    plt.ylim(2.5, max_y * 1.05)
    texts = []
    for x, y, s, k in zip(data_temp["diff_mod_rate_s4_vs_s3"], data_temp["logP"],
                          data_temp["Gene name"], data_temp["kmer"]):
        if (abs(x) > 0.25 and y > y_thresh) or y > 5 or abs(x) > 0.6:
            texts.append(plt.text(x, y + (max_y / 50), f"{s}, {k}", ha="center", va="center", color="k", fontsize=7))
    if True:
        adjust_text(texts)
    plt.text(x=0.92, y=0.92, s=base, fontsize=22, ha="center", va="center", color=color, fontweight="bold",
             transform=plt.gca().transAxes)

    plt.xlabel("Difference in modification rate")
    plt.ylabel("-log(p-value)")
    plt.title("C4 vs. C3")
    plt.savefig(f"direct_RNA_seq/0_PLOTS/XPORE_Vulcano_plots/Vulcano_plot_{base}.svg")
    plt.show()
