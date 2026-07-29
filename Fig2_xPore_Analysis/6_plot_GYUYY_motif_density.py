import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)

# data = pd.read_table("majority_direction_kmer_diffmod_padj_GENE_NAMES.tsv")
# data = data[data["adj_pval_s3_vs_s4"] <= 0.05]

data = pd.read_table("direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS_GENENAMES.tsv")
data = data[data["adj_pval_s3_vs_s4"] <= 0.05]
# data = data[data.adj_pval_s3_vs_s4 <= 0.05]
data = data[abs(data.diff_mod_rate_s3_vs_s4) >= 0.25]

selected_motif = "G[CU]U[UC][UC]"
print(len(data["id"].unique()))

data = data[data.kmer.str.contains(selected_motif, flags=re.IGNORECASE, regex=True)]

print(len(data["id"].unique()))

s3 = [x for x in data.columns if "mod_rate_s3" in x]
s4 = [x for x in data.columns if "mod_rate_s4" in x]
data = data[data[s3].mean(axis=1).between(0, 1)]
data = data[data[s4].mean(axis=1).between(0, 1)]


# plt.pie([25.9, 100 - 25.9])
# plt.show()
# quit()
plt.style.use('ggplot')
plt.figure(figsize=(3, 3))
sns.kdeplot(data[s3].mean(axis=1), label="C3", color="#8DD3C7", fill=True, alpha=0.8, cut=0)
sns.kdeplot(data[s4].mean(axis=1),  label="C4", color="#FDB462", fill=True, alpha=0.8, cut=0)
plt.xlabel("Modification rate")
plt.legend(loc="upper right")
plt.title("GUU-U/C-U/C")
plt.tight_layout()
plt.show()
quit()
# data = data.dropna(subset="transcript_name")
# data = data[data.adj_pval_s3_vs_s4 <= 0.05]
# data = data[abs(data.diff_mod_rate_s3_vs_s4) >= 0.25]
print(len(data[data.transcript_name.str.startswith("FN1-")]))
print(data[data.transcript_name.str.startswith("FN1-")])


data = pd.read_table("majority_direction_kmer_diffmod_before_padj_GENE_NAMES.tsv")
data = data.dropna(subset="transcript_name")
print(len(data[data.transcript_name.str.startswith("FN1-")]))
print(data[data.transcript_name.str.startswith("FN1-")])

quit()
print(len(data))
print(data.columns)
print(data.adj_pval_s3_vs_s4.max())
quit()

data["kmer"] = data["kmer"].apply(lambda x: x.replace("T", "U"))

# a = data.groupby("kmer").count().sort_values("id", ascending=False)
# print(a[a.id > 100])
# quit()

selected_motif = "GUU[UC][UC]"
# selected_motif = "G[U]U[U][UC]"
# selected_motif = "GUUUU"
print(a := len(data.loc[ data["adj_pval_s3_vs_s4"] < 0.05, "id"].unique()))
data = data[data.kmer.str.contains(selected_motif, flags=re.IGNORECASE, regex=True)]
print(b := len(data.loc[ data["adj_pval_s3_vs_s4"] < 0.05, "id"].unique()))
print(b /(a + b) * 100
      )
print(data)

plt.pie([15.8, 100 - 15.8])
plt.show()
quit()
# data = data[data.kmer == "GUUUU"]
# data = data.dropna(subset="transcript_name")
# data["genes"] = data["transcript_name"].apply(lambda x: x.split("-")[0])
#
# print(a := data.groupby("genes").count().sort_values("id", ascending=False))
# print()
#
# for g in a[a["id"] > 1].index:
#     print(g)

s3 = [x for x in data.columns if "mod_rate_s3" in x]
s4 = [x for x in data.columns if "mod_rate_s4" in x]
data = data[data[s3].mean(axis=1).between(0, 1)]
data = data[data[s4].mean(axis=1).between(0, 1)]

plt.style.use('ggplot')
sns.kdeplot(data[s3].mean(axis=1), label="C3", color="#8DD3C7", fill=True, alpha=0.8, cut=0)
sns.kdeplot(data[s4].mean(axis=1),  label="C4", color="#FDB462", fill=True, alpha=0.8, cut=0)
plt.xlabel("Modification rate")
plt.legend(loc="upper right")
plt.title("GUU-U/C-U/C")
plt.tight_layout()
plt.show()

# print(data[data["id"].str.contains("ENST00000640411")])
quit()

data["five_len"] = abs(data["cDNA_end_fiveUTR"] - data["cDNA_start_fiveUTR"])
data["three_len"] = abs(data["cDNA_end_threeUTR"] - data["cDNA_start_threeUTR"])

print("FIVE UTR")
print("MIN", data.five_len.min())
print("MEDIAN", data.five_len.median())
print("MAX", data.five_len.max())

print("\nTHREE UTR")
print("MIN",data.three_len.min())
print("MEDIAN", data.three_len.median())
print("MAX", data.three_len.max())

quit()
data = pd.read_table("majority_direction_kmer_diffmod.table")

print(len(data))
print(len(data["id"].unique()))
