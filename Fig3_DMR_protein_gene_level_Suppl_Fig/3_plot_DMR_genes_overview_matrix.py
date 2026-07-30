import os

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from scipy import stats
import matplotlib as mpl

# mpl.rcParams['axes.spines.left'] = False
# mpl.rcParams['axes.spines.right'] = False
# mpl.rcParams['axes.spines.top'] = False
# mpl.rcParams['axes.spines.bottom'] = False


def calc_ci(data):

    if len(data) == 0:
        return [np.nan]
    return (min(data), max(data))

    if (len(data) == 1):
        return [data[0]]

    if len(data) < 4:
        return (min(data), max(data))
    else:
        return (np.percentile(data, 25), np.percentile(data, 75))

    # mean = np.mean(data)
    # sem = stats.sem(data)               # standard error
    # ci = stats.t.interval(0.95, len(data)-1, loc=mean, scale=sem)
    #
    # return ci


def fill_dict(mod_=["higher"], val="relative_cdna_pos"):
    gene_dict = dict()
    for gene in data.gene_name:
        base_dict = dict()
        print(gene)
        xpore_temp = xpore[xpore["mod_assignment"].isin(mod_)]
        print(temp := xpore_temp.loc[xpore["Gene name"] == gene,
        ["relative_cdna_pos", "z_score_s4_vs_s3",
         "middle_base", "mod_assignment"]].sort_values(by="mod_assignment"))
        for base in bases:
            tbl2 = temp[temp.middle_base == base]
            ci = calc_ci(tbl2[val].values)

            mean_ = tbl2[val].mean()
            # out = [round(x, 2) for x in out if x is not np.nan]
            print(base, ci, mean_)
            base_dict[base] = (mean_, ci, tbl2[val].values, tbl2["mod_assignment"])
        gene_dict[gene] = base_dict
    return gene_dict


def plot_dmr_distr(xpore_dict, plot_pos=True, fig_title=""):
    genes = xpore_dict.keys()
    shift = -.3

    if plot_pos:
        plt.xlim(-0.2, 3.2)
        plt.xticks([])
        plt.text(x=0.5, y=len(genes) + .01, s="5'-UTR", fontweight="regular", color="k", fontsize=8, ha="center")
        plt.text(x=1.5, y=len(genes) + .01, s="CDS", fontweight="regular", color="k", fontsize=8, ha="center")
        plt.text(x=2.5, y=len(genes) + .01, s="3'-UTR", fontweight="regular", color="k", fontsize=8, ha="center")
        plt.axvline(x=1, linestyle="dotted", color="#555555", zorder=2)
        plt.axvline(x=2, linestyle="dotted", color="#555555", zorder=2)
    else:
        plt.axvline(x=0, linestyle="dotted", color="#555555", zorder=2)
        plt.xticks(fontsize=8)
    for i in range(0, len(genes)):
        if i % 2 == 0:
            continue
        plt.axhspan(i - 0.5, i + 0.5, color="#555555", alpha=0.1, zorder=1)

    for base in bases:
        median_temp, vals_temp, xpore_mod_ = [], [], []
        y_idx = 0
        for gene, elem in xpore_dict.items():
            median_temp.append(elem[base][0])
            vals_temp.append(elem[base][2])
            xpore_mod_.append(elem[base][3])
            ci = elem[base][1]
            if len(ci) > 1:
                print(ci)
                print("plot line", y_idx + shift)
                plt.hlines(y=y_idx + shift, xmin=ci[0], xmax=ci[1], color=colors[base], zorder=3,
                           linewidth=1.2)
            y_idx += 1
        print(vals_temp)
        y_temp = []
        x_temp = []
        for idx, e1 in enumerate(vals_temp):
            for e2 in e1:
                y_temp.append(idx)
                x_temp.append(e2)
        xpore_mod_temp = []
        for idx, e1 in enumerate(xpore_mod_):
            for e2 in e1:
                if e2 == "higher":
                    xpore_mod_temp.append("^")
                else:
                    xpore_mod_temp.append("v")
        y_temp = np.array(y_temp) + shift
        for xi, yi, m in zip(x_temp, y_temp, xpore_mod_temp):
            plt.scatter(xi, yi, marker=m, color=colors[base], s=50,
                        edgecolors="white", linewidths=0.1,
                        zorder=4)
        # plt.scatter(y=y_temp, x=x_temp, marker=xpore_mod_temp,color=colors[base], s=16)
        # plt.scatter(y=np.arange(0, len(genes)) + shift, x=median_temp, marker="d", color=colors[base], s=14)
        shift += 0.2

    plt.gca().yaxis.set_inverted(True)  # inverted axis with autoscaling
    plt.ylim(len(genes) - 0.5, -0.5)

    plt.yticks(np.arange(0, len(genes)), genes)
    plt.yticks([])
    plt.title(fig_title, rotation=45, fontsize=9, ha="left")

    return plt.gcf()


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)

data = pd.read_table("direct_RNA_seq/0_TABLES/DMR_MASTER_TABLE_Gene_Level_C3_vs_C4_padj_DMR_025_PAPER.tsv")
xpore = pd.read_table(r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore\direct_RNA_seq\Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4\majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS_GENENAMES.tsv")
xpore["z_score_s4_vs_s3"] = xpore["z_score_s3_vs_s4"] * -1
data = data[data.dmr_count > 1]
tcf21 = ['AQP5', 'NAT8', 'FERMT3', 'RNF166', 'RAB43', 'HSPBP1', 'STK3', 'VAV2', 'LRRC56', 'TMX4', 'TM9SF2', 'TGM2', 'PSMD8', 'CENPT', 'PKDCC', 'FAM3B', 'PDE2A', 'TGFA', 'CMTM1', 'TAS2R41', 'CAMKK1', 'SCAF8', 'WBSCR22', 'PHOX2A', 'PEAR1', 'ITPRIPL1', 'RUVBL2', 'F7', 'TECTB', 'ALDOA', 'AFG3L2', 'RECK', 'HIBCH', 'HRAS', 'ALG2', 'IL1RAP', 'TRAPPC6A', 'AICDA', 'LIG1', 'BRSK1', 'CA2', 'CNGA1', 'FETUB', 'PPP4C', 'HNRNPA3', 'POLR2M', 'SPSB1', 'DNAJC30', 'GYS1', 'PPP1R7', 'S100A3', 'ELMO3', 'PLA2G6', 'EIF4A2', 'HNRNPA2B1', 'DHODH', 'SLC38A7', 'CTU2', 'TSPAN33', 'SGSM1', 'CLDN22', 'DHRS7B', 'TAF6', 'NPRL2', 'TMEM119', 'SLC22A7', 'TIMM8B', 'AKR7A3', 'EHMT2', 'THAP11', 'PASK', 'DDX4', 'A1CF', 'CRYGN', 'SEC61B', 'ARL1', 'VSIG8', 'SKP1', 'KLF6', 'INS', 'MADCAM1', 'ARL9', 'MEOX1', 'BRWD1', 'CD3EAP', 'SLC51A', 'CDX2', 'KRTAP13-2', 'NID2', 'S100A4', 'TMEM95', 'ADCY7', 'DDAH2', 'RAB3C', 'KLK1', 'OTUD5', 'CBX3', 'ITIH4', 'PFKP', 'TSHR', 'PBX3', 'MRPL17', 'LHX6', 'TIPIN', 'RTF1', 'AKAP4', 'FCMR', 'MTF1', 'SCX', 'C4B', 'IL1RAPL1']
tcf21_chea = ['SCARB2', 'HSP90AB1', 'CLTC', 'FHL2', 'SLC2A1', 'ATP2A2', 'LAMC1', 'RND3', 'TSKU', 'FADS3', 'GPR176', 'C1QTNF1', 'EFEMP1', 'CSRP1', 'CREB3L1', 'CHMP1B', 'TIMP2', 'EMILIN1', 'PHLDA1', 'RAB8B', 'CAP1', 'SH3GLB1', 'IGFBP5', 'HEG1', 'IGFBP3', 'MYOF', 'ACSL3', 'UCHL3', 'SRP9', 'EEF1A1', 'ACLY', 'PTP4A2', 'MMP14', 'CCDC80', 'FRMD6', 'LOX', 'COL6A1', 'PXDN', 'MAPKAPK2', 'MXRA8', 'FSCN1', 'ADAM9', 'ANGPTL4', 'KCTD12', 'MXRA7', 'ERGIC1', 'PDXK', 'SHC1', 'LAMA4', 'MVP', 'LTBP2', 'FBLN2', 'THBS1', 'C3', 'CNN2', 'LIMA1', 'HSPH1', 'CALD1', 'S100A16', 'STC2', 'ESD', 'CDIPT', 'AP1M1', 'STAT1', 'SGTA', 'NEK7', 'FN1', 'PARVA', 'DHCR24', 'VASN', 'RAB10', 'COL1A1', 'COL3A1', 'COPS6', 'MYO1C', 'TNIP1', 'MAGED2', 'PKN1', 'P4HB', 'EIF3B', 'FBN1']
show_cols = list(data.columns[:10])
show_cols.append("TCF21")
show_cols.extend(["ont_gene_log2FC_C4vsC3", "protein_log2FC_C4vsC3"])
data["TCF21"] = 0
data["diff_gene_protein"] = abs(data["ont_gene_log2FC_C4vsC3"] - data["protein_log2FC_C4vsC3"])
data = data[data["diff_gene_protein"] > 0.25]

data.loc[data["gene_name"].isin(tcf21) | data["gene_name"].isin(tcf21_chea), "TCF21"] = 1
# data = data[data.TCF21 == 1]
data = data[data.dmr_count > 2]
data = data.sort_values(by=["protein_log2FC_C4vsC3"], ascending=False)
data["ont_gene_log2FC_C4vsC3"] = data["ont_gene_log2FC_C4vsC3"].fillna(0)

data = data[abs(data["protein_log2FC_C4vsC3"]) > 0.5]

xpore["relative_cdna_pos"] = xpore["relative_cdna_pos"].apply(lambda x: round(x, 2))
xpore["z_score_s4_vs_s3"] = xpore["z_score_s4_vs_s3"].apply(lambda x: round(x, 1))
# flip lower
# xpore.loc[xpore["mod_assignment"] == "lower", "z_score_s4_vs_s3"] = xpore.loc[xpore["mod_assignment"] == "lower", "z_score_s4_vs_s3"] * -1

rel_pos_dict = dict()

bases = ["A",
         "C",
         "G",
         "U"]



# plt.hist(data.loc[data["TCF21"] == 0, "dmr_count"], bins=20, alpha=0.5, label="Any")
# plt.hist(data.loc[data["TCF21"] == 1, "dmr_count"], bins=20, label="TCF21")
# plt.legend()
#
# plt.show()
print(show_cols)
old_cols = ['gene_name', 'transcript_count', 'dmr_transcript_count', 'dmr_count', 'dmr_A_count', 'dmr_C_count', 'dmr_G_count', 'dmr_U_count',
            "ont_gene_log2FC_C4vsC3", "protein_log2FC_C4vsC3"
            # 'TCF21'
            ]
new_cols = ['Gene', 'Isoforms', 'DMR+ isoforms', 'DMR count', 'A', 'C', 'G', 'U', 'log$_2$(FC) gene', 'log$_2$(FC) protein'

            # 'TCF21'
            ]


higher = fill_dict()
lower = fill_dict(mod_=["lower"])
both = fill_dict(mod_=["higher",
                       "lower"
                       ]
                 )
zscore = fill_dict(mod_=["higher",
                         "lower"
                         ]
                   , val="z_score_s4_vs_s3")


n_subplots = (1, 4)
data = data[old_cols].rename(columns=dict(zip(old_cols, new_cols)))
data = data.set_index("Gene")
cmap = mcolors.LinearSegmentedColormap.from_list("white_firebrick", ["white", "white", "firebrick"])
cmap2 = mcolors.LinearSegmentedColormap.from_list("blue_firebrick", ["steelblue", "white", "firebrick"])
plt.figure(figsize=(7, 15))
plt.subplot2grid(n_subplots, (0, 0))
sns.heatmap(data[['Isoforms', 'DMR+ isoforms', 'DMR count', 'A', 'C', 'G', 'U']],
            vmax=8, vmin=0, center=0, cmap=cmap, linewidths=0.05, square=False,
            annot=True, linecolor="lightgray", annot_kws={"size": 7},
            xticklabels=True, yticklabels=True,
            cbar=False)
plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position("top")
plt.xticks(rotation=45, ha="left")
plt.ylabel("")

plt.subplot2grid(n_subplots, (0, 1))
sns.heatmap(data[['log$_2$(FC) gene', 'log$_2$(FC) protein']],
            vmax=1.5, vmin=-1.5, center=0, cmap=cmap2, linewidths=0.05, square=False,
            annot=True, linecolor="lightgray", annot_kws={"size": 7},
            xticklabels=True, yticklabels=True,
            cbar=False)
plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position("top")
plt.xticks(rotation=45, ha="left")
plt.ylabel("")
plt.yticks([])

bases = ["A", "U", "C", "G"]
colors = dict(zip(bases, ["#109648", "#D62839", "#255C99", "#F7B32B"]))
plt.subplot2grid(n_subplots, (0, 2))

plot_dmr_distr(both, fig_title="Relative cDNA position")

plt.subplot2grid((1, 4), (0, 3))
plot_dmr_distr(zscore, plot_pos=False,  fig_title="Z-score")

plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='on')
plt.subplots_adjust(top=0.75, left=0.15, right=0.9, bottom=0.05, wspace=0.05)

plt.savefig("DMR_GENES_overview_2.svg")
os.startfile("DMR_GENES_overview_2.svg")
# plt.show()
print(data)
