import os.path
import pandas as pd
import glob
import statsmodels.stats.multitest as multi
from umap import UMAP
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib import colormaps
import seaborn as sns
# from statsmodels.stats.multitest import multipletests
import os
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.markers import MarkerStyle
import matplotlib.patches as patches
from matplotlib.transforms import Affine2D
import json
import re
from Bio import SeqIO
import requests
# from Gene_set_enrichment_analysis import enrichr_api as enrich

from scipy import stats
import squarify
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def make_text_elements(text, x=0.0, y=0.0, width=1.0, height=1.0, color='blue', edgecolor="none",
                       font=FontProperties(family='monospace')):
    # https://github.com/const-ae/motif_plotter
    tp = TextPath((0.0, 0.0), text, size=1, prop=font)
    bbox = tp.get_extents()
    bwidth = bbox.x1 - bbox.x0
    bheight = bbox.y1 - bbox.y0
    trafo = Affine2D()
    trafo.translate(-bbox.x0, -bbox.y0)
    trafo.scale(1 / bwidth * width, 1 / bheight * height)
    trafo.translate(x, y)
    tp = tp.transformed(trafo)
    return patches.PathPatch(tp, facecolor=color, edgecolor=edgecolor)


def plot_motif(data, col_name="kmer", num_positions=5):

    xpore_modi = ["lower", "higher"]
    for m in xpore_modi:
        data_temp_mode = data[data.mod_assignment == m]
        t = 0
        for tag in ["neg_zscore", "pos_zscore"]:
            if tag == "neg_zscore":
                data_z = data_temp_mode[data_temp_mode[z_scor_col] < 0]
            elif tag == "pos_zscore":
                data_z = data_temp_mode[data_temp_mode[z_scor_col] > 0]
            else:
                raise ValueError("Please define filter")
            # prepare tbl

            fig = plt.figure(figsize=(6, 8))
            plt.axis("off")
            plt.title("C" + cond[1] + " vs. C3", fontdict={'weight': 'bold', 'size': 20}, pad=50)
            plt.subplots_adjust(hspace=0.5, top=0.8)
            row = 1
            # if True:
            for middle, temp_data in data_z.groupby('middle_base'):
                # middle, temp_data = "All", tbl
                # temp_data = temp_data[temp_data[col_name].map(temp_data[col_name].value_counts()) > t]
                # temp_data = all_data.copy()
                if len(temp_data) == 0:
                    pass
                    # continue
                num_bases = len(temp_data)
                ratios = []
                for i in range(num_positions):

                    temp_dict = {}
                    temp_data[i] = temp_data[col_name].apply(lambda x: x[i])
                    for base in ["A", "U", "C", "G"]:
                        temp_dict[base] = len(temp_data[temp_data[i] == base]) / num_bases
                    temp_dict = dict(sorted(temp_dict.items(), key=lambda x:x[1]))
                    print(i, temp_dict)
                    ratios.append(temp_dict)

                # print(ratios)
                ax1 = fig.add_subplot(4, 1, row)
                plt.axis("off")
                plt.xticks([])
                plt.yticks([])
                # ax1.spines[['right', 'top', 'bottom']].set_visible(False)
                plt.arrow(x=-0.12, y=0, dx=0, dy=1, width=0.04, length_includes_head=True, facecolor='k',
                          head_length=0.12)
                x_shift = 0.35
                if num_positions > 5:
                    x_shift = 0.6
                plt.text(x=-x_shift, y=0.5, s=str(len(temp_data)), ha='center', va='center', rotation=90,
                         fontsize=19)

                row += 1
                ax1.set_ylim(0, 1)
                ax1.set_xlim(-0.3, num_positions)
                bases = ["A", "U", "C", "G"]
                colors = dict(zip(bases, ["#109648", "#D62839", "#255C99", "#F7B32B"]))
                # sentence_shape = make_text_elements('Hello', x=-0.5, y=0.25, width=0.5, height=0.5)
                for x in range(num_positions):
                    y_start = 0
                    for base, base_height in ratios[x].items():
                        shape = make_text_elements(base, x=x, y=y_start, width=1, height=base_height, color=colors[base])
                        ax1.add_patch(shape)
                        y_start = y_start + base_height
                # plt.title(middle, fontsize=20)
            plt.axis("off")
            plt.savefig(f"direct_RNA_seq/0_PLOTS/Motifs/{col_name}_{cond}_motif_plot_above{t}_{tag}_{m}.svg")
            plt.savefig(f"direct_RNA_seq/0_PLOTS/Motifs/{col_name}_{cond}_motif_plot_above{t}_{tag}_{m}.png")
            plt.show()
            plt.close(fig)



if __name__ == '__main__':
    z_scor_col = "z_score_s3_vs_s4"
    cond = "C4"
    filter_name = f"direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_padj_modrate_025.tsv"  #  majority_filer
    xpore_data = pd.read_table(filter_name)
    xpore_data["kmer"] = xpore_data["kmer"].str.replace("T", "U")
    xpore_data["middle_base"] = xpore_data["middle_base"].str.replace("T", "U")
    plot_motif(xpore_data, col_name="kmer", num_positions=5)














