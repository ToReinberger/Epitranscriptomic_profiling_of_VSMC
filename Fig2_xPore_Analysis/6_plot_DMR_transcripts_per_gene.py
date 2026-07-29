import os.path
import pandas as pd
from pybiomart import Dataset
from typing import List, Dict, Any, Tuple
from collections import defaultdict
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

marker_dict = {
     # '.': 'point',
     # ',': 'pixel',
     'o': 'circle',
     # 'v': 'triangle_down',
     # '^': 'triangle_up',
     '<': 'triangle_left',
     '>': 'triangle_right',
     # '1': 'tri_down',
     # '2': 'tri_up',
     # '3': 'tri_left',
     # '4': 'tri_right',
     '8': 'octagon',
     's': 'square',
     'p': 'pentagon',
     '*': 'star',
     'h': 'hexagon1',
     'H': 'hexagon2',
     'D': 'diamond',
     'd': 'thin_diamond',
     # '|': 'vline',
     # '_': 'hline',
     'P': 'plus_filled',
     'X': 'x_filled',
    # '+': 'plus',
    # 'x': 'x',
     # 0: 'tickleft',
     # 1: 'tickright',
     # 2: 'tickup',
     # 3: 'tickdown',
     # 4: 'caretleft',
     # 5: 'caretright',
     # 6: 'caretup',
     # 7: 'caretdown',
     # 8: 'caretleftbase',
     # 9: 'caretrightbase',
     # 10: 'caretupbase',
     # 11: 'caretdownbase'
}
filled_markers = (
        'o', 'v', '^', '<', '>', '8', 's', 'p', 'h', 'H', 'D', 'd',
        'P', 'X')
marker_style = []
fill_style = ["full", "left", "right", "top", "bottom"]
for style in fill_style:
    for marker in filled_markers:
        marker_style.append((marker, style))




"""def plot_bar(tbl):
    # plot kmers as network
    kmers_counts = tbl.value_counts("kmer")
    kmers_count_map = dict(kmers_counts)
    for base, tbl in tbl.groupby("middle_base"):
        tbl = tbl[abs(tbl.z_score) < 100]
        # plot bar
        print(tbl.nlargest(10, "z_score"))
        print(tbl.nsmallest(10, "z_score"))
        tbl["kmer_4"] = tbl["kmer"].apply(lambda x: x[0:4])
        tbl["kmer_3"] = tbl["kmer"].apply(lambda x: x[1:4])
        # plt.show()
        tbl2 = tbl[["kmer_3", "z_score"]].groupby("kmer_3").quantile(0.3).reset_index("kmer_3")
        plt.figure(figsize=(7, 16))
        plt.title(f"C{cond[1]} vs C3 | {base} #{len(tbl)}", fontsize=18)
        sns.barplot(y="kmer_3", x="z_score", tbl=tbl2, # cut=1
                    )
        plt.savefig(f"Xpore/PLOTS/Barplot_kmer_3_30_quantile_{cond}_{base}.png")
        plt.close()
        print(tbl2)"""


def plot_single_transcripts_per_gene(data, save_img=False, middle_base="All"):
    expr = load_data_sets()

    expr2 = pd.read_table(r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore\direct_RNA_seq\Nanocount\DESeq2_results_Tobias\ONT_transcript_level_C3_C4_DEseq2_sum_stats.tsv")

    data["id"] = data["id"].apply(lambda x: x.split(".")[0])
    expr["transcript"] = expr["transcript"].apply(lambda x: x.split(".")[0])
    print(data.columns, expr.columns)
    data = pd.merge(data, expr2, left_on="id", right_on="transcript_id", how="left")

    # transcript_gene = transcript_gene[["transcript", "strand", "seqid", "start", "end"]]
    # transcript_gene["start"] = transcript_gene["start"] * 10**-6
    # transcript_gene["end"] = transcript_gene["end"] * 10 ** -6
    # data["gene"] = data["transcript_name"].apply(lambda x: "-".join(x.split("-")[:-1]))
    data = data.rename(columns={"Strand": "strand",
                                "Chromosome/scaffold name": "seqid",
                                "Gene start (bp)": "start",
                                "Gene end (bp)": "end"})

    data["start"] = data["start"] * 10**-6
    data["end"] = data["end"] * 10 ** -6

    for gene, gene_tbl in data.groupby("Gene name"):
        if len(gene_tbl["transcript_id"].unique()) <= 2:
            pass
        if gene != "ITGB1":
            pass
        if gene not in selected_genes:
            continue
        print(gene, gene_tbl)
        gene_tbl = gene_tbl.sort_values(by="baseMean", ascending=True)
        transcript_ids = gene_tbl["id"].unique()
        exon_tbl = get_exon_info_from_biomart(transcript_ids)
        out = []
        t_ids = exon_tbl["Transcript stable ID"].unique()

        for t in t_ids:
            tbl = exon_tbl[exon_tbl["Transcript stable ID"] == t]
            tbl["exon_len"] = abs(tbl["Exon region end (bp)"] - tbl["Exon region start (bp)"])
            tbl = tbl.sort_values(by="Exon rank in transcript")
            tbl["exon_start"] = 1
            tbl["exon_end"] = 1
            idx_start = 0
            tbl.loc[tbl["Exon rank in transcript"] == 1, "exon_end"] = tbl.loc[tbl["Exon rank in transcript"] == 1, "exon_len"].values + 1
            for i in tbl["Exon rank in transcript"].values[1:]:
                tbl.loc[tbl["Exon rank in transcript"] == i, "exon_start"] = tbl.loc[
                    tbl["Exon rank in transcript"] == i - 1, "exon_end"].values + 1
                tbl.loc[tbl["Exon rank in transcript"] == i, "exon_end"] = tbl.loc[
                    tbl["Exon rank in transcript"] == i, "exon_start"].values + tbl.loc[
                    tbl["Exon rank in transcript"] == i, "exon_len"].values
            out.append(tbl)
            print(tbl)
        exon_tbl = pd.concat(out, ignore_index=True)

        for col in ["5' UTR start", "5' UTR end",  "3' UTR start",  "3' UTR end",
                     "Exon region start (bp)", "Exon region end (bp)"]:
            exon_tbl[col] = exon_tbl[col] * 10**-6

        # strand = gene_tbl.strand.values[0]
        print(gene_tbl)
        chrom = exon_tbl["Chromosome/scaffold name"].values[0]
        unique_transcripts = len(gene_tbl["transcript_id"].unique())
        if unique_transcripts != 1:
            pass
        if len(gene_tbl) < 8:
            pass
        # gene_tbl = gene_tbl.sort_values(by="middle_base")
        unique_kmers = list(dict.fromkeys(gene_tbl.sort_values(by="middle_base").kmer.values))
        kmer_dict = {"A": [], "C": [], "G": [], "U": []}

        for km in unique_kmers:
            mb = km[2]
            kmer_dict[mb].append(km)

        bases = ["A", "U", "C", "G"]
        colors = dict(zip(bases, ["#109648", "#D62839", "#255C99", "#F7B32B"]))
        for mb in kmer_dict.keys():
            a_temp = kmer_dict[mb]
            kmer_dict[mb] = "\n".join(a_temp)

        if unique_transcripts == 1:
            h = 3
            space = 0.25
            top = 0.7
        elif unique_transcripts < 4:
            h = unique_transcripts + 1
            space = 0.25
            top = 0.7
        else:
            h = unique_transcripts - 1
            space = 0.3
            top = 0.85
        if h > 6:
            h = 6
        plt.figure(figsize=(10, int(h * 0.85)))  # default figsize=(10, int(h * 2.5)
        plt.subplots_adjust(left=0.25, bottom=0.2, right=0.72, top=top)
        start, end = gene_tbl.start.min(), gene_tbl.end.max()
        plot_distance = end - start
        transcript_num = 0
        transcripts = []
        gene_tbl = gene_tbl.sort_values(by="baseMean", ascending=True)
        # print(gene_tbl)
        max_pos = gene_tbl.end.max()
        min_pos = gene_tbl.start.min()

        ####################################################################################
        # plot transcripts
        ####################################################################################
        strand = exon_tbl.Strand.values[0]

        for iso in transcript_ids:
            iso_tbl = gene_tbl[gene_tbl["transcript_id"] == iso].copy()
            if iso_tbl.empty:
                continue
            exon_tbl_for_id = exon_tbl[exon_tbl["Transcript stable ID"] == iso]
            # log_count_c3 = round(iso_tbl["mrna_ont_C3"].values[0], 1)
            # log_count_c4 = round(iso_tbl["mrna_ont_C4"].values[0], 1)
            log2fc = round(iso_tbl["log2FoldChange"].values[0], 2)
            iso_name = iso
            # transcripts.append(iso_name + f"\nC3={log_count_c3} | C4={log_count_c4}")
            transcripts.append(iso_name + f"\nlog$_2$(FC)={log2fc}")
            iso = iso.split(".")[0]
            transcr_start = exon_tbl_for_id["Exon region start (bp)"].min()
            transcr_end = exon_tbl_for_id["Exon region end (bp)"].max()
            plt.hlines(transcript_num, transcr_start, transcr_end, linewidth=2, linestyles="solid", color="lightgray")

            ###########################################################################################################
            # plot UTRs  ##############################################################################################
            ###########################################################################################################
            print("plot UTRs")
            five_prime = exon_tbl_for_id[["5' UTR start", "5' UTR end"]].dropna()
            for s, e in zip(five_prime["5' UTR start"], five_prime["5' UTR end"]):
                plt.hlines(transcript_num, s, e, linewidth=16,
                           linestyles="solid", color="silver", zorder=10)
            print(five_prime)
            three_prime = exon_tbl_for_id[["3' UTR start", "3' UTR end"]].dropna()
            for s, e in zip(three_prime["3' UTR start"], three_prime["3' UTR end"]):
                plt.hlines(transcript_num, s, e, linewidth=16,
                           linestyles="solid", color="silver", zorder=10)

            print(three_prime)

            if strand == 1:
                plt.text(x=five_prime["5' UTR start"].values[0], y=transcript_num, s="5'  ", ha="right",
                         va="center", color="gray")
                plt.text(x=three_prime["3' UTR end"].values[0], y=transcript_num, s="  3'", ha="left", va="center",
                         color="gray")
            if strand == -1:
                plt.text(x=five_prime["5' UTR end"].values[0], y=transcript_num, s="  5'", ha="left", va="center",
                         color="gray")
                plt.text(x=three_prime["3' UTR start"].values[0], y=transcript_num, s="3'  ", ha="right",
                         va="center", color="gray")


            #################################################################################################
            # plot exons
            #################################################################################################
            print("plot exons")
            c = 1
            exon_block = 0
            for exon_start, exon_end in zip(exon_tbl_for_id["Exon region start (bp)"],
                                            exon_tbl_for_id["Exon region end (bp)"]):
                if abs(transcr_end - transcr_start) > 5000:
                    if abs(exon_end - exon_start) < (300 * 10 ** -6):
                        exon_end = exon_start + (300 * 10 ** -6)
                plt.hlines(transcript_num, exon_start, exon_end, linewidth=16,
                           linestyles="solid",
                           color="dimgray", zorder=2)
                exon_block = exon_block + (exon_end - exon_start)
                # plt.scatter(x=(exon_end + exon_start)/2, y=transcript_num)
                c += 1

            ###########################################################################################################
            # prepare_kmer_pos / get genome position
            ###########################################################################################################

            new_genome_pos = []
            for kmer, pos in zip(iso_tbl.kmer, iso_tbl.position):
                exon_tbl_for_id_temp = exon_tbl_for_id[(exon_tbl_for_id["exon_start"] <= pos) &
                                                       (exon_tbl_for_id["exon_end"] >= pos)]
                if exon_tbl_for_id_temp.empty:
                    print("DMR not in exon")
                    new_genome_pos.append(np.nan)
                    continue

                distance = pos - exon_tbl_for_id_temp["exon_start"].values[0]
                exon_start = exon_tbl_for_id_temp["Exon region start (bp)"].values[0] * 10**6
                exon_end = exon_tbl_for_id_temp["Exon region end (bp)"].values[0] * 10 ** 6

                if strand == 1:
                    new_pos = exon_start + distance
                else:
                    new_pos = exon_end - distance
                new_genome_pos.append(new_pos)
            iso_tbl["genome_position"] = new_genome_pos

            print("plot DMRs")
            pos_temps, pos_temps2 = [0], [0]
            for position, tbl3 in iso_tbl.groupby("genome_position"):
                kmer_count_per_pos = tbl3.value_counts("position")
                # print(kmer_count_per_pos)
                for kmer_ in tbl3.kmer.unique():
                    mb = kmer_[2]
                    transcript_num_temp = transcript_num
                    position = position * 10**-6
                    if (abs(pos_temps[-1] - position) < plot_distance * 0.008 and not
                            abs(pos_temps2[-1] - position) < plot_distance * 0.008):
                        transcript_num_temp = transcript_num_temp + space
                        pos_temps2.append(position)
                    elif (abs(pos_temps[-1] - position) < plot_distance * 0.008 and
                            abs(pos_temps2[-1] - position) < plot_distance * 0.008):
                        transcript_num_temp = transcript_num_temp - space
                    pos_temps.append(position)
                    z_score_val = tbl3.z_score_s4_vs_s3.median()
                    if z_score_val < 0:
                        ms = "v"
                    else:
                        ms = "^"
                    size = abs(z_score_val) * 15
                    """marker_style_ = dict(color=color, marker=m[0],
                                        markersize=15, markerfacecoloralt='tab:red')"""
                    plt.scatter(x=position, y=transcript_num_temp,
                                marker=MarkerStyle(ms,   # m[0]
                                    fillstyle="full"),
                                s=size, c=colors[mb], zorder=99,
                                edgecolors="k", linewidths=0.6,
                                )
                    """plt.scatter(x=position, y=transcript_num_temp, marker=MarkerStyle(m[0], fillstyle="left"),
                                s=60, c=m[2], zorder=100,
                                edgecolors="k", linewidths=0.5,
                                )"""
                    """plt.text(x=position, y=transcript_num + .3, s=",".join(tbl3.kmer.values),
                             fontdict={"weight": "regular", "size": 6}, color=color, ha="center", va="bottom")"""
            transcript_num += 1

        # plt.xlim(start, end)
        if not save_img:
            plt.close()
            continue
        # plot legend
        legend_elements = []
        step = (max_pos - min_pos) * 0.09
        x_pos = max_pos + (max_pos - min_pos) * 0.15
        # print(x_pos)
        if transcript_num == 1:
            y_pos1 = 0
            y_pos2 = -0.35
        elif transcript_num < 4:
            y_pos1 = transcript_num / 2
            y_pos2 = y_pos1 * 0.6
        else:
            y_pos1 = transcript_num * 0.67
            y_pos2 = transcript_num * 0.6
        for k, v in kmer_dict.items():
            if len(v) == 0:
                continue
            plt.text(x=x_pos, y=y_pos1, s=k, rotation=0, ha="center", va="center",
                     fontdict={"weight": "regular", "size": 9, # "font": "Courier New"
                               },
                     bbox=dict(facecolor=colors[k], edgecolor='k', boxstyle='round,pad=0.3',
                               alpha=0.7))
            plt.text(x=x_pos, y=y_pos2, s=v, rotation=0, ha="center", va="top",
                     fontsize=6.5)
            x_pos = x_pos + step
            legend_elements.append(plt.Line2D([0], [0],
                                              marker="s", fillstyle="full", linestyle='None', markerfacecolor=colors[k],
                                              mfcalt='white', label=v, mec="k", color='w', mew=0.5, ms=6))

        z_score_legend = [plt.Line2D([0], [0],
                                     marker="v", fillstyle="full", linestyle='None', markerfacecolor="gray",
                                     label="<0", mec='gray', color='w', mew=0.25, ms=6),
                          plt.Line2D([0], [0],
                                     marker="^", fillstyle="full", linestyle='None', markerfacecolor="gray", label=">0",
                                     mec='gray', color='w', mew=0.25, ms=6)
                          ]

        if transcript_num == 1:
            y_pos3 = 0.7
        else:
            y_pos3 = 0.9

        legend2 = plt.legend(handles=z_score_legend, fontsize=7.5,
                             ncol=2, title="z-score",
                             title_fontsize=9, frameon=False, facecolor='white',
                             handletextpad=0.15, labelspacing=0.5, bbox_to_anchor=(1.04, y_pos3),
                             loc='upper left')  # default upper right

        # legend2.get_frame().set_linewidth(0.3)
        legend2._legend_box.align = "left"
        plt.gca().add_artist(legend2)
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.xlabel(f"Chr{chrom} (Mb)")
        plt.tick_params(left=False)
        plt.yticks(np.arange(0, len(transcripts)), transcripts)
        if strand == -1:
            sign = "< (reverse)"
        else:
            sign = "> (forward)"
        plt.title(gene + sign, pad=0)
        plt.gca().spines[['right', 'top', 'left']].set_visible(False)

        plt.ylim(-0.5, transcript_num + 0.5)
        plt.savefig(f"{gene}_xpore_{xpore_mode}.png")
        plt.savefig(f"{gene}_xpore_{xpore_mode}.svg")
        plt.show()


def load_data_sets(load_gff=False):
    # load_gff

    if load_gff:
        exons_ = pd.read_table(r"C:\Users\tobia\PycharmProjects\GENES\Genebank_Ensembl_GFF3\gff_exons.tsv", low_memory=False)
        temp_gff_ = pd.read_table(r"C:\Users\tobia\PycharmProjects\GENES\Genebank_Ensembl_GFF3\transcript_gene.tsv", usecols=[8, 9])
        # temp_gff_["transcript"] = temp_gff_["transcript"].apply(lambda x: x.split(".")[0])
        exons_["transcript"] = exons_["transcript"].apply(lambda x: x.split(".")[0])
        temp_gff_ = pd.merge(exons_, temp_gff_, on="transcript", how="left")
        temp_gff_["exon_len"] = temp_gff_["end"] - temp_gff_["start"]
        transcritps, genes, transcript_len = [], [], []
        print(temp_gff_)
        for t, tbl in temp_gff_.groupby(["transcript", "gene"]):
            transcritps.append(t[0])
            genes.append(t[1])
            transcript_len.append(tbl["exon_len"].sum())
        gff_ = pd.DataFrame(data={"protein": genes,
                                  "transcript": transcritps,
                                  "transcript_len": transcript_len})

    # load_expr_data
    overwrite_expr = False
    expr_temp_file = "plot_transcr_kmer_corr_temp.tsv"
    if not os.path.exists(expr_temp_file) or overwrite_expr:
        expr_ = pd.read_table(r"ONT_transcript_level_C3_C4_deseq2_norm_tpm.tsv")
        expr_["mrna_ont_C4"] = expr_[["C4_R1", "C4_R2", "C4_R3"]].mean(axis=1)
        expr_["mrna_ont_C3"] = expr_[["C3_R1", "C3_R2", "C3_R3"]].mean(axis=1)
        expr_["mrna_ont_C4"] = np.log2(expr_["mrna_ont_C4"] + 1)
        expr_["mrna_ont_C3"] = np.log2(expr_["mrna_ont_C3"] + 1)
        print(expr_[expr_["HGNC symbol"] == "ITGB1"])

        conditions = ["mrna_ont_C3", "mrna_ont_C4"]
        expr_["FC"] = expr_["mrna_ont_C4"] - expr_["mrna_ont_C3"]
        expr_["log_count"] = expr_[conditions].mean(axis=1)
        expr_ = expr_.rename(columns={"transcript_id": "transcript"})
        expr_ = expr_[["transcript", "log_count", "FC", "mrna_ont_C3", "mrna_ont_C4", "HGNC symbol"]]
        expr_.to_csv(expr_temp_file, sep="\t", index=False)
    else:
        expr_ = pd.read_table(expr_temp_file)
    return expr_


def get_exon_info_from_biomart(tids=Any):
    attrs = [
        "ensembl_transcript_id", "ensembl_exon_id", "rank",
        "chromosome_name", "strand",
        "exon_chrom_start", "exon_chrom_end",
        "cds_start", "cds_end",
        "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end",
    ]

    gene_out = []
    for tid in tids:
        df = dataset.query(attributes=attrs,
                           filters={"link_ensembl_transcript_stable_id": tid})
        print(df)
        gene_out.append(df)
    gene_out_final = pd.concat(gene_out)
    return gene_out_final


def get_gene_names(tids=Any):
    attrs = [
        "ensembl_transcript_id",
        "ensembl_gene_id",
        "external_gene_name",
        "chromosome_name",
        "strand",
        "start_position",
        "end_position"
    ]

    gene_out = []
    for tid in tids:
        print(tid)
        df = dataset.query(attributes=attrs,
                           filters={"link_ensembl_transcript_stable_id": tid})
        print(df)

        gene_out.append(df)
    gene_out_final = pd.concat(gene_out)
    print(gene_out_final)
    return gene_out_final


if __name__ == '__main__':

    # load biomart dataset
    dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    print(dataset.filters.keys())
    # quit()
    # load gene level table
    cond = "C3"
    gene_level_data = pd.read_table('DMR_MASTER_TABLE_Gene_Level_C3_vs_C4_padj_DMR_025_PAPER.tsv')
    dmr_count_thresh = 5
    gene_level_data["is_DMR"] = False
    gene_level_data.loc[gene_level_data["dmr_count"] >= dmr_count_thresh, "is_DMR"] = True
    gene_level_data = gene_level_data[
        (gene_level_data[f"gene_level_ont_{cond}"] > 0)
        & (gene_level_data[f"protein_level_{cond}"] > 0)
        ]
    gene_level_data["ont_gene_log2FC_C4vsC3"] = gene_level_data["ont_gene_log2FC_C4vsC3"].fillna(0)
    gene_level_data["protein_log2FC_C4vsC3"] = gene_level_data["protein_log2FC_C4vsC3"].fillna(0)
    selected_genes = gene_level_data.loc[(gene_level_data.protein_log2FC_C4vsC3 > 0)
                                         & gene_level_data.is_DMR, "gene_name"].values

    selected_genes = ["POSTN"]

    # ["FN1", "THBS2", "ITGB1", "HDLBP", "POSTN", "CALD1", "SULF1", "FGF2", "NAMPT", "PHLDA1", "CDV3"]


    xpore_file = r"majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS_GENENAMES.tsv"
   
    # load tables
    xpore = pd.read_table(xpore_file)
    # xpore = xpore[xpore.adj_pval_s3_vs_s4 <= 0.05]
    # xpore = xpore[abs(xpore.diff_mod_rate_s4_vs_s3) >= 0.25]
    # xpore = xpore.dropna(subset="gene_id")
    # xpore["diff_mod_rate_s3_vs_s4"] = xpore["diff_mod_rate_s3_vs_s4"] * -1
    # xpore["z_score_s3_vs_s4"] = xpore["z_score_s3_vs_s4"] * -1
    # xpore = xpore.rename(columns={"diff_mod_rate_s3_vs_s4": "diff_mod_rate_s4_vs_s3",
    #                               "z_score_s3_vs_s4": "z_score_s4_vs_s3"})

    # out = get_gene_names(xpore["id"].unique())
    # # out = out.drop("Gene name", axis=1)
    # for col in ["Gene stable ID", "Chromosome/scaffold name", "Strand", "Gene start (bp)", "Gene end (bp)"]:
    #     xpore[col] = xpore["id"].map(dict(zip(out["Transcript stable ID"], out[col])))
    # # xpore["gene_name"] = xpore["id"].map(dict(zip(out["Transcript stable ID"], out["Gene name"])))
    # xpore.to_csv(xpore_file, sep="\t", index=False)
    # quit()
    xpore_mode = "higher"  # lower  or higher or both
    if xpore_mode != "both":
        xpore = xpore[xpore["mod_assignment"] == xpore_mode]
    xpore["z_score_s4_vs_s3"] = xpore["z_score_s3_vs_s4"] * -1
    plot_single_transcripts_per_gene(xpore, save_img=True, middle_base="All")











