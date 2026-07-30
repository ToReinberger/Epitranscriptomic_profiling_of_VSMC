import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy import stats
import numpy as np
from matplotlib import cm
from enrichr_api import fetch_enrichr_with_background
import os

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)



# temp = data.sort_values(by="gene_level_ont_C3", ascending=False)
# for g in temp.gene_name:
#     print(g, end=",")

"""
This gene list best represents extracellular matrix organization, angiogenesis, and smooth muscle cell phenotypic modulation toward a synthetic, pro-remodeling state.
Extracellular matrix (ECM) & remodeling:
FN1, THBS2, POSTN, VCAN, LAMC1, LAMA4, COL1A1, COL3A1, COL6A1, P4HA1, ADAMTS1, ADAM9, MXRA8 → strong ECM and matrix remodeling signature.

Angiogenesis & vascular signaling:
VEGFA, FGF2, GREM1, IL6ST, CEMIP, HIF1A, HEG1, NPR3, C1QTNF1 → all linked to vascular development, angiogenesis, or vascular inflammation.

Cell adhesion / migration:
ITGB1, PARVA, FERMT2, CALD1, CAP1, CERCAM → adhesion and cytoskeletal regulation.

Stress / plasticity regulators:
PHLDA1, IER3, NAMPT, NR3C1 (glucocorticoid receptor), TPBG, TNIP1, HSPs → stress, survival, inflammatory adaptation.

Protein folding / trafficking / secretion:
GANAB, P4HB, CALU, PDIA3, USO1, RAB10, RABGEF1 → ER/Golgi and secretory machinery, often up in synthetic/ECM-secreting phenotypes.
"""


def run_enrichment(overwrite=False):

    if not overwrite and os.path.exists(file_out):
        return pd.read_excel(file_out)

    gsea_results = fetch_enrichr_with_background(gene_list=dmr_genes,
                                                 background=expressed_genes)

    gsea_results = gsea_results[gsea_results["num_genes"] >= num_overlapping_genes]
    gsea_results = gsea_results[~gsea_results[("Term name"
                                               "")].str.contains("Mouse")]
    out = []
    for lib, tbl in gsea_results.groupby("library"):
        if lib == "ChEA_2022":
           tbl = tbl[tbl["Term name"].str.contains("SMOOTH MUSCLE Human")]
        temp = tbl.sort_values(by="Rank").head(2)
        temp["library"] = lib
        out.append(temp)
    gsea_results = pd.concat(out)

    gsea_results.to_excel(file_out, index=False)
    print(gsea_results)
    return gsea_results


def clean_terms(term):
    if term == "Regulation of IGF Transport and Uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)":
        term = "Insulin-like growth factor transport"
    if term == "Post-translational Protein Phosphorylation":
        term = "Protein Phosphorylation"
    if "ChIP-Seq" in term:
        term = term.split(" ")[0] + " response genes"
    if "Pperoxisome" in term:
        term = "Peroxisome"
    if "UV Response" in term:
        term = term.replace("Dn", "down")
    if " (GO" in term:
        term = term[:term.find(" (GO")]
    elif " R-HSA-" in term:
        term = term[:term.find(" R-HSA")]
    terms = term.split(" ")
    terms_temp = [terms[0]]
    for term in terms[1:]:
        if "NF-" in term or "A/1" in term or "I-" in term:
            terms_temp.append(term)
        else:
            terms_temp.append(term.lower())
    return " ".join(terms_temp)



def plot_enrichment():
    libs = ["MSigDB Hallmark 2020",
            "KEGG 2021 Human",
            "Reactome 2022",
            "GO Biological Process 2023",
            "ChEA_2022",
            "Reactome_Pathways_2024",
            "GWAS_Catalog_2025"
            # "DSigDB",
            # "GWAS_Catalog_2023",
            # "Human_Gene_Atlas",
            # "GTEx_Tissues_V8_2023",
            ]
    libs = dict(zip([x.replace(" ", "_") for x in libs], ["MSig", "K", "R", "GO", "ChEA", "RP", "GWAS"]))

    color_code = {"C4": "#FDB462", "C3": "#8DD3C7", "general": "#555555"}
    color = color_code["general"]

    df = gsea.copy()
    df = df[(~df.library.isin(["GO_Molecular_Function_2023",
                               "TargetScan_microRNA",
                               # "KEGG_2021_Human",
                               "GO_Cellular_Component_2023"]))
            & (df.num_genes > num_overlapping_genes)]

    df = df[(~df["Term name"].str.contains("External Encapsulating Structure Organization"))]

    df = df.sort_values(by="Odds ratio", ascending=True)
    df.reset_index(drop=True, inplace=True)

    print(df)
    df["Term name"] = df["Term name"].apply(clean_terms)
    df = df.sort_values(by="num_genes", ascending=False).drop_duplicates(subset=["Term name"]).sort_values(
        by="num_genes", ascending=True)
    df.reset_index(drop=True, inplace=True)
    terms = df["Term name"].values

    plt.figure(figsize=(4, 3))
    plt.style.use("ggplot")
    plt.xlabel("Odds ratio", fontsize=11)
    plt.ylim(-0.5, len(df) - 0.5)
    # plt.xlim(34, 0)
    plt.grid(False)
    # plt.gca().yaxis.tick_right()
    plt.barh(df.index, df["Odds ratio"], color=color, zorder=1, alpha=0.8)

    # plt.gca().tick_params(axis="y", direction="in", pad=-10, zorder=2)
    plt.yticks(np.arange(0, len(df)), terms, fontsize=10, ha="left", zorder=2)
    plt.yticks([])
    x_pos = 0.05 * df["Odds ratio"].max()
    for y, name, num_genes, lib in zip(np.arange(0, len(df)), terms, df.num_genes, df.library):
        plt.text(x_pos, y, name + f" ({num_genes}, {libs[lib]})", fontsize=9.5, ha="left", va="center")
    plt.title(title_name)
    plt.tight_layout()
    plt.subplots_adjust(right=0.83, bottom=0.25)
    plt.savefig(file_out_plot + ".svg")
    plt.show()


if __name__ == '__main__':
    data = pd.read_table("DMR_MASTER_TABLE_Gene_Level_PAPER.tsv")
    data["median_transcript_len"] = data["median_transcript_len"].fillna(0) + 1
    data["transcript_count"] = data["transcript_count"].fillna(0) + 1
    data["norm_DMR_ratio"] = (data["dmr_count"] / data["median_transcript_len"]) / data["transcript_count"]
    print(data.columns)

    dmr_thresh = 1
    expr = 1
    folder = "direct_RNA_seq/0_TABLES"
    folder2 = "direct_RNA_seq/0_PLOTS/DMR_mRNA_Enrichment"
    file_name = f"/DMR_genes_gsea_results_above_{dmr_thresh}DMRs_protein_sign_up_C4_vs_C3_background_up.xlsx"  # _protein_sign_up_C4_vs_C3
    file_out = folder + file_name
    file_out_plot = folder2 + file_name
    title_name = "GSEA of DMR+ mRNAs"
    print(data)
    num_overlapping_genes = 5

    dmr_genes = data.loc[
        (data.dmr_count > dmr_thresh)
        & (data.protein_log2FC_C4vsC3 > 0)
        & (data.protein_qval <= 0.05)
        & ((data.gene_level_ont_C3 >= expr) | (data.gene_level_ont_C4 >= expr))
        & ((data.protein_level_C3 > 0) | (data.protein_level_C4 > 0)),
        'gene_name']

    # print(dmr_genes)

    expressed_genes = data.loc[
        (data.dmr_count >= 0)
        & (data.protein_log2FC_C4vsC3 > 0)
        # & (data.protein_qval <= 0.05)
        & ((data.gene_level_ont_C3 >= expr) | (data.gene_level_ont_C4 >= expr))
        & ((data.protein_level_C3 > 0) | (data.protein_level_C4 > 0)),
        'gene_name']

    gsea = run_enrichment(overwrite=True)
    print(gsea)
    plot_enrichment()

