import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import requests
import requests
# collaps transcripts/ isoforms info to gene level

# CMRs = sum
# TPMs = sum etc.

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def get_gene_name_from_gene_id(gene_id):
    server = "https://rest.ensembl.org"
    endpoint = f"/lookup/id/{gene_id}"
    headers = {"Content-Type": "application/json"}

    response = requests.get(server + endpoint, headers=headers)

    if not response.ok:
        print(f"Error: {response.status_code} - {response.text}")
        return None

    data = response.json()
    print(data)
    gene_name = data.get("display_name", None)
    print(gene_name)
    gene_id = data.get("id", None)

    if gene_name and gene_id:
        print(f"Gene ID: {gene_id}")
        print(f"Gene Name: {gene_name}")
        return gene_name
    else:
        print(f"Gene name or gene ID not found for transcript ID {gene_id}")
        return "."


xpore_file = r"\majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS_GENENAMES.tsv"
xpore_file = r"direct_RNA_seq\Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4" + xpore_file
mrna_level_file = r"..\Illumina\hSMC_TGFB_PDGF\DATA\results_C3_C4_2025\tables\logcount-matrix\model_condition.logcount-matrix_postprocessed.tsv"
mrna_level_ont_file = r"direct_RNA_seq/Nanocount/DESeq2_results_Tobias/ONT_transcript_level_C3_C4_deseq2_norm_tpm.tsv"
gene_level_ont_file = r"direct_RNA_seq/Nanocount/DESeq2_results_Tobias/ONT_gene_level_C3_C4_deseq2_norm_tpm.tsv"
degs_ont_deseq2_file = r"direct_RNA_seq/Nanocount/DESeq2_results_Tobias/ONT_gene_level_C3_C4_DEseq2_sum_stats.tsv"

# take DEseq2 normalized counts for Master table
degs_file_illumina = r"..\Illumina\hSMC_TGFB_PDGF\DATA\results_C3_C4_2025\tables\diffexp\model_condition.genes-representative.diffexp_postprocessed.tsv"
protein_level_file = r"..\..\Proteomics\hSMC_TGFB_PDGF\SMC_PDGF_IL_vs_TGFB_maxLFQ_protein_intensity_2pept_norm_log2.xlsx"
deps_file = r"..\..\Proteomics\hSMC_TGFB_PDGF\Report_SMC_PDGF_IL_vs TGFB_2_pept.xlsx"
# gene_file = r"../0_REF_TABLES\ensembl_gene_transcript_UTRs_release_113_cDNA_pos.tsv"
# cross_ref_file = "../0_REF_TABLES/Homo_sapiens.GRCh38.114.uniprot.tsv"
# cross_ref_file_cdna_len = "../0_REF_TABLES/cdna_fasta2df.tsv"
# cross_ref_file_ensembl = r"direct_RNA_seq/Nanocount/rep123_Nanocount_tsv/n3_S3_Col_1_Matrix_cdna_FAR26827_dcbb2f3fdb27aa26296f1aedc26052bae70586a6.tsv"

hgnc_ref_file = "../0_REF_TABLES/hgnc_db.tsv"
hgnc_ref = pd.read_table(hgnc_ref_file)
hgnc_ref = hgnc_ref[["uniprot_ids", "ensembl_gene_id", "symbol"]]


# load tables
xpore = pd.read_table(xpore_file)
xpore.middle_base = xpore.middle_base.replace("T", "U")

# xpore = xpore[xpore.adj_pval_s3_vs_s4 <= 0.05]
# xpore = xpore[abs(xpore.diff_mod_rate_s3_vs_s4) >= 0.25]
# xpore["transcript_name"] = xpore["transcript_name"].fillna(".")
# xpore["gene_name"] = xpore["transcript_name"].apply(lambda x: "-".join(x.split("-")[:-1]))
xpore_counts = xpore.groupby("Gene name").agg(# gene_name=("gene_name", "first"),
                                            dmr_transcript_count=("id", "nunique"),
                                            dmr_count=("kmer", "count"),
                                            dmr_A_count=("middle_base", lambda x: (x == "A").sum()),
                                            dmr_C_count=("middle_base", lambda x: (x == "C").sum()),
                                            dmr_G_count=("middle_base", lambda x: (x == "G").sum()),
                                            dmr_U_count=("middle_base", lambda x: (x == "U").sum()),
                                      ).sort_values(by="dmr_count", ascending=False)
xpore["id"] = xpore["id"].apply(lambda x: x.split(".")[0])
print(xpore_counts)

# load mRNA level
mrna_level_illumina = pd.read_table(mrna_level_file)
mrna_level_ont = pd.read_table(mrna_level_ont_file)
mrna_level_ont["mrna_ont_C3"] = np.log2(mrna_level_ont[["C3_R1", "C3_R2", "C3_R3"]].mean(axis=1) + 1)
mrna_level_ont["mrna_ont_C4"] = np.log2(mrna_level_ont[["C4_R1", "C4_R2", "C4_R3"]].mean(axis=1) + 1)
mrna_level_ont = mrna_level_ont.rename(columns={"HGNC symbol": "gene_name"})
mrna_level_ont = mrna_level_ont[["transcript_id", "gene_id", "mrna_ont_C3", "mrna_ont_C4", "gene_name", "transcript_length"]]

# calc transcript_count
mrna_level_ont = mrna_level_ont[(mrna_level_ont["mrna_ont_C3"] >= 1)
                                | (mrna_level_ont["mrna_ont_C3"] >= 1)]
print(mrna_level_ont)
transcript_counts = mrna_level_ont.groupby("gene_name").agg(transcript_count=("transcript_id", "nunique"),
                                                          # gene_name=("gene_name", "first")
                                                          ).sort_values(by="transcript_count", ascending=False)
print(transcript_counts)
# load gene level ont
gene_level_ont = pd.read_table(gene_level_ont_file, index_col=0)  # in TPM !!!
gene_level_ont["gene_level_ont_C3"] = np.log2(gene_level_ont[["C3_R1", "C3_R2", "C3_R3"]].mean(axis=1) + 1)
gene_level_ont["gene_level_ont_C4"] = np.log2(gene_level_ont[["C4_R1", "C4_R2", "C4_R3"]].mean(axis=1) + 1)
gene_level_ont = gene_level_ont.rename(columns={"HGNC symbol": "gene_name"})
gene_level_ont = gene_level_ont[["gene_level_ont_C3", "gene_level_ont_C4", # "gene_name",
                                 "median_transcript_len"]]

# load DEGs (ONT + Illumina)
degs_ont = pd.read_table(degs_ont_deseq2_file, index_col=0, usecols=[0, 2, 5, 6])

print(degs_ont)


degs_illumina = pd.read_table(degs_file_illumina)


# load protein
protein_level = pd.read_excel(protein_level_file)
deps = pd.read_excel(deps_file)
c4 = [x for x in protein_level.columns if "PDGF" in x]
c3 = [x for x in protein_level.columns if "TGF" in x]
protein_level["protein_level_C4"] = protein_level[c4].mean(axis=1)
protein_level["protein_level_C3"] = protein_level[c3].mean(axis=1)
# protein_level = protein_level.drop_duplicates(subset=["Single_Protein_ID"])

protein_level = pd.merge(protein_level, deps[["Single_Protein_ID", "Gene_names",
                                              "PDGF/TGFB|signal_log2_ratio",
                                              "PDGF/TGFB|raw_p_value",
                                              "PDGF/TGFB|adjusted_p_value"]],
                         on="Single_Protein_ID", how="left")
protein_level = protein_level[["Single_Protein_ID",
                               "Gene_names",
                               "protein_level_C3",
                               "protein_level_C4",
                               "PDGF/TGFB|signal_log2_ratio",
                               "PDGF/TGFB|raw_p_value",
                               "PDGF/TGFB|adjusted_p_value"]]

protein_level = pd.merge(protein_level, hgnc_ref, left_on="Single_Protein_ID", right_on="uniprot_ids", how="left")
protein_level = protein_level.drop(["Gene_names", "uniprot_ids"], axis=1).rename(columns={"symbol": "gene_name",
                                                                                    "ensembl_gene_id": "gene_id"})
protein_level = protein_level.dropna(subset=["gene_name"])
protein_level = protein_level.set_index("gene_name", verify_integrity=True)
# cross_ref = pd.read_table(cross_ref_file, usecols=[0, 1, 2, 3])
# cross_ref = cross_ref.drop_duplicates(subset="transcript_stable_id")
# cross_ref["xref"] = cross_ref["xref"].apply(lambda x: x.split("-")[0])
# cross_ref = cross_ref.drop_duplicates(subset="xref")
# protein_level = pd.merge(protein_level, cross_ref, left_on="Single_Protein_ID", right_on="xref", how="left")

print(protein_level)

# protein_level = protein_level.drop_duplicates(subset="gene_stable_id")
# protein_level = protein_level.set_index("gene_stable_id")


so = [transcript_counts,
      xpore_counts,
      gene_level_ont,
      degs_ont,
      protein_level,
      ]


# for tbl in so:
#     print("\n")
#     print(tbl.head(2))

gene_level_data = pd.concat(so, axis=1)
gene_level_data["dmr_transcript_count"] = gene_level_data["dmr_transcript_count"].fillna(0)
gene_level_data["transcript_count"] = gene_level_data["transcript_count"].fillna(0)
gene_level_data["dmr_count"] = gene_level_data["dmr_count"].fillna(0)
gene_level_data["dmr_A_count"] = gene_level_data["dmr_A_count"].fillna(0)
gene_level_data["dmr_C_count"] = gene_level_data["dmr_C_count"].fillna(0)
gene_level_data["dmr_G_count"] = gene_level_data["dmr_G_count"].fillna(0)
gene_level_data["dmr_U_count"] = gene_level_data["dmr_U_count"].fillna(0)

gene_level_data["gene_level_ont_C4"] = gene_level_data["gene_level_ont_C4"].fillna(0)
gene_level_data["gene_level_ont_C3"] = gene_level_data["gene_level_ont_C3"].fillna(0)

gene_level_data["protein_level_C3"] = gene_level_data["protein_level_C3"].fillna(0)
gene_level_data["protein_level_C4"] = gene_level_data["protein_level_C4"].fillna(0)
gene_level_data["PDGF/TGFB|signal_log2_ratio"] = gene_level_data["PDGF/TGFB|signal_log2_ratio"].fillna(0)
gene_level_data["protein_level_C4"] = gene_level_data["protein_level_C4"].fillna(0)

gene_level_data["PDGF/TGFB|raw_p_value"] = gene_level_data["PDGF/TGFB|raw_p_value"].fillna(1)
gene_level_data["PDGF/TGFB|adjusted_p_value"] = gene_level_data["PDGF/TGFB|adjusted_p_value"].fillna(1)
gene_level_data["Single_Protein_ID"] = gene_level_data["Single_Protein_ID"].fillna(".")

# gene_level_data["gene_name"] = gene_level_data["gene_name"].fillna(gene_level_data["Gene_names"])
# gene_level_data["gene_name"] = gene_level_data["gene_name"].fillna(".")

print("#################################\n\n")
print(gene_level_data)
# for col in gene_level_data.columns[:5]:
#     gene_level_data[col] = gene_level_data[col].fillna(0).astype(int)
gene_level_data = gene_level_data.sort_values(by="dmr_count", ascending=False)
gene_level_data = gene_level_data.reset_index()

print(gene_level_data.columns)

old = ['index', 'transcript_count', 'dmr_transcript_count', 'dmr_count', 'dmr_A_count', 'dmr_C_count', 'dmr_G_count', 'dmr_U_count', 'gene_level_ont_C3', 'gene_level_ont_C4', 'median_transcript_len', 'log2FoldChange', "pvalue", 'padj', 'Single_Protein_ID', 'protein_level_C3', 'protein_level_C4', 'PDGF/TGFB|signal_log2_ratio', 'PDGF/TGFB|raw_p_value', 'PDGF/TGFB|adjusted_p_value', 'gene_id']
new = ['gene_name', 'transcript_count', 'dmr_transcript_count', 'dmr_count', 'dmr_A_count', 'dmr_C_count', 'dmr_G_count', 'dmr_U_count', 'gene_level_ont_C3', 'gene_level_ont_C4', 'median_transcript_len', 'ont_gene_log2FC_C4vsC3', 'ont_gene_pval', 'ont_gene_qval', 'protein_id', 'protein_level_C3', 'protein_level_C4', 'protein_log2FC_C4vsC3', 'protein_pval', 'protein_qval', 'gene_id']
gene_level_data = gene_level_data.rename(columns=dict(zip(old, new)))

new_order = ['gene_name', 'gene_id', 'protein_id', 'transcript_count',
             'dmr_transcript_count', 'dmr_count', 'dmr_A_count', 'dmr_C_count', 'dmr_G_count', 'dmr_U_count',
             'gene_level_ont_C3', 'gene_level_ont_C4', 'ont_gene_log2FC_C4vsC3', 'ont_gene_pval', 'ont_gene_qval',
             'protein_level_C3', 'protein_level_C4', 'protein_log2FC_C4vsC3', 'protein_pval', 'protein_qval',
             'median_transcript_len']

gene_level_data = gene_level_data[new_order]
# gene_level_data["protein_id"] = gene_level_data["gene_id"].map(dict(zip(cross_ref["gene_stable_id"],
#                                                                         cross_ref["xref"])))

for col in [x for x in gene_level_data.columns if "count" in x]:
    gene_level_data[col] = gene_level_data[col].astype(int)

# gene = pd.read_table(gene_file, usecols=[8, 9, 10])
# gene = gene.drop_duplicates(subset=["gene_id"])
# gene_level_data.loc[gene_level_data.gene_name == ".", "gene_name"] = gene_level_data.loc[gene_level_data.gene_name == ".", "gene_id"].map(dict(zip(gene["gene_id"], gene["transcript_name"])))
# gene_level_data["gene_name"] = gene_level_data["gene_name"].fillna(".")
# gene_level_data["gene_name"] = gene_level_data["gene_name"].apply(lambda x: x.split("-")[0])
gene_names = pd.read_table("xpore_transcript_id_gene_name.tsv")
# header: Transcript stable ID   Gene stable ID Gene name Source of gene name
gene_level_data.loc[gene_level_data.gene_name == ".", "gene_name"] = gene_level_data.loc[gene_level_data.gene_name == ".", "gene_id"].map(dict(zip(gene_names["Gene stable ID"], gene_names["Gene name"])))

# [["uniprot_ids", "ensembl_gene_id", "symbol"]]
gene_id_map = dict(zip(hgnc_ref["symbol"], hgnc_ref["ensembl_gene_id"]))
protein_id_map = dict(zip(hgnc_ref["symbol"], hgnc_ref["uniprot_ids"]))
gene_level_data.loc[gene_level_data["gene_id"].isna(), "gene_id"] = gene_level_data.loc[gene_level_data["gene_id"].isna(), "gene_name"].map(gene_id_map)
gene_level_data.loc[gene_level_data["protein_id"] == ".", "protein_id"] = gene_level_data.loc[gene_level_data["protein_id"] == ".", "gene_name"].map(protein_id_map)
print(gene_level_data)


# load available xpore site per transcript / gene
xpore_pre_filter_file = r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore\direct_RNA_seq\Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4\majority_direction_kmer_diffmod_GENENAMES.table"
xpore_pre_filter = pd.read_table(xpore_pre_filter_file)
print(xpore_pre_filter.head(2))
xpore_pre_filter = xpore_pre_filter.rename(columns={"Transcript stable ID": "gene_name"})
xpore_sites_per_gene = xpore_pre_filter["Gene name"].value_counts(dropna=False).to_dict()
print(xpore_sites_per_gene)
gene_level_data["xpore_sites_per_gene"] = gene_level_data["gene_name"].map(xpore_sites_per_gene)
# print(gene_level_data.info())
gene_level_data.to_csv("direct_RNA_seq/0_TABLES/DMR_MASTER_TABLE_Gene_Level_PAPER_NEW_GENENAMES.tsv",
                       sep="\t",

                       index=False)