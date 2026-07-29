import pandas as pd

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


folder = r"direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025"
# file_in = "diffmod_s3-s4/majority_direction_kmer_diffmod.table"

file_in = "diffmod_s3-s4/majority_direction_kmer_diffmod_modRate_05.tsv"
file_out = file_in.split(".")[0] + "_GENOME_cDNA_UTR.tsv"

xpore = pd.read_table(f"{folder}/{file_in}")
genome = pd.read_table("../0_REF_TABLES/ensembl_gene_transcript_UTRs_release_113_cDNA_pos_cDNA_seq.tsv",
                       usecols=[0, 3, 4, 6, 8, 9, 10, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24])

# genome = genome.drop(["cdna_seq", "gene_id", "transcript_name"], axis=1)

print(genome.columns)
genome["gene_body_len"] = abs(genome["genomic_start"] - genome["genomic_end"])
xpore["id"] = xpore["id"].apply(lambda x: x.split(".")[0])
xpore = pd.merge(xpore, genome, left_on="id", right_on="transcript_id", how="left").drop("transcript_id", axis=1)
xpore.to_csv(f"{folder}/{file_out}", sep="\t", index=False)
