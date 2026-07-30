import requests
import pandas as pd


def get_gene_name_from_transcript(transcript_id):
    server = "https://rest.ensembl.org"
    endpoint = f"/lookup/id/{transcript_id}"
    headers = {"Content-Type": "application/json"}

    response = requests.get(server + endpoint, headers=headers)

    if not response.ok:
        print(f"Error: {response.status_code} - {response.text}")
        return None

    data = response.json()
    gene_name = data.get("display_name", None)
    gene_id = data.get("Parent", None)

    if gene_name and gene_id:
        print(f"Transcript ID: {transcript_id}")
        print(f"Gene Name: {gene_name}")
        print(f"Gene ID: {gene_id}")
        return gene_name, gene_id
    else:
        print(f"Gene name or gene ID not found for transcript ID {transcript_id}")
        return None, None


if __name__ == '__main__':
    data = pd.read_table(
        r"direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_before_padj.tsv")

    data["id"] = data["id"].apply(lambda x: x.split(".")[0])
    gff_file = r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\0_REF_TABLES\mRNA_cDNA_release_113.tsv"
    gff = pd.read_table(gff_file)
    data["transcript_name"] = data["id"].map(dict(zip(gff.transcript_id, gff.transcript_name)))
    data["gene_id"] = data["id"].map(dict(zip(gff.transcript_id, gff.gene_id)))
    print(gff)
    print(data)
    genes = []
    gene_ids = []
    transcript_ids = []
    for t in data.drop_duplicates(subset="id").loc[data["transcript_name"].isna(), "id"]:
        print(t)
        g, gene_id = get_gene_name_from_transcript(t)
        transcript_ids.append(t)
        genes.append(g)
        gene_ids.append(gene_id)
    data.loc[data["transcript_name"].isna(), "transcript_name"] = data.loc[data["transcript_name"].isna(), "id"].map(dict(zip(transcript_ids, genes)))
    data.loc[data["transcript_name"].isna(), "gene_id"] = data.loc[data["transcript_name"].isna(), "id"].map(
        dict(zip(transcript_ids, gene_ids)))

    # Example usage
    # get_gene_name_from_transcript("ENST00000708566")
    data.to_csv(r"direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_before_padj_GENE_NAMES.tsv",
                sep="\t", index=False)



