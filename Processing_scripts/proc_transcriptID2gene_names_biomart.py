# pip install pybiomart pandas
from pybiomart import Dataset
import pandas as pd
from pyensembl import EnsemblRelease
import numpy as np


def get_gene_name_from_biomart(data):
    attrs = ["ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "external_gene_source"]

    # example: tids = ["ENST00000446046"]  # ENST00000354785  ENST00000715744
    out = []
    batch_size = 100
    num = len(data["id"].unique())
    idx = 0
    t_ids = [x.split(".")[0] for x in data["id"].unique()]
    for x in range(0, len(t_ids)- batch_size, batch_size):
        idx += batch_size
        print(x, " / ", num)
        # tid = tid.split(".")[0]
        # print(tid)
        df = biomart_dataset.query(attributes=attrs,
                                   filters={"link_ensembl_transcript_stable_id": t_ids[x:x + batch_size]}
                                   )
        out.append(df)
        # print(df)
    gene_names = pd.concat(out)
    gene_names = gene_names.drop_duplicates(subset="Transcript stable ID")
    print(gene_names.head())
    gene_names = gene_names.rename(columns={"Transcript stable ID": "id"})
    gene_names.to_csv("xpore_transcript_id_gene_name.tsv", sep="\t", index=False)
    data["id"] = data["id"].apply(lambda x: x.split(".")[0])
    gene_names = pd.merge(data, gene_names, on="id", how="left")
    gene_names.to_csv(file_out_genename, sep="\t", index=False)


if __name__ == '__main__':

    # xpore = pd.read_table(r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore\direct_RNA_seq\Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4\majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS.tsv")
    # file_out_genename = r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore\direct_RNA_seq\Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4\majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS_GENENAMES.tsv"
    xpore = pd.read_table(
        r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore\direct_RNA_seq\Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4\majority_direction_kmer_diffmod.table")
    file_out_genename = r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore\direct_RNA_seq\Xpore\DATA\xpore_cdna.all_group_compare_2025\diffmod_s3-s4\majority_direction_kmer_diffmod_GENENAMES.table"

    biomart_dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    get_gene_name_from_biomart(xpore)
