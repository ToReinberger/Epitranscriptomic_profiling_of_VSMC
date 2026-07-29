# pip install pybiomart pandas
from pybiomart import Dataset
import pybiomart
import pandas as pd
from pyensembl import EnsemblRelease
import numpy as np


print(pybiomart.__version__)
quit()

pd.set_option('display.max_rows', 100)

pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def ensembl_genomic_to_cdna(position, transcript_id):
    #
    # in cmd terminal
    # 1) pip install pyensembl
    # 2) pyensembl install --release 111 --species homo_sapiens

    # Example usage
    """chromosome = "1"
    genomic_position = 65419 + 20  # Replace with your genomic position
    transcript_id = "ENST00000641515"  # Replace with your transcript ID
    cdna_position = genomic_to_cdna(genomic_position, transcript_id)
    print(cdna_position)
    quit()
    """

    # Get the transcript object
    # print(transcript)
    # Convert genomic position to spliced offset
    # Initialize the EnsemblRelease object

    # print(position, transcript_id)
    try:
        transcript = ensembl_data.transcript_by_id(transcript_id)
        cdna_position = transcript.spliced_offset(int(position))
        return int(cdna_position)
    except ValueError:
        return None


def get_utrs_from_biomart():
    attrs = [
        # "ensembl_transcript_id", "ensembl_exon_id", "rank",
        # "chromosome_name", "strand", "cdna",
        "start_position", "end_position",
        "rank", "exon_chrom_start", "exon_chrom_end",
        # "5utr", "3utr",
        "strand", "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end",
    ]

    # quit()

    # example: tids = ["ENST00000446046"]  # ENST00000354785  ENST00000715744
    out = []
    num = len(xpore["id"].unique())
    idx = 0
    for tid in xpore["id"].unique():
        idx += 1

        print(idx, " / ", num)
        tid = tid.split(".")[0]
        print(tid)

        df = biomart_dataset.query(attributes=attrs,
                           filters={"link_ensembl_transcript_stable_id": tid}
                           )
        if df.empty:
            f"{tid} is not in Ensembl"
            continue

        df["exon_len"] = df["Exon region end (bp)"] - df["Exon region start (bp)"] + 1
        df["five_prime"] = False
        df["three_prime"] = False
        df["five_prime_len"] = 0
        df["three_prime_len"] = 0
        df.loc[df["5' UTR start"] > 0, "five_prime"] = True
        df.loc[df["3' UTR start"] > 0, "three_prime"] = True

        df.loc[df["five_prime"], "five_prime_len"] = abs(df.loc[df["five_prime"], "5' UTR start"]
                                                         - df.loc[df["five_prime"], "5' UTR end"]) + 1
        df.loc[df["three_prime"], "three_prime_len"] = abs(df.loc[df["three_prime"], "3' UTR start"]
                                                           - df.loc[df["three_prime"], "3' UTR end"]) + 1
        df["cDNA_start_fiveUTR"] = 0
        df["cDNA_start_threeUTR"] = 0
        df.loc[df["three_prime"], "cDNA_start_threeUTR"] = df["exon_len"].sum() - df.loc[df["three_prime"], "exon_len"]

        # print(df)
        # print(df)
        collapsed = (
            df.agg(
                start=("Gene start (bp)", "min"),
                end=("Gene end (bp)", "max"),
                # "Exon region start (bp)": "min",
                #  "Exon region end (bp)": "max",
                five_prime_len=("five_prime_len", "sum"),
                three_prime_len=("three_prime_len", "sum"),
                cDNA_start_fiveUTR=("cDNA_start_fiveUTR", "sum"),
                cDNA_start_threeUTR=("cDNA_start_threeUTR", "sum"),
                cdna_len=("exon_len", "sum")).T.fillna(0).astype(int)
        )
        collapsed = collapsed.sum(axis=0)
        collapsed["id"] = tid
        # print(collapsed)

        # if len(df.Strand.values) == 0:
        #     continue
        # if df.Strand.values[0] == -1:
        #     print("reverse")
        #     # utrs = ["cDNA_end", "cDNA_start",
        #     #         "cDNA_end_fiveUTR", "cDNA_start_fiveUTR",
        #     #         "cDNA_end_threeUTR", "cDNA_start_threeUTR", "exon_len"]
        # else:
        #     pass

        out.append(collapsed.to_frame().T)

        """
           cDNA_start  cDNA_end  cDNA_end_fiveUTR  cDNA_start_fiveUTR  cDNA_end_threeUTR  cDNA_start_threeUTR
                0        8389         0               265                   0               8389                 7700
        """

    utrs = pd.concat(out)
    utrs["cDNA_end_fiveUTR"] = utrs["five_prime_len"]
    utrs["cDNA_end_threeUTR"] = utrs["cDNA_start_threeUTR"] + utrs["three_prime_len"]
    utr_cols = ["id", "start", "end",
                "cDNA_start_fiveUTR", "cDNA_end_fiveUTR", "five_prime_len",
                "cDNA_start_threeUTR", "cDNA_end_threeUTR", "three_prime_len", "cdna_len"]
    utrs = utrs[utr_cols]
    print(utrs.head())
    utrs.to_csv(file_out, sep="\t", index=False)


def get_cdna_seq_from_biomart():
    attrs = ["ensembl_transcript_id", "cdna", "strand"]

    # quit()

    # example: tids = ["ENST00000446046"]  # ENST00000354785  ENST00000715744
    out = []
    num = len(xpore["id"].unique())
    idx = 0
    for tid in xpore["id"].unique():
        idx += 1
        print(idx, " / ", num)
        tid = tid.split(".")[0]
        print(tid)
        df = biomart_dataset.query(attributes=attrs,
                                   filters={"link_ensembl_transcript_stable_id": tid}
                                   )
        out.append(df)
        # print(df)
    cdna = pd.concat(out)
    print(cdna.head())
    cdna = cdna.rename(columns={"Transcript stable ID": "id"})
    cdna.to_csv(file_out_cdna, sep="\t", index=False)


def get_gene_name_from_biomart():
    attrs = ["ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "external_gene_source"]

    # quit()

    # example: tids = ["ENST00000446046"]  # ENST00000354785  ENST00000715744
    out = []
    num = len(xpore["id"].unique())
    idx = 0
    for tid in xpore["id"].unique():
        idx += 1
        print(idx, " / ", num)
        tid = tid.split(".")[0]
        print(tid)
        df = biomart_dataset.query(attributes=attrs,
                                   filters={"link_ensembl_transcript_stable_id": tid}
                                   )
        out.append(df)
        # print(df)
    gene_names = pd.concat(out)
    print(gene_names.head())
    gene_names = gene_names.rename(columns={"Transcript stable ID": "id"})
    return gene_names

if __name__ == '__main__':

    chroms = list(np.arange(1, 23))
    chroms_str = [str(x) for x in chroms]
    # ensembl_data = EnsemblRelease(release=111)  # Use appropriate release number only 11 or below
    biomart_dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')

    # for f in biomart_dataset.filters.keys():
    #     print(f)
    # for k in biomart_dataset.attributes.keys():
    #     print(k)
    #
    #
    file = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_C3_modRate_05.tsv"
    file_out = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_C3_modRate_05_UTRs.tsv"
    file_out_cdna = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_C3_modRate_05_cDNA.tsv"
    file_out_cdna_utr = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_C3_modRate_05_cDNA_UTR.tsv"

    # file = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_padj_modrate_025.tsv"
    # file_out = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_padj_modrate_025_UTR.tsv"
    # file_out_cdna = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_padj_modrate_025_cDNA.tsv"
    # file_out_cdna_utr = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR.tsv"

    # file = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_C3_modRate_07.tsv"
    # # file_temp = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_C3_modRate_07.tsv"
    # file_out = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_C3_modRate_07_UTRs.tsv"
    # file_out_cdna = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_C3_modRate_07_cDNA.tsv"
    # file_out_cdna_utr = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_C3_modRate_07_cDNA_UTR.tsv"

    # xpore = pd.read_table(file)
    #
    # xpore["id"] = xpore["id"].apply(lambda x: x.split(".")[0])
    #
    # get_utrs_from_biomart()
    # get_cdna_seq_from_biomart()
    # utr_ = pd.read_table(file_out)
    # cdna_ = pd.read_table(file_out_cdna)
    # cdna_.rename(columns={"Transcript stable ID": "id"}, inplace=True)

    # map UTR and cDNA seq
    # xpore = pd.merge(xpore, cdna_, on='id', how='left')
    # xpore = pd.merge(xpore, utr_, on='id', how='left')
    # print(xpore.head())
    # xpore.to_csv(file_out_cdna_utr, sep="\t", index=False)
    # # print(xpore)

    file_for_gene_names = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS.tsv"
    file_out_genename = "direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_padj_modrate_025_cDNA_UTR_REL_POS_GENENAMES.tsv"
    xpore = pd.read_table(file_for_gene_names)
    g = get_gene_name_from_biomart()
    gene_names = pd.merge(xpore, g, on="id", how="left")
    gene_names.to_csv(file_out_genename, sep="\t", index=False)

