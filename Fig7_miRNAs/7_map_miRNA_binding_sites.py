import os.path
import time
import re
import matplotlib.pyplot as plt
import pandas as pd
from pymed import PubMed as PM
import numpy as np
import requests
import sys


sys.path.append(r'C:\Users\tobia\PycharmProjects\LitSearch\SearchSemanticScholar')

# Set up the ProxyGenerator
# pg = ProxyGenerator()
# pg.FreeProxies()


# quit()
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def filter_kmers_in_targetscan(file_name="Resources/miR_Family_Info_kmer_filtered_conservation_score_m1.xlsx", overwrite=False):
    if os.path.exists(file_name) and not overwrite:
        return pd.read_excel(file_name)


    # xpore_data2 = xpore_data2[xpore_data2.transcript.isin(rna_fold.transcript_ID)] #filter for accessable site

    xpore_data2 = xpore_data.dropna(subset=["context_seq"])
    xpore_data2["context_seq"] = xpore_data2["context_seq"].apply(lambda s: s.replace("T", "U"))

    start = 2
    len_motif = 9  # max 7
    xpore_data2["mirna_binding_site"] = xpore_data2["context_seq"].apply(lambda s: s[start:start + len_motif])
    start = 4
    len_motif = 7  # max 7
    xpore_data2["mirna_binding_site_core"] = xpore_data2["context_seq"].apply(lambda s: s[start:start + len_motif])
    kmer_count = xpore_data2.value_counts("mirna_binding_site_core")
    kmer_count = kmer_count[kmer_count >= 2]
    print(kmer_count)


    # print(mirna_info[mirna_info["Seed+m8"].str.contains("CCAGU")])
    # quit()
    # mirna_info["Seed+m8"] = mirna_info["Seed+m8"].apply(lambda x: "U" + x)
    # reverse _complement
    complement = {"A": "U",
                  "U": "A",
                  "G": "C",
                  "C": "G"}
    # consider G-U wobble ????

    mirna_info["Seed+m8_orf_seq"] = mirna_info["Seed+m8"].apply(lambda x: "".join([complement[s] for s in x[::-1]]))
    mirna_info["Seed+m8_extend"] = mirna_info["Mature sequence"].apply(lambda x: x[:11])
    print(mirna_info)
    # https://www.targetscan.org/cgi-bin/targetscan/vert_80/view_gene.cgi?rs=ENST00000289104.4&taxid=9606&showcnc=0&shownc=0&shownc_nc=&showncf1=&showncf2=&subset=1#miR-183-5p.2
    # miRNA quantificaction
    # https://www.qiagen.com/us/product-categories/discovery-and-translational-research/pcr-qpcr-dpcr/qpcr-assays-and-instruments/mirna-qpcr-assay-and-panels?cmpid=PC_PCR_PCR_pcr-sales_1120_SEA_GA&gad_source=1&gclid=Cj0KCQjwkdO0BhDxARIsANkNcrc96nAKGd32Sb0ODNbkHdmer_H84FUn_dCOJUkEg3OJwP5WNbMbhHAaAoAREALw_wcB

    out = []
    for motif, counts in kmer_count.items():
        # print(motif)
        for kmer in xpore_data2.loc[xpore_data2["mirna_binding_site_core"] == motif, "kmer"]:
            kmer_rc = "".join([complement[x] for x in kmer[::-1]])
            temp = mirna_info[mirna_info["Seed+m8_extend"].str.contains(kmer_rc)].copy()
            # check if Seed+m8_orf_seq in context_seq
            unique_seeds = temp["Seed+m8_orf_seq"].unique()
            is_in_binding_site = []
            for seq in unique_seeds:
                if xpore_data2[xpore_data2["mirna_binding_site"].str.contains(seq)].empty:
                    is_in_binding_site.append(False)
                else:
                    is_in_binding_site.append(True)

            temp["is_in_binding_site"] = temp["Seed+m8_orf_seq"].map(dict(zip(unique_seeds, is_in_binding_site)))
            temp = temp[temp["is_in_binding_site"]]
            if not temp.empty:
                # print(motif, counts)
                # print(mirna_info.loc[mirna_info["Seed+m8"].str.match(motif), "miR family"].unique())
                # print(temp[["Seed+m8", "miR family"]])
                # temp["context_seq"] = motif
                temp["kmer"] = kmer
                temp["kmer_rc"] = kmer_rc
                temp["counts"] = counts
                out.append(temp)

    df_out = pd.concat(out, ignore_index=True)
    df_out = df_out.drop_duplicates(subset="MiRBase ID")

    # print(df_out[["kmer", "context_seq", "Seed+m8_orf_seq"]])
    # print(df_out["Family Conservation?"].unique())
    print(df_out)
    df_out = df_out[df_out["Family Conservation?"] >= -2]
    df_out.to_excel(file_name, index=False)
    # print(df_out["miR family"].unique())
    return df_out


def extract_random_mirna(file_name="Resources/miR_Family_Info_RANDOM_filtered.xlsx", overwrite=False):
    if os.path.exists(file_name) and not overwrite:
        return pd.read_excel(file_name)

    data = pd.read_table(r"Resources\miRNA\TargetScan_miR_Family_Info\miR_Family_Info.txt")
    data = data[data["Species ID"] == 9606]
    data = data[data["Family Conservation?"] >= 0]
    df_mirna = pd.read_excel("Resources/miR_Family_Info_kmer_filtered.xlsx")
    data = data[~data["miR family"].isin(df_mirna["miR family"])]
    x = np.random.randint(len(data), size=(len(df_mirna["MiRBase ID"].unique())))
    random_mirna = data.iloc[x]
    print(random_mirna)
    random_mirna.to_excel(file_name, index=False)
    return random_mirna


def search_pubmed(df, file_name, tag=""):

    print(df)
    pubmed = PM(tool="miRNA_search",
                email="tobias.reinberger@uni-luebeck.de")
    # https://pubmed.ncbi.nlm.nih.gov/33378036/
    out = []
    mirnas = df["MiRBase ID"].unique()
    num_publ = []

    for mirna in mirnas:
        mirna_temp = "-".join(mirna.split("-")[1:])
        if "." in mirna_temp:
            mirna_temp = mirna_temp[:-2]

        print(mirna_temp)
        mirna_temp2 = mirna_temp.replace("miR", "MicroRNA")
        q = f"({mirna_temp} OR {mirna_temp2}) AND smooth muscle cell"
        results = pubmed.query(q, max_results=50)
        # AND (atherosclerosis OR coronary artery disease)
        # f"AND (coronary artery disease OR atherosclerosis OR plaque) "
        # f"NOT review",

        idx = 0
        abstracts, titles, dates, pubmed_id, dois = [], [], [], [], []
        for article in results:
            idx += 1
            # Accessing basic information from the article
            print_article = True
            titles.append(article.title)
            if article.abstract is not None:
                abstracts.append(article.abstract.strip())
            else:
                abstracts.append(None)

            dates.append(article.publication_date)
            pubmed_id.append(article.pubmed_id.split("\n")[0])
            if article.doi is not None:
                dois.append(article.doi.split("\n")[0])
            else:
                dois.append(None)
            print_article = False
            if print_article:
                print("Title:", article.title)
                # print("Authors:", ", ".join(author['lastname'] for author in article.authors if 'lastname' in author))
                # print("Journal:", article.journal)
                print("Publication Date:", article.publication_date)
                # print("Abstract:", article.abstract)
                # print("DOI:", article.doi)
                print("PubMed ID:", article.pubmed_id.split("\n")[0])
        temp = pd.DataFrame(data={"title": titles, "year": dates,
                                  "abstract": abstracts,
                                  "pubmed_id": pubmed_id, "DOI": dois})
        temp["miRNA"] = mirna
        out.append(temp)
        print(mirna, idx)
        print("=" * 120 + "\n")
        num_publ.append(idx)
    df_out = pd.concat(out)
    df_out.to_excel(f"Resources/miR_Family_Info_LitResearch_{tag}.xlsx", index=False)
    df["num_publication"] = df["MiRBase ID"].map(dict(zip(mirnas, num_publ)))
    print(df["num_publication"].sum(axis=0))
    df.to_excel(file_name)
    return
    print(df)
    out = []
    mirnas = df["MiRBase ID"].unique()
    num_publ = []

    for mirna in mirnas:
        mirna_temp = "-".join(mirna.split("-")[1:])
        if "." in mirna_temp:
            mirna_temp = mirna_temp[:-2]
        print("search semantic scholar for ", mirna_temp, end="\r")

        idx = 0
        abstracts, titles, dates, urls, journals, ranks = [], [], [], [], [], []
        for a in search_query:
            if a["bib"]["venue"] != "NA":
                # print(a)
                citations = a["num_citations"]
                year = a["bib"]["pub_year"]
                if citations == "NA" or year == "NA":
                    continue
                title = a["bib"]["title"]
                gs_rank = a["gsrank"]
                journal = a["bib"]["venue"]
                url = a["pub_url"]
                snippets = a["bib"]["abstract"]

                if int(citations) > 5 or int(year) >= 2023:
                    idx += 1
                    if idx < max_results:
                        titles.append(title)
                        abstracts.append(snippets)
                        dates.append(year)
                        urls.append(url)
                        journals.append(journal)
                        ranks.append(gs_rank)
                        print("\r")
                        print("Title:", title)
                        print(url)

        temp = pd.DataFrame(data={"miRNA": [mirna] * len(titles),
                                  "title": titles, "year": dates, "gs_rank": ranks,
                                  "abstract": abstracts,
                                  "journal": journals, "URL": urls})
        out.append(temp)
        print(mirna_temp, idx)
        print("=" * 120 + "\n")
        num_publ.append(idx)
    df_out = pd.concat(out)
    df_out.to_excel(f"Resources/miR_Family_Info_GS_Research_{tag}.xlsx", index=False)
    df["num_publication"] = df["MiRBase ID"].map(dict(zip(mirnas, num_publ)))
    print(df["num_publication"].sum(axis=0))
    df.to_excel(file_name)


def search_google_scholar(df=None, file_name="out.xlsx", tag="", max_results=10):

    user_agent = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                  "(KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36")
    headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_5)\
                  AppleWebKit/537.36 (KHTML, like Gecko) Cafari/537.36'}
    # mirna_temp = "miR-199a-5p"
    """googleTrendsUrl = 'https://google.com'
         response = requests.get(googleTrendsUrl)
         if response.status_code == 200:
             g_cookies = response.cookies.get_dict()"""
    # q = f'intitle:"{mirna_temp}" OR (intext:"{mirna_temp}" AROUND (30) intext:"smooth-muscle-cell")'
    # q = q.replace(" ", "%20")
    # q = "https://scholar.google.de/scholar?hl=de&as_sdt=0%2C5&as_ylo=2010&as_yhi=2024&as_vis=1&q=f%27intitle%3A%22miR-199a-5p%22+OR+%28intext%3A%22miR-199a-5p%22+AROUND+%2830%29+intext%3A%22smooth-muscle-cell%22%29%27&btnG="
    g_example = "https://scholar.google.com/citations?hl=en&view_op=search_authors&mauthors=james%20watson"
    resp = requests.get(g_example,
                        # headers=headers,
                        # cookies=g_cookies
                        )
    # print(resp.text)
    if resp.status_code != 200:
        print(f"Request failed with {resp.status_code} because {resp.reason}")

    print(df)
    out = []
    mirnas = df["MiRBase ID"].unique()
    num_publ = []

    for mirna in mirnas:
        mirna_temp = "-".join(mirna.split("-")[1:])
        if "." in mirna_temp:
            mirna_temp = mirna_temp[:-2]
        print("search google scholar for ", mirna_temp, end="\r")
        search_query = scholarly.search_pubs(
            f'intitle:"{mirna_temp}" '
            f'OR (intext:"{mirna_temp}" AROUND (30) intext:"smooth-muscle-cell")',
            patents=False,
            year_low=2014,
            citations=False)
        time.sleep(10)

        idx = 0
        abstracts, titles, dates, urls, journals, ranks = [], [], [], [], [], []
        for a in search_query:
            if a["bib"]["venue"] != "NA":
                # print(a)
                citations = a["num_citations"]
                year = a["bib"]["pub_year"]
                if citations == "NA" or year == "NA":
                    continue
                title = a["bib"]["title"]
                gs_rank = a["gsrank"]
                journal = a["bib"]["venue"]
                url = a["pub_url"]
                snippets = a["bib"]["abstract"]

                if int(citations) > 5 or int(year) >= 2023:
                    idx += 1
                    if idx < max_results:
                        titles.append(title)
                        abstracts.append(snippets)
                        dates.append(year)
                        urls.append(url)
                        journals.append(journal)
                        ranks.append(gs_rank)
                        print("\r")
                        print("Title:", title)
                        print(url)

        temp = pd.DataFrame(data={"miRNA": [mirna] * len(titles),
                                  "title": titles, "year": dates, "gs_rank": ranks,
                                  "abstract": abstracts,
                                  "journal": journals, "URL": urls})
        out.append(temp)
        print(mirna_temp, idx)
        print("=" * 120 + "\n")
        num_publ.append(idx)
    df_out = pd.concat(out)
    df_out.to_excel(f"Resources/miR_Family_Info_GS_Research_{tag}.xlsx", index=False)
    df["num_publication"] = df["MiRBase ID"].map(dict(zip(mirnas, num_publ)))
    print(df["num_publication"].sum(axis=0))
    df.to_excel(file_name)
    return


def search_semantic_scholar(df=None, file_name="out.xlsx", tag="", max_results=10):
    print(df)
    out = []
    mirnas = df["MiRBase ID"].unique()
    num_publ = []

    for mirna in mirnas:
        mirna_temp = "-".join(mirna.split("-")[1:])
        if "." in mirna_temp:
            mirna_temp = mirna_temp[:-2]
        print("search google scholar for ", mirna_temp, end="\r")
        # q = f'"{mirna_temp}" + ("smooth muscle cell" | "coronary artery disease" | "atherosclerosis")'
        q = f'"{mirna_temp}" + "smooth muscle"'
        # print(q)
        search = SearchPaper(query_key=mirna_temp, query=q, limit=1000)
        results = search.find_basis_paper()

        idx = 0
        abstracts, titles, dates, urls, journals, ranks = [], [], [], [], [], []
        for a in results:
            # print(a)

            # print(a)
            citations = a["citationCount"]
            year = a["year"]
            if citations == "NA" or year == "NA":
                pass
            title = a["title"]
            journal = a["journal"]
            url = a["url"]
            snippets = a["abstract"]

            if int(citations) >= 1 or int(year) >= 2000:
                idx += 1
                if idx < max_results:
                    titles.append(title)
                    abstracts.append(snippets)
                    dates.append(year)
                    urls.append(url)
                    journals.append(journal)
                    print(title, url)

        temp = pd.DataFrame(data={"miRNA": [mirna] * len(titles),
                                  "title": titles, "year": dates,
                                  "abstract": abstracts,
                                  "journal": journals, "URL": urls})
        out.append(temp)
        print(mirna_temp, idx)
        print("\n")
        print("=" * 120 + "\n")
        num_publ.append(idx)
    df_out = pd.concat(out)
    df_out.to_excel(f"Resources/miR_Family_Info_SemanticScholar_Research_{tag}.xlsx", index=False)
    df["num_publication"] = df["MiRBase ID"].map(dict(zip(mirnas, num_publ)))
    print(df["num_publication"].sum(axis=0))
    df.to_excel(file_name)


if __name__ == '__main__':
    # https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2014.00023/full
    # Common features of microRNA target prediction tools
    # https://library.acg.edu/how-to-guides/google-scholar/advanced-searching
    # https://serpapi.com/organic-results

    """
    G-U wobble pairs are frequently observed in miRNA seed region binding (nucleotides 2–8 of the miRNA), 
    though perfect Watson-Crick pairing is more common and generally leads to stronger binding.
    Outside the seed region (e.g., at the 3′ end of the miRNA), G-U wobbles are much more tolerated and often present.  
    """

    complement = {"A": "U",
                  "U": "A",
                  "G": "C",
                  "C": "G"}

    # load miRNA expression data
    mirna_file = (r"..\miRNA_Expression/Snakemake_report_alle_Konditionen/results_condition/"
                  r"tables/logcount-matrix/model_X.logcount-matrix.tsv")
    mirna_expr = pd.read_table(mirna_file, index_col="gene")
    mirna_expr = mirna_expr.drop("transcript", axis=1)
    mirna_expr = mirna_expr.fillna(0)

    log_count = 10

    mirna_expr = mirna_expr[mirna_expr.max(axis=1) >= log_count].dropna()
    print(len(mirna_expr))

    # load mirnaBase
    mirna_info = pd.read_table(r"Resources\miRNA\TargetScan_miR_Family_Info\miR_Family_Info.txt")
    mirna_info = mirna_info[mirna_info["Species ID"] == 9606]
    mirna_info = mirna_info[mirna_info["MiRBase ID"].isin(mirna_expr.index)]

    mirna_info_unique_seed = mirna_info.drop_duplicates(subset="Seed+m8")

    # load xpore_data data
    filter_name = "../Nanopore/direct_RNA_seq/Xpore/DATA/xpore_cdna.all_group_compare_2025/diffmod_s3-s4/majority_direction_kmer_diffmod_padj_modrate_025_GENOME_cDNA_UTR.tsv"
    xpore_data = pd.read_table(filter_name)

    cdna_seq = pd.read_table("../0_REF_TABLES/biomart_cdna_seq_transcript_id.tsv")
    xpore_data["cdna_seq"] = xpore_data["id"].map(dict(zip(cdna_seq["Transcript stable ID"],
                                                            cdna_seq["cDNA sequences"])))
    xpore_data = xpore_data.dropna(subset="cdna_seq")
    print(xpore_data[xpore_data.cdna_seq.isna()])

    # filter 3' UTR if needed
    xpore_data = xpore_data[xpore_data.relative_cdna_pos >= 2]

    take_random_positions = False
    if take_random_positions:
        save_file = f"C3_vs_C4_xpore_data_mirna_binding_logcount_{log_count}_RANDOM.tsv"
        if os.path.exists("xpore_random_pos.tsv"):
            xpore_data = pd.read_table("xpore_random_pos.tsv")
        else:
            xpore_data["position"] = xpore_data.apply(lambda x:
                                                      np.random.choice(np.array([
                                                          t_pos for t_pos in range(10, len(x["cdna_seq"]) - 10)
                                                          if t_pos not in list(np.arange(x["position"] - 10,
                                                                                         x["position"] + 10))
                                                      ])),
                                                      axis=1)
            # pos2 = xpore_data["position"].values
            # plt.scatter(pos1, pos2)
            # plt.show()
            print("CAUTION !!!!! >>> random positions created")
            xpore_data[["id", "position", "kmer", "z_score", "mod_assignment", "cdna_seq"]].to_csv("xpore_random_pos.tsv",
                                                                                          sep="\t", index=False)
    else:
        save_file = f"C3_vs_C4_xpore_data_mirna_binding_logcount_{log_count}.tsv"  # _3prime_UTR

    window = 25
    left = int(np.floor(window / 2))
    right = window - left
    xpore_data["context_seq"] = xpore_data.apply(lambda x: x["cdna_seq"][int(x["position"]) - left:
                                                                    int(x["position"]) + right], axis=1)
    xpore_data["context_seq"] = xpore_data["context_seq"].apply(lambda x: x.replace("T", "U"))
    xpore_data["kmer_rc"] = xpore_data["kmer"].apply(lambda x: "".join([complement[b] for b in x][::-1]))
    xpore_data["context_seq_rc"] = xpore_data["context_seq"].apply(lambda x: "".join([complement[b] for b in x][::-1]))
    xpore_data["context_seq_rc_seed"] = xpore_data["context_seq_rc"].apply(lambda x: x[6:-5])

    # modbase in seed 7mer-m8
    out = []
    kmer_binding_mirnas_out = []
    mirnas_temp = []

    for (t, pos, kmer, kmer_rc, seq, seed_seq, context_seq), tbl in xpore_data.groupby(["id", "position", "kmer",
                                                                               "kmer_rc",
                                                                               "context_seq_rc",
                                                                               "context_seq_rc_seed", "context_seq"]):
        # print(t, pos, kmer)
        mirna_found = False
        binding_mirnas = []
        kmer_binding_mirnas = []
        for seed, mirna in zip(mirna_info["Seed+m8"], mirna_info["MiRBase ID"]):
            if seed in seq:
                binding_mirnas.append(mirna)
                mirna_found = True
            if seed in seed_seq:
                kmer_binding_mirnas.append(mirna)
                kmer_binding_mirnas_out.append([t, pos, kmer, kmer_rc, context_seq, seed_seq, mirna, seed])
                mirnas_temp.append(mirna)
        if mirna_found:
            tbl["context_binding_mirnas"] = ", ".join(binding_mirnas)
            tbl["kmer_binding_mirnas"] = ", ".join(kmer_binding_mirnas)
            out.append(tbl)
    xpore_data_mirna = pd.concat(out)
    xpore_data_mirna_dmr = pd.DataFrame(kmer_binding_mirnas_out, columns=["id", "position", "kmer",
                                                                          "kmer_rc",
                                                                          "context_seq",
                                                                          "context_seq_rc_seed",
                                                                          "mirna",
                                                                          "Seed+m8"
                                                                          ])

    xpore_data_mirna = xpore_data_mirna[["id", "position", "kmer", # "z_score",
                                         "mod_assignment", "context_seq",
                                         "kmer_rc", "context_seq_rc", "context_seq_rc_seed",
                                         "context_binding_mirnas", "kmer_binding_mirnas"]]
    print(xpore_data_mirna_dmr.columns)
    xpore_data_mirna_dmr["pos_modbase_in_seed"] = xpore_data_mirna_dmr.apply(lambda x:
                                                                    abs(x["context_seq_rc_seed"].find(x["Seed+m8"]) - 6) + 2,
                                                                                                            axis=1)
    print(xpore_data_mirna_dmr)

    xpore_data_mirna_dmr.to_csv(save_file.replace(".tsv", "_mirna_at_dmr.tsv"), sep="\t", index=False)
    xpore_data_mirna.to_csv(save_file, sep="\t", index=False)
    print(xpore_data_mirna)
    mirnas = np.unique(mirnas_temp)
    counts = {}
    for mi in mirnas:
        counts[mi] = mirnas_temp.count(mi)
    mirna_counts = pd.DataFrame(data={"miRNA": counts.keys(), "Counts": counts.values()})
    print(len(mirnas))
    mirna_counts = mirna_counts[mirna_counts["Counts"] > 5]

    print(mirna_counts.sort_values("Counts", ascending=False))

    quit()
    # rna_fold = pd.read_csv(
    #     r"C:\Users\tobia\PyCharmWetLab\WetLab\RNAseq\Nanopore\direct_RNA_seq\miRNA_binding_RNAfold\RNAfolds\170bp_window\0_sumstats\seq_params_best_5_combined_cDNA.csv",
    #     usecols=[12, 13])
    # rna_fold = rna_fold[rna_fold["kmer_in_structure"].str.contains("hhhhh")]

    df_selected = filter_kmers_in_targetscan(overwrite=True, file_name="expressed_mirna_in_kmer_pos.xlsx")
    print(df_selected)
    quit()

    # check Literature
    df_rand = extract_random_mirna(overwrite=True)
    files = ["Resources/miR_Family_Info_kmer_filtered.xlsx",
             "Resources/miR_Family_Info_RANDOM_filtered.xlsx"]
    search_semantic_scholar(df_selected, file_name=files[0], tag="kmer_select")
    search_semantic_scholar(df_rand, file_name=files[1], tag="RANDOM")
    # search_google_scholar(df_selected, file_name=files[0], tag="kmer_select")
    # search_google_scholar(df_rand, file_name=files[1], tag="RANDOM")
    # search_pubmed(df=df_selected, file_name=files[0], tag="kmer_select")
    # search_pubmed(df=df_rand, file_name=files[1], tag="RANDOM")
