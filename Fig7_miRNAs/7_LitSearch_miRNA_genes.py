import time
import matplotlib.pyplot as plt
from wordcloud import WordCloud
from pymed import PubMed as PM
import pandas as pd
import numpy as np
import json
import os
import requests
import xml.etree.ElementTree as xml
from pprint import pprint



def search_pubmed(df, file_name, tag=""):

    pubmed = PM(tool="ASAPsearch",
                email="tobias.reinberger@uni-luebeck.de")
    # https://pubmed.ncbi.nlm.nih.gov/33378036/
    out = []
    df["Author"] = df["Vorname"] + " " + df["Nachname"]
    num_publ = []

    for author in df["Author"].unique():

        print(author)

        q = f"{author} AND (atherosclerosis OR plaque OR lesion) AND (mouse OR mice OR murine OR rat)"
        results = pubmed.query(q, max_results=100)

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

            if print_article:
                print("Title:", article.title)
                # print("Authors:", ", ".join(author['lastname'] for author in article.authors if 'lastname' in author))
                # print("Journal:", article.journal)
                print("Publication Date:", article.publication_date)
                print("Abstract:", article.abstract)
                # print("DOI:", article.doi)
                print("PubMed ID:", article.pubmed_id.split("\n")[0])
        year = [int(str(x).split("-")[0]) for x in dates]
        if len(year) == 0 or min(year) < 2014:
            continue
        temp = pd.DataFrame(data={"title": titles, "year": dates,
                                  "abstract": abstracts,
                                  "pubmed_id": pubmed_id, "DOI": dois})
        temp["Author"] = author
        out.append(temp)
        print(author, idx)
        print("=" * 120 + "\n")
        num_publ.append(idx)
    df_out = pd.concat(out)
    df_out.to_excel(file_name)
    return


def get_publications_for_terms(terms, file_name):

    results = pubmed.query(terms, max_results=9999)
    idx = 0
    for article in results:
        # print(article)
        # print(article.title)
        artile_json = dict()
        pubmed_id = article.pubmed_id.split("\n")[0]
        json_out = f"MarketScreen/GlobalSearch_New/{pubmed_id}.json"
        if os.path.exists(json_out):
            continue

        # Accessing basic information from the article
        print(idx)
        abstract = article.abstract
        if article.abstract is not None and "patients" in abstract:
            continue

        year = int(str(article.publication_date).split("-")[0])
        if year < 2015:
            continue
        idx += 1
        artile_json["title"] = article.title

        if abstract is not None:
            abstract = abstract.strip()
            artile_json["abstract"] = abstract
        else:
            artile_json["abstract"] = None
        artile_json["publication_date"] = str(article.publication_date)

        artile_json["pubmed_id"] = pubmed_id

        if article.doi is not None:
            doi = article.doi.split("\n")[0]
            artile_json["doi"] = doi
        else:
            artile_json["doi"] = None

        # Get authors and their affiliations
        authors = list(article.authors)
        if len(authors) > 0:
            first_author = authors[0]
            last_author = authors[-1]
            artile_json["first_author"] = first_author
            artile_json["last_author"] = last_author

        if print_article:
            print(artile_json)
        if save_json:
            with open(json_out, "w", encoding="UTF-8") as f:
                json.dump(artile_json, f)
        print("=" * 120 + "\n")

    print("number of articles:", idx)


def xml_to_dict(element):
    if len(element) == 0:  # If no children
        return element.text.strip() if element.text else None
    result = {}
    for child in element:
        if child.tag == "AbstractText":
            result[child.tag] = child.text.strip()
        elif child.tag == "ArticleTitle":
            result[child.tag] = child.text.strip()

        elif child.tag not in result:
            print(child.tag)
            child_dict = xml_to_dict(child)
            result[child.tag] = child_dict
        else:
            # If tag already exists, convert to list
            if type(result[child.tag]) is not list:
                result[child.tag] = [result[child.tag]]
            result[child.tag].append(child_dict)
    # Include attributes
    result.update(element.attrib)
    return result


def get_pmids(term):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&api_key={api_key}&usehistory=y&retmax={limit}&term="

    term = term.replace(" ", "+")
    year_range = "\"2015\"[Date - Publication] : \"3000\"[Date - Publication]"
    year_range = year_range.replace(" ", "+")

    # receive pubmed_ids for search_terms
    r = requests.get(url + term + "+AND+" + year_range, headers={"Content-Type": "xml"})
    response = r.text
    root = xml.fromstring(response)
    # Loop over the articles and construct article objects
    idx = 0
    pubmed_ids = []
    for article in root[5]:
        idx = idx + 1
        pubmed_ids.append(article.text)
    print("Number of publications:", idx)
    return pubmed_ids


def efetch_pubmed(pubmed_ids):

    idx = 1
    for i in range(0, len(pubmed_ids), batch):
        print("batch: ", i)
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&api_key={api_key}&retmode=xml&id="
        ids = ",".join(pubmed_ids[i: i + batch])

        batch_file = f"{folder_out}/XML/batch_{i + 1}.xml"
        if os.path.exists(batch_file) and not overwrite:
            continue
        if not os.path.exists(batch_file):
            r = requests.get(url + ids, headers={"Content-Type": "xml",
                                             "user-agent": "Mozilla/5.0"})

            root = xml.fromstring(r.text)
            tree = xml.ElementTree(root)
            tree.write(batch_file,
                       xml_declaration=True,
                       encoding='utf-8')
        else:
            tree = xml.parse(batch_file)
            root = tree.getroot()

        articles = root.findall('PubmedArticle')
        for article in articles:
            article_json = dict()
            medlinecitation = article.find('MedlineCitation')
            pubmeddata = article.find("PubmedData")
            id_list = pubmeddata.find("ArticleIdList")
            id_list = dict(zip([x.attrib["IdType"] for x in id_list], [x.text for x in id_list]))
            pmid = id_list["pubmed"]
            json_file = f"{folder_out}/json/{pmid}.json"
            if os.path.exists(json_file) and not overwrite:
                continue
            date = None
            for elem in pubmeddata.find("History"):
                if elem.attrib["PubStatus"] == "pubmed":
                    date = dict(zip([d.tag for d in elem], [d.text for d in elem]))
                    break
            coi = medlinecitation.find("CoiStatement")
            coi = coi.text if coi is not None else None

            article_content = medlinecitation.find("Article")
            title = article_content.find("ArticleTitle")
            title = title.text if title is not None else None

            abstract = article_content.findall(".//AbstractText")
            abstract_text = []
            if abstract:
                for elem in abstract:
                    if elem.text is not None:
                        elem = elem.text.strip() if elem is not None else None
                        abstract_text.append(elem)
            abstract_text = " ".join(abstract_text)

            author_list = article_content.find("AuthorList")
            if author_list is None:
                continue
            authors = author_list.findall("Author")

            def get_author_dict(author):
                author_dict = dict(zip([x.tag for x in author], [x.text for x in author]))
                author_affiliation = author.findall(".//Affiliation")
                for a in author_affiliation:
                    a = a.text.strip() if a is not None else None
                    if "Electronic address" in a:
                        mail = a.split("Electronic address: ")[-1][:-1]
                        mail = mail.strip()
                        author_dict["Mail"] = mail
                        break
                    elif "@" in a:
                        temp_mail = a[:a.find("@")]
                        mail_start = temp_mail.rfind(" ")
                        mail = a[mail_start + 1:]
                        mail = mail.strip()
                        author_dict["Mail"] = mail
                        break
                author_affiliation = [x.text.split(" Electronic address")[0][:-1] for x in author_affiliation]

                # clean mail
                temp_array = []
                for temp in author_affiliation:
                    if "@" in temp:
                        temp_mail = temp[:temp.find("@")]
                        mail_start = temp_mail.rfind(" ")
                        temp_array.append(temp[:mail_start])
                    else:
                        temp_array.append(temp)

                # clean brackets
                temp_array2 = []
                for temp in temp_array:
                    if "(" in temp and ")" in temp:
                        first, second = temp.find("("), temp.find(")")
                        temp = temp.replace(temp[first: second + 1], "")
                    temp_array2.append(temp.strip().replace("\xa0", " "))

                author_affiliation = temp_array2
                author_dict["Affiliation"] = author_affiliation
                if "AffiliationInfo" in author_dict.keys():
                    author_dict.pop("AffiliationInfo")
                return author_dict

            first_author_dict = get_author_dict(authors[0])
            last_author_dict = get_author_dict(authors[-1])

            article_json["title"] = title
            article_json["abstract"] = abstract_text
            article_json["publ_date"] = date
            article_json["article_ids"] = id_list
            article_json["coi"] = coi
            article_json["first_author"] = first_author_dict
            article_json["last_author"] = last_author_dict

            with open(json_file, "w", encoding="UTF-8") as f:
                json.dump(article_json, f)

            print("col_idx: ", idx)
            print("title: ", title)
            # print("abstract: ", abstract)
            print("publ_date: ", date)
            #  print("doi: ", id_list["doi"])
            print("pubmed_ids: ", pmid)
            print("coi: ", coi)
            print("first_author: ", first_author_dict)
            print("last_author: ", last_author_dict)
            idx += 1
            # pprint(article)
            print("\n")
            time.sleep(1)


if __name__ == '__main__':

    api_key = "ee534511189e6eca5ec153cdca17278bfc08"

    # apply api_key here: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

    fetch_pub_info = True
    overwrite = True

    term = "Title/Abstract"
    out = {}
    folder_out = "miRNA_binding_PUBMED"
    os.makedirs(folder_out,  exist_ok=True)
    os.makedirs(folder_out + r"/XML", exist_ok=True)
    os.makedirs(folder_out + r"/json", exist_ok=True)
    num_ = []
    mirnas = pd.read_table("../Nanopore/direct_RNA_seq/0_TABLES/C3_vs_C4_xpore_data_mirna_binding_logcount_10_3prime_UTR_mirna_at_dmr_annotated_RNAfold_publ_n.tsv")
    for gene, mirna in zip(mirnas["gene_name"], mirnas["mirna"]):
        mirna = "-".join(mirna.split("-")[1:-1])
        query_string =f"({gene}[Text Word] AND {mirna}[Text Word]) OR ({gene}[Title/Abstract] AND {mirna}[Title/Abstract])"
        # query_string = query_string + " AND 2020:2025[dp]"
        print(query_string)
        # quit()

        limit = 9999
        pubmed_ids_ = get_pmids(query_string)
        print(n_publ := len(pubmed_ids_))
        num_.append(n_publ)
        if not fetch_pub_info:
            continue
        if len(query_string) == 0:
            continue
        batch = 10  # By including an API key, a site can post up to 10 requests per second by default.
        efetch_pubmed(pubmed_ids_)

    mirnas["Num_pub"] = num_
    mirnas.to_excel("miRNA_Gene_PUBMED.xlsx", index=False)
    """ folder = "Zebrafish"
    terms = ["Atherosclerosis[All Fields] OR atherosclerosis[All Fields]) AND (mouse[organism] OR mice[organism] OR murine[organism] OR rat[organism])",
             "Atherosclerosis[All Fields] OR atherosclerosis[All Fields] OR angiogenesis[All Fields]) AND (zebrafish[organism] OR danio rerio[organism])",
             ]"""


