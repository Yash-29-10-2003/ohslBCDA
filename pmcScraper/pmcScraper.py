import requests
import xml.etree.ElementTree as ET      #lightweight and efficient API for parsing and creating XML data
import csv

#fetching total number of results for a particular search term (used to later determine the number of outputs)
def fetch_total_results(query):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pmc",
        "term": query,
        "retmax": 1,
        "mindate": "2020/01/01"
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        root = ET.fromstring(response.content)
        count_elem = root.find('.//Count')
        total_results = int(count_elem.text) if count_elem is not None else 0
        return total_results
    else:
        print("Error fetching total results from PubMed Central")
        return 0

#fetches PubMed Central IDs for articles matching the query, handling pagination if necessary
def fetch_pubmed_central_data(query, total_results, max_results_per_batch=100):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    id_list = []
    for start in range(0, total_results, max_results_per_batch):
        params = {
            "db": "pmc",
            "term": query,
            "retstart": start,
            "retmax": min(max_results_per_batch, total_results - start),
            "mindate": "2020/01/01",
            "sort": "relevance"
        }
        response = requests.get(base_url, params=params)
        if response.status_code == 200:
            root = ET.fromstring(response.content)
            id_list.extend([id_elem.text for id_elem in root.findall('.//Id')])
        else:
            print("Error fetching data from PubMed Central")
            break
    return id_list

#fetches detailed information for articles given a list of PubMed Central IDs
def fetch_article_details(id_list):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pmc",
        "id": ",".join(id_list),
        "retmode": "xml"
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return response.content
    else:
        print("Error fetching article details from PubMed Central")
        return None

#parsing all the details like title , abstract etc. and making sure to not be error prone
def parse_article_details(article_details, id_list):
    root = ET.fromstring(article_details)
    articles = []
    for i, article in enumerate(root.findall(".//article")):
        title_elem = article.find(".//article-title")
        title = title_elem.text.strip() if title_elem is not None and title_elem.text is not None else ""
        abstract_elem = article.find(".//abstract/p")
        abstract = abstract_elem.text.strip() if abstract_elem is not None and abstract_elem.text is not None else ""
        link = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{id_list[i]}"
        article_type = None
        article_meta_elem = article.find(".//article-meta")
        if article_meta_elem is not None:
            pub_type_elem = article_meta_elem.find(".//article-categories/subj-group/subject")
            if pub_type_elem is not None:
                article_type = pub_type_elem.text.strip()
        supplementary_datasets = fetch_supplementary_materials(id_list[i])
        articles.append({"Title": title, "Link": link, "ArticleType": article_type , "Abstract": abstract, "SupplementaryDatasets": supplementary_datasets})
    return articles

#function to check for availability of supp mats
def fetch_supplementary_materials(article_id):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pmc",
        "id": article_id,
        "retmode": "xml"
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return parse_supplementary_materials(response.content, article_id)
    else:
        print("Error fetching supplementary materials")
        return []

#function to make links for supp mats
def parse_supplementary_materials(article_details, article_id):
    root = ET.fromstring(article_details)
    datasets = []
    for article in root.findall(".//sec[@sec-type='supplementary-material']/sec"):
        link_elems = article.findall(".//media")
        for link_elem in link_elems:
            dataset_link = link_elem.get("{http://www.w3.org/1999/xlink}href", "")
            if dataset_link:
                full_link = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{article_id}/bin/{dataset_link.split('/')[-1]}"
                datasets.append(full_link)
    return datasets


#function to save files to csv
def save_to_csv(articles, filename="pmcScraper/pmc_results.csv"):
    with open(filename, "w", newline='', encoding='utf-8') as csvfile:
        fieldnames = ["Title", "Link", "ArticleType", "Abstract", "SupplementaryDatasets"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for article in articles:
            article['SupplementaryDatasets'] = ", ".join(article['SupplementaryDatasets'])
            writer.writerow(article)

def main():
    #Having a max result amount so that the cpu doesnt go haywire processing bad or general search terms
    #Afterwards , adding everything to the csv
    query = input("Enter your query: ")
    max_results = 500

    print(f"Fetching total number of results for query: {query}")
    total_results = fetch_total_results(query)
    print(f"Total results found: {total_results}")

    results_to_fetch = min(max_results, total_results)
    print(f"Fetching {results_to_fetch} results.")
    
    id_list = fetch_pubmed_central_data(query, results_to_fetch)
    print(f"Found {len(id_list)} articles.")

    if id_list:
        article_details = fetch_article_details(id_list)
        if article_details:
            articles = parse_article_details(article_details, id_list)
            save_to_csv(articles)
            print("Results saved to pmc_results.csv")
    else:
        print("No articles found.")

if __name__ == "__main__":
    main()
