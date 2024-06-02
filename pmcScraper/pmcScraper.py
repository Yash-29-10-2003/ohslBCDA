import requests
import xml.etree.ElementTree as ET
import csv
import pandas as pd
import time

# Function to fetch PubMed Central IDs for articles matching the query, handling pagination if necessary
def fetch_pubmed_central_data(query, max_results=500):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    id_list = []
    for start in range(0, max_results, 100):
        params = {
            "db": "pmc",
            "term": query,
            "retstart": start,
            "retmax": min(100, max_results - start),
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

# Function to fetch detailed information for articles given a list of PubMed Central IDs
def fetch_article_details(id_list):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    articles = []
    batch_size = 100
    for i in range(0, len(id_list), batch_size):
        batch_ids = id_list[i:i+batch_size]
        params = {
            "db": "pmc",
            "id": ",".join(batch_ids),
            "retmode": "xml"
        }
        response = requests.get(base_url, params=params)
        if response.status_code == 200:
            articles.append(response.content)
        else:
            print(f"Error fetching article details for batch starting at index {i}")
            time.sleep(5)  # Wait for 5 seconds before retrying
            response = requests.get(base_url, params=params)
            if response.status_code == 200:
                articles.append(response.content)
            else:
                print(f"Failed again fetching article details for batch starting at index {i}")
    return articles

# Function to parse article details and extract necessary information
def parse_article_details(article_details, id_list):
    articles = []
    for content in article_details:
        root = ET.fromstring(content)
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
            articles.append({"Title": title, "Link": link, "ArticleType": article_type, "Abstract": abstract, "SupplementaryDatasets": supplementary_datasets})
    return articles

# Function to check for availability of supplementary materials
def fetch_supplementary_materials(article_id, retries=3):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pmc",
        "id": article_id,
        "retmode": "xml"
    }
    for attempt in range(retries):
        response = requests.get(base_url, params=params)
        if response.status_code == 200:
            return parse_supplementary_materials(response.content, article_id)
        else:
            print(f"Error fetching supplementary materials for article {article_id}, attempt {attempt + 1}")
            time.sleep(2)  # Wait for 5 seconds before retrying
    return []

# Function to parse supplementary materials
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

# Function to save articles to a CSV file
def save_to_csv(articles, filename="pmcScraper/pmc_results.csv"):
    with open(filename, "w", newline='', encoding='utf-8') as csvfile:
        fieldnames = ["Title", "Link", "ArticleType", "Abstract", "SupplementaryDatasets"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for article in articles:
            article['SupplementaryDatasets'] = ", ".join(article['SupplementaryDatasets'])
            writer.writerow(article)

# Main function to handle the workflow
def main():
    query = input("Enter your query: ")
    max_results = 500  # Fetch 500 results

    print(f"Fetching data for query: {query}")
    id_list = fetch_pubmed_central_data(query, max_results)
    print(f"Found {len(id_list)} articles.")

    if id_list:
        article_details = fetch_article_details(id_list)
        if article_details:
            articles = parse_article_details(article_details, id_list)
            save_to_csv(articles)
            print("Results saved to pmc_results.csv")

            df = pd.read_csv("pmcScraper/pmc_results.csv")
            # Filtering out rows with empty SupplementaryDatasets column
            df_filtered = df.dropna(subset=['SupplementaryDatasets'])
            df_filtered.to_csv("pmcScraper/filtered_file.csv", index=False)
            print("Filtered file saved successfully")

    else:
        print("No articles found.")

if __name__ == "__main__":
    main()
