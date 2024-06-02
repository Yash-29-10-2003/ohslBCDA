#this file is just an experimental file to get supplementary files to work 


# test id - 11137110


import requests
import xml.etree.ElementTree as ET
import csv

def fetch_supplementary_materials(article_id):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pmc",
        "id": article_id,
        "retmode": "xml"
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return response.content
    else:
        print("Error fetching supplementary materials")
        return None

def parse_supplementary_materials(article_details):
    root = ET.fromstring(article_details)
    articles = []
    for i, article in enumerate(root.findall(".//sec[@sec-type='supplementary-material']/sec")):
        title_elem = article.find("title")
        title = title_elem.text.strip() if title_elem is not None and title_elem.text is not None else ""  # Extract title text if available
        link_elems = article.findall(".//media")
        datasets = []
        for link_elem in link_elems:
            dataset_link = link_elem.get("{http://www.w3.org/1999/xlink}href", "")
            dataset_id = link_elem.get("id", "")
            datasets.append({"ID": dataset_id, "Link": dataset_link})
        articles.append({"Title": title, "SupplementaryDatasets": datasets})
    return articles

def save_to_csv(articles, filename="pmcScraper/prototypes/supplementary_materials.csv"):
    with open(filename, "w", newline='', encoding='utf-8') as csvfile:
        fieldnames = ["Title", "SupplementaryDatasets"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for article in articles:
            writer.writerow(article)

def main():
    article_id = input("Enter the article ID: ")

    article_details = fetch_supplementary_materials(article_id)
    if article_details:
        articles = parse_supplementary_materials(article_details)
        save_to_csv(articles)
        print("Supplementary materials saved to supplementary_materials.csv")
    else:
        print("No supplementary materials found for the given article ID.")

if __name__ == "__main__":
    main()