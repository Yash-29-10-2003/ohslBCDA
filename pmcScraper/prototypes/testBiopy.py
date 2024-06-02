#this is just a test file to see weather pmc would also work like pubmed , it doesnt D:


import requests

def fetch_paper_content(article_id):
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
        print("Error fetching paper content")
        return None

# Example article ID
article_id = "11137110"

# Fetch the content of the paper
paper_content = fetch_paper_content(article_id)

if paper_content:
    print("Paper content:")
    print(paper_content.decode("utf-8"))  # Decode bytes to string for printing
else:
    print("Failed to fetch paper content")