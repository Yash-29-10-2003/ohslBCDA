import csv
from Bio import Entrez

def search(query):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax=0,  # Set retmax to 0 to get the total count of results
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    total_results = int(results['Count'])
    return total_results

def fetch_details(query, total_results):
    batch_size = 100  # Number of results to fetch per request
    id_list = []
    for start in range(0, total_results, batch_size):
        Entrez.email = 'your.email@example.com'
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax=batch_size,
                                retstart=start,
                                retmode='xml',
                                term=query)
        results = Entrez.read(handle)
        id_list.extend(results['IdList'])
    Entrez.email = 'your.email@example.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=','.join(id_list))
    papers = Entrez.read(handle)
    return papers

def get_abstract(paper):
    abstract = ''
    if 'Abstract' in paper['MedlineCitation']['Article']:
        abstract = paper['MedlineCitation']['Article']['Abstract']['AbstractText']
        if isinstance(abstract, list):
            abstract = ' '.join(abstract)
    return abstract

if __name__ == '__main__':
    # Provide your desired query here
    query = input("Enter the query ! : ")
    
    # Perform the search to get the total number of results
    total_results = search(query)
    
    # Fetch details for all papers
    papers = fetch_details(query, total_results)
    
    # Process the fetched papers and save them to a CSV file
    with open('pubmed_results.csv', 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Title', 'Authors', 'Abstract']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for paper in papers['PubmedArticle']:
            title = paper['MedlineCitation']['Article']['ArticleTitle']
            try:
                author_list = paper['MedlineCitation']['Article']['AuthorList']
                authors = ', '.join([author.get('LastName', '') for author in author_list])
            except KeyError:
                authors = 'N/A'
            abstract = get_abstract(paper)
            
            writer.writerow({'Title': title, 'Authors': authors, 'Abstract': abstract})
