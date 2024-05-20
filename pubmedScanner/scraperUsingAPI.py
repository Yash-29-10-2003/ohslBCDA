import csv
import time
import torch          #to run model on gpu
from Bio import Entrez
from transformers import pipeline

# Checking for GPU availability
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

def search(query):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax=0,  # Requesting total count
                            retmode='xml',
                            term=query,
                            mindate='2020/01/01')  # Data filter
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
                                term=query,
                                mindate='2020/01/01')
        results = Entrez.read(handle)
        id_list.extend(results['IdList'])
    Entrez.email = 'your.email@example.com'
    handle = Entrez.efetch( db='pubmed',
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

def get_publication_type(paper):
    pub_types = []
    if 'PublicationTypeList' in paper['MedlineCitation']['Article']:
        pub_types = [pub_type for pub_type in paper['MedlineCitation']['Article']['PublicationTypeList']]
    return ', '.join(pub_types)

def get_url(paper):
    pmid = paper['MedlineCitation']['PMID']
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

def extract_key_findings(abstract, summarizer):
    summary = summarizer(abstract, max_length=150, min_length=40, do_sample=False)
    return summary[0]['summary_text']

if __name__ == '__main__':
    # Initializing the summarizer
    summarizer = pipeline('summarization', model='facebook/bart-large-cnn', device=0 if torch.cuda.is_available() else -1)
    
    query = input("Enter the query ! : ")
    
    # Performing the search to get the total number of results
    total_results = search(query)
    
    # Fetching details for all papers
    papers = fetch_details(query, total_results)
    
    # Processing the fetched papers and save them to a CSV file
    with open('pubmed_results.csv', 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Title', 'Authors', 'PublicationType', 'URL', 'KeyFindings', 'Abstract']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        total_time = 0
        num_papers = len(papers['PubmedArticle'])
        
        for paper in papers['PubmedArticle']:
            title = paper['MedlineCitation']['Article']['ArticleTitle']
            try:
                author_list = paper['MedlineCitation']['Article']['AuthorList']
                authors = ', '.join([author.get('LastName', '') for author in author_list])
            except KeyError:
                authors = 'N/A'
            abstract = get_abstract(paper)
            publication_type = get_publication_type(paper)
            url = get_url(paper)
            
            start_time = time.time()
            key_findings = extract_key_findings(abstract, summarizer)
            end_time = time.time()
            
            total_time += (end_time - start_time)
            
            writer.writerow({'Title': title, 'Authors': authors, 'PublicationType': publication_type, 'URL': url, 'KeyFindings': key_findings, 'Abstract': abstract})
        
        average_time = total_time / num_papers
        print(f"Processed {num_papers} papers in {total_time:.2f} seconds. Average time per paper: {average_time:.2f} seconds.")
