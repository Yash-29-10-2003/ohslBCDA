# ohslBCDA
All my models/assets for OHSL BCDA project.


# Pubmed Scanner :

Original Base Article the model is based on : https://medium.com/@kliang933/scraping-big-data-from-public-research-repositories-e-g-pubmed-arxiv-2-488666f6f29b

## How to use :
Either download the entire repo as a zip file or fork and clone it.
Steps to forking and cloning : [https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo]

### Once downloaded :
- Create a virtual environment : python -m venv venv

- acitivate the environment : venv\Scripts\activate

- download required dependency : pip install biopython

- run the app : python pubmedScanner/scraperUsingAPI.py

- Input the search parameter : Ex : (((Breast) OR (Breast cancer)) AND (primary prevention)) AND (lifestyle behavior)

## How it works :

#Code follows the following steps :
1. Search: Searches the PubMed database using a specified query and retrieves the total number of results.
2. Fetch Details: Retrieves detailed information for all the papers matching the query.
3. Process: Processes the fetched papers, extracting the title, authors, and abstract for each paper.
4. Save to CSV: Saves the extracted information to a CSV file named "pubmed_results.csv" with columns for title, authors, and abstract.

#Following is an example of how the resulting csv would look after the input : (((Breast) OR (Breast cancer)) AND (primary prevention)) AND ("oxidative stress")

![image](https://github.com/Yash-29-10-2003/ohslBCDA/assets/89728102/e6583a43-92b8-4afd-9c89-83975e95e409)

