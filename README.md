# ohslBCDA
All my models/assets for OHSL BCDA project.


# Pubmed Scanner :

- Original Base Article the model is based on : https://medium.com/@kliang933/scraping-big-data-from-public-research-repositories-e-g-pubmed-arxiv-2-488666f6f29b
- This code is an extension of it having the following extra features :

1. Lets the user add a search term like : (((Breast) OR (Breast cancer)) AND (primary prevention)) AND (oxidative stress) , for personalized results.
2. Saves the result as a csv for better accessibility.
3. Uses pagination to be able to process a large number of papers.
4. Also displays links , publication type and key findings based on the abstract.
5. Limits the output to be of 2020 and afterwards for relevancy.


## How to use :
Either download the entire repo as a zip file or fork and clone it.
Steps to forking and cloning : [https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo]

### Once downloaded :
- Create a virtual environment :
```
python -m venv venv
```
- acitivate the environment :
```
venv\Scripts\activate
```
- download required dependencies :

```
pip install biopython
pip install transformers torch
```

- run the app :
```
python pubmedScanner/scraperUsingAPI.py  
```
[the app may download a ~1.5 gb file to store the summarizing model in local memory]

- Input the search parameter :
```
 Ex : (((Breast) OR (Breast cancer)) AND (primary prevention)) AND (lifestyle behavior)
```
## How it works :

### Code follows the following steps :
1. Search: Searches the PubMed database using a specified query and retrieves the total number of results .
2. Fetch Details: Retrieves detailed information for all the papers matching the query , it only feteches results from 2020 and afterwards.
3. Process: Processes the fetched papers, extracting the title, authors, url to the paper , paper type and abstract for each paper.
4. Runs the abstract through a nlp model (huggingface : Facebook bart)
5. Save to CSV: Saves the extracted information to a CSV file named "pubmed_results.csv" .

- Since the code extracts ALL the related papers from pubmed and sorts them on the basis of relevance , the lower we might go the less relevant the papers might become , filter the data as required !

- Following is an example of how the resulting csv would look after the input : (((Breast) OR (Breast cancer)) AND (primary prevention)) AND (oxidative stress) 

![image](https://github.com/Yash-29-10-2003/ohslBCDA/assets/89728102/a29d883f-d1db-45d1-8095-e1e37c202795)

Everytime we would put in a new query the code would update the pubmed-results.csv so save a copy before running another command .



