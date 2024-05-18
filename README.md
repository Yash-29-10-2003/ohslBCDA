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


