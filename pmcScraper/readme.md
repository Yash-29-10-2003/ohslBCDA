# Pubmed Central Scraper 

I have made this pmc scraper to get associated data from articles using a search term , the script searches pmc for articles based on a search terms and retrieves a csv of title , publication type , abstract and supplementary material. It also generates an additional filtered csv containing only the rows which have associated data in them (for better accessibility).

## Features:

1. Allows the user to add a search term for personalized results, such as (((Breast) OR (Breast cancer)) AND (primary prevention)) AND (oxidative stress).
2. Saves the results as a CSV file for better accessibility.
3. Displays titles, links, publication types and associated data for the search.
4. Limits the output to papers published in 2020 and afterward for relevancy.

NOTE : Since PMC does aggressive rate limiting , I have made the script so that it retries for associated data after getting rate limited at max of 3 times . This ensures the best data collection but in turn takes quite a bit of time even on google colabs resources. On Colab , it takes around 7-8 minutes on the CPU and 5-6 minutes on the t4 gpu.


## How to Use:

1. Clone or download the entire repository as a ZIP file. Guide for forking and cloning : [https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo]
2. Set up a virtual environment:
    ```bash
    python -m venv venv
    ```
3. Activate the environment:
    ```bash
    venv\Scripts\activate
    ```
4. Install the required dependencies:
    ```bash
    pip install pandas
    ```
5. Run the app:
    ```bash
    python pubmedScanner/pmcScraper.py
    ```

6. Input the search parameter when prompted.

Or , if you dont want to use it in your local repository , you can access it on [google colab.](https://colab.research.google.com/drive/19JFpvmGIIu8P7hIYlRAcT83R4OwX1bAl?usp=sharing)

## Miscellaneous:

- Since the code extracts all related papers from PubMed Central and sorts them by relevance, the further down the results, the less relevant the papers might become. Filter the data as required.
- After each query, the resulting CSV files will be updated, so save a copy before running another command.
- Following is an example of how the output might look for the following input :
   ```bash
     (((Breast) OR (Breast cancer)) AND (primary prevention)) AND (oxidative stress) 
   ```
![image](https://github.com/Yash-29-10-2003/ohslBCDA/assets/89728102/f9773a2f-9ae9-42ae-9e20-01968e6ba34f)
