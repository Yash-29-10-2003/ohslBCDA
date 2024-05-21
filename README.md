# OHSL BCDA

This repository contains all the models and assets for the OHSL BCDA project.

# Pubmed Scanner:

- Original Base Article the model is based on: [Scraping Big Data from Public Research Repositories (e.g., PubMed, arXiv)](https://medium.com/@kliang933/scraping-big-data-from-public-research-repositories-e-g-pubmed-arxiv-2-488666f6f29b)
- This code is an extension of it, with the following additional features:

1. Allows the user to add a search term for personalized results, such as (((Breast) OR (Breast cancer)) AND (primary prevention)) AND (oxidative stress).
2. Saves the results as a CSV file for better accessibility.
3. Utilizes pagination to process a large number of papers efficiently.
4. Displays links, publication types, and key findings based on the abstract using huggingface models.
5. Limits the output to papers published in 2020 and afterward for relevancy.

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
    pip install biopython
    pip install transformers torch
    ```
5. Run the app:
    ```bash
    python pubmedScanner/scraperUsingAPI.py
    ```
    Note: The app may download a ~1.5 GB file to store the summarizing model in local memory.

6. Input the search parameter when prompted.

## Improving Processing Time per Request:

- The approximate processing times per abstract as processed on my laptop(r5 5600h , gtx 1650) for different models  are as follows:
    - ~12 seconds per abstract for the Facebook BART-Large-CNN on CPU.
    - ~1.8 seconds per abstract for the Facebook BART-Large-CNN on GPU.
    - ~7.8 seconds per abstract for the DistilBART-CNN-12-6 on CPU (a smaller and faster model).
    - ~1.1 seconds per abstract for the DistilBART-CNN-12-6 on CPU.

### Using CUDA Support for PyTorch [dedicated gpu required]:

1. Check if CUDA support is already available:
    ```bash
    python pubmedScanner/cudaTest.py
    ```
2. If CUDA support is not available, follow these steps:
    - Download CUDA 12.1 from the [NVIDIA website](https://developer.nvidia.com/cuda-12-1-0-download-archive).
    - Download the cuDNN library for CUDA from the [NVIDIA website](https://developer.nvidia.com/cudnn).
    - Set the `CUDA_HOME` environment variable to `C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\bin`.
    - Add `C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\bin` to the system's environment variables.
    - Update the virtual environment's path:
        ```bash
        $env:PATH = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\bin;$env:PATH"
        ```
    - Install PyTorch with CUDA support:
        ```bash
        pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
        ```
3. Check if CUDA support is now available:
    ```bash
    python pubmedScanner/test.py
    ```
Running the scraper after successfuly following these steps would use the gpu for saving time.
If you dont have a GPU and also have a weak CPU and each paper is taking too much time to be summarized , you can use a lighter model by changing the model in ln68 of the scraper (trades of accuracy):
    ```bash
    summarizer = pipeline('summarization', model='sshleifer/distilbart-cnn-12-6', device=0 if torch.cuda.is_available() else -1)
    ```

## Workings:

1. **Search**: Searches the PubMed database using a specified query and retrieves the total number of results.
2. **Fetch Details**: Retrieves detailed information for all the papers matching the query, fetching only those published in 2020 and afterward.
3. **Process**: Processes the fetched papers, extracting the title, authors, URL to the paper, publication type, and abstract for each paper.
4. **Summarize**: Runs the abstract through an NLP model (Hugging Face: Facebook BART) to extract key findings. Uses the computers gpu if followed the steps above to save time .
5. **Save to CSV**: Saves the extracted information to a CSV file named "pubmed_results.csv".

## Miscellaneous:

- Since the code extracts all related papers from PubMed and sorts them by relevance, the further down the results, the less relevant the papers might become. Filter the data as required.
- After each query, the resulting CSV file will be updated, so save a copy before running another command.
- Following is an example of how the output might look for the following input :
   ```bash
     (((Breast) OR (Breast cancer)) AND (primary prevention)) AND (oxidative stress) 
   ```
![image](https://github.com/Yash-29-10-2003/ohslBCDA/assets/89728102/80797f8f-bbad-4a46-a3e1-e144189c07ba)

