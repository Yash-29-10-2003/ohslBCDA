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
## Imoproving the time required per request.
For my laptop (r5 5600h , gtx 1650) the approx times per abstract summary were like this :
- ~ 12 seconds per abstract for the facebook bart-large-cnn on CPU.
- ~ 1.8 seconds per abstract for the facebook bart-large-cnn on GPU.
- ~ 7.8 seconds per abstract for the distilbart-cnn-12-6 on CPU. (A smaller and faster model. Trades of summary quality)
- ~ 1.1 seconds per abstract for the distilbart-cnn-12-6 on CPU. 

### If you have a dedicated GPU , you can get pytorch-cuda support to have faster processing by :
- Check if you already have it by running the following code in your v-env :

```
pubmedScanner/cudaTest.py
```
- If the output is true , you already have it . If not , follow the following steps :
1. Download cuda 12.1 from [https://developer.nvidia.com/cuda-12-1-0-download-archive] for your devices' specification.
2. Download the cudNN library for cuda from [https://developer.nvidia.com/cudnn] . (Recommended for dnn models.)
3. The default address for the downloaded cuda should be C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1 . If its not , find your path and use that in the upcioming steps wherever ive added my path. Add the path C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\bin in your pc's environment variables . [System Properties > Env. Variables ] . Set CUDA_HOME as the variable name and C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\bin as the value.
4. If you're using the virtual environment , also add the path in your virtual environments path . For vsc powershell use the following command :

```
$env:PATH = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\bin;$env:PATH"
```
5. To check youve followed all the steps correctly , put in the following command : nvcc --version in your venv command terminal and ur pc's default command terminal. If it doesnt tell you about your build like following you've made some mistake in following the steps , refer to a forum or make sure youve followed the steps correctly.

```
(venv) PS D:\Coding\ohslBCDA> nvcc --version
nvcc: NVIDIA (R) Cuda compiler driver
Copyright (c) 2005-2023 NVIDIA Corporation
Built on Wed_Feb__8_05:53:42_Coordinated_Universal_Time_2023
Cuda compilation tools, release 12.1, V12.1.66
Build cuda_12.1.r12.1/compiler.32415258_0
```
6. If everything works fine upto this point , just run the following command in your terminal to download the cuda support for pytorch (Might need to change the command a bit for different versions of cuda and pytorch):
```
(venv) PS D:\Coding\ohslBCDA> pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
```
7. After that check weather youve followed all the steps correctly by reruning the cudaTest.py script by :
```
python pubmedScanner/test.py
```
8. When running the pubmed scanner after this , the model should use your gpu to process the abstracts , saving time by a lot !

## How it works :

### Code follows the following steps :
1. Search: Searches the PubMed database using a specified query and retrieves the total number of results .
2. Fetch Details: Retrieves detailed information for all the papers matching the query , it only feteches results from 2020 and afterwards.
3. Process: Processes the fetched papers, extracting the title, authors, url to the paper , paper type and abstract for each paper.
4. Runs the abstract through a nlp model (huggingface : Facebook bart)
5. Save to CSV: Saves the extracted information to a CSV file named "pubmed_results.csv" .

## Miscellaneous :
- Since the code extracts ALL the related papers from pubmed and sorts them on the basis of relevance , the lower we might go the less relevant the papers might become , filter the data as required !

- Following is an example of how the resulting csv would look after the input : (((Breast) OR (Breast cancer)) AND (primary prevention)) AND (oxidative stress) 

![image](https://github.com/Yash-29-10-2003/ohslBCDA/assets/89728102/a3dceaba-2449-4bee-9547-8cd1633fe210)


Everytime we would put in a new query the code would update the pubmed-results.csv so save a copy before running another command .



