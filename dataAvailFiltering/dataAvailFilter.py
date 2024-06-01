#this script removes all rows without supplementary data 

import pandas as pd

# Read the CSV file into a pandas DataFrame
df = pd.read_csv("dataAvailFiltering\pmc_results.csv")

# Filter out rows with empty SupplementaryDatasets column
df_filtered = df.dropna(subset=['SupplementaryDatasets'])

# Save the filtered DataFrame to a new CSV file
df_filtered.to_csv("dataAvailFiltering/filtered_file.csv", index=False)