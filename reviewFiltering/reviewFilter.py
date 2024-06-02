import pandas as pd
import os

print("Current working directory:", os.getcwd())

try:
    df = pd.read_csv('reviewFiltering/file.csv')
    print("File loaded successfully")
except FileNotFoundError:
    print("File not found. Please check the file path.")
    exit()

print("Initial DataFrame:")
print(df.head())

filtered_df = df[~df['PublicationType'].str.contains('Review', na=False)]

print("Filtered DataFrame:")
print(filtered_df.head())

filtered_df.to_csv('reviewFiltering/filtered_file.csv', index=False)
print("Filtered file saved successfully")

review_df = df[df['PublicationType'].str.contains('Review', na=False)]

print("Review DataFrame:")
print(review_df.head())

review_df.to_csv('reviewFiltering/review_file.csv', index=False)
print("Review file saved successfully")