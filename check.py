import pandas as pd 

fpkm_with_metadata =  pd.read_csv(r'C:\Users\Shrinidhi\Desktop\Test\fpkm_with_metadata.csv', sep='\t', index_col=0)
print(fpkm_with_metadata.head())
print(fpkm_with_metadata.columns)
# Check the count of normal vs cancer tissues
cancer_status_counts = fpkm_with_metadata['cancer_status'].value_counts()

# Print the results
print(cancer_status_counts)
