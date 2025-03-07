import pandas as pd
import numpy as np

# FPKM data
fpkm_df = pd.read_csv(r'C:\Users\Shrinidhi\Desktop\Test\TcgaTargetGtex_rsem_gene_fpkm\TcgaTargetGtex_rsem_gene_fpkm',  sep='\t', index_col=0)  # Adjust file path as needed

# Phenotypic data
phenotype_df = pd.read_csv(r'C:\Users\Shrinidhi\Desktop\Test\TcgaTargetGTEX_phenotype.txt\TcgaTargetGTEX_phenotype.txt', 
                             sep='\t')

# Make sure the samples match
common_samples = fpkm_df.columns.intersection(phenotype_df['sample'])

fpkm_df_filtered = fpkm_df[common_samples]
phenotype_df_filtered = phenotype_df[phenotype_df['sample'].isin(common_samples)]

# Cancer status column to label
phenotype_df_filtered['cancer_status'] = phenotype_df_filtered['sample_type'].apply(lambda x: 'cancer' if 'Tumor' in x else 'normal')

print(phenotype_df_filtered['cancer_status'].value_counts())

# Order of samples
phenotype_df_filtered = phenotype_df_filtered.set_index('sample').loc[fpkm_df_filtered.columns]

# Add the column to the RNA sequenced gene expression data
fpkm_with_metadata = fpkm_df_filtered.T
fpkm_with_metadata['cancer_status'] = phenotype_df_filtered['cancer_status'].values

# Save file
fpkm_with_metadata.to_csv('fpkm_with_metadata.csv', sep='\t')

