import numpy as np
import pandas as pd

# Data import
fpkm_data = pd.read_csv(r'C:\Users\Shrinidhi\Desktop\Test\fpkm_with_metadata.csv', sep='\t', index_col=0)

# Print the original shape of the data
print("Original data shape:", fpkm_data.shape)  # Should be (19131, 60499)

# Separate the cancer_status column so it won't be affected by the numeric filtering
cancer_status = fpkm_data['cancer_status']

# Drop the 'cancer_status' column temporarily to focus on numeric gene data
fpkm_data_numeric = fpkm_data.drop(columns=['cancer_status'])

# Transpose the data so that rows represent genes and columns represent samples
transposed_data = fpkm_data_numeric.T  # (60499, 19131)

# Print the transposed data shape
print("Transposed data shape (genes as rows):", transposed_data.shape)  # Should be (60499, 19131)

### 1. Variance-based Gene Selection ###
# Compute variance for each gene (row-wise since genes are now in rows)
variances = transposed_data.var(axis=1)

# Select top 3000 genes based on variance
top_genes_variance = variances.nlargest(3000).index

# Debug: print the selected top genes based on variance to check
print("Selected top genes based on variance:", top_genes_variance)

# Subset the original (non-transposed) data to include only these top genes
fpkm_data_top_variance = fpkm_data_numeric[top_genes_variance]

# Add back the cancer_status column to the filtered dataset
fpkm_data_top_variance['cancer_status'] = cancer_status

# Print the shape of the final dataset
print("Variance-selected data shape (with cancer_status):", fpkm_data_top_variance.shape)  # Should be (19131, 3001)

# Save the variance-selected subset to a file
fpkm_data_top_variance.to_csv('fpkm_top_3000_variance_genes.tsv', sep='\t')


### 2. Mean-based Gene Selection ###
# Compute mean for each gene (row-wise since genes are now in rows)
means = transposed_data.mean(axis=1)

# Select top 3000 genes based on mean expression
top_genes_mean = means.nlargest(3000).index

# Debug: print the selected top genes based on mean to check
print("Selected top genes based on mean:", top_genes_mean)

# Subset the original (non-transposed) data to include only these top genes
fpkm_data_top_mean = fpkm_data_numeric[top_genes_mean]

# Add back the cancer_status column to the filtered dataset
fpkm_data_top_mean['cancer_status'] = cancer_status

# Print the shape of the final dataset
print("Mean-selected data shape (with cancer_status):", fpkm_data_top_mean.shape)  # Should be (19131, 3001)

# Save the mean-selected subset to a file
fpkm_data_top_mean.to_csv('fpkm_top_3000_mean_genes.tsv', sep='\t')


### 3. Coefficient of Variance (CV)-based Gene Selection ###
# Compute standard deviation and mean for each gene
std_dev = transposed_data.std(axis=1)
mean_expr = transposed_data.mean(axis=1)

# Compute coefficient of variance (CV) as std_dev / mean
cv = std_dev / mean_expr

# Select top 3000 genes based on CV
top_genes_cv = cv.nlargest(3000).index

# Debug: print the selected top genes based on CV to check
print("Selected top genes based on CV:", top_genes_cv)

# Subset the original (non-transposed) data to include only these top genes
fpkm_data_top_cv = fpkm_data_numeric[top_genes_cv]

# Add back the cancer_status column to the filtered dataset
fpkm_data_top_cv['cancer_status'] = cancer_status

# Print the shape of the final dataset
print("CV-selected data shape (with cancer_status):", fpkm_data_top_cv.shape)  # Should be (19131, 3001)

# Save the CV-selected subset to a file
fpkm_data_top_cv.to_csv('fpkm_top_3000_cv_genes.tsv', sep='\t')

print("Gene selection completed for variance, mean, and CV.")
