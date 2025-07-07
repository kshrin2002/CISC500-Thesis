import re

# Load the gene list
with open("C://Users//Shrinidhi//Desktop//Thesis//Top 100 genes//LR//top_100_genes_lr_mean.csv", "r") as f:
    lines = f.readlines()

# Remove decimals from ENSG IDs while keeping everything else the same
cleaned_lines = [re.sub(r"(ENSG\d+)\.\d+", r"\1", line) for line in lines]

# Save the cleaned gene list
with open("C://Users//Shrinidhi//Desktop//Thesis//cleaned_genes_lr.txt", "w") as f:
    f.writelines(cleaned_lines)


