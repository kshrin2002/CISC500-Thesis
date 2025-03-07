import re

# Load the gene list
with open("C://Users//Shrinidhi//Desktop//Thesis//david_genes.txt", "r") as f:
    lines = f.readlines()

# Remove decimals from ENSG IDs while keeping everything else the same
cleaned_lines = [re.sub(r"(ENSG\d+)\.\d+", r"\1", line) for line in lines]

# Save the cleaned gene list
with open("C://Users//Shrinidhi//Desktop//Thesis//cleaned_david_genes_dnn.txt", "w") as f:
    f.writelines(cleaned_lines)

print("✅ Cleaned gene list saved as 'cleaned_david_genes_dnn.txt'. Upload this to DAVID!")
