import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the DAVID functional annotation table
file_path = "C://Users//Shrinidhi//Desktop//DAVID_mean_table.txt"
df = pd.read_csv(file_path, sep="\t")

# Select relevant columns
df_filtered = df[["Category", "Term", "PValue"]].dropna()

# Convert PValue to -log10 scale for better visualization
df_filtered["-log10(PValue)"] = -np.log10(df_filtered["PValue"])

# Select top 10 terms for visualization
top_terms = df_filtered.nsmallest(10, "PValue")  

# Plot horizontal bar graph
plt.figure(figsize=(12, 8))  # Adjust size
plt.barh(top_terms["Term"], top_terms["-log10(PValue)"], color="royalblue")
plt.xlabel("Pathways", fontsize=14)  # X-axis label size
plt.ylabel("-log10(P-value)", fontsize=14)  # Y-axis label size
plt.title("DAVID Pathway Enrichment", fontsize=16)  # Title size
plt.xticks(fontsize=12)  # X-axis tick labels
plt.yticks(fontsize=12) 
plt.xticks(rotation=45, ha="right")  # For bar plots  
plt.yticks(rotation=0)  # Ensure y-axis text is readable
plt.title("Top Enriched GO Terms & Pathways from DAVID")
plt.gca().invert_yaxis()  # Flip the order so the smallest p-value is at the top
plt.tight_layout()  # Automatically adjusts spacing
plt.show()
