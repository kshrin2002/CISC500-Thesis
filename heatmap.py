import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import glob

# Set path and sort files for consistent processing
file_paths = glob.glob("C://Users//Shrinidhi//Desktop//Thesis//Top 100 genes//**/*.csv", recursive=True)
file_paths.sort()

# Separate and sort files by type
mean_files = sorted([f for f in file_paths if "_mean" in f])
covar_files = sorted([f for f in file_paths if "covar" in f or "_cv" in f])
var_files = sorted([f for f in file_paths if "_var" in f])

# Preprocess SHAP files with genes as column headers
def preprocess_shap_file(file):
    df = pd.read_csv(file)
    if df.columns[0].startswith("ENSG"):
        df = df.T.reset_index()
        df.rename(columns={"index": "Gene"}, inplace=True)
    if "cancer_status" in df.columns:
        df.drop(columns=["cancer_status"], inplace=True)
    return df

# Stable gene selection + presence matrix
def process_gene_files(files):
    gene_presence = {}
    all_genes = set()

    for file in files:
        df = preprocess_shap_file(file) if "SHAP" in file else pd.read_csv(file)
        gene_column = df.columns[0]
        genes = df[gene_column].unique()
        all_genes.update(genes)

        method_name = file.split("/")[-1].replace(".csv", "")
        for gene in genes:
            gene_presence.setdefault(gene, {})[method_name] = 1

    gene_df = pd.DataFrame.from_dict(gene_presence, orient="index").fillna(0)

    all_genes_df = pd.DataFrame(index=sorted(list(all_genes)))
    gene_df = all_genes_df.join(gene_df, how="left").fillna(0)

    # Stable: sort by frequency, then alphabetically
    gene_df["_sum"] = gene_df.sum(axis=1)
    gene_df["_gene"] = gene_df.index

    # Sort by frequency (_sum), then gene name (_gene) for deterministic output
    gene_df = gene_df.sort_values(by=["_sum", "_gene"], ascending=[False, True])
    gene_df = gene_df.drop(columns=["_sum", "_gene"])

    # Return just the top 20 rows
    top_genes = gene_df.head(20)
    return top_genes

# Clean column labels like "DNN + SHAP"
def simplify_x_labels(columns):
    cleaned = []
    for col in columns:
        col_lower = col.lower()
        model = (
            "SVM" if "svm" in col_lower else
            "DNN" if "dnn" in col_lower else
            "LR" if "lr" in col_lower else
            "Unknown"
        )
        if "shap" in col_lower:
            xai = "SHAP"
        elif "lime" in col_lower:
            xai = "LIME"
        elif "deeplift" in col_lower:
            xai = "DeepLIFT"
        elif "coef" in col_lower:
            xai = "Coefficient"
        elif "pfi" in col_lower or "perm" in col_lower:
            xai = "Permutation"
        else:
            xai = "Other"
        cleaned.append(f"{model} + {xai}")
    return cleaned

# Plot heatmap with a properly aligned vertical arrow

def plot_heatmap(data, title):
    data = data.copy()
    data.index = [f"Sample {i+1}" for i in range(len(data))]
    data.columns = simplify_x_labels(data.columns)

    fig, ax = plt.subplots(figsize=(16, 10))
    sns.heatmap(
        data,
        cmap="Blues",
        annot=True,
        fmt=".0f",
        linewidths=0.5,
        cbar=False,
        ax=ax
    )

    # Add vertical arrow to the far right of the figure, not overlapping heatmap
    fig.subplots_adjust(right=0.85)  # Make space on the right
    arrow_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])  # Position manually outside plot
    arrow_ax.axis("off")
    arrow_ax.annotate(
        "",
        xy=(0.5, 1.0),
        xytext=(0.5, 0.0),
        arrowprops=dict(facecolor="black", edgecolor="black", arrowstyle="-|>", lw=3),
        annotation_clip=False
    )
    arrow_ax.text(
        1.5, 0.5,
        "Higher gene\nimportance",
        va="center",
        ha="left",
        rotation=90,
        fontsize=12,
        fontweight="bold"
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)

    ax.set_xlabel("Model + XAI Method", fontsize=12)
    ax.set_ylabel("Samples", fontsize=12)
    ax.set_title(title, fontsize=14)
    plt.show()

# Process and plot all sets
mean_df = process_gene_files(mean_files)
covar_df = process_gene_files(covar_files)
var_df = process_gene_files(var_files)
print("Top 20 genes (mean-based):")
print(mean_df.index.tolist())

print("\nTop 20 genes (covariance-based):")
print(covar_df.index.tolist())

print("\nTop 20 genes (variance-based):")
print(var_df.index.tolist())

plot_heatmap(mean_df, "Gene Presence Across XAI Methods (Mean-Based Genes)")
plot_heatmap(covar_df, "Gene Presence Across XAI Methods (Covariance-Based Genes)")
plot_heatmap(var_df, "Gene Presence Across XAI Methods (Variance-Based Genes)")
