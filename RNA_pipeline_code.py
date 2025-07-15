import os
import gzip
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# --------------------------
# Step 1: Load Expression Data
# --------------------------
def load_data(file_path):
    """
    Load expression data from supported file types: CSV, TSV, FASTA, or tab-separated .txt.gz.
    Automatically handles malformed .txt.gz files by skipping problematic lines.
    """
    if file_path.endswith(".csv") or file_path.endswith(".tsv"):
        sep = "," if file_path.endswith(".csv") else "\t"
        df = pd.read_csv(file_path, sep=sep, index_col=0)
        print(f"✅ Count matrix loaded: {df.shape[0]} genes × {df.shape[1]} samples")
        return df

    elif file_path.endswith(".fasta") or file_path.endswith(".fa"):
        sequences = list(SeqIO.parse(file_path, "fasta"))
        print("✅ FASTA file loaded:")
        for record in sequences:
            print(f"  - {record.id}: {len(record.seq)} bases")
        return sequences

    elif file_path.endswith(".txt.gz"):
        print(f"✅ Attempting to load .txt.gz file: {file_path}")
        try:
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t', index_col=0, on_bad_lines='skip')
                print(f"✅ Matrix loaded: {df.shape[0]} genes × {df.shape[1]} samples")
                return df
        except Exception as e:
            print("❌ Failed to load .txt.gz as a tab-separated matrix:")
            print(str(e))
            raise ValueError("The .txt.gz file is not in a valid matrix format.")

    else:
        raise ValueError("Unsupported file format. Use .csv, .tsv, .fasta, or .txt.gz")


# --------------------------
# Step 2: Perform DGE
# --------------------------
def perform_dge(counts, sample_groups):
    """
    Perform differential gene expression using Welch's t-test.
    Returns a DataFrame with Log2 Fold Change, p-values, adjusted p-values, and significance flags.
    """

    # Validate sample names match
    if not set(sample_groups.index).issubset(counts.columns):
        print("❌ ERROR: Some sample names in 'sample_groups' are not found in the dataset!")
        print("🧪 Available columns:", list(counts.columns))
        print("📋 Provided group labels:", list(sample_groups.index))
        return pd.DataFrame()  # Return empty safely

    group1, group2 = sample_groups.unique()
    samples1 = sample_groups[sample_groups == group1].index
    samples2 = sample_groups[sample_groups == group2].index

    print(f"✅ Performing DGE between '{group1}' ({len(samples1)} samples) and '{group2}' ({len(samples2)} samples)")

    results = []

    for gene in counts.index:
        try:
            vals1 = counts.loc[gene, samples1].astype(float)
            vals2 = counts.loc[gene, samples2].astype(float)
            t_stat, p_val = ttest_ind(vals1, vals2, equal_var=False)
            log2fc = np.log2(vals2.mean() + 1) - np.log2(vals1.mean() + 1)
            results.append({"Gene": gene, "Log2FC": log2fc, "p-value": p_val})
        except Exception as e:
            continue

    if not results:
        print("⚠️ No valid results found. All rows may have failed during processing.")
        return pd.DataFrame()

    dge_df = pd.DataFrame(results)
    dge_df["adj_p-value"] = multipletests(dge_df["p-value"], method="fdr_bh")[1]
    dge_df["Significant"] = dge_df["adj_p-value"] < 0.05
    return dge_df


# --------------------------
# Step 3: Volcano Plot
# --------------------------
def plot_volcano(dge_df, output_folder="output"):
    """
    Generate a volcano plot from DGE results and save it as an image.
    """
    os.makedirs(output_folder, exist_ok=True)
    dge_df["-log10(p-value)"] = -np.log10(dge_df["p-value"])

    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=dge_df,
        x="Log2FC",
        y="-log10(p-value)",
        hue="Significant",
        alpha=0.7,
        palette={True: "red", False: "grey"}
    )
    plt.axhline(-np.log10(0.05), linestyle="--", color="black")
    plt.axvline(1, linestyle="--", color="blue")
    plt.axvline(-1, linestyle="--", color="blue")
    plt.title("Volcano Plot")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-log10(p-value)")
    plt.tight_layout()

    plot_path = os.path.join(output_folder, "volcano_plot.png")
    plt.savefig(plot_path)
    print(f"📊 Volcano plot saved to {plot_path}")

# --------------------------
# Main Script
# --------------------------
if __name__ == "__main__":
    input_file = "C:/Users/acer/Downloads/RNA_Pipeline/data/GSE255990_genes_counts_edit_240210.csv"
    data = load_data(input_file)

    if isinstance(data, pd.DataFrame):
        # 🔧 Customize your sample groups based on your data
        sample_groups = pd.Series({
        "UoE": "urban",
        "UoF": "urban",  # ✅ was previously written as "Uof"
        "UoH": "urban",
        "UoL": "urban",
        "UoN": "urban",
        "UoO": "urban",
        "HyA": "hybrid",
        "HyB": "hybrid",
        "HyD": "hybrid",
        "HyI": "hybrid",
        "HyK": "hybrid",
        "HyX": "hybrid",
        "UgJ": "rural",
        "UgP": "rural",
        "UgQ": "rural",
        "UgR": "rural",
        "UgS": "rural",
        "UgT": "rural"
    })
        # Filter only the samples from the two groups you want to compare
        sample_groups = sample_groups[sample_groups.isin(["urban", "rural"])]




        dge = perform_dge(data, sample_groups)

        if dge.empty:
            print("⚠️ DGE results are empty — check your sample group names and data format.")
            exit()

        plot_volcano(dge)



        os.makedirs("output", exist_ok=True)
        dge.to_csv("output/dge_results.csv", index=False)
        print("📁 DGE results saved to output/dge_results.csv")

        plot_volcano(dge)

    elif isinstance(data, list):
        print("🔬 FASTA file loaded — expression analysis not applicable.")
