import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import io

st.set_page_config(page_title="RNA-Seq Volcano Plot", layout="centered")

st.title("🧬 RNA-Seq Volcano Plot Generator")

# Upload expression matrix
uploaded_file = st.file_uploader("📁 Upload Expression Matrix (.csv)", type=["csv"])

if uploaded_file:
    try:
        df = pd.read_csv(uploaded_file, index_col=0)
        st.success(f"✅ Matrix loaded: {df.shape[0]} genes × {df.shape[1]} samples")
        st.dataframe(df.head())
    except Exception as e:
        st.error("❌ Failed to read the file. Please ensure it's a valid count matrix.")
        st.stop()

    # Sample name input
    samples = df.columns.tolist()
    st.subheader("🧪 Define Sample Groups")

    group1_samples = st.multiselect("Select samples for Group 1 (e.g., Urban)", samples)
    group2_samples = st.multiselect("Select samples for Group 2 (e.g., Rural)", samples)

    group1_label = st.text_input("Label for Group 1", value="Group1")
    group2_label = st.text_input("Label for Group 2", value="Group2")

    if st.button("Run DGE Analysis"):
        if not group1_samples or not group2_samples:
            st.warning("⚠️ Please select at least one sample per group.")
            st.stop()

        # Create sample group map
        all_samples = group1_samples + group2_samples
        sample_groups = pd.Series(
            {s: group1_label for s in group1_samples} | {s: group2_label for s in group2_samples}
        )

        # Perform DGE
        results = []
        for gene in df.index:
            try:
                vals1 = df.loc[gene, group1_samples].astype(float)
                vals2 = df.loc[gene, group2_samples].astype(float)
                t_stat, p_val = ttest_ind(vals1, vals2, equal_var=False)
                log2fc = np.log2(vals2.mean() + 1) - np.log2(vals1.mean() + 1)
                results.append({"Gene": gene, "Log2FC": log2fc, "p-value": p_val})
            except:
                continue

        if not results:
            st.error("❌ No valid results found.")
            st.stop()

        dge_df = pd.DataFrame(results)
        dge_df["adj_p-value"] = multipletests(dge_df["p-value"], method="fdr_bh")[1]
        dge_df["Significant"] = dge_df["adj_p-value"] < 0.05
        dge_df["-log10(p-value)"] = -np.log10(dge_df["p-value"])

        # Sort and show top up/down-regulated genes
        st.subheader("📈 Top Differentially Expressed Genes")

        top_up = dge_df[dge_df["Log2FC"] > 0].sort_values("adj_p-value").head(10)
        top_down = dge_df[dge_df["Log2FC"] < 0].sort_values("adj_p-value").head(10)

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("🔺 **Top 10 Upregulated Genes**")
            st.dataframe(top_up[["Gene", "Log2FC", "adj_p-value"]])

        with col2:
            st.markdown("🔻 **Top 10 Downregulated Genes**")
            st.dataframe(top_down[["Gene", "Log2FC", "adj_p-value"]])


        # Save DGE results
        csv_buffer = io.StringIO()
        dge_df.to_csv(csv_buffer, index=False)
        csv_download = csv_buffer.getvalue()

                # Show volcano plot
                # Volcano Plot with Top Genes Labeled
        st.subheader("📊 Enhanced Volcano Plot")

        fig, ax = plt.subplots(figsize=(10, 6))
        sns.scatterplot(
            data=dge_df,
            x="Log2FC",
            y="-log10(p-value)",
            hue="Significant",
            palette={True: "red", False: "grey"},
            alpha=0.6,
            ax=ax
        )

        # Highlight and label top genes
        top_genes = pd.concat([top_up, top_down])
        for _, row in top_genes.iterrows():
            ax.text(row["Log2FC"], row["-log10(p-value)"], row["Gene"],
                    fontsize=8, ha='right' if row["Log2FC"] < 0 else 'left',
                    color="black")

        ax.axhline(-np.log10(0.05), linestyle="--", color="black")
        ax.axvline(1, linestyle="--", color="blue")
        ax.axvline(-1, linestyle="--", color="blue")

        ax.set_title("Volcano Plot (Top Genes Highlighted)")
        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("-log10(p-value)")
        st.pyplot(fig)


        # Download
        st.download_button("⬇️ Download DGE Results", data=csv_download, file_name="dge_results.csv", mime="text/csv")

