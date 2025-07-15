import streamlit as st
import pandas as pd
import os
import tempfile
from txtgz_to_csv_converter import GeneExpressionConverter

st.set_page_config(page_title="TXT.GZ to CSV Converter", layout="centered")
st.title("🧬 GEO .txt.gz Expression File Converter")

st.markdown("""
Upload a `.txt.gz` expression matrix (like from NCBI GEO),
and this tool will clean it into a tabular gene × sample format.
""")

uploaded_file = st.file_uploader("📁 Upload your `.txt.gz` file", type=["gz"])

if uploaded_file:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".txt.gz") as temp_input:
        temp_input.write(uploaded_file.read())
        input_path = temp_input.name

    converter = GeneExpressionConverter()

    with st.spinner("🧪 Processing file..."):
        output_path = input_path.replace(".txt.gz", ".tsv")
        result_path = converter.process_file(input_path, output_path)

    if result_path:
        st.success("✅ File converted successfully!")
        df = pd.read_csv(result_path, sep="\t")
        st.write("🔍 Preview of converted file:")
        st.dataframe(df.head())

        with open(result_path, "rb") as f:
            st.download_button(
                label="⬇️ Download cleaned matrix",
                data=f,
                file_name=os.path.basename(result_path),
                mime="text/tab-separated-values"
            )
    else:
        st.error("❌ Conversion failed. Please check the input format.")
