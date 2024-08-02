import numpy as np
import gzip
import streamlit as st
import pandas as pd

@st.cache_data
def read_gene_embeddings(file_path):
    gene_embeddings = {}
    with gzip.open(file_path, 'rt') as file:
        next(file)
        for line in file:
            parts = line.strip().split()
            gene = parts[0]
            embedding = np.array([float(x) for x in parts[1:]])
            gene_embeddings[gene] = embedding
    return gene_embeddings

@st.cache_data
def load_gene_id_mapping(file_path):
    mapping_df = pd.read_csv(file_path, sep='\t', header=None, names=['ensembl_id', 'gene_symbol'])
    ensembl_to_symbol = dict(zip(mapping_df['ensembl_id'], mapping_df['gene_symbol']))
    symbol_to_ensembl = {v: k for k, v in ensembl_to_symbol.items() if v != k}
    return ensembl_to_symbol, symbol_to_ensembl

def get_gene_id(gene, gene_embeddings, ensembl_to_symbol, symbol_to_ensembl):
    if gene in gene_embeddings:
        return gene
    elif gene in symbol_to_ensembl and symbol_to_ensembl[gene] in gene_embeddings:
        return symbol_to_ensembl[gene]
    elif gene in ensembl_to_symbol and ensembl_to_symbol[gene] in gene_embeddings:
        return ensembl_to_symbol[gene]
    else:
        return None

def add_google_analytics(measurement_id):
    ga_script = f"""
    <!-- Global site tag (gtag.js) - Google Analytics -->
    <script async src="https://www.googletagmanager.com/gtag/js?id={measurement_id}"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){{dataLayer.push(arguments);}}
      gtag('js', new Date());

      gtag('config', '{measurement_id}');
    </script>
    """
    st.components.v1.html(ga_script, height=0)