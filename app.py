import streamlit as st
import numpy as np
from scipy.spatial.distance import cosine
from sklearn.manifold import TSNE
import plotly.express as px
import pandas as pd
import umap
import gzip
# Function to read gene embeddings
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

# Function to find similar genes
def find_closest_genes(target_gene, gene_embeddings, n=20):
    target_embedding = gene_embeddings[target_gene]
    similarities = [(gene, 1 - cosine(embedding, target_embedding)) 
                    for gene, embedding in gene_embeddings.items() 
                    if gene != target_gene]
    return sorted(similarities, key=lambda x: x[1], reverse=True)[:n]

# Function to calculate gene similarity
def calculate_similarity(gene1, gene2, gene_embeddings):
    return 1 - cosine(gene_embeddings[gene1], gene_embeddings[gene2])

# Function to perform gene calculation
def gene_calculation(gene_a, gene_b, gene_c, gene_embeddings):
    result_vector = (gene_embeddings[gene_b] - gene_embeddings[gene_a] + gene_embeddings[gene_c])
    similarities = [(gene, 1 - cosine(embedding, result_vector)) 
                    for gene, embedding in gene_embeddings.items()]
    return sorted(similarities, key=lambda x: x[1], reverse=True)[0]

# Load gene embeddings
file_path = "data/GeneRAIN-vec.200d.txt.gz"
gene_embeddings = read_gene_embeddings(file_path)

# Streamlit app
st.title("Gene Embedding Analysis")

# Sidebar for selecting functionality
function = st.sidebar.selectbox("Select Function", 
                                ["Similar Genes", "Visualization", "Calculator", "Computing Similarity"])

if function == "Similar Genes":
    st.header("Find Similar Genes")
    target_gene = st.text_input("Enter a gene name:")
    if target_gene:
        if target_gene in gene_embeddings:
            similar_genes = find_closest_genes(target_gene, gene_embeddings)
            st.write(f"Genes similar to {target_gene}:")
            for gene, similarity in similar_genes:
                st.write(f"{gene}: {similarity:.4f}")
            
            # Visualization
            st.subheader("Visualization of Similar Genes")
            genes_to_plot = [target_gene] + [gene for gene, _ in similar_genes[:19]]
            embeddings_to_plot = np.array([gene_embeddings[gene] for gene in genes_to_plot])
            
            method = st.radio("Select visualization method:", ["t-SNE", "UMAP"])
            
            if method == "t-SNE":
                tsne = TSNE(n_components=2, random_state=42)
                embeddings_2d = tsne.fit_transform(embeddings_to_plot)
            else:
                reducer = umap.UMAP(random_state=42)
                embeddings_2d = reducer.fit_transform(embeddings_to_plot)
            
            df = pd.DataFrame({
                'x': embeddings_2d[:, 0],
                'y': embeddings_2d[:, 1],
                'gene': genes_to_plot,
                'is_target': [gene == target_gene for gene in genes_to_plot]
            })
            
            fig = px.scatter(df, x='x', y='y', text='gene', color='is_target',
                             color_discrete_map={True: 'red', False: 'blue'},
                             title=f"{method} Visualization of Similar Genes")
            fig.update_traces(textposition='top center')
            st.plotly_chart(fig)
        else:
            st.error(f"Gene {target_gene} not found in the embeddings.")

elif function == "Visualization":
    st.header("Gene List Visualization")
    gene_lists = st.text_area("Enter comma-separated gene lists (one or two lists, separate lists with a newline):")
    if gene_lists:
        lists = [list.strip().split(',') for list in gene_lists.split('\n') if list.strip()]
        if len(lists) > 2:
            st.error("Please enter at most two lists of genes.")
        else:
            all_genes = [gene.strip() for sublist in lists for gene in sublist]
            valid_genes = [gene for gene in all_genes if gene in gene_embeddings]
            
            if len(valid_genes) > 1:
                embeddings_to_plot = np.array([gene_embeddings[gene] for gene in valid_genes])
                
                method = st.radio("Select visualization method:", ["t-SNE", "UMAP"])
                
                if method == "t-SNE":
                    tsne = TSNE(n_components=2, random_state=42)
                    embeddings_2d = tsne.fit_transform(embeddings_to_plot)
                else:
                    reducer = umap.UMAP(random_state=42)
                    embeddings_2d = reducer.fit_transform(embeddings_to_plot)
                
                df = pd.DataFrame({
                    'x': embeddings_2d[:, 0],
                    'y': embeddings_2d[:, 1],
                    'gene': valid_genes,
                    'list': ['List 1' if gene in lists[0] else 'List 2' for gene in valid_genes]
                })
                
                fig = px.scatter(df, x='x', y='y', text='gene', color='list',
                                 title=f"{method} Visualization of Gene Lists")
                fig.update_traces(textposition='top center')
                st.plotly_chart(fig)
            else:
                st.error("Not enough valid genes to visualize.")

elif function == "Calculator":
    st.header("Gene Relationship Calculator")
    st.write("Calculate: gene_D is to gene_C as gene_A is to gene_B")
    gene_a = st.text_input("Enter gene A:")
    gene_b = st.text_input("Enter gene B:")
    gene_c = st.text_input("Enter gene C:")
    
    if gene_a and gene_b and gene_c:
        if all(gene in gene_embeddings for gene in [gene_a, gene_b, gene_c]):
            result_gene, similarity = gene_calculation(gene_a, gene_b, gene_c, gene_embeddings)
            st.write(f"The gene D that completes the relationship is: {result_gene}")
            st.write(f"Similarity score: {similarity:.4f}")
        else:
            st.error("One or more genes not found in the embeddings.")

elif function == "Computing Similarity":
    st.header("Compute Gene Similarity")
    genes = st.text_input("Enter two space-separated genes:")
    if genes:
        gene1, gene2 = genes.split()
        if gene1 in gene_embeddings and gene2 in gene_embeddings:
            similarity = calculate_similarity(gene1, gene2, gene_embeddings)
            st.write(f"Similarity between {gene1} and {gene2}: {similarity:.4f}")
        else:
            st.error("One or both genes not found in the embeddings.")

st.sidebar.info("This app analyzes gene embeddings using various techniques.")