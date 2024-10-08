import streamlit as st
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import plotly.express as px
import plotly.graph_objects as go
from utils.data_utils import get_gene_id

# Global colorblind-friendly color palette
COLORBLIND_COLORS = {
    'blue': '#3182bd',
    'orange': '#e6550d',
    'green': '#31a354',
    'red': '#de2d26',
    'purple': '#756bb1',
    'brown': '#8c6d31',
    'pink': '#fd8d3c',
    'gray': '#969696'
}

def plot_gene_embeddings(embeddings, genes, method, target_gene=None, lists=None, group_labels=None):
    """
    Plots gene embeddings using PCA, t-SNE, or UMAP.

    Args:
        embeddings (list): List of embedding vectors.
        genes (list): List of gene identifiers.
        method (str): Dimensionality reduction method ("PCA", "t-SNE", "UMAP").
        target_gene (str, optional): Gene identifier to highlight.
        lists (list of lists, optional): Lists of genes for coloring.
        group_labels (dict, optional): Custom labels for groups when target_gene is provided.
    """
    if len(embeddings) < 3:
        st.error("Not enough genes to visualize. Please provide at least 3 genes.")
        return

    if method == "PCA":
        reducer = PCA(n_components=2, random_state=42)
    elif method == "t-SNE":
        n_samples = len(embeddings)
        perplexity = min(30, max(5, n_samples // 5))
        reducer = TSNE(n_components=2, random_state=42, perplexity=perplexity, n_iter=1000, learning_rate='auto')
    else:
        n_neighbors = min(15, len(embeddings) - 1)
        reducer = umap.UMAP(n_neighbors=n_neighbors, random_state=42, n_components=2)

    embeddings_2d = reducer.fit_transform(np.array(embeddings))
    
    df = pd.DataFrame({
        'x': embeddings_2d[:, 0],
        'y': embeddings_2d[:, 1],
        'gene': genes
    })
    
    if target_gene:
        # Create a boolean column for the target gene
        df['is_input'] = df['gene'] == target_gene
        
        if group_labels:
            # Map boolean to custom labels
            df['group'] = df['is_input'].map(group_labels)
            # Define color mapping based on labels
            fig = px.scatter(
                df, 
                x='x', 
                y='y', 
                text='gene', 
                color='group',
                color_discrete_map={
                    group_labels[True]: COLORBLIND_COLORS['red'], 
                    group_labels[False]: COLORBLIND_COLORS['blue']
                },
                title=f"{method} Visualization of Similar Genes"
            )
        else:
            # Default labeling
            fig = px.scatter(
                df, 
                x='x', 
                y='y', 
                text='gene', 
                color='is_input',
                color_discrete_map={True: COLORBLIND_COLORS['red'], False: COLORBLIND_COLORS['blue']},
                title=f"{method} Visualization of Similar Genes"
            )

    elif lists:
        color_map = {'List 1': COLORBLIND_COLORS['blue'], 
                     'List 2': COLORBLIND_COLORS['orange'], 
                     'List 3': COLORBLIND_COLORS['green']}
        df['list'] = ['List 1' if gene in lists[0] else 'List 2' if len(lists) > 1 and gene in lists[1] else 'List 3' for gene in genes]
        fig = px.scatter(df, x='x', y='y', text='gene', color='list',
                        color_discrete_map=color_map,
                        title=f"{method} Visualization of Gene Lists")
    else:
        fig = px.scatter(df, x='x', y='y', text='gene',
                         title=f"{method} Visualization of Gene Embeddings")
    
    fig.update_traces(textposition='top center')
    fig.update_layout(
        height=600,
        plot_bgcolor='rgba(240, 240, 240, 0.8)',  # Light gray background
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        font=dict(family="Arial, sans-serif"),  # Elegant font
        shapes=[
            dict(
                type="rect",
                xref="paper", yref="paper",
                x0=0, y0=0, x1=1, y1=1,
                line=dict(color="rgba(0,0,0,0)", width=0),
                fillcolor="rgba(255, 255, 255, 0)"
            )
        ]
    )
    st.plotly_chart(fig, use_container_width=True)

def plot_gene_relationship(gene_a, gene_b, gene_c, gene_d, gene_embeddings):
    # Get embeddings for the four genes
    embeddings = [gene_embeddings[gene] for gene in [gene_a, gene_b, gene_c, gene_d]]
    
    # Use PCA to reduce to 2D for visualization
    pca = PCA(n_components=2)
    embeddings_2d = pca.fit_transform(embeddings)
    
    # Create a DataFrame for the plot
    df = pd.DataFrame({
        'x': embeddings_2d[:, 0],
        'y': embeddings_2d[:, 1],
        'gene': [gene_a, gene_b, gene_c, gene_d]
    })
    
    # Create the plot
    fig = go.Figure()
    
    # Add points
    fig.add_trace(go.Scatter(
        x=df['x'], y=df['y'], text=df['gene'],
        mode='markers+text', textposition="top center",
        marker=dict(size=10, color=[COLORBLIND_COLORS['red'], COLORBLIND_COLORS['blue'], 
                                    COLORBLIND_COLORS['green'], COLORBLIND_COLORS['purple']])
    ))
    
    # Add arrows
    fig.add_annotation(
        x=df.loc[df['gene'] == gene_b, 'x'].iloc[0],
        y=df.loc[df['gene'] == gene_b, 'y'].iloc[0],
        ax=df.loc[df['gene'] == gene_a, 'x'].iloc[0],
        ay=df.loc[df['gene'] == gene_a, 'y'].iloc[0],
        xref="x", yref="y", axref="x", ayref="y",
        showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2,
        arrowcolor=COLORBLIND_COLORS['red']
    )
    fig.add_annotation(
        x=df.loc[df['gene'] == gene_d, 'x'].iloc[0],
        y=df.loc[df['gene'] == gene_d, 'y'].iloc[0],
        ax=df.loc[df['gene'] == gene_c, 'x'].iloc[0],
        ay=df.loc[df['gene'] == gene_c, 'y'].iloc[0],
        xref="x", yref="y", axref="x", ayref="y",
        showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2,
        arrowcolor=COLORBLIND_COLORS['green']
    )
    
    fig.update_layout(
        title="Gene Relationship Visualization",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        showlegend=False,
        height=500,
        width=700
    )
    
    return fig