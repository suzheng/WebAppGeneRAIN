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


def plot_calculator_plus(terms, result_vector, gene_embeddings, ensembl_to_symbol, result_gene=None):
    COLORBLIND_COLORS = [
        '#377eb8',  # Blue
        '#ff7f00',  # Orange
        '#4daf4a',  # Green
        '#f781bf',  # Pink
        '#a65628',  # Brown
        '#984ea3',  # Purple
        '#999999',  # Gray
        '#e41a1c',  # Red
        '#dede00'   # Yellow
    ]
    """
    Visualizes the linear combination of gene embeddings as connected vectors in a 2D space.
    
    Each gene is represented as an arrow starting from the end of the previous vector.
    The first vector starts at the origin (0,0). The resultant vector is plotted from origin to the end of the last vector.
    If a result_gene is provided, its embedding is visualized as an additional vector starting from the end of the last term vector.
    
    Args:
        terms (list of tuples): Each tuple contains (coefficient, gene_id).
        result_vector (numpy array): The resulting vector from the linear combination.
        gene_embeddings (dict): Dictionary of gene embeddings.
        result_gene (str, optional): The Ensembl ID of the resulting gene.
    """
    # Collect all vectors: term vectors
    term_vectors = []
    gene_names = []
    coefficients = []
    
    for coeff, gene in terms:
        embedding = gene_embeddings[gene]
        term_vectors.append(coeff * embedding)
        gene_names.append(gene)
        coefficients.append(coeff)
    
    # If result_gene is provided, get its embedding
    if result_gene and result_gene in gene_embeddings:
        closest_gene_embedding = gene_embeddings[result_gene]
        all_vectors = term_vectors + [result_vector, closest_gene_embedding]
        include_closest = True
    else:
        all_vectors = term_vectors + [result_vector]
        include_closest = False

    # Perform PCA on all vectors to reduce to 2D
    pca = PCA(n_components=2)
    vectors_2d = pca.fit_transform(all_vectors)
    
    # Separate PCA-transformed vectors
    term_vectors_2d = vectors_2d[:-2] if include_closest else vectors_2d[:-1]
    resultant_vector_2d = vectors_2d[-2] if include_closest else vectors_2d[-1]
    if include_closest:
        closest_gene_vector_2d = vectors_2d[-1]
    
    # Initialize starting point
    x_start, y_start = 0, 0
    
    # Initialize Plotly figure
    fig = go.Figure()
    
    # Plot each term as a vector
    for idx, (vec, gene, coeff) in enumerate(zip(term_vectors_2d, gene_names, coefficients)):
        x_end, y_end = vec
        color = COLORBLIND_COLORS[idx % len(COLORBLIND_COLORS)]
        
        # Format the label with coefficient and sign
        if coeff == 1:
            label = f"+ {ensembl_to_symbol.get(gene, gene)}"
        elif coeff == -1:
            label = f"- {ensembl_to_symbol.get(gene, gene)}"
        elif coeff > 0:
            label = f"{coeff} * {ensembl_to_symbol.get(gene, gene)}"
        else:
            label = f"{coeff} * {ensembl_to_symbol.get(gene, gene)}"
        
        # Add arrow for the vector
        fig.add_trace(go.Scatter(
            x=[x_start, x_end],
            y=[y_start, y_end],
            mode='lines',
            line=dict(color=color, width=3),
            showlegend=False,
            hoverinfo='none'
        ))
        
        # Add marker at the end of the vector
        fig.add_trace(go.Scatter(
            x=[x_end],
            y=[y_end],
            mode='markers',
            marker=dict(color=color, size=8),
            showlegend=False,
            hoverinfo='none'
        ))
        
        # Add text label at the midpoint of the vector
        fig.add_trace(go.Scatter(
            x=[(x_start + x_end) / 2],
            y=[(y_start + y_end) / 2],
            mode='text',
            text=[label],
            textposition='middle center',
            showlegend=False,
            hoverinfo='none',
            textfont=dict(color=color, size=12)
        ))
        
        # Update starting point for the next vector
        x_start, y_start = x_end, y_end
    
    # Plot the resultant vector from origin to the end of the last vector
    fig.add_trace(go.Scatter(
        x=[0, x_start],
        y=[0, y_start],
        mode='lines',
        line=dict(color='black', width=4, dash='dash'),
        name='Resultant Vector',
        hoverinfo='none'
    ))
    
    # Add marker for the resultant vector end
    fig.add_trace(go.Scatter(
        x=[x_start],
        y=[y_start],
        mode='markers',
        marker=dict(color='black', size=10),
        showlegend=False,
        hoverinfo='none'
    ))
    
    # Add text label for the resultant vector at its midpoint
    fig.add_trace(go.Scatter(
        x=[x_start / 2],
        y=[y_start / 2],
        mode='text',
        text=['Resultant'],
        textposition='middle center',
        showlegend=False,
        hoverinfo='none',
        textfont=dict(color='black', size=12)
    ))
    
    # If closest gene is to be included, plot its embedding
    if include_closest:
        # Draw an arrow from the end of the last term vector to the closest gene's embedding
        fig.add_trace(go.Scatter(
            x=[x_start, closest_gene_vector_2d[0]],
            y=[y_start, closest_gene_vector_2d[1]],
            mode='lines',
            line=dict(color='green', width=3, dash='dot'),
            name='Closest Gene',
            hoverinfo='none'
        ))
        
        # Add marker for the closest gene
        fig.add_trace(go.Scatter(
            x=[closest_gene_vector_2d[0]],
            y=[closest_gene_vector_2d[1]],
            mode='markers',
            marker=dict(color='green', size=10),
            showlegend=False,
            hoverinfo='none'
        ))
        
        # Add text label for the closest gene at its position
        fig.add_trace(go.Scatter(
            x=[closest_gene_vector_2d[0]],
            y=[closest_gene_vector_2d[1]],
            mode='text',
            text=[f"Closest: {ensembl_to_symbol.get(result_gene, result_gene)}"],
            textposition='top right',
            showlegend=False,
            hoverinfo='none',
            textfont=dict(color='green', size=12)
        ))
    
    # Update layout for better visualization
    fig.update_layout(
        title="Calculator Plus Visualization",
        xaxis_title="PCA Component 1",
        yaxis_title="PCA Component 2",
        showlegend=True,
        height=600,
        width=700,
        template='plotly_white',
        xaxis=dict(scaleanchor="y", scaleratio=1),
        yaxis=dict(scaleanchor="x", scaleratio=1),
        shapes=[
            dict(
                type="line",
                x0=0, y0=0, x1=0, y1=0,
                line=dict(color="rgba(0,0,0,0)", width=0)
            )
        ],
        annotations=[
            dict(
                x=0, y=0,
                xref="x", yref="y",
                text="Origin",
                showarrow=False,
                font=dict(color="black")
            )
        ]
    )
    
    # Determine plot ranges based on vectors
    all_x = vectors_2d[:, 0]
    all_y = vectors_2d[:, 1]
    buffer = 0.1 * max(np.max(np.abs(all_x)), np.max(np.abs(all_y)))
    fig.update_xaxes(range=[min(all_x) - buffer, max(all_x) + buffer])
    fig.update_yaxes(range=[min(all_y) - buffer, max(all_y) + buffer])
    
    # Display the plot in Streamlit
    st.plotly_chart(fig, use_container_width=True)