import numpy as np
from scipy.spatial.distance import cosine

def find_closest_genes(target_gene, gene_embeddings, n=20):
    target_embedding = gene_embeddings[target_gene]
    similarities = [(gene, 1 - cosine(embedding, target_embedding)) 
                    for gene, embedding in gene_embeddings.items() 
                    if gene != target_gene]
    return sorted(similarities, key=lambda x: x[1], reverse=True)[:n]

def calculate_similarity(gene1, gene2, gene_embeddings):
    similarity = 1 - cosine(gene_embeddings[gene1], gene_embeddings[gene2])
    
    # Calculate rank and quantile for gene1
    similarities_1 = [1 - cosine(gene_embeddings[gene1], emb) for emb in gene_embeddings.values()]
    rank_1 = sum(s > similarity for s in similarities_1)
    quantile_1 = np.percentile([s for s in similarities_1 if s < 1], (1 - rank_1 / len(similarities_1)) * 100)
    
    # Calculate rank and quantile for gene2
    similarities_2 = [1 - cosine(gene_embeddings[gene2], emb) for emb in gene_embeddings.values()]
    rank_2 = sum(s > similarity for s in similarities_2)
    quantile_2 = np.percentile([s for s in similarities_2 if s < 1], (1 - rank_2 / len(similarities_2)) * 100)
    
    return similarity, rank_1, quantile_1, rank_2, quantile_2

def gene_calculation(gene_a, gene_b, gene_c, gene_embeddings):
    result_vector = (gene_embeddings[gene_b] - gene_embeddings[gene_a] + gene_embeddings[gene_c])
    similarities = [(gene, 1 - cosine(embedding, result_vector)) 
                    for gene, embedding in gene_embeddings.items()
                    if gene not in [gene_a, gene_b, gene_c]]  # Exclude input genes
    return sorted(similarities, key=lambda x: x[1], reverse=True)[0]
