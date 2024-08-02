import streamlit as st
from utils.data_utils import read_gene_embeddings, load_gene_id_mapping, get_gene_id, add_google_analytics
from utils.gene_utils import find_closest_genes, calculate_similarity, gene_calculation
from utils.visualization_utils import plot_gene_embeddings, plot_gene_relationship
import pandas as pd
from utils.data_utils import get_gene_id
from PIL import Image
# Load gene embeddings and ID mapping
file_path = "data/GeneRAIN-vec.200d.txt.gz"
mapping_file_path = "data/gencode.v43.Ensembl_ID_gene_symbol_mapping.GeneRAIN.txt"
gene_embeddings = read_gene_embeddings(file_path)
ensembl_to_symbol, symbol_to_ensembl = load_gene_id_mapping(mapping_file_path)

# Sidebar for navigation
st.sidebar.title("GeneRAIN-vec")

# Create a container for the navigation buttons
nav_container = st.sidebar.container()

# Create buttons that look like tabs
if nav_container.button("Home"):
    st.session_state.page = "Home"
if nav_container.button("Similar Genes"):
    st.session_state.page = "Similar Genes"
if nav_container.button("Visualization"):
    st.session_state.page = "Visualization"
if nav_container.button("Calculator"):
    st.session_state.page = "Calculator"
if nav_container.button("Computing Similarity"):
    st.session_state.page = "Computing Similarity"
st.markdown("""
<style>
    /* Existing styles for navigation buttons */
    .stButton>button {
        width: 100%;
        border-radius: 0;
        border: none;
        border-bottom: 1px solid #e6e6e6;
        text-align: left;
        padding: 0.5rem 1rem;
    }
    .stButton>button:hover {
        background-color: #f0f2f6;
    }
    .stButton>button:focus {
        box-shadow: none;
        background-color: #e6e6e6;
    }
    /* New styles for primary buttons */
    .stButton.primary-button>button {
        background-color: white !important;
        color: #4CAF50 !important;
        border: 2px solid #4CAF50 !important;
        border-radius: 20px !important;
        text-align: center !important;
        transition: all 0.3s ease !important;
    }
    .stButton.primary-button>button:hover {
        background-color: #4CAF50 !important;
        color: white !important;
    }
</style>
""", unsafe_allow_html=True)

# Home page
def home():
    try:
        measurement_id = st.secrets['Measurement_Id']
        add_google_analytics()
        # st.write(f"Google Analytics has been integrated successfully, ID {measurement_id}.")
    except FileNotFoundError:
        st.write("Secrets file not found. Google Analytics integration skipped.")
    except KeyError:
        st.write("Google Analytics Measurement ID is not available in the secrets.")

    image = Image.open("data/crop2.jpg")
    st.image(image, use_column_width=True)
    # st.title("GeneRAIN-vec Gene Embedding Analysis Tools")
    st.markdown("""
    ## GeneRAIN-vec Gene Embedding Analysis Tools

    This web application allows you to explore and analyze gene embeddings derived from the GeneRAIN model, a state-of-the-art deep learning approach for understanding gene relationships and functions.

    ### What are Gene Embeddings?
    Gene embeddings are vector representations of genes in a high-dimensional space. These embeddings capture complex relationships between genes based on their expression patterns. In our case, each gene is represented by a 200-dimensional vector.

    ### About the GeneRAIN Models
    GeneRAIN are transformer-based models trained on a large dataset of 410K human bulk RNA-seq samples. These embeddings are derived from the GPT protein-coding+lncRNA model, which uses a novel 'Binning-By-Gene' normalization method and a GPT (Generative Pre-trained Transformer) architecture to learn multifaceted representations of genes.

    ### How to Use the Tools
    Use the sidebar to select from four main functions:
    1. Find Similar Genes
    2. Visualize Gene Lists
    3. Gene Relationship Calculator
    4. Compute Gene Similarity

    Each function provides unique insights into gene relationships and properties.

    For more information, please refer to our paper: 
                
    Zheng Su, Mingyan Fang, Andrei Smolnikov, Marcel E. Dinger, Emily Oates, Fatemeh Vafaee. Multifaceted Representation of Genes via Deep Learning of Gene Expression Networks, bioRxiv 2024.03.07.583777, [doi: https://doi.org/10.1101/2024.03.07.583777](https://www.biorxiv.org/content/10.1101/2024.03.07.583777)

    The complete gene embedding matrix is avaiable in [our Zenodo repo](https://zenodo.org/records/10408775).
                
    Contact: zheng.su1@unsw.edu.au
    """)

def similar_genes():
    st.header("Find Similar Genes")
    st.markdown("""
    This function allows you to find genes that are similar to a given gene based on their embeddings.
    The similarity is calculated using cosine similarity in the 200-dimensional embedding space.
    
    ### How to use:
    1. Enter a gene symbol or Ensembl ID in the text box below.
    2. Select the number of similar genes to display.
    3. The app will display the most similar genes and their similarity scores.
    4. A visualization of these genes will be shown using PCA, t-SNE, or UMAP for dimensionality reduction.
    """)
    
    col1, col2 = st.columns([3, 1])
    with col1:
        target_gene = st.text_input("Enter a gene name or Ensembl ID (e.g. TP53 or ENSG00000141510):", key="target_gene_input")
        num_genes = st.selectbox("Number of similar genes to return:", options=[10, 20, 50, 100], index=0)
        submit_button = st.button("Submit", key="submit_button", use_container_width=True, help="Click to find similar genes")
        st.markdown('<div class="primary-button"></div>', unsafe_allow_html=True)       


    if submit_button or target_gene:
        gene_id = get_gene_id(target_gene, gene_embeddings, ensembl_to_symbol, symbol_to_ensembl)
        if gene_id:
            similar_genes = find_closest_genes(gene_id, gene_embeddings, n=num_genes)
            
            # Prepare data for the table
            table_data = []
            for gene, similarity in similar_genes:
                symbol = ensembl_to_symbol.get(gene, gene)
                ensembl = symbol_to_ensembl.get(gene, gene)
                if symbol == ensembl:
                    display_name = symbol
                else:
                    display_name = f"{symbol} ({ensembl})"
                table_data.append({"Gene Symbol": symbol, "Ensembl ID": ensembl, "Similarity Score": f"{similarity:.4f}"})
            
            # Display the table with pagination
            st.write(f"Genes similar to {target_gene}:")
            df = pd.DataFrame(table_data)
            page_size = 10
            num_pages = (len(df) + page_size - 1) // page_size
            
            # Page selector
            col1, col2 = st.columns([1, 3])
            with col1:
                page_number = st.selectbox("Page", options=range(1, num_pages + 1))
            
            # Display table based on page selection
            start_idx = (page_number - 1) * page_size
            end_idx = min(start_idx + page_size, len(df))
            st.table(df.iloc[start_idx:end_idx])
            
            # Download button
            csv = df.to_csv(index=False)
            st.download_button(
                label="Download CSV",
                data=csv,
                file_name=f"similar_genes_to_{target_gene}.csv",
                mime="text/csv",
            )
            
            # Visualization
            st.subheader("Visualization of Similar Genes")
            # genes_to_plot = [gene_id] + [gene for gene, _ in similar_genes[:min(19, num_genes-1)]]
            genes_to_plot = [gene_id] + [gene for gene, _ in similar_genes]
            embeddings_to_plot = [gene_embeddings[gene] for gene in genes_to_plot]
            
            method = st.radio("Select visualization method:", ["PCA", "t-SNE", "UMAP"])
            plot_gene_embeddings(embeddings_to_plot, genes_to_plot, method, gene_id)
        else:
            st.error(f"Gene {target_gene} not found in the embeddings.")

def visualization():
    st.header("Gene List Visualization")
    st.markdown("""
    This function allows you to visualize the relationships between multiple genes in a 2D space.
    
    ### How to use:
    1. Enter up to three comma-separated lists of genes (symbols or Ensembl IDs) in the text area below.
    2. If entering multiple lists, separate them with a newline.
    3. The app will create a 2D projection of these genes using PCA, t-SNE, or UMAP.
    4. Genes from different lists will be colored differently.

    """)

    def update_input():
        st.session_state.gene_lists = "TP53,MDM2,CDKN1A,BAX\nMYC,CCND1,CDK4,E2F1\nBRCA1,BRCA2,ATM,CHEK2"

    # Initialize gene_lists in session state if it doesn't exist
    if 'gene_lists' not in st.session_state:
        st.session_state.gene_lists = ""

    gene_lists = st.text_area("Enter comma-separated gene lists (up to three lists, separate lists with a newline):", st.session_state.gene_lists)
    
    # Add the example link
    st.button("Input example genes", key="viz_example_button", on_click=update_input, help="Click to input example genes")

    method = st.radio("Select visualization method:", ["PCA", "t-SNE", "UMAP"])
    submit_button = st.button("Submit", key="viz_submit_button", use_container_width=True, help="Click to visualize genes")
    st.markdown('<div class="primary-button"></div>', unsafe_allow_html=True)

    if submit_button and gene_lists:
        # Store the input in session state
        st.session_state.gene_lists = gene_lists


        # Parse input
        if ',' in gene_lists:
            lists = [list(map(str.strip, gene_group.split(','))) for gene_group in gene_lists.split('\n') if gene_group.strip()]
        else:
            lists = [gene_lists.split()]

        if len(lists) > 3:
            st.error("Please enter at most three lists of genes.")
        else:
            all_genes = [gene.strip() for sublist in lists for gene in sublist]
            valid_genes = []
            invalid_genes = []
            for gene in all_genes:
                gene_id = get_gene_id(gene, gene_embeddings, ensembl_to_symbol, symbol_to_ensembl)
                if gene_id:
                    valid_genes.append(gene_id)
                else:
                    invalid_genes.append(gene)
            
            if invalid_genes:
                st.warning(f"The following genes were not found and will be excluded: {', '.join(invalid_genes)}")
            
            if len(valid_genes) > 1:
                embeddings_to_plot = [gene_embeddings[gene] for gene in valid_genes]
                plot_gene_embeddings(embeddings_to_plot, valid_genes, method, lists=lists)
            else:
                st.error("Not enough valid genes to visualize.")

# Update the calculator function
def calculator():
    st.header("Gene Relationship Calculator")
    st.markdown("""
    This function allows you to explore gene relationships using vector arithmetic on gene embeddings.
    It's similar to word analogies in natural language processing.
    
    ### How to use:
    1. Enter three gene symbols or Ensembl IDs for A, B, and C.
    2. The app will find gene D such that the relationship A:B is similar to C:D.
    3. The result shows the most likely gene D and its similarity score.
    4. A visualization of the four genes and their relationships will be displayed.
    
    """)
    
    def update_input():
        st.session_state.gene_a = "BRCA1"
        st.session_state.gene_b = "BRCA2"
        st.session_state.gene_c = "TP53"

    st.write("Calculate: gene_D is to gene_C as gene_A is to gene_B")
    gene_a = st.text_input("Enter gene A:", key="gene_a")
    gene_b = st.text_input("Enter gene B:", key="gene_b")
    gene_c = st.text_input("Enter gene C:", key="gene_c")
    
    # Add the example link
    st.button("Input example genes", key="calc_example_button", on_click=update_input, help="Click to input example genes")

    if st.button("Calculate", key="calc_button", use_container_width=True, help="Click to calculate gene relationship"):
        st.markdown('<div class="primary-button"></div>', unsafe_allow_html=True)
    
        if gene_a and gene_b and gene_c:
            gene_a_id = get_gene_id(gene_a, gene_embeddings, ensembl_to_symbol, symbol_to_ensembl)
            gene_b_id = get_gene_id(gene_b, gene_embeddings, ensembl_to_symbol, symbol_to_ensembl)
            gene_c_id = get_gene_id(gene_c, gene_embeddings, ensembl_to_symbol, symbol_to_ensembl)
            
            if all([gene_a_id, gene_b_id, gene_c_id]):
                result_gene, similarity = gene_calculation(gene_a_id, gene_b_id, gene_c_id, gene_embeddings)
                display_name = ensembl_to_symbol.get(result_gene, result_gene)
                st.write(f"The gene D that completes the relationship is: {display_name} ({result_gene})")
                st.write(f"Similarity score: {similarity:.4f}")
                
                # Visualize the gene relationship
                fig = plot_gene_relationship(gene_a_id, gene_b_id, gene_c_id, result_gene, gene_embeddings)
                st.plotly_chart(fig)
                
            else:
                invalid_genes = [gene for gene, gene_id in zip([gene_a, gene_b, gene_c], [gene_a_id, gene_b_id, gene_c_id]) if not gene_id]
                st.error(f"The following genes were not found: {', '.join(invalid_genes)}")
        else:
            st.error("Please enter all three genes.")

def computing_similarity():
    st.header("Compute Gene Similarity")
    st.markdown("""
    This function calculates the similarity between two genes based on their embeddings.
    
    ### How to use:
    1. Enter two gene symbols or Ensembl IDs, separated by a space.
    2. Click the 'Calculate' button.
    3. The app will display the similarity between the genes and how this similarity ranks among all other genes.
    """)

    # Function to update the input field
    def update_input():
        st.session_state.gene_input = "TP53 MDM2"

    # Display the input field
    genes = st.text_input("Enter two space-separated genes:", key="gene_input")
    
    # Add the example link
    st.button("Input example genes", key="example_button", on_click=update_input, help="Click to input example genes")


    if st.button("Calculate", key="sim_calc_button", use_container_width=True, help="Click to calculate gene embedding similarity"):
        st.markdown('<div class="primary-button"></div>', unsafe_allow_html=True)
        if genes:
            gene1, gene2 = genes.split()
            gene1_id = get_gene_id(gene1, gene_embeddings, ensembl_to_symbol, symbol_to_ensembl)
            gene2_id = get_gene_id(gene2, gene_embeddings, ensembl_to_symbol, symbol_to_ensembl)
            
            if gene1_id and gene2_id:
                similarity, rank_1, _, rank_2, _ = calculate_similarity(gene1_id, gene2_id, gene_embeddings)
                total_genes = len(gene_embeddings)
                
                st.write(f"Similarity between {gene1}  and {gene2} : {similarity:.4f}")
                st.write(f"\nCompared to {gene2}, {rank_1}/{total_genes} other genes ({rank_1/total_genes*100:.2f}%) have higher similarity with {gene1}.")
                st.write(f"Compared to {gene1}, {rank_2}/{total_genes} other genes ({rank_2/total_genes*100:.2f}%) have higher similarity with {gene2}.")
                
                st.write("\nNote: The rankings may differ for each gene due to the nature of the embedding space.")
            else:
                invalid_genes = [gene for gene, gene_id in zip([gene1, gene2], [gene1_id, gene2_id]) if not gene_id]
                st.error(f"The following genes were not found: {', '.join(invalid_genes)}")
        else:
            st.error("Please enter two genes.")



# Initialize the page state if it doesn't exist
if 'page' not in st.session_state:
    st.session_state.page = "Home"

# Rest of your code remains the same, but replace the if-elif statements with:
if st.session_state.page == "Home":
    home()
elif st.session_state.page == "Similar Genes":
    similar_genes()
elif st.session_state.page == "Visualization":
    visualization()
elif st.session_state.page == "Calculator":
    calculator()
elif st.session_state.page == "Computing Similarity":
    computing_similarity()

# Footer
st.sidebar.markdown("---")
st.sidebar.info("""
This app analyzes gene embeddings using various techniques. 
The embeddings are derived from the GeneRAIN model, which was trained on a large dataset of human bulk RNA-seq samples.
For more details, please refer to our paper.
""")

# Add license statement
st.sidebar.markdown("---")
st.sidebar.markdown("""
<small>
This work is licensed under a 
<a href="https://creativecommons.org/licenses/by-nc/4.0/" target="_blank">
CC BY-NC 4.0 License</a>.
</small>
""", unsafe_allow_html=True)