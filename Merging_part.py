import argparse
import requests
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import plotly.express as px

def fetch_and_process_data(api_url, file_path, region_name, genes, cell_type):
    print("Fetching Data...")
    # Fetch the data from the API link and save it to a local file
    response = requests.get(api_url)
    response.raise_for_status()  # Ensure the request was successful
    with open(file_path, "wb") as f:
        f.write(response.content)

    adata = sc.read_h5ad(file_path)

    # Filter the specified cell type
    cell_type_data = adata[adata.obs['cell_type'] == cell_type]
    print("Filtering data for specified cell type...")

    # Select only the relevant NMDA genes
    valid_genes = [gene for gene in genes if gene in adata.var_names]

    if not valid_genes:
        print("No NMDA receptor genes found in the dataset.")

    # Subset the data to only NMDA receptor genes
    nmda_data = cell_type_data[:, valid_genes]

    # Add region information
    nmda_data.obs['region'] = region_name

    return nmda_data

def plot_gene_expression(combined_data, genes, plot_type):
    df = combined_data.to_df()
    df['region'] = combined_data.obs['region']

    if plot_type == "dot":
        sc.pl.dotplot(combined_data, var_names=genes, groupby='region', standard_scale='var', show=False)
        print("Creating dot plot...")

        # Get the current figure
        fig = plt.gcf()

        # Rotate x-axis labels
        for ax in fig.axes:
            plt.sca(ax)
            plt.xticks(rotation=90)

        # Add label to color bar
        for ax in fig.axes:
            if ax.get_ylabel() == 'Mean expression\nin group':
                cbar = ax.collections[0].colorbar
                cbar.set_label('Mean expression in group')

        # Adjust spacing between subplots to make sure labels are not cut off
        plt.subplots_adjust(top=0.85, bottom=0.25, left=0.2, right=0.95, hspace=0.3, wspace=0.4)

        # Set title at the top
        plt.suptitle('Dot Plot of Selected Genes', fontsize=24, y=0.98)

        # Show the plot
        plt.show()

        print("Dot plot created.")

    elif plot_type == "cluster_profile":
        print("Creating cluster profiles...")
        cluster_means = combined_data[:, genes].to_df().groupby(combined_data.obs["region"]).mean()
        fig = px.line(cluster_means.T, labels={"value": "Expression Level", "index": "Genes", "variable": "region"},
                      title="Cluster Profiles")
        fig.show()
        print("Cluster profiles created.")

def find_top_genes(combined_data, genes):
    regions = combined_data.obs['region'].unique()
    top_genes = {}

    for region in regions:
        region_data = combined_data[combined_data.obs['region'] == region]
        mean_expression = np.mean(region_data[:, genes].X.toarray(), axis=0)
        top_gene = genes[np.argmax(mean_expression)]
        top_genes[region] = (top_gene, mean_expression[np.argmax(mean_expression)])

    return top_genes

def find_overall_top_gene(combined_data, genes):
    mean_expression = np.mean(combined_data[:, genes].X.toarray(), axis=0)
    top_gene = genes[np.argmax(mean_expression)]
    return top_gene, mean_expression[np.argmax(mean_expression)]

def plot_top_genes(top_genes, overall_top_gene=None):
    regions = list(top_genes.keys())
    genes = [top_genes[region][0] for region in regions]
    expressions = [top_genes[region][1] for region in regions]

    plt.figure(figsize=(10, 6))
    plt.barh(regions, expressions, color='skyblue')
    for index, value in enumerate(expressions):
        plt.text(value, index, f'{genes[index]}: {value:.2f}')
    plt.xlabel('Expression Level')
    plt.ylabel('Hippocampus Region')
    plt.title('Top Expressed Genes in Each Hippocampus Region')

    if overall_top_gene:
        top_gene, expression = overall_top_gene
        plt.suptitle(f'Top Expressed Gene Overall: {top_gene} with Expression Level {expression:.2f}', y=0.95,
                     fontsize=14)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

def plot_combined_top_gene(combined_data, genes):
    mean_expression = np.mean(combined_data[:, genes].X.toarray(), axis=0)
    top_gene = genes[np.argmax(mean_expression)]
    expression_level = mean_expression[np.argmax(mean_expression)]

    plt.figure(figsize=(6, 6))
    plt.bar([top_gene], [expression_level], color='skyblue')
    plt.xlabel('Gene')
    plt.ylabel('Expression Level')
    plt.title(f'Top Expressed Gene Overall: {top_gene}')
    plt.show()

def coexpression_analysis(combined_data, genes):
    print("Performing co-expression analysis...")

    # Extract the unique regions
    regions = combined_data.obs['region'].unique()

    for region in regions:
        region_data = combined_data[combined_data.obs['region'] == region]

        # Extract the expression data for the selected genes
        expression_data = region_data[:, genes].X.toarray()

        # Calculate the correlation matrix
        correlation_matrix = np.corrcoef(expression_data, rowvar=False)

        # Create a DataFrame for the correlation matrix
        corr_df = pd.DataFrame(correlation_matrix, index=genes, columns=genes)

        # Plot the correlation matrix as a heatmap using seaborn
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_df, annot=True, cmap="coolwarm", xticklabels=True, yticklabels=True)
        plt.title(f"Gene Co-expression Heatmap for {region}")
        plt.show()

    # Combined co-expression analysis for all regions
    combined_expression_data = combined_data[:, genes].X.toarray()
    combined_correlation_matrix = np.corrcoef(combined_expression_data, rowvar=False)
    combined_corr_df = pd.DataFrame(combined_correlation_matrix, index=genes, columns=genes)

    plt.figure(figsize=(10, 8))
    sns.heatmap(combined_corr_df, annot=True, cmap="coolwarm", xticklabels=True, yticklabels=True)
    plt.title("Combined Gene Co-expression Heatmap for All Regions")
    plt.show()

    print("Co-expression analysis completed.")

def main(api_files_regions, genes, plot_type, cell_type):
    combined_data_list = []
    region_names = []

    for api_url, file_path, region_name in api_files_regions:
        region_data = fetch_and_process_data(api_url, file_path, region_name, genes, cell_type)
        combined_data_list.append(region_data)
        region_names.append(region_name)

    # Combine all the data
    combined_data = sc.concat(combined_data_list, join='outer', label='region', keys=region_names)

    # Manually add var information
    combined_data.var = combined_data_list[0].var.copy()

    # Plotting
    if plot_type in ['dot', 'cluster_profile']:
        plot_gene_expression(combined_data, genes, plot_type)
    elif plot_type == "coexpression":
        coexpression_analysis(combined_data, genes)
    elif plot_type == "top_genes":
        top_genes = find_top_genes(combined_data, genes)
        overall_top_gene = find_overall_top_gene(combined_data, genes)
        plot_top_genes(top_genes, overall_top_gene)
        plot_combined_top_gene(combined_data, genes)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and plot NMDA gene expression data.")
    parser.add_argument("--api_files_regions", nargs='+', type=str,
                        help="List of API URLs, file paths, and region names.")
    parser.add_argument("--genes", nargs='+', type=str, help="List of NMDA genes.")
    parser.add_argument("--plot_type", type=str,
                        choices=['dot', 'cluster_profile', 'coexpression', 'top_genes'],
                        help="Type of plot to generate.")
    parser.add_argument("--cell_type", type=str, help="Type of cell to filter the data.")

    args = parser.parse_args()

    # Ensure the necessary arguments are provided
    if not args.api_files_regions or not args.genes or not args.plot_type or not args.cell_type:
        parser.print_help()
        exit(1)

    # Convert input strings to tuples for api_files_regions
    api_files_regions = [tuple(item.split(',')) for item in args.api_files_regions]

    main(api_files_regions, args.genes, args.plot_type, args.cell_type)
