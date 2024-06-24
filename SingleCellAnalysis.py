import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd
import networkx as nx
import plotly.graph_objs as go
import scipy.sparse as sp
from sklearn.decomposition import PCA

class SingleCellAnalysis:
    def __init__(self, file_path, genes, groupby='cell_type', preprocess=True, n_neighbors=10, n_pcs=20,
                 subset_fraction=0.1, plots='all'):
        """
        Initialize the SingleCellAnalysis class with the provided parameters.

        Args:
            file_path (str): Path to the input h5ad file containing the single-cell RNA sequencing data.
            genes (str): Comma-separated list of genes to analyze.
            groupby (str): Category to group cells by. Default is 'cell_type'.
            preprocess (bool): Whether to preprocess the data. Default is True.
            n_neighbors (int): Number of neighbors for UMAP. Default is 10.
            n_pcs (int): Number of principal components for PCA. Default is 20.
            subset_fraction (float): Fraction of data to subset. Default is 0.1.
            plots (str): Comma-separated list of plots to generate. Default is 'all'.
        """
        self.file_path = file_path  # Veri dosyasının yolu
        self.genes = genes.split(",")  # Analiz edilecek genlerin listesi
        self.groupby = groupby  # Hücreleri gruplayacak kategori
        self.adata = None  # Anndata objesi, veriyi tutacak
        self.preprocess = preprocess  # Veriyi ön işlemden geçirip geçirmeme seçeneği
        self.n_neighbors = n_neighbors  # Komşu sayısı (neighbors) parametresi
        self.n_pcs = n_pcs  # PCA bileşen sayısı
        self.subset_fraction = subset_fraction  # Veri alt küme oranı
        self.plots = plots.split(",")  # List of plots to generate

    def load_data(self):
        """
        Load the single-cell RNA sequencing data from the specified h5ad file.
        """
        try:
            print("Loading data...")
            self.adata = sc.read_h5ad(self.file_path)  # Veriyi yükleme
            print("Data loaded.")
        except Exception as e:
            print(f"An error occurred while loading the data: {e}")

    def preprocess_data(self):
        """
        Preprocess the loaded data by filtering cells and genes, normalizing counts, and performing PCA.
        """
        try:
            print("Preprocessing data...")
            self.adata = self.adata[
                np.random.choice(self.adata.shape[0], int(self.adata.shape[0] * self.subset_fraction),
                                 replace=False)]  # Veriyi alt kümeleme

            sc.pp.filter_cells(self.adata, min_genes=200)  # Minimum gen sayısına göre hücreleri filtreleme
            sc.pp.filter_genes(self.adata, min_cells=3)  # Minimum hücre sayısına göre genleri filtreleme

            total_counts = np.array(self.adata.X.sum(axis=1)).flatten()  # Toplam gen ekspresyon sayısını hesaplama
            self.adata.obs['total_counts'] = total_counts

            mt_gene_mask = self.adata.var_names.str.startswith('MT-')  # Mitokondri genlerini belirleme
            mt_gene_counts = np.array(
                self.adata[:, mt_gene_mask].X.sum(axis=1)).flatten()  # Mitokondri gen sayısını hesaplama
            self.adata.obs['pct_counts_mt'] = (
                mt_gene_counts / total_counts) * 100  # Mitokondri gen yüzdesini hesaplama

            sc.pp.normalize_total(self.adata, target_sum=1e4)  # Toplam gen sayısını normalize etme
            sc.pp.log1p(self.adata)  # Logaritmik ölçekleme
            sc.pp.highly_variable_genes(self.adata, min_mean=0.0125, max_mean=3,
                                        min_disp=0.5)  # Yüksek değişkenlik gösteren genleri belirleme
            self.adata = self.adata[:,
                         self.adata.var.highly_variable]  # Sadece yüksek değişkenlik gösteren genleri tutma

            sc.pp.regress_out(self.adata, ['total_counts',
                                           'pct_counts_mt'])  # Toplam gen sayısı ve mitokondri yüzdesini regresyonla çıkarma
            sc.pp.scale(self.adata, max_value=10)  # Veriyi ölçekleme
            print("Data preprocessed.")
        except Exception as e:
            print(f"An error occurred during preprocessing: {e}")

    def check_genes(self):
        """
        Check the genes in the data to include only those specified by the user.
        """
        print("Checking genes...")
        # Assuming 'feature_name' column has the gene names
        if 'feature_name' in self.adata.var.columns:
            name_to_ensembl = {name: idx for idx, name in zip(self.adata.var.index, self.adata.var["feature_name"])}
            valid_genes = [gene for gene in self.genes if gene in name_to_ensembl]
            valid_gene_ids = [name_to_ensembl[gene] for gene in valid_genes]
        else:
            print("'feature_name' column not found in adata.var.")
            valid_genes = []
            valid_gene_ids = []

        if not valid_genes:
            print("None of the genes of interest are present in the dataset.")
        else:
            print(f"Valid genes found in dataset: {valid_genes}")
        return valid_genes, valid_gene_ids

    def plot_pca_and_scree(self):
        """
        Plot PCA results and a scree plot to visualize the explained variance.
        """
        pca = PCA(n_components=10)
        pca_result = pca.fit_transform(self.adata.X.toarray())
        self.adata.obs['PC1'] = pca_result[:, 0]
        self.adata.obs['PC2'] = pca_result[:, 1]

        cell_types = self.adata.obs[self.groupby].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types)))
        cell_colors = dict(zip(cell_types, colors))

        plt.figure(figsize=(10, 8))
        for cell_type, color in cell_colors.items():
            mask = (self.adata.obs[self.groupby] == cell_type)
            plt.scatter(self.adata.obs.loc[mask, 'PC1'],
                        self.adata.obs.loc[mask, 'PC2'],
                        label=cell_type,
                        color=color,
                        alpha=0.6,
                        edgecolors='w',
                        s=50)

        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)')
        plt.title('PCA Analysis')
        plt.legend(title=self.groupby, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Legend için daha fazla yer açma
        plt.show()

        plt.figure(figsize=(10, 8))
        plt.plot(range(1, len(pca.explained_variance_ratio_) + 1), pca.explained_variance_ratio_, marker='o',
                 linestyle='--')
        plt.xlabel('Principal Component')
        plt.ylabel('Variance Explained')
        plt.title('Scree Plot')
        plt.grid()
        plt.show()

    def plot_violin(self, valid_genes, valid_gene_ids):
        """
        Generate violin plots for the specified genes.

        Args:
        valid_genes (list): List of valid genes to plot.
        valid_gene_ids (list): List of valid gene IDs to plot.
        """
        if valid_genes:
            print("Creating violin plot...")
            gene_name_to_id = {gene: gene_id for gene, gene_id in zip(valid_genes, valid_gene_ids)}
            sc.pl.violin(self.adata, keys=list(gene_name_to_id.values()), groupby=self.groupby, jitter=0.2,
                         multi_panel=True, rotation=90, show=False)

            fig = plt.gcf()
            for ax, gene_name in zip(fig.axes, valid_genes):
                ax.set_title(gene_name)

            plt.suptitle('Violin Plot of Selected Genes', y=1.0, fontsize=16)
            plt.show()
            print("Violin plot created.")
        else:
            print("No valid genes to plot for violin plot.")

    def plot_umap(self, valid_genes, valid_gene_ids):
        """
        Generate UMAP plots for the specified genes.

        Args:
        valid_genes (list): List of valid genes to plot.
        valid_gene_ids (list): List of valid gene IDs to plot.
        """
        if valid_genes:
            print("Creating UMAP plots...")
            sc.tl.pca(self.adata, svd_solver='arpack', n_comps=self.n_pcs)
            sc.pp.neighbors(self.adata, n_neighbors=self.n_neighbors, n_pcs=self.n_pcs)
            sc.tl.umap(self.adata)
            for gene_id, gene_name in zip(valid_gene_ids, valid_genes):
                plt.figure(figsize=(9, 5))
                sc.pl.umap(self.adata, color=gene_id, show=False)
                plt.title(f'UMAP Plot of {gene_name}', fontsize=15)
                plt.show()
            print("UMAP plots created.")
        else:
            print("No valid genes to plot for UMAP plot.")

    def plot_scatter(self, valid_genes, valid_gene_ids):
        """
        Generate scatter plots for the specified genes.

        Args:
        valid_genes (list): List of valid genes to plot.
        valid_gene_ids (list): List of valid gene IDs to plot.
        """
        if len(valid_genes) >= 2:
            print("Creating scatter plots...")
            cell_types = self.adata.obs[self.groupby].unique()
            colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types)))
            for i in range(len(valid_gene_ids) - 1):
                for j in range(i + 1, len(valid_gene_ids)):
                    x = self.adata[:, valid_gene_ids[i]].X.toarray().flatten()
                    y = self.adata[:, valid_gene_ids[j]].X.toarray().flatten()
                    cell_labels = self.adata.obs[self.groupby]
                    plt.figure(figsize=(10, 8))
                    for cell_type, color in zip(cell_types, colors):
                        mask = (cell_labels == cell_type)
                        plt.scatter(x[mask], y[mask], label=cell_type, color=color, alpha=0.6, edgecolors='w', s=50)
                    plt.xlabel(valid_genes[i])
                    plt.ylabel(valid_genes[j])
                    plt.title(f'Scatter Plot of {valid_genes[i]} vs {valid_genes[j]}')
                    plt.legend(title=self.groupby, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
                    plt.subplots_adjust(top=0.9, bottom=0.1, right=0.7)
                    plt.show()
                    print(f'Scatter plot for {valid_genes[i]} vs {valid_genes[j]} created.')
        else:
            print("Need at least 2 valid genes for scatter plot.")

    def plot_dotplot(self, valid_genes, valid_gene_ids):
        """
        Generate dot plots for the specified genes.

        Args:
        valid_genes (list): List of valid genes to plot.
        valid_gene_ids (list): List of valid gene IDs to plot.
        """
        if valid_genes:
            print("Creating dot plot...")
            dotplot = sc.pl.dotplot(self.adata, var_names=valid_gene_ids, groupby=self.groupby, show=False)
            fig = plt.gcf()
            # Adjusting the gene names on the x-axis
            ax = fig.axes[0]
            ax.set_xticklabels(valid_genes, rotation=90)
            fig.suptitle('Dot Plot of Selected Genes', y=0.90, fontsize=24)
            plt.subplots_adjust(top=0.85, bottom=0.3, left=0.2, right=0.95, hspace=0.3, wspace=0.4)
            for ax in fig.axes:
                if ax.get_ylabel() == 'Fraction of cells\nin group (%)':
                    ax.yaxis.set_label_position("left")
                    ax.yaxis.tick_left()
                    for label in ax.get_yticklabels():
                        label.set_horizontalalignment('left')
                        label.set_x(-0.1)
            plt.show()
            print("Dot plot created.")
        else:
            print("No valid genes to plot for dot plot.")

    def plot_interactive_pca(self, valid_genes, valid_gene_ids):
        """
        Generate an interactive PCA plot for the specified genes.

        Args:
        valid_genes (list): List of valid genes to plot.
        valid_gene_ids (list): List of valid gene IDs to plot.
        """
        if valid_genes:
            print("Creating interactive PCA plot...")
            sc.tl.pca(self.adata, n_comps=3)
            pca_df = pd.DataFrame(self.adata.obsm['X_pca'][:, :3], columns=['PC1', 'PC2', 'PC3'])
            pca_df[self.groupby] = self.adata.obs[self.groupby].values
            fig = px.scatter_3d(pca_df, x='PC1', y='PC2', z='PC3', color=self.groupby, title='Interactive PCA Plot')
            fig.show()
            print("Interactive PCA plot created.")
        else:
            print("No valid genes to plot for interactive PCA plot.")

    def run_analysis(self):
        """
        Run the entire analysis workflow: loading data, preprocessing, filtering genes, and generating plots.
        """
        self.load_data()
        if self.preprocess:
            self.preprocess_data()
        valid_genes, valid_gene_ids = self.check_genes()

        if 'pca_and_scree' in self.plots or 'all' in self.plots:
            self.plot_pca_and_scree()
        if 'violin' in self.plots or 'all' in self.plots:
            self.plot_violin(valid_genes, valid_gene_ids)
        if 'umap' in self.plots or 'all' in self.plots:
            self.plot_umap(valid_genes, valid_gene_ids)
        if 'scatter' in self.plots or 'all' in self.plots:
            self.plot_scatter(valid_genes, valid_gene_ids)
        if 'dotplot' in self.plots or 'all' in self.plots:
            self.plot_dotplot(valid_genes, valid_gene_ids)
        if 'interactive_heatmap' in self.plots or 'all' in self.plots:
            self.plot_interactive_heatmap(valid_genes, valid_gene_ids)
        if 'interactive_pca' in self.plots or 'all' in self.plots:
            self.plot_interactive_pca(valid_genes, valid_gene_ids)
        plt.show()

def main():
    """
    Parse command-line arguments and run the single-cell analysis.
    """
    parser = argparse.ArgumentParser(description='Single Cell Analysis Tool')
    parser.add_argument('file_path', type=str, help='Path to the input h5ad file')
    parser.add_argument('genes', type=str, help='Comma-separated list of genes to analyze')
    parser.add_argument('--groupby', type=str, default='cell_type', help='Category to group cells by')
    parser.add_argument('--preprocess', action='store_true', help='Preprocess the data')
    parser.add_argument('--n_neighbors', type=int, default=10, help='Number of neighbors for UMAP')
    parser.add_argument('--n_pcs', type=int, default=20, help='Number of principal components for PCA')
    parser.add_argument('--subset_fraction', type=float, default=0.1, help='Fraction of data to subset')
    parser.add_argument('--plots', type=str, default='all', help='Comma-separated list of plots to generate (options: pca_and_scree, violin, umap, scatter, dotplot, interactive_heatmap, interactive_pca)')
    args = parser.parse_args()

    analysis = SingleCellAnalysis(args.file_path, args.genes, args.groupby, args.preprocess, args.n_neighbors,
                                  args.n_pcs, args.subset_fraction, args.plots)
    analysis.run_analysis()

if __name__ == "__main__":
    main()

