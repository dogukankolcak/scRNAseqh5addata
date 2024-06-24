^# Single-Cell RNA Sequencing Analysis of NMDAR Expression in the Human Hippocampus
## Overview
This Python project focuses on analyzing single-cell RNA sequencing (scRNA-seq) data to investigate the expression of N-methyl-D-aspartate receptors (NMDARs) in the human hippocampus. The analysis is divided into two main parts:

1.SingleCellAnalysis: Analyzes scRNA-seq data from a specific region of the human hippocampus.
2.Data Merging: Combines scRNA-seq data from multiple regions of the human hippocampus to provide a comprehensive view of NMDAR expression.
Project Aim
The main goal of this project is to analyze differential NMDAR expression across different cell types in the hippocampus using scRNA-seq data. Additionally, the project aims to merge scRNA-seq data from various regions of the hippocampus to visualize NMDAR expression in specific cell types, such as neurons, across the entire hippocampus.

## Background
NMDARs are crucial for synaptic plasticity, learning, and memory in the vertebrate nervous system. These receptors are heteromeric complexes consisting of various subunits (NR1, NR2A-D, NR3A-B) and are important for processes like long-term potentiation (LTP) and long-term depression (LTD). Single-cell RNA sequencing allows scientists to examine gene expression patterns in individual cells, providing insights into cellular heterogeneity and function.

## Methods
### SingleCellAnalysis
This part of the project focuses on analyzing scRNA-seq data from a specific region of the hippocampus. The analysis involves the following steps:

1. **Data Loading:** Load scRNA-seq data in h5ad format.
2. **Gene Expression Analysis:** Check if the provided NMDAR genes are significantly expressed in the cells.
3. **Data Preprocessing:** Optionally preprocess the data (e.g., filtering out mitochondrial counts or duplicates).
4. **Visualization:** Generate various plots (e.g., PCA, Violin plots, UMAP, Scatter plots, Dot plots, Interactive PCA) to visualize gene expression.

Usage:

python SingleCellAnalysis.py rostralca3.h5ad GRIN1,GRIN2A,GRIN2B,GRIN2C,GRIN2D,GRIN3A,GRIN3B --groupby cell_type --plots all

## Data Merging
This part of the project combines scRNA-seq data from multiple regions of the hippocampus to analyze NMDAR expression across the entire hippocampus.

Usage:

sh
Kodu kopyala
python Merging_Code.py --api_files_regions "<API_LINK_1>" "<API_LINK_2>" ... --genes ENSG00000176884 ENSG00000183454 ENSG00000273079 ENSG00000161509 ENSG00000105464 ENSG00000198785 ENSG00000116032 --plot_type top_genes --cell_type neuron

## Results
The results include visualizations of NMDAR gene expression across different cell types in the hippocampus. The analysis provides insights into the distribution and expression levels of NMDARs in specific cell types, aiding in understanding their roles in synaptic function and neurological disorders.

## Discussion
The project initially aimed to analyze NMDAR expression in the entire human brain, but due to data size constraints, the focus was shifted to the hippocampus. The project involved collaboration among team members, with each contributing to different parts of the code. The use of optional preprocessing and various visualization tools helped optimize the analysis and provide comprehensive insights into NMDAR expression.

## References
For a detailed description of the methods and results, please refer to the project report included in this repository.

## How to Run
Install required dependencies:
sh
Kodu kopyala
pip install -r requirements.txt
Run the SingleCellAnalysis:
sh
Kodu kopyala
python SingleCellAnalysis.py <h5ad_file> <gene_list> --groupby <groupby> --plots <plots>
Run the Data Merging script:
sh
Kodu kopyala
python Merging_Code.py --api_files_regions "<API_LINK_1>" "<API_LINK_2>" ... --genes <gene_ids> --plot_type <plot_type> --cell_type <cell_type>

## Contact
For any questions or contributions, please contact Doğukan Kolçak.
