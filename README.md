##Running the Analysis Pipeline

The following Python script performs the analysis steps. Make sure to adjust the paths to your data files as needed.

# Import Libraries
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# Activate the pandas2ri conversion
pandas2ri.activate()

# Import Seurat if using it
seurat = importr('Seurat')

# Load the data (assuming 10X Genomics format)
adata = sc.read_10x_mtx('path_to_filtered_feature_bc_matrix')

# Perform basic QC filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter cells based on QC metrics
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.total_counts < 10000, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# Scale the data
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

# Compute the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# UMAP
sc.tl.umap(adata)
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])

# Clustering
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden'])

# Rank genes by groups
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# Save processed data
adata.write('processed_data.h5ad')

# Generate QC report
qc_report = adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']]
qc_report.to_csv('qc_report.csv')

# Save clustering results
clustering_results = adata.obs['leiden']
clustering_results.to_csv('clustering_results.csv')

# Save differential expression results
de_results = sc.get.rank_genes_groups_df(adata, group=None)
de_results.to_csv('de_results.csv')

# Save UMAP plot
sc.pl.umap(adata, color=['leiden'], save='umap_plot.png')


#Pipeline Overview
1)Quality Control and Preprocessing: Filters cells and genes, calculates QC metrics, and filters cells based on QC criteria.

2)Normalization and Scaling: Normalizes the data, identifies highly variable genes, and scales the data.

3)Clustering and Dimensionality Reduction: Performs PCA, computes the neighborhood graph, runs UMAP, and clusters cells.

4)Differential Expression Analysis: Identifies differentially expressed genes between clusters.

5)Save Results and Generate Reports: Saves the processed data, QC report, clustering results, differential expression results, and UMAP plots.

##Output
The pipeline generates several output files:

processed_data.h5ad: Processed AnnData object.
qc_report.csv: Quality control report.
clustering_results.csv: Clustering results.
de_results.csv: Differential expression results.
umap_plot.png: UMAP plot.