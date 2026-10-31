# Single-cell RNA-seq clustering of mouse liver fibrosis data

R/Seurat workflow on a public dataset: Fabre et al., Broad Single Cell Portal
SCP2053 (10x 3' scRNA-seq, CCl4 and high-fat diet mouse models).

Steps: QC plots, normalization, 2,000 variable genes, PCA, clustering (30 clusters),
UMAP, macrophage marker plots (Cd14, Cd68, Fcgr3).

Data: download matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz from
https://singlecell.broadinstitute.org/single_cell/study/SCP2053 and place them in
the working directory.

Code: compbioscript.R
