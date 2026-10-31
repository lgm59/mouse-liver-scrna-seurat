# Single-cell RNA-seq clustering of mouse liver fibrosis data (R / Seurat)
# Data: Fabre et al., Broad Single Cell Portal SCP2053
# Place matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz in the working directory.

### Run once, then leave commented out
# install.packages("Seurat")
# install.packages("dplyr")
# install.packages("ggplot2")
library(Seurat)
library(dplyr)
library(ggplot2)
###

data <- Read10X(data.dir = ".")
seurat_obj <- CreateSeuratObject(counts = data,
                                 project = "CompBioProj",
                                 min.cells = 3,
                                 min.features = 200)

seurat_obj

# Mouse mitochondrial genes start with "mt-"
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# Run PCA using the 2000 variable features
seurat_obj <- RunPCA(seurat_obj,
                     features = VariableFeatures(object = seurat_obj)) # NORMALLY features = all.genes
# Visualize top genes driving each PC
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
# Visualize cells in PC space
DimPlot(seurat_obj, reduction = "pca")
# Elbow plot - look for the "elbow"
ElbowPlot(seurat_obj, ndims = 50)

# Find neighbors in PC space (using first 20 PCs)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
table(seurat_obj$seurat_clusters)

# Run UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
# Plot clusters
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

FeaturePlot(seurat_obj, features = c("Cd14", "Cd68", "Fcgr3"))
