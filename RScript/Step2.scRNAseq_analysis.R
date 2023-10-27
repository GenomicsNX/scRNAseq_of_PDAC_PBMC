library(Seurat)

# Remove genes that have zero values across all cells
counts.m=counts.m[rowSums(counts.m)>0,]

# Normalize data, select highly variable genes (HVGs), and perform clustering with visualization using the Seurat package
num_hvgs=2000
num_principal_components=30

pbmc_seurat_obj = CreateSeuratObject(counts=counts.m)
pbmc_seurat_obj = pbmc_seurat_obj %>% NormalizeData() %>% FindVariableFeatures(nfeatures=num_hvgs) %>% ScaleData() %>% RunPCA(npcs=num_principal_components) 
pbmc_seurat_obj = pbmc_seurat_obj %>% FindNeighbors(dims=1:num_principal_components) %>% FindClusters() %>% RunUMAP(dims = 1:num_principal_components) %>% RunTSNE(dims = 1:num_principal_components)

"""
# Perform batch correction using the integration workflow within the Seurat package
# This process referred to the "https://satijalab.org/seurat/archive/v3.0/integration.html"
"""

# Divide the dataset into a list of sample-level data
split_obj <- SplitObject(pbmc_seurat_obj, split.by = "sample")

# Independently perform normalization and identification of variable features for each sample dataset
split_obj <- lapply(X=split_obj, FUN=function(x) {
    x = x %>% NormalizeData(x) %>% FindVariableFeatures(nfeatures=num_hvgs)
})

anchors <- FindIntegrationAnchors(object.list=split_obj, dims=1:num_principal_components)
integrated_obj <- IntegrateData(anchorset=anchors, dims=1:num_principal_components)

# Perform the standard workflow for data visualization and clustering analysis within the Seurat pipeline
integrated_obj <- integrated_obj %>% ScaleData() %>% RunPCA(npcs=num_principal_components) %>% RunUMAP(reduction="pca", dims=1:num_principal_components) %>% RunTSNE(reduction="pca", dims=1:num_principal_components)
