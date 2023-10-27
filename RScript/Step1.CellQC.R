library(DropletUtils)
library(scater)
library(dplyr)

setwd(path) # path: Includes subdirectories containing Cell Ranger count outputs for each sample
sample.list=dir(path)

# Set a threshold for cell quality control (QC) to remove empty droplets and low-quality cells
emptydropfilter.threshold=0.01
threshold.log10UMIs=3
threshold.pct_mito=10

counts.m=NULL
i=0

for (sample in sample.list){
  setwd(path)
  data.path=paste(path,sample,sep="/")
  
  # Load data
  umi.counts <- Seurat::Read10X_h5(filename = paste(data.path,"raw_feature_bc_matrix.h5",sep="/"), use.names = TRUE)
  colnames(umi.counts)=paste(sample,colnames(umi.counts),sep="-")
  sce <- SingleCellExperiment(list(counts=umi.counts))
  e.out <- emptyDrops(counts(sce))
  filt.sce <- sce[,which(e.out$FDR<emptydropfilter.threshold)]
  
  # Cell QC: Calculate the ratio of mitochondrially-mapped reads in the genome
  is.mito <- grepl("^MT-", rownames(filt.sce))
  colData(filt.sce) <- perCellQCMetrics(filt.sce, subsets=list(Mt=is.mito))
  filt.sce$log10UMIs=log10(filt.sce$sum)
  filt.sce$log10nGenes=log10(filt.sce$detected)
   
  colnames(filt.sce) <- paste(sample, colnames(filt.sce), sep="-")
  filt.sce2=filt.sce[,filt.sce$log10UMIs>=threshold.log10UMIs & filt.sce$subsets_Mt_percent<=threshold.pct_mito]
  
  # Save the results
  dir.create(output_dir) # output_dir: Directory for saving the results of QC data for each sample
  setwd(output_dir)
  saveRDS(sce,file="sce.rds")
  saveRDS(filt.sce,file="filt.sce.rds")
  saveRDS(filt.sce2,file="filt.sce2.rds")
  
  if (i==0){
    counts.m=counts(filt.sce2)
  } else {counts.m=RowMergeSparseMatrices(counts.m,counts(filt.sce2))}
  i=i+1
  print (i)
}
