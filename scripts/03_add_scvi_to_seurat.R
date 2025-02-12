## Set library path
.libPaths("/ihome/acillo/arc85/Rlibs_Oct_2024")

## Load reticulate
Sys.setenv(RETICULATE_PYTHON_ENV = "ret_250101")
library(reticulate)

## Load packages
library(Seurat)
library(SeuratData)
library(anndata)
library(hdf5r)
library(tidyverse)

## Parse input
cli <- commandArgs(trailingOnly=T)
args <- strsplit(cli,"=",fixed=T)
input_ser_file <- args[[1]][2]
input_adata_file <- args[[2]][2]
output_ser_file <- args[[3]][2]

## Load data
ser <- readRDS(input_ser_file)
ad <- read_h5ad(input_adata_file)

## Add new scvi results to old object
X_scvi <- ad$obsm['X_scVI'][[1]]
colnames(X_scvi) <- paste("Xscvi_,",seq(1,10,1),sep="")
rownames(X_scvi) <- rownames(ad$layers['counts'])
scvi <- CreateDimReducObject(embeddings=X_scvi,key="scvi_",assay="RNA")
ser[["x_SCVI"]] <- scvi

## Save Seurat object with scvi embeddings
saveRDS(ser,file=output_ser_file)

