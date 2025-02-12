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
input_file <- args[[1]][2]
output_file <- args[[2]][2]

## Load data
ser <- readRDS(input_file)

## Save as anndata
mat <- GetAssayData(ser,layer="counts",assay="RNA") %>% t()
meta <- ser@meta.data
pca <- Embeddings(ser,"pca")
umap <- Embeddings(ser,"umap")

ad <- AnnData(
  X = mat,
  obs = meta,
  var = data.frame(gene=rownames(ser)),
  layers = list(
  ),
  obsm = list(
      pca = pca,
      umap = umap
  ),
  varm = list(
  ),
  uns = list(
  )
)

## Write out anndata
write_h5ad(anndata=ad,filename=output_file)
