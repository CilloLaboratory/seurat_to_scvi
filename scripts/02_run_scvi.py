## Import packages
import sys
import os
import numpy as np
import pandas as pd
import scvi
import scanpy as sc
import torch
import torch.nn as nn
from scipy.sparse import csr_matrix
from sklearn_ann.kneighbors.annoy import AnnoyTransformer

## Parse input variables
input_file = sys.argv[1].split("=")[1]
integration_var = sys.argv[2].split("=")[1]
model_path_save = sys.argv[3].split("=")[1]
output_file = sys.argv[4].split("=")[1]

## Check GPU availability
torch.cuda.is_available()
torch.cuda.device_count()

## Load data
adata = sc.read(input_file)

## Pre-process data
sc.pp.filter_genes(adata, min_counts=3)
adata.layers['counts'] = csr_matrix(adata.X.copy())
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

## Identify highly variable genes 
sc.pp.highly_variable_genes(
    adata,
    layer='counts',
    n_top_genes=2000,
    subset=True,
    flavor="seurat_v3",
    batch_key="Method")

## Set up model
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=[integration_var]
)
model = scvi.model.SCVI(adata)

## Train the model
model.train()

## Save the model 
model.save(model_path_save)

## Add scvi normalized data to object
SCVI_NORMALIZED_KEY = "scvi_normalized"
adata.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=10e4)

## Add latent space to object 
SCVI_LATENT_KEY = "X_scVI"
latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent

## Save object as new anndata
adata.write(output_file)
