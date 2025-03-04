## Workflow for scvi integration of RNAseq data

## Pre-requisites 
- Reticulate environment installed in R
- Conda environment with scvi tools and other packages installed, which can be created using the enviornment.yml file provided in this repository. The following code will create a python environment called 'scvi-env' using conda:
	```
	conda env create -f environment.yml
	```
- Preprocessed Seurat object with a column in meta.data representing the batch effect to correct for

## Start an sbatch job on the GPU cluster
```
salloc --time=0-02:00:00 --gres=gpu:1 --cluster=gpu --partition=A100
ssh arc85@gpu-n35 ## Assumes the job is assigned to gpu-n35
```

## Load R modules
```
module load python/3.7.0 gcc/12.2.0 r/4.4.0
```

## Create anndata from pre-processed Seurat object
```
Rscript --vanilla 01_create_anndata_from_seurat.R \
	input_file=/ix1/acillo/arc85/03_workspace_241010/06_scvi_integration/03_output/pbmcsca_ser_obj.rds \
	output_file=/ix1/acillo/arc85/03_workspace_241010/06_scvi_integration/03_output/pbmcsca_adata_obj_250212.h5ad
```

## Run scvi integration
```
conda activate scvi-env

python3 02_run_scvi.py \
	input_file=/ix1/acillo/arc85/03_workspace_241010/06_scvi_integration/03_output/pbmcsca_adata_obj_250212.h5ad \
	integration_var=Method \
	model_path_save=/ix1/acillo/arc85/03_workspace_241010/06_scvi_integration/03_output/pbmcsca_scvi_model_obj_250212 \
	output_file=/ix1/acillo/arc85/03_workspace_241010/06_scvi_integration/03_output/pbmcsca_adata_scvi_250212.h5ad

conda deactivate
```

## Add scvi latent space back to Seurat object 
```
Rscript --vanilla 03_add_scvi_to_seurat.R \
	input_ser_file=/ix1/acillo/arc85/03_workspace_241010/06_scvi_integration/03_output/pbmcsca_ser_obj.rds \
	input_adata_file=/ix1/acillo/arc85/03_workspace_241010/06_scvi_integration/03_output/pbmcsca_adata_scvi_250212.h5ad \
	output_ser_file=/ix1/acillo/arc85/03_workspace_241010/06_scvi_integration/03_output/pbmcsca_ser_scvi_250212.rds
```