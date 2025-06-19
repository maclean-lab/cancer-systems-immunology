# Filename: convert_to_anndata.R
# Description: Convert IMC data object (SingleCellExperiment) to AnnData
# Author: Xiaojun Wu
# Email: xiaojunw@usc.edu

library(SingleCellExperiment)
library(zellkonverter)

# load data
imc.data <- readRDS("data/imc_cells_subclustered.rds")

# add rownames of colData as a column
colData(imc.data)$cell_id <- rownames(colData(imc.data))

# write data to h5ad
writeH5AD(imc.data, "data/imc_cells_subclustered.h5ad")
