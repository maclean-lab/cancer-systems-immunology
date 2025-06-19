# Filename: convert_to_anndata.R
# Description: Merge two IMC data objects into one. Multiple data fields are
# added/updated/corrected, including marker names, patient id, time point
# (group), batch, publication ID, ROI, tissue type, and metacluster names.
# Author: Xiaojun Wu
# Email: xiaojunw@usc.edu

# load libraries
library(dplyr)
library(SingleCellExperiment)

# load data
cat("Loading data...\n")
# batch 1; original file name: RoussosTorres_IMC_Cells.rds
imc.data.1 <- readRDS("data/imc_cells_1.rds")
# batch 2; original file name: breast_mets_data_Cells.rds
imc.data.2 <- readRDS("data/imc_cells_2.rds")
cat("Data loaded.\n")

# unify marker names
cat("Updating marker names...\n")
rownames(imc.data.1) <- recode(
  rownames(imc.data.1),
  "Ruthenium" = "Ruthenium1",
  "DCSIGN" = "DC-SIGN",
  "DCLAMP" = "DC-LAMP",
  "ToxTox2" = "Tox/Tox2",
  "191Ir" = "Ir191",
  "193Ir" = "Ir193",
)
rownames(imc.data.2) <- recode(
  rownames(imc.data.2),
  "Tox_Tox2" = "Tox/Tox2"
)
common.markers <- intersect(rownames(imc.data.1), rownames(imc.data.2))
cat("Marker names updated.\n")

# remove all columns in rowData for both data.1 and data.2
cat("Removing all columns in rowData...\n")
rowData(imc.data.1) <- rowData(imc.data.1)[, 0]
rowData(imc.data.2) <- rowData(imc.data.2)[, 0]
cat("Columns in rowData removed.\n")

# update patient and group (time point) info both batches
cat("Updating patient and group (time point) info...\n")
colData(imc.data.1)$group <- recode(colData(imc.data.1)$group, "Wk8" = "WK8", "Pre" = "BL")
colData(imc.data.2)$Patient <- sapply(strsplit(colData(imc.data.2)$sample_id, "_"), "[", 1)
colData(imc.data.2)$group <- sapply(strsplit(colData(imc.data.2)$sample_id, "_"), "[", 2)
colData(imc.data.2)$Patient <- recode(colData(imc.data.2)$Patient, "42" = "042")
cat("Patient and group (time point) info updated.\n")

# add batch info
cat("Adding batch info...\n")
colData(imc.data.1)$batch <- 1
colData(imc.data.2)$batch <- 2
cat("Batch info added.\n")

# merge data
cat("Merging batch 1 and batch 2 data into a single data object...\n")
imc.data <- cbind(imc.data.1[common.markers, ], imc.data.2[common.markers, ])
cat("Data merged.\n")

# unify metacluster names
cat("Updating metacluster names...\n")
cluster.col.index <- which(colnames(colData(imc.data)) == "Clusters")
colnames(colData(imc.data))[cluster.col.index] <- "cluster"
colData(imc.data)$cluster <- recode(
  colData(imc.data)$cluster,
  "B cells" = "B cell",
  "CD4+ T cells" = "CD4+ T cell",
  "CD8+ T cells" = "CD8+ T cell",
  "Macrophages" = "macrophage",
  "Neutrophils" = "classical monocyte",
  "NK cells" = "NK cell",
  "Epithelial cells" = "epithelial cell",
  "Stromal cells" = "stromal cell",
  "MDSCs" = "MDSC",
  "Dendritic cells" = "dendritic cell",
  "CD4+ T-cells" = "CD4+ T cell",
  "CD8+ T-cells" = "CD8+ T cell"
)
cat("Metacluster names updated.\n")

# create a new patient ID column and a publication ID column
cat("Adding numeric patient IDs to colData...\n")
colData(imc.data)$patient_id <- gsub("-", "_", colData(imc.data)$Patient)
# keep patient ID after last "_"
colData(imc.data)$patient_id <- sub(".*_", "", colData(imc.data)$patient_id)
colData(imc.data)$patient_id <- as.numeric(colData(imc.data)$patient_id)
cat("Numeric patient IDs added.\n")

# map patient IDs to publication IDs
cat("Adding publication IDs to colData\n")
subject.ids <- read.csv("subject_ids.csv")
subject.ids$Sample <- as.numeric(sub(".*-", "", subject.ids$Sample))
subject.ids$Publication.ID <- as.numeric(sub("R-", "", subject.ids$Publication.ID))
sample.to.pub <- setNames(subject.ids$Publication.ID, subject.ids$Sample)
colData(imc.data)$pub_id <- recode(colData(imc.data)$patient_id,
                                   !!!sample.to.pub, .default = 0)
cat("Publication IDs added.\n")

# extract ROI info
cat("Adding ROI to colData...\n")
colData(imc.data)$roi <- as.numeric(sub(".*_", "", colData(imc.data)$sample_id))
cat("ROI added.\n")

# add a new column for sample IDs based publication ID
cat("Creating new sample IDs based on publication ID...\n")
sample.id.col.index <- which(colnames(colData(imc.data)) == "sample_id")
colnames(colData(imc.data))[sample.id.col.index] <- "sample_id_old"
colData(imc.data)$sample_id <- paste0(
  "R-", sprintf("%02d", colData(imc.data)$pub_id), "_",
  colData(imc.data)$group, "_", sprintf("%03d", colData(imc.data)$roi)
)
cat("New sample IDs created.\n")

# add tissue type
cat("Adding tissue types (normal vs tumor) to colData...\n")
normal.tissues <- data.frame(
  patient_id = c(27, 36, 42, 44, 44, 44),
  group = c("C1D1", "C1D1", "BL", "BL", "C1D1", "C1D1"),
  ROI = c(6, 4, 2, 3, 1, 2)
)

colData(imc.data)$tissue_type <- apply(colData(imc.data), 1, function(row) {
  if (any(normal.tissues$patient_id == row["patient_id"] &
            normal.tissues$group == row["group"] &
            normal.tissues$ROI == row["roi"])) {
    "normal"
  } else {
    "tumor"
  }
})
cat("Tissue types added.\n")

# Add responder vs non-responder info
cat("Adding responder vs non-responder info to colData...\n")
subject.info <- read.csv("data/subject_info.csv", row.names = 1)
responder.info <- subject.info %>% distinct(pub_id, responder)
rownames(responder.info) <- responder.info$pub_id
colData(imc.data)$responder <- responder.info[as.character(colData(imc.data)$pub_id), "responder"]
cat("Responder vs non-responder info added.\n")

# save merged data
cat("Saving merged data...\n")
saveRDS(imc.data, "data/imc_cells_merged.rds")
cat("Data saved.\n")
