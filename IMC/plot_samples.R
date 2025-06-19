# Filename: plot_samples.R
# Description: Plot IMC images for select immune clusters
# Author: Xiaojun Wu
# Email: xiaojunw@usc.edu

# load libraries ####
library(dplyr)
library(SingleCellExperiment)
library(cytomapper)

# load data ####
cat("Loading data...\n")
imc.data <- readRDS("data/imc_cells_subclustered.rds")
imc.metadata <- as.data.frame(colData(imc.data))
clusters <- c("B cell", "CD8+ T cell", "CD4+ T cell", "dendritic cell",
              "plasma cell")
imc.data <- imc.data[, imc.metadata$cluster %in% clusters]
imc.metadata <- as.data.frame(colData(imc.data))  # update metadata

subject.info <- read.csv("data/subject_info.csv", row.names = 1)
subject.info <- subject.info %>% distinct() %>% arrange(pub_id)
rownames(subject.info) <- sprintf("R-%02d", subject.info$pub_id)

input.dir <- "data/masks-fig-s7h"
output.dir <- "imc-images-fig-s7h"
imc.masks <- loadImages(input.dir, as.is = TRUE)
mcols(imc.masks) <- data.frame(sample_id = names(imc.masks))

# plot samples ####
cat("Plotting samples...\n")
scale.bar.spec <- list(length = 250, label = "", lwidth = 5)
output.path <- file.path(output.dir, "image.png")
imc.images <- plotCells(imc.masks, object = imc.data, cell_id = "ObjectNumber",
                        img_id = "sample_id", colour_by = "cluster",
                        outline_by = "sample_id", missing_colour = "black",
                        return_plot = TRUE, display = "single",
                        scale_bar = scale.bar.spec, image_title = NULL,
                        save_plot = list(filename = output.path))

# rename output files
cat("Renaming output files...\n")
imc.image.ids <- names(imc.images[["plot"]])
imc.image.ids <- imc.image.ids[2:length(imc.image.ids)]
for (i in seq_along(imc.image.ids)) {
  old.filename <- sprintf(file.path(output.dir, "image_%d.png"), i)
  new.filename <- sprintf(file.path(output.dir, "%s.png"), imc.image.ids[i])

  if (file.exists(old.filename)) {
    file.rename(old.filename, new.filename)
  }
}
