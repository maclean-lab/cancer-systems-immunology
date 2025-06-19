# load libraries ####
library(dplyr)
library(CATALYST)
library(circlize)
library(ComplexHeatmap)

# load data ####
imc.data <- readRDS("data/imc_cells_merged.rds")
rng.seed <- 88
set.seed(rng.seed)

# define cell subtypes and their markers
subtypes <- list(
  "MDSC" = c("G-MDSC", "M-MDSC"),
  "NK cell" = c("immature NK cell", "cytotoxic NK cell"),
  "CD4+ T cell" = c("Treg", "naive CD4+ T cell", "memory T helper cell"),
  "CD8+ T cell" = c("naive CD8+ T cell", "memory CD8+ T cell",
                    "exhausted CD8+ T cell", "effector CD8+ T cell",
                    "proliferating CD8+ T cell"),
  "B cell" = "B cell",
  "macrophage" = c("suppressive macrophage", "proliferating macrophage"),
  "dendritic cell" = c("immature DC", "mature DC"),
  "epithelial cell" = "epithelial cell",
  "stromal cell" = "stromal cell",
  "classical monocyte" = "classical monocyte"
)
subtype.markers <- list(
  "G-MDSC" = c("CD45", "CD33", "CD15", "PDL1", "PD1", "Tox/Tox2", "pSTAT3", "ARG1"),
  "M-MDSC" = c("CD45", "CD33", "CD14", "PDL1", "PD1", "Tox/Tox2", "pSTAT3", "ARG1"),
  "immature NK cell" = c("CD45", "CD57", "PD1", "pSTAT3"),
  "cytotoxic NK cell" = c("CD45", "CD57", "GZMB", "CD16", "PD1", "Tox/Tox2", "pSTAT3", "LAG3", "CD137"),
  "Treg" = c("CD45", "CD3", "CD4", "FOXP3", "PDL1", "DC-LAMP", "PD1", "Tox/Tox2", "pSTAT3", "LAG3", "CD137"),
  "naive CD4+ T cell" = c("CD45", "CD3", "CD4", "CD45RA", "DC-LAMP"),
  "memory T helper cell" = c("CD45", "CD3", "CD4", "CD45RO", "CD57", "DC-LAMP", "PD1", "FOXP3", "LAG3", "CD137", "CD86"),
  "naive CD8+ T cell" = c("CD45", "CD3", "CD8", "CD45RA", "DC-LAMP"),
  "memory CD8+ T cell" = c("CD45", "CD3", "CD8", "CD45RO", "CD57", "DC-LAMP", "PD1", "FOXP3", "LAG3", "CD137", "CD86"),
  "exhausted CD8+ T cell" = c("CD45", "CD3", "CD8", "LAG3", "Tox/Tox2", "PD1", "PDL1", "CD57", "FOXP3", "pSTAT3"),
  "effector CD8+ T cell" = c("CD45", "CD3", "CD8", "CD137", "CD57", "PD1", "pSTAT3", "GZMB", "LAG3", "CD86"),
  "proliferating CD8+ T cell" = c("CD45", "CD3", "CD8", "KI67", "CD137", "PD1"),
  "plasma cell" = c("CD45", "CD138", "HLADR"),
  "suppressive macrophage" = c("CD45", "CD68", "CD163", "CD14", "CD4", "CD33", "HLADR", "CD16", "CD15", "ARG1", "PDL1", "PD1", "Tox/Tox2", "FOXP3", "pSTAT3", "ARG1", "CD86", "CD137", "CD86"),
  "proliferating macrophage" = c("CD45", "CD68", "CD163", "CD14", "HLADR", "KI67", "ARG1", "CD137", "CD86"),
  "immature DC" = c("CD45", "DC-SIGN", "HLADR", "PDL1", "PD1", "FOXP3", "pSTAT3"),
  "mature DC" = c("CD45", "DC-SIGN", "HLADR", "DC-LAMP", "CD86", "PDL1", "CD4", "CD68", "PD1", "FOXP3", "pSTAT3", "CD86", "ARG1", "LAG3", "CD137"),
  "B cell" = c("CD45", "CD20", "HLADR", "CD86", "KI67", "PDL1", "PD1", "LAG3", "CD137", "CD45RA", "CD45RO"),
  "epithelial cell" = c("Ecadherin", "CK", "Podoplanin", "DC-LAMP", "HLADR", "PDL1", "pSTAT3"),
  "stromal cell" = c("SMA", "VIM", "Podoplanin", "Collagen", "CD4", "CD15", "HLADR", "PDL1", "pSTAT3", "ARG1"),
  "classical monocyte" = c("CD45", "CD14", "CD4", "CD163", "CD68", "CD33", "HLADR", "CD16", "CD15", "DC-SIGN", "PDL1")
)

# create container for sub-clustering results
subcluster.data <- list()

# define functions to for sub-clustering ####
run.subclustering <- function(metacluster) {
  cat(sprintf("Running subclustering for %s...\n", metacluster))

  markers <- unique(unlist(subtype.markers[subtypes[[metacluster]]]))
  CATALYST::cluster(
    imc.data[, imc.data$cluster == metacluster],
    features = markers, xdim = 3, ydim = 3, maxK = 5, seed = rng.seed
  )
}

plot.subcluster.heatmap <- function(metacluster, cluster.labels = NULL,
                                    markers = NULL, draw.plot = TRUE,
                                    save.plot = FALSE) {
  cat(sprintf("Plotting subcluster heatmap for %s...\n", metacluster))

  # get marker expressions for cells in the metacluster
  curr.subcluster.data <- subcluster.data[[metacluster]]
  curr.subcluster.metadata <- as.data.frame(colData(curr.subcluster.data))
  if (is.null(markers)) {
    markers <- unique(unlist(subtype.markers[subtypes[[metacluster]]]))
  }
  # NOTE: marker.exprs is a data frame with markers as columns and cells as
  # rows. cluster_id is also added as a column for grouping purpose
  marker.exprs <- data.frame(t(assay(curr.subcluster.data[markers, ], "exprs"))) %>%
    bind_cols(curr.subcluster.metadata %>% select(cluster_id)) %>%
    group_by(cluster_id)

  # computer median expressions of markers in each sub-cluster
  marker.medians <- marker.exprs %>% summarise_all(median) %>% select(-cluster_id)
  cluster.sizes <- marker.exprs %>% dplyr::count(cluster_id)

  # make heatmap for median expression of markers
  heatmap.data <- as.matrix(marker.medians)
  # add cluster ids to heatmap.data, which are lost when cast to matrix
  rownames(heatmap.data) <- rownames(marker.medians)
  plot.title <- sprintf("Subclustering of %s", metacluster)
  cluster.size.annotation <- HeatmapAnnotation(
    "Cluster size" = anno_text(
      format(cluster.sizes$n, nsmall = 0), which = "row", just = "right",
      height = unit(2, "cm"), location = 0.9,
      gp = gpar(fontsize=14, col = "black", border = "black")
    ),
    which = "row"
  )
  cluster.heatmap <- Heatmap(
    heatmap.data, colorRamp2(c(0, 2), c("white", "red")),
    "Median expression", column_title = plot.title,
    left_annotation = cluster.size.annotation
  )

  if (draw.plot) {
    draw(cluster.heatmap)
  }

  if (save.plot) {
    save.path <- paste(metacluster, "subclusters.png", sep = "_")
    save.path <- gsub(" ", "_", save.path)
    save.path <- file.path("subclusters", save.path)
    png(save.path, width = 720)
    draw(cluster.heatmap)
    dev.off()
  }
}

# run sub-clustering for MDSCs ####
subcluster.data[["MDSC"]] <- run.subclustering("MDSC")
# markers.to.plot <- unique(unlist(subtype.markers[subtypes[["MDSC"]]]))
# markers.to.plot <- c(markers.to.plot, "HLADR", "CD68", "CD163", "CD16")
markers.to.plot <- rownames(imc.data)
plot.subcluster.heatmap("MDSC", markers = markers.to.plot, draw.plot = TRUE, save.plot = FALSE)

# run sub-clustering for NK cells ####
subcluster.data[["NK cell"]] <- run.subclustering("NK cell")
plot.subcluster.heatmap("NK cell", markers = markers.to.plot, draw.plot = TRUE, save.plot = FALSE)

# run sub-clustering for CD4+ T cells ####
subcluster.data[["CD4+ T cell"]] <- run.subclustering("CD4+ T cell")
plot.subcluster.heatmap("CD4+ T cell", markers = markers.to.plot, draw.plot = TRUE, save.plot = FALSE)

# run sub-clustering for CD8+ T cells ####
subcluster.data[["CD8+ T cell"]] <- run.subclustering("CD8+ T cell")
plot.subcluster.heatmap("CD8+ T cell", markers = markers.to.plot, draw.plot = TRUE, save.plot = FALSE)

# run sub-clustering for macrophages ####
subcluster.data[["macrophage"]] <- run.subclustering("macrophage")
plot.subcluster.heatmap("macrophage", markers = markers.to.plot, draw.plot = TRUE, save.plot = FALSE)

# run sub-clustering for dendritic cells ####
subcluster.data[["dendritic cell"]] <- run.subclustering("dendritic cell")
plot.subcluster.heatmap("dendritic cell", markers = markers.to.plot, draw.plot = TRUE, save.plot = FALSE)

# run sub-clustering for B cells ####
subcluster.data[["B cell"]] <- run.subclustering("B cell")
plot.subcluster.heatmap("B cell", markers = markers.to.plot, draw.plot = TRUE, save.plot = FALSE)

# run sub-clustering for epithelial cells ####
subcluster.data[["epithelial cell"]] <- run.subclustering("epithelial cell")
plot.subcluster.heatmap("epithelial cell", markers = markers.to.plot, draw.plot = TRUE, save.plot = FALSE)

# run sub-clustering for stromal cells ####
subcluster.data[["stromal cell"]] <- run.subclustering("stromal cell")
plot.subcluster.heatmap("stromal cell", markers = markers.to.plot, draw.plot = TRUE, save.plot = FALSE)

# run sub-clustering for classical monocyte ####
subcluster.data[["classical monocyte"]] <- run.subclustering("classical monocyte")
plot.subcluster.heatmap("classical monocyte", markers = markers.to.plot, draw.plot = TRUE, save.plot = FALSE)

# add annotations back to the original data ####
cat("Adding annotations to subclusters...\n")
subcluster.annotations <- read.csv("data/subcluster_annotations.csv", row.names = 1)
colData(imc.data)$subcluster <- NA

for (metacluster in names(subcluster.data)) {
  # get subcluster metadata
  curr.subcluster.data <- subcluster.data[[metacluster]]
  curr.subcluster.metadata <- as.data.frame(colData(curr.subcluster.data))

  # add annotations
  cell.ids <- rownames(curr.subcluster.metadata)
  subcluster.ids <- paste0("C", curr.subcluster.metadata$cluster_id)
  colData(imc.data)[cell.ids, "subcluster"] <- t(subcluster.annotations[metacluster, subcluster.ids])
}

# fix metacluster labels according to subcluster annotations ####
metacluster.labels <- colData(imc.data)$cluster
subcluster.labels <- colData(imc.data)$subcluster
colData(imc.data)$cluster <- case_when(
  subcluster.labels == "neutrophil" ~ "neutrophil",
  subcluster.labels == "non-classical monocyte" ~ "monocyte",
  subcluster.labels == "pro-tumor mac" ~ "macrophage",
  subcluster.labels == "pro-inflammatory mac" ~ "macrophage",
  subcluster.labels == "Treg" ~ "CD4+ T cell",
  subcluster.labels == "G-MDSC" ~ "MDSC",
  subcluster.labels == "B cell" ~ "B cell",
  subcluster.labels == "naive B cell" ~ "B cell",
  subcluster.labels == "proliferating B cell" ~ "B cell",
  subcluster.labels == "activated B cell" ~ "B cell",
  subcluster.labels == "classical monocyte" ~ "monocyte",
  subcluster.labels == "non-classical monocyte" ~ "monocyte",
  subcluster.labels == "fibroblast" ~ "stromal cell",
  subcluster.labels == "plasma cell" ~ "plasma cell",
  subcluster.labels == "endothelial cell" ~ "endothelial cell",
  subcluster.labels == "non-cell" ~ "non-cell",
  subcluster.labels == "other" ~ "other",
  .default = metacluster.labels
)

# use metacluster labels for cells without subcluster labels ####
subcluster.labels <- colData(imc.data)$subcluster
metacluster.labels <- colData(imc.data)$cluster
colData(imc.data)$subcluster <- ifelse(is.na(subcluster.labels),
                                       metacluster.labels, subcluster.labels)

# save annotated data ####
cat("Saving annotated data...\n")
saveRDS(imc.data, "data/imc_cells_subclustered.rds")
