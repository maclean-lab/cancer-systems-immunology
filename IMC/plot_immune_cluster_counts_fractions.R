# Filename: plot_immune_cluster_counts_fractions.R
# Description: Plot cell proportions for all meta/subclusters (computed in
# compare_immune_cluster_fractions.R)
# Author: Xiaojun Wu
# Email: xiaojunw@usc.edu

# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(paletteer)
library(colorspace)
library(SingleCellExperiment)

# load data
cat("Loading data...\n")
imc.data <- readRDS("data/imc_cells_subclustered.rds")
imc.metadata <- as.data.frame(colData(imc.data))
imc.metadata$pub_id_formatted <- sprintf("R-%02d", imc.metadata$pub_id)
subject.info <- read.csv("data/subject_info.csv", row.names = 1)
subject.info <- subject.info %>% distinct() %>% arrange(pub_id)
rownames(subject.info) <- sprintf("R-%02d", subject.info$pub_id)
output.dir <- file.path(".", "immune-cluster-counts-fractions")

# exclude normal samples and include specific metaclusters only
metaclusters <- c("CD4+ T cell", "CD8+ T cell", "NK cell", "MDSC",
                  "macrophage", "dendritic cell", "B cell", "plasma cell",
                  "neutrophil", "monocyte")
metacluster.colors <- c("CD4+ T cell" = "#b4b637",
                        "CD8+ T cell" = "#c22b2c",
                        "NK cell" = "#c970a1",
                        "MDSC" = "#866099",
                        "macrophage" = "#7e5146",
                        "dendritic cell" = "#2b9746",
                        "B cell" = "#286ca0",
                        "plasma cell" = "#7ec6bb",
                        "neutrophil" = "#e6762b",
                        "monocyte" = "#f3b183")  # manually assigned colors
imc.metadata <- imc.metadata %>%
  filter(tissue_type == "tumor", cluster %in% metaclusters)

# define function to compute cell counts and fractions
get.cell.counts <- function(metadata, cluster.col = "cluster") {
  cluster.col.sym <- rlang::sym(cluster.col)
  metadata %>%
    dplyr::count(pub_id, group, !!cluster.col.sym) %>%
    tidyr::pivot_wider(names_from = !!cluster.col.sym, values_from = n,
                       values_fill = 0) %>%
    mutate(pub_id = sprintf("R-%02d", pub_id))
}

get.cell.fracs <- function(counts) {
  counts %>%
    rowwise() %>%
    mutate(across(3:ncol(counts), ~ .x / sum(c_across(3:ncol(counts))))) %>%
    ungroup()
}

# count and plot immune cells at metacluster level (one page per patient) ####
# count immune cells per patient and time point
cat("Counting cells in metaclusters...\n")

metacluster.counts <- get.cell.counts(imc.metadata)
metacluster.fracs <- get.cell.fracs(metacluster.counts)

write.csv(metacluster.counts, file.path(output.dir, "metacluster_counts.csv"))
write.csv(metacluster.fracs, file.path(output.dir, "metacluster_fractions.csv"))

# plot fractions of immune metacluster per patient at each time point
cat("Plotting immune metacluster fractions...\n")
# NOTE: manually assigned colors are used for this plot
figure.path <- file.path(output.dir, "metacluster_fractions.pdf")
pdf(figure.path, width = 6, height = 4)

for (patient in rownames(subject.info)) {
  pub.id <- subject.info[patient, "pub_id"]
  patient.metadata <- imc.metadata %>% filter(pub_id == pub.id)

  tissue.type <- subject.info[patient, "tissue"]
  if (is.na(subject.info[patient, "responder"])) {
    responder.str <- "response unknown"
  } else {
    responder.str <- ifelse(subject.info[patient, "responder"],
                            "responder", "non-responder")
  }
  plot.title <- paste0("Patient ", patient, "\n", tissue.type, ", ",
                       responder.str)

  p <- ggplot(patient.metadata, aes(group, fill = cluster)) +
    geom_bar(position = "fill", width = 0.8) +
    scale_fill_manual(values = metacluster.colors) +
    labs(title = plot.title, x = "time point", y = "fraction of cells") +
    theme_minimal()


  print(p)
}

dev.off()

# update colors and order of metaclusters and subclusters ####
cat("Generating colors for metaclusters and subclusters from pallette...\n")
metacluster.labels <- imc.metadata %>%
  distinct(cluster) %>%
  arrange(cluster) %>%
  select(cluster) %>%
  pull()
imc.metadata$cluster <- factor(imc.metadata$cluster, levels = metacluster.labels)

subcluster.labels <- imc.metadata %>%
  filter(cluster %in% metaclusters) %>%
  distinct(cluster, subcluster) %>%
  arrange(cluster) %>%
  select(subcluster) %>%
  pull()
imc.metadata$subcluster <- factor(imc.metadata$subcluster, levels = subcluster.labels)

# create color scheme for subclusters
metacluster.hues <- paletteer_d("ggthemes::Classic_10")
metacluster.hues <- setNames(metacluster.hues, metacluster.labels)

generate.subcluster.colors <- function(mc, scs) {
  base.color <- metacluster.hues[mc]

  if (length(scs) == 1) {
    return(setNames(base.color, scs))
  }

  alpha.vals <- seq(1, 0.4, length.out = length(scs))
  shades <- sapply(alpha.vals, function(a) alpha(base.color, a))
  setNames(shades, scs)
}

subcluster.colors <- unlist(
  lapply(metaclusters, function(mc) {
    scs <- imc.metadata %>%
      filter(cluster == mc) %>%
      distinct(subcluster) %>%
      pull(subcluster)
    generate.subcluster.colors(mc, scs)
  })
)
subcluster.colors <- subcluster.colors[subcluster.labels]

# count and plot immune cells at subcluster level (one page per patient) ####
# count immune cells per patient and time point
cat("Counting cells in immune subclusters...\n")

subcluster.counts <- get.cell.counts(imc.metadata, cluster.col = "subcluster")
subcluster.fracs <- get.cell.fracs(subcluster.counts)

write.csv(subcluster.counts, file.path(output.dir, "subcluster_counts.csv"))
write.csv(subcluster.fracs, file.path(output.dir, "subcluster_fractions.csv"))

figure.path <- file.path(output.dir, "subcluster_fractions.pdf")
pdf(figure.path, width = 6, height = 4)

# plot fractions of immune subcluster per patient at each time point
cat("Plotting subcluster fractions ...\n")
for (patient in rownames(subject.info)) {
  pub.id <- subject.info[patient, "pub_id"]
  patient.metadata <- imc.metadata %>% filter(pub_id == pub.id)

  tissue.type <- subject.info[patient, "tissue"]
  if (is.na(subject.info[patient, "responder"])) {
    responder.str <- "response unknown"
  } else {
    responder.str <- ifelse(subject.info[patient, "responder"],
                            "responder", "non-responder")
  }
  plot.title <- paste0("Patient ", patient, "\n", tissue.type, ", ",
                       responder.str)

  p <- ggplot(patient.metadata, aes(group, fill = subcluster)) +
    geom_bar(position = "fill", width = 0.8) +
    scale_fill_manual(values = subcluster.colors,
                      guide = guide_legend(title = "Subcluster", ncol = 2)) +
    labs(title = plot.title, x = "time point", y = "fraction of cells") +
    theme_minimal()

  print(p)
}

dev.off()

# group patients by response in metadata for combined plots ####
cat("Grouping patients by response in metadata (for combined plots)...\n")
pub.ids.ordered <- rownames(subject.info %>% arrange(responder, pub_id))
imc.metadata <- imc.metadata %>%
  mutate(pub_id_formatted = factor(pub_id_formatted, levels = pub.ids.ordered))

# Create a combined plot of metacluster fractions for all patients ####
cat("Creating combined plot of metacluster fractions for all patients...\n")
p <- ggplot(imc.metadata, aes(group, fill = cluster)) +
  geom_bar(position = "fill", width = 0.8) +
  scale_fill_manual(values = metacluster.hues,
                    guide = guide_legend(title = "Metacluster")) +
  labs(x = "Time Point", y = "Fraction of cells") +
  theme_minimal() +
  facet_wrap(~pub_id_formatted, scales = "free_x", nrow = 1)

figure.path <- file.path(output.dir, "metacluster_fractions_combined.pdf")
ggsave(figure.path, p, width = 16, height = 8)

# Create a combined plot of subcluster fractions for all patients ####
cat("Creating combined plot of subcluster fractions for all patients...\n")
p <- ggplot(imc.metadata, aes(group, fill = subcluster)) +
  geom_bar(position = "fill", width = 0.8) +
  scale_fill_manual(values = subcluster.colors,
                    guide = guide_legend(title = "Subcluster", ncol = 2)) +
  labs(x = "Time Point", y = "Fraction of cells") +
  theme_minimal() +
  facet_wrap(~pub_id_formatted, scales = "free_x", nrow = 1)

figure.path <- file.path(output.dir, "subcluster_fractions_combined.pdf")
ggsave(figure.path, p, width = 16, height = 8)
