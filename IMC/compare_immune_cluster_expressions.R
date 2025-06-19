# Filename: compare_immune_cluster_expressions.R
# Description: Compare IMC marker expression for immune cells using t-test
# Author: Xiaojun Wu
# Email: xiaojunw@usc.edu

# load libraries ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(SingleCellExperiment)

rng.seed <- 2025

# load data ####
cat("Loading data...\n")
subject.info <- read.csv("data/subject_info.csv", row.names = 1)
subject.info <- subject.info %>% distinct() %>% arrange(pub_id)
rownames(subject.info) <- sprintf("R-%02d", subject.info$pub_id)

imc.data <- readRDS("data/imc_cells_subclustered.rds")
imc.data <- imc.data[, colData(imc.data)$tissue_type == "tumor"]
imc.metadata <- as.data.frame(colData(imc.data))
imc.metadata$pub_id_formatted <- sprintf("R-%02d", imc.metadata$pub_id)
pub.ids.ordered <- rownames(subject.info %>% arrange(responder, pub_id))
imc.metadata <- imc.metadata %>%
  mutate(pub_id_formatted = factor(pub_id_formatted, levels = pub.ids.ordered))

output.dir <- file.path(".", "immune-cluster-expressions")

# specify clusters to compare for each marker
marker.clusters <- list(
  PDL1 = c("all", "Treg", "macrophage", "MDSC", "immature DC"),
  PD1 = c("T cell", "CD8+ T cell", "Treg"),
  "Tox/Tox2" = c("T cell", "CD8+ T cell", "NK cell"),
  pSTAT3 = c("macrophage", "pro-tumor mac", "MDSC", "Treg", "dendritic cell",
             "NK cell"),
  ARG1 = c("Treg", "M-MDSC", "mature B cell", "proliferating B cell",
           "immature DC", "exhausted CD8+ T cell"),
  GZMB = c("NK cell", "CD8+ T cell", "naive CD4+ T cell",
           "memory CD4+ T cell"),
  LAG3 = c("T cell", "CD8+ T cell", "exhausted CD8+ T cell", "Treg",
           "NK cell", "mature B cell"),
  CD137 = c("T cell", "CD8+ T cell", "B cell", "macrophage",
            "pro-inflammatory mac"),
  HLADR = c("dendritic cell", "mature B cell", "macrophage",
            "pro-inflammatory mac")
)

# map metaclusters to subclusters
meta.to.subs <- list()
meta.to.subs[["T cell"]] <- imc.metadata %>%
  filter(cluster %in% c("CD4+ T cell", "CD8+ T cell")) %>%
  pull("subcluster") %>%
  unique()
meta.to.subs[["mature B cell"]] <- c("B cell", "activated B cell")
for (mc in intersect(unlist(marker.clusters), unique(imc.metadata$cluster))) {
  meta.to.subs[[mc]] <- imc.metadata %>%
    filter(cluster == mc) %>%
    pull("subcluster") %>%
    unique()
}

# define function to filter data by marker and cluster
filter.data <- function(marker, cluster) {
  if (cluster == "all") {
    curr.data <- imc.data[marker, ]
    curr.metadata <- imc.metadata
  } else {
    if (cluster %in% names(meta.to.subs)) {
      curr.data <- imc.data[marker, imc.metadata$subcluster %in% meta.to.subs[[cluster]]]
      curr.metadata <- imc.metadata %>% filter(subcluster %in% meta.to.subs[[!!cluster]])
    } else {
      curr.data <- imc.data[marker, imc.metadata$subcluster == cluster]
      curr.metadata <- imc.metadata %>% filter(subcluster == !!cluster)
    }
  }

  list(data = curr.data, metadata = curr.metadata)
}

# define colors for time points
group.colors <- list(BL = "#a6cee3", C1D1 = "#1f78b4", WK8 = "#08306b")

# define functions to compare marker expressions for individual patients ####
# function to initialize an empty data.frame for t-test results between time
# points for individual patients
init.t.test.patient.results <- function() {
  data.frame(
    marker = character(),
    cluster = character(),
    groups = character(),
    pub_id = character(),
    group_1_size = integer(),
    group_2_size = integer(),
    hypothesis = character(),
    statistic = numeric(),
    p_value = numeric()
  )
}

# function to compare marker expressions between time points for individual
# patients
compare.patient.expressions <- function(marker.data, marker, cluster, groups) {
  t.hypotheses <- c("two.sided", "less", "greater")
  marker.metadata <- as.data.frame(colData(marker.data))
  results <- init.t.test.patient.results()
  set.seed(rng.seed)

  for (pub.id in sort(unique(marker.metadata$pub_id))) {
    pub.id.formatted <- sprintf("R-%02d", pub.id)
    cat(sprintf("Running t-tests for: %s, %s, %s vs %s, %s\n", marker, cluster,
                groups[1], groups[2], pub.id.formatted))

    group.1.data <- marker.data[, marker.metadata$pub_id == pub.id &
                                  marker.metadata$group == groups[1]]
    group.1.exprs <- assay(group.1.data, "exprs")
    group.2.data <- marker.data[, marker.metadata$pub_id == pub.id &
                                  marker.metadata$group == groups[2]]
    group.2.exprs <- assay(group.2.data, "exprs")

    if (length(group.1.exprs) < 2 || length(group.2.exprs) < 2) {
      next
    }

    for (hypothesis in t.hypotheses) {
      t.test.output <- tryCatch(
        {
          t.test(group.1.exprs, group.2.exprs, alternative = hypothesis)
        },
        error = function(e) {
          NULL
        }
      )

      if (is.null(t.test.output)) {
        next
      }

      results <- rbind(
        results,
        data.frame(
          marker = marker,
          cluster = cluster,
          groups = paste(groups, collapse = " vs "),
          pub_id = pub.id.formatted,
          group_1_size = length(group.1.exprs),
          group_2_size = length(group.2.exprs),
          hypothesis = hypothesis,
          statistic = t.test.output$statistic,
          p_value = t.test.output$p.value
        )
      )
    }
  }

  results
}

# compare marker expressions between time points for individual patients ####
cat("Comparing marker expressions between time points for individual patients...\n")
t.test.patient.results <- init.t.test.patient.results()
pdf(file.path(output.dir, "violin_indiv_patients.pdf"), width = 12, height = 6)

for (marker in names(marker.clusters)) {
  for (cluster in marker.clusters[[marker]]) {
    filtered.data <- filter.data(marker, cluster)
    curr.data <- filtered.data$data
    curr.metadata <- filtered.data$metadata

    # run t-test for each pair of time points
    t.test.patient.results <- rbind(
      t.test.patient.results,
      compare.patient.expressions(curr.data, marker, cluster, c("BL", "C1D1"))
    )
    t.test.patient.results <- rbind(
      t.test.patient.results,
      compare.patient.expressions(curr.data, marker, cluster, c("BL", "WK8"))
    )
    t.test.patient.results <- rbind(
      t.test.patient.results,
      compare.patient.expressions(curr.data, marker, cluster, c("C1D1", "WK8"))
    )

    # make combined violin plots for current marker and cluster for all patients
    # and time points
    cat(sprintf("Making violin plots for: %s, %s\n", marker, cluster))
    p <- ggplot(
      curr.metadata,
      aes(x = group, y = as.vector(assay(curr.data, "exprs")), fill = group)
    ) +
      geom_violin(trim = FALSE) +
      labs(
        title = sprintf("Marker: %s, cluster: %s", marker, cluster),
        x = "Patient - time point",
        y = "Expression",
        fill = "Time point"
      ) +
      theme_minimal() +
      facet_wrap(~pub_id_formatted, scales = "free_x", nrow = 1) +
      scale_x_discrete(limits = c("BL", "C1D1", "WK8"))
    print(p)
  }
}

dev.off()
write.csv(t.test.patient.results,
          file.path(output.dir, "t_test_indiv_patients.csv"),
          row.names = FALSE)

# define functions to compare marker expressions for responders or non-responders ####
# function to initialize an empty data.frame for t-test results between time
# points for either responders or non-responders
init.t.test.response.results <- function() {
  data.frame(
    marker = character(),
    cluster = character(),
    groups = character(),
    responder = logical(),
    group_1_size = integer(),
    group_2_size = integer(),
    hypothesis = character(),
    statistic = numeric(),
    p_value = numeric(),
    p_value_adj = numeric()
  )
}

# function to compare marker expressions between time points for responders or
# non-responders
compare.response.expressions <- function(marker.data, marker, cluster, groups) {
  t.hypotheses <- c("two.sided", "less", "greater")
  marker.metadata <- as.data.frame(colData(marker.data))
  results <- init.t.test.response.results()
  set.seed(rng.seed)

  for (is.responder in c(TRUE, FALSE)) {
    cat(sprintf("Running t-tests for: %s, %s, %s vs %s, %s\n", marker, cluster,
                groups[1], groups[2],
                ifelse(is.responder, "responders", "non-responders")))

    is.cell.selected <- replace_na(marker.metadata$responder == is.responder,
                                   FALSE)  # select cell by response status
    group.1.data <- marker.data[, marker.metadata$group == groups[1] &
                                  is.cell.selected]
    group.1.exprs <- assay(group.1.data, "exprs")
    group.2.data <- marker.data[, marker.metadata$group == groups[2] &
                                  is.cell.selected]
    group.2.exprs <- assay(group.2.data, "exprs")

    if (length(group.1.exprs) < 2 || length(group.2.exprs) < 2) {
      next
    }

    for (hypothesis in t.hypotheses) {
      t.test.output <- tryCatch(
        {
          t.test(group.1.exprs, group.2.exprs, alternative = hypothesis)
        },
        error = function(e) {
          NULL
        }
      )

      if (is.null(t.test.output)) {
        next
      }

      results <- rbind(
        results,
        data.frame(
          marker = marker,
          cluster = cluster,
          groups = paste(groups, collapse = " vs "),
          responder = is.responder,
          group_1_size = length(group.1.exprs),
          group_2_size = length(group.2.exprs),
          hypothesis = hypothesis,
          statistic = t.test.output$statistic,
          p_value = t.test.output$p.value,
          p_value_adj = NA
        )
      )
    }
  }

  results
}

# compare marker expressions between time points for responders and non-responders ####
cat("Comparing marker expressions between time points for responders and non-responders...\n")
t.test.response.results <- init.t.test.response.results()
pdf(file.path(output.dir, "violin_response_status.pdf"), width = 12, height = 6)
test.groups <- list(c("BL", "C1D1"), c("BL", "WK8"), c("C1D1", "WK8"))

for (marker in names(marker.clusters)) {
  for (cluster in marker.clusters[[marker]]) {
    filtered.data <- filter.data(marker, cluster)
    curr.data <- filtered.data$data
    curr.metadata <- filtered.data$metadata

    # run t-test for each pair of time points
    curr.results <- do.call(
      rbind,
      lapply(
        test.groups,
        function(groups) {
          compare.response.expressions(curr.data, marker, cluster, groups)
        }
      )
    )
    # adjust p-values for multiple testing
    curr.results$p_value_adj <- p.adjust(curr.results$p_value,
                                         method = "bonferroni")
    t.test.response.results <- rbind(t.test.response.results, curr.results)

    # make combined violin plots for current marker and cluster for all
    # response status and time points
    cat(sprintf("Making violin plots for: %s, %s\n", marker, cluster))
    p <- ggplot(
      curr.metadata %>% mutate(responder = factor(
        ifelse(is.na(responder), "unknown",
               ifelse(responder, "responder", "non-responder")),
        levels = c("responder", "non-responder", "unknown")
      )),
      aes(x = group, y = as.vector(assay(curr.data, "exprs")), fill = group)
    ) +
      geom_violin(trim = FALSE) +
      labs(
        title = sprintf("Marker: %s, cluster: %s", marker, cluster),
        x = "Patient - time point",
        y = "Expression",
        fill = "Time point"
      ) +
      theme_minimal() +
      facet_wrap(~responder, scales = "free_x", nrow = 1) +
      scale_x_discrete(limits = c("BL", "C1D1", "WK8")) +
      scale_fill_manual(values = group.colors)
    print(p)
  }
}

dev.off()
write.csv(t.test.response.results,
          file.path(output.dir, "t_test_response_status.csv"),
          row.names = FALSE)

# define functions to compare marker expressions between responders and non-responders ####
# function to initialize an empty data.frame for t-test results between
# responders and non-responders at each time point
init.t.test.between.responses.results <- function() {
  data.frame(
    marker = character(),
    cluster = character(),
    group = character(),
    responder_size = integer(),
    non_responder_size = integer(),
    hypothesis = character(),
    statistic = numeric(),
    p_value = numeric(),
    p_value_adj = numeric()
  )
}

# function to compare marker expressions between responders and non-responders
compare.expressions.between.responses <- function(marker.data, marker, cluster) {
  t.hypotheses <- c("two.sided", "less", "greater")
  marker.metadata <- as.data.frame(colData(marker.data))
  results <- init.t.test.between.responses.results()
  set.seed(rng.seed)

  # filter out patients without known response status
  marker.data <- marker.data[, !is.na(marker.metadata$responder)]
  marker.metadata <- marker.metadata[!is.na(marker.metadata$responder), ]

  for (group in sort(unique(marker.metadata$group))) {
    cat(sprintf("Running t-tests for: %s, %s, responders vs non-responders, at %s\n",
                marker, cluster, group))

    responder.data <- marker.data[, marker.metadata$group == group &
                                    marker.metadata$responder == TRUE]
    responder.exprs <- assay(responder.data, "exprs")
    non.responder.data <- marker.data[, marker.metadata$group == group &
                                        marker.metadata$responder == FALSE]
    non.responder.exprs <- assay(non.responder.data, "exprs")

    if (length(responder.exprs) < 2 || length(non.responder.exprs) < 2) {
      next
    }

    for (hypothesis in t.hypotheses) {
      t.test.output <- tryCatch(
        {
          t.test(responder.exprs, non.responder.exprs, alternative = hypothesis)
        },
        error = function(e) {
          NULL
        }
      )

      if (is.null(t.test.output)) {
        next
      }

      results <- rbind(
        results,
        data.frame(
          marker = marker,
          cluster = cluster,
          group = group,
          responder_size = length(responder.exprs),
          non_responder_size = length(non.responder.exprs),
          hypothesis = hypothesis,
          statistic = t.test.output$statistic,
          p_value = t.test.output$p.value,
          p_value_adj = NA
        )
      )
    }
  }

  # adjust p-values for multiple testing
  results$p_value_adj <- p.adjust(results$p_value, method = "bonferroni")

  results
}

# compare marker expressions between responders and non-responders at each time point ####
cat("Comparing marker expressions between responders and non-responders at each time point...\n")
t.test.between.responses.results <- init.t.test.between.responses.results()

for (marker in names(marker.clusters)) {
  for (cluster in marker.clusters[[marker]]) {
    filtered.data <- filter.data(marker, cluster)
    curr.data <- filtered.data$data

    # run t-test for current marker and cluster
    t.test.between.responses.results <- rbind(
      t.test.between.responses.results,
      compare.expressions.between.responses(curr.data, marker, cluster)
    )
  }
}

write.csv(t.test.between.responses.results,
          file.path(output.dir, "t_test_between_responses.csv"),
          row.names = FALSE)
