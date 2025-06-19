# Filename: compare_immune_cluster_fractions.R
# Description: Compute cell proportions for all meta/subclusters
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
subject.info$pub_id <- rownames(subject.info)
output.dir <- file.path(".", "immune-cluster-counts-fractions")

# load cluster and subcluster counts
metacluster.counts <- read.csv(file.path(output.dir, "metacluster_counts.csv"),
                               row.names = 1, check.names = FALSE)
subcluster.counts <- read.csv(file.path(output.dir, "subcluster_counts.csv"),
                              row.names = 1, check.names = FALSE)

# compute cell fractions from counts
get.cell.fracs <- function(counts) {
  counts %>%
    rowwise() %>%
    mutate(across(3:ncol(counts), ~ .x / sum(c_across(3:ncol(counts))))) %>%
    ungroup()
}

metacluster.fracs <- get.cell.fracs(metacluster.counts)
metacluster.fracs <- metacluster.fracs %>%
  left_join(subject.info %>% select(pub_id, responder), by = "pub_id")
subcluster.fracs <- get.cell.fracs(subcluster.counts)
subcluster.fracs <- subcluster.fracs %>%
  left_join(subject.info %>% select(pub_id, responder), by = "pub_id")

# define function to initialize an empty data.frame for t-test results between time points ####
init.t.test.group.results <- function() {
  data.frame(
    groups = character(),
    cluster = character(),
    responder = logical(),
    group_1_size = integer(),
    group_2_size = integer(),
    hypothesis = character(),
    statistic = numeric(),
    p_value = numeric()
  )
}

# define function to analyze cell fractions between two time point ####
# (separately for responders and non-responders)
compare.groups <- function(fracs, clusters, group.1, group.2, plot.bar = FALSE) {
  results <- init.t.test.group.results()
  t.hypotheses <- c("two.sided", "less", "greater")
  set.seed(rng.seed)

  for (cluster in clusters) {
    for (is.responder in c(TRUE, FALSE)) {
      group.1.fracs <- fracs %>%
        filter(responder == is.responder, group == group.1) %>%
        pull({{ cluster }})
      group.2.fracs <- fracs %>%
        filter(responder == is.responder, group == group.2) %>%
        pull({{ cluster }})

      if (length(group.1.fracs) < 2 || length(group.2.fracs) < 2) {
        next
      }

      # perform t-test
      for (hypothesis in t.hypotheses) {
        t.test.output <- t.test(group.1.fracs, group.2.fracs,
                                alternative = hypothesis)
        results <- bind_rows(results, list(
          groups = sprintf("%s vs %s", group.1, group.2),
          cluster = cluster,
          responder = is.responder,
          group_1_size = length(group.1.fracs),
          group_2_size = length(group.2.fracs),
          hypothesis = hypothesis,
          statistic = t.test.output$statistic,
          p_value = t.test.output$p.value
        ))
      }

      # make a plot for current comparison, with a bar for mean cell fraction
      # of each time point and points for individual subjects
      if (!plot.bar) {
        next
      }

      plot.data <- fracs %>%
        filter(responder == is.responder, group %in% c(group.1, group.2)) %>%
        select(pub_id, group, {{ cluster }}) %>%
        dplyr::rename(cell_fraction = {{ cluster }}, time.point = group)
      if (is.responder) {
        plot.title <- sprintf("Responders, %s, %s vs %s", cluster, group.1,
                              group.2)
      } else {
        plot.title <- sprintf("Non-responders, %s, %s vs %s", cluster, group.1,
                              group.2)
      }
      p <- ggplot(plot.data, aes(x = time.point, y = cell_fraction,
                  fill = time.point)) +
        geom_bar(stat = "summary", fun = "mean", width = 0.5) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                                   dodge.width = 0.9),
                   size = 2) +
        labs(x = "Time point", y = "Cell fraction", fill = "Time point",
             title = plot.title) +
        theme_minimal() +
        theme(legend.position = "right")
      print(p)
    }
  }

  results
}

# analyze naive CD4+ T cells, naive CD8+ T cells, MDSCs and macrophages (C1D1/WK8 vs BL, responders or non-responders only) ####
cat("Analyzing naive CD4+ T cells, naive CD8+ T cells, MDSCs and macrophages by response status...\n")
test.subclusters <- c("naive CD4+ T cell", "naive CD8+ T cell", "effector CD8+ T cell")
test.clusters <- c("MDSC", "macrophage")

t.test.group.results <- init.t.test.group.results()
pdf(file.path(output.dir, "cell_fractions_time_points.pdf"), width = 6,
    height = 4.5)

for (time.point in c("C1D1", "WK8")) {
  t.test.group.results <- rbind(
    t.test.group.results,
    compare.groups(subcluster.fracs, test.subclusters, "BL", time.point,
                   plot.bar = TRUE)
  )

  t.test.group.results <- rbind(
    t.test.group.results,
    compare.groups(metacluster.fracs, test.clusters, "BL", time.point,
                   plot.bar = TRUE)
  )
}

dev.off()
write.csv(t.test.group.results, file.path(output.dir, "t_test_time_points.csv"),
          row.names = FALSE)

# define function to initialize an empty data.frame for t-test results between responders and non-responders ####
init.t.test.response.results <- function() {
  data.frame(
    cluster = character(),
    group = character(),
    responder_size = integer(),
    non_responder_size = integer(),
    hypothesis = character(),
    statistic = numeric(),
    p_value = numeric()
  )
}

# define function to analyze cell fractions between responders and non-responders ####
# (at specific time point for each cluster)
compare.responses <- function(fracs, clusters, groups, plot.bar = FALSE) {
  results <- init.t.test.response.results()
  t.hypotheses <- c("two.sided", "less", "greater")
  set.seed(rng.seed)

  for (cluster in clusters) {
    for (group in groups) {
      responder.fracs <- fracs %>%
        filter(responder == TRUE, group == !!group) %>%
        pull({{ cluster }})
      non.responder.fracs <- fracs %>%
        filter(responder == FALSE, group == !!group) %>%
        pull({{ cluster }})

      if (length(responder.fracs) < 2 || length(non.responder.fracs) < 2) {
        next
      }

      for (hypothesis in t.hypotheses) {
        t.test.output <- t.test(responder.fracs, non.responder.fracs,
                                alternative = hypothesis)
        results <- rbind(results, data.frame(
          cluster = cluster,
          group = group,
          responder_size = length(responder.fracs),
          non_responder_size = length(non.responder.fracs),
          hypothesis = hypothesis,
          statistic = t.test.output$statistic,
          p_value = t.test.output$p.value
        ))
      }

      # make a plot for current comparison, with a bar for mean cell fraction
      # for each response group and points for individual subjects
      if (!plot.bar) {
        next
      }

      plot.data <- fracs %>%
        filter(group == !!group, !is.na(responder)) %>%
        select(pub_id, responder, {{ cluster }}) %>%
        dplyr::rename(cell_fraction = {{ cluster }}) %>%
        mutate(responder = factor(ifelse(responder, "responder", "non-responder"),
                                  levels = c("responder", "non-responder")))
      plot.title <- sprintf("%s, %s, responders vs non-responders", cluster,
                            group)
      p <- ggplot(plot.data, aes(x = responder, y = cell_fraction,
                  fill = responder)) +
        geom_bar(stat = "summary", fun = "mean", width = 0.5) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                                   dodge.width = 0.9),
                   size = 2) +
        labs(x = "Response group", y = "Cell fraction", fill = "Response group",
             title = plot.title) +
        theme_minimal() +
        theme(legend.position = "right")
      print(p)
    }
  }

  results
}

# analyze macrophages and NK cells for responders vs non-responders at BL ####
cat("Analyzing macrophages and NK cells by response status at BL...\n")
test.clusters <- c("macrophage", "NK cell")
test.groups <- c("BL")

pdf(file.path(output.dir, "cell_fractions_response.pdf"), width = 6,
    height = 4.5)
t.test.response.results <- compare.responses(metacluster.fracs, test.clusters,
                                             test.groups, plot.bar = TRUE)
dev.off()
write.csv(t.test.response.results, file.path(output.dir, "t_test_response.csv"),
          row.names = FALSE)
