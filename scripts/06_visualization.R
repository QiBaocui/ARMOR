#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_visualization.R
# Description: Generate publication-quality figures for AMR prediction results
#              Based on manuscript Figure 3, 4, 5
# Author: ARMOR Pipeline
# Date: 2026-04-10
# ==============================================================================

library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(pROC)

# ------------------------------------------------------------------------------
# Command line options
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Directory containing model results [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory for figures [required]"),
  make_option(c("--dpi"), type = "integer", default = 300,
              help = "Resolution for PNG output [default: 300]"),
  make_option(c("--theme"), type = "character", default = "bw",
              help = "Plot theme: bw/minimal/classic [default: bw]"),
  make_option(c("--prefix"), type = "character", default = "",
              help = "Prefix for output files")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Generate publication-quality figures")
opt <- parse_args(opt_parser)

# Validate arguments
if (is.null(opt$input) || is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Both --input and --outdir are required!")
}

# Create output directory
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# Set theme
if (opt$theme == "bw") {
  plot_theme <- theme_bw()
} else if (opt$theme == "minimal") {
  plot_theme <- theme_minimal()
} else {
  plot_theme <- theme_classic()
}

# ------------------------------------------------------------------------------
# Figure 1: ROC Curve Comparison
# Corresponds to manuscript Figure 3
# ------------------------------------------------------------------------------

#' Generate ROC curve comparison plot
#' @param input_dir Directory with ROC data files
#' @param output_file Output file path
#' @return ggplot object
generate_roc_comparison <- function(input_dir, output_file) {
  message("Generating ROC curve comparison...")

  # Collect ROC data from all models
  roc_data_list <- list()

  roc_files <- list.files(input_dir, pattern = "_rocData.tsv$", full.names = TRUE)

  for (roc_file in roc_files) {
    base <- basename(roc_file)
    # Parse: {Feature}_{Algorithm}_rocData.tsv
    parts <- strsplit(gsub("_rocData.tsv$", "", base), "_")[[1]]

    if (length(parts) >= 2) {
      feature <- parts[1]
      algorithm <- parts[2]

      data <- fread(roc_file, sep = "\t", header = TRUE)

      if ("fpr" %in% colnames(data) && "tpr" %in% colnames(data)) {
        data$Feature <- feature
        data$Algorithm <- algorithm
        data$Label <- paste0(feature, "_", algorithm)
        roc_data_list[[length(roc_data_list) + 1]] <- data
      }
    }
  }

  if (length(roc_data_list) == 0) {
    message("No ROC data found!")
    return(NULL)
  }

  roc_data <- bind_rows(roc_data_list)

  # Get AUC values for labels
  auc_labels <- roc_data %>%
    group_by(Feature, Algorithm, Label) %>%
    summarise(AUC = max(tpr, na.rm = TRUE), .groups = "drop") %>%
    mutate(Label = paste0(Algorithm, " (AUC = ", round(AUC, 3), ")"))

  # Plot ROC curves
  p <- ggplot(roc_data, aes(x = fpr, y = tpr, color = Algorithm)) +
    geom_line(size = 1) +
    facet_wrap(~Feature, ncol = 2) +
    scale_color_brewer(palette = "Set1") +
    plot_theme +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("False Positive Rate (1 - Specificity)") +
    ylab("True Positive Rate (Sensitivity)") +
    ggtitle("ROC Curves by Feature Type and Algorithm")

  ggsave(output_file, p, width = 10, height = 8, dpi = opt$dpi)
  message(sprintf("Saved: %s", output_file))

  return(p)
}

# ------------------------------------------------------------------------------
# Figure 2: Performance Metrics Comparison (Multi-panel)
# Corresponds to manuscript Figure 4
# ------------------------------------------------------------------------------

#' Generate multi-metric comparison plot
#' @param input_dir Directory with performance files
#' @param output_file Output file path
#' @return ggplot object
generate_performance_comparison <- function(input_dir, output_file) {
  message("Generating performance comparison...")

  # Collect all performance metrics
  perf_data <- list()

  perf_files <- list.files(input_dir, pattern = "_cutoffpermance.tsv$", full.names = TRUE)

  for (perf_file in perf_files) {
    base <- basename(perf_file)
    parts <- strsplit(gsub("_cutoffpermance.tsv$", "", base), "_")[[1]]

    if (length(parts) >= 2) {
      feature <- parts[1]
      algorithm <- parts[2]

      data <- fread(perf_file, sep = "\t", header = TRUE)

      if (nrow(data) > 0) {
        data$Feature <- feature
        data$Algorithm <- algorithm
        perf_data[[length(perf_data) + 1]] <- data
      }
    }
  }

  if (length(perf_data) == 0) {
    message("No performance data found!")
    return(NULL)
  }

  perf_df <- bind_rows(perf_data)

  # Calculate metrics if not present
  if ("sensitivity" %in% colnames(perf_df)) {
    perf_df$Specificity <- 1 - perf_df$fpr
  }

  # Reshape for faceted plot
  metrics_to_plot <- c("auc", "sensitivity", "specificity", "accuracy", "f1_score")
  available_metrics <- intersect(metrics_to_plot, colnames(perf_df))

  if (length(available_metrics) < 2) {
    message("Insufficient metrics for comparison plot")
    return(NULL)
  }

  plot_data <- perf_df %>%
    select(Feature, Algorithm, all_of(available_metrics)) %>%
    pivot_longer(cols = all_of(available_metrics),
                 names_to = "Metric", values_to = "Value")

  # Metric labels for display
  metric_labels <- c(
    auc = "AUC",
    sensitivity = "Sensitivity",
    specificity = "Specificity",
    accuracy = "Accuracy",
    f1_score = "F1 Score"
  )

  plot_data$Metric <- factor(plot_data$Metric, levels = available_metrics)

  # Create faceted bar plot
  p <- ggplot(plot_data, aes(x = Algorithm, y = Value, fill = Feature)) +
    geom_bar(stat = "identity", position = position.dodge(width = 0.7), width = 0.6) +
    facet_wrap(~Metric, scales = "free_y", labeller = labeller(Metric = metric_labels)) +
    scale_fill_brewer(palette = "Set2") +
    plot_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    ylab("Score") +
    xlab("") +
    ggtitle("Model Performance Metrics Comparison")

  ggsave(output_file, p, width = 12, height = 8, dpi = opt$dpi)
  message(sprintf("Saved: %s", output_file))

  return(p)
}

# ------------------------------------------------------------------------------
# Figure 3: Feature Importance Visualization
# Corresponds to manuscript Figure 5
# ------------------------------------------------------------------------------

#' Generate feature importance heatmap
#' @param input_dir Directory with variable importance files
#' @param gene_annotation_file Optional gene annotation file
#' @param output_file Output file path
#' @return ggplot object
generate_feature_importance <- function(input_dir, gene_annotation_file = NULL,
                                        output_file) {
  message("Generating feature importance heatmap...")

  # Collect variable importance from all models
  varimp_data <- list()

  varimp_files <- list.files(input_dir, pattern = "_varimp.tsv$", full.names = TRUE)

  for (varimp_file in varimp_files) {
    base <- basename(varimp_file)
    parts <- strsplit(gsub("_varimp.tsv$", "", base), "_")[[1]]

    if (length(parts) >= 2) {
      feature_type <- parts[1]
      algorithm <- parts[2]

      # Only process Gene-based models for this visualization
      if (feature_type == "Gene") {
        data <- fread(varimp_file, sep = "\t", header = TRUE)

        if (nrow(data) > 0 && "variable" %in% colnames(data)) {
          # Normalize importance scores
          if ("scaled_importance" %in% colnames(data)) {
            importance_col <- "scaled_importance"
          } else if ("importance" %in% colnames(data)) {
            importance_col <- "importance"
          } else {
            next
          }

          data$FeatureType <- feature_type
          data$Algorithm <- algorithm
          varimp_data[[length(varimp_data) + 1]] <- data[, c("variable", importance_col, "FeatureType", "Algorithm")]
          colnames(data)[1:2] <- c("Gene", "Importance")
        }
      }
    }
  }

  if (length(varimp_data) == 0) {
    message("No variable importance data found!")
    return(NULL)
  }

  varimp_df <- bind_rows(varimp_data)

  # Get top genes across all models
  top_genes <- varimp_df %>%
    group_by(variable) %>%
    summarise(mean_importance = mean(get(colnames(varimp_df)[2]), na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(mean_importance)) %>%
    head(40) %>%
    pull(variable)

  # Filter data for top genes
  plot_data <- varimp_df %>%
    filter(variable %in% top_genes)

  # Calculate mean importance per gene per antibiotic (Algorithm as proxy)
  gene_heatmap <- plot_data %>%
    group_by(variable, Algorithm) %>%
    summarise(Importance = mean(get(colnames(plot_data)[2]), na.rm = TRUE),
              .groups = "drop")

  # Create heatmap
  p <- ggplot(gene_heatmap, aes(x = Algorithm, y = variable, fill = Importance)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "steelblue", na.value = "gray90") +
    plot_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    xlab("Algorithm") +
    ylab("Gene") +
    ggtitle("Top 40 Resistance-Associated Genes")

  ggsave(output_file, p, width = 10, height = 12, dpi = opt$dpi)
  message(sprintf("Saved: %s", output_file))

  return(p)
}

# ------------------------------------------------------------------------------
# Additional visualizations
# ------------------------------------------------------------------------------

#' Generate sample size distribution plot
generate_sample_distribution <- function(group_file, output_file) {
  if (!file.exists(group_file)) {
    message("Group file not found, skipping sample distribution plot")
    return(NULL)
  }

  group_data <- read.table(group_file, sep = "\t", header = TRUE)

  # Count by group
  dist_data <- group_data %>%
    count(Group) %>%
    mutate(Percentage = n / sum(n) * 100)

  p <- ggplot(dist_data, aes(x = Group, y = n, fill = Group)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = paste0(n, "\n(", round(Percentage, 1), "%)")),
              vjust = -0.5, size = 4) +
    plot_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    xlab("Resistance Phenotype") +
    ylab("Number of Isolates") +
    ggtitle("Sample Distribution by Resistance Phenotype")

  ggsave(output_file, p, width = 6, height = 5, dpi = opt$dpi)
  message(sprintf("Saved: %s", output_file))

  return(p)
}

# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

main <- function() {
  message("=== ARMOR Visualization Pipeline ===")
  message(sprintf("Input directory: %s", opt$input))
  message(sprintf("Output directory: %s", opt$outdir))

  # Figure 1: ROC Curves
  message("\n--- Generating Figure 1: ROC Curves ---")
  generate_roc_comparison(
    opt$input,
    file.path(opt$outdir, paste0(opt$prefix, "figure1_roc_curves.png"))
  )

  # Figure 2: Performance Comparison
  message("\n--- Generating Figure 2: Performance Comparison ---")
  generate_performance_comparison(
    opt$input,
    file.path(opt$outdir, paste0(opt$prefix, "figure2_performance_comparison.png"))
  )

  # Figure 3: Feature Importance
  message("\n--- Generating Figure 3: Feature Importance ---")
  generate_feature_importance(
    opt$input,
    NULL,
    file.path(opt$outdir, paste0(opt$prefix, "figure3_feature_importance.png"))
  )

  message("\n=== Visualization Complete ===")
}

main()
