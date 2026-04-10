#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_model_evaluation.R
# Description: Comprehensive model evaluation and comparison
#              Calculate and compare AUC, sensitivity, specificity, accuracy, F1
# Author: ARMOR Pipeline
# Date: 2026-04-10
# ==============================================================================

library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pROC)

# ------------------------------------------------------------------------------
# Command line options
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Directory containing model results [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory for evaluation results [required]"),
  make_option(c("-p", "--pattern"), type = "character", default = "*_cutoffpermance.tsv",
              help = "File pattern for performance metrics [default: '*_cutoffpermance.tsv']"),
  make_option(c("--feature-types"), type = "character", default = "Gene,SNP,Kmer,All",
              help = "Feature types to analyze [default: Gene,SNP,Kmer,All]"),
  make_option(c("--algorithms"), type = "character", default = "GBM,GLM,RF,DL,stacked",
              help = "Algorithms to analyze [default: GBM,GLM,RF,DL,stacked]")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Model evaluation and comparison")
opt <- parse_args(opt_parser)

# Validate arguments
if (is.null(opt$input) || is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Both --input and --outdir are required!")
}

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

#' Calculate performance metrics from confusion matrix values
#' @param tp True positives
#' @param fp False positives
#' @param tn True negatives
#' @param fn False negatives
#' @return Named vector with metrics
calculate_metrics <- function(tp, fp, tn, fn) {
  # Sensitivity (Recall, TPR)
  sensitivity <- tp / (tp + fn)

  # Specificity (TNR)
  specificity <- tn / (tn + fp)

  # Accuracy
  accuracy <- (tp + tn) / (tp + tn + fp + fn)

  # Precision (PPV)
  precision <- tp / (tp + fp)

  # F1 Score
  f1 <- 2 * (precision * sensitivity) / (precision + sensitivity)

  # Balanced Accuracy
  bal_accuracy <- (sensitivity + specificity) / 2

  return(c(
    sensitivity = sensitivity,
    specificity = specificity,
    accuracy = accuracy,
    precision = precision,
    f1_score = f1,
    balanced_accuracy = bal_accuracy
  ))
}

#' Parse performance file
#' @param file_path Path to cutoff performance TSV
#' @return Data frame with metrics
parse_performance <- function(file_path) {
  if (!file.exists(file_path)) {
    return(NULL)
  }

  perf <- fread(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  # Extract metric names from column names
  # Expected columns: tpr, fpr, threshold, etc.

  # Calculate metrics if raw values available
  if ("tpr" %in% colnames(perf) && "fpr" %in% colnames(perf)) {
    perf$sensitivity <- perf$tpr
    perf$specificity <- 1 - perf$fpr
  }

  return(perf)
}

#' Extract model info from filename
#' @param filename Performance file name
#' @return List with feature type and algorithm
parse_filename <- function(filename) {
  # Expected format: {FeatureType}_{Algorithm}_cutoffpermance.tsv
  base <- basename(filename)
  base <- gsub("_cutoffpermance.tsv$", "", base)

  parts <- strsplit(base, "_")[[1]]

  if (length(parts) >= 2) {
    feature_type <- parts[1]
    algorithm <- parts[2]
  } else {
    feature_type <- "Unknown"
    algorithm <- "Unknown"
  }

  return(list(feature = feature_type, algorithm = algorithm))
}

# ------------------------------------------------------------------------------
# Main analysis functions
# ------------------------------------------------------------------------------

#' Collect all performance metrics
#' @param input_dir Directory with results
#' @param feature_types Feature types to include
#' @param algorithms Algorithms to include
#' @return Combined data frame
collect_all_metrics <- function(input_dir, feature_types, algorithms) {
  all_metrics <- list()

  for (feature in feature_types) {
    for (algo in algorithms) {
      perf_file <- file.path(input_dir,
                             sprintf("%s_%s_cutoffpermance.tsv", feature, algo))

      if (file.exists(perf_file)) {
        perf <- parse_performance(perf_file)

        if (!is.null(perf) && nrow(perf) > 0) {
          perf$FeatureType <- feature
          perf$Algorithm <- algo
          all_metrics[[length(all_metrics) + 1]] <- perf

          message(sprintf("Loaded: %s_%s (AUC: %.4f)",
                         feature, algo,
                         ifelse("auc" %in% colnames(perf), perf$auc[1], NA)))
        }
      }
    }
  }

  if (length(all_metrics) == 0) {
    stop("No performance files found!")
  }

  return(bind_rows(all_metrics))
}

#' Create comparison plots
#' @param metrics_df Combined metrics data frame
#' @param outdir Output directory
#' @return List of ggplot objects
create_comparison_plots <- function(metrics_df, outdir) {
  plots <- list()

  # ---------------------------------------------------------------------------
  # Plot 1: AUC comparison by feature type and algorithm
  # ---------------------------------------------------------------------------
  if ("auc" %in% colnames(metrics_df)) {
    p1 <- ggplot(metrics_df, aes(x = FeatureType, y = auc, fill = Algorithm)) +
      geom_bar(stat = "identity", position = position.dodge(width = 0.8), width = 0.7) +
      geom_errorbar(aes(ymin = auc - sd(auc, na.rm = TRUE) / sqrt(n()),
                        ymax = auc + sd(auc, na.rm = TRUE) / sqrt(n())),
                    position = position.dodge(width = 0.8), width = 0.2) +
      scale_y_continuous(limits = c(0.5, 1.0), expand = c(0, 0)) +
      scale_fill_brewer(palette = "Set2") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        legend.title = element_blank()
      ) +
      xlab("Feature Type") +
      ylab("AUC") +
      ggtitle("Model Performance Comparison: AUC by Feature Type")

    plots[["auc_comparison"]] <- p1
    ggsave(file.path(outdir, "auc_comparison.png"), p1, width = 10, height = 6, dpi = 150)
    ggsave(file.path(outdir, "auc_comparison.pdf"), p1, width = 10, height = 6)
  }

  # ---------------------------------------------------------------------------
  # Plot 2: Heatmap of all metrics
  # ---------------------------------------------------------------------------
  # Reshape for heatmap
  metrics_long <- metrics_df %>%
    select(FeatureType, Algorithm, auc, sensitivity, specificity, f1_score) %>%
    pivot_longer(cols = c(auc, sensitivity, specificity, f1_score),
                 names_to = "Metric", values_to = "Value")

  p2 <- ggplot(metrics_long, aes(x = FeatureType, y = Algorithm, fill = Value)) +
    geom_tile(color = "white") +
    facet_grid(~Metric) +
    scale_fill_gradient2(low = "#fee0d2", mid = "#fc8d59", high = "#d7301f",
                         midpoint = 0.75, limits = c(0, 1)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank()
    ) +
    xlab("Feature Type") +
    ylab("Algorithm") +
    ggtitle("Performance Heatmap")

  plots[["heatmap"]] <- p2
  ggsave(file.path(outdir, "performance_heatmap.png"), p2, width = 10, height = 8, dpi = 150)
  ggsave(file.path(outdir, "performance_heatmap.pdf"), p2, width = 10, height = 8)

  # ---------------------------------------------------------------------------
  # Plot 3: Radar chart for best models
  # ---------------------------------------------------------------------------
  # Select best model from each feature type
  best_models <- metrics_df %>%
    group_by(FeatureType) %>%
    slice_max(order_by = auc, n = 1) %>%
    ungroup()

  plots[["best_models"]] <- NULL  # Placeholder for radar chart

  return(plots)
}

#' Generate summary statistics
#' @param metrics_df Combined metrics data frame
#' @return Summary data frame
generate_summary <- function(metrics_df) {
  summary_stats <- metrics_df %>%
    group_by(FeatureType, Algorithm) %>%
    summarise(
      AUC_mean = mean(auc, na.rm = TRUE),
      AUC_sd = sd(auc, na.rm = TRUE),
      Sensitivity_mean = mean(sensitivity, na.rm = TRUE),
      Specificity_mean = mean(specificity, na.rm = TRUE),
      F1_mean = mean(f1_score, na.rm = TRUE),
      Accuracy_mean = mean(accuracy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(AUC_mean))

  return(summary_stats)
}

# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

main <- function() {
  # Create output directory
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }

  # Parse feature types and algorithms
  feature_types <- strsplit(opt$`feature-types`, ",")[[1]]
  feature_types <- trimws(feature_types)

  algorithms <- strsplit(opt$algorithms, ",")[[1]]
  algorithms <- trimws(algorithms)

  message("=== Collecting Performance Metrics ===")
  message(sprintf("Feature types: %s", paste(feature_types, collapse = ", ")))
  message(sprintf("Algorithms: %s", paste(algorithms, collapse = ", ")))

  # Collect all metrics
  metrics_df <- collect_all_metrics(opt$input, feature_types, algorithms)

  message(sprintf("\nTotal performance records: %d", nrow(metrics_df)))

  # Generate summary statistics
  message("\n=== Generating Summary Statistics ===")
  summary_stats <- generate_summary(metrics_df)

  # Save summary
  summary_file <- file.path(opt$outdir, "model_comparison_summary.tsv")
  fwrite(summary_stats, summary_file, sep = "\t")
  message(sprintf("Summary saved: %s", summary_file))

  # Create visualizations
  message("\n=== Creating Visualizations ===")
  plots <- create_comparison_plots(metrics_df, opt$outdir)

  for (plot_name in names(plots)) {
    if (!is.null(plots[[plot_name]])) {
      message(sprintf("  Created: %s", plot_name))
    }
  }

  # Save raw data
  raw_file <- file.path(opt$outdir, "all_metrics_raw.tsv")
  fwrite(metrics_df, raw_file, sep = "\t")
  message(sprintf("Raw data saved: %s", raw_file))

  message("\n=== Model Evaluation Complete ===")
}

main()
