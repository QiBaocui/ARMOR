#!/usr/bin/env Rscript

# ==============================================================================
# Script: 08_compare_models.R
# Description: Compare model performance across antibiotics and feature types
#              Generate summary tables and statistical comparisons
# Author: ARMOR Pipeline
# Date: 2026-04-10
# ==============================================================================

library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------------------------
# Command line options
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Directory containing results from multiple antibiotics [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory [required]"),
  make_option(c("--pattern"), type = "character", default = "*_cutoffpermance.tsv",
              help = "File pattern to match [default: '*_cutoffpermance.tsv']"),
  make_option(c("--output-format"), type = "character", default = "png,pdf,tsv",
              help = "Output formats (comma-separated) [default: png,pdf,tsv]")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Compare model performance across conditions")
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Both --input and --outdir are required!")
}

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# Parse output formats
output_formats <- strsplit(opt$`output-format`, ",")[[1]]
output_formats <- trimws(output_formats)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

#' Collect performance data from multiple result directories
#' @param base_dir Base directory containing subdirectories for each antibiotic
#' @return Combined data frame
collect_performance_data <- function(base_dir) {
  all_data <- list()

  # Get subdirectories (each represents an antibiotic condition)
  subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

  message(sprintf("Found %d antibiotic conditions: %s",
                  length(subdirs), paste(subdirs, collapse = ", ")))

  for (antibiotic in subdirs) {
    antibiotic_dir <- file.path(base_dir, antibiotic)

    # Find all cutoff performance files
    perf_files <- list.files(antibiotic_dir, pattern = "_cutoffpermance.tsv$",
                             full.names = TRUE)

    for (perf_file in perf_files) {
      base <- basename(perf_file)
      # Parse: {Feature}_{Algorithm}_cutoffpermance.tsv
      parts <- strsplit(gsub("_cutoffpermance.tsv$", "", base), "_")[[1]]

      if (length(parts) >= 2) {
        feature_type <- parts[1]
        algorithm <- parts[2]

        data <- fread(perf_file, sep = "\t", header = TRUE)

        if (nrow(data) > 0) {
          data$Antibiotic <- antibiotic
          data$FeatureType <- feature_type
          data$Algorithm <- algorithm
          all_data[[length(all_data) + 1]] <- data
        }
      }
    }
  }

  if (length(all_data) == 0) {
    stop("No performance data found!")
  }

  return(bind_rows(all_data))
}

#' Generate comparison summary table
#' @param perf_df Combined performance data frame
#' @return Summary data frame
generate_summary_table <- function(perf_df) {
  summary_df <- perf_df %>%
    group_by(Antibiotic, FeatureType, Algorithm) %>%
    summarise(
      AUC = mean(auc, na.rm = TRUE),
      AUC_SD = sd(auc, na.rm = TRUE),
      Sensitivity = mean(sensitivity, na.rm = TRUE),
      Specificity = mean(specificity, na.rm = TRUE),
      Accuracy = mean(accuracy, na.rm = TRUE),
      F1 = mean(f1_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Antibiotic, FeatureType, desc(AUC))

  return(summary_df)
}

#' Generate comparison plots
#' @param perf_df Combined performance data frame
#' @param outdir Output directory
generate_comparison_plots <- function(perf_df, outdir) {
  # ---------------------------------------------------------------------------
  # Plot 1: AUC comparison across antibiotics (grouped bar plot)
  # ---------------------------------------------------------------------------
  p1 <- ggplot(perf_df, aes(x = Antibiotic, y = auc, fill = FeatureType)) +
    geom_bar(stat = "summary", fun = "mean", position = position.dodge(width = 0.7)) +
    geom_errorbar(stat = "summary", fun.data = "mean_se",
                  position = position.dodge(width = 0.7), width = 0.2) +
    facet_wrap(~Algorithm, ncol = 2) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    ylab("AUC (mean ± SE)") +
    xlab("") +
    ggtitle("Model Performance Across Antibiotics") +
    ylim(0.5, 1.0)

  ggsave(file.path(outdir, "auc_comparison_across_antibiotics.png"),
         p1, width = 12, height = 8, dpi = 150)

  if ("pdf" %in% output_formats) {
    ggsave(file.path(outdir, "auc_comparison_across_antibiotics.pdf"),
           p1, width = 12, height = 8)
  }

  message("  Created: AUC comparison across antibiotics")

  # ---------------------------------------------------------------------------
  # Plot 2: Feature type comparison (violin + box plot)
  # ---------------------------------------------------------------------------
  p2 <- ggplot(perf_df, aes(x = FeatureType, y = auc, fill = FeatureType)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.15, alpha = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      legend.position = "none"
    ) +
    xlab("Feature Type") +
    ylab("AUC") +
    ggtitle("AUC Distribution by Feature Type")

  ggsave(file.path(outdir, "auc_by_feature_type.png"),
         p2, width = 8, height = 6, dpi = 150)

  if ("pdf" %in% output_formats) {
    ggsave(file.path(outdir, "auc_by_feature_type.pdf"),
           p2, width = 8, height = 6)
  }

  message("  Created: AUC by feature type")

  # ---------------------------------------------------------------------------
  # Plot 3: Algorithm comparison
  # ---------------------------------------------------------------------------
  p3 <- ggplot(perf_df, aes(x = Algorithm, y = auc, fill = Algorithm)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, alpha = 0.3, size = 2) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none"
    ) +
    xlab("Algorithm") +
    ylab("AUC") +
    ggtitle("AUC Distribution by Algorithm")

  ggsave(file.path(outdir, "auc_by_algorithm.png"),
         p3, width = 8, height = 6, dpi = 150)

  if ("pdf" %in% output_formats) {
    ggsave(file.path(outdir, "auc_by_algorithm.pdf"),
           p3, width = 8, height = 6)
  }

  message("  Created: AUC by algorithm")

  # ---------------------------------------------------------------------------
  # Plot 4: Heatmap of best AUC per antibiotic
  # ---------------------------------------------------------------------------
  best_auc <- perf_df %>%
    group_by(Antibiotic, FeatureType, Algorithm) %>%
    summarise(AUC = mean(auc, na.rm = TRUE), .groups = "drop") %>%
    group_by(Antibiotic, FeatureType) %>%
    slice_max(order_by = AUC, n = 1) %>%
    ungroup()

  p4 <- ggplot(best_auc, aes(x = FeatureType, y = Antibiotic, fill = AUC)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "#fee0d2", high = "#d7301f", limits = c(0.6, 1.0)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank()
    ) +
    xlab("Feature Type") +
    ylab("Antibiotic") +
    ggtitle("Best Model AUC per Antibiotic")

  ggsave(file.path(outdir, "best_auc_heatmap.png"),
         p4, width = 8, height = 6, dpi = 150)

  if ("pdf" %in% output_formats) {
    ggsave(file.path(outdir, "best_auc_heatmap.pdf"),
           p4, width = 8, height = 6)
  }

  message("  Created: Best AUC heatmap")
}

#' Statistical comparison
#' @param perf_df Combined performance data frame
#' @return Statistical test results
perform_statistical_comparison <- function(perf_df) {
  # Compare feature types using ANOVA
  message("\n=== Statistical Comparison ===")

  # ANOVA for feature types
  anova_formula <- auc ~ FeatureType
  anova_result <- aov(anova_formula, data = perf_df)
  anova_summary <- summary(anova_result)

  message("ANOVA: AUC ~ Feature Type")
  message(sprintf("  F-value: %.2f", anova_summary[[1]]$`F value`[1]))
  message(sprintf("  P-value: %.2e", anova_summary[[1]]$`Pr(>F)`[1]))

  # Pairwise t-tests
  pairwise_result <- pairwise.t.test(perf_df$auc, perf_df$FeatureType,
                                     p.adjust.method = "BH")

  stats_results <- list(
    anova = anova_summary,
    pairwise = pairwise_result
  )

  return(stats_results)
}

# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

main <- function() {
  message("=== Model Comparison Pipeline ===")
  message(sprintf("Input directory: %s", opt$input))
  message(sprintf("Output directory: %s", opt$outdir))

  # Collect data
  message("\n=== Collecting Performance Data ===")
  perf_df <- collect_performance_data(opt$input)

  message(sprintf("Total records: %d", nrow(perf_df)))
  message(sprintf("Antibiotics: %d", length(unique(perf_df$Antibiotic))))
  message(sprintf("Feature types: %s",
                  paste(unique(perf_df$FeatureType), collapse = ", ")))
  message(sprintf("Algorithms: %s",
                  paste(unique(perf_df$Algorithm), collapse = ", ")))

  # Generate summary table
  message("\n=== Generating Summary Table ===")
  summary_df <- generate_summary_table(perf_df)

  summary_file <- file.path(opt$outdir, "performance_summary.tsv")
  fwrite(summary_df, summary_file, sep = "\t")
  message(sprintf("Summary saved: %s", summary_file))

  # Generate plots
  message("\n=== Generating Comparison Plots ===")
  generate_comparison_plots(perf_df, opt$outdir)

  # Statistical comparison
  message("\n=== Statistical Analysis ===")
  stats_results <- perform_statistical_comparison(perf_df)

  # Save statistical results
  stats_file <- file.path(opt$outdir, "statistical_comparison.txt")
  sink(stats_file)
  cat("=== Model Comparison Statistical Analysis ===\n\n")
  cat("ANOVA: AUC ~ Feature Type\n")
  print(stats_results$anova)
  cat("\n\nPairwise Comparisons (BH adjusted p-values):\n")
  print(stats_results$pairwise)
  sink()
  message(sprintf("Statistical results saved: %s", stats_file))

  message("\n=== Comparison Complete ===")
}

main()
