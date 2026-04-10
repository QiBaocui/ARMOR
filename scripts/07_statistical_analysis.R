#!/usr/bin/env Rscript

# ==============================================================================
# Script: 07_statistical_analysis.R
# Description: Statistical analysis of features and model comparisons
#              Chi-square tests, Odds Ratios, and significance testing
# Author: ARMOR Pipeline
# Date: 2026-04-10
# ==============================================================================

library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)
library(broom)

# ------------------------------------------------------------------------------
# Command line options
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Feature matrix file [required]"),
  make_option(c("-g", "--group"), type = "character", default = NULL,
              help = "Group/phenotype file [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory [required]"),
  make_option(c("--p-adjust"), type = "character", default = "BH",
              help = "P-value adjustment method: BH/Bonferroni/Holm [default: BH]"),
  make_option(c("--or-threshold"), type = "double", default = 2.0,
              help = "Odds ratio threshold for significant associations [default: 2.0]"),
  make_option(c("--top-n"), type = "integer", default = 50,
              help = "Number of top features to report [default: 50]")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Statistical analysis of feature associations")
opt <- parse_args(opt_parser)

# Validate arguments
if (is.null(opt$input) || is.null(opt$group) || is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Arguments --input, --group, and --outdir are required!")
}

# Create output directory
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Statistical analysis functions
# ------------------------------------------------------------------------------

#' Perform chi-square test and calculate odds ratio
#' @param feature_vector Binary feature vector
#' @param group_vector Group labels
#' @return List with test statistics
perform_statistical_test <- function(feature_vector, group_vector) {
  # Create contingency table
  # Convert to binary
  feature_binary <- ifelse(feature_vector > 0, 1, 0)

  # Handle missing values
  valid_idx <- !is.na(feature_binary) & !is.na(group_vector)
  if (sum(valid_idx) < 10) {
    return(list(
      p_value = NA,
      odds_ratio = NA,
      ci_low = NA,
      ci_high = NA,
      n = sum(valid_idx)
    ))
  }

  feature_clean <- feature_binary[valid_idx]
  group_clean <- group_vector[valid_idx]

  contingency <- table(feature_clean, group_clean)

  # Ensure 2x2 table
  if (nrow(contingency) < 2 || ncol(contingency) < 2) {
    return(list(
      p_value = NA,
      odds_ratio = NA,
      ci_low = NA,
      ci_high = NA,
      n = sum(valid_idx)
    ))
  }

  # Chi-square test
  chi_result <- tryCatch({
    chisq.test(contingency, correct = TRUE)
  }, warning = function(w) {
    # Use Fisher's exact test if expected counts are low
    fisher.test(contingency)
  }, error = function(e) {
    NULL
  })

  if (is.null(chi_result)) {
    return(list(
      p_value = NA,
      odds_ratio = NA,
      ci_low = NA,
      ci_high = NA,
      n = sum(valid_idx)
    ))
  }

  # Calculate odds ratio
  # Table structure:
  #           Group1  Group2
  # feat=1      a       b
  # feat=0      c       d
  a <- as.numeric(contingency[2, 2])
  b <- as.numeric(contingency[2, 1])
  c <- as.numeric(contingency[1, 2])
  d <- as.numeric(contingency[1, 1])

  # Add pseudo-counts if needed
  if (b == 0 || c == 0) {
    a <- a + 0.5
    b <- b + 0.5
    c <- c + 0.5
    d <- d + 0.5
  }

  odds_ratio <- (a * d) / (b * c)

  # Calculate 95% CI for OR
  se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)
  log_or <- log(odds_ratio)
  ci_low <- exp(log_or - 1.96 * se_log_or)
  ci_high <- exp(log_or + 1.96 * se_log_or)

  return(list(
    p_value = chi_result$p.value,
    odds_ratio = odds_ratio,
    ci_low = ci_low,
    ci_high = ci_high,
    n = sum(valid_idx)
  ))
}

#' Analyze all features
#' @param feature_matrix Data frame with features
#' @param group_data Group labels
#' @param p_adjust_method P-value adjustment method
#' @return Data frame with statistics
analyze_all_features <- function(feature_matrix, group_data,
                                  p_adjust_method = "BH") {
  # Prepare data
  if ("ID" %in% colnames(feature_matrix)) {
    rownames(feature_matrix) <- feature_matrix$ID
    feature_matrix <- feature_matrix[, -which(colnames(feature_matrix) == "ID")]
  }

  if ("ID" %in% colnames(group_data)) {
    rownames(group_data) <- group_data$ID
    group_data <- group_data[, -which(colnames(group_data) == "ID")]
  }

  # Get group vector
  if ("Group" %in% colnames(group_data)) {
    group_vector <- group_data$Group
  } else {
    group_vector <- group_data[, 1]
  }

  # Match samples
  common_samples <- intersect(rownames(feature_matrix), names(group_vector))
  message(sprintf("Analyzing %d features in %d samples...",
                  ncol(feature_matrix), length(common_samples)))

  feature_matrix <- feature_matrix[common_samples, , drop = FALSE]
  group_vector <- group_vector[common_samples]

  # Results storage
  results <- data.frame(
    Feature = colnames(feature_matrix),
    P_value = numeric(ncol(feature_matrix)),
    P_adjusted = numeric(ncol(feature_matrix)),
    Odds_Ratio = numeric(ncol(feature_matrix)),
    CI_Low = numeric(ncol(feature_matrix)),
    CI_High = numeric(ncol(feature_matrix)),
    Significant = logical(ncol(feature_matrix)),
    stringsAsFactors = FALSE
  )

  # Analyze each feature
  for (i in seq_along(colnames(feature_matrix))) {
    if (i %% 500 == 0) {
      message(sprintf("  Processing feature %d/%d...", i, ncol(feature_matrix)))
    }

    test_result <- perform_statistical_test(
      feature_vector = feature_matrix[, i],
      group_vector = group_vector
    )

    results$P_value[i] <- test_result$p_value
    results$Odds_Ratio[i] <- test_result$odds_ratio
    results$CI_Low[i] <- test_result$ci_low
    results$CI_High[i] <- test_result$ci_high
  }

  # Adjust p-values
  results$P_adjusted <- p.adjust(results$P_value, method = p_adjust_method)

  # Define significance
  results$Significant <- results$P_adjusted < 0.05 &
                         !is.na(results$Odds_Ratio) &
                         results$Odds_Ratio > 2

  return(results)
}

# ------------------------------------------------------------------------------
# Visualization functions
# ------------------------------------------------------------------------------

#' Generate volcano plot
#' @param results Statistical results data frame
#' @param output_file Output file path
generate_volcano_plot <- function(results, output_file) {
  # Prepare data for volcano plot
  plot_data <- results %>%
    mutate(
      log10_p = -log10(P_adjusted + 1e-10),
      log2_or = log2(Odds_Ratio + 1e-10),
      Significant = Significant
    )

  # Create volcano plot
  p <- ggplot(plot_data, aes(x = log2_or, y = log10_p, color = Significant)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("gray", "red")) +
    theme_bw() +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    ) +
    xlab("log2(Odds Ratio)") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle("Volcano Plot: Feature Associations with Resistance") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = log2(2), linetype = "dashed", color = "blue")

  ggsave(output_file, p, width = 8, height = 6, dpi = 150)
  message(sprintf("Volcano plot saved: %s", output_file))

  return(p)
}

#' Generate Manhattan-style plot
#' @param results Statistical results data frame
#' @param output_file Output file path
generate_manhattan_plot <- function(results, output_file) {
  plot_data <- results %>%
    mutate(
      Feature_Index = 1:n(),
      log10_p = -log10(P_adjusted + 1e-10)
    )

  p <- ggplot(plot_data, aes(x = Feature_Index, y = log10_p)) +
    geom_point(alpha = 0.5, color = "steelblue", size = 1) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank()
    ) +
    xlab("Feature Index") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle("Manhattan Plot: Feature Association Statistics") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

  ggsave(output_file, p, width = 10, height = 5, dpi = 150)
  message(sprintf("Manhattan plot saved: %s", output_file))

  return(p)
}

#' Generate top features bar plot
#' @param results Statistical results data frame
#' @param n Number of top features to show
#' @param output_file Output file path
generate_top_features_plot <- function(results, n = 50, output_file) {
  top_features <- results %>%
    filter(Significant) %>%
    arrange(desc(Odds_Ratio)) %>%
    head(n)

  if (nrow(top_features) == 0) {
    message("No significant features found for plotting")
    return(NULL)
  }

  p <- ggplot(top_features, aes(x = reorder(Feature, Odds_Ratio), y = Odds_Ratio)) +
    geom_bar(stat = "identity", fill = "coral") +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.margin = margin(1, 2, 1, 1, "cm")
    ) +
    xlab("Feature") +
    ylab("Odds Ratio") +
    ggtitle(sprintf("Top %d Resistance-Associated Features", nrow(top_features)))

  ggsave(output_file, p, width = 10, height = 8, dpi = 150)
  message(sprintf("Top features plot saved: %s", output_file))

  return(p)
}

# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

main <- function() {
  message("=== Statistical Analysis Pipeline ===")

  # Load data
  message("Loading feature matrix...")
  feature_matrix <- fread(opt$input, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE)

  message("Loading group file...")
  group_data <- read.table(opt$group, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE)

  # Perform statistical analysis
  message("\n=== Performing Statistical Tests ===")
  results <- analyze_all_features(
    feature_matrix = feature_matrix,
    group_data = group_data,
    p_adjust_method = opt$`p-adjust`
  )

  # Save full results
  results_file <- file.path(opt$outdir, "statistical_analysis_full.tsv")
  fwrite(results, results_file, sep = "\t")
  message(sprintf("\nFull results saved: %s", results_file))

  # Summary statistics
  n_significant <- sum(results$Significant, na.rm = TRUE)
  n_total <- nrow(results)
  n_sig_or <- sum(results$Odds_Ratio > opt$`or-threshold`, na.rm = TRUE)
  n_sig_p <- sum(results$P_adjusted < 0.05, na.rm = TRUE)

  message("\n=== Summary ===")
  message(sprintf("Total features analyzed: %d", n_total))
  message(sprintf("Significant (adjusted p < 0.05): %d", n_sig_p))
  message(sprintf("OR > %.1f: %d", opt$`or-threshold`, n_sig_or))
  message(sprintf("Both criteria met: %d", n_significant))

  # Save summary
  summary_lines <- c(
    "=== Statistical Analysis Summary ===",
    sprintf("Analysis date: %s", Sys.time()),
    sprintf("Input file: %s", opt$input),
    sprintf("Group file: %s", opt$group),
    sprintf("P-value adjustment: %s", opt$`p-adjust`),
    sprintf("OR threshold: %.1f", opt$`or-threshold`),
    "",
    "=== Results ===",
    sprintf("Total features: %d", n_total),
    sprintf("Significant (adj. p < 0.05): %d (%.1f%%)",
            n_sig_p, n_sig_p / n_total * 100),
    sprintf("OR > %.1f: %d (%.1f%%)",
            opt$`or-threshold`, n_sig_or, n_sig_or / n_total * 100),
    sprintf("Both criteria: %d (%.1f%%)",
            n_significant, n_significant / n_total * 100),
    "",
    "=== Top 10 Features by Odds Ratio ==="
  )

  top_10 <- results %>%
    filter(Significant) %>%
    arrange(desc(Odds_Ratio)) %>%
    head(10)

  if (nrow(top_10) > 0) {
    for (i in 1:nrow(top_10)) {
      summary_lines <- c(summary_lines,
        sprintf("  %d. %s (OR=%.2f, p=%.2e)",
                i, top_10$Feature[i], top_10$Odds_Ratio[i], top_10$P_adjusted[i])
      )
    }
  }

  summary_file <- file.path(opt$outdir, "statistical_analysis_summary.txt")
  writeLines(summary_lines, summary_file)
  message(sprintf("Summary saved: %s", summary_file))

  # Generate visualizations
  message("\n=== Generating Visualizations ===")

  generate_volcano_plot(results, file.path(opt$outdir, "volcano_plot.png"))
  generate_manhattan_plot(results, file.path(opt$outdir, "manhattan_plot.png"))
  generate_top_features_plot(results, opt$top_n,
                             file.path(opt$outdir, "top_features_barplot.png"))

  # Save significant features list
  sig_features <- results %>%
    filter(Significant) %>%
    arrange(desc(Odds_Ratio)) %>%
    select(Feature, Odds_Ratio, CI_Low, CI_High, P_adjusted)

  sig_file <- file.path(opt$outdir, "significant_features.tsv")
  fwrite(sig_features, sig_file, sep = "\t")
  message(sprintf("Significant features saved: %s", sig_file))

  message("\n=== Statistical Analysis Complete ===")
}

main()
