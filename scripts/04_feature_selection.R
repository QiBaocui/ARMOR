#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_feature_selection.R
# Description: Feature selection using between-group difference strategy
#              Based on manuscript methods: |PR - PS| ranking
# Author: ARMOR Pipeline
# Date: 2026-04-10
# ==============================================================================

#' Feature Selection for AMR Prediction
#'
#' This script performs feature selection based on the between-group difference
#' strategy as described in the manuscript:
#'
#' For each candidate feature:
#' - Calculate presence proportion in Resistant (R) group: PR
#' - Calculate presence proportion in Susceptible (S) group: PS
#' - Discriminative score = |PR - PS|
#' - Rank features by score in descending order
#' - Select top N features (default: 1000)
#'
#' Additionally performs statistical tests:
#' - Chi-square test for association (p < 0.05)
#' - Odds Ratio calculation (OR > 2)

library(optparse)
library(data.table)
library(dplyr)

# ------------------------------------------------------------------------------
# Command line options
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input feature matrix file [required]"),
  make_option(c("-g", "--group"), type = "character", default = NULL,
              help = "Sample group/phenotype file [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory for selected features [required]"),
  make_option(c("-n", "--top-n"), type = "integer", default = 1000,
              help = "Number of top features to select [default: 1000]"),
  make_option(c("-t", "--type"), type = "character", default = "auto",
              help = "Feature type: gene/snp/kmer/auto [default: auto]"),
  make_option(c("--p-threshold"), type = "double", default = 0.05,
              help = "P-value threshold for chi-square test [default: 0.05]"),
  make_option(c("--or-threshold"), type = "double", default = 2.0,
              help = "Odds ratio threshold [default: 2.0]"),
  make_option(c("--min-presence"), type = "double", default = 0.01,
              help = "Minimum presence rate to retain feature [default: 0.01]"),
  make_option(c("--max-presence"), type = "double", default = 0.99,
              help = "Maximum presence rate to retain feature [default: 0.99]"),
  make_option(c("--output-all"), action = "store_true", default = FALSE,
              help = "Output all features with scores (not just top N)")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Feature selection using between-group difference strategy")
opt <- parse_args(opt_parser)

# Validate arguments
required_args <- c("input", "group", "outdir")
for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(sprintf("Argument --%s is required!", arg))
  }
}

# ------------------------------------------------------------------------------
# Core functions
# ------------------------------------------------------------------------------

#' Calculate between-group difference score
#' @param feature_vector Numeric vector of feature values (0/1 or continuous)
#' @param group_vector Factor vector of group labels (R/S or 0/1)
#' @return List with PR, PS, score, p_value, odds_ratio
calculate_feature_score <- function(feature_vector, group_vector,
                                   min_presence = 0.01, max_presence = 0.99) {

  # Ensure binary groups
  groups <- unique(group_vector)
  if (length(groups) != 2) {
    return(list(
      PR = NA, PS = NA, score = NA,
      p_value = NA, odds_ratio = NA,
      passed_filter = FALSE
    ))
  }

  # Define resistant and susceptible groups
  # Assume "R" or "1" or "Resistant" is resistant group
  resistant_labels <- c("R", "1", "Resistant", "resistant", "RES")
  if (groups[1] %in% resistant_labels || groups[2] %in% resistant_labels) {
    resistant_group <- groups[groups %in% resistant_labels][1]
  } else {
    # Use first group as resistant by default
    resistant_group <- groups[1]
  }
  susceptible_group <- groups[groups != resistant_group][1]

  # Calculate presence proportions
  r_idx <- which(group_vector == resistant_group)
  s_idx <- which(group_vector == susceptible_group)

  if (length(r_idx) == 0 || length(s_idx) == 0) {
    return(list(
      PR = NA, PS = NA, score = NA,
      p_value = NA, odds_ratio = NA,
      passed_filter = FALSE
    ))
  }

  # For binary features (gene presence/absence, SNP, k-mer)
  # Convert to binary if needed
  feature_binary <- ifelse(feature_vector > 0, 1, 0)

  # Presence proportion in R group
  PR <- sum(feature_binary[r_idx] > 0) / length(r_idx)

  # Presence proportion in S group
  PS <- sum(feature_binary[s_idx] > 0) / length(s_idx)

  # Check presence rate filters
  overall_presence <- sum(feature_binary > 0) / length(feature_binary)
  passed_filter <- overall_presence >= min_presence & overall_presence <= max_presence

  # Calculate discriminative score
  score <- abs(PR - PS)

  # Chi-square test
  contingency_table <- table(feature_binary, group_vector)

  # Ensure table has proper dimensions
  if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2) {
    p_value <- NA
    odds_ratio <- NA
  } else {
    # Chi-square test with Yates correction for 2x2
    tryCatch({
      chi_test <- chisq.test(contingency_table, correct = TRUE)
      p_value <- chi_test$p.value
    }, warning = function(w) {
      # If expected counts are low, use Fisher's exact test
      tryCatch({
        fisher_test <- fisher.test(contingency_table)
        p_value <<- fisher_test$p.value
      }, error = function(e) {
        p_value <<- NA
      })
    }, error = function(e) {
      p_value <<- NA
    })

    # Calculate odds ratio
    tryCatch({
      # OR = (a*d) / (b*c) for table:
      #           R    S
      # feat=1    a    b
      # feat=0    c    d
      a <- contingency_table[2, 2]  # feature=1, resistant
      b <- contingency_table[2, 1]  # feature=1, susceptible
      c <- contingency_table[1, 2]  # feature=0, resistant
      d <- contingency_table[1, 1]  # feature=0, susceptible

      # Add pseudo-count to avoid division by zero
      if (b == 0 || c == 0) {
        a <- a + 0.5
        b <- b + 0.5
        c <- c + 0.5
        d <- d + 0.5
      }

      odds_ratio <- (a * d) / (b * c)
    }, error = function(e) {
      odds_ratio <<- NA
    })
  }

  return(list(
    PR = PR,
    PS = PS,
    score = score,
    p_value = p_value,
    odds_ratio = odds_ratio,
    passed_filter = passed_filter
  ))
}

#' Perform feature selection on a matrix
#' @param feature_matrix Data frame with samples as rows, features as columns
#' @param group_data Data frame with sample IDs and group labels
#' @param top_n Number of top features to select
#' @param p_threshold P-value threshold
#' @param or_threshold Odds ratio threshold
#' @return List with selected features and statistics
perform_feature_selection <- function(feature_matrix, group_data,
                                       top_n = 1000,
                                       p_threshold = 0.05,
                                       or_threshold = 2.0,
                                       min_presence = 0.01,
                                       max_presence = 0.99,
                                       output_all = FALSE) {

  # Merge feature matrix with group data
  if ("ID" %in% colnames(feature_matrix)) {
    rownames(feature_matrix) <- feature_matrix$ID
    feature_matrix <- feature_matrix[, -which(colnames(feature_matrix) == "ID")]
  }

  if ("ID" %in% colnames(group_data)) {
    rownames(group_data) <- group_data$ID
    group_data <- group_data[, -which(colnames(group_data) == "ID")]
  }

  # Get group column (assume first column or column named "Group")
  if ("Group" %in% colnames(group_data)) {
    group_vector <- group_data$Group
  } else {
    group_vector <- group_data[, 1]
  }

  # Match samples
  common_samples <- intersect(rownames(feature_matrix), names(group_vector))
  message(sprintf("Common samples: %d", length(common_samples)))

  feature_matrix <- feature_matrix[common_samples, , drop = FALSE]
  group_vector <- group_vector[common_samples]

  # Get feature columns (exclude non-feature columns)
  feature_cols <- setdiff(colnames(feature_matrix), c("Group", "ID", "Sample", "Label"))

  message(sprintf("Analyzing %d features...", length(feature_cols)))

  # Calculate scores for all features
  results <- data.frame(
    Feature = feature_cols,
    PR = numeric(length(feature_cols)),
    PS = numeric(length(feature_cols)),
    Score = numeric(length(feature_cols)),
    P_value = numeric(length(feature_cols)),
    Odds_Ratio = numeric(length(feature_cols)),
    Passed_Filter = logical(length(feature_cols)),
    stringsAsFactors = FALSE
  )

  # Process features
  for (i in seq_along(feature_cols)) {
    feat <- feature_cols[i]

    if (i %% 1000 == 0) {
      message(sprintf("  Processing feature %d/%d...", i, length(feature_cols)))
    }

    feature_vec <- feature_matrix[, feat]

    # Handle missing values
    valid_idx <- !is.na(feature_vec) & !is.na(group_vector)
    if (sum(valid_idx) < 10) {
      results[i, ] <- c(feat, NA, NA, NA, NA, NA, FALSE)
      next
    }

    score_result <- calculate_feature_score(
      feature_vector = feature_vec[valid_idx],
      group_vector = group_vector[valid_idx],
      min_presence = min_presence,
      max_presence = max_presence
    )

    results[i, ] <- c(
      feat,
      score_result$PR,
      score_result$PS,
      score_result$score,
      score_result$p_value,
      score_result$odds_ratio,
      score_result$passed_filter
    )
  }

  # Convert types
  numeric_cols <- c("PR", "PS", "Score", "P_value", "Odds_Ratio")
  results[, numeric_cols] <- lapply(results[, numeric_cols], as.numeric)
  results$Passed_Filter <- as.logical(results$Passed_Filter)

  # Sort by score (descending)
  results <- results[order(-results$Score, na.last = TRUE), ]

  # Add rank
  results$Rank <- 1:nrow(results)

  # Filter by statistical significance
  results$Significant <- !is.na(results$P_value) &
                         results$P_value < p_threshold &
                         !is.na(results$Odds_Ratio) &
                         results$Odds_Ratio > or_threshold

  # Select top N features that pass filters
  filtered_results <- results[results$Passed_Filter, ]

  if (output_all) {
    selected_features <- filtered_results$Feature
  } else {
    top_n_actual <- min(top_n, nrow(filtered_results))
    selected_features <- head(filtered_results$Feature, top_n_actual)
  }

  message(sprintf("Features passing filters: %d", sum(results$Passed_Filter, na.rm = TRUE)))
  message(sprintf("Significant features (p<%.3f, OR>%.1f): %d",
                  p_threshold, or_threshold, sum(results$Significant, na.rm = TRUE)))
  message(sprintf("Selected top %d features", length(selected_features)))

  return(list(
    all_results = results,
    selected_features = selected_features,
    statistics = list(
      total_features = nrow(results),
      passed_filter = sum(results$Passed_Filter, na.rm = TRUE),
      significant = sum(results$Significant, na.rm = TRUE),
      top_n = length(selected_features)
    )
  ))
}

#' Extract feature names from selected matrix
#' @param feature_matrix Full feature matrix
#' @param selected_features Vector of selected feature names
#' @return Subset feature matrix
extract_selected_features <- function(feature_matrix, selected_features) {
  feature_cols <- intersect(selected_features, colnames(feature_matrix))

  if ("ID" %in% colnames(feature_matrix)) {
    output <- feature_matrix[, c("ID", feature_cols), drop = FALSE]
  } else {
    output <- feature_matrix[, feature_cols, drop = FALSE]
  }

  return(output)
}

# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

main <- function() {
  # Create output directory
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }

  # Load input data
  message("=== Loading Input Data ===")

  message(sprintf("Loading feature matrix: %s", opt$input))
  feature_matrix <- fread(opt$input, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE)

  message(sprintf("Loading group file: %s", opt$group))
  group_data <- read.table(opt$group, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE)

  message(sprintf("Feature matrix: %d samples x %d features",
                  nrow(feature_matrix), ncol(feature_matrix) - 1))
  message(sprintf("Group file: %d samples", nrow(group_data)))

  # Determine feature type
  if (opt$type == "auto") {
    if (grepl("gene", basename(opt$input), ignore.case = TRUE)) {
      opt$type <- "gene"
    } else if (grepl("snp", basename(opt$input), ignore.case = TRUE)) {
      opt$type <- "snp"
    } else if (grepl("kmer", basename(opt$input), ignore.case = TRUE)) {
      opt$type <- "kmer"
    } else {
      opt$type <- "unknown"
    }
  }
  message(sprintf("Feature type: %s", opt$type))

  # Perform feature selection
  message("=== Performing Feature Selection ===")

  selection_results <- perform_feature_selection(
    feature_matrix = feature_matrix,
    group_data = group_data,
    top_n = opt$top_n,
    p_threshold = opt$p_threshold,
    or_threshold = opt$or_threshold,
    min_presence = opt$min_presence,
    max_presence = opt$max_presence,
    output_all = opt$output_all
  )

  # Save results
  message("=== Saving Results ===")

  # Save all feature statistics
  all_stats_file <- file.path(opt$outdir, "feature_selection_all_stats.txt")
  fwrite(selection_results$all_results, all_stats_file,
         sep = "\t", col.names = TRUE, row.names = FALSE)
  message(sprintf("All feature statistics: %s", all_stats_file))

  # Save selected feature names
  selected_file <- file.path(opt$outdir, "selected_features.txt")
  writeLines(selection_results$selected_features, selected_file)
  message(sprintf("Selected features: %s", selected_file))

  # Save selected feature matrix
  selected_matrix <- extract_selected_features(feature_matrix,
                                                selection_results$selected_features)
  selected_matrix_file <- file.path(opt$outdir, "selected_feature_matrix.txt")
  fwrite(selected_matrix, selected_matrix_file,
         sep = "\t", col.names = TRUE, row.names = FALSE)
  message(sprintf("Selected feature matrix: %s", selected_matrix_file))

  # Save summary
  summary_file <- file.path(opt$outdir, "feature_selection_summary.txt")
  summary_lines <- c(
    "=== Feature Selection Summary ===",
    sprintf("Input features: %d", selection_results$statistics$total_features),
    sprintf("Passed presence filter: %d", selection_results$statistics$passed_filter),
    sprintf("Statistically significant: %d", selection_results$statistics$significant),
    sprintf("Top N selected: %d", selection_results$statistics$top_n),
    sprintf("P-value threshold: %.3f", opt$p_threshold),
    sprintf("Odds ratio threshold: %.1f", opt$or_threshold),
    sprintf("Min presence: %.2f", opt$min_presence),
    sprintf("Max presence: %.2f", opt$max_presence),
    "",
    "Top 10 features by score:",
    paste("  ", head(selection_results$all_results$Feature, 10), collapse = "\n")
  )
  writeLines(summary_lines, summary_file)
  message(sprintf("Summary: %s", summary_file))

  # Create visualization of score distribution
  message("=== Generating Visualization ===")

  tryCatch({
    library(ggplot2)

    # Score distribution
    p1 <- ggplot(selection_results$all_results, aes(x = Score)) +
      geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
      theme_bw() +
      xlab("Discriminative Score (|PR - PS|)") +
      ylab("Number of Features") +
      ggtitle("Distribution of Feature Discriminative Scores")

    ggsave(file.path(opt$outdir, "score_distribution.png"), p1,
           width = 8, height = 6, dpi = 150)

    # Top features bar plot
    top_features <- head(selection_results$all_results, 20)
    p2 <- ggplot(top_features, aes(x = reorder(Feature, Score), y = Score)) +
      geom_bar(stat = "identity", fill = "coral") +
      coord_flip() +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8)) +
      xlab("Feature") +
      ylab("Discriminative Score") +
      ggtitle("Top 20 Features by Discriminative Score")

    ggsave(file.path(opt$outdir, "top_features.png"), p2,
           width = 10, height = 8, dpi = 150)

    message("Visualizations saved to output directory")
  }, error = function(e) {
    message(sprintf("Could not generate visualization: %s", e$message))
  })

  message("=== Feature Selection Complete ===")
}

# Run main function
main()
