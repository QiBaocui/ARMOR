#!/usr/bin/env Rscript

# ==============================================================================
# Script: 09_predict_new_samples.R
# Description: Apply trained models to predict resistance in new samples
# Author: ARMOR Pipeline
# Date: 2026-04-10
# ==============================================================================

library(optparse)
library(data.table)
library(dplyr)
library(h2o)

# ------------------------------------------------------------------------------
# Command line options
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Feature matrix file for new samples [required]"),
  make_option(c("-m", "--model"), type = "character", default = NULL,
              help = "Path to trained H2O model directory [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory for predictions [required]"),
  make_option(c("-t", "--type"), type = "character", default = "Gene",
              help = "Feature type: Gene/SNP/Kmer/All [default: Gene]"),
  make_option(c("-p", "--port"), type = "integer", default = 54321,
              help = "H2O cluster port [default: 54321]"),
  make_option(c("--threshold"), type = "double", default = NULL,
              help = "Custom prediction threshold (default: use Youden optimal)")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Predict resistance in new samples using trained models")
opt <- parse_args(opt_parser)

# Validate arguments
required_args <- c("input", "model", "outdir")
for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(sprintf("Argument --%s is required!", arg))
  }
}

if (!file.exists(opt$input)) {
  stop(sprintf("Input file not found: %s", opt$input))
}

if (!dir.exists(opt$model)) {
  stop(sprintf("Model directory not found: %s", opt$model))
}

# Create output directory
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

#' Load trained H2O model
#' @param model_dir Directory containing saved model
#' @return H2O model object
load_trained_model <- function(model_dir) {
  # Find model file
  model_files <- list.files(model_dir, pattern = "\\.zip$", full.names = TRUE)

  if (length(model_files) == 0) {
    # Try to find model directory
    model_dirs <- list.dirs(model_dir, recursive = FALSE)

    if (length(model_dirs) == 0) {
      stop(sprintf("No model found in: %s", model_dir))
    }

    model_files <- list.files(model_dirs[1], pattern = "\\.zip$", full.names = TRUE)
  }

  if (length(model_files) == 0) {
    stop("No model file (.zip) found!")
  }

  message(sprintf("Loading model: %s", model_files[1]))
  model <- h2o.loadModel(model_files[1])

  return(model)
}

#' Prepare feature matrix for prediction
#' @param input_file Feature matrix file
#' @param model Trained H2O model
#' @return H2O frame ready for prediction
prepare_data <- function(input_file, model) {
  message(sprintf("Loading feature matrix: %s", input_file))

  # Load data
  data <- fread(input_file, sep = "\t", header = TRUE,
                stringsAsFactors = FALSE, check.names = FALSE)

  # Store sample IDs
  if ("ID" %in% colnames(data)) {
    sample_ids <- data$ID
    data <- data[, -which(colnames(data) == "ID")]
  } else {
    sample_ids <- rownames(data)
    rownames(data) <- NULL
  }

  message(sprintf("Loaded %d samples with %d features",
                  nrow(data), ncol(data)))

  # Get required features from model
  # Note: This depends on the model structure
  # For now, assume all columns are features

  # Convert to H2O frame
  data.h2o <- as.h2o(data)

  return(list(data = data, data.h2o = data.h2o, sample_ids = sample_ids))
}

#' Generate predictions
#' @param model Trained H2O model
#' @param data.h2o H2O frame with features
#' @param sample_ids Sample identifiers
#' @param threshold Prediction threshold
#' @return Data frame with predictions
generate_predictions <- function(model, data.h2o, sample_ids, threshold = NULL) {
  message("Generating predictions...")

  # Get raw predictions
  pred <- h2o.predict(object = model, newdata = data.h2o) %>%
    as.data.frame()

  # Add sample IDs
  pred$SampleID <- sample_ids

  # Reorder columns
  if ("p0" %in% colnames(pred) && "p1" %in% colnames(pred)) {
    # Binary classification: p0 = S, p1 = R
    pred <- pred[, c("SampleID", "p0", "p1", "predict")]
    colnames(pred) <- c("SampleID", "Prob_S", "Prob_R", "Prediction")
  } else {
    colnames(pred)[1:2] <- c("Prob_S", "Prob_R")
  }

  # Apply threshold if provided
  if (!is.null(threshold)) {
    pred$Prediction <- ifelse(pred$Prob_R >= threshold, "R", "S")
    message(sprintf("Using custom threshold: %.3f", threshold))
  }

  return(pred)
}

# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

main <- function() {
  message("=== ARMOR Prediction Pipeline ===")
  message(sprintf("Input: %s", opt$input))
  message(sprintf("Model: %s", opt$model))
  message(sprintf("Output: %s", opt$outdir))

  # Initialize H2O
  message("\n=== Initializing H2O ===")
  h2o.init(ip = "localhost", port = opt$port)

  # Load model
  message("\n=== Loading Trained Model ===")
  model <- load_trained_model(opt$model)
  message(sprintf("Model type: %s", class(model)[1]))

  # Prepare data
  message("\n=== Preparing Data ===")
  data_info <- prepare_data(opt$input, model)

  # Generate predictions
  message("\n=== Generating Predictions ===")
  predictions <- generate_predictions(
    model = model,
    data.h2o = data_info$data.h2o,
    sample_ids = data_info$sample_ids,
    threshold = opt$threshold
  )

  # Save predictions
  message("\n=== Saving Results ===")

  # Main predictions file
  pred_file <- file.path(opt$outdir, "predictions.tsv")
  fwrite(predictions, pred_file, sep = "\t")
  message(sprintf("Predictions saved: %s", pred_file))

  # Summary statistics
  summary_stats <- predictions %>%
    count(Prediction) %>%
    mutate(Percentage = n / sum(n) * 100)

  summary_file <- file.path(opt$outdir, "prediction_summary.tsv")
  fwrite(summary_stats, summary_file, sep = "\t")
  message(sprintf("Summary saved: %s", summary_file))

  message("\nPrediction Summary:")
  print(summary_stats)

  # Cleanup
  h2o.shutdown(prompt = FALSE)

  message("\n=== Prediction Complete ===")
}

main()
