#!/usr/bin/env Rscript

# ==============================================================================
# File: basic_analyse.R
# Description: Main analysis pipeline for AMR prediction using WGS data
#              This is the primary entry point for running the full analysis
#
# Workflow:
#   1. Load input data (gene matrix, SNP matrix, k-mer matrix, group file)
#   2. Split samples into training and testing sets
#   3. Merge features with group labels
#   4. Initialize H2O cluster
#   5. Train models for each feature type (Gene, SNP, Kmer, All)
#   6. Save results and trained models
#
# Usage:
#   Rscript basic_analyse.R -o outdir -g group.txt -G gene.txt -s snp.txt -k kmer.txt
#
# Author: ARMOR Pipeline
# Reference: Manuscript.docx - Methods section
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
source("model.R")  # Load modeling functions

library(dplyr)
library(data.table)
library(ggplot2)
library(optparse)
library(caret)
library(h2o)

# ------------------------------------------------------------------------------
# Command line options
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = NULL,
    help = "Output directory for results [required]"
  ),
  make_option(
    c("-g", "--group"),
    type = "character",
    default = NULL,
    help = "Sample group/phenotype file [required]"
  ),
  make_option(
    c("-G", "--gene"),
    type = "character",
    default = NULL,
    help = "Gene presence/absence matrix file [required]"
  ),
  make_option(
    c("-s", "--snp"),
    type = "character",
    default = NULL,
    help = "SNP matrix file [required]"
  ),
  make_option(
    c("-k", "--kmer"),
    type = "character",
    default = NULL,
    help = "K-mer frequency matrix file [required]"
  ),
  make_option(
    c("-t", "--trainNum"),
    type = "numeric",
    default = 0.6,
    help = "Training set proportion [default: 0.6]"
  ),
  make_option(
    c("-S", "--split_type"),
    type = "character",
    default = "proportional",
    help = "Split method: 'proportional' or 'sample_size' [default: proportional]"
  ),
  make_option(
    c("-p", "--ip"),
    type = "integer",
    default = 54321,
    help = "H2O cluster port [default: 54321]"
  ),
  make_option(
    c("--nthreads"),
    type = "integer",
    default = -1,
    help = "Number of CPU threads for H2O (-1 = all) [default: -1]"
  ),
  make_option(
    c("--nfolds"),
    type = "integer",
    default = 10,
    help = "Number of cross-validation folds [default: 10]"
  ),
  make_option(
    c("--algorithms"),
    type = "character",
    default = "GBM,GLM,RF,DL",
    help = "ML algorithms to train (comma-separated) [default: GBM,GLM,RF,DL]"
  ),
  make_option(
    c("--skip-snp"),
    action = "store_true",
    default = FALSE,
    help = "Skip SNP model training"
  ),
  make_option(
    c("--skip-kmer"),
    action = "store_true",
    default = FALSE,
    help = "Skip K-mer model training"
  ),
  make_option(
    c("--gene-only"),
    action = "store_true",
    default = FALSE,
    help = "Train gene models only"
  ),
  make_option(
    c("--seed"),
    type = "integer",
    default = 42,
    help = "Random seed for reproducibility [default: 42]"
  )
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "ARMOR: AMR Prediction Pipeline using WGS data",
  epilogue = "
Example usage:
  # Full analysis with all feature types
  Rscript basic_analyse.R -o results -g group.txt -G gene_matrix.txt \\
    -s SNP_matrix.txt -k kmer_matrix.txt

  # Gene-only analysis
  Rscript basic_analyse.R -o results -g group.txt -G gene_matrix.txt \\
    -s SNP_matrix.txt -k kmer_matrix.txt --gene-only

  # Custom algorithms
  Rscript basic_analyse.R -o results -g group.txt -G gene_matrix.txt \\
    -s SNP_matrix.txt -k kmer_matrix.txt --algorithms GBM,RF,stacked
"
)

opt <- parse_args(opt_parser)

# ------------------------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------------------------
validate_args <- function() {
  required_args <- c("outdir", "group", "gene", "snp", "kmer")
  missing <- character()

  for (arg in required_args) {
    if (is.null(opt[[arg]])) {
      missing <- c(missing, arg)
    } else if (!file.exists(opt[[arg]])) {
      stop(sprintf("File not found: %s = %s", arg, opt[[arg]]))
    }
  }

  if (length(missing) > 0) {
    print_help(opt_parser)
    stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))
  }

  # Validate split_type
  if (!opt$split_type %in% c("proportional", "sample_size")) {
    stop("split_type must be 'proportional' or 'sample_size'")
  }

  # Parse algorithms
  opt$algorithms <- strsplit(opt$algorithms, ",")[[1]]
  opt$algorithms <- trimws(opt$algorithms)

  valid_algos <- c("GBM", "GLM", "RF", "DL", "xGboost", "Bayes", "stacked")
  invalid <- setdiff(opt$algorithms, valid_algos)
  if (length(invalid) > 0) {
    stop(sprintf("Invalid algorithms: %s. Valid options: %s",
                 paste(invalid, collapse = ", "),
                 paste(valid_algos, collapse = ", ")))
  }
}

# ------------------------------------------------------------------------------
# Main analysis function
# ------------------------------------------------------------------------------
run_analysis <- function() {
  # Set random seed for reproducibility
  set.seed(opt$seed)

  # Create output directory
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }

  # Create log directory
  log_dir <- file.path(opt$outdir, "logs")
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }

  # Start logging
  log_file <- file.path(log_dir, sprintf("analysis_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
  sink(log_file)
  sink(log_file, type = "message")

  message("========================================")
  message("ARMOR Pipeline Started")
  message("========================================")
  message(sprintf("Start time: %s", Sys.time()))
  message(sprintf("Output directory: %s", opt$outdir))
  message(sprintf("Random seed: %d", opt$seed))
  message("")

  # ---------------------------------------------------------------------------
  # Step 1: Load input data
  # ---------------------------------------------------------------------------
  message("=== Step 1: Loading Input Data ===")

  message(sprintf("Loading group file: %s", opt$group))
  group <- read.table(opt$group, sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)
  message(sprintf("  Samples: %d", nrow(group)))
  message(sprintf("  Groups: %s", paste(unique(group$Group), collapse = ", ")))

  # Handle gene-only or full analysis
  if (opt$gene_only) {
    opt$skip_snp <- TRUE
    opt$skip_kmer <- TRUE
    message("Gene-only mode enabled")
  }

  # Load feature matrices
  message(sprintf("Loading gene matrix: %s", opt$gene))
  gene_dat <- fread(opt$gene, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    as.data.frame()
  message(sprintf("  Dimensions: %d samples x %d genes",
                  nrow(gene_dat), ncol(gene_dat) - 1))

  if (!opt$skip_snp) {
    message(sprintf("Loading SNP matrix: %s", opt$snp))
    snp_dat <- fread(opt$snp, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
      as.data.frame()
    message(sprintf("  Dimensions: %d samples x %d SNPs",
                    nrow(snp_dat), ncol(snp_dat) - 1))
  }

  if (!opt$skip_kmer) {
    message(sprintf("Loading k-mer matrix: %s", opt$kmer))
    kmer_dat <- fread(opt$kmer, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
      as.data.frame()
    message(sprintf("  Dimensions: %d samples x %d k-mers",
                    nrow(kmer_dat), ncol(kmer_dat) - 1))
  }

  # ---------------------------------------------------------------------------
  # Step 2: Split samples into training and testing sets
  # ---------------------------------------------------------------------------
  message("")
  message("=== Step 2: Splitting Samples ===")

  message(sprintf("Split type: %s", opt$split_type))
  if (opt$split_type == "proportional") {
    message(sprintf("Training proportion: %.1f%%", opt$trainNum * 100))
  } else {
    message(sprintf("Sample size method: %d total samples", opt$trainNum))
  }

  group_split <- splitSample(
    group = group,
    split_type = opt$split_type,
    ratio = opt$trainNum,
    num = opt$trainNum
  )

  train_group <- group_split$train_group
  test_group <- group_split$test_group

  message(sprintf("Training set: %d samples", nrow(train_group)))
  message(sprintf("  Resistant: %d (%.1f%%)",
                  sum(train_group$Group == "R"),
                  sum(train_group$Group == "R") / nrow(train_group) * 100))
  message(sprintf("  Susceptible: %d (%.1f%%)",
                  sum(train_group$Group == "S"),
                  sum(train_group$Group == "S") / nrow(train_group) * 100))

  message(sprintf("Test set: %d samples", nrow(test_group)))
  message(sprintf("  Resistant: %d (%.1f%%)",
                  sum(test_group$Group == "R"),
                  sum(test_group$Group == "R") / nrow(test_group) * 100))
  message(sprintf("  Susceptible: %d (%.1f%%)",
                  sum(test_group$Group == "S"),
                  sum(test_group$Group == "S") / nrow(test_group) * 100))

  # ---------------------------------------------------------------------------
  # Step 3: Merge features with group labels
  # ---------------------------------------------------------------------------
  message("")
  message("=== Step 3: Merging Features ===")

  # Gene features
  gene_train <- left_join(train_group, gene_dat, by = "ID")
  gene_test <- left_join(test_group, gene_dat, by = "ID")
  message(sprintf("Gene train: %d samples x %d features",
                  nrow(gene_train), ncol(gene_train) - 1))

  # SNP features
  if (!opt$skip_snp) {
    snp_train <- left_join(train_group, snp_dat, by = "ID")
    snp_test <- left_join(test_group, snp_dat, by = "ID")
    message(sprintf("SNP train: %d samples x %d features",
                    nrow(snp_train), ncol(snp_train) - 1))
  }

  # K-mer features
  if (!opt$skip_kmer) {
    kmer_train <- left_join(train_group, kmer_dat, by = "ID")
    kmer_test <- left_join(test_group, kmer_dat, by = "ID")
    message(sprintf("K-mer train: %d samples x %d features",
                    nrow(kmer_train), ncol(kmer_train) - 1))
  }

  # Combined features (All)
  train <- gene_train
  test <- gene_test

  if (!opt$skip_snp) {
    train <- left_join(train, snp_train[, -which(colnames(snp_train) == "ID")], by = "ID")
    test <- left_join(test, snp_test[, -which(colnames(snp_test) == "ID")], by = "ID")
  }

  if (!opt$skip_kmer) {
    train <- left_join(train, kmer_train[, -which(colnames(kmer_train) == "ID")], by = "ID")
    test <- left_join(test, kmer_test[, -which(colnames(kmer_test) == "ID")], by = "ID")
  }

  message(sprintf("Combined train: %d samples x %d features",
                  nrow(train), ncol(train) - 1))

  # ---------------------------------------------------------------------------
  # Step 4: Prepare data for H2O
  # ---------------------------------------------------------------------------
  message("")
  message("=== Step 4: Preparing Data Frames ===")

  # Set row names for R data frames
  prepare_data <- function(df) {
    rownames(df) <- df[, 1]
    df <- df[, -1]
    return(df)
  }

  gene_train <- prepare_data(gene_train)
  gene_test <- prepare_data(gene_test)

  if (!opt$skip_snp) {
    snp_train <- prepare_data(snp_train)
    snp_test <- prepare_data(snp_test)
  }

  if (!opt$skip_kmer) {
    kmer_train <- prepare_data(kmer_train)
    kmer_test <- prepare_data(kmer_test)
  }

  train <- prepare_data(train)
  test <- prepare_data(test)

  # Save processed data
  dat <- list(
    train_group = train_group,
    test_group = test_group,
    train = train,
    test = test,
    gene_train = gene_train,
    gene_test = gene_test
  )

  if (!opt$skip_snp) {
    dat$snp_train <- snp_train
    dat$snp_test <- snp_test
  }

  if (!opt$skip_kmer) {
    dat$kmer_train <- kmer_train
    dat$kmer_test <- kmer_test
  }

  save(dat, file = paste0(opt$outdir, "/train_data.Rdata"))
  message(sprintf("Saved training data: %s/train_data.Rdata", opt$outdir))

  # ---------------------------------------------------------------------------
  # Step 5: Initialize H2O cluster
  # ---------------------------------------------------------------------------
  message("")
  message("=== Step 5: Initializing H2O Cluster ===")

  h2o.init(
    ip = "localhost",
    port = opt$ip,
    nthreads = opt$nthreads,
    min_mem_size = "8G",
    max_mem_size = "32G"
  )

  message(sprintf("H2O version: %s", h2o.getVersion()))
  message(sprintf("Cluster status: %s", h2o.clusterStatus()))

  # ---------------------------------------------------------------------------
  # Step 6: Train models
  # ---------------------------------------------------------------------------
  message("")
  message("=== Step 6: Training Models ===")

  # Gene-based models
  message("")
  message("--- Gene-based Models ---")
  modelFun(
    train.R = gene_train,
    test.R = gene_test,
    outdir = opt$outdir,
    feature = "Gene",
    feature_type_list = opt$algorithms,
    nfolds = opt$nfolds
  )

  # SNP-based models
  if (!opt$skip_snp) {
    message("")
    message("--- SNP-based Models ---")
    modelFun(
      train.R = snp_train,
      test.R = snp_test,
      outdir = opt$outdir,
      feature = "SNP",
      feature_type_list = opt$algorithms,
      nfolds = opt$nfolds
    )
  }

  # K-mer-based models
  if (!opt$skip_kmer) {
    message("")
    message("--- K-mer-based Models ---")
    modelFun(
      train.R = kmer_train,
      test.R = kmer_test,
      outdir = opt$outdir,
      feature = "Kmer",
      feature_type_list = opt$algorithms,
      nfolds = opt$nfolds
    )
  }

  # Combined models (All features)
  if (!opt$gene_only && (!opt$skip_snp || !opt$skip_kmer)) {
    message("")
    message("--- Combined Models (All Features) ---")
    modelFun(
      train.R = train,
      test.R = test,
      outdir = opt$outdir,
      feature = "All",
      feature_type_list = opt$algorithms,
      nfolds = opt$nfolds
    )
  }

  # ---------------------------------------------------------------------------
  # Step 7: Cleanup
  # ---------------------------------------------------------------------------
  message("")
  message("=== Step 7: Cleanup ===")

  h2o.shutdown(prompt = FALSE)
  message("H2O cluster shut down")

  # Restore output
  sink()
  sink(type = "message")

  message("")
  message("========================================")
  message("ARMOR Pipeline Completed")
  message("========================================")
  message(sprintf("End time: %s", Sys.time()))
  message(sprintf("Results saved to: %s", opt$outdir))
}

# ------------------------------------------------------------------------------
# Run main analysis
# ------------------------------------------------------------------------------
validate_args()
run_analysis()
