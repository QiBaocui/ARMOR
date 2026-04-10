#!/usr/bin/env Rscript

# ==============================================================================
# Script: 03_kmer_feature_generation.R
# Description: Generate k-mer frequency matrix from WGS data using Jellyfish
#              Based on manuscript methods: 11-mer canonical counting with filtering
# Author: ARMOR Pipeline
# Date: 2026-04-10
# ==============================================================================

library(optparse)
library(data.table)

# ------------------------------------------------------------------------------
# Command line options
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Directory containing FASTQ files [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory for k-mer matrix [required]"),
  make_option(c("-k", "--kmer-size"), type = "integer", default = 11,
              help = "K-mer size [default: 11]"),
  make_option(c("-t", "--threads"), type = "integer", default = 8,
              help = "Number of threads [default: 8]"),
  make_option(c("-s", "--sketch-size"), type = "character", default = "100M",
              help = "Jellyfish hash size (e.g., 100M, 1G) [default: 100M]"),
  make_option(c("--min-count"), type = "integer", default = 2,
              help = "Minimum k-mer count threshold [default: 2]"),
  make_option(c("--min-presence"), type = "double", default = 0.01,
              help = "Minimum presence rate across samples [default: 0.01]"),
  make_option(c("--max-presence"), type = "double", default = 0.99,
              help = "Maximum presence rate across samples [default: 0.99]"),
  make_option(c("--normalize"), action = "store_true", default = FALSE,
              help = "Normalize k-mer counts by total counts per sample")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Generate k-mer frequency matrix using Jellyfish")
opt <- parse_args(opt_parser)

# Validate arguments
if (is.null(opt$input) || is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Both --input and --outdir are required!")
}

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

#' Check if command exists
check_command <- function(cmd) {
  result <- system2("which", cmd, stdout = NULL, stderr = NULL)
  return(result == 0)
}

#' Run Jellyfish count for a single sample
#' @param input_file Path to FASTQ file
#' @param out_prefix Output prefix
#' @param kmer_size K-mer size
#' @param sketch_size Hash size
#' @param threads Number of threads
#' @return Path to jellyfish database
run_jellyfish_count <- function(input_file, out_prefix, kmer_size = 11,
                                 sketch_size = "100M", threads = 8) {
  # Check if paired-end
  is_paired <- grepl("_R1|_r1|_1\\.", input_file, ignore.case = TRUE)

  if (is_paired) {
    r1_file <- input_file
    r2_file <- gsub("_R1|_r1|_1\\.", "_R2_", r1_file)
    if (!file.exists(r2_file)) {
      r2_file <- gsub("_1\\.", "_2.", r1_file)
    }

    # Count k-mers from both reads
    cmd <- sprintf("jellyfish count -C -m %d -s %s -t %d -o %s.jf %s %s",
                   kmer_size, sketch_size, threads, out_prefix, r1_file, r2_file)
  } else {
    cmd <- sprintf("jellyfish count -C -m %d -s %s -t %d -o %s.jf %s",
                   kmer_size, sketch_size, threads, out_prefix, input_file)
  }

  system(cmd)

  return(paste0(out_prefix, ".jf"))
}

#' Dump Jellyfish database to text format
#' @param jf_file Path to Jellyfish database
#' @param out_file Output text file
dump_jellyfish <- function(jf_file, out_file) {
  cmd <- sprintf("jellyfish dump -o %s %s", out_file, jf_file)
  system(cmd)
  return(out_file)
}

#' Convert Jellyfish output to tabular format
#' @param dump_file Path to Jellyfish dump file
#' @param out_file Output tabular file
convert_to_tsv <- function(dump_file, out_file) {
  # Read FASTA-like format and convert to TSV
  lines <- readLines(dump_file)

  kmer_data <- data.frame(
    kmer = character(),
    count = integer(),
    stringsAsFactors = FALSE
  )

  for (i in seq(1, length(lines), by = 2)) {
    if (i + 1 <= length(lines)) {
      kmer <- gsub(">", "", lines[i])
      count <- lines[i + 1]
      kmer_data <- rbind(kmer_data, data.frame(kmer = kmer, count = as.integer(count)))
    }
  }

  fwrite(kmer_data, out_file, sep = "\t", col.names = TRUE)
  return(out_file)
}

#' Merge k-mer counts from multiple samples
#' @param sample_files Named list of sample k-mer count files
#' @param out_file Output merged matrix file
#' @param min_count Minimum count threshold
#' @param min_presence Minimum presence rate
#' @param max_presence Maximum presence rate
#' @param normalize Whether to normalize counts
merge_kmer_counts <- function(sample_files, out_file, min_count = 2,
                               min_presence = 0.01, max_presence = 0.99,
                               normalize = FALSE) {

  message("Merging k-mer counts across samples...")

  # Read all sample files
  sample_data <- list()
  total_counts <- list()

  for (sample_name in names(sample_files)) {
    message(sprintf("  Reading %s...", sample_name))

    # Read k-mer counts
    kmer_df <- fread(sample_files[[sample_name]], sep = "\t",
                     header = FALSE, col.names = c("kmer", "count"))

    # Filter by minimum count
    kmer_df <- kmer_df[kmer_df$count >= min_count, ]

    # Store as named vector
    sample_data[[sample_name]] <- setNames(kmer_df$count, kmer_df$kmer)

    # Store total for normalization
    total_counts[[sample_name]] <- sum(kmer_df$count)
  }

  # Get all unique k-mers
  all_kmers <- unique(unlist(lapply(sample_data, names)))
  message(sprintf("Total unique k-mers: %d", length(all_kmers)))

  # Create matrix
  kmer_matrix <- matrix(0, nrow = length(all_kmers), ncol = length(sample_data))
  rownames(kmer_matrix) <- all_kmers
  colnames(kmer_matrix) <- names(sample_data)

  # Fill matrix
  for (sample_name in names(sample_data)) {
    kmer_vec <- sample_data[[sample_name]]
    matching_kmers <- intersect(names(kmer_vec), all_kmers)
    kmer_matrix[matching_kmers, sample_name] <- kmer_vec[matching_kmers]
  }

  # Filter by presence rate
  presence_rate <- rowSums(kmer_matrix > 0) / ncol(kmer_matrix)

  message(sprintf("Filtering k-mers by presence rate (%.1f%% - %.1f%%)...",
                  min_presence * 100, max_presence * 100))

  keep_kmers <- presence_rate >= min_presence & presence_rate <= max_presence
  message(sprintf("  K-mers before filtering: %d", nrow(kmer_matrix)))
  message(sprintf("  K-mers after filtering: %d", sum(keep_kmers)))

  kmer_matrix <- kmer_matrix[keep_kmers, ]

  # Normalize if requested
  if (normalize) {
    message("Normalizing k-mer counts by total counts per sample...")
    for (sample_name in colnames(kmer_matrix)) {
      kmer_matrix[, sample_name] <- kmer_matrix[, sample_name] / total_counts[[sample_name]] * 1e6
    }
  }

  # Transpose to have samples as rows
  kmer_matrix <- t(kmer_matrix)

  # Convert to data frame with ID column
  kmer_df <- as.data.frame(kmer_matrix)
  kmer_df <- cbind(ID = rownames(kmer_df), kmer_df)

  # Write output
  message(sprintf("Writing k-mer matrix to: %s", out_file))
  fwrite(kmer_df, out_file, sep = "\t", col.names = TRUE)

  message(sprintf("Final dimensions: %d samples x %d k-mers",
                  nrow(kmer_df), ncol(kmer_df) - 1))

  return(kmer_df)
}

#' Alternative: Use KMC for efficient k-mer counting
run_kmc <- function(input_file, out_prefix, kmer_size = 11, threads = 8) {
  # KMC is faster and more memory-efficient for large datasets

  is_paired <- grepl("_R1|_r1|_1\\.", input_file, ignore.case = TRUE)

  if (is_paired) {
    r1_file <- input_file
    r2_file <- gsub("_R1|_r1|_1\\.", "_R2_", r1_file)
    if (!file.exists(r2_file)) {
      r2_file <- gsub("_1\\.", "_2.", r1_file)
    }

    # Create file list for KMC
    file_list <- paste0(out_prefix, "_filelist.txt")
    writeLines(c(r1_file, r2_file), file_list)

    cmd <- sprintf("kmc -k%d -t%d -m64 -ci2 -cs10000 @%s %s %s",
                   kmer_size, threads, file_list, out_prefix, tempdir())
  } else {
    cmd <- sprintf("kmc -k%d -t%d -m64 -ci2 -cs10000 %s %s %s",
                   kmer_size, threads, input_file, out_prefix, tempdir())
  }

  system(cmd)

  # Dump results
  dump_cmd <- sprintf("kmc_tools transform %s dump %s.kmc.txt",
                      out_prefix, out_prefix)
  system(dump_cmd)

  return(paste0(out_prefix, ".kmc.txt"))
}

# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

main <- function() {
  # Create output directory
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }

  # Create temp directory for intermediate files
  temp_dir <- file.path(opt$outdir, "temp")
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }

  # Get list of input files
  input_patterns <- c("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
  input_files <- c()
  for (pattern in input_patterns) {
    files <- list.files(opt$input, pattern = pattern,
                        full.names = TRUE, ignore.case = TRUE)
    input_files <- c(input_files, files)
  }

  # Only keep R1 files for paired-end
  input_files <- input_files[!grepl("_R2|_r2|_2\\.", input_files, ignore.case = TRUE)]

  if (length(input_files) == 0) {
    stop(sprintf("No FASTQ files found in %s", opt$input))
  }

  message(sprintf("Found %d FASTQ files", length(input_files)))

  # Check for k-mer counting tools
  use_kmc <- FALSE
  if (check_command("kmc")) {
    message("Using KMC for k-mer counting (recommended)")
    use_kmc <- TRUE
  } else if (check_command("jellyfish")) {
    message("Using Jellyfish for k-mer counting")
  } else {
    stop("Neither KMC nor Jellyfish found. Please install one of them.")
  }

  # Step 1: Count k-mers for each sample
  message("=== Step 1: K-mer Counting ===")

  sample_files <- list()

  for (input_file in input_files) {
    sample_id <- tools::file_path_sans_ext(basename(input_file))
    sample_id <- gsub("_R1$", "", sample_id)

    message(sprintf("Processing %s...", sample_id))

    out_prefix <- file.path(temp_dir, sample_id)

    if (use_kmc) {
      kmer_file <- run_kmc(input_file, out_prefix, opt$kmer_size, opt$threads)
    } else {
      jf_file <- run_jellyfish_count(input_file, out_prefix,
                                      opt$kmer_size, opt$sketch_size, opt$threads)
      dump_file <- dump_jellyfish(jf_file, paste0(out_prefix, ".dump"))
      kmer_file <- convert_to_tsv(dump_file, paste0(out_prefix, ".tsv"))
    }

    sample_files[[sample_id]] <- kmer_file
  }

  # Step 2: Merge and filter k-mers
  message("=== Step 2: Merge and Filter K-mers ===")

  output_file <- file.path(opt$outdir, "kmer_matrix.txt")

  kmer_matrix <- merge_kmer_counts(
    sample_files = sample_files,
    out_file = output_file,
    min_count = opt$min_count,
    min_presence = opt$min_presence,
    max_presence = opt$max_presence,
    normalize = opt$normalize
  )

  # Clean up temp files
  message("Cleaning up temporary files...")
  unlink(temp_dir, recursive = TRUE)

  message("=== K-mer Feature Generation Complete ===")
  message(sprintf("Output: %s", output_file))
}

# Run main function
main()
