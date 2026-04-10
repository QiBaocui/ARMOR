#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_gene_feature_generation.R
# Description: Generate gene presence/absence matrix from WGS data
#              Based on manuscript methods: Bakta annotation + Panaroo pangenome
# Author: ARMOR Pipeline
# Date: 2026-04-10
# ==============================================================================

library(optparse)
library(system2)

# ------------------------------------------------------------------------------
# Command line options
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Directory containing assembled genome FASTA files [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory for gene matrix [required]"),
  make_option(c("-a", "--annotator"), type = "character", default = "bakta",
              help = "Annotation tool (bakta/prokka) [default: bakta]"),
  make_option(c("-t", "--threads"), type = "integer", default = 8,
              help = "Number of threads for annotation [default: 8]"),
  make_option(c("--skip-annotation"), action = "store_true", default = FALSE,
              help = "Skip annotation step if GFF files already exist"),
  make_option(c("--pangenome-only"), action = "store_true", default = FALSE,
              help = "Only run pangenome analysis (expect GFF files in input dir)")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Generate gene presence/absence matrix from WGS assemblies")
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
#' @param cmd Command name to check
#' @return TRUE if command exists, FALSE otherwise
check_command <- function(cmd) {
  result <- system2("which", cmd, stdout = NULL, stderr = NULL)
  return(result == 0)
}

#' Run Bakta annotation
#' @param fasta_path Path to input FASTA file
#' @param out_dir Output directory for annotation results
#' @param threads Number of threads to use
#' @return Path to generated GFF file
run_bakta <- function(fasta_path, out_dir, threads = 8) {
  sample_id <- tools::file_path_sans_ext(basename(fasta_path))
  sample_out <- file.path(out_dir, sample_id)

  cmd <- sprintf("bakta --db db --output %s --threads %d --prefix %s %s",
                 sample_out, threads, sample_id, fasta_path)

  message(sprintf("[%s] Running Bakta annotation...", sample_id))
  system(cmd)

  gff_file <- file.path(sample_out, paste0(sample_id, ".gff3"))
  return(gff_file)
}

#' Run Panaroo pangenome analysis
#' @param gff_files Vector of paths to GFF files
#' @param out_dir Output directory for pangenome results
#' @param threads Number of threads to use
#' @return Path to gene presence/absence matrix
run_panaroo <- function(gff_files, out_dir, threads = 8) {
  gff_list <- paste(gff_files, collapse = ",")
  out_prefix <- file.path(out_dir, "panaroo")

  cmd <- sprintf("panaroo -i %s -o %s --clean-mode strict -t %d --algos gene",
                 gff_list, basename(out_prefix), threads)

  message("Running Panaroo pangenome analysis...")
  system(cmd)

  # Panaroo output files
  pa_matrix <- file.path(out_dir, "panaroo_gene_presence_absence.CSV")
  return(pa_matrix)
}

#' Convert Panaroo output to binary matrix format
#' @param panaroo_file Path to Panaroo gene presence/absence CSV
#' @param group_file Optional path to sample group/phenotype file
#' @param out_file Output path for formatted matrix
#' @return Data frame of gene presence/absence matrix
convert_panaroo_output <- function(panaroo_file, group_file = NULL, out_file) {
  message("Converting Panaroo output to binary matrix...")

  # Read Panaroo CSV (tab-separated with CSV-like quoting)
  pa_data <- read.csv(panaroo_file, sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)

  # Extract gene names and sample columns
  # Panaroo output: Gene, Annotation, ... (sample columns)
  gene_names <- pa_data[, 1]

  # Sample columns start from column 9 (after metadata columns)
  sample_cols <- 9:ncol(pa_data)
  sample_names <- colnames(pa_data)[sample_cols]

  # Create binary matrix (0/1 for absence/presence)
  binary_matrix <- pa_data[, sample_cols]

  # Convert to binary: any non-empty value = 1, empty = 0
  binary_matrix <- apply(binary_matrix, c(1, 2), function(x) {
    if (is.na(x) || x == "" || x == "-") return(0)
    return(1)
  })
  binary_matrix <- as.data.frame(binary_matrix)
  rownames(binary_matrix) <- gene_names
  colnames(binary_matrix) <- sample_names

  # Transpose to have samples as rows and genes as columns
  binary_matrix <- t(binary_matrix)
  binary_matrix <- as.data.frame(binary_matrix)

  # Add sample ID column
  binary_matrix <- cbind(ID = rownames(binary_matrix), binary_matrix)

  # Merge with group file if provided
  if (!is.null(group_file) && file.exists(group_file)) {
    group_data <- read.table(group_file, sep = "\t", header = TRUE,
                             stringsAsFactors = FALSE, check.names = FALSE)
    binary_matrix <- merge(group_data, binary_matrix, by = "ID", all.x = TRUE)
  }

  # Write output
  write.table(binary_matrix, out_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE, na = "")

  message(sprintf("Gene matrix saved to: %s", out_file))
  message(sprintf("Dimensions: %d samples x %d genes",
                  nrow(binary_matrix), ncol(binary_matrix) - 1))

  return(binary_matrix)
}

# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

main <- function() {
  # Create output directory
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }

  # Get list of FASTA files
  fasta_files <- list.files(opt$input, pattern = "\\.fna$|\\.fasta$|\\.fa$",
                            full.names = TRUE, ignore.case = TRUE)

  if (length(fasta_files) == 0) {
    stop(sprintf("No FASTA files found in %s", opt$input))
  }

  message(sprintf("Found %d genome assemblies", length(fasta_files)))

  # Step 1: Annotation (if not skipped)
  gff_files <- c()

  if (!opt$skip_annotation && !opt$pangenome_only) {
    message("=== Step 1: Genome Annotation ===")

    if (opt$annotator == "bakta") {
      if (!check_command("bakta")) {
        stop("Bakta not found. Please install Bakta or use --annotator prokka")
      }

      for (fasta in fasta_files) {
        gff <- run_bakta(fasta, opt$outdir, opt$threads)
        gff_files <- c(gff_files, gff)
      }
    } else if (opt$annotator == "prokka") {
      if (!check_command("prokka")) {
        stop("Prokka not found. Please install Prokka or use --annotator bakta")
      }

      for (fasta in fasta_files) {
        sample_id <- tools::file_path_sans_ext(basename(fasta))
        sample_out <- file.path(opt$outdir, sample_id)

        cmd <- sprintf("prokka --outdir %s --cpus %d --prefix %s %s",
                       sample_out, opt$threads, sample_id, fasta)
        message(sprintf("[%s] Running Prokka annotation...", sample_id))
        system(cmd)

        gff <- file.path(sample_out, paste0(sample_id, ".gff"))
        gff_files <- c(gff_files, gff)
      }
    }
  } else {
    # Look for existing GFF files
    gff_files <- list.files(opt$input, pattern = "\\.gff$|\\.gff3$",
                            full.names = TRUE, ignore.case = TRUE)
    if (length(gff_files) == 0) {
      stop("No GFF files found for pangenome analysis")
    }
    message(sprintf("Using %d existing GFF files", length(gff_files)))
  }

  # Step 2: Pangenome analysis
  message("=== Step 2: Pangenome Analysis ===")

  if (!check_command("panaroo")) {
    stop("Panaroo not found. Please install Panaroo (pip install panaroo)")
  }

  pa_matrix_file <- run_panaroo(gff_files, opt$outdir, opt$threads)

  # Step 3: Convert to binary matrix
  message("=== Step 3: Generate Binary Matrix ===")

  output_file <- file.path(opt$outdir, "gene_matrix.txt")
  gene_matrix <- convert_panaroo_output(pa_matrix_file, NULL, output_file)

  message("=== Gene Feature Generation Complete ===")
  message(sprintf("Output: %s", output_file))
}

# Run main function
main()
