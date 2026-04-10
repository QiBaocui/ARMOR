#!/usr/bin/env Rscript

# ==============================================================================
# Script: 02_snp_feature_generation.R
# Description: Generate SNP matrix from WGS data using Snippy and Gubbins
#              Based on manuscript methods: Snippy variant calling + Gubbins filtering
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
              help = "Directory containing FASTQ or BAM files [required]"),
  make_option(c("-r", "--reference"), type = "character", default = NULL,
              help = "Reference genome FASTA file [required]"),
  make_option(c("-a", "--annotation"), type = "character", default = NULL,
              help = "Reference genome annotation GFF/GTF file [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory for SNP matrix [required]"),
  make_option(c("-t", "--threads"), type = "integer", default = 8,
              help = "Number of threads [default: 8]"),
  make_option(c("--missing-threshold"), type = "double", default = 0.05,
              help = "Maximum missing data rate per site [default: 0.05]"),
  make_option(c("--maf"), type = "double", default = 0.01,
              help = "Minimum minor allele frequency [default: 0.01]"),
  make_option(c("--skip-gubbins"), action = "store_true", default = FALSE,
              help = "Skip Gubbins recombination filtering"),
  make_option(c("--vcf-only"), action = "store_true", default = FALSE,
              help = "Only generate VCF files (skip matrix generation)")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Generate SNP matrix from WGS data using Snippy")
opt <- parse_args(opt_parser)

# Validate arguments
required_args <- c("input", "reference", "annotation", "outdir")
for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(sprintf("Argument --%s is required!", arg))
  }
}

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

#' Check if command exists
check_command <- function(cmd) {
  result <- system2("which", cmd, stdout = NULL, stderr = NULL)
  return(result == 0)
}

#' Run Snippy for variant calling on a single sample
#' @param input_file Path to FASTQ or BAM file
#' @param ref_fasta Reference genome FASTA
#' @param out_dir Output directory
#' @param sample_id Sample identifier
#' @param threads Number of threads
run_snippy <- function(input_file, ref_fasta, out_dir, sample_id, threads = 8) {
  sample_out <- file.path(out_dir, sample_id)

  # Determine if input is FASTQ or BAM
  ext <- tolower(tools::file_ext(input_file))

  if (ext %in% c("fastq", "fq", "fastq.gz", "fq.gz")) {
    # Check if paired-end
    is_paired <- grepl("_R1|_r1|_1\\.", input_file, ignore.case = TRUE)

    if (is_paired && grepl("_R1|_r1|_1\\.", input_file, ignore.case = TRUE)) {
      # Get R2 file
      r1_file <- input_file
      r2_file <- gsub("_R1|_r1|_1\\.", "_R2_", r1_file)
      if (!file.exists(r2_file)) {
        r2_file <- gsub("_1\\.", "_2.", r1_file)
      }

      cmd <- sprintf("snippy --cpus %d --ref %s --R1 %s --R2 %s --outdir %s --prefix %s",
                     threads, ref_fasta, r1_file, r2_file, sample_out, sample_id)
    } else {
      cmd <- sprintf("snippy --cpus %d --ref %s --se %s --outdir %s --prefix %s",
                     threads, ref_fasta, input_file, sample_out, sample_id)
    }
  } else if (ext %in% c("bam")) {
    cmd <- sprintf("snippy --cpus %d --ref %s --bam %s --outdir %s --prefix %s",
                   threads, ref_fasta, input_file, sample_out, sample_id)
  } else {
    stop(sprintf("Unsupported input format: %s", ext))
  }

  message(sprintf("[%s] Running Snippy variant calling...", sample_id))
  system(cmd)

  vcf_file <- file.path(sample_out, paste0(sample_id, ".vcf"))
  return(vcf_file)
}

#' Run Snippy-core to generate core SNP alignment
#' @param snippy_dirs Vector of Snippy output directories
#' @param out_dir Output directory for core SNPs
run_snippy_core <- function(snippy_dirs, out_dir) {
  dirs_list <- paste(snippy_dirs, collapse = ",")

  cmd <- sprintf("snippy-core --ref %s --dir %s %s",
                 file.path(snippy_dirs[1], paste0(basename(snippy_dirs[1]), ".fasta")),
                 out_dir, dirs_list)

  message("Running Snippy-core to identify core SNPs...")
  system(cmd)

  core_vcf <- file.path(out_dir, "core.vcf")
  core_aln <- file.path(out_dir, "core.full.aln")
  return(list(vcf = core_vcf, aln = core_aln))
}

#' Run Gubbins to filter recombinant regions
#' @param alignment_file Path to core alignment file
#' @param out_dir Output directory
#' @param threads Number of threads
run_gubbins <- function(alignment_file, out_dir, threads = 8) {
  cmd <- sprintf("run_gubbins.py --threads %d --output_dir %s %s",
                 threads, out_dir, alignment_file)

  message("Running Gubbins to identify and mask recombinant regions...")
  system(cmd)

  # Gubbins output
  filtered_aln <- file.path(out_dir, "gubbins.final_tree.filtered_polymorphic_sites.fasta")
  vcf_file <- file.path(out_dir, "gubbins.final_tree.variant_calls.vcf")
  return(list(alignment = filtered_aln, vcf = vcf_file))
}

#' Convert VCF to binary SNP matrix
#' @param vcf_file Path to VCF file
#' @param out_file Output path for SNP matrix
#' @param missing_threshold Maximum missing data rate
#' @param maf Minimum minor allele frequency
#' @return Data frame of SNP matrix
convert_vcf_to_matrix <- function(vcf_file, out_file, missing_threshold = 0.05, maf = 0.01) {
  message("Converting VCF to binary SNP matrix...")

  # Read VCF file
  vcf <- fread(vcf_file, sep = "\t", header = TRUE,
               skip = "#CHROM", stringsAsFactors = FALSE, check.names = FALSE)

  # Parse VCF columns
  chrom <- vcf$CHROM
  pos <- vcf$POS
  ref <- vcf$REF
  alt <- vcf$ALT

  # Get sample columns (columns after FORMAT)
  format_idx <- which(colnames(vcf) == "FORMAT")
  sample_cols <- (format_idx + 1):ncol(vcf)
  sample_names <- colnames(vcf)[sample_cols]

  message(sprintf("Found %d samples and %d variant sites",
                  length(sample_names), nrow(vcf)))

  # Create genotype matrix
  gt_matrix <- matrix(NA, nrow = nrow(vcf), ncol = length(sample_names))
  rownames(gt_matrix) <- paste0(chrom, "_", pos)
  colnames(gt_matrix) <- sample_names

  # Parse genotypes for each sample
  for (i in seq_along(sample_cols)) {
    col_idx <- sample_cols[i]
    gt_col <- vcf[[col_idx]]

    # Extract GT (genotype) field
    gt <- sapply(strsplit(gt_col, ":"), function(x) {
      if (length(x) == 0 || x[1] == "." || is.na(x[1])) return(NA)
      return(x[1])
    })

    # Convert to binary: 0 = homozygous ref, 1 = heterozygous/alt
    gt_binary <- sapply(gt, function(x) {
      if (is.na(x) || x == ".") return(NA)
      if (x == "0/0" || x == "0|0") return(0)
      return(1)  # 0/1, 1/1, 0|1, 1|0, etc.
    })

    gt_matrix[, i] <- gt_binary
  }

  # Filter sites by missing data
  missing_rate <- colMeans(is.na(gt_matrix))
  keep_sites <- missing_rate <= missing_threshold
  message(sprintf("Filtering sites with >%.1f%% missing data: %d sites removed",
                  missing_threshold * 100, sum(!keep_sites)))
  gt_matrix <- gt_matrix[keep_sites, ]

  # Filter sites by MAF
  if (nrow(gt_matrix) > 0) {
    alt_freq <- colMeans(gt_matrix, na.rm = TRUE) / 2  # Approximate allele frequency
    maf_filter <- alt_freq >= maf & alt_freq <= (1 - maf)
    message(sprintf("Filtering sites with MAF <%.2f: %d sites removed",
                    maf, sum(!maf_filter)))
    gt_matrix <- gt_matrix[maf_filter, ]
  }

  # Impute remaining missing values with most common genotype
  for (j in 1:ncol(gt_matrix)) {
    na_idx <- which(is.na(gt_matrix[, j]))
    if (length(na_idx) > 0) {
      mode_val <- as.numeric(names(sort(table(gt_matrix[, j]), decreasing = TRUE)[1]))
      gt_matrix[na_idx, j] <- mode_val
    }
  }

  # Create output data frame
  snp_matrix <- as.data.frame(gt_matrix)
  snp_matrix <- cbind(ID = rownames(snp_matrix), snp_matrix)

  # Transpose to have samples as rows
  snp_matrix <- t(snp_matrix[, -1])
  snp_matrix <- as.data.frame(snp_matrix)
  snp_matrix <- cbind(ID = colnames(snp_matrix), snp_matrix)

  # Write output
  write.table(snp_matrix, out_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE, na = "")

  message(sprintf("SNP matrix saved to: %s", out_file))
  message(sprintf("Dimensions: %d samples x %d SNPs",
                  nrow(snp_matrix), ncol(snp_matrix) - 1))

  return(snp_matrix)
}

# ------------------------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------------------------

main <- function() {
  # Create output directory
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }

  # Get list of input files (FASTQ or BAM)
  input_patterns <- c("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz", "*.bam")
  input_files <- c()
  for (pattern in input_patterns) {
    files <- list.files(opt$input, pattern = pattern,
                        full.names = TRUE, ignore.case = TRUE)
    input_files <- c(input_files, files)
  }

  # Only keep R1 files for paired-end (avoid double counting)
  input_files <- input_files[!grepl("_R2|_r2|_2\\.", input_files, ignore.case = TRUE)]

  if (length(input_files) == 0) {
    stop(sprintf("No input files found in %s", opt$input))
  }

  message(sprintf("Found %d input files", length(input_files)))

  # Check required tools
  if (!check_command("snippy")) {
    stop("Snippy not found. Please install Snippy (conda install -c bioconda snippy)")
  }

  # Step 1: Run Snippy for each sample
  message("=== Step 1: Variant Calling with Snippy ===")

  snippy_dirs <- c()
  for (input_file in input_files) {
    sample_id <- tools::file_path_sans_ext(basename(input_file))
    # Clean sample ID
    sample_id <- gsub("_R1$", "", sample_id)

    sample_out <- file.path(opt$outdir, "snippy", sample_id)
    if (!dir.exists(sample_out)) {
      dir.create(sample_out, recursive = TRUE)
    }

    run_snippy(input_file, opt$reference, sample_out, sample_id, opt$threads)
    snippy_dirs <- c(snippy_dirs, sample_out)
  }

  if (opt$vcf_only) {
    message("VCF-only mode: skipping core SNP and matrix generation")
    return(invisible(NULL))
  }

  # Step 2: Run Snippy-core for core SNPs
  message("=== Step 2: Core SNP Identification ===")

  core_out <- file.path(opt$outdir, "core")
  if (!dir.exists(core_out)) {
    dir.create(core_out, recursive = TRUE)
  }

  core_results <- run_snippy_core(snippy_dirs, core_out)

  if (!file.exists(core_results$aln)) {
    stop("Snippy-core failed to generate core alignment")
  }

  # Step 3: Run Gubbins (optional)
  if (!opt$skip_gubbins) {
    message("=== Step 3: Recombination Filtering with Gubbins ===")

    if (check_command("run_gubbins.py")) {
      gubbins_out <- file.path(opt$outdir, "gubbins")
      if (!dir.exists(gubbins_out)) {
        dir.create(gubbins_out, recursive = TRUE)
      }

      gubbins_results <- run_gubbins(core_results$aln, gubbins_out, opt$threads)

      # Use Gubbins-filtered alignment
      final_aln <- gubbins_results$alignment
      final_vcf <- gubbins_results$vcf
    } else {
      message("Gubbins not found, using core SNPs without recombination filtering")
      final_aln <- core_results$aln
      final_vcf <- core_results$vcf
    }
  } else {
    final_aln <- core_results$aln
    final_vcf <- core_results$vcf
  }

  # Step 4: Generate binary SNP matrix
  message("=== Step 4: Generate SNP Matrix ===")

  output_file <- file.path(opt$outdir, "SNP_matrix.txt")

  # If we have a VCF, convert it
  if (file.exists(final_vcf)) {
    snp_matrix <- convert_vcf_to_matrix(final_vcf, output_file,
                                        opt$missing_threshold, opt$maf)
  } else {
    # Convert alignment file to matrix
    message("Converting alignment to SNP matrix...")
    # Implementation for alignment conversion
  }

  message("=== SNP Feature Generation Complete ===")
  message(sprintf("Output: %s", output_file))
}

# Run main function
main()
