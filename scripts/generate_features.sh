#!/bin/bash

# ==============================================================================
# ARMOR Feature Generation Script
# ==============================================================================
# Usage:
#   bash scripts/generate_features.sh [options]
#
# Options:
#   --genomes         Genome FASTA directory
#   --fastq           FASTQ directory
#   --reference       Reference genome FASTA
#   --annotation      Reference annotation GFF
#   --output          Output directory (default: features)
#   --threads         Number of threads (default: 8)
#   --kmer-size       K-mer size (default: 11)
#   --skip-gene       Skip gene feature generation
#   --skip-snp        Skip SNP feature generation
#   --skip-kmer       Skip K-mer feature generation
#   --dry-run         Show commands without executing
# ==============================================================================

set -e

# Default values
OUTPUT_DIR="features"
THREADS=8
KMER_SIZE=11
DRY_RUN=false

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

usage() {
    cat << EOF
ARMOR Feature Generation

Usage: $0 [options]

Required (depending on feature type):
  --genomes         Genome FASTA directory (for gene features)
  --fastq           FASTQ directory (for SNP and K-mer features)
  --reference       Reference genome FASTA (for SNP features)
  --annotation      Reference annotation GFF (for SNP features)

Optional:
  --output          Output directory (default: features)
  --threads         Number of threads (default: 8)
  --kmer-size       K-mer size (default: 11)
  --skip-gene       Skip gene feature generation
  --skip-snp        Skip SNP feature generation
  --skip-kmer       Skip K-mer feature generation
  --dry-run         Show commands without executing
  -h, --help        Show this help

Examples:
  # Generate all features
  $0 --genomes genomes/ --fastq fastq/ --reference ref.fasta --annotation ref.gff

  # Generate only gene features
  $0 --genomes genomes/ --skip-snp --skip-kmer

  # Generate only K-mer features
  $0 --fastq fastq/ --skip-gene --skip-snp --kmer-size 15

EOF
    exit 1
}

# Parse arguments
GENOME_DIR=""
FASTQ_DIR=""
REFERENCE=""
ANNOTATION=""
SKIP_GENE=false
SKIP_SNP=false
SKIP_KMER=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --genomes)
            GENOME_DIR="$2"
            shift 2
            ;;
        --fastq)
            FASTQ_DIR="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --annotation)
            ANNOTATION="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --kmer-size)
            KMER_SIZE="$2"
            shift 2
            ;;
        --skip-gene)
            SKIP_GENE=true
            shift
            ;;
        --skip-snp)
            SKIP_SNP=true
            shift
            ;;
        --skip-kmer)
            SKIP_KMER=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate arguments
if [[ "$SKIP_GENE" == false ]] && [[ -z "$GENOME_DIR" ]]; then
    log_error "Genome directory is required for gene feature generation!"
    usage
fi

if [[ "$SKIP_SNP" == false ]] && [[ -z "$FASTQ_DIR" ]]; then
    log_error "FASTQ directory is required for SNP feature generation!"
    usage
fi

if [[ "$SKIP_SNP" == false ]] && [[ -z "$REFERENCE" ]]; then
    log_error "Reference genome is required for SNP feature generation!"
    usage
fi

if [[ "$SKIP_KMER" == false ]] && [[ -z "$FASTQ_DIR" ]]; then
    log_error "FASTQ directory is required for K-mer feature generation!"
    usage
fi

# Create output directory
log_info "Creating output directory: $OUTPUT_DIR"
if [[ "$DRY_RUN" == false ]]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Check dependencies
check_dependencies() {
    log_info "Checking dependencies..."

    if [[ "$SKIP_GENE" == false ]]; then
        if ! command -v bakta &> /dev/null; then
            log_warning "Bakta not found. Install with: conda install -c bioconda bakta"
        fi
        if ! command -v panaroo &> /dev/null; then
            log_warning "Panaroo not found. Install with: pip install panaroo"
        fi
    fi

    if [[ "$SKIP_SNP" == false ]]; then
        if ! command -v snippy &> /dev/null; then
            log_warning "Snippy not found. Install with: conda install -c bioconda snippy"
        fi
    fi

    if [[ "$SKIP_KMER" == false ]]; then
        if command -v kmc &> /dev/null; then
            log_info "Using KMC for k-mer counting"
        elif command -v jellyfish &> /dev/null; then
            log_info "Using Jellyfish for k-mer counting"
        else
            log_warning "Neither KMC nor Jellyfish found!"
        fi
    fi
}

# Run feature generation
run_gene_generation() {
    log_info "=== Generating Gene Features ==="

    CMD="Rscript scripts/01_gene_feature_generation.R"
    CMD+=" -i $GENOME_DIR"
    CMD+=" -o $OUTPUT_DIR/gene"
    CMD+=" -t $THREADS"

    log_info "Running: $CMD"

    if [[ "$DRY_RUN" == false ]]; then
        mkdir -p "$OUTPUT_DIR/gene"
        $CMD
        log_success "Gene features generated: $OUTPUT_DIR/gene/gene_matrix.txt"
    fi
}

run_snp_generation() {
    log_info "=== Generating SNP Features ==="

    CMD="Rscript scripts/02_snp_feature_generation.R"
    CMD+=" -i $FASTQ_DIR"
    CMD+=" -r $REFERENCE"
    CMD+=" -a $ANNOTATION"
    CMD+=" -o $OUTPUT_DIR/snp"
    CMD+=" -t $THREADS"

    log_info "Running: $CMD"

    if [[ "$DRY_RUN" == false ]]; then
        mkdir -p "$OUTPUT_DIR/snp"
        $CMD
        log_success "SNP features generated: $OUTPUT_DIR/snp/SNP_matrix.txt"
    fi
}

run_kmer_generation() {
    log_info "=== Generating K-mer Features ==="

    CMD="Rscript scripts/03_kmer_feature_generation.R"
    CMD+=" -i $FASTQ_DIR"
    CMD+=" -o $OUTPUT_DIR/kmer"
    CMD+=" -k $KMER_SIZE"
    CMD+=" -t $THREADS"

    log_info "Running: $CMD"

    if [[ "$DRY_RUN" == false ]]; then
        mkdir -p "$OUTPUT_DIR/kmer"
        $CMD
        log_success "K-mer features generated: $OUTPUT_DIR/kmer/kmer_matrix.txt"
    fi
}

# Main execution
check_dependencies

if [[ "$SKIP_GENE" == false ]]; then
    run_gene_generation
fi

if [[ "$SKIP_SNP" == false ]]; then
    run_snp_generation
fi

if [[ "$SKIP_KMER" == false ]]; then
    run_kmer_generation
fi

log_success "Feature generation completed!"
log_info "Output directory: $OUTPUT_DIR"
