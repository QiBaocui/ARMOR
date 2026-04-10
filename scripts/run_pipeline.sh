#!/bin/bash

# ==============================================================================
# ARMOR Pipeline Runner
# ==============================================================================
# Usage:
#   bash scripts/run_pipeline.sh [options]
#
# Options:
#   -g, --group       Group/phenotype file (required)
#   -G, --gene        Gene matrix file (required)
#   -s, --snp         SNP matrix file
#   -k, --kmer        K-mer matrix file
#   -o, --output      Output directory (default: results)
#   -t, --threads     Number of threads (default: 8)
#   -n, --nfolds      Cross-validation folds (default: 10)
#   --seed            Random seed (default: 42)
#   --algorithms      ML algorithms (default: GBM,GLM,RF,DL)
#   --gene-only       Use gene features only
#   --dry-run         Show commands without executing
# ==============================================================================

set -e

# Default values
OUTPUT_DIR="results"
THREADS=8
NFOLDS=10
SEED=42
ALGORITHMS="GBM,GLM,RF,DL"
GENE_ONLY=false
DRY_RUN=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Print usage
usage() {
    cat << EOF
ARMOR Pipeline - Antibiotic Resistance Prediction

Usage: $0 -g group.txt -G gene_matrix.txt [options]

Required:
  -g, --group       Group/phenotype file
  -G, --gene        Gene matrix file

Optional:
  -s, --snp         SNP matrix file
  -k, --kmer        K-mer matrix file
  -o, --output      Output directory (default: results)
  -t, --threads     Number of threads (default: 8)
  -n, --nfolds      Cross-validation folds (default: 10)
  --seed            Random seed (default: 42)
  --algorithms      ML algorithms (default: GBM,GLM,RF,DL)
  --gene-only       Use gene features only
  --dry-run         Show commands without executing
  -h, --help        Show this help message

Example:
  $0 -g group_data.txt -G gene_matrix.txt -s SNP_matrix.txt -k kmer_matrix.txt -o results

EOF
    exit 1
}

# Parse arguments
GROUP_FILE=""
GENE_FILE=""
SNP_FILE=""
KMER_FILE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -g|--group)
            GROUP_FILE="$2"
            shift 2
            ;;
        -G|--gene)
            GENE_FILE="$2"
            shift 2
            ;;
        -s|--snp)
            SNP_FILE="$2"
            shift 2
            ;;
        -k|--kmer)
            KMER_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -n|--nfolds)
            NFOLDS="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
            shift 2
            ;;
        --algorithms)
            ALGORITHMS="$2"
            shift 2
            ;;
        --gene-only)
            GENE_ONLY=true
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

# Validate required arguments
if [[ -z "$GROUP_FILE" ]] || [[ -z "$GENE_FILE" ]]; then
    log_error "Group file (-g) and Gene matrix file (-G) are required!"
    usage
fi

# Check if files exist
log_info "Checking input files..."
for file in "$GROUP_FILE" "$GENE_FILE"; do
    if [[ ! -f "$file" ]]; then
        log_error "File not found: $file"
        exit 1
    fi
done

if [[ -n "$SNP_FILE" ]] && [[ ! -f "$SNP_FILE" ]]; then
    log_error "SNP matrix file not found: $SNP_FILE"
    exit 1
fi

if [[ -n "$KMER_FILE" ]] && [[ ! -f "$KMER_FILE" ]]; then
    log_error "K-mer matrix file not found: $KMER_FILE"
    exit 1
fi

# Create output directory
log_info "Creating output directory: $OUTPUT_DIR"
if [[ "$DRY_RUN" == false ]]; then
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR/logs"
    mkdir -p "$OUTPUT_DIR/figures"
fi

# Build command
CMD="Rscript basic_analyse.R"
CMD+=" -o $OUTPUT_DIR"
CMD+=" -g $GROUP_FILE"
CMD+=" -G $GENE_FILE"

if [[ -n "$SNP_FILE" ]]; then
    CMD+=" -s $SNP_FILE"
fi

if [[ -n "$KMER_FILE" ]]; then
    CMD+=" -k $KMER_FILE"
fi

if [[ "$GENE_ONLY" == true ]]; then
    CMD+=" --gene-only"
fi

CMD+=" -t 0.6"
CMD+=" --nfolds $NFOLDS"
CMD+=" --seed $SEED"
CMD+=" --algorithms $ALGORITHMS"
CMD+=" --nthreads $THREADS"

# Print command
log_info "Running command:"
echo "  $CMD"

if [[ "$DRY_RUN" == true ]]; then
    log_info "Dry run - not executing"
    exit 0
fi

# Run main analysis
log_info "Starting ARMOR pipeline..."
START_TIME=$(date +%s)

$CMD 2>&1 | tee "$OUTPUT_DIR/logs/pipeline.log"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_success "Pipeline completed in $((DURATION / 60)) minutes $((DURATION % 60)) seconds"

# Generate visualizations
if [[ -d "$OUTPUT_DIR" ]]; then
    log_info "Generating visualizations..."

    Rscript scripts/06_visualization.R \
        -i "$OUTPUT_DIR" \
        -o "$OUTPUT_DIR/figures" \
        --dpi 300 \
        2>&1 | tee -a "$OUTPUT_DIR/logs/pipeline.log"

    log_success "Visualizations saved to $OUTPUT_DIR/figures"
fi

# Print summary
log_info "=== Pipeline Summary ==="
log_info "Output directory: $OUTPUT_DIR"
log_info "Results: $OUTPUT_DIR/*_cutoffpermance.tsv"
log_info "Figures: $OUTPUT_DIR/figures/"
log_info "Log file: $OUTPUT_DIR/logs/pipeline.log"

log_success "ARMOR pipeline completed successfully!"
