# ==============================================================================
# ARMOR Pipeline Makefile
# ==============================================================================
# Usage:
#   make help              - Show available targets
#   make all               - Run complete pipeline (requires input files)
#   make features          - Generate all feature types
#   make models            - Train models
#   make figures           - Generate publication figures
#   make clean             - Remove output files
# ==============================================================================

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

# Input files (modify these paths for your data)
GROUP_FILE = group_data.txt
GENE_MATRIX = gene_matrix.txt
SNP_MATRIX = SNP_matrix.txt
KMER_MATRIX = kmer_matrix.txt

# Directories
GENOME_DIR = genomes
FASTQ_DIR = fastq
REFERENCE = reference/GCF_000005845.2.fasta
ANNOTATION = reference/GCF_000005845.2.gff

# Output directories
OUTPUT_DIR = results
FEATURES_DIR = features
FIGURES_DIR = figures

# Parameters
TRAIN_RATIO = 0.6
NFOLDS = 10
SEED = 42
THREADS = 8
KMER_SIZE = 11

# R executable
R = Rscript

# ------------------------------------------------------------------------------
# Targets
# ------------------------------------------------------------------------------

.PHONY: all help clean features models figures stats predict

# Default target
all: features models figures stats

# Help message
help:
	@echo "ARMOR Pipeline - Available Targets:"
	@echo "===================================="
	@echo ""
	@echo "Main Targets:"
	@echo "  make all          - Run complete pipeline"
	@echo "  make features     - Generate all feature types"
	@echo "  make models       - Train machine learning models"
	@echo "  make figures      - Generate publication figures"
	@echo "  make stats        - Run statistical analysis"
	@echo "  make predict      - Predict new samples"
	@echo "  make clean        - Remove output files"
	@echo ""
	@echo "Feature Generation:"
	@echo "  make gene         - Generate gene features"
	@echo "  make snp          - Generate SNP features"
	@echo "  make kmer         - Generate k-mer features"
	@echo ""
	@echo "Configuration:"
	@echo "  Edit Makefile variables to customize paths and parameters"
	@echo ""

# ------------------------------------------------------------------------------
# Feature Generation
# ------------------------------------------------------------------------------

features: gene snp kmer

gene: $(FEATURES_DIR)/gene/gene_matrix.txt
$(FEATURES_DIR)/gene/gene_matrix.txt: $(GENOME_DIR)
	@echo "=== Generating Gene Features ==="
	mkdir -p $(FEATURES_DIR)/gene
	$(R) scripts/01_gene_feature_generation.R \
		-i $(GENOME_DIR) \
		-o $(FEATURES_DIR)/gene \
		-t $(THREADS)

snp: $(FEATURES_DIR)/snp/SNP_matrix.txt
$(FEATURES_DIR)/snp/SNP_matrix.txt: $(FASTQ_DIR) $(REFERENCE) $(ANNOTATION)
	@echo "=== Generating SNP Features ==="
	mkdir -p $(FEATURES_DIR)/snp
	$(R) scripts/02_snp_feature_generation.R \
		-i $(FASTQ_DIR) \
		-r $(REFERENCE) \
		-a $(ANNOTATION) \
		-o $(FEATURES_DIR)/snp \
		-t $(THREADS)

kmer: $(FEATURES_DIR)/kmer/kmer_matrix.txt
$(FEATURES_DIR)/kmer/kmer_matrix.txt: $(FASTQ_DIR)
	@echo "=== Generating K-mer Features ==="
	mkdir -p $(FEATURES_DIR)/kmer
	$(R) scripts/03_kmer_feature_generation.R \
		-i $(FASTQ_DIR) \
		-o $(FEATURES_DIR)/kmer \
		-k $(KMER_SIZE) \
		-t $(THREADS)

# ------------------------------------------------------------------------------
# Model Training
# ------------------------------------------------------------------------------

models: $(OUTPUT_DIR)/train_data.Rdata

$(OUTPUT_DIR)/train_data.Rdata: $(GENE_MATRIX) $(SNP_MATRIX) $(KMER_MATRIX) $(GROUP_FILE)
	@echo "=== Training Machine Learning Models ==="
	mkdir -p $(OUTPUT_DIR)
	$(R) scripts/basic_analyse.R \
		-o $(OUTPUT_DIR) \
		-g $(GROUP_FILE) \
		-G $(GENE_MATRIX) \
		-s $(SNP_MATRIX) \
		-k $(KMER_MATRIX) \
		-t $(TRAIN_RATIO) \
		-S proportional \
		--nfolds $(NFOLDS) \
		--seed $(SEED)

# Quick run (gene-only, fewer algorithms)
models-quick:
	@echo "=== Training Quick Models ==="
	mkdir -p $(OUTPUT_DIR)_quick
	$(R) scripts/basic_analyse.R \
		-o $(OUTPUT_DIR)_quick \
		-g $(GROUP_FILE) \
		-G $(GENE_MATRIX) \
		-s $(SNP_MATRIX) \
		-k $(KMER_MATRIX) \
		--gene-only \
		--algorithms GBM,RF \
		--nfolds 5 \
		--seed $(SEED)

# ------------------------------------------------------------------------------
# Evaluation and Visualization
# ------------------------------------------------------------------------------

figures: $(FIGURES_DIR)/figure1_roc_curves.png

$(FIGURES_DIR)/figure1_roc_curves.png: $(OUTPUT_DIR)
	@echo "=== Generating Figures ==="
	mkdir -p $(FIGURES_DIR)
	$(R) scripts/06_visualization.R \
		-i $(OUTPUT_DIR) \
		-o $(FIGURES_DIR) \
		--dpi 300

stats: $(FIGURES_DIR)/volcano_plot.png
$(FIGURES_DIR)/volcano_plot.png: $(GENE_MATRIX) $(GROUP_FILE)
	@echo "=== Running Statistical Analysis ==="
	mkdir -p $(FIGURES_DIR)/stats
	$(R) scripts/07_statistical_analysis.R \
		-i $(GENE_MATRIX) \
		-g $(GROUP_FILE) \
		-o $(FIGURES_DIR)/stats \
		--top-n 50

evaluate:
	@echo "=== Evaluating Models ==="
	mkdir -p $(OUTPUT_DIR)/evaluation
	$(R) scripts/05_model_evaluation.R \
		-i $(OUTPUT_DIR) \
		-o $(OUTPUT_DIR)/evaluation

# ------------------------------------------------------------------------------
# Prediction
# ------------------------------------------------------------------------------

predict:
	@echo "=== Predicting New Samples ==="
	@echo "Usage: make predict NEW_SAMPLES=<matrix> MODEL=<model_dir> OUT=<outdir>"
	@if [ -n "$(NEW_SAMPLES)" ] && [ -n "$(MODEL)" ] && [ -n "$(OUT)" ]; then \
		$(R) scripts/09_predict_new_samples.R \
			-i $(NEW_SAMPLES) \
			-m $(MODEL) \
			-o $(OUT); \
	else \
		echo "Please provide NEW_SAMPLES, MODEL, and OUT parameters"; \
	fi

# ------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------

clean:
	@echo "=== Cleaning Output Files ==="
	rm -rf $(OUTPUT_DIR)
	rm -rf $(OUTPUT_DIR)_quick
	rm -rf $(FEATURES_DIR)
	rm -rf $(FIGURES_DIR)
	rm -rf *.log *.Rdata
	@echo "Cleanup complete. Input files preserved."

# Clean specific components
clean-features:
	rm -rf $(FEATURES_DIR)

clean-models:
	rm -rf $(OUTPUT_DIR)

clean-figures:
	rm -rf $(FIGURES_DIR)

# ------------------------------------------------------------------------------
# Utility
# ------------------------------------------------------------------------------

# Check dependencies
check-deps:
	@echo "=== Checking Dependencies ==="
	@echo -n "R: "
	@which Rscript > /dev/null && echo "OK" || echo "NOT FOUND"
	@echo -n "H2O: "
	$(R) -e "library(h2o)" 2>/dev/null && echo "OK" || echo "NOT FOUND"
	@echo -n "data.table: "
	$(R) -e "library(data.table)" 2>/dev/null && echo "OK" || echo "NOT FOUND"
	@echo -n "ggplot2: "
	$(R) -e "library(ggplot2)" 2>/dev/null && echo "OK" || echo "NOT FOUND"
	@echo -n "caret: "
	$(R) -e "library(caret)" 2>/dev/null && echo "OK" || echo "NOT FOUND"
	@echo ""
	@echo "External tools:"
	@echo -n "Bakta: "
	@which bakta > /dev/null && echo "OK" || echo "NOT FOUND"
	@echo -n "Panaroo: "
	@which panaroo > /dev/null && echo "OK" || echo "NOT FOUND"
	@echo -n "Snippy: "
	@which snippy > /dev/null && echo "OK" || echo "NOT FOUND"
	@echo -n "Jellyfish: "
	@which jellyfish > /dev/null && echo "OK" || echo "NOT FOUND"

# Show configuration
show-config:
	@echo "=== ARMOR Pipeline Configuration ==="
	@echo "Input Files:"
	@echo "  Group:        $(GROUP_FILE)"
	@echo "  Gene Matrix:  $(GENE_MATRIX)"
	@echo "  SNP Matrix:   $(SNP_MATRIX)"
	@echo "  K-mer Matrix: $(KMER_MATRIX)"
	@echo ""
	@echo "Directories:"
	@echo "  Genomes:      $(GENOME_DIR)"
	@echo "  FastQ:        $(FASTQ_DIR)"
	@echo "  Reference:    $(REFERENCE)"
	@echo ""
	@echo "Output:"
	@echo "  Results:      $(OUTPUT_DIR)"
	@echo "  Features:     $(FEATURES_DIR)"
	@echo "  Figures:      $(FIGURES_DIR)"
	@echo ""
	@echo "Parameters:"
	@echo "  Train Ratio:  $(TRAIN_RATIO)"
	@echo "  CV Folds:     $(NFOLDS)"
	@echo "  Seed:         $(SEED)"
	@echo "  Threads:      $(THREADS)"
	@echo "  K-mer Size:   $(KMER_SIZE)"
