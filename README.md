# ARMOR

A comprehensive pipeline for predicting antibiotic resistance phenotypes in bacteria using whole-genome sequencing (WGS) data and machine learning.

## 

## Features

- **Multi-Feature Integration**: Combines gene presence/absence, SNPs, and k-mers for comprehensive genomic analysis
- **Multiple ML Algorithms**: Supports GBM, GLM, Random Forest, Deep Learning, and Stacked Ensemble
- **Automated Pipeline**: From raw WGS data to resistance predictions with minimal manual intervention
- **Publication-Quality Visualizations**: ROC curves, performance comparisons, and feature importance heatmaps
- **Statistical Analysis**: Chi-square tests, odds ratios, and significance testing for feature associations
- **Pre-trained Models**: Save and load models for predicting new samples

## Installation

```bash
# Clone the repository
git clone https://github.com/your-org/ARMOR.git
cd ARMOR

# Create conda environment
conda env create -f environment.yml
conda activate armor

# Install Python dependencies
pip install -r requirements.txt

# Install R dependencies
Rscript -e "install.packages(c('data.table', 'dplyr', 'ggplot2', 'optparse', 'caret', 'tidyr', 'pROC', 'cowplot'))"
Rscript -e "install.packages('h2o', repos='http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R')"
```

## Quick Start

### 1. Generate Feature Matrices

```bash
# Generate all feature types from raw data
bash scripts/generate_features.sh \
    --genomes genomes/ \
    --fastq fastq/ \
    --reference reference/GCF_000005845.2.fasta \
    --annotation reference/GCF_000005845.2.gff \
    --output features/ \
    --threads 16
```

### 2. Train Models

```bash
# Full analysis with all feature types
Rscript scripts/basic_analyse.R \
    -o results \
    -g group_data.txt \
    -G features/gene_matrix.txt \
    -s features/snp/SNP_matrix.txt \
    -k features/kmer/kmer_matrix.txt \
    -t 0.6 \
    --nfolds 10 \
    --seed 42
```

### 3. Generate Visualizations

```bash
# Generate publication-quality figures
Rscript scripts/06_visualization.R \
    -i results \
    -o figures \
    --dpi 300
```

### 4. Predict New Samples

```bash
# Use trained model to predict new samples
Rscript scripts/09_predict_new_samples.R \
    -i new_samples_matrix.txt \
    -m results/Gene_stacked_model/ \
    -o predictions/
```

## Complete Pipeline

```bash
# Run the complete pipeline with one command
bash scripts/run_pipeline.sh \
    -g group_data.txt \
    -G gene_matrix.txt \
    -s SNP_matrix.txt \
    -k kmer_matrix.txt \
    -o results \
    --nfolds 10 \
    --seed 42
```

## Project Structure

```
ARMOR/
├── README.md
├── LICENSE
├── pyproject.toml
├── environment.yml
├── requirements.txt
├── Makefile
├── example_config.txt           # Configuration examples
├── config/
│   └── pipeline_config.yaml     # Pipeline configuration
├── src/
│   ├── __init__.py
│   ├── feature_generation.py    # Feature generation (Gene/SNP/K-mer)
│   ├── feature_selection.py     # Feature selection
│   ├── model_trainer.py         # Model training with H2O
│   ├── evaluate.py              # Model evaluation
│   ├── visualize.py             # Visualization generation
│   ├── predict.py               # Prediction for new samples
│   └── utils/
│       ├── __init__.py
│       ├── io.py                # I/O utilities
│       └── metrics.py           # Evaluation metrics
├── scripts/
│   ├── basic_analyse.R          # Main R analysis pipeline (entry point)
│   ├── model.R                  # R modeling functions
│   ├── run_pipeline.sh          # Complete pipeline runner
│   ├── generate_features.sh     # Feature generation runner
│   ├── 01_gene_feature_generation.R
│   ├── 02_snp_feature_generation.R
│   ├── 03_kmer_feature_generation.R
│   ├── 04_feature_selection.R
│   ├── 05_model_evaluation.R
│   ├── 06_visualization.R
│   ├── 07_statistical_analysis.R
│   ├── 08_compare_models.R
│   └── 09_predict_new_samples.R
├── docs/
│   └── USAGE.md                 # Detailed usage guide
├── data/
│   ├── reference/               # Reference genomes
│   └── knowledge/               # Knowledge base documents
├── results/
│   ├── figures/                 # Output figures
│   ├── tables/                  # Output tables
│   └── metrics/                 # Evaluation metrics
└── tests/
    ├── test_feature_selection.py
    └── test_model_trainer.py
```

## ARMOR Framework Overview

ARMOR consists of four major components:

1. **Feature Generation**
   
   - Gene presence/absence from pangenome analysis (Bakta + Panaroo)
   - Core genome SNPs (Snippy + Gubbins)
   - K-mer frequencies (Jellyfish/KMC)

2. **Feature Selection**
   
   - Between-group difference strategy: Score = |PR - PS|
   - Statistical filtering: Chi-square p < 0.05, OR > 2
   - Top N feature selection (default: 1000)

3. **Model Training**
   
   - Multiple algorithms: GBM, GLM, RF, DL
   - Cross-validation (default: 10-fold)
   - Stacked ensemble for optimal performance

4. **Evaluation & Visualization**
   
   - AUC, sensitivity, specificity, accuracy, F1 score
   - ROC curves and performance heatmaps
   - Feature importance analysis

## Feature Types

| Feature Type | Method           | Dimensions            | Description               |
| ------------ | ---------------- | --------------------- | ------------------------- |
| Gene         | Bakta + Panaroo  | ~3000-5000 genes      | Presence/absence of genes |
| SNP          | Snippy + Gubbins | ~10000-50000 SNPs     | Core genome polymorphisms |
| K-mer        | Jellyfish/KMC    | ~100000-500000 k-mers | 11-mer frequency profiles |

## Machine Learning Algorithms

| Algorithm | H2O Function            | Description                            |
| --------- | ----------------------- | -------------------------------------- |
| GBM       | `h2o.gbm()`             | Gradient Boosting Machine              |
| GLM       | `h2o.glm()`             | Generalized Linear Model (Elastic Net) |
| RF        | `h2o.randomForest()`    | Random Forest                          |
| DL        | `h2o.deeplearning()`    | Deep Neural Network                    |
| Stacked   | `h2o.stackedEnsemble()` | Ensemble of all above                  |

## Evaluation Metrics

### Read-level (Classification)

- **Precision**: TP / (TP + FP)
- **Recall (Sensitivity)**: TP / (TP + FN)
- **F1-score**: Harmonic mean of Precision and Recall

### Profile-level (Prediction)

- **AUC**: Area under ROC curve
- **Specificity**: TN / (TN + FP)
- **Accuracy**: (TP + TN) / (TP + TN + FP + FN)

## Key Findings

Based on the manuscript results:

1. **Gene-based models** show superior performance (AUC: 0.8936-0.9787)
2. **Integrated multi-feature models** maintain high performance (AUC: 0.8888-0.9879)
3. **Aminoglycosides (GEN, TOB) and Ciprofloxacin (CIP)** achieve near-ceiling performance
4. **Beta-lactams** show heterogeneity, with AMX/CLA having lowest AUC
5. **40 core resistance genes** identified with high concordance to known mechanisms

## 
