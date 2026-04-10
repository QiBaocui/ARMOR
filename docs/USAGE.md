# ARMOR 使用指南

## 目录

1. [安装](#安装)
2. [快速开始](#快速开始)
3. [详细配置](#详细配置)
4. [命令行工具](#命令行工具)
5. [Python API](#python-api)
6. [结果解释](#结果解释)
7. [故障排除](#故障排除)

---

## 安装

### 系统要求

- R >= 4.0.0
- Python >= 3.8
- Linux/macOS 操作系统
- 最低 16GB RAM（推荐 32GB+）
- 100GB 可用磁盘空间

### 安装步骤

#### 1. 克隆仓库

```bash
git clone https://github.com/your-org/ARMOR.git
cd ARMOR
```

#### 2. 创建 Conda 环境（推荐）

```bash
conda env create -f environment.yml
conda activate armor
```

#### 3. 或使用 pip 安装

```bash
pip install -r requirements.txt
```

#### 4. 安装 R 依赖

```r
# 安装 CRAN 包
install.packages(c("data.table", "dplyr", "ggplot2", "optparse",
                   "caret", "tidyr", "broom", "pROC", "cowplot"))

# 安装 H2O
install.packages("h2o", repos = "http://h2o-release.s3.amazonaws.com/h2o/latest_stable_R")
```

#### 5. 安装生物信息学工具

```bash
# 使用 conda 安装
conda install -c bioconda bakta panaroo snippy gubbins jellyfish kmc
```

---

## 快速开始

### 完整分析流程

```bash
# 1. 生成特征矩阵
bash scripts/generate_features.sh \
    --genomes genomes/ \
    --fastq fastq/ \
    --reference reference/ \
    --output features/

# 2. 训练模型
Rscript scripts/basic_analyse.R \
    -o results \
    -g group_data.txt \
    -G gene_matrix.txt \
    -s SNP_matrix.txt \
    -k kmer_matrix.txt

# 3. 生成可视化
Rscript scripts/06_visualization.R \
    -i results \
    -o figures
```

### 分步执行

#### 1. 生成基因特征

```bash
Rscript scripts/01_gene_feature_generation.R \
    -i genomes/ \
    -o features/gene \
    -t 16
```

#### 2. 生成 SNP 特征

```bash
Rscript scripts/02_snp_feature_generation.R \
    -i fastq/ \
    -r reference/GCF_000005845.2.fasta \
    -a reference/GCF_000005845.2.gff \
    -o features/snp \
    -t 16
```

#### 3. 生成 K-mer 特征

```bash
Rscript scripts/03_kmer_feature_generation.R \
    -i fastq/ \
    -o features/kmer \
    -k 11 \
    -t 16
```

#### 4. 特征选择

```bash
Rscript scripts/04_feature_selection.R \
    -i features/gene/gene_matrix.txt \
    -g group_data.txt \
    -o features/selected \
    -n 1000
```

#### 5. 训练模型

```bash
Rscript scripts/basic_analyse.R \
    -o results \
    -g group_data.txt \
    -G features/gene/gene_matrix.txt \
    -s features/snp/SNP_matrix.txt \
    -k features/kmer/kmer_matrix.txt \
    -t 0.6 \
    --nfolds 10 \
    --algorithms GBM,GLM,RF,DL
```

#### 6. 评估模型

```bash
Rscript scripts/05_model_evaluation.R \
    -i results \
    -o evaluation
```

#### 7. 生成图表

```bash
Rscript scripts/06_visualization.R \
    -i results \
    -o figures \
    --dpi 300
```

#### 8. 统计分析

```bash
Rscript scripts/07_statistical_analysis.R \
    -i gene_matrix.txt \
    -g group_data.txt \
    -o stats
```

---

## 详细配置

### 配置文件 (config/pipeline_config.yaml)

```yaml
# 数据路径
data:
  group_file: group_data.txt
  gene_matrix: gene_matrix.txt
  snp_matrix: SNP_matrix.txt
  kmer_matrix: kmer_matrix.txt

# 特征生成
feature_generation:
  threads: 8
  kmer_size: 11
  min_presence: 0.01
  max_presence: 0.99

# 特征选择
feature_selection:
  top_n: 1000
  p_threshold: 0.05
  or_threshold: 2.0

# 模型训练
training:
  train_ratio: 0.6
  split_type: proportional
  seed: 42
  nfolds: 10

  algorithms:
    - GBM
    - GLM
    - RF
    - DL
```

### 样本量计算

对于不平衡数据集，建议使用 `sample_size` 分割方法：

```bash
Rscript scripts/basic_analyse.R \
    -o results \
    -g group_data.txt \
    -G gene_matrix.txt \
    -s SNP_matrix.txt \
    -k kmer_matrix.txt \
    -S sample_size \
    -t 200  # 每组 200 个样本
```

---

## 命令行工具

### armor-generate

特征生成工具。

```bash
# 基因特征
armor-generate gene \
    --input genomes/ \
    --output features/gene \
    --threads 16

# SNP 特征
armor-generate snp \
    --input fastq/ \
    --reference reference.fasta \
    --output features/snp \
    --threads 16

# K-mer 特征
armor-generate kmer \
    --input fastq/ \
    --output features/kmer \
    --kmer-size 11 \
    --threads 16
```

### armor-select

特征选择工具。

```bash
armor-select \
    --input gene_matrix.txt \
    --group group_data.txt \
    --output selected_features/ \
    --top-n 1000 \
    --p-threshold 0.05 \
    --or-threshold 2.0
```

### armor-train

模型训练工具。

```bash
armor-train \
    --input gene_matrix.txt \
    --group group_data.txt \
    --output results/ \
    --feature-type Gene \
    --algorithms GBM,RF,DL \
    --nfolds 10 \
    --seed 42
```

### armor-evaluate

模型评估工具。

```bash
armor-evaluate \
    --input results/ \
    --output evaluation/
```

### armor-visualize

可视化工具。

```bash
armor-visualize \
    --input results/ \
    --output figures/ \
    --dpi 300
```

### armor-predict

新样本预测工具。

```bash
armor-predict \
    --input new_samples_matrix.txt \
    --model results/Gene_stacked_model/ \
    --output predictions/
```

---

## Python API

### 特征生成

```python
from src.feature_generation import FeatureGenerator

# 初始化
generator = FeatureGenerator(
    reference_fasta="reference.fasta",
    reference_gff="annotation.gff",
    threads=16,
)

# 生成基因矩阵
gene_matrix = generator.generate_gene_matrix(
    genome_dir="genomes/",
    output_dir="features/gene/",
)

# 生成 SNP 矩阵
snp_matrix = generator.generate_snp_matrix(
    fastq_dir="fastq/",
    output_dir="features/snp/",
)

# 生成 K-mer 矩阵
kmer_matrix = generator.generate_kmer_matrix(
    fastq_dir="fastq/",
    output_dir="features/kmer/",
    kmer_size=11,
)
```

### 特征选择

```python
from src.feature_selection import FeatureSelector
import pandas as pd

# 加载数据
feature_matrix = pd.read_csv("gene_matrix.txt", sep="\t")
group_data = pd.read_csv("group_data.txt", sep="\t")

# 初始化选择器
selector = FeatureSelector(
    top_n=1000,
    p_threshold=0.05,
    or_threshold=2.0,
)

# 执行特征选择
selected_matrix, stats = selector.select_features(
    feature_matrix, group_data
)

# 保存结果
selector.save_results(stats, "selected_features/")
```

### 模型训练

```python
from src.model_trainer import ModelTrainer
import pandas as pd

# 初始化训练器
trainer = ModelTrainer(
    nfolds=10,
    seed=42,
    h2o_port=54321,
)

# 初始化 H2O
trainer.init_h2o()

try:
    # 加载数据
    feature_matrix = pd.read_csv("gene_matrix.txt", sep="\t")
    group_data = pd.read_csv("group_data.txt", sep="\t")

    # 准备数据
    train_h2o, test_h2o, features = trainer.prepare_data(
        feature_matrix, group_data
    )

    # 训练模型
    results = trainer.train_models(
        train_h2o, test_h2o, features,
        algorithms=["GBM", "GLM", "RF", "DL"]
    )

    # 保存模型
    trainer.save_models(results["models"], "results/", "Gene")

finally:
    trainer.shutdown_h2o()
```

### 模型评估

```python
from src.evaluate import Evaluator

evaluator = Evaluator()

# 收集结果
perf_df = evaluator.collect_results("results/")

# 生成摘要
summary = evaluator.generate_summary(perf_df)

# 比较特征类型
comparisons = evaluator.compare_feature_types(perf_df)

# 保存结果
evaluator.save_results(perf_df, summary, comparisons, "evaluation/")
```

### 可视化

```python
from src.visualize import Visualizer

visualizer = Visualizer(
    output_dir="figures/",
    dpi=300,
)

# 生成 ROC 曲线
visualizer.generate_roc_curves("results/")

# 生成性能比较
visualizer.generate_performance_comparison("results/")

# 生成特征重要性
visualizer.generate_feature_importance("results/", top_n=40)
```

### 预测新样本

```python
from src.predict import Predictor
import pandas as pd

# 初始化预测器
predictor = Predictor(
    model_path="results/Gene_stacked_model/",
    h2o_port=54321,
)

try:
    # 加载模型
    predictor.load_model()

    # 加载新样本特征
    new_samples = pd.read_csv("new_samples_matrix.txt", sep="\t")

    # 预测
    results = predictor.predict(new_samples)

    # 保存结果
    results.to_csv("predictions.tsv", sep="\t", index=False)

finally:
    predictor.close()
```

---

## 结果解释

### 性能指标

| 指标 | 说明 | 取值范围 |
|------|------|----------|
| AUC | ROC 曲线下面积 | 0-1 |
| Sensitivity (TPR) | 真阳性率 | 0-1 |
| Specificity (TNR) | 真阴性率 | 0-1 |
| Accuracy | 准确率 | 0-1 |
| F1 Score | F1 分数 | 0-1 |
| PPV | 阳性预测值 | 0-1 |
| NPV | 阴性预测值 | 0-1 |

### 显著性标记

- `***`: p < 0.001
- `**`: p < 0.01
- `*`: p < 0.05
- `.`: p < 0.1
- `ns`: 不显著 (p ≥ 0.05)

### 输出文件说明

```
results/
├── Gene_GBM_pred.tsv              # 预测结果
├── Gene_GBM_ROC.png               # ROC 曲线
├── Gene_GBM_cutoffpermance.tsv    # 最佳阈值性能
├── Gene_GBM_varimp.tsv            # 变量重要性
├── Gene_stacked_model/            # 保存的模型
└── ...
```

---

## 故障排除

### 常见问题

#### 1. H2O 初始化失败

**错误**: `H2OConnectionError`

**解决**:
- 检查端口是否被占用：`netstat -tlnp | grep 54321`
- 尝试不同端口：`--port 54322`
- 增加内存：`--min_mem_size 16G --max_mem_size 64G`

#### 2. 内存不足

**错误**: `MemoryError` 或模型训练崩溃

**解决**:
- 减少特征数量（使用特征选择）
- 减少交叉验证折数：`--nfolds 5`
- 增加 H2O 内存限制

#### 3. 特征矩阵过大

**错误**: 无法加载大型特征矩阵

**解决**:
```r
# 使用 data.table 快速读取
library(data.table)
feature_matrix <- fread("large_matrix.txt")

# 或使用 Python
import pandas as pd
feature_matrix = pd.read_csv("large_matrix.txt", sep="\t")
```

#### 4. 生物信息学工具未找到

**错误**: `bakta: command not found`

**解决**:
```bash
# 使用 conda 安装
conda install -c bioconda bakta panaroo snippy
```

### 获取帮助

- GitHub Issues: https://github.com/your-org/ARMOR/issues
- 文档：https://github.com/your-org/ARMOR/docs
