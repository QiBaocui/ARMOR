#!/usr/bin/env python3
"""
Visualization Module for ARMOR Pipeline

Generate publication-quality figures for AMR prediction results.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)

# Set style
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)


class Visualizer:
    """Generate visualizations for AMR prediction results."""

    def __init__(
        self,
        output_dir: str,
        dpi: int = 300,
        theme: str = "bw",
    ):
        """
        Initialize visualizer.

        Args:
            output_dir: Output directory for figures
            dpi: Resolution for raster images
            theme: Plot theme
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.dpi = dpi
        self.theme = theme

    def generate_roc_curves(
        self,
        results_dir: str,
        output_file: str = None,
    ) -> None:
        """
        Generate ROC curve comparison plot.

        Args:
            results_dir: Directory containing ROC data files
            output_file: Output file path
        """
        results_path = Path(results_dir)
        roc_data_list = []

        # Collect ROC data
        for roc_file in results_path.glob("*_rocData.tsv"):
            parts = roc_file.stem.split("_")
            if len(parts) >= 2:
                feature_type = parts[0]
                algorithm = parts[1]

                data = pd.read_csv(roc_file, sep="\t")

                if "fpr" in data.columns and "tpr" in data.columns:
                    data["FeatureType"] = feature_type
                    data["Algorithm"] = algorithm
                    roc_data_list.append(data)

        if len(roc_data_list) == 0:
            logger.warning("No ROC data found")
            return

        roc_df = pd.concat(roc_data_list, ignore_index=True)

        # Create plot
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()

        feature_types = roc_df["FeatureType"].unique()

        for idx, feature_type in enumerate(feature_types[:4]):
            ax = axes[idx]
            subset = roc_df[roc_df["FeatureType"] == feature_type]

            for algo in subset["Algorithm"].unique():
                algo_data = subset[subset["Algorithm"] == algo]
                auc = algo_data["label"].iloc[0] if "label" in algo_data.columns else algo
                ax.plot(algo_data["fpr"], algo_data["tpr"], label=algo, linewidth=2)

            ax.plot([0, 1], [0, 1], "k--", linewidth=1, alpha=0.5)
            ax.set_xlabel("False Positive Rate")
            ax.set_ylabel("True Positive Rate")
            ax.set_title(f"{feature_type} Features")
            ax.legend(loc="lower right", fontsize=8)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)

        plt.tight_layout()

        if output_file is None:
            output_file = self.output_dir / "roc_curves.png"
        else:
            output_file = Path(output_file)

        plt.savefig(output_file, dpi=self.dpi, bbox_inches="tight")
        plt.close()

        logger.info(f"ROC curves saved to {output_file}")

    def generate_performance_comparison(
        self,
        results_dir: str,
        output_file: str = None,
    ) -> None:
        """
        Generate performance metrics comparison plot.

        Args:
            results_dir: Directory containing performance files
            output_file: Output file path
        """
        results_path = Path(results_dir)
        all_metrics = []

        for perf_file in results_path.glob("*_cutoffpermance.tsv"):
            parts = perf_file.stem.split("_")
            if len(parts) >= 2:
                feature_type = parts[0]
                algorithm = parts[1]

                perf = pd.read_csv(perf_file, sep="\t")

                if len(perf) > 0:
                    perf["FeatureType"] = feature_type
                    perf["Algorithm"] = algorithm
                    all_metrics.append(perf)

        if len(all_metrics) == 0:
            logger.warning("No performance data found")
            return

        perf_df = pd.concat(all_metrics, ignore_index=True)

        # Metrics to plot
        metrics_map = {
            "auc": "AUC",
            "tpr": "Sensitivity",
            "tnr": "Specificity",
        }

        fig, axes = plt.subplots(1, 3, figsize=(14, 5))

        for idx, (metric, label) in enumerate(metrics_map.items()):
            if metric in perf_df.columns:
                ax = axes[idx]

                # Prepare data for plotting
                plot_data = perf_df.pivot_table(
                    values=metric,
                    index="Algorithm",
                    columns="FeatureType",
                )

                # Create heatmap
                sns.heatmap(
                    plot_data,
                    annot=True,
                    fmt=".3f",
                    cmap="YlOrRd",
                    vmin=0.5,
                    vmax=1.0,
                    ax=ax,
                    cbar_kws={"label": label},
                )

                ax.set_title(label)
                ax.set_xlabel("Feature Type")
                ax.set_ylabel("Algorithm")

        plt.tight_layout()

        if output_file is None:
            output_file = self.output_dir / "performance_comparison.png"
        else:
            output_file = Path(output_file)

        plt.savefig(output_file, dpi=self.dpi, bbox_inches="tight")
        plt.close()

        logger.info(f"Performance comparison saved to {output_file}")

    def generate_feature_importance(
        self,
        results_dir: str,
        top_n: int = 40,
        output_file: str = None,
    ) -> None:
        """
        Generate feature importance heatmap.

        Args:
            results_dir: Directory containing variable importance files
            top_n: Number of top features to show
            output_file: Output file path
        """
        results_path = Path(results_dir)
        varimp_data = []

        for varimp_file in results_path.glob("*_varimp.tsv"):
            parts = varimp_file.stem.split("_")
            if len(parts) >= 2:
                feature_type = parts[0]
                algorithm = parts[1]

                if feature_type == "Gene":  # Only gene features for this plot
                    data = pd.read_csv(varimp_file, sep="\t")

                    if "variable" in data.columns:
                        # Get importance column
                        imp_col = "scaled_importance" if "scaled_importance" in data.columns else "importance"
                        if imp_col in data.columns:
                            data["FeatureType"] = feature_type
                            data["Algorithm"] = algorithm
                            data["Importance"] = data[imp_col]
                            varimp_data.append(data[["variable", "Algorithm", "Importance"]])

        if len(varimp_data) == 0:
            logger.warning("No variable importance data found")
            return

        varimp_df = pd.concat(varimp_data, ignore_index=True)

        # Get top genes
        top_genes = (
            varimp_df.groupby("variable")["Importance"]
            .mean()
            .nlargest(top_n)
            .index.tolist()
        )

        # Filter data
        plot_data = varimp_df[varimp_df["variable"].isin(top_genes)]

        # Create pivot table
        pivot = plot_data.pivot_table(
            values="Importance",
            index="variable",
            columns="Algorithm",
            aggfunc="mean",
        )

        # Create heatmap
        plt.figure(figsize=(10, 12))
        sns.heatmap(
            pivot,
            cmap="Blues",
            cbar_kws={"label": "Importance"},
            linewidths=0.5,
        )

        plt.title(f"Top {top_n} Resistance-Associated Genes")
        plt.xlabel("Algorithm")
        plt.ylabel("Gene")

        plt.tight_layout()

        if output_file is None:
            output_file = self.output_dir / "feature_importance.png"
        else:
            output_file = Path(output_file)

        plt.savefig(output_file, dpi=self.dpi, bbox_inches="tight")
        plt.close()

        logger.info(f"Feature importance saved to {output_file}")

    def generate_all_figures(self, metrics_file: str) -> None:
        """
        Generate all standard figures.

        Args:
            metrics_file: Path to combined metrics file
        """
        metrics = pd.read_csv(metrics_file)

        # Figure 1: AUC comparison
        self.plot_auc_comparison(metrics)

        # Figure 2: Performance by feature type
        self.plot_feature_type_performance(metrics)

        # Figure 3: Algorithm comparison
        self.plot_algorithm_comparison(metrics)

        logger.info("All figures generated")

    def plot_auc_comparison(self, metrics: pd.DataFrame) -> None:
        """Plot AUC comparison across feature types."""
        plt.figure(figsize=(10, 6))

        sns.barplot(
            data=metrics,
            x="FeatureType",
            y="auc",
            hue="Algorithm",
            palette="Set2",
        )

        plt.xlabel("Feature Type")
        plt.ylabel("AUC")
        plt.title("Model Performance Comparison: AUC by Feature Type")
        plt.ylim(0.5, 1.0)
        plt.legend(title="Algorithm", bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()
        plt.savefig(self.output_dir / "auc_comparison.png", dpi=self.dpi)
        plt.close()

    def plot_feature_type_performance(self, metrics: pd.DataFrame) -> None:
        """Plot performance distribution by feature type."""
        plt.figure(figsize=(8, 6))

        sns.violinplot(
            data=metrics,
            x="FeatureType",
            y="auc",
            fill=True,
            palette="Set2",
        )

        sns.boxplot(
            data=metrics,
            x="FeatureType",
            y="auc",
            width=0.1,
            showcaps=True,
            boxprops={"facecolor": "white", "alpha": 0.5},
            showfliers=False,
        )

        plt.xlabel("Feature Type")
        plt.ylabel("AUC")
        plt.title("AUC Distribution by Feature Type")

        plt.tight_layout()
        plt.savefig(self.output_dir / "auc_by_feature_type.png", dpi=self.dpi)
        plt.close()

    def plot_algorithm_comparison(self, metrics: pd.DataFrame) -> None:
        """Plot performance distribution by algorithm."""
        plt.figure(figsize=(8, 6))

        sns.violinplot(
            data=metrics,
            x="Algorithm",
            y="auc",
            fill=True,
            palette="Set1",
        )

        plt.xlabel("Algorithm")
        plt.ylabel("AUC")
        plt.title("AUC Distribution by Algorithm")
        plt.xticks(rotation=45)

        plt.tight_layout()
        plt.savefig(self.output_dir / "auc_by_algorithm.png", dpi=self.dpi)
        plt.close()


def main():
    """Main entry point for visualization."""
    import click

    @click.command()
    @click.option("--input", "-i", required=True, help="Results directory")
    @click.option("--output", "-o", required=True, help="Output directory")
    @click.option("--dpi", default=300, help="Output DPI")
    def cli(input, output, dpi):
        """Generate visualizations."""
        visualizer = Visualizer(output, dpi=dpi)

        # Generate figures
        visualizer.generate_roc_curves(input)
        visualizer.generate_performance_comparison(input)
        visualizer.generate_feature_importance(input)

    cli()


if __name__ == "__main__":
    main()
