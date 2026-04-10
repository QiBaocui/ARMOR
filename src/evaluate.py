#!/usr/bin/env python3
"""
Evaluation Module for ARMOR Pipeline

Comprehensive model evaluation and comparison.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class Evaluator:
    """Evaluate and compare AMR prediction models."""

    def __init__(self, metrics: List[str] = None):
        """
        Initialize evaluator.

        Args:
            metrics: List of metrics to calculate
        """
        if metrics is None:
            metrics = ["AUC", "Sensitivity", "Specificity", "Accuracy", "F1"]
        self.metrics = metrics

    def collect_results(self, results_dir: str) -> pd.DataFrame:
        """
        Collect performance metrics from multiple model results.

        Args:
            results_dir: Directory containing model results

        Returns:
            Combined performance DataFrame
        """
        results_path = Path(results_dir)
        all_metrics = []

        # Find all cutoff performance files
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
            raise ValueError("No performance files found!")

        return pd.concat(all_metrics, ignore_index=True)

    def generate_summary(self, perf_df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate summary statistics.

        Args:
            perf_df: Combined performance DataFrame

        Returns:
            Summary statistics DataFrame
        """
        summary = perf_df.groupby(["FeatureType", "Algorithm"]).agg({
            "auc": ["mean", "std"],
            "tpr": ["mean", "std"],  # Sensitivity
            "tnr": ["mean", "std"],  # Specificity
        }).round(4)

        # Flatten column names
        summary.columns = ["_".join(col).strip() for col in summary.columns]
        summary = summary.reset_index()

        return summary

    def compare_feature_types(self, perf_df: pd.DataFrame) -> Dict:
        """
        Compare performance across feature types.

        Args:
            perf_df: Combined performance DataFrame

        Returns:
            Dictionary with comparison statistics
        """
        from scipy import stats

        feature_types = perf_df["FeatureType"].unique()
        comparisons = {}

        for i, ft1 in enumerate(feature_types):
            for ft2 in feature_types[i+1:]:
                auc1 = perf_df[perf_df["FeatureType"] == ft1]["auc"]
                auc2 = perf_df[perf_df["FeatureType"] == ft2]["auc"]

                # T-test
                t_stat, p_value = stats.ttest_ind(auc1, auc2, equal_var=False)

                key = f"{ft1}_vs_{ft2}"
                comparisons[key] = {
                    "mean_1": auc1.mean(),
                    "mean_2": auc2.mean(),
                    "t_stat": t_stat,
                    "p_value": p_value,
                }

        return comparisons

    def save_results(
        self,
        perf_df: pd.DataFrame,
        summary_df: pd.DataFrame,
        comparisons: Dict,
        output_dir: str,
    ) -> None:
        """Save evaluation results."""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Save raw data
        perf_df.to_csv(output_path / "all_metrics_raw.tsv", sep="\t", index=False)

        # Save summary
        summary_df.to_csv(output_path / "model_comparison_summary.tsv", sep="\t", index=False)

        # Save comparisons
        comp_df = pd.DataFrame(comparisons).T
        comp_df.to_csv(output_path / "feature_type_comparisons.tsv", sep="\t")

        logger.info(f"Evaluation results saved to {output_dir}")


def main():
    """Main entry point for evaluation."""
    import click

    @click.command()
    @click.option("--input", "-i", required=True, help="Results directory")
    @click.option("--output", "-o", required=True, help="Output directory")
    def cli(input, output):
        """Evaluate model performance."""
        evaluator = Evaluator()

        # Collect results
        perf_df = evaluator.collect_results(input)
        logger.info(f"Collected {len(perf_df)} performance records")

        # Generate summary
        summary_df = evaluator.generate_summary(perf_df)

        # Compare feature types
        comparisons = evaluator.compare_feature_types(perf_df)

        # Save results
        evaluator.save_results(perf_df, summary_df, comparisons, output)

    cli()


if __name__ == "__main__":
    main()
