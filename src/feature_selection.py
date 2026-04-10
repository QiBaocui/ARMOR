#!/usr/bin/env python3
"""
Feature Selection Module for ARMOR Pipeline

Implements between-group difference-based feature selection strategy
as described in the manuscript.
"""

import logging
from typing import Dict, List, Optional, Tuple
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

logger = logging.getLogger(__name__)


class FeatureSelector:
    """
    Feature selection for AMR prediction.

    Uses between-group difference strategy:
    - Calculate presence proportion in Resistant (R) and Susceptible (S) groups
    - Discriminative score = |PR - PS|
    - Rank features by score and select top N
    """

    def __init__(
        self,
        top_n: int = 1000,
        p_threshold: float = 0.05,
        or_threshold: float = 2.0,
        min_presence: float = 0.01,
        max_presence: float = 0.99,
    ):
        """
        Initialize feature selector.

        Args:
            top_n: Number of top features to select
            p_threshold: P-value threshold for chi-square test
            or_threshold: Odds ratio threshold
            min_presence: Minimum presence rate to retain feature
            max_presence: Maximum presence rate to retain feature
        """
        self.top_n = top_n
        self.p_threshold = p_threshold
        self.or_threshold = or_threshold
        self.min_presence = min_presence
        self.max_presence = max_presence

    def select_features(
        self,
        feature_matrix: pd.DataFrame,
        group_data: pd.DataFrame,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Perform feature selection.

        Args:
            feature_matrix: Feature matrix (samples x features)
            group_data: Group/phenotype data (samples x 1)

        Returns:
            Tuple of (selected features matrix, feature statistics)
        """
        logger.info("Starting feature selection...")

        # Merge data
        if "ID" in feature_matrix.columns:
            feature_matrix = feature_matrix.set_index("ID")
        if "ID" in group_data.columns:
            group_data = group_data.set_index("ID")

        # Align samples
        common_samples = feature_matrix.index.intersection(group_data.index)
        logger.info(f"Common samples: {len(common_samples)}")

        X = feature_matrix.loc[common_samples]
        y = group_data.loc[common_samples]

        # Get group column
        if "Group" in y.columns:
            y = y["Group"]
        else:
            y = y.iloc[:, 0]

        # Calculate scores for all features
        stats_results = []

        for feature in tqdm(X.columns, desc="Calculating feature scores"):
            result = self._calculate_feature_score(X[feature].values, y.values)
            result["Feature"] = feature
            stats_results.append(result)

        stats_df = pd.DataFrame(stats_results)

        # Sort by score
        stats_df = stats_df.sort_values("Score", ascending=False)
        stats_df["Rank"] = range(1, len(stats_df) + 1)

        # Filter features
        stats_df["Passed_Filter"] = (
            (stats_df["Overall_Presence"] >= self.min_presence) &
            (stats_df["Overall_Presence"] <= self.max_presence)
        )

        stats_df["Significant"] = (
            (stats_df["P_value"] < self.p_threshold) &
            (stats_df["Odds_Ratio"] > self.or_threshold)
        )

        # Select top N features
        filtered = stats_df[stats_df["Passed_Filter"]]
        selected_features = filtered.head(self.top_n)["Feature"].tolist()

        logger.info(f"Features passing filters: {filtered.shape[0]}")
        logger.info(f"Significant features: {stats_df['Significant'].sum()}")
        logger.info(f"Selected top {len(selected_features)} features")

        # Extract selected feature matrix
        selected_matrix = X[selected_features].copy()
        selected_matrix = selected_matrix.reset_index()

        return selected_matrix, stats_df

    def _calculate_feature_score(
        self,
        feature_vector: np.ndarray,
        group_vector: np.ndarray,
    ) -> Dict:
        """
        Calculate discriminative score for a single feature.

        Args:
            feature_vector: Feature values (binary or continuous)
            group_vector: Group labels (R/S)

        Returns:
            Dictionary with score statistics
        """
        # Convert to binary
        feature_binary = (feature_vector > 0).astype(int)

        # Handle missing values
        valid_mask = ~np.isnan(feature_binary) & ~pd.isna(group_vector)
        feature_clean = feature_binary[valid_mask]
        group_clean = group_vector[valid_mask]

        if len(group_clean) < 10:
            return {
                "PR": np.nan,
                "PS": np.nan,
                "Score": np.nan,
                "P_value": np.nan,
                "Odds_Ratio": np.nan,
                "Overall_Presence": np.nan,
            }

        # Identify groups
        groups = np.unique(group_clean)
        if len(groups) != 2:
            return {
                "PR": np.nan,
                "PS": np.nan,
                "Score": np.nan,
                "P_value": np.nan,
                "Odds_Ratio": np.nan,
                "Overall_Presence": np.nan,
            }

        # Assume "R" or first group is resistant
        if "R" in groups:
            resistant_label = "R"
        else:
            resistant_label = groups[0]

        r_mask = group_clean == resistant_label
        s_mask = group_clean != resistant_label

        # Calculate presence proportions
        PR = feature_clean[r_mask].sum() / r_mask.sum()
        PS = feature_clean[s_mask].sum() / s_mask.sum()

        # Overall presence
        overall_presence = feature_clean.sum() / len(feature_clean)

        # Discriminative score
        score = abs(PR - PS)

        # Chi-square test
        contingency = pd.crosstab(feature_clean, group_clean)

        try:
            if contingency.shape == (2, 2):
                chi2, p_value, _, _ = stats.chi2_contingency(contingency, correction=True)

                # Calculate odds ratio
                a, b = contingency.iloc[1, 1], contingency.iloc[1, 0]
                c, d = contingency.iloc[0, 1], contingency.iloc[0, 0]

                # Add pseudo-counts if needed
                if b == 0 or c == 0:
                    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5

                odds_ratio = (a * d) / (b * c)
            else:
                chi2, p_value, _, _ = stats.chi2_contingency(contingency, correction=True)
                odds_ratio = np.nan
        except Exception:
            p_value = np.nan
            odds_ratio = np.nan

        return {
            "PR": PR,
            "PS": PS,
            "Score": score,
            "P_value": p_value,
            "Odds_Ratio": odds_ratio,
            "Overall_Presence": overall_presence,
        }

    def save_results(
        self,
        stats_df: pd.DataFrame,
        output_dir: str,
    ) -> None:
        """Save feature selection results."""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Save all statistics
        stats_df.to_csv(output_path / "feature_selection_all_stats.tsv", sep="\t", index=False)

        # Save selected features
        selected = stats_df.head(self.top_n)
        selected["Feature"].to_csv(
            output_path / "selected_features.txt",
            sep="\t",
            index=False,
            header=False,
        )

        # Save summary
        summary_lines = [
            "=== Feature Selection Summary ===",
            f"Total features: {len(stats_df)}",
            f"Passed presence filter: {(stats_df['Overall_Presence'] >= self.min_presence).sum()}",
            f"Statistically significant: {stats_df['Significant'].sum()}",
            f"Top N selected: {self.top_n}",
            f"P-value threshold: {self.p_threshold}",
            f"Odds ratio threshold: {self.or_threshold}",
            "",
            "Top 10 features by score:",
        ]

        for i, row in selected.head(10).iterrows():
            summary_lines.append(f"  {row['Feature']}: Score={row['Score']:.4f}, OR={row['Odds_Ratio']:.2f}")

        with open(output_path / "feature_selection_summary.txt", "w") as f:
            f.write("\n".join(summary_lines))

        logger.info(f"Results saved to {output_dir}")


def main():
    """Main entry point for feature selection."""
    import click

    @click.command()
    @click.option("--input", "-i", required=True, help="Feature matrix file")
    @click.option("--group", "-g", required=True, help="Group/phenotype file")
    @click.option("--output", "-o", required=True, help="Output directory")
    @click.option("--top-n", "-n", default=1000, help="Number of features to select")
    @click.option("--p-threshold", default=0.05, help="P-value threshold")
    @click.option("--or-threshold", default=2.0, help="Odds ratio threshold")
    def cli(input, group, output, top_n, p_threshold, or_threshold):
        """Perform feature selection."""
        # Load data
        feature_matrix = pd.read_csv(input, sep="\t")
        group_data = pd.read_csv(group, sep="\t")

        # Select features
        selector = FeatureSelector(
            top_n=top_n,
            p_threshold=p_threshold,
            or_threshold=or_threshold,
        )

        selected_matrix, stats_df = selector.select_features(feature_matrix, group_data)
        selector.save_results(stats_df, output)

        # Save selected matrix
        selected_matrix.to_csv(Path(output) / "selected_feature_matrix.tsv", sep="\t", index=False)

    cli()


if __name__ == "__main__":
    main()
