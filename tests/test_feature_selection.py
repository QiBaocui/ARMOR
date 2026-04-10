#!/usr/bin/env python3
"""
Tests for Feature Selection Module
"""

import unittest
import tempfile
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from src.feature_selection import FeatureSelector


class TestFeatureSelector(unittest.TestCase):
    """Test feature selection functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

        # Create sample data
        np.random.seed(42)
        n_samples = 100
        n_features = 50

        # Create feature matrix
        self.feature_matrix = pd.DataFrame({
            "ID": [f"SAMPLE_{i:03d}" for i in range(n_samples)],
            **{
                f"Feature_{i}": np.random.binomial(1, 0.3, n_samples)
                for i in range(n_features)
            }
        })

        # Create group data (balanced)
        self.group_data = pd.DataFrame({
            "ID": [f"SAMPLE_{i:03d}" for i in range(n_samples)],
            "Group": ["R"] * (n_samples // 2) + ["S"] * (n_samples // 2)
        })

    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.temp_dir)

    def test_initialization(self):
        """Test FeatureSelector initialization."""
        selector = FeatureSelector(
            top_n=100,
            p_threshold=0.05,
            or_threshold=2.0,
        )

        self.assertEqual(selector.top_n, 100)
        self.assertEqual(selector.p_threshold, 0.05)
        self.assertEqual(selector.or_threshold, 2.0)

    def test_calculate_feature_score(self):
        """Test feature score calculation."""
        selector = FeatureSelector()

        # Create a feature strongly associated with resistance
        feature_vector = np.array([1] * 40 + [0] * 10)  # Mostly in R
        group_vector = np.array(["R"] * 50 + ["S"] * 50)

        result = selector._calculate_feature_score(feature_vector, group_vector)

        self.assertIn("PR", result)
        self.assertIn("PS", result)
        self.assertIn("Score", result)
        self.assertGreater(result["PR"], result["PS"])  # Should be higher in R
        self.assertGreater(result["Score"], 0)

    def test_select_features(self):
        """Test feature selection."""
        selector = FeatureSelector(top_n=10)

        selected_matrix, stats_df = selector.select_features(
            self.feature_matrix, self.group_data
        )

        # Check output types
        self.assertIsInstance(selected_matrix, pd.DataFrame)
        self.assertIsInstance(stats_df, pd.DataFrame)

        # Check number of selected features
        n_features = len([c for c in selected_matrix.columns if c != "ID"])
        self.assertLessEqual(n_features, 10)

        # Check stats columns
        expected_cols = ["Feature", "PR", "PS", "Score", "P_value", "Odds_Ratio", "Rank"]
        for col in expected_cols:
            self.assertIn(col, stats_df.columns)

    def test_save_results(self):
        """Test saving feature selection results."""
        selector = FeatureSelector(top_n=10)

        _, stats_df = selector.select_features(
            self.feature_matrix, self.group_data
        )

        selector.save_results(stats_df, self.temp_dir)

        # Check output files exist
        self.assertTrue((Path(self.temp_dir) / "feature_selection_all_stats.tsv").exists())
        self.assertTrue((Path(self.temp_dir) / "selected_features.txt").exists())
        self.assertTrue((Path(self.temp_dir) / "feature_selection_summary.txt").exists())


class TestFeatureSelectionIntegration(unittest.TestCase):
    """Integration tests for feature selection."""

    def setUp(self):
        """Set up with known significant features."""
        self.temp_dir = tempfile.mkdtemp()

        np.random.seed(42)
        n_samples = 200

        # Create features with known associations
        # Feature 1: Strongly associated with R
        feature_1 = np.array([1] * 80 + [0] * 20 + [0] * 100)  # 80% in R, 0% in S

        # Feature 2: Not associated
        feature_2 = np.random.binomial(1, 0.5, n_samples)

        # Feature 3: Moderately associated
        feature_3 = np.array([1] * 60 + [0] * 40 + [0.2] * 100)  # Higher in R

        self.feature_matrix = pd.DataFrame({
            "ID": [f"S_{i:03d}" for i in range(n_samples)],
            "Strong_Feature": feature_1,
            "Noise_Feature": feature_2,
            "Moderate_Feature": feature_3,
        })

        self.group_data = pd.DataFrame({
            "ID": [f"S_{i:03d}" for i in range(n_samples)],
            "Group": ["R"] * 100 + ["S"] * 100
        })

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_strong_feature_ranked_high(self):
        """Test that strongly associated features are ranked high."""
        selector = FeatureSelector(top_n=10)

        _, stats_df = selector.select_features(
            self.feature_matrix, self.group_data
        )

        # Strong feature should be ranked first
        strong_row = stats_df[stats_df["Feature"] == "Strong_Feature"]
        self.assertEqual(strong_row["Rank"].values[0], 1)

        # Should have high score
        self.assertGreater(strong_row["Score"].values[0], 0.5)


if __name__ == "__main__":
    unittest.main()
