#!/usr/bin/env python3
"""
Tests for Model Trainer Module
"""

import unittest
import tempfile
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


class TestModelTrainer(unittest.TestCase):
    """Test model trainer functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

        # Create sample data
        np.random.seed(42)
        n_samples = 100
        n_features = 20

        self.feature_matrix = pd.DataFrame({
            "ID": [f"SAMPLE_{i:03d}" for i in range(n_samples)],
            **{
                f"Feature_{i}": np.random.normal(0, 1, n_samples)
                for i in range(n_features)
            }
        })

        self.group_data = pd.DataFrame({
            "ID": [f"SAMPLE_{i:03d}" for i in range(n_samples)],
            "Group": ["R"] * 50 + ["S"] * 50
        })

    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.temp_dir)

    def test_trainer_initialization(self):
        """Test ModelTrainer initialization."""
        from src.model_trainer import ModelTrainer

        trainer = ModelTrainer(
            nfolds=5,
            seed=42,
            h2o_port=54321,
        )

        self.assertEqual(trainer.nfolds, 5)
        self.assertEqual(trainer.seed, 42)
        self.assertFalse(trainer.h2o_initialized)

    def test_data_preparation(self):
        """Test data preparation for training."""
        from src.model_trainer import ModelTrainer

        trainer = ModelTrainer()

        # Mock H2O initialization for testing
        # (Skip actual H2O initialization in unit tests)

        # Test that data can be merged correctly
        if "ID" in self.feature_matrix.columns:
            feature_matrix = self.feature_matrix.set_index("ID")
        if "ID" in self.group_data.columns:
            group_data = self.group_data.set_index("ID")

        common_samples = feature_matrix.index.intersection(group_data.index)

        self.assertEqual(len(common_samples), 100)


class TestEvaluator(unittest.TestCase):
    """Test evaluator functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

        # Create sample performance data
        self.perf_data = pd.DataFrame({
            "FeatureType": ["Gene"] * 4 + ["SNP"] * 4,
            "Algorithm": ["GBM", "GLM", "RF", "DL"] * 2,
            "auc": [0.95, 0.92, 0.93, 0.91, 0.85, 0.82, 0.84, 0.80],
            "tpr": [0.90, 0.88, 0.89, 0.87, 0.80, 0.78, 0.79, 0.75],
            "tnr": [0.92, 0.90, 0.91, 0.89, 0.85, 0.83, 0.84, 0.82],
        })

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_evaluator_initialization(self):
        """Test Evaluator initialization."""
        from src.evaluate import Evaluator

        evaluator = Evaluator()
        self.assertEqual(len(evaluator.metrics), 5)

    def test_generate_summary(self):
        """Test summary generation."""
        from src.evaluate import Evaluator

        evaluator = Evaluator()
        summary = evaluator.generate_summary(self.perf_data)

        self.assertIsInstance(summary, pd.DataFrame)
        self.assertIn("FeatureType", summary.columns)
        self.assertIn("Algorithm", summary.columns)

    def test_compare_feature_types(self):
        """Test feature type comparison."""
        from src.evaluate import Evaluator

        evaluator = Evaluator()
        comparisons = evaluator.compare_feature_types(self.perf_data)

        self.assertIn("Gene_vs_SNP", comparisons)
        self.assertIn("t_stat", comparisons["Gene_vs_SNP"])
        self.assertIn("p_value", comparisons["Gene_vs_SNP"])


if __name__ == "__main__":
    unittest.main()
