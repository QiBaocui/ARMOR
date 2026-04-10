#!/usr/bin/env python3
"""
Model Training Module for ARMOR Pipeline

Train machine learning models using H2O framework:
- GBM (Gradient Boosting Machine)
- GLM (Generalized Linear Model)
- Random Forest
- Deep Learning
- Stacked Ensemble
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class ModelTrainer:
    """Train and evaluate ML models for AMR prediction."""

    def __init__(
        self,
        nfolds: int = 10,
        seed: int = 42,
        h2o_port: int = 54321,
        nthreads: int = -1,
    ):
        """
        Initialize model trainer.

        Args:
            nfolds: Number of cross-validation folds
            seed: Random seed for reproducibility
            h2o_port: H2O cluster port
            nthreads: Number of CPU threads (-1 for all)
        """
        self.nfolds = nfolds
        self.seed = seed
        self.h2o_port = h2o_port
        self.nthreads = nthreads
        self.h2o_initialized = False

    def init_h2o(self) -> None:
        """Initialize H2O cluster."""
        try:
            import h2o
            h2o.init(
                ip="localhost",
                port=self.h2o_port,
                nthreads=self.nthreads,
                min_mem_size="8G",
                max_mem_size="32G",
            )
            self.h2o_initialized = True
            logger.info(f"H2O initialized: {h2o.cluster_status()}")
        except Exception as e:
            logger.error(f"Failed to initialize H2O: {e}")
            raise

    def shutdown_h2o(self) -> None:
        """Shutdown H2O cluster."""
        if self.h2o_initialized:
            import h2o
            h2o.shutdown(prompt=False)
            self.h2o_initialized = False

    def prepare_data(
        self,
        feature_matrix: pd.DataFrame,
        group_data: pd.DataFrame,
    ) -> tuple:
        """
        Prepare data for H2O modeling.

        Args:
            feature_matrix: Feature matrix
            group_data: Group/phenotype data

        Returns:
            Tuple of (train_h2o, test_h2o, feature_names)
        """
        import h2o
        from sklearn.model_selection import train_test_split

        # Merge data
        if "ID" in feature_matrix.columns:
            feature_matrix = feature_matrix.set_index("ID")
        if "ID" in group_data.columns:
            group_data = group_data.set_index("ID")

        # Align samples
        common_samples = feature_matrix.index.intersection(group_data.index)
        X = feature_matrix.loc[common_samples]
        y = group_data.loc[common_samples]

        # Get group labels
        if "Group" in y.columns:
            y = y["Group"]
        else:
            y = y.iloc[:, 0]

        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.4, random_state=self.seed, stratify=y
        )

        logger.info(f"Training samples: {len(X_train)}, Test samples: {len(X_test)}")

        # Convert to H2O frames
        train_df = X_train.copy()
        train_df["Group"] = y_train.values
        test_df = X_test.copy()
        test_df["Group"] = y_test.values

        train_h2o = h2o.H2OFrame(train_df)
        test_h2o = h2o.H2OFrame(test_df)

        # Convert response to factor
        train_h2o["Group"] = train_h2o["Group"].asfactor()
        test_h2o["Group"] = test_h2o["Group"].asfactor()

        feature_names = [c for c in X.columns if c not in ["ID", "Group"]]

        return train_h2o, test_h2o, feature_names

    def train_models(
        self,
        train_h2o,
        test_h2o,
        feature_names: List[str],
        algorithms: List[str] = None,
    ) -> Dict:
        """
        Train multiple models.

        Args:
            train_h2o: Training H2O frame
            test_h2o: Test H2O frame
            feature_names: List of feature column names
            algorithms: List of algorithms to train

        Returns:
            Dictionary of trained models and results
        """
        import h2o

        if algorithms is None:
            algorithms = ["GBM", "GLM", "RF", "DL"]

        y = "Group"
        x = feature_names

        models = {}
        results = {}

        for algo in algorithms:
            logger.info(f"Training {algo} model...")

            if algo == "GBM":
                model = h2o.gbm(
                    x=x, y=y,
                    training_frame=train_h2o,
                    nfolds=self.nfolds,
                    keep_cross_validation_predictions=True,
                    fold_assignment="Stratified",
                    seed=self.seed,
                )
            elif algo == "GLM":
                model = h2o.glm(
                    x=x, y=y,
                    training_frame=train_h2o,
                    nfolds=self.nfolds,
                    keep_cross_validation_predictions=True,
                    fold_assignment="Stratified",
                    seed=self.seed,
                )
            elif algo == "RF":
                model = h2o.randomForest(
                    x=x, y=y,
                    training_frame=train_h2o,
                    nfolds=self.nfolds,
                    keep_cross_validation_predictions=True,
                    fold_assignment="Stratified",
                    seed=self.seed,
                )
            elif algo == "DL":
                model = h2o.deeplearning(
                    x=x, y=y,
                    training_frame=train_h2o,
                    nfolds=self.nfolds,
                    keep_cross_validation_predictions=True,
                    fold_assignment="Stratified",
                    seed=self.seed,
                )
            else:
                logger.warning(f"Unknown algorithm: {algo}")
                continue

            models[algo] = model

            # Evaluate
            perf = model.model_performance(test_h2o)
            auc = perf.auc()

            results[algo] = {
                "AUC": auc,
                "model": model,
            }

            logger.info(f"{algo} AUC: {auc:.4f}")

        # Train stacked ensemble if multiple models
        if len(models) >= 2:
            logger.info("Training Stacked Ensemble...")
            stacked = h2o.stackedEnsemble(
                x=x, y=y,
                training_frame=train_h2o,
                base_models=list(models.values()),
            )
            models["Stacked"] = stacked

            perf = stacked.model_performance(test_h2o)
            results["Stacked"] = {
                "AUC": perf.auc(),
                "model": stacked,
            }
            logger.info(f"Stacked Ensemble AUC: {perf.auc():.4f}")

        return {"models": models, "results": results}

    def get_optimal_cutoff(self, perf) -> float:
        """
        Get optimal probability cutoff using Youden index.

        Youden Index = Sensitivity + Specificity - 1 = TPR - FPR
        """
        # Extract thresholds and metrics
        metrics = perf.thresholds_and_metric_scores

        # Calculate Youden index
        youden = metrics["tpr"] - metrics["fpr"]

        # Find optimal cutoff
        max_idx = np.argmax(youden)
        optimal_cutoff = metrics["threshold"][max_idx]

        return optimal_cutoff

    def save_models(
        self,
        models: Dict,
        output_dir: str,
        feature_type: str,
    ) -> None:
        """Save trained models."""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        for name, model in models.items():
            model_path = output_path / f"{feature_type}_{name}_model"
            model.save_model(str(model_path))
            logger.info(f"Saved model: {model_path}")

    def save_predictions(
        self,
        model,
        test_h2o,
        test_df: pd.DataFrame,
        output_dir: str,
        feature_type: str,
        model_name: str,
    ) -> None:
        """Save predictions and evaluation metrics."""
        import h2o

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Generate predictions
        pred = model.predict(test_h2o).as_data_frame()
        pred_df = test_df[["Group"]].copy()
        pred_df["Prob_S"] = pred["p0"].values
        pred_df["Prob_R"] = pred["p1"].values
        pred_df["Prediction"] = pred["predict"].values

        # Save predictions
        pred_file = output_path / f"{feature_type}_{model_name}_pred.tsv"
        pred_df.to_csv(pred_file, sep="\t")

        # Calculate performance metrics
        perf = model.model_performance(test_h2o)
        auc = perf.auc()

        # Get optimal cutoff
        optimal_cutoff = self.get_optimal_cutoff(perf)

        # Extract metrics at optimal cutoff
        metrics = perf.thresholds_and_metric_scores
        cutoff_idx = np.argmin(np.abs(metrics["threshold"] - optimal_cutoff))

        cutoff_metrics = {
            "threshold": float(metrics["threshold"][cutoff_idx]),
            "tpr": float(metrics["tpr"][cutoff_idx]),
            "fpr": float(metrics["fpr"][cutoff_idx]),
            "tnr": float(metrics["tnr"][cutoff_idx]),
            "ppv": float(metrics["ppv"][cutoff_idx]),
            "npv": float(metrics["npv"][cutoff_idx]),
            "auc": float(auc),
            "model": model_name,
        }

        # Save metrics
        metrics_file = output_path / f"{feature_type}_{model_name}_cutoffpermance.tsv"
        pd.DataFrame([cutoff_metrics]).to_csv(metrics_file, sep="\t", index=False)

        # Save ROC data
        roc_data = pd.DataFrame({
            "tpr": metrics["tpr"],
            "fpr": metrics["fpr"],
            "threshold": metrics["threshold"],
        })
        roc_data = pd.concat([
            pd.DataFrame({"tpr": [0], "fpr": [0], "threshold": [1]}),
            roc_data
        ])
        roc_data = roc_data.sort_values("fpr")
        roc_data["label"] = f"{model_name}_AUC: {auc:.3f}"

        roc_file = output_path / f"{feature_type}_{model_name}_rocData.tsv"
        roc_data.to_csv(roc_file, sep="\t", index=False)

        # Save variable importance
        try:
            if model_name == "Stacked":
                varimp = h2o.varimp(model).as_data_frame()
            else:
                varimp = h2o.varimp(model).as_data_frame()
            varimp_file = output_path / f"{feature_type}_{model_name}_varimp.tsv"
            varimp.to_csv(varimp_file, sep="\t", index=False)
        except Exception as e:
            logger.warning(f"Could not save variable importance: {e}")

        logger.info(f"Saved predictions and metrics for {model_name}")


def main():
    """Main entry point for model training."""
    import click

    @click.command()
    @click.option("--input", "-i", required=True, help="Feature matrix file")
    @click.option("--group", "-g", required=True, help="Group file")
    @click.option("--output", "-o", required=True, help="Output directory")
    @click.option("--feature-type", "-f", default="Gene", help="Feature type")
    @click.option("--algorithms", default="GBM,GLM,RF,DL", help="Algorithms to train")
    @click.option("--nfolds", default=10, help="Cross-validation folds")
    @click.option("--seed", default=42, help="Random seed")
    @click.option("--port", default=54321, help="H2O port")
    def cli(input, group, output, feature_type, algorithms, nfolds, seed, port):
        """Train ML models."""
        # Load data
        feature_matrix = pd.read_csv(input, sep="\t")
        group_data = pd.read_csv(group, sep="\t")

        # Initialize trainer
        trainer = ModelTrainer(
            nfolds=nfolds,
            seed=seed,
            h2o_port=port,
        )

        # Initialize H2O
        trainer.init_h2o()

        try:
            # Prepare data
            train_h2o, test_h2o, feature_names = trainer.prepare_data(
                feature_matrix, group_data
            )

            # Parse algorithms
            algo_list = [a.strip() for a in algorithms.split(",")]

            # Train models
            results = trainer.train_models(
                train_h2o, test_h2o, feature_names, algo_list
            )

            # Save models
            trainer.save_models(results["models"], output, feature_type)

            # Save predictions
            test_df = test_h2o.as_data_frame()
            for name in results["models"]:
                trainer.save_predictions(
                    results["models"][name],
                    test_h2o,
                    test_df,
                    output,
                    feature_type,
                    name,
                )

        finally:
            trainer.shutdown_h2o()

    cli()


if __name__ == "__main__":
    main()
