#!/usr/bin/env python3
"""
Prediction Module for ARMOR Pipeline

Apply trained models to predict resistance in new samples.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class Predictor:
    """Predict antibiotic resistance using trained models."""

    def __init__(
        self,
        model_path: str,
        h2o_port: int = 54321,
    ):
        """
        Initialize predictor.

        Args:
            model_path: Path to trained H2O model
            h2o_port: H2O cluster port
        """
        self.model_path = Path(model_path)
        self.h2o_port = h2o_port
        self.model = None
        self.h2o_initialized = False

    def load_model(self) -> None:
        """Load trained H2O model."""
        import h2o

        if not self.h2o_initialized:
            h2o.init(ip="localhost", port=self.h2o_port)
            self.h2o_initialized = True

        # Find model file
        model_files = list(self.model_path.glob("*.zip"))

        if len(model_files) == 0:
            # Try subdirectories
            for subdir in self.model_path.iterdir():
                if subdir.is_dir():
                    model_files = list(subdir.glob("*.zip"))
                    if len(model_files) > 0:
                        break

        if len(model_files) == 0:
            raise FileNotFoundError(f"No model file found in {self.model_path}")

        logger.info(f"Loading model from {model_files[0]}")
        self.model = h2o.load_model(str(model_files[0]))

    def predict(
        self,
        feature_matrix: pd.DataFrame,
        threshold: float = None,
    ) -> pd.DataFrame:
        """
        Predict resistance for new samples.

        Args:
            feature_matrix: Feature matrix for new samples
            threshold: Custom prediction threshold (None for default)

        Returns:
            DataFrame with predictions and probabilities
        """
        import h2o

        if self.model is None:
            raise RuntimeError("Model not loaded. Call load_model() first.")

        # Prepare data
        if "ID" in feature_matrix.columns:
            sample_ids = feature_matrix["ID"]
            feature_matrix = feature_matrix.drop(columns=["ID"])
        else:
            sample_ids = feature_matrix.index

        # Convert to H2O frame
        data_h2o = h2o.H2OFrame(feature_matrix)

        # Generate predictions
        pred = self.model.predict(data_h2o).as_data_frame()

        # Create output DataFrame
        results = pd.DataFrame({
            "SampleID": sample_ids,
            "Prob_S": pred["p0"].values,
            "Prob_R": pred["p1"].values,
            "Prediction": pred["predict"].values,
        })

        # Apply custom threshold if provided
        if threshold is not None:
            results["Prediction"] = np.where(
                results["Prob_R"] >= threshold, "R", "S"
            )

        return results

    def predict_batch(
        self,
        input_files: List[str],
        output_dir: str,
        threshold: float = None,
    ) -> None:
        """
        Predict resistance for multiple sample files.

        Args:
            input_files: List of feature matrix files
            output_dir: Output directory for predictions
            threshold: Custom prediction threshold
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        for input_file in input_files:
            logger.info(f"Processing {input_file}")

            # Load feature matrix
            feature_matrix = pd.read_csv(input_file, sep="\t")

            # Predict
            results = self.predict(feature_matrix, threshold)

            # Save results
            sample_name = Path(input_file).stem
            output_file = output_path / f"{sample_name}_predictions.tsv"
            results.to_csv(output_file, sep="\t", index=False)

            # Save summary
            summary = results["Prediction"].value_counts()
            summary_file = output_path / f"{sample_name}_summary.tsv"
            summary.to_csv(summary_file, sep="\t", header=["Count"])

            logger.info(f"Predictions saved to {output_file}")

    def close(self) -> None:
        """Close H2O connection."""
        if self.h2o_initialized:
            import h2o
            h2o.shutdown(prompt=False)
            self.h2o_initialized = False


def main():
    """Main entry point for prediction."""
    import click

    @click.command()
    @click.option("--input", "-i", required=True, help="Feature matrix file")
    @click.option("--model", "-m", required=True, help="Model directory")
    @click.option("--output", "-o", required=True, help="Output directory")
    @click.option("--threshold", "-t", type=float, default=None, help="Prediction threshold")
    @click.option("--port", "-p", default=54321, help="H2O port")
    def cli(input, model, output, threshold, port):
        """Predict resistance in new samples."""
        # Load feature matrix
        feature_matrix = pd.read_csv(input, sep="\t")

        # Initialize predictor
        predictor = Predictor(model, h2o_port=port)

        try:
            # Load model
            predictor.load_model()

            # Predict
            results = predictor.predict(feature_matrix, threshold)

            # Save results
            output_path = Path(output)
            output_path.mkdir(parents=True, exist_ok=True)

            results.to_csv(output_path / "predictions.tsv", sep="\t", index=False)

            # Save summary
            summary = results["Prediction"].value_counts()
            summary.to_csv(output_path / "summary.tsv", sep="\t", header=["Count"])

            logger.info(f"Predictions saved to {output_path}")

        finally:
            predictor.close()

    cli()


if __name__ == "__main__":
    main()
