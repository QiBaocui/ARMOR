#!/usr/bin/env python3
"""
I/O utilities for ARMOR Pipeline
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def load_feature_matrix(
    file_path: str,
    feature_type: str = None,
) -> pd.DataFrame:
    """
    Load feature matrix from file.

    Args:
        file_path: Path to feature matrix file
        feature_type: Type of features (gene/snp/kmer)

    Returns:
        DataFrame with samples as rows and features as columns
    """
    logger.info(f"Loading feature matrix: {file_path}")

    # Detect separator
    if file_path.endswith(".csv"):
        sep = ","
    else:
        sep = "\t"

    df = pd.read_csv(file_path, sep=sep)

    # Validate
    if df.shape[1] < 2:
        raise ValueError("Feature matrix must have at least ID and one feature column")

    logger.info(f"Loaded {df.shape[0]} samples with {df.shape[1] - 1} features")

    return df


def load_group_data(file_path: str) -> pd.DataFrame:
    """
    Load group/phenotype data from file.

    Args:
        file_path: Path to group file

    Returns:
        DataFrame with sample IDs and group labels
    """
    logger.info(f"Loading group data: {file_path}")

    # Detect separator
    if file_path.endswith(".csv"):
        sep = ","
    else:
        sep = "\t"

    df = pd.read_csv(file_path, sep=sep)

    # Validate columns
    if "ID" not in df.columns:
        # Assume first column is ID
        df.columns = ["ID"] + list(df.columns[1:])

    if "Group" not in df.columns:
        # Assume second column is Group
        if df.shape[1] >= 2:
            df.columns = ["ID", "Group"] + list(df.columns[2:])
        else:
            raise ValueError("Group file must have at least ID and Group columns")

    logger.info(f"Loaded {df.shape[0]} samples")
    logger.info(f"Group distribution: {df['Group'].value_counts().to_dict()}")

    return df


def save_results(
    data: pd.DataFrame,
    output_path: str,
    sep: str = "\t",
    index: bool = False,
) -> None:
    """
    Save results to file.

    Args:
        data: DataFrame to save
        output_path: Output file path
        sep: Field separator
        index: Whether to include index
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    data.to_csv(output, sep=sep, index=index)
    logger.info(f"Results saved to {output}")


def save_json(
    data: Dict,
    output_path: str,
    indent: int = 2,
) -> None:
    """
    Save results to JSON file.

    Args:
        data: Dictionary to save
        output_path: Output file path
        indent: JSON indentation
    """
    import json

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(output, "w") as f:
        json.dump(data, f, indent=indent)

    logger.info(f"JSON saved to {output}")


def load_config(config_path: str) -> Dict:
    """
    Load configuration from YAML file.

    Args:
        config_path: Path to configuration file

    Returns:
        Configuration dictionary
    """
    import yaml

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    logger.info(f"Configuration loaded from {config_path}")

    return config


def ensure_directories(dirs: List[str]) -> None:
    """
    Ensure directories exist.

    Args:
        dirs: List of directory paths
    """
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)
        logger.debug(f"Directory ensured: {d}")
