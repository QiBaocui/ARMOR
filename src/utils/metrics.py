#!/usr/bin/env python3
"""
Evaluation metrics for ARMOR Pipeline
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def calculate_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_prob: np.ndarray = None,
) -> Dict:
    """
    Calculate classification metrics.

    Args:
        y_true: True labels
        y_pred: Predicted labels
        y_prob: Predicted probabilities (optional)

    Returns:
        Dictionary with metrics
    """
    from sklearn.metrics import (
        accuracy_score,
        precision_score,
        recall_score,
        f1_score,
        roc_auc_score,
        confusion_matrix,
    )

    # Convert to binary if needed
    if len(np.unique(y_true)) == 2:
        pos_label = np.unique(y_true)[-1]
    else:
        pos_label = 1

    metrics = {
        "accuracy": accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, pos_label=pos_label, zero_division=0),
        "recall": recall_score(y_true, y_pred, pos_label=pos_label, zero_division=0),
        "f1": f1_score(y_true, y_pred, pos_label=pos_label, zero_division=0),
    }

    # Calculate AUC if probabilities provided
    if y_prob is not None:
        try:
            metrics["auc"] = roc_auc_score(y_true, y_prob)
        except Exception:
            metrics["auc"] = np.nan

    # Confusion matrix
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    metrics["true_positive"] = int(tp)
    metrics["true_negative"] = int(tn)
    metrics["false_positive"] = int(fp)
    metrics["false_negative"] = int(fn)

    # Additional metrics
    if tp + fp > 0:
        metrics["ppv"] = tp / (tp + fp)  # Positive predictive value
    else:
        metrics["ppv"] = np.nan

    if tn + fn > 0:
        metrics["npv"] = tn / (tn + fn)  # Negative predictive value
    else:
        metrics["npv"] = np.nan

    if tp + fn > 0:
        metrics["sensitivity"] = tp / (tp + fn)
    else:
        metrics["sensitivity"] = np.nan

    if tn + fp > 0:
        metrics["specificity"] = tn / (tn + fp)
    else:
        metrics["specificity"] = np.nan

    return metrics


def calculate_youden_index(
    sensitivity: np.ndarray,
    specificity: np.ndarray,
) -> np.ndarray:
    """
    Calculate Youden index for each threshold.

    Youden Index = Sensitivity + Specificity - 1

    Args:
        sensitivity: Array of sensitivity values
        specificity: Array of specificity values

    Returns:
        Array of Youden indices
    """
    return sensitivity + specificity - 1


def find_optimal_cutoff(
    thresholds: np.ndarray,
    sensitivity: np.ndarray,
    specificity: np.ndarray,
) -> Tuple[float, float]:
    """
    Find optimal cutoff using Youden index.

    Args:
        thresholds: Array of probability thresholds
        sensitivity: Array of sensitivity values
        specificity: Array of specificity values

    Returns:
        Tuple of (optimal_threshold, youden_index)
    """
    youden = calculate_youden_index(sensitivity, specificity)
    optimal_idx = np.argmax(youden)

    return thresholds[optimal_idx], youden[optimal_idx]


def calculate_roc_metrics(
    y_true: np.ndarray,
    y_prob: np.ndarray,
) -> Dict:
    """
    Calculate ROC-related metrics.

    Args:
        y_true: True labels
        y_prob: Predicted probabilities

    Returns:
        Dictionary with ROC metrics
    """
    from sklearn.metrics import roc_curve, auc

    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    roc_auc = auc(fpr, tpr)

    # Find optimal cutoff
    youden = tpr - fpr
    optimal_idx = np.argmax(youden)
    optimal_cutoff = thresholds[optimal_idx]

    return {
        "fpr": fpr,
        "tpr": tpr,
        "thresholds": thresholds,
        "auc": roc_auc,
        "optimal_cutoff": optimal_cutoff,
        "max_youden": youden[optimal_idx],
    }


def calculate_feature_importance(
    X: pd.DataFrame,
    y: pd.Series,
    method: str = "chi2",
) -> pd.DataFrame:
    """
    Calculate feature importance scores.

    Args:
        X: Feature matrix
        y: Target variable
        method: Method for importance calculation

    Returns:
        DataFrame with feature importance scores
    """
    from scipy import stats

    n_features = X.shape[1]
    importance = pd.DataFrame({
        "Feature": X.columns,
        "Score": np.zeros(n_features),
        "P_value": np.ones(n_features),
    })

    if method == "chi2":
        for i, feature in enumerate(X.columns):
            # Binary conversion
            feature_binary = (X[feature] > 0).astype(int)

            # Contingency table
            contingency = pd.crosstab(feature_binary, y)

            try:
                chi2, p_value, _, _ = stats.chi2_contingency(contingency)
                importance.loc[i, "Score"] = chi2
                importance.loc[i, "P_value"] = p_value
            except Exception:
                importance.loc[i, "Score"] = 0
                importance.loc[i, "P_value"] = 1

    elif method == "mutual_info":
        from sklearn.feature_selection import mutual_info_classif

        scores = mutual_info_classif(X, y)
        importance["Score"] = scores

    # Sort by score
    importance = importance.sort_values("Score", ascending=False)

    return importance


def compare_models(
    metrics_df: pd.DataFrame,
    group_by: str = "FeatureType",
) -> pd.DataFrame:
    """
    Compare model performance across groups.

    Args:
        metrics_df: DataFrame with performance metrics
        group_by: Column to group by

    Returns:
        DataFrame with comparison statistics
    """
    from scipy import stats

    groups = metrics_df[group_by].unique()
    comparisons = []

    for i, g1 in enumerate(groups):
        for g2 in groups[i + 1:]:
            auc1 = metrics_df[metrics_df[group_by] == g1]["auc"]
            auc2 = metrics_df[metrics_df[group_by] == g2]["auc"]

            # T-test
            t_stat, p_value = stats.ttest_ind(auc1, auc2, equal_var=False)

            comparisons.append({
                "Group1": g1,
                "Group2": g2,
                "Mean1": auc1.mean(),
                "Mean2": auc2.mean(),
                "T_statistic": t_stat,
                "P_value": p_value,
                "Significant": "Yes" if p_value < 0.05 else "No",
            })

    return pd.DataFrame(comparisons)
