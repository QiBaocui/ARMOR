"""
ARMOR: Antibiotic Resistance prediction using Machine learning with Omics data

A comprehensive pipeline for predicting antibiotic resistance phenotypes
in Escherichia coli using whole-genome sequencing (WGS) data.
"""

__version__ = "1.0.0"
__author__ = "ARMOR Team"
__email__ = "armor@example.com"

from .feature_generation import FeatureGenerator
from .feature_selection import FeatureSelector
from .model_trainer import ModelTrainer
from .evaluate import Evaluator
from .visualize import Visualizer
from .predict import Predictor

__all__ = [
    "FeatureGenerator",
    "FeatureSelector",
    "ModelTrainer",
    "Evaluator",
    "Visualizer",
    "Predictor",
]
