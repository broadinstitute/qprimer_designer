"""Model loading and inference utilities."""

import warnings
from importlib.resources import files
from pathlib import Path
from typing import Tuple

import joblib
import torch


def get_model_path(filename: str) -> Path:
    """
    Get path to a bundled model file.

    Args:
        filename: Name of the model file (e.g., 'combined_classifier.pth')

    Returns:
        Path to the model file
    """
    return files('qprimer_designer.data').joinpath(filename)


def load_scaler():
    """
    Load the pre-trained StandardScaler for feature normalization.

    Returns:
        sklearn StandardScaler object
    """
    scaler_path = get_model_path('standard_scaler.joblib')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return joblib.load(scaler_path)


def load_classifier(device: str = 'cpu'):
    """
    Load the pre-trained binary classifier model.

    Args:
        device: Device to load the model on ('cpu' or 'cuda')

    Returns:
        CombinedModelClassifier in eval mode
    """
    model_path = get_model_path('combined_classifier.pth')
    model = torch.load(model_path, map_location=device, weights_only=False)
    model.to(device)
    model.eval()
    return model


def load_regressor(device: str = 'cpu'):
    """
    Load the pre-trained regressor model.

    Args:
        device: Device to load the model on ('cpu' or 'cuda')

    Returns:
        CombinedModel in eval mode
    """
    model_path = get_model_path('combined_regressor.pth')
    model = torch.load(model_path, map_location=device, weights_only=False)
    model.to(device)
    model.eval()
    return model


def load_models(device: str = None) -> Tuple:
    """
    Load all models and scaler for inference.

    Args:
        device: Device to load models on. If None, auto-detects CUDA.

    Returns:
        Tuple of (scaler, classifier, regressor, device_str)
    """
    if device is None:
        device = 'cuda' if torch.cuda.is_available() else 'cpu'

    scaler = load_scaler()
    classifier = load_classifier(device)
    regressor = load_regressor(device)

    return scaler, classifier, regressor, device


# Feature columns expected by the models
FEATURE_COLUMNS = [
    "len_f", "Tm_f", "GC_f", "indel_f", "mm_f",
    "len_r", "Tm_r", "GC_r", "indel_r", "mm_r",
    "prod_len", "prod_Tm",
]

# Sequence columns needed for encoding
SEQUENCE_COLUMNS = ["pseq_f", "tseq_f", "pseq_r", "tseq_r"]
