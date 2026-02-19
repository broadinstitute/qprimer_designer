"""Tests for model loading and inference."""

from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest
import torch

from qprimer_designer.models.inference import (
    get_model_path,
    load_scaler,
    load_classifier,
    load_regressor,
    load_models,
)


class TestLoadScaler:
    """Tests for load_scaler."""

    def test_loads_scaler(self):
        """Scaler should load successfully and have a transform method."""
        scaler = load_scaler()
        assert hasattr(scaler, "transform")
        assert hasattr(scaler, "mean_")


class TestLoadClassifier:
    """Tests for load_classifier."""

    def test_loads_on_cpu(self):
        """Classifier should load on CPU and be in eval mode."""
        model = load_classifier(device="cpu")
        assert not model.training  # eval mode
        # Should be a CombinedModelClassifier
        assert hasattr(model, "mlp")
        assert hasattr(model, "ssm")
        assert hasattr(model, "combiner")


class TestLoadRegressor:
    """Tests for load_regressor."""

    def test_loads_on_cpu(self):
        """Regressor should load on CPU and be in eval mode."""
        model = load_regressor(device="cpu")
        assert not model.training
        assert hasattr(model, "mlp")
        assert hasattr(model, "combiner")


class TestLoadModels:
    """Tests for load_models."""

    def test_loads_all(self):
        """Should load scaler, classifier, regressor, and report device."""
        scaler, classifier, regressor, device = load_models(device="cpu")
        assert hasattr(scaler, "transform")
        assert not classifier.training
        assert not regressor.training
        assert device == "cpu"

    def test_auto_device(self):
        """With device=None, should auto-detect (CPU in test env)."""
        scaler, classifier, regressor, device = load_models(device=None)
        assert device in ("cpu", "cuda")
