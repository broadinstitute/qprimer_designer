"""Tests for ML model architectures and inference utilities."""

import numpy as np
import pytest
import torch
import torch.nn as nn

from qprimer_designer.models.architectures import (
    PGC,
    DropoutNd,
    S4DKernel,
    S4D,
    Janus,
    MLP,
    CombinedModel,
    CombinedModelClassifier,
    PcrDataset,
)
from qprimer_designer.models.inference import (
    get_model_path,
    FEATURE_COLUMNS,
    SEQUENCE_COLUMNS,
    _register_model_classes_for_pickle,
)


# --- Architecture tests ---


class TestPGC:
    """Tests for Position-wise Gated Convolution block."""

    def test_output_shape(self):
        model = PGC(model_dim=32)
        x = torch.randn(2, 10, 32)
        out = model(x)
        assert out.shape == (2, 10, 32)

    def test_expansion_factor(self):
        model = PGC(model_dim=32, expansion_factor=2.0)
        x = torch.randn(2, 10, 32)
        out = model(x)
        assert out.shape == (2, 10, 32)

    def test_with_dropout(self):
        model = PGC(model_dim=16, dropout=0.1)
        x = torch.randn(1, 5, 16)
        out = model(x)
        assert out.shape == (1, 5, 16)


class TestDropoutNd:
    """Tests for N-dimensional dropout."""

    def test_valid_probability(self):
        DropoutNd(p=0.0)
        DropoutNd(p=0.5)

    def test_invalid_probability(self):
        with pytest.raises(ValueError, match="dropout probability"):
            DropoutNd(p=1.0)
        with pytest.raises(ValueError, match="dropout probability"):
            DropoutNd(p=-0.1)

    def test_eval_mode_passthrough(self):
        """In eval mode, dropout should not modify input."""
        model = DropoutNd(p=0.5)
        model.eval()
        x = torch.randn(2, 4, 10)
        out = model(x)
        torch.testing.assert_close(out, x)

    def test_train_mode_changes_output(self):
        """In train mode with p>0, output should differ from input (probabilistically)."""
        model = DropoutNd(p=0.5)
        model.train()
        torch.manual_seed(42)
        x = torch.ones(2, 32, 100)
        out = model(x)
        # With p=0.5 on a large tensor, some values should be zero
        assert (out == 0).any()

    def test_shape_preserved(self):
        model = DropoutNd(p=0.3)
        x = torch.randn(2, 4, 10)
        out = model(x)
        assert out.shape == x.shape


class TestS4DKernel:
    """Tests for S4D kernel."""

    def test_output_shape(self):
        kernel = S4DKernel(model_dim=16, N=32)
        K = kernel(L=20)
        assert K.shape == (16, 20)

    def test_register_buffer(self):
        kernel = S4DKernel(model_dim=8, N=16, lr=0.0)
        # With lr=0.0, tensors are registered as buffers
        assert hasattr(kernel, 'log_dt')

    def test_register_parameter(self):
        kernel = S4DKernel(model_dim=8, N=16)
        # Without lr=0.0, tensors are parameters
        assert isinstance(kernel.log_dt, nn.Parameter)


class TestS4D:
    """Tests for S4D layer."""

    def test_output_shape_transposed(self):
        model = S4D(model_dim=16, state_dim=32, transposed=True)
        x = torch.randn(2, 16, 10)  # (batch, dim, length)
        out = model(x)
        assert out.shape == (2, 16, 10)

    def test_output_shape_not_transposed(self):
        model = S4D(model_dim=16, state_dim=32, transposed=False)
        x = torch.randn(2, 10, 16)  # (batch, length, dim)
        out = model(x)
        assert out.shape == (2, 10, 16)

    def test_no_dropout(self):
        model = S4D(model_dim=16, dropout=0.0)
        x = torch.randn(1, 16, 5)
        out = model(x)
        assert out.shape == (1, 16, 5)


class TestJanus:
    """Tests for Janus architecture."""

    def test_output_shape(self):
        model = Janus(input_dim=10, output_dim=8, model_dim=32)
        x = torch.randn(2, 56, 10)
        out = model(x)
        assert out.shape == (2, 8)

    def test_single_sample(self):
        model = Janus(input_dim=10, output_dim=4, model_dim=16)
        x = torch.randn(1, 20, 10)
        out = model(x)
        assert out.shape == (1, 4)


class TestMLP:
    """Tests for MLP."""

    def test_output_shape(self):
        model = MLP(input_dim=10, output_dim=5, hidden_dim=20)
        x = torch.randn(3, 10)
        out = model(x)
        assert out.shape == (3, 5)

    def test_single_sample(self):
        model = MLP(input_dim=8, output_dim=1, hidden_dim=16)
        x = torch.randn(1, 8)
        out = model(x)
        assert out.shape == (1, 1)

    def test_forward_deterministic(self):
        """MLP without dropout should be deterministic."""
        model = MLP(input_dim=4, output_dim=2, hidden_dim=8)
        model.eval()
        x = torch.randn(1, 4)
        out1 = model(x)
        out2 = model(x)
        torch.testing.assert_close(out1, out2)


class TestCombinedModel:
    """Tests for CombinedModel (regression)."""

    def test_output_shape(self):
        model = CombinedModel(
            DL=Janus,
            mlp_dims=(12, 8, 16),
            dl_dims=(10, 8, 32),
            combined_hidden=16,
            final_output=1,
        )
        mlp_input = torch.randn(2, 12)
        dl_input = torch.randn(2, 56, 10)
        out = model(mlp_input, dl_input)
        assert out.shape == (2, 1)


class TestCombinedModelClassifier:
    """Tests for CombinedModelClassifier (classification)."""

    def test_output_shape(self):
        model = CombinedModelClassifier(
            mlp_dims=(12, 8, 16),
            ssm_dims=(10, 8, 32),
            combined_hidden=16,
            final_output=1,
        )
        mlp_input = torch.randn(2, 12)
        ssm_input = torch.randn(2, 56, 10)
        out = model(mlp_input, ssm_input)
        assert out.shape == (2, 1)

    def test_sigmoid_output_range(self):
        """Classifier output should be in [0, 1] due to sigmoid."""
        model = CombinedModelClassifier(
            mlp_dims=(12, 8, 16),
            ssm_dims=(10, 8, 32),
            combined_hidden=16,
            final_output=1,
        )
        model.eval()
        mlp_input = torch.randn(5, 12)
        ssm_input = torch.randn(5, 56, 10)
        with torch.no_grad():
            out = model(mlp_input, ssm_input)
        assert (out >= 0).all()
        assert (out <= 1).all()


class TestPcrDataset:
    """Tests for PcrDataset."""

    def test_len(self):
        ds = PcrDataset(
            encoded_input=torch.randn(10, 56, 10),
            custom_features=np.zeros((10, 12)),
            scores=np.zeros(10),
        )
        assert len(ds) == 10

    def test_getitem(self):
        encoded = torch.randn(5, 56, 10)
        features = np.ones((5, 12))
        scores = np.arange(5, dtype=float)
        ds = PcrDataset(encoded, features, scores)

        enc, feat, score = ds[2]
        torch.testing.assert_close(enc, encoded[2])
        assert score == 2.0

    def test_empty(self):
        ds = PcrDataset(
            encoded_input=torch.randn(0, 56, 10),
            custom_features=np.zeros((0, 12)),
            scores=np.zeros(0),
        )
        assert len(ds) == 0


# --- Inference utility tests ---


class TestGetModelPath:
    """Tests for get_model_path."""

    def test_returns_path(self):
        path = get_model_path("combined_classifier.pth")
        assert "combined_classifier.pth" in str(path)

    def test_returns_path_for_scaler(self):
        path = get_model_path("standard_scaler.joblib")
        assert "standard_scaler.joblib" in str(path)


class TestConstants:
    """Tests for inference module constants."""

    def test_feature_columns_count(self):
        assert len(FEATURE_COLUMNS) == 12

    def test_feature_columns_contents(self):
        assert "len_f" in FEATURE_COLUMNS
        assert "Tm_f" in FEATURE_COLUMNS
        assert "GC_f" in FEATURE_COLUMNS
        assert "prod_len" in FEATURE_COLUMNS
        assert "prod_Tm" in FEATURE_COLUMNS

    def test_sequence_columns(self):
        assert SEQUENCE_COLUMNS == ["pseq_f", "tseq_f", "pseq_r", "tseq_r"]


class TestRegisterModelClasses:
    """Tests for _register_model_classes_for_pickle."""

    def test_registers_classes(self):
        import __main__
        _register_model_classes_for_pickle()
        assert hasattr(__main__, 'CombinedModelClassifier')
        assert hasattr(__main__, 'MLP')
        assert hasattr(__main__, 'Janus')
        assert hasattr(__main__, 'PGC')
        assert hasattr(__main__, 'S4D')
