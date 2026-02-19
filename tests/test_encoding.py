"""Tests for sequence encoding utilities."""

import numpy as np
import pandas as pd
import pytest

from qprimer_designer.utils.encoding import one_hot_encode, encode_primer_pair, encode_batch_parallel


class TestOneHotEncode:
    """Tests for one_hot_encode."""

    def test_basic_encoding(self):
        """Test encoding of standard nucleotides."""
        result = one_hot_encode("ATCG", length=4)
        expected = np.array([
            [1, 0, 0, 0, 0],  # A
            [0, 1, 0, 0, 0],  # T
            [0, 0, 1, 0, 0],  # C
            [0, 0, 0, 1, 0],  # G
        ])
        np.testing.assert_array_equal(result, expected)

    def test_shape(self):
        """Test output shape with default length."""
        result = one_hot_encode("ATCG")
        assert result.shape == (28, 5)

    def test_custom_length(self):
        """Test output shape with custom length."""
        result = one_hot_encode("ATCG", length=10)
        assert result.shape == (10, 5)

    def test_gap_encoding(self):
        """Test encoding of gap character."""
        result = one_hot_encode("-", length=1)
        expected = np.array([[0, 0, 0, 0, 1]])
        np.testing.assert_array_equal(result, expected)

    def test_n_encoding(self):
        """N should encode as all zeros."""
        result = one_hot_encode("N", length=1)
        expected = np.array([[0, 0, 0, 0, 0]])
        np.testing.assert_array_equal(result, expected)

    def test_right_padding_short_sequence(self):
        """Short sequences are right-padded with N (all zeros)."""
        result = one_hot_encode("AT", length=4)
        # AT followed by two N positions (all zeros)
        assert result.shape == (4, 5)
        np.testing.assert_array_equal(result[0], [1, 0, 0, 0, 0])  # A
        np.testing.assert_array_equal(result[1], [0, 1, 0, 0, 0])  # T
        np.testing.assert_array_equal(result[2], [0, 0, 0, 0, 0])  # N pad
        np.testing.assert_array_equal(result[3], [0, 0, 0, 0, 0])  # N pad

    def test_left_truncation_long_sequence(self):
        """Long sequences keep the last N characters."""
        result = one_hot_encode("ATCGATCG", length=4)
        # Last 4 characters: ATCG
        assert result.shape == (4, 5)
        np.testing.assert_array_equal(result[0], [1, 0, 0, 0, 0])  # A
        np.testing.assert_array_equal(result[1], [0, 1, 0, 0, 0])  # T
        np.testing.assert_array_equal(result[2], [0, 0, 1, 0, 0])  # C
        np.testing.assert_array_equal(result[3], [0, 0, 0, 1, 0])  # G

    def test_exact_length(self):
        """Exact-length sequence should be unchanged."""
        result = one_hot_encode("ATCG", length=4)
        assert result.shape == (4, 5)
        np.testing.assert_array_equal(result[0], [1, 0, 0, 0, 0])

    def test_lowercase_input(self):
        """Lowercase input should be handled (uppercased)."""
        result = one_hot_encode("atcg", length=4)
        expected = one_hot_encode("ATCG", length=4)
        np.testing.assert_array_equal(result, expected)


class TestEncodePrimerPair:
    """Tests for encode_primer_pair."""

    def test_output_shape(self):
        """Output should be (56, 10)."""
        row = {
            "pseq_f": "ATCGATCGATCGATCGATCG",
            "tseq_f": "ATCGATCGATCGATCGATCG",
            "pseq_r": "CGATCGATCGATCGATCGAT",
            "tseq_r": "CGATCGATCGATCGATCGAT",
        }
        result = encode_primer_pair(row)
        assert result.shape == (56, 10)

    def test_dtype(self):
        """Output should be numeric numpy array."""
        row = {
            "pseq_f": "ATCG",
            "tseq_f": "ATCG",
            "pseq_r": "GCTA",
            "tseq_r": "GCTA",
        }
        result = encode_primer_pair(row)
        assert isinstance(result, np.ndarray)
        assert result.dtype in (np.int64, np.float64, np.int32, np.float32)

    def test_structure(self):
        """First 28 rows are forward, last 28 are reverse; first 5 cols are target, last 5 are primer."""
        row = {
            "pseq_f": "A" * 28,
            "tseq_f": "T" * 28,
            "pseq_r": "C" * 28,
            "tseq_r": "G" * 28,
        }
        result = encode_primer_pair(row)
        # Target encoding is columns 0-4, primer encoding is columns 5-9
        # Forward target (tseq_f) is rows 0-27, reverse target (tseq_r) is rows 28-55
        assert result.shape == (56, 10)


class TestEncodeBatchParallel:
    """Tests for encode_batch_parallel."""

    def test_output_shape(self):
        """Test batch encoding produces correct tensor shape."""
        df = pd.DataFrame({
            "pseq_f": ["ATCGATCGATCGATCGATCG", "GCTAGCTAGCTAGCTAGCTA"],
            "tseq_f": ["ATCGATCGATCGATCGATCG", "GCTAGCTAGCTAGCTAGCTA"],
            "pseq_r": ["CGATCGATCGATCGATCGAT", "TAGCTAGCTAGCTAGCTAGC"],
            "tseq_r": ["CGATCGATCGATCGATCGAT", "TAGCTAGCTAGCTAGCTAGC"],
        })
        result = encode_batch_parallel(df, threads=1)
        assert result.shape == (2, 56, 10)

    def test_output_dtype(self):
        """Output should be float32 tensor."""
        import torch
        df = pd.DataFrame({
            "pseq_f": ["ATCG"],
            "tseq_f": ["ATCG"],
            "pseq_r": ["GCTA"],
            "tseq_r": ["GCTA"],
        })
        result = encode_batch_parallel(df, threads=1)
        assert isinstance(result, torch.Tensor)
        assert result.dtype == torch.float32

    def test_single_sample(self):
        """Test with a single sample."""
        df = pd.DataFrame({
            "pseq_f": ["ATCG"],
            "tseq_f": ["ATCG"],
            "pseq_r": ["GCTA"],
            "tseq_r": ["GCTA"],
        })
        result = encode_batch_parallel(df, threads=1)
        assert result.shape == (1, 56, 10)
