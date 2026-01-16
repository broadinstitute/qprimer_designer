"""Pytest configuration and fixtures."""

import pytest


@pytest.fixture
def sample_dna_sequence():
    """Return a sample DNA sequence for testing."""
    return "ATCGATCGATCGATCGATCG"


@pytest.fixture
def sample_primer_pair():
    """Return a sample primer pair for testing."""
    return {
        "forward": "ATCGATCGATCGATCGATCG",
        "reverse": "CGATCGATCGATCGATCGAT",
    }
