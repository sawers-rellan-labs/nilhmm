"""Tests for nilhmm core functions."""

import pytest
import numpy as np
import pandas as pd
from nilhmm.core import introgression_hmm, call_introgressions
from nilhmm.io import read_vcf


def test_introgression_hmm_basic():
    """Test basic HMM functionality with minimal data."""
    # Create minimal test data
    n_samples, n_markers = 10, 100
    geno_matrix = np.random.randint(0, 4, size=(n_samples, n_markers))

    # Create marker dictionary
    marker_dict = {1: list(range(n_markers))}

    # Run HMM
    calls = introgression_hmm(
        geno=geno_matrix,
        marker_dict=marker_dict,
        return_calls=True
    )

    # Check output shape
    assert calls.shape == (n_samples, n_markers)

    # Check output values are valid
    assert np.all(np.isin(calls, [0, 1, 2]))


def test_hmm_parameters():
    """Test HMM with different parameter values."""
    n_samples, n_markers = 5, 50
    geno_matrix = np.random.randint(0, 4, size=(n_samples, n_markers))
    marker_dict = {1: list(range(n_markers))}

    # Test with different parameters
    calls1 = introgression_hmm(
        geno=geno_matrix,
        marker_dict=marker_dict,
        nir=0.01,
        germ=0.05,
        return_calls=True
    )

    calls2 = introgression_hmm(
        geno=geno_matrix,
        marker_dict=marker_dict,
        nir=0.1,
        germ=0.1,
        return_calls=True
    )

    # Results should be different with different parameters
    assert not np.array_equal(calls1, calls2)


def test_multiple_chromosomes():
    """Test HMM with multiple chromosomes."""
    n_samples = 5
    n_markers_chr1, n_markers_chr2 = 30, 20
    total_markers = n_markers_chr1 + n_markers_chr2

    geno_matrix = np.random.randint(0, 4, size=(n_samples, total_markers))

    marker_dict = {
        1: list(range(n_markers_chr1)),
        2: list(range(n_markers_chr1, total_markers))
    }

    calls = introgression_hmm(
        geno=geno_matrix,
        marker_dict=marker_dict,
        return_calls=True
    )

    assert calls.shape == (n_samples, total_markers)
