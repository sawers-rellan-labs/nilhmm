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


def test_introgression_hmm_counts_detects_donor_block():
    """Count-based HMM should recover a donor-homozygous block from sparse AD."""
    from nilhmm.core import introgression_hmm_counts

    rng = np.random.default_rng(0)
    n_samples, n_markers = 4, 600
    ref = np.zeros((n_samples, n_markers), dtype=int)
    alt = np.zeros((n_samples, n_markers), dtype=int)

    # ~30% of sites covered with a single read (depth-1 skim regime)
    covered = rng.random((n_samples, n_markers)) < 0.3
    block = slice(200, 400)  # true donor-homozygous segment

    for i in range(n_samples):
        for j in range(n_markers):
            if not covered[i, j]:
                continue
            # donor-hom block emits ALT reads; elsewhere B73 emits REF reads
            if block.start <= j < block.stop:
                alt[i, j] = 1
            else:
                ref[i, j] = 1

    marker_dict = {1: list(range(n_markers))}
    calls = introgression_hmm_counts(ref, alt, marker_dict, r=0.01,
                                     f_1=0.0625, f_2=0.0938, return_calls=True)
    assert calls.shape == (n_samples, n_markers)
    assert np.all(np.isin(calls, [0, 1, 2]))
    # the block should be called donor (state 2) the large majority of the time
    block_calls = calls[:, block]
    assert (block_calls == 2).mean() > 0.8
    # flanks should be mostly B73 (state 0)
    flank = np.c_[calls[:, :200], calls[:, 400:]]
    assert (flank == 0).mean() > 0.8


def test_read_vcf_counts(tmp_path):
    """read_vcf_counts parses FORMAT/AD into ref/alt matrices."""
    from nilhmm.io import read_vcf_counts

    vcf = tmp_path / "t.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
        "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT:AD\t0/0:5,0\t./.:.",
        "chr1\t200\trs2\tA\tT\t.\tPASS\t.\tGT:AD\t1/1:0,3\t0/1:2,2",
    ]
    vcf.write_text("\n".join(lines) + "\n")
    ref, alt, mdict, samples, minfo = read_vcf_counts(str(vcf))
    assert samples == ["S1", "S2"]
    assert ref.shape == (2, 2) and alt.shape == (2, 2)
    assert list(ref[0]) == [5, 0] and list(alt[0]) == [0, 3]   # S1
    assert list(ref[1]) == [0, 2] and list(alt[1]) == [0, 2]   # S2 (missing AD -> 0,0)
    assert mdict[1] == [0, 1]
