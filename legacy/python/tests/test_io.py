"""Tests for nilhmm I/O functions."""

import pytest
import numpy as np
import pandas as pd
import tempfile
import os
from nilhmm.io import parse_vcf_header, write_results


def test_write_results():
    """Test result writing functionality."""
    # Create test results
    n_samples, n_markers = 5, 20
    calls = np.random.randint(0, 3, size=(n_samples, n_markers))

    sample_names = [f"sample_{i}" for i in range(n_samples)]

    marker_info = pd.DataFrame({
        'CHROM': [1] * n_markers,
        'POS': range(1000, 1000 + n_markers),
        'ID': [f"marker_{i}" for i in range(n_markers)],
        'REF': ['A'] * n_markers,
        'ALT': ['T'] * n_markers
    })

    parameters = {"nir": 0.01, "germ": 0.05}

    results = {
        "calls": calls,
        "sample_names": sample_names,
        "marker_info": marker_info,
        "parameters": parameters
    }

    # Test writing to temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        output_prefix = os.path.join(tmpdir, "test_results")

        write_results(results, output_prefix)

        # Check that files were created
        expected_files = [
            "test_results_introgression_calls.txt",
            "test_results_introgression_calls.csv",
            "test_results_introgression_summary.csv",
            "test_results_marker_info.csv",
            "test_results_parameters.csv"
        ]

        for filename in expected_files:
            filepath = os.path.join(tmpdir, filename)
            assert os.path.exists(filepath), f"Missing file: {filename}"

            # Check file is not empty
            assert os.path.getsize(filepath) > 0, f"Empty file: {filename}"


def test_parse_vcf_header():
    """Test VCF header parsing with mock data."""
    # Create a mock VCF file
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=1,length=307041717>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
1	1000	.	A	T	60	PASS	.	GT	0/0	0/1	1/1
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_content)
        f.flush()

        try:
            samples = parse_vcf_header(f.name)
            expected_samples = ["sample1", "sample2", "sample3"]
            assert samples == expected_samples
        finally:
            os.unlink(f.name)
