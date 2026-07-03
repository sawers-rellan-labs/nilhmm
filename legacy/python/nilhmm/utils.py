"""Utility functions for nilhmm package."""

import numpy as np
import pandas as pd
from typing import Dict, List
import logging


def estimate_data_parameters(
    geno_matrix: np.ndarray,
    expected_maf: float = 0.0156
) -> Dict[str, float]:
    """
    Estimate data-specific parameters from genotype matrix.

    Parameters:
    -----------
    geno_matrix : np.ndarray
        Genotype matrix (samples x markers)
    expected_maf : float
        Expected minor allele frequency for population

    Returns:
    --------
    Dict[str, float]
        Estimated parameters
    """
    # Convert missing values to NaN for calculation
    geno_clean = geno_matrix.astype(float)
    geno_clean[geno_clean == 3] = np.nan

    # Estimate missing data rate
    missing_rate = np.mean(geno_matrix == 3)

    # Estimate minor allele frequency
    observed_maf = np.nanmean(geno_clean, axis=0) * 0.5
    mean_maf = np.nanmean(observed_maf)

    # Estimate non-informative rate
    nir = max((expected_maf - mean_maf) / expected_maf, 0.001)

    return {
        "missing_rate": missing_rate,
        "observed_maf": mean_maf,
        "estimated_nir": nir,
        "expected_maf": expected_maf
    }


def calculate_recombination_rate(
    marker_positions: pd.DataFrame,
    total_map_length: float = 1500.0,
    generations: int = 2
) -> float:
    """
    Calculate average recombination rate between markers.

    Parameters:
    -----------
    marker_positions : pd.DataFrame
        DataFrame with marker positions
    total_map_length : float
        Total genetic map length in cM (default: 1500 for maize)
    generations : int
        Effective number of meioses (default: 2 for BC2S3)

    Returns:
    --------
    float
        Average recombination rate between adjacent markers
    """
    n_markers = len(marker_positions)
    avg_r = generations * total_map_length / (100 * n_markers)
    return avg_r


def validate_genotype_matrix(geno_matrix: np.ndarray) -> bool:
    """
    Validate genotype matrix format.

    Parameters:
    -----------
    geno_matrix : np.ndarray
        Genotype matrix to validate

    Returns:
    --------
    bool
        True if valid, raises ValueError if invalid
    """
    if not isinstance(geno_matrix, np.ndarray):
        raise ValueError("Genotype matrix must be numpy array")

    if len(geno_matrix.shape) != 2:
        raise ValueError("Genotype matrix must be 2D")

    valid_values = {0, 1, 2, 3}
    unique_values = set(np.unique(geno_matrix))

    if not unique_values.issubset(valid_values):
        raise ValueError(f"Invalid genotype values: {unique_values - valid_values}")

    logging.info(f"Genotype matrix validation passed: {geno_matrix.shape}")
    return True


def setup_logging(level: str = "INFO") -> None:
    """
    Set up logging for nilhmm package.

    Parameters:
    -----------
    level : str
        Logging level ("DEBUG", "INFO", "WARNING", "ERROR")
    """
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {level}')

    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    logging.info(f"nilhmm logging initialized at {level} level")
