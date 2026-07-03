"""HMM parameter optimization using grid search."""

import numpy as np
import pandas as pd
import itertools
from typing import Dict, List, Optional
import logging
from .core import introgression_hmm


def optimize_parameters(
    geno_matrix: np.ndarray,
    marker_dict: Dict[int, List[int]],
    nir_values: List[float] = [0.001, 0.01, 0.1, 0.3, 0.5],
    germ_values: List[float] = [0.001, 0.01, 0.05],
    gert_values: List[float] = [0.001, 0.01, 0.05],
    p_values: List[float] = [0.25, 0.5, 0.75],
    r_multipliers: List[float] = [0.5, 1.0, 2.0],
    base_r: float = 0.01,
    output_file: Optional[str] = None
) -> pd.DataFrame:
    """
    Perform grid search over HMM parameters.

    Parameters:
    -----------
    geno_matrix : np.ndarray
        Genotype matrix (samples x markers)
    marker_dict : Dict[int, List[int]]
        Dictionary mapping chromosomes to marker indices
    nir_values : List[float]
        Non-informative rate values to test
    germ_values : List[float]
        Genotyping error rate values (homozygotes) to test
    gert_values : List[float]
        Genotyping error rate values (heterozygotes) to test
    p_values : List[float]
        Proportion of homozygous errors as heterozygous calls to test
    r_multipliers : List[float]
        Multipliers for base recombination rate
    base_r : float
        Base recombination rate between markers
    output_file : Optional[str]
        Path to save results CSV file

    Returns:
    --------
    pd.DataFrame
        Results for all parameter combinations
    """
    logging.info("Starting HMM parameter grid search")

    results = []
    total_combinations = (len(nir_values) * len(germ_values) *
                         len(gert_values) * len(p_values) * len(r_multipliers))

    logging.info(f"Testing {total_combinations} parameter combinations")

    current = 0
    for nir, germ, gert, p, r_mult in itertools.product(
        nir_values, germ_values, gert_values, p_values, r_multipliers
    ):
        current += 1
        r = base_r * r_mult

        if current % 10 == 0:
            logging.info(f"Progress: {current}/{total_combinations}")

        try:
            # Run HMM with current parameters
            calls = introgression_hmm(
                geno=geno_matrix,
                marker_dict=marker_dict,
                nir=nir,
                germ=germ,
                gert=gert,
                p=p,
                r=r,
                return_calls=True
            )

            # Calculate quality metrics
            metrics = calculate_quality_metrics(calls)

            # Store results
            result = {
                "nir": nir,
                "germ": germ,
                "gert": gert,
                "p": p,
                "r": r,
                **metrics
            }
            results.append(result)

        except Exception as e:
            logging.warning(f"Failed for nir={nir}, germ={germ}, gert={gert}, p={p}, r={r}: {e}")
            continue

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    if output_file:
        results_df.to_csv(output_file, index=False)
        logging.info(f"Grid search results saved to {output_file}")

    logging.info("Grid search completed")
    return results_df


def calculate_quality_metrics(calls: np.ndarray) -> Dict[str, float]:
    """
    Calculate quality metrics for introgression calls.

    Parameters:
    -----------
    calls : np.ndarray
        Introgression calls matrix (samples x markers)

    Returns:
    --------
    Dict[str, float]
        Quality metrics
    """
    n_samples, n_markers = calls.shape

    # Calculate percentages of each call type
    pct_b73 = np.mean(calls == 0) * 100
    pct_het = np.mean(calls == 1) * 100
    pct_donor = np.mean(calls == 2) * 100

    # Calculate per-sample statistics
    sample_het_rates = np.mean(calls == 1, axis=1)
    sample_donor_rates = np.mean(calls == 2, axis=1)

    # Calculate number of samples with no introgressions
    samples_no_introg = np.sum(np.all(calls == 0, axis=1))
    pct_samples_no_introg = (samples_no_introg / n_samples) * 100

    return {
        "pct_b73_calls": pct_b73,
        "pct_het_calls": pct_het,
        "pct_donor_calls": pct_donor,
        "mean_het_rate": np.mean(sample_het_rates),
        "mean_donor_rate": np.mean(sample_donor_rates),
        "std_het_rate": np.std(sample_het_rates),
        "std_donor_rate": np.std(sample_donor_rates),
        "pct_samples_no_introg": pct_samples_no_introg
    }


def select_best_parameters(
    results_df: pd.DataFrame,
    criteria: str = "donor_rate"
) -> Dict[str, float]:
    """
    Select best parameters from grid search results.

    Parameters:
    -----------
    results_df : pd.DataFrame
        Grid search results
    criteria : str
        Selection criteria ("donor_rate", "het_rate", or "balanced")

    Returns:
    --------
    Dict[str, float]
        Best parameter combination
    """
    if criteria == "donor_rate":
        # Maximize donor introgression rate
        best_idx = results_df["mean_donor_rate"].idxmax()
    elif criteria == "het_rate":
        # Maximize heterozygous rate (for F2 populations)
        best_idx = results_df["mean_het_rate"].idxmax()
    elif criteria == "balanced":
        # Balance donor rate and minimize samples with no introgressions
        score = (results_df["mean_donor_rate"] * 0.7 -
                results_df["pct_samples_no_introg"] * 0.3)
        best_idx = score.idxmax()
    else:
        raise ValueError(f"Unknown criteria: {criteria}")

    best_row = results_df.loc[best_idx]

    return {
        "nir": best_row["nir"],
        "germ": best_row["germ"],
        "gert": best_row["gert"],
        "p": best_row["p"],
        "r": best_row["r"]
    }
