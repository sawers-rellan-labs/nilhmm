"""Core HMM functions for introgression calling in NIL populations."""

import numpy as np
from hmmlearn import hmm
from typing import Dict, List, Tuple, Optional
import logging


def introgression_hmm(
    geno: np.ndarray,
    marker_dict: Dict[int, List[int]],
    nir: float = 0.01,
    germ: float = 0.05,
    gert: float = 0.10,
    p: float = 0.5,
    mr: float = 0.15,
    r: float = 0.01,
    f_1: float = 0.25,
    f_2: float = 0.05,
    return_calls: bool = True
) -> Optional[np.ndarray]:
    """
    Call introgressions using HMM approach.

    Refactored from File_S11_callIntrogressions.py to be more modular.

    Parameters:
    -----------
    geno : np.ndarray
        Genotype matrix (individuals x markers), 0/1/2/3 encoding
    marker_dict : Dict[int, List[int]]
        Dictionary mapping chromosomes to marker indices
    nir : float
        Non-informative rate
    germ : float
        SNP calling error rate on true introgression homozygotes
    gert : float
        SNP calling error rate on true introgression heterozygotes
    p : float
        Proportion of homozygous SNP call errors resulting in het call
    mr : float
        Missing call rate
    r : float
        Expected recombination rate between adjacent markers
    f_1 : float
        Expected frequency of heterozygotes
    f_2 : float
        Expected frequency of homozygous introgressions
    return_calls : bool
        Whether to return Viterbi best calls

    Returns:
    --------
    np.ndarray or None
        Introgression calls matrix if return_calls=True
    """
    # Starting probabilities
    f_0 = 1 - f_1 - f_2  # B73 homozygous frequency

    states = [0, 1, 2]  # "B73", "het", "donor"

    # Transition probabilities
    p00 = 1 - r
    p01 = r * (f_1 / (f_1 + f_2))
    p02 = r * (f_2 / (f_1 + f_2))
    p10 = r * f_0 / (f_0 + f_2)
    p11 = 1 - r
    p12 = r * f_2 / (f_0 + f_2)
    p20 = r * f_0 / (f_0 + f_1)
    p21 = r * f_1 / (f_0 + f_1)
    p22 = 1 - r

    tmat = np.array([
        [p00, p01, p02],
        [p10, p11, p12],
        [p20, p21, p22]
    ])

    # Emission probabilities
    emimat = np.array([
        [(1-germ)*(1-mr), p*germ*(1-mr), (1-p)*germ*(1-mr), mr],
        [(((1-nir)*0.5*gert) + nir*(1-germ))*(1-mr),
         (((1-nir)*(1-gert)) + (nir*germ*p))*(1-mr),
         (((1-nir)*0.5*gert) + nir*germ*(1-p))*(1-mr), mr],
        [((1-nir)*germ*(1-p) + (nir*(1-germ)))*(1-mr),
         germ*p*(1-mr),
         ((1-nir)*(1-germ) + (nir*germ*(1-p)))*(1-mr), mr]
    ])

    # Set up the HMM
    model = hmm.MultinomialHMM(
        n_components=len(states),
        n_trials=1,
        init_params=''
    )

    model.n_features = len(states) + 1  # Additional NA 'feature'
    model.startprob_ = np.array([f_0, f_1, f_2])
    model.transmat_ = tmat
    model.emissionprob_ = emimat

    # Process each chromosome
    results = {}
    for chrom in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
        results[chrom] = []
        geno_current_chr = geno[:, marker_dict[chrom]]

        logging.info(f"Processing chromosome {chrom}")

        for i in range(geno.shape[0]):
            nil_i = np.identity(4)[np.nan_to_num(geno_current_chr[i, :]).astype(int)]
            preds = model.predict(nil_i).reshape(1, -1)
            results[chrom].append(preds)

    # Pack results back into single array
    results_by_chrom = {}
    for chrom in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
        results_by_chrom[chrom] = np.concatenate(results[chrom], axis=0)

    # Merge chromosomes into one array
    nil_calls = np.concatenate(list(results_by_chrom.values()), axis=1)

    if return_calls:
        return nil_calls

    return None


def call_introgressions(
    vcf_file: str,
    output_prefix: str = "introgressions",
    coverage_level: str = "low",
    **hmm_params
) -> Dict:
    """
    Main function to call introgressions from VCF file.

    Parameters:
    -----------
    vcf_file : str
        Path to input VCF file
    output_prefix : str
        Prefix for output files
    coverage_level : str
        Coverage level ("low", "medium", "high") for parameter defaults
    **hmm_params
        Additional HMM parameters to override defaults

    Returns:
    --------
    Dict
        Results dictionary with calls and metadata
    """
    from .io import read_vcf, write_results

    # Load data
    logging.info(f"Reading VCF file: {vcf_file}")
    geno_matrix, marker_dict, sample_names, marker_info = read_vcf(vcf_file)

    # Set coverage-specific defaults
    coverage_defaults = {
        "low": {"nir": 0.02, "germ": 0.08, "gert": 0.15, "mr": 0.20},
        "medium": {"nir": 0.01, "germ": 0.05, "gert": 0.10, "mr": 0.10},
        "high": {"nir": 0.005, "germ": 0.02, "gert": 0.05, "mr": 0.05}
    }

    # Merge defaults with user parameters
    params = coverage_defaults.get(coverage_level, coverage_defaults["low"])
    params.update(hmm_params)

    logging.info(f"Using parameters: {params}")

    # Call introgressions
    calls = introgression_hmm(
        geno=geno_matrix,
        marker_dict=marker_dict,
        **params
    )

    # Prepare results
    results = {
        "calls": calls,
        "sample_names": sample_names,
        "marker_info": marker_info,
        "parameters": params
    }

    # Write results
    write_results(results, output_prefix)

    logging.info(f"Analysis complete. Results saved with prefix: {output_prefix}")

    return results
