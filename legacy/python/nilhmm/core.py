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

    # Process each chromosome present in the marker map (sorted for stable
    # ordering). Avoids assuming all 10 maize chromosomes are present, e.g.
    # for subset/filtered VCFs or test fixtures.
    chroms = sorted(marker_dict.keys())
    results = {}
    for chrom in chroms:
        results[chrom] = []
        geno_current_chr = geno[:, marker_dict[chrom]]

        logging.info(f"Processing chromosome {chrom}")

        for i in range(geno.shape[0]):
            nil_i = np.identity(4)[np.nan_to_num(geno_current_chr[i, :]).astype(int)]
            preds = model.predict(nil_i).reshape(1, -1)
            results[chrom].append(preds)

    # Pack results back into single array
    results_by_chrom = {}
    for chrom in chroms:
        results_by_chrom[chrom] = np.concatenate(results[chrom], axis=0)

    # Merge chromosomes into one array
    nil_calls = np.concatenate(list(results_by_chrom.values()), axis=1)

    if return_calls:
        return nil_calls

    return None


def _build_transition(r: float, f_1: float, f_2: float) -> Tuple[np.ndarray, np.ndarray]:
    """Build (startprob, transmat) from recombination r and state freqs.

    Identical parameterization to ``introgression_hmm`` so the count-based and
    GT-based callers share the same transition/init structure.
    """
    f_0 = 1 - f_1 - f_2
    p01 = r * (f_1 / (f_1 + f_2)); p02 = r * (f_2 / (f_1 + f_2))
    p10 = r * f_0 / (f_0 + f_2);   p12 = r * f_2 / (f_0 + f_2)
    p20 = r * f_0 / (f_0 + f_1);   p21 = r * f_1 / (f_0 + f_1)
    tmat = np.array([
        [1 - r, p01,   p02],
        [p10,   1 - r, p12],
        [p20,   p21,   1 - r],
    ])
    startprob = np.array([f_0, f_1, f_2])
    return startprob, tmat


def _log_viterbi(log_startprob: np.ndarray, log_transmat: np.ndarray,
                 log_emission: np.ndarray) -> np.ndarray:
    """Standard log-space Viterbi. log_emission: (T, K). Returns path (T,)."""
    T, K = log_emission.shape
    delta = np.empty((T, K))
    psi = np.empty((T, K), dtype=int)
    delta[0] = log_startprob + log_emission[0]
    for t in range(1, T):
        scores = delta[t - 1][:, None] + log_transmat   # (i -> j)
        psi[t] = np.argmax(scores, axis=0)
        delta[t] = scores[psi[t], np.arange(K)] + log_emission[t]
    path = np.empty(T, dtype=int)
    path[-1] = int(np.argmax(delta[-1]))
    for t in range(T - 2, -1, -1):
        path[t] = psi[t + 1, path[t + 1]]
    return path


def introgression_hmm_counts(
    ref: np.ndarray,
    alt: np.ndarray,
    marker_dict: Dict[int, List[int]],
    err: float = 0.01,
    conc: float = 20.0,
    r: float = 0.01,
    f_1: float = 0.0625,
    f_2: float = 0.0938,
    return_calls: bool = True,
) -> Optional[np.ndarray]:
    """
    Call introgressions from allelic depth (ref/alt) counts via a beta-binomial
    emission HMM. Designed for low-coverage data where per-site GT is
    uninformative but the HMM pools single-read observations along a segment.

    States 0/1/2 = B73-hom / het / donor-hom, with expected alt fractions
    theta = [err, 0.5, 1-err]. Emission for marker i, state s is
    BetaBinomial(alt_i | n_i, theta_s*conc, (1-theta_s)*conc); markers with zero
    depth contribute a flat (uninformative) emission. ``conc`` controls
    overdispersion (conc -> inf approaches a Binomial). Transition/init reuse the
    same r / f_1 / f_2 parameterization as the GT caller.

    Implementation: emissions are evaluated once per DISTINCT (n, a) count pair
    (sparse memoization via an integer-key np.unique; cost scales with distinct
    pairs ~ coverage, not samples x markers, and is depth-robust in memory), and
    the Viterbi is BATCHED across all samples of a chromosome (vectorized over the
    sample axis). This is the same structure as RTIGER's getlogpsi memoization;
    see Implementation.md section 7. Mathematically identical to a per-sample loop:
    depth-0 maps to BetaBinomial(0|0,.) = 1 (log 0), a flat/uninformative emission.
    """
    from scipy.stats import betabinom

    startprob, tmat = _build_transition(r, f_1, f_2)
    log_start = np.log(startprob)
    log_trans = np.log(tmat)

    theta = np.array([err, 0.5, 1.0 - err])
    a_s = theta * conc
    b_s = (1.0 - theta) * conc

    n_samples, n_markers = ref.shape
    total = ref + alt

    # --- emission: one BetaBinomial eval per distinct (n, a) pair ---
    base = int(total.max()) + 1 if total.size else 1          # a <= n < base
    key = (total.astype(np.int64) * base + alt).ravel()
    uniq, inv = np.unique(key, return_inverse=True)           # distinct pairs, cell->pair map
    un = uniq // base
    ua = uniq % base
    em_uniq = np.empty((uniq.size, 3))
    for s in range(3):
        em_uniq[:, s] = betabinom.logpmf(ua, un, a_s[s], b_s[s])   # logpmf(0|0,.) = 0
    inv = inv.reshape(n_samples, n_markers)

    # --- batched Viterbi per chromosome (vectorized over samples) ---
    calls = np.zeros((n_samples, n_markers), dtype=np.int8)
    for chrom in sorted(marker_dict.keys()):
        idx = np.asarray(marker_dict[chrom])
        if idx.size == 0:
            continue
        logging.info(f"Processing chromosome {chrom} ({idx.size} markers)")
        E = em_uniq[inv[:, idx]]                              # (S, T, 3)
        S, T = E.shape[0], E.shape[1]
        delta = log_start[None, :] + E[:, 0, :]              # (S, 3)
        psi = np.empty((T, S, 3), dtype=np.int8)
        for t in range(1, T):
            scores = delta[:, :, None] + log_trans[None, :, :]        # (S, prev, cur)
            best = np.argmax(scores, axis=1)                          # (S, cur)
            psi[t] = best
            delta = np.take_along_axis(scores, best[:, None, :], axis=1)[:, 0, :] + E[:, t, :]
        path = np.empty((T, S), dtype=np.int8)
        path[T - 1] = np.argmax(delta, axis=1)
        for t in range(T - 2, -1, -1):
            path[t] = np.take_along_axis(psi[t + 1], path[t + 1][:, None], axis=1)[:, 0]
        calls[:, idx] = path.T

    if return_calls:
        return calls.astype(int)
    return None


def call_introgressions_counts(
    vcf_file: str,
    output_prefix: str = "introgressions",
    **hmm_params,
) -> Dict:
    """Call introgressions from VCF FORMAT/AD counts and write standard outputs."""
    from .io import read_vcf_counts, write_results

    logging.info(f"Reading VCF (AD counts): {vcf_file}")
    ref, alt, marker_dict, sample_names, marker_info = read_vcf_counts(vcf_file)

    defaults = {"err": 0.01, "conc": 20.0, "r": 0.01, "f_1": 0.0625, "f_2": 0.0938}
    params = {**defaults, **hmm_params}
    logging.info(f"Using count-HMM parameters: {params}")

    calls = introgression_hmm_counts(ref=ref, alt=alt, marker_dict=marker_dict, **params)

    results = {"calls": calls, "sample_names": sample_names,
               "marker_info": marker_info, "parameters": params}
    write_results(results, output_prefix)
    logging.info(f"Analysis complete. Results saved with prefix: {output_prefix}")
    return results


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
