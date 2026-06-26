# NILHMM: Hidden Markov Model for NIL Introgression Analysis

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

NILHMM implements Hidden Markov Model approaches for calling introgressions in Near Isogenic Lines (NILs) from sequencing data. The package is optimized for maize genetics applications, particularly low-coverage sequencing of populations like BZea (B73 × Teosinte).

Based upon Jim Holland's methodology:
- Repository: https://github.com/ncsumaize/nNIL
- Preprint: Zhong, T. et al. (2025). A maize near-isogenic line population designed for gene discovery and characterization of allelic effects. bioRxiv. https://doi.org/10.1101/2025.01.29.635337

## Two callers

NILHMM provides **two emission models** over the same 3-state hidden chain
(B73-homozygote / heterozygote / donor-homozygote), chosen by what the data can support:

| caller | reads | emission | use when |
|---|---|---|---|
| **GT** (original) | VCF `GT` field | 3×4 categorical with explicit genotyping-error model (`germ/gert/nir/mr`) | genotypes are called/dense (≥~0.5× or imputed) |
| **counts** (added) | VCF `FORMAT/AD` | BetaBinomial over (ref, alt) depths, per-state alt fraction θ=[`err`, 0.5, 1−`err`], concentration `conc` | ultra-low coverage where per-site GT is mostly missing / heterozygotes invisible (e.g. ~0.4× skim) |

Both share the same transition/init parameterization (recombination `r`, state
frequencies `f_1`, `f_2`) and decode per chromosome by Viterbi. The **count caller** exists
because at ~0.4× the `GT` field is mostly missing and the donor-homozygous state is nearly
invisible (a single read can't show zygosity); reading `AD` lets the HMM pool single-read
observations along a segment — the same principle RTIGER uses. See
[`Implementation.md`](Implementation.md) for the full rationale, and
[`calibration/`](calibration/) for the BZea SNP50K calibration.

## Key Features

- **Two emission models** (called GT, or allelic-depth counts) over a shared 3-state chain
- **Individual-marker resolution** with explicit genotyping-error modeling (GT caller)
- **BetaBinomial count model** for ultra-low coverage / sparse GT (count caller)
- **VCF processing** (GT and `AD`) with support for compressed files
- **Vectorized count caller** — emissions memoized over distinct (depth, alt) pairs + Viterbi
  batched across samples (cost scales with coverage, not sample×marker count)
- **Parameter optimization** / KS-against-simulation calibration tools

## Installation

### From Source
```bash
git clone https://github.com/sawers-rellan-labs/nilhmm.git
cd nilhmm
pip install -e .
```

### Dependencies
```bash
pip install numpy pandas scipy hmmlearn scikit-learn
```
(`scipy` is required by the count caller's BetaBinomial emission.)

## Quick Start

### Command Line (GT caller)
```bash
# Basic introgression calling from called genotypes
call-bzea-introgressions your_data.vcf -o results

# With custom parameters for your data
call-bzea-introgressions your_data.vcf -o results --mr 0.2 --germ 0.08
```

### Python API — GT caller
```python
import nilhmm

results = nilhmm.call_introgressions(
    vcf_file="your_data.vcf",
    output_prefix="results",
    coverage_level="low",   # low/medium/high parameter presets
)
```

### Python API — count caller (allelic depths)
```python
import nilhmm

# Reads FORMAT/AD; BetaBinomial emission tuned for ultra-low coverage
results = nilhmm.call_introgressions_counts(
    vcf_file="your_data.vcf",     # must carry FORMAT/AD
    output_prefix="results",
    r=3e-5, err=0.01, conc=20.0,  # transition + emission params (calibrate per dataset)
    f_1=0.0625, f_2=0.0938,       # state frequencies (here: BC2S2)
)

# Lower-level access
ref, alt, marker_dict, samples, markers = nilhmm.read_vcf_counts("your_data.vcf")
calls = nilhmm.introgression_hmm_counts(ref, alt, marker_dict, r=3e-5)
```

## Model Architecture

Three hidden states — B73 homozygote, heterozygote, donor homozygote — with
recombination/population-genetics transition probabilities (parameters `r`, `f_1`, `f_2`).
The emission is one of:

- **GT caller:** four observable genotype states (`0/0`, `0/1`, `1/1`, missing) via a 3×4
  categorical matrix that models genotyping-error sources (`germ`, `gert`, `nir`, `mr`).
- **Count caller:** allelic depths (ref, alt) per marker via a BetaBinomial with per-state
  expected alt fraction θ = [`err`, 0.5, 1−`err`] and concentration `conc`; zero-depth markers
  contribute a flat (uninformative) emission. `conc → ∞` approaches a Binomial.

Decoding is per-chromosome Viterbi. The count caller memoizes the BetaBinomial over the
distinct (depth, alt) pairs and batches the Viterbi across samples, so its cost scales with
sequencing depth (number of distinct count pairs), not the sample×marker product.

## Applications

- **QTL mapping** in introgression populations
- **Adaptive introgression** studies
- **Breeding program** donor segment tracking
- **Population genomics** introgression pattern analysis

## Package Structure

```
nilhmm/
├── nilhmm/           # Core package (GT + count callers, IO)
├── scripts/          # Analysis scripts
├── tests/            # Unit tests
├── calibration/      # BZea SNP50K calibration (sweeps, params, staging)
├── Implementation.md # Count-path design, calibration, results
├── docs/             # Documentation
└── setup.py          # Package metadata
```

## Citation

If you use this software in your research, please cite:

```
Zhong, T., Mullens, A., Morales, L., Swarts, K.L., Stafstrom, W.C., He, Y.,
Sermons, S.M., Yang, Q., Lopez-Zuniga, L.O., Rucker, E., Thomason, W.E.,
Nelson, R.J., Jamann, T., Balint-Kurti, P.J., & Holland, J.B. (2025).
A maize near-isogenic line population designed for gene discovery and
characterization of allelic effects. bioRxiv.
https://doi.org/10.1101/2025.01.29.635337
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome! Please submit issues and pull requests on GitHub.
