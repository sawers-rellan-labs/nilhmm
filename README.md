# NILHMM: Hidden Markov Model for NIL Introgression Analysis

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

NILHMM implements Hidden Markov Model approaches for calling introgressions in Near Isogenic Lines (NILs) from sequencing data. The package is optimized for maize genetics applications, particularly low-coverage sequencing of populations like BZea (B73 × Teosinte).

Based upon Jim Holland's methodology:
- Repository: https://github.com/ncsumaize/nNIL
- Preprint: Zhong, T. et al. (2025). A maize near-isogenic line population designed for gene discovery and characterization of allelic effects. bioRxiv. https://doi.org/10.1101/2025.01.29.635337

## Key Features

- **Individual SNP resolution** with sophisticated genotyping error modeling
- **VCF file processing** with support for compressed formats
- **Parameter optimization** tools for different sequencing conditions
- **Comprehensive output** including calls, summaries, and marker information
- **Optimized for low coverage** sequencing (0.5-1.0× coverage)

## Installation

### From Source
```bash
git clone https://github.com/sawers-rellan-labs/nilhmm.git
cd nilhmm
pip install -e .
```

### Dependencies
```bash
pip install numpy pandas hmmlearn scikit-learn
```

## Quick Start

### Command Line
```bash
# Basic introgression calling
call-bzea-introgressions your_data.vcf -o results

# With custom parameters for your data
call-bzea-introgressions your_data.vcf -o results --mr 0.2 --germ 0.08
```

### Python API
```python
import nilhmm

# Load VCF and call introgressions
results = nilhmm.call_introgressions(
    vcf_file="your_data.vcf",
    output_prefix="results",
    coverage_level="low"
)

# Access results
calls = results.introgression_calls
summary = results.summary_stats
```

## Model Architecture

The HMM uses three hidden states (B73 homozygote, heterozygote, donor homozygote) and four observable genotype states (0/0, 0/1, 1/1, missing). Transition probabilities incorporate recombination rates and population genetics expectations, while emission probabilities model various sources of genotyping error.

## Applications

- **QTL mapping** in introgression populations
- **Adaptive introgression** studies
- **Breeding program** donor segment tracking
- **Population genomics** introgression pattern analysis

## Package Structure

```
nilhmm/
├── nilhmm/           # Core package
├── scripts/          # Analysis scripts
├── tests/            # Unit tests
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
