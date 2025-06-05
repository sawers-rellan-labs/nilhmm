# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NIL HMM is a Hidden Markov Model implementation for calling introgressions in Near Isogenic Lines (NILs) from low-coverage sequencing data. The project is specifically designed for maize genetics research, particularly the BZea population (B73 × Teosinte crosses).

## Current Project Structure

### Python Package (`nilhmm/`)

The main implementation is organized as a Python package:

- **`nilhmm/core.py`**: Core HMM algorithms and state inference
  - `introgression_hmm()`: Main HMM implementation with complex transition/emission matrices
  - `call_introgressions()`: High-level interface with coverage-specific parameter defaults
  
- **`nilhmm/io.py`**: VCF file parsing and data input/output utilities
  - `read_vcf()`: Converts VCF files to genotype matrices
  - `write_results()`: Outputs multiple file formats for analysis
  
- **`nilhmm/utils.py`**: Helper functions and data validation
  - Logging setup and utility functions
  
- **`nilhmm/grid_search.py`**: Parameter optimization routines
  - `optimize_parameters()`: Grid search for HMM parameter tuning

### Command-Line Scripts (`scripts/`)

- **`scripts/call_bzea_introgressions.py`**: Main command-line interface
  - Handles argument parsing and parameter configuration
  - Provides coverage-level parameter presets (low/medium/high)
  - Calls `nilhmm.call_introgressions()` with user parameters
  
- **`scripts/parameter_tuning.py`**: Grid search parameter optimization
- **`scripts/preprocess_vcf.py`**: VCF preprocessing utilities

### Documentation (`docs/`)

- **`docs/README.md`**: Detailed project documentation
- **`docs/jim_hmm_diagram.Rmd`**: R Markdown file for HMM architecture diagrams
- **`docs/package_structure.md`**: Package organization and design decisions

### Environment and Configuration

- **`envs/nilhmm.yml`**: Conda environment specification
- **`requirements.txt`**: Python package dependencies
- **`setup.py`**: Package installation configuration
- **`codemcp.toml`**: Development tooling configuration

### Additional Files

- **`python_hmm_diagram.Rmd`**: Python HMM implementation details and comparison with R version
- **`File_S11_callIntrogressions.py`**: Original HMM implementation (basis for `nilhmm/core.py`)

## Core HMM Architecture

The HMM uses a 3-state model with sophisticated mathematical foundations:

### Hidden States
- **State 0**: B73 homozygote (recurrent parent)
- **State 1**: Heterozygote (B73/Teosinte) 
- **State 2**: Donor homozygote (Teosinte)

### Observable States
- **0**: Major allele homozygote (minor allele count = 0)
- **1**: Heterozygote (minor allele count = 1) 
- **2**: Minor allele homozygote (minor allele count = 2)
- **3**: Missing genotype call

### Key Mathematical Components

**Transition Matrix**: Uses genetics-based probabilities incorporating recombination
- No recombination: P(S_t = i | S_{t-1} = i) = 1-r
- Recombination transitions weighted by expected state frequencies

**Emission Matrix**: Complex 3×4 matrix incorporating multiple error sources
- Genotyping errors on homozygotes vs heterozygotes
- Non-informative marker effects
- Missing data handling

## Key Parameters

The HMM uses 8 critical parameters that must be carefully tuned:

- **`nir`**: Non-informative marker rate (proportion of markers that don't distinguish ancestry)
- **`germ`**: SNP calling error rate on true homozygotes  
- **`gert`**: SNP calling error rate on true heterozygotes
- **`p`**: Proportion of homozygous errors resulting in heterozygous calls
- **`mr`**: Missing genotype call rate
- **`r`**: Recombination rate between adjacent markers
- **`f_1`**: Expected frequency of heterozygotes in population
- **`f_2`**: Expected frequency of donor homozygotes

### Coverage-Specific Parameter Defaults

The package provides optimized defaults for different sequencing coverage levels:

- **Low coverage (0.5-1x)**: `nir=0.02, germ=0.08, gert=0.15, mr=0.20`
- **Medium coverage (2-5x)**: `nir=0.01, germ=0.05, gert=0.10, mr=0.10`
- **High coverage (10x+)**: `nir=0.005, germ=0.02, gert=0.05, mr=0.05`

## Data Flow

1. **Input**: VCF files with GT (genotype) field
2. **Processing**: Convert to numeric matrix (samples × markers) via `nilhmm.io.read_vcf()`
3. **Chromosome splitting**: Process each chromosome 1-10 separately 
4. **HMM inference**: Apply Viterbi algorithm per sample per chromosome via `nilhmm.core.introgression_hmm()`
5. **Output**: Concatenate results across chromosomes and save via `nilhmm.io.write_results()`

## Common Development Commands

### Environment Setup
```bash
# Create conda environment
conda env create -f envs/nilhmm.yml
conda activate nilhmm

# Or install via pip
pip install -r requirements.txt
```

### Package Installation
```bash
# Install in development mode
pip install -e .
```

### Running the Main Pipeline
```bash
# Basic usage with automatic coverage detection
python scripts/call_bzea_introgressions.py input_data.vcf

# Specify coverage level for parameter defaults
python scripts/call_bzea_introgressions.py input_data.vcf --coverage low

# Custom parameters
python scripts/call_bzea_introgressions.py input_data.vcf \
    --nir 0.01 --germ 0.05 --gert 0.10 --p 0.5 --mr 0.15 --r 0.01

# Custom output prefix
python scripts/call_bzea_introgressions.py input_data.vcf -o custom_results
```

### Development Tools (via codemcp)
```bash
# Format code
codemcp format

# Run tests with coverage
codemcp test-coverage

# Lint code
codemcp lint

# Check git status
codemcp status
```

### Using the Package Directly
```python
import nilhmm

# High-level interface
results = nilhmm.call_introgressions(
    vcf_file="data.vcf",
    output_prefix="results",
    coverage_level="low"
)

# Lower-level access
geno_matrix, marker_dict, samples, markers = nilhmm.read_vcf("data.vcf")
calls = nilhmm.introgression_hmm(geno_matrix, marker_dict)
```

## Input Data Requirements

VCF files must have:
- Standard VCF format (v4.0+)
- GT (genotype) field in FORMAT column
- Chromosomes 1-10 (standard maize chromosomes)
- Sample names in header line
- Markers sorted by chromosome and position

## Output Files

The pipeline generates four standardized output files:
1. **`*_introgression_calls.txt`**: Raw numeric matrix
2. **`*_introgression_calls.csv`**: Labeled calls with headers
3. **`*_introgression_summary.csv`**: Per-sample summary statistics  
4. **`*_marker_info.csv`**: Marker details (chr, pos, alleles)

## Important Design Constraints

- **Chromosome Processing**: Each chromosome (1-10) must be processed separately due to linkage assumptions
- **Marker Ordering**: Input markers MUST be sorted by chromosome and position for HMM to work correctly
- **Missing Data Encoding**: Uses integer 3 for missing data, not NaN
- **Low Coverage Optimization**: Parameters specifically tuned for 0.5-1× coverage sequencing data

## Integration Notes

When extending this codebase:
- The HMM relies on `hmmlearn.MultinomialHMM` with custom transition/emission matrices
- Parameter values significantly impact results - use grid search for optimization
- All genotype matrices use samples as rows, markers as columns
- The `marker_dict` maps chromosome numbers to column indices in the genotype matrix
- Core functions are available through the package interface: `nilhmm.call_introgressions()`, `nilhmm.introgression_hmm()`

## Project Configuration

The project uses `codemcp.toml` for configuration management, which defines:
- Code formatting and linting standards (black, flake8, mypy)
- Testing commands (pytest with coverage)
- Git automation settings
- Python version requirements (>=3.8)
- File ignore patterns

The configuration enables automatic commits (`auto_commit = true`) and provides convenient development commands for maintaining code quality.
