# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NIL HMM is a Hidden Markov Model implementation for calling introgressions in Near Isogenic Lines (NILs) from low-coverage sequencing data. The project is specifically designed for maize genetics research, particularly the BZea population (B73 × Teosinte crosses).

## Core Architecture

### HMM Implementation (`File_S11_callIntrogressions.py`)

The core HMM uses a 3-state model with sophisticated mathematical foundations:

**Hidden States:**
- State 0: B73 homozygote (recurrent parent)
- State 1: Heterozygote (B73/Teosinte) 
- State 2: Donor homozygote (Teosinte)

**Observable States:**
- 0: Major allele homozygote (minor allele count = 0)
- 1: Heterozygote (minor allele count = 1) 
- 2: Minor allele homozygote (minor allele count = 2)
- 3: Missing genotype call

**Key Functions:**
- `call_intros()`: Main HMM implementation with complex transition/emission matrices
- `call_intros_one_chrom()`: Processes individual chromosomes using Viterbi algorithm

### VCF Processing Pipeline (`bzea_vcf_introgression_caller.py`)

Main command-line interface that:
- Parses VCF files (compressed or uncompressed)
- Converts genotype calls to numeric matrix format (0,1,2,3)
- Calls `call_intros()` function with user-specified parameters
- Outputs multiple file formats for downstream analysis

### Parameter Optimization (`File_S14.Step_2_nNIL_introgression_calls_from_chip_data.py`)

Grid search implementation for parameter tuning using high-quality chip data as ground truth. Demonstrates the full parameter optimization workflow.

## Key Parameters

The HMM uses 8 critical parameters that must be carefully tuned:

- `nir`: Non-informative marker rate (proportion of markers that don't distinguish ancestry)
- `germ`: SNP calling error rate on true homozygotes  
- `gert`: SNP calling error rate on true heterozygotes
- `p`: Proportion of homozygous errors resulting in heterozygous calls
- `mr`: Missing genotype call rate
- `r`: Recombination rate between adjacent markers
- `f_1`: Expected frequency of heterozygotes in population
- `f_2`: Expected frequency of donor homozygotes

## Data Flow

1. **Input**: VCF files with GT (genotype) field
2. **Processing**: Convert to numeric matrix (samples × markers)
3. **Chromosome splitting**: Process each chromosome 1-10 separately 
4. **HMM inference**: Apply Viterbi algorithm per sample per chromosome
5. **Output**: Concatenate results across chromosomes

## Common Development Commands

### Running the Main Pipeline
```bash
# Basic usage
python bzea_vcf_introgression_caller.py input_data.vcf

# With custom parameters
python bzea_vcf_introgression_caller.py input_data.vcf \
    --nir 0.01 --germ 0.05 --gert 0.10 --p 0.5 --mr 0.15 --r 0.01

# With custom output prefix
python bzea_vcf_introgression_caller.py input_data.vcf -o custom_results
```

### Dependencies
Install required packages:
```bash
pip install numpy pandas hmmlearn argparse pathlib
```

### Testing Data Format
VCF files must have:
- Standard VCF format (v4.0+)
- GT (genotype) field in FORMAT column
- Chromosomes 1-10 (standard maize chromosomes)
- Sample names in header line

## Output Files

The pipeline generates four standardized output files:
1. `*_introgression_calls.txt`: Raw numeric matrix
2. `*_introgression_calls.csv`: Labeled calls with headers
3. `*_introgression_summary.csv`: Per-sample summary statistics  
4. `*_marker_info.csv`: Marker details (chr, pos, alleles)

## Mathematical Implementation Notes

### Transition Matrix
Uses genetics-based probabilities incorporating recombination:
- No recombination: P(S_t = i | S_{t-1} = i) = 1-r
- Recombination transitions weighted by expected state frequencies

### Emission Matrix  
Complex 3×4 matrix incorporating multiple error sources:
- Genotyping errors on homozygotes vs heterozygotes
- Non-informative marker effects
- Missing data handling

## Important Design Constraints

- **Chromosome Processing**: Each chromosome (1-10) must be processed separately due to linkage assumptions
- **Marker Ordering**: Input markers MUST be sorted by chromosome and position for HMM to work correctly
- **Missing Data Encoding**: Uses integer 3 for missing data, not NaN
- **Low Coverage Optimization**: Parameters specifically tuned for 0.8× coverage sequencing data

## Integration Notes

When extending this codebase:
- The HMM relies on `hmmlearn.MultinomialHMM` with custom transition/emission matrices
- Parameter values significantly impact results - use grid search for optimization
- All genotype matrices use samples as rows, markers as columns
- The `marker_dict` maps chromosome numbers to column indices in the genotype matrix