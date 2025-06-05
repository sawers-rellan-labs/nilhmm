"""I/O functions for VCF reading and result writing."""

import numpy as np
import pandas as pd
import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import logging


def parse_vcf_header(vcf_file: str) -> List[str]:
    """Parse VCF header to extract sample names."""
    samples = []

    opener = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'

    with opener(vcf_file, mode) as f:
        for line in f:
            if line.startswith('#CHROM'):
                header_fields = line.strip().split('\t')
                samples = header_fields[9:]  # Sample names start from column 10
                break
            elif not line.startswith('#'):
                raise ValueError("No header line found in VCF file")

    return samples


def read_vcf(
    vcf_file: str,
    chromosomes: Optional[List[int]] = None
) -> Tuple[np.ndarray, Dict[int, List[int]], List[str], pd.DataFrame]:
    """
    Convert VCF file to genotype matrix suitable for HMM.

    Parameters:
    -----------
    vcf_file : str
        Path to VCF file
    chromosomes : List[int], optional
        List of chromosomes to process (default: 1-10 for maize)

    Returns:
    --------
    Tuple containing:
        - geno_matrix : np.ndarray (samples x markers)
        - marker_dict : Dict[int, List[int]] (chromosome -> marker indices)
        - sample_names : List[str]
        - marker_info : pd.DataFrame
    """
    if chromosomes is None:
        chromosomes = list(range(1, 11))  # Maize chromosomes 1-10

    # Get sample names from header
    sample_names = parse_vcf_header(vcf_file)

    # Lists to store data
    genotypes = []
    marker_info = []

    opener = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'

    logging.info(f"Reading VCF file: {vcf_file}")
    logging.info(f"Found {len(sample_names)} samples")

    with opener(vcf_file, mode) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')

            chrom = fields[0]
            pos = int(fields[1])
            marker_id = fields[2]
            ref = fields[3]
            alt = fields[4]

            # Skip if chromosome not in target list
            try:
                chrom_num = int(chrom.replace('chr', ''))
                if chrom_num not in chromosomes:
                    continue
            except ValueError:
                continue

            # Store marker information
            marker_info.append({
                'CHROM': chrom_num,
                'POS': pos,
                'ID': marker_id,
                'REF': ref,
                'ALT': alt
            })

            # Parse genotypes
            format_field = fields[8]
            gt_index = format_field.split(':').index('GT')

            sample_genos = []
            for sample_field in fields[9:]:
                gt_call = sample_field.split(':')[gt_index]

                # Convert genotype to numeric
                if gt_call in ['./.', '.', '.|.']:
                    geno_code = 3  # Missing
                elif gt_call in ['0/0', '0|0']:
                    geno_code = 0  # Homozygous reference
                elif gt_call in ['0/1', '1/0', '0|1', '1|0']:
                    geno_code = 1  # Heterozygous
                elif gt_call in ['1/1', '1|1']:
                    geno_code = 2  # Homozygous alternate
                else:
                    geno_code = 3  # Treat unknown as missing

                sample_genos.append(geno_code)

            genotypes.append(sample_genos)

            if line_num % 10000 == 0:
                logging.info(f"Processed {line_num} variants...")

    # Convert to numpy array (samples x markers)
    geno_matrix = np.array(genotypes).T
    marker_df = pd.DataFrame(marker_info)

    logging.info(f"Final genotype matrix shape: {geno_matrix.shape}")
    logging.info(f"Samples: {geno_matrix.shape[0]}, Markers: {geno_matrix.shape[1]}")

    # Create marker dictionary by chromosome
    marker_dict = {}
    for chrom in chromosomes:
        chrom_indices = marker_df[marker_df['CHROM'] == chrom].index.tolist()
        marker_dict[chrom] = chrom_indices
        logging.info(f"Chromosome {chrom}: {len(chrom_indices)} markers")

    return geno_matrix, marker_dict, sample_names, marker_df


def write_results(results: Dict, output_prefix: str) -> None:
    """
    Write introgression calling results to files.

    Parameters:
    -----------
    results : Dict
        Results dictionary from call_introgressions
    output_prefix : str
        Prefix for output files
    """
    calls = results["calls"]
    sample_names = results["sample_names"]
    marker_info = results["marker_info"]
    parameters = results["parameters"]

    # Save raw calls matrix
    np.savetxt(
        f"{output_prefix}_introgression_calls.txt",
        calls,
        fmt='%d',
        delimiter='\t'
    )

    # Create and save labeled DataFrame
    marker_names = [
        f"{row['CHROM']}_{row['POS']}"
        for _, row in marker_info.iterrows()
    ]

    calls_df = pd.DataFrame(
        calls,
        index=sample_names,
        columns=marker_names
    )
    calls_df.to_csv(f"{output_prefix}_introgression_calls.csv")

    # Calculate summary statistics
    summary_stats = []
    for i, sample in enumerate(sample_names):
        sample_calls = calls[i, :]
        n_total = len(sample_calls)
        n_b73 = np.sum(sample_calls == 0)
        n_het = np.sum(sample_calls == 1)
        n_donor = np.sum(sample_calls == 2)

        summary_stats.append({
            'Sample': sample,
            'Total_markers': n_total,
            'B73_homoz': n_b73,
            'Heterozygous': n_het,
            'Donor_homoz': n_donor,
            'Pct_B73': (n_b73 / n_total) * 100,
            'Pct_Het': (n_het / n_total) * 100,
            'Pct_Donor': (n_donor / n_total) * 100
        })

    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(f"{output_prefix}_introgression_summary.csv", index=False)

    # Save marker information
    marker_info.to_csv(f"{output_prefix}_marker_info.csv", index=False)

    # Save parameters used
    params_df = pd.DataFrame([parameters])
    params_df.to_csv(f"{output_prefix}_parameters.csv", index=False)

    logging.info(f"Results saved:")
    logging.info(f"  - {output_prefix}_introgression_calls.txt (raw matrix)")
    logging.info(f"  - {output_prefix}_introgression_calls.csv (labeled)")
    logging.info(f"  - {output_prefix}_introgression_summary.csv (summary stats)")
    logging.info(f"  - {output_prefix}_marker_info.csv (marker details)")
    logging.info(f"  - {output_prefix}_parameters.csv (HMM parameters)")
