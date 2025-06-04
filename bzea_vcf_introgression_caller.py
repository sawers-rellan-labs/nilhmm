#!/usr/bin/env python3
"""
VCF Introgression Caller for BZea Population

Reads VCF file from low coverage sequencing (0.8x) of BZea maize NIL population
and calls introgressions using HMM-based approach.

Author: Based on File_S11_callIntrogressions.py by jholland
"""

import numpy as np
import pandas as pd
import argparse
import sys
from pathlib import Path
import gzip
from File_S11_callIntrogressions import call_intros


def parse_vcf_header(vcf_file):
    """Parse VCF header to extract sample names and metadata."""
    samples = []
    
    if vcf_file.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(vcf_file, mode) as f:
        for line in f:
            if line.startswith('#CHROM'):
                header_fields = line.strip().split('\t')
                samples = header_fields[9:]  # Sample names start from column 10
                break
            elif not line.startswith('#'):
                raise ValueError("No header line found in VCF file")
    
    return samples


def vcf_to_genotype_matrix(vcf_file, chromosomes=None):
    """
    Convert VCF file to genotype matrix suitable for call_intros function.
    
    Parameters:
    -----------
    vcf_file : str
        Path to VCF file
    chromosomes : list, optional
        List of chromosomes to process (default: 1-10 for maize)
    
    Returns:
    --------
    geno_matrix : numpy.ndarray
        Genotype matrix (samples x markers) with values:
        0 = homozygous reference
        1 = heterozygous
        2 = homozygous alternate
        3 = missing data
    marker_dict : dict
        Dictionary mapping chromosome to marker indices
    sample_names : list
        List of sample names
    marker_info : pandas.DataFrame
        DataFrame with marker information (CHROM, POS, ID, REF, ALT)
    """
    if chromosomes is None:
        chromosomes = list(range(1, 11))  # Maize chromosomes 1-10
    
    # Get sample names from header
    sample_names = parse_vcf_header(vcf_file)
    
    # Lists to store data
    genotypes = []
    marker_info = []
    
    if vcf_file.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    print(f"Reading VCF file: {vcf_file}")
    print(f"Found {len(sample_names)} samples")
    
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
                print(f"Processed {line_num} variants...")
    
    # Convert to numpy array (samples x markers)
    geno_matrix = np.array(genotypes).T
    marker_df = pd.DataFrame(marker_info)
    
    print(f"Final genotype matrix shape: {geno_matrix.shape}")
    print(f"Samples: {geno_matrix.shape[0]}, Markers: {geno_matrix.shape[1]}")
    
    # Create marker dictionary by chromosome
    marker_dict = {}
    for chrom in chromosomes:
        chrom_indices = marker_df[marker_df['CHROM'] == chrom].index.tolist()
        marker_dict[chrom] = chrom_indices
        print(f"Chromosome {chrom}: {len(chrom_indices)} markers")
    
    return geno_matrix, marker_dict, sample_names, marker_df


def call_bzea_introgressions(vcf_file, output_prefix, 
                           nir=0.01, germ=0.05, gert=0.10, p=0.5, mr=0.15, 
                           r=0.01, f_1=0.25, f_2=0.05):
    """
    Call introgressions in BZea population from VCF file.
    
    Parameters optimized for low coverage (0.8x) sequencing data.
    
    Parameters:
    -----------
    vcf_file : str
        Path to input VCF file
    output_prefix : str
        Prefix for output files
    nir : float
        Non-informative rate (default: 0.01)
    germ : float
        SNP calling error rate on true introgression homozygotes (default: 0.05)
    gert : float
        SNP calling error rate on true introgression heterozygotes (default: 0.10)
    p : float
        Proportion of homozygous SNP call errors resulting in het call (default: 0.5)
    mr : float
        Missing call rate (default: 0.15, higher for low coverage)
    r : float
        Expected recombination rate between adjacent markers (default: 0.01)
    f_1 : float
        Expected frequency of heterozygotes (default: 0.25)
    f_2 : float
        Expected frequency of homozygous introgressions (default: 0.05)
    """
    print("BZea VCF Introgression Caller")
    print("=" * 40)
    
    # Parse VCF file
    geno_matrix, marker_dict, sample_names, marker_df = vcf_to_genotype_matrix(vcf_file)
    
    # Call introgressions using HMM
    print("\nCalling introgressions using HMM...")
    print(f"Parameters: nir={nir}, germ={germ}, gert={gert}, p={p}")
    print(f"            mr={mr}, r={r}, f_1={f_1}, f_2={f_2}")
    
    introgression_calls = call_intros(
        geno=geno_matrix,
        marker_dict=marker_dict,
        nir=nir,
        germ=germ,
        gert=gert,
        p=p,
        mr=mr,
        r=r,
        f_1=f_1,
        f_2=f_2,
        return_calls=True
    )
    
    print(f"Introgression calls shape: {introgression_calls.shape}")
    
    # Save results
    print("\nSaving results...")
    
    # Save raw calls matrix
    np.savetxt(f"{output_prefix}_introgression_calls.txt", 
               introgression_calls, fmt='%d', delimiter='\t')
    
    # Create and save summary DataFrame
    calls_df = pd.DataFrame(
        introgression_calls,
        index=sample_names,
        columns=[f"{marker_df.iloc[i]['CHROM']}_{marker_df.iloc[i]['POS']}" 
                for i in range(len(marker_df))]
    )
    calls_df.to_csv(f"{output_prefix}_introgression_calls.csv")
    
    # Calculate summary statistics
    summary_stats = []
    for i, sample in enumerate(sample_names):
        sample_calls = introgression_calls[i, :]
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
    marker_df.to_csv(f"{output_prefix}_marker_info.csv", index=False)
    
    print(f"Results saved to:")
    print(f"  - {output_prefix}_introgression_calls.txt (raw matrix)")
    print(f"  - {output_prefix}_introgression_calls.csv (labeled)")
    print(f"  - {output_prefix}_introgression_summary.csv (summary stats)")
    print(f"  - {output_prefix}_marker_info.csv (marker details)")
    
    return introgression_calls, summary_df, marker_df


def main():
    parser = argparse.ArgumentParser(
        description="Call introgressions in BZea NIL population from VCF file"
    )
    
    parser.add_argument("vcf_file", help="Input VCF file path")
    parser.add_argument("-o", "--output", default="bzea_introgressions",
                       help="Output prefix (default: bzea_introgressions)")
    
    # HMM parameters
    parser.add_argument("--nir", type=float, default=0.01,
                       help="Non-informative rate (default: 0.01)")
    parser.add_argument("--germ", type=float, default=0.05,
                       help="Error rate on true homozygotes (default: 0.05)")
    parser.add_argument("--gert", type=float, default=0.10,
                       help="Error rate on true heterozygotes (default: 0.10)")
    parser.add_argument("--p", type=float, default=0.5,
                       help="Prop. of homoz errors as het calls (default: 0.5)")
    parser.add_argument("--mr", type=float, default=0.15,
                       help="Missing call rate (default: 0.15)")
    parser.add_argument("--r", type=float, default=0.01,
                       help="Recombination rate between markers (default: 0.01)")
    parser.add_argument("--f1", type=float, default=0.25,
                       help="Expected frequency of heterozygotes (default: 0.25)")
    parser.add_argument("--f2", type=float, default=0.05,
                       help="Expected freq. of homoz introgressions (default: 0.05)")
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not Path(args.vcf_file).exists():
        print(f"Error: VCF file '{args.vcf_file}' not found")
        sys.exit(1)
    
    # Run introgression calling
    try:
        call_bzea_introgressions(
            vcf_file=args.vcf_file,
            output_prefix=args.output,
            nir=args.nir,
            germ=args.germ,
            gert=args.gert,
            p=args.p,
            mr=args.mr,
            r=args.r,
            f_1=args.f1,
            f_2=args.f2
        )
        print("\nIntrogression calling completed successfully!")
        
    except Exception as e:
        print(f"Error during introgression calling: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
