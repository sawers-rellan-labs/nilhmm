#!/usr/bin/env python3
"""
Parameter optimization script for nilhmm.
"""

import argparse
import sys
from pathlib import Path
import nilhmm
from nilhmm.utils import setup_logging


def main():
    parser = argparse.ArgumentParser(
        description="Optimize HMM parameters using grid search"
    )

    parser.add_argument(
        "vcf_file",
        help="Input VCF file path"
    )

    parser.add_argument(
        "-o", "--output",
        default="parameter_optimization",
        help="Output prefix (default: parameter_optimization)"
    )

    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO)"
    )

    args = parser.parse_args()

    # Set up logging
    setup_logging(args.log_level)

    # Check if input file exists
    if not Path(args.vcf_file).exists():
        print(f"Error: VCF file '{args.vcf_file}' not found")
        sys.exit(1)

    try:
        print(f"Loading VCF file: {args.vcf_file}")

        # Read VCF data
        geno_matrix, marker_dict, sample_names, marker_info = nilhmm.read_vcf(args.vcf_file)

        print(f"Loaded {len(sample_names)} samples and {len(marker_info)} markers")
        print("Starting parameter optimization...")

        # Run grid search
        results_df = nilhmm.optimize_parameters(
            geno_matrix=geno_matrix,
            marker_dict=marker_dict,
            output_file=f"{args.output}_grid_search_results.csv"
        )

        # Select best parameters
        best_params = nilhmm.grid_search.select_best_parameters(
            results_df,
            criteria="balanced"
        )

        print(f"\nOptimization completed!")
        print(f"Best parameters: {best_params}")
        print(f"Results saved to: {args.output}_grid_search_results.csv")

    except Exception as e:
        print(f"Error during optimization: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
