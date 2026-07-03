#!/usr/bin/env python3
"""
Main script for calling introgressions in Bzea NIL population.
"""

import argparse
import sys
from pathlib import Path
import nilhmm
from nilhmm.utils import setup_logging


def main():
    parser = argparse.ArgumentParser(
        description="Call introgressions in Bzea NIL population using HMM"
    )

    # Required arguments
    parser.add_argument(
        "vcf_file",
        help="Input VCF file path"
    )

    # Optional arguments
    parser.add_argument(
        "-o", "--output",
        default="bzea_introgressions",
        help="Output prefix (default: bzea_introgressions)"
    )

    parser.add_argument(
        "--coverage",
        default="low",
        choices=["low", "medium", "high"],
        help="Sequencing coverage level for parameter defaults (default: low)"
    )

    # HMM parameters (optional overrides)
    parser.add_argument(
        "--nir",
        type=float,
        help="Non-informative rate (default: coverage-dependent)"
    )

    parser.add_argument(
        "--germ",
        type=float,
        help="Error rate on homozygotes (default: coverage-dependent)"
    )

    parser.add_argument(
        "--gert",
        type=float,
        help="Error rate on heterozygotes (default: coverage-dependent)"
    )

    parser.add_argument(
        "--p",
        type=float,
        default=0.5,
        help="Proportion of homozygous errors as het calls (default: 0.5)"
    )

    parser.add_argument(
        "--mr",
        type=float,
        help="Missing rate (default: coverage-dependent)"
    )

    parser.add_argument(
        "--r",
        type=float,
        default=0.01,
        help="Recombination rate between markers (default: 0.01)"
    )

    parser.add_argument(
        "--f1",
        type=float,
        default=0.25,
        help="Expected frequency of heterozygotes (default: 0.25)"
    )

    parser.add_argument(
        "--f2",
        type=float,
        default=0.05,
        help="Expected frequency of donor homozygotes (default: 0.05)"
    )

    # Other options
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

    # Prepare HMM parameters
    hmm_params = {}

    # Add user-specified parameters
    if args.nir is not None:
        hmm_params["nir"] = args.nir
    if args.germ is not None:
        hmm_params["germ"] = args.germ
    if args.gert is not None:
        hmm_params["gert"] = args.gert
    if args.mr is not None:
        hmm_params["mr"] = args.mr

    hmm_params.update({
        "p": args.p,
        "r": args.r,
        "f_1": args.f1,
        "f_2": args.f2
    })

    try:
        # Run introgression calling
        print(f"Processing VCF file: {args.vcf_file}")
        print(f"Coverage level: {args.coverage}")
        print(f"Output prefix: {args.output}")

        results = nilhmm.call_introgressions(
            vcf_file=args.vcf_file,
            output_prefix=args.output,
            coverage_level=args.coverage,
            **hmm_params
        )

        # Print summary
        n_samples, n_markers = results["calls"].shape
        print(f"\nAnalysis completed successfully!")
        print(f"Processed {n_samples} samples and {n_markers} markers")
        print(f"Results saved with prefix: {args.output}")

    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
