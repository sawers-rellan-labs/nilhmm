#!/usr/bin/env python3
"""
VCF preprocessing script for nilhmm.
"""

import argparse
import sys
from pathlib import Path
import subprocess
import logging


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess VCF files for nilhmm analysis"
    )

    parser.add_argument(
        "input_vcf",
        help="Input VCF file path"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output VCF file path"
    )

    parser.add_argument(
        "--min-qual",
        type=float,
        default=20.0,
        help="Minimum variant quality (default: 20.0)"
    )

    parser.add_argument(
        "--min-depth",
        type=int,
        default=5,
        help="Minimum read depth (default: 5)"
    )

    parser.add_argument(
        "--max-missing",
        type=float,
        default=0.5,
        help="Maximum missing data rate (default: 0.5)"
    )

    parser.add_argument(
        "--chromosomes",
        nargs="+",
        type=int,
        default=list(range(1, 11)),
        help="Chromosomes to include (default: 1-10)"
    )

    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(level=logging.INFO)

    # Check if input file exists
    if not Path(args.input_vcf).exists():
        print(f"Error: VCF file '{args.input_vcf}' not found")
        sys.exit(1)

    try:
        print(f"Preprocessing VCF file: {args.input_vcf}")

        # Build bcftools filter command
        chrom_filter = ",".join([str(c) for c in args.chromosomes])

        cmd = [
            "bcftools", "view",
            "-r", chrom_filter,  # Select chromosomes
            "-q", str(args.min_qual),  # Quality filter
            "-e", f"F_MISSING > {args.max_missing}",  # Missing data filter
            "-O", "z",  # Output compressed VCF
            "-o", args.output,
            args.input_vcf
        ]

        print(f"Running: {' '.join(cmd)}")

        # Run bcftools
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Error running bcftools: {result.stderr}")
            sys.exit(1)

        # Index the output file
        index_cmd = ["bcftools", "index", args.output]
        subprocess.run(index_cmd, check=True)

        print(f"Preprocessing completed!")
        print(f"Filtered VCF saved to: {args.output}")

    except subprocess.CalledProcessError as e:
        print(f"Error during preprocessing: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
