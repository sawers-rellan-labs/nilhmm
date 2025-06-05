# NILHMM Python Package Structure (R-like)

```
nilhmm/
├── nilhmm/                     # Main package directory (like R/)
│   ├── __init__.py            # Package imports (like NAMESPACE)
│   ├── core.py                # Main HMM functions
│   ├── io.py                  # VCF reading/writing
│   ├── grid_search.py         # Parameter optimization
│   └── utils.py               # Helper functions
├── scripts/                   # Executable scripts (like inst/scripts/)
│   ├── call_bzea_introgressions.py
│   ├── preprocess_vcf.py
│   └── parameter_tuning.py
├── tests/                     # Unit tests
│   ├── test_core.py
│   ├── test_io.py
│   └── test_data/
├── docs/                      # Documentation
│   └── README.md
├── setup.py                   # Package metadata (like DESCRIPTION)
├── requirements.txt           # Dependencies (like Imports:)
├── README.md                  # Package description
├── LICENSE                    # License file
└── .gitignore                 # Git ignore patterns

```

## What `call_bzea_introgressions.py` would be:

This would be your **main analysis script** (equivalent to an R script that calls your package functions):

```python
#!/usr/bin/env python3
"""
Main script for calling introgressions in Bzea population.
Equivalent to a analysis script that uses your R package.
"""

import argparse
from pathlib import Path
import nilhmm

def main():
    parser = argparse.ArgumentParser(description="Call introgressions in Bzea NIL population")
    parser.add_argument("vcf_file", help="Input VCF file")
    parser.add_argument("-o", "--output", default="bzea_results", help="Output prefix")
    parser.add_argument("--coverage", default="low", choices=["low", "medium", "high"])
    
    args = parser.parse_args()
    
    # Use your package functions
    results = nilhmm.call_introgressions(
        vcf_file=args.vcf_file,
        output_prefix=args.output,
        coverage_level=args.coverage
    )
    
    print(f"Analysis complete. Results saved to {args.output}_*")

if __name__ == "__main__":
    main()
```

## Python Development Workflow (devtools equivalent):

### 1. **Package Setup** (like `usethis::create_package()`)
```bash
# Create the basic structure
mkdir nilhmm
cd nilhmm
touch setup.py requirements.txt README.md LICENSE .gitignore
mkdir nilhmm scripts tests docs
touch nilhmm/__init__.py
```

### 2. **Development Installation** (like `devtools::load_all()`)
```bash
# Install in development mode - changes are immediately available
pip install -e .
```

### 3. **Testing** (like `devtools::test()`)
```bash
# Run tests
python -m pytest tests/
```

### 4. **Package Check** (like `devtools::check()`)
```bash
# Check package can be built
python setup.py check
python setup.py sdist  # builds source distribution
```

### 5. **Documentation** (like `devtools::document()`)
```bash
# Generate documentation (if using sphinx)
sphinx-build docs/ docs/_build/
```

### 6. **Installation from GitHub** (like `devtools::install_github()`)
```bash
pip install git+https://github.com/sawers-rellan-labs/nilhmm.git
```

## Key Files You Need:

### `setup.py` (equivalent to DESCRIPTION):
```python
from setuptools import setup, find_packages

setup(
    name="nilhmm",
    version="0.1.0",
    author="Fausto Rodriguez-Zapata",
    author_email="faustovrz@gmail.com", 
    description="NIL HMM for maize genetics introgression analysis",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.18.0",
        "pandas>=1.0.0", 
        "hmmlearn>=0.2.0",
        "scikit-learn"
    ],
    python_requires=">=3.8",
    entry_points={
        'console_scripts': [
            'call-bzea-introgressions=scripts.call_bzea_introgressions:main',
        ],
    },
)
```

### `nilhmm/__init__.py` (equivalent to NAMESPACE):
```python
"""NILHMM: HMM-based introgression calling for NIL populations"""

from .core import call_introgressions, IntrogressionHMM
from .io import read_vcf, write_results  
from .grid_search import optimize_parameters

__version__ = "0.1.0"
__all__ = ["call_introgressions", "IntrogressionHMM", "read_vcf", "write_results", "optimize_parameters"]
```

## Workflow Summary:

1. **Development**: `pip install -e .` (changes immediately available)
2. **Testing**: `python -m pytest`
3. **Check**: `python setup.py check`
4. **Install from GitHub**: `pip install git+https://github.com/user/nilhmm.git`

This is much simpler than my initial suggestion and closely mirrors your R package development workflow. The `scripts/` directory contains your analysis scripts that use the package, just like R scripts that load and use your R package functions.

Would you like me to help set up the basic files to get started with this structure?
