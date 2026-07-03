# Legacy Python implementation (frozen)

This directory holds the **original Python implementation** of nilHMM. It is
**frozen** and kept for provenance and reproducibility only — it receives no
further development. All active work happens in the R + Rcpp package at the
repository root.

## Why it's here

nilHMM began as a Python package built on Jim Holland's nNIL introgression-calling
methodology (<https://github.com/ncsumaize/nNIL>). It has since been reimplemented
as a unified R + Rcpp package — a single duration-aware 3-state (REF/HET/ALT) HMM
engine with swappable emission and duration layers that expresses the `nnil`,
`rtiger`, `binhmm`, and `atlas` callers. See the root [`README.md`](../../README.md).

The R port reproduces the Python count caller bit-for-bit on the golden-slice
fixtures, so this code is no longer on the critical path. It is preserved rather
than deleted so the original algorithm and its calibration remain retrievable.

## What's here

```
legacy/python/
├── nilhmm/        # the Python package (core.py, io.py, grid_search.py, utils.py)
├── scripts/       # command-line entry points (call_bzea_introgressions.py, ...)
├── calibration/   # Python KS-sweep / calibration scripts (CSV outputs stay in ../../calibration/)
├── envs/          # conda environment for the Python package
├── tests/         # Python unit tests (pytest)
├── setup.py
└── requirements.txt
```

## Retrieving the last active state

The final commit of the Python implementation before retirement is tagged
[`python-final`](https://github.com/sawers-rellan-labs/nilhmm/releases/tag/python-final):

```bash
git checkout python-final     # the tree as it was, Python at the repo root
```

## Origin

- **nNIL methodology:** Jim Holland, <https://github.com/ncsumaize/nNIL>
- **Preprint:** Zhong, T. *et al.* (2025). *A maize near-isogenic line population
  designed for gene discovery and characterization of allelic effects.* bioRxiv.
  <https://doi.org/10.1101/2025.01.29.635337>
