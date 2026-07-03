# nilHMM: duration-aware HMM ancestry caller

A unified R + Rcpp engine for calling REF/HET/ALT ancestry in NILs from
sequencing data. One 3-state HMM with swappable emission
([`emission_count()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_count.md),
[`emission_gt()`](https://sawers-rellan-labs.github.io/nilhmm/reference/emission_gt.md))
and duration
([`duration_geometric()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_geometric.md),
[`duration_rigidity()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_rigidity.md),
[`duration_hsmm()`](https://sawers-rellan-labs.github.io/nilhmm/reference/duration_hsmm.md))
layers expresses the `nnil` and `rtiger` callers (plus the `binhmm`
bin/cluster caller). See `REFACTOR_R_PACKAGE.md` for the design.

## Details

The package is **data-agnostic**: no hardcoded paths, sample lists, or
mount locations. Functions take `(data, params)` and return calls;
pipeline scripts (not part of the package) own which files/samples and
where outputs go.

## See also

Useful links:

- <https://github.com/sawers-rellan-labs/nilhmm>

- <https://sawers-rellan-labs.github.io/nilhmm/>

- Report bugs at <https://github.com/sawers-rellan-labs/nilhmm/issues>

## Author

**Maintainer**: Fausto Rodriguez Zapata <thegemmalab@ncsu.edu>

Authors:

- Fausto Rodriguez Zapata <thegemmalab@ncsu.edu>
