# Write imputed per-marker genotypes as a VCF (LB-Impute-style deliverable)

Optional, decoupled from the engine: turns a per-marker state table
(from
[`call_states()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_states.md),
typically `caller = "lbimpute"`) into an imputed biallelic VCF –
LB-Impute's native output form (missing filled, false homozygotes
corrected), as opposed to nilHMM's default segment schema. State -\>
genotype is `0 -> 0/0` (REF-hom), `1 -> 0/1` (het), `2 -> 1/1`
(ALT-hom); any marker x sample absent from `states` is written `./.`.
This never re-reads the input VCF; supply per-marker `REF`/`ALT` alleles
via `markers` if you want the real alleles rather than the `N`
placeholder.

## Usage

``` r
write_vcf_impute(states, path, markers = NULL, ref = "N", alt = "N")
```

## Arguments

- states:

  A per-marker state table with columns `name, chr, pos, state` (0/1/2),
  e.g. from `call_states(..., caller = "lbimpute")`.

- path:

  Output `.vcf` path.

- markers:

  Optional data.frame keyed by `chr, pos` carrying `ref`, `alt` (and
  optional `id`) alleles to emit; when `NULL`, REF/ALT default to `N`.

- ref, alt:

  Fallback single-character REF/ALT alleles when `markers` is `NULL`
  (default `"N"`).

## Value

`path`, invisibly.

## Examples

``` r
st <- data.frame(name = c("NIL1", "NIL1", "NIL2", "NIL2"),
                 chr = 1L, pos = c(1e5, 2e5, 1e5, 2e5),
                 state = c(0L, 2L, 1L, 0L))
write_vcf_impute(st, tempfile(fileext = ".vcf"))
```
