# Per-site genotype caller with a swappable prior

A closed-form GATK-style diploid genotype call from allelic read counts,
decided independently per (marker, sample) – no linkage model. ALT =
donor (teosinte) = dosage 2. Given per-read error `error` and counts
`(n_ref, n_alt)`, the three genotype log-likelihoods (GL) are \$\$L(0) =
(1-\epsilon)^{n\_{ref}}\\\epsilon^{n\_{alt}}\$\$ \$\$L(2) =
\epsilon^{n\_{ref}}\\(1-\epsilon)^{n\_{alt}}\$\$ \$\$L(1) =
0.5^{\\n\_{ref}+n\_{alt}}\$\$ and the call is `argmax_g L(g) * pi(g)`
over the prior `pi` (the MAP / argmax-GP estimate). All arithmetic is in
log space and vectorized over the whole matrix (no per-cell loop).
Markers with zero total depth carry no signal and are returned as `NA`.

## Usage

``` r
call_gt(
  n_ref,
  n_alt,
  prior = "hwe",
  af = NULL,
  error = 0.01,
  return = c("call", "dosage", "post")
)
```

## Arguments

- n_ref, n_alt:

  Reference / alternate read counts. Integer (or numeric) matrices
  `markers x samples`, or plain vectors (treated as one sample =
  `markers x 1`). Same shape required.

- prior:

  The single-locus prior: `"hwe"` (default), `"flat"`, or a length-3
  numeric vector `c(f_REF, f_HET, f_ALT)` (fixed genome-wide,
  renormalized).

- af:

  For `prior = "hwe"`: per-marker ALT allele frequency `p`, length
  `nrow`. `NULL` (default) self-estimates `p` from the reads
  (discouraged).

- error:

  Per-read error rate `epsilon` in `(0, 0.5)` (default `0.01`).

- return:

  `"call"` (default): integer 0/1/2 hard call, `NA` at zero depth.
  `"dosage"`: the same hard call cast to numeric (convenience for
  [`interpolate_genotype()`](https://sawers-rellan-labs.github.io/nilhmm/reference/interpolate_genotype.md))
  – this is the hardcall as a double, **not** `E[G | reads]`; for the
  true posterior-mean dosage take `return = "post"` and compute
  `gp[,,2]*1 + gp[,,3]*2`. `"post"`: the normalized genotype-posterior
  (GP) array `markers x samples x 3` (`NA` rows at zero depth). "GP" is
  the VCF FORMAT field for `P(G | reads)`; the hard `"call"` is its MAP
  (argmax-GP).

## Value

For `"call"`/`"dosage"`: a matrix `markers x samples` (a plain vector if
the inputs were vectors). For `"post"`: the genotype-posterior (GP)
array, `markers x samples x 3` (slices ordered REF/HET/ALT = dosage
0/1/2).

## Details

The prior *is* the argument. `prior` is polymorphic:

- `"flat"`:

  Uniform prior -\> pure argmax-GL (the ML call). Het-blind at depth 1
  (a single ALT read is called hom-ALT), the pessimistic naive caller.

- `"hwe"`:

  Hardy-Weinberg from the per-marker ALT (donor/teosinte) allele
  frequency `p`: `pi = {(1-p)^2, 2p(1-p), p^2}`. This is the
  **het-excess** control – a single ALT read is called HET. `af`
  supplies `p` per marker; if `NULL`, `p` is self-estimated per marker
  from the reads (`rowSums(n_alt) / rowSums(n_ref + n_alt)`) –
  discouraged, as the estimate is circular at low depth.

- numeric `c(f_REF, f_HET, f_ALT)`:

  A **fixed** genome-wide prior, renormalized to sum 1. Covers both the
  **design** prior (a vector derived from the cross – see
  [`design_prior()`](https://sawers-rellan-labs.github.io/nilhmm/reference/design_prior.md))
  and an arbitrary **custom** prior; the code path is identical.

A length-3 numeric vector is genome-wide fixed – it cannot express a
per-marker custom prior. (An `M x 3` matrix could be accepted for that
later; out of scope now.)

## See also

[`design_prior()`](https://sawers-rellan-labs.github.io/nilhmm/reference/design_prior.md),
[`interpolate_genotype()`](https://sawers-rellan-labs.github.io/nilhmm/reference/interpolate_genotype.md),
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)

## Examples

``` r
# A single ALT read decided under four priors:
call_gt(0, 1, prior = "flat")                 # 2 (hom-ALT, het-blind: argmax-GL)
#> [1] 2
call_gt(0, 1, prior = "hwe", af = 0.30)        # 1 (HET, the het-excess control)
#> [1] 1
call_gt(0, 1, prior = design_prior("BC2S3"))   # 2 (design prior resists the het flip)
#> [1] 2
call_gt(0, 1, prior = c(.98, .01, .01))        # custom fixed prior
#> [1] 2
# High depth: the data dominates the prior.
call_gt(0, 10, prior = "hwe", af = 0.05)       # 2
#> [1] 2
```
