# Per-site genotype-likelihood (GL) genotype caller with a swappable prior

A closed-form GATK-style diploid genotype call from allelic read counts,
decided independently per (marker, sample) – no linkage model. ALT =
donor (teosinte) = dosage 2. Given per-read error `error` and counts
`(n_ref, n_alt)`, the three genotype log-likelihoods are \$\$L(0) =
(1-\epsilon)^{n\_{ref}}\\\epsilon^{n\_{alt}}\$\$ \$\$L(2) =
\epsilon^{n\_{ref}}\\(1-\epsilon)^{n\_{alt}}\$\$ \$\$L(1) =
0.5^{\\n\_{ref}+n\_{alt}}\$\$ and the call is `argmax_g L(g) * pi(g)`
over the prior `pi`. All arithmetic is in log space and vectorized over
the whole matrix (no per-cell loop). Markers with zero total depth carry
no signal and are returned as `NA`.

## Usage

``` r
call_gl(
  n_ref,
  n_alt,
  prior = c("hwe", "flat", "breeding"),
  af = NULL,
  error = 0.01,
  f = NULL,
  return = c("call", "dosage", "post")
)
```

## Arguments

- n_ref, n_alt:

  Reference / alternate read counts. Integer (or numeric) matrices
  `markers x samples`, or plain vectors (treated as one sample =
  `markers x 1`). Same shape required.

- prior:

  One of `"hwe"` (default), `"flat"`, `"breeding"`.

- af:

  For `prior = "hwe"`: per-marker ALT allele frequency `p`, length
  `nrow`. `NULL` (default) estimates `p` from the reads.

- error:

  Per-read error rate `epsilon` in `(0, 0.5)` (default `0.01`).

- f:

  For `prior = "breeding"`: length-3 genotype frequencies
  `c(f_REF, f_HET, f_ALT)` (renormalized).

- return:

  `"call"` (default): integer 0/1/2 hard call, `NA` at zero depth.
  `"dosage"`: the same hard call as numeric (convenience for
  [`interpolate_genotype()`](https://sawers-rellan-labs.github.io/nilhmm/reference/interpolate_genotype.md)).
  `"post"`: the normalized posterior array `markers x samples x 3` (`NA`
  rows at zero depth).

## Value

For `"call"`/`"dosage"`: a matrix `markers x samples` (a plain vector if
the inputs were vectors). For `"post"`: a `markers x samples x 3` array.

## Details

Priors (`prior`):

- `"hwe"`:

  Hardy-Weinberg from the per-marker ALT (donor/teosinte) allele
  frequency `p`: `pi = {(1-p)^2, 2p(1-p), p^2}`. This is the
  **het-excess** control – a single ALT read is called HET. `af`
  supplies `p` per marker; if `NULL`, `p` is estimated per marker from
  the reads (`rowSums(n_alt) / rowSums(n_ref + n_alt)`).

- `"flat"`:

  Uniform prior -\> pure argmax-GL. Het-blind at depth 1 (a single ALT
  read is called hom-ALT), the pessimistic naive caller.

- `"breeding"`:

  Fixed design frequencies `f = c(f_REF, f_HET, f_ALT)` (e.g. BC1S4 ~
  `c(0.84, 0.06, 0.10)`); renormalized to sum 1.

## See also

[`interpolate_genotype()`](https://sawers-rellan-labs.github.io/nilhmm/reference/interpolate_genotype.md),
[`call_ancestry()`](https://sawers-rellan-labs.github.io/nilhmm/reference/call_ancestry.md)

## Examples

``` r
# A single ALT read: het-blind under flat, HET under HWE (het-excess control).
call_gl(0, 1, prior = "flat")               # 2 (hom-ALT)
#> [1] 2
call_gl(0, 1, prior = "hwe", af = 0.3)       # 1 (HET)
#> [1] 1
# High depth: the data dominates the prior.
call_gl(0, 10, prior = "hwe", af = 0.05)     # 2
#> [1] 2
```
