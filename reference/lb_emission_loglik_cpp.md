# LB-Impute coverage-aware emission (log-space, REF/HET/ALT)

Per-marker emission from allelic read depths, the model of
`ImputeOffspring.getprobabilities2`. Raw state likelihoods are
\\E\_{REF}\propto(1-err)^{n\_{ref}}err^{n\_{alt}}\\,
\\E\_{ALT}\propto(1-err)^{n\_{alt}}err^{n\_{ref}}\\, \\E\_{HET}\propto
0.5^{n\_{ref}+n\_{alt}}\\; each is divided by the per-marker maximum,
scaled by \\1-2\\err_g\\ and offset by \\err_g\\, bounding every
emission to \\\[err_g,\\1-err_g\]\\ so a single artifactual marker
cannot dominate the path. A zero-coverage marker emits flat (all states
equal).

## Usage

``` r
lb_emission_loglik_cpp(nref, nalt, err, errg)
```

## Arguments

- nref, nalt:

  Integer per-marker reference / alternate read counts.

- err:

  Per-read sequencing-error probability (LB-Impute `readerr`).

- errg:

  Coverage-independent genotyping-error probability (LB-Impute
  `genotypeerr`); sets the emission floor/ceiling.

## Value

A T x 3 matrix of log emission probabilities (columns REF/HET/ALT).
