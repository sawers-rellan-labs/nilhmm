# RTIGER windowed emission product (productpsi)

PSI(i,t) = sum of psi over the window (t-r+1 .. t) (sliding window of
r), as a running cumulative sum. PSI is s x (T+r): the tail columns
T+1..T+r-1 carry the trailing partial sums and column T+r is left 0
(exactly as the fork).

## Usage

``` r
rtiger_productpsi_cpp(logpsi, r)
```

## Arguments

- logpsi:

  s x T matrix from rtiger_getlogpsi_cpp.

- r:

  Rigidity.

## Value

s x (T+r) matrix.
