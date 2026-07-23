# TODO (later): consolidate the general Viterbi kernels into one

Scope-B cleanup, deferred. Scope A (a unified sweep kernel, `viterbi_sweep_cpp`)
is being built first because it is what the emission x transition parameter grid
(`caller_grid()`, the nnil/Holland calibration, MolBreeding) actually needs, and
its parity gate is contained (2 kernels). This note captures the bigger cleanup so
it can be picked up cleanly.

## The mess: six Viterbi kernels, five of them the same recursion

| kernel | file | transition | batch axis | threaded |
|---|---|---|---|---|
| `viterbi_log_cpp` | src/viterbi.cpp | const matrix | none (1 decode) | no |
| `viterbi_batch_cpp` | src/viterbi_batch.cpp | const (shared) | samples (memoized emission) | no |
| `viterbi_batch_par_cpp` | src/viterbi_batch_par.cpp | const (shared) | samples | yes (RcppParallel) |
| `lb_viterbi_cpp` | src/lbimpute.cpp | distance + drp | none (1 decode) | no |
| `lb_viterbi_sweep_cpp` | src/lbimpute.cpp | distance + drp | transition grid | no |
| `rtiger_viterbi_cpp` | src/rtiger.cpp | rigidity (rigid-duration max-product) | - | no |

The top five run the identical 3-op max-plus recursion + traceback; they differ
only on three axes: transition model (const vs distance+drp), what is batched
(nothing / sample axis + memoized emission / a parameter grid), and threading.

## Scope A (in progress, separate): `viterbi_sweep_cpp`

One kernel with `mode = "const" | "distance"` that sweeps a TRANSITION grid with a
fixed emission. Retires `lb_viterbi_sweep_cpp` (distance mode == it, byte-parity
gated) and `lb_viterbi_cpp` (grid-of-1). The emission axis of a full grid is
driven at the R level (`caller_grid()`: build `log_emit` once per emission combo,
one sweep call per combo). See tests/testthat/test-viterbi-sweep-parity.R.

## Scope B (this TODO): fold the sample-decode kernels in too

Extend `viterbi_sweep_cpp` (or a shared core it calls) to also carry the SAMPLE
batch axis + memoized emission + optional RcppParallel threading, then retire:
`viterbi_log_cpp`, `viterbi_batch_cpp`, `viterbi_batch_par_cpp` (and the already
retired `lb_viterbi_*`). Net end state: ONE general max-plus Viterbi kernel
(axes: samples x emission-combos x transition-values; transition mode const or
distance+drp; threads flag), plus the two genuinely different structures kept
separate: `rtiger_viterbi_cpp` (rigidity) and the FSFHap 5-state EM kernel.

### Gates (do NOT delete before these pass)
1. **Parity matrix, byte-for-byte** against each retired kernel on real fixtures:
   - `viterbi_log_cpp` (single, const), `viterbi_batch_cpp` (sample batch),
     `viterbi_batch_par_cpp` (threaded == serial, the risky one),
     `lb_viterbi_cpp`, `lb_viterbi_sweep_cpp`.
   - Use tests/fixtures/baseline_pre_refactor/ (SHA256SUMS) as the regression anchor.
2. **Keep the transition representation compact** -- const = one shared matrix,
   distance = tpos + param + drp built per-gap in C++. Never materialize T-1 dense
   S x S matrices, or LB-Impute's speed/memory regresses.
3. **Preserve tie-break semantics** (`viterbi_log_cpp` tie_break 0 = first, 1 =
   incumbent; matters for categorical/GT emission-degenerate boundaries).
4. Grep every consumer (engine.R decode path uses `viterbi_batch_cpp`/`_par`;
   rtiger.R, lbimpute.R) + the RcppParallel/libtbb load-order note in
   nilHMM-package.R before removing the threaded kernel.

Risk note: `viterbi_batch_par_cpp` is the engine's hot decode path; collapsing it
mid-calibration is where a silent regression would hurt most. Do Scope B as its
own PR, gated on the full parity matrix, when no calibration run depends on it.
