# Pre-refactor regression baseline (frozen 2026-06-30)

Frozen outputs of the **trusted Python `nilhmm`** (and the RTIGER rigidity reference), captured
*before* the R+Rcpp reimplementation begins. These are the "previous-to-refactor" baselines of
`REFACTOR_R_PACKAGE.md` §9.2. The R engine must reproduce them within tolerance (§9.4) —
**distributional** equivalence (call concordance, segment-size KS), **not** bit-identical
(EM init / Viterbi tie-breaking differ across languages).

**Immutable.** Do not regenerate from the new R code — that would be circular. Verify integrity
with `shasum -c SHA256SUMS` (run from this directory).

Inputs are large and intentionally **not vendored**; they stay in the zealtiger analysis project
and are referenced by path below. Paths are relative to the zealtiger project root (parent of
`agent/nilhmm`).

---

## What is frozen

### `skim_nilhmm/` — SNP50K, count caller (BC2S2)
- `calls_common_schema.csv`, `calls_checks_common_schema.csv`
- Source: `results/sim_calibration/nilhmm_calls/`
- Regime: skim ~0.4×, `count` emission, **BC2S2** priors (bulked S3; memory `bzea-skim-bulked-bc2s2`).
- Input counts (not vendored): `data/rtiger_50K/counts` (~2.9 GB; 50K cohort AD counts, not imputed).

### `brb_nilhmm/` — BRB-seq, count caller (BC2S3)
- `calls_common_schema.csv`, `calls_checks_common_schema.csv` — the calls.
- `calibration_params.csv` — per-taxon `best_r / D / median_mb / donor_frac / b73_control_dosage /
  at_grid_edge` (was `summary.csv`). **All five species rail to the rigid grid edge (`r=1e-8`,
  `at_grid_edge=True`)** — expected: donor *rate* is r-invariant on BRB (findings §"Calibration
  note"), so KS-on-block-size is a poor objective here.
- `Z{d,h,l,v,x}_ks_sweep.csv` — per-species KS-vs-sim `r` sweeps.
- `per_sample_fracs_nilhmm.csv`, `het_tail_comparison.csv` — HET-tail diagnostics (the science:
  nilHMM drops the BrB-RTI high-HET artifact tail; see `BRB_run_findings.md`).
- Source: `results/sim_calibration/nilhmm_brbseq/`
- Regime: BRB ~2.8×, `count` emission, **BC2S3** priors (`f_1=0.0312, f_2=0.1094`).
- Input counts (not vendored): `results/sim_calibration/brbseq_ks_wideseq/counts`
  (~178 MB; 379 samples × 21,665 wideseq-thinned markers).
- B73 control: clean BrB rep **`PN3_SID213`** (memory `b73-control-impostors`).
- Producing scripts: `agent/brb_nilhmm_counts.py` (TSV→counts adapter),
  `agent/ks_sweep_nilhmm_brbseq.py`, `agent/brb_nilhmm_het_compare.py`.

### `sim_truth/` — KS calibration targets / Null
- `bc2s3_nil_segments.csv` — BC2S3 simulated NIL segments (`results/nil_segments.csv`); the
  grey-Null / KS target for fragment-size calibration.
- `brbseq_benchmark_truth_segments.csv` + `brbseq_benchmark_params.json` — the BRB benchmark sim
  truth and its recipe (`results/sim_calibration/brbseq_benchmark/`). Recipe: seed `20260620`,
  `lambda_mean=2.81`, Stahl `m=10, p=0`, simcross `0.8`, realized mean coverage 2.81×.

### `rtiger_rigidity_ref/` — RTIGER cross-check for `rigidity` mode (§7)
- `calls_taxa_r5.csv` — the existing Skim-RTI calls (`data/rtiger_50K/calls_taxa_r5.csv`). The
  engine's reimplemented `rigidity` mode validates against these, NOT against nilHMM.

---

## Known caveat carried into the refactor (do NOT treat as a regression)

The Python collapses **ALT→HET on BRB** (ALT mean ~0.0003) because its **fixed** emission mean
`θ_ALT = 1−err = 0.99` cannot represent BRB donor-hom markers that read at alt-fraction ~0.5
(reference-mapping bias). The R engine deliberately fixes this with **fittable / bias-corrected
emission means** (§10) — so the R ALT calls are *expected to differ* from this baseline. Validate
BRB on the HET-tail behavior and the count-caller distribution, and use the skim baseline (where
fixed means are adequate) as the strict concordance gate. See `../../../BRB_run_findings.md`.