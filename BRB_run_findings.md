# nilHMM on BRB-seq — run & HET-tail findings (2026-06-26)

Ran the nilHMM count caller on the **wideseq-thinned BRB counts**
(`results/sim_calibration/brbseq_ks_wideseq/counts`, 379 samples × 21,665 markers) with
**BC2S3 priors** (`f_1=0.0312, f_2=0.1094`), per teosinte species, calibrating `r` by KS of the
donor-block size distribution vs the **BC2S3 sim** (`results/nil_segments.csv`). Control = clean
BrB rep `PN3_SID213`. Scripts: `agent/brb_nilhmm_counts.py` (TSV→counts adapter),
`agent/ks_sweep_nilhmm_brbseq.py`, `agent/brb_nilhmm_het_compare.py`. Outputs:
`results/sim_calibration/nilhmm_brbseq/`.

## Q (the hypothesis): does nilHMM drop the high-HET tail BrB-RTI shows? — **YES**

Per-sample HET genome fraction (length-weighted), 341 NILs:

| method | median | p90 | mean | frac HET=0 | **frac HET>0.20** |
|---|---|---|---|---|---|
| BrB-RTI (r=7) | 0.00 | **0.44** | 0.093 | 61.9% | **12.9%** |
| nilHMM (BC2S3) | 0.033 | **0.091** | 0.042 | 2.3% | **0.9%** |

BrB-RTI is **bimodal** — 62% HET=0 plus a heavy artifact tail (12.9% of samples >0.20, p90=0.44 =
contamination/mislabel/RNA noise, `b73-control-impostors`). **nilHMM suppresses that tail**
(0.9% >0.20, p90=0.09) and is unimodal/moderate. Hypothesis confirmed.

## But two caveats came with it

1. **nilHMM under-calls donor broadly.** Donor (HET+ALT) genome fraction: nilHMM **0.042** vs
   BrB-RTI **0.174** vs BC2S3 expected ~0.14. So it doesn't cleanly separate noise from signal —
   it's *globally conservative* (BC2S3 low-HET prior + explicit error model), trading RTIGER's
   noisy over-calling for systematic under-calling (~⅓ of expected donor). The tail is gone, but
   so is sensitivity.
2. **nilHMM zeroes ALT** (mean 0.0003; all donor → HET). This is **not** the depth-1 degeneracy
   (BRB ~2.8×). Diagnosis: among BRB markers with any alt read, the **alt-fraction is median 0.500
   / mean 0.563** (only 41% > 0.8) — donor-hom signal sits near 0.5, not 1 (reference-mapping bias
   and/or genuine het). nilHMM's **fixed** emission mean `θ_ALT = 1−err = 0.99` cannot represent a
   donor-hom that reads at ~0.5, so those markers land on the HET state. RTIGER recovers ALT
   (0.081) because its emissions are **EM-fit** and adapt to the observed alt-fraction.

## `conc` is NOT the knob (the refactor plan's suspect was wrong)

Probed Zx at `r=1e-8`, `conc ∈ {3,8,20,50}`: ALT stays ≈0 (0.0001→0.0010), donor_frac ~0.05–0.06
throughout. Overdispersion doesn't recover ALT — because the problem is the emission **mean**, not
its spread. → For RNA / reference-biased data the engine needs **EM-fit or bias-corrected emission
means**, not the fixed `[err, 0.5, 1−err]`.

## Calibration note (r-edge)

The KS-vs-sim `r` sweep hits the **rigid grid edge** for every species (argmin r=1e-8; D 0.12–0.25;
median block 4–9 Mb < sim 10.74). Donor *rate* is nearly **r-invariant** (Zx donor_frac
0.051→0.049 across r=1e-8→1e-4) — set by the prior+emission, not by `r`. So `r` only fragments;
the limiting factor is the donor-calling rate. **KS-on-block-size is a poor calibration target when
the caller under-calls donor.** (This is why the HET comparison is robust to the `r` choice — the
donor rate doesn't move with `r`.)

## Implications for the refactor (update open items)

- **Emission means must be fittable**, not fixed — the BRB ALT collapse is a fixed-`θ` artifact
  under reference bias. EM-fit emissions (RTIGER-style) or an explicit reference-bias term. This is
  the concrete requirement BRB surfaced; supersedes the earlier "calibrate conc" guess.
- **`conc` is a near-no-op here** — depth>1 did not make it matter (the emission *mean* mismatch
  dominates).
- **Per-platform calibration objective:** KS-on-donor-block-size assumes the caller detects the
  sim's donor footprint; on BRB nilHMM doesn't, so the objective rails to the edge. A
  rate-aware or footprint-aware objective (or fixing the emission first) is needed for BRB.

## Bottom line

nilHMM-on-BRB **does drop the pathological high-HET tail** — your hunch was right — but by being
**conservative across the board** (under-calls donor, collapses ALT→HET) rather than by clean
noise/signal separation. The ALT collapse is a fixed-emission-mean artifact under BRB's
reference-biased ~0.5 alt-fraction, fixable only by EM-fit/bias-aware emissions (a refactor
requirement), not by `conc`.
