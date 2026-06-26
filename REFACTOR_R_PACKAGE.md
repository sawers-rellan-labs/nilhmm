# nilHMM → R-package refactor plan (2026-06-26)

Design spec for turning `nilhmm` (currently Python: `hmmlearn` + the custom counts path) into a
single **R package, `nilHMM`**, housing all functions required for ancestry calls. Companion to
`Implementation.md` (the Python design this refactors). Cross-repo paths (`data/…`, `agent/…`,
`results/…`, `docs/…`) are relative to the **`zealtiger`** analysis project (parent of
`agent/nilhmm`), not this repo — same convention as `Implementation.md`.

> **Sequencing decision (do this first):** run nilHMM on BRB-seq with the **current Python**
> implementation *before* refactoring — see §9. This note is written down now so the design is
> fixed before the BRB work; the refactor itself starts after BRB.

---

## 1. Goal & motivation

Make `nilHMM` an R package with all ancestry-calling functions, called from R **pipeline
scripts**. Drivers:

- **RTIGER is unmaintained (4+ years).** Its Julia + JuliaCall stack is a maintenance liability
  (Julia/version drift, bridge breakage; the zealtiger handover is full of `JULIA_HOME`/`juliaup`
  friction), not a live dependency. Reimplementation is warranted.
- **R is the lab language** (Bioconductor/GenomicRanges already in use). Consolidating the
  current polyglot toolchain — RTIGER (R+Julia), nilHMM (Python), Skim-BIN (R) — into **one
  R + Rcpp package** removes two language runtimes and unifies the interface.
- **The callers are already one model.** nilHMM's counts path is ≈ RTIGER's core (BetaBinomial-HMM
  on the same reads); Skim-BIN is the same HMM with a Gaussian emission. They differ in
  *emission* and *duration*, not in kind → one engine expresses all three.

## 2. Naming (resolve the collision)

- **Package = `nilHMM`** (keeps continuity with the nNIL-of-Holland origin; packages outgrow
  their origin name).
- **Callers = explicit methods**, never "nilHMM": `caller = "nnil"` (Holland's nNIL),
  `"rtiger"` (rigidity mode), `"skimbin"`. So "nilHMM" denotes the *package*; the callers are
  named methods inside it. This kills the package-vs-caller ambiguity.

## 3. Architecture — three layers

1. **Engine** (core) — the duration-aware HMM (§4).
2. **Presets** (on top of the engine) — emission-by-depth regimes (§5) and breeding-design
   references (§6). Convenience/defaults; not the engine.
3. **Pipeline scripts** (consumers) — `library(nilHMM)`, call its functions on actual data,
   write common-schema outputs. **Not part of the package.** They own the paths/sample lists.

**Hard rule:** the package is **data-agnostic — no hardcoded paths, sample lists, or mount
locations.** Functions take `(data, params) → calls`; scripts own "which files, which samples,
where outputs go." This is what lets one package serve skim / BRB / MolBreeding / future data.

## 4. The engine

- **States:** 3 hidden — REF / HET / ALT.
- **Swappable emission** (the key abstraction; the axis that distinguishes callers):
  - `count` — **BetaBinomial** on `(n_ref, n_alt)`, state means `θ = {err, 0.5, 1−err}`,
    depth-0 = flat. (nilHMM counts path, RTIGER.)
  - `gt` — **categorical** over `{0,1,2,missing}` + genotype-error matrix. (Holland's native GT
    path; MolBreeding regime.)
  - `dosage` — **Gaussian/Beta** centred at `0/1/2`, variance = imputation uncertainty. (Skim-BIN
    style; imputed-dosage sources.)
- **Duration layer** (transition):
  - `geometric` — self-transition rate `r` (memoryless). (nilHMM.)
  - `rigidity` / phase-type — hard minimum run length of `r` markers via state expansion (a chain
    of `r` sub-states). **This is RTIGER's rigidity** — and the same machinery as an explicit
    Erlang duration.
  - `hsmm` (optional, later) — explicit-duration sojourn distribution. **Not** for baking in the
    fitted Gamma (see §7) — reserved if a genuine duration model is ever needed.
- **Fit/decode:** Baum-Welch EM for emission params where applicable; Viterbi for the path; **Rcpp
  for the hot loops** (Viterbi, emission, EM) — recovers Julia-class speed within an R package
  using the standard, CRAN/Bioconductor-friendly toolchain.
- **Agnostic core:** no generation prior in the transition (RTIGER-style). Generation enters via
  presets (priors) and calibration, never hardcoded in the engine.

## 5. Emission-by-depth presets (the regime axis)

Emission choice is **depth-driven** (see zealtiger `docs/emission_by_depth_regime`):

| Source | Depth | Regime | Emission |
|---|---|---|---|
| Skim (SNP50K) | ~0.4× | degenerate (HET invisible at depth 1) | `count` (best available; still degenerate) |
| BRB-seq | ~2.8× | intermediate — counts pay off, `conc` matters | `count` |
| MolBreeding (GBTS) | ≥~20× | **saturated** — BetaBinomial → delta = hard call | `gt` (equivalent + cheaper) |

Selector rule: `depth-saturated (≥~20×) → gt; intermediate (~1–20×) → count; imputed → dosage`.
Cost basis: BetaBinomial cost ∝ #distinct `(n,k)` ∝ coverage → cap ~20–30× → above the cap it's a
hard call (memory `rtiger-betabinomial-cost`).

## 6. Breeding-design presets (the generation axis) + bundled data

**Orthogonal** to the emission axis: a dataset carries *both* a depth-regime (→ emission) and a
breeding design (→ generation-derived quantities). E.g. **skim = `count` × BC2S2**, **BRB =
`count` × BC2S3+** (memory `bzea-skim-bulked-bc2s2`).

Generation determines: the single-locus genotype-freq priors (`f_1, f_2`) and the **expected
fragment-size law** (the fitted Gamma) used as **Null / calibration target**.

**Bundled reference data:**
- `data/maize_map_v5` — consensus map (marker → chr, cM, bp), **version-namespaced**, overridable;
  powers position files (`c_m`), cM↔Mb, optional map-aware transitions.
- `data/breeding_designs` — per design (`BC1S1`, `BC2S2`, `BC2S3`, …): `g`, `f_1`, `f_2`, Gamma
  `(k, λ)`, mean cM, provenance.

**Accessors:** `design_priors("BC2S2")`, `expected_fragment_dist("BC2S3")` (params + density/CDF
closure for Nulls), `fit_design_gamma(sim_segments)` (extensible), `cm_to_mb(seg, map)`.

**Regeneration (`data-raw/`, not shipped):** the bundled objects are *derived*; `data-raw/` holds
the **scripts** that regenerate them — `make_breeding_designs.R` runs simcross from a **pinned
seed + recipe** (pedigree, map cM lengths, Stahl `m=10, p=0`, `n=1500`) → fits the Gamma →
`use_data()`. No raw simulation CSV is vendored; the seed + recipe regenerate it. **Pin the seed
*and* the simcross/R versions + `RNGkind()`** — a seed reproduces the exact shipped `(k,λ)` only
under a fixed environment. The stored `data/` value is authoritative at runtime; the seeded script
is for verification/extension. Same pattern for the map.

## 7. Key boundaries / decisions

- **Reimplement RTIGER as the engine's `rigidity` mode — do NOT wrap it.** It's unmaintained;
  reimplementing the *subset we use* (BetaBinomial emission + rigidity + Viterbi + KS-calibrated
  `r`) is tractable and kills the Julia dependency. Skip `optimize_R` (we don't use it; it
  over-rigidifies — memory `optimize-R-overrigidifies-nil-sim`). Validate against the existing
  Skim-RTI calls (`data/rtiger_50K/calls_taxa_r5.csv`). Cite the RTIGER paper as the algorithm.
- **The fitted Gammas are calibration/Null only — NEVER engine priors.** Baking the design's Gamma
  into the transition/duration would make the caller reproduce it by construction and **mask the
  open "callers run longer than the model" signal** (recombination suppression). Presets feed
  *plotting* (the grey Null) and *calibration* (KS-vs-sim target for tuning `r`); the engine stays
  agnostic. (The Task-2 circularity trap — zealtiger `docs/16`.)
- **`r` is not MLE'd** (memory `rigidity-not-mle`): it's a resolution hyperparameter set by
  KS-vs-simulated-truth, not by maximizing a likelihood. The generative parameter is the Gamma's
  `λ`; `r` is the caller knob tuned to it.
- **Map is version-namespaced; warn on assembly mismatch.** cM-space Gammas are assembly-robust;
  anything bp/Mb is tied to B73 v5.

## 8. Package structure

```
nilHMM/
  R/    engine (fit/decode), emissions {count|gt|dosage}, duration {geometric|rigidity|hsmm},
        callers {nnil|rtiger|skimbin}, presets_regime (emission by depth),
        presets_design (priors + expected Gamma), map utils, calibrate (KS-vs-sim), plot, io
  src/  Rcpp: Viterbi + EM hot loops
  data/        maize_map_v5, breeding_designs           (lazy-loaded, overridable)
  data-raw/    make_maize_map.R, make_breeding_designs.R (seed + recipe; .Rbuildignore'd)
  tests/       regression vs frozen baselines (§9); preset-value checks
  vignettes/   the zealtiger regime notebooks → package docs
```

Top-level API: `call_ancestry(data, caller = c("nnil","rtiger","skimbin"), source/design presets, …)`
wrapping the engine `fit()/decode()`; emission and duration as pluggable interfaces.

## 9. Sequencing & validation (decided)

**Run BRB on the current Python nilHMM FIRST, then refactor under a regression net.** Rationale:
answers the science now (independent of the refactor); produces a BRB regression fixture in a
*new regime* (depth ~2.8×, BC2S3 priors, `conc` active); keeps new-science and new-code separate
(never debug biology + port at once); and surfaces engine requirements before the API is frozen.

1. **Run BRB on current Python** (Task 1): TSV→counts adapter for `read_vcf_counts` (reuse the
   wideseq-thinned counts `results/sim_calibration/brbseq_ks_wideseq/counts`); **BC2S3 priors**
   (`f_1=0.0312, f_2=0.1094`, NOT the BC2S2 skim values); calibrate `r`/`err`/**`conc`** per
   teosinte species vs the BC2S3 sim (`results/nil_segments.csv`); B73 control = clean BrB rep
   `PN3_SID213`. → answers the HET-tail hypothesis (does nilHMM suppress the RNA-noise HET tail
   RTIGER shows, while contamination/impostors persist?).
2. **Freeze fixtures:** the skim per-taxon and BRB common-schema call CSVs **plus** the calibration
   params (`r/err/conc` per taxon) — the "previous-to-refactor" baselines.
3. **Refactor** to the R package (§§3–8).
4. **Validate:** the R engine reproduces both fixtures within tolerance — **distributional**
   equivalence (call concordance, segment-size KS), not bit-identical (EM init / Viterbi
   tie-breaking differ across languages; same bar as the counts-optimization element-wise check
   in `Implementation.md`).

## 10. Open items — BRB run settled some of these (2026-06-26, see `BRB_run_findings.md`)

- **Emission means must be FITTABLE, not fixed** — **(SETTLED by BRB; new requirement).** The BRB
  run showed nilHMM collapses ALT→HET (ALT mean 0.0003) because BRB donor signal reads at
  alt-fraction ~0.5 (median 0.500 among alt-bearing markers — reference-mapping bias and/or het),
  which the fixed `θ_ALT = 1−err = 0.99` cannot represent. RTIGER recovers ALT (0.081) via EM-fit
  emissions. → the engine needs **EM-fit or reference-bias-corrected emission means**, not the
  fixed `[err, 0.5, 1−err]`, for RNA / ref-biased data.
- **`conc` is a near-no-op — (SETTLED, was the wrong suspect).** Probed `conc ∈ {3,8,20,50}` on BRB:
  ALT stays ≈0, donor_frac ~0.05 throughout. Overdispersion is not the knob; the emission *mean*
  mismatch dominates. The earlier "calibrate conc at depth>1" guess is superseded.
- **Per-platform calibration objective** — **(SETTLED: KS-on-block-size is insufficient for BRB).**
  The `r` sweep rails to the rigid grid edge because donor *rate* is r-invariant (set by
  prior+emission); KS-on-donor-block-size assumes the caller detects the sim footprint, which
  nilHMM doesn't on BRB. Need a rate/footprint-aware objective, or fix the emission first.
- **BC2S3 design-preset values** — confirm `g, f_1, f_2` (used: 0.0312/0.1094), and the fitted
  BC2S3 Gamma for `breeding_designs`.
- **Map-aware (position-dependent) transitions** — still open. The real gain is map-awareness
  (cM→Mb / LogNormal), not the small cM-Gamma shape (`k≈1.15`).
- **GT path resurrection** — MolBreeding (saturated) is the dataset that exercises the `gt` emission
  skim couldn't; plan a consistency check (gt vs count on MolBreeding's real depths).

## 11. References

- zealtiger docs: `docs/emission_by_depth_regime.html`, `docs/generative_vs_discriminative_callers.html`,
  `docs/elai_admixture_mapping_plan.html`, `docs/16-parameter-search-and-rigidity-mle.md`,
  `fragment_size_cm_validation`.
- This repo: `Implementation.md` (Python design + the counts extension + calibration).
- Memories: `rtiger-betabinomial-cost`, `nilhmm-counts-path`, `bzea-skim-bulked-bc2s2`,
  `b73-control-impostors`, `optimize-R-overrigidifies-nil-sim`, `rigidity-not-mle`,
  `molbreeding-45k-source5`.
