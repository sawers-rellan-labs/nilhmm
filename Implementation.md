# nilhmm on BZea SNP50K skim — implementation notes (2026-06-25)

This document records the work done to bring `nilhmm` onto the BZea SNP50K ("skim")
dataset as a genotyping caller, the data problems found, the count-based extension added,
the calibration sweeps and their results, and a set of suggestions and open questions.

It spans two repositories:
- **this repo** (`nilhmm`) — the caller code that was extended.
- the **`zealtiger`** analysis project (parent of `agent/nilhmm`) — data staging, simulation
  references, sweep scripts, and the cross-caller comparison docs. Paths like
  `data/nilhmm/…`, `agent/…`, `results/…` below are relative to that project, not this repo.

---

## 1. Intent (restated)

The goal is to add **`nilHMM`** as a **4th ancestry/introgression caller** alongside the three
already used in the BZea NIL analysis:

| caller | reads | model |
|---|---|---|
| Skim-RTI | SNP50K ref/alt **counts** | RTIGER BetaBinomial HMM (rigidity `r`) |
| Skim-BIN | SNP50K per-1 Mb **alt fraction** | GMM-HMM |
| BrB-RTI | BRB-seq 3′ RNA counts | RTIGER BetaBinomial HMM |
| **nilHMM** | called **GT** (per the upstream design) | genotype-HMM + explicit error model |

Key constraints and decisions that frame the work, as I understood them:

1. **Start small.** Don't run the whole 1,439-sample cohort first; calibrate on **one small,
   single-divergence group** (one teosinte taxon) before scaling out. This mirrors the existing
   RTIGER scheme, which fits **per teosinte species (taxon)**, not as a mixed cohort.
2. **Skim is BC2S2.** The SNP50K skim genotypes are bulks of ~6 S3 plants per S2 ear, so the
   bulk recovers the **S2** parent → calibrate skim against the **BC2S2** generation, not BC2S3.
   (BRB-seq is a different, later generation — out of scope here.)
3. **Calibrate against the simulation.** Use the established KS-on-segment-size framework: tune
   the caller so the **donor-block size distribution** matches the BC2S2 simcross simulation.
4. **B73 controls must stay clean.** Any parameter choice has to keep the B73 negative controls
   at ≈0 donor.
5. **Per-taxon frequencies / grouping**, like RTIGER (species differ in divergence from B73).
6. **Local commits, hold the push.** Work is committed locally in this clone; pushing to the
   shared lab repo is a separate, explicit decision.

---

## 2. Environment & sanity (Phase 0)

- Conda env `nilhmm` = clone of the existing `base` env (numpy/pandas/scipy/sklearn already
  present) **+ `hmmlearn` 0.3.3 + `pip install -e .`** of this repo. We avoided rebuilding from
  `envs/nilhmm.yml` because `base` already had the heavy stack; only `hmmlearn` was missing.
- Chromosome naming: the SNP50K cohort VCF uses `chr1..chr10` (verified via `tabix -l`).
  `io.read_vcf` strips the `chr` prefix to an int natively, so **no renaming preprocess** is
  needed; non-numeric contigs are safely skipped.
- **Bug fixed:** `core.introgression_hmm` (the GT caller) hardcoded `for chrom in [1..10]` and
  unconditionally indexed `marker_dict[chrom]`, raising `KeyError` whenever the input lacked all
  10 chromosomes (its own test fixtures, or any subset/filtered VCF). Changed both loops to
  iterate `sorted(marker_dict.keys())`. Behavior is identical for full-genome input.
  Tests went from 3 failed / 2 passed → **5 passed**.

---

## 3. The data problem (Phase 1): GT is unusable at skim coverage

`nilhmm` as shipped reads **only the `GT` field**. At ~0.4× skim, per-site GT is mostly missing
and, where present, single-read (so heterozygotes look homozygous). Measured on the smallest
taxon, **Zh (n=61, 49,002 markers)**:

| file | content | missing GT | het (of called) | **donor-hom** |
|---|---|---|---|---|
| `bzea_50K_cohort.vcf.gz` | GT only | **73%** | 3.9% | **0.00%** |
| `bzea_50K_cohort_ref.vcf.gz` | GT only (+217 ref-panel samples) | 73% | 3.9% | 0.00% |

The `README.md` in the joint-genotyping output confirms these are raw `bcftools mpileup+call`
genotypes — **never imputed** — and no imputed/phased VCF exists anywhere on the mount. The
consequence is fundamental: **the donor-homozygous state — the thing nilHMM exists to detect —
is essentially invisible in the GT (0.00%)**, because at 0.4× you almost never see two donor
reads at one site.

One file does carry counts: **`cohort.vcf.gz` (1,439 samples) has `GT` + `AD` + `PL`.** It is the
only viable input. Its allele-depth coverage is the same ~27% of sites, mean depth ≈ **1.25
reads** where covered — no per-site genotype resolution, but enough for an HMM that **pools**
single-read observations along a segment (this is exactly why RTIGER works on these data).

**Decision:** rather than imputing GT or shelving nilHMM, extend it to read **AD counts**
directly. (Accepted trade-off: a count-based nilHMM converges toward what RTIGER already does;
its remaining distinctiveness is the explicit per-state error model + per-taxon priors.)

---

## 4. The count-based extension

Added to this repo (GT path untouched; tests now 7 passed):

- **`io.read_vcf_counts(vcf, chromosomes=None)`** — parses `FORMAT/AD` into `ref` and `alt`
  integer matrices (samples × markers), biallelic SNPs only; missing/absent AD → `(0,0)`.
  Returns the same `marker_dict` / `marker_info` / sample layout as `read_vcf`.
- **`core.introgression_hmm_counts(ref, alt, marker_dict, err, conc, r, f_1, f_2)`** — a
  beta-binomial emission HMM:
  - states `0/1/2` = B73-hom / het / donor-hom, with expected alt fractions
    `theta = [err, 0.5, 1-err]`.
  - emission for marker *i*, state *s* = `BetaBinomial(alt_i | n_i, theta_s·conc, (1-theta_s)·conc)`.
  - markers with depth 0 contribute a flat (uninformative) emission.
  - transition/init reuse the **same `r`/`f_1`/`f_2` parameterization as the GT caller**
    (`_build_transition`), with a log-space Viterbi (`_log_viterbi`) per chromosome.
- **`core.call_introgressions_counts(vcf, output_prefix, **params)`** — CLI-equivalent wrapper
  that writes the standard nilhmm output files.

New tests: `test_read_vcf_counts` (AD parsing) and `test_introgression_hmm_counts_detects_donor_block`
(a synthetic depth-1 donor block is recovered as state 2).

### First validation (uncalibrated defaults, r=0.01)

| group | n | REF / HET / ALT | donor dosage (mean / median) |
|---|---|---|---|
| Zh teosinte NILs | 61 | 91.9 / 8.0 / 0.1 | 0.041 / 0.039 |
| B73 controls | 12 | 99.36 / 0.64 / 0.00 | 0.0032 / 0.0000 |

The count path **recovers introgression signal where GT gave nothing** (0% → real donor calls),
and the B73 controls behave as a clean negative control (`PN10_SID893` = 0.0000).

---

## 5. Subsetting (per taxon), mirroring RTIGER

`agent/make_nilhmm_taxon_subset.sh <GROUP>` builds a per-group counts subset from
`cohort.vcf.gz` (the AD source), biallelic SNPs, keeping GT+AD →
`data/nilhmm/<GROUP>_{samples.txt,counts.vcf.gz}`.

- **Teosinte taxa** (`Zx/Zv/Zd/Zl/Zh`): sample → `donor_id` via
  `data/sample_metadata_master.csv` (`project==bzea`, not a check, `donor_id!=NA`),
  `taxon = substr(donor_id, 1, 2)`.
- **Check groups** (`B73`, `Purple`): checks have `donor_id=NA`, so selected by `sample_group`.

Group sizes within the 1,439-sample cohort VCF:

| group | species | n |
|---|---|---|
| Zx | ssp. *mexicana* | 580 |
| Zv | ssp. *parviglumis* | 386 |
| Zd | *Z. diploperennis* | 255 |
| Zl | *Z. luxurians* | 121 |
| Zh | ssp. *huehuetenangensis* | 61 |
| B73 (check) | — | 12 |
| Purple (check) | — | 17 |

Zh (smallest taxon) was used as the calibration starter.

---

## 6. Calibration (Phase 2) on Zh

**Objective:** KS distance `D` between the called **donor-block sizes (Mb)** and the **BC2S2**
simulation (`agent/bc2s2_segments.csv`, 15,567 segments, median 12.41 Mb). A "donor block" is a
contiguous run of `state>0` (HET ∪ ALT, matching the simulation's donor-union segments), with
block boundaries placed at midpoints between flanking opposite-state markers. Blocks pooled
across all 61 Zh samples. Script: `agent/ks_sweep_nilhmm_zh.py`.

### 6.1 Recombination `r` sweep (err=0.01, conc=20, BC2S2 freqs f_1=0.0625, f_2=0.0938)

| r | D | median block (Mb) | donor frac |
|---|---|---|---|
| 1e-5 | 0.067 | 15.26 | 0.087 |
| **3e-5** | **0.043** | **12.89** | 0.086 |
| 7e-5 | 0.056 | 11.84 | 0.086 |
| 1e-4 | 0.067 | 11.31 | 0.086 |
| 2e-4 | 0.076 | 10.16 | 0.086 |
| 3e-4 | 0.085 | 9.24 | 0.086 |
| 5e-4 | 0.109 | 8.13 | 0.086 |
| 1e-3 | 0.142 | 7.18 | 0.085 |
| 3e-3 | 0.201 | 5.54 | 0.084 |
| 1e-2 | 0.306 | 3.36 | 0.081 |
| 3e-2 | 0.472 | 1.36 | 0.075 |
| 1e-1 | 0.678 | 0.39 | 0.061 |
| 3e-1 | 0.906 | 0.06 | 0.026 |

**Argmin at `r = 3e-5`, D = 0.043** — a clean U-shaped minimum (not pinned to a grid edge),
median block 12.89 Mb vs the sim's 12.41. For reference, the BrB-RTI calibration's best was
D ≈ 0.103, so this is a notably tighter fit. The total donor fraction (~8.6%) is essentially
flat across `r` — i.e. `r` sets segment **length**, not how much donor is called.

### 6.2 Emission `err` sweep (at r=3e-5)

| err | D (seg-size) | HET% | ALT% | dosage | HET:ALT | B73 max dosage |
|---|---|---|---|---|---|---|
| **0.01** | **0.043** | 8.60 | 0.014 | 0.043 | 620 | 0.032 |
| 0.02 | 0.081 | 7.95 | 0.015 | 0.040 | 520 | 0.031 |
| 0.05 | 0.121 | 6.93 | 0.027 | 0.035 | 252 | 0.017 |
| 0.10 | 0.193 | 5.86 | 0.064 | 0.030 | 92 | 0.008 |
| 0.20 | 0.341 | 4.22 | 0.245 | 0.024 | 17 | 0.000 |
| 0.30 | 0.353 | 2.02 | 0.765 | 0.018 | 2.6 | 0.000 |

**Tension:** raising `err` recovers ALT toward skim-RTI's ~2:1 HET:ALT ratio, **but it destroys
the segment-size fit** (D 0.043 → 0.35) and shrinks total donor. Since KS-on-segment-size is the
blessed calibration objective, `err` is kept low. Note `conc` is irrelevant at depth 1
(`BetaBinomial(a | n=1) = theta` regardless of concentration); it would only matter at the rare
deeper sites.

### 6.3 Calibrated Zh parameters

```
r = 3e-5,  err = 0.01,  conc = 20,  f_1 = 0.0625,  f_2 = 0.0938   (skim = BC2S2)
→ KS D = 0.043,  median donor block 12.89 Mb (sim 12.41),  donor frac 0.086
→ B73 control guard: dosage mean 0.0026, median 0.0000, max 0.0317   (clean)
```

Saved to `results/sim_calibration/nilhmm_zh/calibrated_params.csv`.

**Honest limitation to carry forward:** at ~1× effective depth the **het/hom split is not
identifiable per site** — nilHMM-counts is strongly het-biased (HET:ALT ≈ 600). The reliable
outputs are the **donor footprint** and its **block-size distribution**; the het-vs-hom call
within a donor block is not. This is the count-path version of the original "GT-at-0.4× is the
make-or-break" caveat: the counts recover the *where*, not the *zygosity*.

---

## 7. Performance note: emission lookup + batched Viterbi (reverted, to redo)

A vectorization was attempted (precompute emissions in a `(n, a)` lookup table; batch the
Viterbi across samples) to make the larger taxa tractable. It was committed in a **wrecked,
unverified state and then reverted** (commits `8c1312e` → `d0f4377`); the current code is the
verified per-sample implementation. The reasoning, for when it's redone properly:

The two pieces are different:

- **Batched Viterbi (across samples)** — depth-independent; always valid and helpful.
- **Emission lookup table** — *correct* at any depth, but only *worthwhile* at low depth.

Why the table is a low-depth optimization, with the actual Zh numbers
(3.0M `(n,a)` cells, 61 × 49,002):

| total reads per cell | count |
|---|---|
| 0 | 2,194,103 |
| 1 | 624,770 |
| 2 | 138,552 |
| 3 | 25,493 |
| 4 | 4,868 |
| 5 | 981 |
| max depth | **36** |

- **Dedup payoff.** The speedup is computing `betabinom.logpmf` once per *distinct* `(n,a)`
  instead of once per cell × state. At depth ~1 there are only **77 distinct `(n,a)` pairs across
  3.0M cells** (~40,000× dedup). At high depth nearly every cell is unique and the table saves
  nothing.
- **Table size is `O(nmax²)`.** A dense grid is `(nmax+1)² × 3` entries: nmax=36 → ~4,100
  (trivial); nmax=1,000 → ~3M (24 MB); nmax=10,000 (deep WGS) → ~300M ≈ 2.4 GB. It blows up
  quadratically *because depth is the side length of the table.*

**Robust redo:** tabulate only the distinct pairs that actually occur (`np.unique` → those 77
rows) and index back via the inverse map. That keeps the dedup speedup, drops the quadratic
memory, and degrades gracefully at any depth. And **verify it reproduces D=0.043 + passes pytest
before committing.**

### 7.1 Relation to RTIGER's `getlogpsi` (same computation, opposite coverage regime)

This is **not a new idea** — it is exactly the structure RTIGER already uses. RTIGER fits a
BetaBinomial-HMM whose cost scales with the **number of distinct `(n, k)` = (total depth,
alt count) pairs, not the marker count**, and it memoizes the BetaBinomial over those distinct
pairs:

- **`getlogpsi`** (E-step) — memoized BetaBinomial log-pdf over the distinct `(n, k)` pairs.
  This is the direct analog of our emission lookup table.
- **`emissionUpdateState`** (M-step) — re-optimizes the dispersion τ by summing over the distinct
  pairs, **per state, every EM iteration**.

The difference is which end of the same `cost ∝ distinct pairs ∝ coverage` curve each sits on:

| | RTIGER runs (seqcapture/target sequencing) | nilHMM-counts (skim) |
|---|---|---|
| coverage | ~110× | ~1× (mean 1.25) |
| distinct `(n,k)` pairs | 7,555 @110× → 817 @30× → 66 @3× | **77** (Zh) |
| fit time | 240 s @110×, 18 s @30× | ~15 s (one Viterbi pass) |
| problem | too **many** pairs | none — already minimal |
| fix | **downsample / cap to ~20–30×** | nothing needed (already ≪ 20×) |

So the lab's RTIGER **20× cap** and this document's **lookup table** are two responses to the
*same* fact. RTIGER had too much coverage and capped it down to bound the pair count; skim is
~1×, so the pair count is tiny and the memoization is essentially free.

**One structural difference:** RTIGER runs full **EM**, so it pays the distinct-pairs sum in
*both* the E-step (`getlogpsi`) and the M-step (`emissionUpdateState`, every iteration × per
state) — which is the larger reason high coverage hurt it. nilHMM-counts here uses **fixed
emissions** (θ=[err,0.5,1−err], fixed `conc`) and a **single Viterbi pass — no M-step**, so it
pays the emission evaluation only once. If nilHMM-counts is ever fed deeper-than-skim data, the
same cap logic applies, and per the lab's architecture rule the thinning belongs in the **caller**
(zealtiger), not in this package. See the `zealtiger` memory `rtiger-betabinomial-cost`.

---

## 8. What's done vs. next

**Done:** env + Phase 0 bugfix; Phase 1 data assessment; count caller (`read_vcf_counts`,
`introgression_hmm_counts`, wrapper) + tests; per-taxon/check subset builder; Phase 2
calibration on Zh; B73 control validated.

**Next:**
1. **Per-taxon sweeps** for Zx / Zv / Zd / Zl — each has different B73 divergence, so re-fit `r`
   (and re-check effective `err`) per taxon, as RTIGER does. (Faster after the optimization redo.)
2. **Phase 3 — common-schema staging.** Run-length-encode per-marker calls → segments
   `name, chr, start_bp, end_bp, state(REF=0/HET=1/ALT=2)` →
   `results/sim_calibration/nilhmm_calls/` (mirror `stage_skimbin_calls.R`).
3. **Phase 4 — integrate** nilHMM into the cross-caller comparison
   (`skim_brbseq_skimbin.qmd`): painting lanes, confusion, Jaccard, intro-size histogram.
4. (Optional) **Phase 5 — donor founder matching** (which teosinte donor), using the teosinte
   reference panel — a capability the other callers lack.

---

## 9. Suggestions

- **Make the simulation reference the single source of truth for every `r`.** nilHMM is now
  calibrated against the BC2S2 sim; the other skim callers should target the same blessed
  reference so the four callers are compared on a common footing.
- **Report the het/hom split as low-confidence.** Given §6.3, downstream summaries (dosage,
  HET:ALT) should treat the donor *footprint* as the trustworthy quantity and flag the
  het-vs-hom partition as depth-limited. Consider reporting "donor fraction" (state>0) as the
  primary metric rather than separate HET/ALT.
- **Per-taxon, not per-accession.** Keep the RTIGER grouping (pool within a taxon's divergence;
  avoid mixing taxa in one emission fit). Singletons / tiny accessions don't get their own fit.
- **Redo the optimization defensively** (sparse-unique table + batched Viterbi) with a
  reproduce-D=0.043 gate, so the large taxa don't each take ~15–20 min/sweep.
- **Keep a clean B73 guard in every calibration**, as a hard constraint, not just a diagnostic.

## 10. Open questions

1. **Is the count path distinct enough from RTIGER to be worth a 4th caller?** Both are now
   BetaBinomial-HMMs on the same skim counts. The differences are nilHMM's explicit per-state
   error model and per-taxon priors. Is that a meaningful methodological contrast, or should
   nilHMM be positioned specifically for a **higher-coverage / imputed GT** dataset where its
   GT/error-model design is its real strength?
2. **Should we impute after all?** A reference panel (teosinte + founders, Grzybowski 2023) sits
   next to the cohort at the same sites. Imputing GT would let nilHMM run in its *native* GT mode
   and make it genuinely different from the count callers. Worth the extra pipeline step?
3. **Per-taxon `f_1`/`f_2`.** You mentioned having frequencies computed per donor species. Should
   the calibration use **observed per-taxon** `f_1`/`f_2`, or the **theoretical BC2S2** values
   (0.0625 / 0.0938)? The HMM is fairly prior-insensitive at this `r`, but observed freqs would
   bake in the real donor depletion / het excess.
4. **Does `r` need to be per-taxon at all,** or is the segment-size distribution similar enough
   across taxa that one `r` suffices? The sweep on each taxon will answer this — if the argmin is
   ~3e-5 everywhere, a single `r` is simpler and more defensible.
5. **Block-boundary definition.** Donor blocks currently use midpoints between flanking markers.
   Is that the right convention to compare against the simulation's true segment boundaries, or
   should we adopt whatever convention `stage_skimbin_calls.R` / the RTIGER segment objects use,
   for consistency across callers?
6. **Should the `err` ↔ segment-size tension be resolved differently** — e.g. a two-component
   objective (segment size **and** genotype frequency), or accepting that footprint is the only
   well-determined output? Right now we optimize size only and report the het-bias as a caveat.
