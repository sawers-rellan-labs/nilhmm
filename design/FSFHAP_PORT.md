# Native port of FSFHap → R/Rcpp (planning)

Porting **FSFHap** (Full-Sib Family Haplotype Imputation; Swarts et al. 2014,
*The Plant Genome* 7(3), doi:10.3835/plantgenome2014.05.0023) into the nilHMM
engine as a new caller. FSFHap is *not* a NIL introgression caller in the nNIL
sense — it is a **full-sib-family parental-haplotype reconstructor + imputer**.
It fits this package because TeoNAM and the nNILs are the **same BCₙSₙ lineage**
(donor × B73, backcrossed then selfed), differing only in generation depth, and
both organize into per-donor families that FSFHap pools.

**Method (chosen):** native port from spec + the TASSEL 5 Java source — **no git
subtree / no vendoring** (matches the `lbimpute` precedent). Java read over the
bitbucket API; no clone needed.

## Reference sources
- Paper: `agent/swarts2014.txt` (full text). Algorithm = M&M §"Viterbi Algorithm"
  / "Initializing the Matrices"; follows Andolfatto et al. 2011 / King et al. 2012.
- Java: `bitbucket.org/tasseladmin/tassel-5-source` →
  `src/net/maizegenetics/analysis/imputation/`. **Classes that actually carry the
  algorithm** (verified by reading, not just names):
  - **`FSFHapImputationPlugin`** — user entry point (`-FSFHapImputationPlugin`); wires
    `CallParentAllelesPlugin` (stage 1) → `ViterbiAlgorithmPlugin` (stage 3).
  - `CallParentAllelesPlugin` — driver: reads pedigree, filters to family,
    dispatches one parent-calling route.
  - **`ViterbiAlgorithmPlugin`** — FSFHap's imputation step: calls
    `imputeUsingViterbiFiveState` + `fillGapsInAlignment`.
  - `NucleotideImputationUtils` — **core**: `callParentAllelesByWindowForBackcrosses`
    (stage 1, BC), `callParentAllelesUsingClusters`, `callParentAllelesByWindow`,
    `imputeUsingViterbiFiveState` (stage 3, 5-state EM), `fillGapsInAlignment`,
    `HaplotypeClusterer`.
  - `TransitionProbability` (basic) / `TransitionProbabilityWithVariableRecombination`
    (distance-scaled); `ViterbiAlgorithm` (+ `…VariableStateNumber`) — decode.
  - **NOT on the FSFHap path** (driven by the separate `-ImputeProgenyStatesPlugin`,
    which needs an external parent-haplotype map): `BiparentalHaplotypeFinder` /
    `SelfedHaplotypeFinder`, `ImputeCrossProgeny`, `CrossProgenyEmissionMatrix`,
    `RephaseParents`.
  - `RephaseParents` — the single EM-like re-estimation pass.
  - `PopulationData` — pedigree-file parsing, family/contribution/F.
- Golden baseline for parity: `zealhmm agent/teonam_fsfhap.R` (installed
  `/Applications/TASSEL 5` FSFHap on TeoNAM families).

## Applicability & the reference-panel question
- **No donor/founder reference panel required.** Parental haplotypes are
  reconstructed *from the full-sib progeny themselves*: "Because FILLIN and FSFHap
  do not require known parental genotypes, these algorithms provide pedigree
  independent imputation." Consistent with nNIL's donor-agnostic philosophy — the
  donor is reconstructed, never looked up. (Resolves the standing concern: **no panel.**)
- **Hard input requirement: family grouping** (which lines share a donor cross), so
  their scattered introgressions collectively phase the donor haplotype.
- **TASSEL default is backcross-specific** (`useBCFilter = true` →
  `callParentAllelesByWindowForBackcrosses`, keying on sites segregating at the
  **0.25 ratio**) — a dedicated mode matching our BCₙSₙ design. Deeper BC depth is
  carried by the pedigree `contribution` columns (below).
- **TeoNAM** — clean fit: family membership known and correct.
- **nNIL** — applicable, caveat: ~1/3 of donor assignments disagreed with pedigree,
  so pedigree-based grouping is unreliable; run on confidently-grouped lines or
  expect a contaminated family. Likely why Jim used a *per-line* donor-agnostic HMM
  for nNILs and FSFHap for the family-structured TeoNAM.

### When `BiparentalHaplotypeFinder` is (NOT) the right tool — biological applicability
The `biparental` route reconstructs **two parental haplotypes** by clustering full-sib
progeny into two balanced groups per window. Two assumptions are load-bearing, and
both fail for NILs with heterozygous/outbred founders — use the **donor-agnostic HMM
(`nnil`/`emission_gt`) instead** in that case.

1. **Homozygous (inbred) founders.** Each parent must contribute ONE defined haplotype.
   A heterozygous/outbred founder contributes two (up to 4 haplotypes segregate), so the
   two-haplotype model breaks. Swarts 2014: *"identify exactly two parental haplotypes,
   ignoring any sites heterozygous in either parent … has NOT been tested for … outbred,
   heterozygous parents."* The finder discards het-in-parent sites → for a truly outbred
   donor (e.g. teosinte) that discards most informative sites. **Not useful.**
2. **Balanced full-sib segregation.** Windows need two substantial clusters (one per
   parental haplotype) and a seed window whose two majorities differ by `≥ 2·window−4`.
   **NILs violate this**: >90% of the genome is homozygous recurrent parent across all
   lines (one cluster → no seed), and the donor appears only as small, scattered
   introgressions in different places per line (no coherent donor haplotype to phase).

**The F2 cascade (task #5) is the mild preview of failure #1**: ~50% het fragments
window clusters → one candidate → a parent list empties → cascade. NILs would hit this
almost everywhere. **Validated behavior:** bit-exact vs TASSEL on RIL/inbred (its target)
and BC1; degrades on het-heavy F2; inappropriate for NILs.

| population | founders | segregation | right caller |
|---|---|---|---|
| RIL / inbred full-sib | homozygous | balanced | `fsfhap` (biparental) — bit-exact |
| BC1 (TeoNAM) | (recurrent majority) | backcross | `fsfhap` (BC route) — bit-exact |
| **NILs, esp. heterozygous/outbred donor** | heterozygous/diverse | recurrent-dominated | **`nnil` donor-agnostic HMM — NOT `fsfhap`** |

This is why nNIL/BzeaSeq (teosinte NILs) route through the donor-agnostic HMM: it needs
neither a reconstructed donor haplotype nor an inbred founder (the `nir` non-informative
rate absorbs het/diverse founders statistically), whereas `BiparentalHaplotypeFinder`
requires both.

**Primary-source confirmation (Swarts 2014, Test Datasets).** FSFHap is benchmarked on
**exactly one population: the NAM RILs** — 25 biparental full-sib families, ~200 **F6**
progeny each, *"highly inbred, average heterozygosity per line ≈ 0.001"*, with *"for each
family we expect MAF of 0.5"* (balanced 1:1 segregation). NILs are **never** benchmarked.
The paper's other datasets (diverse inbreds; outbred **diverse landraces**, het 0.052) are
FILLIN/Beagle comparisons, and the abstract states *"Beagle v. 4 is still preferable for
diverse heterozygous populations."* It even filters residual het out of the RIL benchmark
(subpop 1:1-segregation check + LD-with-neighbors + **removing RILs > 30% heterozygous**).
So: RIL/inbred + balanced segregation is FSFHap's validated domain (matches our bit-exact
RIL gold test); het populations and NILs are out of scope by the authors' own design —
use the donor-agnostic HMM for teosinte NILs.

## Pipeline (what `CallParentAllelesPlugin.performFunction` actually does)
Per family, chromosome-separate, markers sorted by position:

1. **Read pedigree → filter to family members** (`PopulationData.readPedigreeFile`,
   `FilterGenotypeTable`); optionally drop hets (`useHets`, default true).
2. **Call parent alleles** — one of three routes (mutually exclusive; default = BC):
   - **`callParentAllelesByWindowForBackcrosses`** *(default, `useBCFilter=true`)* —
     find sites segregating at ratio 0.25 (`whichSitesSegregateCorrectly`), drop
     same-tag SNPs (`whichSnpsAreFromSameTag`, R²>0.8), optional LD filter (minR>0),
     assign `alleleA`/`alleleC`, drop taxa <200 gametes coverage, recode to A/C/het.
   - **`callParentAllelesUsingClusters`** — filter by MAF/missing/het (maxHet 0.06);
     seed an LD window with ≤2 haplotype clusters (`HaplotypeClusterer`,
     `Haplotype.distanceFrom`, maxDiff ≤4); `extendClusters` bidirectionally; link
     windows via overlap similarity 0.8; subpopulations (NAM/monolithic); score
     sites by flanking consistency (halfWindowSize 20); refine at cutoffs
     [0.7, 0.8, 0.9].
   - **`callParentAllelesByWindow`** — window-LD variant.
   Then **`imputeUsingViterbiFiveState`** smooths the parent-origin track (5-state,
   below).
3. **Reconstruct/assemble parental haplotypes** (`BiparentalHaplotypeFinder` /
   `SelfedHaplotypeFinder`). **NOT on the backcross path — see route dispatch below.**

> **ROUTE DISPATCH (verified from `CallParentAllelesPlugin.performFunction`).**
> The four parent-calling routes are mutually exclusive, in this order:
> ```
> if (useWindowLD)                                          -> callParentAllelesByWindow
> else if (useClusterAlgorithm)                             -> callParentAllelesUsingClusters
> else if (useBCFilter && contribution1 ∈ {0.75, 0.25})     -> callParentAllelesByWindowForBackcrosses  [stage 1a/1b]
> else if (useMultipleBCFilter)                             -> callParentAllelesByWindowForMultipleBC
> else if (assignHaplotypesFromParents)                     -> UseParentHaplotypes
> else                                                      -> BiparentalHaplotypeFinder  [DEFAULT; stage 2]
> ```
> - The BC filter is gated on `contribution1 == 0.75 || 0.25` — i.e. **BC1 only**.
>   **TeoNAM (BC1S4, contribution 0.75)** takes the BC route → parent-calling is
>   **stage 1a/1b (done)**; `BiparentalHaplotypeFinder` is **never reached**, and the
>   critical path goes straight to stage 4 (`ImputeCrossProgeny`). In a BC1 the
>   recurrent parent is the majority at every segregating site, so `major=A(recurrent)`
>   / `minor=C(donor)` needs no windowed clustering.
> - **nNIL does NOT use FSFHap at all** (resolved 2026-07-07 from the pipeline): the
>   only FSFHap invocation in zealhmm is `agent/teonam_fsfhap.R` (TeoNAM); nNIL is
>   called by the **donor-agnostic HMM** (the `nnil`/`emission_gt` caller, already in
>   the package). So the BC2 routing question (`multipleBC` vs default
>   `BiparentalHaplotypeFinder`) is **moot for the current pipeline** — hypothetical
>   future work only.
>
> **Net: the sole real FSFHap consumer is TeoNAM (BC1) → BC route → stage 1a/1b (done)
> → stage 3 (`ViterbiAlgorithmPlugin` / `imputeUsingViterbiFiveState`, see step 4 below).
> `BiparentalHaplotypeFinder` (stage 2) AND `ImputeCrossProgeny` (4-state) are off every
> current critical path.** For reference, the BC2 fall-through route
> `callParentAllelesByWindowForMultipleBC` is NOT windowed clustering either — it is
> per-site `alleleA=major`/`alleleC=minor` assignment (like stage 1b) with a
> `whichSitesArePolymorphic(minMinorAlleleCount)` site filter instead of the 0.25
> binomial test. So even hypothetical nNIL-via-FSFHap would reuse stage 1's shape, not
> stage 2.
4. **Impute** (`ViterbiAlgorithmPlugin` → `imputeUsingViterbiFiveState`).
   **CORRECTION (verified 2026-07-07):** the FSFHap plugin does NOT use
   `ImputeCrossProgeny`/`CrossProgenyEmissionMatrix` (the 4-state, depth-dependent,
   phased-haplotype path — that belongs to a *different* plugin and would NPE on the
   GT-only TeoNAM HapMap). `FSFHapImputationPlugin.performFunction` is, unconditionally:
   ```java
   cpaResult = CallParentAllelesPlugin.performFunction(input);      // stage 1
   vap = new ViterbiAlgorithmPlugin(null);
   vap.setFillGapsInAlignment(fillgaps); vap.setProbHeterozygous(phet);
   vapResult = vap.performFunction(cpaResult);                       // stage 3 (below)
   ```
   and `ViterbiAlgorithmPlugin` does:
   ```java
   phet = (family.inbredCoef in [0,1]) ? (1 - family.inbredCoef)/2 : probHeterozygous; // DESIGN-DERIVED from pedigree F
   family.imputed = NucleotideImputationUtils.imputeUsingViterbiFiveState(tba, phet, family.name, useVariableTransition);
   if (fillGapsInAlignment) NucleotideImputationUtils.fillGapsInAlignment(family);
   ```
   - **`imputeUsingViterbiFiveState`** — **5-state Viterbi-training EM**: obs
     `AA→0 / AC(het)→1 / CC→2` on non-missing sites; init emission = Table 1 (below),
     init transition diag `{.999,.0001,.0003,.0001,.0005}` w/ `setAverageSegmentLength(chrLen/nSites)`,
     init `pTrue = {phom, .25·phet, .5·phet, .25·phet, phom}`, `phom=(1-phet)/2`.
     Loop (≤ **50** iters): Viterbi-decode all taxa → accumulate `transitionCounts[5][5]`
     (consecutive states) → `tp.setTransitionCounts(counts, chrLen, ntaxa)` and
     `emissionCounts[5][3]` (state×obs) → **re-estimate** `emissionProb[r][c] = counts[r][c]/rowsum`
     (**no pseudocount**); **converge** when the emission-count matrix equals the prior
     iteration's exactly. Output map **0→AA, 1-3→AC, 4→CC**. (Matches the parity log:
     "Iteration 0..3", "Imputation counts, rows=states, columns=observations".)
   - **`fillGapsInAlignment`** — per taxon forward-fill: consecutive non-`N` calls with
     the same value fill the intervening `N`; reset on a value change.
   - **`phet` is design-derived** from the pedigree `F` (`(1-F)/2`; BC1S4 F=0.9375 →
     0.03125), reinforcing the [design-routing dispatcher] — not a magic constant.
   - **END-TO-END TASSEL PARITY IS AVAILABLE HERE:** the imputed genotype HapMap
     (`imputed_genotypes_Chr*`) is exported, so unlike stages 1a/1b/2 (unloggable/
     unexportable intermediates), stage 3's per-cell output can be diffed against TASSEL.

> `ImputeCrossProgeny` / `CrossProgenyEmissionMatrix` (4-state, depth-aware) +
> `RephaseParents` are driven by a **different CLI tool — `-ImputeProgenyStatesPlugin`**
> (button "Progeny States"), VERIFIED: it does `new ImputeCrossProgeny(); icp.setHaplotypeMap(parentHaplotypeFilename); icp.improveImputedProgenyStates(); icp.imputeAll();`.
> That tool imputes progeny from an **externally supplied parental-haplotype map** +
> parentage file (hence its depth use, which would NPE on a GT-only HapMap). **No
> FSFHap option reaches it** — `-FSFHapImputationPlugin` uses `CallParentAllelesPlugin`
> + `ViterbiAlgorithmPlugin` (5-state EM). Do NOT port `ImputeCrossProgeny` for the
> FSFHap target. (This doc previously mis-attributed it to FSFHap by class-name
> inference — corrected 2026-07-07 by reading the plugin dispatch.)

## State spaces
> **CORRECTION:** an earlier version framed this as "5-vs-4-vs-3, reproduce all". That
> was wrong. The FSFHap path uses only the **5-state** imputation (`imputeUsingViterbiFiveState`,
> stage 3) projected to the **3-state** REF/HET/ALT output schema. The **4-state**
> cross-progeny model is `-ImputeProgenyStatesPlugin`'s, NOT FSFHap's — it is not
> reproduced. The 5-state and 3-state entries below are current; the 4-state entry is
> retained only as a contrast/reference.
- **Parent-allele smoother — 5 states** (`imputeUsingViterbiFiveState`), residual-het
  / bulked frequency classes `{all A, 3A:1C, 1A:1C, 1A:3C, all C}`. Emission (Table 1,
  rows = state, cols = obs {homA, het, homC}), **verbatim**:
  ```
  {{.998,.001,.001},
   {.6, .2, .2},
   {.4, .2, .4},
   {.2, .2, .6},
   {.001,.001,.998}}
  ```
- **Cross-progeny phasing — 4 states** (`CrossProgenyEmissionMatrix`), parental
  haplotype combinations `{hap0/hap0, hap0/hap1, hap1/hap0, hap1/hap1}`;
  `missingState = 4`. Emission is **coverage-aware**: het-called-as-hom
  `probHetAsHom = (depth==1) ? 0.99 : pow(0.5, totalDepth−1)`; homozygous states
  `0.998` match / `0.001` otherwise; het-mismatch floor `0.0001`.
- **Common schema — 3 states** (REF/HET/ALT) is only the **output projection**:
  `A/homA → REF(0)`, hets/intermediate → `HET(1)`, `C/homC → ALT(2)` with A = B73
  (recurrent), C = donor. Output = standard segment schema
  `source, donor, name, chr, start_bp, end_bp, state`.

> **Design consequence:** faithful reproduction ports the **5-state smoother AND the
> 4-state phaser**; reusing the engine's **3-state emission** would skip parental
> reconstruction and is a *different* algorithm. Reuse only the **Viterbi kernel**
> (`viterbi.cpp` / `viterbi_batch.cpp`) as a primitive, driven by FSFHap's own
> emission/transition. `forward_backward.cpp` is **not** on the critical path.

## Transition probability
- **Basic** (`TransitionProbability`, used in `ImputeCrossProgeny`): 4×4 with
  `avgTP = 1/(nNonMissing − 1)` off-diagonal scale; diagonal = 1 − Σ off-diag.
- **Distance-scaled** (`TransitionProbabilityWithVariableRecombination`, when a recomb
  map is supplied), per interval `node`:
  ```
  adjustedProbability[i][j] = segmentLength · rrr · transitionCounts[i][j] / numberOfTaxa
  adjustedProbability[i][i] = 1 − Σ_{j≠i} adjustedProbability[i][j]
  ```
  `segmentLength = positions[node] − positions[node−1]`; `rrr` = mean relative
  recombination rate over the interval; `transitionCounts` = **empirical base matrix
  from NAM** (chromosome-specific); `numberOfTaxa = 195` normalization.
- Same family as `lbimpute`'s distance transition (linear in distance × local recomb
  rate) — reuse that `unit`/`recombdist` parameterization; the empirical NAM
  `transitionCounts` is the one asset we can't reuse and must decide whether to port,
  re-derive from our map, or replace with a parametric form.

## Input formats & IO separation (design principle — REQUIRED)
**The FSFHap algorithm implementation is independent of input format.** Format
knowledge (TASSEL, PLINK, VCF, pedigree TSV) lives *only* in a thin IO/adapter
layer that **normalizes to the engine's canonical data structures**; the core
(`src/fsfhap.cpp`, `R/fsfhap.R`) never parses a file and sees only those normalized
structures — exactly like every other caller. This mirrors the package's existing
pattern: `read_counts()` folds `tsv` / `gatk_table` / `vcf_ad` into **one** long
table, and the emission/engine code is oblivious to which format it came from.

**Canonical internal representation (the "normalizing convention" — reuse, don't
invent):** the same tables/matrices the rest of the machinery uses.
- Genotype calls: the `read_vcf_gt()` long table **`(name, chr, pos, g)`**,
  `g ∈ {0 REF-hom, 1 het, 2 ALT-hom, 3 missing}` (A/H/B → 0/1/2/NA).
- Optional allelic depth: the `read_counts()` long table
  **`(name, chr, pos, n_ref, n_alt)`** — FSFHap's coverage-aware emission
  (`probHetAsHom = pow(0.5, depth−1)`) and `RephaseParents.rephaseUsingAlleleDepth`
  consume depth when present; without it, fall back to a fixed het-as-hom rate.
- Working matrix: **one chromosome at a time, markers sorted by position**, as a
  per-**family** `taxa × markers` matrix — the same 2D per-chromosome layout the
  batched decoder already uses (cf. the count caller's `ref[n_samples, n_markers]`
  and `viterbi_batch.cpp`). FSFHap's family-level stages (parent-allele calling,
  haplotype clustering) simply need the *family* slice of that matrix, not a new
  structure.
- Grouping + priors: a per-line **`family`** key and design **`priors`** (see below).

**Supported inputs, both via adapters (never inside the algorithm):**
- **TASSEL** — HapMap / VCF genotypes + FSFHap pedigree TSV.
- **PLINK** — `.ped`/`.bed` (+`.bim`/`.fam`) genotypes; `.fam` FID as the family key.

**Decision: adapters live in `io.R`, minimal surface.** Prefer *extending existing
readers* (e.g. a new `format =` on `read_vcf_gt()` / `read_counts()`, plus a single
`read_pedigree()` for family/priors) over adding one function per format. The
pipeline still owns paths/sample lists per the data-agnostic rule; `io.R` owns only
format → canonical-table conversion.

### Pedigree / family-grouping normalization
FSFHap's own pedigree file (`PopulationData.readPedigreeFile`) is a **tab-delimited,
header-first TSV** — **not PLINK `.fam`**:

| col | field | notes |
|---|---|---|
| 1 | family name | must be non-empty, not "NA"; rows accumulate by family |
| 2 | individual ID | |
| 3 | parent 1 | |
| 4 | parent 2 | |
| 5 | contribution 1 | double (expected genome fraction from parent 1) |
| 6 | contribution 2 | double |
| 7 | inbreeding coef F | double, optional |

- Columns 1–4 overlap PLINK `.fam` (family, individual, two parents) but **5–7 do
  not** (PLINK cols 5–6 are sex/phenotype). The `contribution` + `F` columns encode
  the **backcross depth and selfing level** — i.e. exactly the `design_priors`
  (`f_1`/`f_2`) inputs. FSFHap's pedigree is *richer* than `.fam` and carries the
  priors directly.
- **The algorithm never reads any of these files.** Contract:
  `call_ancestry(caller = "fsfhap", data, family = <vec/col>, priors = ...)`
  — normalized `data`, a per-line `family` key, and contribution/F-derived `priors`.
  The IO layer (**zealhmm** pipeline, or `io.R` reader) builds `family`/`priors` from
  an FSFHap pedigree TSV, a PLINK `.fam` FID, or the `PN<n>` sample-name prefix
  (`PN1_SID37 → Zx.0120_P1_P2_P2.1.1.1.1`; `PN` = family, `Zx/Zd/Zv/Zl` = donor
  species, trailing `.1.1.1.1` = selfing path).
- Note: no `.fam`/`.ped` is checked into either repo; `genetic_maps/*.map` are
  genetic maps (`marker_id`, cM) and `K_snp.rel.id` is a PLINK2 `--make-rel` kinship
  ID sidecar — neither is a parentage pedigree.

## Stage 2–3 internals (locked from Java)

### `HaplotypeClusterer` (util building block for stage 2)
- A **Haplotype** = byte array of allele calls; **distance** = count of allele
  differences, **excluding sites missing (`N`) in either** haplotype.
- **Greedy/sequential** clustering: a cluster groups haplotypes at **0 pairwise
  distance**; because of missing data a haplotype may join multiple clusters, each
  membership scored `1 / (#clusters it belongs to)`.
- **Merge** (`getMergedClusters`): sequential — if two clusters merge, the merged
  one becomes the head for testing the rest; merge when **max pairwise difference ≤
  `maxdiff`**. `moveAllHaplotypesToBiggestCluster(maxdiff)` consolidates outliers;
  `sortClusters()` orders by cluster score; consensus via `getMajorityHaplotype()` /
  `getUnanimousHaplotype()`.
- Port: a small standalone kernel; the metric is trivial, the greedy merge is the
  only subtlety (order-sensitive → replicate the sequential head rule for parity).

### `BiparentalHaplotypeFinder` (stage 2 — parental haplotype reconstruction)
- **Seed**: slide windows (`window = 100` SNPs, step `= window − overlap = 75`,
  `overlap = 25`) until a window has **exactly two major clusters**, any third
  cluster `< minClusterSize (3)`, and the two haplotypes differ by **≥ `2·window − 4`**
  sites. Seed `h0/h1` = the two clusters' majority haplotypes → `parentHaplotypes`
  (two lists, one per parent).
- **Per window**: `clusterWindow()` keeps taxa with ≥ `minNotMissingProportion
  (0.2)` non-missing → cluster → merge if ≤ `maxDifferenceScore (0)` →
  `moveAllHaplotypesToBiggestCluster` → `removeHeterozygousClusters`.
- **Assign** candidate haplotypes to a parent via `doesOverlapMatch()` (overlap
  region matches with ≤1 mismatch); else nearest parent by min distance.
- **Extend** forward then backward from the seed; `updatePopulationDataAlleles()`
  fills the non-overlap portion and resolves ambiguous/shared-allele sites in the
  overlap; phase is carried by overlap continuity.
- Filters: `minR2 0.2` (50-site LD window), `minMaf 0.05`, `minCoverage 0.2`,
  `maxHetDeviation 5`.
- **Key finding: no BC/F₂-specific logic here** — this class assumes a biparental
  cross with two distinct parental haplotypes. The backcross handling lives
  *upstream* in `callParentAllelesByWindowForBackcrosses` (allele A/C assignment on
  0.25-segregating sites). **Implication for our BCₙSₙ data:** by the time we reach
  haplotype finding, one "parent" track is ~the B73 recurrent consensus and the
  other is the donor mosaic being reconstructed — an asymmetry we can exploit, and
  a risk (a rare donor haplotype may fail the seed's two-cluster / `2·window−4`
  divergence test at deep BC — ties to the min-family-size open question).

### `RephaseParents.rephaseUsingAlleleDepth` (stage 3 re-estimation — the "EM" pass)
- **Inputs**: progeny imputed states `byte[]` (0–3, `4`=missing), allele depths
  `origGeno.depthForAlleles(taxon, site) → int[6]`, prior parental haplotypes
  (`startingParents`).
- **Accumulate**: `stateAlleleDepths[4][6]` (by progeny state) and
  `parentStateAlleleDepths[2][6]` (by focal-parent state).
- **Update**: Bayesian over the 4 parental-haplotype-pair combinations, each term a
  `BinomialDistribution(...).probability(observed_minor_count)` weighted by
  major/minor allele frequency; posterior
  `haplotypeProbability[0][s] = ph00major / (ph00major + ph00minor)`.
- **Runs once** per invocation (no internal loop); feedback into
  `imputeCrossFromParentsUsingProbabilities()` is external (the driver's single
  improvement iteration).
- **Depends on allelic depth** → only meaningful when the canonical `n_ref/n_alt`
  table is present; **skip rephase when depth is absent** (consistent with the
  "optional depth" IO note). This is the only depth-dependent stage.

### Consensus & overlap helpers (locked)
- **`HaplotypeCluster.getMajorityHaplotype()`** — per site, `Multiset<Byte>` of
  non-`N` alleles across members; consensus = **max-count allele**; **ties → `N`**;
  all-`N` → `N`. (Seed haplotypes come from this.)
- **`getUnanimousHaplotype()`** — per site, distinct non-`N` alleles; assign only if
  **exactly one** (`nNucleotidesObserved == 1`); polymorphic or all-`N` → `N`.
- **`BiparentalHaplotypeFinder.doesOverlapMatch(h0,h1,overlap,forward)`** (static):
  compares `overlap`-length ends — forward = end(h0) vs start(h1); backward =
  start(h0) vs end(h1). Non-`ACGT-+` chars → `N`; **mismatches involving `N` don't
  count**; match iff **mismatches < 2** (≤1 real mismatch).
- **`updatePopulationDataAlleles(parentHaplotypes, positions, start, length)`** —
  where the two parent tracks agree (`seqA[ptr] == seqC[ptr]`), set
  `alleleA = alleleC = NN` (missing). **This is FSFHap's non-informative-site
  handling** (parent-of-origin undeterminable where parents share the allele) — the
  structural analog of nNIL's `nir`.

### Decode kernels & segregation test (locked)
- **`imputeUsingViterbiFiveState`** — obs encoding `{AA→0, AC/CA→1 (het), CC→2}`;
  5×3 emission = Table 1 (verbatim above); transition = basic `TransitionProbability`
  with `tp.setAverageSegmentLength(chrLength / nSites)` and diagonal seed
  `{.999,.0001,.0003,.0001,.0005}`; **initial dist**
  `phom = (1−pHet)/2; pTrue = {phom, .25·pHet, .5·pHet, .25·pHet, phom}` (`pHet` =
  the `-phet` flag, **0.03125** in the TeoNAM run); decodes with the **standard
  `ViterbiAlgorithm`**; output map **state 0→A, states 1–3→het, state 4→C**.
- **`ViterbiAlgorithmVariableStateNumber`** — log-space max-product; per-node state
  count via `transitionMatrix.setNode(node)` → `getNumberOfStates()`.
  Init `distance[i] = lnEmission(i,obs[0],0) + ln pTrue[i]`; step
  `cand[i][j] = distance[i] + lnTransition(i,j) + lnEmission(j,obs[node],node)`;
  backtrace via a per-node `history`. **Tie-break: strict `>` in argmax ⇒ lowest
  state index wins.** *Parity: **confirmed** — our `viterbi.cpp` uses strict `>` in
  both the transition argmax (`:48`) and final-state pick (`:59`), so it already
  matches the earliest-index rule.* (The variable-state-number aspect — per-node
  `getNumberOfStates()` — is unused by our fixed 4/5-state chains; reuse `viterbi.cpp`.)
- **Both Viterbis decode the non-missing observation subsequence only**
  (`nonMissingObs`); missing sites are **not** modeled in the chain — they are
  restored afterward by `fillgaps`. Port accordingly: filter → decode → gap-fill.
- **`whichSitesSegregateCorrectly`** — no chi-square; a **binomial dominance test**
  (log-binomial via `GammaFunction.lnGamma`). Per site compute
  `pmono = binomial(Mj+Mn, Mn, 0.002)` (monomorphic/error model, err const 0.002),
  `pquarter` (ratio 0.25, backcross), `phalf` (ratio 0.5, F2). Keep a site iff
  biallelic (`freq[1].length > 1`) and `pMissing ≤ maxMissing` and:
  **backcross → `pquarter > phalf && pquarter > pmono`**; F2 → `phalf/(pmono+pquarter) > 2`.
  Deterministic ⇒ a tight per-stage parity target.

**Constants to carry into `src/fsfhap.cpp`:** `window=100`, `overlap=25`,
`step=75`, `minClusterSize=3`, `minNotMissingProportion=0.2`, `maxDifferenceScore=0`,
seed-divergence `≥2·window−4`, `minR2=0.2`, `minMaf=0.05`, `maxHetDeviation=5`;
five-state: emission Table 1, `pHet=0.03125` default, diagonal seed
`{.999,.0001,.0003,.0001,.0005}`, `avgSegmentLength = chrLen/nSites`; segregation:
error const `0.002`, ratios `0.25` (BC) / `0.5` (F2); overlap-match `< 2` mismatches
(N-tolerant); majority-consensus ties → `N`.

## Engine integration
- **`caller_spec("fsfhap", …)`** (`R/callers.R`): needs its own `R/fsfhap.R` driver
  (like `rtiger`/`lbimpute`), not the shared `emission × duration` combinator — it
  carries multi-stage state spaces, a required `family` grouping, and haplotype
  reconstruction the combinator can't express.
- **`design_priors()`** ← pedigree contribution/F → expected het/donor freqs; TeoNAM
  vs nNIL is a preset difference, not new code.
- **Reuse:** `viterbi.cpp` / `viterbi_batch.cpp` (decode primitive); `lbimpute`
  distance-transition parameterization.
- **New C++ (`src/fsfhap.cpp`) — the FSFHap path (done):** BC-ratio parent-allele
  calling (stage 1a/1b) + the 5-state EM Viterbi imputation + `fillGapsInAlignment`
  (stage 3). The 4-state cross-progeny emission / rephase and biparental haplotype
  assembly are **NOT** on the FSFHap path (see route dispatch + stage-3 correction);
  `HaplotypeClusterer` (stage 2a) is built but off-path.

## Defaults observed in source (for the caller signature)
`windowSize=50`, `minRforSnps=0.2`, `maxMissing=0.9`, `minMinorAlleleFrequency=-1`
(negative ⇒ use expected segregation ratio; cluster route falls back to 0.05),
`maxHetDev=25`, `minUsedClusterSize=5`, `maxDifference=0` (driver) / haplotype
distance `maxDiff=4` (cluster route), `useBCFilter=true`, `useMultipleBCFilter=false`,
`useClusterAlgorithm=false`, `useWindowLD=false`, `useHets=true`; cluster route:
`maxHet=0.06`, overlap similarity `0.8`, tag-SNP R² `0.8`, refine cutoffs
`[0.7,0.8,0.9]`, `halfWindowSize=20`, `minMatchScore=0.9`; BC route: `minGametes=200`.

## Validation plan (tiered)

Reference TASSEL command (from `zealhmm agent/teonam_fsfhap.R`) — the frozen
invocation every gate compares against:
```
-FSFHapImputationPlugin -pedigrees <ped> -bc true -phet 0.03125 -fillgaps true \
  -maxMissing 1.0 -minMaf 0.0 -minR 0.0 -exportType Hapmap
```

**Scaling principle — never gate on hour-long runs.** Every tier runs on a
**small fixture first** (one family, one or a few chromosomes, thinned markers),
which is what CI/`devtools::test()` executes; the full-family / full-genotype run
is a **separate, manual extension** (a `\dontrun` harness in
`agent/`, like `nilhmm_validate_full_cohort.R`), not part of the routine test
suite. All TASSEL reference outputs are **precomputed once and frozen** under
`tests/fixtures/`, so validation never re-invokes the JVM.

### Tier 0 — Speed benchmark (do FIRST)
Head-to-head wall-clock vs TASSEL FSFHap on the **same input**, since "speed/
control" is a stated motivation and must be quantified before correctness work.
- Metric: seconds per family, scaled small → large (1 chr thinned → 1 chr full →
  all chr → full family set); report throughput (taxa·markers/s) and peak memory.
- TASSEL side timed as in `teonam_fsfhap.R` (`system.time` around `run_pipeline.pl`,
  fixed `-Xmx`); port side timed identically. Record the JVM startup overhead
  separately (fixed cost TASSEL pays per invocation, the port avoids).
- Output: a small table in this doc / `benchmarks/benchmark_fsfhap.R`, refreshed as the
  port matures.

**Head-to-head captured 2026-07-07** (`benchmarks/benchmark_fsfhap.R`, synthetic BC1S4
single family; TASSEL 5.2.96, Java 17, Apple-silicon 10-core) — TASSEL vs the
nilHMM `fsfhap` caller (BC route), **now with cross-family/chromosome parallelism**
(`threads=`; `.fsfhap_states` fans out the independent family × chromosome units
via `parallel::mclapply`). Port timed both single-threaded and at the effective
unit count (= n_chr for this single-family ladder):

| rung | lines × markers | cells | TASSEL total | TASSEL compute | port 1-thread | port N-thread | port÷compute (best) |
|---|---|---|---|---|---|---|---|
| xs | 20 × 500 (1 chr) | 10k | 0.37 s | 0.07 s | **0.02 s** | — | 3.2× |
| s | 50 × 5k (1 chr) | 250k | 1.02 s | 0.66 s | **0.60 s** | — | 1.1× |
| m | 100 × 25k (5 chr) | 2.5M | 8.6 s | 8.1 s | 5.20 s | **2.31 s** (5t) | **3.5×** |
| full | 100 × 51k (10 chr) | 5.1M | 16.8 s | 16.1 s | 10.9 s | **3.85 s** (10t) | **4.2×** |

**Reading:** the full-genome gap is not just closed but reversed — the port now
runs full genome in **3.85 s vs TASSEL's 16.1 s compute (4.2× faster)**. Two
changes did it, both landed here:

1. **A serial `order()` hotspot removed.** Profiling the (initially poor) parallel
   scaling exposed that a character multi-key `order()` — in `.fsfhap_states`'s
   final sort *and* in the shared `to_segments()` RLE — was ~half of full-genome
   wall-clock. Switching both to radix on integer factor codes (the same fix
   `.lbimpute_prep` already carried) cut full-genome **single-threaded** time from
   34 s → **10.9 s** on its own — so the port beats TASSEL even serial. (This sped
   up every caller's `to_segments`, not just fsfhap.)
2. **Coarse-grained parallelism.** The family × chromosome units are fully
   independent (chromosome-separate processing, family pooling), so they map
   cleanly onto `mclapply`. Full genome scales **2.8× on 10 units**; the remaining
   gap from linear is Amdahl's serial tail (data prep + the radix sort + mclapply
   IPC), not load imbalance (equal-sized chromosomes). Results are **bit-identical
   to `threads=1L`** (verified: `identical(serial, parallel)` + full suite green).

Real populations widen the win further: TASSEL pays JVM startup **per family**
(~0.5 s each; TeoNAM ~5, nNIL ~18) while the port loads once, and a multi-family
cohort gives `mclapply` more independent units to fill all cores (here a single
family caps units at n_chr = 10). Memory: TASSEL ~0.5 GB peak; the port is lighter
(no JVM heap; `mclapply` forks are copy-on-write). **Net: the port is now faster
than TASSEL both single- and multi-threaded**, on top of dependency removal +
common-schema integration.

**Real-data confirmation (2026-07-07) — TeoNAM TIL01 family, full 51K genotypes.**
Run through the env-var path `FSFHAP_HAPMAP=…/geno_gwas_51k.hmp.txt
FSFHAP_PED=…/pedigree_TIL01.txt Rscript benchmarks/benchmark_fsfhap.R` (the harness now
loads the real genotypes via `read_hapmap`/`read_pedigree`, filtered to the
pedigree's taxa, and derives the design from the pedigree contribution/F →
**BC1S4**; TASSEL runs on the same HapMap + pedigree — a true head-to-head):

| dataset | taxa × markers (chr) | units | cells | TASSEL compute | port 1-thread | port 10-thread | port ÷ TASSEL |
|---|---|---|---|---|---|---|---|
| TIL01 (1 family) | 220 × 51,004 (10) | 10 | 11.2M | 24.3 s (462k c/s) | 12.3 s | **4.96 s (2.26M c/s)** | **4.9×** |
| all 5 families | 1,237 × 51,004 (10) | 50 | 63.1M | 218.5 s (289k c/s) | 78.9 s | **28.6 s (2.21M c/s)** | **7.6×** |

The real numbers confirm the synthetic conclusion and are **stronger**. Single
family (TIL01): the port is **2.0× faster single-threaded** and **4.9× on 10
threads**, 2.5× over its own serial across 10 chromosome units. Full cohort (all
five TeoNAM families, grouped by taxon-name prefix `^([A-Za-z]+[0-9]+)` → TIL01/
03/11/14/25): the margin **widens to 7.6×** (28.6 s vs 218.5 s) — 50 independent
family × chromosome units fill the cores better (2.8× self-scaling) *and* TASSEL's
throughput degrades at scale (289k vs 462k cells/s) while the port holds steady
(~2.2M cells/s at both scales). `read_hapmap` IO is 6–17 s (one-time, excluded from
compute). Peak RSS: TASSEL 0.7–1.1 GB, port lighter. Reproduce via the env-var
path (pedigree optional — falls back to taxon-prefix grouping + a synthesized
TASSEL pedigree):

```
# single family (pedigree-driven; design derived from contribution/F)
FSFHAP_HAPMAP=…/geno_gwas_51k.hmp.txt FSFHAP_PED=…/pedigree_TIL01.txt \
  Rscript benchmarks/benchmark_fsfhap.R
# all families (taxon-prefix grouping; design = FSFHAP_DESIGN, default BC1S4)
FSFHAP_HAPMAP=…/geno_gwas_51k.hmp.txt Rscript benchmarks/benchmark_fsfhap.R
```

### Tier 1 — Port faithfulness: stock TASSEL FSFHap parity (acceptance gate)
"Does the port reproduce the reference on identical input?" — analogous to `rtiger`
vs `calls_taxa_r5.csv` and `nnil` vs the frozen Python baselines.
- **Golden baseline** = `teonam_fsfhap.R` output: `imputed_<fam>.hmp.txt` + the
  parent-call HapMap, per family, with the flags above.
- **Fixture**: skeleton HapMap + 7-col pedigree + TASSEL outputs for a **small**
  family slice (then 2–3 families spanning Zx/Zd/Zv/Zl for the extended run),
  frozen under `tests/fixtures/` with SHA256SUMS.
- **Compare (FSFHap BC path)**: stage-1 segregating-site set + same-tag cardinalities
  (count parity, done) → then the **imputed parent calls** end-to-end (stage 1 →
  5-state EM → fillgaps) vs TASSEL's `imputed_parents` HapMap.
- **Bar**: **ACHIEVED — per-cell concordance 1.00000** on the e2e fixture
  (`test-fsfhap-e2e-parity.R`); the BC-ratio filter, same-tag R², recode, EM Viterbi,
  and gap-fill are all bit-exact. (No 4-state/rephase — off the FSFHap path.)

### Tier 2 — Imputation accuracy: masked-genotype recovery (self-contained)
Swarts 2014's own validation; catches "port == TASSEL but both wrong". Mask a
fraction of high-confidence observed calls on real TeoNAM, impute, measure error on
the held-out calls. Run **port vs TASSEL on the identical mask** → head-to-head
accuracy, not just mutual agreement. No external truth needed; small-fixture-first.

### Tier 3 — Biological ground truth: simulated full-sib family
Only truth that scores breakpoint *placement*. **Gap:** existing `sim_truth/`
(`bc2s3_nil_segments.csv`) is independent NIL segments, **not** a segregating
full-sib family — FSFHap needs parents + progeny with a shared, phaseable donor
haplotype and known breakpoints. Requires extending the simulation (new fixture,
tracked below). Score: donor-fragment Dice, breakpoint distance, marker recall
(the `simulation-calibration.qmd` objective).

### Tier 4 — optional: high-density chip
If any TeoNAM lines have chip genotypes, use as orthogonal truth; likely nNIL-only,
so Tiers 0–3 stand alone.

## Port status (all TBD)
- [x] Read + lock stage-2/3 internals: `BiparentalHaplotypeFinder`,
      `HaplotypeClusterer`, `RephaseParents` (see "Stage 2–3 internals" section).
- [x] Read + lock decode/segregation/consensus bodies: `imputeUsingViterbiFiveState`,
      `ViterbiAlgorithmVariableStateNumber`, `whichSitesSegregateCorrectly`,
      `getMajorityHaplotype`/`getUnanimousHaplotype`, `doesOverlapMatch`,
      `updatePopulationDataAlleles` (see "Decode kernels & segregation test" and
      "Consensus & overlap helpers"). **Java reading complete — ready for C++.**
- [~] `src/fsfhap.cpp` — parent-allele calling, BC-ratio route:
  - [x] **stage 1a — segregating-site test** (`fsfhap_segregating_sites_cpp`, faithful
        `whichSitesSegregateCorrectly` + explicit-lnGamma binomial; `R/fsfhap.R`
        `.fsfhap_segregating`; `tests/testthat/test-fsfhap.R`, 11 checks, MAF-parameterized).
  - [x] **stage 1b — BC parent-allele calling** (`fsfhap_same_tag_keep_cpp`, faithful
        `whichSnpsAreFromSameTag` + `calculateRSqr` on presence contingencies;
        `R/fsfhap.R` `.fsfhap_call_parents_bc` orchestrates seg∧same-tag → A/C
        major/minor assignment → <=200-gamete coverage filter → A/C/het recode into
        the parent-origin frame). +10 tests (21 total). **`ldfilter` (min_r>0)
        deferred** — TeoNAM run uses min_r=0.
  - [x] **TASSEL parity for 1a/1b** (`agent/fsfhap_tassel_parity.R`): matches TASSEL's
        logged cardinalities on identical input — `polybits` (seg), `filteredBits`
        (seg ∧ same-tag), `called snps`. **Exact match** on a 60×580 fixture with the
        same-tag filter firing (569→421). Frozen to `tests/fixtures/fsfhap_tassel/`
        + `test-fsfhap-parity.R` as a no-JVM CI gate. NOTE: per-site matrix diff is
        impossible via CLI (`CallParentAllelesPlugin` emits an unexportable
        `PopulationData`); counts pin the seg-binomial + presence-R² logic. Coverage
        filter / recode values remain gated on unit tests only.
- [~] `src/fsfhap.cpp` — `HaplotypeClusterer` + biparental haplotype assembly (**hardest**).
      **REAL FSFHap route, DEPRIORITIZED — off every current critical path.**
      `BiparentalHaplotypeFinder` IS reachable from FSFHap: it's the DEFAULT `else` branch
      of `CallParentAllelesPlugin` (stage 1), reached for non-BC1 families (contribution ∉
      {0.75,0.25}) with no cluster/windowLD/multipleBC flag. NOT a misattribution (unlike
      `ImputeCrossProgeny`, which is `-ImputeProgenyStatesPlugin`-only and off FSFHap
      entirely). It's just off OUR path: TeoNAM = BC1 (BC route, done); nNIL ≠ FSFHap.
      Stage 2a is banked (a genuine dependency of `BiparentalHaplotypeFinder` AND
      `callParentAllelesUsingClusters`); 2b needed only to support non-BC1 biparental (F2/RIL).
  - [x] **stage 2a — clustering primitive** (`fsfhap_cluster_window_cpp`): faithful
        `makeClusters` (two-pass 0-distance clustering, multi-membership via missing,
        fractional `1/count` scores), `Haplotype.getDistance` g-code metric
        (het/hom=1, hom/hom=2, N=0), majority/unanimous consensus, score↓/size↓ sort.
        `test-fsfhap-cluster.R` (20 checks: multi-membership, het-distance, missing
        tolerance, consensus). majority==unanimous for 0-distance clusters (diverge post-merge).
  - [x] **stage 2a-complete — clusterer ops** (`fsfhap_cluster_window_cpp` flags
        `merge`/`move_biggest`/`max_het`): faithful `mergeClusters` (sequential head-merge,
        `clusterDistanceMaxPairDiff` on non-shared cross pairs, `mergeTwoClusters` union),
        `moveAllHaplotypesToBiggestCluster` (unanimous-distance ≤ maxdiff),
        `removeHeterozygousClusters` (`countHeterozygousSites` > maxHet: majority-het OR
        2nd-allele-count > 1), `recalculateScores`. +8 tests (merge/move at maxdiff, het drop).
  - [x] stage 2b — **all bodies read verbatim** (spec captured, ready to port):
        `assignHaplotyes` seed + fwd/bwd extend loops; `clusterWindow`; `getParentHaplotypes`
        (doesOverlapMatch, else nearest-parent by overlap distance; equidistant → dropped);
        `mergeMajorHaplotypes` (`getCensoredMajorityHaplotype(maxMaf=0.2, maxMinorCount=2)`:
        major unless `maf>maxMaf && minorCount>maxMinorCount` → N; merge clusters within
        maxDistance=4, minClusterSize=3); `updatePopulationDataAlleles` (nhap==1/1: seqA==seqC→NN,
        else non-het allele each; multi-hap: TreeSet size logic, mostly → NN).
  - [x] stage 2b — `preFilterSites` step 1: **`filterSnpsByTag` + `computeRForMissingness`**
        (`fsfhap_filter_snps_by_tag_cpp`): same-tag thinning by presence phi-correlation
        (`<0.7` within 64 bp) + MAF/missing/het gates, head advances per kept site.
        `test-fsfhap-prefilter.R` (7 checks). DECISION 2026-07-07: **TASSEL-faithful** LD
        (port the subsystem, NOT substitute `fastindep` — different criterion, breaks parity).
  - [x] stage 2b — **`preFilterSites` COMPLETE** (`fsfhap_prefilter_sites_cpp`): all 4 steps
        faithful — filterSnpsByTag → het-deviation (`mean + maxHetDeviation·sd`, sample sd) →
        biallelic → LD (per-site avg r² over ±50-filtered-site window, `HetTreatment.Homozygous`
        = hets→missing via `ld_rsqr_hom` reusing `calc_rsqr`, reject `< minR2`; NaN skipped,
        all-NaN not a rejection). `test-fsfhap-prefilter.R` (+6 checks). Suite 629/0.
  - [x] stage 2b — **`BiparentalHaplotypeFinder.assignHaplotyes` DONE** (`fsfhap_biparental_alleles_cpp`):
        seed scan → fwd/bwd extension maintaining `parentHaplotypes`, per window
        `cluster_window` → `merge_major_haplotypes` (`getCensoredMajorityHaplotype`) →
        `get_parent_haplotypes` (`doesOverlapMatch` + nearest-parent) → `update_pop_data_alleles`
        → per-site `alleleA`/`alleleC`. `test-fsfhap-cluster.R` covers the clusterer ops.
  - [x] stage 2b — **non-BC1 route WIRED + F2 END-TO-END PARITY**: `.fsfhap_biparental_call`
        (preFilterSites → biparental alleles → recode → coverage) dispatched by design
        (BC1→bc, else→biparental). **EXACT site-set match (500/500 markers) AND per-cell
        concordance 0.99997 vs TASSEL** on an F2 (1 cell / 39,939 — a lone Viterbi tie).
        The site-set check caught+fixed a real bug: the SEED window's alleles must be written
        (`updatePopulationDataAlleles(...,0,window)`) before extension — was missing, left the
        first window NN. Frozen `biparental_fixture.rds` + `test-fsfhap-biparental-parity.R`
        (asserts site-set match + concordance; no-JVM CI gate). CodeRabbit-clean.
        NOTE: that parity ran LD OFF both sides (`-minR 0`), so it validated preFilterSites'
        non-LD steps + the finder on that fixture.
  - [x] stage 2b — **LD-on site-set parity: VALIDATED (bit-exact)**. `hapFinder.minR2 = minRforSnps`
        (the FSFHap `-minR` flag), so `-minR 0.2` turns the LD filter on. On an F2 with
        interspersed noise markers (random hom, uncorrelated → low windowed r², isolating the
        LD filter), `fsfhap_prefilter_sites_cpp(min_r2=0.2)` keeps **420/500, dropping exactly
        the 80 noise markers — identical to TASSEL** (420, 80/80). `agent/fsfhap_biparental_ld_parity.R`.
  - [x] **RESOLVED (task #5) — finder is bit-exact on RIL, its target population**. Root cause
        (found by instrumenting the forward loop): in a ~50%-het **F2**, het-heavy windows
        fragment the clusters so `mergeMajorHaplotypes` returns a single candidate → one parent
        list empties → the nearest-parent fallback then sends every later candidate to the
        surviving parent (permanent cascade) → only the seed window is written. This is an
        **out-of-scope stress case**: `BiparentalHaplotypeFinder` targets mostly-homozygous
        RIL/inbred full-sib families. On a **RIL** (inbred, ~0 het) the finder is robust and
        **BIT-EXACT vs TASSEL — site sets 500/500 (Jaccard 1.0) + per-cell concordance 1.00000**
        (`agent/fsfhap_ril_parity.R`, frozen `ril_fixture.rds`, `test-fsfhap-biparental-parity.R`).
        The seed=9 F2 fixture (0.99997) is kept as a stress-tolerance check.
  - [x] **Carry-forward guard added** (`fsfhap_biparental_alleles_cpp`, both extension loops;
        NON-TASSEL robustness extension): a degenerate window that would leave a parent list
        EMPTY is skipped, retaining the previous parents so neither can vanish and drive the
        one-sided cascade. **No-op when both parents are filled → RIL/BC1 parity UNCHANGED
        (verified: RIL site-set match + concordance 1.00000).** Mitigates the F2 cascade
        (seed=21 @420: 100→270 assigned; total collapse gone) but does not make F2 fully
        assign — F2 remains an out-of-scope stress case (het windows + post-skip overlap
        staleness). The guard only diverges from TASSEL on windows TASSEL would itself cascade
        on, so it cannot affect the validated RIL/BC1 outputs.
  - [ ] stage 2b — wire the **non-BC1 route** into the dispatcher (contribution ∉
        {0.75,0.25} → `BiparentalHaplotypeFinder`) + **F2 end-to-end parity**
        (`imputed_parents` per-cell diff, same gate that made BC1 bit-exact).
- [x] **Caller wiring + design-routing dispatcher** — `caller = "fsfhap"` in
      `call_states`/`call_ancestry` (`R/fsfhap.R` `.fsfhap_states`): pools **per
      family × chromosome** (stage 1 → stage 3), emits the common REF/HET/ALT segment
      schema (`source,donor,name,chr,start_bp,end_bp,state`). Dispatcher derives
      `phet=(1-F)/2` and the route from `design` (`BC{n}S{m}`); **BC1 only** —
      non-BC1 designs error (BiparentalHaplotypeFinder/multipleBC not ported).
      `family` via a `data` column or a vector arg. `test-fsfhap-caller.R` (14 checks:
      schema, phet derivation, routing guards, missing-column guards, family-vector).
      NOTE: `min_gametes=200` (TASSEL) is hardcoded in `.fsfhap_states` → chromosomes
      with <~100 kept sites/taxon drop all taxa; expose if real data needs it.
- [x] **Cross-family/chromosome parallelism** (`threads=`, closes the full-genome gap):
      `.fsfhap_states` builds one independent work unit per family × chromosome and
      fans them out via `parallel::mclapply` (unix; serial `lapply` elsewhere) — the
      per-unit stages 1+3 are untouched, so results are **bit-identical to serial**
      (`identical(serial, parallel)` verified). `threads` is plumbed from
      `call_states`/`call_ancestry`. Profiling the scaling also exposed a serial
      character `order()` hotspot (`.fsfhap_states` final sort + shared `to_segments`
      RLE); both switched to radix on integer factor codes (cf. `.lbimpute_prep`),
      which alone cut full-genome **single-threaded** 34 s → 10.9 s. Full-genome
      head-to-head now **4.2× faster than TASSEL** (3.85 s vs 16.1 s compute, 10 chr /
      10 threads) — see the Tier-0 table above. Suite green (incl. all parity gates +
      every `to_segments` consumer).
- [ ] ~~Design-routing dispatcher~~ (folded into the caller wiring above) — reads the design
      (pedigree `contribution1` / the `BC{n}S{m}` token) and selects the route +
      segregation ratio, mirroring TASSEL's dispatch:
      `contribution1 ∈ {0.75,0.25}` → BC route (ratio 0.25, stage 1a/1b); else
      `multipleBC` / `BiparentalHaplotypeFinder`. **This is where design→ratio
      derivation belongs** — the leaf routes (`.fsfhap_call_parents_bc`,
      `fsfhap_segregating_sites_cpp`) never default or validate `ratio`; the
      dispatcher guarantees a valid design-derived value. (Rejected CodeRabbit's
      leaf-level `Rcpp::stop` guard as wrong-layer: it treats a design parameter as
      untrusted input. `.fsfhap_segregating` now requires `ratio` — no default.)
- [x] **stage 3 — `imputeUsingViterbiFiveState`** (`fsfhap_impute_five_state_cpp`):
      5-state Viterbi-training EM (≤50 iters; re-estimate emission `counts/rowsum`
      no-pseudocount + transition `setTransitionCounts` row-normalize; distance-scaled
      per-node transition via the Haldane map fn `(1-exp(-2m))/2`; converge on stable
      emission-count matrix), init emission Table 1 + transition FILLIN Table 2, `pTrue`
      from design-derived `phet` ([.fsfhap_phet], `(1-F)/2`), output `0→AA,1-3→AC,4→CC`,
      earliest-index tie-break. **All prior spec gaps (setTransitionCounts, init matrices,
      recode frame) validated by the bit-exact e2e parity below.**
- [x] **stage 3 — `fillGapsInAlignment`** (`fsfhap_fill_gaps_cpp`): per-taxon
      forward-fill between EQUAL flanking non-N calls; `.fsfhap_impute` chains EM+fill.
      `test-fsfhap-impute.R` (18 checks: error-smoothing, phet derivation, gap-fill).
- [x] **stage 3 — END-TO-END TASSEL PARITY (bit-exact)**: full pipeline
      call→g → `.fsfhap_call_parents_bc` (stage 1) → `.fsfhap_impute` (stage 3) vs
      TASSEL's `imputed_parents` HapMap. **Per-cell concordance = 1.00000** over 23,429
      cells (391 markers × 60 taxa; our kept set == TASSEL's). `agent/fsfhap_e2e_parity.R`;
      frozen to `tests/fixtures/fsfhap_tassel/e2e_fixture.rds` + `test-fsfhap-e2e-parity.R`
      (no-JVM CI gate for stages 1+3 combined). This supersedes the stage-1 count-only gate.
- [~] ~~4-state cross-progeny / rephase~~ — **NOT the FSFHap path** (separate
      `ImputeProgenyStatesPlugin` lineage); dropped from the TeoNAM roadmap.
- [ ] `R/fsfhap.R` driver + `caller_spec("fsfhap")` + `family`/priors contract;
      5/4→3 output projection; `design_priors` wiring. **Core takes only the
      canonical normalized matrices — no file parsing.**
- [x] **IO/adapter layer in `io.R`** (minimal surface, separate from the algorithm):
      `read_hapmap()` (TASSEL HapMap — FSFHap's native input; single-char IUPAC or
      two-char diploid → canonical `(name,chr,pos,g)`), `read_plink()` (SNP-major
      `.bed`/`.bim`/`.fam` → `(name,chr,pos,g)`, g = A1 dosage), and `read_pedigree()`
      (FSFHap 7-col TSV / PLINK `.fam` → `taxon,family,parent1,parent2,contribution1/2,F`;
      phet via `.fsfhap_phet(F)`). Chosen as clearly-named sibling readers rather than a
      `format=` overload on the VCF-specific `read_vcf_gt()`. `test-io-adapters.R` (20
      checks: IUPAC + diploid HapMap, hand-built `.bed` round-trip, pedigree both formats,
      end-to-end `read_hapmap → call_ancestry(caller="fsfhap")`). Exported + documented.
- [x] `benchmarks/benchmark_fsfhap.R` — Tier-0 speed harness scaffolded: small→large
      ladder (`xs/s/m/full`), synthetic segregating BC family generator (+true
      breakpoints, reusable as the Tier-3 seed), TASSEL Hapmap/pedigree writer,
      JVM-startup baseline separated, peak-RSS (macOS), throughput + speedup table
      → `benchmarks/benchmark_fsfhap_out/`. TASSEL baseline runnable now
      (`Rscript benchmarks/benchmark_fsfhap.R xs 1`); port side self-skips until the
      caller exists. Override real data via `FSFHAP_HAPMAP`/`FSFHAP_PED`.
- [ ] Freeze fixtures: small family slice (CI) + 2–3 families (extended), skeleton
      HapMap + pedigree + TASSEL outputs, under `tests/fixtures/` + SHA256SUMS.
- [ ] Tier 1 — per-stage + end-to-end TASSEL parity on the small fixture (CI gate);
      extend to full family/full-genotype via an `agent/` `\dontrun` harness.
- [ ] Tier 2 — masked-recovery accuracy (port vs TASSEL on identical mask).
- [ ] Tier 3 — extend the simulator to emit a full-sib family with known breakpoints
      (new `sim_truth/` fixture); score donor-fragment Dice + breakpoint distance.

## Open questions
1. **Which parent-calling route to port first** — the **BC route is the default and
   matches our design**; port it first, add cluster/window routes only if TeoNAM
   parity needs them.
2. **Empirical NAM `transitionCounts`** — port TASSEL's baked matrix, re-derive from
   our B73v5 map, or replace with a parametric distance transition (`lbimpute`-style)?
3. **Iteration count** — the paper says "to convergence" but the code runs a single
   rephase pass; match TASSEL's actual behavior, not the paper's prose.
4. **5-state necessity for single-plant lines** — the residual-het classes target
   bulked DNA; confirm whether TeoNAM/nNIL lines are bulked before deciding if the
   3A:1C / 1A:3C classes are ever active.
5. **nNIL grouping QC** — run only confidently-assigned lines, or add a grouping-QC
   step given the ~1/3 pedigree mismatch?
6. **Worth it vs TASSEL?** Motivation = common-schema integration + speed/control (no
   JVM). **Tier-0 baseline (above) shows speed is NOT the driver** — TASSEL is
   ~18 s/full-family, linear, <1 GB. So justify the port on integration/dependency
   removal + parity, and confirm it clears the TASSEL parity bar (Tiers 1–3) before
   retiring the `teonam_fsfhap.R` TASSEL call.
