# Terminology & naming conventions (nilHMM package)

Coding and naming conventions for the nilHMM caller engine — the vocabulary a
developer or downstream analyst needs to name objects, functions, and states
consistently. The governing idea is a hard wall between **observed genotypes**
(inputs) and **inferred ancestry** (our output).

> **Scope.** This file is nilHMM's *coding/naming* terminology. Paper-facing and
> general/biological terms for the ZEAL/BZea work — phenotype definitions, taxa,
> datasets, population-structure covariates — live in `zealhmm/TERMINOLOGY.md`, not
> here.

## The two layers: genotype vs. ancestry mosaic

| | **Genotype** | **Ancestry mosaic** |
|-|--------------|---------------------|
| Role | **input / evidence** | **our inference / output** |
| Source | actual observations, or calls from a prior caller | inferred by us from the evidence |
| Method | per-site, **no linkage / no HMM** | HMM **across** sites (recombination + design prior) |
| Verb / fn | `call_gt()` | `call_ancestry()` / `call_states()` |
| Output suffix | **`_gt`** | **`_mosaic`** |
| States 0/1/2 | genotype dosage: 0 = REF, 1 = HET, 2 = ALT (donor) | ancestry: 0 = recurrent, 1 = het, 2 = donor |

**Why the wall matters:** a mosaic *overrides* the genotype — a donor-ancestry block
reports state `2` even at an invariant site (allele identical to the recurrent
parent), and it collapses observed hets into an ancestry state. A mosaic is **not** a
genotype and must never be filed, labeled, or served as one. "Imputed genotype" is
also banned: it is ambiguous and blurs the wall. **Never call the output
"genotypes"** — genotypes are inputs; the mosaic is what we infer.

## Callers

A **caller** is a method = **(emission × duration + priors)** over the shared 3-state
(REF/HET/ALT) HMM engine (see `architecture.md`, `REFACTOR_R_PACKAGE.md`). The
`caller` names:

| `caller` | engine |
|----------|--------|
| `nnil` | categorical GT emission × geometric duration (Holland nNIL, hard calls) |
| `bbnil` | count / BetaBinomial emission × geometric (low-coverage count extension of nNIL) |
| `catiger` | categorical GT emission × rigidity |
| `rtiger` | count / BetaBinomial × rigidity (Julia-free RTIGER port) |
| `binhmm` | 3-state Gaussian on ~1 Mb bins |
| `googa` | categorical GT × geometric, GOOGA competitive-alignment thresholds (faithful GOOGA/Veltsos) |
| `atlas` | categorical GT × rigidity, GOOGA thresholds (this work's transcript caller) |
| `lbimpute` | LB-Impute port (external baseline) |
| `fsfhap` | FSFHap port (external baseline; per-family) |
| `pedigree` | family-coupled belief propagation over the pedigree × genome grid |

The four grid callers are the pure (emission × duration) cells: `nnil` (gt +
geometric), `bbnil` (count + geometric), `catiger` (gt + rigidity), `rtiger`
(count + rigidity).

### Casing: the lowercase caller vs. its namesake

Caller identifiers are **always lowercase** — in code, and set in **fixed-width font**
in the paper (`\texttt{}`): `nnil`, `bbnil`, `catiger`, `rtiger`, `googa`, `atlas`,
`binhmm`, `lbimpute`, `fsfhap`, `pedigree` (and the genotype caller `call_gt`).
Lowercase always means **the method as implemented in this package** (this engine's
emission × duration cell), never anything external.

A capitalized / display form refers to something **different** — the external namesake
or the population, not our caller:

- `nnil` = our caller (the ported Holland nested-NIL HMM). **nNIL** = the *population*
  (Holland's nested near-isogenic lines). Never write `nNIL` for the caller.
- `rtiger` = the caller as implemented here (the Julia-free port). **RTIGER** = the
  *original algorithm and its Julia implementation* (Rowan et al.).
- The same rule extends to any caller with an upstream namesake, e.g. `googa` (our port)
  vs. GOOGA (the original Flagel/Veltsos method). Callers introduced by this work
  (`bbnil`, `catiger`, `atlas`) have no separate namesake and are only ever lowercase.

So `nnil`/`rtiger`/`googa` name what runs in this package; `nNIL` names a population and
`RTIGER`/GOOGA name upstream algorithms.

`call_gt()` (flat prior = maximum likelihood; HWE prior = MAP / het-excess) is a
**genotype** caller, **not** an ancestry caller — it is not in this table and not
dispatchable through `call_ancestry()`. The genotype-vs-mosaic wall below is why:
a caller produces an ancestry *mosaic*; `call_gt()` produces *genotypes*.

- **`mosaic` is a noun** (the ancestry-state matrix), **never a caller name.** A caller
  is the *method*; the mosaic is its *output*.
- Emission is the swappable observation channel: **count** (`emission_count`,
  BetaBinomial over read depths) vs. **gt** (`emission_gt`, categorical
  genotype-confusion over hard calls). Depth selects it (saturated → gt; intermediate
  → count).

## GL / GP / MAP

- **GL / PL** — genotype *likelihood* `P(reads | G)` (phred-scaled = PL). Argmax GL = ML.
- **GP** — genotype *posterior* `P(G | reads)` (VCF ≥ 4.3 FORMAT field).
- **MAP** — reserved for the *estimator* only: argmax over GP. Not a token, not a
  filename. Kept distinct from `map` (the genetic/recombination map) and `Map()` (the
  function) — hence GP, not MAP, in `_gt` tokens. `bcftools call -m` (and
  `call_gt(prior = "hwe")`) apply an HWE prior → the reported `GT` is the MAP
  ("HWE-posterior").

## Verbs vs. nouns (quick reference)

- **caller** = a method (engine + params). **genotype** = input. **mosaic** = output.
- `call_gt()` — verb, genotype layer → a per-site **genotype** (`_gt`).
- `call_ancestry()` = `call_states()` |> `to_segments()` — verb, mosaic layer → an
  **ancestry mosaic** (`_mosaic`).