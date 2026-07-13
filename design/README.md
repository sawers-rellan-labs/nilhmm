# Design & development notes

Internal design documents for the nilHMM R package. These are **not** part of the
built package (the `design/` folder is `.Rbuildignore`d); they record the reasoning
behind the implementation. User-facing documentation is the top-level `README.md`,
the roxygen help (`?call_ancestry`), and the pkgdown site.

- **architecture.md** — the 3-state HMM math (hidden states, emission/transition
  model) with the architecture diagram.
- **REFACTOR_R_PACKAGE.md** — the Python → R + Rcpp refactor plan (engine =
  one HMM × swappable emission × duration; the S-numbered design sections).
- **RTIGER_PORT.md** — porting the RTIGER rigidity HMM off Julia; kernel-by-kernel
  checklist and bit-for-bit validation notes.
- **VALIDATION.md** — full-cohort validation of the R callers against the frozen
  Python baselines (concordance, segment-size KS).
- **Implementation.md** — the count-path (BetaBinomial) design, calibration, and
  results.
- **BRB_run_findings.md** — the 3′ RNA-seq (BrB) run: reference-bias / ALT-collapse
  findings that motivated EM-fit emission means.
- **PEDIGREE_HMM.md** — pedigree-aware ancestry *refinement* (`refine_ancestry`):
  structured loopy BP over the pedigree × genome grid that couples relatives to
  correct per-individual `call_states()` calls; design decisions, C++ BP kernel,
  and the R/C++ boundary.
- **package_structure.md** — package organization and design decisions.
- **Zv_RTIGER_divergence.md** — the Zv (parviglumis) low-divergence identifiability
  investigation.
- **jim_hmm_diagram.Rmd**, **jim_hmm.png** — the original HMM architecture diagram.
