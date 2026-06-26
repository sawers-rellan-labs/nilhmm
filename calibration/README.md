# Calibration scripts (BZea SNP50K skim)

Caller-side scripts that calibrate `nilhmm`'s count caller against the BC2S2
simulation, plus the resulting parameters. See `../Implementation.md` for the full
write-up (data problem, sweeps, tables, results).

**These are snapshots from the `zealtiger` analysis project and run from its root**,
not from this package — they reference `zealtiger` paths (`data/nilhmm/…`,
`agent/bc2s2_segments.csv`, `results/…`). They are kept here so the nilHMM
calibration story is self-contained and reproducible; they are not part of the
installable package.

| file | purpose |
|---|---|
| `make_nilhmm_taxon_subset.sh` | build a per-taxon (or B73/Purple check) counts subset (GT+AD) from `cohort.vcf.gz` |
| `ks_sweep_nilhmm_zh.py` | r sweep on Zh, KS vs BC2S2 sim donor-block sizes; B73 guard |
| `err_sweep_zh.py` | emission `err` sweep at the best r (segment-size vs HET:ALT tradeoff) |
| `ks_sweep_nilhmm_taxa.py` | per-taxon r sweeps (Zx/Zv/Zd/Zl); writes per-taxon CSVs + summary |
| `Zh_calibrated_params.csv` | calibrated Zh parameters (r=3e-5, err=0.01, conc=20, BC2S2 freqs; KS D=0.043) |

Run (from the zealtiger project root, in the `nilhmm` conda env):

```bash
bash agent/make_nilhmm_taxon_subset.sh Zh        # and Zx Zv Zd Zl, B73, Purple
conda run -n nilhmm python agent/ks_sweep_nilhmm_zh.py
conda run -n nilhmm python agent/ks_sweep_nilhmm_taxa.py Zx Zv Zd Zl
```

Reference distribution: `agent/bc2s2_segments.csv` (n=1500 simcross BC2S2 donor-union
segments; skim = bulked BC2S2). Calibration target = KS distance of called donor-block
sizes (Mb) vs that simulation.

## Per-taxon calibration results

Calibrated per teosinte taxon (KS-on-donor-block-size vs the BC2S2 sim, median 12.41 Mb;
err=0.01, conc=20, f_1=0.0625, f_2=0.0938 fixed). All minima are interior (not at a grid
edge); the B73 control stays clean at every taxon's optimum. Full sweeps in `results/`;
parameters in `calibrated_params_all_taxa.csv`.

| taxon | species | n | best r | KS D | median block (Mb) | donor frac | B73 max dosage |
|---|---|---|---|---|---|---|---|
| Zh | ssp. *huehuetenangensis* | 61 | 3e-5 | 0.043 | 12.89 | 0.086 | 0.032 |
| Zx | ssp. *mexicana* | 580 | 7e-6 | 0.062 | 12.63 | 0.079 | 0.032 |
| Zv | ssp. *parviglumis* | 386 | 3e-6 | 0.072 | 12.22 | 0.065 | 0.032 |
| Zd | *Z. diploperennis* | 255 | 7e-6 | 0.060 | 12.26 | 0.055 | 0.032 |
| Zl | *Z. luxurians* | 121 | 3e-6 | 0.089 | 11.18 | 0.054 | 0.032 |

**Findings:**
- **`r` is taxon-dependent** (3e-6 … 3e-5, ~10× spread) — per-taxon fitting was warranted,
  not a single global `r`. All optima are far below the physical per-marker recombination
  rate (~3e-4), i.e. the calls need strong smoothing to match BC2S2 block sizes (consistent
  with low-coverage noise and the "callers run longer than the model" pattern).
- **Donor fraction tracks the species** (Zh 8.6% → Zl 5.4%), a real per-taxon difference in
  introgression content, not a calibration artifact (donor_frac is ~r-stable within a taxon).
- **B73 control clean everywhere** (max dosage 0.032), so each taxon's `r` respects the
  negative-control constraint.
- Same depth-1 caveat as Zh applies: the donor **footprint** and **block-size distribution**
  are the reliable outputs; the het/hom split within blocks is not identifiable at ~1× and is
  strongly het-biased (see `../Implementation.md` §6.3).
