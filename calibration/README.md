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
