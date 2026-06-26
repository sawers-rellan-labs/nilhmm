#!/usr/bin/env python
"""Per-taxon Phase-2 calibration of the nilHMM count caller (skim -> BC2S2).

For each teosinte taxon, sweep the per-adjacent-marker recombination r and pick the
argmin KS distance of the called donor-block size distribution (Mb) vs the BC2S2 sim
(agent/bc2s2_segments.csv). emission err/conc fixed at the Zh-calibrated values; freqs
= BC2S2 single-locus. Guard: re-run B73 controls at each taxon's argmin r.

Usage: conda run -n nilhmm python agent/ks_sweep_nilhmm_taxa.py [Zx Zv Zd Zl ...]
"""
import sys, time
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

sys.path.insert(0, "agent/nilhmm")
from nilhmm.io import read_vcf_counts
from nilhmm.core import introgression_hmm_counts

REF_CSV = "agent/bc2s2_segments.csv"
OUT = "results/sim_calibration/nilhmm_taxa"
ERR, CONC, F1, F2 = 0.01, 20.0, 0.0625, 0.0938
RS = [1e-5, 2e-5, 3e-5, 5e-5, 7e-5, 1e-4, 2e-4, 3e-4, 5e-4, 1e-3]

def donor_blocks_mb(calls, chrom_of, pos):
    sizes = []
    for i in range(calls.shape[0]):
        for c in np.unique(chrom_of):
            m = chrom_of == c
            st = calls[i, m] > 0
            p = pos[m]
            if not st.any():
                continue
            mid = (p[:-1] + p[1:]) / 2.0
            lo = np.concatenate(([p[0]], mid)); hi = np.concatenate((mid, [p[-1]]))
            d = np.diff(np.concatenate(([0], st.astype(int), [0])))
            s = np.where(d == 1)[0]; e = np.where(d == -1)[0] - 1
            sizes.extend((hi[e] - lo[s]) / 1e6)
    return np.array(sizes)

def main():
    import os, logging
    os.makedirs(OUT, exist_ok=True)
    logging.basicConfig(level=logging.ERROR)
    taxa = sys.argv[1:] or ["Zx", "Zv", "Zd", "Zl"]

    sim = pd.read_csv(REF_CSV)["mb"].values
    print(f"BC2S2 ref: {len(sim)} segs, median {np.median(sim):.2f} Mb", flush=True)

    # B73 control once
    rb, ab, mdb, sb, mib = read_vcf_counts("data/nilhmm/B73_counts.vcf.gz")

    summary = []
    for tx in taxa:
        ref, alt, md, samp, mi = read_vcf_counts(f"data/nilhmm/{tx}_counts.vcf.gz")
        chrom = mi["CHROM"].values; pos = mi["POS"].values.astype(float)
        print(f"\n=== {tx}: {ref.shape[0]} samples x {ref.shape[1]} markers ===", flush=True)
        rows = []
        for r in RS:
            t0 = time.time()
            calls = introgression_hmm_counts(ref, alt, md, err=ERR, conc=CONC,
                                             r=r, f_1=F1, f_2=F2)
            db = donor_blocks_mb(calls, chrom, pos)
            D = ks_2samp(db, sim).statistic
            rows.append(dict(taxon=tx, r=r, D=D, n_blocks=len(db),
                             median_mb=np.median(db), donor_frac=(calls > 0).mean()))
            print(f"  r={r:<7g} D={D:.4f} n={len(db):5d} median={np.median(db):6.2f}Mb "
                  f"frac={(calls>0).mean():.4f} ({time.time()-t0:.1f}s)", flush=True)
        df = pd.DataFrame(rows).sort_values("r")
        df.to_csv(f"{OUT}/{tx}_ks_sweep.csv", index=False)
        best = df.loc[df["D"].idxmin()]
        # B73 guard at this taxon's argmin r
        cb = introgression_hmm_counts(rb, ab, mdb, err=ERR, conc=CONC,
                                      r=float(best.r), f_1=F1, f_2=F2)
        bdos = ((cb == 1).sum(1) + 2*(cb == 2).sum(1)) / (2*cb.shape[1])
        edge = "EDGE!" if best.r in (RS[0], RS[-1]) else ""
        print(f"  -> argmin r={best.r:g} D={best.D:.4f} median={best.median_mb:.2f}Mb "
              f"B73max={bdos.max():.4f} {edge}", flush=True)
        summary.append(dict(taxon=tx, n_samples=ref.shape[0], best_r=best.r, D=best.D,
                            median_mb=best.median_mb, donor_frac=best.donor_frac,
                            b73_max_dosage=bdos.max(), at_grid_edge=bool(edge)))

    s = pd.DataFrame(summary)
    s.to_csv(f"{OUT}/summary.csv", index=False)
    print("\n=== SUMMARY ===", flush=True)
    print(s.to_string(index=False), flush=True)

if __name__ == "__main__":
    main()
