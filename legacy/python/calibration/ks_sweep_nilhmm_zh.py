#!/usr/bin/env python
"""Phase-2 calibration of the nilHMM count caller on Zh (skim -> BC2S2).

Calibrate the per-adjacent-marker recombination r (the segment-size driver) by KS
distance of the called donor-block size distribution (Mb) vs the BC2S2 simulation
(agent/bc2s2_segments.csv). Donor block = contiguous run of state>0 (HET union ALT),
matching the sim's donor-union segments. Blocks pooled across all Zh samples.

Guard: re-run B73 controls at the chosen r and confirm donor dosage stays ~0.

Usage: conda run -n nilhmm python agent/ks_sweep_nilhmm_zh.py
"""
import sys, time
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

sys.path.insert(0, "agent/nilhmm")
from nilhmm.io import read_vcf_counts
from nilhmm.core import introgression_hmm_counts

ZH = "data/nilhmm/Zh_counts.vcf.gz"
B73 = "data/nilhmm/B73_counts.vcf.gz"
REF_CSV = "agent/bc2s2_segments.csv"     # BC2S2 (g=4) donor-union sim segments
OUT = "results/sim_calibration/nilhmm_zh"

# fixed emission params for the r sweep; refined after
ERR, CONC = 0.01, 20.0
F1, F2 = 0.0625, 0.0938                   # BC2S2 single-locus freqs

def donor_block_sizes_mb(calls, chrom_of, pos):
    """Pool donor-block sizes (Mb) across samples. Boundaries at midpoints
    between flanking opposite-state markers (avoids 0-Mb singletons)."""
    sizes = []
    chroms = np.unique(chrom_of)
    for i in range(calls.shape[0]):
        for c in chroms:
            m = chrom_of == c
            st = calls[i, m] > 0
            p = pos[m]
            if not st.any():
                continue
            # boundaries between consecutive markers
            mid = (p[:-1] + p[1:]) / 2.0
            lo = np.concatenate(([p[0]], mid))      # left edge of each marker's cell
            hi = np.concatenate((mid, [p[-1]]))     # right edge
            # runs of donor
            d = np.diff(np.concatenate(([0], st.astype(int), [0])))
            starts = np.where(d == 1)[0]
            ends = np.where(d == -1)[0] - 1
            for s, e in zip(starts, ends):
                sizes.append((hi[e] - lo[s]) / 1e6)
    return np.array(sizes)

def run(ref, alt, mdict, chrom_of, pos, r):
    calls = introgression_hmm_counts(ref, alt, mdict, err=ERR, conc=CONC,
                                     r=r, f_1=F1, f_2=F2, return_calls=True)
    return calls

def main():
    import os; os.makedirs(OUT, exist_ok=True)
    import logging; logging.basicConfig(level=logging.ERROR)

    sim = pd.read_csv(REF_CSV)["mb"].values
    print(f"BC2S2 ref segments: {len(sim)} (median {np.median(sim):.2f} Mb, mean {sim.mean():.2f})")

    ref, alt, mdict, samples, minfo = read_vcf_counts(ZH)
    chrom_of = minfo["CHROM"].values
    pos = minfo["POS"].values.astype(float)
    print(f"Zh: {ref.shape[0]} samples x {ref.shape[1]} markers")

    # time one run
    t0 = time.time()
    c0 = run(ref, alt, mdict, chrom_of, pos, 0.01)
    print(f"single run: {time.time()-t0:.1f}s")

    rs = [1e-5, 3e-5, 7e-5, 1e-4, 2e-4, 3e-4, 5e-4, 1e-3]
    rows = []
    for r in rs:
        t0 = time.time()
        calls = run(ref, alt, mdict, chrom_of, pos, r)
        db = donor_block_sizes_mb(calls, chrom_of, pos)
        if len(db) == 0:
            print(f"r={r:<7g} no donor blocks"); continue
        D = ks_2samp(db, sim).statistic
        donor_frac = (calls > 0).mean()
        rows.append(dict(r=r, D=D, n_blocks=len(db), median_mb=np.median(db),
                         mean_mb=db.mean(), donor_frac=donor_frac,
                         blocks_per_sample=len(db)/ref.shape[0]))
        print(f"r={r:<7g} D={D:.4f} n_blocks={len(db):5d} median={np.median(db):6.2f}Mb "
              f"donor_frac={donor_frac:.4f} ({time.time()-t0:.1f}s)")

    out = pd.DataFrame(rows).sort_values("r")
    out.to_csv(f"{OUT}/zh_ks_sweep.csv", index=False)
    best = out.loc[out["D"].idxmin()]
    print(f"\nargmin D: r={best.r:g} (D={best.D:.4f}, median {best.median_mb:.2f} Mb vs sim {np.median(sim):.2f})")

    # B73 guard at best r
    rb, ab, mdb, sb, mib = read_vcf_counts(B73)
    cb = introgression_hmm_counts(rb, ab, mdb, err=ERR, conc=CONC, r=float(best.r),
                                  f_1=F1, f_2=F2, return_calls=True)
    het = (cb == 1).sum(1); altc = (cb == 2).sum(1); n = cb.shape[1]
    dos = (het + 2*altc) / (2*n)
    print(f"B73 guard @ r={best.r:g}: dosage mean={dos.mean():.4f} median={np.median(dos):.4f} max={dos.max():.4f}")

if __name__ == "__main__":
    main()
