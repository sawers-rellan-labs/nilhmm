#!/usr/bin/env python
"""Phase 3 — stage nilHMM count-caller calls into the common segment schema.

For each teosinte taxon, run the calibrated count caller (per-taxon r from
calibration/calibrated_params_all_taxa.csv; err=0.01, conc=20, BC2S2 freqs),
run-length-encode per sample x chromosome into segments, and emit the common schema
shared with Skim-RTI / Skim-BIN:

    source, donor, name, chr(int), start_bp, end_bp, state(REF=0/HET=1/ALT=2)

Segment boundaries = first/last marker position of each contiguous state run (same
convention as the RTIGER segments). Outputs:
    results/sim_calibration/nilhmm_calls/calls_common_schema.csv         (5 taxa)
    results/sim_calibration/nilhmm_calls/calls_checks_common_schema.csv  (B73, Purple)

Usage: conda run -n nilhmm python agent/stage_nilhmm_calls.py
"""
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "agent/nilhmm")
from nilhmm.io import read_vcf_counts
from nilhmm.core import introgression_hmm_counts

OUT = "results/sim_calibration/nilhmm_calls"
SOURCE = "nilHMM_SNP50K_counts"
ERR, CONC, F1, F2 = 0.01, 20.0, 0.0625, 0.0938
PARAMS = "agent/nilhmm/calibration/calibrated_params_all_taxa.csv"
CHECK_R = 3e-5   # checks have no introgression structure; B73 is ~flat in r (dosage ~0.03 max)

def rle_segments(calls, chrom_of, pos, name_of, donor):
    """Run-length-encode per sample x chrom -> common-schema rows."""
    rows = []
    chroms = np.unique(chrom_of)
    for i in range(calls.shape[0]):
        nm = name_of[i]
        for c in chroms:
            m = chrom_of == c
            s = calls[i, m]
            p = pos[m]
            if s.size == 0:
                continue
            # boundaries where state changes
            brk = np.where(np.diff(s) != 0)[0]
            starts = np.concatenate(([0], brk + 1))
            ends = np.concatenate((brk, [s.size - 1]))
            for a, b in zip(starts, ends):
                rows.append((SOURCE, donor, nm, int(c), int(p[a]), int(p[b]), int(s[a])))
    return rows

def call_taxon(tx, r):
    ref, alt, md, samples, mi = read_vcf_counts(f"data/nilhmm/{tx}_counts.vcf.gz")
    calls = introgression_hmm_counts(ref, alt, md, err=ERR, conc=CONC, r=r,
                                     f_1=F1, f_2=F2)
    return rle_segments(calls, mi["CHROM"].values, mi["POS"].values.astype(float),
                        np.array(samples), donor=tx)

def main():
    import os, logging
    os.makedirs(OUT, exist_ok=True)
    logging.basicConfig(level=logging.ERROR)
    cols = ["source", "donor", "name", "chr", "start_bp", "end_bp", "state"]

    par = pd.read_csv(PARAMS).set_index("taxon")["r"].to_dict()
    taxa = ["Zh", "Zx", "Zv", "Zd", "Zl"]

    rows = []
    for tx in taxa:
        r = float(par[tx])
        trows = call_taxon(tx, r)
        rows += trows
        print(f"{tx}: r={r:g}  {len({x[2] for x in trows})} samples, {len(trows)} segments", flush=True)
    df = pd.DataFrame(rows, columns=cols).sort_values(["donor", "name", "chr", "start_bp"])
    df.to_csv(f"{OUT}/calls_common_schema.csv", index=False)
    st = df["state"].value_counts().sort_index().to_dict()
    print(f"\nTAXA -> {OUT}/calls_common_schema.csv : "
          f"{df['name'].nunique()} samples, {len(df)} segments, states {st}", flush=True)

    # checks (B73, Purple) at a nominal r
    crows = []
    for grp in ["B73", "Purple"]:
        ref, alt, md, samples, mi = read_vcf_counts(f"data/nilhmm/{grp}_counts.vcf.gz")
        calls = introgression_hmm_counts(ref, alt, md, err=ERR, conc=CONC, r=CHECK_R,
                                         f_1=F1, f_2=F2)
        g = rle_segments(calls, mi["CHROM"].values, mi["POS"].values.astype(float),
                         np.array(samples), donor=grp)
        crows += g
        print(f"{grp}: r={CHECK_R:g}  {len(samples)} samples, {len(g)} segments", flush=True)
    cdf = pd.DataFrame(crows, columns=cols).sort_values(["donor", "name", "chr", "start_bp"])
    cdf.to_csv(f"{OUT}/calls_checks_common_schema.csv", index=False)
    print(f"CHECKS -> {OUT}/calls_checks_common_schema.csv : "
          f"{cdf['name'].nunique()} samples, {len(cdf)} segments", flush=True)

if __name__ == "__main__":
    main()
