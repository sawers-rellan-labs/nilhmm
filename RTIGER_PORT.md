# Faithful port of RTIGER → R/Rcpp (in progress)

Reproducing RTIGER by **direct port of the user's fork** `faustovrz/RTIGER`
(commit `649cbf6`, branch `optimize-julia-core`, = installed `RTIGER` 2.1.0), file
`julia/rHMM_methods.jl` (1708 lines). This REPLACES the earlier invented
`duration_rigidity`/`p_switch` rtiger caller (that was a fabrication, not RTIGER).

**Method:** port the **serial path** (fork marks it bit-identical to upstream);
verify every kernel against the real Julia via `agent/rtref.jl`
(`RT=<pkg> julia rtref.jl`). Julia has Distributions + Optim. Reproduction bar is
distributional (RTIGER's init is `runif`-randomized → not bit-identical end-to-end).

## Data convention (generateObject)
Input file cols: chr,pos,P1allele,**k=P1.Allele.Count**,P2allele,P2count; **n = k + P2count**.
For our skim TSVs (chr,pos,ref,nref,alt,nalt): **k = nref, n = nref+nalt** → ref-allele fraction.
States ordered **[pat (ref≈0.95), het (0.5), mat/donor (ref≈0.05)]** — note this is the
REVERSE of nilHMM's [REF,HET,ALT]; remap on output to common schema.

## Initialization (generate_params, randomize=TRUE)
- transition `A`: `0.1*ones(s,s) + 10*I` (+ runif(0,.01)), row-normalized → diag≈0.98.
- emission `α=[20,20,1]*ps`, `β=[1,20,20]*ps`, `ps≈1±0.5`.
- `pi`: `1 + runif(-1,3)` per state, normalized.

## Algorithm (per chain = sample×chr; pooled M-step over all chains)
E-step: `getlogpsi`(BetaBinomial on k,n) → `productpsi`(windowed log-sum, s×(T+r))
→ `forward`/`backward` (rigidity: only diagonal "stay" each step + r-step "enter/leave")
→ `zeta` (T×s×s, normalized by max) → `gamma`.
M-step: transition `A[i,j]=Σ_band zeta / rowsum`; start `=mean gamma[:,1]`;
emission per state: `m=Σγk/Σγn` clamp[0.01,0.99], `tau`=Brent argmax Σ w·logpdf(BB(n,tau·m,tau·(1-m)),k)
(init tau=α/m+β/(1-m), cap 100), `a=tau·m, b=tau·(1-m)`.
Converge: `max(|Δα|,|Δβ|) ≤ eps` (R default eps=0.01, max.iter=50).
Decode: rigidity `viterbi` (max-product; backtrace fills r-1 on a switch) on ALL samples.
Post: `postprocessing` shifts segment borders via `findsplit` (cumsum of emission, argmax).

## Port status (src/rtiger.cpp, R/ driver TBD)
- [x] `rtiger_getlogpsi_cpp`  — verified vs Julia (max|Δ| 4e-9)
- [x] `rtiger_productpsi_cpp` — verified vs Julia (max|Δ| 5e-9)
- [x] `rtiger_forward_cpp` / `rtiger_backward_cpp` (L95-196) — verified vs Julia (max|Δ| 5e-9)
- [x] `rtiger_zeta_cpp` / `rtiger_gamma_cpp` (L218-322) — verified vs Julia (max|Δ| 5e-9); E-step complete
- [x] `rtiger_viterbi_cpp` (L585-646) — verified vs Julia (exact path incl. switches)
- [ ] emission M-step `emissionUpdateState` (Brent; L459-517) — use R optimize/C++ Brent
- [ ] transition/start M-step (L353-410) ; EM loop + convergence (fit, L1412-1582)
- [ ] R driver: build obs (k,n) per sample×chr, init params, run EM, Viterbi, postproc,
      remap states→common schema. New caller `caller="rtiger"` (replace invented path).
- [ ] validate vs frozen `rtiger_rigidity_ref/calls_taxa_r5.csv` (distributional + concordance)

Reference harness: `agent/rtref.jl` (extend to dump forward/backward/zeta/gamma/viterbi
for the same tiny fixed input to verify each kernel as ported).
