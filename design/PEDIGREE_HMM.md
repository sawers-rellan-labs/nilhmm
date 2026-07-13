# Pedigree-Aware Ancestry Inference for BC2S3 Selfing Families
## Full Implementation Specification (C++, nilHMM-compatible)

## 0. Scope and design decisions (locked)

- **Population**: BC2S3 selfing families. Donor (`Zd.xxxx`) and recurrent
  backcross parents are pure lines and non-segregating.
- **Data**: joint-called, teosinte-informative SNPs only (~20k), low depth
  (mean λ ≈ 2.8×). Depth modeling happens **upstream** in `call_states()` (the
  count/BetaBinomial emission over donor-allele count vs. depth). This refine
  stage consumes only its output — the per-marker **hard state calls**, no read
  depth (see Role below).
- **Intermediate generations (P2/P3/P4) are NOT genotyped.** They are latent
  nodes with emission ≡ 1; their beliefs are driven entirely by messages.
  They are retained because they impose chromosome continuity on the consensus.
- **Method**: structured loopy belief propagation over the pedigree × genome
  grid. Each individual's chromosome chain is solved exactly (forward–backward);
  pedigree edges are the loopy part handled by damped message passing.
  This is NOT exact inference — the grid is loopy (4-cycles between adjacent
  markers across a parent–child edge). Report it as exact-within-chain-FB,
  converged-message-passing, linear-cost approximation. Note the qualifier
  "within-chain FB": the per-site cavity (§9) removes only the site-`i` value of
  the excluded message, not that message's influence on `belief[i]` propagated
  through the chain from neighboring sites, so the cross-chain step is
  approximate even though each chain solve is exact (see §9 note).
- **Role in nilHMM**: this is a **refinement stage** on top of the per-individual
  engine. It consumes a family's **hard per-marker state calls** (the
  `call_states()` mosaic — one 0/1/2 call per leaf per marker, **no depth**),
  builds a depth-blind genotyping-error emission from them via `emission_gt()`,
  couples the leaves through the pedigree, and returns refined per-marker states.
  It does not re-read counts and does not replace the engine's emission/duration
  axes — it corrects per-individual calls using relatives.
- **Language / boundary**: the BP kernel is **pure, data-agnostic C++**. It
  receives a fully parsed pedigree (edges, `meioses`, `hasData`), the per-leaf
  emission matrices, per-node marker-0 priors, and the per-interval `r` vector —
  all constructed in R. **R owns all IO**: `call_states()` produces the input
  mosaic, `emission_gt()` turns its hard states into emissions, `design_priors()`
  the priors, `read_pedigree()` the forest, `to_segments()` the output schema.
  The user-facing verb is `refine_ancestry(mosaic, pedigree = "pop.ped", ...)`
  where `mosaic` is the `call_states()` table (§17). No paths, sample-name
  grammars, or design literals live in the C++; project-specific name→pedigree
  conversion lives in `zealhmm` (see §12, §17).

## 1. State space

Per-individual **donor dosage** x ∈ {0, 1, 2}:

    0 = RP/RP  (homozygous recurrent)
    1 = RP/D   (heterozygous)          <- the only segregating state
    2 = D/D    (homozygous donor)

Rationale: het (1) vs. fixed-donor (2) is exactly what selfing segregation
acts on and what carries the pedigree information. The phased-haplotype model
({0,1}^2 + copy-switch latent) is the "correct" generative model but forces a
non-separable transmission factor; deferred to v2.

## 2. Selfing transmission kernel (per locus)

Row = parent dosage, column = child dosage:

                 child=0   child=1   child=2
    parent=0     1.00      0.00      0.00
    parent=1     0.25      0.50      0.25     <- only segregating row
    parent=2     0.00      0.00      1.00

Two rows deterministic: a homozygous locus stays homozygous in all descendants
forever. Heterozygosity halves each selfing generation.

**Approximation stated for methods**: using this as a per-marker factor drops
*transmission linkage* (a het parent transmits a donor segment as a coherent
block, not marker-by-marker). Each child's own chain FB re-imposes block
structure within that child, but the *cross-sib* block coherence — the mechanism
by which many genotyped siblings sharing an IBD block pool their depth over that
block — is exactly what this factorization discards. What survives is
marker-by-marker marginal pooling (many sibs het at marker `i` raises the
parent's het belief there, which feeds back). That is real but weaker than
block-coherent sharing, so temper the "many genotyped sibs" expected-gain claim
(§16) accordingly and report the sib-count gain explicitly in validation. If
breakpoint localization or the sib-count gain plateaus above the per-individual
baseline, dropped linkage is the first suspect and the phased copy-switch kernel
(§18) is the fix.

## 3. Chromosome transition (the ONLY change to the HMM engine)

The spatial (along-chromosome) transition must do two things: (a) create blocks —
adjacent markers usually share ancestry — and (b) relax toward a stationary that
does **not** re-inject a genotype-frequency prior already carried by the
transmission messages. A relax-to-stationary transition does both:

    A[x][y] = (1 - r_eff) * δ_xy  +  r_eff * π[y]

- **Stationary = `π` exactly** (`π A = π`) for any `π`. **Which `π`?** — see the
  correction below. The founder uses `π = π_0`; **non-root nodes use `π = uniform`**
  (marginal-neutral smoothing).
- `r_eff` is the **meiosis-compounded** recombination fraction (block-length
  control), the probability of an odd number of crossovers over `meioses`
  meioses:

      r_eff = (1 - (1 - 2*r_interval)^meioses) / 2

  which saturates correctly at the 0.5 ceiling. (Naive `1-(1-r)^m` — P(≥1
  crossover) — agrees to first order for tight SNPs but wrongly saturates at 1.0
  and overstates mixing across large gaps.)

**Why not `lump(B ⊗ B)`** (the earlier design, kept here for the record): two
independent per-haplotype chains sharing `B` force a *panmictic* stationary
`((1-φ)², 2φ(1-φ), φ²)` — Hardy–Weinberg at `φ=1/8`, i.e. inbreeding `F=0` —
*regardless of generation*. For BC2S3 that over-weights het ~7× at every interior
marker (0.219 vs. ~0.031). Independent haplotypes structurally cannot represent
selfing-induced homozygosity; that needs correlated haplotypes (the copy-switch
kernel, §18). So `B ⊗ B` is replaced by relax-to-`π_t`.

**CORRECTION (from the implementation + §16 validation).** An earlier draft set
the stationary to the node's *own* generation marginal `π_t` at **every** node, to
"cure heterozygote inflation." That **double-counts**: applying `Tsib` to `π_0`
already yields `π_1`, `π_2`, ... , so the transmission messages propagate the
generation marginal downward on their own. Injecting `π_t` *again* as each node's
chain stationary multiplies the prior against the transmission message and
**over-suppresses conditionally-correct heterozygosity** (a selfing child of a
HET parent is genuinely 50% het at that locus; `π_t` wrongly pulls it toward the
3% population marginal). The smoke test confirmed the skew, and the fix restored
the correct `(0.25,0.5,0.25)` behavior. So:

- **Founder (root)**: `rho = π_0`, chain stationary `π = π_0` (its only prior; §4).
- **Non-root nodes**: `rho = uniform`, chain stationary `π = uniform`
  (marginal-neutral). Genotype frequencies flow in through `Tsib` transmission,
  not the chain. This is what [refine_ancestry()] implements.

Heterozygote inflation is then handled correctly by the founder prior +
transmission (not a per-node stationary), and is further checked by sibling
consensus. The §16 run bears this out: refine *reduces* the HET false-call rate.

**Tradeoff** (state in methods): relax-to-`π` permits a direct 0→2 step at rate
`r_eff·π(2)`, i.e. it does **not** suppress double recombination the way
`B ⊗ B` did (which required both haplotypes to flip). At the tiny per-interval
`r_eff` of 20k SNPs this is negligible; the phased kernel (§18) recovers exact
double-recomb suppression **and** the inbred stationary together.

## 4. Priors and roots

All of the quantities below are **derived in R and passed in** — no literals in
the C++. This keeps the design assumptions in one auditable place
(`design_priors()`) instead of scattered as magic numbers.

- Data is a **forest**: one tree per founder. Grouping into trees and the
  parent–child edges are produced by the R name-parser (§12) and handed to the
  kernel as a pedigree; the kernel never inspects names. Families are fully
  separable → process (and parallelize) independently.
- **Root** of each tree = the founder (NOT the donor).
- **Founder prior `π_0`** (the only genotype-frequency prior in the model; §3
  correction). With donor *genome* fraction `q = f_2 + f_1/2` derived from
  `design_priors()` (BC2S3 → `q = 0.1094 + 0.0312/2 = 0.125`; selfing preserves
  the mean), the founder is heterozygous for all its donor content (at BC2 one
  whole haplotype is fixed RP), so founder het `H_0 = 2q` and

      π_0 = { 1 - 2q ,   2q ,   0 }        # (P0, P1, P2); BC2 → {0.75, 0.25, 0}

  Do **not** hard-code `{0.875, 0.125, 0}` (that conflates genome fraction with
  het). The deeper `π_t` marginals (`≈{0.859,0.031,0.109}` for a BC2S3 leaf) are
  **emergent**, not injected: applying `Tsib` to `π_0` gives `π_1, π_2, ...`, so
  the transmission messages carry them downward automatically (§3 correction).
- **Per-node stationary** `nd.pi`: `π_0` for the founder, **uniform** for every
  non-root node (marginal-neutral; §3). **`meioses`** per node = meioses
  accumulated to the founder (backcross generations + the F1 cross; BC2 → 3),
  +1 per descending selfing edge. Both derived in R (`.founder_prior()` /
  depth-from-root in [refine_ancestry()]), passed in per node.
- **Marker-0 prior** `rho`: founder uses `rho_root = π_0 = {0.75, 0.25, 0}`;
  non-root nodes use `rho = uniform` (the inherited prior arrives through the
  marker-0 downward message, so an informative `rho` there would double-count it —
  the same reason non-root chains are marginal-neutral, §3). `rho` touches only
  `alpha[0]`.

## 5. Data types

```cpp
constexpr int    K     = 3;
constexpr double FLOOR = 1e-12;      // cavity-division stability
using Col = std::array<double, K>;
using Mat = std::vector<Col>;        // M markers x K

static const double Tsib[K][K] = {
    {1.00, 0.00, 0.00},
    {0.25, 0.50, 0.25},
    {0.00, 0.00, 1.00}
};

inline void normalize(Col& c) {
    double s = c[0] + c[1] + c[2];
    if (s <= 0) { c = {1.0/3, 1.0/3, 1.0/3}; return; }  // see note below
    for (double& v : c) v /= s;
}
```

The `s <= 0` branch is a **deliberate** guard, not a silent bug-hider: an
all-floored column (every entry driven to `FLOOR` by cavity division against a
deterministic `Tsib` row, §9/§14) is a legitimate transient state, and falling
back to uniform is the correct recovery — erroring there would abort valid loopy
BP mid-sweep. It is *not* a place to absorb corrupt data: non-finite/negative
inputs are rejected at the R boundary (§17) before the kernel runs, so `s` here
is only ever a small non-negative sum, never `NaN`.

## 6. Transition construction

```cpp
// Meiosis-compounded recombination fraction: P(odd # crossovers) over `meioses`
// meioses, capped at the 0.5 ceiling. Controls block length only.
inline double rEff(double r_interval, int meioses) {
    int m = std::max(1, meioses);
    return 0.5 * (1.0 - std::pow(1.0 - 2.0 * r_interval, m));
}

// Relax-to-stationary dosage transition. Stationary == pi EXACTLY for any pi
// (row-stochastic: (1-r) + r*sum(pi) = 1; pi*A = pi). `pi` is the node's
// generation-aware genotype marginal from design_priors() (see sec 4).
inline std::array<std::array<double,K>,K>
makeA(double r_interval, int meioses, const Col& pi) {
    double r = rEff(r_interval, meioses);
    std::array<std::array<double,K>,K> A;
    for (int x = 0; x < K; ++x)
        for (int y = 0; y < K; ++y)
            A[x][y] = (x == y ? 1.0 - r : 0.0) + r * pi[y];
    return A;
}
```

## 7. Node / Family structures

Messages live at the sending node. `upMsg` is over the parent's state space;
`downMsg[k]` is over child k's state space. Both coincide with the recipient's
own state space.

```cpp
struct Node {
    int parent = -1;
    std::vector<int> kids;
    int  meioses = 0;
    bool hasData = false;

    Mat  emit;                   // M x K from emission_gt() over hard states (depth-blind); empty if !hasData
    Col  rho;                    // marker-0 prior      = pi_t (this node's generation) [sec 4]
    Col  pi;                     // generation stationary = pi_t; the transition relaxes to it [sec 3,4]

    Mat              upMsg;      // v -> parent      (parent states)
    std::vector<Mat> downMsg;    // v -> each child  (child states)
    Mat              belief;     // gamma_v, refreshed each sweep
};

struct Family {
    std::vector<Node> nodes;
    int root;
    int M;                       // markers (one chromosome; see §18)
    std::vector<double> r;       // per-interval recombination fraction, size M-1
};
```

## 8. Tilted forward–backward (only FB change: emission multiplier)

Effective emission = emit ⊙ tilt, where `emit` is the depth-blind `emission_gt()`
likelihood of the leaf's hard call at marker `i` (state `s` → ≈`1-err` on `s`,
`err` spread over the other two), and tilt = product of incoming messages.
Column-wise renormalization of alpha/beta is stable over 20k markers and
cancels in gamma = normalize(alpha ⊙ beta).

```cpp
Mat tiltedFB(const Node& nd, const Family& fam, const Mat& tilt) {
    const int M = fam.M;
    auto E = [&](int i, int x) {
        double e = nd.hasData ? nd.emit[i][x] : 1.0;
        return e * tilt[i][x];
    };
    Mat alpha(M), beta(M);

    for (int x = 0; x < K; ++x) alpha[0][x] = nd.rho[x] * E(0, x);
    normalize(alpha[0]);
    for (int i = 1; i < M; ++i) {
        auto A = makeA(fam.r[i-1], nd.meioses, nd.pi);
        for (int x = 0; x < K; ++x) {
            double s = 0;
            for (int y = 0; y < K; ++y) s += alpha[i-1][y] * A[y][x];
            alpha[i][x] = s * E(i, x);
        }
        normalize(alpha[i]);
    }
    beta[M-1] = {1,1,1}; normalize(beta[M-1]);
    for (int i = M-2; i >= 0; --i) {
        auto A = makeA(fam.r[i], nd.meioses, nd.pi);
        for (int x = 0; x < K; ++x) {
            double s = 0;
            for (int y = 0; y < K; ++y) s += A[x][y] * E(i+1, y) * beta[i+1][y];
            beta[i][x] = s;
        }
        normalize(beta[i]);
    }
    Mat gamma(M);
    for (int i = 0; i < M; ++i) {
        for (int x = 0; x < K; ++x) gamma[i][x] = alpha[i][x] * beta[i][x];
        normalize(gamma[i]);
    }
    return gamma;
}
```

## 9. Incoming product, cavities, damping

**Cavity is exact only at the site, not along the chain.** `belief[i]` is
proportional to `tilt[i][x]`, so `belief[i] / excludedMsg[i]` cleanly removes the
excluded message's *site-`i`* contribution. But that same message entered the FB
as a tilt at every other site `j ≠ i` and propagated into `belief[i]` through the
chain transitions; the per-site division does not remove that indirect
influence, so a neighbor's message leaks back around the 4-cycle. Leakage grows
with chain coupling and is non-trivial at 20k dense SNPs. Damping (λ ≈ 0.5)
controls it in practice; this is the standard structured-BP heuristic. If
accuracy plateaus, this leakage is the second suspect after dropped transmission
linkage (§2).

```cpp
Mat incoming(const Family& fam, int v) {
    const Node& nd = fam.nodes[v];
    Mat prod(fam.M, Col{1,1,1});
    for (int child : nd.kids) {                       // upMsg from children
        const Mat& m = fam.nodes[child].upMsg;
        for (int i = 0; i < fam.M; ++i)
            for (int x = 0; x < K; ++x) prod[i][x] *= m[i][x];
    }
    if (nd.parent >= 0) {                              // downMsg from parent
        const Node& p = fam.nodes[nd.parent];
        int slot = 0; while (p.kids[slot] != v) ++slot;
        const Mat& m = p.downMsg[slot];
        for (int i = 0; i < fam.M; ++i)
            for (int x = 0; x < K; ++x) prod[i][x] *= m[i][x];
    }
    return prod;
}

inline Col cavity(const Col& belief, const Col& excluded) {
    Col c;
    for (int x = 0; x < K; ++x) c[x] = belief[x] / std::max(excluded[x], FLOOR);
    normalize(c);
    return c;
}

double dampInto(Mat& cur, const Mat& tgt, double lambda) {
    if (cur.empty()) cur.assign(tgt.size(), Col{1,1,1});
    double maxd = 0;
    for (size_t i = 0; i < cur.size(); ++i) {
        Col nx;
        for (int x = 0; x < K; ++x)
            nx[x] = std::pow(std::max(cur[i][x], FLOOR), 1-lambda)
                  * std::pow(std::max(tgt[i][x], FLOOR),   lambda);
        normalize(nx);
        for (int x = 0; x < K; ++x)
            maxd = std::max(maxd, std::fabs(nx[x] - cur[i][x]));
        cur[i] = nx;
    }
    return maxd;
}
```

## 10. Message updates

Directions differ only by which index of Tsib is contracted:
- UP (child→parent): sum over child states, result over parent states.
- DOWN (parent→child): sum over parent states, result over child states.

```cpp
double sendUp(Family& fam, int v, double lambda) {
    Node& nd = fam.nodes[v];
    Node& p  = fam.nodes[nd.parent];
    int slot = 0; while (p.kids[slot] != v) ++slot;

    Mat tgt(fam.M);
    for (int i = 0; i < fam.M; ++i) {
        Col cav = cavity(nd.belief[i], p.downMsg[slot][i]);  // exclude parent msg
        Col out{0,0,0};
        for (int xp = 0; xp < K; ++xp)
            for (int xc = 0; xc < K; ++xc) out[xp] += Tsib[xp][xc] * cav[xc];
        normalize(out); tgt[i] = out;
    }
    return dampInto(nd.upMsg, tgt, lambda);
}

double sendDown(Family& fam, int v, int k, double lambda) {
    Node& nd = fam.nodes[v];
    Mat tgt(fam.M);
    for (int i = 0; i < fam.M; ++i) {
        Col cav = cavity(nd.belief[i], fam.nodes[nd.kids[k]].upMsg[i]); // excl child k
        Col out{0,0,0};
        for (int xc = 0; xc < K; ++xc)
            for (int xp = 0; xp < K; ++xp) out[xc] += Tsib[xp][xc] * cav[xp];
        normalize(out); tgt[i] = out;
    }
    return dampInto(nd.downMsg[k], tgt, lambda);
}
```

## 11. Driver: post-order up, pre-order down, iterate

```cpp
void runFamilyBP(Family& fam, int maxIters = 30,
                 double tol = 1e-4, double lambda = 0.5) {
    std::vector<int> post;
    std::function<void(int)> dfs = [&](int v){
        for (int c : fam.nodes[v].kids) dfs(c);
        post.push_back(v);
    };
    dfs(fam.root);
    std::vector<int> pre(post.rbegin(), post.rend());

    // Initialize ALL messages to uniform M x K. Must be non-empty before the
    // first upward sweep: incoming() and sendUp() read the parent's downMsg
    // (and children's upMsg) before any send has populated them; empty Mats
    // would be indexed out of bounds. dampInto lazily inits `cur`, but it runs
    // only inside the sends, which is too late for these reads.
    const Mat uniformMsg(fam.M, Col{1.0/3, 1.0/3, 1.0/3});
    for (Node& nd : fam.nodes) {
        nd.upMsg = uniformMsg;
        nd.downMsg.assign(nd.kids.size(), uniformMsg);
    }

    // init: independent HMM per node (flat tilt)
    for (int v : post) {
        Mat flat(fam.M, Col{1,1,1});
        fam.nodes[v].belief = tiltedFB(fam.nodes[v], fam, flat);
    }

    for (int it = 0; it < maxIters; ++it) {
        double maxd = 0;
        for (int v : post) {                          // upward sweep
            fam.nodes[v].belief = tiltedFB(fam.nodes[v], fam, incoming(fam, v));
            if (fam.nodes[v].parent >= 0)
                maxd = std::max(maxd, sendUp(fam, v, lambda));
        }
        for (int v : pre) {                           // downward sweep
            fam.nodes[v].belief = tiltedFB(fam.nodes[v], fam, incoming(fam, v));
            for (int k = 0; k < (int)fam.nodes[v].kids.size(); ++k)
                maxd = std::max(maxd, sendDown(fam, v, k, lambda));
        }
        if (maxd < tol) break;
    }
    for (int v : pre)                                 // final beliefs
        fam.nodes[v].belief = tiltedFB(fam.nodes[v], fam, incoming(fam, v));
}
```

## 12. Pedigree input: `read_pedigree()`, not a name grammar

The kernel is data-agnostic and so is `refine_ancestry()`: the pedigree enters
through the existing `read_pedigree()` adapter, which returns one row per taxon
with `taxon, family, parent1, parent2, contribution1, contribution2, F`. This
already encodes the forest — `family` keys the tree, `parent1`/`parent2` give the
edges, and latent ungenotyped ancestors appear as a `taxon` referenced as a
parent but with no genotype column in `mosaic` (`hasData=false`).
`refine_ancestry()` calls `read_pedigree()` (or accepts its data.frame), joins
`taxon` to the sample columns of `mosaic`, and builds the `Family` structs of §7
(node list, `kids`, `meioses`, `hasData`, per-node `emit`/`rho`, per-interval
`r`). No sample-name grammar lives in the package.

**Where the pedigree file comes from.** Turning project-specific sample names
into a `read_pedigree()`-readable file (PLINK `.fam`/`.ped` or FSFHap TSV) is a
`zealhmm` pipeline step, not package code. For the current data the grammar is
`Zd.<donor>_P<b2>_P<s1>_P<s2>.<a>.<b>.<c>`:
- `Zd.xxxx_Pn` prefix (donor + founder) keys the tree / `family` (FID).
- Each subsequent `_Pn` is an internal (latent) pedigree node.
- The trailing `.a.b.c` are the final branching path to the sampled leaf; each
  dotted level is one selfing generation → one edge (+1 meiosis).
- Longest-common-prefix over names reconstructs shared ancestors →
  `parent1`/`parent2`; genotyped samples become leaves, inferred ancestors become
  latent rows. `meioses` = design backcross count at the founder, +1 per
  descending generation.

A different data source supplies its own pedigree file (however produced) and
the kernel is untouched. Only leaves present in `mosaic` get `hasData=true` +
`emit` from `emission_gt()` over their hard per-marker `state` (depth-blind;
error `err`). When the FSFHap `F` column is present it pins the founder het
fraction directly (`phet = (1 - F)/2`), overriding the design default.

## 13. Complexity and resources

- Time:  O(iterations · |V| · M · K^2),  K^2 = 9 → linear in pedigree and
  markers. Converges in a few sweeps (near-deterministic kernel + damping).
- Memory: O(Σ_v (1 + |kids_v|) · M · K) doubles per family; hundreds of MB at
  most for a few hundred leaves at 20k SNPs.
- Parallelism: families are fully separable → runFamilyBP is embarrassingly
  parallel across the forest.

## 14. Numerical guards

- Keep FLOOR + max(...) in cavity division: the deterministic Tsib rows create
  hard zeros. Do NOT remove floors to "clean up."
- Column-wise renormalization in FB; damping λ ≈ 0.5 for loop stability.
- If cavity BP is unstable, fall back to non-cavity (use belief directly) +
  stronger damping: cruder, mild over-counting bias, usually fine.

## 15. Free correctness / QC checks

1. **Mendelian consistency** holds approximately after convergence (up to the
   deterministic-kernel enforcement and the per-site cavity leakage of §9). A
   strong residual violation (child P(d=0)≈1 under parent P(d=2)≈1) flags either
   non-convergence OR a pedigree-string error from the R name-parser. Run as a QC
   pass over the parsed pedigree before trusting output.
2. **Isolated d=2 spikes** at a leaf (not a coherent block traceable to a
   selfing-fixation event) suggest the emission error rate `err` is set too low.
   Distinguish from genuine rare double-recombinants by requiring the spike to
   span fewer than `min_run` consecutive markers before flagging.

## 16. Validation protocol

Using the known pedigree as ground truth:
- Simulate selfing families with a gamete simulator (known ancestry mosaics).
- Down-sample reads to real depth (λ ≈ 2.8×).
- Metrics: (a) dosage accuracy, (b) breakpoint localization error.
- Compare per-individual HMM vs. pedigree-BP, stratified by depth and by number
  of genotyped sibs.
- Expected gains largest at low depth, near breakpoints, and for leaves with
  many genotyped siblings.

## 17. R / C++ boundary, decode → segments, and API surface

The kernel returns beliefs; R turns them into the package's deliverable and owns
everything data-specific.

- **C++ (pure kernel)**: `runFamilyBP(Family&, ...)` mutates `belief` in place.
  Exposed via one `// [[Rcpp::export]]` entry (`pedigree_bp_cpp`) taking the
  `Family` fields as plain R vectors/lists (edges, `meioses`, `hasData`, `emit`,
  `rho`, `pi`, `r`, `maxIters`, `tol`, `lambda`) and returning the per-node
  `belief` matrices. Regenerate with `Rcpp::compileAttributes()` +
  `devtools::document()`; never hand-edit `RcppExports.*`.
  - **Boundary validation (fail fast in R, before the kernel runs)**: the R
    wrapper `stopifnot`-checks `M >= 1`, `length(r) == M-1`, every `emit`/message
    matrix is `M × 3`, all `parent`/`kids` indices are in range and acyclic, and
    every `rho`/`pi` row is finite, non-negative, and sums to > 0. The kernel
    trusts its inputs (tight inner loops over 20k × |V|); malformed R input must
    never reach it, so all shape/finiteness/index guards live at this boundary.
- **Decode → segments (R)**: per genotyped leaf, `argmax` over `belief[i]` gives
  dosage 0/1/2, mapped to `REF`/`HET`/`ALT`, then RLE to the common schema via
  the existing `to_segments()` (`segments.cpp`). Only `hasData` leaves are
  emitted; latent internal nodes are discarded. Optionally also return the
  posterior `belief` for downstream soft use.
- **Chromosome handling**: `Family` holds one chromosome (single `M`, single `r`
  of length `M-1`), consistent with nilHMM's chromosome-separate rule and sorted
  markers. R loops families × chromosomes and concatenates segments; the founder
  marker-0 `rho` is applied per chromosome.
- **API surface**: a standalone verb, not a subcase of the per-individual
  `call_ancestry(data, params)` contract, because the unit of work is a
  family/forest. The input is the **per-marker hard-state mosaic** from
  `call_states()` — NOT `call_ancestry()`'s RLE segments. Reason: BP runs on the
  marker grid, and the per-marker table carries it (`chr, pos, state`), whereas
  RLE segments (`start_bp, end_bp`) have dropped the grid. No depth, no counts
  needed. **Recombination source**: `pos` alone does *not* give the per-interval
  `r` — bp gaps become recombination fractions only through a rate or a genetic
  map. `refine_ancestry` derives the size-`M-1` `r` vector exactly as the engine
  does: a uniform `rrate` (default from the engine, same units as `call_states`)
  applied to `pos` deltas, or a genetic `map` (cM) when supplied (Haldane), via
  `load_map()`/`cm_to_mb`. State which was used in methods.

  Concrete end-to-end (real function names, no pseudocode):

  ```r
  # 1. per-individual calls -- this IS the count/depth modeling; it happens here, once.
  mosaic <- call_states(counts, caller = "nnil", design = "BC2S3", err = 0.01)
  #    mosaic: source, donor, name, chr, pos, state   (one 0/1/2 call per marker/leaf)

  # 2. refine those calls over the pedigree forest (depth-blind from here on).
  refined <- refine_ancestry(
    mosaic,                       # call_states() per-marker table above
    pedigree = "pop.ped",         # path (dispatched via read_pedigree) OR a read_pedigree() df
    design   = "BC2S3",           # -> per-node pi_t (rho/pi) + meioses via design_priors()
    err      = 0.02,              # emission_gt() call-error rate (see caveat below)
    rrate    = 0.01,              # recomb rate applied to pos deltas (or map = <cM> instead)
    maxiter  = 30, tol = 1e-4, lambda = 0.5
  )                               # -> refined per-marker states, same columns as `mosaic`

  calls <- to_segments(refined)   # RLE to the common segment schema when you want intervals
  ```

  Wiring details:
  - `mosaic` is the `call_states()` per-marker table being refined.
    `refine_ancestry` builds each leaf's `emit` via `emission_gt()` from its
    `state` column (depth-blind), couples leaves through the pedigree, and
    returns a refined table of the same shape. Call `to_segments()` yourself for
    intervals (mirrors `call_ancestry() = call_states() |> to_segments()`).
  - **Depth caveat (state it in methods).** Because emission is over hard calls,
    every marker is equally confident regardless of the depth behind it: a
    low-depth miscall looks as trustworthy as a high-depth correct one. `err`
    must reflect the *upstream caller's* per-marker error rate, not raw
    sequencing error, and per-marker depth weighting is unrecoverable at this
    stage. If that loss matters, the depth-aware alternative is to run the same
    BP directly on counts with the BetaBinomial `emission_count()` (deferred
    alongside the phased kernel, §18) — but that is a new caller, not "refine".
  - `pedigree` reuses `read_pedigree()` (`format = "fam"` for PLINK-style
    `FID IID PID MID`, or `"fsfhap"`). Its `taxon` column joins to `mosaic$name`;
    `parent1`/`parent2` build the forest, **including latent ungenotyped
    ancestors** (a `taxon` used as a parent but absent from `mosaic` →
    `hasData=false`). `meioses` derives from generation depth + design backcross
    count; the FSFHap `F` column, when present, pins the founder het fraction
    directly (`phet = (1 - F)/2`).
  - `design` feeds `design_priors()` for each node's `rho`/`pi` (the generation
    marginal `π_t`, §4) and `meioses`.
  - The `caller_spec` table stays unchanged; `refine_ancestry` reuses the
    engine's per-interval `r` and adds only the pedigree-BP layer.
- **Generation-aware stationary (het-inflation fix, now in v1)**: each node's
  chain relaxes toward its own generation marginal `π_t` (§3, §4), not a shared
  panmictic stationary, so deep low-coverage leaves are pulled toward the correct
  high-homozygosity distribution instead of Hardy–Weinberg. Residual het-inflation
  risk is now limited to the relax-to-`π_t` double-recomb permissiveness (§3
  tradeoff) and the cavity leakage (§9); validation (§16) should still stratify
  het accuracy by depth and selfing depth to confirm.

## 18. Deferred to v2

**Priority order** (from the design review): (1) generation-aware genotype
priors/stationary — *done in v1* (§3, §4); (2) depth-aware pedigree emission;
(3) phased copy-switch transmission; (4) more-exact cavity inference. Improve the
biological model (2, 3) before the inference approximation (4).

- **Depth-aware pedigree caller** (a new caller, not "refine"): run the same BP
  directly on read counts with the BetaBinomial `emission_count()` instead of the
  depth-blind `emission_gt()` over hard calls, so pedigree information and
  per-marker depth confidence combine in one pass. Recovers the depth weighting
  the refine stage discards (see §17 depth caveat); benchmark against
  refine-over-`call_states()` to decide whether the extra coupling earns its cost.
- Phased copy-switch transmission kernel (preserves transmission linkage).
- Breakpoint-consensus approximation (Implementation 2): pass per-interval
  recombination posteriors r_v(i) = Σ_{x≠x'} xi_v(i,x,x') and do coincidence-
  scored consensus instead of full-posterior messages. Faster, lower-dim,
  heuristic transmission semantics; benchmark against Implementation 1 once the
  accuracy ceiling is established.