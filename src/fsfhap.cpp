// FSFHap port (design/FSFHAP_PORT.md) — stage 1a: segregating-site test.
//
// Faithful port of TASSEL's NucleotideImputationUtils.whichSitesSegregateCorrectly
// (+ its private binomialProbability). Per site, a binomial "dominance" test picks
// the segregation model best explaining the minor-allele count: monomorphic/error
// (p=0.002), backcross (p=0.25), or F2/Mendelian (p=0.5). A site is kept when the
// design's model wins. Deterministic — a tight per-stage TASSEL-parity target.
//
// binomialProbability is reproduced with the SAME explicit lnGamma expansion TASSEL
// uses (GammaFunction.lnGamma), not R::dbinom, so the arithmetic matches term for
// term (Java lnGamma vs R::lgammafn agree to machine precision; the model
// comparisons are relative and robust to any residual ulp difference / underflow —
// TASSEL's Math.exp(logprob) underflows identically to std::exp here).
//
// Genotype coding is the engine's canonical g in {0 REF-hom, 1 het, 2 ALT-hom,
// 3 missing}; allele counts follow: refCount = 2*nRefHom + nHet, altCount likewise.
// This is inherently biallelic, matching FSFHap's biallelic (freq[1].length>1) gate.

#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

// TASSEL binomialProbability(trials, successes, pSuccess), verbatim formula:
//   lnGamma(n+1) - lnGamma(k+1) - lnGamma(n-k+1) + k*ln(p) + (n-k)*ln(1-p), exp'd.
static inline double binom_prob(double n, double k, double p) {
  double logprob = R::lgammafn(n + 1.0) - R::lgammafn(k + 1.0)
                 - R::lgammafn(n - k + 1.0)
                 + k * std::log(p) + (n - k) * std::log(1.0 - p);
  return std::exp(logprob);
}

//' FSFHap segregating-site test (whichSitesSegregateCorrectly), faithful port
//'
//' Per site, keep it if its minor-allele count fits the design's segregation
//' model better than the alternatives. Backcross (`ratio` 0.25 or 0.75):
//' keep iff `pquarter > phalf && pquarter > pmono`. F2 (`ratio` 0.5): keep iff
//' `phalf / (pmono + pquarter) > 2`. Sites that are monomorphic (single allele)
//' or exceed `max_missing` are dropped.
//'
//' @param G Integer matrix, taxa x sites, values 0/1/2/3 = REF-hom / het /
//'   ALT-hom / missing (engine canonical `g`). One family, one chromosome,
//'   sites sorted by position.
//' @param max_missing Max missing-genotype proportion for a site to be tested.
//' @param ratio Expected minor-allele frequency: 0.25/0.75 backcross, 0.5 F2.
//' @return List: `seg` (logical, kept sites), `Mj`/`Mn` (major/minor allele
//'   counts), `p_missing`, and the model probs `pmono`/`pquarter`/`phalf`
//'   (NA at untested sites) — for per-stage TASSEL-parity checks.
//' @keywords internal
// [[Rcpp::export]]
List fsfhap_segregating_sites_cpp(IntegerMatrix G, double max_missing, double ratio) {
  const int ntaxa  = G.nrow();
  const int nsites = G.ncol();
  LogicalVector seg(nsites);
  IntegerVector Mj(nsites), Mn(nsites);
  NumericVector p_missing(nsites), pmono(nsites), pquarter(nsites), phalf(nsites);

  for (int s = 0; s < nsites; ++s) {
    int nRef = 0, nHet = 0, nAlt = 0, nMiss = 0;
    for (int t = 0; t < ntaxa; ++t) {
      switch (G(t, s)) {
        case 0:  ++nRef;  break;
        case 1:  ++nHet;  break;
        case 2:  ++nAlt;  break;
        default: ++nMiss; break;   // 3 or any other value = missing
      }
    }
    const int refCount = 2 * nRef + nHet;      // allele counts (diploid)
    const int altCount = 2 * nAlt + nHet;
    const int mj = std::max(refCount, altCount);
    const int mn = std::min(refCount, altCount);
    const double pM = ntaxa > 0 ? (double) nMiss / (double) ntaxa : 1.0;
    Mj[s] = mj; Mn[s] = mn; p_missing[s] = pM;

    // FSFHap gate: biallelic (both alleles present) and missingness acceptable.
    const bool biallelic = (refCount > 0 && altCount > 0);
    double pmo = NA_REAL, pq = NA_REAL, ph = NA_REAL;
    bool keep = false;
    if (biallelic && pM <= max_missing) {
      const double N = mj + mn;
      pmo = binom_prob(N, mn, 0.002);
      pq  = binom_prob(N, mn, 0.25);
      ph  = binom_prob(N, mn, 0.5);
      // `ratio` is design-selected UPSTREAM (the routing dispatcher derives it
      // from the pedigree contribution, as TASSEL does) — this leaf does not guess
      // or validate it, and mirrors TASSEL's permissive else (any non-BC ratio -> F2).
      if (ratio == 0.25 || ratio == 0.75) keep = (pq > ph && pq > pmo);
      else                                keep = (ph / (pmo + pq) > 2.0);
    }
    seg[s] = keep; pmono[s] = pmo; pquarter[s] = pq; phalf[s] = ph;
  }

  return List::create(_["seg"] = seg, _["Mj"] = Mj, _["Mn"] = Mn,
                      _["p_missing"] = p_missing, _["pmono"] = pmono,
                      _["pquarter"] = pquarter, _["phalf"] = phalf);
}

// TASSEL calculateRSqr on a 2x2 allele-PRESENCE contingency (hets count toward
// both alleles). Returns NaN when under-sampled (< minTaxa) or a marginal is
// monomorphic — matching Double.NaN, which the caller treats as "keep".
static inline double calc_rsqr(int countAB, int countAb, int countaB, int countab,
                               int minTaxa) {
  const double N = (double)(countAB + countAb + countaB + countab);
  if (N < minTaxa) return std::numeric_limits<double>::quiet_NaN();
  const double freqA = (double)(countAB + countAb) / N;
  const double freqB = (double)(countAB + countaB) / N;
  if (freqA == 0.0 || freqB == 0.0 || freqA == 1.0 || freqB == 1.0)
    return std::numeric_limits<double>::quiet_NaN();
  double r = ((double)countAB / N) * ((double)countab / N)
           - ((double)countaB / N) * ((double)countAb / N);
  r *= r;
  r /= freqA * (1.0 - freqA) * freqB * (1.0 - freqB);
  return r;
}

//' FSFHap same-tag SNP filter (whichSnpsAreFromSameTag), faithful port
//'
//' Within a 64 bp window, consecutive SNPs whose presence-based R^2 to the
//' window's anchor is >= `min_rsq` are treated as coming from the same GBS tag
//' and dropped, keeping only the anchor; a SNP that is >= 64 bp away, or has
//' R^2 < `min_rsq` (or NaN), becomes the next kept anchor. **One chromosome per
//' call** (positions only; the chromosome-equality guard is implicit), markers
//' sorted by position.
//'
//' @param G Integer matrix, taxa x sites, canonical `g` in {0,1,2,3}.
//' @param pos Integer marker positions (bp), length = ncol(G), sorted ascending.
//' @param major_is_ref Logical, length = ncol(G): TRUE if REF is the major
//'   allele at that site (so allele presence is defined per the site's own
//'   major/minor, as TASSEL's `majorAllele`/`minorAllele` do).
//' @param min_rsq R^2 threshold for "same tag" (TASSEL default 0.8).
//' @return Logical vector, length = ncol(G): TRUE = keep (a distinct tag).
//' @keywords internal
// [[Rcpp::export]]
LogicalVector fsfhap_same_tag_keep_cpp(IntegerMatrix G, IntegerVector pos,
                                       LogicalVector major_is_ref, double min_rsq) {
  const int ntaxa  = G.nrow();
  const int nsites = G.ncol();
  if ((int) pos.size() != nsites)
    stop("fsfhap_same_tag_keep_cpp(): length(pos) must equal ncol(G)");
  if ((int) major_is_ref.size() != nsites)
    stop("fsfhap_same_tag_keep_cpp(): length(major_is_ref) must equal ncol(G)");
  for (int s = 1; s < nsites; ++s)                     // 64 bp spans assume sorted pos
    if (pos[s] < pos[s - 1]) stop("fsfhap_same_tag_keep_cpp(): pos must be sorted ascending");
  LogicalVector isSelected(nsites);            // all FALSE
  if (nsites == 0) return isSelected;
  isSelected[0] = true;

  // allele presence at (t, s), keyed on the site's own major/minor (hets carry
  // both; missing g==3 carries neither) — mirrors allelePresenceForAllTaxa.
  auto majPresent = [&](int t, int s) -> bool {
    const int g = G(t, s);
    return major_is_ref[s] ? (g == 0 || g == 1) : (g == 1 || g == 2);
  };
  auto minPresent = [&](int t, int s) -> bool {
    const int g = G(t, s);
    return major_is_ref[s] ? (g == 1 || g == 2) : (g == 0 || g == 1);
  };

  int firstSite = 0;
  int firstPos  = pos[0];
  while (firstSite < nsites - 1) {
    int nextSite = firstSite + 1;
    int nextPos  = pos[nextSite];
    while (nextPos - firstPos < 64) {          // single chromosome: chr guard implicit
      int AB = 0, Ab = 0, aB = 0, ab = 0;
      for (int t = 0; t < ntaxa; ++t) {
        const bool mr = majPresent(t, firstSite), nr = minPresent(t, firstSite);
        const bool mc = majPresent(t, nextSite),  nc = minPresent(t, nextSite);
        if (mr && mc) ++AB;
        if (mr && nc) ++Ab;
        if (nr && mc) ++aB;
        if (nr && nc) ++ab;
      }
      const double rsq = calc_rsqr(AB, Ab, aB, ab, 2);
      if (std::isnan(rsq) || rsq < min_rsq) { isSelected[nextSite] = true; break; }
      ++nextSite;
      if (nextSite >= nsites) break;
      nextPos = pos[nextSite];
    }
    firstSite = nextSite;
    if (firstSite < nsites) { firstPos = pos[firstSite]; isSelected[firstSite] = true; }
  }
  return isSelected;
}

// ==== Stage 2a: HaplotypeClusterer / HaplotypeCluster (faithful port) ========
// Operates on our canonical g in {0=A-hom, 1=het, 2=C-hom, 3=missing}. A
// "haplotype" is a taxon's row over the window. Only makeClusters + consensus +
// sort are ported here (the pieces BiparentalHaplotypeFinder's seed needs);
// getMergedClusters / moveAllHaplotypesToBiggestCluster / removeHeterozygousClusters
// land in stage 2b with the finder's clusterWindow orchestration.

// Haplotype.getDistance on g-codes: N at either -> 0; equal -> 0;
// het(1) vs hom -> 1; hom vs different hom -> 2.
static int seq_dist(const std::vector<int>& a, const std::vector<int>& b) {
  const int n = (int) a.size();
  int d = 0;
  for (int s = 0; s < n; ++s) {
    const int b0 = a[s], b1 = b[s];
    if (b0 == b1 || b0 == 3 || b1 == 3) continue;
    const bool h0 = (b0 == 1), h1 = (b1 == 1);
    if (h0)      { if (!h1) ++d; }   // het vs hom
    else if (h1) ++d;                // hom vs het
    else         d += 2;             // hom vs different hom
  }
  return d;
}

struct HCluster { std::vector<int> members; double score; };

// HaplotypeClusterer.makeClusters(): 0-distance clustering, multi-membership via
// missing, fractional score 1/count. Two passes, exactly as TASSEL:
//   pass 1 — build cluster reps (haplotypes >0 from every existing rep);
//   pass 2 — add each remaining haplotype to every cluster it is 0-distance to
//            (checking all current members), splitting score 1/count.
static std::vector<HCluster> make_clusters(const std::vector<std::vector<int> >& pool) {
  std::vector<HCluster> clusters;
  const int ntaxa = (int) pool.size();
  if (ntaxa == 0) return clusters;
  HCluster c0; c0.members.push_back(0); c0.score = 1.0; clusters.push_back(c0);

  std::vector<int> leftover;
  for (int t = 1; t < ntaxa; ++t) {
    bool inCluster = false;
    for (size_t c = 0; c < clusters.size(); ++c)
      if (seq_dist(pool[clusters[c].members[0]], pool[t]) == 0) { inCluster = true; break; }
    if (!inCluster) { HCluster nc; nc.members.push_back(t); nc.score = 1.0; clusters.push_back(nc); }
    else leftover.push_back(t);
  }

  const int nclusters = (int) clusters.size();   // fixed before pass 2 (no new clusters)
  for (size_t li = 0; li < leftover.size(); ++li) {
    const int t = leftover[li];
    std::vector<char> incl(nclusters, 1);
    int count = 0;
    for (int c = 0; c < nclusters; ++c) {
      for (size_t mi = 0; mi < clusters[c].members.size(); ++mi)
        if (seq_dist(pool[clusters[c].members[mi]], pool[t]) > 0) { incl[c] = 0; break; }
      if (incl[c]) { ++count; clusters[c].members.push_back(t); }
    }
    const double hs = count > 0 ? 1.0 / (double) count : 0.0;
    for (int c = 0; c < nclusters; ++c) if (incl[c]) clusters[c].score += hs;
  }
  return clusters;
}

// HaplotypeCluster.getMajorityHaplotype(): per site, majority g among non-missing
// members; tie or all-missing -> N(3). (== unanimous for 0-distance clusters.)
static std::vector<int> cluster_majority(const HCluster& c,
    const std::vector<std::vector<int> >& pool, int nsites) {
  std::vector<int> maj(nsites, 3);
  for (int s = 0; s < nsites; ++s) {
    int cnt[3] = {0, 0, 0};
    for (size_t mi = 0; mi < c.members.size(); ++mi) {
      const int g = pool[c.members[mi]][s];
      if (g >= 0 && g <= 2) ++cnt[g];
    }
    int maxc = 0; for (int a = 0; a < 3; ++a) if (cnt[a] > maxc) maxc = cnt[a];
    if (maxc == 0) { maj[s] = 3; continue; }
    int nmax = 0, which = -1;
    for (int a = 0; a < 3; ++a) if (cnt[a] == maxc) { ++nmax; which = a; }
    maj[s] = (nmax == 1) ? which : 3;
  }
  return maj;
}

// HaplotypeCluster.getUnanimousHaplotype(): per site, the single non-missing g if
// exactly one distinct value is present, else N(3).
static std::vector<int> cluster_unanimous(const HCluster& c,
    const std::vector<std::vector<int> >& pool, int nsites) {
  std::vector<int> un(nsites, 3);
  for (int s = 0; s < nsites; ++s) {
    int seen = -1; bool multi = false;
    for (size_t mi = 0; mi < c.members.size(); ++mi) {
      const int g = pool[c.members[mi]][s];
      if (g == 3 || g < 0 || g > 2) continue;
      if (seen == -1) seen = g; else if (seen != g) { multi = true; break; }
    }
    un[s] = (seen != -1 && !multi) ? seen : 3;
  }
  return un;
}

// --- clusterer operations (HaplotypeClusterer merge/move/removeHet) ----------
// sort order: score desc, then size desc (HaplotypeCluster.compareTo)
static bool cluster_less(const HCluster& a, const HCluster& b) {
  if (a.score != b.score) return a.score > b.score;
  return a.members.size() > b.members.size();
}
// mergeTwoClusters: c0 absorbs c1's unique members (dedup by taxon index)
static void merge_two(HCluster& c0, const HCluster& c1) {
  for (size_t i = 0; i < c1.members.size(); ++i)
    if (std::find(c0.members.begin(), c0.members.end(), c1.members[i]) == c0.members.end())
      c0.members.push_back(c1.members[i]);
}
// clusterDistanceMaxPairDiff: max distance over NON-shared cross pairs
static int cluster_maxpair_diff(const HCluster& c0, const HCluster& c1,
                                const std::vector<std::vector<int> >& pool) {
  std::vector<int> s0, s1;
  for (size_t i = 0; i < c0.members.size(); ++i)
    if (std::find(c1.members.begin(), c1.members.end(), c0.members[i]) == c1.members.end()) s0.push_back(c0.members[i]);
  for (size_t i = 0; i < c1.members.size(); ++i)
    if (std::find(c0.members.begin(), c0.members.end(), c1.members[i]) == c0.members.end()) s1.push_back(c1.members[i]);
  int mx = 0;
  for (size_t a = 0; a < s0.size(); ++a) for (size_t b = 0; b < s1.size(); ++b)
    mx = std::max(mx, seq_dist(pool[s0[a]], pool[s1[b]]));
  return mx;
}
// recalculateScores: score = Σ_hap 1/(#clusters containing hap); drop score-0
static void recalc_scores(std::vector<HCluster>& cl, int ntaxa) {
  for (size_t c = 0; c < cl.size(); ++c) cl[c].score = 0.0;
  for (int t = 0; t < ntaxa; ++t) {
    int cnt = 0;
    for (size_t c = 0; c < cl.size(); ++c)
      if (std::find(cl[c].members.begin(), cl[c].members.end(), t) != cl[c].members.end()) ++cnt;
    if (cnt == 0) continue;
    const double add = 1.0 / (double) cnt;
    for (size_t c = 0; c < cl.size(); ++c)
      if (std::find(cl[c].members.begin(), cl[c].members.end(), t) != cl[c].members.end()) cl[c].score += add;
  }
  std::vector<HCluster> keep;
  for (size_t c = 0; c < cl.size(); ++c) if (cl[c].score != 0.0) keep.push_back(cl[c]);
  cl.swap(keep);
}
// mergeClusters(maxdiff): sequential merge (head absorbs ≤maxdiff neighbours), recalc
static void merge_clusters(std::vector<HCluster>& cl, int maxdiff,
                           const std::vector<std::vector<int> >& pool, int ntaxa) {
  std::stable_sort(cl.begin(), cl.end(), cluster_less);
  std::vector<HCluster> cand(cl.begin(), cl.end()), merged;
  size_t head_i = 0;
  while (head_i < cand.size()) {
    HCluster head = cand[head_i];
    std::vector<HCluster> rest;
    for (size_t i = head_i + 1; i < cand.size(); ++i) {
      if (cluster_maxpair_diff(head, cand[i], pool) <= maxdiff) merge_two(head, cand[i]);
      else rest.push_back(cand[i]);
    }
    merged.push_back(head);
    cand = rest; head_i = 0;
  }
  cl = merged;
  recalc_scores(cl, ntaxa);
}
// moveAllPossibleHaplotypesToCluster(idx, higherOnly, maxdiff): move haps within
// maxdiff of cluster idx's unanimous haplotype up from higher-indexed clusters
static void move_to_cluster(std::vector<HCluster>& cl, int idx, int maxdiff,
                            const std::vector<std::vector<int> >& pool, int nsites) {
  if (cl[idx].members.empty()) return;
  const std::vector<int> chap = cluster_unanimous(cl[idx], pool, nsites);
  for (int c = idx + 1; c < (int) cl.size(); ++c) {
    std::vector<int> newmem;
    for (size_t i = 0; i < cl[c].members.size(); ++i) {
      const int m = cl[c].members[i];
      if (seq_dist(chap, pool[m]) <= maxdiff) {
        if (std::find(cl[idx].members.begin(), cl[idx].members.end(), m) == cl[idx].members.end())
          cl[idx].members.push_back(m);
      } else newmem.push_back(m);
    }
    cl[c].members.swap(newmem);
  }
}
static void move_to_biggest(std::vector<HCluster>& cl, int maxdiff,
                            const std::vector<std::vector<int> >& pool, int nsites, int ntaxa) {
  std::stable_sort(cl.begin(), cl.end(), cluster_less);
  const int n = (int) cl.size();
  for (int i = 0; i < n; ++i) move_to_cluster(cl, i, maxdiff, pool, nsites);
  std::vector<HCluster> keep;
  for (size_t c = 0; c < cl.size(); ++c) if (!cl[c].members.empty()) keep.push_back(cl[c]);
  cl.swap(keep);
  recalc_scores(cl, ntaxa);
  std::stable_sort(cl.begin(), cl.end(), cluster_less);
}
// countHeterozygousSites: per site, het if the most-common g is het(1), OR the
// 2nd-most-frequent g count > 1 (HaplotypeCluster.countHeterozygousSites).
static int count_het_sites(const HCluster& c, const std::vector<std::vector<int> >& pool, int nsites) {
  int nhet = 0;
  for (int s = 0; s < nsites; ++s) {
    int cnt[3] = {0, 0, 0};
    for (size_t i = 0; i < c.members.size(); ++i) { const int g = pool[c.members[i]][s]; if (g >= 0 && g <= 2) ++cnt[g]; }
    int maxc = 0, maxv = -1;
    for (int a = 0; a < 3; ++a) if (cnt[a] > maxc) { maxc = cnt[a]; maxv = a; }
    if (maxv == 1) { ++nhet; continue; }                 // most-common allele is het
    int second = 0;
    for (int a = 0; a < 3; ++a) if (a != maxv && cnt[a] > second) second = cnt[a];
    if (second > 1) ++nhet;                              // a real 2nd allele (count > 1)
  }
  return nhet;
}
static void remove_het_clusters(std::vector<HCluster>& cl, int maxHet,
                                const std::vector<std::vector<int> >& pool, int nsites) {
  std::vector<HCluster> keep;
  for (size_t c = 0; c < cl.size(); ++c) if (count_het_sites(cl[c], pool, nsites) <= maxHet) keep.push_back(cl[c]);
  cl.swap(keep);
}

//' FSFHap stage 2a: cluster a window of parent-called haplotypes
//'
//' Faithful port of `HaplotypeClusterer.makeClusters` + `HaplotypeCluster`
//' consensus, on the parent-origin frame (`g` in {0 A-hom, 1 het, 2 C-hom,
//' 3 missing}). Clusters group taxa whose window haplotypes are 0-distance
//' (identical modulo missing); a haplotype 0-distance to members of several
//' clusters joins all of them with fractional score `1/count`. Clusters are
//' returned sorted by score (desc) then size (desc) — `HaplotypeCluster.compareTo`.
//'
//' @param Gw Integer matrix, taxa x window-sites, canonical `g` in {0,1,2,3}.
//' @param maxdiff Distance threshold for `merge`/`move_biggest` (TASSEL
//'   `maxDifferenceScore`, 0 on the BC/finder path).
//' @param merge Apply `mergeClusters(maxdiff)` after `makeClusters`
//'   (clusterWindow does this only when `maxdiff > 0`).
//' @param move_biggest Apply `moveAllHaplotypesToBiggestCluster(maxdiff)`.
//' @param max_het If `>= 0`, drop clusters with more than `max_het` heterozygous
//'   sites (`removeHeterozygousClusters`; the finder passes `maxdiff + 5`).
//' @return List: `size`, `score` (per cluster); `majority`, `unanimous`
//'   (clusters x sites consensus, `3` = N); `members` (list of 1-based taxon
//'   indices per cluster). For 0-distance clusters `majority == unanimous`;
//'   they diverge only after merges.
//' @keywords internal
// [[Rcpp::export]]
List fsfhap_cluster_window_cpp(IntegerMatrix Gw, int maxdiff = 0, bool merge = false,
                               bool move_biggest = false, int max_het = -1) {
  const int ntaxa = Gw.nrow(), nsites = Gw.ncol();
  std::vector<std::vector<int> > pool(ntaxa, std::vector<int>(nsites));
  for (int t = 0; t < ntaxa; ++t)
    for (int s = 0; s < nsites; ++s) pool[t][s] = Gw(t, s);

  std::vector<HCluster> clusters = make_clusters(pool);
  if (merge) merge_clusters(clusters, maxdiff, pool, ntaxa);        // clusterWindow (maxdif>0)
  std::stable_sort(clusters.begin(), clusters.end(), cluster_less); // sortClusters
  if (move_biggest) move_to_biggest(clusters, maxdiff, pool, nsites, ntaxa);
  if (max_het >= 0) remove_het_clusters(clusters, max_het, pool, nsites);
  std::stable_sort(clusters.begin(), clusters.end(), cluster_less);

  const int nc = (int) clusters.size();
  IntegerVector size(nc); NumericVector score(nc);
  IntegerMatrix majority(nc, nsites), unanimous(nc, nsites);
  List members(nc);
  for (int c = 0; c < nc; ++c) {
    size[c] = (int) clusters[c].members.size();
    score[c] = clusters[c].score;
    std::vector<int> maj = cluster_majority(clusters[c], pool, nsites);
    std::vector<int> un  = cluster_unanimous(clusters[c], pool, nsites);
    for (int s = 0; s < nsites; ++s) { majority(c, s) = maj[s]; unanimous(c, s) = un[s]; }
    IntegerVector mem(clusters[c].members.size());
    for (size_t mi = 0; mi < clusters[c].members.size(); ++mi) mem[mi] = clusters[c].members[mi] + 1; // 1-based
    members[c] = mem;
  }
  return List::create(_["size"] = size, _["score"] = score,
                      _["majority"] = majority, _["unanimous"] = unanimous,
                      _["members"] = members);
}

// ==== Stage 3: imputeUsingViterbiFiveState (the REAL FSFHap imputation) =======
// 5-state Viterbi-training EM (via ViterbiAlgorithmPlugin). States = allele-freq
// classes {AA, 3A:1C, 1A:1C, 1A:3C, CC}. obs from parent-called g: 0=A-hom, 1=het,
// 2=C-hom (missing skipped). Emission init = Table 1; transition init = FILLIN
// Table 2, distance-rescaled per node (Haldane map fn), both re-estimated each
// iteration from Viterbi state counts; converge when the 5x3 emission-count matrix
// is unchanged (<=50 iters). Output: state 0->A(0), 1-3->het(1), 4->C(2).

// per-node distance-scaled transition (TransitionProbability.setNode), log-space.
static void setnode_logtrans(const double base[5][5], double segLen, double avgSeg,
                             double out[5][5]) {
  for (int r = 0; r < 5; ++r) {
    double offdiag = 0.0;
    for (int c = 0; c < 5; ++c) if (c != r) {
      double b = base[r][c];
      if (b < 0.0) b = 0.0; if (b > 0.4999) b = 0.4999;   // guard log(1-2b)
      const double m = -std::log(1.0 - 2.0 * b) * segLen / avgSeg / 2.0;
      double p = (1.0 - std::exp(-2.0 * m)) / 2.0;
      out[r][c] = p; offdiag += p;
    }
    out[r][r] = 1.0 - offdiag;
  }
  for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c)
    out[r][c] = out[r][c] > 0.0 ? std::log(out[r][c]) : R_NegInf;
}

// 5-state log-space Viterbi for one taxon; distance-scaled transition per node.
// Tie-break: strict > => lowest state index wins (matches TASSEL ViterbiAlgorithm).
static void viterbi5(const std::vector<int>& obs, const std::vector<int>& pos,
                     const double emitLog[5][3], const double base[5][5],
                     double avgSeg, const double pTrueLog[5], std::vector<int>& out) {
  const int T = (int) obs.size();
  out.assign(T, 0);
  if (T == 0) return;
  std::vector<std::array<double,5> > delta(T);
  std::vector<std::array<int,5> > psi(T);
  for (int k = 0; k < 5; ++k) delta[0][k] = pTrueLog[k] + emitLog[k][obs[0]];
  double lt[5][5];
  for (int t = 1; t < T; ++t) {
    const double segLen = std::abs((double)(pos[t] - pos[t - 1]));
    setnode_logtrans(base, segLen, avgSeg, lt);
    for (int k = 0; k < 5; ++k) {
      double best = R_NegInf; int arg = 0;
      for (int j = 0; j < 5; ++j) {
        const double v = delta[t - 1][j] + lt[j][k];
        if (v > best) { best = v; arg = j; }
      }
      delta[t][k] = best + emitLog[k][obs[t]];
      psi[t][k] = arg;
    }
  }
  double best = R_NegInf; int last = 0;
  for (int k = 0; k < 5; ++k) if (delta[T - 1][k] > best) { best = delta[T - 1][k]; last = k; }
  out[T - 1] = last;
  for (int t = T - 1; t > 0; --t) out[t - 1] = psi[t][out[t]];
}

//' FSFHap stage 3: 5-state Viterbi-training EM imputation (imputeUsingViterbiFiveState)
//'
//' The real FSFHap imputation step (via `ViterbiAlgorithmPlugin`). Faithful port
//' of `imputeUsingViterbiFiveState`: per-taxon 5-state Viterbi on non-missing
//' parent calls with a distance-scaled transition, EM re-estimating emission +
//' transition from state-count matrices until the emission-count matrix stabilizes.
//'
//' @param G Integer matrix, taxa x sites, parent-called `g` in {0 A-hom, 1 het,
//'   2 C-hom, 3 missing} (stage-1b output); one family, one chromosome, sorted.
//' @param pos Integer marker positions (bp), length = ncol(G).
//' @param phet Design-derived expected heterozygosity (`(1-F)/2`); sets the
//'   initial state distribution `{phom, .25 phet, .5 phet, .25 phet, phom}`.
//' @param max_iter EM iteration cap (TASSEL 50).
//' @return List: `imputed` (taxa x sites, `0`=A / `1`=het / `2`=C / `3`=missing on
//'   undecoded sites), `iters` (EM iterations run), `emission` (final 5x3).
//' @keywords internal
// [[Rcpp::export]]
List fsfhap_impute_five_state_cpp(IntegerMatrix G, IntegerVector pos,
                                  double phet, int max_iter = 50) {
  const int ntaxa = G.nrow(), nsites = G.ncol();
  if ((int) pos.size() != nsites)
    stop("fsfhap_impute_five_state_cpp(): length(pos) must equal ncol(G)");
  for (int s = 1; s < nsites; ++s)
    if (pos[s] < pos[s - 1]) stop("fsfhap_impute_five_state_cpp(): pos must be sorted ascending");
  // per-taxon observation + position sequences (non-missing sites only)
  std::vector<std::vector<int> > obs(ntaxa), opos(ntaxa), osite(ntaxa);
  for (int t = 0; t < ntaxa; ++t)
    for (int s = 0; s < nsites; ++s) {
      const int g = G(t, s);
      if (g == 0 || g == 1 || g == 2) { obs[t].push_back(g); opos[t].push_back(pos[s]); osite[t].push_back(s); }
    }
  // guard zero-span (all positions equal) so avgSeg is never 0 in setnode_logtrans
  const double rawChrLen = nsites > 1 ? (double)(pos[nsites - 1] - pos[0]) : 1.0;
  const double chrLen = rawChrLen > 0.0 ? rawChrLen : 1.0;

  double emission[5][3] = {{.998,.001,.001},{.6,.2,.2},{.4,.2,.4},{.2,.2,.6},{.001,.001,.998}};
  double base[5][5] = {                                   // FILLIN Table 2
    {.999,.0001,.0003,.0001,.0005},{.0002,.999,.00005,.00005,.0002},
    {.0002,.00005,.999,.00005,.0002},{.0002,.00005,.00005,.999,.0002},
    {.0005,.0001,.0003,.0001,.999}};
  double avgSeg = nsites > 0 ? chrLen / (double) nsites : 1.0;
  const double phom = (1.0 - phet) / 2.0;
  double pTrueLog[5];
  { double pt[5] = {phom, .25*phet, .5*phet, .25*phet, phom};
    for (int k = 0; k < 5; ++k) pTrueLog[k] = pt[k] > 0 ? std::log(pt[k]) : R_NegInf; }

  std::vector<std::vector<int> > states(ntaxa);
  int prevEmC[5][3]; for (int r=0;r<5;++r) for (int c=0;c<3;++c) prevEmC[r][c] = -1;
  int iter = 0;
  for (; iter < max_iter; ++iter) {
    double emitLog[5][3];
    for (int r=0;r<5;++r) for (int c=0;c<3;++c) emitLog[r][c] = emission[r][c] > 0 ? std::log(emission[r][c]) : R_NegInf;
    int transC[5][5] = {{0}}, emC[5][3] = {{0}};
    for (int t = 0; t < ntaxa; ++t) {
      viterbi5(obs[t], opos[t], emitLog, base, avgSeg, pTrueLog, states[t]);
      const int T = (int) states[t].size();
      for (int s = 0; s < T; ++s) emC[states[t][s]][obs[t][s]]++;
      for (int s = 1; s < T; ++s) transC[states[t][s-1]][states[t][s]]++;
    }
    bool converged = true;
    for (int r=0;r<5;++r) for (int c=0;c<3;++c) if (emC[r][c] != prevEmC[r][c]) { converged = false; prevEmC[r][c] = emC[r][c]; }
    // re-estimate transition (setTransitionCounts): base = row-normalized counts;
    // avgSeg = chrLen*ntaxa/totalTransitions
    long total = 0; int rowSum5[5] = {0};
    for (int r=0;r<5;++r) for (int c=0;c<5;++c) { total += transC[r][c]; rowSum5[r] += transC[r][c]; }
    if (total > 0) avgSeg = chrLen * (double) ntaxa / (double) total;
    for (int r=0;r<5;++r) if (rowSum5[r] > 0) for (int c=0;c<5;++c) base[r][c] = (double) transC[r][c] / (double) rowSum5[r];
    // re-estimate emission: row-normalized counts (keep prior row if unused)
    for (int r=0;r<5;++r) { int rs=0; for (int c=0;c<3;++c) rs += emC[r][c];
      if (rs > 0) for (int c=0;c<3;++c) emission[r][c] = (double) emC[r][c] / (double) rs; }
    if (converged) { ++iter; break; }
  }

  // output: map final states to A/het/C on decoded sites; 3 (missing) elsewhere
  IntegerMatrix imputed(ntaxa, nsites);
  std::fill(imputed.begin(), imputed.end(), 3);
  const int st2geno[5] = {0, 1, 1, 1, 2};   // 0->A, 1-3->het, 4->C
  for (int t = 0; t < ntaxa; ++t)
    for (size_t s = 0; s < states[t].size(); ++s)
      imputed(t, osite[t][s]) = st2geno[states[t][s]];
  NumericMatrix emOut(5, 3);
  for (int r=0;r<5;++r) for (int c=0;c<3;++c) emOut(r,c) = emission[r][c];
  return List::create(_["imputed"] = imputed, _["iters"] = iter, _["emission"] = emOut);
}

//' FSFHap stage 3: forward-fill gaps (fillGapsInAlignment)
//'
//' Per taxon (row), across sites in order: when two non-missing calls flanking a
//' run of missing (`3`) are EQUAL, fill the run with that value; a differing
//' non-missing call resets the anchor (no fill). Faithful to `fillGapsInAlignment`.
//'
//' @param G Integer matrix, taxa x sites, `0`/`1`/`2`/`3` (`3` = missing).
//' @return `G` with eligible missing runs forward-filled.
//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix fsfhap_fill_gaps_cpp(IntegerMatrix G) {
  const int ntaxa = G.nrow(), nsites = G.ncol();
  IntegerMatrix out = clone(G);
  for (int t = 0; t < ntaxa; ++t) {
    int segStart = -1, segVal = 3;
    for (int s = 0; s < nsites; ++s) {
      const int g = out(t, s);
      if (g == 3) continue;                       // missing: skip
      if (segStart >= 0 && g == segVal)           // equal flanks -> fill the gap
        for (int i = segStart + 1; i < s; ++i) out(t, i) = segVal;
      segStart = s; segVal = g;                   // advance/reset anchor
    }
  }
  return out;
}

// ==== Stage 2b (preFilterSites): filterSnpsByTag + computeRForMissingness ======
// phi-correlation of presence (non-missing) patterns between two sites — TASSEL's
// computeRForMissingness (presence = major∪minor allele present = g != 3).
static double r_for_missingness(const IntegerMatrix& G, int s1, int s2, int ntaxa) {
  long c1 = 0, c2 = 0, prod = 0;
  for (int t = 0; t < ntaxa; ++t) {
    const bool p1 = G(t, s1) != 3, p2 = G(t, s2) != 3;
    c1 += p1; c2 += p2; if (p1 && p2) ++prod;
  }
  const double N = (double) ntaxa;
  const double num = (double) prod - (double) c1 * (double) c2 / N;
  const double den = ((double) c1 * (N - c1) / N) * ((double) c2 * (N - c2) / N);
  if (den == 0.0) return NA_REAL;
  return num / std::sqrt(den);
}

//' FSFHap preFilterSites step 1: filterSnpsByTag (faithful port)
//'
//' Thins SNPs from the same GBS tag and applies per-site quality gates. Finds the
//' first site passing MAF/missing/het thresholds (the head), then keeps each later
//' site that (is >= 64 bp from the head OR has presence-correlation `< 0.7` to it)
//' AND passes the gates; the head advances to each kept site. Faithful to TASSEL's
//' `filterSnpsByTag(a, minMaf, maxMissing, maxHet)`.
//'
//' @param G Integer matrix, taxa x sites, canonical `g` in {0,1,2,3}; one
//'   chromosome, sites sorted by position.
//' @param pos Integer marker positions (bp), length = ncol(G).
//' @param min_maf,max_missing,max_het Per-site gates (minor-allele freq,
//'   missing proportion, heterozygous proportion). preFilterSites calls this with
//'   `(minMaf, 1 - minCoverage, 1.0)`.
//' @return Logical vector, length = ncol(G): TRUE = kept site.
//' @keywords internal
// [[Rcpp::export]]
LogicalVector fsfhap_filter_snps_by_tag_cpp(IntegerMatrix G, IntegerVector pos,
                                            double min_maf, double max_missing, double max_het) {
  const int ntaxa = G.nrow(), nsites = G.ncol();
  if ((int) pos.size() != nsites)
    stop("fsfhap_filter_snps_by_tag_cpp(): length(pos) must equal ncol(G)");
  for (int s = 1; s < nsites; ++s)                     // 64 bp spans assume sorted pos
    if (pos[s] < pos[s - 1]) stop("fsfhap_filter_snps_by_tag_cpp(): pos must be sorted ascending");
  LogicalVector keep(nsites);
  if (nsites == 0) return keep;

  // per-site MAF / missing / het proportions
  std::vector<double> maf(nsites), pmiss(nsites), phet(nsites);
  for (int s = 0; s < nsites; ++s) {
    int nRef = 0, nHet = 0, nAlt = 0, nMiss = 0;
    for (int t = 0; t < ntaxa; ++t) switch (G(t, s)) {
      case 0: ++nRef; break; case 1: ++nHet; break; case 2: ++nAlt; break; default: ++nMiss; }
    const int refC = 2 * nRef + nHet, altC = 2 * nAlt + nHet, tot = refC + altC;
    const int npresent = ntaxa - nMiss;
    maf[s]   = tot > 0 ? (double) std::min(refC, altC) / (double) tot : 0.0;
    pmiss[s] = ntaxa > 0 ? (double) nMiss / (double) ntaxa : 1.0;
    phet[s]  = npresent > 0 ? (double) nHet / (double) npresent : 0.0;
  }
  auto ok = [&](int s) { return maf[s] >= min_maf && pmiss[s] <= max_missing && phet[s] <= max_het; };

  int head = -1;
  do { ++head; if (head >= nsites) return keep; } while (!ok(head));  // first quality site
  keep[head] = true;
  // Faithful to TASSEL (`for (int s = 1; ...)`): sites before `head` are the first
  // quality site's predecessors and all fail `ok(s)`, so `head` never regresses —
  // starting at 1 is result-identical to head+1 and matches the source verbatim.
  for (int s = 1; s < nsites; ++s) {
    const int dist = pos[s] - pos[head];
    if ((dist >= 64 || r_for_missingness(G, head, s, ntaxa) < 0.7) && ok(s)) {
      keep[s] = true; head = s;
    }
  }
  return keep;
}

// LD r^2 (Hill-Robertson, = calc_rsqr) with HetTreatment.Homozygous: hets are
// treated as MISSING, so only homozygotes contribute to the presence contingency.
static double ld_rsqr_hom(const IntegerMatrix& G, int s1, int s2,
                          bool mir1, bool mir2, int ntaxa) {
  int AB = 0, Ab = 0, aB = 0, ab = 0;
  for (int t = 0; t < ntaxa; ++t) {
    const int g1 = G(t, s1), g2 = G(t, s2);
    const bool maj1 = mir1 ? (g1 == 0) : (g1 == 2);   // hom-major only (het/missing excluded)
    const bool min1 = mir1 ? (g1 == 2) : (g1 == 0);
    const bool maj2 = mir2 ? (g2 == 0) : (g2 == 2);
    const bool min2 = mir2 ? (g2 == 2) : (g2 == 0);
    if (maj1 && maj2) ++AB; if (maj1 && min2) ++Ab;
    if (min1 && maj2) ++aB; if (min1 && min2) ++ab;
  }
  return calc_rsqr(AB, Ab, aB, ab, 2);
}

//' FSFHap preFilterSites (full): filterSnpsByTag -> het-deviation -> biallelic -> LD
//'
//' Faithful port of `BiparentalHaplotypeFinder.preFilterSites`. Step 1 =
//' [fsfhap_filter_snps_by_tag_cpp()] with `(minMaf, 1 - minCoverage, 1.0)`. On the
//' result: (2) drop sites whose het fraction exceeds `mean + maxHetDeviation * sd`
//' (sample sd over the filtered sites); (3) drop monomorphic (non-biallelic) sites;
//' (4) if `minR2 > 0`, drop sites whose average `r^2` over a ±50-filtered-site window
//' (hets→missing, [calc_rsqr]) is `< minR2` (NaN pairs skipped; all-NaN average is
//' not a rejection, matching TASSEL's `NaN < minR2` = false).
//'
//' @param G Integer matrix, taxa x sites, canonical `g` in {0,1,2,3}; one chromosome.
//' @param pos Integer marker positions (bp), length = ncol(G).
//' @param min_maf,min_coverage,max_het_deviation,min_r2 TASSEL fields
//'   (0.05 / 0.2 / 5 / 0.2).
//' @return Logical vector, length = ncol(G): TRUE = site kept by preFilterSites.
//' @keywords internal
// [[Rcpp::export]]
LogicalVector fsfhap_prefilter_sites_cpp(IntegerMatrix G, IntegerVector pos,
    double min_maf, double min_coverage, double max_het_deviation, double min_r2) {
  const int ntaxa = G.nrow(), nsites = G.ncol();
  if ((int) pos.size() != nsites)
    stop("fsfhap_prefilter_sites_cpp(): length(pos) must equal ncol(G)");
  LogicalVector out(nsites);
  if (nsites == 0) return out;

  // step 1: filterSnpsByTag(minMaf, 1 - minCoverage, 1.0)
  LogicalVector keep1 = fsfhap_filter_snps_by_tag_cpp(G, pos, min_maf, 1.0 - min_coverage, 1.0);
  std::vector<int> f;                                   // filtered original-site indices
  for (int s = 0; s < nsites; ++s) if (keep1[s]) f.push_back(s);
  const int nf = (int) f.size();
  if (nf == 0) return out;

  // per filtered-site: het fraction, major-is-ref, biallelic
  std::vector<double> pHet(nf);
  std::vector<char> mir(nf), sel(nf, 1);
  for (int j = 0; j < nf; ++j) {
    const int s = f[j];
    int nRef = 0, nHet = 0, nAlt = 0, nMiss = 0;
    for (int t = 0; t < ntaxa; ++t) switch (G(t, s)) {
      case 0: ++nRef; break; case 1: ++nHet; break; case 2: ++nAlt; break; default: ++nMiss; }
    const int notmiss = ntaxa - nMiss, refC = 2 * nRef + nHet, altC = 2 * nAlt + nHet;
    pHet[j] = notmiss > 0 ? (double) nHet / (double) notmiss : 0.0;
    mir[j]  = refC >= altC;
    if (refC == 0 || altC == 0) sel[j] = 0;             // step 3: biallelic (drop monomorphic)
  }

  // step 2: het-deviation — drop sites with pHet > mean + maxHetDeviation * sd
  if (nf >= 2) {
    double mean = 0.0; for (int j = 0; j < nf; ++j) mean += pHet[j]; mean /= nf;
    double var = 0.0; for (int j = 0; j < nf; ++j) { const double d = pHet[j] - mean; var += d * d; }
    var /= (nf - 1);                                    // StatUtils.variance (sample, n-1)
    const double maxPhet = mean + max_het_deviation * std::sqrt(var);
    for (int j = 0; j < nf; ++j) if (pHet[j] > maxPhet) sel[j] = 0;
  }

  // step 4: LD filter — average r^2 over ±50-filtered-site window (hom-only)
  if (min_r2 > 0.0) {
    const int w = 50;
    for (int j = 0; j < nf; ++j) if (sel[j]) {
      double sum = 0.0; int cnt = 0;
      const int lo = std::max(0, j - w), hi = std::min(nf - 1, j + w);
      for (int i = lo; i <= hi; ++i) if (i != j) {
        const double r2 = ld_rsqr_hom(G, f[j], f[i], mir[j], mir[i], ntaxa);
        if (!std::isnan(r2)) { sum += r2; ++cnt; }
      }
      const double avg = cnt > 0 ? sum / cnt : NA_REAL;
      if (avg < min_r2) sel[j] = 0;                     // NaN < min_r2 is false -> not rejected
    }
  }

  for (int j = 0; j < nf; ++j) if (sel[j]) out[f[j]] = true;
  return out;
}

// ==== Stage 2b: BiparentalHaplotypeFinder.assignHaplotyes (faithful port) =====
// Reconstructs two parental haplotypes per window and chains them across a
// bidirectional window scan, writing per-site alleleA/alleleC (g-frame: 0=REF-hom,
// 2=ALT-hom, 3=NN). Operates on the ALREADY preFiltered genotype matrix.
// window=100, overlap=25, startIncr=75, minClusterSize=3, maxDifferenceScore=0.
static const int BHF_WINDOW = 100, BHF_OVERLAP = 25, BHF_MINCLUST = 3;

// HaplotypeCluster.getCensoredMajorityHaplotype(maxMaf, maxMinorCount): per site the
// majority g, unless (minorFreq > maxMaf AND minorCount > maxMinorCount) -> N(3).
static std::vector<int> cluster_censored_majority(const HCluster& c,
    const std::vector<std::vector<int> >& pool, int nsites, double maxMaf, int maxMinorCount) {
  std::vector<int> h(nsites, 3);
  for (int s = 0; s < nsites; ++s) {
    int cnt[3] = {0,0,0};
    for (size_t i = 0; i < c.members.size(); ++i) { int g = pool[c.members[i]][s]; if (g>=0&&g<=2) ++cnt[g]; }
    int total = cnt[0]+cnt[1]+cnt[2];
    if (total == 0) continue;
    int maxc = 0, maxv = -1; for (int a=0;a<3;++a) if (cnt[a] > maxc) { maxc = cnt[a]; maxv = a; }
    int minorcount = total - maxc;
    double maf = (double) minorcount / (double) total;
    h[s] = (maf <= maxMaf || minorcount <= maxMinorCount) ? maxv : 3;
  }
  return h;
}

// doesOverlapMatch: compare overlap-length ends; N-tolerant; match iff < 2 mismatches.
static bool does_overlap_match(const std::vector<int>& a, const std::vector<int>& b, int overlap, bool forward) {
  const int la = (int)a.size(), lb = (int)b.size(); int mism = 0;
  for (int i = 0; i < overlap; ++i) {
    const int av = forward ? a[la-overlap+i] : a[i];
    const int bv = forward ? b[i] : b[lb-overlap+i];
    if (av != 3 && bv != 3 && av != bv) { if (++mism >= 2) return false; }
  }
  return mism < 2;
}
// distance over the overlap region (hapstart vs parentEnd), seq_dist metric
static int overlap_dist(const std::vector<int>& hap, const std::vector<int>& ph, int overlap, bool forward) {
  const int lh = (int)hap.size(), lp = (int)ph.size(); int d = 0;
  for (int i = 0; i < overlap; ++i) {
    const int hv = forward ? hap[i] : hap[lh-overlap+i];
    const int pv = forward ? ph[lp-overlap+i] : ph[i];
    if (hv == pv || hv == 3 || pv == 3) continue;
    const bool h=(hv==1), p=(pv==1);
    if (h) { if(!p) ++d; } else if (p) ++d; else d += 2;
  }
  return d;
}

struct WinClust { std::vector<HCluster> clusters; std::vector<std::vector<int> > pool; };
// clusterWindow + sortClusters + moveAllHaplotypesToBiggestCluster + removeHeterozygousClusters
static WinClust cluster_window(const IntegerMatrix& Gf, int start, int len,
                               int maxdif, int minNotMissing, int maxHet) {
  const int ntaxa = Gf.nrow();
  WinClust wc;
  for (int t = 0; t < ntaxa; ++t) {
    std::vector<int> seq(len); int nm = 0;
    for (int p = 0; p < len; ++p) { const int g = Gf(t, start+p); seq[p] = g; if (g>=0&&g<=2) ++nm; }
    if (nm >= minNotMissing) wc.pool.push_back(seq);
  }
  if (wc.pool.empty()) return wc;
  const int np = (int) wc.pool.size();
  wc.clusters = make_clusters(wc.pool);
  if (maxdif > 0) merge_clusters(wc.clusters, maxdif, wc.pool, np);
  std::stable_sort(wc.clusters.begin(), wc.clusters.end(), cluster_less);
  move_to_biggest(wc.clusters, maxdif, wc.pool, len, np);
  if (maxHet >= 0) remove_het_clusters(wc.clusters, maxHet, wc.pool, len);
  return wc;
}

// mergeMajorHaplotypes: merge similar clusters (<=4 censored-majority distance),
// return the censored-majority haplotypes of the surviving big-enough clusters.
static std::vector<std::vector<int> > merge_major_haplotypes(WinClust& wc, int len, int minClusterSize) {
  const double maxMaf = 0.2; const int maxMinorCount = 2, maxDistance = 4;
  std::vector<HCluster>& cl = wc.clusters;
  int cc = 1;
  while (cc < (int)cl.size() && (int)cl[cc].members.size() >= minClusterSize) {
    std::vector<int> comp = cluster_censored_majority(cl[cc], wc.pool, len, maxMaf, maxMinorCount);
    for (int i = 0; i < cc; ++i) {
      std::vector<int> head = cluster_censored_majority(cl[i], wc.pool, len, maxMaf, maxMinorCount);
      if (seq_dist(head, comp) <= maxDistance) {
        merge_two(cl[i], cl[cc]); cl.erase(cl.begin()+cc); --cc; break;
      }
    }
    ++cc;
  }
  std::vector<std::vector<int> > out;
  for (int i = 0; i < cc && i < (int)cl.size(); ++i)
    out.push_back(cluster_censored_majority(cl[i], wc.pool, len, maxMaf, maxMinorCount));
  return out;
}

// getParentHaplotypes: assign each candidate to parent 0/1 via overlap match, else
// nearest parent by overlap distance (equidistant -> dropped).
static void get_parent_haplotypes(const std::vector<std::vector<int> > prev[2],
    const std::vector<std::vector<int> >& cands, int overlap, bool forward,
    std::vector<std::vector<int> > out[2]) {
  out[0].clear(); out[1].clear();
  const int n0 = (int)prev[0].size(), n1 = (int)prev[1].size(), maxN = std::max(n0,n1);
  for (size_t k = 0; k < cands.size(); ++k) {
    const std::vector<int>& hap = cands[k];
    bool match = false;
    for (int i = 0; i < maxN && !match; ++i) {
      if (i < n0 && does_overlap_match(prev[0][i], hap, overlap, forward)) { out[0].push_back(hap); match = true; }
      else if (i < n1 && does_overlap_match(prev[1][i], hap, overlap, forward)) { out[1].push_back(hap); match = true; }
    }
    if (!match) {
      int p0 = 1000000, p1 = 1000000;
      for (size_t i = 0; i < prev[0].size(); ++i) p0 = std::min(p0, overlap_dist(hap, prev[0][i], overlap, forward));
      for (size_t i = 0; i < prev[1].size(); ++i) p1 = std::min(p1, overlap_dist(hap, prev[1][i], overlap, forward));
      if (p0 < p1) out[0].push_back(hap); else if (p0 > p1) out[1].push_back(hap);  // tie -> dropped
    }
  }
}

// updatePopulationDataAlleles: fill alleleA/alleleC (g-frame 0/2/3) over the window's
// non-overlap positions from the parent haplotypes.
static void update_pop_data_alleles(const std::vector<std::vector<int> > phap[2],
    IntegerVector& alleleA, IntegerVector& alleleC, int winStart, int ptrStart, int length, int nf) {
  const int nhap0 = (int)phap[0].size(), nhap1 = (int)phap[1].size();
  if (nhap0 == 0 || nhap1 == 0) return;
  if (nhap0 == 1 && nhap1 == 1) {
    const std::vector<int>& seqA = phap[0][0]; const std::vector<int>& seqC = phap[1][0];
    for (int ptr = ptrStart; ptr < ptrStart + length; ++ptr) {
      const int ap = winStart + ptr; if (ap < 0 || ap >= nf) continue;
      const int a = seqA[ptr], c = seqC[ptr];
      if (a == c) { alleleA[ap] = 3; alleleC[ap] = 3; }
      else { alleleA[ap] = (a == 1) ? 3 : a; alleleC[ap] = (c == 1) ? 3 : c; }
    }
  } else {
    for (int ptr = ptrStart; ptr < ptrStart + length; ++ptr) {
      const int ap = winStart + ptr; if (ap < 0 || ap >= nf) continue;
      bool p0[2] = {false,false}, p1[2] = {false,false};   // alleles present: [REF, ALT]
      for (int h = 0; h < nhap0; ++h) { int v = phap[0][h][ptr]; if (v==0){p0[0]=true;} else if (v==2){p0[1]=true;} else if (v==1){p0[0]=p0[1]=true;} }
      for (int h = 0; h < nhap1; ++h) { int v = phap[1][h][ptr]; if (v==0){p1[0]=true;} else if (v==2){p1[1]=true;} else if (v==1){p1[0]=p1[1]=true;} }
      const int n0 = p0[0]+p0[1], n1 = p1[0]+p1[1];
      auto homOf = [](int allele){ return allele == 0 ? 0 : 2; };  // REF->0-hom, ALT->2-hom
      int A = 3, C = 3;
      if (n0 == 0) { A = 3; C = (n1 == 1) ? homOf(p1[0]?0:1) : 3; }
      else if (n0 == 1) {
        int Aallele = p0[0] ? 0 : 1;
        if (n1 == 0) { A = homOf(Aallele); C = 3; }
        else if (n1 == 1) { int Callele = p1[0] ? 0 : 1;
          if (Aallele == Callele) { A = 3; C = 3; } else { A = homOf(Aallele); C = homOf(Callele); } }
        else { A = 3; C = 3; }
      } else { A = 3; C = 3; }
      alleleA[ap] = A; alleleC[ap] = C;
    }
  }
}

//' FSFHap stage 2b: BiparentalHaplotypeFinder.assignHaplotyes (faithful port)
//'
//' Reconstruct two parental haplotypes across a bidirectional window scan on the
//' preFiltered genotype matrix, producing per-site parent alleles. Seed: first
//' window with exactly two clusters (3rd `< minClusterSize`) whose majority
//' haplotypes differ by `>= 2*window-4`; then extend forward and backward,
//' matching each window's candidate haplotypes to the running parents and writing
//' the non-overlap alleles.
//'
//' @param Gf Integer matrix, taxa x sites, canonical `g` in {0,1,2,3}; the
//'   ALREADY preFiltered chromosome (see [fsfhap_prefilter_sites_cpp()]).
//' @return List: `alleleA`, `alleleC` (length ncol(Gf); g-frame `0`=REF-hom /
//'   `2`=ALT-hom / `3`=NN) parent-of-origin alleles per site; `seeded` (logical).
//' @keywords internal
// [[Rcpp::export]]
List fsfhap_biparental_alleles_cpp(IntegerMatrix Gf) {
  const int nf = Gf.ncol();
  IntegerVector alleleA(nf, 3), alleleC(nf, 3);
  const int window = BHF_WINDOW, overlap = BHF_OVERLAP, startIncr = window - overlap;
  const int diff = 0;   // maxDifferenceScore
  if (nf < window) return List::create(_["alleleA"]=alleleA, _["alleleC"]=alleleC, _["seeded"]=false);

  // ---- seed: first window with two well-separated clusters --------------------
  int initialStart = 0; const int maxStart = nf - window; bool seeded = false;
  std::vector<int> h0, h1;
  const int minNotMissing = (int)(window * 0.2);
  // strict `<` is faithful to TASSEL (`initialStart < maxStart`): a chromosome with
  // exactly one window (nf == window) is not seeded, matching the reference.
  while (!seeded && initialStart < maxStart) {
    WinClust wc = cluster_window(Gf, initialStart, window, diff, minNotMissing, 5 + diff);
    if ((int)wc.clusters.size() >= 2) {
      h0 = cluster_majority(wc.clusters[0], wc.pool, window);
      h1 = cluster_majority(wc.clusters[1], wc.pool, window);
      int c2 = (int)wc.clusters.size() > 2 ? (int)wc.clusters[2].members.size() : 0;
      if (c2 < BHF_MINCLUST && seq_dist(h0, h1) >= 2 * window - 4) { seeded = true; break; }
    }
    initialStart += window;
  }
  if (!seeded) return List::create(_["alleleA"]=alleleA, _["alleleC"]=alleleC, _["seeded"]=false);

  std::vector<std::vector<int> > parent[2];
  parent[0].push_back(h0); parent[1].push_back(h1);
  // write the SEED window's alleles (ptrStart=0, length=window) before extending
  update_pop_data_alleles(parent, alleleA, alleleC, initialStart, 0, window, nf);

  // ---- forward extension ------------------------------------------------------
  for (int start = initialStart + startIncr; start < nf - overlap; start += startIncr) {
    int windowSize = window; if (start + window > nf) windowSize = nf - start;
    const int mnm = (int)(windowSize * 0.2);
    WinClust wc = cluster_window(Gf, start, windowSize, diff, mnm, 5 + diff);
    std::vector<std::vector<int> > cands = merge_major_haplotypes(wc, windowSize, BHF_MINCLUST);
    std::vector<std::vector<int> > np[2]; get_parent_haplotypes(parent, cands, overlap, true, np);
    // CARRY-FORWARD GUARD (extension beyond faithful TASSEL): a degenerate window
    // that would leave a parent list EMPTY (het-heavy F2 → single candidate) is
    // skipped, retaining the previous parents so neither can vanish and trigger the
    // one-sided cascade. No-op when both parents are filled → RIL/BC1 paths (and
    // TASSEL parity there) are unchanged. See design/FSFHAP_PORT.md (task #5).
    if (np[0].empty() || np[1].empty()) continue;
    parent[0] = np[0]; parent[1] = np[1];
    update_pop_data_alleles(parent, alleleA, alleleC, start, overlap, windowSize - overlap, nf);
  }

  // ---- backward extension (reset to seed) -------------------------------------
  parent[0].assign(1, h0); parent[1].assign(1, h1);
  for (int start = initialStart - startIncr; start > -startIncr; start -= startIncr) {
    int s0 = start, windowSize = window;
    if (s0 < 0) { windowSize = window + s0; s0 = 0; }        // end = start+window fixed; start clamped
    if (windowSize <= overlap) continue;
    const int mnm = (int)(windowSize * 0.2);
    WinClust wc = cluster_window(Gf, s0, windowSize, diff, mnm, 5 + diff);
    std::vector<std::vector<int> > cands = merge_major_haplotypes(wc, windowSize, BHF_MINCLUST);
    std::vector<std::vector<int> > np[2]; get_parent_haplotypes(parent, cands, overlap, false, np);
    if (np[0].empty() || np[1].empty()) continue;             // carry-forward guard (see forward loop)
    parent[0] = np[0]; parent[1] = np[1];
    update_pop_data_alleles(parent, alleleA, alleleC, s0, 0, windowSize - overlap, nf);
  }
  return List::create(_["alleleA"]=alleleA, _["alleleC"]=alleleC, _["seeded"]=true);
}
