// Pedigree-aware ancestry refinement: structured loopy belief propagation over
// one family's (pedigree x genome) grid. Design of record: design/PEDIGREE_HMM.md.
//
// Each individual's chromosome chain is solved exactly (tilted forward-backward);
// parent-child selfing edges are the loopy part, handled by damped cavity message
// passing. Exact-within-chain-FB, converged-message-passing, linear-cost
// approximation. R owns all IO and builds the emission matrices; this kernel is
// data-agnostic and processes ONE (family, chromosome) per call (markers sorted).

#include <Rcpp.h>
#include <vector>
#include <array>
#include <cmath>
#include <functional>
using namespace Rcpp;

namespace {

const int    K      = 3;
const double PFLOOR = 1e-12;                 // cavity-division stability (sec 14)
typedef std::array<double, K> Col;
typedef std::vector<Col>      Mat;           // M markers x K

// Selfing transmission kernel (row = parent dosage, col = child dosage).
const double Tsib[K][K] = {
  {1.00, 0.00, 0.00},
  {0.25, 0.50, 0.25},
  {0.00, 0.00, 1.00}
};

inline void pnormalize(Col& c) {
  double s = c[0] + c[1] + c[2];
  if (s <= 0) { c[0] = c[1] = c[2] = 1.0 / 3; return; }  // deliberate guard (sec 5/14)
  c[0] /= s; c[1] /= s; c[2] /= s;
}

// Relax-to-stationary dosage transition: A[x][y] = (1-r_eff) d_xy + r_eff pi[y].
// Stationary == pi exactly. r_eff is the meiosis-compounded recombination
// fraction, P(odd # crossovers) over `meioses` meioses, capped at 0.5.
inline void makeA(double A[K][K], double r_interval, int meioses, const Col& pi) {
  int mm = meioses < 1 ? 1 : meioses;
  double r = 0.5 * (1.0 - std::pow(1.0 - 2.0 * r_interval, mm));
  for (int x = 0; x < K; ++x)
    for (int y = 0; y < K; ++y)
      A[x][y] = (x == y ? 1.0 - r : 0.0) + r * pi[y];
}

struct Node {
  int parent = -1;
  std::vector<int> kids;
  int  meioses = 1;
  bool hasData = false;
  Mat  emit;                   // M x K (only when hasData)
  Col  rho;                    // marker-0 prior
  Col  pi;                     // generation stationary
  Mat  upMsg;                  // v -> parent   (parent states)
  std::vector<Mat> downMsg;    // v -> each kid (kid states)
  Mat  belief;                 // gamma_v
};

struct Fam {
  std::vector<Node> nodes;
  int root = 0;
  int M = 0;
  std::vector<double> r;       // per-interval recomb fraction, size M-1
};

// Tilted forward-backward: effective emission = emit (or 1 for latent) x tilt,
// tilt = product of incoming messages. Column-wise renormalized; returns gamma.
Mat tiltedFB(const Node& nd, const Fam& fam, const Mat& tilt) {
  const int M = fam.M;
  Mat alpha(M), beta(M), gamma(M);
  double A[K][K];

  for (int x = 0; x < K; ++x) {
    double e = nd.hasData ? nd.emit[0][x] : 1.0;
    alpha[0][x] = nd.rho[x] * e * tilt[0][x];
  }
  pnormalize(alpha[0]);
  for (int i = 1; i < M; ++i) {
    makeA(A, fam.r[i - 1], nd.meioses, nd.pi);
    for (int x = 0; x < K; ++x) {
      double s = 0;
      for (int y = 0; y < K; ++y) s += alpha[i - 1][y] * A[y][x];
      double e = nd.hasData ? nd.emit[i][x] : 1.0;
      alpha[i][x] = s * e * tilt[i][x];
    }
    pnormalize(alpha[i]);
  }
  beta[M - 1] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
  for (int i = M - 2; i >= 0; --i) {
    makeA(A, fam.r[i], nd.meioses, nd.pi);
    for (int x = 0; x < K; ++x) {
      double s = 0;
      for (int y = 0; y < K; ++y) {
        double e = nd.hasData ? nd.emit[i + 1][y] : 1.0;
        s += A[x][y] * e * tilt[i + 1][y] * beta[i + 1][y];
      }
      beta[i][x] = s;
    }
    pnormalize(beta[i]);
  }
  for (int i = 0; i < M; ++i) {
    for (int x = 0; x < K; ++x) gamma[i][x] = alpha[i][x] * beta[i][x];
    pnormalize(gamma[i]);
  }
  return gamma;
}

// Product of all incoming messages at node v (children upMsg + parent downMsg).
Mat incoming(const Fam& fam, int v) {
  const Node& nd = fam.nodes[v];
  Mat prod(fam.M, Col{1, 1, 1});
  for (int child : nd.kids) {
    const Mat& mchild = fam.nodes[child].upMsg;
    for (int i = 0; i < fam.M; ++i)
      for (int x = 0; x < K; ++x) prod[i][x] *= mchild[i][x];
  }
  if (nd.parent >= 0) {
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
  for (int x = 0; x < K; ++x)
    c[x] = belief[x] / (excluded[x] > PFLOOR ? excluded[x] : PFLOOR);
  pnormalize(c);
  return c;
}

// Damped geometric update cur <- cur^(1-l) * tgt^l; returns max abs change.
double dampInto(Mat& cur, const Mat& tgt, double lambda) {
  double maxd = 0;
  for (size_t i = 0; i < cur.size(); ++i) {
    Col nx;
    for (int x = 0; x < K; ++x)
      nx[x] = std::pow(cur[i][x] > PFLOOR ? cur[i][x] : PFLOOR, 1 - lambda)
            * std::pow(tgt[i][x] > PFLOOR ? tgt[i][x] : PFLOOR, lambda);
    pnormalize(nx);
    for (int x = 0; x < K; ++x) maxd = std::max(maxd, std::fabs(nx[x] - cur[i][x]));
    cur[i] = nx;
  }
  return maxd;
}

// child v -> parent: contract Tsib over child states (exclude parent's message).
double sendUp(Fam& fam, int v, double lambda) {
  Node& nd = fam.nodes[v];
  Node& p  = fam.nodes[nd.parent];
  int slot = 0; while (p.kids[slot] != v) ++slot;
  Mat tgt(fam.M);
  for (int i = 0; i < fam.M; ++i) {
    Col cav = cavity(nd.belief[i], p.downMsg[slot][i]);
    Col out{0, 0, 0};
    for (int xp = 0; xp < K; ++xp)
      for (int xc = 0; xc < K; ++xc) out[xp] += Tsib[xp][xc] * cav[xc];
    pnormalize(out); tgt[i] = out;
  }
  return dampInto(nd.upMsg, tgt, lambda);
}

// parent v -> child k: contract Tsib over parent states (exclude child k's message).
double sendDown(Fam& fam, int v, int k, double lambda) {
  Node& nd = fam.nodes[v];
  int child = nd.kids[k];
  Mat tgt(fam.M);
  for (int i = 0; i < fam.M; ++i) {
    Col cav = cavity(nd.belief[i], fam.nodes[child].upMsg[i]);
    Col out{0, 0, 0};
    for (int xc = 0; xc < K; ++xc)
      for (int xp = 0; xp < K; ++xp) out[xc] += Tsib[xp][xc] * cav[xp];
    pnormalize(out); tgt[i] = out;
  }
  return dampInto(nd.downMsg[k], tgt, lambda);
}

int runFamilyBP(Fam& fam, int maxIters, double tol, double lambda) {
  // post-order (children before parents) and its reverse (pre-order).
  std::vector<int> post;
  std::function<void(int)> dfs = [&](int v) {
    for (int c : fam.nodes[v].kids) dfs(c);
    post.push_back(v);
  };
  dfs(fam.root);
  std::vector<int> pre(post.rbegin(), post.rend());

  // Initialize ALL messages to uniform M x K before any sweep: incoming()/sendUp()
  // read parent downMsg (and children upMsg) before a send populates them.
  Mat uni(fam.M, Col{1.0 / 3, 1.0 / 3, 1.0 / 3});
  for (Node& nd : fam.nodes) {
    nd.upMsg = uni;
    nd.downMsg.assign(nd.kids.size(), uni);
  }
  // init beliefs: independent HMM per node (flat tilt)
  Mat flat(fam.M, Col{1, 1, 1});
  for (int v : post) fam.nodes[v].belief = tiltedFB(fam.nodes[v], fam, flat);

  int it = 0;
  for (; it < maxIters; ++it) {
    double maxd = 0;
    for (int v : post) {                        // upward sweep
      fam.nodes[v].belief = tiltedFB(fam.nodes[v], fam, incoming(fam, v));
      if (fam.nodes[v].parent >= 0)
        maxd = std::max(maxd, sendUp(fam, v, lambda));
    }
    for (int v : pre) {                          // downward sweep
      fam.nodes[v].belief = tiltedFB(fam.nodes[v], fam, incoming(fam, v));
      for (int k = 0; k < (int)fam.nodes[v].kids.size(); ++k)
        maxd = std::max(maxd, sendDown(fam, v, k, lambda));
    }
    if (maxd < tol) { ++it; break; }
  }
  for (int v : pre)                              // final beliefs
    fam.nodes[v].belief = tiltedFB(fam.nodes[v], fam, incoming(fam, v));
  return it;
}

}  // namespace

//' Loopy belief propagation over one family's pedigree x genome grid
//'
//' Data-agnostic BP kernel for [refine_ancestry()]; processes ONE (family,
//' chromosome). Node fields are pre-built in R. See design/PEDIGREE_HMM.md.
//'
//' @param M Number of markers (chromosome length).
//' @param parent 0-based parent index per node (`-1` for the root).
//' @param meioses Per-node accumulated meiosis count (transition block length).
//' @param hasData Per-node logical; `TRUE` for genotyped leaves.
//' @param emit Length-V list; for a `hasData` node an `M x 3` emission matrix
//'   (P(obs | state)), otherwise ignored (latent nodes emit 1).
//' @param rho V x 3 marker-0 prior per node.
//' @param pimat V x 3 generation stationary per node (transition relaxes to it).
//' @param r Length `M-1` per-interval recombination fraction (base; meiosis
//'   compounding is applied per node).
//' @param root 0-based index of the family root (founder).
//' @param maxIters,tol,lambda Message-passing iterations, convergence tolerance,
//'   damping.
//' @return Length-V list of `M x 3` posterior belief matrices; attribute `iters`
//'   records sweeps run.
//' @keywords internal
// [[Rcpp::export]]
List pedigree_bp_cpp(int M, IntegerVector parent, IntegerVector meioses,
                     LogicalVector hasData, List emit,
                     NumericMatrix rho, NumericMatrix pimat,
                     NumericVector r, int root,
                     int maxIters, double tol, double lambda) {
  const int V = parent.size();
  if (M < 1) stop("pedigree_bp_cpp(): M must be >= 1");
  if ((int)r.size() != M - 1) stop("pedigree_bp_cpp(): length(r) must be M-1");
  if (rho.nrow() != V || pimat.nrow() != V)
    stop("pedigree_bp_cpp(): rho/pimat must have one row per node");

  Fam fam;
  fam.root = root;
  fam.M = M;
  fam.r.assign(r.begin(), r.end());
  fam.nodes.resize(V);
  for (int v = 0; v < V; ++v) {
    Node& nd = fam.nodes[v];
    nd.parent  = parent[v];
    nd.meioses = meioses[v];
    nd.hasData = hasData[v];
    for (int x = 0; x < K; ++x) { nd.rho[x] = rho(v, x); nd.pi[x] = pimat(v, x); }
    if (nd.parent >= 0) fam.nodes[nd.parent].kids.push_back(v);
    if (nd.hasData) {
      NumericMatrix e = emit[v];
      if (e.nrow() != M || e.ncol() != K)
        stop("pedigree_bp_cpp(): emit[[%d]] must be M x 3", v + 1);
      nd.emit.resize(M);
      for (int i = 0; i < M; ++i)
        for (int x = 0; x < K; ++x) nd.emit[i][x] = e(i, x);
    }
  }

  int iters = runFamilyBP(fam, maxIters, tol, lambda);

  List out(V);
  for (int v = 0; v < V; ++v) {
    NumericMatrix b(M, K);
    for (int i = 0; i < M; ++i)
      for (int x = 0; x < K; ++x) b(i, x) = fam.nodes[v].belief[i][x];
    out[v] = b;
  }
  out.attr("iters") = iters;
  return out;
}
