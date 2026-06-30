# data-raw/ holds the SCRIPTS that regenerate the bundled `breeding_designs`
# object; the raw simulation CSV is NOT vendored -- the seed + recipe regenerate
# it (REFACTOR_R_PACKAGE.md §6). The stored data/ value is authoritative at
# runtime; this script is for verification/extension.
#
# Pin the seed AND the simcross/R versions + RNGkind() -- a seed reproduces the
# exact shipped (k, lambda) only under a fixed environment.
#
# Recipe per design (BC1S1, BC2S2, BC2S3, ...):
#   simcross from pinned (pedigree, map cM lengths, Stahl m = 10, p = 0,
#   n = 1500) -> segment sizes (cM) -> fit_design_gamma() -> (k, lambda),
#   plus g, f_1, f_2, mean cM, provenance. Then usethis::use_data().
#
# Confirm the BC2S3 values surfaced by the BRB run (§10): g, f_1 = 0.0312,
# f_2 = 0.1094, and the fitted BC2S3 Gamma.

stop("data-raw/make_breeding_designs.R is a scaffold stub (Task 4 / §6).")