# Extracted from test-scaffold.R:29

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "nilHMM", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_error(caller_spec("skimbin"), NA)
d <- data.frame(name = "s", chr = 1L, pos = 1:3, n_ref = c(1L,0L,1L),
                  n_alt = c(0L,1L,0L), donor = "Zx")
expect_error(call_ancestry(d, caller = "skimbin", design = "BC2S2", r = 0.01),
               "only the count emission")
