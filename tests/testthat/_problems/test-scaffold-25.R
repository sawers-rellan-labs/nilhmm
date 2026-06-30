# Extracted from test-scaffold.R:25

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "nilHMM", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_error(fit(NULL, emission_count(), duration_geometric(), list()),
               "not yet implemented")
