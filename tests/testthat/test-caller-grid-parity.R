# caller_grid must reproduce a cold call_ancestry() at every (emission, rrate)
# grid point -- the "exact per grid point" promise. Small synthetic nnil cohort.

test_that("caller_grid(nnil) equals call_ancestry per grid point", {
  set.seed(7)
  mk <- do.call(rbind, lapply(1:2, function(ch) {
    data.frame(chr = ch, pos = sort(sample.int(1e7, 120)))
  }))
  data <- do.call(rbind, lapply(paste0("S", 1:3), function(nm) {
    d <- mk
    d$name <- nm
    d$g <- sample(0:3, nrow(d), replace = TRUE, prob = c(0.7, 0.1, 0.1, 0.1))
    d
  }))[, c("name", "chr", "pos", "g")]

  eg <- expand.grid(nir = c(0.01, 0.3, 0.9), germ = c(1e-3, 1e-2))
  rr <- c(2e-4, 6e-4)
  grid <- as.data.frame(caller_grid(data,
    caller = "nnil",
    emission_grid = eg, rrate = rr, design = "BC2S3"
  ))

  seg_cols <- c("source", "donor", "name", "chr", "start_bp", "end_bp", "state")
  ord <- function(x) {
    x <- x[, seg_cols, drop = FALSE]
    x[order(x$name, x$chr, x$start_bp), , drop = FALSE]
  }
  for (i in seq_len(nrow(eg))) {
    for (r in rr) {
      ref <- as.data.frame(call_ancestry(data,
        caller = "nnil", design = "BC2S3",
        rrate = r, germ = eg$germ[i], nir = eg$nir[i]
      ))
      sub <- grid[grid$nir == eg$nir[i] & grid$germ == eg$germ[i] & grid$rrate == r, , drop = FALSE]
      expect_equal(ord(sub), ord(ref), ignore_attr = TRUE)
    }
  }
})
