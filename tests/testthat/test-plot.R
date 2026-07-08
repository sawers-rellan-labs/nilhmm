test_that("paint_calls returns a ggplot and validates its input", {
  skip_if_not_installed("ggplot2")
  calls <- data.frame(
    name = c("A", "A", "B"), donor = "B", chr = 1L,
    start_bp = c(1e6, 4e6, 1e6), end_bp = c(4e6, 9e6, 9e6),
    state = c(0L, 2L, 0L), method = c("m1", "m1", "m2"),
    stringsAsFactors = FALSE)

  expect_s3_class(paint_calls(calls), "ggplot")                 # one band per sample
  expect_s3_class(paint_calls(calls, track = "method"), "ggplot")  # stacked tracks
  expect_s3_class(paint_calls(calls, samples = "A"), "ggplot")

  # a factor track is accepted (its level order sets the band order)
  calls$method <- factor(calls$method, levels = c("m2", "m1"))
  expect_s3_class(paint_calls(calls, track = "method"), "ggplot")

  expect_error(paint_calls(calls[, c("name", "chr")]), "needs columns")
  expect_error(paint_calls(calls, track = "nope"), "not in")
  expect_error(paint_calls(calls[0, ]), "no segments")
  expect_error(paint_calls(calls, palette = c("a", "b")), "length-3")
})
