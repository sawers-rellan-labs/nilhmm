# Plotting (S8). Fragment-size distributions against the grey Null (the design's
# expected Gamma from expected_fragment_dist() -- calibration/Null only, S7),
# and chromosome-painting of called segments.

# aes() below maps precomputed helper columns by bare name; declare them so R CMD
# check's global-variable analysis stays quiet.
utils::globalVariables(c(".xmb0", ".xmb1", ".yb0", ".yb1", ".statef", ".chrf"))

#' Plot called fragment sizes against the expected Null
#'
#' @param calls Segment calls (common schema).
#' @param design Design key for the grey Null (see [expected_fragment_dist()]).
#' @param space One of `"cM"`, `"Mb"`.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # Planned (Task 4): called fragment sizes vs the expected Null.
#' calls <- call_ancestry(read_counts("counts/"), caller = "nnil", design = "BC2S2")
#' plot_fragment_sizes(calls, design = "BC2S2", space = "Mb")
#' }
#' @export
plot_fragment_sizes <- function(calls, design, space = c("Mb", "cM")) {
  space <- match.arg(space)
  stop("nilHMM::plot_fragment_sizes() not yet implemented (Task 4)")
}

#' Chromosome-paint called ancestry segments
#'
#' Paint the common-schema segments as coloured rectangles along each chromosome
#' --- REF / HET / ALT (`state` 0/1/2). Each **sample** (`name`) is a facet row and
#' each chromosome a facet column. When `track` is supplied, its levels are stacked
#' as horizontal bands *within* every cell, so several callers (or data sources)
#' can be overlaid on identical axes to check whether they agree on where the donor
#' blocks land --- the "sanity paint" comparison. Build the multi-track input by
#' `rbind()`-ing each caller's [call_ancestry()] output with a column naming the
#' track, then pass that column name as `track`.
#'
#' @param calls Segment calls in the common schema: columns `name`, `chr`,
#'   `start_bp`, `end_bp`, `state` (`0`/`1`/`2`). Extra columns are ignored.
#' @param track Optional column name whose levels are stacked as bands within each
#'   `name` x `chr` cell (e.g. `"method"` for a caller comparison, or `"donor"`).
#'   `NULL` (default) paints one band per sample. A factor `track` keeps its level
#'   order (top band = first level).
#' @param samples Optional character vector restricting to a subset of `name`
#'   values (keeps busy figures legible).
#' @param palette Named length-3 fill for `c(REF, HET, ALT)`.
#' @return A \pkg{ggplot2} object (requires the suggested \pkg{ggplot2}).
#' @examples
#' calls <- data.frame(
#'   name = "NIL01", donor = "B", chr = 1L,
#'   start_bp = c(1e6, 3e6, 5e6), end_bp = c(3e6, 5e6, 9e6),
#'   state = c(0L, 2L, 0L))
#' if (requireNamespace("ggplot2", quietly = TRUE)) paint_calls(calls)
#' @export
paint_calls <- function(calls, track = NULL, samples = NULL,
                        palette = c(REF = "gold", HET = "springgreen4", ALT = "purple4")) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("paint_calls() needs the 'ggplot2' package (add it to your library)")
  calls <- as.data.frame(calls, stringsAsFactors = FALSE)
  need <- c("name", "chr", "start_bp", "end_bp", "state")
  if (!all(need %in% names(calls)))
    stop("paint_calls(): `calls` needs columns ", paste(need, collapse = ", "))
  if (!is.null(track) && !track %in% names(calls))
    stop("paint_calls(): track column '", track, "' not in `calls`")
  if (length(palette) != 3L || is.null(names(palette)))
    stop("paint_calls(): `palette` must be a named length-3 vector (REF, HET, ALT)")
  if (!is.null(samples)) calls <- calls[calls$name %in% samples, , drop = FALSE]
  if (!nrow(calls)) stop("paint_calls(): no segments to plot")

  # stack the tracks as horizontal bands within each (name, chr) facet cell
  if (is.null(track)) {
    trk <- factor(rep("", nrow(calls)))
  } else {
    tv <- calls[[track]]
    trk <- if (is.factor(tv)) tv else factor(tv, levels = unique(tv))
  }
  M <- nlevels(trk); k <- M - as.integer(trk); gap <- if (M > 1L) 0.015 else 0
  calls$.yb0 <- k / M + gap
  calls$.yb1 <- (k + 1L) / M - gap
  calls$.xmb0 <- calls$start_bp / 1e6
  calls$.xmb1 <- calls$end_bp / 1e6
  calls$.statef <- factor(calls$state, levels = 0:2, labels = names(palette))
  chr_lv <- sort(unique(calls$chr))
  calls$.chrf <- factor(calls$chr, levels = chr_lv, labels = paste0("chr", chr_lv))

  p <- ggplot2::ggplot(calls, ggplot2::aes(xmin = .xmb0, xmax = .xmb1,
                                           ymin = .yb0, ymax = .yb1, fill = .statef)) +
    ggplot2::geom_rect() +
    ggplot2::facet_grid(name ~ .chrf, scales = "free_x", space = "free_x", switch = "y") +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE, name = "State") +
    ggplot2::labs(x = "Position (Mb)") +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(legend.position = "bottom",
                   strip.background = ggplot2::element_blank(),
                   strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
                   panel.spacing.x = ggplot2::unit(0.1, "lines"),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  if (M > 1L)
    p <- p + ggplot2::scale_y_continuous(
      breaks = (M - seq_len(M) + 0.5) / M, labels = levels(trk),
      position = "right", expand = ggplot2::expansion(mult = 0.01))
  else
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                            axis.ticks.y = ggplot2::element_blank())
  p
}