# Callers = named (emission x duration) combinations (S2, S8). "nilHMM" is the
# PACKAGE; the callers are explicit methods inside it. This kills the
# package-vs-caller ambiguity.
#
# The named cells of the grid. Two things vary: the emission x duration axes, and
# (on the gt row) the hard-call thresholding lineage -- Holland's 1/3-2/3 dosage
# cutoffs vs GOOGA's competitive-alignment fraction thresholds:
#
#                              geometric        rigidity
#   gt, Holland thresholds     nnil             catiger
#   gt, GOOGA thresholds       googa            atlas
#   count / BetaBinomial       bbnil            rtiger
#
# `googa` is the faithful GOOGA/Veltsos reproduction (categorical + geometric =
# recombination-fraction F2 HMM); `atlas` is this work's transcript caller (the
# same GOOGA thresholding with the rigidity duration). `nnil`/`googa` share the
# emission x duration spec (gt + geometric) and differ only in the hard-call rule;
# likewise `catiger`/`atlas` (gt + rigidity). The GOOGA thresholding itself is
# applied in [call_ancestry()], not here.

#' Resolve a named caller into emission + duration specs
#'
#' The named cells of the shared engine's grid (emission x duration; the gt row
#' also varies by hard-call thresholding lineage):
#' - `nnil`    : categorical `gt` emission + geometric duration -- Holland's
#'   original NIL caller on hard genotype calls (1/3-2/3 dosage cutoffs).
#' - `bbnil`   : count/BetaBinomial emission + geometric duration -- the
#'   low-coverage read-count extension of nNIL; the self-transition is the smoother.
#' - `catiger` : categorical `gt` emission + rigidity duration (S7) -- the
#'   Holland-threshold categorical + minimum-run-length caller.
#' - `rtiger`  : count/BetaBinomial emission + rigidity duration (S7).
#' - `googa`   : categorical `gt` emission + geometric duration, with GOOGA
#'   competitive-alignment thresholds -- the faithful GOOGA/Veltsos reproduction.
#' - `atlas`   : categorical `gt` emission + rigidity duration, with GOOGA
#'   thresholds -- this work's transcript caller (rigidity variant of `googa`).
#'   The GOOGA thresholding is applied in [call_ancestry()].
#'
#' @param caller One of `"nnil"`, `"bbnil"`, `"catiger"`, `"rtiger"`, `"googa"`,
#'   `"atlas"`.
#' @param rrate Geometric callers (`nnil`/`bbnil`/`googa`): expected per-marker
#'   recombination rate (self-stay = `1 - rrate`). Holland's nNIL sets it to
#'   `2 * total_cM / (100 * n_markers)`.
#' @param rigidity Rigidity callers (`catiger`/`rtiger`/`atlas`): integer minimum
#'   run length (e.g. `5`).
#' @param err Count-emission baseline error.
#' @param conc Count-emission BetaBinomial concentration.
#' @param fit_means EM-fit emission means (count emission; S10).
#' @param xrate Exit rate of the **rigidity duration** ([duration_rigidity()]):
#'   free-state (post-minimum-run) switch probability. A nilHMM construct, not a
#'   RTIGER parameter.
#' @param germ,gert,p,mr,nir Genotype-error rates for the gt emission (the gt
#'   callers `nnil`/`catiger`/`googa`/`atlas`; Holland's nNIL error model): hom
#'   error, het error, hom-error->het fraction, missing rate,
#'   non-informative-marker rate.
#' @param ... Ignored extra args (e.g. `f_1`/`f_2` consumed by [call_ancestry()]).
#' @return `list(emission, duration)`.
#' @examples
#' caller_spec("nnil", rrate = 1e-4)      # gt emission + geometric duration
#' caller_spec("bbnil", rrate = 1e-4)     # count emission + geometric duration
#' caller_spec("catiger", rigidity = 5)   # gt emission + rigidity duration
#' caller_spec("rtiger", rigidity = 5)    # count emission + rigidity duration
#' str(caller_spec("googa"))              # gt emission + geometric (faithful GOOGA)
#' str(caller_spec("atlas", rigidity = 5))# gt emission + rigidity (this-work transcript)
#' @export
caller_spec <- function(caller = c("nnil", "bbnil", "catiger", "rtiger", "googa", "atlas"),
                        rrate = 0.01, rigidity = NULL, err = 0.01, conc = 20,
                        fit_means = FALSE, xrate = 0.01,
                        germ = 0.05, gert = 0.10, p = 0.5,
                        mr = 0.10, nir = 0.01, ...) {
  caller <- match.arg(caller)
  rig <- if (is.null(rigidity)) 5L else as.integer(rigidity)   # minimum run length
  switch(caller,
    nnil    = list(emission = emission_gt(germ, gert, p, mr, nir),
                   duration = duration_geometric(rrate)),
    googa   = list(emission = emission_gt(germ, gert, p, mr, nir),
                   duration = duration_geometric(rrate)),
    bbnil   = list(emission = emission_count(err, conc, fit_means),
                   duration = duration_geometric(rrate)),
    catiger = list(emission = emission_gt(germ, gert, p, mr, nir),
                   duration = duration_rigidity(rig, xrate)),
    atlas   = list(emission = emission_gt(germ, gert, p, mr, nir),
                   duration = duration_rigidity(rig, xrate)),
    rtiger  = list(emission = emission_count(err, conc, fit_means),
                   duration = duration_rigidity(rig, xrate))
  )
}
