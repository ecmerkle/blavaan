### exported entry point for posterior predictive p-value (ppp).
###
### For single-level ordinal/mixed target="stan"/"cmdstan" models, ppp is
### no longer computed by default (see R/blavaan.R) because the
### replacement statistic (R/blav_limited_info.R) is real post-hoc R
### computation, not free like the old Stan-side value. This is the
### "second step" the fit-time NOTE points users to: it returns the
### already-computed ppp when one exists, and computes it on demand
### otherwise.

ppp <- function(object, ...) {
  lavopt <- lavInspect(object, "options")

  if (lavopt$target %in% c("stan", "cmdstan")) {
    if (blavInspect(object, "nlevels") > 1) {
      if (lavopt$test == "none") {
        stop("blavaan ERROR: on-demand ppp computation for two-level models is ",
             "not yet supported; refit with test != \"none\" if you need ppp.")
      }
      return(fitMeasures(object, "ppp"))
    }

    if (lavopt$test != "none") {
      return(fitMeasures(object, "ppp"))
    }

    ## test == "none": nothing was computed at fit time. test="none" is
    ## only the *implicit* default for ordinal/mixed models -- for
    ## continuous-only models it only happens by explicit user choice,
    ## and pp_limited_info() was never scoped to handle that case
    ## (continuous-only stays on the old Stan mechanism by design).
    if (lavInspect(object, "categorical")) {
      return(pp_limited_info(object, ...)$ppval)
    } else {
      ## no Stan-side value was ever computed to fall back on (do_test
      ## was 0 during sampling), and object@external$samplls was never
      ## computed either (R/blavaan.R only computes it when test !=
      ## "none"), so postpred()'s default measure="logl" fast path
      ## (which indexes into samplls) is unusable here. Match ppmc()'s
      ## own precedent (R/ctr_ppmc.R) of requesting a measure that
      ## forces postpred()'s self-contained get_ll()-based path instead.
      return(postpred(object@external$samplls, object,
                      measure = c("chisq", "chisq.scaled"))$ppval)
    }
  } else {
    return(postpred(object@external$samplls, object)$ppval)
  }
}
