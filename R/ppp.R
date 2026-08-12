### exported entry point for posterior predictive p-value (ppp).
###
### For single-level ordinal/mixed and two-level (multilevel)
### target="stan"/"cmdstan" models, ppp is no longer computed by default
### (see R/blavaan.R) because the replacement statistics
### (R/blav_limited_info.R, R/blav_twolevel_ppp.R) are real post-hoc R
### computation, not free like the old Stan-side value. This is the
### "second step" the fit-time NOTE points users to: it returns the
### already-computed ppp when one exists, and computes it on demand
### otherwise.

ppp <- function(object, ...) {
  lavopt <- lavInspect(object, "options")

  if (lavopt$target %in% c("stan", "cmdstan")) {
    if (blavInspect(object, "nlevels") > 1) {
      if (lavopt$test != "none") {
        ## already computed at fit time via R/blav_test.R's pp_twolevel()
        ## dispatch
        return(fitMeasures(object, "ppp"))
      }
      return(pp_twolevel(object, ...)$ppval)
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
