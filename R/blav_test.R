blav_model_test <- function(lavmodel       = NULL, 
                            lavpartable    = NULL, 
                            lavsamplestats = NULL, 
                            lavoptions     = NULL, 
                            x              = NULL, 
                            VCOV           = NULL, 
                            lavcache       = NULL,
                            lavdata        = NULL,
                            lavjags        = NULL,
                            lavobject      = NULL,
                            samplls        = NULL,
                            jagextra       = NULL,
                            stansumm       = NULL,
                            domll          = NULL,
                            control        = list()) {


  TEST <- list()

  ## marginal log-likelihood approximation
  ## needs original partable with rhos
  if("syntax" %in% names(jagextra)) {
    warning("blavaan WARNING: Marginal log-likelihood cannot be approximated when there is additional JAGS syntax.", call. = FALSE)
    mll <- NA
  } else if(domll) {
    mll <- try(margloglik(lavpartable, lavmodel, lavoptions,
                          lavsamplestats, lavdata, lavcache,
                          lavjags, VCOV, x, stansumm),
               silent=TRUE)
    if(inherits(mll, "try-error")) mll <- NA
  } else {
    mll <- NA # not tested, priors may cause problems
  }

  if(lavoptions$target %in% c("stan", "cmdstan")) {
    has_ordinal <- length(lavdata@ordered) > 0
    if (isTRUE(lavoptions$.multilevel)) {
      ## two-level: always computed in R now (R/blav_twolevel_ppp.R) -- the
      ## Stan-side EM/saturated-model machinery for this was removed (see
      ## inst/stan/stanmarg.stan and R/blavaan.R's do_test computation), so
      ## stansumm's 'ppp'/'log_lik_sat' values are no longer meaningful for
      ## multilevel fits. Reached only when the user explicitly overrides
      ## test="standard"/"ppp" on a two-level fit, since test defaults to
      ## "none" there (R/blavaan.R) -- see R/ppp.R for the on-demand
      ## second-step path.
      ppp <- pp_twolevel(lavobject)$ppval
    } else if (!has_ordinal) {
      ## continuous-only (existing mechanism already appropriate) models
      ## keep the Stan-computed value; same compiled Stan program either
      ## way, so stansumm already has the correctly-dispatched-on-
      ## missingness ppp value regardless of which rstan/cmdstanr summary
      ## path built it
      ppp <- stansumm['ppp', 'mean']
    } else {
      ## single-level ordinal/mixed: limited-information pairwise ppp,
      ## computed in R from saved posterior draws (see R/blav_limited_info.R).
      ## Reached only when the user explicitly overrides test="standard"
      ## on an ordinal model, since test defaults to "none" there
      ## (R/blavaan.R) -- see R/ppp.R for the on-demand second-step path.
      ppp <- pp_limited_info(lavobject)$ppval
    }
  } else {
    ppp <- postpred(samplls, lavobject)$ppval
  }
        
  TEST[[1]] <- list(test="mloglik",
                    stat=as.numeric(mll),
                    stat.group=as.numeric(NA),
                    df=as.integer(NA),
                    refdistr="NA",
                    pvalue=as.numeric(NA))

  TEST[[2]] <- list(test="ppp",
                    ## DIC: 2*ll(theta_hat) - 4*mean(ll(theta_samp))
                    stat=as.numeric(ppp),
                    stat.group=as.numeric(NA),
                    df=as.integer(NA),
                    refdistr="NA",
                    pvalue=as.numeric(NA))

  TEST
}
