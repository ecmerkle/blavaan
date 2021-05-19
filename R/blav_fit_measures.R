setMethod("fitMeasures", signature(object = "blavaan"),
function(object, fit.measures = "all", baseline.model = NULL) {
    blav_fit_measures(object = object, fit.measures = fit.measures,
                     baseline.model = baseline.model)
})

# lowercase 'm'
setMethod("fitmeasures", signature(object = "blavaan"),
function(object, fit.measures = "all", baseline.model = NULL) {
    blav_fit_measures(object = object, fit.measures = fit.measures,
                     baseline.model = baseline.model)
})


blav_fit_measures <- function(object, fit.measures = "all", 
                              baseline.model = NULL) {

    # has the model converged?
    if(object@Fit@npar > 0L && !object@optim$converged &&
       class(object@external$mcmcout) != "NULL") {
        warning("blavaan WARNING: the chains may not have converged.", call. = FALSE)
    }

    # do we have a test statistic?
    bopts <- blavInspect(object, "options")
    if(bopts$test == "none") {
        if(bopts$target != "stan") {
            stop("blavaan ERROR: fit measures cannot be obtained when test=\"none\"")
        } else {
            warning("blavaan WARNING: not all fit measures are available when test=\"none\"",
                    call. = FALSE)
        }
    }

    if(bopts$prisamp) {
        warning("blavaan WARNING: These metrics are based on prior samples so may be meaningless.",
                call. = FALSE)
    }
  
    if("all" %in% fit.measures) {
        class.flag <- TRUE
    } else {
        class.flag <- FALSE
    }

    # collect info from the lavaan slots
    GLIST <- object@Model@GLIST

    # N versus N-1
    # this affects BIC, RMSEA, cn_01/05, MFI and ECVI
    # Changed 0.5-15: suggestion by Mark Seeto
    if(object@Options$estimator %in% c("ML","PML","FML") && 
       object@Options$likelihood == "normal") {
        N <- object@SampleStats@ntotal
    } else {
        N <- object@SampleStats@ntotal - object@SampleStats@ngroups
    }

    # Change 0.5-13: take into account explicit equality constraints!!
    # reported by Mark L. Taper (affects AIC and BIC)
    npar <- object@Fit@npar
    if(nrow(object@Model@con.jac) > 0L) {
        ceq.idx <- attr(object@Model@con.jac, "ceq.idx")
        if(length(ceq.idx) > 0L) {
            neq <- qr(object@Model@con.jac[ceq.idx,,drop=FALSE])$rank
            npar <- npar - neq
        }
    }

    fx <- object@Fit@fx
    fx.group <- object@Fit@fx.group
    meanstructure <- object@Model@meanstructure
    categorical   <- object@Model@categorical
    multigroup    <- object@Data@ngroups > 1L
    estimator     <- "ML" #object@Options$estimator
    test          <- object@Options$test
    G <- object@Data@ngroups  # number of groups
    ## TODO get these from something like get_ll
    ##X2 <- object@Fit@test[[1]]$stat
    ##df <- object@Fit@test[[1]]$df

    # define 'sets' of fit measures:
    fit.always <- c("npar")

    # logl at posterior means
    fit.logl <- "logl"
    
    # posterior predictive p
    fit.chisq <- "ppp" #, "chisq.pi.lower", "chisq.pi.upper")

    # information criteria
    fit.ic <- c("bic", "dic", "p_dic")
    if("sampkls" %in% names(object@external)){
        fit.ic <- c(fit.ic, "dic_jags", "p_dic_jags")
    }
    if("csamplls" %in% names(object@external)){
        fit.ic <- c(fit.ic, "dic_cond", "p_dic_cond")
    }
    if("csampkls" %in% names(object@external)){
        fit.ic <- c(fit.ic, "dic_cond_j", "p_dic_cond_j")
    }
    fit.ic <- c(fit.ic, "waic", "p_waic", "se_waic",
                "looic", "p_loo", "se_loo")
    if("csamplls" %in% names(object@external)){
        fit.ic <- c(fit.ic, "waic_cond", "p_waic_cond", "se_waic_cond",
                    "looic_cond", "p_loo_cond", "se_loo_cond")
    }
  
    fit.rmsea <- "rmsea"
    
    # marginal logl
    fit.mll <- "margloglik"

    # lower case
    fit.measures <- tolower(fit.measures)

    # select 'default' fit measures
    if(length(fit.measures) == 1L) {
        if(fit.measures == "default" | fit.measures == "all") {
            if(estimator == "ML") {
                fit.measures <- c(fit.always, fit.logl, fit.chisq, fit.ic, 
                                  fit.rmsea, fit.mll)
            } else {
                stop("blavaan ERROR: Invalid estimator.")

                #if(estimator == "MML") {
                #  fit.measures <- c(fit.always, fit.logl)
            }

        }
    }

    # main container
    indices <- list()

    if("npar" %in% fit.measures) {
        indices["npar"] <- npar
    }

    if("logl" %in% fit.measures) {
        indices["logl"] <- object@Fit@fx
    }
    
    # posterior predictive p
    if("ppp" %in% fit.measures & bopts$test != "none") {
        indices["ppp"] <- object@Fit@test[[2]]$stat
    }
    if(any(c("bic", "dic", "p_dic") %in% fit.measures & bopts$test != "none")) {
        df <- 2*(object@Fit@fx - mean(as.numeric(object@external$samplls[,,1])))
        indices["bic"] <- -2*object@Fit@fx + npar*log(N)
        indices["dic"] <- -2*object@Fit@fx + 2*df
        indices["p_dic"] <- df

        if("sampkls" %in% names(object@external)){
          dfj <- mean(object@external$sampkls)/2
          indices["dic_jags"] <- -2*object@Fit@fx + 2*dfj
          indices["p_dic_jags"] <- dfj
        }

        if("csamplls" %in% names(object@external)){
          if(!is.na(as.numeric(object@external$csamplls)[1])){
            cllmn <- mean(as.numeric(object@external$csamplls[,,1]))
          } else {
            cllmn <- NA
          }
          dfc <- 2*(object@external$cfx - cllmn)
          indices["dic_cond"] <- -2*object@external$cfx + 2*dfc
          indices["p_dic_cond"] <- dfc
        }

        if("csampkls" %in% names(object@external)){
          dfjc <- mean(object@external$csampkls)/2
          indices["dic_cond_j"] <- -2*object@external$cfx + 2*dfjc
          indices["p_dic_cond_j"] <- dfjc
        }
    }
    if(any(c("waic", "p_waic", "looic", "p_loo") %in% fit.measures)) {
        lavopt <- object@Options
        lavopt$estimator <- "ML"
        if(lavopt$target == "stan"){
          casells <- loo::extract_log_lik(object@external$mcmcout)
        } else {
          casells <- case_lls(object@external$mcmcout, object@Model,
                              object@ParTable, object@SampleStats,
                              lavopt, object@Cache,
                              object@Data, make_mcmc(object@external$mcmcout))
        }

        fitres <- waic(casells)
        fitse <- fitres$estimates[,'SE']
        fitres <- fitres$estimates[,'Estimate']
        indices["waic"] <- fitres[["waic"]]
        indices["p_waic"] <- fitres[["p_waic"]]
        indices["se_waic"] <- fitse[["waic"]]
        nchain <- blavInspect(object, "n.chains")
        ref <- relative_eff(casells, chain_id =
                            rep(1:nchain, each = nrow(casells)/nchain))
        fitres <- loo(casells, r_eff = ref)
        fitse <- fitres$estimates[,'SE']
        fitres <- fitres$estimates[,'Estimate']
        indices["looic"] <- fitres[["looic"]]
        indices["p_loo"] <- fitres[["p_loo"]]
        indices["se_loo"] <- fitse[["looic"]]

        if("csamplls" %in% names(object@external)){
            if("stanlvs" %in% names(object@external)){
                samps <- make_mcmc(object@external$mcmcout, object@external$stanlvs)
            } else {
                samps <- make_mcmc(object@external$mcmcout)
            }
            casells <- case_lls(object@external$mcmcout, object@Model,
                                object@ParTable, object@SampleStats,
                                lavopt, object@Cache,
                                object@Data, samps,
                                lavobject = object, conditional = TRUE)

            fitres <- waic(casells)
            fitse <- fitres$estimates[,'SE']
            fitres <- fitres$estimates[,'Estimate']
            indices["waic_cond"] <- fitres[["waic"]]
            indices["p_waic_cond"] <- fitres[["p_waic"]]
            indices["se_waic_cond"] <- fitse[["waic"]]
            ref <- relative_eff(casells, chain_id =
                                rep(1:nchain, each = nrow(casells)/nchain))
            fitres <- loo(casells, r_eff = ref)
            fitse <- fitres$estimates[,'SE']
            fitres <- fitres$estimates[,'Estimate']
            indices["looic_cond"] <- fitres[["looic"]]
            indices["p_loo_cond"] <- fitres[["p_loo"]]
            indices["se_loo_cond"] <- fitse[["looic"]]
        }
    }
    if("margloglik" %in% fit.measures & test != "none") {
        indices["margloglik"] <- object@test[[1]]$stat
    }
    
    out <- unlist(indices[fit.measures])

    if(length(out) > 0L) {
        class(out) <- c("lavaan.vector", "numeric")
    } else {
        return( invisible(numeric(0)) )
    }

    out
}
