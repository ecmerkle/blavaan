## calculate model log-likelihood given some sampled parameter values
## (with lvs integrated out by default)
get_ll <- function(postsamp       = NULL, # one posterior sample
                   lavmodel       = NULL,
                   lavpartable    = NULL,
                   lavsamplestats = NULL,
                   lavoptions     = NULL,
                   lavcache       = NULL,
                   lavdata        = NULL,
                   lavobject      = NULL, # just to use lavPredict()
                   measure        = "logl",
                   casewise       = FALSE,
                   conditional    = FALSE){

    if(length(postsamp) > 0){
        lavmodel <- fill_params(postsamp, lavmodel, lavpartable)
    }

    if(conditional){
      implied <- cond_moments(postsamp, lavmodel, lavpartable, lavsamplestats,
                              lavdata, lavobject)
    } else {
      implied <- lav_model_implied(lavmodel)
    }
  
    ## check for missing, to see if we can easily get baseline ll for chisq
    mis <- FALSE
    if(any(is.na(unlist(lavdata@X)))){
      mis <- TRUE
      lavmvh1 <- getFromNamespace("lav_mvnorm_missing_h1_estimate_moments", "lavaan")
      lavmvll <- getFromNamespace("lav_mvnorm_missing_loglik_data", "lavaan")
    }

    if(measure[1] %in% c("logl", "chisq") & !mis){
        if(casewise){
            ll.samp <- rep(NA, sum(unlist(lavdata@nobs)))
        } else if(conditional){
            ll.samp <- 0
        } else {
            ## logl + baseline logl
            ll.samp <- c(0,0)
        }

        for(g in 1:length(implied$cov)){
            if(conditional){
              mnvec <- implied$mean[[g]]
            } else {
              mnvec <- as.numeric(implied$mean[[g]])
            }

            ## ensure symmetric:
            cmat <- (implied$cov[[g]] + t(implied$cov[[g]]))/2
            tmpll <- try(dmnorm(lavdata@X[[g]], mnvec, cmat, log=TRUE))
            if(inherits(tmpll, "try-error")) tmpll <- NA

            ## subtract logl.X
            x.idx <- lavsamplestats@x.idx[[g]]
            if(!is.null(x.idx) && length(x.idx) > 0L){
                Mu.X <- lavsamplestats@mean.x[[g]]
                Sigma.X <- lavsamplestats@cov.x[[g]]
                if(is.null(Mu.X)){
                    Mu.X <- mnvec[x.idx]
                }
                if(is.null(Sigma.X)){
                    Sigma.X <- cmat[x.idx, x.idx, drop=FALSE]
                }
                tmpll.x <- try(dmnorm(lavdata@X[[g]][,x.idx], Mu.X, Sigma.X, log=TRUE))
                if(inherits(tmpll.x, "try-error")) tmpll.x <- NA
                tmpll <- tmpll - tmpll.x
            }

            if(!conditional){
                sampmn <- apply(lavdata@X[[g]], 2, mean, na.rm=TRUE)
                sampcov <- ((lavdata@nobs[[g]]-1)/(lavdata@nobs[[g]]))*cov(lavdata@X[[g]])

                basell <- try(dmnorm(lavdata@X[[g]], sampmn, sampcov, log=TRUE))
                if(inherits(basell, "try-error")){
                  basell <- NA
                  warning("blavaan WARNING: sample covariance matrix is not positive definite.", call. = FALSE)
                }

                if(!is.null(x.idx) && length(x.idx) > 0L){
                    Mu.X <- lavsamplestats@mean.x[[g]]
                    Sigma.X <- lavsamplestats@cov.x[[g]]
                    if(is.null(Mu.X)){
                        Mu.X <- sampmn[x.idx]
                    }
                    if(is.null(Sigma.X)){
                        Sigma.X <- sampcov[x.idx, x.idx, drop=FALSE]
                    }
                    tmpll.x <- try(dmnorm(lavdata@X[[g]][,x.idx], Mu.X, Sigma.X, log=TRUE))
                    if(inherits(tmpll.x, "try-error")) tmpll.x <- NA
                    basell <- basell - tmpll.x
                }
            }

            if(casewise){
                ll.samp[lavdata@case.idx[[g]]] <- tmpll
            } else if(conditional){
                ll.samp <- ll.samp + sum(tmpll)
            } else {
                tmpll <- c(sum(tmpll), sum(basell))
                #tmpll <- -2*(sum(tmpll) - sum(basell))
                ll.samp <- ll.samp + tmpll
            }

            ## handling missing data, but need EM to
            ## obtain saturated LL
            ## mis <- which(is.na(tmpll))
            ## nmis <- length(mis)
            ## if(nmis > 0){
            ##     ## do this by missingness pattern;
            ##     ## see nonnest2::llcont()
            ##     for(i in 1:nmis){
            ##         obs <- which(!is.na(lavdata@X[[g]][mis[i],]))
            ##         tmpll[mis[i]] <- dmnorm(lavdata@X[[g]][mis[i],obs],
            ##                                 mnvec[obs], cmat[obs,obs], log=TRUE)
            ##     }
            ## }

        }
    } else if(measure[1] %in% c("logl", "chisq") & !casewise) {
        if(lavoptions$target == "stan" | lavoptions$categorical){ ## FIXME for categorical
            tmpll <- NA # we'll get it from stan
        } else {
            tmpobj <- lavobject
            tmpobj@implied <- lav_model_implied(lavmodel)
            tmpll <- sum(llcont(tmpobj))
        }
        tmpsat <- 0
        for(g in 1:length(implied$cov)){
          ## high tolerance speed boost
          satmod <- lavmvh1(Y = lavdata@X[[g]], Mp = lavdata@Mp[[g]],
                            #max.iter = 20,
                            tol = 1e-2)
          tmpsat <- tmpsat + lavmvll(Y = lavdata@X[[g]], Mu = satmod$Mu, Sigma = satmod$Sigma,
                                     x.idx = lavsamplestats@x.idx[[g]])
        }
        ll.samp <- c(tmpll, tmpsat)
    } else {
        ## other measures require us to run lavaan
        lavoptions$se <- "none"
        lavoptions$test <- "standard"
        lavoptions$estimator <- "ML"
        ## control() is part of lavmodel (for now)
        lavoptions$optim.method <- "none"
        lavoptions$check.gradient <- FALSE
        if("control" %in% slotNames(lavmodel)){
            lavmodel@control <- list(optim.method="none")
        }
        ## FIXME: not sure when 'free' becomes numeric,
        ## but S4 doesn't like it in the lavaan call
        ## (only sometimes; related to extra monitors?)
        lavpartable$free <- as.integer(lavpartable$free)
        fit.samp <- try(lavaan(slotParTable = lavpartable,
                               slotModel = lavmodel,
                               slotOptions = lavoptions,
                               slotSampleStats = lavsamplestats,
                               slotData = lavdata,
                               slotCache = lavcache), silent=TRUE)
        if(!inherits(fit.samp, "try-error")){
            fit.samp@Options$se <- "standard" # for nonnest2
            fit.samp@test[[1]]$test <- "standard" # for do.fit=FALSE

            if(casewise){
                ll.samp <- llcont(fit.samp)
            } else if(measure[1] == "logl"){
                ll.samp <- fitMeasures(fit.samp,
                                       c("logl", "unrestricted.logl"))
            } else {
                ll.samp <- fitMeasures(fit.samp, measure)
            }
        } else {
            ll.samp <- NA
        }
    }

    ##TDJ: preserve names of other fitMeasures(), when requested by postpred()
    ATTR <- attributes(ll.samp)
    ll.samp <- as.numeric(ll.samp)
    attributes(ll.samp) <- ATTR
    ll.samp
}

## get log-likelihoods for each sampled parameter
samp_lls <- function(lavjags        = NULL,
                     lavmodel       = NULL,
                     lavpartable    = NULL,
                     lavsamplestats = NULL,
                     lavoptions     = NULL,
                     lavcache       = NULL,
                     lavdata        = NULL,
                     lavmcmc        = NULL,
                     lavobject      = NULL,
                     thin           = 1,
                     conditional    = FALSE){
    itnums <- sampnums(lavjags, thin = thin)
    nsamps <- length(itnums)

    nchain <- length(lavmcmc)

    if(lavoptions$target != "stan" | conditional) {
      loop.args <- list(X = 1:nsamps, future.seed = TRUE, FUN = function(i){
        tmpmat <- matrix(NA, nchain, 2)
        for(j in 1:nchain){
          tmpmat[j,1:2] <- get_ll(lavmcmc[[j]][itnums[i],],
                                  lavmodel,
                                  lavpartable,
                                  lavsamplestats,
                                  lavoptions,
                                  lavcache,
                                  lavdata,
                                  lavobject,
                                  conditional = conditional)
        }
        tmpmat})

      llmat <- do.call("future_lapply", loop.args)
      llmat <- array(unlist(llmat), c(nchain, 2, nsamps)) ## logl + baseline logl
      llmat <- aperm(llmat, c(3,1,2))
    } else {
      ## the model log-likelihoods have already been computed in stan
      llmat <- array(NA, c(nsamps, nchain, 2))
      lls <- loo::extract_log_lik(lavjags)
      llsat <- loo::extract_log_lik(lavjags, parameter_name = "log_lik_sat")
      for(j in 1:nchain){
        idx <- (j-1)*nsamps + itnums
        llmat[itnums,j,1] <- rowSums(lls[idx,])
        llmat[itnums,j,2] <- rowSums(llsat[idx,]) + llmat[itnums,j,1]
      }
    }
    llmat
}

## casewise log-likelihoods
case_lls <- function(lavjags        = NULL,
                     lavmodel       = NULL,
                     lavpartable    = NULL,
                     lavsamplestats = NULL,
                     lavoptions     = NULL,
                     lavcache       = NULL,
                     lavdata        = NULL,
                     lavmcmc        = NULL,
                     lavobject      = NULL,
                     conditional    = FALSE,
                     thin           = 1){

    ## mcmc draws always in list
    itnums <- sampnums(lavjags, thin=thin)
    nsamps <- length(itnums)
    nchain <- length(lavmcmc)

    ntot <- sum(unlist(lavdata@nobs))
    llmat <- matrix(NA, nchain*nsamps, sum(unlist(lavdata@nobs)))

    tmpres <- vector("list", nchain)
    for(j in 1:nchain){
      loop.args <- list(X = 1:nsamps, FUN = function(i, j){
        get_ll(lavmcmc[[j]][itnums[i],],
               lavmodel,
               lavpartable,
               lavsamplestats,
               lavoptions,
               lavcache,
               lavdata,
               lavobject,
               casewise = TRUE,
               conditional = conditional)}, j = j)
      tmpres[[j]] <- t(do.call("future_sapply", loop.args))
    }

    llmat <- do.call("rbind", tmpres)
      
    llmat
}
