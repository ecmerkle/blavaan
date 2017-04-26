## utility functions for blavaan

## calculate model log-likelihood given some sampled parameter
## values (with lvs integrated out by default)
get_ll <- function(postsamp       = NULL, # one posterior sample
                   lavmodel       = NULL, 
                   lavpartable    = NULL, 
                   lavsamplestats = NULL, 
                   lavoptions     = NULL, 
                   lavcache       = NULL,
                   lavdata        = NULL,
                   measure        = "logl",
                   casewise       = FALSE,
                   conditional    = FALSE){

    if(length(postsamp) > 0){
        lavmodel <- fill_params(postsamp, lavmodel, lavpartable)
    }

    if(conditional){
      stop("blavaan ERROR: conditional log-likelihoods currently unavailable.")
      eta <- fill_eta(postsamp, lavpartable, lavsamplestats, lavdata)

      ## implied meanvec + covmat
      ## TODO replace with lav_predict_yhat and lavInspect?
      ## (lav_predict_yhat unavailable from lavPredict with custom ETA)
      mnvec <- NULL
      covmat <- NULL
      #mnvec <- lavaan:::computeYHAT(lavmodel, lavmodel@GLIST,
      #                              lavsamplestats, ETA = eta)
      #covmat <- lavaan:::computeTHETA(lavmodel, lavmodel@GLIST)

      ngroups <- lavsamplestats@ngroups
      implied <- list(cov = covmat, mean = mnvec,
                      slopes = vector("list", ngroups),
                      th = vector("list", ngroups),
                      group.w = vector("list", ngroups))
    } else {
      implied <- lav_model_implied(lavmodel)
    }

    ## check for missing, to see if we can easily get baseline ll for chisq
    mis <- FALSE
    if(any(is.na(unlist(lavdata@X)))) mis <- TRUE

    if(measure %in% c("logl", "chisq") & !mis){
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

            if(!conditional){
                sampmn <- apply(lavdata@X[[g]], 2, mean, na.rm=TRUE)
                sampcov <- ((lavdata@nobs[[g]]-1)/(lavdata@nobs[[g]]))*cov(lavdata@X[[g]])

                basell <- dmnorm(lavdata@X[[g]], sampmn, sampcov, log=TRUE)
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
    } else {
        ## other measures require us to run lavaan
        lavoptions$se <- "none"
        lavoptions$test <- "standard"
        lavoptions$estimator <- "ML"
        ## control() is part of lavmodel (for now)
        lavoptions$optim.method <- "none"
        if("control" %in% slotNames(lavmodel)){
            lavmodel@control <- list(optim.method="none")
        }

        fit.samp <- try(lavaan(slotParTable = lavpartable,
                               slotModel = lavmodel,
                               slotOptions = lavoptions,
                               slotSampleStats = lavsamplestats,
                               slotData = lavdata,
                               slotCache = lavcache), silent=TRUE)
        if(!inherits(fit.samp, "try-error")){
            fit.samp@Options$se <- "standard" # for nonnest2

            if(casewise){
                ll.samp <- llcont(fit.samp)
            } else if(measure == "logl"){
                ll.samp <- c(fitMeasures(fit.samp, "logl"),
                             fitMeasures(fit.samp, "unrestricted.logl"))
            } else {
                ll.samp <- fitMeasures(fit.samp, measure)
            }
        } else {
            ll.samp <- NA
        }
    }

    as.numeric(ll.samp)
}

## get log-likelihoods for each sampled parameter
samp_lls <- function(lavjags        = NULL,
                     lavmodel       = NULL, 
                     lavpartable    = NULL, 
                     lavsamplestats = NULL, 
                     lavoptions     = NULL, 
                     lavcache       = NULL,
                     lavdata        = NULL,
                     thin           = 5,
                     conditional    = FALSE){
    itnums <- sampnums(lavjags, thin = thin)
    nsamps <- length(itnums)

    lavmcmc <- make_mcmc(lavjags)
    
    nchain <- length(lavmcmc)
    llmat <- array(NA, c(nsamps, nchain, 2)) ## logl + baseline logl

    for(i in 1:nsamps){
        for(j in 1:nchain){
            llmat[i,j,1:2] <- get_ll(lavmcmc[[j]][itnums[i],],
                                     lavmodel,
                                     lavpartable, 
                                     lavsamplestats, 
                                     lavoptions, 
                                     lavcache,
                                     lavdata,
                                     conditional = conditional)
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
                     thin           = 5){

    itnums <- sampnums(lavjags, thin=5)
    nsamps <- length(itnums)

    ## mcmc draws always in list
    lavmcmc <- make_mcmc(lavjags)
  
    nchain <- length(lavmcmc)

    llmat <- matrix(NA, nchain*nsamps, sum(unlist(lavdata@nobs)))

    for(i in 1:nsamps){
        for(j in 1:nchain){
            llmat[(i-1)*nchain + j,] <- get_ll(lavmcmc[[j]][itnums[i],],
                                               lavmodel,
                                               lavpartable, 
                                               lavsamplestats, 
                                               lavoptions, 
                                               lavcache,
                                               lavdata,
                                               casewise = TRUE)
        }
    }

    llmat
}

## fill in lavmodel with new parameter values (from a posterior sample)
fill_params <- function(postsamp      = NULL,
                        lavmodel      = NULL,
                        lavpartable   = NULL){
  if("jagpnum" %in% names(lavpartable)){
    filled <- lav_model_set_parameters(lavmodel,
                                       x = postsamp[lavpartable$jagpnum[!is.na(lavpartable$jagpnum)]])
  } else {
    filled <- lav_model_set_parameters(lavmodel,
                                       x = postsamp[lavpartable$stanpnum[lavpartable$free > 0][order(lavpartable$free[lavpartable$free > 0])]])
  }
  filled
}

## re-arrange columns of parameter samples to match that of blavaan object
rearr_params <- function(mcmc         = NULL,
                         lavpartable  = NULL){
    if(class(mcmc) == "list" & length(mcmc) > 1){
        fullmat <- do.call("rbind", mcmc)
    }
    fullmat <- mcmc[[1]]
    if(length(mcmc) > 1){
        for(i in 2:length(mcmc)){
            fullmat <- rbind(fullmat, mcmc[[i]])
        }
    }
    if("jagpnum" %in% names(lavpartable)){
        fullmat[,lavpartable$jagpnum[lavpartable$free > 0]]
    } else {
        fullmat[,lavpartable$stanpnum[lavpartable$free > 0][order(lavpartable$free[lavpartable$free > 0])]]
    }
}   

## iteration numbers for samp_lls and postpred
sampnums <- function(lavmcmc, thin){
    if("mcmc" %in% names(lavmcmc)){
        niter <- nrow(lavmcmc$mcmc[[1]])
    } else {
        niter <- dim(as.array(lavmcmc))[1]
    }
    nsamps <- min(1000,floor(niter/thin))
    psamp <- seq(1, niter, length.out=nsamps)

    psamp
}

## define blocks for multivariate normal priors
set_blocks <- function(partable){
    loadblks <- unique(partable$lhs[partable$op == "=~"])
    regintblks <- unique(partable$lhs[partable$op %in% c("~", "~1")])
    blks <- c(loadblks, regintblks)
    blks <- blks[!duplicated(blks)]

    blknum <- 0
    partable$blk <- rep(NA, length(partable$lhs))
    for(i in 1:length(blks)){
        parrows <- which(partable$lhs == blks[i] &
                         grepl("dnorm", partable$prior) &
                         !grepl("T\\(", partable$prior) &
                         !grepl("I\\(", partable$prior))

        if(length(parrows) < 2) next

        blknum <- blknum + 1

        partable$blk[parrows] <- blknum
    }
    partable
}

## evaluate prior density using result of jagsdist2r
eval_prior <- function(pricom, thetstar, pxname){
    ## check for truncation and [sd]/[var] modifiers
    trun <- which(pricom == "T")
    sdvar <- which(pricom %in% c("[sd]","[var]"))

    if(length(trun) > 0 | length(sdvar) > 0){
        snip <- min(c(trun, sdvar))
        dpars <- as.numeric(pricom[2:(snip-1)])
        if(length(c(trun, sdvar)) > 1){
            ## assume [sd],[var] come after truncation
            trunend <- sdvar - 1
        } else{
            trunend <- length(pricom)
        }
    } else {
        dpars <- as.numeric(pricom[2:length(pricom)])
    }

    ## thetstar modifications:
    ## convert to precision or sd, vs variance (depending on prior)
    if(grepl("theta", pxname) | grepl("psi", pxname)){
        ## FIXME assumes correlation prior under srs is dbeta
        if(length(sdvar) == 0 & pricom[1] != "dbeta") thetstar <- 1/thetstar
        if(any(grepl("\\[sd", pricom))) thetstar <- sqrt(thetstar)
    }
    ## dt() in R assumes mean=0, precision=1
    if(pricom[1] == "dt"){
        tmn <- dpars[2]
        tprec <- dpars[3]
        dpars <- dpars[1]
        thetstar <- (thetstar - as.numeric(tmn))*sqrt(tprec)
    }

    ## for truncated distributions:
    support.prob <- 1
    ## is prior truncated
    if(length(trun) > 0){
        ## FIXME deal with censored priors
        ## warning("blavaan WARNING: Cannot yet handle censored priors in marginal log-likelihood computation.\nMarginal log-likelihood and Bayes factor approximations may be poor.\n")
        cdf.fun <- gsub("^d", "p", pricom[1])
        if(trunend - trun == 1){
            ## FIXME assumes truncation from below, cannot
            ## handle truncation from above (without below)
            support.prob <- 1 - do.call(cdf.fun, c(as.numeric(pricom[trunend]), as.list(dpars)))
        } else {
            support.prob <- do.call(cdf.fun, c(as.numeric(pricom[trunend]), as.list(dpars))) - do.call(cdf.fun, c(as.numeric(pricom[(trun+1)]), as.list(dpars)))
        }
    }

    dens <- do.call(pricom[1], c(thetstar, as.list(dpars), log=TRUE)) - log(support.prob)

    dens
}

## add extra monitors from jagextra to parameter table as defined parameters
add_monitors <- function(lavpartable, lavjags, jagextra){
    monres <- vector("list", length(jagextra$monitor))
    for(i in 1:length(jagextra$monitor)){
        tmploc <- grep(paste("^", jagextra$monitor[i], sep=""), rownames(lavjags$summaries))
        tmploc2 <- grep(paste("^", jagextra$monitor[i], sep=""), rownames(lavjags$psrf$psrf))
        monres[[i]] <- list(xlocs = tmploc, psrfloc = tmploc2,
                            xnms = rownames(lavjags$summaries)[tmploc],
                            nvars = length(tmploc))
    }
    xlocs <- sapply(monres, function(x) x$xlocs)
    psrflocs <- sapply(monres, function(x) x$psrfloc)
    xnms <- sapply(monres, function(x) x$xnms)
    nvars <- sapply(monres, function(x) x$nvars)
        
    nrs <- length(lavpartable$id)
    lavpartable$id <- c(lavpartable$id, (nrs + 1):(nrs + sum(nvars)))
    lavpartable$lhs <- c(lavpartable$lhs, xnms)
    lavpartable$op <- c(lavpartable$op, rep(":=", sum(nvars)))
    lavpartable$rhs <- c(lavpartable$rhs, rep("", sum(nvars)))
    lavpartable$user <- c(lavpartable$user, rep(2, sum(nvars)))
    lavpartable$group <- c(lavpartable$group, rep(1, sum(nvars)))
    lavpartable$block <- c(lavpartable$block, rep(1, sum(nvars)))
    lavpartable$free <- c(lavpartable$free, rep(0, sum(nvars)))
    lavpartable$ustart <- c(lavpartable$ustart, rep(NA, sum(nvars)))
    lavpartable$exo <- c(lavpartable$exo, rep(0, sum(nvars)))
    lavpartable$label <- c(lavpartable$label, rep("", sum(nvars)))
    lavpartable$plabel <- c(lavpartable$plabel, rep("", sum(nvars)))
    lavpartable$start <- c(lavpartable$start, rep(0, sum(nvars)))
    lavpartable$est <- c(lavpartable$est, lavjags$summaries[xlocs,'Mean'])
    lavpartable$se <- c(lavpartable$se, lavjags$summaries[xlocs,'SD'])
    lavpartable$prior <- c(lavpartable$prior, rep("", sum(nvars)))
    lavpartable$psrf <- c(lavpartable$psrf, lavjags$psrf$psrf[psrflocs,1])
    lavpartable$pxnames <- c(lavpartable$pxnames, xnms)
    lavpartable$jagpnum <- c(lavpartable$jagpnum, xlocs)
    lavpartable$logBF <- c(lavpartable$logBF, rep(NA, sum(nvars)))

    lavpartable
}

## forbidden variable names (don't confuse with parameter names)
namecheck <- function(ov.names){
    forbidden <- c("mu", "invthetstar", "invtheta", "nu", "lambda", "eta",
                   "mu_eta", "invpsistar", "invpsi", "alpha", "beta",
                   "rho", "theta", "psi", "rstar", "cov", "ibpsi", "bpsi", "iden", "yvec",
                   paste(".phant", 1:100, sep=""), "def")

    forbid.idx <- which(ov.names %in% forbidden)
    
    if(length(forbid.idx) > 0L){
        stop("blavaan ERROR: the following variable names must be changed:\n",
             "                   ", paste(forbidden[forbid.idx], collapse = " "))
    }
}

## compute undirected K-L divergence across all draws
## (each draw paired with one from another chain)
samp_kls <- function(draws          = NULL, # all chains in 1 matrix
                     lavmodel       = NULL, 
                     lavpartable    = NULL, 
                     lavsamplestats = NULL, 
                     lavoptions     = NULL, 
                     lavcache       = NULL,
                     lavdata        = NULL,
                     conditional    = FALSE){

    ## need to implement plummer's approach of generating y_rep
    ##mis <- FALSE
    ##if(any(is.na(unlist(lavdata@X)))) mis <- TRUE
    ##if(mis | lavoptions$categorical) stop("blavaan ERROR: K-L divergence not implemented for missing data or ordinal variables.")
  
    ndraws <- nrow(draws)
    halfdraws <- floor(ndraws/2)
    ngroups <- lavsamplestats@ngroups

    klres <- rep(NA, halfdraws)
    for(i in 1:halfdraws){
        lavmodel0 <- fill_params(draws[i,], lavmodel, lavpartable)
        lavmodel1 <- fill_params(draws[(halfdraws + i),], lavmodel,
                                 lavpartable)

        if(conditional){
            stop("blavaan ERROR: conditional kl-distance unavailable.")
            eta0 <- fill_eta(draws[i,], lavpartable, lavsamplestats,
                             lavdata)
            eta1 <- fill_eta(draws[(halfdraws + i),], lavpartable,
                             lavsamplestats, lavdata)

            ## mnvec0 <- lavaan:::computeYHAT(lavmodel0,
            ##                                lavmodel0@GLIST,
            ##                                lavsamplestats,
            ##                                ETA = eta0)
            ## cmat0 <- lavaan:::computeTHETA(lavmodel0,
            ##                                lavmodel0@GLIST)
            ## mnvec1 <- lavaan:::computeYHAT(lavmodel1,
            ##                                lavmodel1@GLIST,
            ##                                lavsamplestats,
            ##                                ETA = eta1)
            ## cmat1 <- lavaan:::computeTHETA(lavmodel1,
            ##                                lavmodel1@GLIST)
            implied0 <- list(cov = cmat0, mean = mnvec0,
                             slopes = vector("list", ngroups),
                             th = vector("list", ngroups),
                             group.w = vector("list", ngroups))
            implied1 <- list(cov = cmat1, mean = mnvec1,
                             slopes = vector("list", ngroups),
                             th = vector("list", ngroups),
                             group.w = vector("list", ngroups))
        } else {
            implied0 <- lav_model_implied(lavmodel0)
            implied1 <- lav_model_implied(lavmodel1)
        }

        tmpkl <- 0
        for(g in 1:lavsamplestats@ngroups){
            ## ensure symmetric:
            cmat0 <- (implied0$cov[[g]] + t(implied0$cov[[g]]))/2
            invcmat0 <- solve(cmat0)
            det0 <- det(cmat0)
            cmat1 <- (implied1$cov[[g]] + t(implied1$cov[[g]]))/2
            invcmat1 <- solve(cmat1)
            det1 <- det(cmat1)
            if(conditional){
                mnvec0 <- implied0$mean[[g]]
                mnvec1 <- implied1$mean[[g]]

                for(j in 1:nrow(mnvec0)){
                  tmpkl <- tmpkl + kl_und(mnvec0[j,], mnvec1[j,],
                                          cmat0, invcmat0, cmat1,
                                          invcmat1, det0, det1)
                }
            } else {
                mnvec0 <- as.numeric(implied0$mean[[g]])
                mnvec1 <- as.numeric(implied1$mean[[g]])

                tmpkl <- tmpkl + lavsamplestats@nobs[[g]] *
                         kl_und(mnvec0, mnvec1, cmat0, invcmat0,
                                cmat1, invcmat1, det0, det1)
            }
        }
        klres[i] <- tmpkl
    }
    klres
}        

## fill in eta matrices (1 per group, in list)
fill_eta <- function(postsamp, lavpartable, lavsamplestats, lavdata){
    ## FIXME: deal with phantom lvs
    nlv <- length(lav_partable_attributes(lavpartable)$vnames$lv[[1]])
    etapars <- grepl("^eta", names(postsamp))
    etamat <- matrix(postsamp[etapars], lavsamplestats@ntotal, nlv)
    ngroups <- lavsamplestats@ngroups

    eta <- vector("list", ngroups)
    for(g in 1:ngroups){
        eta[[g]] <- etamat[lavdata@case.idx[[g]], , drop = FALSE]
    }

    eta
}
        
## compute undirected K-L divergence between two normal distributions
kl_und <- function(mn0, mn1, cov0, invcov0, cov1, invcov1,
                   det0, det1){
  k <- nrow(cov0)

  kl01 <- sum(diag(invcov1 %*% cov0)) +
    t(mn1 - mn0) %*% invcov1 %*% (mn1 - mn0) -
    k + log(det1/det0)

  kl10 <- sum(diag(invcov0 %*% cov1)) +
    t(mn0 - mn1) %*% invcov0 %*% (mn0 - mn1) -
    k + log(det0/det1)

  (1/2) * (kl01 + kl10)
}

make_mcmc <- function(mcmcout){
  ## extract mcmc draws from jags/stan object
  if(class(mcmcout) == "runjags"){
    lavmcmc <- mcmcout$mcmc
  } else {
    lavmcmc <- as.array(mcmcout)
    lavmcmc <- lapply(seq(dim(lavmcmc)[2]), function(x) lavmcmc[,x,])
  }
  lavmcmc
}
