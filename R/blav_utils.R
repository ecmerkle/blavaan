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
                   lavobject      = NULL, # just to use lavPredict()
                   measure        = "logl",
                   casewise       = FALSE,
                   conditional    = FALSE){

    if(length(postsamp) > 0){
        lavmodel <- fill_params(postsamp, lavmodel, lavpartable)
    }

    if(conditional){
      eta <- fill_eta(postsamp, lavmodel, lavpartable,
                      lavsamplestats, lavdata)

      ## implied meanvec + covmat
      #mnvec <- lavaan:::computeYHAT(lavmodel, lavmodel@GLIST,
      #                              lavsamplestats, ETA = eta)
      lavobject@Model <- lavmodel
      mnvec <- lavPredict(lavobject, type="ov", ETA = eta)
      if(inherits(mnvec, "matrix")) mnvec <- list(mnvec)

      #covmat <- lavaan:::computeTHETA(lavmodel, lavmodel@GLIST)
      covmat <- lavInspect(lavobject, 'theta')
      if(inherits(covmat, "matrix")) covmat <- list(covmat)
      ## to avoid warnings from mnormt::pd.solve
      covmat <- lapply(covmat, function(x){
        class(x) <- "matrix"
        x})

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
                  warning("blavaan WARNING: sample covariance matrix is not positive definite")
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

    loop.args <- list(X = 1:nsamps, FUN = function(i){
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

## fill in lavmodel with new parameter values (from a posterior sample)
fill_params <- function(postsamp      = NULL,
                        lavmodel      = NULL,
                        lavpartable   = NULL){
  if("jagpnum" %in% names(lavpartable)){
    filled <- lav_model_set_parameters(lavmodel,
                                       x = postsamp[lavpartable$jagpnum[!is.na(lavpartable$jagpnum)]])
  } else {
    filled <- lav_model_set_parameters(lavmodel,
                                       x = postsamp[lavpartable$stansumnum[lavpartable$free > 0][order(lavpartable$free[lavpartable$free > 0])]])
  }
  filled
}

## re-arrange columns of parameter samples to match that of blavaan object
rearr_params <- function(mcmc         = NULL,
                         lavpartable  = NULL){
    if(inherits(mcmc, "list") & length(mcmc) > 1){
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
        fullmat[,lavpartable$stansumnum[lavpartable$free > 0][order(lavpartable$free[lavpartable$free > 0])]]
    }
}

## iteration numbers for samp_lls and postpred
sampnums <- function(lavmcmc, thin, lout = FALSE){
    if("mcmc" %in% names(lavmcmc)){
        niter <- nrow(lavmcmc$mcmc[[1]])
    } else {
        niter <- dim(as.array(lavmcmc))[1]
    }
    #nsamps <- min(1000,floor(niter/thin))
    if(lout){
        psamp <- seq(1, niter, length.out = thin)
    } else {
        psamp <- seq(1, niter, thin)
    }

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
        dpars <- suppressWarnings(as.numeric(pricom[2:length(pricom)]))
        ## handle text like sqrt, ^, etc
        if(any(is.na(dpars))){
            nas <- which(is.na(dpars))
            for(i in nas){
                dpars[i] <- eval(parse(text=pricom[i+1]))
            }
        }
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

dist2r <- function(priors, target){
    ## convert jags/stan priors to R distribution
    ## return lists where distribution + parameters (+ truncation)
    ## appear as character vectors.

    ## explicitly change sqrt() to ^.5, because it may often be
    ## used to express sd's
    gsub("sqrt\\((.*)\\)\\).*", "\\1^.5\\)", priors)

    if(target == "jags"){
        out <- jagsdist2r(priors)
    } else if(target == "stan"){
        ## TODO need exported, or reverse rstan::lookup()
        #rosetta <- rstan:::rosetta
        ## alternate way to possibly get around export
        rloc <- paste0(system.file("R", package="rstan"), "/sysdata")
        lazyLoad(rloc)
        rosetta <- rosetta
        prisplit <- strsplit(priors, "[, ()]+")
        pridist <- sapply(prisplit, function(x) x[1])
        newdist <- rosetta$RFunction[match(pridist, rosetta$StanFunction)]
        for(i in 1:length(newdist)){
            prisplit[[i]][1] <- newdist[i]
        }

        out <- prisplit
    }

    out
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
    lavpartable$user <- c(lavpartable$user, rep(2L, sum(nvars)))
    lavpartable$group <- c(lavpartable$group, rep(1L, sum(nvars)))
    lavpartable$block <- c(lavpartable$block, rep(1L, sum(nvars)))
    lavpartable$free <- c(lavpartable$free, rep(0L, sum(nvars)))
    lavpartable$ustart <- c(lavpartable$ustart, rep(NA, sum(nvars)))
    lavpartable$exo <- c(lavpartable$exo, rep(0L, sum(nvars)))
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
             "                   ", paste(ov.names[forbid.idx], collapse = " "))
    }
}

## compute undirected K-L divergence across all draws
## (each draw paired with one from another chain)
samp_kls <- function(lavjags        = NULL,
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

    ## need to implement plummer's approach of generating y_rep
    ##mis <- FALSE
    ##if(any(is.na(unlist(lavdata@X)))) mis <- TRUE
    ##if(mis | lavoptions$categorical) stop("blavaan ERROR: K-L divergence not implemented for missing data or ordinal variables.")

    itnums <- sampnums(lavjags, thin = thin)
    lavmcmc <- lapply(lavmcmc, function(x) x[itnums,])
    draws <- do.call("rbind", lavmcmc)

    ndraws <- nrow(draws)
    halfdraws <- floor(ndraws/2)
    ngroups <- lavsamplestats@ngroups

    klres <- rep(NA, halfdraws)
    for(i in 1:halfdraws){
        lavmodel0 <- fill_params(draws[i,], lavmodel, lavpartable)
        lavmodel1 <- fill_params(draws[(halfdraws + i),], lavmodel,
                                 lavpartable)

        if(conditional){
            eta0 <- fill_eta(draws[i,], lavmodel, lavpartable,
                             lavsamplestats, lavdata)
            eta1 <- fill_eta(draws[(halfdraws + i),], lavmodel,
                             lavpartable, lavsamplestats, lavdata)

            lavobject@Model <- lavmodel0
            mnvec0 <- lavPredict(lavobject, type="ov", ETA=eta0)
            if(inherits(mnvec0, "matrix")) mnvec0 <- list(mnvec0)
            cmat0 <- lavInspect(lavobject, 'theta')
            if(inherits(cmat0, "matrix")) cmat0 <- list(cmat0)

            lavobject@Model <- lavmodel1
            mnvec1 <- lavPredict(lavobject, type="ov", ETA=eta1)
            if(inherits(mnvec1, "matrix")) mnvec1 <- list(mnvec1)
            cmat1 <- lavInspect(lavobject, 'theta')
            if(inherits(cmat1, "matrix")) cmat1 <- list(cmat1)

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
fill_eta <- function(postsamp, lavmodel, lavpartable, lavsamplestats, lavdata){
    nlv <- length(lavmodel@GLIST$alpha)
    etapars <- grepl("^eta", names(postsamp))
    cnums <- strsplit(names(postsamp)[etapars], "\\[|,|\\]")
    cnums <- sapply(cnums, function(x) as.numeric(x[3]))
    etavec <- postsamp[etapars][order(cnums)]

    ## need to worry about (1) excluding phantom lvs
    ## and (2) including dummy lvs
    foundlvs <- sum(etapars)/lavsamplestats@ntotal
    etamat <- matrix(etavec, lavsamplestats@ntotal, foundlvs)
    if(foundlvs < nlv) etamat <- cbind(etamat, matrix(0, lavsamplestats@ntotal, (nlv - foundlvs)))

    ngroups <- lavsamplestats@ngroups

    eta <- vector("list", ngroups)
    for(g in 1:ngroups){
      eta[[g]] <- etamat[lavdata@case.idx[[g]], 1:nlv, drop = FALSE]

      ## fill in eta with dummys, if needed
      dummyov <- c(lavmodel@ov.x.dummy.ov.idx[[g]], lavmodel@ov.y.dummy.ov.idx[[g]])
      dummylv <- c(lavmodel@ov.x.dummy.lv.idx[[g]], lavmodel@ov.y.dummy.lv.idx[[g]])
      if(length(dummyov) > 0){
        eta[[g]][, dummylv] <- lavdata@X[[g]][, dummyov]
      }
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

## now defunct:
## get various fit metrics from a fitted model for each
## posterior draw
samp_idx <- function(lavjags        = NULL,
                     lavmodel       = NULL,
                     lavpartable    = NULL,
                     lavsamplestats = NULL,
                     lavoptions     = NULL,
                     lavcache       = NULL,
                     lavdata        = NULL,
                     lavmcmc        = NULL,
                     thin           = 1,
                     measure        = "logl"){
    itnums <- sampnums(lavjags, thin = thin)
    nsamps <- length(itnums)
    lavmcmc <- make_mcmc(lavjags)

    nchain <- length(lavmcmc)
    idxmat <- matrix(NA, nsamps, nchain)

    for(i in 1:nsamps){
        for(j in 1:nchain){
            idxmat[i,j] <- get_ll(lavmcmc[[j]][itnums[i],],
                                  lavmodel,
                                  lavpartable,
                                  lavsamplestats,
                                  lavoptions,
                                  lavcache,
                                  lavdata,
                                  measure)[1]
        }
    }

    idxmat
}

make_mcmc <- function(mcmcout){
  ## extract mcmc draws from jags/stan object
  if(inherits(mcmcout, "runjags")){
    lavmcmc <- mcmcout$mcmc
  } else {
    ## for stan: as.array() gives parameters in a different order from summary()
    ##           so reorder
    tmpsumm <- rstan::summary(mcmcout)
    lavmcmc <- as.array(mcmcout)
    lavmcmc <- lapply(seq(dim(lavmcmc)[2]), function(x) lavmcmc[,x,])
    reord <- match(rownames(tmpsumm$summary), colnames(lavmcmc[[1]]))
    lavmcmc <- lapply(lavmcmc, function(x) x[,reord])
  }
  lavmcmc
}

## check that a package is installed via requireNamespace
pkgcheck <- function(x){
  suppressMessages(requireNamespace(x, quietly = TRUE))
}

pkgload <- function(x){
  try(suppressMessages(attachNamespace(x)), silent = TRUE)
}

## get plabels that have "wiggle"
wiglabels <- function(lavpartable, wiggle, wiggle.sd, target = "stan"){
  ## allowable group.equal names
  gqnames <- c("loadings", "intercepts", "regressions", "means", "thresholds")
  gqops <- c("=~", "~1", "~", "~1", "|")
  lv.names <- unique(unlist(lav_partable_attributes(lavpartable, pta=NULL)$vnames$lv))
  lpt <- lavpartable[lavpartable$label != "" & !is.na(lavpartable$label),]

  tmplabs <- lapply(wiggle, function(x){
    if(any(grepl(x, lpt$label))){
      if(any(lpt$op[lpt$label == x] == "~~")){
        stop("blavaan ERROR: wiggle cannot be used on variance parameters.")
      }
      lpt$plabel[lpt$label == x]
    } else if(x %in% gqnames){
      wname <- which(gqnames == x)
      tmppt <- lpt[lpt$op == gqops[wname],]
      if(x == 'intercepts'){
        tmppt <- tmppt[!(tmppt$lhs %in% lv.names),]
      }
      if(x == 'means'){
        tmppt <- tmppt[tmppt$lhs %in% lv.names,]
      }
      if(NROW(tmppt) == 0L) stop(paste0("blavaan ERROR: use of wiggle='", x, "' also requires group.equal='", x, "'."))

      out <- lapply(unique(tmppt$label), function(y){
        tmppt$plabel[tmppt$label == y]
        })
      out
    } else {
      stop("blavaan ERROR: poorly-specified wiggle argument (cannot be used on variances).")
    }
  })

  ## fix list nesting, in case group.equal was used
  outlist <- NULL
  if(length(tmplabs) == 1 & inherits(tmplabs[[1]], "list")){
    tmplabs <- tmplabs[[1]]
  }

  for(i in 1:length(tmplabs)){
    if(inherits(tmplabs[[i]], "list")){
      tmpelem <- tmplabs[[i]]
    } else {
      tmpelem <- list(tmplabs[[i]])
    }
    if(length(tmpelem[[1]]) > 1){
      outlist <- c(outlist, tmpelem)
    }
  }

  ## prior for partable
  if(!("prior" %in% names(lavpartable))) lavpartable$prior <- rep("", length(lavpartable$lhs))
  for(i in 1:length(outlist)){
    tmprows <- which(lavpartable$plabel %in% outlist[[i]])
    eqrows <- NULL
    if(target == "stan"){
      parname <- with(lavpartable, paste0(mat[tmprows[1]], "[", group[tmprows[1]], ",",
                                          row[tmprows[1]], ",", col[tmprows[1]], "]"))
      wigpri <- paste0("normal(", parname, ",", wiggle.sd, ")")
    } else {
      dname <- ifelse(grepl("stan", target), "normal(", "dnorm(")
      wigsc <- ifelse(grepl("stan", target), wiggle.sd, wiggle.sd^(-2))
      parname <- with(lavpartable, paste0(mat[tmprows[1]], "[", row[tmprows[1]], ",",
                                          col[tmprows[1]], ",", group[tmprows[1]], "]"))
      wigpri <- paste0(dname, parname, ",", wigsc, ")")

      ## nuke == rows
      eqrows <- with(lavpartable, which(op == "==" & (rhs %in% plabel[tmprows])))
    }
    lavpartable$prior[tmprows] <- c(lavpartable$prior[tmprows[1]],
                                    rep(wigpri, length(tmprows) - 1))
    if(length(eqrows) > 0){
      lavpartable <- lavpartable[-eqrows,]
    }
  }

  list(outlist = outlist, lavpartable = lavpartable)
}
  
