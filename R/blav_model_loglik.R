## calculate model log-likelihood given some sampled parameter values
## (with lvs integrated out by default)
## NB: column ordering of postsamp is obtained from make_mcmc()
get_ll <- function(postsamp       = NULL, # one posterior sample
                   lavobject      = NULL,
                   measure        = "logl",
                   casewise       = FALSE,
                   conditional    = FALSE,
                   standata       = NULL){

    if(lavobject@Options$categorical){
      ll.samp <- get_ll_ord(postsamp, lavobject, measure, casewise, conditional, standata)
    } else {
      ll.samp <- get_ll_cont(postsamp, lavobject, measure, casewise, conditional)
    }

    ll.samp
}

get_ll_cont <- function(postsamp       = NULL, # one posterior sample
                        lavobject      = NULL,
                        measure        = "logl",
                        casewise       = FALSE,
                        conditional    = FALSE){

  lavmodel <- lavobject@Model
  lavpartable <- lavobject@ParTable
  lavsamplestats <- lavobject@SampleStats
  lavoptions <- lavobject@Options
  lavcache <- lavobject@Cache
  lavdata <- lavobject@Data
  
  if(length(postsamp) > 0){
    lavmodel <- fill_params(postsamp, lavmodel, lavpartable)
  }

  if(conditional){
    implied <- cond_moments(postsamp, lavmodel, lavpartable, lavsamplestats,
                            lavdata, lavobject)
  } else {
    implied <- lav_model_implied(lavmodel, delta = (lavmodel@parameterization == "delta"))
  }
  
  ## check for missing, to see if we can easily get baseline ll for chisq
  mis <- FALSE
  if(any(is.na(unlist(lavdata@X)))){
    mis <- TRUE
    lavmvh1 <- getFromNamespace("lav_mvnorm_missing_h1_estimate_moments", "lavaan")
    lavmvll <- getFromNamespace("lav_mvnorm_missing_loglik_data", "lavaan")
  }

  if(measure[1] %in% c("logl", "chisq") & !mis & length(measure) == 1){
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
  } else if(measure[1] %in% c("logl", "chisq") & length(measure) == 1) {
    if(lavoptions$target == "stan"){
      tmpll <- NA # we'll get it from stan
    } else {
      tmpobj <- lavobject
      tmpobj@implied <- lav_model_implied(lavmodel, delta = (lavmodel@parameterization == "delta"))
      tmpobj@Options$estimator <- "ML"
      tmpll <- llcont(tmpobj)
    }
    if(casewise) {
      ll.samp <- tmpll
    } else {
      tmpsat <- 0
      for(g in 1:length(implied$cov)){
        ## high tolerance speed boost
        satmod <- lavmvh1(Y = lavdata@X[[g]], Mp = lavdata@Mp[[g]],
                          #max.iter = 20,
                          tol = 1e-2)
        tmpsat <- tmpsat + lavmvll(Y = lavdata@X[[g]], Mu = satmod$Mu, Sigma = satmod$Sigma,
                                   x.idx = lavsamplestats@x.idx[[g]])
      }
      ll.samp <- c(sum(tmpll), tmpsat)
    }
  } else {
    ## other measures require us to run lavaan
    lavoptions$se <- "none"
    lavoptions$test <- "standard"
    lavoptions$estimator <- "ML"
    if(lavoptions$categorical) lavoptions$estimator <- "DWLS"
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

get_ll_ord <- function(postsamp       = NULL, # one posterior sample
                       lavobject      = NULL,
                       measure        = "logl",
                       casewise       = FALSE,
                       conditional    = FALSE,
                       standata       = NULL){

  lavmodel <- lavobject@Model
  lavpartable <- lavobject@ParTable
  lavsamplestats <- lavobject@SampleStats
  lavoptions <- lavobject@Options
  lavcache <- lavobject@Cache
  lavdata <- lavobject@Data
  
  if(length(postsamp) > 0){
    lavmodel <- fill_params(postsamp, lavmodel, lavpartable)
  }
  
  if(conditional){
    lavmodel@GLIST$delta <- NULL # or else predictions will be rescaled
    implied <- cond_moments(postsamp, lavmodel, lavpartable, lavsamplestats,
                            lavdata, lavobject)
  } else {
    implied <- lav_model_implied(lavmodel, delta = (lavmodel@parameterization == "delta"))
  }

  if(is.null(standata)){
    if("mcmcdata" %in% names(lavobject@external)){
      standata <- lavobject@external$mcmcdata
    } else {
      stop("blavaan ERROR: Problem approximating ordinal log-likelihood.")
    }
  }

  ## continuous and ordinal indicators
  num.idx <- lavmodel@num.idx
  th.idx <- lavmodel@th.idx

  if("llnsamp" %in% names(lavoptions)){
    llnsamp <- lavoptions$llnsamp
  } else {
    ## we need to use tmvnsim for more than 20 dimensions
    if(length(unique(th.idx[[1]][th.idx[[1]] > 0])) > 20) {
      llnsamp <- lavoptions$llnsamp <- 100
    }
  }

  if(measure[1] %in% c("logl", "chisq") & length(measure) == 1){
    if(casewise){
      ll.samp <- rep(NA, sum(unlist(lavdata@nobs)))
    } else {
      ## logl + baseline logl
      ll.samp <- c(0,0)
    }

    Np <- standata$Np
    grpnum <- standata$grpnum
    Obsvar <- standata$Obsvar
    Nobs <- standata$Nobs
    startrow <- standata$startrow
    endrow <- standata$endrow

    ## full data
    YX <- matrix(NA, NROW(standata$YX), NCOL(standata$YX) + NCOL(standata$YXo))
    YX[, standata$contidx] <- standata$YX
    YX[, standata$ordidx] <- standata$YXo
        
    for(mm in 1:Np) {
      g <- grpnum[mm]
      r1 <- startrow[mm]
      r2 <- endrow[mm]
      obsidx <- Obsvar[mm, ]

      mnvec <- implied$mean[[g]]
      if(!conditional) {
        mnvec <- as.numeric(mnvec)
      }

      ## ensure symmetric:
      cmat <- (implied$cov[[g]] + t(implied$cov[[g]]))/2

      obsnum <- num.idx[[g]][num.idx[[g]] %in% obsidx[1:Nobs[mm]]]
      if(length(obsnum) > 0){
        tmpll <- try(dmnorm(YX[r1:r2, obsnum, drop=FALSE],
                            mnvec[obsnum],
                            cmat[obsnum, obsnum], log=TRUE))
        if(inherits(tmpll, "try-error")) tmpll <- NA

        ## subtract logl.X; assume always observed
        ## TODO: handle ordinal x.idx?
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
          tmpll.x <- try(dmnorm(YX[r1:r2, x.idx, drop=FALSE], Mu.X, Sigma.X, log=TRUE))
          if(inherits(tmpll.x, "try-error")) tmpll.x <- NA
          tmpll <- tmpll - tmpll.x
        }
      }
      ## baseline logliks only needed for ppp, which is now handled in Stan.
      basell <- 0L

      ## now condition on continuous, get ordinal ll by mnormt::sadmvn() or tmvnsim()
      TH.idx <- th.idx[[g]][th.idx[[g]] > 0]
      ord.idx <- unique(TH.idx)
      ord.idx <- ord.idx[ord.idx %in% obsidx[1:Nobs[mm]]]
      nord <- length(ord.idx)
      if(length(obsnum) > 0){
        s12s22i <- cmat[ord.idx, obsnum] %*% chol2inv(chol(cmat[obsnum, obsnum]))
        cov.ord <- cmat[ord.idx, ord.idx] - s12s22i %*% cmat[obsnum, ord.idx]
      } else {
        tmpll <- rep(0, r2 - r1 + 1)
        cov.ord <- cmat[ord.idx, ord.idx, drop=FALSE]
        mu.ord <- mnvec[ord.idx]
      }

      mm.in.group <- 1:lavmodel@nmat[g] + cumsum(c(0,lavmodel@nmat))[g]
      mms <- lavmodel@GLIST[mm.in.group]
      tau <- mms$tau

      ## thresholds for all cases
      lowtau <- hitau <- matrix(NA, NROW(YX[r1:r2,]), length(ord.idx))
      for(j in seq_len(nord)){
        tmptau <- c(-Inf, tau[TH.idx == ord.idx[j]], Inf)
        lowtau[,j] <- tmptau[YX[r1:r2,ord.idx[j]]]
        hitau[,j] <- tmptau[YX[r1:r2,ord.idx[j]] + 1]
      }

      for(i in r1:r2){
        llidx <- i - r1 + 1

        if(conditional){
          catprob <- pnorm(hitau[llidx,], mean = mnvec[i,ord.idx], sd = sqrt(diag(cmat)[ord.idx])) -
            pnorm(lowtau[llidx,], mean = mnvec[i,ord.idx], sd = sqrt(diag(cmat)[ord.idx]))
          lsigi <- sum( dbinom(1, size = 1, prob = catprob, log = TRUE) )
          tmpll[llidx] <- tmpll[llidx] + lsigi
        } else {
          if(length(obsnum) > 0){
            mu.ord <- mnvec[ord.idx] + s12s22i %*% (YX[i,obsnum] - mnvec[obsnum])
          }

          if("llnsamp" %in% names(lavoptions)){
            ## run tmvnsim to approximate marginal logl
            lsigi <- try(tmvnsim::tmvnsim(llnsamp, nord,
                                          lower = lowtau[llidx,], upper = hitau[llidx,],
                                          means = mu.ord, sigma = cov.ord), silent = TRUE)
            if(!inherits(lsigi, 'try-error')) lsigi <- mean(lsigi$wts)
          } else {
            lsigi <- try(mnormt::sadmvn(lowtau[llidx,], hitau[llidx,], mean = mu.ord, varcov = cov.ord, abseps = 1e-2))
          }

          if(inherits(lsigi, 'try-error')){
            tmpll[llidx] <- NA
          } else {
            tmpll[llidx] <- tmpll[llidx] + log(lsigi)
          }
        }
      }
            
      if(casewise){
        ll.samp[r1:r2] <- tmpll
      } else {
        tmpll <- c(sum(tmpll), sum(basell))
        #tmpll <- -2*(sum(tmpll) - sum(basell))
        ll.samp <- ll.samp + tmpll
      }
    }

    ## reorder to match original data
    if(casewise){
      rorig <- sapply(lavdata@Mp, function(x) unlist(x$case.idx))
      if(inherits(rorig, "list")){
        ## multiple groups
        for(ii in length(rorig)){
          rorig[[ii]] <- lavdata@case.idx[[ii]][rorig[[ii]]]
        }
        rorig <- unlist(rorig)
      }
      ll.samp[rorig] <- ll.samp
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
                     lavmcmc        = NULL,
                     lavobject      = NULL,
                     thin           = 1,
                     conditional    = FALSE,
                     standata       = NULL){
  lavoptions <- lavobject@Options
  itnums <- sampnums(lavjags, thin = thin)
  nsamps <- length(itnums)

  nchain <- length(lavmcmc)

  if(lavoptions$target != "stan" | conditional | lavoptions$categorical) {
    loop.args <- list(X = 1:nsamps, future.seed = TRUE, FUN = function(i){
      tmpmat <- matrix(NA, nchain, 2)
      for(j in 1:nchain){
        tmpmat[j,1:2] <- get_ll(lavmcmc[[j]][itnums[i],],
                                lavobject, conditional = conditional, standata = standata)
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
                     lavmcmc        = NULL,
                     lavobject      = NULL,
                     conditional    = FALSE,
                     thin           = 1){

  lavdata <- lavobject@Data
  
  ## mcmc draws always in list
  itnums <- sampnums(lavjags, thin=thin)
  nsamps <- length(itnums)
  nchain <- length(lavmcmc)

  ntot <- sum(unlist(lavdata@nobs))
  llmat <- matrix(NA, nchain*nsamps, sum(unlist(lavdata@nobs)))

  tmpres <- vector("list", nchain)
  for(j in 1:nchain){
    loop.args <- list(X = 1:nsamps, FUN = function(i, j){
      get_ll(lavmcmc[[j]][itnums[i],], lavobject,
             casewise = TRUE, conditional = conditional)},
      j = j, future.seed = TRUE)
    tmpres[[j]] <- t(do.call("future_sapply", loop.args))
  }

  llmat <- do.call("rbind", tmpres)
      
  llmat
}
