margloglik <- function(lavpartable, lavmodel, lavoptions, 
                       lavsamplestats, lavdata, lavcache, lavjags,
                       VCOV, xest, stansumm) {
  ## compute marginal log-likelihood given model output
  ## the lavs are created via blavaan()
  bayesout <- lavjags
  target <- lavoptions$target

  ## unique parameters (remove equality constraints and
  ## cov parameters under srs priors)
  eqpars <- lavpartable$rhs[lavpartable$op == "=="]
  fixpars <- which(lavpartable$free == 0 | lavpartable$prior == "")
  urows <- 1:length(lavpartable$pxnames)
  urows <- urows[!is.na(lavpartable$pxnames) &
                 !(lavpartable$plabel %in% eqpars) &
                 !(urows %in% fixpars)]

  q <- length(urows)

  ## re-arrange crosscorr rows/columns to match partable,
  ## then remove redundant rows/columns
  rearr <- match(lavpartable$pxnames[urows], rownames(VCOV))

  Jinv <- VCOV[rearr,rearr]

  if(target == "jags"){
    summstats <- bayesout$summary$statistics
  } else if(target == "stan"){
    summstats <- stansumm
  }
  cmatch <- match(lavpartable$pxnames[urows],
                  rownames(summstats),
                  nomatch=0)

  ## this is potentially under srs parameterization (unless stanmarg)
  if(target == "jags"){
    thetstar <- summstats[cmatch,"Mean"]
  } else if(target == "stanclassic"){
    thetstar <- summstats[cmatch,"mean"]
  } else {
    thetstar <- lavpartable$est[urows]
  }
  names(thetstar) <- NULL

  ## convert prior strings to commands + parameters
  pri <- lavpartable$prior[urows]
  pricom <- dist2r(pri, target = target)

  ## warn about fa priors
  if(lavoptions$cp == "fa") warning("blavaan WARNING: marginal log-likelihoods under cp='fa' may be unstable.")

  ## evaluate each prior
  priloglik <- 0

  ## first deal with any wishart priors
  wps <- which(sapply(pricom, function(x) x[1] == "dwish"))
  if(length(wps) > 0){
    ngroups <- max(lavpartable$group)
    for(k in 1:ngroups){
      varpars <- which(grepl("dwish", lavpartable$prior) &
                       lavpartable$group == k &
                       lavpartable$lhs == lavpartable$rhs)
      dimen <- length(varpars)
      tmpmat <- diag(lavpartable$est[varpars])
      ## TODO? ensure that covpars are ordered the same as varpars?
      covpars <- which(grepl("dwish", lavpartable$prior) &
                       lavpartable$group == k &
                       lavpartable$lhs != lavpartable$rhs)
      if(length(covpars) > 0){
        tmpmat[lower.tri(tmpmat)] <- lavpartable$est[covpars]
      }
      tmpmat <- tmpmat + t(tmpmat)
      diag(tmpmat) <- diag(tmpmat)/2

      ## TODO do we really need MCMCpack, or should we just
      ## compute the log density ourselves?
      ## NB wishart on precision matrix, so need to invert:
      priloglik <- priloglik + log(MCMCpack::dwish(solve(tmpmat), (dimen+1), diag(dimen)))
    }
  } else {
    ## just put a 0 here to avoid the error
    wps <- 0
  }

  ## pick out var/cov parameters under fa priors
  ocpcovs <- which(lavpartable$lhs[urows] %in% unlist(lavdata@ov.names) &
                   lavpartable$rhs[urows] %in% unlist(lavdata@ov.names) &
                   lavpartable$op[urows] == "~~" &
                   lavpartable$lhs[urows] != lavpartable$rhs[urows])

  lv.names <- lav_partable_attributes(partable = lavpartable, pta = NULL)$vnames$lv
  lcpcovs <- which(lavpartable$lhs[urows] %in% unlist(lv.names) &
                   lavpartable$rhs[urows] %in% unlist(lv.names) &
                   lavpartable$op[urows] == "~~" &
                   lavpartable$lhs[urows] != lavpartable$rhs[urows])
  
  if(lavoptions$cp == "fa"){
    ocpvars <- which(lavpartable$lhs[urows] %in% unlist(lavdata@ov.names) &
                     lavpartable$rhs[urows] %in% unlist(lavdata@ov.names) &
                     lavpartable$op[urows] == "~~" &
                     lavpartable$lhs[urows] == lavpartable$rhs[urows])

    lcpvars <- which(lavpartable$lhs[urows] %in% unlist(lv.names) &
                     lavpartable$rhs[urows] %in% unlist(lv.names) &
                     lavpartable$op[urows] == "~~" &
                     lavpartable$lhs[urows] == lavpartable$rhs[urows])
  } else {
    ocpvars <- ""
    lcpvars <- ""
  }

  for(i in 1:q){
    if(i %in% wps) next # already got it above

    if(i %in% c(ocpcovs, lcpcovs)){
      if(lavoptions$cp == "srs"){
        ## for srs, convert to correlation and use eval_prior()
        var1 <- which(lavpartable$lhs == lavpartable$lhs[urows][i] &
                      lavpartable$lhs == lavpartable$rhs &
                      lavpartable$op == "~~" &
                      lavpartable$group == lavpartable$group[urows][i])
        var2 <- which(lavpartable$lhs == lavpartable$rhs[urows][i] &
                      lavpartable$lhs == lavpartable$rhs &
                      lavpartable$op == "~~" &
                      lavpartable$group == lavpartable$group[urows][i])
        tstar <- thetstar[i]/sqrt(lavpartable$est[var1] * lavpartable$est[var2])
        # convert to support on 0,1
        tstar <- (tstar + 1)/2

        tmpdens <- eval_prior(pricom[[i]], tstar, lavpartable$pxnames[urows][i])
      } else {
        ## fa; TODO? use dp in kernel density to get correct priors?
        ## deal with covariances under fa parameterization using
        ## kernel density estimation of covariance parameter's prior
        tmpdist <- rnorm(1e5,sd=sqrt(1/1e-4))*rnorm(1e5,sd=sqrt(1/1e-4))/rgamma(1e5,1,.5)
        tmpkd <- density(tmpdist)
        tmpdens <- log(approx(tmpkd$x, tmpkd$y, thetstar[i])$y)
      }
    } else if(i %in% c(ocpvars, lcpvars)){
      ## kernel density estimation of fa variance's prior
      partype <- ifelse(i %in% ocpvars, "itheta", "ipsi")
      varpri <- jagsdist2r(dpriors()[[partype]])
      tmpdist <- (rnorm(1e5,sd=sqrt(1/1e-4))*rnorm(1e5,sd=sqrt(1/1e-4)))/rgamma(1e5,1,.5) +
                 1/rgamma(1e5, as.numeric(varpri[[1]][2]), as.numeric(varpri[[1]][3]))
      tmpkd <- density(tmpdist)
      tmpdens <- log(approx(tmpkd$x, tmpkd$y, thetstar[i])$y)
    } else {
      ## we have an explicit prior;
      ## convert to R parameterization of distributions & evaluate
      tmpdens <- eval_prior(pricom[[i]], thetstar[i], lavpartable$pxnames[urows][i])
    }
    priloglik <- priloglik + tmpdens
  }

  loglik <- attr(xest, "fx")
  ## have lavaan calculate the log-likelihood
  ## switch off se, force test = "standard"
  ## lavoptions$se <- "none"
  ## lavoptions$test <- "standard"
  ## lavoptions$estimator <- "ML"
  ## ## control() is part of lavmodel (for now)
  ## lavoptions$optim.method <- "none"
  ## if("control" %in% slotNames(lavmodel)){
  ##   lavmodel@control <- list(optim.method="none")
  ## }

  ## fit.new <- try(lavaan(slotParTable = lavpartable,
  ##                       slotModel = lavmodel,
  ##                       slotOptions = lavoptions,
  ##                       slotSampleStats = lavsamplestats,
  ##                       slotData = lavdata,
  ##                       slotCache = lavcache), silent=TRUE)

  ## if(!inherits(fit.new, "try-error")){
  ##   loglik <- fitMeasures(fit.new, "logl")
  ## } else {
  ##   loglik <- NA
  ## }
  #print(c(priloglik, loglik, log(det(Jinv))))
  
  margloglik <- (q/2)*log(2*pi) + determinant(Jinv, logarithm=TRUE)$modulus/2 +
                priloglik + loglik
  names(margloglik) <- ""
  margloglik
}
