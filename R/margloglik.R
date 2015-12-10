margloglik <- function(lavpartable, lavmodel, lavoptions, 
                       lavsamplestats, lavdata, lavcache, lavjags) {
  ## compute marginal log-likelihood given model output
  ## the lavs are created via blavaan()

  bayesout <- lavjags
  
  ## unique parameters (remove equality constraints and
  ## cov parameters under srs priors)
  urows <- 1:length(lavpartable$jlabel)
  ## FIXME If there are inequality constraints or complex
  ##       equality constraints, we may have problems here:
  eqpars <- lavpartable$rhs[lavpartable$op == "=="]
  urows <- urows[!(lavpartable$plabel %in% eqpars) &
                 !grepl("@", lavpartable$plabel) &
                 !lavpartable$jlabel == ""]

  q <- length(urows)

  ## re-arrange crosscorr rows/columns to match partable,
  ## then remove redundant rows/columns
  rearr <- match(lavpartable$jlabel[urows], rownames(bayesout$crosscorr))
  rearr <- rearr[!is.na(rearr)] # assume NAs only come at end (for deviance, etc)

  Jinv <- diag(bayesout$summary$statistics[rearr,"SD"]) %*%
          bayesout$crosscorr[rearr,rearr] %*%
          diag(bayesout$summary$statistics[rearr,"SD"])

  cmatch <- match(lavpartable$jlabel[urows],
                  rownames(bayesout$summary$statistics),
                  nomatch=0)

  ## this is potentially under srs parameterization,
  ## need to change below to use lavaan's log-likelihood
  thetstar <- bayesout$summary$statistics[cmatch,"Mean"]
  names(thetstar) <- NULL

  ## convert prior strings to commands + parameters
  pri <- lavpartable$prior[urows]
  pricom <- jagsdist2r(pri)

  ## warn about fa priors
  if(any((sapply(pricom, length) == 0) & grepl("cov", lavpartable$jlabel[urows]))) warning("blavaan WARNING: marginal log-likelihoods under ov.cp=fa or lv.cp=fa may be unstable.")

  ## evaluate each prior
  priloglik <- 0

  ## first deal with any wishart priors
  wps <- which(sapply(pricom, function(x) x[1] == "dwish"))
  if(length(wps) > 0){
    ngroups <- max(lavpartable$group)
    dimen <- sum(grepl("dwish", lavpartable$prior) &
                 grepl("psi", lavpartable$jlabel) &
                 grepl(",1]", lavpartable$jlabel))
    for(k in 1:ngroups){
      tmpmat <- diag(lavpartable$est[grepl("psi", lavpartable$jlabel) &
                                     grepl(paste(",", ngroups, "]", sep=""), lavpartable$jlabel)])
      covlocs <- grepl("cov", lavpartable$jlabel) &
                                                   grepl("dwish", lavpartable$prior) &
                                                   grepl(paste(",", ngroups, "]", sep=""), lavpartable$jlabel)
      if(sum(covlocs) > 0){
        tmpmat[lower.tri(tmpmat)] <- lavpartable$est[covlocs]
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

  ## pick out variance parameters under fa priors
  if(lavoptions$ov.cp == "fa"){
    ocpcovs <- which(lavpartable$lhs[urows] %in% lavdata@ov.names &
                     lavpartable$rhs[urows] %in% lavdata@ov.names &
                     lavpartable$op[urows] == "~~" &
                     lavpartable$lhs[urows] != lavpartable$rhs[urows])
    ocpvars <- unique(c(lavpartable$lhs[urows][ocpcovs], lavpartable$rhs[urows][ocpcovs]))
  } else {
    ocpvars <- ""
  }
  if(lavoptions$lv.cp == "fa"){
    lcpcovs <- which(!(lavpartable$lhs[urows] %in% lavdata@ov.names) &
                     !(lavpartable$rhs[urows] %in% lavdata@ov.names) &
                     lavpartable$op[urows] == "~~" &
                     lavpartable$lhs[urows] != lavpartable$rhs[urows])
    lcpvars <- unique(c(lavpartable$lhs[urows][lcpcovs], lavpartable$rhs[urows][lcpcovs]))
  } else {
    lcpvars <- ""
  }

  for(i in 1:q){
    if(i %in% wps) next # already got it above
    ## TODO do we need to use dp in kernel density to get correct priors?
    if(length(pricom[[i]])==0 & grepl("cov", lavpartable$jlabel[urows][i])){
      ## deal with covariances under fa parameterization using
      ## kernel density estimation of covariance parameter's prior
      tmpdist <- rnorm(1e5,sd=sqrt(1/1e-4))*rnorm(1e5,sd=sqrt(1/1e-4))/rgamma(1e5,1,.5)
      tmpkd <- density(tmpdist)
      tmpdens <- log(approx(tmpkd$x, tmpkd$y, thetstar[i])$y)
    } else if(lavpartable$lhs[urows][i] %in% c(ocpvars, lcpvars) &
              lavpartable$op[urows][i] == "~~" &
              lavpartable$rhs[urows][i] == lavpartable$lhs[urows][i]){
      ## kernel density estimation of fa variance's prior
      tmpdist <- (rnorm(1e5,sd=sqrt(1/1e-4))*rnorm(1e5,sd=sqrt(1/1e-4)))/rgamma(1e5,1,.5) +
                 1/rgamma(1e5, as.numeric(pricom[[i]][2]), as.numeric(pricom[[i]][3]))
      tmpkd <- density(tmpdist)
      tmpdens <- log(approx(tmpkd$x, tmpkd$y, thetstar[i])$y)
    } else {
      ## we have an explicit prior;
      ## convert to R parameterization of distributions & evaluate
      tmpdens <- eval_prior(pricom[[i]], thetstar[i], lavpartable$jlabel[urows][i])
    }
    priloglik <- priloglik + tmpdens
  }

  ## have lavaan calculate the log-likelihood
  ## switch off se, force test = "standard"
  lavoptions$se <- "none"
  lavoptions$test <- "standard"
  lavoptions$estimator <- "ML"
  ## control() is part of lavmodel (for now)
  lavmodel@control <- list(optim.method="none")

  ## to compute likelihood, we need to remove the rho rows
  ## from lavpartable
  rhos <- grep("rho", lavpartable$jlabel)
  if(length(rhos) > 0) lavpartable <- lapply(lavpartable, function(x) x[-rhos])
  lavpartable$plabel <- sapply(lavpartable$plabel, function(x) strsplit(x, "@")[[1]][1])

  fit.new <- try(lavaan(slotParTable = lavpartable,
                        slotModel = lavmodel,
                        slotOptions = lavoptions,
                        slotSampleStats = lavsamplestats,
                        slotData = lavdata,
                        slotCache = lavcache), silent=TRUE)

  if(!inherits(fit.new, "try-error")){
    loglik <- fitMeasures(fit.new, "logl")
  } else {
    loglik <- NA
  }
  ##print(c(priloglik, loglik, log(det(Jinv))))
  margloglik <- (q/2)*log(2*pi) + log(det(Jinv))/2 +
                priloglik + loglik
  names(margloglik) <- ""
  margloglik
}
