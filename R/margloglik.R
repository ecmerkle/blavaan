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

  ## convert prior strings to commands + parameters
  pri <- lavpartable$prior[urows]
  pricom <- strsplit(pri, "[, ()]+")

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
      ## convert to R parameterization of distributions
      par1 <- as.numeric(pricom[[i]][2])
      par2 <- as.numeric(pricom[[i]][3])
      ## sd parameterization of dnorm
      if(pricom[[i]][1] == "dnorm") par2 <- sqrt(1/par2)
      ## convert beta with (-1,1) support to beta with (0,1)
      if(grepl("rho", lavpartable$jlabel[urows][i])) thetstar[i] <- (thetstar[i]+1)/2
      ## convert to precision, vs variance (priors are on precisions)
      if(grepl("theta", lavpartable$jlabel[urows][i]) | grepl("psi", lavpartable$jlabel[urows][i])) thetstar[i] <- 1/thetstar[i]
      ## for truncated/censored distributions:
      support.prob <- 1
      ## is prior on sd or variance, is it truncated
      if(length(pricom[[i]]) > 3){
        prisd <- grepl("sd", pricom[[i]][4])
        privar <- grepl("var", pricom[[i]][4])
        if(prisd | privar){
          thetstar[i] <- 1/thetstar[i]
          if(prisd) thetstar[i] <- sqrt(thetstar[i])
        } else if(pricom[[i]][4] != "T"){
          warning("blavaan WARNING: Cannot yet handle censored priors in marginal log-likelihood computation.\nMarginal log-likelihood and Bayes factor approximations may be poor.\n")
        } else {
          cdf.fun <- gsub("^d", "p", pricom[[i]][1])
          if(length(pricom[[i]]) == 5){
            support.prob <- 1 - do.call(cdf.fun, list(as.numeric(pricom[[i]][5]), par1, par2))
          } else {
            support.prob <- do.call(cdf.fun, list(as.numeric(pricom[[i]][6]), par1, par2)) - do.call(cdf.fun, list(as.numeric(pricom[[i]][5]), par1, par2))
          }
        }
      }
      
      tmpdens <- do.call(pricom[[i]][1], list(thetstar[i], par1, par2, log=TRUE)) - log(support.prob)
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
