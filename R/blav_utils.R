## utility functions for blavaan

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
    sdvar <- which(pricom %in% c("[sd]","[var]","[prec]"))

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
        ## FIXME? assumes correlation prior under srs is beta or unif
        if(length(sdvar) == 0 & !(grepl("beta", pricom[1])) &
         !(grepl("unif", pricom[1]))) thetstar <- 1/thetstar
        if(any(grepl("\\[sd", pricom))) thetstar <- sqrt(thetstar)
        if(any(grepl("\\[prec", pricom))) thetstar <- 1/thetstar
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
        pridist <- paste0(pridist, "_lpdf")
        newdist <- rosetta$RFunction[match(pridist, rosetta$StanFunction)]
        for(i in 1:length(newdist)){
            if(!is.na(newdist[i])) prisplit[[i]][1] <- newdist[i]
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

make_mcmc <- function(mcmcout, stanlvs = NULL){
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
    lavmcmc <- lapply(lavmcmc, function(x) x[, reord, drop = FALSE])

    if(is.array(stanlvs)){
      stanlvs <- lapply(seq(dim(stanlvs)[2]), function(x) stanlvs[,x,])
      lavmcmc <- lapply(1:length(lavmcmc), function(x) cbind(lavmcmc[[x]], stanlvs[[x]]))
    }
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

  badnames <- c("residuals", "residual.covariances", "lv.variances", "lv.covariances")
  if(any(wiggle %in% badnames)){
    stop("blavaan ERROR: wiggle cannot be used on (co-)variance parameters.")
  }
  
  if(all(wiggle %in% gqnames)){
    ## ensure we have these things in this level of the model
    rmvars <- NULL
    for(i in 1:length(wiggle)){
      relop <- gqops[match(wiggle[i], gqnames)]
      if(!any(relop %in% lavpartable$op)) rmvars <- c(rmvars, i)
    }
    if(length(rmvars) > 0) wiggle <- wiggle[-rmvars]
  }

  if(!any(wiggle %in% lavpartable$label) && !any(wiggle %in% gqnames)) return( list(outlist = NULL, lavpartable = list(prior = NULL)) )
  
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
      stop("blavaan ERROR: poorly-specified wiggle argument.")
    }
  })

  ## fix list nesting, in case group.equal was used
  outlist <- NULL
  outsd <- NULL
  if(length(tmplabs) == 1 & inherits(tmplabs[[1]], "list")){
    tmplabs <- tmplabs[[1]]
  }

  if(length(wiggle.sd) == 1) wiggle.sd <- rep(wiggle.sd, length(tmplabs))
  
  for(i in 1:length(tmplabs)){
    if(inherits(tmplabs[[i]], "list")){
      tmpelem <- tmplabs[[i]]
    } else {
      tmpelem <- list(tmplabs[[i]])
    }
    multpars <- which(sapply(tmpelem, length) > 1)
    if(length(multpars) > 0){
      outlist <- c(outlist, tmpelem[multpars])
      outsd <- c(outsd, rep(wiggle.sd[i], length(tmpelem[multpars])))
    }
  }

  ## prior for partable
  if(!("prior" %in% names(lavpartable))) lavpartable$prior <- rep("", length(lavpartable$lhs))
  stanpris <- lavpartable$prior
  for(i in 1:length(outlist)){
    tmprows <- which(lavpartable$plabel %in% outlist[[i]])
    eqrows <- NULL
    if(target == "stan"){
      parname <- with(lavpartable, paste0(mat[tmprows[1]], "[", group[tmprows[1]], ",",
                                          row[tmprows[1]], ",", col[tmprows[1]], "]"))
      wigpri <- paste0("normal(", parname, ",", outsd[i], ")")
      spri <- paste0("normal(0,", outsd[i], ")")
    } else {
      dname <- ifelse(grepl("stan", target), "normal(", "dnorm(")
      wigsc <- ifelse(grepl("stan", target), wiggle.sd, wiggle.sd^(-2))
      parname <- with(lavpartable, paste0(mat[tmprows[1]], "[", row[tmprows[1]], ",",
                                          col[tmprows[1]], ",", group[tmprows[1]], "]"))
      wigpri <- spri <- paste0(dname, parname, ",", wigsc, ")")

      ## nuke == rows
      eqrows <- with(lavpartable, which(op == "==" & (rhs %in% plabel[tmprows])))
    }
    lavpartable$prior[tmprows] <- c(lavpartable$prior[tmprows[1]],
                                    rep(wigpri, length(tmprows) - 1))
    stanpris[tmprows] <- c(stanpris[tmprows[1]], rep(spri, length(tmprows) - 1))
    if(length(eqrows) > 0){
      lavpartable <- lavpartable[-eqrows,]
    }
  }

  list(outlist = outlist, lavpartable = lavpartable, stanpris = stanpris)
}
  
## wishart log-density
ldwish <- function(S, df, V){
  p <- NROW(S)

  logdetS <- log(det(S))

  chV <- chol(V)
  invV <- chol2inv(chV)
  logdetV <- 2*sum(log(diag(chV)))
  lmvgamma <- sum(log(gamma((df + 1 - 1:p)/2)))
  
  return(.5 * ((df - p - 1) * logdetS - sum(diag(invV %*% S)) - df*p*log(2) - df * logdetV - log(pi) * p * (p - 1)/2) - lmvgamma)
}

## for missing data, index row in data matrix based on index in Mp
## tricky because some cases can be fully missing, so they are excluded from the model
## Mp is @Data@Mp; case.idx is @Data@case.idx
## if exclude.empty, indexes are decreased to go from 1 to sum(lavInspect(,"nobs"))
## if !exclude.empty, indexes can go from 1 to sum(lavInspect(,"norig"))
Mp2dataidx <- function(Mp, case.idx, exclude.empty = TRUE){
  Ng <- length(Mp)

  for(g in 1:Ng){
    for(j in 1:length(Mp[[g]]$case.idx)){
      Mp[[g]]$case.idx[[j]] <- case.idx[[g]][Mp[[g]]$case.idx[[j]]]
    }
  }

  out <- unlist(sapply(Mp, function(x) unlist(x$case.idx)))
  if(exclude.empty){
    out <- rank(out)
  }

  out
}

## check for restricted model covariance matrices, which causes problems for priors
checkcovs <- function(lavobject){
  free <- lavInspect(lavobject, 'free')

  ## ensure each list entry is one group
  if (inherits(free[[1]], "matrix")) free <- list(free)

  if (nrow(free[[1]]$psi) > 0) {
    psis <- lapply(free, function(x) x$psi)
    psinums <- sapply(psis, function(x) x[lower.tri(x)])
    diagpsi <- all(unlist(psinums) == 0L, na.rm = TRUE)
    fullpsi <- all(unlist(psinums) > 0L, na.rm = TRUE) & (anyDuplicated(unlist(psinums), MARGIN = 0) == 0L)
    ## check for blocks of free covariances that have no impact on each other
    psiblk <- sapply(psis, function(x) {
      x[!lower.tri(x)] <- 0
      frnums <- which(x > 0, arr.ind = TRUE)
      if (nrow(frnums) > 0) {
        blk <- all(!duplicated(as.numeric(frnums)))
      } else {
        blk <- TRUE
      }
      blk} )
    blkp <- all(psiblk)
  } else {
    diagpsi <- FALSE
    fullpsi <- TRUE
    blkp <- TRUE
  }

  if (nrow(free[[1]]$theta) > 0) {
    thets <- lapply(free, function(x) x$theta)
    thetnums <- sapply(thets, function(x) x[lower.tri(x)])
    diagthet <- all(unlist(thetnums) == 0L, na.rm = TRUE)
    ## surprising if this happens:
    fullthet <- all(unlist(thetnums) > 0L, na.rm = TRUE) & (anyDuplicated(unlist(thetnums), MARGIN = 0) == 0L)
    ## check for blocks of free covariances that have no impact on each other
    ## FIXME
    thetblk <- sapply(thets, function(x) {
      x[!lower.tri(x)] <- 0
      frnums <- which(x > 0, arr.ind = TRUE)
      if (nrow(frnums) > 0) {
        blk <- !duplicated(as.numeric(frnums))
      } else {
        blk <- TRUE
      }
      blk} )
    blkt <- all(unlist(thetblk))
  } else {
    diagthet <- FALSE
    fullthet <- TRUE
    blkt <- TRUE
  }

  list(diagpsi = diagpsi, fullpsi = fullpsi, diagthet = diagthet, fullthet = fullthet, dobf = (blkp && blkt))
}

## check whether model cov matrix is block diagonal
## uses matrices from lavInspect(, "free")
## eqcon is attributes(lavInspect(, "free"))$header
blkdiag <- function(mat, eqcon = NULL) {
  isblk <- TRUE
  cnum <- 1L
  nblks <- 0
  matrows <- NROW(mat)
  blkse <- matrix(0, matrows, 3) # start/end of each block, is it a fully unrestricted block
  while (cnum <= matrows) {
    ## ending row of this potential block
    currend <- which(mat[, cnum] > 0)
    if (length(currend) > 0) {
      currend <- max(currend)
    } else {
      ## a fixed diagonal entry, which is its own block
      currend <- cnum
    }

    nblks <- nblks + 1
    if (currend > cnum) {
      ## check that columns cnum+1 to max() also equal max()
      othend <- sapply((cnum + 1):currend, function(i) {
        nzents <- which(mat[,i] > 0)
        if (length(nzents) > 0) {
          out <- max(nzents, i)
        } else {
          out <- i
        }
        out})

      if (all(othend == currend)) {
        ## is this entire submatrix free?
        submat <- mat[cnum:currend, cnum:currend]
        submat <- submat[lower.tri(submat)]
        if (all(submat > 0) & all(!(submat %in% as.numeric(eqcon$rhs)))) {
          blkse[nblks,] <- c(cnum, currend, TRUE)
        } else {
          isblk <- FALSE
          blkse[nblks,] <- c(cnum, currend, FALSE)
        }
      } else {
        isblk <- FALSE
        blkse[nblks,] <- c(cnum, currend, FALSE)
      }
    } else if (currend == cnum) {
      ## 1x1 block
      blkse[nblks,] <- c(cnum, cnum, TRUE)
    } else {
      isblk <- FALSE
    }
    cnum <- currend + 1
  }

  if (nrow(blkse) > 0) blkse <- blkse[1:nblks, , drop = FALSE]
  
  list(isblk = isblk, nblks = nblks, blkse = blkse)
}

## level labels for two-level models. this is taken from the similar lavaan function
## with some small modifications.
blav_partable_level_values <- function(partable) {
    if(is.null(partable$level)) {
        level.values <- 1L
    } else if(is.numeric(partable$level)) {
        tmp <- partable$level[  partable$level > 0L &
                               !partable$op %in% c("==", "<", ">", ":=") ]
        level.values <- unique(na.omit(tmp))
    } else { # character
        tmp <- partable$level[ nchar(partable$level) > 0L &
                              !partable$op %in% c("==", "<", ">", ":=") ]
        level.values <- unique(na.omit(tmp))
    }

    level.values
}

## Approximate posterior modes
modeapprox <- function(draws) {
  if(!inherits(draws, "matrix")) draws <- as.matrix(draws)
  
  ## can the modeest package be used?
  out <- rep(NA, ncol(draws))
  if (suppressMessages(requireNamespace("modeest", quietly = TRUE))) {
    tryMode <- try(apply(draws, 2, modeest::mlv, method = "kernel", na.rm = TRUE),
                   silent = TRUE)
    if (!inherits(tryMode, "try-error") && is.numeric(tryMode)) {
      out <- tryMode
    }
  }

  if (all(is.na(out))) {
    ## if not, use the quick-and-dirty way
    dd <- try(apply(draws, 2, density, na.rm = TRUE), silent = TRUE)

    if (!inherits(dd, "try-error")) out <- sapply(dd, function(z) z$x[which.max(z$y)])
  }

  out
}
