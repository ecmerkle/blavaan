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
    lavmcmc <- lapply(lavmcmc, function(x) x[,reord])

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
  
