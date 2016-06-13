## utility functions for blavaan
mergejag <- function(lavpartable, coefvec){
    ## add results of jags translation to the partable
    ## TODO plus starting values
    plmatch <- sapply(coefvec$plabel, function(x) strsplit(x, "@")[[1]][1])
    pldups <- duplicated(plmatch)

    ## distinguish rhos from others
    rhos <- grepl("rho", coefvec$jlabel)

    ## match existing parameters with partable
    parlocs <- match(plmatch[!rhos & !pldups], lavpartable$plabel)
    lavpartable$plabel[parlocs] <- coefvec$plabel[!rhos & !pldups]

    ## add priors to partable
    nrlpt <- length(lavpartable$plabel)
    if(!("prior" %in% names(lavpartable))) lavpartable$prior <- rep("", nrlpt)
    lavpartable$prior[parlocs] <- coefvec$prior[!rhos & !pldups]

    ## add jags labels (jlabel)
    lavpartable$jlabel <- rep("", nrlpt)
    lavpartable$jlabel[parlocs] <- coefvec$jlabel[!rhos & !pldups]

    ## create new rho rows (that's a pun!) if needed and add info
    nrhos <- sum(rhos)
    if(nrhos > 0){
        lavpartable <- lapply(lavpartable, function(x) c(x, rep(NA, nrhos)))
        rhof <- nrlpt + 1
        rhol <- nrlpt + nrhos

        ## fill in most stuff the same as covariances
        covpars <- which(grepl("cov", lavpartable$jlabel) & !(grepl("dwish", lavpartable$prior)) & grepl("@", lavpartable$plabel))

        samevals <- c("id", "lhs", "op", "rhs", "user", "group", "ustart",
                      "exo", "label", "prior")

        for(i in 1:length(samevals)){
            nameloc <- which(names(lavpartable) == samevals[i])
            lavpartable[[nameloc]][rhof:rhol] <- lavpartable[[nameloc]][covpars]
        }

        ## remove priors from covariance parameters, because they technically have none?
        lavpartable$prior[covpars] <- ""
        lavpartable$free[rhof:rhol] <- 0L # so that parameterEstimates() works

        lavpartable$plabel[rhof:rhol] <- coefvec$plabel[rhos]
        lavpartable$jlabel[rhof:rhol] <- coefvec$jlabel[rhos]
    }
    
    lavpartable
}

## calculate model log-likelihood given some sampled parameter
## values (with lvs integrated out)
get_ll <- function(postsamp       = NULL, # one posterior sample
                   lavmodel       = NULL, 
                   lavpartable    = NULL, 
                   lavsamplestats = NULL, 
                   lavoptions     = NULL, 
                   lavcache       = NULL,
                   lavdata        = NULL,
                   measure        = "logl",
                   casewise       = FALSE){

    if(length(postsamp) > 0){
        lavmodel <- fill_params(postsamp, lavmodel, lavpartable)
    }

    implied <- lav_model_implied(lavmodel)

    ## check for missing, to see if we can easily get baseline ll for chisq
    mis <- FALSE
    if(any(is.na(unlist(lavdata@X)))) mis <- TRUE

    if(measure %in% c("logl", "chisq") & !mis){
        if(casewise){
            ll.samp <- rep(NA, sum(unlist(lavdata@nobs)))
        } else {
            ## logl + baseline logl
            ll.samp <- c(0,0)
        }
        
        for(g in 1:length(implied$cov)){
            mnvec <- as.numeric(implied$mean[[g]])
            ## ensure symmetric:
            cmat <- (implied$cov[[g]] + t(implied$cov[[g]]))/2
            tmpll <- dmnorm(lavdata@X[[g]], mnvec, cmat, log=TRUE)

            sampmn <- apply(lavdata@X[[g]], 2, mean, na.rm=TRUE)
            sampcov <- ((lavdata@nobs[[g]]-1)/(lavdata@nobs[[g]]))*cov(lavdata@X[[g]])

            basell <- dmnorm(lavdata@X[[g]], sampmn, sampcov, log=TRUE)

            if(casewise){
                ll.samp[lavdata@case.idx[[g]]] <- tmpll
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
        lavmodel@control <- list(optim.method="none")

        fit.samp <- try(lavaan(slotParTable = lavpartable,
                               slotModel = lavmodel,
                               slotOptions = lavoptions,
                               slotSampleStats = lavsamplestats,
                               slotData = lavdata,
                               slotCache = lavcache), silent=TRUE)

        if(casewise){
            ll.samp <- llcont(fit.samp)
        } else if(measure == "logl"){
            ll.samp <- c(fitMeasures(fit.samp, "logl"),
                         fitMeasures(fit.samp, "unrestricted.logl"))
        } else {
            ll.samp <- fitMeasures(fit.samp, measure)
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
                     thin           = 5){
    itnums <- sampnums(lavjags, thin=5)
    nsamps <- length(itnums)
    nchain <- length(lavjags$mcmc)
    llmat <- array(NA, c(nsamps, nchain, 2)) ## logl + baseline logl

    for(i in 1:nsamps){
        for(j in 1:nchain){
            llmat[i,j,1:2] <- get_ll(lavjags$mcmc[[j]][itnums[i],],
                                     lavmodel,
                                     lavpartable, 
                                     lavsamplestats, 
                                     lavoptions, 
                                     lavcache,
                                     lavdata)
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
    nchain <- length(lavjags$mcmc)

    llmat <- matrix(NA, nchain*nsamps, sum(unlist(lavdata@nobs)))

    for(i in 1:nsamps){
        for(j in 1:nchain){
            llmat[(i-1)*nchain + j,] <- get_ll(lavjags$mcmc[[j]][itnums[i],],
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
    ## this code related to coeffun():
    cnames <- lavpartable$jlabel
    rhos <- grep("rho", cnames)
    nn <- c(rhos, which(cnames == ""))
    if(length(nn) > 0) {
        cnames <- cnames[-nn]
    }
    cmatch <- match(cnames, names(postsamp), nomatch=0)

    lav_model_set_parameters(lavmodel, x = postsamp[cmatch])
}

## re-arrange columns of parameter samples to match that of blavaan object
rearr_params <- function(mcmc         = NULL,
                         lavpartable  = NULL){
    fullmat <- mcmc[[1]]
    if(length(mcmc) > 1){
        for(i in 2:length(mcmc)){
            fullmat <- rbind(fullmat, mcmc[[i]])
        }
    }
    
    cnames <- lavpartable$jlabel
    rhos <- grep("rho", cnames)
    nn <- c(rhos, which(cnames == ""))
    if(length(nn) > 0) {
        cnames <- cnames[-nn]
    }
    cmatch <- match(cnames, colnames(fullmat), nomatch=0)

    fullmat[,cmatch]
}   

## iteration numbers for samp_lls and postpred
sampnums <- function(lavjags, thin){
    niter <- nrow(lavjags$mcmc[[1]])
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
eval_prior <- function(pricom, thetstar, jlabel){
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
    ## convert beta with (-1,1) support to beta with (0,1)
    if(grepl("rho", jlabel)) thetstar <- (thetstar+1)/2
    ## convert to precision or sd, vs variance (depending on prior)
    if(jlabel %in% c("theta", "psi")){
        if(length(sdvar) == 0) thetstar <- 1/thetstar
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
        monres[[i]] <- list(xlocs = tmploc, xnms = rownames(lavjags$summaries)[tmploc],
                            nvars = length(tmploc))
    }
    xlocs <- sapply(monres, function(x) x$xlocs)
    xnms <- sapply(monres, function(x) x$xnms)
    nvars <- sapply(monres, function(x) x$nvars)
        
    nrs <- length(lavpartable$id)
    lavpartable$id <- c(lavpartable$id, (nrs + 1):(nrs + sum(nvars)))
    lavpartable$lhs <- c(lavpartable$lhs, xnms)
    lavpartable$op <- c(lavpartable$op, rep(":=", sum(nvars)))
    lavpartable$rhs <- c(lavpartable$rhs, rep("", sum(nvars)))
    lavpartable$user <- c(lavpartable$user, rep(2, sum(nvars)))
    lavpartable$group <- c(lavpartable$group, rep(1, sum(nvars)))
    lavpartable$free <- c(lavpartable$free, rep(0, sum(nvars)))
    lavpartable$ustart <- c(lavpartable$ustart, rep(NA, sum(nvars)))
    lavpartable$exo <- c(lavpartable$exo, rep(0, sum(nvars)))
    lavpartable$label <- c(lavpartable$label, rep("", sum(nvars)))
    lavpartable$plabel <- c(lavpartable$plabel, rep("", sum(nvars)))
    lavpartable$start <- c(lavpartable$start, rep(0, sum(nvars)))
    lavpartable$est <- c(lavpartable$est, lavjags$summaries[xlocs,'Mean'])
    lavpartable$se <- c(lavpartable$se, lavjags$summaries[xlocs,'SD'])
    lavpartable$prior <- c(lavpartable$prior, rep("", sum(nvars)))
    lavpartable$jlabel <- c(lavpartable$jlabel, xnms)
    lavpartable$logBF <- c(lavpartable$logBF, rep(NA, sum(nvars)))

    lavpartable
}
