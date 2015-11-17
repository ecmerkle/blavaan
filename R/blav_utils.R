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
    for(i in 2:length(mcmc)){
        fullmat <- rbind(fullmat, mcmc[[i]])
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
