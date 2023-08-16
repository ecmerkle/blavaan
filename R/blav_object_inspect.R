## inspect blavaan object (wrapper around lavInspect with
## some additions)
blavTech <- function(blavobject, what, ...) {

    blavInspect(blavobject, what, ...)
}

## use lavInspect everywhere we can:
blavInspect <- function(blavobject, what, ...) {

    stopifnot(inherits(blavobject, "blavaan"))

    what <- tolower(what)
  
    dotdotdot <- list(...)
    dotNames <- names(dotdotdot)
    add.labels <- TRUE
    if(any(dotNames == "add.labels")) add.labels <- dotdotdot$add.labels

    jagtarget <- lavInspect(blavobject, "options")$target == "jags"
  
    ## whats unique to blavaan
    blavwhats <- c("start", "starting.values", "inits", "psrf",
                   "ac.10", "neff", "mcmc", "draws", "samples",
                   "n.chains", "cp", "dp", "postmode", "postmean",
                   "postmedian", "hpd", "jagnames", "stannames",
                   "fscores", "lvs", "fsmeans", "lvmeans", "mcobj",
                   "rhat", "n_eff", "nchain", "nchains")

    ## blavwhats that don't require do.fit
    blavnofit <- c("start", "starting.values", "inits", "n.chains", "cp", "dp",
                   "jagnames", "stannames", "nchain", "nchains")
  
    ## whats that are not handled
    nowhats <- c("mi", "modindices", "modification.indices",
                 "wls.est", "wls.obs", "wls.v")

    if(what %in% blavwhats){
        if(!(what %in% blavnofit) & !blavobject@Options$do.fit){
            stop(paste0("blavaan ERROR: ", what, " does not exist when do.fit = FALSE"))
        }
        if(jagtarget){
            idx <- blavobject@ParTable$jagpnum
            idx <- idx[!is.na(idx)]
        } else {
            idx <- blavobject@ParTable$stansumnum
            if("pxnames" %in% names(blavobject@ParTable)){
              drows <- grepl("^def", blavobject@ParTable$pxnames)
            } else {
              drows <- grepl("def", blavobject@ParTable$mat)
            }
            idx <- idx[blavobject@ParTable$free > 0 | drows]
        }
        labs <- lav_partable_labels(blavobject@ParTable, type = "free")
        if(what %in% c("start", "starting.values", "inits")){
            blavobject@external$inits
        } else if(what %in% c("psrf", "ac.10", "neff", "rhat", "n_eff")){
            if(jagtarget){
                mcmcsumm <- blavobject@external$mcmcout$summaries
            } else {
                mcmcsumm <- rstan::summary(blavobject@external$mcmcout)$summary
            }
            if(what %in% c("psrf", "rhat")){
                if(jagtarget){
                    OUT <- mcmcsumm[idx,'psrf']
                } else {
                    OUT <- mcmcsumm[idx,'Rhat']
                }
            }else if(what == "ac.10"){
                if(jagtarget){
                    OUT <- mcmcsumm[idx,'AC.10']
                } else {
                    stop("blavaan ERROR: autocorrelation stat currently unavailable for Stan.")
                }
            } else {
                if(jagtarget){
                    OUT <- mcmcsumm[idx,'SSeff']
                } else {
                    OUT <- mcmcsumm[idx,'n_eff']
                }
            }
            if(add.labels) names(OUT) <- labs
            OUT
        } else if(what %in% c("mcmc", "draws", "samples", "hpd")){
            ## add defined parameters to labels
            pt <- blavobject@ParTable
            pt$free[pt$op == ":="] <- max(pt$free, na.rm = TRUE) + 1:sum(pt$op == ":=")
            labs <- lav_partable_labels(pt, type = "free")
            draws <- make_mcmc(blavobject@external$mcmcout)
            draws <- lapply(draws, function(x) mcmc(x[,idx]))
            draws <- mcmc.list(draws)
            if(what == "hpd"){
                pct <- .95
                if("level" %in% dotNames) pct <- dotdotdot$level
                draws <- mcmc(do.call("rbind", draws))
                draws <- HPDinterval(draws, pct)
                if(add.labels) rownames(draws) <- labs
            }
            draws
        } else if(what == "mcobj"){
            blavobject@external$mcmcout
        } else if(what %in% c("fscores","lvs","fsmeans","lvmeans")){
            if(jagtarget){
                etas <- any(blavobject@external$mcmcout$monitor == "eta")
            } else {
                etas <- any(grepl("^eta", rownames(blavobject@external$stansumm)))
            }

            ## how many lvs, excluding phantoms
            lvmn <- lavInspect(blavobject, "mean.lv")
            if(!inherits(lvmn, "list")){
                lvmn <- list(lvmn)
            }
            nlv <- length(lvmn[[1]])
            nlv2 <- 0
            if(length(lvmn) > 1) nlv2 <- length(lvmn[[2]])

            if(nlv == 0) stop("blavaan ERROR: no latent variables are in the model")
            if(!etas) stop("blavaan ERROR: factor scores not saved; set save.lvs=TRUE")
            
            nsamp <- sum(lavInspect(blavobject, "nobs"))
            nclus <- lavInspect(blavobject, "nclusters")

            draws <- make_mcmc(blavobject@external$mcmcout, blavobject@external$stanlvs)
            drawcols <- grep("^eta\\[", colnames(draws[[1]]))
            drawcols2 <- grep("^eta_b", colnames(draws[[1]]))

            if(jagtarget){
                ## remove phantoms
                drawcols <- drawcols[1:(nlv * nsamp)]
            } else {
                nfound <- length(drawcols)/nsamp
                drawcols <- drawcols[as.numeric(matrix(1:length(drawcols),
                                                       nsamp, nfound,
                                                       byrow=TRUE)[,1:nlv])]

                if(any(nclus > 1)){
                  nfound2 <- length(drawcols2)/sum(nclus)
                  drawcols2 <- drawcols2[as.numeric(matrix(1:length(drawcols2),
                                                           sum(nclus), nfound2,
                                                           byrow=TRUE)[,1:nlv2])]
                  drawcols <- c(drawcols, drawcols2)
                }
            }
            draws <- lapply(draws, function(x) mcmc(x[,drawcols]))

            ## for target="stan" + missing, use @Data@Mp to reorder rows to correspond
            ## to original data
            mis <- any(is.na(unlist(blavobject@Data@X)))
            Mp <- blavobject@Data@Mp
            if(blavobject@Options$target == "stan" & mis){
                rorig <- sapply(Mp, function(x) unlist(x$case.idx))
                empties <- sapply(Mp, function(x) x$empty.idx)
                if(inherits(rorig, "list")){
                    ## multiple groups
                    for(ii in length(rorig)){
                        rorig[[ii]] <- blavobject@Data@case.idx[[ii]][rorig[[ii]]]
                    }
                    rorig <- unlist(rorig)
                }
                cids <- Mp2dataidx(Mp, blavobject@Data@case.idx)

                ## reordering for lvs:
                nfit <- sum(lavInspect(blavobject, 'nobs'))
                rsamps <- rep(NA, nlv*nfit)
                for(j in 1:nlv){
                    rsamps[((j-1)*nfit + 1):(j*nfit)] <- (j-1)*nfit + cids
                }

                for(j in 1:length(draws)){
                    draws[[j]][,rsamps] <- draws[[j]]
                }
            }
            draws <- mcmc.list(draws)

            if((what %in% c("lvmeans", "fsmeans")) | ("means" %in% dotdotdot)){
                br <- TRUE
                if(jagtarget){
                    summ <- blavobject@external$mcmcout$summaries
                    summname <- "Mean"
                    br <- FALSE
                } else {
                    summ <- blavobject@external$stansumm
                    summname <- "mean"
                }
                mnrows <- grep("^eta\\[", rownames(summ))
                mnrows2 <- grep("^eta_b", rownames(summ))

                draws <- matrix(summ[mnrows,summname], nsamp,
                                length(mnrows)/nsamp, byrow=br)[,1:nlv,drop=FALSE]
                ## FIXME multiple groups?
                colnames(draws) <- names(lvmn[[1]])

                if(any(nclus > 1)){
                  draws2 <- matrix(summ[mnrows2,summname], sum(nclus),
                                   length(mnrows2)/sum(nclus), byrow=br)[,1:nlv2,drop=FALSE]
                  colnames(draws2) <- names(lvmn[[2]])

                  draws <- list(draws, draws2)
                }

                if(blavobject@Options$target == "stan" & mis){
                    draws[rank(rorig),] <- draws
                }
            }
            draws
        } else if(what %in% c("n.chains", "nchain", "nchains")){
            draws <- make_mcmc(blavobject@external$mcmcout)
            length(draws)
        } else if(what == "cp"){
            blavobject@Options$cp
        } else if(what == "dp"){
            blavobject@Options$dp
        } else if(what %in% c("postmode", "postmean", "postmedian")){
            if(jagtarget){
                mcmcsumm <- blavobject@external$mcmcout$summaries
            } else {
                mcmcsumm <- rstan::summary(blavobject@external$mcmcout)$summary
            }

            if(what == "postmean"){
                if(jagtarget){
                    OUT <- mcmcsumm[idx,'Mean']
                } else {
                    OUT <- mcmcsumm[idx,'mean']
                }
            }else if(what == "postmedian"){
                if(jagtarget){
                    OUT <- mcmcsumm[idx,'Median']
                } else {
                    OUT <- mcmcsumm[idx,'50%']
                }
            } else {
                if(jagtarget){
                    OUT <- mcmcsumm[idx,'Mode']
                } else {
                    stop("blavaan ERROR: Modes unavailable for Stan.")
                }
            }
            if(add.labels) names(OUT) <- labs
            OUT
        } else if(what == "jagnames"){
            if(!jagtarget) stop("blavaan ERROR: JAGS was not used for model estimation.")
            OUT <- blavobject@ParTable$pxnames[blavobject@ParTable$free > 0]
            OUT <- OUT[order(blavobject@ParTable$free[blavobject@ParTable$free > 0])]
            if(add.labels) names(OUT) <- labs
            OUT
        } else if(what == "stannames"){
            if(jagtarget) stop("blavaan ERROR: Stan was not used for model estimation.")
            mcmcsumm <- rstan::summary(blavobject@external$mcmcout)$summary
            OUT <- rownames(mcmcsumm)[idx]
            if(add.labels) names(OUT) <- labs
            OUT
        }
    } else if(what %in% nowhats){
        stop(paste("blavaan ERROR: argument", what,
                   "not available for Bayesian models."))
    } else {
        ## we can use lavInspect
        lavargs <- c(dotdotdot, list(object = blavobject, what = what))
        do.call("lavInspect", lavargs)
    }
}
