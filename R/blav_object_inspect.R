## inspect blavaan object (wrapper around lavInspect with
## some additions)
blavTech <- function(blavobject, what, ...) {

    blavInspect(blavobject, what, ...)
}

## use lavInspect everywhere we can:
blavInspect <- function(blavobject, what, ...) {

    stopifnot(inherits(blavobject, "blavaan"))
  
    dotdotdot <- list(...)
    dotNames <- names(dotdotdot)
    add.labels <- TRUE
    if(any(dotNames == "add.labels")) add.labels <- dotdotdot$add.labels
  
    ## whats unique to blavaan
    blavwhats <- c("start", "starting.values", "inits", "psrf",
                   "ac.10", "neff", "mcmc", "draws", "samples",
                   "n.chains", "cp", "dp", "postmode", "postmean",
                   "postmedian", "hpd")

    ## whats that are not handled
    nowhats <- c("mi", "modindices", "modification.indices",
                 "wls.est", "wls.obs", "wls.v")

    if(what %in% blavwhats){
        idx <- blavobject@ParTable$jagpnum
        idx <- idx[!is.na(idx)]
        labs <- lav_partable_labels(blavobject@ParTable, type = "free")
        if(what %in% c("start", "starting.values", "inits")){
            blavobject@external$runjags$inits
        } else if(what %in% c("psrf", "ac.10", "neff")){
            if(what == "psrf"){
                OUT <- blavobject@external$runjags$summaries[idx,'psrf']
                ## blavobject@ParTable$psrf[!is.na(blavobject@ParTable$psrf)]
            }else if(what == "ac.10"){
                OUT <- blavobject@external$runjags$summaries[idx,'AC.10']
            } else{
                OUT <- blavobject@external$runjags$summaries[idx,'SSeff']
            }
            if(add.labels) names(OUT) <- labs
            OUT
        } else if(what %in% c("mcmc", "draws", "samples", "hpd")){
            draws <- blavobject@external$runjags$mcmc
            draws <- lapply(draws, function(x) x[,idx])
            draws <- mcmc.list(draws)
            if(what == "hpd"){
                pct <- .95
                if("level" %in% names(dotdotdot)) pct <- dotdotdot$level
                draws <- mcmc(do.call("rbind", draws))
                draws <- HPDinterval(draws, pct)
                if(add.labels) rownames(draws) <- labs
            }
            draws
        } else if(what == "n.chains"){
            length(blavobject@external$runjags$mcmc)
        } else if(what == "cp"){
            blavobject@Options$cp
        } else if(what == "dp"){
            blavobject@Options$dp
        } else if(what %in% c("postmode", "postmean", "postmedian")){
            if(what == "postmean"){
                OUT <- blavobject@external$runjags$summaries[idx,'Mean']
            }else if(what == "postmedian"){
                OUT <- blavobject@external$runjags$summaries[idx,'Median']
            } else{
                OUT <- blavobject@external$runjags$summaries[idx,'Mode']
            }
            if(add.labels) names(OUT) <- labs
            OUT
        }
    } else if(what %in% nowhats){
        stop(paste("blavaan ERROR: argument", what,
                   "not available for Bayesian models."))
    } else {
        ## we can use lavInspect
        lavargs <- c(dotdotdot, list(lavobject = blavobject, what = what))
        do.call("lavInspect", lavargs)
    }
}
