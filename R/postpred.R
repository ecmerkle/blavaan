### Ed Merkle
### posterior predictive model checking
###   - returns PPP value to store in @test slot
###   - returns posterior predictive distribution (PP-Dist) for additional
###     discrepancy functions evaluated on observed and replicated data
### Last updated: 24 May 2018
### Mauricio Garnier-Villarreal:
###   - update to return "chisq" PP-Dist from observed and replicated data
### Terrence D. Jorgensen:
###   - added "discFUN" argument to evaluate custom discrepancy function(s)
###     on observed and replicated data.  This is distinct from the "measure"
###     argument, which only returns values from fitMeasures().
###   - made notes with "FIXME TDJ" in places where postpred() could be updated

postpred <- function(lavpartable, lavmodel, lavoptions, 
                     lavsamplestats, lavdata, lavcache, lavjags,
                     samplls, measure = "logl", thin = 1, discFUN = NULL) {
    ## check custom discrepancy function(s)
    if (!is.null(discFUN)) {
      allFuncs <- if (is.list(discFUN)) all(sapply(discFUN, is.function)) else FALSE
      if (!(is.function(discFUN) || allFuncs)) stop('The "discFUN" argument must',
                                                    ' be a (list of) function(s).')
    }
    discFUN <- NULL #FIXME: Not implemented yet
  
    ## run through lavjags$mcmc, generate data from various posterior
    ## samples. thin like we do in samp_lls
    lavmcmc <- make_mcmc(lavjags)
    samp.indices <- sampnums(lavjags, thin=thin)
    n.chains <- length(lavmcmc)
    psamp <- length(samp.indices)
    
    ## parallel across chains if we can
    ncores <- NA
    loop.comm <- "lapply"
    if(.Platform$OS.type != "windows" & requireNamespace("parallel", quietly = TRUE)){
      ncores <- min(n.chains, parallel::detectCores())
      loop.comm <- "mclapply"
    }
    
    origlavmodel <- lavmodel
    origlavdata <- lavdata

    ## check for missing, to see if we can easily get baseline ll for chisq
    mis <- FALSE
    if(any(is.na(unlist(lavdata@X)))) mis <- TRUE
  
    loop.args <- list(X = 1:n.chains, FUN = function(j){
      ind <- csdist <- csboots <- rep(NA, psamp)
      for(i in 1:psamp){
        ## translate each posterior sample to a model-implied mean vector +
        ## cov matrix.
        lavmodel <- fill_params(lavmcmc[[j]][samp.indices[i],],
                                origlavmodel, lavpartable)

        ## generate data (some code from lav_bootstrap.R)
        implied <- lav_model_implied(lavmodel)
        Sigma.hat <- implied$cov
        Mu.hat <- implied$mean
        dataeXo <- lavdata@eXo

        dataX <- origlavdata@X
        for(g in 1:lavsamplestats@ngroups) {
          x.idx <- lavsamplestats@x.idx[[g]]
          nox <- (1:nrow(Mu.hat[[g]]))[-x.idx]
          if(!is.null(x.idx) && length(x.idx) > 0L){
            ## for fixed.x, generate the other ovs
            ## conditional on the x values
            if(!mis){
              tm1 <- Sigma.hat[[g]][nox,x.idx] %*% solve(Sigma.hat[[g]][x.idx,x.idx])
              cmu <- Mu.hat[[g]][nox,] +
                tm1 %*% apply(origlavdata@X[[g]][,x.idx,drop=FALSE], 1,
                              function(x) (x - Mu.hat[[g]][x.idx,]))
              csig <- Sigma.hat[[g]][nox,nox] - tm1 %*% Sigma.hat[[g]][x.idx,nox]
              sigchol <- chol(csig)
            
              dataX[[g]][,nox] <- t(apply(cmu, 2, function(x) mnormt::rmnorm(n=1,
                                                                      sqrt=sigchol,
                                                                      mean=x)))
            } else {
              ## condition only on observed x values; to be re-examined after
              ## lav_mvnorm_missing_loglik_samplestats handles x.idx
              M <- lavsamplestats@missing[[g]]
              Mp <- lavdata@Mp[[g]]
              for(p in 1:length(M)){
                var.idx <- M[[p]][["var.idx"]]
                obsx <- x.idx[var.idx[x.idx]]

                ## could also generate missing x's, but has no
                ## impact
                ##misx <- x.idx[!var.idx[x.idx]]
                ##if(length(misx) > 0){
                ##  dataX[[g]][Mp$case.idx[[p]],misx] <- mnormt::rmnorm(n = M[[p]]$freq,
                ##                                         varcov = Sigma.hat[[g]][misx,misx],
                ##                                         mean = Mu.hat[[g]][misx,])
                ##  obsx <- x.idx
                ##}
                  
                if(length(obsx) > 0){
                  xp.idx <- obsx
                  tm1 <- Sigma.hat[[g]][nox,xp.idx] %*% solve(Sigma.hat[[g]][xp.idx,xp.idx])
                  cmu <- Mu.hat[[g]][nox,] +
                    tm1 %*% apply(origlavdata@X[[g]][Mp$case.idx[[p]],xp.idx,drop=FALSE], 1,
                                  function(x) (x - Mu.hat[[g]][xp.idx,]))
                  csig <- Sigma.hat[[g]][nox,nox] - tm1 %*% Sigma.hat[[g]][xp.idx,nox]
                  sigchol <- chol(csig)

                  dataX[[g]][Mp$case.idx[[p]],nox] <- t(apply(cmu, 2,
                                                              function(x) mnormt::rmnorm(n=1,
                                                                                         sqrt=sigchol,
                                                                                         mean=x)))
                } else {
                  cmu <- Mu.hat[[g]][nox,]
                  csig <- Sigma.hat[[g]][nox,nox]

                  dataX[[g]][Mp$case.idx[[p]],nox] <- as.matrix(mnormt::rmnorm(n = M[[p]]$freq,
                                                                               varcov = csig,
                                                                               mean = cmu))
                }
              }
            }
          } else {
            nox <- 1:nrow(Mu.hat[[g]])
            cmu <- Mu.hat[[g]]
            csig <- Sigma.hat[[g]]

            dataX[[g]] <- as.matrix(mnormt::rmnorm(n = nrow(dataX[[g]]),
                                    varcov = csig, mean = cmu))
          }

          dataX[[g]][is.na(origlavdata@X[[g]])] <- NA

          ## get rid of completely missing
          if(length(origlavdata@Mp[[g]]$empty.idx) > 0){
            dataX[[g]] <- dataX[[g]][-origlavdata@Mp[[g]]$empty.idx,,drop=FALSE]
          }
        }
        
        ## compute (i) X2 of generated data and model-implied
        ## moments, along with (ii) X2 of real data and model-implied
        ## moments.
        chisq.obs <- -2*(samplls[i, j, 1] -
                         samplls[i, j, 2])
        
        #FIXME TDJ: Apply custom "discFUN" here.
        #           Need to create a lavaan object using original data,
        #           like the hack below does for generated data.
  
        if(!mis & length(discFUN) == 0){ #TDJ: if discFUN is supplied, we go right to "else"
          lavdata@X <- dataX
          chisq.boot <- 2*diff(get_ll(lavmodel = lavmodel,
                                      lavsamplestats = lavsamplestats,
                                      lavdata = lavdata,
                                      measure = measure))
          #FIXME TDJ: no way to apply custom "discFUN" here. Use hack below?
        } else {
          ## we need lavaan to get the saturated log-l for missing data (EM)
                                         
          # YR: ugly hack to avoid lav_samplestats_from_data:
          # reconstruct data + call lavaan()
          # ed: if we need lavaan() anyway, might as well
          #     get the chisq while we're here:
          # TDJ: this also enables us to apply custom "discFUN" argument to
          #     fitted lavaan object -- also use this hack when !is.null(discFUN)?
          #     ed: yes, probably want to pull necessary stuff from "out" object below
          DATA.X <- do.call("rbind", dataX)
          colnames(DATA.X) <- lavdata@ov.names[[1L]]
          DATA.eXo <- do.call("rbind", dataeXo)
          empties <- any(sapply(lavdata@Mp, function(x) length(x$empty.idx)) > 0)
          if(empties){
            empties <- as.numeric(unlist(sapply(lavdata@Mp, function(x) x$empty.idx)))
            if(!any(is.na(empties))){
              DATA.eXo <- DATA.eXo[-empties, , drop = FALSE]
            }
          }
          if(!is.null(DATA.eXo)) {
            colnames(DATA.eXo) <- lavdata@ov.names.x[[1L]]
            DATA <- cbind(DATA.X, DATA.eXo)
          } else {
            DATA <- DATA.X
          }
          DATA <- as.data.frame(DATA)

          lavoptions2 <- lavoptions
          lavoptions2$verbose <- FALSE
          lavoptions2$estimator <- "ML"
          lavoptions2$se <- "none"
          lavoptions2$test <- "standard"
          lavoptions2$optim.method <- "none"
          lavmodel2 <- lavmodel
          if("control" %in% slotNames(lavmodel2)){
            lavmodel2@control <- list(optim.method="none")
          }
          if(lavsamplestats@ngroups > 1L) {
            DATA$.g. <- rep(1:lavdata@ngroups, 
                            times = unlist(lavdata@nobs))
            out <- lavaan(slotOptions = lavoptions2, 
                          slotParTable = lavpartable,
                          slotSampleStats = NULL, slotData = NULL, 
                          slotModel = lavmodel2, slotCache = lavcache, 
                          data = DATA, group = ".g.")
          } else {
            out <- lavaan(slotOptions = lavoptions2, 
                          slotParTable = lavpartable,
                          slotSampleStats = NULL, slotData = NULL, 
                          slotModel = lavmodel2, slotCache = lavcache, 
                          data = DATA)
          }
          # bootSampleStats <- out@SampleStats
          # end of ugly hack
  
          if(measure %in% c("logl", "chisq")){
            chisq.boot <- fitMeasures(out, "chisq")
          } else {
            chisq.boot <- fitMeasures(out, measure)
          }
  
          ## see lines 286-298 of lav_bootstrap to avoid fixed.x errors?
          ## chisq.boot <- 2*diff(get_ll(lavmodel = lavmodel,
          ##                             lavpartable = lavpartable,
          ##                             lavsamplestats = bootSampleStats,
          ##                             lavoptions = lavoptions,
          ##                             lavcache = lavcache,
          ##                             lavdata = lavdata,
          ##                             measure = measure))
        }
        ## record whether observed value is larger
        ind[i] <- chisq.obs < chisq.boot
        csdist[i] <- chisq.obs
        csboots[i] <- chisq.boot
        #FIXME TDJ: extract and organize custom "discFUN" output here
        
      } # i
      
      result <- list(ind = ind, csdist = csdist, csboots = csboots)
      # if (!is.null(discFUN)) result <- c(result, discFUN_results)
      result
    })
  
    if(loop.comm == "mclapply"){
        loop.args <- c(loop.args, list(mc.cores = ncores))
        res <- do.call(parallel::mclapply, loop.args)
    } else {
        res <- do.call(lapply, loop.args)
    }
  
    ind <- unlist(lapply(res, function(x) x$ind))
    csdist <- unlist(lapply(res, function(x) x$csdist))
    csboots <- unlist(lapply(res, function(x) x$csboots))
    #FIXME TDJ: extract custom "discFUN" output here
    
    ppval <- mean(as.numeric(ind))
    cspi <- quantile(as.numeric(csdist), c(.025,.975))
    
    #FIXME TDJ: check whether to add custom "discFUN" output to returned list
    list(ppval = ppval, cspi = cspi, chisqs = cbind(obs = csdist, reps = csboots))
}


