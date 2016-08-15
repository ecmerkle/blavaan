postpred <- function(lavpartable, lavmodel, lavoptions, 
                     lavsamplestats, lavdata, lavcache, lavjags,
                     measure = "logl") {

    ## run through lavjags$mcmc, generate data from various posterior
    ## samples. thin like we do in samp_lls

    samp.indices <- sampnums(lavjags, thin=5)
    n.chains <- length(lavjags$mcmc)
    psamp <- length(samp.indices)
  
    ## parallel across chains if we can
    ncores <- NA
    loop.comm <- "lapply"
    if(requireNamespace("parallel", quietly = TRUE)){
      ncores <- min(n.chains, parallel::detectCores())
      loop.comm <- "mclapply"
    }
  
    origlavmodel <- lavmodel
    origlavdata <- lavdata

    loop.args <- list(X = 1:n.chains, FUN = function(j){
      ind <- csdist <- rep(NA, psamp)
      for(i in 1:psamp){
        ## translate each posterior sample to a model-implied mean vector +
        ## cov matrix.
        lavmodel <- fill_params(lavjags$mcmc[[j]][samp.indices[i],],
                                origlavmodel, lavpartable)

        ## generate data (some code from lav_bootstrap.R)
        implied <- lav_model_implied(lavmodel)
        Sigma.hat <- implied$cov
        Mu.hat <- implied$mean
        dataeXo <- lavdata@eXo

        ## TODO? this generates complete cases; maybe we want missing
        ## observations to stay missing in the generated data:
        dataX <- vector("list", length=lavdata@ngroups)
        for(g in 1:lavsamplestats@ngroups) {
          dataX[[g]] <- MASS::mvrnorm(n     = lavsamplestats@nobs[[g]],
                                      Sigma = Sigma.hat[[g]],
                                      mu    = Mu.hat[[g]])
          dataX[[g]][is.na(origlavdata@X[[g]])] <- NA
        }

        ## compute (i) X2 of generated data and model-implied
        ## moments, along with (ii) X2 of real data and model-implied
        ## moments.
        chisq.obs <- -2*(lavjags$samplls[i, j, 1] -
                         lavjags$samplls[i, j, 2])
                             #get_ll(lavmodel = lavmodel,
                             #    lavpartable = lavpartable,
                             #    lavsamplestats = lavsamplestats,
                             #    lavoptions = lavoptions,
                             #    lavcache = lavcache,
                             #    lavdata = origlavdata,
                             #    measure = measure)

        ## check for missing, to see if we can easily get baseline ll for chisq
        mis <- FALSE
        if(any(is.na(unlist(lavdata@X)))) mis <- TRUE

        if(!mis){
          lavdata@X <- dataX
          
          chisq.boot <- 2*diff(get_ll(lavmodel = lavmodel,
                                      lavdata = lavdata,
                                      measure = measure))
        } else {
          ## we need lavaan to get the saturated log-l for missing data (EM)
                                         
          # YR: ugly hack to avoid lav_samplestats_from_data:
          # reconstruct data + call lavaan()
          # ed: if we need lavaan() anyway, might as well
          # get the chisq while we're here:
          DATA.X <- do.call("rbind", dataX)
          colnames(DATA.X) <- lavdata@ov.names[[1L]]
          DATA.eXo <- do.call("rbind", dataeXo)
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
          lavmodel2 <- lavmodel
          lavmodel2@control <- list(optim.method="none")
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
      } # i
      list(ind = ind, csdist = csdist)
    })

    if(loop.comm == "mclapply"){
        loop.args <- c(loop.args, list(mc.cores = ncores))
        res <- do.call(parallel::mclapply, loop.args)
    } else {
        res <- do.call(lapply, loop.args)
    }

    ind <- unlist(lapply(res, function(x) x$ind))
    csdist <- unlist(lapply(res, function(x) x$csdist))

    ppval <- mean(as.numeric(ind))
    cspi <- quantile(as.numeric(csdist), c(.025,.975))
    
    list(ppval=ppval, cspi=cspi)
}
