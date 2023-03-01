### posterior predictive data generation and model checking
###   - returns PPP value to store in @test slot
###   - returns posterior predictive distribution (PP-Dist) for additional
###     discrepancy functions evaluated on observed and replicated data

### authors:
### Ed Merkle
###   - data generation and ppp computations
### Yves Rosseel
###   - "dirty hack" to manufacture a lavaan object from 1 posterior sample
### Mauricio Garnier-Villarreal:
###   - update to return "chisq" PP-Dist from observed and replicated data
### Terrence D. Jorgensen:
###   - added "discFUN" argument to evaluate custom discrepancy function(s)
###     on observed and replicated data.  This is distinct from the "measure"
###     argument, which only returns values from fitMeasures().
###   - added "probs" argument to optionally return different quantiles, and
###     changed the returned name from "cspi" to "quantiles" (now a named list)
###   - updated if-else logic to first check for custom discrepancy (only used
###     in public ppmc() function), then check for complete-data logl/chisq,
###     otherwise use hack for (in)complete data using ANY fit.measures from
###     fitMeasures().  Note that the latter could also be a custom discFUN.
###   - for custom discFUN, updated simulated data to retain grouping variable
###     name and the group labels

postpred <- function(samplls = NULL, lavobject = NULL, measure = "logl", thin = 1,
                     discFUN = NULL, probs = c(.025, .975)) {

  lavpartable <- lavobject@ParTable
  lavmodel <- lavobject@Model
  lavoptions <- lavobject@Options
  catmod <- lavInspect(lavobject, "categorical")
  lavsamplestats <- lavobject@SampleStats
  lavdata <- lavobject@Data
  lavcache <- lavobject@Cache
  lavjags <- lavobject@external$mcmcout
  
  ## check custom discrepancy function(s)
  if (!is.null(discFUN)) {
    if (!is.list(discFUN)) discFUN <- list(discFUN)
    if (!all(sapply(discFUN, is.function)))
      stop('blavaan ERROR: The "discFUN" argument must be a (list of) function(s).')
  }

  ## verify that probs is a numeric vector
  probs <- as.numeric(probs)

  ## add sampled lvs to posterior samples, if they exist
  ## and add stan data, too
  save_lvs <- FALSE
  standata <- NULL
  if (!is.null(lavobject)) {    
    save_lvs <- ("stanlvs" %in% names(lavobject@external)) & length(discFUN)
    standata <- lavobject@external$mcmcdata
  }

  if (save_lvs) {
    lavmcmc <- make_mcmc(lavjags, lavobject@external$stanlvs)
  } else {
    lavmcmc <- make_mcmc(lavjags)
  }
  n.chains <- length(lavmcmc)
  samp.indices <- sampnums(lavjags, thin=thin)
  psamp <- length(samp.indices)
  lmorig <- lavmodel
  ldorig <- lavdata

  ## check for missing, to see if we can easily get baseline ll for chisq
  mis <- FALSE
  if (any(is.na(unlist(lavdata@X)))) mis <- TRUE

  loop.args <- list(X = 1:psamp, FUN = function(i){
    ## lists to store output (reduced after concatenated across chains)
    if (length(discFUN)) {
      ## Nested lists (iterations within discrepancy functions)
      ind <- csdist <- csboots <- vector("list", length(discFUN))
      for (d in seq_along(discFUN)) {
        ind[[d]] <- csdist[[d]] <- csboots[[d]] <- vector("list", n.chains)
      }

    } else {
      ## scalar or vecctor from fitMeasures()
      ind <- csdist <- csboots <- vector("list", n.chains)
    }

    for(j in 1:n.chains){
      ## compute (i) X2 of generated data and model-implied moments,
      ## along with (ii) X2 of real data and model-implied moments.
      fit <- lavobject

      ## REAL DATA
      if (length(discFUN)) {
        ## update options
        lavoptions1 <- lavoptions
        lavoptions1$verbose <- FALSE
        lavoptions1$estimator <- "ML"
        lavoptions1$se <- "none"
        lavoptions1$test <- "standard"
        lavoptions1$optim.method <- "none"
        lavmodel1 <- fill_params(lavmcmc[[j]][i,], # update parameters
                                 lmorig, lavpartable)
        if ("control" %in% slotNames(lavmodel1)) {
          lavmodel1@control <- list(optim.method = "none")
        }
        fit <- lavaan(slotOptions = lavoptions1, slotModel = lavmodel1, # updated slots
                      slotSampleStats = lavsamplestats, slotData = ldorig,
                      slotParTable = lavpartable, slotCache = lavcache)

        if (save_lvs) {
          eta <- fill_eta(lavmcmc[[j]][i,], lavmodel1, lavpartable, lavsamplestats, lavdata)
          implied <- cond_moments(lavmcmc[[j]][i,], lavmodel1, lavpartable, lavsamplestats,
                                  lavdata, fit)
          fit@external$eta <- eta
          fit@external$implied <- implied
        }
        
        ## Apply custom "discFUN" to observed data
        chisq.obs <- vector("list", length(discFUN))
        for (d in seq_along(discFUN)) {
          chisq.obs[[d]] <- discFUN[[d]](fit)
          if (mode(chisq.obs[[d]]) != "numeric")
            stop('Each function in discFUN must return an object with numeric ',
                 'mode (e.g., vector, matrix, or array)')
        }

      } else if (measure[1] %in% c("logl", "chisq") & length(measure) == 1) {
        ## standard chi-squared-based PPP
        chisq.obs <- -2*(samplls[samp.indices[i], j, 1] -
                         samplls[samp.indices[i], j, 2])

      } else {
        ## save alternative fit.measures=measure for observed data
        chisq.obs <- get_ll(postsamp = lavmcmc[[j]][i,], lavobject = lavobject,
                            measure = measure, standata = standata)
      }

      ## supply extra args to postdata so that we only generate
      ## a single dataset
      dataX <- postdata(samp.indices = samp.indices[i],
                        chain.num = j,
                        lavmodel = lmorig, lavdata = ldorig,
                        lavjags = lavjags, lavpartable = lavpartable,
                        lavsamplestats = lavsamplestats,
                        lavobject = fit,
                        conditional = save_lvs)
      lavmodel <- dataX$lavmod[[1]]
      dataX <- dataX[[1]]
      dataeXo <- lavdata@eXo

      simstandata <- NULL
      if(!is.null(standata)) {
        fit@Data@X <- dataX
        simstandata <- lav2standata(fit)
      }

      ## SIMULATED DATA
      if (!mis & !length(discFUN) & measure[1] %in% c("logl", "chisq") & length(measure) == 1) {
        lavobject@Data@X <- dataX
        chisq.boot <- 2*diff(get_ll(lavobject = lavobject,
                                    measure = measure[1], standata = simstandata))

      ## chi-squared for (in)complete data
      } else if (!length(discFUN) & measure[1] %in% c("logl", "chisq") & length(measure) == 1) {
        lavoptions$target <- "jags" ## because we need both h0 and h1 ll
        lavobject@Data@X <- dataX
        chisq.boot <- 2*diff(get_ll(lavobject = lavobject,
                                    measure = measure[1], standata = simstandata))
      } else {
        ## any other measure/discFUN

        # YR: ugly hack to avoid lav_samplestats_from_data:
        # reconstruct data + call lavaan()
        # TDJ: this also enables us to apply custom "discFUN" argument to
        #     fitted lavaan object -- also use this hack when !is.null(discFUN)?
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
        if(catmod) lavoptions2$estimator <- "DWLS"
        if(!lavoptions2$meanstructure) lavoptions2$missing <- 'listwise' ## to avoid "missing" ll calculations, possibly from a bug in lavaan 0.6-13
        lavoptions2$se <- "none"
        lavoptions2$test <- "standard"
        lavoptions2$optim.method <- "none"
        lavmodel2 <- lavmodel
        if("control" %in% slotNames(lavmodel2)){
          lavmodel2@control <- list(optim.method="none")
        }
        if(lavsamplestats@ngroups > 1L) {
          DATA[ , lavdata@group] <- rep(lavdata@group.label,
                                        times = unlist(lavdata@nobs))
          out <- lavaan(slotOptions = lavoptions2,
                        slotParTable = lavpartable,
                        slotSampleStats = NULL, slotData = NULL,
                        slotModel = lavmodel2, slotCache = lavcache,
                        data = DATA, group = lavdata@group)
        } else {
          out <- lavaan(slotOptions = lavoptions2,
                        slotParTable = lavpartable,
                        slotSampleStats = NULL, slotData = NULL,
                        slotModel = lavmodel2, slotCache = lavcache,
                        data = DATA)
        }
        # bootSampleStats <- out@SampleStats
        if (save_lvs) {
          out@external$eta <- eta
          out@external$implied <- implied
        }
        # end of ugly hack

        ## prioritize custom discrepancy function, so "logl" can remain default
        if (length(discFUN)) {
          #TODO TDJ: Apply custom "discFUN" to simulated data
          chisq.boot <- vector("list", length(discFUN))
          for (d in seq_along(discFUN)) {
            chisq.boot[[d]] <- discFUN[[d]](out)
            # if (mode(chisq.boot[[d]]) != "numeric") # already checked for chisq.obs
            #   stop('Each function in discFUN must return an object with numeric ',
            #        'mode (e.g., vector, matrix, or array)')
          }

        } else if (measure[1] %in% c("logl", "chisq") & length(measure) == 1) {
          ## standard chi-squared-based PPP
          chisq.boot <- fitMeasures(out, "chisq")

        } else {
          ## save alternative fit.measures=measure for simulated data
          chisq.boot <- fitMeasures(out, measure)
        }

        ## see lines 286-298 of lav_bootstrap to avoid fixed.x errors?
        ## probably only needed for missing='ml.x'?
        ## chisq.boot <- 2*diff(get_ll(lavmodel = lavmodel,
        ##                             lavpartable = lavpartable,
        ##                             lavsamplestats = bootSampleStats,
        ##                             lavoptions = lavoptions,
        ##                             lavcache = lavcache,
        ##                             lavdata = lavdata,
        ##                             measure = measure))
      }


      ## COMPARE REAL v. SIMULATED, record whether observed value is larger
      if (length(discFUN)) {
        ## loop over each "discFUN" to extract output in a list
        for (d in seq_along(discFUN)) {
          csdist[[d]][[j]] <- chisq.obs[[d]]
          csboots[[d]][[j]] <- chisq.boot[[d]]
          ind[[d]][[j]] <- chisq.obs[[d]] < chisq.boot[[d]]
          ## pass long names, etc.
          attributes(ind[[d]][[j]]) <- attributes(chisq.obs[[d]])
        }

      } else {
        ## either a scalar or vector from fitMeasures()
        ind[[j]] <- chisq.obs < chisq.boot
        ## pass long names, etc.
        if (length(chisq.obs) > 1L) attributes(ind[[j]]) <- attributes(chisq.obs)
        csdist[[j]] <- chisq.obs
        csboots[[j]] <- chisq.boot
      }

    } # j

    result <- list(ind = ind, csdist = csdist, csboots = csboots)
    result
  }, future.seed = TRUE)

  res <- do.call("future_lapply", loop.args)
  
  ## extract PPP and posterior (realized & predictive) distributions
  if (length(discFUN)) {
    ## store in list per discrepancy function
    ind <- csdist <- csboots <- ppval <- vector("list", length(discFUN))

    for (d in seq_along(discFUN)) {
      ## concatenate lists from each chain
      ind[[d]]     <- do.call("c", lapply(res, function(x) x$ind[[d]]     ))
      csdist[[d]]  <- do.call("c", lapply(res, function(x) x$csdist[[d]]  ))
      csboots[[d]] <- do.call("c", lapply(res, function(x) x$csboots[[d]] ))
      ## mean for any k-dim array;
      ## FIXME? Reduce() is very slow here 
      ppval[[d]]   <- Reduce("+", ind[[d]]) / length(ind[[d]])
      attributes(ppval[[d]]) <- attributes(res[[1]]$ind[[d]][[1]])
      ## Assign names, if they exist
      if (!is.null(names(discFUN))) {
        names(ppval)   <- names(discFUN)
        names(csdist)  <- names(discFUN)
        names(csboots) <- names(discFUN)
      }
    }

  } else {
    ## either scalar or vector from fitMeasures()
    ## was originally based on a list of length n.chains,
    ## instead of psamp
    ind <- do.call("c", lapply(res, function(x) x$ind))
    csdist <- do.call("c", lapply(1:n.chains, function(j) lapply(res, function(x) x$csdist[[j]])))
    csboots <- do.call("c", lapply(1:n.chains, function(j) lapply(res, function(x) x$csboots[[j]])))

    ## mean vector
    ppval <- colMeans(do.call("rbind", ind))
    attributes(ppval) <- attributes(res[[1]]$ind[[1]])
  }

  list(ppval = ppval,
       ppdist = list(obs = csdist, reps = csboots))
}


## generate data from posterior predictive dist
postdata <- function(object = NULL, nrep = 50L, conditional = FALSE, type = "response", ...){

  ## parallel=TRUE can be sent here, but it slows everything down if
  ## used within postpred() because it leads to nested parallelization.
  ddd <- list(...)

  ## users can supply object; postpred() will supply the slots:
  if(length(object) == 0L){
    if(!all(c('lavmodel', 'lavdata', 'lavjags', 'lavpartable',
              'lavsamplestats') %in% names(ddd))){
      stop("blavaan ERROR: either object or lav-model/data/jags/partable/samplestats must be supplied.")
    }
    lavmodel <- ddd$lavmodel
    lavdata <- ddd$lavdata
    lavjags <- ddd$lavjags
    lavpartable <- ddd$lavpartable
    lavsamplestats <- ddd$lavsamplestats
    lavobject <- ddd$lavobject
  } else {
    lavmodel <- object@Model
    lavdata <- object@Data
    lavjags <- object@external$mcmcout
    lavpartable <- object@ParTable
    lavsamplestats <- object@SampleStats
    lavobject <- object
  }

  if(conditional){
    if(!(("eta" %in% names(lavobject@external) & "implied" %in% names(lavobject@external)))){
      stop("blavaan ERROR: could not find lvs.")
    }
  }
  
  lavmcmc <- make_mcmc(lavjags)
  n.chains <- length(lavmcmc)
  chnums <- 1:n.chains

  if(all(c("samp.indices", "chain.num") %in% names(ddd))){
    ## 1 sample, for postpred
    samp.indices <- ddd$samp.indices
    chnums <- ddd$chain.num
  } else {
    nper <- ceiling(nrep/n.chains)
    samp.indices <- sampnums(lavjags, thin=nper, lout=TRUE)
  }
  psamp <- length(samp.indices)

  origlavmodel <- lavmodel
  origlavdata <- lavdata

  ## check for missing, to see if we can easily generate new data
  mis <- FALSE
  if(any(is.na(unlist(lavdata@X)))) mis <- TRUE

  ## check for ordinal, which requires to chop the data
  ordresp <- FALSE
  if(length(lavdata@ordered) > 0 && type == "response"){
    ordresp <- TRUE
    ordidx <- which(lavdata@ov$type == "ordered")
    thidx <- lavmodel@th.idx
  }
  
  postdat <- vector("list", psamp)
  lavmod <- vector("list", psamp)

  loop.args <- list(X = chnums, FUN = function(j){
    ind <- csdist <- csboots <- rep(NA, psamp)
    for(i in 1:psamp){
      ## translate each posterior sample to a model-implied mean vector +
      ## cov matrix.
      lavmodel <- fill_params(lavmcmc[[j]][samp.indices[i],],
                              origlavmodel, lavpartable)
      lavmod[[i]] <- lavmodel

      ## generate data (some code from lav_bootstrap.R)
      if(conditional){
        implied <- lavobject@external$implied
      } else {
        implied <- lav_model_implied(lavmodel, delta = (lavmodel@parameterization == "delta"))
      }
      Sigma.hat <- implied$cov
      Mu.hat <- implied$mean
      if(is.null(Mu.hat[[1]])) Mu.hat <- lapply(lavsamplestats@mean, matrix)
      dataeXo <- lavdata@eXo

      dataX <- origlavdata@X
      for(g in 1:lavsamplestats@ngroups) {
        x.idx <- lavsamplestats@x.idx[[g]]
        nox <- (1:nrow(Mu.hat[[g]]))[-x.idx]
        if(!is.null(x.idx) && length(x.idx) > 0L){
          if(conditional) stop("blavaan ERROR: conditional=TRUE does not work for models with fixed.x")
          
          ## for fixed.x, generate the other ovs
          ## conditional on the x values.
          ## can use approach similar to
          ## lav_mvnorm_missing_impute_pattern (lav_mvnorm_missing.R)
          if(!mis){
            tm1 <- Sigma.hat[[g]][nox,x.idx] %*% solve(Sigma.hat[[g]][x.idx,x.idx])
            csig <- Sigma.hat[[g]][nox,nox] - tm1 %*% Sigma.hat[[g]][x.idx,nox]
            sigchol <- chol(csig)

            cmu <- Mu.hat[[g]][nox,] +
              tm1 %*% apply(origlavdata@X[[g]][,x.idx,drop=FALSE], 1,
                            function(x) (x - Mu.hat[[g]][x.idx,]))

            dataX[[g]][,nox] <- t(apply(cmu, 2, function(x) mnormt::rmnorm(n=1,
                                                                           sqrt=sigchol,
                                                                           mean=x)))
          } else {
            ## condition only on observed x values;
            ## this is only needed for missing = "ml.x",
            ## which for awhile was missing = "ml" in lavaan.
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
          } # mis
        } else {
          cmu <- Mu.hat[[g]]
          csig <- Sigma.hat[[g]]

          if(!conditional) {
            dataX[[g]] <- as.matrix(mnormt::rmnorm(n = nrow(dataX[[g]]),
                                                   varcov = csig, mean = cmu))
          } else {
            sigchol <- chol(csig)
            dataX[[g]] <- t(sapply(1:NROW(cmu), function(i) mnormt::rmnorm(n=1,
                                                                           sqrt=sigchol,
                                                                           mean = cmu[i,])))
          }
        }

        ## convert to ordinal if requested
        if(ordresp){
          for(oj in ordidx){
            tmpth <- implied$mean[[g]][oj] + implied$th[[g]][thidx[[g]] == oj] * sqrt(implied$cov[[g]][oj, oj])
            tmpth <- c(-Inf, tmpth, Inf)
            dataX[[g]][, oj] <- cut(dataX[[g]][, oj], breaks = tmpth, labels = FALSE)
          }
        }
        
        dataX[[g]][is.na(origlavdata@X[[g]])] <- NA

        ## get rid of completely missing
        if(length(origlavdata@Mp[[g]]$empty.idx) > 0){
          dataX[[g]] <- dataX[[g]][-origlavdata@Mp[[g]]$empty.idx,,drop=FALSE]
        }
      } # g

      postdat[[i]] <- dataX
    }
    list(postdat = postdat, lavmod = lavmod)})

  lapcomm <- "lapply"
  if("parallel" %in% names(ddd)){
    if(ddd$parallel){
      lapcomm <- "future_lapply"
      loop.args <- c(loop.args, future.seed = TRUE)
    }
  }
  res <- do.call(lapcomm, loop.args)

  lavmod <- sapply(res, function(x) x$lavmod)
  res <- lapply(res, function(x) x$postdat)

  ## undo list over chains
  res <- do.call("c", res)

  ## remove extras when nrep/nchains is non-integer
  if(length(object) > 0L & nrep < length(res)){
    res <- res[1:nrep]
  }

  ## return lavmodel when called by postpred
  if(length(object) == 0L){
    res <- c(res, list(lavmod = lavmod))
  }

  ## list structure: generated data/groups within generated data
  res
}
