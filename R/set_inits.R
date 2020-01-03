set_inits <- function(partable, ov.cp, lv.cp, n.chains, inits){
  ## Generate initial values for each chain
  ## TODO write start values to new columns of coefvec, so can include in partable
  initvals <- vector("list", n.chains)
  names(initvals) <- paste("c", 1:n.chains, sep="")
  pveclen <- max(partable$freeparnums[partable$mat != ""], na.rm = TRUE)

  for(i in 1:n.chains){
    initvals[[i]] <- list(parvec = rep(NA, pveclen))
  }

  ## find parameters arising from wishart
  ## TODO this currently skips over priors that have been placed on
  ## variances/sds; could instead set them with some extra handling
  wps <- grep("dwish", partable$prior)

  ## handle wishart inits separately
  if(length(wps) > 0){
    wdimen <- sum(grepl("dwish", partable$prior) &
                  partable$group == 1 &
                  partable$lhs == partable$rhs)
    ngroups <- length(wps)/(wdimen*(wdimen+1)/2)
    ## generate values
    for(i in 1:n.chains){
      ## get something close to an identity matrix, otherwise
      ## the chains can go crazy places
      wvals <- rWishart(ngroups, wdimen*500, diag(wdimen))/(wdimen*500)
      initvals[[i]] <- c(initvals[[i]], list(ibpsi = wvals))
    }
  }

  for(i in 1:nrow(partable)){
    eqcons <- which(partable$lhs == partable$label[i] &
                    partable$op %in% c("==", ":=", ">", "<"))
    if((i %in% wps) | partable$free[i] == 0 | partable$prior[i] == "") next
    ## next unless it is a simple equality constraint:
    if(length(eqcons) > 0 & !grepl('^\\.p', partable$rhs[eqcons[1]])) next
    
    tmppri <- partable$prior[i]
      
    pricom <- unlist(strsplit(tmppri, "[, ()]+"))
    
    if(inits == "prior"){
      ## Try to set sensible starting values, using some of the
      ## prior information
      if(grepl("dnorm", pricom[1])){
        pricom[3] <- "1"
        ## keep loadings/regressions on one side
        if(grepl("lambda", partable$mat[i]) | grepl("beta", partable$prior[i])){
          pricom[1] <- "dunif"
          pricom[2] <- ".75"
          pricom[3] <- "2"
        }
      }
      ## Extreme correlations lead to errors, so keep them close to 0
      if(grepl("dbeta", pricom[1])){
        pricom[2] <- "100"
        pricom[3] <- "100"
      }

      ## Switch to r instead of d for random inits
      pricom[1] <- gsub("^d", "r", pricom[1])

      ## Generate initial values
      ## FIXME do something smarter upon failure
      ivs <- try(do.call(pricom[1], list(n.chains, as.numeric(pricom[2]),
                                     as.numeric(pricom[3]))), silent = TRUE)
      if(inherits(ivs, "try-error")) ivs <- rep(partable$start[i], n.chains)
    } else {
      ivs <- rep(partable$start[i], n.chains)
    }

    ## now (try to) ensure the jittered values won't crash on us
    ## and converge
    if(grepl("\\[sd\\]", partable$prior[i]) |
       grepl("\\[var\\]", partable$prior[i])){
      powval <- ifelse(grepl("\\[sd\\]", partable$prior[i]), -.5, -1)
      ivs <- ivs^powval
      ivs[ivs <= 0] <- -ivs[ivs <= 0]
    }
    if(grepl("dbeta", partable$prior[i])){
      ivs <- rep(.5, n.chains)
    }

    ## extract matrix, dimensions
    for(j in 1:n.chains){
      initvals[[j]][["parvec"]][partable$freeparnums[i]] <- ivs[j]
    }
  }

  ## if an entire init matrix/array is NA, remove it.  This can happen when
  ## all lvs are covered by dmnorm/dwish (vs invpsi) or when phantom lvs
  ## are used in tandem with dmnorm/dwish.
  all.na <- FALSE
  if(all(is.na(as.numeric(initvals[[1]][["parvec"]])))) all.na <- TRUE

  if(all.na){
    for(j in 1:n.chains){
      initvals[[j]][["parvec"]] <- NULL
    }
  }
  initvals
}

set_inits_stan <- function(partable, nfree, n.chains, inits,
                           ntot = NULL, nlvno0 = 0){
  ## Generate initial values for each chain
  initvals <- vector("list", n.chains)
  names(initvals) <- paste("c", 1:n.chains, sep="")
  pveclen <- nfree[nfree > 0]

  initmats <- list()
  for(i in 1:length(pveclen)){
    initmats <- c(initmats, list(array(NA, dim=pveclen[i])))
    ## if(pveclen[i] == 1){
    ##   initmats <- c(initmats, list(as.vector(NA)))
    ## } else {
    ##   initmats <- c(initmats, list(rep(NA, pveclen[i])))
    ## }
  }
  names(initmats) <- paste0(names(pveclen), "free")
  if(nlvno0 > 0){
    initmats <- c(initmats, list(etafree = array(1, dim = c(ntot, nlvno0))))
  }

  for(i in 1:n.chains){
    initvals[[i]] <- initmats
  }

  partable$freeparnums[is.na(partable$freeparnums)] <- 0
  freepartable <- partable[partable$freeparnums > 0,]
  if("rhoidx" %in% names(freepartable)){
      rhorows <- which(!is.na(freepartable$rhoidx) &
                       freepartable$free > 0 &
                       freepartable$mat == "rho")
      if(length(rhorows) > 0){
          freepartable$freeparnums[rhorows] <- 1:length(rhorows)
      }
      lvrhorows <- which(!is.na(freepartable$rhoidx) &
                         freepartable$free > 0 &
                         freepartable$mat == "lvrho")
      if(length(rhorows) > 0){
          freepartable$freeparnums[lvrhorows] <- 1:length(lvrhorows)
      }
  }

  ## TODO need exported, or reverse rstan::lookup()
  ## rosetta <- rstan:::rosetta
  ## alternate way to possibly get around export
  rloc <- paste0(system.file("R", package="rstan"), "/sysdata")
  lazyLoad(rloc)
  rosetta <- rosetta

  prilist <- dist2r(freepartable$prior, target = "stan")
  for(i in 1:nrow(freepartable)){
    if(inits == "prior"){
      ## Try to set sensible starting values, using some of the
      ## prior information
      pricom <- prilist[[i]]
      if(grepl("dnorm", pricom[1])){
        pricom[3] <- "1"
        ## keep loadings/regressions on one side
        if(grepl("lambda", freepartable$mat[i]) | grepl("beta", freepartable$prior[i])){
          pricom[1] <- "dunif"
          pricom[2] <- ".75"
          pricom[3] <- "2"
        }
      }
      ## Extreme correlations lead to errors, so keep them close to 0
      if(grepl("dbeta", pricom[1])){
        pricom[2] <- "100"
        pricom[3] <- "100"
      }

      ## Switch to r instead of d for random inits
      pricom[1] <- gsub("^d", "r", pricom[1])

      ## Generate initial values
      ## FIXME do something smarter upon failure
      ivs <- try(do.call(pricom[1], list(n.chains, as.numeric(pricom[2]),
                                         as.numeric(pricom[3]))), silent = TRUE)

      if(inherits(ivs, "try-error")){
        ivs <- rep(1, n.chains)
      } else if(pricom[1] == "rgamma" & !grepl("[sd]", freepartable$prior[i], fixed = TRUE) &
                !grepl("[var]", freepartable$prior[i], fixed = TRUE)){
        ## free parameter is a precision, not a variance/sd:
        ivs <- 1/ivs
      }
    } else {
      ivs <- rep(freepartable$start[i], n.chains)
    }

    ## now (try to) ensure the jittered values won't crash on us
    ## and converge
    ## if(grepl("\\[sd\\]", freepartable$prior[i]) |
    ##    grepl("\\[var\\]", freepartable$prior[i])){
    ##   powval <- ifelse(grepl("\\[sd\\]", freepartable$prior[i]), -.5, -1)
    ##   ivs <- ivs^powval
    ##   ivs[ivs <= 0] <- -ivs[ivs <= 0]
    ## }
    if(grepl("beta", freepartable$prior[i])){
      ivs <- rep(.5, n.chains)
    }

    ## extract matrix, dimensions
    for(j in 1:n.chains){
      matidx <- which(names(initvals[[j]]) == paste0(freepartable$mat[i], "free"))
      initvals[[j]][[matidx]][freepartable$freeparnums[i]] <- ivs[j]
    }
  }

  initvals
}
