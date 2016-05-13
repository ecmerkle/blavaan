set_inits <- function(partable, coefvec, matdims, ov.cp, lv.cp, n.chains, inits){
  ## Generate initial values for each chain
  ## TODO write start values to new columns of coefvec, so can include in partable
    
  initvals <- vector("list", n.chains)
  names(initvals) <- paste("c", 1:n.chains, sep="")

  ## replace variance parameter names with precisions  
  if(ov.cp == "srs"){
    coefvec[,1] <- gsub("theta", "invtheta", coefvec[,1])
    rownames(matdims)[1] <- "invtheta"
  } else {
    coefvec[,1] <- gsub("theta", "invthetstar", coefvec[,1])
    rownames(matdims)[1] <- "invthetstar"
  }
  if(lv.cp == "srs"){
    coefvec[,1] <- gsub("psi", "invpsi", coefvec[,1])
  } else {
    coefvec[,1] <- gsub("psi", "invpsistar", coefvec[,1])
    rownames(matdims)[2] <- "invpsistar"
  }

  ## invtheta <- matrix(NA, matdims[grep("invthet", rownam,1], matdims["invthet",2])
  ## invpsi <- matrix(NA, matdims["invpsi",1], matdims["invpsi",2])

  ## only include parameter matrices with some free parameters
  parnames <- rownames(matdims)
  for(i in 1:length(parnames)){
    if(parnames[i] == "ibpsi"){
      if(!is.na(matdims[i,1])){
        tmpmat <- array(NA, dim=c(matdims[parnames[i],1], matdims[parnames[i],1],
                                  matdims[parnames[i],2]))

        for(j in 1:n.chains){
          initvals[[j]] <- c(initvals[[j]], list(tmpmat))
          names(initvals[[j]])[length(initvals[[j]])] <- parnames[i]
        }
      }
    } else {
      tmpmat <- matrix(NA, matdims[parnames[i],1], matdims[parnames[i],2])
      tmpname <- parnames[i]
      if(tmpname == "rstar") tmpname <- "rho"
      if(any(grepl(tmpname, coefvec[,1]))){
        for(j in 1:n.chains){
          initvals[[j]] <- c(initvals[[j]], list(tmpmat))
          names(initvals[[j]])[length(initvals[[j]])] <- parnames[i]
        }
      }
    }
  }

  ## find redundant parameters and parameters arising from wishart
  ## TODO this currently skips over priors that have been placed on
  ## variances/sds; could instead set them with some extra handling
  rps <- which(coefvec[,3] == "" | grepl("@", coefvec[,2]) |
               grepl(" <-", coefvec[,1]) | grepl("\\[", coefvec[,3]))
  wps <- grep("dwish", coefvec[,3])

  ## handle wishart inits separately
  if("ibpsi" %in% names(initvals[[1]])){
    wdimen <- sum(grepl("dwish", coefvec[,3]) & grepl("invpsi", coefvec[,1]) & grepl(",1]", coefvec[,1]))
    ngroups <- length(wps)/(wdimen*(wdimen+1)/2)
    ## generate values
    for(i in 1:n.chains){
      ## get something close to an identity matrix, otherwise
      ## the chains can go crazy places
      wvals <- rWishart(ngroups, wdimen*500, diag(wdimen))/(wdimen*500)
      initvals[[i]][["ibpsi"]] <- wvals
    }
  }

  for(i in 1:nrow(coefvec)){
    if((i %in% rps) | (i %in% wps)) next

    tmppri <- coefvec[i,3]
      
    pricom <- unlist(strsplit(tmppri, "[, ()]+"))
    
    if(inits == "prior"){
      ## Try to set sensible starting values, using some of the
      ## prior information
      if(grepl("dnorm", pricom[1])){
        pricom[3] <- "1"
        ## keep loadings/regressions on one side
        if(grepl("lambda", coefvec[i,1]) | grepl("beta", coefvec[i,1])){
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
      if(inherits(ivs, "try-error")) ivs <- rep(NA, n.chains)
    } else {
      ptrow <- which(partable$plabel == coefvec[i,2])[1]
      ivs <- rep(partable$start[ptrow], n.chains)

      ## now (try to) ensure the jittered values won't crash on us
      if(grepl("thet", coefvec[i,1]) | grepl("psi", coefvec[i,1])){
        ivs <- 1/ivs
        ivs[ivs <= 0] <- -ivs[ivs <= 0]
      }
      if(grepl("rho", coefvec[i,1]) | grepl("rstar", coefvec[i,1])){
        ivs <- rep(.5, n.chains)
      }
    }

    ## extract matrix, dimensions
    mat <- unlist(strsplit(coefvec[i,1], "[, \\[^\\]]+", perl=T))
    if(mat[1] == "rho") mat[1] <- "rstar"

    for(j in 1:n.chains){
      initvals[[j]][[mat[1]]][as.numeric(mat[2]), as.numeric(mat[3])] <- ivs[j]
    }
  }

  ## if an entire init matrix/array is NA, remove it.  This can happen when
  ## all lvs are covered by dmnorm/dwish (vs invpsi) or when phantom lvs
  ## are used in tandem with dmnorm/dwish.
  all.na <- NULL
  n.mats <- length(initvals[[1]])
  for(i in 1:n.mats){
    if(all(is.na(as.numeric(initvals[[1]][[i]])))) all.na <- c(all.na, i)
  }
  if(length(all.na) > 0){
    for(j in 1:n.chains){
      ## backwards because the list indexes change when we go forwards:
      for(i in length(all.na):1){
        initvals[[j]][[all.na[i]]] <- NULL
      }
    }
  }

  initvals
}
