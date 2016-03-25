lav2jags <- function(model, lavdata = NULL, ov.cp = "srs", lv.cp = "srs", lv.x.wish = FALSE, dp = dpriors(), n.chains = 1, jagextra = "", inits = "prior", pta = NULL) {
  ## lots of code is taken from lav_export_bugs.R

  if(class(model)[1]=="lavaan"){
    partable <- parTable(model)
  } else {
    partable <- model
    # we assume it is a data.frame further on!
    if(!is.data.frame(partable)) {
        partable <- as.data.frame(partable, stringsAsFactors = FALSE)
    }
  }    

  ## get names of ovs before we add phantom variables
  old.pta <- lav_partable_attributes(partable = partable, pta = pta)
  old.vnames <- old.pta$vnames
  ngroups <- old.pta$ngroups
  nparam <- sum(partable$free > 0)
  orig.ov.names <- old.vnames$ov[[1]]; nov <- length(orig.ov.names)
  orig.lv.names <- old.vnames$lv[[1]]; orig.lv.names.x <- old.vnames$lv.x[[1]]
  ## so ordering stays consistent:
  orig.lv.names <- c(orig.lv.names[orig.lv.names %in% orig.lv.names.x],
                     orig.lv.names[!(orig.lv.names %in% orig.lv.names.x)])
  orig.ov.names.x <- old.vnames$ov.x[[1]]
  nlvx <- length(orig.lv.names.x)
  
  ## if lv.x.wish and default prior, change df parameter for this model
  if(lv.x.wish & nlvx > 1 & dp[["ibpsi"]] == dpriors()[["ibpsi"]]){
    dp[["ibpsi"]] <- paste("dwish(iden,", length(orig.lv.names.x) + 1, ")", sep="")
  }
  
  ## add prior column if it doesn't exist
  if(is.na(match("prior", names(partable)))) partable$prior <- rep("", length(partable$id))

  ## decide whether we need px on ovs by searching
  ## for covariances among mvs:
  mvcovs <- length(which(partable$lhs != partable$rhs &
                         partable$op == "~~" &
                         (partable$lhs %in% orig.ov.names |
                          partable$rhs %in% orig.ov.names)))
  tvname <- ifelse(mvcovs > 0, "invthetstar", "invtheta")

  ## add necessary phantom lvs/mvs to model:
  partable <- set_phantoms(partable, orig.ov.names, orig.lv.names, orig.ov.names.x, orig.lv.names.x, ov.cp, lv.cp, lv.x.wish, ngroups)
  ## ensure group parameters are in order, for parameter indexing:
  partable <- partable[order(partable$group),]
  ## get parameter table attributes 
  pta <- lav_partable_attributes(partable = partable, pta = pta)
  vnames <- pta$vnames; nvar <- pta$nvar; nfac <- pta$nfac
  ov.names.nox <- vnames$ov.nox[[1]]; nov.nox <- length(ov.names.nox)
  ov.names.x <- vnames$ov.x[[1]]; nov.x <- length(ov.names.x)
  ov.ord <- vnames$ov.ord[[1]]
  if(length(ov.ord) > 0){
    ## figure out how many categories are in the ordered variables
    ## TODO seems like this is already hidden somewhere in lavaan...
    ncats <- sapply(ov.ord, function(x) length(grep(x, vnames$th[[1]])) + 1)
  }
    
  lv.nox <- vnames$lv.nox[[1]]
  lv.names <- vnames$lv[[1]]
  ## ensure that lv.x names always come first (so we can possibly use dmnorm)
  lv.names <- c(lv.names[lv.names %in% orig.lv.names.x],
                lv.names[!(lv.names %in% orig.lv.names.x)])

  ## check that variables are the same in all groups:
  for(g in 1:ngroups){
    if(!identical(orig.ov.names, old.vnames$ov[[g]])){
      stop("blavaan ERROR: observed variables are not the same in each group.")
    }
    if(!identical(lv.names, vnames$lv[[g]])){
      if(all(lv.names %in% vnames$lv[[g]]) & length(lv.names) == length(vnames$lv[[g]])){
        next
      } else {
        stop("blavaan ERROR: latent variables are not the same in each group.")
      }
    }
  }

  ## add phantom lvs if not already there
  phnames <- unique(partable$lhs[grep(".phant", partable$lhs)])
  lv.names <- c(lv.names, phnames[!(phnames %in% lv.names)])

  ## tabs
  t1 <- paste(rep(" ", 2L), collapse="")
  t2 <- paste(rep(" ", 4L), collapse="")
  t3 <- paste(rep(" ", 6L), collapse="")
  t4 <- paste(rep(" ", 8L), collapse="")
  
  ## TXT header
  TXT <- paste("model {\n", sep="")

  TXT <- paste(TXT, t1,
               "for(i in 1:N) {\n", sep="")

  ## Second object for priors/constraints
  TXT2 <- "\n"
  ## Matrix that keeps track of parameter ordering, priors,
  ## and starting values
  coefvec <- matrix(NA, nparam, 3)
  ## Combine these to be passed to set_priors()
  priorres <- list(TXT2=TXT2, coefvec=coefvec)

  ## Third object for variance/covariance parameters
  ## (taking into account phantoms)
  TXT3 <- paste("\n", t1, "# variances & covariances\n", sep="")

  ## Decide whether we need to model exogenous x's
  if(length(ov.names.x) > 0){ # & !is.na(ov.names.x)){
    ## FIXME: this NA catch is related to filling in the exo column
    ##        in set_phantoms()
    if(!any(is.na(ov.names.x))){
      exotab <- partable[which(partable$lhs %in% old.vnames$ov.x[[1]]),]
      if(all(exotab$free==0)){
        nmvs <- nov.nox
        ov.names <- ov.names.nox
      } else {
        nmvs <- nov
        ov.names <- orig.ov.names
      }
    } else {
      ## only get here due to NA catch
      nmvs <- nov.nox
      ov.names <- orig.ov.names[orig.ov.names %in% ov.names.nox]
    }
  } else {
    nmvs <- nov.nox
    ov.names <- orig.ov.names[orig.ov.names %in% ov.names.nox]
    #ov.names <- ov.names.nox
  }

  ## Keep track of jags matrix dimensions
  ## (for initial values)
  matdims <- matrix(NA, 8, 2)
  rownames(matdims) <- c("invthet", "invpsi", "rstar",
                         "nu", "alpha", "lambda", "beta", "ibpsi")
  ## NB: don't add priors on precisions here because they mess up initial values,
  ##     and they wouldn't be included in blocks anyway
  mvv <- (partable$lhs == partable$rhs &
          partable$op == "~~" &
          partable$lhs %in% ov.names)
  #partable$prior[mvv & partable$prior == ""] <- dp[["itheta"]]
  mvvdim <- sum(mvv)/ngroups
  matdims[1,] <- c(mvvdim, ngroups)
  lvv <- (partable$lhs == partable$rhs &
          partable$op == "~~" &
          partable$lhs %in% lv.names)
  #partable$prior[lvv & partable$prior == ""] <- dp[["ipsi"]]
  lvvdim <- sum(lvv)/ngroups
  matdims[2,] <- c(lvvdim, ngroups)
  ## only find covariances under srs; fa parameterization
  ## is covered by the phantom variables and wishart is covered separately:
  if(!lv.x.wish){
    covs <- (partable$lhs != partable$rhs & partable$op == "~~")
  } else {
    covs <- (partable$lhs != partable$rhs &
             partable$op == "~~" &
             !(partable$lhs %in% orig.lv.names.x) &
             !(partable$rhs %in% orig.lv.names.x))
  }
  cov.eq <- which(covs & partable$free == 0) #partable$op == "==" & partable$rhs %in% partable$plabel[covs])
  ##partable$prior[covs & partable$prior==""] <- dp[["rho"]]
  covdim <- sum(covs)/ngroups
  #covdim <- (sum(covs) - length(cov.eq))/ngroups
  matdims[3,] <- c(covdim, ngroups)

  eqlabs <- partable$rhs[partable$op == "=="]

  ## TODO add thresholds here
  ##      we only need scale factors ~*~ in delta parameterization...
  ovi <- (partable$op == "~1" & partable$lhs %in% ov.names)
  ovifree <- (ovi & (partable$label == "" | !duplicated(partable$label)) &
              partable$free > 0 & !(partable$label %in% eqlabs))
  partable$prior[ovifree & partable$prior==""] <- dp[["nu"]]
  lvi <- (partable$op == "~1" & partable$lhs %in% lv.names)
  lvifree <- (lvi & (partable$label == "" | !duplicated(partable$label)) &
              partable$free > 0 & !(partable$label %in% eqlabs))
  partable$prior[lvifree & partable$prior==""] <- dp[["alpha"]]
  load <- (partable$op == "=~")
  loadfree <- (load & (partable$label == "" | !duplicated(partable$label)) &
               partable$free > 0 & !(partable$label %in% eqlabs))
  partable$prior[loadfree & partable$prior==""] <- dp[["lambda"]]
  reg <- (partable$op == "~" &
          partable$lhs %in% c(orig.ov.names, lv.nox) &
          (partable$rhs %in% lv.names |
           partable$rhs %in% orig.ov.names))
  regfree <- (reg & (partable$label == "" | !duplicated(partable$label)) &
              partable$free > 0 & !(partable$label %in% eqlabs))
  partable$prior[regfree & partable$prior==""] <- dp[["beta"]]

  ## block mvn priors for efficiency
  #partable <- set_blocks(partable)
  partable$blk <- rep(NA, length(partable$lhs))

  ## TODO We now attach equality constraints to these tables, could
  ##      also deal with inequality constraints.
  ## Smaller partables for different parameter types +
  ## dimensions of parameter matrices (for initial values)
  ovintercepts <- partable[ovi,]
  matdims[4,] <- c(nrow(ovintercepts)/ngroups, ngroups)
  ovintercepts <- rbind(ovintercepts, partable[which(partable$op == "=="),])
  lvintercepts <- partable[lvi,]
  matdims[5,] <- c(nrow(lvintercepts)/ngroups, ngroups)
  lvintercepts <- rbind(lvintercepts, partable[which(partable$op == "=="),])
  loadings <- partable[load,]
  matdims[6,] <- c(nrow(loadings)/ngroups, ngroups)
  loadings <- rbind(loadings, partable[which(partable$op == "=="),])
  regressions <- partable[reg,]
  matdims[7,] <- c(nrow(regressions)/ngroups, ngroups)
  regressions <- rbind(regressions, partable[which(partable$op == "=="),])
  if(lv.x.wish & nlvx > 1){
    matdims[8,] <- c(nlvx, ngroups)
  }
  
  ## Define univariate distributions of each observed variable
  ## Loop if everything is continuous
  if(length(ov.ord) == 0){
    TXT <- paste(TXT, t2,
                 "for(j in 1:", nmvs, ") {\n", sep="")
    TXT <- paste(TXT, t3,
                 "y[i,j] ~ dnorm(mu[i,j],", tvname, "[j,g[i]])\n", sep="")
  
    TXT <- paste(TXT, t2, "}\n", sep="")
  } else {
    for(j in 1:nmvs){
      if(ov.names[j] %in% ov.ord){
        ord.num <- match(ov.names[j], ov.ord)
        TXT <- paste(TXT, t2, "ones[i,", ord.num, "] ~ dbern(probs[i,", ord.num,
                     ",y[i,", j, "]])\n", sep="")

        ## category probs from pnorm
        TXT <- paste(TXT, t2, "probs[i,", ord.num, ",1] <- pnorm(tau[", ord.num,
                     ",1,g[i]], mu[i,", j, "], invthetstar[", j, ",g[i]])\n", sep="")
        if(ncats[ord.num] > 2){
          for(k in 2:(ncats[ord.num] - 1)){
            TXT <- paste(TXT, t2, "probs[i,", ord.num, ",", k, "] <- pnorm(tau[",
                         ord.num, ",", k, ",g[i]], mu[i,", j, "], invthetstar[",
                         j, ",g[i]]) - sum(probs[i,", ord.num, ",1:", (k-1), "])\n",
                         sep="")
          }
        }
        TXT <- paste(TXT, t2, "probs[i,", ord.num, ",", ncats[ord.num],
                     "] <- 1 - sum(probs[i,", ord.num, ",1:", ncats[ord.num]-1, "])\n",
                     sep="")
      } else {
        TXT <- paste(TXT, t2, "y[i,j] ~ dnorm(mu[i,j],", tvname, "[j,g[i]])\n", sep="")
      }
    }
  }

  ## Define mean of each observed variable
  ## This assumes that the data matrix passed to jags
  ## is ordered in the same way as ov.names.nox.
  ## data would be cbind(ov.names.nox, ov.names.x)
  for(i in 1:nmvs) {
    ov.idx <- i
    TXT <- paste(TXT, "\n", t2,
                 "mu[i,", ov.idx, "] <- ", sep="")

    ## find rhs for this observed variable
    ## 1. intercept?
    
    ## Always include intercept parameters, fix to zero
    ## if they are not desired
    TXT <- paste(TXT, "nu[", ov.idx, ",g[i]]", sep="")
    int.idx <- which(partable$op == "~1" &
                     partable$lhs == ov.names[i] &
                     partable$group == 1)

    ## Now deal with intercept constraints/priors:
    if(length(int.idx) == 0L) {
      for(j in 1:ngroups){
        priorres$TXT2 <- paste(priorres$TXT2, t1, "nu[", ov.idx, ",", j, "] <- 0\n", sep="")
      }
    } else {
      priorres <- set_priors(priorres, ovintercepts, ov.idx, ov.names, ngroups, "int", dp, is.na(partable$blk[int.idx]))
    }

    ## 2. factor loading? 
    lam.idx <- which(loadings$op == "=~" &
                     loadings$rhs == ov.names[i] &
                     loadings$group == 1)
    if(length(lam.idx) > 0){
      for(j in 1:length(lam.idx)) {
        TXT <- paste(TXT, " + ",
                     "lambda[", lam.idx[j], ",g[i]]*eta[i,", 
                     match(loadings$lhs[lam.idx[j]], lv.names)
                     , "]", sep="")

        ## Now assign priors/constraints
        priorres <- set_priors(priorres, loadings, i, ov.names, ngroups, "loadings", dp, is.na(loadings$blk[lam.idx[j]]), j=j)
      } # end j loop
    }

    ## 3. regression?
    r.idx <- which(regressions$lhs == ov.names[i] &
                   regressions$group == 1)
    for(j in r.idx) {
      ## what is the rhs?
      rhs <- regressions$rhs[j]
      if(rhs %in% lv.names) {
        RHS <- paste("eta[i,",
                     match(rhs, lv.names), "]", sep="")
      } else if(rhs %in% orig.ov.names) {
        RHS <- paste("y[i,",
                     match(rhs, orig.ov.names), "]", sep="")
      }
      
      ## deal with fixed later
      TXT <- paste(TXT, " + ",
                   "beta[", j, ",g[i]]*", RHS, sep="")
      ## 25Sept15 was orig.ov.names; changed to ov.names
      priorres <- set_priors(priorres, regressions, i, ov.names, ngroups, "regressions", dp, is.na(regressions$blk[j]), j=j)
    }

    ## 4. residual variance (with phantoms)
    p.idx <- which(loadings$rhs == ov.names[i] &
                   loadings$op == "=~" &
                   grepl(".phant", loadings$lhs) &
                   loadings$group == 1)
    for(k in 1:ngroups){
      if(ov.cp == "srs" & mvcovs > 0){
        TXT3 <- paste(TXT3, t1, "invthetstar[", i, ",", k,
                      "] <- 1/(theta[", i, ",", k, "]", sep="")
        for(j in p.idx){
          var.idx <- match(loadings$lhs[j], lv.names)
          TXT3 <- paste(TXT3,
                        " - (lambda[", j, ",", k, "]^2/invpsi[",
                        var.idx, ",", k, "])", sep="")
        }
        TXT3 <- paste(TXT3, ")\n", sep="")
      } else if (mvcovs > 0){
        TXT3 <- paste(TXT3, t1, "invtheta[", i, ",", k,
                      "] <- 1/(1/invthetstar[", i, ",", k,
                      "]", sep="")
        for(j in p.idx){
          var.idx <- match(loadings$lhs[j], lv.names)
          TXT3 <- paste(TXT3,
                        " + (lambda[", j, ",", k, "]^2/invpsi[",
                        var.idx, ",", k, "])", sep="")
        }
        TXT3 <- paste(TXT3, ")\n", sep="")
      }
    }

    ## 5. TODO thresholds
    if(length(ov.ord) > 0){

    }
  }
  TXT3 <- paste(TXT3, "\n", sep="")

  ## lvs
  if(length(lv.names) > 0L) {
    TXT <- paste(TXT, "\n\n", t2,
                 "# lvs", sep="")
    TXT2 <- paste(TXT2, "\n", sep="")

    ## for skipping over interactions in mu.eta and invpsi:
    lv.ind <- rbind(c(0,0))

    lvstart <- 1
    if(lv.x.wish & nlvx > 1){
      lvstart <- nlvx + 1
      lv.ind <- rbind(lv.ind, cbind(1:nlvx, 1:nlvx))

      TXT <- paste(TXT, "\n", t2,
                   "eta[i,1:", nlvx, "] ~ dmnorm(mu.eta[i,1:",
                   nlvx, "], ibpsi[1:", nlvx,",1:", nlvx, ",g[i]])", sep="")
    }
    
    nlv <- length(lv.names)

    eta.eq <- rep(NA, nlv)

    if(nlv >= lvstart){
      for(j in lvstart:nlv) {
        psi.free.idx <- which(partable$group == 1 &
                              partable$op == "~~" &
                              partable$lhs == partable$rhs &
                              partable$lhs == lv.names[j])

        if(length(psi.free.idx) != 1L) {
          stop("lavaan ERROR: parameter for residual variance ",
               lv.names[j], " not found")
        }

        ## check for latent interaction
        tmp.eq <- which(partable$op == "==" &
                        partable$lhs == lv.names[j] &
                        grepl(":", partable$rhs))
        if(length(tmp.eq) == 0L){
          lv.ind <- rbind(lv.ind, c(j, lv.ind[nrow(lv.ind),2] + 1))
          mu.ind <- lv.ind[nrow(lv.ind),2]
          ## TODO see whether we need invpsistar?
          
          lv.var <- which(partable$lhs == lv.names[mu.ind] &
                          partable$rhs == lv.names[mu.ind] &
                          partable$op == "~~")
          ## if lv variance is 0, don't assign a distribution!
          if(partable$free[lv.var] == 0 & partable$ustart[lv.var] == 0){
            TXT <- paste(TXT, "\n", t2,
                         "eta[i,", j, "] <- mu.eta[i,", mu.ind, "]", sep="")
            ## now change ustart to 1000 so no divide by 0 in jags
            partable$ustart[lv.var] <- 1000
          } else {
            TXT <- paste(TXT, "\n", t2,
                         ## TODO check for alternative distribution?
                         "eta[i,", j, "] ~ dnorm(mu.eta[i,", 
                         mu.ind, "], invpsistar[", mu.ind, ",g[i]])", sep="")
          }
        } else {
          ## latent interaction:
          eta.eq[j] <- tmp.eq
          lv.terms <- strsplit(partable$rhs[tmp.eq], ":")[[1]]
          lvs <- match(lv.terms, lv.names)
          lvs <- lvs[!is.na(lvs)]
          ## in case there is interaction with exogenous mv
          mvs <- match(lv.terms, orig.ov.names)
          mvs <- mvs[!is.na(mvs)]
          TXT <- paste(TXT, "\n", t2,
                       "eta[i,", j, "] <- eta[i,", lvs[1], "]*",
                       sep="")
          if(length(lvs) == 2){
            TXT <- paste(TXT, "eta[i,", lvs[2], "]", sep="")
          } else {
            if(length(mvs) == 0) stop("Problem with lv interaction")
            TXT <- paste(TXT, "y[i,", mvs, "]", sep="")
          }
        }
      } # j
    } # if

    ## After lv distributions are defined, now define means/regressions:
    for(j in 1:nlv) {
      ## this is for latent interactions
      if(!is.na(eta.eq[j])) next
      mu.ind <- lv.ind[which(lv.ind[,1]==j),2]
      TXT <- paste(TXT, "\n", t2,
                   ## TODO check for alternative distribution
                   ## in parameter table.
                   "mu.eta[i,", mu.ind, "] <- ", sep="")

      ## lhs elements regression
      ## 1. intercept? (even exogenous can have an intercept)
      int.idx <- which(lvintercepts$group == 1 &
                       lvintercepts$op == "~1" &
                       lvintercepts$lhs == lv.names[j])
      if(length(int.idx) == 1L) {
        ## fixed or free?
        TXT <- paste(TXT, "alpha[", mu.ind, ",", "g[i]]", sep="")

        priorres <- set_priors(priorres, lvintercepts, i, lv.names, ngroups, "lv.nox.int", dp, is.na(lvintercepts$blk[int.idx]), j=mu.ind)
      } else { # no intercept, say '0', so we always have rhs
        TXT <- paste(TXT, "0", sep="")
      }

      ## FIXME!! Possibility of loadings here
      if(lv.names[j] %in% lv.nox){
        ## 2. loadings?
        lam.idx <- which(loadings$op == "=~" &
                         loadings$rhs == lv.names[j] &
                         loadings$group == 1)
        if(length(lam.idx) > 0){
          for(k in 1:length(lam.idx)){
            TXT <- paste(TXT, " + lambda[", lam.idx[k],
                         ",g[i]]*eta[i,",
                         match(loadings$lhs[lam.idx[k]], lv.names),
                         "]", sep="")

            ## priors/constraints
            priorres <- set_priors(priorres, loadings, j, lv.names,
                                   ngroups, "loadings", dp,
                                   is.na(loadings$blk[lam.idx[k]]),
                                   j=k)
          } # end k loop
        }                         
        
        ## 3. regressions?
        rhs.idx <- which(regressions$lhs == lv.names[j] &
                         regressions$group == 1)
        np <- length(rhs.idx)
        if(np > 0){ # there could be none if we have higher-order factors
          for(p in 1:np){
            TXT <- paste(TXT, " + ",
                         "beta[", rhs.idx[p], ",g[i]]", sep="")

            ## Is the rhs an lv or ov?
            lvmatch <- match(regressions$rhs[rhs.idx[p]], lv.names)
            if(is.na(lvmatch)){
              TXT <- paste(TXT, "*y[i,", match(regressions$rhs[rhs.idx[p]], orig.ov.names), "]", sep="")
            } else {
              TXT <- paste(TXT, "*eta[i,", lvmatch, "]", sep="")
            }

            ## Now assign priors/constraints
            priorres <- set_priors(priorres, regressions, i, lv.names, ngroups, "lv.nox.reg",
                                   dp, is.na(regressions$blk[rhs.idx[p]]), j=lv.names[j], p=rhs.idx[p])
          } # end p loop
        } # end if np
      } # end if

      ## 4. lv variances (with phantoms)
      p.idx <- which(regressions$lhs == lv.names[j] &
                     regressions$op == "~" &
                     grepl(".phant", regressions$rhs) &
                     regressions$group == 1)
      ## TODO decide whether we need invpsistar?
      for(k in 1:ngroups){
        if(lv.cp == "srs"){
          TXT3 <- paste(TXT3, t1, "invpsistar[", mu.ind, ",", k,
                        "] <- 1/(psi[", mu.ind, ",", k, "]", sep="")
          for(p in p.idx){
            tmp.idx <- match(regressions$rhs[p], lv.names)
            var.idx <- lv.ind[which(lv.ind[,1]==tmp.idx),2]

            TXT3 <- paste(TXT3,
                          " - (beta[", p, ",", k, "]^2/invpsi[",
                          var.idx, ",", k, "])", sep="")
          }
        } else {
          TXT3 <- paste(TXT3, t1, "invpsi[", mu.ind, ",", k,
                        "] <- 1/(1/invpsistar[", mu.ind, ",", k,
                        "]", sep="")
          for(p in p.idx){
            tmp.idx <- match(regressions$rhs[p], lv.names)
            var.idx <- lv.ind[which(lv.ind[,1]==tmp.idx),2]

            TXT3 <- paste(TXT3,
                          " + (beta[", p, ",", k, "]^2/invpsi[",
                          var.idx, ",", k, "])", sep="")
          }
        }
        TXT3 <- paste(TXT3, ")\n", sep="")
      } # k
    } # end j loop
    if(any(lv.names %in% lv.nox)) TXT2 <- paste(TXT2, "\n", sep="")
  } # end if length(lv.names)
  
  ## end of main model specification (still need priors + equality constraints)
  TXT <- paste(TXT, "\n", t1,
               "}", sep="")

  ## now get priors for residual variance parameters of non-exo mvs and lvs
  priorres <- set_priors(priorres, partable, 1, c(ov.names, lv.names), ngroups,
                         type="vars", dp=dp, blk=TRUE, nov=nmvs, lv.names.x=orig.lv.names.x, ov.cp=ov.cp, lv.cp=lv.cp, lv.x.wish=lv.x.wish, mvcovs=mvcovs)
  ## and for correlation parameters
  priorres <- set_priors(priorres, partable, 1, c(ov.names, lv.names), ngroups,
                         type="covs", dp=dp, blk=TRUE, nov=nmvs, lv.names.x=orig.lv.names.x, ov.cp=ov.cp, lv.cp=lv.cp, lv.x.wish=lv.x.wish)

  ## and for blocked multivariate normal parameters, now that
  ## priorres$coefvec should have all the jags labels
  if(any(!is.na(partable$blk))) priorres <- block_priors(priorres, partable)
    
  ## covariances resulting from phantoms go in TXT3
  TXT3 <- paste(TXT3, "\n", sep="")
  covtable <- partable[which(partable$op == "~~" &
                             partable$lhs %in% phnames &
                             partable$group == 1),]
  if(nrow(covtable) > 0){
    for(j in 1:nrow(covtable)){
      phname <- covtable$lhs[j]
      tmp.ind <- match(phname, lv.names)
      var.ind <- lv.ind[which(lv.ind[,1]==tmp.ind),2]
      for(k in 1:ngroups){
        TXT3 <- paste(TXT3, t1, "cov[", j, ",", k, "] <- ",
                      "psi[", var.ind, ",", k, "]", sep="")
      
        tmp.ov <- partable$rhs[partable$lhs == phname &
                               partable$op == "=~" &
                               partable$group == k]
        ## find parameters from loadings, multiply by
        ## 1/invtheta
        if(length(tmp.ov) > 0L){
          for(p in 1:length(tmp.ov)){
            lam.idx <- which(loadings$op == "=~" &
                             loadings$lhs == phname &
                             loadings$rhs == tmp.ov[p] &
                             loadings$group == 1)
            TXT3 <- paste(TXT3, "*lambda[", lam.idx, ",", k, "]", sep="")
          }
        }
      
        tmp.lv <- partable$lhs[partable$rhs == phname &
                               partable$op == "~" &
                               partable$group == k]
        ## find parameters from regressions, multiply by
        ## 1/invpsi
        if(length(tmp.lv) > 0L){
          for(p in 1:length(tmp.lv)){
            bet.idx <- which(regressions$op == "~" &
                             regressions$lhs == tmp.lv[p] &
                             regressions$rhs == phname &
                             regressions$group == 1)
            TXT3 <- paste(TXT3, "*beta[", bet.idx, ",", k, "]", sep="")
          }
        }
        TXT3 <- paste(TXT3, "\n", sep="")
      } # end k
    } # end j
  } # end if

  ## end of model
  TXT <- paste(TXT, "\n\n", t1, "# Priors/constraints", priorres$TXT2, sep="")
  TXT <- paste(TXT, TXT3, sep="")
  ## extra stuff from the user, formatted to look nice-ish
  if("syntax" %in% names(jagextra)){
    jagextra <- unlist(strsplit(jagextra$syntax, "\n"))
    jagextra <- gsub("^\\s+|\\s+$", "", jagextra)
    jagextra <- paste(t1, jagextra, sep="", collapse="\n")
    TXT <- paste(TXT, "\n", jagextra, "\n", sep="")
  }
  TXT <- paste(TXT, "\n", "} # End of model\n", sep="")

  out <- TXT
  class(out) <- c("lavaan.character", "character")

  priorres$coefvec <- data.frame(priorres$coefvec, stringsAsFactors = FALSE)
  names(priorres$coefvec) <- c("jlabel", "plabel", "prior")
  out <- list(model = out, coefvec = priorres$coefvec, inits = NA)
    
  ## Initial values
  if(inits != "jags"){
      inits <- set_inits(partable, priorres$coefvec, matdims, ov.cp, lv.cp, n.chains, inits)
      out$inits <- inits
  }
        
  ## Now add data for jags if we have it
  if(!is.null(lavdata) | class(model)=="lavaan"){
    if(class(model) == "lavaan") lavdata <- model@Data
    ntot <- sum(unlist(lavdata@norig))
    y <- matrix(NA, ntot, length(orig.ov.names))
    g <- rep(NA, ntot)
    for(k in 1:ngroups){
      y[lavdata@case.idx[[k]],] <- lavdata@X[[k]]
      g[lavdata@case.idx[[k]]] <- k
    }
    ## remove deleted rows
    nas <- which(apply(is.na(y), 1, sum) == length(orig.ov.names))
    if(length(nas) > 0){
        jagsdata <- list(y=y[-nas,], g=g[-nas], N=sum(unlist(lavdata@nobs)))
    } else {
        jagsdata <- list(y=y, g=g, N=ntot)
    }

    ## identity matrix for wishart prior
    ## TODO allow user to specify this matrix
    if(lv.x.wish & length(orig.lv.names.x) > 1){
      iden <- diag(length(orig.lv.names.x))
      jagsdata <- c(jagsdata, list(iden=iden))
    }
    
    out <- c(out, list(data=jagsdata))
  }
  
  out
}

coeffun <- function(lavpartable, rjob, fun = "mean") {
  ## Extract posterior means from coda.samples() object.
  ## jagmod is the result of lav2jags;
  ## rjob is the result of run.jags() applied to
  ## the jags model.

  ## remove any rhos or unnamed parameters
  cnames <- lavpartable$jlabel
  rhos <- grep("rho", cnames)
  nn <- c(rhos, which(cnames == ""))
  if(length(nn) > 0) {
      cnames <- cnames[-nn]
  }

  ## posterior means:
  if(fun == "mean"){
    b.est <- rjob$summary$statistics[,"Mean"]
  } else if(fun == "median"){
    b.est <- rjob$hpd[,"Median"]
  }

  cmatch <- match(cnames, names(b.est), nomatch=0)
  ptmatch <- match(names(b.est), lavpartable$jlabel, nomatch=0)
  lavpartable$est[ptmatch[ptmatch != 0]] <- b.est[ptmatch != 0]

  list(x = b.est[cmatch], lavpartable = lavpartable,
       vcorr = rjob$crosscorr[cmatch,cmatch],
       sd = rjob$summary$statistics[cmatch,2])
}
