lav2stan <- function(model, lavdata = NULL, dp = NULL, n.chains = 1, mcmcextra = "", inits = "prior") {
  ## lots of code is taken from lav_export_bugs.R

  if(class(model)[1]=="lavaan"){
    partable <- parTable(model)
  } else {
    stop("blavaan ERROR: model must be class lavaan")
  }
  
  eqop <- "="
  commop <- "// "
  eolop <- ";"
  if(length(dp) == 0) dp <- dpriors(target = "stan")
  
  ## get names of ovs before we add phantom variables
  old.pta <- lav_partable_attributes(partable = partable, pta = NULL)
  old.vnames <- old.pta$vnames
  ngroups <- old.pta$ngroups
  orig.ov.names <- old.vnames$ov[[1]]; nov <- length(orig.ov.names)
  orig.lv.names <- old.vnames$lv[[1]]; orig.lv.names.x <- old.vnames$lv.x[[1]]
  ## so ordering stays consistent:
  orig.lv.names <- c(orig.lv.names[orig.lv.names %in% orig.lv.names.x],
                     orig.lv.names[!(orig.lv.names %in% orig.lv.names.x)])
  orig.ov.names.x <- old.vnames$ov.x[[1]]
  nlvx <- length(orig.lv.names.x)

  ## set up mvs with fixed 0 variances (single indicators of lvs)
  partable <- set_mv0(partable, orig.ov.names, ngroups)
  ## convert covariances to corr * sd1 * sd2
  partable <- set_stancovs(partable)

  ## ensure group parameters are in order, for parameter indexing:
  partable <- partable[order(partable$group),]
  ## get parameter table attributes 
  pta <- lav_partable_attributes(partable = partable, pta = NULL)
  vnames <- pta$vnames; nvar <- pta$nvar; nfac <- pta$nfac
  ov.names.nox <- vnames$ov.nox[[1]]; nov.nox <- length(ov.names.nox)
  ov.names.x <- vnames$ov.x[[1]]
  ## lavaan FIXME? if no x, ov.names.x is sometimes length 0,
  ## sometimes NA
  if(length(ov.names.x) > 0){
    if(all(is.na(ov.names.x))){
      nov.x <- 0
    } else {
      nov.x <- length(ov.names.x)
    }
  } else {
    nov.x <- 0
  }
  ov.ord <- vnames$ov.ord[[1]]
  lv.nox <- vnames$lv.nox[[1]]
  lv.names <- vnames$lv[[1]]
  ## ensure that lv.x names always come first (so we can possibly use dmnorm)
  #lv.names <- c(lv.names[lv.names %in% orig.lv.names.x],
  #              lv.names[!(lv.names %in% orig.lv.names.x)])
  nlv <- length(lv.names)
  
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

  ## tabs
  t1 <- paste(rep(" ", 2L), collapse="")
  t2 <- paste(rep(" ", 4L), collapse="")
  t3 <- paste(rep(" ", 6L), collapse="")
  t4 <- paste(rep(" ", 8L), collapse="")
  
  ## stan blocks
  datblk <- parblk <- TPS <- TXT <- ""
  ## hold priors to put at bottom of model block
  TXT2 <- ""
  TXT <- paste("model {\n", sep="")

  ## lvs with variances fixed to 0
  lv0 <- which(partable$op == "~~" &
               partable$lhs %in% lv.names &
               partable$rhs == partable$lhs &
               partable$group == 1 &
               partable$free == 0 &
               partable$ustart == 0)
  if(length(lv0) > 0){
    lv0.names <- partable$lhs[lv0]
    lv0.idx <- which(lv.names %in% lv0.names)
    nlvno0 <- nlv - length(lv0.idx)
  } else {
    lv0.names <- NULL
    lv0.idx <- NULL
    nlvno0 <- nlv
  }

  nmvs <- nov
  ov.names <- orig.ov.names

  eqlabs <- partable$rhs[partable$op %in% c("==", ":=")]
  eqplabs <- partable$lhs[partable$op %in% c("==", ":=")]
  eqplabs <- eqplabs[eqplabs %in% partable$label]
  eqlabs <- c(eqlabs, eqplabs)

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

  loadings <- partable[load,]
  loadings <- rbind(loadings, partable[which(partable$op %in% c("==", ":=")),])
  regressions <- partable[reg,]
  regressions <- rbind(regressions, partable[which(partable$op %in% c("==", ":=")),])

  ## number of free parameters per type, for stan parameter vectors
  ## (need separated so can use "lower" and "upper")
  parmats <- lavInspect(model)
  parmattable <- lavInspect(model, 'est')
  parconst <- attr(parmats, "header")
  gamind <- "gamma" %in% names(parmats[[1]])

  ## so it is always a list of lists
  if(ngroups == 1){
    parmats <- list(g1 = parmats)
    parmattable <- list(g1 = parmattable)
  }

  ## decide whether psi is diagonal and whether beta is
  ## lower/upper triangular, for faster matrix computations
  ## in stan
  diagpsi <- 0L
  if("psi" %in% names(parmattable[[1]])){
    tmppsi <- parmattable[[1]]$psi
    tmppsi <- tmppsi[lower.tri(tmppsi)]
    if(all(tmppsi == 0)) diagpsi <- 1L
  }
  fullbeta <- 1L

  if("beta" %in% names(parmattable[[1]])){
    tmpbeta <- parmattable[[1]]$beta
    if(all(tmpbeta[lower.tri(tmpbeta)] == 0) |
       all(tmpbeta[upper.tri(tmpbeta)] == 0)) fullbeta <- 0L
  }

  nfree <- sapply(parmats, sapply, function(x){
    if(class(x)[1] == "lavaan.matrix.symmetric"){
      # off-diagonals handled via rho parameters, unless they
      # are both ov.names.x
      if(FALSE){ #rownames(x)[1] %in% c(lv.names, ov.names.x)){
        covpars <- which(partable$op == "~~" &
                         partable$lhs != partable$rhs &
                         partable$free > 0L &
                         partable$lhs %in% ov.names.x)
        length(covpars) + sum(diag(x) > 0)
      } else {
        sum(diag(x) > 0)
      }
    } else {
      sum(x > 0)
    }})
  if(length(parconst) > 0){
    nfix <- sapply(parmats, sapply, function(x){
      if(class(x)[1] == "lavaan.matrix.symmetric"){
        sum(diag(x) %in% parconst$rhs)
      } else {
        sum(x %in% parconst$rhs)
      }})
  } else {
    nfix <- 0
  }
  nfree <- apply(nfree - nfix, 1, sum)

  parblk <- paste0(parblk, "parameters{\n")
  for(i in 1:length(nfree)){
    if(nfree[i] == 0) next

    parnm <- names(nfree)[i]
    parblk <- paste0(parblk, t1, "vector")
    if(parnm %in% c("theta", "psi")){
      parblk <- paste0(parblk, "<lower=0>")        
    }
    parblk <- paste0(parblk, "[", nfree[i], "]")
    parblk <- paste0(parblk, " ", parnm, "free", eolop, "\n")
  }

  if(any(partable$mat == "rho")){
    #nrhofix <- sum(sapply(parmats, function(x){
    #  sum(x$theta[lower.tri(x$theta)] %in% parconst$rhs)
    #}))
    nrho <- sum(partable$mat == "rho" &
                partable$free > 0 &
                !is.na(partable$rhoidx), na.rm = TRUE)# - nrhofix
    nfree <- c(nfree, rho = nrho)
    parblk <- paste0(parblk, t1, "vector<lower=0,upper=1>[",
                     nrho, "] rhofree;\n")
  }
  if(any(partable$mat == "lvrho")){
    #nlrhofix <- sum(sapply(parmats, function(x){
    #  sum(x$psi[lower.tri(x$psi)] %in% parconst$rhs)
    #}))
    nlrho <- sum(partable$mat == "lvrho" &
                 partable$free > 0 &
                 !is.na(partable$rhoidx), na.rm = TRUE)# - nlrhofix
    nfree <- c(nfree, lvrho = nlrho)
    parblk <- paste0(parblk, t1, "vector<lower=0,upper=1>[",
                     nlrho, "] lvrhofree;\n")
  }

  if(nlvno0 > 0){
    parblk <- paste0(parblk, t1, "matrix[N, ", nlvno0, "] etavec",
                     eolop, "\n")
  }
  parblk <- paste0(parblk, "}\n\n")                     

  psi.ov <- which(partable$lhs %in% ov.names &
                  partable$op == "~~" &
                  partable$lhs == partable$rhs &
                  partable$group == 1 &
                  partable$mat == "psi")
  n.psi.ov <- length(psi.ov)
  ny <- nov - n.psi.ov
  if(n.psi.ov > 0){
    psi.ov.names <- partable$lhs[psi.ov]
    thet.ov.names <- ov.names[!(ov.names %in% psi.ov.names)]
  } else {
    psi.ov.names <- ""
    thet.ov.names <- ov.names
  }
  yind <- which(ov.names %in% thet.ov.names)
  xind <- which(ov.names %in% psi.ov.names)
  ov.dummy.idx <- c(model@Model@ov.y.dummy.ov.idx[[1]],
                    model@Model@ov.x.dummy.ov.idx[[1]])
  lv.dummy.idx <- c(model@Model@ov.y.dummy.lv.idx[[1]],
                    model@Model@ov.x.dummy.lv.idx[[1]])
  dumov <- 0L
  if(length(ov.dummy.idx) > 0) dumov <- 1L

  ## FIXME? see .internal_get_ALPHA from lav_representation_lisrel.R
  ## for alternative (better) way to handle this than eqs.x  
  if(nov.x > 0 | length(vnames$eqs.x[[1]]) > 0){
    xnames <- c(ov.names.x, vnames$eqs.x[[1]])
    exoind <- which(ov.names[xind] %in% xnames)
    regind <- which(!(ov.names[xind] %in% xnames))
    etaind <- 1:nlv
    if(nlv > 0 & length(lv0.idx) < nlv){
      if(length(lv0.idx) > 0){
        nlvno0 <- nlv - length(lv0.idx)
        regind <- c((1:nlv)[-lv0.idx], (nlvno0+regind))
        exoind <- nlvno0 + exoind
        etaind <- etaind[-lv0.idx]
      } else {
        nlvno0 <- nlv
        regind <- c(1:nlv, (nlv+regind))
        exoind <- nlv + exoind
      }
    }
    lvindall <- c(regind, exoind)
  } else {
    regind <- xind
    exoind <- rep(0,length(xind))
    lvindall <- regind
    etaind <- exoind
    if(nlv > 0){
      if(length(lv0.idx) < nlv & length(lv0.idx) > 0){
        nlvno0 <- nlv - length(lv0.idx)
        etaind <- (1:nlv)[-lv0.idx]
      } else {
        nlvno0 <- nlv
        etaind <- 1:nlv
      }
    }
  }
  
  ## missingness of ovs split by whether or not they appear
  ## in psi
  missflag <- FALSE
  miss.psi <- FALSE
  if(n.psi.ov > 0){
    for(k in 1:ngroups){
      miss.psi <- (miss.psi | any(is.na(lavdata@X[[k]][,xind])))
    }
  }
  if(length(yind) > 0){
    for(k in 1:ngroups){
      missflag <- (missflag | any(is.na(lavdata@X[[k]][,yind])))
    }
  }
  
  TXT <- paste0(TXT, t1, "for(i in 1:N) {\n")

  if(ny > 0){
    if(missflag){
      TXT <- paste0(TXT, t2,
                    "segment(y[i], 1, nseen[i]) ~ ",
                    "multi_normal_cholesky(",
                    "to_vector(mu[i])[obsvar[i,1:nseen[i]]],",
                    "thetld[g[i],obsvar[i,1:nseen[i]],",
                    "obsvar[i,1:nseen[i]]]);\n")
    } else {
      TXT <- paste0(TXT, t2,
                    "y[i] ~ multi_normal_cholesky(",
                    "to_vector(mu[i,1:", (nov - n.psi.ov),
                    "]), thetld[g[i]]);\n")
    }
  }


  TXT <- paste0(TXT, t1, "}\n\n")

  if((nlv + n.psi.ov) > 0){
    TXT <- paste0(TXT, t1, "eta ~ ")

    if(miss.psi){
      TXT <- paste0(TXT, "sem_lv_missing_lpdf(")
    } else {
      TXT <- paste0(TXT, "sem_lv_lpdf(")
    }

    TXT <- paste0(TXT, "alpha, beta, psi, ")
    TXT <- paste0(TXT, ifelse(gamind, "gamma", "beta"), ", ")
    TXT <- paste0(TXT, as.numeric(gamind), ", meanx, ")
    TXT <- paste0(TXT, "g, ", (nlv + n.psi.ov), ", N, ",
                  ngroups, ", ", diagpsi, ", ", fullbeta, ", ", nlv,
                  ", etaind, ", nlvno0)
    if(miss.psi){
      TXT <- paste0(TXT, ", nseenx, obsvarx")
    }
    TXT <- paste0(TXT, ");\n")
  }

  ## for missing=="fi", to model variables on rhs of regression
  ovreg <- unique(regressions$rhs[regressions$rhs %in% ov.names])
  ovcol <- which(ov.names %in% ovreg)

  if(nlvno0 < nlv){
    TPS <- paste0(TPS, t1, "for(i in 1:N) {\n")
    TPS <- paste0(TPS, t2, "eta[i,etaind] = etavec[i];\n", t1, "}\n")
    
    TPS <- paste0(TPS, t1, "mueta = sem_mean_eta(alpha, eta, ",
                  "beta, ", ifelse(gamind, "gamma", "beta"),
                  ", g, ", (nlv + n.psi.ov),
                  ", N, ", ngroups, ", lvind);\n")
  }
  
  ## Define mean of each observed variable
  ## This assumes that the data matrix passed to jags
  ## is ordered in the same way as ov.names.nox.
  ## data would be cbind(ov.names.nox, ov.names.x)
  TPS <- paste(TPS, t1, commop, "mu definitions\n", t1,
               "for(i in 1:N) {\n", sep="")

  if(nlvno0 < nlv){
    #TPS <- paste0(TPS, t2, "eta[i,etaind] = etavec[i];\n")
    TPS <- paste0(TPS, t2, "eta[i,eta0ind] = mueta[i,eta0ind]';\n")
  } else if(nlv > 0){
    TPS <- paste0(TPS, t2, "eta[i,1:", nlv, "] = etavec[i];\n")
  }

  if(n.psi.ov > 0){
    TPS <- paste0(TPS, t2, "eta[i,", (nlv+1), ":", (nlv + n.psi.ov),
                  "] = x[i]';\n")
  }
  
  if(ny > 0) {
    for(i in 1:ny) {
      ov.idx <- i
      if(i > 1) TPS <- paste(TPS, eolop, sep="")
      TPS <- paste(TPS, "\n", t2, "mu[i,", ov.idx, "] ", eqop, " ",
                   sep="")

      ## find rhs for this observed variable
      ## 1. intercept?
    
      ## Always include intercept parameters, fix to zero
      ## if they are not desired
      int.idx <- which(partable$op == "~1" &
                       partable$lhs == thet.ov.names[i] &
                       partable$group == 1)

      ## Now deal with intercept constraints/priors:
      if(length(int.idx) == 0L) {
        TPS <- paste(TPS, "nu[", ov.idx, ",1,g[i]]", sep="")
      } else {
        TPS <- paste(TPS, partable$mat[int.idx], "[", partable$row[int.idx], ",",
                     partable$col[int.idx], ",g[i]]", sep="")
      }

      ## 2. factor loading? 
      lam.idx <- which(loadings$op == "=~" &
                       loadings$rhs == thet.ov.names[i] &
                       loadings$group == 1)
      if(length(lam.idx) > 0){
        for(j in 1:length(lam.idx)) {
          TPS <- paste(TPS, " + ", loadings$mat[lam.idx[j]], "[",
                       loadings$row[lam.idx[j]], ",", loadings$col[lam.idx[j]],
                       ",g[i]]*eta[i,", match(loadings$lhs[lam.idx[j]], lv.names)
                     , "]", sep="")
        } # end j loop
      }

      ## 3. regression?
      r.idx <- which(regressions$lhs == thet.ov.names[i] &
                     regressions$group == 1)
      for(j in r.idx) {
        ## what is the rhs?
        rhs <- regressions$rhs[j]
        if(rhs %in% lv.names) {
          RHS <- paste("etamat[i,", match(rhs, lv.names), "]", sep="")
        } else if(rhs %in% thet.ov.names) {
          RHS <- paste("etamat[i,", nlv + match(rhs, thet.ov.names), "]", sep="")
        } else if(rhs %in% psi.ov.names) {
          RHS <- paste("x[i,", match(rhs, psi.ov.names), "]", sep="")
        }
      
        ## deal with fixed later
        TPS <- paste(TPS, " + ", regressions$mat[j], "[",
                     regressions$row[j], ",", regressions$col[j],
                     ",g[i]]*", RHS, sep="")
      }
    }
  }
  TPS <- paste(TPS, eolop, "\n", sep="")

  ## priors/constraints
  TXT2 <- set_stanpars(TXT2, partable, nfree, dp, orig.lv.names.x)
  partable$prior <- TXT2$partable$prior
  partable$freeparnums <- TXT2$partable$freeparnums
  TXT3 <- TXT2$TXT3
  TXT2 <- TXT2$TXT2
  ## end of main model specification

  ## extra stuff from the user, formatted to look nice-ish
  if("syntax" %in% names(mcmcextra)){
    mcmcextra <- unlist(strsplit(mcmcextra$syntax, "\n"))
    mcmcextra <- gsub("^\\s+|\\s+$", "", mcmcextra)
    mcmcextra <- paste(t1, mcmcextra, sep="", collapse="\n")
    TXT <- paste(TXT, "\n", mcmcextra, "\n", sep="")
  }
  
  out <- TXT
  out <- paste0(out, TXT3, "\n}")
  class(out) <- c("lavaan.character", "character")
  out <- list(model = out, inits = NA)

  ## Initial values
  inits <- set_inits_stan(partable, nfree, n.chains, inits)
  out$inits <- inits

  ## Now add data if we have it
  datablk <- paste0("data{\n", t1, "int N;\n", t1, "int g[N];\n",
                    t1, "int lvind[", length(lvindall), "];\n",
                    t1, "int etaind[", length(etaind), "];\n")
  if(length(lv0.idx) > 0){
    datablk <- paste0(datablk, t1, "int eta0ind[", length(lv0.idx),
                      "];\n")
  }

  ## NB: if meanx is empty, we won't use it. so just
  ## set meanx to smean for stan.
  smean <- do.call("cbind", model@SampleStats@mean)
  if(length(model@SampleStats@mean.x[[1]]) > 0){
    meanx <- do.call("cbind", model@SampleStats@mean.x)
  } else {
    meanx <- smean
  }
  datablk <- paste0(datablk, t1, "real sampmean[", nrow(smean), ",",
                    ncol(smean), "];\n", t1,
                    "real meanx[", nrow(meanx), ",", ncol(meanx),
                    "];\n")

  if(length(ov.dummy.idx) == 0){
    ov.dummy.idx <- rep(0,2)
    lv.dummy.idx <- rep(0,2)
  }
  datablk <- paste0(datablk, t1, "int dummyov[",
                    length(ov.dummy.idx), "];\n", t1,
                    "int dummylv[", length(lv.dummy.idx), "];\n")

  if(!is.null(lavdata) | class(model)[1]=="lavaan"){
    if(class(model)[1] == "lavaan") lavdata <- model@Data
    ntot <- sum(unlist(lavdata@norig))

    ## exogenous x's go in their own matrix
    y <- matrix(NA, ntot, ny) #lapply(1:tmpnmvs, function(x) rep(NA,ntot))
    if(n.psi.ov > 0) x <- matrix(NA, ntot, n.psi.ov)

    if(ny > 0){
      obsvar <- matrix(-999, ntot, ny)
      nseen <- rep(NA, ntot)
      misvar <- matrix(-999, ntot, ny)
      nmis <- rep(NA, ntot)
    }
    if(n.psi.ov > 0){
      obsvarx <- matrix(-999, ntot, n.psi.ov)
      nseenx <- rep(NA, ntot)
      misvarx <- matrix(-999, ntot, n.psi.ov)
      nmisx <- rep(NA, ntot)
    }
    if(length(exoind) > 0){
      obsexo <- matrix(-999, ntot, length(exoind))
      nseenexo <- rep(NA, ntot)
    } else {
      obsexo <- matrix(1, ntot, 1)
      nseenexo <- rep(1, ntot)
    }
    
    g <- rep(NA, ntot)

    for(k in 1:ngroups){
      if(ny > 0){
        for(j in 1:ny){
          y[lavdata@case.idx[[k]],j] <- lavdata@X[[k]][,yind[j]]
        }
      }
      if(n.psi.ov > 0){
        for(j in 1:n.psi.ov){
          x[lavdata@case.idx[[k]],j] <- lavdata@X[[k]][,xind[j]]
        }
      }
      g[lavdata@case.idx[[k]]] <- k
      
      ## missingness patterns
      npatt <- lavdata@Mp[[k]]$npatterns
      for(m in 1:npatt){
        if(ny > 0){
          tmpobs <- which(lavdata@Mp[[k]]$pat[m,])
          tmpobs <- tmpobs[tmpobs %in% yind]
          tmpmis <- which(!lavdata@Mp[[k]]$pat[m,])
          tmpmis <- tmpmis[tmpmis %in% yind]
          tmpidx <- lavdata@Mp[[k]]$case.idx[[m]]
          nseen[tmpidx] <- length(tmpobs)
          if(length(tmpobs) > 0){
            tmpobs <- matrix(tmpobs, length(tmpidx), length(tmpobs),
                             byrow=TRUE)
            obsvar[tmpidx,1:nseen[tmpidx[1]]] <- tmpobs
          }
          nmis[tmpidx] <- length(tmpmis)
          if(length(tmpmis) > 0){
            tmpmis <- matrix(tmpmis, length(tmpidx), length(tmpmis),
                             byrow=TRUE)
            misvar[tmpidx,1:nmis[tmpidx[1]]] <- tmpmis
          }
        }

        if(n.psi.ov > 0){
          ## now for x
          tmpobsx <- which(lavdata@Mp[[k]]$pat[m,])
          tmpobs <- tmpobsx[tmpobsx %in% xind]
          tmpmis <- which(!lavdata@Mp[[k]]$pat[m,])
          tmpmis <- tmpmis[tmpmis %in% xind]
          tmpidx <- lavdata@Mp[[k]]$case.idx[[m]]
          nseenx[tmpidx] <- length(tmpobs)
          tmpobsexo <- which(tmpobs %in% exoind)
          nseenexo[tmpidx] <- length(tmpobsexo)
          if(length(tmpobs) > 0){
            tmpobs <- matrix(tmpobs, length(tmpidx), length(tmpobs),
                             byrow=TRUE)
            obsvarx[tmpidx,1:nseenx[tmpidx[1]]] <- tmpobs
          }
          if(length(tmpobsexo) > 0){
            tmpobsexo <- matrix(tmpobsexo, length(tmpidx),
                                length(tmpobsexo), byrow=TRUE)
            obsexo[tmpidx,1:nseenexo[tmpidx[1]]] <- tmpobsexo
          }
          nmisx[tmpidx] <- length(tmpmis)
          if(length(tmpmis) > 0){
            tmpmis <- matrix(tmpmis, length(tmpidx), length(tmpmis),
                             byrow=TRUE)
            misvarx[tmpidx,1:nmisx[tmpidx[1]]] <- tmpmis
          }
        }
      }
    }
    if(ny > 0) colnames(y) <- ov.names[yind]
    if(n.psi.ov > 0) colnames(x) <- ov.names[xind]
    
    ## remove fully deleted rows
    yna <- rep(FALSE, ntot)
    xna <- rep(FALSE, ntot)
    if(ny > 0) yna <- apply(is.na(y), 1, sum) == ny
    if(n.psi.ov > 0) xna <- apply(is.na(x), 1, sum) == n.psi.ov
    nas <- which(yna | xna)

    if(length(nas) > 0){
      if(ny > 0) y <- y[-nas,]
      if(n.psi.ov > 0) x <- x[-nas]
      g <- g[-nas]
      if(ny > 0){
        obsvar <- obsvar[-nas,]
        misvar <- misvar[-nas,]
        nseen <- nseen[-nas]
        nmis <- nmis[-nas]
      }
      if(n.psi.ov > 0){
        obsvarx <- obsvarx[-nas,]
        misvarx <- misvarx[-nas,]
        obsexo <- obsexo[-nas,]
        nseenx <- nseenx[-nas]
        nmisx <- nmisx[-nas]
        nseenexo <- nseenexo[-nas]
      }
      ntot <- sum(unlist(lavdata@nobs))
    }

    ## move observed variables all to the left, because stan only
    ## allows segment() to be contiguous
    if(missflag){
      if(ny > 0){
        for(i in 1:nrow(y)){
          y[i,1:nseen[i]] <- y[i,obsvar[i,1:nseen[i]]]
          if(ny - nseen[i] > 0){
            y[i,(nseen[i]+1):ny] <- -999
          }
        }
      }
    }
    if(miss.psi){
      if(n.psi.ov > 0){
        for(i in 1:nrow(x)){
          x[i,1:nseenx[i]] <- x[i,obsvarx[i,1:nseenx[i]]]
          if(n.psi.ov - nseenx[i] > 0){
            x[i,(nseenx[i]+1):n.psi.ov] <- -999
          }
        }
      }
    }

    standata <- list(g=g, N=ntot, regind=array(regind),
                     exoind=array(exoind), lvind=array(lvindall),
                     etaind=array(etaind))
    if(length(lv0.idx) > 0){
      standata <- c(standata, list(eta0ind=array(lv0.idx)))
    }
    standata <- c(standata, list(dummyov=array(ov.dummy.idx),
                                 dummylv=array(lv.dummy.idx),
                                 sampmean=smean,
                                 meanx=meanx))

    if(ny > 0) standata <- c(standata, list(y=y))
    if(n.psi.ov > 0) standata <- c(standata, list(x=x))
    if(missflag){
      if(ny > 0){
        standata <- c(standata, list(obsvar=obsvar, misvar=misvar,
                                     nseen=nseen, nmis=nmis))
        standata$y[is.na(standata$y)] <- -999
      }
    }
    if(miss.psi){
      if(n.psi.ov > 0){
        standata <- c(standata, list(obsvarx=obsvarx,
                                     misvarx=misvarx,
                                     obsexo=obsexo,
                                     nseenx=nseenx, nmisx=nmisx,
                                     nseenexo=nseenexo))
        standata$x[is.na(standata$x)] <- -999
      }
    }
    ## TODO needed?
    if(any(partable$op == "|")){
      standata <- c(standata, list(ones = matrix(1, ntot, nov.nox)))
    }

    ## stan data block
    if(ny > 0){
      datablk <- paste0(datablk, t1, "vector[", ny, "] y[N];\n")
    }
    if(n.psi.ov > 0){
      datablk <- paste0(datablk, t1, "vector[", n.psi.ov,
                        "] x[N];\n")
    }
    if(missflag){
      if(ny > 0){
        datablk <- paste0(datablk, t1, "int obsvar[N,", ny,
                          "];\n", t1, "int misvar[N,", ny,
                          "];\n", t1, "int nseen[N];\n",
                          t1, "int nmis[N];\n")
      }
    }
    if(miss.psi){
      if(n.psi.ov > 0){
        datablk <- paste0(datablk, t1, "int obsvarx[N,", n.psi.ov,
                          "];\n", t1, "int misvarx[N,", n.psi.ov,
                          "];\n", t1, "int obsexo[N,",
                          ncol(obsexo), "];\n", t1,
                          "int nseenx[N];\n", t1,
                          "int nmisx[N];\n", t1,
                          "int nseenexo[N];\n")
      }
    }

    ## parameter matrices/vectors
    matrows <- sapply(parmats[[1]], nrow)
    matcols <- sapply(parmats[[1]], ncol)
    if("rho" %in% names(nfree)){
      matrows <- c(matrows, rho = matrows[["theta"]])
      matcols <- c(matcols, rho = matcols[["theta"]])
    }
    if("lvrho" %in% names(nfree)){
      matrows <- c(matrows, lvrho = matrows[["psi"]])
      matcols <- c(matcols, lvrho = matcols[["psi"]])
    }

    pmats <- vector("list", length(matrows))
    for(i in 1:length(pmats)){
        if(names(matrows)[i] == "lambda"){
            tmpmat <- model@Model@GLIST$lambda
            pmats[[i]] <- array(tmpmat,
                                c(nrow(tmpmat), ncol(tmpmat), ngroups))
        } else {
            pmats[[i]] <- array(0, c(matrows[i], matcols[i], ngroups))
        }
    }
    names(pmats) <- names(matrows)

    ## monitored parameters
    monitors <- with(partable[partable$mat != "",], unique(mat))

    ## these are passed in as data in stan, so are the "frames"
    tpnames <- names(pmats)
    names(pmats) <- paste0(names(pmats), "frame")

    ## declare data variables and defined params
    datdecs <- tpdecs <- tpeqs <- ""
    for(i in 1:length(tpnames)){
      tmpdim <- dim(pmats[[i]])
      datdecs <- paste0(datdecs, t1, "real ",
                        names(pmats)[i], "[", tmpdim[1],
                        ",", tmpdim[2], ",", tmpdim[3], "];\n")
      tpdecs <- paste0(tpdecs, t1, "real ",
                       tpnames[i], "[", tmpdim[1],
                       ",", tmpdim[2], ",", tmpdim[3], "];\n")
      if(tpnames[i] == "theta"){
        tpdecs <- paste0(tpdecs, t1, "matrix[", ny,
                         ",", ny, "] thetld[", tmpdim[3],
                         "];\n")
      }
      tpeqs <- paste0(tpeqs, t1, tpnames[i], " = ",
                      names(pmats)[i], ";\n")
    }
    tpdecs <- paste0(tpdecs, t1, "real mu[N,", nov, "];\n")

    if(any(partable$mat == "def")){
      ndecs <- sum(partable$mat == "def" &
                   partable$group == 1)
      tpdecs <- paste0(tpdecs, "\n", t1, "real def[", ndecs, ", 1, ",
                       ngroups, "];\n")
    }
    
    if(nlv + n.psi.ov > 0){
      tpdecs <- paste0(tpdecs, t1, "matrix[N,", (nlv + n.psi.ov), "] eta;\n")
      if(nlvno0 < nlv){
        tpdecs <- paste0(tpdecs, t1, "vector[", (nlv + n.psi.ov),
                         "] mueta[N];\n")
      }
      tpdecs <- paste0(tpdecs, "\n", t1,
                       "eta = rep_matrix(0, N, ", (nlv + n.psi.ov),
                       ");\n")
    }
    
    ## if no beta, define it as 0 matrix
    if(!("beta" %in% tpnames)){
        matrows <- c(matrows, beta = matrows[["psi"]])
        matcols <- c(matcols, beta = matcols[["psi"]])
        pmats <- c(pmats, list(beta = array(0, c(matrows[["psi"]],
                                                 matcols[["psi"]],
                                                 ngroups))))
        datdecs <- paste0(datdecs, t1, "real beta[",
                          matrows[["psi"]], ",", matcols[["psi"]],
                          ",", ngroups, "];\n")
    }

    ## add cholesky decomp of theta matrix (and psi for nov.x);
    ## non-eXo ov vars sometimes show up in psi, so handle that as well.
    TPS <- paste0(TPS, t1, "}\n\n")
    TPS <- paste0(TPS, t1, "for(j in 1:", ngroups, "){\n")
    if(any(partable$mat == "theta")){
      if(n.psi.ov > 0){
        for(i in 1:length(yind)){
          for(j in i:length(yind)){
            TPS <- paste0(TPS, t2, "thetld[j,", i, ",", j, "] = ",
                          "theta[", yind[i], ",", yind[j], ",j];\n")
          }
        }
        TPS <- paste0(TPS, t2, "thetld[j] = fill_lower(thetld[j]);\n")
      } else {
        TPS <- paste0(TPS, t2, "thetld[j] = fill_lower(to_matrix(",
                      "theta[,,j]));\n")
      }
      TPS <- paste0(TPS, t2, "thetld[j] = cholesky_decompose(",
                    "thetld[j]);\n")
    }
    if(dumov & !model@Options$fixed.x &
       !all(parmattable$lambda == diag(nrow(parmattable$lambda)))){
      TPS <- paste0(TPS, t2, "alpha[dummylv,1,j] = to_array_1d(inverse((to_matrix(lambda[,,j]) * inverse(diag_matrix(rep_vector(1.0, ", (nlv + n.psi.ov), ")) - to_matrix(beta[,,j])))[dummyov,dummylv]) * to_vector(to_array_1d(alpha[dummylv,1,j])")
      TPS <- paste0(TPS, "));\n")
    }
    
    TPS <- paste0(TPS, t1, "}\n")

    TPS <- paste0("transformed parameters{\n", tpdecs, "\n",
                  tpeqs, TXT2, "\n\n", TPS,
                  "\n}\n\n")
    datablk <- paste0(datablk, datdecs, "}\n\n")

    standata <- c(standata, pmats)
        
    out <- c(out, list(data=standata))
  }

  funblk <- "functions{\n"
  if((nlv + n.psi.ov) > 0){
    funblk <- paste0(funblk, t1, "#include 'sem_mean.stan' \n")
    if(nlvno0 < nlv){
      funblk <- paste0(funblk, t1, "#include 'sem_mean_eta.stan' \n")
    }
    if(miss.psi){
      funblk <- paste0(funblk, t1, "#include 'sem_lv_missing.stan' \n")
    } else {
      funblk <- paste0(funblk, t1, "#include 'sem_lv.stan' \n")
    }
  }
  funblk <- paste0(funblk, t1, "#include 'fill_lower.stan' \n")
  ## could insert other functions as needed
  funblk <- paste0(funblk, "}\n\n")

  fullmodel <- paste0(funblk, datablk, parblk, TPS, out$model, "\n")

  ## insert function files, similar to brms approach:
  tmp <- tempfile(fileext = ".stan")
  cat(fullmodel, file = tmp)
  isystem <- system.file("stanfuns", package = "blavaan")
  out$model <- rstan::stanc_builder(file = tmp, isystem = isystem,
                                    obfuscate_model_name = TRUE)$model_code

  out <- c(out, list(monitors = monitors, pxpartable = partable))

  out
}

coeffun_stan <- function(lavpartable, rsob, fun = "mean") {
  ## Extract posterior means from coda.samples() object.
  ## rsob is the result of rstan().
  rssumm <- rstan::summary(rsob)
  rsmcmc <- as.array(rsob)

  ## posterior means:
  if(fun == "mean"){
    b.est <- rssumm$summary[,"mean"]
  } else if(fun == "median"){
    b.est <- rssumm$summary[,"50%"]
  }

  ## move "free" parameters from rho to theta
  rhopars <- grep("rho", lavpartable$mat)
  if(length(rhopars) > 0){
    for(i in 1:length(rhopars)){
      idx <- rhopars[i]
      matname <- ifelse(lavpartable$mat[idx] == "rho", "theta", "psi")
      newidx <- which(lavpartable$mat == matname &
                      lavpartable$row == lavpartable$row[idx] &
                      lavpartable$col == lavpartable$col[idx] &
                      lavpartable$group == lavpartable$group[idx])

      tmpfree <- lavpartable$free[idx]

      lavpartable$free[idx] <- 0L
      lavpartable$free[newidx] <- tmpfree
    }
  }
  lavord <- order(lavpartable$id)
  lavpartable <- lapply(lavpartable, function(x) x[lavord])
  
  ## from stan to partable
  ptnames <- with(lavpartable, paste0(mat, "[", row, ",", col, ",",
                                      group, "]"))
  cmatch <- match(ptnames, names(b.est), nomatch=0)
  lavpartable$est[cmatch > 0] <- b.est[cmatch]
  lavpartable$psrf[cmatch > 0] <- rssumm$summary[cmatch,"Rhat"]

  ## NB: order of parameters in mcmc array differs from order
  ##     of parameters in summary()
  lavpartable$stanpnum <- match(ptnames, names(rsmcmc[1,1,]), nomatch=0)
  lavpartable$stansumnum <- match(ptnames, rownames(rssumm$summary), nomatch=0)

  sdvec <- rssumm$summary[cmatch, "sd"]

  ## vcorr
  draw_mat <- as.matrix(rsob)
  cmatch <- match(ptnames[lavpartable$free > 0][order(lavpartable$free[lavpartable$free > 0])], colnames(draw_mat))
  vcorr <- cor(draw_mat[,cmatch])

  svmatch <- match(colnames(vcorr), names(sdvec), nomatch = 0)
  sdvec <- sdvec[svmatch]
  
  ## convert to list
  lavpartable <- as.list(lavpartable, seq(ncol(lavpartable)))

  list(x = lavpartable$est[lavpartable$free > 0][order(lavpartable$free[lavpartable$free > 0])],
       lavpartable = lavpartable,
       vcorr = vcorr,
       sd = sdvec, stansumm = rssumm$summary)
}
