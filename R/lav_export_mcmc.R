lav2mcmc <- function(model, lavdata = NULL, cp = "srs", lv.x.wish = FALSE, dp = NULL, n.chains = 1, mcmcextra = "", inits = "prior", blavmis = "da", pta = NULL, target = "stan") {
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
  
  eqop <- ifelse(target == "stan", "=", "<-")
  commop <- ifelse(target == "stan", "// ", "# ")
  eolop <- ifelse(target == "stan", ";", "")
  if(length(dp) == 0) dp <- dpriors(target = target)
  
  ## get names of ovs before we add phantom variables
  old.pta <- lav_partable_attributes(partable = partable, pta = pta)
  old.vnames <- old.pta$vnames
  ngroups <- old.pta$ngroups
  orig.ov.names <- old.vnames$ov[[1]]; nov <- length(orig.ov.names)
  orig.lv.names <- old.vnames$lv[[1]]; orig.lv.names.x <- old.vnames$lv.x[[1]]
  ## so ordering stays consistent:
  orig.lv.names <- c(orig.lv.names[orig.lv.names %in% orig.lv.names.x],
                     orig.lv.names[!(orig.lv.names %in% orig.lv.names.x)])
  orig.ov.names.x <- old.vnames$ov.x[[1]]
  nlvx <- length(orig.lv.names.x)
  
  ## if lv.x.wish and default prior, change df parameter for this model
  if(lv.x.wish & nlvx > 1 & dp[["ibpsi"]] == dpriors(target = target)[["ibpsi"]]){
    if(target == "jags"){
      dp[["ibpsi"]] <- paste("dwish(iden,", length(orig.lv.names.x) + 1, ")", sep="")
    } else {
      dp[["ibpsi"]] <- paste("wishart(",
                             length(orig.lv.names.x) + 1,
                             ",iden)", sep="")
    }
  }

  ## set up mvs with fixed 0 variances (single indicators of lvs)
  partable <- set_mv0(partable, orig.ov.names, ngroups)
  ## add necessary phantom lvs/mvs to model:
  partable <- set_phantoms(partable, orig.ov.names, orig.lv.names, orig.ov.names.x, orig.lv.names.x, cp, cp, lv.x.wish, ngroups)
  facovs <- partable$facovs
  partable <- partable$partable
  ## set equality constraints for phantom variances
  partable <- set_phanvars(partable, orig.ov.names, orig.lv.names, cp, cp, ngroups)
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
  TPS <- ""

  if(blavmis == "da"){
    TXT <- paste(TXT, t1,
                 "for(i in 1:N) {\n", sep="")
  } else {
    TXT <- paste(TXT, t1, "for(i in 1:nrows) {\n", sep="")
  }

  ## Second object for priors/constraints
  TXT2 <- ""

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
  eqlabs <- partable$rhs[partable$op %in% c("==", ":=")]
  eqplabs <- partable$lhs[partable$op %in% c("==", ":=")]
  eqplabs <- eqplabs[eqplabs %in% partable$label]
  eqlabs <- c(eqlabs, eqplabs)

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

  ## TODO We now attach equality constraints to these tables, could
  ##      also deal with inequality constraints.
  ## Smaller partables for different parameter types +
  ## dimensions of parameter matrices (for initial values)
  ovintercepts <- partable[ovi,]
  ovintercepts <- rbind(ovintercepts, partable[which(partable$op %in% c("==", ":=")),])
  lvintercepts <- partable[lvi,]
  lvintercepts <- rbind(lvintercepts, partable[which(partable$op %in% c("==", ":=")),])
  loadings <- partable[load,]
  loadings <- rbind(loadings, partable[which(partable$op %in% c("==", ":=")),])
  regressions <- partable[reg,]
  regressions <- rbind(regressions, partable[which(partable$op %in% c("==", ":=")),])

  ## for missing=="fi", to model variables on rhs of regression
  ovreg <- unique(regressions$rhs[regressions$rhs %in% ov.names])
  ovcol <- which(ov.names %in% ovreg)

  ## Define univariate distributions of each observed variable
  ## Loop if everything is continuous
  if(length(ov.ord) == 0){
    if(blavmis == "da"){
      for(j in 1:nmvs){
        ## decide whether we need px on ovs/lvs by searching
        ## for covariances:
        mvvar <- which(partable$op == "~~" &
                       partable$lhs == ov.names[j] &
                       partable$lhs == partable$rhs &
                       !(grepl("star", partable$mat)) &
                       partable$group == 1)
        tvname <- partable$mat[mvvar]
        mvcovs <- length(which(grepl(".phant", partable$lhs) &
                               partable$op == "=~" &
                               partable$rhs == ov.names[j]))

        if(mvcovs > 0){
          tvname <- paste(tvname, "star", sep="")
        }
          
        TXT <- paste(TXT, t2, ov.names[j], sep="")
        if(target == "stan"){
          TXT <- paste(TXT, "[i] ~ normal(mu[i,", j, "], sqrt(",
                       sep="")
        } else {
          TXT <- paste(TXT, "[i] ~ dnorm(mu[i,", j, "], 1/",
                       sep="")
        }
        TXT <- paste(TXT, tvname, "[", partable$row[mvvar], ",",
                     partable$col[mvvar], ",g[i]])", sep="")
        if(target == "stan"){
          TXT <- paste(TXT, ")", sep="")
        }
        TXT <- paste(TXT, eolop, "\n", sep="")
      }
    } else {
      if(target == "stan") stop("blavaan ERROR: missing data in stan not yet available.")
      TXT <- paste(TXT, t2, "yvec[i] ~ dnorm(mu[sub[i], mv[i]], 1/",
                   tvname, "[mv[i], mv[i], g[i]])\n", sep="")
      # now close this loop and start the usual one
      TXT <- paste(TXT, t1, "}\n\n", sep="")
      TXT <- paste(TXT, t1, "for(i in 1:N) {\n", sep="")
      # now model regressions on rhs
      if(length(ovreg) > 0){
        colidx <- match(ovreg, ov.names)
        for(j in 1:length(ovreg)){
          TXT <- paste(TXT, t2, ovreg[j],
                       "[i] ~ dnorm(mu[i,", colidx[j], "], 1/", tvname,
                       "[", colidx[j], ",g[i]])", eolop, "\n", sep="")
        }
      }
    }
  } else {
    if(target == "stan") stop("blavaan ERROR: ordinal models cannot yet be exported to stan")
    # TODO revisit "fi" approach to missing data
    if(blavmis == "fi") stop("blavaan ERROR: missing='fi' not yet supported for ordinal data")
    for(j in 1:nmvs){
      if(ov.names[j] %in% ov.ord){
        tvvar <- which(partable$op == "~~" &
                       partable$lhs == ov.names[j] &
                       partable$lhs == partable$rhs &
                       !(grepl("star", partable$mat)) & # needed?
                       partable$group == 1)
        tvname <- partable$mat[tvvar]
        ord.num <- match(ov.names[j], ov.ord)
        TXT <- paste(TXT, t2, "ones[i,", ord.num, "] ~ dbern(probs[i,", ord.num,
                     ",", ov.names[j], "[i]])\n", sep="")

        ## category probs from pnorm
        taus <- which(partable$lhs == ov.names[j] &
                      partable$op == "|" &
                      partable$group == 1)
        TXT <- paste(TXT, t2, "probs[i,", ord.num, ",1] ", eqop,
                     " pnorm(tau[", partable$row[taus[1]],
                     ",", partable$col[taus[1]], ",g[i]], mu[i,", j, "], 1/", tvname, "[", j, ",", j, ",g[i]])\n", sep="")
        if(ncats[ord.num] > 2){
          for(k in 2:(ncats[ord.num] - 1)){
            TXT <- paste(TXT, t2, "probs[i,", ord.num, ",", k, "] ",
                         eqop, " pnorm(tau[",
                         partable$row[taus[k]], ",", partable$col[taus[k]], ",g[i]], mu[i,", j, "], 1/", tvname, "[", j, ",",
                         j, ",g[i]]) - sum(probs[i,", ord.num, ",1:", (k-1), "])\n",
                         sep="")
          }
        }
        TXT <- paste(TXT, t2, "probs[i,", ord.num, ",", ncats[ord.num],
                     "] ", eqop, " 1 - sum(probs[i,", ord.num, ",1:", ncats[ord.num]-1, "])\n",
                     sep="")
      } else {
        TXT <- paste(TXT, t2, ov.names[j], "[i] ~ dnorm(mu[i,", j, "],", tvname, "[", partable$row[tvvar], ",", partable$col[tvvar], ",g[i]])\n", sep="")
      }
    }
  }

  ## Define mean of each observed variable
  ## This assumes that the data matrix passed to jags
  ## is ordered in the same way as ov.names.nox.
  ## data would be cbind(ov.names.nox, ov.names.x)
  TPS <- paste(TPS, t1, commop, "mu definitions\n", t1,
               "for(i in 1:N) {", sep="")
  for(i in 1:nmvs) {
    ov.idx <- i
    if(i > 1) TPS <- paste(TPS, eolop, sep="")
    TPS <- paste(TPS, "\n", t2, "mu[i,", ov.idx, "] ", eqop, " ",
                 sep="")

    ## find rhs for this observed variable
    ## 1. intercept?
    
    ## Always include intercept parameters, fix to zero
    ## if they are not desired
    int.idx <- which(partable$op == "~1" &
                     partable$lhs == ov.names[i] &
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
                     loadings$rhs == ov.names[i] &
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
    r.idx <- which(regressions$lhs == ov.names[i] &
                   regressions$group == 1)
    for(j in r.idx) {
      ## what is the rhs?
      rhs <- regressions$rhs[j]
      if(rhs %in% lv.names) {
        RHS <- paste("eta[i,", match(rhs, lv.names), "]", sep="")
      } else if(rhs %in% orig.ov.names) {
        RHS <- paste(rhs, "[i]", sep="")
      }
      
      ## deal with fixed later
      TPS <- paste(TPS, " + ", regressions$mat[j], "[",
                   regressions$row[j], ",", regressions$col[j],
                   ",g[i]]*", RHS, sep="")
    }

    ## 4. residual variances now handled separately

    ## 5. TODO thresholds
    if(length(ov.ord) > 0){

    }
  }
  TPS <- paste(TPS, eolop, "\n", sep="")

  ## lvs
  if(length(lv.names) > 0L) {
    TXT <- paste(TXT, "\n", t2,
                 commop, "lvs", sep="")

    lvstart <- 1
    if(lv.x.wish & nlvx > 1){
      lvstart <- nlvx + 1

      TXT <- paste(TXT, "\n", t2,
                   "eta[i,1:", nlvx, "] ~ ", sep="")
      if(target == "stan"){
        TXT <- paste(TXT, "multi_normal_prec(", sep="")
      } else {
        TXT <- paste(TXT, "dmnorm(", sep="")
      }
      TXT <- paste(TXT, "mu_eta[i,1:", nlvx, "], ibpsi[1:",
                   nlvx,",1:", nlvx, ",g[i]])", eolop, sep="")
    }
    
    nlv <- length(lv.names)

    eta.eq <- rep(NA, nlv)

    if(nlv >= lvstart){
      for(j in lvstart:nlv) {
        psi.free.idx <- which(partable$group == 1 &
                              partable$op == "~~" &
                              partable$lhs == partable$rhs &
                              partable$lhs == lv.names[j] &
                              !grepl("star", partable$mat))

        if(length(psi.free.idx) != 1L) {
          stop("blavaan ERROR: parameter for residual variance ",
               lv.names[j], " not found")
        }

        lv.var <- which(partable$lhs == lv.names[j] &
                        partable$rhs == lv.names[j] &
                        partable$op == "~~")
        if(any(partable$free[lv.var] == 0 & partable$ustart[lv.var] == 0)){
          TXT <- paste(TXT, "\n", t2, "eta[i,", j, "] ", eqop,
                       " mu_eta[i,", j, "]", eolop, sep="")
          ## now change ustart to 1000 so no divide by 0 in jags
          partable$ustart[lv.var] <- 1000
        } else {
          lvcovs <- length(which(partable$lhs == lv.names[j] &
                                 grepl(".phant", partable$rhs) &
                                 partable$op == "~"))

          pvname <- ifelse(lvcovs > 0, "psistar", "psi")
            
          TXT <- paste(TXT, "\n", t2,
                       ## TODO check for alternative distribution?
                       "eta[i,", j, "] ~ ", sep="")
          if(target == "stan"){
            TXT <- paste(TXT, "normal(mu_eta[i,", j, "], sqrt(",
                         sep="")
          } else {
            TXT <- paste(TXT, "dnorm(mu_eta[i,", j, "], 1/", sep="")
          }
          TXT <- paste(TXT, pvname, "[",
                       partable$row[psi.free.idx], ",", partable$col[psi.free.idx],
                       ",g[i]])", sep="")
          if(target == "stan"){
            TXT <- paste(TXT, ")", eolop, sep="")
          }
        }
      } # j
    } # if

    ## After lv distributions are defined, now define means/regressions:
    for(j in 1:nlv) {
      if(j > 1) TPS <- paste(TPS, eolop, sep="")

      TPS <- paste(TPS, "\n", t2,
                   ## TODO check for alternative distribution
                   ## in parameter table.
                   "mu_eta[i,", j, "] ", eqop, " ", sep="")

      ## lhs elements regression
      ## 1. intercept? (even exogenous can have an intercept)
      int.idx <- which(lvintercepts$group == 1 &
                       lvintercepts$op == "~1" &
                       lvintercepts$lhs == lv.names[j])
      if(length(int.idx) == 1L) {
        ## fixed or free?
        TPS <- paste(TPS, lvintercepts$mat[int.idx], "[", lvintercepts$row[int.idx],
                     ",", lvintercepts$col[int.idx], ",g[i]]", sep="")
      } else { # no intercept, say '0', so we always have rhs
        TPS <- paste(TPS, "0", sep="")
      }

      if(lv.names[j] %in% lv.nox){
        ## 2. loadings?
        lam.idx <- which(loadings$op == "=~" &
                         loadings$rhs == lv.names[j] &
                         loadings$group == 1)
        if(length(lam.idx) > 0){
          for(k in lam.idx){
            TPS <- paste(TPS, " + ", loadings$mat[k], "[", loadings$row[k],
                         ",", loadings$col[k], ",g[i]]*eta[i,",
                         match(loadings$lhs[k], lv.names), "]", sep="")
          } # end k loop
        }                         
        
        ## 3. regressions?
        rhs.idx <- which(regressions$lhs == lv.names[j] &
                         regressions$op == "~" &
                         regressions$group == 1)
        np <- length(rhs.idx)
        if(np > 0){ # there could be none if we have higher-order factors
          for(p in rhs.idx){
            TPS <- paste(TPS, " + ", regressions$mat[p], "[", regressions$row[p],
                         ",", regressions$col[p], ",g[i]]", sep="")

            ## Is the rhs an lv or ov?
            lvmatch <- match(regressions$rhs[p], lv.names)
            if(is.na(lvmatch)){
              TPS <- paste(TPS, "*", regressions$rhs[p], "[i]", sep="")
            } else {
              TPS <- paste(TPS, "*eta[i,", lvmatch, "]", sep="")
            }
          } # end p loop
        } # end if np
      } # end if

      ## 4. lv variances now handled separately
    } # end j loop
    TPS <- paste(TPS, eolop, sep="")
  } # end if length(lv.names)

  ## end of main model specification (still need priors + equality constraints)
  TXT <- paste(TXT, "\n", t1, "}", sep="")
  TPS <- paste(TPS, "\n", t1, "}", sep="")

  ## priors/constraints
  TXT2 <- set_parvec(TXT2, partable, dp, cp, lv.x.wish, orig.lv.names.x, target)
  partable$prior <- TXT2$partable$prior
  partable$freeparnums <- TXT2$partable$freeparnums
  TXT3 <- TXT2$TXT3
  TXT2 <- TXT2$TXT2

  ## end of model
  if(target == "jags"){
    TXT <- paste(TXT, "\n\n", TPS, "\n\n", t1, "# Assignments from parameter vector & equality constraints", TXT2, TXT3, sep="")
  }

  ## extra stuff from the user, formatted to look nice-ish
  if("syntax" %in% names(mcmcextra)){
    mcmcextra <- unlist(strsplit(mcmcextra$syntax, "\n"))
    mcmcextra <- gsub("^\\s+|\\s+$", "", mcmcextra)
    mcmcextra <- paste(t1, mcmcextra, sep="", collapse="\n")
    TXT <- paste(TXT, "\n", mcmcextra, "\n", sep="")
  }
  TXT <- paste(TXT, "\n", sep="")
  if(target == "jags") TXT <- paste0(TXT, "}\n")
  
  out <- TXT
  if(target == "stan") out <- paste0(out, TXT3, "\n}")
  class(out) <- c("lavaan.character", "character")
  out <- list(model = out, inits = NA)
    
  ## Initial values
  if(inits != "jags"){
      if(target == "stan") stop("blavaan ERROR: random inits not yet available for stan")
      inits <- set_inits(partable, cp, cp, n.chains, inits)
      out$inits <- inits
  }
        
  ## Now add data for jags if we have it
  datablk <- paste0("data{\n", t1, "int N;\n", t1, "int g[N];\n")
  if(!is.null(lavdata) | class(model)[1]=="lavaan"){
    if(class(model)[1] == "lavaan") lavdata <- model@Data
    ntot <- sum(unlist(lavdata@norig))
    ## pick up exogenous x's
    tmpnmvs <- length(orig.ov.names)

    y <- lapply(1:tmpnmvs, function(x) rep(NA,ntot))
    g <- rep(NA, ntot)
    nX <- ncol(lavdata@X[[1]])
    for(k in 1:ngroups){
      for(j in 1:nX){
        y[[j]][lavdata@case.idx[[k]]] <- lavdata@X[[k]][,j]
      }
      if(tmpnmvs > nX){
        for(j in 1:(tmpnmvs - nX)){
          y[[j + nX]][lavdata@case.idx[[k]]] <- lavdata@eXo[[k]][,j]
        }
      }
      g[lavdata@case.idx[[k]]] <- k
    }
    names(y) <- orig.ov.names

    ## stan data block
    for(j in 1:tmpnmvs){
      datablk <- paste0(datablk, t1, "vector[N] ",
                        orig.ov.names[j], ";\n")
    }
    
    ## remove fully deleted rows
    ymat <- matrix(unlist(y), ntot, tmpnmvs)
    nas <- which(apply(is.na(ymat), 1, sum) == tmpnmvs)

    if(length(nas) > 0){
      y <- lapply(y, function(x) x[-nas])
      g <- g[-nas]
      ntot <- sum(unlist(lavdata@nobs))
      ymat <- ymat[-nas,]
    }

    if(blavmis == "fi"){
      ## variables on rhs of regression, in case missing
      y <- y[names(y) %in% ovreg]
    }

    jagsdata <- c(y, list(g=g, N=ntot))
    if(any(partable$op == "|")){
      jagsdata <- c(jagsdata, list(ones = matrix(1, ntot, tmpnmvs)))
    }

    if(blavmis == "fi"){
      ## keep only modeled y's not on rhs of regression
      matvars <- which(orig.ov.names %in% ov.names &
                       !(orig.ov.names %in% ovreg))
      ymat <- ymat[,matvars]
      nmvs <- length(matvars)
  
      ydf <- data.frame(y=as.numeric(ymat),
                        g=rep(g, nmvs),
                        sub=rep(1:ntot, nmvs),
                        mv=rep(matvars, each=ntot))
      
      ydf <- subset(ydf, !is.na(ydf$y) & !(ydf$mv %in% ovcol))
      ## sub index excluding completely missing observations
      ydf$sub <- as.numeric(as.factor(ydf$sub))

      jagsdata <- c(jagsdata, list(yvec=ydf$y, sub=ydf$sub,
                                   mv=ydf$mv, nrows=nrow(ydf)))
      jagsdata$g <- ydf$g
      jagsdata$N <- max(ydf$sub)
    }

    ## parameter matrices/vectors
    matrows <- with(partable[partable$mat != "",], tapply(row, mat, max, na.rm=TRUE))
    matcols <- with(partable[partable$mat != "",], tapply(col, mat, max, na.rm=TRUE))
    ngrp <- max(partable$group, na.rm=TRUE)

    pmats <- vector("list", length(matrows))
    for(i in 1:length(pmats)){
        pmats[[i]] <- array(0, c(matrows[i], matcols[i], ngrp))
    }
    names(pmats) <- names(matrows)

    ## replace parameter entries with NA to please jags
    if(target == "jags"){
      for(i in 1:nrow(partable)){
        if(partable$mat[i] == "") next

        wmat <- match(partable$mat[i], names(pmats))
        pmats[[wmat]][partable$row[i], partable$col[i], partable$group[i]] <- NA
      }
    }

    ## monitored parameters
    monitors <- with(partable[partable$mat != "",],
                     paste(mat, "[", row, ",", col, ",", group, "]",
                           sep=""))
    
    ## inferential covariances under fa priors
    if(cp == "fa" & length(facovs) > 0){
        if(nrow(facovs) > 0){
            for(i in 1:nrow(facovs)){
                wmat <- match(facovs$mat[i], names(pmats))

                if(target == "jags"){
                  pmats[[wmat]][facovs$row[i], facovs$col[i], facovs$group[i]] <- NA
                }

                monitors <- c(monitors, paste(facovs$mat[i], "[", facovs$row[i], ",",
                                              facovs$col[i], ",",
                                              facovs$group[i], "]", sep=""))
            }
            ## re-add fa covariances, for sending back to lavaan
            partable <- partable[,match(names(partable), names(facovs), nomatch=0)]
            partable <- rbind(partable, facovs)
        }
    }

    ## these are passed in as data in stan, so are the "frames"
    if(target == "stan"){
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
        tpeqs <- paste0(tpeqs, t1, tpnames[i], " = ",
                        names(pmats)[i], ";\n")
      }
      tpdecs <- paste0(tpdecs, t1, "real mu[N,", tmpnmvs, "];\n")
      if(length(lv.names) > 0){
        tpdecs <- paste0(tpdecs, t1, "real mu_eta[N,",
                         length(lv.names), "];\n")
      }

      ## define rho matrices as needed
      if(any(grepl("thetastar", tpnames))){
        thetloc <- which(names(pmats) == "thetastarframe")
        thetdim <- dim(pmats[[thetloc]])
        rhoframe <- array(0, thetdim)
        pmats <- c(pmats, list(rhoframe = rhoframe))
        datdecs <- paste0(datdecs, t1, "real ",
                          "rhoframe[", thetdim[1], ",",
                          thetdim[2], ",", thetdim[3], "];\n")
        tpdecs <- paste0(tpdecs, t1, "real ",
                          "rho[", thetdim[1], ",",
                          thetdim[2], ",", thetdim[3], "];\n")
        tpeqs <- paste0(tpeqs, t1, "rho = rhoframe;\n")
      }

      if(any(grepl("psistar", tpnames))){
        psiloc <- which(names(pmats) == "psistarframe")
        psidim <- dim(pmats[[psiloc]])
        lvrhoframe <- array(0, psidim)
        pmats <- c(pmats, list(lvrhoframe = lvrhoframe))
        datdecs <- paste0(datdecs, t1, "real ",
                          "lvrhoframe[", psidim[1], ",",
                          psidim[2], ",", psidim[3], "];\n")
        tpdecs <- paste0(tpdecs, t1, "real ",
                          "lvrho[", psidim[1], ",",
                          psidim[2], ",", psidim[3], "];\n")
        tpeqs <- paste0(tpeqs, t1, "lvrho = lvrhoframe;\n")
      }

      TPS <- paste0("transformed parameters{\n", tpdecs, "\n",
                    tpeqs, TXT2, "\n\n", TPS,
                    "\n}\n\n")
      datablk <- paste0(datablk, datdecs, "}\n\n")
    }

    jagsdata <- c(jagsdata, pmats)
    
    ## identity matrix for wishart prior
    ## TODO allow user to specify this matrix
    ## or could be specified via mcmcextra argument?
    if(lv.x.wish & length(orig.lv.names.x) > 1){
      iden <- diag(length(orig.lv.names.x))
      jagsdata <- c(jagsdata, list(iden=iden))
    }
    
    out <- c(out, list(data=jagsdata))
  }

  if(target == "stan"){
    nparms <- max(partable$freeparnums, na.rm = TRUE)
    parmblk <- paste0("parameters{\n", t1, "vector[",
                      nparms, "] parvec;\n")
    if(length(lv.names) > 0){
      parmblk <- paste0(parmblk, t1,
                      "real eta[N,", length(lv.names),
                      "];\n")
    }
    parmblk <- paste0(parmblk, "}\n\n")    
    out$model <- paste0(datablk, parmblk, TPS, out$model)
  }
  
  out <- c(out, list(monitors = monitors, pxpartable = partable))
    
  out
}

coeffun <- function(lavpartable, pxpartable, rjob, fun = "mean") {
  ## Extract posterior means from coda.samples() object.
  ## jagmod is the result of lav2jags;
  ## rjob is the result of run.jags() applied to
  ## the jags model.

  ## remove any rhos or unnamed parameters
  pxpartable <- pxpartable[!is.na(pxpartable$id) & pxpartable$op != "==",]
  pxnames <- paste(pxpartable$mat, "[", pxpartable$row, ",", pxpartable$col,
                   ",", pxpartable$group, "]", sep="")
  pxpartable$pxnames <- pxnames

  ## posterior means:
  if(fun == "mean"){
    b.est <- rjob$summary$statistics[,"Mean"]
  } else if(fun == "median"){
    b.est <- rjob$hpd[,"Median"]
  }

  ## from jags to pxpartable
  cmatch <- match(pxnames, names(b.est), nomatch=0)
  pxpartable$est <- b.est[cmatch]
  pxpartable$psrf <- rep(NA, length(pxpartable$free))
  if(length(rjob$mcmc) > 1){
    psrfmatch <- match(pxnames, rownames(rjob$psrf$psrf))
    pxpartable$psrf <- rjob$psrf$psrf[psrfmatch,1]
  }
  pxpartable$jagpnum <- cmatch

  ## from pxpartable to lavpartable
  ## first check for px parameters with "free" labels (fa priors)
  pxmats <- c("theta", "psi")
  for (j in 1:length(pxmats)){
    stars <- which(grepl(paste(pxmats[j], "star", sep=""),
                          pxpartable$mat) &
                    pxpartable$free > 0)
    if(length(stars) > 0){
      for(i in 1:length(stars)){
        infpar <- which(pxpartable$mat == pxmats[j] &
                        pxpartable$row == pxpartable$row[stars[i]] &
                        pxpartable$col == pxpartable$col[stars[i]] &
                        pxpartable$group == pxpartable$group[stars[i]])
        pxpartable$free[infpar] <- pxpartable$free[stars[i]]
        pxpartable$free[stars[i]] <- 0
      }
    }
  }

  ptmatch <- match(lavpartable$free[lavpartable$free > 0], pxpartable$free)
  if("est" %in% names(pxpartable)){
    ## to handle do.fit = FALSE
    lavpartable$est[lavpartable$free > 0] <- pxpartable$est[ptmatch]
  }
  lavpartable$psrf <- rep(NA, length(lavpartable$free))
  lavpartable$psrf[lavpartable$free > 0] <- pxpartable$psrf[ptmatch]
  lavpartable$prior[lavpartable$free > 0] <- pxpartable$prior[ptmatch]
  lavpartable$pxnames[lavpartable$free > 0] <- pxpartable$pxnames[ptmatch]
  lavpartable$jagpnum[lavpartable$free > 0] <- pxpartable$jagpnum[ptmatch]

  ## defined variables
  defmatch <- which(pxpartable$op == ":=")
  if(length(defmatch) > 0){
    lavpartable$est[lavpartable$op == ":="] <- pxpartable$est[defmatch]
    lavpartable$psrf[lavpartable$op == ":="] <- pxpartable$psrf[defmatch]
    lavpartable$pxnames[lavpartable$op == ":="] <- pxpartable$pxnames[defmatch]
    lavpartable$jagpnum[lavpartable$op == ":="] <- pxpartable$jagpnum[defmatch]
  }
  
  ## NB this automatically removes fixed parameters, just
  ##    like the psrf
  lmatch <- match(lavpartable$pxnames[lavpartable$free > 0],
                  rownames(rjob$crosscorr))
  vcorr <- rjob$crosscorr[lmatch, lmatch]
  smatch <- match(lavpartable$pxnames[lavpartable$free > 0 | lavpartable$op == ":="],
                  rownames(rjob$summary$statistics),
                  nomatch=0)
  sdvec <- rjob$summary$statistics[smatch, "SD"]

  list(x = lavpartable$est[lavpartable$free > 0],
       lavpartable = lavpartable,
       vcorr = vcorr,
       sd = sdvec)
}
