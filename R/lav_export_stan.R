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
  partable <- set_stancovs(partable, orig.ov.names, orig.ov.names.x,
                           dp)

  ## ensure group parameters are in order, for parameter indexing:
  partable <- partable[order(partable$group),]
  ## get parameter table attributes 
  pta <- lav_partable_attributes(partable = partable, pta = NULL)
  vnames <- pta$vnames; nvar <- pta$nvar; nfac <- pta$nfac
  ov.names.nox <- vnames$ov.nox[[1]]; nov.nox <- length(ov.names.nox)
  ov.names.x <- vnames$ov.x[[1]]; nov.x <- length(ov.names.x)
  ov.ord <- vnames$ov.ord[[1]]
  lv.nox <- vnames$lv.nox[[1]]
  lv.names <- vnames$lv[[1]]
  ## ensure that lv.x names always come first (so we can possibly use dmnorm)
  lv.names <- c(lv.names[lv.names %in% orig.lv.names.x],
                lv.names[!(lv.names %in% orig.lv.names.x)])
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

  ## Decide whether we need to model exogenous x's
  if(length(ov.names.x) > 0){ # & !is.na(ov.names.x)){
    exotab <- partable[which(partable$lhs %in% old.vnames$ov.x[[1]]),]
    if(all(exotab$free==0)){
      nmvs <- nov.nox
      ov.names <- ov.names.nox
    } else {
      nmvs <- nov
      ov.names <- orig.ov.names
    }
  } else {
    nmvs <- nov.nox
    ov.names <- orig.ov.names[orig.ov.names %in% ov.names.nox]
  }
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

  ## NOTE We attach equality constraints to these tables, could
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

  ## number of free parameters per type, for stan parameter vectors
  ## (need separated so can use "lower" and "upper")
  parmats <- lavInspect(model)
  parconst <- attr(parmats, "header")
  
  ## so it is always a list of lists
  if(model@Data@ngroups == 1) parmats <- list(g1 = parmats)

  ## decide whether psi is diagonal, for faster matrix computations
  ## in stan
  diagpsi <- 0L
  if("psi" %in% names(parmats[[1]])){
    tmppsi <- parmats[[1]]$psi
    tmppsi <- tmppsi[lower.tri(tmppsi)]
    if(all(tmppsi == 0)) diagpsi <- 1L
  }
  
  nfree <- sapply(parmats, sapply, function(x){
    if(class(x)[1] == "lavaan.matrix.symmetric"){
      # off-diagonals handled via rho parameters!
      sum(diag(x) > 0)
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
    nrho <- max(partable$rhoidx[partable$mat == "rho"], na.rm = TRUE)# - nrhofix
    nfree <- c(nfree, rho = nrho)
    parblk <- paste0(parblk, t1, "vector<lower=0,upper=1>[",
                     nrho, "] rhofree;\n")
  }
  if(any(partable$mat == "lvrho")){
    #nlrhofix <- sum(sapply(parmats, function(x){
    #  sum(x$psi[lower.tri(x$psi)] %in% parconst$rhs)
    #}))
    nlrho <- max(partable$rhoidx[partable$mat == "lvrho"], na.rm = TRUE)# - nlrhofix
    nfree <- c(nfree, lvrho = nlrho)
    parblk <- paste0(parblk, t1, "vector<lower=0,upper=1>[",
                     nlrho, "] lvrhofree;\n")
  }

  if(nlv > 0){
    parblk <- paste0(parblk, t1, "matrix[N, ", nlv, "] eta", eolop,
                     "\n")
  }
  parblk <- paste0(parblk, "}\n\n")                     

  if(any(grepl("0", model@Data@Mp[[1]]$id))){
    stop("blavaan ERROR: missing data not yet handled in stan.")
    ## start off as if blavmis=="fi"
    TXT <- paste(TXT, t1, "for(i in 1:nrows) {\n", sep="")
  } else {
    TXT <- paste0(TXT, t1, "for(i in 1:N) {\n", t2,
                 "y[i] ~ multi_normal_cholesky(",
                 "to_vector(mu[i]), thetld[g[i]]);\n", t1,
                 "}\n\n", t1,
                 "eta ~ sem_lv_lpdf(alpha, beta, psi, g, ",
                 nlv, ", N, ", ngroups, ", ", diagpsi, ");\n")
  }

  
  ## for missing=="fi", to model variables on rhs of regression
  ovreg <- unique(regressions$rhs[regressions$rhs %in% ov.names])
  ovcol <- which(ov.names %in% ovreg)

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
  datablk <- paste0("data{\n", t1, "int N;\n", t1, "int g[N];\n")
  if(!is.null(lavdata) | class(model)[1]=="lavaan"){
    if(class(model)[1] == "lavaan") lavdata <- model@Data
    ntot <- sum(unlist(lavdata@norig))
    ## pick up exogenous x's
    tmpnmvs <- length(orig.ov.names)

    y <- matrix(NA, ntot, tmpnmvs) #lapply(1:tmpnmvs, function(x) rep(NA,ntot))
    g <- rep(NA, ntot)
    nX <- ncol(lavdata@X[[1]])
    for(k in 1:ngroups){
      for(j in 1:nX){
        y[lavdata@case.idx[[k]],j] <- lavdata@X[[k]][,j]
      }
      ## TODO? separate eXo variables?
      if(tmpnmvs > nX){
        for(j in 1:(tmpnmvs - nX)){
          y[lavdata@case.idx[[k]],(j+nX)] <- lavdata@eXo[[k]][,j]
        }
      }
      g[lavdata@case.idx[[k]]] <- k
    }
    colnames(y) <- orig.ov.names

    ## stan data block
    datablk <- paste0(datablk, t1, "vector[", tmpnmvs,
                      "] y[N];\n")
    
    ## remove fully deleted rows
    nas <- which(apply(is.na(y), 1, sum) == tmpnmvs)

    if(length(nas) > 0){
      y <- y[-nas,]
      g <- g[-nas]
      ntot <- sum(unlist(lavdata@nobs))
    }

    if(FALSE){ #blavmis == "fi"){
      ## variables on rhs of regression, in case missing
      y <- y[names(y) %in% ovreg]
    }

    standata <- list(y=y, g=g, N=ntot)
    ## TODO needed?
    if(any(partable$op == "|")){
      standata <- c(standata, list(ones = matrix(1, ntot, tmpnmvs)))
    }

    if(FALSE){ #blavmis == "fi"){
      ## keep only modeled y's not on rhs of regression
      matvars <- which(orig.ov.names %in% ov.names &
                       !(orig.ov.names %in% ovreg))
      ## NB y is now a matrix, so no need for ymat
      ymat <- ymat[,matvars]
      nmvs <- length(matvars)
  
      ydf <- data.frame(y=as.numeric(ymat),
                        g=rep(g, nmvs),
                        sub=rep(1:ntot, nmvs),
                        mv=rep(matvars, each=ntot))
      
      ydf <- subset(ydf, !is.na(ydf$y) & !(ydf$mv %in% ovcol))
      ## sub index excluding completely missing observations
      ydf$sub <- as.numeric(as.factor(ydf$sub))

      standata <- c(standata, list(yvec=ydf$y, sub=ydf$sub,
                                   mv=ydf$mv, nrows=nrow(ydf)))
      standata$g <- ydf$g
      standata$N <- max(ydf$sub)
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
        pmats[[i]] <- array(0, c(matrows[i], matcols[i], ngroups))
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
        tpdecs <- paste0(tpdecs, t1, "matrix[", tmpdim[1],
                         ",", tmpdim[2], "] thetld[", tmpdim[3],
                         "];\n")
      }
      tpeqs <- paste0(tpeqs, t1, tpnames[i], " = ",
                      names(pmats)[i], ";\n")
    }
    tpdecs <- paste0(tpdecs, t1, "real mu[N,", tmpnmvs, "];\n")

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
    
    ## add cholesky decomp of theta matrix
    TPS <- paste0(TPS, t1, "}\n\n")
    TPS <- paste0(TPS, t1, "for(j in 1:", ngroups, "){\n")
    TPS <- paste0(TPS, t2, "thetld[j] = fill_lower(to_matrix(theta[,,j]));\n")
    TPS <- paste0(TPS, t2, "thetld[j] = cholesky_decompose(",
                  "thetld[j]);\n", t1, "}\n")
    
    TPS <- paste0("transformed parameters{\n", tpdecs, "\n",
                  tpeqs, TXT2, "\n\n", TPS,
                  "\n}\n\n")
    datablk <- paste0(datablk, datdecs, "}\n\n")

    standata <- c(standata, pmats)
        
    out <- c(out, list(data=standata))
  }

  funblk <- paste0("functions{\n", t1, "#include 'sem_lv.stan' \n",
                   t1, "#include 'fill_lower.stan' \n")
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

      lavpartable$free[idx] <- 0
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
       sd = sdvec)
}
