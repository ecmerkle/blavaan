lvgqs <- function(modmats, standata, eeta = NULL, getlvs = TRUE) {
  Ng <- length(modmats)

  ## stan data
  ## FIXME for getlvs=FALSE, YX only contains continuous data so dimensions
  ## are wrong below when ordinal data are involved
  p <- standata$p
  q <- standata$q
  m <- standata$m
  Np <- standata$Np
  grpnum <- standata$grpnum
  Obsvar <- standata$Obsvar
  Nobs <- standata$Nobs
  Ntot <- standata$Ntot
  YX <- standata$YX
  usepsi <- standata$usepsi
  nopsi <- standata$nopsi
  startrow <- standata$startrow
  endrow <- standata$endrow
  w9no <- standata$w9no
  w9use <- standata$w9use

  I <- diag(m)
  precision <- matrix(0, (p+q), (p+q))
  if (getlvs) {
    eta <- with(standata, matrix(NA, Ntot, w9use + w9no))
  } else {
    YXimp <- matrix(NA, NROW(YX), NCOL(YX))
  }

  ## new matrices
  ovmean <- lvimpmean <- vector("list", Ng)
  for (g in 1:Ng) {
    if ("nu" %in% names(modmats[[g]])){
      ovmean[[g]] <- modmats[[g]]$nu
    } else {
      ovmean[[g]] <- t(standata$YXbar[g, , drop = FALSE])
    }

    if(!("beta" %in% names(modmats[[g]]))) modmats[[g]]$beta <- matrix(0, standata$m, standata$m)
    if (!("alpha" %in% names(modmats[[g]]))) modmats[[g]]$alpha <- matrix(0, m, 1)
    if (p > 0){
      lvimpmean[[g]] <-  solve(I - modmats[[g]]$beta) %*% modmats[[g]]$alpha[1:m,]

      if (standata$m > 0 && "alpha" %in% names(modmats[[g]])){
        ovmean[[g]][1:p] <- ovmean[[g]][1:p,] + modmats[[g]]$lambda %*% lvimpmean[[g]]
      }
    }
  }

  if (is.null(eeta)) {
    eeta <- vector("list", Ng)
    for (g in 1:Ng) {
      eeta[[g]] <- rep(0, standata$m)
    }
  }

  if ((w9use + w9no) > 0 | !getlvs) {
    if (!getlvs) {
      ## find row/col of missing observations that we are monitoring
      for (mm in 1:Np) {
        obsidx <- Obsvar[mm, ]
        r1 <- startrow[mm]
        r2 <- endrow[mm]
        YXimp[r1:r2, obsidx[1:Nobs[mm]]] <- YX[r1:r2, obsidx[1:Nobs[mm]]]
      }
      misvals <- is.na(YXimp)
    }
    
    for (mm in 1:Np) {
      grpidx <- grpnum[mm]
      
      A <- solve(I - modmats[[grpidx]]$beta)
      #total_eta_eta <- A - I
      #indirect_eta_eta <- total_eta_eta - modmats[[grpidx]]$beta
      #total_eta_y <- modmats[[grpidx]]$lambda %*% A
      #indirect_eta_y <- total_eta_y - modmats[[grpidx]]$lambda

      Psi_star <- A %*% modmats[[grpidx]]$psi %*% t(A)
      L_Yt <- t(modmats[[grpidx]]$lambda)

      cov_eta <- Psi_star
      top_left <- modmats[[grpidx]]$lambda %*% cov_eta %*% L_Yt + modmats[[grpidx]]$theta

      ## FIXME?? what if obsidx also extends to x variables?
      obsidx <- Obsvar[mm, ]
      misidx <- (1:NROW(top_left))[-obsidx[1:Nobs[mm]]]
      anymis <- Nobs[mm] < NROW(top_left)
      
      precision[1:Nobs[mm],1:Nobs[mm]] <- solve(top_left[obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]], drop=FALSE])

      if (getlvs) {
        corner <- cov_eta %*% L_Yt
        bottom_right <- cov_eta

        ## NB SEM-specific expressions for this matrix also exist
        L <- bottom_right[usepsi,usepsi,drop=FALSE] - (corner[,obsidx[1:Nobs[mm]],drop=FALSE] %*% precision[1:Nobs[mm],1:Nobs[mm]] %*% t(corner[,obsidx[1:Nobs[mm]],drop=FALSE]))[usepsi,usepsi,drop=FALSE]
        L <- try(chol(L))
        if (inherits(L, 'try-error')) {
          ## occasionally negative variance
          L <- matrix(0, nrow=length(usepsi), ncol=length(usepsi))
        }
        if (anymis) {
          corner <- corner[,obsidx[1:Nobs[mm]]]
        }
      } else if (anymis) {
        ## impute missing observed values
        corner <- top_left[misidx, obsidx[1:Nobs[mm]], drop=FALSE]
        bottom_right <- top_left[misidx, misidx, drop=FALSE]
        L <- chol(bottom_right - (corner %*% precision[1:Nobs[mm],1:Nobs[mm],drop=FALSE] %*% t(corner)))
      }

      if (getlvs | anymis) {
        beta <- corner %*% precision[1:Nobs[mm], 1:Nobs[mm], drop=FALSE]
      }

      r1 <- startrow[mm]
      r2 <- endrow[mm]
      for (idx in r1:r2){
        if (getlvs) {
          lvmean <- modmats[[grpidx]]$alpha + beta[, 1:Nobs[mm], drop=FALSE] %*% (YX[idx, obsidx[1:Nobs[mm]]] - ovmean[[grpidx]][obsidx[1:Nobs[mm]]])
          eta[idx,usepsi] <- t(rmnorm(1, lvmean[usepsi], sqrt = L) + eeta[[grpidx]][usepsi])
          if (w9no > 0) {
            eta[idx,nopsi] <- eta[idx,usepsi,drop=FALSE] %*% t(A[nopsi,usepsi,drop=FALSE])
          }
        } else if (anymis) {
          ovreg <- ovmean[[grpidx]][misidx] + beta[, 1:Nobs[mm], drop=FALSE] %*% (YX[idx, obsidx[1:Nobs[mm]]] - ovmean[[grpidx]][obsidx[1:Nobs[mm]]])
          YXimp[idx, obsidx[1:Nobs[mm]]] <- YX[idx, obsidx[1:Nobs[mm]]]
          YXimp[idx, misidx] <- t(rmnorm(1, ovreg, sqrt = L))
        }
      }
    }
  }

  if (getlvs) {
    out <- eta
  } else {
    out <- YXimp[misvals]
  }
  out
}

samp_lvs <- function(mcobj, lavmodel, lavpartable, standata, eeta, categorical, thin = 1) {
  lavmcmc <- make_mcmc(mcobj)
  itnums <- sampnums(mcobj, thin = thin)
  nsamps <- length(itnums)
  nchain <- length(lavmcmc)

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks

  loop.args <- list(X = 1:nsamps, FUN = function(i){
      tmpmat <- array(NA, dim=c(nchain, standata$Ntot, standata$w9use + standata$w9no))
      for(j in 1:nchain){
        lavmodel <- fill_params(lavmcmc[[j]][i,], lavmodel, lavpartable)

        modmats <- vector("list", nblocks)
        for(g in 1:nblocks) {
          ## which mm belong to group g?
          mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]

          modmats[[g]] <- lavmodel@GLIST[ mm.in.group ]
        }

        standata2 <- standata
        if(categorical) {
          ## use standata and lavmcmc[[j]][i,] to create YX matrix that has both
          ## continuous and ordinal
          tmpsamp <- lavmcmc[[j]][i,]
          YXo <- matrix(tmpsamp[grep("YXostar", names(tmpsamp))], standata$Ntot, standata$Nord,
                        byrow = TRUE)
          YXstar <- matrix(0, standata$Ntot, standata$p + standata$q)

          YXstar[, standata$ordidx] <- YXo

          if(with(standata, p + q - Nord) > 0){
            YXstar[, standata$contidx] <- standata$YX
          }

          standata2$YX <- YXstar
        }
        
        tmpmat[j,,] <- lvgqs(modmats, standata2, eeta)
      }
      tmpmat}, future.seed = TRUE)

  etasamps <- do.call("future_lapply", loop.args)
  etasamps <- array(unlist(etasamps), with(standata, c(nchain, Ntot, w9use + w9no, nsamps)))
  etasamps <- aperm(etasamps, c(4,1,3,2))
  dim(etasamps) <- with(standata, c(nsamps, nchain, Ntot * (w9use + w9no)))
  if((standata$w9use + standata$w9no) > 0){
    dimnames(etasamps)[[3]] <- with(standata, paste0("eta[", rep(1:Ntot, each=(w9use + w9no)), ",",
                                                     rep(1:(w9use + w9no), Ntot), "]"))
  }
  
  etasamps
}

samp_data <- function(mcobj, lavmodel, lavpartable, standata, lavdata, thin = 1) {
  lavmcmc <- make_mcmc(mcobj)
  itnums <- sampnums(mcobj, thin = thin)
  nsamps <- length(itnums)
  nchain <- length(lavmcmc)

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks

  loop.args <- list(X = 1:nsamps, future.seed = TRUE, FUN = function(i){
      tmplist <- vector("list", nchain)
      for(j in 1:nchain){
        lavmodel <- fill_params(lavmcmc[[j]][i,], lavmodel, lavpartable)

        modmats <- vector("list", nblocks)
        for (g in 1:nblocks) {
          ## which mm belong to group g?
          mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]

          modmats[[g]] <- lavmodel@GLIST[ mm.in.group ]
        }

        tmplist[[j]] <- lvgqs(modmats, standata, getlvs = FALSE)
      }
      tmplist})

  missamps <- do.call("future_lapply", loop.args)
  missamps <- array(unlist(missamps), with(standata, c(NROW(missamps[[1]][[1]]), nchain, nsamps)))
  missamps <- aperm(missamps, c(3,2,1))

  ## reorder to correspond to original data, vs to missingness pattern
  mp <- lavdata@Mp
  idx <- matrix(NA, dim(missamps)[3], 2)
  vnm <- rep(NA, dim(missamps)[3])

  srow <- 0L
  for (i in 1:length(mp)){
    for (j in 1:mp[[i]]$npatterns){
      misvars <- which(!mp[[i]]$pat[j,])
      cidx <- mp[[i]]$case.idx[[j]]
      nr <- length(misvars) * length(cidx)

      if(nr > 0){
        idx[(srow + 1):(srow + nr),1] <- rep(cidx, each = length(misvars))
        idx[(srow + 1):(srow + nr),2] <- rep(misvars, length(cidx))
        vnm[(srow + 1):(srow + nr)] <- rep(lavdata@ov.names[[i]][misvars], length(cidx))
      }
      srow <- srow + nr
    }
  }

  vnm <- vnm[order(idx[,1])]
  missamps[,,order(idx[,1])] <- missamps
  idx <- idx[order(idx[,1]),]
  dim(missamps) <- c(prod(dim(missamps)[1:2]), dim(missamps)[3])
  dimnames(missamps)[[2]] <- paste0(vnm, "[", idx[,1], "]")

  missamps
}

samp_lvs_2lev <- function(mcobj, lavmodel, lavsamplestats, lavdata, lavpartable, standata, eeta, thin = 1, debug = FALSE) {
  lavmcmc <- make_mcmc(mcobj)
  itnums <- sampnums(mcobj, thin = thin)
  nsamps <- length(itnums)
  nchain <- length(lavmcmc)

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  stanorig <- standata

  lav_implied22l <- getFromNamespace("lav_mvnorm_cluster_implied22l", "lavaan")
  lav_estep <- getFromNamespace("lav_mvnorm_cluster_em_estep_ranef", "lavaan")
  
  loop.args <- list(X = 1:nsamps, FUN = function(i){
      tmpmat <- array(NA, dim=c(nchain, standata$Ntot, standata$w9use + standata$w9no))
      tmpmat2 <- array(NA, dim=c(nchain, sum(standata$nclus[,2]), standata$w9use_c + standata$w9no_c))
      for(j in 1:nchain){
        lavmodel <- fill_params(lavmcmc[[j]][i,], lavmodel, lavpartable)

        ## get model-implied matrices from this lavmodel
        modmats <- vector("list", length = lavmodel@nblocks)
        nmat <- lavmodel@nmat
        for(b in seq_len(lavmodel@nblocks)) {
            mm.in.group <- 1:nmat[b] + cumsum(c(0,nmat))[b]

            modmats[[b]] <- lavmodel@GLIST[mm.in.group]
        }
        modmat2 <- modmats[2 * (1:standata$Ng)]
        clusmns <- vector("list", length(modmat2))
        modimp <- lav_model_implied(lavmodel) ## for all groups

        for(g in 1:length(modmat2)){
          if(!("beta" %in% names(modmat2[[g]]))) modmat2[[g]]$beta <- matrix(0, standata$m_c, standata$m_c)

          out <- lav_implied22l(lavdata@Lp[[g]], lapply(modimp, function(x) x[(2*g - 1):(2*g)]))
          clusmns[[g]] <- lav_estep(YLp = lavsamplestats@YLp[[g]], Lp = lavdata@Lp[[g]],
                                    sigma.w = out$sigma.w, sigma.b = out$sigma.b,
                                    sigma.zz = out$sigma.zz, sigma.yz = out$sigma.yz,
                                    mu.z = out$mu.z, mu.w = out$mu.w, mu.b = out$mu.b, se = FALSE)
        }
        clusmns <- do.call("rbind", clusmns)
        YX.B <- matrix(0, nrow = nrow(clusmns), ncol = ncol(stanorig$YX))
        Lp <- lavdata@Lp[[1]]
        YX.B[, Lp$ov.idx[[1]]] <- clusmns
        between.idx <- Lp$between.idx[[2]]

        if(length(between.idx) > 0L){
          YX.B[, between.idx] <- stanorig$YX[!duplicated(Lp$cluster.idx[[2]]), between.idx]
        }

        ## manipulations to reuse existing lvgqs code
        standata$p <- standata$p_c
        standata$q <- 0
        standata$m <- standata$m_c
        standata$usepsi <- standata$usepsi_c
        standata$nopsi <- standata$nopsi_c
        standata$w9use <- standata$w9use_c
        standata$w9no <- standata$w9no_c
        standata$endrow <- cumsum(standata$nclus[,2])
        standata$startrow <- c(1, standata$endrow[-length(standata$endrow)] + 1)
        standata$YX <- YX.B[, lavdata@Lp[[1]]$ov.idx[[2]], drop = FALSE]
        standata$Ntot <- sum(standata$nclus[,2])
        standata$Nobs <- with(standata, rep(N_between + N_both, Np))
        standata$Obsvar <- with(standata, matrix(1:standata$Nobs[1], Np, N_between + N_both, byrow = TRUE))
        tmpmat2[j,,] <- lvgqs(modmat2, standata, eeta[2*(1:standata$Ng)])

        ## now level 1
        standata <- stanorig
        clusidx <- rep(1:length(standata$cluster_size), standata$cluster_size)
        standata$YX <- with(standata, YX[, between_idx[(N_between + 1):p_tilde]]) - clusmns[clusidx,]
        modmat1 <- modmats[2 * (1:standata$Ng) - 2 + 1]
        for(g in 1:length(modmat1)){
          if(!("beta" %in% names(modmat1[[g]]))) modmat1[[g]]$beta <- matrix(0, standata$m, standata$m)
        }

        tmpmat[j,,] <- lvgqs(modmat1, standata, eeta[2 * (1:standata$Ng) - 1])
      }
      list(tmpmat, tmpmat2)})
  if(!debug) {
    loop.args <- c(loop.args, future.seed = TRUE)
    funcall <- "future_lapply"
  } else {
    funcall <- "lapply"
  }

  etasamps <- do.call(funcall, loop.args)

  etaout <- vector("list", 2)
  for (i in 1:2) {
    tmpeta <- lapply(etasamps, function(x) x[[i]])
    tmpN <- ifelse(i==1, standata$Ntot, sum(standata$nclus[,2]))
    tmpw9 <- ifelse(i==1, standata$w9use + standata$w9no, standata$w9use_c + standata$w9no_c)

    tmpeta <- array(unlist(tmpeta), with(standata, c(nchain, tmpN, tmpw9, nsamps)))

    tmpeta <- aperm(tmpeta, c(4,1,3,2))
    dim(tmpeta) <- with(standata, c(nsamps, nchain, tmpN * tmpw9))
    if(tmpw9 > 0) {
      dimnames(tmpeta)[[3]] <- with(standata, paste0("eta", ifelse(i==2, "_b", ""), "[",
                                                     rep(1:tmpN, each=tmpw9), ",",
                                                     rep(1:tmpw9, tmpN), "]"))
    }
    etaout[[i]] <- tmpeta
  }
  
  etaout
}

if(FALSE){
  model <- ' 
       # latent variable definitions
         ind60 =~ x1 + x2 + x3
         dem60 =~ y1 + a*y2 + b*y3 + c*y4
         dem65 =~ y5 + a*y6 + b*y7 + c*y8
     
       # regressions
         dem60 ~ ind60
         dem65 ~ ind60 + dem60
     
       # residual correlations
         y1 ~~ y5
         y2 ~~ y4 + y6
         y3 ~~ y7
         y4 ~~ y8
         y6 ~~ y8
     '

  fit <- bsem(model, data=PoliticalDemocracy, target=mytarg, burnin=100, sample=100,
              mcmcfile = TRUE)
  load('lavExport/semstan.rda')

  samp_lvs(fit@external$mcmcout, fit@Model, fit@ParTable, stantrans$data)
}
