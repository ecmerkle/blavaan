lvgqs <- function(modmats, standata, getlvs = TRUE) {
  Ng <- length(modmats)

  ## stan data
  p <- standata$p
  q <- standata$q
  m <- standata$m
  Np <- standata$Np
  grpnum <- standata$grpnum
  Obsvar <- standata$Obsvar
  Nobs <- standata$Nobs
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
  ovmean <- L_Y_A <- vector("list", Ng)
  for (g in 1:Ng){
    ovmean[[g]] <- modmats[[g]]$nu

    if (p > 0){
      L_Y_A[[g]] <- modmats[[g]]$lambda %*% solve(I - modmats[[g]]$beta)

      if (standata$m > 0){
        ovmean[[g]][1:p] <- ovmean[[g]][1:p,] + L_Y_A[[g]] %*% modmats[[g]]$alpha[1:m,]
      }
    }
  }

  if ((w9use + w9no) > 0 | !getlvs) {
    if (!getlvs) {
      ## find row/col of missing observations that we are monitoring
      for (mm in 1:Np) {
        obsidx <- Obsvar[mm, ]
        r1 <- startrow[mm]
        r2 <- endrow[mm]
        YXimp[r1:r2, obsidx[1:Nobs[mm]]] <- YX[r1:r2, 1:Nobs[mm]]
      }
      misvals <- is.na(YXimp)
    }
    
    for (mm in 1:Np) {
      grpidx <- grpnum[mm]
      
      A <- solve(I - modmats[[grpidx]]$beta)
      total_eta_eta <- A - I
      indirect_eta_eta <- total_eta_eta - modmats[[grpidx]]$beta
      total_eta_y <- modmats[[grpidx]]$lambda %*% A
      indirect_eta_y <- total_eta_y - modmats[[grpidx]]$lambda

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
          lvmean <- modmats[[grpidx]]$alpha + beta[, 1:Nobs[mm], drop=FALSE] %*% (YX[idx, 1:Nobs[mm]] - ovmean[[grpidx]][obsidx[1:Nobs[mm]]])
          eta[idx,usepsi] <- t(rmnorm(1, lvmean[usepsi], sqrt = L))
          if (w9no > 0) {
            eta[idx,nopsi] <- eta[idx,usepsi,drop=FALSE] %*% t(A[nopsi,usepsi,drop=FALSE])
          }
        } else if (anymis) {
          ovreg <- ovmean[[grpidx]][misidx] + beta[, 1:Nobs[mm], drop=FALSE] %*% (YX[idx, 1:Nobs[mm]] - ovmean[[grpidx]][obsidx[1:Nobs[mm]]])
          YXimp[idx, obsidx[1:Nobs[mm]]] <- YX[idx, 1:Nobs[mm]]
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

samp_lvs <- function(mcobj, lavmodel, lavpartable, standata, thin = 1) {
  lavmcmc <- make_mcmc(mcobj)
  itnums <- sampnums(mcobj, thin = thin)
  nsamps <- length(itnums)
  nchain <- length(lavmcmc)

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks

  loop.args <- list(X = 1:nsamps, future.seed = TRUE, FUN = function(i){
      tmpmat <- array(NA, dim=c(nchain, standata$Ntot, standata$w9use + standata$w9no))
      for(j in 1:nchain){
        lavmodel <- fill_params(lavmcmc[[j]][i,], lavmodel, lavpartable)

        modmats <- vector("list", nblocks)
        for (g in 1:nblocks) {
          ## which mm belong to group g?
          mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]

          modmats[[g]] <- lavmodel@GLIST[ mm.in.group ]
          if(!("beta" %in% names(modmats[[g]]))) modmats[[g]]$beta <- matrix(0, standata$m, standata$m)
        }

        tmpmat[j,,] <- lvgqs(modmats, standata)
      }
      tmpmat})

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
          if(!("beta" %in% names(modmats[[g]]))) modmats[[g]]$beta <- matrix(0, standata$m, standata$m)
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

  blavaan:::samp_lvs(fit@external$mcmcout, fit@Model, fit@ParTable, stantrans$data)
}
