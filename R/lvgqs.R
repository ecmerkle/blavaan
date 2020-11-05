lvgqs <- function(modmats, standata) {
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
  eta <- with(standata, matrix(NA, Ntot, w9use + w9no))
  
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

  if ((w9use + w9no) > 0) {
    for (mm in 1:Np){
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

      corner <- cov_eta %*% L_Yt
      bottom_right <- cov_eta

      ## FIXME?? what if obsidx also extends to x variables?
      obsidx <- Obsvar[mm, ]
      precision[1:Nobs[mm],1:Nobs[mm]] <- solve(top_left[obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]], drop=FALSE])
      L <- bottom_right[usepsi,usepsi,drop=FALSE] - (corner[,obsidx[1:Nobs[mm]],drop=FALSE] %*% precision[1:Nobs[mm],1:Nobs[mm]] %*% t(corner[,obsidx[1:Nobs[mm]],drop=FALSE]))[usepsi,usepsi,drop=FALSE]
      if (all(round(L, 6) == 0)) {
        L <- matrix(0, nrow=NROW(L), ncol=NCOL(L))
      } else {
        L <- chol(L)
      }
      beta <- corner[, obsidx[1:Nobs[mm]], drop=FALSE] %*% precision[1:Nobs[mm], 1:Nobs[mm], drop=FALSE]

      r1 <- startrow[mm]
      r2 <- endrow[mm]

      for (idx in r1:r2){
        lvmean <- modmats[[grpidx]]$alpha + beta[, 1:Nobs[mm], drop=FALSE] %*% (YX[idx, 1:Nobs[mm]] - ovmean[[grpidx]][obsidx[1:Nobs[mm]]])
        eta[idx,usepsi] <- t(rmnorm(1, lvmean[usepsi], sqrt = L))
        if (w9no > 0) {
          eta[idx,nopsi] <- eta[idx,usepsi,drop=FALSE] %*% t(A[nopsi,usepsi,drop=FALSE])
        }
      }
    }
  }
  
  eta
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
