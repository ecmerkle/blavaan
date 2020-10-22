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

  I <- diag(m)
  precision <- matrix(0, (p+q), (p+q))
  
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
    precision[1:Nobs[mm],1:Nobs[mm]] <- solve(top_left[obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]])
    ## FIXME? is chol the same as cholesky_decompose?
    browser()
    L <- chol(bottom_right[usepsi,usepsi] - (corner[,obsidx[1:Nobs[mm]]] %*% precision[1:Nobs[mm],1:Nobs[mm]] %*% t(corner[,obsidx[1:Nobs[mm]]]))[usepsi,usepsi])
    beta[, 1:Nobs[mm]] <- corner[, obsidx[1:Nobs[mm]]] %*% precision[1:Nobs[mm], 1:Nobs[mm]]

    r1 <- startrow[mm]
    r2 <- endrow[mm]

    for (idx in r1:r2){
      lvmean <- modmats[[grpidx]]$alpha + beta[, 1:Nobs[mm]] %*% (YX[idx, 1:Nobs[mm]] - ovmean[grpidx, obsidx[1:Nobs[mm]]])
      eta[idx,usepsi] <- t(rmnorm(1, lvmean[usepsi], sqrt = L));
      if (w9no > 0) {
        eta[idx,nopsi] <- eta[idx,usepsi] %*% t(A[nopsi,usepsi]);
      }
    }
  }
  
  eta
}

matsetup <- function(object, standata, thin = 1) {
  lavmcmc <- blavInspect(object, 'mcmc')
  itnums <- sampnums(blavInspect(object, 'mcobj'), thin = thin)

  ## TODO change indexing:
  lavmodel <- fill_params(lavmcmc[[1]][5,], object@Model, object@ParTable)

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks

  modmats <- vector("list", nblocks)
  for (g in 1:nblocks) {
    ## which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]

    modmats[[g]] <- lavmodel@GLIST[ mm.in.group ]
  }
browser()
  eta <- lvgqs(modmats, standata)
  
  modmats
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

  blavaan:::matsetup(fit, stantrans$data)
}
