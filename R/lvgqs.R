lvgqs <- function() {
  for (g in 1:N){
    ovmean[[g]] <- Nu[[g]]

    if (p > 0){
      L_Y_A[[g]] <- L_Y[[g]] %*% solve(I - Bet[[g]])

      if (m > 0){
        ovmean[[g]][1:p] <- L_Y_A[[g]] %*% Alpha[[g]][1:m]
      }
    }
  }

  for (mm in 1:Np){
    grpidx <- grpnum[mm]

    A <- solve(I - Bet[[grpidx]])
    total_eta_eta <- A - I
    indirect_eta_eta <- total_eta_eta - Bet[[grpidx]]
    total_eta_y <- L_Y[[grpidx]] %*% A
    indirect_eta_y <- total_eta_y - L_Y[[grpidx]]

    Psi_star <- A %*% PS[[grpidx]] %*% t(A)
    Pi_t <- t(total_xi_eta)
    L_Yt <- t(L_Y[[grpidx]])

    cov_eta <- Psi_star
    top_left <- L_Y[[grpidx]] %*% cov_eta %*% L_Yt + Theta[[grpidx]]

    corner <- cov_eta %*% L_Yt
    bottom_right <- cov_eta

    ## FIXME?? what if obsidx also extends to x variables?
    obsidx <- Obsvar[mm, ]
    precision[1:Nobs[mm],1:Nobs[mm]] <- solve(top_left[obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]])
    ## FIXME? is chol the same as cholesky_decompose?
    L <- chol(bottom_right[usepsi,usepsi] - (corner[,obsidx[1:Nobs[mm]]] %*% precision[1:Nobs[mm],1:Nobs[mm]] %*% transpose(corner[,obsidx[1:Nobs[mm]]]))[usepsi,usepsi])
    beta[, 1:Nobs[mm]] <- corner[, obsidx[1:Nobs[mm]]] %*% precision[1:Nobs[mm], 1:Nobs[mm]]

    r1 <- startrow[mm]
    r2 <- endrow[mm]

    for (idx in r1:r2){
      lvmean <- Alpha[[grpidx]] + beta[, 1:Nobs[mm]] %*% (YX[idx, 1:Nobs[mm]] - ovmean[grpidx, obsidx[1:Nobs[mm]]])
      eta[idx,usepsi] <- t(rmnorm(1, lvmean[usepsi], sqrt = L));
      if (w9no > 0) {
        eta[idx,nopsi] <- eta[idx,usepsi] %*% t(A[nopsi,usepsi]);
      }
    }
  }
}
