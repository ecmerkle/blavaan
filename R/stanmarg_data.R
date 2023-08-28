# This code is based on code from the LERSIL package, authored by Ben Goodrich

# Convert a skeleton matrix in R to the internal format used by Stan
#
# @param skeleton A matrix that indicates the restrictions placed on
#   its elements. If `NA`, then the element is unrestricted. If
#   `Inf` or `-Inf`, the element is unrestricted but is constrained
#   to be positive or negative respectively. Otherwise, the element is 
#   fixed to the specified number, which is often zero but can be any finite 
#   value.
# @return A list containing the sparse representation of the input matrix
make_sparse_skeleton <- function(skeleton) {
  stopifnot(is.matrix(skeleton))
  #skeleton[is.na(skeleton)] <- 1L
  vals <- c(t(skeleton)) # vals needs to be in row-major order
  ## addresses change to Matrix clashing with extract_sparse_parts
  spmat <- Matrix::Matrix(!(skeleton==0L), doDiag=FALSE, sparse=TRUE)
  parts <- rstan::extract_sparse_parts(spmat) #is.na(skeleton))
  parts$w <- as.array(vals[!(vals==0L)]) #is.na(vals)])
  parts$v <- as.array(parts$v)
  parts$u <- as.array(parts$u)
  return(parts)
}

group_sparse_skeleton <- function(skeleton) {
  stopifnot(length(dim(skeleton)) == 3)

  Ng <- dim(skeleton)[1]

  g_len <- u_len <- array(NA, Ng)
  tmpu <- tmpw <- tmpv <- list()

  for (g in 1:Ng) {
    ## removing first dimension and maintaining dims 2+3 is tricky:
    parts <- make_sparse_skeleton(array(skeleton[g,,],
                                        dim = dim(skeleton)[2:3]))
    wlen <- length(parts$w)
    g_len[g] <- wlen
    u_len[g] <- length(parts$u)
    tmpu <- c(tmpu, list(parts$u))
    tmpw <- c(tmpw, list(parts$w))
    tmpv <- c(tmpv, list(parts$v))
  }

  vdat <- wdat <- matrix(0, Ng, max(g_len))
  udat <- matrix(0, Ng, max(u_len))

  for (g in 1:Ng) {
    if (g_len[g] > 0) {
      vdat[g, 1:g_len[g]] <- tmpv[[g]]
      wdat[g, 1:g_len[g]] <- tmpw[[g]]
    }
    if (u_len[g] > 0) {
      udat[g, 1:u_len[g]] <- tmpu[[g]]
    }
  }
  
  out <- list(g_len = g_len, v = vdat, w = wdat, u = udat)
  return(out)
}

# Get prior parameters in a manner that Stan will like
#
# @param lavpartable A lavaan partable with "priors" column
# @param mat The matrix for which we are obtaining priors
# @return A list containing the prior parameters
format_priors <- function(lavpartable, level = 1L) {
  ## parameter matrices are filled in by row, so need to make
  ## sure we get parameters in the right order!
  if ("group" %in% names(lavpartable)) {
    lavpartable <- lavpartable[order(lavpartable$group, lavpartable$col, lavpartable$row),]
  } else {
    lavpartable <- lavpartable[order(lavpartable$col, lavpartable$row),]
  }

  transtab <- list(c('lambda_y_mn', 'lambda_y_sd', 'len_lam_y'),
                   c('lambda_x_mn', 'lambda_x_sd', 'len_lam_x'),
                   c('gamma_mn', 'gamma_sd', 'len_gam'),
                   c('b_mn', 'b_sd', 'len_b'),
                   c('theta_sd_shape', 'theta_sd_rate', 'len_thet_sd', 'theta_pow'),
                   c('theta_x_sd_shape', 'theta_x_sd_rate', 'len_thet_x_sd', 'theta_x_pow'),
                   c('theta_r_alpha', 'theta_r_beta', 'len_thet_r'),
                   c('theta_x_r_alpha', 'theta_x_r_beta', 'len_thet_x_r'),
                   c('psi_sd_shape', 'psi_sd_rate', 'len_psi_sd', 'psi_pow'),
                   c('psi_r_alpha', 'psi_r_beta', 'len_psi_r'),
                   c('phi_sd_shape', 'phi_sd_rate', 'len_phi_sd', 'phi_pow'),
                   c('phi_r_alpha', 'phi_r_beta', 'len_phi_r'),
                   c('nu_mn', 'nu_sd', 'len_nu'),
                   c('alpha_mn', 'alpha_sd', 'len_alph'),
                   c('tau_mn', 'tau_sd', 'len_tau'))

  mats <- c('lambda', 'lambda_x', 'gamma', 'beta', 'thetavar', 'cov.xvar', 'thetaoff',
            'cov.xoff', 'psivar', 'psioff', 'phivar', 'phioff', 'nu', 'alpha', 'tau')
  if (level == 2L) {
    newmats <- c('lambda', 'beta', 'thetavar', 'thetaoff', 'psivar', 'psioff', 'nu', 'alpha')
    subloc <- match(newmats, mats)
    mats <- newmats
    transtab <- transtab[subloc]
    transtab <- lapply(transtab, function(x) paste0(x, '_c'))
  }

  out <- list()
  
  for (i in 1:length(mats)) {
    mat <- origmat <- mats[i]
    
    if (grepl("var", mat)) {
      mat <- gsub("var", "", mat)
      prisel <- lavpartable$row == lavpartable$col
    } else if (grepl("psioff", mat)) {
      mat <- "lvrho"
      prisel <- lavpartable$row != lavpartable$col
    } else if (grepl("thetaoff", mat)) {
      mat <- "rho"
      prisel <- lavpartable$row != lavpartable$col
    } else {
      prisel <- rep(TRUE, length(lavpartable$row))
    }

    if (mat == "nu") {
      prisel <- prisel & (lavpartable$mat %in% c(mat, "mean.x"))
    } else {
      prisel <- prisel & (lavpartable$mat == mat)
    }
  
    prisel <- prisel & (lavpartable$free > 0)
    thepris <- lavpartable$prior[prisel]

    if (length(thepris) > 0) {
      textpris <- thepris[thepris != ""]

      prisplit <- strsplit(textpris, "[, ()]+")
      
      param1 <- sapply(prisplit, function(x) x[2])

      if (!grepl("\\[", prisplit[[1]][3])) {
        param2 <- sapply(prisplit, function(x) x[3])
      } else {
        param2 <- rep(NA, length(param1))
      }

      ## check that var/sd/prec is same for all
      powargs <- sapply(prisplit, tail, 1)
      if (any(grepl("\\[", powargs))) {
        if (length(unique(powargs)) > 1) stop(paste0("blavaan ERROR: In matrix ", mat, ", all priors must be placed on either vars, sds, or precisions."))
      }
      powpar <- 1
      powarg <- powargs[1]
      if (grepl("\\[var\\]", powarg)) {
        powpar <- 2
      } else if (!grepl("\\[sd\\]", powarg)) {
        powpar <- -2
      }
    
      param1 <- array(as.numeric(param1), length(param1))
      param2 <- array(as.numeric(param2), length(param2))
    } else {
      param1 <- array(0, 0)
      param2 <- array(0, 0)
      powpar <- 1
    }

    out[[ transtab[[i]][1] ]] <- param1
    out[[ transtab[[i]][2] ]] <- param2
    out[[ transtab[[i]][3] ]] <- length(param1)

    if (origmat %in% c('thetavar', 'cov.xvar', 'psivar', 'phivar')) {
      out[[ transtab[[i]][4] ]] <- powpar
    }
  } # mats
  
  return(out)
}

# Check that priors match what is in the stan file
#
# @param lavpartable A lavaan partable with "priors" column
# @param mat The matrix for which we are obtaining priors
# @return nothing
check_priors <- function(lavpartable) {
  right_pris <- sapply(dpriors(target = "stan"), function(x) strsplit(x, "[, ()]+")[[1]][1])
  pt_pris <- sapply(lavpartable$prior[lavpartable$prior != ""], function(x) strsplit(x, "[, ()]+")[[1]][1])
  names(pt_pris) <- lavpartable$mat[lavpartable$prior != ""]

  primatch <- match(names(pt_pris), names(right_pris))

  badpris <- which(pt_pris != right_pris[primatch])

  if (length(badpris) > 0) {
    badtxt <- unique(paste(names(pt_pris)[badpris], pt_pris[badpris]))
    stop(paste0("blavaan ERROR: For target='stan', the following priors are not allowed:\n", paste(badtxt, collapse = "\n"), "\n See dpriors(target='stan') for each parameter's required distribution."))
  }
}


#' Obtain data list for stanmarg.
stanmarg_data <- function(YX = NULL, S = NULL, YXo = NULL, N, Ng, grpnum, # data
                          miss, Np, Nobs, Obsvar, # missing
                          ord, Nord, ordidx, contidx, nlevs,
                          Noent, Nordobs, OrdObsvar, # ordinal
                          Xvar, Xdatvar, Xbetvar, Nx, Nx_between, # fixed.x
                          startrow, endrow, save_lvs = FALSE, do_test = TRUE,
                          Lambda_y_skeleton, # skeleton matrices
                          Lambda_x_skeleton, Gamma_skeleton, B_skeleton,
                          Theta_skeleton, Theta_r_skeleton,
                          Theta_x_skeleton, Theta_x_r_skeleton,
                          Psi_skeleton, Psi_r_skeleton, Phi_skeleton,
                          Phi_r_skeleton, Nu_skeleton, Alpha_skeleton, Tau_skeleton,
                          w1skel, w2skel, w3skel, # eq constraint matrices
                          w4skel, w5skel, w6skel, w7skel, w8skel,
                          w9skel, w10skel, w11skel, w12skel, w13skel,
                          w14skel, w15skel, emiter,
                          lam_y_sign, lam_x_sign, # sign constraint matrices
                          gam_sign, b_sign, psi_r_sign, fullpsi, phi_r_sign,
                          lavpartable = NULL, # for priors
                          dumlv = NULL, # for sampling lvs
                          wigind = NULL, # wiggle indicator
                          pri_only = FALSE, # prior predictive sampling
                          do_reg = FALSE, # regression sampling
                          multilev, mean_d, cov_w, log_lik_x, cov_d, nclus, cluster_size, # level 2 data
                          ncluster_sizes, mean_d_full, cov_d_full, log_lik_x_full, xbar_w, xbar_b,
                          cov_b, gs, cluster_sizes, cluster_size_ns, between_idx, N_between,
                          within_idx, N_within, both_idx, N_both, ov_idx1, ov_idx2, p_tilde, N_lev,
                          Lambda_y_skeleton_c = NULL, # level 2 matrices
                          B_skeleton_c = NULL, Theta_skeleton_c = NULL, Theta_r_skeleton_c = NULL,
                          Psi_skeleton_c = NULL, Psi_r_skeleton_c = NULL, Nu_skeleton_c = NULL,
                          Alpha_skeleton_c = NULL,
                          w1skel_c = NULL, w4skel_c = NULL, w5skel_c = NULL, w7skel_c = NULL,
                          w9skel_c = NULL, w10skel_c = NULL, w13skel_c = NULL, w14skel_c = NULL,
                          lam_y_sign_c = NULL, b_sign_c = NULL, psi_r_sign_c = NULL,
                          phi_r_sign_c = NULL, fullpsi_c = NULL, dumlv_c = NULL, wigind_c = NULL,
                          ...) {
  
  dat <- list()

  dat$Ng <- Ng
  dat$N <- array(N, length(N))
  dat$Ntot <- sum(dat$N)
  dat$grpnum <- grpnum
  dat$missing <- miss
  dat$Np <- Np
  dat$Nobs <- array(Nobs, length(Nobs))
  dat$Obsvar <- Obsvar
  dat$ord <- ord
  dat$Nord <- Nord
  dat$Noent <- Noent
  dat$Nordobs <- Nordobs
  dat$OrdObsvar <- OrdObsvar
  dat$ordidx <- ordidx
  dat$contidx <- contidx
  dat$nlevs <- nlevs
  dat$startrow <- startrow
  dat$endrow <- endrow
  dat$Xvar <- Xvar
  dat$Xdatvar <- Xdatvar
  dat$Xbetvar <- Xbetvar
  dat$Nx <- array(Nx, length(Nx))
  dat$Nx_between <- array(Nx_between, length(Nx_between))
  dat$emiter <- emiter
  dat$do_reg <- do_reg
  dat$pri_only <- pri_only
  dat$multilev <- multilev
  
  dat$YX <- YX
  dat$YXo <- YXo

  dat$use_suff <- 1L
  if (ord | multilev) dat$use_suff <- 0L

  dat$has_data <- 0L
  dat$use_cov <- 0L
  if (pri_only) {
    dat$use_suff <- 0L
    if (dim(Nu_skeleton)[2] == 0L) dat$use_cov <- 1L
    tmparr <- array(dim = c(dat$Ng, ncol(YX) + 1, ncol(YX) + 1))
    for (i in 1:Ng) {
      tmparr[i,,] <- diag(nrow=ncol(YX) + 1)
    }
    dat$S <- tmparr
    dat$YXbar <- array(0, dim = c(dat$Ng, ncol(YX)))
  } else {
    if (missing(YX)) {
      dat$has_data <- 0L
      if (is.null(S)) stop("S must be specified if YX is missing")
      if (!is.list(S)) stop("S must be a list")
      if (is.null(N)) stop("N must be specified if YX is missing")
      if (miss) stop("blavaan ERROR: missingness requires raw data.")
      ## TODO if YXbar missing, set use_cov = 1 otherwise use_cov = 0 (but need new ppp)
      dat$use_cov <- 1L
      dat$S <- array(0, dim = c(dat$Ng, nrow(S[[1]]) + 1, nrow(S[[1]]) + 1))
      for (i in 1:Ng) {
        dat$S[i, 1:nrow(S[[1]]), 1:nrow(S[[1]])] <- S[[i]]
      }
      dat$YX <- array(0, dim = c(dat$Ntot, ncol(S[[1]])))
      dat$YXo <- array(0, dim = c(dat$Ntot, ncol(YXo)))
      dat$YXbar <- array(0, dim = c(dat$Ng, ncol(S[[1]])))
    } else {
      if (NROW(YX) != dat$Ntot) stop("blavaan ERROR: nrow(YX) != Ntot.")

      dat$YXbar <- array(0, dim=c(Np, NCOL(YX)))
      dat$S <- array(1, dim=c(Np, NCOL(YX) + 1, NCOL(YX) + 1))
      dat$has_data <- 1L
      if (dim(Nu_skeleton)[2] == 0L) dat$use_cov <- 1L
      if (length(contidx) > 0) {
        for (i in 1:dat$Np) {
          tmpyxbar <- colMeans(YX[(startrow[i] : endrow[i]), , drop = FALSE])
          tmpyxbar[is.na(tmpyxbar)] <- 0
          dat$YXbar[i,] <- tmpyxbar

          tmpN <- endrow[i] - startrow[i] + 1
          if(tmpN > 1) {
            tmpS <- cov(YX[(startrow[i] : endrow[i]), , drop = FALSE]) * (tmpN - 1) / tmpN
          } else {
            tmpS <- matrix(0, NCOL(YX), NCOL(YX))
          }
          dat$S[i,1:NCOL(YX),1:NCOL(YX)] <- tmpS
        }
      }
    }
  }

  dat$save_lvs <- save_lvs
  dat$do_test <- do_test

  ## level 2 data
  dat$mean_d <- mean_d
  dat$cov_w <- cov_w
  dat$log_lik_x <- log_lik_x
  dat$cov_d <- cov_d
  dat$mean_d_full <- mean_d_full
  dat$cov_d_full <- cov_d_full
  dat$log_lik_x_full <- log_lik_x_full
  dat$xbar_w <- xbar_w
  dat$xbar_b <- xbar_b
  dat$cov_b <- cov_b
  dat$gs <- gs
  dat$nclus <- nclus
  dat$cluster_size <- cluster_size
  dat$ncluster_sizes <- ncluster_sizes
  dat$cluster_sizes <- cluster_sizes
  dat$cluster_size_ns <- cluster_size_ns
  dat$between_idx <- between_idx; dat$N_between <- N_between
  dat$within_idx <- within_idx; dat$N_within <- N_within
  dat$both_idx <- both_idx; dat$N_both <- N_both
  dat$ov_idx1 <- ov_idx1; dat$ov_idx2 <- ov_idx2
  dat$p_tilde <- p_tilde
  dat$N_lev <- N_lev

  ## level 1 matrix info
  dat <- c(dat, stanmarg_matdata(dat, Lambda_y_skeleton, Lambda_x_skeleton,
                                 Gamma_skeleton, B_skeleton, Theta_skeleton,
                                 Theta_r_skeleton, Theta_x_skeleton, Theta_x_r_skeleton,
                                 Psi_skeleton, Psi_r_skeleton, Phi_skeleton, Phi_r_skeleton,
                                 Nu_skeleton, Alpha_skeleton, Tau_skeleton, dumlv, level = 1L))
  dat$lam_y_sign <- lam_y_sign
  dat$w1skel <- w1skel
  #dat$lam_x_sign <- lam_x_sign
  dat$w3skel <- w3skel
  dat$gam_sign <- gam_sign
  dat$w4skel <- w4skel
  dat$b_sign <- b_sign
  dat$w5skel <- w5skel
  dat$w6skel <- w6skel
  dat$w7skel <- w7skel
  dat$w8skel <- w8skel
  dat$w9skel <- w9skel
  dat$w10skel <- w10skel
  dat$psi_r_sign <- psi_r_sign
  dat$fullpsi <- fullpsi

  #dat$w11skel <- w11skel
  #dat$w12skel <- w12skel
  #dat$phi_r_sign <- phi_r_sign
  dat$w13skel <- w13skel
  dat$w14skel <- w14skel
  dat$w15skel <- w15skel

  ## level 2 matrices
  dat <- c(dat, stanmarg_matdata(dat, Lambda_y_skeleton = Lambda_y_skeleton_c,
                                 B_skeleton = B_skeleton_c, Theta_skeleton = Theta_skeleton_c,
                                 Theta_r_skeleton = Theta_r_skeleton_c,
                                 Psi_skeleton = Psi_skeleton_c, Psi_r_skeleton = Psi_r_skeleton_c,
                                 Nu_skeleton = Nu_skeleton_c, Alpha_skeleton = Alpha_skeleton_c,
                                 dumlv = dumlv_c, level = 2L))
  dat$lam_y_sign_c <- lam_y_sign_c
  dat$w1skel_c <- w1skel_c
  dat$w4skel_c <- w4skel_c
  dat$b_sign_c <- b_sign_c
  dat$w5skel_c <- w5skel_c
  dat$w7skel_c <- w7skel_c
  dat$w9skel_c <- w9skel_c
  dat$w10skel_c <- w10skel_c
  dat$psi_r_sign_c <- psi_r_sign_c
  dat$fullpsi_c <- fullpsi_c
  dat$w13skel_c <- w13skel_c
  dat$w14skel_c <- w14skel_c
  

  ## priors, first making sure they match what is in the stan file
  check_priors(lavpartable)
  dat$wigind <- wigind
  dat$wigind_c <- wigind_c

  if (!dat$multilev) {
    dat <- c(dat, format_priors(lavpartable))
    dat <- c(dat, format_priors(lavpartable[0,], level = 2L))
  } else {
    dat <- c(dat, format_priors(lavpartable[lavpartable$level == "within",]))
    dat <- c(dat, format_priors(lavpartable[lavpartable$level == "between",], level = 2L))
  }
  
  return(dat)
}


## create data on model matrices
stanmarg_matdata <- function(indat, Lambda_y_skeleton, Lambda_x_skeleton = NULL,
                             Gamma_skeleton = NULL, B_skeleton, Theta_skeleton,
                             Theta_r_skeleton, Theta_x_skeleton = NULL,
                             Theta_x_r_skeleton = NULL, Psi_skeleton,
                             Psi_r_skeleton, Phi_skeleton = NULL,
                             Phi_r_skeleton = NULL, Nu_skeleton, Alpha_skeleton,
                             Tau_skeleton = NULL, dumlv = NULL, level = 1L) {

  dat <- list()
  dat$p <- dim(Lambda_y_skeleton)[2]
  dat$m <- dim(Lambda_y_skeleton)[3]
  tmpres <- group_sparse_skeleton(Lambda_y_skeleton)
  dat$len_w1 <- max(tmpres$g_len)
  dat$w1 <- tmpres$w
  dat$v1 <- tmpres$v
  dat$u1 <- tmpres$u
  dat$wg1 <- array(tmpres$g_len, length(tmpres$g_len))

  if (level == 1L) {
    dat$q <- dim(Lambda_x_skeleton)[2]
    dat$n <- dim(Lambda_x_skeleton)[3]

    #tmpres <- group_sparse_skeleton(Lambda_x_skeleton)
    #dat$len_w2 <- max(tmpres$g_len)
    #dat$w2 <- tmpres$w
    #dat$v2 <- tmpres$v
    #dat$u2 <- tmpres$u
    #dat$wg2 <- array(tmpres$g_len, length(tmpres$g_len))
    #dat$w2skel <- w2skel


    tmpres <- group_sparse_skeleton(Gamma_skeleton)
    dat$len_w3 <- max(tmpres$g_len)
    dat$w3 <- tmpres$w
    dat$v3 <- tmpres$v
    dat$u3 <- tmpres$u
    dat$wg3 <- array(tmpres$g_len, length(tmpres$g_len))
  }
    
  tmpres <- group_sparse_skeleton(B_skeleton)
  dat$len_w4 <- max(tmpres$g_len)
  dat$w4 <- tmpres$w
  dat$v4 <- tmpres$v
  dat$u4 <- tmpres$u
  dat$wg4 <- array(tmpres$g_len, length(tmpres$g_len))

  dThet <- Theta_skeleton
  for (g in 1:indat$Ng) {
    tmpmat <- as.matrix(dThet[g,,])
    tmpmat[lower.tri(tmpmat)] <- tmpmat[upper.tri(tmpmat)] <- 0L
    dThet[g,,] <- tmpmat
  }
  tmpres <- group_sparse_skeleton(dThet)
  dat$len_w5 <- max(tmpres$g_len)
  dat$w5 <- sqrt(tmpres$w) # because we do SRS in the model
  dat$v5 <- tmpres$v
  dat$u5 <- tmpres$u
  dat$wg5 <- array(tmpres$g_len, length(tmpres$g_len))

  ## for lv sampling
  usethet <- array(which(diag(as.matrix(Theta_skeleton[1,,])) != 0))
  dat$w5use <- length(usethet)
  dat$usethet <- usethet

  if (level == 1L) {
    dThetx <- Theta_x_skeleton
    for (g in 1:indat$Ng) {
      tmpmat <- as.matrix(dThetx[g,,])
      tmpmat[lower.tri(tmpmat)] <- tmpmat[upper.tri(tmpmat)] <- 0L
      dThetx[g,,] <- tmpmat
    }
    tmpres <- group_sparse_skeleton(dThetx)
    dat$len_w6 <- max(tmpres$g_len)
    dat$w6 <- tmpres$w
    dat$v6 <- tmpres$v
    dat$u6 <- tmpres$u
    dat$wg6 <- array(tmpres$g_len, length(tmpres$g_len))
  }

  tmpres <- group_sparse_skeleton(Theta_r_skeleton)
  dat$len_w7 <- max(tmpres$g_len)
  dat$w7 <- tmpres$w
  dat$v7 <- tmpres$v
  dat$u7 <- tmpres$u
  dat$wg7 <- array(tmpres$g_len, length(tmpres$g_len))

  if (level == 1L) {
    tmpres <- group_sparse_skeleton(Theta_x_r_skeleton)
    dat$len_w8 <- max(tmpres$g_len)
    dat$w8 <- tmpres$w
    dat$v8 <- tmpres$v
    dat$u8 <- tmpres$u
    dat$wg8 <- array(tmpres$g_len, length(tmpres$g_len))
  }

  dPsi <- Psi_skeleton
  for (g in 1:indat$Ng) {
    tmpmat <- as.matrix(dPsi[g,,])
    tmpmat[lower.tri(tmpmat)] <- tmpmat[upper.tri(tmpmat)] <- 0L
    dPsi[g,,] <- tmpmat
  }
  tmpres <- group_sparse_skeleton(dPsi)
  dat$len_w9 <- max(tmpres$g_len)
  dat$w9 <- sqrt(tmpres$w) ## because we do SRS in the model
  dat$v9 <- tmpres$v
  dat$u9 <- tmpres$u
  dat$wg9 <- array(tmpres$g_len, length(tmpres$g_len))

  ## for lv sampling
  usepsi <- useorig <- array(which(diag(as.matrix(Psi_skeleton[1,,])) != 0))
  if (length(dumlv) > 0) {
    dums <- match(dumlv, usepsi)
    usepsi <- array(usepsi[-dums[!is.na(dums)]])
  }
  dat$w9use <- length(usepsi)
  dat$usepsi <- usepsi
  dat$nopsi <- array((1:dim(Psi_skeleton)[2])[-useorig])
  dat$w9no <- length(dat$nopsi)  

  tmpres <- group_sparse_skeleton(Psi_r_skeleton)
  dat$len_w10 <- max(tmpres$g_len)
  dat$w10 <- tmpres$w
  dat$v10 <- tmpres$v
  dat$u10 <- tmpres$u
  dat$wg10 <- array(tmpres$g_len, length(tmpres$g_len))

  if (level == 1L) {
    dPhi <- Phi_skeleton
    for (g in 1:indat$Ng) {
      tmpmat <- as.matrix(dPhi[g,,])
      tmpmat[lower.tri(tmpmat)] <- tmpmat[upper.tri(tmpmat)] <- 0L
      dPhi[g,,] <- tmpmat
    }
    #tmpres <- group_sparse_skeleton(dPhi)
    #dat$len_w11 <- max(tmpres$g_len)
    #dat$w11 <- tmpres$w
    #dat$v11 <- tmpres$v
    #dat$u11 <- tmpres$u
    #dat$wg11 <- array(tmpres$g_len, length(tmpres$g_len))

    #tmpres <- group_sparse_skeleton(Phi_r_skeleton)
    #dat$len_w12 <- max(tmpres$g_len)
    #dat$w12 <- tmpres$w
    #dat$v12 <- tmpres$v
    #dat$u12 <- tmpres$u
    #dat$wg12 <- array(tmpres$g_len, length(tmpres$g_len))
  }

  if (indat$has_data & is.null(Nu_skeleton)) stop("blavaan ERROR: Nu_skeleton not provided")
  tmpres <- group_sparse_skeleton(Nu_skeleton)
  dat$len_w13 <- max(tmpres$g_len)
  dat$w13 <- tmpres$w
  dat$v13 <- tmpres$v
  dat$u13 <- tmpres$u
  dat$wg13 <- array(tmpres$g_len, length(tmpres$g_len))

  tmpres <- group_sparse_skeleton(Alpha_skeleton)
  dat$len_w14 <- max(tmpres$g_len)
  dat$w14 <- tmpres$w
  dat$v14 <- tmpres$v
  dat$u14 <- tmpres$u
  dat$wg14 <- array(tmpres$g_len, length(tmpres$g_len))

  if (level == 1L) {
    tmpres <- group_sparse_skeleton(Tau_skeleton)
    dat$len_w15 <- max(tmpres$g_len)
    dat$w15 <- tmpres$w
    dat$v15 <- tmpres$v
    dat$u15 <- tmpres$u
    dat$wg15 <- array(tmpres$g_len, length(tmpres$g_len))
  }

  ## level 2 mats get "_c" suffix
  if (level == 2L) names(dat) <- paste0(names(dat), "_c")
  
  return(dat)
}
