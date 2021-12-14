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
format_priors <- function(lavpartable, mat) {
  ## parameter matrices are filled in by row, so need to make
  ## sure we get parameters in the right order!
  lavpartable <- lavpartable[order(lavpartable$group, lavpartable$col, lavpartable$row),]

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
    param1 <- array(0,0)
    param2 <- array(0,0)
    powpar <- 1
  }
  
  return(list(p1=param1, p2=param2, powpar=powpar))
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
#'
#' @export
#' @param S If `YX` is not supplied, then you must supply
#'   `S`, which is a covariance matrix among observed variables (divided by N-1)
#' @param N If `lavobject` is not supplied, then you must provide
#'   `N`, which is an integer indicating the number of observations
#' @param If `lavobject` is not supplied, then you must provide number of groups (for multi-group analysis)
#' @param Lambda_y_skeleton,Lambda_x_skeleton,Gamma_skeleton,B_skeleton
#'   Matrices indicating the restrictions placed on the elements using
#'   the parameterization of the LISREL software. If `NA`, then the
#'   element is unrestricted but presumably not too far from zero. If
#'   `Inf` or `-Inf`, the element is unrestricted but is constrained
#'   to be positive or negative respectively and is presumably far from
#'   zero. Otherwise, the element is fixed to the specified number, 
#'   which is often zero but can be any finite value.
#' @param theta_sd_rate,theta_x_sd_rate,psi_sd_rate,phi_sd_rate Vectors (possibly of length one)
#'   containing the rate parameter under independent exponential priors
#'   on the standard deviations of the measurement errors in `Y` and `X`
#'   respectively. If either of these are of length one, this value is
#'   recylcled to the number of columns in `Y` or `X` respectively
#' @param ... Further arguments
#' @return A list of data
#' @details Explain this
stanmarg_data <- function(YX = NULL, S = NULL, YXo = NULL, N, Ng, grpnum, # data
                          miss, Np, Nobs, Obsvar, # missing
                          ord, Nord, ordidx, contidx, nlevs, neach, # ordinal
                          Xvar, Xdatvar, Nx, # fixed.x
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
  dat$ordidx <- ordidx
  dat$contidx <- contidx
  dat$nlevs <- nlevs
  dat$neach <- neach
  dat$startrow <- startrow
  dat$endrow <- endrow
  dat$Xvar <- Xvar
  dat$Xdatvar <- Xdatvar
  dat$Nx <- array(Nx, length(Nx))
  dat$emiter <- emiter

  dat$YX <- YX
  dat$YXo <- YXo
  stopifnot(nrow(dat$YX) == dat$Ntot)

  dat$has_data <- dat$has_cov <- 0L
  if (pri_only) {
    tmparr <- array(dim = c(dat$Ng, ncol(YX) + 1, ncol(YX) + 1))
    for (i in 1:Ng) {
      tmparr[i,,] <- diag(nrow=ncol(YX) + 1)
    }
    dat$S <- tmparr
  } else {
    if (missing(YX)) {
      dat$has_data <- 0L
      if (is.null(S)) stop("S must be specified if YX is missing")
      if (!is.list(S)) stop("S must be a list")
      if (is.null(N)) stop("N must be specified if YX is missing")
      if (miss) stop("blavaan ERROR: missingness requires raw data.")
      dat$S <- array(NA, dim=c(Ng, nrow(S), nrow(S)))
      for (i in 1:Ng) {
        dat$S[i,,] <- (N[i] - 1) * S[[i]]
      }
      dat$YX <- array(NA_real_, dim = c(dat$Ntot, ncol(S)))
      dat$YXo <- array(NA_real_, dim = c(0, ncol(YXo)))
    } else {
      dat$has_data <- 1L

      if (NROW(YX) != dat$Ntot) stop("blavaan ERROR: nrow(YX) != Ntot.")
    
      dat$S <- array(NA, dim=c(Ng, NCOL(YX) + 1, NCOL(YX) + 1))
      for (i in 1:Ng) {
        dat$S[i,,] <- diag(1, NCOL(YX) + 1)
        ## not added because, if not pd, stan fails: (dat$N[i] - 1) * cov(YX[(startrow[i] : endrow[i]), , drop = FALSE]) # NB!! this multiplication is needed to use wishart_lpdf
      }
    }
  }
  dat$save_lvs <- save_lvs
  dat$do_test <- do_test
  
  dat$p <- dim(Lambda_y_skeleton)[2]
  dat$m <- dim(Lambda_y_skeleton)[3]
  tmpres <- group_sparse_skeleton(Lambda_y_skeleton)
  dat$len_w1 <- max(tmpres$g_len)
  dat$w1 <- tmpres$w
  dat$v1 <- tmpres$v
  dat$u1 <- tmpres$u
  dat$wg1 <- array(tmpres$g_len, length(tmpres$g_len))
  dat$w1skel <- w1skel
  dat$lam_y_sign <- lam_y_sign

  dat$q <- dim(Lambda_x_skeleton)[2]
  dat$n <- dim(Lambda_x_skeleton)[3]
  tmpres <- group_sparse_skeleton(Lambda_x_skeleton)
  #dat$len_w2 <- max(tmpres$g_len)
  #dat$w2 <- tmpres$w
  #dat$v2 <- tmpres$v
  #dat$u2 <- tmpres$u
  #dat$wg2 <- array(tmpres$g_len, length(tmpres$g_len))
  #dat$w2skel <- w2skel
  #dat$lam_x_sign <- lam_x_sign

  tmpres <- group_sparse_skeleton(Gamma_skeleton)
  dat$len_w3 <- max(tmpres$g_len)
  dat$w3 <- tmpres$w
  dat$v3 <- tmpres$v
  dat$u3 <- tmpres$u
  dat$wg3 <- array(tmpres$g_len, length(tmpres$g_len))
  dat$w3skel <- w3skel
  dat$gam_sign <- gam_sign

  tmpres <- group_sparse_skeleton(B_skeleton)
  dat$len_w4 <- max(tmpres$g_len)
  dat$w4 <- tmpres$w
  dat$v4 <- tmpres$v
  dat$u4 <- tmpres$u
  dat$wg4 <- array(tmpres$g_len, length(tmpres$g_len))
  dat$w4skel <- w4skel
  dat$b_sign <- b_sign

  dThet <- Theta_skeleton
  for (g in 1:Ng) {
    tmpmat <- as.matrix(dThet[g,,])
    tmpmat[lower.tri(tmpmat)] <- tmpmat[upper.tri(tmpmat)] <- 0L
    dThet[g,,] <- tmpmat
  }
  tmpres <- group_sparse_skeleton(dThet)
  dat$len_w5 <- max(tmpres$g_len)
  dat$w5 <- tmpres$w
  dat$v5 <- tmpres$v
  dat$u5 <- tmpres$u
  dat$wg5 <- array(tmpres$g_len, length(tmpres$g_len))
  dat$w5skel <- w5skel
  ## for lv sampling
  usethet <- array(which(diag(as.matrix(Theta_skeleton[1,,])) != 0))
  dat$w5use <- length(usethet)
  dat$usethet <- usethet
  
  dThetx <- Theta_x_skeleton
  for (g in 1:Ng) {
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
  dat$w6skel <- w6skel

  tmpres <- group_sparse_skeleton(Theta_r_skeleton)
  dat$len_w7 <- max(tmpres$g_len)
  dat$w7 <- tmpres$w
  dat$v7 <- tmpres$v
  dat$u7 <- tmpres$u
  dat$wg7 <- array(tmpres$g_len, length(tmpres$g_len))
  dat$w7skel <- w7skel

  tmpres <- group_sparse_skeleton(Theta_x_r_skeleton)
  dat$len_w8 <- max(tmpres$g_len)
  dat$w8 <- tmpres$w
  dat$v8 <- tmpres$v
  dat$u8 <- tmpres$u
  dat$wg8 <- array(tmpres$g_len, length(tmpres$g_len))
  dat$w8skel <- w8skel

  dPsi <- Psi_skeleton
  for (g in 1:Ng) {
    tmpmat <- as.matrix(dPsi[g,,])
    tmpmat[lower.tri(tmpmat)] <- tmpmat[upper.tri(tmpmat)] <- 0L
    dPsi[g,,] <- tmpmat
  }
  tmpres <- group_sparse_skeleton(dPsi)
  dat$len_w9 <- max(tmpres$g_len)
  dat$w9 <- tmpres$w
  dat$v9 <- tmpres$v
  dat$u9 <- tmpres$u
  dat$wg9 <- array(tmpres$g_len, length(tmpres$g_len))
  dat$w9skel <- w9skel
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
  dat$w10skel <- w10skel
  dat$psi_r_sign <- psi_r_sign
  dat$fullpsi <- fullpsi

  dPhi <- Phi_skeleton
  for (g in 1:Ng) {
    tmpmat <- as.matrix(dPhi[g,,])
    tmpmat[lower.tri(tmpmat)] <- tmpmat[upper.tri(tmpmat)] <- 0L
    dPhi[g,,] <- tmpmat
  }
  tmpres <- group_sparse_skeleton(dPhi)
  #dat$len_w11 <- max(tmpres$g_len)
  #dat$w11 <- tmpres$w
  #dat$v11 <- tmpres$v
  #dat$u11 <- tmpres$u
  #dat$wg11 <- array(tmpres$g_len, length(tmpres$g_len))
  #dat$w11skel <- w11skel

  tmpres <- group_sparse_skeleton(Phi_r_skeleton)
  #dat$len_w12 <- max(tmpres$g_len)
  #dat$w12 <- tmpres$w
  #dat$v12 <- tmpres$v
  #dat$u12 <- tmpres$u
  #dat$wg12 <- array(tmpres$g_len, length(tmpres$g_len))
  #dat$w12skel <- w12skel
  #dat$phi_r_sign <- phi_r_sign

  if(dat$has_data & is.null(Nu_skeleton)) stop("blavaan ERROR: Nu_skeleton not provided")
  tmpres <- group_sparse_skeleton(Nu_skeleton)
  dat$len_w13 <- max(tmpres$g_len)
  dat$w13 <- tmpres$w
  dat$v13 <- tmpres$v
  dat$u13 <- tmpres$u
  dat$wg13 <- array(tmpres$g_len, length(tmpres$g_len))
  dat$w13skel <- w13skel

  tmpres <- group_sparse_skeleton(Alpha_skeleton)
  dat$len_w14 <- max(tmpres$g_len)
  dat$w14 <- tmpres$w
  dat$v14 <- tmpres$v
  dat$u14 <- tmpres$u
  dat$wg14 <- array(tmpres$g_len, length(tmpres$g_len))
  dat$w14skel <- w14skel

  tmpres <- group_sparse_skeleton(Tau_skeleton)
  dat$len_w15 <- max(tmpres$g_len)
  dat$w15 <- tmpres$w
  dat$v15 <- tmpres$v
  dat$u15 <- tmpres$u
  dat$wg15 <- array(tmpres$g_len, length(tmpres$g_len))
  dat$w15skel <- w15skel
  
  ## priors; first make sure they match what is in the stan file
  check_priors(lavpartable)
  dat$wigind <- wigind
  
  pris <- format_priors(lavpartable, "lambda")
  dat$lambda_y_mn <- pris[['p1']]; dat$lambda_y_sd <- pris[['p2']]
  dat$len_lam_y <- length(dat$lambda_y_mn)
  
  pris <- format_priors(lavpartable, "lambda_x")
  dat$lambda_x_mn <- pris[['p1']]; dat$lambda_x_sd <- pris[['p2']]
  dat$len_lam_x <- length(dat$lambda_x_mn)

  pris <- format_priors(lavpartable, "gamma")
  dat$gamma_mn <- pris[['p1']]; dat$gamma_sd <- pris[['p2']]
  dat$len_gam <- length(dat$gamma_mn)
  
  pris <- format_priors(lavpartable, "beta")
  dat$b_mn <- pris[['p1']]; dat$b_sd <- pris[['p2']]
  dat$len_b <- length(dat$b_mn)

  pris <- format_priors(lavpartable, "thetavar")
  dat$theta_sd_shape <- pris[['p1']]
  dat$theta_sd_rate <- pris[['p2']]
  dat$len_thet_sd <- length(dat$theta_sd_rate)
  dat$theta_pow <- pris[['powpar']]

  pris <- format_priors(lavpartable, "cov.xvar")
  dat$theta_x_sd_shape <- pris[['p1']]
  dat$theta_x_sd_rate <- pris[['p2']]
  dat$len_thet_x_sd <- length(dat$theta_x_sd_rate)
  dat$theta_x_pow <- pris[['powpar']]
  
  pris <- format_priors(lavpartable, "thetaoff")
  dat$theta_r_alpha <- pris[['p1']]; dat$theta_r_beta <- pris[['p2']]
  dat$len_thet_r <- length(dat$theta_r_alpha)
  
  pris <- format_priors(lavpartable, "cov.xoff")
  dat$theta_x_r_alpha <- pris[['p1']]; dat$theta_x_r_beta <- pris[['p2']]
  dat$len_thet_x_r <- length(dat$theta_x_r_alpha)

  pris <- format_priors(lavpartable, "psivar")
  dat$psi_sd_shape <- pris[['p1']]
  dat$psi_sd_rate <- pris[['p2']]
  dat$len_psi_sd <- length(dat$psi_sd_rate)
  dat$psi_pow <- pris[['powpar']]

  pris <- format_priors(lavpartable, "psioff")
  dat$psi_r_alpha <- pris[['p1']]; dat$psi_r_beta <- pris[['p2']]
  dat$len_psi_r <- length(dat$psi_r_alpha)

  pris <- format_priors(lavpartable, "phivar")
  dat$phi_sd_shape <- pris[['p1']]
  dat$phi_sd_rate <- pris[['p2']]
  dat$len_phi_sd <- length(dat$phi_sd_rate)
  dat$phi_pow <- pris[['powpar']]
  
  pris <- format_priors(lavpartable, "phioff")
  dat$phi_r_alpha <- pris[['p1']]; dat$phi_r_beta <- pris[['p2']]
  dat$len_phi_r <- length(dat$phi_r_alpha)
  
  pris <- format_priors(lavpartable, "nu")
  dat$nu_mn <- pris[['p1']]; dat$nu_sd <- pris[['p2']]
  dat$len_nu <- length(dat$nu_mn)

  pris <- format_priors(lavpartable, "alpha")
  dat$alpha_mn <- pris[['p1']]; dat$alpha_sd <- pris[['p2']]
  dat$len_alph <- length(dat$alpha_mn)

  pris <- format_priors(lavpartable, "tau")
  dat$tau_mn <- pris[['p1']]; dat$tau_sd <- pris[['p2']]
  dat$len_tau <- length(dat$tau_mn)
  
  return(dat)
}
