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
  parts <- rstan::extract_sparse_parts(!(skeleton==0L)) #is.na(skeleton))
  parts$w <- as.array(vals[!(vals==0L)]) #is.na(vals)])
  parts$v <- as.array(parts$v)
  parts$u <- as.array(parts$u)
  return(parts)
}

# Get prior parameters in a manner that Stan will like
#
# @param skeleton A matrix that indicates the restrictions placed on
#   its elements. If `NA`, then the element is unrestricted. If
#   `Inf` or `-Inf`, the element is unrestricted but is constrained
#   to be positive or negative respectively. Otherwise, the element is 
#   fixed to the specified number, which is often zero but can be any finite 
#   value.
# @param Ng Number of groups
# @param wskel Matrix of equality constraints
# @param param1 first parameter of prior distribution (possibly supplied by user)
# @param param2 second parameter of prior distribution (possibly supplied by user)
# @return A list containing the prior parameters
format_priors <- function(skeleton, Ng, wskel, param1, param2) {
  nfree <- Ng * sum(!is.finite(skeleton)) - sum(wskel[,1] == 1L)
  stopifnot(length(param1) == 1L || length(param1) == nfree)
  if (length(param1) == 1L) param1 <- array(param1, nfree)
  else param1 <- as.numeric(param1)

  stopifnot(length(param2) == 1L || length(param2) == nfree)
  if (length(param2) == 1L) param2 <- array(param2, nfree)
  else param2 <- as.numeric(param2)

  return(list(p1=param1, p2=param2))
}


#' Bayesian Structural Equation Models via Stan
#'
#' Obtain data list for LERSIL().
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
#' 
#' @examples
#' data("Bollen", package = "sem")
#' # fill in
stanmarg_data <- function(YX = NULL, S = NULL, N, Ng, grpnum, # data
                          miss, Np, Nobs, Obsvar, # missing
                          startrow, endrow, save_lvs = FALSE, 
                          Lambda_y_skeleton, # skeleton matrices
                          Lambda_x_skeleton, Gamma_skeleton, B_skeleton,
                          Theta_skeleton, Theta_r_skeleton,
                          Theta_x_skeleton, Theta_x_r_skeleton,
                          Psi_skeleton, Psi_r_skeleton, Phi_skeleton,
                          Phi_r_skeleton, Nu_skeleton, Alpha_skeleton,
                          w1skel, w2skel, w3skel, # eq constraint matrices
                          w4skel, w5skel, w6skel, w7skel, w8skel,
                          w9skel, w10skel, w11skel, w12skel, w13skel,
                          w14skel,
                          lam_y_sign, lam_x_sign, # sign constraint matrices
                          gam_sign, b_sign, psi_r_sign, phi_r_sign,
                          lambda_y_mn = 0, lambda_y_sd = 10, # prior settings
                          lambda_x_mn = 0, lambda_x_sd = 10,
                          gamma_mn = 0, gamma_sd = 10,
                          b_mn = 0, b_sd = 10, theta_r_alpha = 1,
                          theta_r_beta = 1, theta_x_r_alpha = 1,
                          theta_x_r_beta = 1, psi_r_alpha = 1, psi_r_beta = 1,
                          phi_r_alpha = 1, phi_r_beta = 1, nu_mn = 0,
                          nu_sd = 50, alpha_mn = 0, alpha_sd = 10,
                          theta_sd_rate = .5, theta_x_sd_rate = .5,
                          psi_sd_rate = .5, phi_sd_rate = .5,
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
  dat$startrow <- startrow
  dat$endrow <- endrow

  dat$YX <- YX
  stopifnot(nrow(dat$YX) == dat$Ntot)
  
  if (missing(YX)) {
    dat$has_data <- 0L
    if (is.null(S)) stop("S must be specified if YX is missing")
    if (!is.list(S)) stop("S must be a list")
    if (is.null(N)) stop("N must be specified if YX is missing")
    if (miss) stop("blavaan ERROR: missingness requires raw data.")
    dat$S <- array(NA, dim=c(Ng, nrow(S), nrow(S)))
    for(i in 1:Ng) {
      dat$S[i,,] <- (N[i] - 1) * S[[i]]
    }
    dat$YX <- array(NA_real_, dim = c(dat$Ntot, ncol(S)))
  } else {
    dat$has_data <- 1L

    if (NROW(YX) != dat$Ntot) stop("blavaan ERROR: nrow(YX) != Ntot.")
    
    dat$S <- array(NA, dim=c(Ng, NCOL(YX), NCOL(YX)))
    for(i in 1:Ng) {
      dat$S[i,,] <- (dat$N[i] - 1) * cov(YX[(startrow[i] : endrow[i]),]) # NB!! this multiplication is needed to use
                                        # wishart_lpdf
    }
  }

  if(dat$missing & save_lvs) stop("blavaan ERROR: lvs cannot currently be saved when data are missing.")
  dat$save_lvs <- save_lvs
  
  parts <- make_sparse_skeleton(Lambda_y_skeleton)
  dat$p <- nrow(Lambda_y_skeleton)
  dat$m <- ncol(Lambda_y_skeleton)
  dat$len_w1 <- length(parts$w)
  dat$w1 <- parts$w
  dat$v1 <- parts$v
  dat$u1 <- parts$u
  vals <- Lambda_y_skeleton[!is.finite(Lambda_y_skeleton)]
  dat$small_w1 <- as.array(which(is.na(vals)))
  dat$len_small_w1 <- length(dat$small_w1)
  dat$w1skel <- w1skel
  dat$lam_y_sign <- lam_y_sign
  
  parts <- make_sparse_skeleton(Lambda_x_skeleton)
  dat$q <- nrow(Lambda_x_skeleton)
  dat$n <- ncol(Lambda_x_skeleton)
  dat$len_w2 <- length(parts$w)
  dat$w2 <- parts$w
  dat$v2 <- parts$v
  dat$u2 <- parts$u
  vals <- Lambda_x_skeleton[!is.finite(Lambda_x_skeleton)]
  dat$small_w2 <- as.array(which(is.na(vals)))
  dat$len_small_w2 <- length(dat$small_w2)
  dat$w2skel <- w2skel
  dat$lam_x_sign <- lam_x_sign
  
  parts <- make_sparse_skeleton(Gamma_skeleton)
  dat$len_w3 <- length(parts$w)
  dat$w3 <- parts$w
  dat$v3 <- parts$v
  dat$u3 <- parts$u
  vals <- Gamma_skeleton[!is.finite(Gamma_skeleton)]
  dat$small_w3 <- as.array(which(is.na(vals)))
  dat$len_small_w3 <- length(dat$small_w3)
  dat$w3skel <- w3skel
  dat$gam_sign <- gam_sign

  parts <- make_sparse_skeleton(B_skeleton)
  dat$len_w4 <- length(parts$w)
  dat$w4 <- parts$w
  dat$v4 <- parts$v
  dat$u4 <- parts$u
  vals <- B_skeleton[!is.finite(B_skeleton)]
  dat$small_w4 <- as.array(which(is.na(vals)))
  dat$len_small_w4 <- length(dat$small_w4)
  dat$w4skel <- w4skel
  dat$b_sign <- b_sign

  dThet <- Theta_skeleton
  dThet[lower.tri(dThet)] <- dThet[upper.tri(dThet)] <- 0L
  parts <- make_sparse_skeleton(dThet)
  dat$len_w5 <- length(parts$w)
  dat$w5 <- parts$w
  dat$v5 <- parts$v
  dat$u5 <- parts$u
  vals <- dThet[!is.finite(dThet)]
  dat$small_w5 <- as.array(which(is.na(vals)))
  dat$len_small_w5 <- length(dat$small_w5)
  dat$w5skel <- w5skel

  dThetx <- Theta_x_skeleton
  dThetx[lower.tri(dThetx)] <- dThetx[upper.tri(dThetx)] <- 0L
  parts <- make_sparse_skeleton(dThetx)
  dat$len_w6 <- length(parts$w)
  dat$w6 <- parts$w
  dat$v6 <- parts$v
  dat$u6 <- parts$u
  vals <- dThetx[!is.finite(dThetx)]
  dat$small_w6 <- as.array(which(is.na(vals)))
  dat$len_small_w6 <- length(dat$small_w6)
  dat$w6skel <- w6skel

  parts <- make_sparse_skeleton(Theta_r_skeleton)
  dat$len_w7 <- length(parts$w)
  dat$w7 <- parts$w
  dat$v7 <- parts$v
  dat$u7 <- parts$u  
  vals <- Theta_r_skeleton[!is.finite(Theta_r_skeleton)]
  dat$small_w7 <- as.array(which(is.na(vals)))
  dat$len_small_w7 <- length(dat$small_w7)
  dat$w7skel <- w7skel

  parts <- make_sparse_skeleton(Theta_x_r_skeleton)
  dat$len_w8 <- length(parts$w)
  dat$w8 <- parts$w
  dat$v8 <- parts$v
  dat$u8 <- parts$u
  vals <- Theta_x_r_skeleton[!is.finite(Theta_x_r_skeleton)]
  dat$small_w8 <- as.array(which(is.na(vals)))
  dat$len_small_w8 <- length(dat$small_w8)
  dat$w8skel <- w8skel
  
  dPsi <- Psi_skeleton
  dPsi[lower.tri(dPsi)] <- dPsi[upper.tri(dPsi)] <- 0L
  parts <- make_sparse_skeleton(dPsi)
  dat$len_w9 <- length(parts$w)
  dat$w9 <- parts$w
  dat$v9 <- parts$v
  dat$u9 <- parts$u
  vals <- dPsi[!is.finite(dPsi)]
  dat$small_w9 <- as.array(which(is.na(vals)))
  dat$len_small_w9 <- length(dat$small_w9)
  dat$w9skel <- w9skel

  parts <- make_sparse_skeleton(Psi_r_skeleton)
  dat$len_w10 <- length(parts$w)
  dat$w10 <- parts$w
  dat$v10 <- parts$v
  dat$u10 <- parts$u
  vals <- Psi_r_skeleton[!is.finite(Psi_r_skeleton)]
  dat$small_w10 <- as.array(which(is.na(vals)))
  dat$len_small_w10 <- length(dat$small_w10)
  dat$w10skel <- w10skel
  dat$psi_r_sign <- psi_r_sign

  dPhi <- Phi_skeleton
  dPhi[lower.tri(dPhi)] <- dPhi[upper.tri(dPhi)] <- 0L
  parts <- make_sparse_skeleton(dPhi)
  dat$len_w11 <- length(parts$w)
  dat$w11 <- parts$w
  dat$v11 <- parts$v
  dat$u11 <- parts$u
  vals <- dPhi[!is.finite(dPhi)]
  dat$small_w11 <- as.array(which(is.na(vals)))
  dat$len_small_w11 <- length(dat$small_w11)
  dat$w11skel <- w11skel

  parts <- make_sparse_skeleton(Phi_r_skeleton)
  dat$len_w12 <- length(parts$w)
  dat$w12 <- parts$w
  dat$v12 <- parts$v
  dat$u12 <- parts$u
  vals <- Phi_r_skeleton[!is.finite(Phi_r_skeleton)]
  dat$small_w12 <- as.array(which(is.na(vals)))
  dat$len_small_w12 <- length(dat$small_w12)
  dat$w12skel <- w12skel
  dat$phi_r_sign <- phi_r_sign

  if(dat$has_data & is.null(Nu_skeleton)) stop("blavaan ERROR: Nu_skeleton not provided")
  parts <- make_sparse_skeleton(Nu_skeleton)
  dat$len_w13 <- length(parts$w)
  dat$w13 <- parts$w
  dat$v13 <- parts$v
  dat$u13 <- parts$u
  vals <- Nu_skeleton[!is.finite(Nu_skeleton)]
  dat$small_w13 <- as.array(which(is.na(vals)))
  dat$len_small_w13 <- length(dat$small_w13)
  dat$w13skel <- w13skel

  parts <- make_sparse_skeleton(Alpha_skeleton)
  dat$len_w14 <- length(parts$w)
  dat$w14 <- parts$w
  dat$v14 <- parts$v
  dat$u14 <- parts$u
  vals <- Alpha_skeleton[!is.finite(Alpha_skeleton)]
  dat$small_w14 <- as.array(which(is.na(vals)))
  dat$len_small_w14 <- length(dat$small_w14)
  dat$w14skel <- w14skel

  ## priors
  pris <- format_priors(Lambda_y_skeleton, dat$Ng, dat$w1skel, lambda_y_mn, lambda_y_sd)
  dat$lambda_y_mn <- pris[['p1']]; dat$lambda_y_sd <- pris[['p2']]
  dat$len_lam_y <- length(dat$lambda_y_mn)
  
  pris <- format_priors(Lambda_x_skeleton, dat$Ng, dat$w2skel, lambda_x_mn, lambda_x_sd)
  dat$lambda_x_mn <- pris[['p1']]; dat$lambda_x_sd <- pris[['p2']]
  dat$len_lam_x <- length(dat$lambda_x_mn)

  pris <- format_priors(Gamma_skeleton, dat$Ng, dat$w3skel, gamma_mn, gamma_sd)
  dat$gamma_mn <- pris[['p1']]; dat$gamma_sd <- pris[['p2']]
  dat$len_gam <- length(dat$gamma_mn)
  
  pris <- format_priors(B_skeleton, dat$Ng, dat$w4skel, b_mn, b_sd)
  dat$b_mn <- pris[['p1']]; dat$b_sd <- pris[['p2']]
  dat$len_b <- length(dat$b_mn)

  pris <- format_priors(dThet, dat$Ng, dat$w5skel, theta_sd_rate, theta_sd_rate)
  dat$theta_sd_rate <- pris[['p1']] #; dat$b_sd <- pris[['p2']]
  dat$len_thet_sd <- length(dat$theta_sd_rate)

  pris <- format_priors(dThetx, dat$Ng, dat$w6skel, theta_x_sd_rate, theta_x_sd_rate)
  dat$theta_x_sd_rate <- pris[['p1']] #; dat$b_sd <- pris[['p2']]
  dat$len_thet_x_sd <- length(dat$theta_x_sd_rate)
  
  pris <- format_priors(Theta_r_skeleton, dat$Ng, dat$w7skel, theta_r_alpha, theta_r_beta)
  dat$theta_r_alpha <- pris[['p1']]; dat$theta_r_beta <- pris[['p2']]
  dat$len_thet_r <- length(dat$theta_r_alpha)
  
  pris <- format_priors(Theta_x_r_skeleton, dat$Ng, dat$w8skel, theta_x_r_alpha, theta_x_r_beta)
  dat$theta_x_r_alpha <- pris[['p1']]; dat$theta_x_r_beta <- pris[['p2']]
  dat$len_thet_x_r <- length(dat$theta_x_r_alpha)

  pris <- format_priors(dPsi, dat$Ng, dat$w9skel, psi_sd_rate, psi_sd_rate)
  dat$psi_sd_rate <- pris[['p1']] #; dat$b_sd <- pris[['p2']]
  dat$len_psi_sd <- length(dat$psi_sd_rate)
  
  pris <- format_priors(Psi_r_skeleton, dat$Ng, dat$w10skel, psi_r_alpha, psi_r_beta)
  dat$psi_r_alpha <- pris[['p1']]; dat$psi_r_beta <- pris[['p2']]
  dat$len_psi_r <- length(dat$psi_r_alpha)

  pris <- format_priors(dPhi, dat$Ng, dat$w11skel, phi_sd_rate, phi_sd_rate)
  dat$phi_sd_rate <- pris[['p1']] #; dat$b_sd <- pris[['p2']]
  dat$len_phi_sd <- length(dat$phi_sd_rate)

  pris <- format_priors(Phi_r_skeleton, dat$Ng, dat$w12skel, phi_r_alpha, phi_r_beta)
  dat$phi_r_alpha <- pris[['p1']]; dat$phi_r_beta <- pris[['p2']]
  dat$len_phi_r <- length(dat$phi_r_alpha)
  
  pris <- format_priors(Nu_skeleton, dat$Ng, dat$w13skel, nu_mn, nu_sd)
  dat$nu_mn <- pris[['p1']]; dat$nu_sd <- pris[['p2']]
  dat$len_nu <- length(dat$nu_mn)

  pris <- format_priors(Alpha_skeleton, dat$Ng, dat$w14skel, alpha_mn, alpha_sd)
  dat$alpha_mn <- pris[['p1']]; dat$alpha_sd <- pris[['p2']]
  dat$len_alph <- length(dat$alpha_mn)

  stopifnot(length(theta_sd_rate) == 1L || length(theta_sd_rate) == (dat$Ng * dat$p - sum(dat$w5skel[,1] == 1L)))
  if (length(theta_sd_rate) == 1L) dat$theta_sd_rate <- rep(theta_sd_rate, dat$Ng * dat$p - sum(dat$w5skel[,1] == 1L))
  else dat$theta_sd_rate <- as.numeric(theta_sd_rate)

  stopifnot(length(theta_x_sd_rate) == 1L   || length(theta_x_sd_rate) == (dat$Ng * dat$q - sum(dat$w6skel[,1] == 1L)))
  if (length(theta_x_sd_rate) == 1L) dat$theta_x_sd_rate <- rep(theta_x_sd_rate, dat$Ng * dat$q - sum(dat$w6skel[,1] == 1L))
  else dat$theta_x_sd_rate <- as.numeric(theta_x_sd_rate)
  
  return(dat)
}
