## post-sampling, "adapted" gauss-hermite quadrature
adapted_ghq <- function(fit, ngq, samprow = NULL) {
  samps <- do.call("rbind", make_mcmc(fit@external$mcmcout, fit@external$stanlvs))

  lavmodel <- fill_params(samps[samprow, , drop = FALSE], fit@Model, fit@ParTable)
  GLIST <- lavmodel@GLIST
  if (any(GLIST$theta[lower.tri(GLIST$theta)] != 0L)) stop("blavaan ERROR: The quadrature method cannot be used with non-diagonal theta matrix.")
  alphas <- GLIST[which(names(GLIST) == "alpha")]
  psis <- GLIST[which(names(GLIST) == "psi")]
  N <- fit@SampleStats@ntotal

  grpidx <- rep(1, N)
  if(fit@Data@ngroups > 1) grpidx <- fit@Data@group

  ## compute mean and cov of each case i
  etamncov <- eta_moments(samps, N)
  
  ## get weights
  aws <- adapted_weights(samps, ngq, alphas, psis, grpidx, etamncov[[1]], etamncov[[2]], N)
  x.star.list <- aws[[1]]; w.star.list <- aws[[2]]

  ## loop thru quadrature points
  nqpt <- NROW(x.star.list[[1]])
  out <- matrix(NA, length(etamncov[[1]]), nqpt)
  
  for(i in 1:nqpt){
    samps[samprow,grep("^eta", colnames(samps))] <- as.numeric( sapply(1:length(etamncov[[1]]), function(k) x.star.list[[k]][i,]) )

    out[,i] <- blavaan:::get_ll(postsamp = samps[samprow,], fit,
                                casewise = TRUE, conditional = TRUE)
  }

  for(i in 1:nrow(out)){
    out[i,] <- exp(out[i,]) * w.star.list[[i]]
  }

  return( log( rowSums(out) ) )
}
  
adapted_weights <- function(samps, ngq, alphas, psis, grpidx, etamns, etacovs, N) {
  ## adapt gh nodes/weights to each case
  ndim <- NROW(alphas[[1]])
  XW <- lavaan:::lav_integration_gauss_hermite(n = ngq, ndim = ndim, dnorm = TRUE)
  eXWxcp <- exp(0.5 * apply(XW$x, 1, crossprod))

  x.star.list <- vector("list", length(etamns))
  w.star.list <- vector("list", length(etamns))
  XW2pi <- XW$w * (2*pi)^(ndim/2)
  
  for(i in 1:N) {
    C <- t(chol(etacovs[[j]]))
    tmpmn <- as.numeric(etamns[[i]])

    x.star.list[[i]] <- t(as.matrix(C %*% t(XW$x)) + tmpmn)
    w.star.list[[i]] <- XW2pi * eXWxcp * prod(diag(C)) * ## = det(C) for triangular matrix
      lavaan:::lav_mvnorm_dmvnorm(x.star.list[[i]], Mu = alphas[[grpidx[i]]],
                                  Sigma = psis[[grpidx[i]]], log = FALSE)
  }
  
  list(x.star.list, w.star.list)
}

eta_moments <- function(samps, N) {
  ## columns containing etas
  etasamps <- samps[, grep("^eta", colnames(samps))]

  etamns <- etacovs <- vector("list", N)
  
  for (i in 1:N) {
    tmpcol <- grep(paste0("eta[", i, ","), colnames(etasamps), fixed = TRUE)

    etamns[[i]] <- colMeans(etasamps[, tmpcol])

    etacovs[[i]] <- cov(etasamps[, tmpcol])
  }

  list(etamns, etacovs)
}
