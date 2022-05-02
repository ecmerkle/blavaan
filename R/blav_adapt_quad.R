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

## fixed gauss-hermite quadrature, to reuse quadrature points across cases
## FIXME this is unfinished!!
fixed_ghq <- function(fit, ngq, samprow = NULL) {
  GLIST <- fit@Model@GLIST
  if (any(GLIST$theta[lower.tri(GLIST$theta)] != 0L)) stop("blavaan ERROR: The quadrature method cannot be used with non-diagonal theta matrix.")
  ndim <- NROW(GLIST$alpha)
  
  samps <- do.call("rbind", make_mcmc(fit@external$mcmcout, fit@external$stanlvs))
  if(length(samprow) > 1) samps <- samps[samprow, , drop = FALSE]

  XW <- lavaan:::lav_integration_gauss_hermite(n = ngq, ndim = ndim, dnorm = TRUE)
  x.star <- XW$x
  x.star.eval <- apply(XW$x, 2, unique)
  w.star <- XW$w

  ## response patterns
  standata <- fit@external$mcmcdata
  YX <- matrix(NA, NROW(standata$YX), NCOL(standata$YX) + NCOL(standata$YXo))
  YX[, standata$contidx] <- standata$YX
  YX[, standata$ordidx] <- standata$YXo
  rpatts <- apply(standata$YXo, 1, paste0)
  upatts <- as.numeric(as.factor(rpatts))
  YXou <- standata$YXo[!duplicated(upatts), , drop = FALSE]
  deltas <- which(names(fit@Model@GLIST) == "delta")

  origlm <- fit@Model

  out <- matrix(NA, NROW(samps), NROW(YX))
  
  for(i in 1:NROW(samps)) {
    lavmodel <- fill_params(samps[i, , drop = FALSE], origlm, fit@ParTable)
    lavmodel@GLIST[[deltas]] <- NULL
    fit@Model <- lavmodel
    mnvec <- lavPredict(fit, type = "ov", ETA = x.star.eval)
    if(inherits(mnvec, "matrix")) mnvec <- list(mnvec)

    ## for each entry in mnvec, do line 345 of model_loglik for each set of thresholds
    ## a matrix per column of mnvec: number of rows in x.star.eval by number of ordered categories

    ## check for continuous data and throw error for now

    ## for each response pattern, use x.star to pull values out of the above matrices and sum
    qpt.uniq <- matrix(NA, NROW(YXou), NROW(x.star))



    qpt.uniq <- sweep(exp(qpt.uniq), 2, w.star, FUN = "*")
      
    ## assign values to full data matrix, for each response pattern

    out[i,] <- above
  }

  out
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
