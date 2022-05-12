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

    out[,i] <- get_ll(postsamp = samps[samprow,], fit,
                      casewise = TRUE, conditional = TRUE)
  }

  out <- exp(out) * do.call("rbind", w.star.list)
  
  return( log( rowSums(out) ) )
}

## fixed gauss-hermite quadrature, to reuse quadrature points across cases
fixed_ghq <- function(fit, ngq, samprow = NULL) {
  GLIST <- fit@Model@GLIST
  if (any(GLIST$theta[lower.tri(GLIST$theta)] != 0L)) stop("blavaan ERROR: The quadrature method cannot be used with non-diagonal theta matrix.")
  ndim <- NROW(GLIST$alpha)

  if (blavInspect(fit, 'ngroups') > 1) stop("blavaan ERROR: The quadrature method currently does not support multiple groups.")
  
  samps <- do.call("rbind", make_mcmc(fit@external$mcmcout, fit@external$stanlvs))
  if(length(samprow) > 0) samps <- samps[samprow, , drop = FALSE]

  lavigh <- getFromNamespace("lav_integration_gauss_hermite", "lavaan")
  XW <- lavigh(n = ngq, ndim = ndim, dnorm = TRUE)
  x.star <- XW$x
  x.star.eval <- apply(XW$x, 2, unique)
  w.star <- XW$w

  ## response patterns
  standata <- fit@external$mcmcdata
  if(length(standata$YX) > 0) stop("blavaan ERROR: The fixed quadrature method cannot handle mixes of continuous variables yet.")
  YX <- matrix(NA, NROW(standata$YX), NCOL(standata$YX) + NCOL(standata$YXo))
  YX[, standata$contidx] <- standata$YX
  YX[, standata$ordidx] <- standata$YXo
  rpatts <- apply(standata$YXo, 1, paste0, collapse = "")
  upatts <- as.numeric(as.factor(rpatts))
  ulocs <- which(!duplicated(upatts))
  ## FIXME: also need to consider Ng > 1 in response patterns:
  YXou <- standata$YXo[!duplicated(upatts), , drop = FALSE]
  deltas <- which(names(fit@Model@GLIST) == "delta")
  th.idx <- fit@Model@th.idx
  Ng <- blavInspect(fit, 'ngroups')
  TH.idx <- lapply(1:Ng, function(g) th.idx[[g]][th.idx[[g]] > 0])

  origlm <- fit@Model
  out <- matrix(NA, NROW(samps), NROW(YX))
  
  for(i in 1:NROW(samps)) {
    lavmodel <- fill_params(samps[i, , drop = FALSE], origlm, fit@ParTable)
    lavmodel@GLIST[[deltas]] <- NULL
    ## fit@Model <- lavmodel
    ## mnvec <- lavPredict(fit, type = "ov", newdata = fakedat,
    ##                     ETA = x.star.eval)
    ## if(inherits(mnvec, "matrix")) mnvec <- list(mnvec)

    ## for each entry in mnvec, compute univariate likelihoods for each set of thresholds
    ## a matrix per column of mnvec: number of rows in x.star.eval by number of ordered categories
    likevals <- array(NA, dim = c(NROW(x.star.eval), max(standata$YXo), NCOL(standata$YXo), Ng))

    for(g in 1:Ng) {
      mm.in.group <- 1:lavmodel@nmat[g] + cumsum(c(0,lavmodel@nmat[g]))[g]
      mms <- lavmodel@GLIST[mm.in.group]
      mnvec <- mms$lambda %*% t(x.star.eval)
      mnvec <- sweep(mnvec, 1, mms$nu, FUN = "+")

      for(j in 1:NCOL(standata$YXo)) {
        tmpidx <- unique(TH.idx[[g]])[j]
        tau <- c(-Inf, mms$tau[TH.idx[[g]] == tmpidx], Inf)
        utau <- tau[2:length(tau)]
        ltau <- tau[1:(length(tau) - 1)]

        for(k in 1:max(standata$YXo[,tmpidx])) {
          tmpprob <- pnorm(utau[k], mean = mnvec[tmpidx,], sd = sqrt(mms$theta[tmpidx, tmpidx])) -
            pnorm(ltau[k], mean = mnvec[tmpidx,], sd = sqrt(mms$theta[tmpidx, tmpidx]))

          #tmpprob[tmpprob == 0] <- 1e-300

          likevals[, k, j, g] <- log(tmpprob)
        }
      }
    }

    ## for each response pattern, use x.star to pull values out of the above matrices and sum
    qpt.uniq <- matrix(NA, NROW(YXou), NROW(x.star))
    diment <- apply(mms$lambda != 0, 1, which) ## FIXME only works for no cross-loadings
    tmpmatch <- sapply(1:ndim, function(j) match(x.star[,j], x.star.eval[,j]))

    ## all entries we need, which could replace the loop. but summing the right entries
    ## might take just as long.
    ##tmpent <- cbind(as.numeric(t(tmpmatch[,diment])), rep(as.numeric(t(YXou)), nrow(x.star)),
    ##                rep(1:9, nrow(x.star) * nrow(YXou)), 1)
    ## tmpeval <- likevals[tmpent]
    for(p in 1:NROW(x.star)) {
      tmpeval <- t(sapply(1:NROW(YXou), function(ii) likevals[cbind(tmpmatch[p,diment], YXou[ii,], 1:NCOL(YXou), 1)])) ## FIXME last index is currently fixed at 1 for group

      qpt.uniq[,p] <- rowSums(tmpeval)
    }

    qpt.uniq <- sweep(exp(qpt.uniq), 2, w.star, FUN = "*")

    ## FIXME deal with continuous data here
    
    ## assign values to full data matrix, for each response pattern
    full.lik <- rep(NA, NROW(YX))
    for(j in 1:length(ulocs)) {
      tmpidx <- which(rpatts == rpatts[ulocs[j]])
      full.lik[tmpidx] <- log(sum(qpt.uniq[j,]))
    }

    out[i,] <- full.lik
  }

  out
}

adapted_weights <- function(samps, ngq, alphas, psis, grpidx, etamns, etacovs, N) {
  ## adapt gh nodes/weights to each case
  ndim <- NROW(alphas[[1]])
  lavigh <- getFromNamespace("lav_integration_gauss_hermite", "lavaan")
  lavdmvnorm <- getFromNamespace("lav_mvnorm_dmvnorm", "lavaan")
  
  XW <- lavigh(n = ngq, ndim = ndim, dnorm = TRUE)
  eXWxcp <- exp(0.5 * apply(XW$x, 1, crossprod))

  x.star.list <- vector("list", length(etamns))
  w.star.list <- vector("list", length(etamns))
  XW2pi <- XW$w * (2*pi)^(ndim/2)

  for(i in 1:N) {
    C <- t(chol(etacovs[[i]]))
    tmpmn <- as.numeric(etamns[[i]])

    x.star.list[[i]] <- t(as.matrix(C %*% t(XW$x)) + tmpmn)
    w.star.list[[i]] <- XW2pi * eXWxcp * prod(diag(C)) * ## = det(C) for triangular matrix
      lavdmvnorm(x.star.list[[i]], Mu = alphas[[grpidx[i]]],
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
