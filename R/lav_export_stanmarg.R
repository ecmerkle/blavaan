matattr <- function(free, est, constraint, mat, Ng, std.lv, ...) {
  ## obtain matrix attributes, equality constraint,
  ## sign constraint information.
  ## free: a specific list of matrices from lavInspect(, 'free')
  ## est: a specific list of matrices from lavInspect(, 'est')
  ## constraint: data.frame of parameter constraints
  ## mat: name of matrix
  ## Ng: number of groups
  ## ...: additional arguments (free2 and sign matrix for
  ##                            associated lambda matrix)

  ddd <- list(...)

  matskel <- est[[1]]
  matskel[free[[1]] != 0L] <- Inf

  ## for corrmats, only use lower triangle
  if (grepl("_r", mat)) {
    matskel[upper.tri(matskel)] <- 0    
    for (i in 1:Ng) {
      free[[i]][upper.tri(free[[i]])] <- 0
    }
  }
  
  ## re-number free parameters just for this matrix
  start <- 1L
  free2 <- free
  for (i in 1:Ng) {
    tm <- free[[i]]
    np <- sum(tm != 0)
    tm[tm != 0] <- start:(start + np - 1)
    free2[[i]] <- tm
    
    start <- start + np
  }

  ## only done here to ensure matrices are right size for
  ## stanmarg_data; we use the results of make_sparse_skeleton
  ## in stanmarg_data instead of here so that stanmarg_data does not
  ## depend on lavaan objects.
  parts <- make_sparse_skeleton(matskel)
  len <- length(parts$w)

  wskel <- matrix(0, Ng * len, 2)
  if (NROW(constraint) > 0) {
    freemat <- do.call("rbind", free)
    free2mat <- do.call("rbind", free2)
    for (i in 1:NROW(constraint)) {
      usecon <- c(any(freemat == constraint$lhs[i]),
                  any(freemat == constraint$rhs[i]))
      if (sum(usecon) == 2) {
        lhsnum <- free2mat[freemat == constraint$lhs[i]]
        rhsnum <- free2mat[freemat == constraint$rhs[i]]
        constraint$rhs2 <- rhsnum ## for sign constraints below

        wskel[rhsnum, 1:2] <- c(1L, lhsnum)
      } else if (sum(usecon) == 1) {
        stop("blavaan ERROR: cross-matrix equality constraints not supported.")
      }
    }

    ## wskel[,2] must refer to only free parameters, skipping over
    ## constrained parameters
    freepars <- cumsum(wskel[,1] == 0)
    wskel[wskel[,1]==1,2] <- freepars[wskel[wskel[,1]==1,2]]
  }
  
  lvmat <- mat %in% c('Gamma', 'B', 'Psi_r')
  lammat <- grepl('Lambda', mat)
  sign <- matrix(0, Ng * len, 2 + lvmat)
  if (std.lv & (lvmat | lammat)) {
    if (lvmat) {
      lamfree <- ddd$free2
      lamsign <- ddd$sign
      
      for (i in 1:length(free2)) {
        fpar <- which(free2[[i]] != 0, arr.ind = TRUE)
        if (nrow(fpar) > 0) {
          for (j in 1:nrow(fpar)) {
            ## find sign-constrained loadings of the two lvs
            lampar1 <- lamfree[[i]][,fpar[j,2]]
            l1 <- min(lampar1[lampar1 != 0L])
            if (lamsign[l1,1] == 1) l1 <- lamsign[l1,2]

            lampar2 <- lamfree[[i]][,fpar[j,1]]
            l2 <- min(lampar2[lampar2 != 0])
            if (lamsign[l2,1] == 1) l2 <- lamsign[l2,2]

            rowloc <- free2[[i]][fpar[j,1], fpar[j,2]]
            sign[rowloc, 2:3] <- c(l1, l2)
          }
        }
      }
    } else {
      for (i in 1:length(free2)) {
        for (j in 1:NCOL(free2[[i]])) {
          col <- free2[[i]][,j]
          parnums <- col[col != 0L]
          psign <- min(parnums)
          ## if equality constraint, sign must involve the
          ## "free" parameter
          if (wskel[psign,1] == 1L) psign <- wskel[psign,2]
          sign[parnums, 2] <- psign
        }
      }
    }
  }
  
  out <- list(matskel = matskel, free2 = free2, free = free, wskel = wskel,
              sign = sign)

  return(out)
}

lav2stanmarg <- function(lavobject, dp) {
  ## extract model and data characteristics from lavaan object
  dat <- list()
  opts <- lavInspect(lavobject, 'options')

  ## data
  Ng <- dat$Ng <- lavInspect(lavobject, 'ngroups')
  YX <- lavobject@Data@X
  nvar <- ncol(YX[[1]])
  dat$N <- lavInspect(lavobject, 'nobs')

  ## lavobject@SampleStats@missing.flag is TRUE when missing='ml',
  ## regardless of whether data are missing
  misflag <- any(sapply(lavobject@Data@X, function(x) any(is.na(x))))
  if (misflag) {
    dat$miss <- 1L
    Mp <- lavobject@Data@Mp
    cases <- lapply(Mp, function(x) do.call("c", x$case.idx))
    misgrps <- lapply(Mp, function(x) x$freq)
#   if (length(misgrps) > 1) {
    misgrps <- do.call("c", misgrps)
#   }
    misgrps <- rep(1:length(misgrps), misgrps)
    npatt <- sapply(Mp, function(x) NROW(x$pat))
    dat$startrow <- tapply(1:NROW(misgrps), misgrps, head, 1)
    dat$endrow <- tapply(1:NROW(misgrps), misgrps, tail, 1)
    dat$grpnum <- rep(1:dat$Ng, npatt)

    dat$Nobs <- do.call("c", lapply(Mp, function(x) rowSums(x$pat)))
    Obsvar <- do.call("c", lapply(Mp, function(x) apply(x$pat, 1, which)))
    
    dat$Np <- length(unique(misgrps))
    dat$Ntot <- sum(dat$N)
    dat$Obsvar <- matrix(0, dat$Np, nvar)
    
    for (i in 1:dat$Np) {
      dat$Obsvar[i, 1:dat$Nobs[i]] <- Obsvar[[i]]
    }

    for (g in 1:dat$Ng) {
      YX[[g]] <- YX[[g]][cases[[g]],]
      YX[[g]] <- t(apply(YX[[g]], 1, function(x) c(x[!is.na(x)], rep(0, sum(is.na(x))))))
    }
  } else {
    dat$miss <- 0L
    dat$grpnum <- rep(1:dat$Ng, dat$N)
    dat$startrow <- tapply(1:NROW(dat$grpnum), dat$grpnum, head, 1)
    dat$endrow <- tapply(1:NROW(dat$grpnum), dat$grpnum, tail, 1)
    dat$grpnum <- unique(dat$grpnum)
    dat$Np <- dat$Ng
    dat$Nobs <- array(nvar, dat$Np)
    dat$Obsvar <- matrix(1:nvar, dat$Np, nvar, byrow=TRUE)
  }
  dat$YX <- do.call("rbind", YX)
  dat$grpnum <- array(dat$grpnum, length(dat$grpnum))

  ## model
  freemats <- lavInspect(lavobject, 'free')
  constrain <- attr(freemats, 'header')
  if (any(constrain$op != "==")){
    ops <- unique(constrain$op)
    ops <- ops[ops != "=="]
    stop(paste("blavaan ERROR: cannot handle constraints with ",
               paste(ops, collapse=" ")))
  }
  estmats <- lavInspect(lavobject, 'est')
  if (Ng == 1) {
    freemats <- list(freemats)
    estmats <- list(estmats)
  }
  free2 <- list()
  
  ## 1. Lambda_y
  if ("lambda" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x) x$lambda)
    es <- lapply(estmats, function(x) x$lambda)
    res <- matattr(fr, es, constrain, mat = "Lambda_y", Ng, opts$std.lv)

    dat$Lambda_y_skeleton <- res$matskel
    dat$w1skel <- res$wskel
    dat$lam_y_sign <- res$sign
    lyfree2 <- res$free2
    free2 <- c(free2, list(lambda = res$free))
  } else {
    dat$Lambda_y_skeleton <- matrix(0, 0, 0)
    dat$w1skel <- matrix(0, 0, 2)
    dat$lam_y_sign <- matrix(0, 0, 2)
  }

  ## 2. Lambda_x; never used because x only pops up in
  ##    the conditional case.
  dat$Lambda_x_skeleton <- matrix(0, 0, 0)
  dat$w2skel <- matrix(0, 0, 2)
  dat$lam_x_sign <- matrix(0, 0, 2)

  ## 3. Gamma
  if ("gamma" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x) x$gamma)
    es <- lapply(estmats, function(x) x$gamma)
    ## Note: if Lambda_x were in use, then FALSE below needs to
    ##       be opts$std.lv:
    res <- matattr(fr, es, constrain, mat = "Gamma", Ng, FALSE)

    dat$Gamma_skeleton <- res$matskel
    dat$w3skel <- res$wskel
    dat$gam_sign <- res$sign
    free2 <- c(free2, list(gamma = res$free))
  } else {
    dat$Gamma_skeleton <- matrix(0, ncol(dat$Lambda_y_skeleton), 0)
    dat$w3skel <- matrix(0, 0, 2)
    dat$gam_sign <- matrix(0, 0, 3)
  }

  ## 4. Beta
  if ("beta" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x) x$beta)
    es <- lapply(estmats, function(x) x$beta)
    res <- matattr(fr, es, constrain, mat = "B", Ng, opts$std.lv,
                   free2 = lyfree2, sign = dat$lam_y_sign)

    dat$B_skeleton <- res$matskel
    dat$w4skel <- res$wskel
    dat$b_sign <- res$sign
    free2 <- c(free2, list(beta = res$free))
  } else {
    dat$B_skeleton <- matrix(0, ncol(dat$Lambda_y_skeleton), 0)
    dat$w4skel <- matrix(0, 0, 2)
    dat$b_sign <- matrix(0, 0, 3)
  }
  
  ## 5. diag(Theta)
  if ("theta" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x){
      dmat <- x$theta
      dmat[lower.tri(dmat)] <- dmat[upper.tri(dmat)] <- 0
      dmat}
      )
    
    es <- lapply(estmats, function(x){
      dmat <- x$theta
      dmat[lower.tri(dmat)] <- dmat[upper.tri(dmat)] <- 0
      dmat}
      )
    
    res <- matattr(fr, es, constrain, mat = "Theta", Ng, opts$std.lv)

    dat$Theta_skeleton <- res$matskel
    dat$w5skel <- res$wskel
    free2 <- c(free2, list(dtheta = res$free))
  } else {
    dat$Theta_skeleton <- matrix(0, 0, 0)
    dat$w5skel <- matrix(0, 0, 2)
  }

  ## 6. diag(Theta_x)
  if ("cov.x" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x){
      dmat <- x$cov.x
      dmat[lower.tri(dmat)] <- dmat[upper.tri(dmat)] <- 0
      dmat}
      )
    
    es <- lapply(estmats, function(x){
      dmat <- x$cov.x
      dmat[lower.tri(dmat)] <- dmat[upper.tri(dmat)] <- 0
      dmat}
      )
    
    res <- matattr(fr, es, constrain, mat = "Theta_x", Ng, opts$std.lv)

    dat$Theta_x_skeleton <- res$matskel
    dat$w6skel <- res$wskel
    free2 <- c(free2, list(dthetax = res$free))
  } else {
    dat$Theta_x_skeleton <- matrix(0, 0, 0)
    dat$w6skel <- matrix(0, 0, 2)
  }

  ## 7. Theta_r
  if ("theta" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x){
      dmat <- x$theta
      diag(dmat) <- 0L
      dmat}
      )
    
    es <- lapply(estmats, function(x){
      dmat <- x$theta
      diag(dmat) <- 1L
      dmat}
      )
    
    res <- matattr(fr, es, constrain, mat = "Theta_r", Ng, opts$std.lv)

    dat$Theta_r_skeleton <- res$matskel
    dat$w7skel <- res$wskel
    free2 <- c(free2, list(rtheta = res$free))
  } else {
    dat$Theta_r_skeleton <- matrix(0, 0, 0)
    dat$w7skel <- matrix(0, 0, 2)
  }

  ## 8. Theta_x_r
  if ("cov.x" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x){
      dmat <- x$cov.x
      diag(dmat) <- 0L
      dmat}
      )
    
    es <- lapply(estmats, function(x){
      dmat <- x$cov.x
      diag(dmat) <- 1L
      dmat}
      )
    
    res <- matattr(fr, es, constrain, mat = "Theta_x_r", Ng, opts$std.lv)

    dat$Theta_x_r_skeleton <- res$matskel
    dat$w8skel <- res$wskel
    free2 <- c(free2, list(cov.x = res$free))
  } else {
    dat$Theta_x_r_skeleton <- matrix(0, 0, 0)
    dat$w8skel <- matrix(0, 0, 2)
  }

  ## 9. diag(Psi)
  if ("psi" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x){
      dmat <- x$psi
      dmat[lower.tri(dmat)] <- dmat[upper.tri(dmat)] <- 0
      dmat}
      )
    
    es <- lapply(estmats, function(x){
      dmat <- x$psi
      dmat[lower.tri(dmat)] <- dmat[upper.tri(dmat)] <- 0
      dmat}
      )

    ## std.lv only matters for off-diagonals
    res <- matattr(fr, es, constrain, mat = "Psi", Ng, FALSE)

    dat$Psi_skeleton <- res$matskel
    dat$w9skel <- res$wskel
    free2 <- c(free2, list(dpsi = res$free))
  } else {
    dat$Psi_skeleton <- matrix(0, 0, 0)
    dat$w9skel <- matrix(0, 0, 2)
  }

  ## 10. Psi_r
  if ("psi" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x){
      dmat <- x$psi
      diag(dmat) <- 0L
      dmat}
      )
    
    es <- lapply(estmats, function(x){
      dmat <- x$psi
      diag(dmat) <- 1L
      dmat}
      )
    
    res <- matattr(fr, es, constrain, mat = "Psi_r", Ng, opts$std.lv,
                   free2 = lyfree2, sign = dat$lam_y_sign)

    dat$Psi_r_skeleton <- res$matskel
    dat$w10skel <- res$wskel
    dat$psi_r_sign <- res$sign
    free2 <- c(free2, list(rpsi = res$free))
  } else {
    dat$Psi_r_skeleton <- matrix(0, 0, 0)
    dat$w10skel <- matrix(0, 0, 2)
    dat$psi_r_sign <- matrix(0, 0, 3)
  }

  ## 11. Phi unused
  dat$Phi_skeleton <- matrix(0, 0, 0)
  dat$w11skel <- matrix(0, 0, 2)

  ## 12. Phi_r unused
  dat$Phi_r_skeleton <- matrix(0, 0, 0)
  dat$w12skel <- matrix(0, 0, 2)
  dat$phi_r_sign <- matrix(0, 0, 3)

  ## 13. Nu NB: unlike lavaan, we paste mean.x to end!!
  if ("nu" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x) {
      out <- x$nu
      if ("mean.x" %in% names(x)) out <- rbind(out, x$mean.x)
      out}
      )
    es <- lapply(estmats, function(x) {
      out <- x$nu
      if ("mean.x" %in% names(x)) out <- rbind(out, x$mean.x)
      out}
      )
    res <- matattr(fr, es, constrain, mat = "Nu", Ng, opts$std.lv)

    dat$Nu_skeleton <- res$matskel
    dat$w13skel <- res$wskel
    free2 <- c(free2, list(nu = res$free))
  } else {
    dat$Nu_skeleton <- matrix(0, 0, 0)
    dat$w13skel <- matrix(0, 0, 2)
  }

  ## 14. Alpha
  if ("alpha" %in% names(freemats[[1]])) {
    fr <- lapply(freemats, function(x) x$alpha)
    es <- lapply(estmats, function(x) x$alpha)
    res <- matattr(fr, es, constrain, mat = "Alpha", Ng, opts$std.lv)

    dat$Alpha_skeleton <- res$matskel
    dat$w14skel <- res$wskel
    free2 <- c(free2, list(alpha = res$free))
  } else {
    dat$Alpha_skeleton <- matrix(0, 0, 0)
    dat$w14skel <- matrix(0, 0, 2)
  }

  lavpartable <- parTable(lavobject)

  ## add priors by using set_stanpars() from classic approach
  lavpartable <- lavMatrixRepresentation(lavpartable, add.attributes = TRUE)
  lavpartable$rhoidx <- rep(0, length(lavpartable$mat))
  if (!("prior" %in% names(lavpartable))) lavpartable$prior <- rep("", length(lavpartable$mat))
  covpars <- with(lavpartable, mat %in% c("theta", "psi") &
                               row != col)
  lavpartable$mat[covpars] <- "rho"
  
  stanprires <- set_stanpars("", lavpartable, free2, dp, "")
  lavpartable$prior <- stanprires$partable$prior
  
  return(list(dat = dat, free2 = free2, lavpartable = lavpartable))
}


coeffun_stanmarg <- function(lavpartable, lavfree, free2, lersdat, rsob, fun = "mean") {
  ## Extract posterior means from marginal stan model.
  ## free2 comes from lav2lers().
  ## lersdat is data passed to sem stan code.
  ## rsob is the result of sampling().
  stanfit <- !is.null(rsob)
  if(stanfit){
    rssumm <- rstan::summary(rsob)
    rsmcmc <- as.array(rsob)

    ## posterior means:
    if(fun == "mean"){
      b.est <- rssumm$summary[,"mean"]
    } else if(fun == "median"){
      b.est <- rssumm$summary[,"50%"]
    }
    sd.est <- rssumm$summary[,"sd"]
  }

  draw_mat <- as.matrix(rsob)
  
  ## lavaan pars to stan par vectors
  mapping <- c(ly_sign = "lambda", g_sign = "gamma",
               bet_sign = "beta", Theta_cov = "theta",
               Theta_var = "theta", Theta_x_cov = "cov.x",
               Theta_x_var = "cov.x", Psi_cov = "psi",
               Psi_var = "psi", Nu_free = "nu", ## includes mean.x!
               Alpha_free = "alpha")

  ## lavaan pars to w?skel (for equality constraints)
  mapping2 <- c("lambda", "gamma", "beta", "theta",
                "theta", "cov.x", "cov.x", "psi",
                "psi", "nu", "alpha")
  names(mapping2) <- as.character(c(1, 3, 4, 7, 5, 8, 6, 10, 9,
                                    13, 14))

  ## stan pars to free2 pars
  mapping3 <- c(lambda = "ly_sign", gamma = "g_sign",
                beta = "bet_sign", rtheta = "Theta_cov",
                dtheta = "Theta_var", rtheta_x = "Theta_x_cov",
                dtheta_x = "Theta_x_var", rpsi = "Psi_cov",
                dpsi = "Psi_var", nu = "Nu_free",
                alpha = "Alpha_free")
  
  ## check names in lavfree
  if(!all(names(lavfree) %in% mapping)){
    ## multiple groups?
    if(!all(names(lavfree[[1]]) %in% mapping)){
      stop("blavaan ERROR: unrecognized lavaan model matrix.")
    }
    ngrp <- length(lavfree)
  } else {
    ## one group, make a list out of it
    lavfree <- list(lavfree)
    ngrp <- 1
  }

  freeidx <- lapply(lavfree, function(x) lapply(x, function(y) which(y > 0, arr.ind = TRUE)))
  freenums <- lapply(free2, function(x) lapply(x, function(y) y[y > 0]))
  nfree <- max(sapply(lavfree, function(x) sapply(x, max)))

  freevec <- rep(NA, nfree)
  rowidx <- rowidx2 <- rep(NA, nfree) # row index of stan est and summary containing the parameters (for vcorr)

  ## 1. get free par vector
  ## 2. expand it using w?skel for eq constraints
  ## 3. fill "x" in using lavaan free
  ## 4. record freeidx, double-counting free parameters
  est <- sdvec <- rep(NA, nfree)

  for(m in 1:length(freeidx[[1]])){
    stanvec <- names(mapping)[mapping == names(freeidx[[1]])[m]]
    wskel <- names(mapping2)[mapping == names(freeidx[[1]])[m]]
    wvec <- paste0("w", wskel)
    wskel <- paste0(wvec, "skel")

    ## 2 for cov/var vectors, 1 otherwise
    if(length(stanvec) > 2) stop("blavaan ERROR: problem with mapping from stan to lavaan")

    for(j in 1:length(stanvec)){
      freename <- names(mapping3)[mapping3 == stanvec[j]]
      parnums <- do.call("c", freenums[[freename]])      
      tmpw <- lersdat[[wskel[j]]]

      if(is.na(parnums[1])) next

      if(any(!is.finite(lersdat[[wvec[j]]]))){
        tmpw <- tmpw[1:(lersdat[['Ng']] * sum(!is.finite(lersdat[[wvec[j]]]))), , drop=FALSE]
      } else {
        tmpw <- NULL
      }

      if(NROW(tmpw) > 0){
        ## need rowvec & rowvec2 because stan summary rows
        ## ordered differently from stan draws rows
        parvec <- tmpsd <- rowvec <- rowvec2 <- rep(NA, NROW(tmpw))
        rowvec[tmpw[,1] == 0] <- grep(stanvec[j], names(b.est))
        rowvec2[tmpw[,1] == 0] <- grep(stanvec[j], colnames(draw_mat))
        parvec[tmpw[,1] == 0] <- b.est[rowvec[tmpw[,1] == 0]]
        tmpsd[tmpw[,1] == 0] <- sd.est[rowvec[tmpw[,1] == 0]]

        eqconst <- tmpw[,2][tmpw[,1] == 1]
        rowvec[tmpw[,1] == 1] <- rowvec[tmpw[,1] == 0][eqconst]
        rowvec2[tmpw[,1] == 1] <- rowvec2[tmpw[,1] == 0][eqconst]
        parvec[tmpw[,1] == 1] <- parvec[eqconst]
        tmpsd[tmpw[,1] == 1] <- tmpsd[eqconst]

        rowidx[parnums] <- rowvec
        rowidx2[parnums] <- rowvec2
        est[parnums] <- parvec
        sdvec[parnums] <- tmpsd
      }
    }
  }

  vcorr <- cor(draw_mat[,rowidx2])

  names(sdvec) <- colnames(vcorr)

  ## indexing of stan objects
  lavpartable$stanpnum <- rep(NA, length(lavpartable$est))
  lavpartable$stansumnum <- rep(NA, length(lavpartable$est))
  lavpartable$stanpnum[lavpartable$free > 0] <- rowidx
  lavpartable$stansumnum[lavpartable$free > 0] <- rowidx2

  ## add matrices, which are used in blavInspect()
  lavpartable <- lavMatrixRepresentation(lavpartable, add.attributes = TRUE, as.data.frame. = FALSE)
  
  list(x = est, vcorr = vcorr, sd = sdvec,
       stansumm = rssumm$summary, lavpartable = lavpartable)
}
