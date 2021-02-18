## compute undirected K-L divergence across all draws
## (each draw paired with one from another chain)
samp_kls <- function(lavjags        = NULL,
                     lavmodel       = NULL,
                     lavpartable    = NULL,
                     lavsamplestats = NULL,
                     lavoptions     = NULL,
                     lavcache       = NULL,
                     lavdata        = NULL,
                     lavmcmc        = NULL,
                     lavobject      = NULL,
                     thin           = 1,
                     conditional    = FALSE){

    ## need to implement plummer's approach of generating y_rep
    ##mis <- FALSE
    ##if(any(is.na(unlist(lavdata@X)))) mis <- TRUE
    ##if(mis | lavoptions$categorical) stop("blavaan ERROR: K-L divergence not implemented for missing data or ordinal variables.")

    itnums <- sampnums(lavjags, thin = thin)
    lavmcmc <- lapply(lavmcmc, function(x) x[itnums,])
    draws <- do.call("rbind", lavmcmc)

    ndraws <- nrow(draws)
    halfdraws <- floor(ndraws/2)
    ngroups <- lavsamplestats@ngroups

    klres <- rep(NA, halfdraws)
    for(i in 1:halfdraws){
        lavmodel0 <- fill_params(draws[i,], lavmodel, lavpartable)
        lavmodel1 <- fill_params(draws[(halfdraws + i),], lavmodel,
                                 lavpartable)

        if(conditional){
            eta0 <- fill_eta(draws[i,], lavmodel, lavpartable,
                             lavsamplestats, lavdata)
            eta1 <- fill_eta(draws[(halfdraws + i),], lavmodel,
                             lavpartable, lavsamplestats, lavdata)

            lavobject@Model <- lavmodel0
            mnvec0 <- lavPredict(lavobject, type="ov", ETA=eta0)
            if(inherits(mnvec0, "matrix")) mnvec0 <- list(mnvec0)
            cmat0 <- lavInspect(lavobject, 'theta')
            if(inherits(cmat0, "matrix")) cmat0 <- list(cmat0)

            lavobject@Model <- lavmodel1
            mnvec1 <- lavPredict(lavobject, type="ov", ETA=eta1)
            if(inherits(mnvec1, "matrix")) mnvec1 <- list(mnvec1)
            cmat1 <- lavInspect(lavobject, 'theta')
            if(inherits(cmat1, "matrix")) cmat1 <- list(cmat1)

            implied0 <- list(cov = cmat0, mean = mnvec0,
                             slopes = vector("list", ngroups),
                             th = vector("list", ngroups),
                             group.w = vector("list", ngroups))
            implied1 <- list(cov = cmat1, mean = mnvec1,
                             slopes = vector("list", ngroups),
                             th = vector("list", ngroups),
                             group.w = vector("list", ngroups))
        } else {
            implied0 <- lav_model_implied(lavmodel0)
            implied1 <- lav_model_implied(lavmodel1)
        }

        tmpkl <- 0
        for(g in 1:lavsamplestats@ngroups){
            ## ensure symmetric:
            cmat0 <- (implied0$cov[[g]] + t(implied0$cov[[g]]))/2
            invcmat0 <- solve(cmat0)
            det0 <- det(cmat0)
            cmat1 <- (implied1$cov[[g]] + t(implied1$cov[[g]]))/2
            invcmat1 <- solve(cmat1)
            det1 <- det(cmat1)
            if(conditional){
                mnvec0 <- implied0$mean[[g]]
                mnvec1 <- implied1$mean[[g]]

                for(j in 1:nrow(mnvec0)){
                  tmpkl <- tmpkl + kl_und(mnvec0[j,], mnvec1[j,],
                                          cmat0, invcmat0, cmat1,
                                          invcmat1, det0, det1)
                }
            } else {
                mnvec0 <- as.numeric(implied0$mean[[g]])
                mnvec1 <- as.numeric(implied1$mean[[g]])

                tmpkl <- tmpkl + lavsamplestats@nobs[[g]] *
                         kl_und(mnvec0, mnvec1, cmat0, invcmat0,
                                cmat1, invcmat1, det0, det1)
            }
        }
        klres[i] <- tmpkl
    }
    klres
}

## fill in eta matrices (1 per group, in list)
fill_eta <- function(postsamp, lavmodel, lavpartable, lavsamplestats, lavdata){
    nlv <- length(lavmodel@GLIST$alpha)
    etapars <- grepl("^eta", names(postsamp))
    cnums <- strsplit(names(postsamp)[etapars], "\\[|,|\\]")
    cnums <- sapply(cnums, function(x) as.numeric(x[3]))
    etavec <- postsamp[etapars][order(cnums)]

    ## need to worry about (1) excluding phantom lvs
    ## and (2) including dummy lvs
    foundlvs <- sum(etapars)/lavsamplestats@ntotal
    etamat <- matrix(etavec, lavsamplestats@ntotal, foundlvs)
    if(foundlvs < nlv) etamat <- cbind(etamat, matrix(0, lavsamplestats@ntotal, (nlv - foundlvs)))

    ## fulleta needs to have rows for any excluded cases
    if(sum(unlist(lavdata@norig)) > lavsamplestats@ntotal){
      fulleta <- matrix(NA, sum(unlist(lavdata@norig)), ncol(etamat))
      empties <- as.numeric(sapply(lavdata@Mp, function(x) x$empty.idx))
      fulleta[-empties,] <- etamat
    } else {
      fulleta <- etamat
    }

    ngroups <- lavsamplestats@ngroups
    eta <- vector("list", ngroups)
    for(g in 1:ngroups){
      eta[[g]] <- fulleta[lavdata@case.idx[[g]], 1:nlv, drop = FALSE]

      ## fill in eta with dummys, if needed
      dummyov <- c(lavmodel@ov.x.dummy.ov.idx[[g]], lavmodel@ov.y.dummy.ov.idx[[g]])
      dummylv <- c(lavmodel@ov.x.dummy.lv.idx[[g]], lavmodel@ov.y.dummy.lv.idx[[g]])
      if(length(dummyov) > 0){
        eta[[g]][, dummylv] <- lavdata@X[[g]][, dummyov]
      }
    }

    eta
}

## compute undirected K-L divergence between two normal distributions
kl_und <- function(mn0, mn1, cov0, invcov0, cov1, invcov1,
                   det0, det1){
  k <- nrow(cov0)

  kl01 <- sum(diag(invcov1 %*% cov0)) +
    t(mn1 - mn0) %*% invcov1 %*% (mn1 - mn0) -
    k + log(det1/det0)

  kl10 <- sum(diag(invcov0 %*% cov1)) +
    t(mn0 - mn1) %*% invcov0 %*% (mn0 - mn1) -
    k + log(det0/det1)

  (1/2) * (kl01 + kl10)
}

cond_moments <- function(postsamp, lavmodel, lavpartable, lavsamplestats, lavdata, lavobject){
  eta <- fill_eta(postsamp, lavmodel, lavpartable, lavsamplestats, lavdata)

  ## implied meanvec + covmat
  ##mnvec <- lavaan:::computeYHAT(lavmodel, lavmodel@GLIST,
  ##                              lavsamplestats, ETA = eta)
  lavobject@Model <- lavmodel
  mnvec <- lavPredict(lavobject, type="ov", ETA = eta)
  if(inherits(mnvec, "matrix")) mnvec <- list(mnvec)

  ##covmat <- lavaan:::computeTHETA(lavmodel, lavmodel@GLIST)
  covmat <- lavInspect(lavobject, 'theta')
  if(inherits(covmat, "matrix")) covmat <- list(covmat)
  ## to avoid warnings from mnormt::pd.solve
  covmat <- lapply(covmat, function(x){
    class(x) <- "matrix"
    zvar <- which(diag(x) == 0L)
    if(length(zvar) > 0) diag(x)[zvar] <- 1e-4
    x})

  ngroups <- lavsamplestats@ngroups
  implied <- list(cov = covmat, mean = mnvec,
                  slopes = vector("list", ngroups),
                  th = vector("list", ngroups),
                  group.w = vector("list", ngroups))

  implied
}
