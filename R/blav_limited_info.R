### R-side limited-information (pairwise) posterior predictive check for
### models with ordinal indicators.
###
### Motivation: the default Stan-computed ppp compares the fitted model to
### a saturated model on the continuous data-augmented y* used internally
### for ordinal likelihoods. Because that augmentation already respects the
### fitted model's own thresholds, the comparison understates ordinal lack
### of fit. This file instead compares observed vs model-implied pairwise
### summaries (contingency tables for ordinal-ordinal pairs, a conditional/
### polyserial statistic for ordinal-continuous pairs, correlations for
### continuous-continuous pairs), separately for the real data and for a
### posterior-predictive replicate, per retained MCMC draw.
###
### Scope: single-level (single- or multi-group) target="stan"/"cmdstan"
### models only; see R/ppp.R for the dispatcher that keeps this out of the
### two-level and continuous-only paths.

## classify every pair of observed variables (assumes structurally
## identical ov.names/types across groups, as blavaan requires) as
## oo (ordinal-ordinal), oc (ordinal-continuous), or cc (continuous-continuous)
li_pairtypes <- function(lavdata, lavsamplestats) {
  ov.names <- lavdata@ov.names[[1]]
  nvar <- length(ov.names)
  ov.type <- lavdata@ov$type[match(ov.names, lavdata@ov$name)]
  is.ord <- ov.type == "ordered"

  ## a cc pair is dropped only if fixed.x in *every* group (conservative:
  ## a pair that is fixed.x in some but not all groups is kept)
  ngroups <- lavsamplestats@ngroups
  fixedx <- rep(TRUE, nvar)
  for (g in 1:ngroups) {
    x.idx <- lavsamplestats@x.idx[[g]]
    thisfx <- rep(FALSE, nvar)
    if (length(x.idx) > 0) thisfx[x.idx] <- TRUE
    fixedx <- fixedx & thisfx
  }

  oo <- oc <- cc <- matrix(numeric(0), 0, 2)
  if (nvar >= 2) {
    allpairs <- combn(nvar, 2)
    for (r in seq_len(ncol(allpairs))) {
      i <- allpairs[1, r]; j <- allpairs[2, r]
      if (is.ord[i] && is.ord[j]) {
        oo <- rbind(oo, c(i, j))
      } else if (is.ord[i] || is.ord[j]) {
        ordv <- if (is.ord[i]) i else j
        contv <- if (is.ord[i]) j else i
        oc <- rbind(oc, c(ordv, contv))
      } else if (!(fixedx[i] && fixedx[j])) {
        cc <- rbind(cc, c(i, j))
      }
    }
  }

  list(ov.names = ov.names, is.ord = is.ord, oo = oo, oc = oc, cc = cc)
}

## fixed (non-draw-dependent) observed-data summaries for oc and cc pairs.
## oo pairs need no bespoke summary here -- handled by lavTables()
## directly off lavdata@X inside li_oo_stat().
li_obs_summary <- function(lavdata, pairs) {
  ngroups <- lavdata@ngroups
  oc_obs <- vector("list", ngroups)
  cc_obs <- vector("list", ngroups)
  for (g in 1:ngroups) {
    X <- lavdata@X[[g]]
    if (nrow(pairs$oc) > 0) {
      oc_obs[[g]] <- lapply(seq_len(nrow(pairs$oc)), function(r) {
        ordv <- pairs$oc[r, 1]; contv <- pairs$oc[r, 2]
        ok <- complete.cases(X[, c(ordv, contv)])
        list(ord = X[ok, ordv], cont = X[ok, contv])
      })
    }
    if (nrow(pairs$cc) > 0) {
      cc_obs[[g]] <- lapply(seq_len(nrow(pairs$cc)), function(r) {
        i <- pairs$cc[r, 1]; j <- pairs$cc[r, 2]
        ok <- complete.cases(X[, c(i, j)])
        list(r = cor(X[ok, i], X[ok, j]), n = sum(ok))
      })
    }
  }
  list(oc = oc_obs, cc = cc_obs)
}

## build a non-refit ("fake fitted") lavaan object from one posterior
## draw's filled-in lavmodel and a data list (either the original observed
## data, or a posterior-predictive replicate from postdata()). Mirrors the
## existing hack in postpred.R's "other measures" branch (ugly hack to
## avoid lav_samplestats_from_data: reconstruct data + call lavaan()).
li_build_fake_lavobj <- function(datalist, lavmodel_i, lavpartable, lavoptions,
                                 lavcache, lavdata) {
  DATA.X <- do.call("rbind", datalist)
  colnames(DATA.X) <- lavdata@ov.names[[1L]]
  DATA <- as.data.frame(DATA.X)

  lavoptions2 <- lavoptions
  lavoptions2$verbose <- FALSE
  lavoptions2$estimator <- if (lavmodel_i@categorical) "DWLS" else "ML"
  lavoptions2$se <- "none"
  lavoptions2$test <- "standard"
  lavoptions2$optim.method <- "none"
  if ("control" %in% slotNames(lavmodel_i)) {
    lavmodel_i@control <- list(optim.method = "none")
  }

  ngroups <- length(datalist)
  if (ngroups > 1L) {
    DATA[, lavdata@group] <- rep(lavdata@group.label,
                                 times = sapply(datalist, nrow))
    out <- lavaan(slotOptions = lavoptions2, slotParTable = lavpartable,
                          slotSampleStats = NULL, slotData = NULL,
                          slotModel = lavmodel_i, slotCache = lavcache,
                          data = DATA, group = lavdata@group)
  } else {
    out <- lavaan(slotOptions = lavoptions2, slotParTable = lavpartable,
                          slotSampleStats = NULL, slotData = NULL,
                          slotModel = lavmodel_i, slotCache = lavcache,
                          data = DATA)
  }
  ## lavaan() with optim.method="none" leaves @implied empty
  implied <- lav_model_implied(lavmodel_i,
                               delta = (lavmodel_i@parameterization == "delta"))
  ## blavaan forces theta parameterization for ordinal models
  ## (R/blavaan.R), so implied$cov's diagonal is not 1 for ordinal
  ## variables' underlying y* -- but implied$th is already on the
  ## standardized (unit-variance) scale (lav_tables_pairwise_model_pi()
  ## uses th directly, with no rescaling of its own). Convert cov to a
  ## genuine correlation matrix so the two are on a consistent scale;
  ## otherwise lavTables()/pbivnorm sees spurious "correlations" > 1.
  implied$cov <- lapply(implied$cov, cov2cor)
  out@implied <- implied
  out
}

## thin wrapper around lavTables()'s pairwise (2-way) X2 table,
## summed across all ordinal-ordinal pairs. lavTables() enumerates ordinal
## pairs itself off fake_lavobj@Data@ov$type -- no pair list needed here.
## A draw with an out-of-[-1,1] implied correlation (rare, but possible for
## a poorly-mixed/divergent draw) makes lavTables()'s pbivnorm call error
## hard rather than degrade gracefully -- treat that draw's oo contribution
## as unusable (NA), mirroring get_ll_ord()'s try()-wrapped NA convention
## for its own mnormt::sadmvn() calls (R/blav_model_loglik.R).
li_oo_stat <- function(fake_lavobj) {
  tf <- try(lavTables(fake_lavobj, dimension = 2L, type = "table",
                      statistic = "X2"), silent = TRUE)
  if (inherits(tf, "try-error")) return(NA_real_)
  if (is.null(tf) || nrow(tf) == 0) return(0)
  sum(tf$X2, na.rm = TRUE)
}

## polyserial-style discrepancy for one ordinal-continuous pair (one group).
## No natural tabulation exists once one variable is continuous, so this is
## a per-case deviance: -2*sum(log(conditional probability of the observed
## category)), given the conditional distribution of the ordinal item's
## underlying y* conditional on the paired continuous variable's value
## (standard bivariate-normal conditioning). A squared-Pearson-residual
## form (obs-p)^2/p was tried first but rejected: it blows up whenever any
## single case's observed category has a small conditional probability
## (routine in the tails, even for a correctly-specified model, once
## multiplied out over many cases) -- a per-case deviance only grows
## logarithmically as p -> 0, so one unusual case can't dominate the sum.
li_oc_stat <- function(ord, cont, mu_o, mu_c, s_oo, s_cc, s_oc, tau_raw,
                       floor = 1e-8) {
  n <- length(ord)
  if (n == 0) return(0)
  rho <- s_oc / sqrt(s_oo * s_cc)
  mu_cond <- mu_o + rho * sqrt(s_oo / s_cc) * (cont - mu_c)
  sd_cond <- sqrt(s_oo * (1 - rho^2))
  cuts <- c(-Inf, tau_raw, Inf)

  upper <- cuts[ord + 1L]
  lower <- cuts[ord]
  catprob_obs <- pnorm((upper - mu_cond) / sd_cond) -
    pnorm((lower - mu_cond) / sd_cond)
  catprob_obs <- pmax(catprob_obs, floor)

  -2 * sum(log(catprob_obs))
}

## continuous-continuous pair discrepancy: squared standardized residual
## between the sample and model-implied correlation.
li_cc_stat <- function(r_data, r_model, n) {
  if (n == 0 || is.na(r_data)) return(0)
  n * (r_data - r_model)^2
}

## per-draw discrepancy, observed vs replicated data, summed across all
## oo/oc/cc pairs.
li_draw_stat <- function(postsamp, lavobject, lavmodel_orig, lavdata,
                         lavsamplestats, lavpartable, lavoptions, lavcache,
                         samp.index, chain.num, pairs, obs_summary, has_oo) {
  lavmodel_i <- fill_params(postsamp, lavmodel_orig, lavpartable)
  implied <- lav_model_implied(lavmodel_i,
                                       delta = (lavmodel_i@parameterization == "delta"))

  gen <- postdata(samp.indices = samp.index, chain.num = chain.num,
                  lavmodel = lavmodel_orig, lavdata = lavdata,
                  lavjags = lavobject@external$mcmcout, lavpartable = lavpartable,
                  lavsamplestats = lavsamplestats, lavobject = lavobject)
  dataX_rep <- gen[[1]]

  ngroups <- lavdata@ngroups
  T_obs <- 0
  T_rep <- 0

  if (has_oo) {
    obs_fake <- li_build_fake_lavobj(lavdata@X, lavmodel_i, lavpartable,
                                     lavoptions, lavcache, lavdata)
    T_obs <- T_obs + li_oo_stat(obs_fake)

    rep_fake <- li_build_fake_lavobj(dataX_rep, lavmodel_i, lavpartable,
                                     lavoptions, lavcache, lavdata)
    T_rep <- T_rep + li_oo_stat(rep_fake)
  }

  for (g in 1:ngroups) {
    Xg_rep <- dataX_rep[[g]]
    mu <- as.numeric(implied$mean[[g]])
    Sigma <- (implied$cov[[g]] + t(implied$cov[[g]])) / 2
    th <- implied$th[[g]]
    th.idx <- lavmodel_i@th.idx[[g]]

    if (nrow(pairs$oc) > 0) {
      for (r in seq_len(nrow(pairs$oc))) {
        ordv <- pairs$oc[r, 1]; contv <- pairs$oc[r, 2]
        tau_raw <- mu[ordv] + th[th.idx == ordv] * sqrt(Sigma[ordv, ordv])

        oo_obs <- obs_summary$oc[[g]][[r]]
        T_obs <- T_obs + li_oc_stat(oo_obs$ord, oo_obs$cont,
                                    mu[ordv], mu[contv],
                                    Sigma[ordv, ordv], Sigma[contv, contv],
                                    Sigma[ordv, contv], tau_raw)

        ok <- complete.cases(Xg_rep[, c(ordv, contv)])
        T_rep <- T_rep + li_oc_stat(Xg_rep[ok, ordv], Xg_rep[ok, contv],
                                    mu[ordv], mu[contv],
                                    Sigma[ordv, ordv], Sigma[contv, contv],
                                    Sigma[ordv, contv], tau_raw)
      }
    }

    if (nrow(pairs$cc) > 0) {
      for (r in seq_len(nrow(pairs$cc))) {
        iv <- pairs$cc[r, 1]; jv <- pairs$cc[r, 2]
        r_model <- Sigma[iv, jv] / sqrt(Sigma[iv, iv] * Sigma[jv, jv])

        cc_obs <- obs_summary$cc[[g]][[r]]
        T_obs <- T_obs + li_cc_stat(cc_obs$r, r_model, cc_obs$n)

        ok <- complete.cases(Xg_rep[, c(iv, jv)])
        r_rep <- cor(Xg_rep[ok, iv], Xg_rep[ok, jv])
        T_rep <- T_rep + li_cc_stat(r_rep, r_model, sum(ok))
      }
    }
  }

  c(T_obs = T_obs, T_rep = T_rep)
}

## top-level entry point: limited-information posterior predictive p-value
## for single-level (single- or multi-group) ordinal/mixed stan/cmdstan
## fits. Mirrors postpred()'s scaffolding and return shape.
pp_limited_info <- function(lavobject, thin = 1, parallel = TRUE) {
  lavpartable <- lavobject@ParTable
  lavmodel <- lavobject@Model
  lavoptions <- lavobject@Options
  lavsamplestats <- lavobject@SampleStats
  lavdata <- lavobject@Data
  lavcache <- lavobject@Cache
  lavjags <- lavobject@external$mcmcout

  lavmcmc <- make_mcmc(lavjags)
  n.chains <- length(lavmcmc)
  samp.indices <- sampnums(lavjags, thin = thin)
  psamp <- length(samp.indices)

  pairs <- li_pairtypes(lavdata, lavsamplestats)
  has_oo <- nrow(pairs$oo) > 0
  obs_summary <- li_obs_summary(lavdata, pairs)

  loop.args <- list(X = 1:n.chains, FUN = function(j) {
    sapply(1:psamp, function(i) {
      li_draw_stat(lavmcmc[[j]][samp.indices[i], ], lavobject, lavmodel,
                  lavdata, lavsamplestats, lavpartable, lavoptions, lavcache,
                  samp.indices[i], j, pairs, obs_summary, has_oo)
    })
  })

  if (parallel) {
    loop.args <- c(loop.args, future.seed = TRUE)
    res <- do.call("future_lapply", loop.args)
  } else {
    res <- do.call("lapply", loop.args)
  }

  T_obs <- unlist(lapply(res, function(x) x["T_obs", ]))
  T_rep <- unlist(lapply(res, function(x) x["T_rep", ]))

  ## draws with an out-of-range implied correlation (see li_oo_stat()) are
  ## NA and excluded from the ppp average, matching how other blavaan
  ## per-draw computations (e.g. get_ll_ord()) drop try()-error draws
  nna <- sum(is.na(T_obs) | is.na(T_rep))
  if (nna > 0) {
    warning("blavaan WARNING: ", nna, " of ", length(T_obs),
            " posterior draws produced an out-of-range model-implied ",
            "correlation and were excluded from the ppp computation.",
            call. = FALSE)
  }

  list(ppval = mean(T_rep > T_obs, na.rm = TRUE),
       ppdist = list(obs = T_obs, reps = T_rep))
}
