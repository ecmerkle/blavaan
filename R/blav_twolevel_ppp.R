### R-side posterior predictive p-value (ppp) for two-level (multilevel)
### target="stan"/"cmdstan" models.
###
### Motivation: two-level ppp used to be computed inside stanmarg.stan via a
### custom Stan port (twolevel_em_step()) of lavaan's saturated-model EM
### algorithm. That machinery has been removed from the Stan model entirely
### (see R/blavaan.R's do_test computation); this file reimplements the same
### statistic -- a saturated-model likelihood-ratio chi-square ppp,
### T = -2*(ll_fit - ll_sat), compared for real vs. posterior-predictive-
### replicated data -- in R, calling lavaan's own lav_mvn_cl_em_sat()/
### lav_mvn_cl_mi_em_sat() (via R/blav_model_loglik.R's get_ll_2l_sat())
### rather than re-deriving the EM algorithm a second time.
###
### Scope: two-level target="stan"/"cmdstan" models, including fixed.x
### variables at the within level, the between level, or both (in different
### variables) -- see tl_replicate_data()/tl_cond_sim() for how fixed.x
### columns are held fixed at their observed values (never resimulated) and
### the endogenous columns are drawn conditional on them, mirroring
### R/postpred.R's single-level fixed.x approach. A single variable that is
### fixed.x at BOTH levels simultaneously is not supported (lavaan itself
### has an open FIXME for this combination -- see outstanding_issues.md item
### 3 in the repo root at the time this was written), nor is a between-level
### fixed.x variable with real missing data (blocked further upstream, at
### fit time, by a known lavaan crash -- R/blavaan.R already refuses to fit
### such a model at all, so pp_twolevel() never sees this combination). See
### R/ppp.R for the dispatcher that reaches this file.
###
### Why the fixed.x correction lavaan itself needs for its FIML likelihood
### (and the old Stan calc_log_lik_x()/calc_log_lik_x_missing() this
### replaces) turns out to be unnecessary here: T = -2*(ll_fit - ll_sat)
### only needs ll_fit and ll_sat to each include a fixed.x contribution that
### CANCELS in the difference, not to actually exclude it. get_ll_2l()/
### get_ll_2l_sat() (R/blav_model_loglik.R) both already evaluate the FULL
### joint density (endogenous + fixed.x jointly, via the complete implied/
### saturated Sigma) with loglik_x=0, i.e. no explicit correction. Verified
### empirically (both exactly, to floating-point precision, for complete
### data, and to EM convergence tolerance for missing data) that the
### fitted-model's implied Sigma's x-block and the saturated EM fit's
### x-block agree on the SAME data -- because both are, by construction,
### genuine (joint-)MLE-consistent estimates of x's own marginal
### distribution, which fixed.x holds structurally unconstrained in both
### the null and saturated model. So the joint-vs-conditional distinction is
### inert to T as long as x is held at its observed values (not
### resimulated) in the replicate, which is exactly what tl_cond_sim()
### below does.

## conditional multivariate-normal simulation for one level (within or
## between), holding a set of "x" columns fixed at their observed values
## (never resimulated) and drawing the remaining "y" columns from their
## conditional distribution given x. Handles per-row missingness in the x
## columns by grouping rows into observed-x patterns (mirrors
## R/postpred.R's per-missing-pattern fixed.x handling) -- a row with some
## (or all) x columns unobserved conditions on whatever subset of x it does
## have observed, degenerating to the unconditional draw when none are
## observed. mu/sigma are already restricted to this level's own
## active-variable ordering; x_local/y_local are positions within that
## ordering (not full-data column indices).
tl_cond_sim <- function(mu, sigma, x_local, y_local, Xobs, n) {
  if (length(x_local) == 0L) {
    return(matrix(rmnorm(n, mean = mu[y_local], varcov = sigma[y_local, y_local, drop = FALSE]),
                  nrow = n, ncol = length(y_local)))
  }

  out <- matrix(NA_real_, n, length(y_local))
  obs_mask <- !is.na(Xobs)
  ## group rows by which x columns are observed, so each group's conditional
  ## mean/covariance is only computed once
  patt_key <- apply(obs_mask, 1, function(r) paste(which(r), collapse = ","))

  for (patt in unique(patt_key)) {
    rows <- which(patt_key == patt)
    obs_cols <- which(obs_mask[rows[1], ])

    if (length(obs_cols) == 0L) {
      cmu <- matrix(mu[y_local], length(rows), length(y_local), byrow = TRUE)
      csig <- sigma[y_local, y_local, drop = FALSE]
    } else {
      xo_local <- x_local[obs_cols]
      Sxx_inv <- solve(sigma[xo_local, xo_local, drop = FALSE])
      Syx <- sigma[y_local, xo_local, drop = FALSE]
      csig <- sigma[y_local, y_local, drop = FALSE] - Syx %*% Sxx_inv %*% t(Syx)
      resid <- sweep(Xobs[rows, obs_cols, drop = FALSE], 2, mu[xo_local], "-")
      cmu <- matrix(mu[y_local], length(rows), length(y_local), byrow = TRUE) +
        resid %*% t(Syx %*% Sxx_inv)
    }
    csig <- (csig + t(csig)) / 2

    sims <- lapply(seq_along(rows), function(i) rmnorm(1, mean = cmu[i, ], varcov = csig))
    out[rows, ] <- matrix(unlist(sims), nrow = length(rows), ncol = length(y_local), byrow = TRUE)
  }

  out
}

## simulate one two-level posterior-predictive replicate dataset from one
## posterior draw's implied within/between moments. Mirrors what Stan used
## to do in generated quantities (multi_normal_cholesky_rng on the implied
## Sigma.W/Sigma.B), but in R via mnormt::rmnorm; fixed.x columns are copied
## unchanged from the observed data instead of being resimulated (see
## tl_cond_sim()).
tl_replicate_data <- function(lavdata, lavsamplestats, implied) {
  Ng <- lavsamplestats@ngroups
  dataX_rep <- vector("list", Ng)

  for (g in 1:Ng) {
    lp <- lavdata@Lp[[g]]
    X_g <- lavdata@X[[g]]
    n_g <- nrow(X_g)
    p <- ncol(X_g)

    cluster_idx <- lp$cluster.idx[[2]]
    n_clusters <- lp$nclusters[[2]]

    Mu.W <- as.numeric(implied$mean[[2 * g - 1]])
    Sigma.W <- implied$cov[[2 * g - 1]]
    Mu.B <- as.numeric(implied$mean[[2 * g]])
    Sigma.B <- implied$cov[[2 * g]]

    ## implied$mean[[.]]/implied$cov[[.]] are already sized to each level's
    ## OWN active-variable set (within: both + within-only; between: both +
    ## between-only) -- NOT embedded in the full "tilde" (all-variables)
    ## space with structural zeros. lp$ov.idx[[1]]/[[2]] give the
    ## corresponding full-data column for each row/col of Sigma.W/Sigma.B,
    ## in order (same convention lavaan's own lav_mvn_cl_implied22l() uses
    ## to scatter mu_w/mu_b into a p_tilde-length vector via
    ## mu_w_tilde[ov_idx[[1]]] <- mu_w). Do NOT subset Sigma.W/Sigma.B
    ## themselves with within.idx/between.idx/both.idx (those are tilde-
    ## universe indices, meant for p_tilde x p_tilde matrices, and index
    ## out of bounds here whenever within-only/between-only variables are
    ## present alongside "both" variables -- Sigma.B, e.g., has no row/col
    ## for a within-only variable at all).
    ov_idx1 <- lp$ov.idx[[1]]
    ov_idx2 <- lp$ov.idx[[2]]

    Sigma.W <- (Sigma.W + t(Sigma.W)) / 2
    Sigma.B <- (Sigma.B + t(Sigma.B)) / 2

    ## fixed.x columns (full-data-column indices); NULL when there are none
    xidx_w <- lp$ov.x.idx[[1]]
    if (is.null(xidx_w)) xidx_w <- integer(0)
    xidx_b <- lp$ov.x.idx[[2]]
    if (is.null(xidx_b)) xidx_b <- integer(0)

    w_local_x <- match(xidx_w, ov_idx1)
    w_local_y <- setdiff(seq_along(ov_idx1), w_local_x)
    Xw_obs <- if (length(xidx_w) > 0L) X_g[, xidx_w, drop = FALSE] else NULL
    Wy <- tl_cond_sim(Mu.W, Sigma.W, w_local_x, w_local_y, Xw_obs, n_g)
    W <- matrix(0, n_g, length(ov_idx1))
    W[, w_local_y] <- Wy
    if (length(xidx_w) > 0L) W[, w_local_x] <- Xw_obs

    b_local_x <- match(xidx_b, ov_idx2)
    b_local_y <- setdiff(seq_along(ov_idx2), b_local_x)
    if (length(xidx_b) > 0L) {
      ## between-level fixed.x is constant within a cluster and never has
      ## real missing data reaching here (blocked at fit time upstream) --
      ## one representative (first) row per cluster gives its value
      Xb_obs <- X_g[match(seq_len(n_clusters), cluster_idx), xidx_b, drop = FALSE]
    } else {
      Xb_obs <- NULL
    }
    By <- tl_cond_sim(Mu.B, Sigma.B, b_local_x, b_local_y, Xb_obs, n_clusters)
    B <- matrix(0, n_clusters, length(ov_idx2))
    B[, b_local_y] <- By
    if (length(xidx_b) > 0L) B[, b_local_x] <- Xb_obs

    Y_rep_g <- matrix(0, n_g, p)
    colnames(Y_rep_g) <- colnames(X_g)
    ## within-only and "both" columns start out as the within-level draw...
    Y_rep_g[, ov_idx1] <- W
    ## ...and "both"/between-only columns additionally get the cluster's
    ## between-level draw added in, so "both" columns end up W + B while
    ## between-only columns end up B alone (they were never touched by the
    ## line above). Fixed.x columns are unaffected either way, since they
    ## were copied directly above rather than drawn from W/B.
    Y_rep_g[, ov_idx2] <- Y_rep_g[, ov_idx2, drop = FALSE] +
      B[cluster_idx, , drop = FALSE]
    ## between-level fixed.x columns are "both.idx"-like only in the sense
    ## that they get column-summed above; since B's own fixed.x columns are
    ## exact observed copies (not draws) and within-level W never touches
    ## those same columns (a variable cannot be within-only AND
    ## between-fixed.x simultaneously), this still ends up exactly right:
    ## Y_rep_g[, xidx_b] == B[cluster_idx, b_local_x] == the observed value

    ## preserve the observed missingness pattern/positions exactly, so the
    ## replicate's Mp stays valid (unchanged) against lavdata@Mp[[g]]
    Y_rep_g[is.na(X_g)] <- NA

    dataX_rep[[g]] <- Y_rep_g
  }

  dataX_rep
}

## top-level entry point: two-level saturated-model-comparison ppp for
## target="stan"/"cmdstan" fits. Mirrors postpred()'s/pp_limited_info()'s
## scaffolding and return shape.
pp_twolevel <- function(lavobject, thin = 1, parallel = FALSE,
                        em_control = list(tol = 1e-3, max_iter = 50L,
                                          acceleration = "squarem")) {
  lavpartable <- lavobject@ParTable
  lavmodel <- lavobject@Model
  lavsamplestats <- lavobject@SampleStats
  lavdata <- lavobject@Data
  lavjags <- lavobject@external$mcmcout

  ## a variable that is fixed.x at BOTH levels simultaneously is not
  ## supported by tl_replicate_data() (lavaan itself has an open FIXME for
  ## this combination -- see the file header); fail clearly up front rather
  ## than mis-simulating. Ordinary within-only or between-only fixed.x is
  ## fine. Between-level fixed.x + real missing data can't reach here at
  ## all (R/blavaan.R already refuses to fit such a model, per a known
  ## upstream lavaan crash), so no separate check is needed for that.
  for (g in seq_len(lavsamplestats@ngroups)) {
    xw <- lavdata@Lp[[g]]$ov.x.idx[[1]]
    xb <- lavdata@Lp[[g]]$ov.x.idx[[2]]
    if (length(intersect(xw, xb)) > 0) {
      stop("blavaan ERROR: on-demand ppp computation for two-level models does ",
           "not yet support a variable that is fixed.x at both levels simultaneously.")
    }
  }

  lavmcmc <- make_mcmc(lavjags)
  n.chains <- length(lavmcmc)
  samp.indices <- sampnums(lavjags, thin = thin)
  psamp <- length(samp.indices)

  ## observed-data saturated-model logl does not depend on the posterior
  ## draw (lav_mvn_cl_em_sat()/lav_mvn_cl_mi_em_sat() take only data), so
  ## compute it once, with full (default, not em_control-loosened)
  ## convergence settings -- mirrors Stan's own one-time-per-chain EM fit
  ## for the observed data, vs. a cheaper per-draw refit for each replicate
  ll_sat_obs <- get_ll_2l_sat(dataX = NULL, lavobject = lavobject)

  tl_draw_stat <- function(postsamp) {
    lavmodel_i <- fill_params(postsamp, lavmodel, lavpartable)
    implied <- lav_model_implied(lavmodel_i,
                                 delta = (lavmodel_i@parameterization == "delta"))

    ll_fit_obs <- get_ll_2l(postsamp, lavobject, standata = NULL)

    dataX_rep <- tl_replicate_data(lavdata, lavsamplestats, implied)

    ll_fit_rep <- get_ll_2l(postsamp, lavobject, standata = NULL, dataX = dataX_rep)
    ll_sat_rep <- do.call(get_ll_2l_sat,
                          c(list(dataX = dataX_rep, lavobject = lavobject), em_control))

    c(T_obs = -2 * (ll_fit_obs - ll_sat_obs),
      T_rep = -2 * (ll_fit_rep - ll_sat_rep))
  }

  loop.args <- list(X = 1:n.chains, FUN = function(j) {
    sapply(1:psamp, function(i) {
      tl_draw_stat(lavmcmc[[j]][samp.indices[i], ])
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

  list(ppval = mean(T_rep > T_obs, na.rm = TRUE),
       ppdist = list(obs = T_obs, reps = T_rep))
}
