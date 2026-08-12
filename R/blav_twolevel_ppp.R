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
### Scope: two-level target="stan"/"cmdstan" models without fixed.x
### variables. fixed.x is not yet supported (see tl_replicate_data()); see
### R/ppp.R for the dispatcher that reaches this file.

## simulate one two-level posterior-predictive replicate dataset from one
## posterior draw's implied within/between moments. Mirrors what Stan used
## to do in generated quantities (multi_normal_cholesky_rng on the implied
## Sigma.W/Sigma.B), but in R via mnormt::rmnorm.
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

    W <- matrix(rmnorm(n_g, mean = Mu.W, varcov = Sigma.W),
               nrow = n_g, ncol = length(ov_idx1))
    B <- matrix(rmnorm(n_clusters, mean = Mu.B, varcov = Sigma.B),
               nrow = n_clusters, ncol = length(ov_idx2))

    Y_rep_g <- matrix(0, n_g, p)
    colnames(Y_rep_g) <- colnames(X_g)
    ## within-only and "both" columns start out as the within-level draw...
    Y_rep_g[, ov_idx1] <- W
    ## ...and "both"/between-only columns additionally get the cluster's
    ## between-level draw added in, so "both" columns end up W + B while
    ## between-only columns end up B alone (they were never touched by the
    ## line above)
    Y_rep_g[, ov_idx2] <- Y_rep_g[, ov_idx2, drop = FALSE] +
      B[cluster_idx, , drop = FALSE]

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

  ## fixed.x is not yet supported by tl_replicate_data()'s replicate
  ## generator (see file header); fail clearly up front rather than
  ## mis-simulating
  if (any(sapply(lavsamplestats@x.idx, length) > 0)) {
    stop("blavaan ERROR: on-demand ppp computation for two-level models with ",
         "fixed.x variables is not yet supported.")
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
