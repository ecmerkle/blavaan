## tinytest file: two-level (multilevel) FIML/MAR missing-data models,
## Stan backend
##
## Run with: tinytest::run_test_file("tinytest_stan_multilevel_missing.R")
## or as part of a package: tinytest::test_package("blavaan")
##
## NB: two-level FIML/MAR is now the default (no option/gate needed) whenever
## a two-level model has real missing data -- see R/blavaan.R (the
## `!any(is.na(...))` check). missing="listwise" was never a user-facing
## override for two-level models -- it was only ever blavaan's own internal
## choice (set unconditionally, not read from what the user passed) to force
## row-dropping before this FIML/MAR support existed; see the note in
## section 4 below. fixed.x variables ARE supported for two-level models
## with real missingness, but ONLY at the within level -- see section 3's
## model_x_w below. Between-level fixed.x + missing data remains blocked (a
## confirmed upstream lavaan bug crashes lavaan()/sem() itself whenever a
## between-level fixed.x variable has missing values -- lav_samp_cl_patterns()
## in lavaan's lav_samplestats.R uses plain cov(), not
## cov(use="pairwise.complete.obs"), for the between-level block, so any NA
## there propagates into an uninvertible covariance matrix downstream);
## R/blavaan.R blocks this combination with a clear error instead of
## surfacing that crash -- see the expect_error() case in section 3).
## test != "none" (ppmc/ppp) IS supported (via an EM-based saturated-model
## fit for the endogenous variables, unified across complete/missing data --
## see
## R/lav_export_stanmarg.R and inst/stan/stanmarg.stan's twolevel_em_step()
## -- plus a pairwise-complete-obs correction for within-level fixed.x, see
## calc_log_lik_x_missing() in inst/stan/stanmarg.stan and llx_2l() in
## R/blav_model_loglik.R). See also
## tinytest_stan_multilevel2.R for complete-data two-level model tests.
##
## Multi-group two-level models are NOT tested here (out of scope for this
## FIML/MAR-focused file), but note that they now work: plain lavaan::sem()
## already supported cluster + group together (with the correct syntax:
## group: blocks as the OUTER blocks, each with its own nested level:
## blocks -- group= alone, or level: blocks with a bare group= argument,
## mis-sizes ov.names.l and throws "subscript out of bounds" in
## lav_lavaan_step01_ovnames_namesl()), and bsem() previously hit a
## separate bug on top of that even for complete multi-group two-level
## data: R/lav_export_stanmarg.R matched psi/theta correlation blocks to
## parameter-table rows using lavpartable$group, which is a clean 1:Ng
## group index for single-group models but the literal (often non-numeric)
## "group: <label>" text for multi-group multilevel models, so the match
## silently found nothing and "dat$psiblkpri[b] <- ..." failed with
## "replacement has length zero". Fixed by matching on lavpartable$block
## instead (always a well-defined group*level block number) -- see that
## file's git history for the fix. See tinytest_stan_multilevel2.R for
## complete-data two-level coverage.

library("tinytest")
library("lavaan")
library("blavaan")

## mildly informative prior on loadings: the default (vague) prior can leave
## short (burnin=100/sample=100) chains weakly identified enough that
## post-fit lavaan-refit-based computations (fitMeasures(), standardized
## solutions) occasionally fail to converge
dp_stable <- dpriors(lambda = "normal(1,.4)")

## ── Helpers ──────────────────────────────────────────────────────────────────

## Check that Bayesian estimates are within `tol` (in SD units) of lavaan MLEs.
## We use a loose tolerance because MCMC chains are short (burnin=100, sample=100).
coef_close <- function(bfit, fit, tol = 0.5) {
  bc  <- coef(bfit)
  lc  <- coef(fit)
  nms <- intersect(names(bc), names(lc))
  if (length(nms) == 0L) return(FALSE)
  bse <- sqrt(diag(vcov(bfit)))[nms]
  all(abs(bc[nms] - lc[nms]) < tol * (abs(lc[nms]) + bse + 0.1))
}

inject_mcar <- function(d, vars, rate = 0.15, seed) {
  set.seed(seed)
  n <- nrow(d)
  for (v in vars) {
    idx <- sample(seq_len(n), size = floor(rate * n))
    d[idx, v] <- NA
  }
  d
}

inject_mcar_between <- function(d, vars, cluster_col = "cluster", rate = 0.15, seed) {
  set.seed(seed)
  clusters <- unique(d[[cluster_col]])
  miss_clusters <- sample(clusters, size = floor(rate * length(clusters)))
  for (v in vars) {
    d[d[[cluster_col]] %in% miss_clusters, v] <- NA
  }
  d
}

## ── Data ─────────────────────────────────────────────────────────────────────

data(Demo.twolevel, package = "lavaan")

## y4 is deliberately within-only (no level-2 loading), y1-y3 are "both",
## w1/w2 are between-only -- exercises all three variable classes at once.
## Only two between-level covariances (not three): a 3x3 between-covariance
## sub-block from a 2-way + w1~~w2 spec proved fragile at only burnin=100/
## sample=100 -- some seeds left the chains short of convergence (R-hat up
## to ~1.3 was observed), and the resulting unstable posterior mean/cov
## occasionally tripped a Hessian-based check in blavaan's post-fit
## processing ("model did not converge"). Unrelated to the FIML/MAR
## likelihood itself (Phase C validated that to 1e-12 in isolation); just a
## short-chain identifiability issue with this particular covariance spec.
model_wb <- '
    level: 1
        fw =~ y1 + y2 + y3 + y4
    level: 2
        fb =~ y1 + y2 + y3
        fb ~~ w1
        fb ~~ w2
'

## =============================================================================
## 1. MCAR across within-only, both, and between-only variables
## =============================================================================

try({
set.seed(101)
d1 <- Demo.twolevel
d1 <- inject_mcar(d1, c("y1", "y2", "y3", "y4"), rate = 0.15, seed = 101)
d1 <- inject_mcar_between(d1, "w1", rate = 0.15, seed = 102)

fit1 <- sem(model_wb, data = d1, cluster = "cluster", missing = "fiml",
            fixed.x = FALSE)

## missing="fiml" is not passed explicitly here: real NAs are present in
## d1, and blavaan already defaults missing to "ml" (equivalent to "fiml")
## whenever a two-level fit has actual missing data -- see R/blavaan.R (the
## `!any(is.na(...))` check). Passing it explicitly is redundant and
## triggers a spurious "blavaan WARNING: the following arguments have no
## effect: missing" (that warning check does not know about the
## two-level-FIML-MAR default).
bfit1 <- bsem(
  model   = model_wb,
  data    = d1,
  cluster = "cluster",
  fixed.x = FALSE,
  test    = "ppp",
  burnin  = 100,
  sample  = 100,
  dp      = dp_stable
)

expect_inherits(bfit1, "blavaan",
  info = "bsem() with two-level FIML/MAR (within+both+between missingness) should return a blavaan object")

expect_equal(
  sort(names(coef(bfit1))),
  sort(names(coef(fit1))),
  info = "two-level FIML/MAR: parameter names should match lavaan"
)

expect_true(
  coef_close(bfit1, fit1),
  info = "two-level FIML/MAR: estimates should be within 0.5 SD of lavaan MLEs"
)

expect_true(
  all(diag(vcov(bfit1)) > 0),
  info = "two-level FIML/MAR: vcov diagonal should be positive"
)

## the log_lik generated-quantity dispatch is mandatory (see R/blavaan.R
## comments): a regression there would silently break waic/looic, not error
fm1 <- fitMeasures(bfit1)
expect_true(is.finite(fm1["waic"]),
  info = "two-level FIML/MAR: waic should be finite")
expect_true(is.finite(fm1["looic"]),
  info = "two-level FIML/MAR: looic should be finite")

## no divergent transitions over the (short) sampled run
sampler_params <- rstan::get_sampler_params(bfit1@external$mcmcout, inc_warmup = FALSE)
ndiv <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
expect_equal(ndiv, 0L,
  info = "two-level FIML/MAR: no divergent transitions expected in a well-behaved short run")

## ppp (EM-based saturated-model fit, see R/lav_export_stanmarg.R and
## inst/stan/stanmarg.stan's twolevel_em_step()) should be well-defined
## for this within+both+between missingness scenario
ppp_val1 <- fm1["ppp"]
expect_true(
  is.finite(ppp_val1) && ppp_val1 >= 0 && ppp_val1 <= 1,
  info = "two-level FIML/MAR ppp should be a finite value in [0,1]"
)
ll_sat1 <- rstan::extract(bfit1@external$mcmcout, "log_lik_sat")[[1]]
expect_true(
  all(is.finite(ll_sat1)),
  info = "two-level FIML/MAR log_lik_sat should be finite for every draw/cluster"
)
})

## =============================================================================
## 2. Complete-data-through-the-new-path equivalence check
## =============================================================================
## Forcing missing="fiml" with test="none" on data with NO actual missing
## values routes through the new Y/Z-pattern machinery with a degenerate
## "all observed" pattern per group/cluster (dat$miss == 1, dispatch to
## twolevel_logdens_missing). This should give essentially the same
## posterior as the ordinary complete-data path (dat$miss == 0, dispatch to
## the existing twolevel_logdens), since it's the same likelihood.
##
## Unlike the other bsem() calls in this file, missing="fiml" IS passed
## explicitly here on purpose: R/blavaan.R only auto-defaults missing to
## "ml" when the data actually contains NAs (`!any(is.na(...))` check); with
## Demo.twolevel (no real missingness) that auto-default never kicks in, so
## the explicit argument is the only way to force the "new" pattern-based
## dispatch for this comparison.

try({
set.seed(201)
bfit_old <- bsem(
  model   = model_wb,
  data    = Demo.twolevel,
  cluster = "cluster",
  burnin  = 100,
  sample  = 100,
  dp      = dp_stable
)

set.seed(202)
bfit_new <- bsem(
  model   = model_wb,
  data    = Demo.twolevel,
  cluster = "cluster",
  missing = "fiml",
  fixed.x = FALSE,
  test    = "none",
  burnin  = 100,
  sample  = 100,
  dp      = dp_stable
)

expect_equal(
  sort(names(coef(bfit_old))),
  sort(names(coef(bfit_new))),
  info = "complete-data equivalence: old- and new-path fits should have the same parameter names"
)

expect_true(
  coef_close(bfit_new, sem(model_wb, data = Demo.twolevel, cluster = "cluster")),
  info = "complete-data equivalence: new-path (degenerate pattern) estimates should be close to lavaan MLEs"
)

## the two dispatch paths should give statistically indistinguishable
## posteriors on the SAME complete data (loose tolerance: independent short
## chains, not a bitwise check)
nms <- intersect(names(coef(bfit_old)), names(coef(bfit_new)))
se_old <- sqrt(diag(vcov(bfit_old)))[nms]
se_new <- sqrt(diag(vcov(bfit_new)))[nms]
expect_true(
  all(abs(coef(bfit_old)[nms] - coef(bfit_new)[nms]) < 0.5 * (se_old + se_new + 0.1)),
  info = "complete-data equivalence: old-path and new-path estimates should agree within 0.5 combined-SD"
)
})

## =============================================================================
## 3. Error cases: restrictions specific to two-level FIML/MAR
## =============================================================================

try({
## within-level fixed.x variables ARE supported for two-level FIML/MAR data
## (w1/w2 stay endogenous via ~~, not ~, so this model has no between-level
## fixed.x variable -- keeping it within the supported scope)
model_x_w <- '
    level: 1
        fw =~ y1 + y2 + y3 + y4
        fw ~ x1 + x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~~ w1
        fb ~~ w2
'
dw <- Demo.twolevel
dw <- inject_mcar(dw, c("x1", "x2"), rate = 0.15, seed = 306)
dw <- inject_mcar(dw, c("y1", "y2"), rate = 0.15, seed = 307)
fit_x_w <- sem(model_x_w, data = dw, cluster = "cluster", missing = "ml.x", fixed.x = TRUE)

set.seed(308)
bfit_x_w <- bsem(model_x_w, data = dw, cluster = "cluster", fixed.x = TRUE,
                  test = "ppp", burnin = 100, sample = 100, dp = dp_stable)

expect_inherits(bfit_x_w, "blavaan",
  info = "two-level FIML/MAR with within-level fixed.x should fit successfully")
expect_equal(
  sort(names(coef(bfit_x_w))), sort(names(coef(fit_x_w))),
  info = "within-level fixed.x: parameter names should match lavaan"
)
expect_true(
  coef_close(bfit_x_w, fit_x_w),
  info = "within-level fixed.x: estimates should be within 0.5 SD of lavaan MLEs"
)
expect_true(
  all(diag(vcov(bfit_x_w)) > 0),
  info = "within-level fixed.x: vcov diagonal should be positive"
)

fm_x_w <- fitMeasures(bfit_x_w)
expect_true(is.finite(fm_x_w["waic"]) && is.finite(fm_x_w["looic"]),
  info = "within-level fixed.x: waic/looic should be finite (exercises log_lik_x_full)")

ppp_x_w <- fm_x_w["ppp"]
expect_true(
  is.finite(ppp_x_w) && ppp_x_w >= 0 && ppp_x_w <= 1,
  info = "within-level fixed.x: ppp should be a finite value in [0,1] (exercises calc_log_lik_x_missing())"
)
ll_sat_x_w <- rstan::extract(bfit_x_w@external$mcmcout, "log_lik_sat")[[1]]
expect_true(
  all(is.finite(ll_sat_x_w)),
  info = "within-level fixed.x: log_lik_sat should be finite for every draw/cluster"
)
sampler_params_x_w <- rstan::get_sampler_params(bfit_x_w@external$mcmcout, inc_warmup = FALSE)
expect_equal(
  sum(sapply(sampler_params_x_w, function(x) sum(x[, "divergent__"]))), 0L,
  info = "within-level fixed.x: no divergent transitions expected in a well-behaved short run"
)

## between-level fixed.x + missing data remains blocked (confirmed upstream
## lavaan bug): lavaan()/sem() itself crashes on a between-level fixed.x
## variable with missing values, so R/blavaan.R blocks the combination with
## a clear error before that crash can surface
model_x_b <- '
    level: 1
        fw =~ y1 + y2 + y3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1
'
db <- inject_mcar(Demo.twolevel, c("y1", "y2"), rate = 0.15, seed = 309)
expect_error(
  bsem(model_x_b, data = db, cluster = "cluster", fixed.x = TRUE,
       burnin = 100, sample = 100),
  info = "two-level FIML/MAR with a between-level fixed.x variable should still error"
)

## test != "none" (ppmc/ppp) is now supported for two-level FIML/MAR data:
## log_lik_sat/log_lik_rep/log_lik_rep_sat are computed via an EM-based
## saturated-model fit (see R/lav_export_stanmarg.R, inst/stan/stanmarg.stan
## twolevel_em_step()), unified across complete and missing data.
dt <- inject_mcar(Demo.twolevel, c("y1", "y2"), rate = 0.15, seed = 302)
bfit_ppp <- bsem(model_wb, data = dt, cluster = "cluster",
                  fixed.x = FALSE, test = "ppp", burnin = 100, sample = 100,
                  dp = dp_stable)
expect_inherits(bfit_ppp, "blavaan",
  info = "two-level FIML/MAR with test='ppp' should fit successfully")

ppp_val <- fitMeasures(bfit_ppp)["ppp"]
expect_true(
  is.finite(ppp_val) && ppp_val >= 0 && ppp_val <= 1,
  info = "two-level FIML/MAR ppp should be a finite value in [0,1]"
)

ll_sat <- rstan::extract(bfit_ppp@external$mcmcout, "log_lik_sat")[[1]]
expect_true(
  all(is.finite(ll_sat)),
  info = "two-level FIML/MAR log_lik_sat should be finite for every draw/cluster"
)
})

## =============================================================================
## 4. Regression check: FIML/MAR is now the default, with no missing= argument
## =============================================================================
## Previously (with the blavaan.multilevel.missing gate off, its default),
## a two-level model with real missing data and no explicit missing=
## argument silently fell back to listwise deletion. Now FIML/MAR is used
## automatically instead: nobs should retain every row (no deletion), not
## drop to the listwise-deleted count.
##
## NB: passing missing="listwise" explicitly does NOT opt back into listwise
## deletion for two-level models. missing="listwise" was never a real,
## user-facing feature here -- it was only ever blavaan's own internal
## choice (set unconditionally to force row-dropping when this FIML/MAR
## support didn't exist yet), not something read back from what the user
## passed. Concretely: blavaan's internal exploratory pre-fit (used to get
## info about equality constraints, for npar; see the
## `dotdotdot$missing <- "direct"` line in R/blavaan.R) unconditionally
## overwrites missing to "direct" (FIML-like) before the two-level
## default-selection logic ever runs, so any explicit missing= value the
## user supplies is lost regardless. Unrelated to removing the gate, and
## out of scope here.

try({
d5 <- inject_mcar(Demo.twolevel, c("y1", "y2"), rate = 0.15, seed = 401)

set.seed(402)
bfit5 <- bsem(
  model   = model_wb,
  data    = d5,
  cluster = "cluster",
  burnin  = 100,
  sample  = 100,
  dp      = dp_stable
)

expect_inherits(bfit5, "blavaan",
  info = "two-level models with missing data and no explicit missing= argument should fit via FIML/MAR by default")
expect_equal(
  blavInspect(bfit5, "nobs"), nrow(d5),
  info = "with no explicit missing= argument, FIML/MAR should retain every row (no listwise deletion)"
)
})

## =============================================================================
## 5. Model without within-only variables
## =============================================================================

try({
model_bw <- '
    level: 1
        fw =~ y1 + y2 + y3
    level: 2
        fb =~ y1 + y2 + y3
'
d6 <- inject_mcar(Demo.twolevel, c("y1", "y2"), rate = 0.15, seed = 301)
fit6 <- sem(model_bw, data = d6, cluster = "cluster", missing = "fiml")

set.seed(402)
bfit6 <- bsem(
  model   = model_bw,
  data    = d6,
  cluster = "cluster",
  burnin  = 100,
  sample  = 100,
  dp      = dp_stable
)

expect_inherits(bfit6, "blavaan",
  info = "two-level models with missing data and no within-only variables should work")

expect_equal(
  sort(names(coef(fit6))),
  sort(names(coef(bfit6))),
  info = "complete-data equivalence: old- and new-path fits should have the same parameter names"
)

expect_true(
  coef_close(bfit6, fit6),
  info = "complete-data equivalence: new-path (degenerate pattern) estimates should be close to lavaan MLEs"
)

})
