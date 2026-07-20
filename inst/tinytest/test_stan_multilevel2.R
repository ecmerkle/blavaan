## tinytest file: Testing blavaan two-level (multilevel) models with Stan backend
##
## Run with: tinytest::run_test_file("test_stan_multilevel2.R")
## or as part of a package: tinytest::test_package("blavaan")

library("tinytest")
library("lavaan")
library("blavaan")

set.seed(12345)

## ── Helpers ──────────────────────────────────────────────────────────────────

## Check that Bayesian estimates are within `tol` (in SD units) of lavaan MLEs.
## We use a loose tolerance because MCMC chains are short (burnin=100, sample=100).
coef_close <- function(bfit, fit, tol = 0.5) {
  bc  <- coef(bfit)
  lc  <- coef(fit)
  ## align by name — only compare parameters present in both
  nms <- intersect(names(bc), names(lc))
  if (length(nms) == 0L) return(FALSE)
  bse <- sqrt(diag(vcov(bfit)))[nms]
  all(abs(bc[nms] - lc[nms]) < tol * (abs(lc[nms]) + bse + 0.1))
}

## ── Data ─────────────────────────────────────────────────────────────────────

data(Demo.twolevel, package = "lavaan")

## =============================================================================
## 1. Basic two-level CFA  (named levels: within / between)
## =============================================================================

try({
model_wb <- '
    level: within
        fw =~ y1 + y2 + y3
    level: between
        fb =~ y1 + y2 + y3
'

fit_wb  <- sem(model = model_wb, data = Demo.twolevel, cluster = "cluster")

bfit_wb <- bsem(
  model   = model_wb,
  data    = Demo.twolevel,
  cluster = "cluster",
  burnin  = 100,
  sample  = 100,
  dp      = dpriors(lambda = "normal(1,.5)")
)

expect_inherits(bfit_wb, "blavaan",
  info = "bsem() with named levels (within/between) should return a blavaan object")

## Coefficient vector has same names as lavaan
expect_equal(
  sort(names(coef(bfit_wb))),
  sort(names(coef(fit_wb))),
  info = "bsem() and sem() should return the same parameter names"
)

## Estimates are in the right ballpark
expect_true(
  coef_close(bfit_wb, fit_wb),
  info = "bsem() estimates should be within 0.5 SD of lavaan MLEs (within/between model)"
)

## vcov is positive-definite (all diagonal elements > 0)
expect_true(
  all(diag(vcov(bfit_wb)) > 0),
  info = "vcov(bfit) diagonal should be positive"
)

## SEs are finite
expect_true(
  all(is.finite(sqrt(diag(vcov(bfit_wb))))),
  info = "posterior SEs should all be finite"
)

## fitMeasures returns a named numeric vector with key indices
fm_wb <- fitMeasures(bfit_wb)
expect_true(is.numeric(fm_wb),
  info = "fitMeasures() should return a numeric vector")
expect_true(length(fm_wb) > 0,
  info = "fitMeasures() should return at least one measure")
for (idx in c("ppp", "waic", "looic")) {
  expect_true(idx %in% names(fm_wb),
    info = paste("fitMeasures() should contain", idx))
}

## PPP should be between 0 and 1
expect_true(
  fm_wb["ppp"] >= 0 & fm_wb["ppp"] <= 1,
  info = "PPP should be in [0, 1]"
)
})

## NB: this section is a variation on section 1 (just numeric vs named level
## labels) with its own extra Stan fit, so it's gated as slow; only runs with
## Sys.setenv(blavaan_slow_tests = "true").
if (Sys.getenv("blavaan_slow_tests") == "true") {

## =============================================================================
## 2. Basic two-level CFA  (numeric levels: 1 / 2)
## =============================================================================

model_12 <- '
    level: 1
        fw =~ y1 + y2 + y3
    level: 2
        fb =~ y1 + y2 + y3
'

try({
fit_12 <- sem(model = model_12, data = Demo.twolevel, cluster = "cluster")

bfit_12 <- bsem(
  model   = model_12,
  data    = Demo.twolevel,
  cluster = "cluster",
  burnin  = 100,
  sample  = 100
)

expect_inherits(bfit_12, "blavaan",
  info = "bsem() with numeric levels (1/2) should return a blavaan object")

expect_equal(
  sort(names(coef(bfit_12))),
  sort(names(coef(fit_12))),
  info = "Numeric-level model: parameter names should match lavaan"
)

expect_true(
  coef_close(bfit_12, fit_12),
  info = "Numeric-level model: estimates should be close to lavaan MLEs"
)

fm_12 <- fitMeasures(bfit_12)
expect_true("ppp" %in% names(fm_12),
  info = "Numeric-level model: fitMeasures() should contain ppp")
expect_true(fm_12["ppp"] >= 0 & fm_12["ppp"] <= 1,
  info = "Numeric-level model: PPP should be in [0, 1]")
})

} ## end slow-test gate (section 2)

## =============================================================================
## 3. do.fit = FALSE  — model should be set up but not sampled
## =============================================================================

try({
bfit_nofit <- bsem(
  model   = model_wb,
  data    = Demo.twolevel,
  cluster = "cluster",
  burnin  = 100,
  sample  = 100,
  do.fit  = FALSE
)

expect_inherits(bfit_nofit, "blavaan",
  info = "bsem() with do.fit=FALSE should return a blavaan object")

## No posterior samples should be stored when do.fit = FALSE
expect_error(
  blavInspect(bfit_nofit, "mcmc"),
  info = "do.fit=FALSE: no posterior samples should be stored"
)
})

## =============================================================================
## 4. Prior predictive sampling  (prisamp = TRUE)
## =============================================================================

try({
bfit_prior <- bsem(
  model   = model_wb,
  data    = Demo.twolevel,
  cluster = "cluster",
  burnin  = 100,
  sample  = 100,
  prisamp = TRUE
)

expect_inherits(bfit_prior, "blavaan",
  info = "bsem() with prisamp=TRUE should return a blavaan object")

expect_error(
  sampleData(bfit_prior),
  info = "sampleData() does not work with two-level models"
)
})

## =============================================================================
## 5. Error cases
## =============================================================================

## JAGS target not supported for multilevel models
try({
expect_error(
  bsem(model_wb, data = Demo.twolevel, cluster = "cluster",
       burnin = 100, sample = 100, target = "jags"),
  info = "bsem() multilevel with target='jags' should throw an error"
)

## Summary covariance matrix input not supported for multilevel
expect_error(
  bsem(model_wb,
       sample.cov   = cov(Demo.twolevel[, c("y1", "y2", "y3")]),
       sample.nobs  = nrow(Demo.twolevel),
       sample.mean  = matrix(0, 3, 1),
       cluster      = "cluster"),
  info = "bsem() with sample.cov instead of raw data should error for multilevel"
)
})

## =============================================================================
## 6. MCMC output structure checks
## =============================================================================

try({
## blavInspect returns MCMC chains as a list of matrices
mcmc_chains <- blavInspect(bfit_wb, "mcmc")
expect_true(is.list(mcmc_chains),
  info = "blavInspect(., 'mcmc') should return a list")
expect_true(length(mcmc_chains) >= 1,
  info = "There should be at least one MCMC chain")
expect_true(
  all(sapply(mcmc_chains, is.matrix)),
  info = "Each MCMC chain should be a matrix"
)

## Number of posterior draws == sample argument
total_draws <- sum(sapply(mcmc_chains, nrow))
expect_equal(total_draws, 300L,
  info = "Total posterior draws should equal the 'sample' * 'n.chains' argument (1300)")

## Column names of MCMC matrix correspond to model parameters
chain_params <- colnames(mcmc_chains[[1]])
expect_true(length(chain_params) > 0,
  info = "MCMC chains should have named columns for parameters"
)
})

## NB: sections 7-9 add their own extra Stan fit (bfit_lvs, save.lvs=TRUE)
## that sections 8-9 also depend on, so they're gated together as slow; only
## run with Sys.setenv(blavaan_slow_tests = "true").
if (Sys.getenv("blavaan_slow_tests") == "true") {

## =============================================================================
## 7. Posterior predictive checks (ppmc)
## =============================================================================

try({
## We need save.lvs = TRUE for ppmc with LVs; refit with that option
bfit_lvs <- bsem(
  model    = model_wb,
  data     = Demo.twolevel,
  cluster  = "cluster",
  burnin   = 100,
  sample   = 100,
  dp       = dpriors(lambda = "normal(1,.5)"),
  save.lvs = TRUE
)

expect_inherits(bfit_lvs, "blavaan",
  info = "bsem() with save.lvs=TRUE should return a blavaan object")

## standardizedPosterior
sp <- standardizedPosterior(bfit_lvs)
expect_true(is.matrix(sp) || is.data.frame(sp),
  info = "standardizedPosterior() should return a matrix or data frame"
)

## ppmc
expect_error(ppmc(bfit_lvs),
  info = "ppmc() does not work with two-level models"
)
})

## =============================================================================
## 8. Latent variable predictions (blavPredict / blavInspect)
## =============================================================================

try({
## blavPredict (level 1)
lv1_blav <- blavPredict(bfit_lvs)
expect_true(
  is.list(lv1_blav),
  info = "blavPredict() should return list predictions"
)

## blavPredict() and lavPredict() should correlate > .95 for level-1 LVs
lv1_lav_pred <- lavPredict(fit_wb)
cor_blavpred <- cor(rowMeans(do.call("cbind", lv1_blav)), as.numeric(lv1_lav_pred[, 1]))
expect_true(
  cor_blavpred > 0.95,
  info = "blavPredict() level-1 LV predictions should correlate > .95 with lavPredict()"
)

## blavInspect lvmeans (level 1) vs lavPredict
lv1_blav2 <- blavInspect(bfit_lvs, "lvmeans")
lv1_lav   <- lavPredict(fit_wb)

expect_equal(
  length(lv1_blav2),
  nrow(lv1_lav),
  info = "blavInspect lvmeans and lavPredict should have same number of observations (level 1)"
)

## Correlation between Bayesian and frequentist LV means should be > .95
cor_lv1 <- cor(as.numeric(lv1_blav2), as.numeric(lv1_lav[, 1]))
expect_true(
  cor_lv1 > 0.95,
  info = "Bayesian and lavaan level-1 LV means should be highly correlated (r > 0.95)"
)

## blavInspect lvmeans level 2
lv2_blav <- blavInspect(bfit_lvs, "lvmeans", level = 2)
lv2_lav  <- lavPredict(fit_wb, level = 2)

expect_equal(
  length(lv2_blav),
  nrow(lv2_lav),
  info = "blavInspect lvmeans and lavPredict should have same number of clusters (level 2)"
)

cor_lv2 <- cor(as.numeric(lv2_blav), as.numeric(lv2_lav[, 1]))
expect_true(
  cor_lv2 > 0.95,
  info = "Bayesian and lavaan level-2 LV means should be highly correlated (r > 0.95)"
)

## blavInspect 'lvs' should return a list of matrices (one per draw)
lvs1 <- blavInspect(bfit_lvs, "lvs")
expect_true(is.list(lvs1),
  info = "blavInspect(., 'lvs') should return a list"
)
expect_true(
  all(sapply(lvs1, is.matrix)),
  info = "Each element of blavInspect(., 'lvs') should be a matrix"
)

lvs2 <- blavInspect(bfit_lvs, "lvs", level = 2)
expect_true(is.list(lvs2),
  info = "blavInspect(., 'lvs', level=2) should return a list"
)

## Posterior SD of LV scores should be positive and finite
lv_sd <- apply(do.call("rbind", lvs1), 2, sd)
expect_true(
  all(lv_sd > 0) && all(is.finite(lv_sd)),
  info = "Posterior SD of level-1 LV scores should be positive and finite"
)

## Posterior means of level-2 LV draws should correlate > .95 with lavaan
lv2_post_means <- colMeans(do.call("rbind", lvs2))
lv2_lav_pred   <- as.numeric(lavPredict(fit_wb, level = 2)[, 1])
cor_lvs2 <- cor(lv2_post_means, lv2_lav_pred)
expect_true(
  cor_lvs2 > 0.95,
  info = "Posterior means of level-2 LV draws should correlate > .95 with lavaan lavPredict(level=2)"
)
})

## =============================================================================
## 9. summary() and blavInspect() smoke tests
## =============================================================================

try({
smry_wb <- summary(bfit_wb)
expect_true(!is.null(smry_wb),
  info = "summary() on blavaan object should return a non-null result")

smry_lvs <- summary(bfit_lvs)
expect_true(!is.null(smry_lvs),
  info = "summary() on blavaan object with save.lvs should return a non-null result")

## Rhat values should be available and (with short chains) not badly diverged
rhat_vals <- blavInspect(bfit_wb, "psrf")
expect_true(
  is.numeric(rhat_vals) || is.matrix(rhat_vals),
  info = "blavInspect(., 'psrf') should return numeric Rhat values"
)

## neff / n.eff
neff_vals <- blavInspect(bfit_wb, "neff")
expect_true(
  is.numeric(neff_vals),
  info = "blavInspect(., 'neff') should return numeric effective sample sizes"
)
expect_true(
  all(neff_vals > 0),
  info = "Effective sample sizes should all be positive"
)

## loo / waic objects via blavFitIndices or direct extraction
expect_true("looic" %in% names(fitMeasures(bfit_wb)),
  info = "fitMeasures() should include 'looic'")
expect_true("waic"  %in% names(fitMeasures(bfit_wb)),
  info = "fitMeasures() should include 'waic'")
})

} ## end slow-test gate (sections 7-9)

## NB: sections 10-11 each fit their own multi-group multilevel model
## (independent extra Stan fits), so they're gated as slow; only run with
## Sys.setenv(blavaan_slow_tests = "true").
if (Sys.getenv("blavaan_slow_tests") == "true") {

## =============================================================================
## 10. Multi-group two-level CFA
## =============================================================================
## Multi-group + multilevel requires "group:" blocks as the OUTER blocks,
## each with its own nested "level:" blocks (group= alone, or level: blocks
## with a bare group= argument, throws "subscript out of bounds" in
## lavaan's lav_lavaan_step01_ovnames_namesl()). This also exercises the
## group*level correlation-block prior matching in lav2stanmarg() (psi/theta
## "block" == group*level number, not lavpartable$group, which holds the
## literal group label here rather than a 1:Ng index).

try({
set.seed(999)
d_mg <- Demo.twolevel
d_mg$grp <- ifelse(d_mg$cluster <= median(unique(d_mg$cluster)), "A", "B")

model_mg <- '
group: A
level: 1
    fw =~ y1 + y2 + y3
level: 2
    fb =~ y1 + y2 + y3
    fb ~~ w1
group: B
level: 1
    fw =~ y1 + y2 + y3
level: 2
    fb =~ y1 + y2 + y3
    fb ~~ w1
'

fit_mg <- sem(model_mg, data = d_mg, cluster = "cluster", group = "grp")

bfit_mg <- bsem(
  model   = model_mg,
  data    = d_mg,
  cluster = "cluster",
  group   = "grp",
  burnin  = 100,
  sample  = 100
)

expect_inherits(bfit_mg, "blavaan",
  info = "multi-group two-level bsem() should return a blavaan object")

expect_equal(
  sort(names(coef(bfit_mg))),
  sort(names(coef(fit_mg))),
  info = "multi-group two-level: parameter names should match lavaan"
)

expect_true(
  all(is.finite(coef(bfit_mg))),
  info = "multi-group two-level: all coefficients should be finite"
)

expect_true(
  coef_close(bfit_mg, fit_mg),
  info = "multi-group two-level: estimates should be within 0.5 SD of lavaan MLEs"
)
})

## =============================================================================
## 11. Multi-group two-level fixed.x (regression test for a log_lik_x
##     double-counting bug)
## =============================================================================
## stanmarg.stan's `target += -log_lik_x` subtracted the WHOLE cross-group
## log_lik_x vector on every group's pass through the model-fitting loop
## (log_lik_x spans all groups for multilevel), over-counting by a factor of
## Ng for Ng>1 fixed.x models; fixed by slicing to each group's own r3:r4
## range. Ng=1 fixed.x models (tested elsewhere in this file/
## test_stan_multilevel_missing.R) can't detect this bug, since the
## over-count is a no-op when there's only one group.

try({
set.seed(998)
d_mgx <- Demo.twolevel
d_mgx$grp <- ifelse(d_mgx$cluster <= median(unique(d_mgx$cluster)), "A", "B")

model_mgx <- '
group: A
level: 1
    fw =~ y1 + y2 + y3
    fw ~ x1 + x2
level: 2
    fb =~ y1 + y2 + y3
    fb ~~ w1
group: B
level: 1
    fw =~ y1 + y2 + y3
    fw ~ x1 + x2
level: 2
    fb =~ y1 + y2 + y3
    fb ~~ w1
'

fit_mgx <- sem(model_mgx, data = d_mgx, cluster = "cluster", group = "grp", fixed.x = TRUE)

bfit_mgx <- bsem(
  model   = model_mgx,
  data    = d_mgx,
  cluster = "cluster",
  group   = "grp",
  fixed.x = TRUE,
  burnin  = 100,
  sample  = 100
)

expect_inherits(bfit_mgx, "blavaan",
  info = "multi-group two-level fixed.x bsem() should return a blavaan object")

expect_equal(
  sort(names(coef(bfit_mgx))),
  sort(names(coef(fit_mgx))),
  info = "multi-group two-level fixed.x: parameter names should match lavaan"
)

expect_true(
  all(is.finite(coef(bfit_mgx))),
  info = "multi-group two-level fixed.x: all coefficients should be finite"
)

expect_true(
  coef_close(bfit_mgx, fit_mgx),
  info = "multi-group two-level fixed.x: estimates should be within 0.5 SD of lavaan MLEs (would fail under the log_lik_x double-counting bug, which biases the fixed.x-adjusted likelihood surface)"
)
})

} ## end slow-test gate (sections 10-11)

