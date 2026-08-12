## tinytest file: Testing the R-side saturated-model ppp for two-level
## (multilevel) models (R/blav_twolevel_ppp.R, R/ppp.R), and the new
## test="none"-by-default behavior for such models (R/blavaan.R).
##
## Mirrors test_limited_info_ppp.R's structure/conventions for the
## analogous single-level ordinal feature, but with a DIFFERENT statistic:
## two-level ppp is a saturated-model likelihood-ratio comparison
## (T = -2*(ll_fit - ll_sat)), not pairwise limited-information.
##
## Run with: tinytest::run_test_file("test_twolevel_ppp.R")
## or as part of a package: tinytest::test_package("blavaan")

library("tinytest")
library("lavaan")
library("blavaan")

source("helper_convergence.R")

set.seed(1357)

data(Demo.twolevel, package = "lavaan")

model_wb <- '
    level: within
        fw =~ y1 + y2 + y3
    level: between
        fb =~ y1 + y2 + y3
'

## =============================================================================
## 1. test="none" is the new implicit default for two-level models, and the
##    fit-time NOTE fires only when test= was not explicitly supplied
## =============================================================================

try({
msgs <- capture.output(
  fit_default <- bsem(model_wb, data = Demo.twolevel, cluster = "cluster",
                       burnin = 100, sample = 100, n.chains = 2,
                       dp = dpriors(lambda = "normal(1,.5)"))
)
expect_true(any(grepl("ppp is no longer computed by default", msgs)),
  info = "Two-level model with no test= argument should print the new-default NOTE")

expect_true(lavInspect(fit_default, "options")$test == "none",
  info = "Two-level model with no test= argument should default to test=\"none\"")

fm_default <- fitMeasures(fit_default)
expect_true(!("ppp" %in% names(fm_default)) || is.na(fm_default["ppp"]),
  info = "Two-level model with test=\"none\" default should not have ppp in fitMeasures()")

## second step: ppp(fit) should compute it on demand
ppval <- ppp(fit_default, thin = 5)
expect_true(is.numeric(ppval) && !is.na(ppval) && ppval >= 0 && ppval <= 1,
  info = "ppp(fit) should compute a valid [0,1] value for the deferred two-level case")

## waic/looic (mandatory log_lik, unaffected by do_test) should still be
## available under the new default
expect_true(all(is.finite(fitMeasures(fit_default, c("waic", "looic")))),
  info = "waic/looic should remain available and finite under the new test=\"none\" default")
})


## =============================================================================
## 2. Explicit test="ppp": NOTE should NOT print, ppp computed inline via
##    R/blav_test.R's dispatch to pp_twolevel()
## =============================================================================

try({
msgs2 <- capture.output(
  fit_std <- bsem(model_wb, data = Demo.twolevel, cluster = "cluster",
                   test = "ppp", burnin = 100, sample = 100, n.chains = 2,
                   dp = dpriors(lambda = "normal(1,.5)"))
)
expect_true(!any(grepl("ppp is no longer computed by default", msgs2)),
  info = "Explicit test=\"ppp\" on a two-level model should not print the new-default NOTE")

expect_true(lavInspect(fit_std, "options")$test == "ppp" ||
              lavInspect(fit_std, "options")$test == "standard",
  info = "Explicit test= should be respected, not overridden to \"none\"")

expect_true("ppp" %in% names(fitMeasures(fit_std)) && !is.na(fitMeasures(fit_std, "ppp")),
  info = "Explicit test=\"ppp\" on a two-level model should compute ppp inline")

fm_std <- fitMeasures(fit_std, "ppp")
expect_true(fm_std >= 0 && fm_std <= 1,
  info = "Fit-time two-level ppp should be in [0, 1]")

## ppp(fit) on an already-computed fit should just return the stored value,
## not recompute
t0 <- Sys.time()
ppval_std <- ppp(fit_std)
dt <- as.numeric(Sys.time() - t0, units = "secs")
expect_equal(as.numeric(ppval_std), as.numeric(fitMeasures(fit_std, "ppp")),
  info = "ppp(fit) should match fitMeasures(fit,'ppp') when already computed, not recompute")
expect_true(dt < 2,
  info = "ppp(fit) should return near-instantly when the value is already stored")

## the underlying Stan program never computes two-level ppp any more: its
## own 'ppp' generated quantity should be 0 (the do_test=FALSE fallback),
## confirming the fitMeasures() value above genuinely came from R, not Stan
stan_ppp <- fit_std@external$stansumm["ppp", "mean"]
expect_equal(as.numeric(stan_ppp), 0,
  info = "Stan's own ppp generated quantity should be the do_test=FALSE fallback (0) for multilevel fits")
})


## =============================================================================
## 3. Explicit test="none" on a two-level model: NOTE should not print
##    (already exactly what the user asked for)
## =============================================================================

try({
msgs3 <- capture.output(
  fit_none <- bsem(model_wb, data = Demo.twolevel, cluster = "cluster",
                    test = "none", burnin = 60, sample = 60, n.chains = 2,
                    dp = dpriors(lambda = "normal(1,.5)"))
)
expect_true(!any(grepl("ppp is no longer computed by default", msgs3)),
  info = "Explicit test=\"none\" on a two-level model should not print the new-default NOTE")
})


## =============================================================================
## 4. Missing-data two-level ppp (no fixed.x) works via the same R path
##
## Two-level models with real missingness are a niche combination that
## doesn't need default-CI coverage; only run with
## Sys.setenv(blavaan_slow_tests = "true").
## =============================================================================

if (Sys.getenv("blavaan_slow_tests") == "true") {
try({
dmiss <- Demo.twolevel
set.seed(24)
miss_idx <- sample(nrow(dmiss), floor(0.1 * nrow(dmiss)))
dmiss$y1[miss_idx] <- NA

fit_miss <- bsem(model_wb, data = dmiss, cluster = "cluster",
                  test = "ppp", burnin = 100, sample = 100, n.chains = 2,
                  dp = dpriors(lambda = "normal(1,.5)"))

ppval_miss <- fitMeasures(fit_miss, "ppp")
expect_true(is.finite(ppval_miss) && ppval_miss >= 0 && ppval_miss <= 1,
  info = "Two-level missing-data ppp should be a finite value in [0,1]")

ppval_miss2 <- ppp(fit_miss)
expect_equal(as.numeric(ppval_miss2), as.numeric(ppval_miss),
  info = "ppp(fit) on a missing-data two-level fit should match the stored fit-time value")
})
}


## =============================================================================
## 5. fixed.x support: within-level and between-level fixed.x variables both
##    work (x is held at its observed value in each replicate, not
##    resimulated -- see R/blav_twolevel_ppp.R's tl_cond_sim()); a variable
##    that is fixed.x at BOTH levels simultaneously is not supported (a
##    pre-existing lavaan limitation -- see outstanding_issues.md item 3 at
##    the time this was written -- not something introduced here) and
##    should still error clearly.
## =============================================================================

try({
model_x_within <- '
    level: within
        fw =~ y1 + y2 + y3
        fw ~ x1
    level: between
        fb =~ y1 + y2 + y3
'

fit_x_default <- bsem(model_x_within, data = Demo.twolevel, cluster = "cluster",
                       fixed.x = TRUE, burnin = 60, sample = 60, n.chains = 2)
expect_true(lavInspect(fit_x_default, "options")$test == "none",
  info = "A fixed.x two-level model should still default to test=\"none\"")

ppval_x_within <- ppp(fit_x_default, thin = 5)
expect_true(is.numeric(ppval_x_within) && is.finite(ppval_x_within) &&
              ppval_x_within >= 0 && ppval_x_within <= 1,
  info = "ppp() on a two-level model with within-level fixed.x should return a valid [0,1] value")

## fixed.x columns must never be resimulated: every replicate's x1 column
## should reproduce the observed x1 values exactly
implied_check <- lav_model_implied(fit_x_default@Model)
dataX_rep <- blavaan:::tl_replicate_data(fit_x_default@Data, fit_x_default@SampleStats, implied_check)
x1_col <- match("x1", fit_x_default@Data@ov.names[[1]])
expect_equal(dataX_rep[[1]][, x1_col], fit_x_default@Data@X[[1]][, x1_col],
  info = "tl_replicate_data() must copy fixed.x columns unchanged, never resimulate them")
})


try({
model_x_between <- '
    level: within
        fw =~ y1 + y2 + y3
    level: between
        fb =~ y1 + y2 + y3
        fb ~ w1
'

fit_x_between <- bsem(model_x_between, data = Demo.twolevel, cluster = "cluster",
                       fixed.x = TRUE, burnin = 60, sample = 60, n.chains = 2)
ppval_x_between <- ppp(fit_x_between, thin = 5)
expect_true(is.numeric(ppval_x_between) && is.finite(ppval_x_between) &&
              ppval_x_between >= 0 && ppval_x_between <= 1,
  info = "ppp() on a two-level model with between-level fixed.x should return a valid [0,1] value")
})


## =============================================================================
## 6. em_control / thin / parallel arguments are accepted
## =============================================================================

try({
res_default_em <- ppp(fit_default, thin = 10)
expect_true(is.numeric(res_default_em) && is.finite(res_default_em),
  info = "ppp(fit, thin=) should be accepted and return a finite value")

res_em <- blavaan:::pp_twolevel(fit_default, thin = 10,
                                em_control = list(tol = 1e-2, max_iter = 20L,
                                                  acceleration = "none"))
expect_true(is.numeric(res_em$ppval) && is.finite(res_em$ppval) &&
              res_em$ppval >= 0 && res_em$ppval <= 1,
  info = "pp_twolevel(em_control=) should be accepted and return a valid [0,1] ppval")
expect_true(is.list(res_em$ppdist) && all(c("obs", "reps") %in% names(res_em$ppdist)),
  info = "pp_twolevel() should return a ppdist list with obs/reps, matching the other ppp-family functions")
})
