## tinytest file: Testing the R-side limited-information ppp for ordinal/mixed
## single-level models (R/blav_limited_info.R, R/ppp.R), and the new
## test="none"-by-default behavior for such models (R/blavaan.R).
##
## Run with: tinytest::run_test_file("test_limited_info_ppp.R")
## or as part of a package: tinytest::test_package("blavaan")

library("tinytest")
library("lavaan")
library("blavaan")

source("helper_convergence.R")

set.seed(2468)

mytarg <- "stan"

## make ordinal data from continuous (same helper convention as test_ord2.R)
makeord <- function(Data, vars = NULL, ncat = 3) {
  if (length(vars) == 0) vars <- 1:NCOL(Data)
  if (length(ncat) != 1 && length(ncat) != length(vars)) stop("bad ncat")
  Data <- rbind(Data, NA)
  Data[nrow(Data), vars] <- ncat
  Data[, vars] <- apply(Data[, vars, drop = FALSE], 2, function(x) {
    nc   <- tail(x, 1)
    tmpp <- (1/nc) + runif(1, -.1/nc, .1/nc)
    brks <- c(min(x, na.rm = TRUE) - .1,
              seq(quantile(x, tmpp, na.rm = TRUE),
                  quantile(x, 1 - tmpp, na.rm = TRUE),
                  length.out = (nc - 1)),
              max(x, na.rm = TRUE) + .1)
    cut(x, breaks = brks, labels = FALSE)
  })
  Data[-nrow(Data), ]
}

## =============================================================================
## 1. test="none" is the new implicit default for single-level ordinal models,
##    and the fit-time NOTE fires only when test= was not explicitly supplied
## =============================================================================

try({
pop.model <- ' f1 =~ x1 + x2 + x3 + x4 '
Data <- simulateData(pop.model, sample.nobs = 300)
Data <- makeord(Data, vars = c(1, 2), ncat = 3)
model <- ' f1 =~ x1 + x2 + x3 + x4 '

## default test=: NOTE should print, ppp absent from fitMeasures()
msgs <- capture.output(
  fit_default <- bcfa(model, data = Data, ordered = c("x1", "x2"),
                       burnin = 100, sample = 100, target = mytarg)
)
expect_true(any(grepl("ppp is no longer computed by default", msgs)),
  info = "Ordinal model with no test= argument should print the new-default NOTE")

expect_true(lavInspect(fit_default, "options")$test == "none",
  info = "Ordinal model with no test= argument should default to test=\"none\"")

fm_default <- fitMeasures(fit_default)
expect_true(!("ppp" %in% names(fm_default)) || is.na(fm_default["ppp"]),
  info = "Ordinal model with test=\"none\" default should not have ppp in fitMeasures()")

## second step: ppp(fit) should compute it on demand
ppval <- ppp(fit_default)
expect_true(is.numeric(ppval) && !is.na(ppval) && ppval >= 0 && ppval <= 1,
  info = "ppp(fit) should compute a valid [0,1] value for the deferred ordinal case")


## explicit test="standard": NOTE should NOT print, ppp computed inline
msgs2 <- capture.output(
  fit_std <- bcfa(model, data = Data, ordered = c("x1", "x2"), test = "standard",
                   burnin = 100, sample = 100, target = mytarg)
)
expect_true(!any(grepl("ppp is no longer computed by default", msgs2)),
  info = "Explicit test=\"standard\" on an ordinal model should not print the new-default NOTE")

expect_true(lavInspect(fit_std, "options")$test == "standard",
  info = "Explicit test=\"standard\" should be respected, not overridden to \"none\"")

expect_true("ppp" %in% names(fitMeasures(fit_std)) && !is.na(fitMeasures(fit_std, "ppp")),
  info = "Explicit test=\"standard\" on an ordinal model should compute ppp inline")

## ppp(fit) on an already-computed fit should just return the stored value,
## not recompute -- check it returns quickly and matches fitMeasures()
t0 <- Sys.time()
ppval_std <- ppp(fit_std)
dt <- as.numeric(Sys.time() - t0, units = "secs")
expect_equal(as.numeric(ppval_std), as.numeric(fitMeasures(fit_std, "ppp")),
  info = "ppp(fit) should match fitMeasures(fit,'ppp') when already computed, not recompute")
expect_true(dt < 2,
  info = "ppp(fit) should return near-instantly when the value is already stored")

## explicit test="none" on an ordinal model: NOTE should not print (already
## exactly what the user asked for)
msgs3 <- capture.output(
  fit_none <- bcfa(model, data = Data, ordered = c("x1", "x2"), test = "none",
                    burnin = 100, sample = 100, target = mytarg)
)
expect_true(!any(grepl("ppp is no longer computed by default", msgs3)),
  info = "Explicit test=\"none\" on an ordinal model should not print the new-default NOTE")


## old Stan-computed ppp stays retrievable alongside the new default, when
## test="standard" is used (both values exist and generally differ)
old_val <- fit_std@external$stansumm["ppp", "mean"]
expect_true(is.numeric(old_val) && !is.na(old_val),
  info = "Old Stan-computed ppp should remain present in stansumm when test=\"standard\"")
})


## =============================================================================
## 2. Continuous-only model: default test stays "standard" (unaffected); an
##    explicit test="none" on such a model should still work via ppp(), which
##    must route through postpred(), not pp_limited_info() (regression test
##    for a dispatcher bug caught during implementation)
## =============================================================================

try({
pop.model2 <- ' f1 =~ x1 + x2 + x3 + x4 '
Data2 <- simulateData(pop.model2, sample.nobs = 300)
model2 <- ' f1 =~ x1 + x2 + x3 + x4 '

msgs4 <- capture.output(
  fit_cont <- bcfa(model2, data = Data2, burnin = 100, sample = 100, target = mytarg)
)
expect_true(!any(grepl("ppp is no longer computed by default", msgs4)),
  info = "Continuous-only model should never print the ordinal-only NOTE")
expect_true(lavInspect(fit_cont, "options")$test == "standard",
  info = "Continuous-only model should keep test=\"standard\" as its default")
expect_true("ppp" %in% names(fitMeasures(fit_cont)),
  info = "Continuous-only model should have ppp computed by default, as before")

fit_cont_none <- bcfa(model2, data = Data2, test = "none",
                       burnin = 100, sample = 100, target = mytarg)
ppval_cont <- tryCatch(ppp(fit_cont_none), error = function(e) e)
expect_true(is.numeric(ppval_cont) && length(ppval_cont) >= 1,
  info = "ppp() on a continuous-only model with explicit test=\"none\" should succeed via postpred(), not error")
})


## =============================================================================
## 3. Discrimination: a clearly misspecified ordinal/mixed model should show
##    a lower ppp than a correctly-specified one (the key behavioral check
##    the pre-existing bounds-only assertions never provided)
##
## Fits two full models just to compare relative ppp ordering; only run with
## Sys.setenv(blavaan_slow_tests = "true").
## =============================================================================

if (Sys.getenv("blavaan_slow_tests") == "true") {
try({
pop.model3 <- '
  f1 =~ 0.8*x1 + 0.8*x2 + 0.8*x3 + 0.1*x4 + 0.1*x5 + 0.1*x6
  f2 =~ 0.1*x1 + 0.1*x2 + 0.1*x3 + 0.8*x4 + 0.8*x5 + 0.8*x6
  f1 ~~ 1*f1
  f2 ~~ 1*f2
  f1 ~~ 0.3*f2
  x1 ~~ 0.36*x1
  x2 ~~ 0.36*x2
  x3 ~~ 0.36*x3
  x4 ~~ 0.36*x4
  x5 ~~ 0.36*x5
  x6 ~~ 0.36*x6
'
Data3 <- simulateData(pop.model3, sample.nobs = 400)
Data3 <- makeord(Data3, vars = c(1, 4), ncat = 3)

correct_model <- '
  f1 =~ x1 + x2 + x3
  f2 =~ x4 + x5 + x6
  f1 ~~ f2
'
wrong_model <- '
  f1 =~ x1 + x2 + x3 + x4 + x5 + x6
'

res_correct <- robust_fit(bcfa, correct_model, data = Data3, ordered = c("x1", "x4"),
                          test = "standard", burnin = 200, sample = 200, target = mytarg)
res_wrong <- robust_fit(bcfa, wrong_model, data = Data3, ordered = c("x1", "x4"),
                        test = "standard", burnin = 200, sample = 200, target = mytarg)

if (res_correct$converged && res_wrong$converged) {
  ppp_correct <- fitMeasures(res_correct$fit, "ppp")
  ppp_wrong <- fitMeasures(res_wrong$fit, "ppp")
  expect_true(ppp_wrong <= ppp_correct,
    info = paste0("Misspecified model's ppp (", round(ppp_wrong, 3),
                  ") should not exceed the correctly-specified model's (",
                  round(ppp_correct, 3), ")"))
} else {
  cat("SKIPPED discrimination assertion: fit(s) did not converge\n")
}
})
}


## =============================================================================
## 4. Sparse-cell robustness: small N with a near-empty joint ordinal cell
## =============================================================================

try({
set.seed(13)
pop.model4 <- ' f1 =~ x1 + x2 + x3 '
Data4 <- simulateData(pop.model4, sample.nobs = 60)
Data4 <- makeord(Data4, vars = c(1, 2), ncat = 4)
## force sparsity: collapse most of one category away
Data4$x1[Data4$x1 == 4 & seq_len(nrow(Data4)) %% 3 != 0] <- 3

model4 <- ' f1 =~ x1 + x2 + x3 '
fit_sparse <- bcfa(model4, data = Data4, ordered = c("x1", "x2"),
                    burnin = 100, sample = 100, target = mytarg)

ppval_sparse <- tryCatch(ppp(fit_sparse), error = function(e) e)
expect_true(is.numeric(ppval_sparse) && is.finite(ppval_sparse) &&
              ppval_sparse >= 0 && ppval_sparse <= 1,
  info = "ppp() should return a finite [0,1] value even with sparse ordinal cells")
})


## =============================================================================
## 5. Multi-group ordinal
##
## Only run with Sys.setenv(blavaan_slow_tests = "true").
## =============================================================================

if (Sys.getenv("blavaan_slow_tests") == "true") {
try({
pop.model5 <- ' f1 =~ x1 + x2 + x3 + x4 '
D1 <- simulateData(pop.model5, sample.nobs = 150)
D2 <- simulateData(pop.model5, sample.nobs = 150)
D1$grp <- "a"; D2$grp <- "b"
Data5 <- rbind(D1, D2)
Data5 <- makeord(Data5, vars = c(1, 2), ncat = 3)

model5 <- ' f1 =~ x1 + x2 + x3 + x4 '
fit_mg5 <- bcfa(model5, data = Data5, group = "grp", ordered = c("x1", "x2"),
                burnin = 100, sample = 100, target = mytarg)

ppval_mg <- tryCatch(ppp(fit_mg5), error = function(e) e)
expect_true(is.numeric(ppval_mg) && is.finite(ppval_mg) &&
              ppval_mg >= 0 && ppval_mg <= 1,
  info = "ppp() should return a valid [0,1] value for a multi-group ordinal model")
})
}


## =============================================================================
## 6. Two-level models: pp_limited_info() is NOT the dispatch target
##
## Two-level models now get their own analogous treatment (default
## test="none", on-demand ppp() via R/blav_twolevel_ppp.R's pp_twolevel()),
## but with a DIFFERENT statistic (saturated-model likelihood-ratio, not
## pairwise limited-information) -- see test_twolevel_ppp.R for full
## coverage of that. This section only checks that the .multilevel branch
## in R/ppp.R/R/blav_test.R correctly routes to pp_twolevel(), not
## pp_limited_info() (which was never scoped to handle two-level data).
## =============================================================================

try({
data(Demo.twolevel, package = "lavaan")
## subset to the first 30 clusters for speed. Deliberately continuous-only
## (an ordinal two-level fit hits an unrelated, pre-existing blavaan
## constraint -- missing="pairwise" is rejected for two-level models, and
## is not overridable via bcfa()'s missing= argument for cluster-based fits
## -- not something introduced by this change).
Data6 <- Demo.twolevel[Demo.twolevel$cluster <= 30, ]

model6 <- '
  level: within
    fw =~ y1 + y2 + y3
  level: between
    fb =~ y1 + y2 + y3
'
fit2l <- bcfa(model6, data = Data6, cluster = "cluster",
              burnin = 100, sample = 100, target = mytarg)

expect_true(lavInspect(fit2l, "options")$test == "none",
  info = "Two-level model should default to test=\"none\", same as ordinal (R/blavaan.R)")

ppval_2l <- tryCatch(ppp(fit2l, thin = 5), error = function(e) e)
expect_true(is.numeric(ppval_2l) && is.finite(ppval_2l) &&
              ppval_2l >= 0 && ppval_2l <= 1,
  info = "ppp(fit) on a two-level model should compute via pp_twolevel(), not pp_limited_info()")
})


## =============================================================================
## 7. Mixed continuous+ordinal: all three pair types (oo/oc/cc) present
##
## Only run with Sys.setenv(blavaan_slow_tests = "true").
## =============================================================================

if (Sys.getenv("blavaan_slow_tests") == "true") {
try({
pop.model7 <- ' f1 =~ x1 + x2 + x3 + x4 + x5 '
Data7 <- simulateData(pop.model7, sample.nobs = 300)
Data7 <- makeord(Data7, vars = c(1, 2), ncat = 3)  ## x1,x2 ordinal (oo);
                                                     ## x1/x2-x3/x4/x5 give oc;
                                                     ## x3,x4,x5 give cc

model7 <- ' f1 =~ x1 + x2 + x3 + x4 + x5 '
fit_mixed <- bcfa(model7, data = Data7, ordered = c("x1", "x2"),
                   burnin = 100, sample = 100, target = mytarg)

ppval_mixed <- tryCatch(ppp(fit_mixed), error = function(e) e)
expect_true(is.numeric(ppval_mixed) && is.finite(ppval_mixed) &&
              ppval_mixed >= 0 && ppval_mixed <= 1,
  info = "ppp() should return a valid [0,1] value for a mixed continuous+ordinal model (oo+oc+cc)")
})
}
