## tinytest file: Testing blavaan ordinal (and mixed continuous/ordinal) models
##
## Run with: tinytest::run_test_file("tinytest_ord2.R")
## or as part of a package: tinytest::test_package("blavaan")

## NB these tests are slow! they will only be run with
## Sys.setenv(blavaan_slow_tests = "true")
if (Sys.getenv("blavaan_slow_tests") == "true") {

library("tinytest")
library("lavaan")
library("blavaan")

set.seed(12345)

mytarg <- "stan"

## ── Helpers ──────────────────────────────────────────────────────────────────

## make ordinal data from continuous
makeord <- function(Data, vars = NULL, ncat = 2) {
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

## Check that Bayesian estimates are reasonably close to lavaan MLEs.
coef_close <- function(bfit, fit, tol = 0.75) {
  bc  <- coef(bfit)
  lc  <- coef(fit)
  nms <- intersect(names(bc), names(lc))
  if (length(nms) == 0L) return(FALSE)
  bse <- sqrt(diag(vcov(bfit)))[nms]
  all(abs(bc[nms] - lc[nms]) < tol * (abs(lc[nms]) + bse + 0.1))
}

## =============================================================================
## 1. Mixed continuous/ordinal SEM
## =============================================================================

pop.model <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + a*y2 + b*y3 + c*y4
  dem65 =~ y5 + a*y6 + b*y7 + c*y8
  dem60 ~ 1*ind60
  dem65 ~ 1*ind60 + 1*dem60
  y1 ~~ y5
  y2 ~~ y4 + y6
  y3 ~~ y7
  y4 ~~ y8
  y6 ~~ y8
'

try({
Data <- simulateData(pop.model, sample.nobs = 200)
Data <- makeord(Data, vars = c(3, 7, 11), ncat = 3)
Data$x3[Data$x3 == 3L] <- 2L   # unequal number of categories

model <- '
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + a*y2 + b*y3 + c*y4
  dem65 =~ y5 + a*y6 + b*y7 + c*y8
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
  y1 ~~ y5
  y3 ~~ y7
  y2 ~~ y4
'

fit <- sem(model, data = Data, meanstructure = TRUE,
           ordered = c("x3", "y4", "y8"), parameterization = "theta")
fit@optim$converged <- TRUE

bfit <- bsem(model, data = Data, ordered = c("x3", "y4", "y8"),
             burnin = 100, sample = 100,
             dp = dpriors(lambda = "normal(1,.5)", psi = "gamma(2,2)[sd]"),
             save.lvs = TRUE)

expect_inherits(bfit, "blavaan",
  info = "Mixed cont/ord SEM: bsem() should return a blavaan object")

expect_equal(
  sort(names(coef(bfit))),
  sort(names(coef(fit))),
  info = "Mixed cont/ord SEM: parameter names should match lavaan"
)

expect_true(
  coef_close(bfit, fit),
  info = "Ordinal CFA: estimates should be close to lavaan"
)

expect_true(
  all(diag(vcov(bfit)) > 0),
  info = "Mixed cont/ord SEM: vcov diagonal should be positive"
)

expect_true(
  all(is.finite(sqrt(diag(vcov(bfit))))),
  info = "Mixed cont/ord SEM: posterior SEs should be finite"
)

## lvmeans vs lavPredict: correlation should be > .95
tmp_lv  <- blavPredict(bfit, type = "lvmeans")
tmp2_lv <- lavPredict(fit)

expect_true(
  is.matrix(tmp_lv),
  info = "blavPredict(type='lvmeans') should return numeric output"
)

## correlate each LV column
for (j in seq_len(ncol(tmp2_lv))) {
  r <- cor(as.numeric(tmp_lv[, j]), as.numeric(tmp2_lv[, j]))
  expect_true(r > 0.95,
    info = paste0("blavPredict lvmeans col ", j,
                  " should correlate > .95 with lavPredict (r=", round(r, 3), ")"))
}

## other blavPredict types should return without error and sensible output
tmp_lv2  <- blavPredict(bfit, type = "lv")
tmp_ov   <- blavPredict(bfit, type = "ov")
tmp_yhat <- blavPredict(bfit, type = "yhat")
expect_error(blavPredict(bfit, type = "ymis"))

expect_true(
  is.matrix(tmp_lv2) || is.list(tmp_lv2),
  info = "blavPredict(type='lv') should return a matrix or list"
)
expect_true(
  is.matrix(tmp_ov) || is.list(tmp_ov),
  info = "blavPredict(type='ov') should return a matrix or list"
)
expect_true(
  is.matrix(tmp_yhat) || is.list(tmp_yhat),
  info = "blavPredict(type='yhat') should return a matrix or list"
)
})

## =============================================================================
## 2. Missing data with ordinal indicators
## =============================================================================

try({
set.seed(9618)
mis <- matrix(rbinom(prod(dim(Data)), 1, .9), nrow(Data), ncol(Data))
pd  <- Data * mis
pd[pd == 0] <- NA

fitm <- bsem(model, data = pd, burnin = 150, sample = 100,
             ordered = c("x3", "y4", "y8"),
             dp = dpriors(lambda = "normal(1,1)"), save.lvs = TRUE)

expect_inherits(fitm, "blavaan",
  info = "Ordinal SEM with missing data: bsem() should return a blavaan object")

expect_true(
  all(is.finite(coef(fitm))),
  info = "Ordinal SEM with missing data: all coefficients should be finite"
)

expect_true(
  all(diag(vcov(fitm)) > 0),
  info = "Ordinal SEM with missing data: vcov diagonal should be positive"
)

## lvmeans vs lavPredict: correlation should be > .95
tmp_lv  <- blavPredict(fitm, type = "lvmeans")

expect_true(
  is.matrix(tmp_lv),
  info = "blavPredict(type='lvmeans') should return numeric output"
)

## correlate each LV column
for (j in seq_len(ncol(tmp2_lv))) {
  r <- cor(as.numeric(tmp_lv[, j]), as.numeric(tmp2_lv[, j]))
  expect_true(r > 0.95,
    info = paste0("blavPredict lvmeans col ", j,
                  " should correlate > .95 with lavPredict (r=", round(r, 3), ")"))
}
})


## =============================================================================
## 3. Ordinal CFA
## =============================================================================

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

try({
dcont <- simulateData(HS.model, sample.nobs = 200)
Data  <- makeord(dcont, vars = c(3, 6, 9), ncat = 3)

fit2 <- cfa(HS.model, data = Data, meanstructure = TRUE,
            ordered = c("x3", "x6", "x9"), parameterization = "theta")

bfit2 <- bcfa(HS.model, data = Data, sample = 100, burnin = 100,
              ordered = c("x3", "x6", "x9"),
              dp = dpriors(lambda = "normal(1,1)", tau = "normal(0,1)"))

expect_inherits(bfit2, "blavaan",
  info = "Ordinal CFA: bcfa() should return a blavaan object")

expect_equal(
  sort(names(coef(bfit2))),
  sort(names(coef(fit2))),
  info = "Ordinal CFA: parameter names should match lavaan"
)

expect_true(
  coef_close(bfit2, fit2),
  info = "Ordinal CFA: estimates should be close to lavaan"
)

expect_true(
  all(diag(vcov(bfit2)) > 0),
  info = "Ordinal CFA: vcov diagonal should be positive"
)

fm2 <- fitMeasures(bfit2)
expect_true(is.numeric(fm2) && length(fm2) > 0,
  info = "Ordinal CFA: fitMeasures() should return a non-empty numeric vector")
expect_true(fm2["ppp"] >= 0 & fm2["ppp"] <= 1,
  info = "Ordinal CFA: PPP should be in [0, 1]")
})


## =============================================================================
## 4. Ordinal CFA with std.lv=TRUE
## =============================================================================

try({
Data  <- makeord(dcont, ncat = 5)

fit2_stdlv <- bcfa(HS.model, data = Data, sample = 100, burnin = 100,
                   std.lv = TRUE,
                   dp = dpriors(psi = "gamma(1,1)", lambda = "normal(1,.5)"),
                   ordered = TRUE)

fit2b_stdlv <- cfa(HS.model, data = Data, meanstructure = TRUE,
                   std.lv = TRUE, ordered = TRUE, parameterization = "theta")

expect_inherits(fit2_stdlv, "blavaan",
  info = "Ordinal CFA std.lv: bcfa() should return a blavaan object")

expect_equal(
  sort(names(coef(fit2_stdlv))),
  sort(names(coef(fit2b_stdlv))),
  info = "Ordinal CFA std.lv: parameter names should match lavaan"
)

expect_true(
  coef_close(fit2_stdlv, fit2b_stdlv),
  info = "Ordinal CFA std.lv: parameter estimates should be close to lavaan"
)

expect_true(
  all(diag(vcov(fit2_stdlv)) > 0),
  info = "Ordinal CFA std.lv: vcov diagonal should be positive"
)

fm_stdlv <- fitMeasures(fit2_stdlv)
expect_true(fm_stdlv["ppp"] >= 0 & fm_stdlv["ppp"] <= 1,
  info = "Ordinal CFA std.lv: PPP should be in [0, 1]")
})


## =============================================================================
## 5. allow.empty.cell — multi-group ordinal CFA
## =============================================================================

try({
Data$grp  <- rep(1:2, each = nrow(Data)/2)
Data$x3   <- as.factor(Data$x3)
Data$x3[Data$grp == 2 & (Data$x3 == 2 | Data$x3 == 3)] <- 1

## without allow.empty.cell should error
expect_error(
  bcfa(HS.model, data = Data, sample = 100, burnin = 100,
       std.lv = TRUE, ordered = c("x3", "x6", "x9"), group = "grp"),
  info = "bcfa() with empty ordinal cell and no allow.empty.cell should error"
)

## with allow.empty.cell should succeed
fit_ec <- bcfa(HS.model, data = Data, sample = 100, burnin = 100,
               std.lv = TRUE, ordered = c("x3", "x6", "x9"),
               group = "grp", allow.empty.cell = TRUE)

expect_inherits(fit_ec, "blavaan",
  info = "allow.empty.cell=TRUE: bcfa() should return a blavaan object")

smry_ec <- summary(fit_ec)
expect_true(!is.null(smry_ec),
  info = "allow.empty.cell fit: summary() should return non-null")

fm_ec <- fitMeasures(fit_ec)
expect_true(is.numeric(fm_ec) && length(fm_ec) > 0,
  info = "allow.empty.cell fit: fitMeasures() should return a non-empty numeric vector")
})


## =============================================================================
## 6. Fixed thresholds in ordinal CFA
## =============================================================================

hs1 <- ' visual  =~ x1 + x2 + x3
         visual ~~ prior("gamma(1,2.5)") * visual
         x1 | -1*t1 + 1*t2
         x2 | -1.5*t1 + .5*t2 + .6*t3 + 1.5*t4
         x3 | -.8*t1 + .8*t2 '

try({
Data <- makeord(dcont, ncat = c(3, 5, 3, rep(3, 6)))

fit_th <- cfa(hs1, data = Data, meanstructure = TRUE,
              ordered = TRUE, parameterization = "theta", parser = "old")

bfit_th <- bcfa(hs1, data = Data, sample = 100, burnin = 100, ordered = TRUE,
                dp = dpriors(lambda = "normal(1,1)", tau = "normal(0,1)"))

expect_inherits(bfit_th, "blavaan",
  info = "Fixed thresholds CFA: bcfa() should return a blavaan object")

expect_equal(
  sort(names(coef(bfit_th))),
  sort(names(coef(fit_th))),
  info = "Fixed thresholds CFA: parameter names should match lavaan"
)

expect_true(any(grepl("gamma(1,2.5)", parTable(bfit_th)$prior, fixed = TRUE)))

expect_true(coef_close(bfit_th, fit_th),
            info = "Fixed thresholds CFA: parameters should be close to lavaan")

expect_true(
  all(diag(vcov(bfit_th)) > 0),
  info = "Fixed thresholds CFA: vcov diagonal should be positive"
)

fm_th <- fitMeasures(bfit_th)
expect_true("ppp" %in% names(fm_th),
  info = "Fixed thresholds CFA: fitMeasures() should contain ppp")

expect_true(fm_th["ppp"] >= 0 & fm_th["ppp"] <= 1,
  info = "Fixed thresholds CFA: PPP should be in [0, 1]")
})


## =============================================================================
## 7. Multi-group CFA with group.equal="loadings" + wiggle
## =============================================================================
hs39 <- HolzingerSwineford1939[,7:15]
hs39o <- makeord(hs39, ncat=3)
hs39o$school <- HolzingerSwineford1939$school
HS.model <- ' visual =~ x1 + x2 + x3 '

try({
fit_mg <- bcfa(HS.model, data = hs39o, sample = 100, burnin = 100,
               group = "school", group.equal = c("thresholds", "loadings"),
               wiggle = "thresholds", wiggle.sd = .1,
               ordered = TRUE,
               dp = dpriors(lambda = "normal(1,.5)", tau = "normal(0,.5)"))

fit_mg_lav <- cfa(HS.model, data = hs39o, group = "school",
                  group.equal = c("thresholds", "loadings"),
                  parameterization = "theta", ordered = TRUE)

expect_inherits(fit_mg, "blavaan",
  info = "Multi-group ordinal CFA: bcfa() should return a blavaan object")

expect_equal(
  sort(names(coef(fit_mg))),
  sort(names(coef(fit_mg_lav))),
  info = "Multi-group ordinal CFA: parameter names should match lavaan"
)

expect_true(coef_close(fit_mg, fit_mg_lav),
            info = "Fixed thresholds CFA: parameters should be close to lavaan")

expect_true(
  all(diag(vcov(fit_mg)) > 0),
  info = "Multi-group ordinal CFA: vcov diagonal should be positive"
)

fm_mg <- fitMeasures(fit_mg)
expect_true(fm_mg["ppp"] >= 0 & fm_mg["ppp"] <= 1,
  info = "Multi-group ordinal CFA: PPP should be in [0, 1]")
})


## =============================================================================
## 8. Ordinal growth curve model
## =============================================================================

dg <- makeord(Demo.growth, ncat = 3)

model_growth <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
                  s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4 '

try({
fit6 <- bgrowth(model_growth, data = dg, burnin = 100, sample = 100,
                ordered = TRUE)

fit6b <- growth(model_growth, data = dg, ordered = TRUE,
                parameterization = "theta", meanstructure = TRUE)

expect_inherits(fit6, "blavaan",
  info = "Ordinal growth model: bgrowth() should return a blavaan object")

expect_equal(
  sort(names(coef(fit6))),
  sort(names(coef(fit6b))),
  info = "Ordinal growth model: parameter names should match lavaan"
)

expect_true(coef_close(fit6, fit6b),
            info = "Ordinal growth model: parameters should be close to lavaan")

expect_true(
  all(diag(vcov(fit6)) > 0),
  info = "Ordinal growth model: vcov diagonal should be positive"
)
})


## =============================================================================
## 9. One-factor + exogenous covariates, missing data, fixed.x=TRUE
## =============================================================================

set.seed(1234)
pop.model_xo <- ' f =~ 0.7*y1 + 0.7*y2 + 0.7*y3 + 0.7*y4 + 0.7*y5
                  f ~ (-2.3)*x1 + 0.8*x2
                  y1 ~ 0.2*x2
                  y3 ~ 0.7*x1 '
Data_xo <- simulateData(pop.model_xo, sample.nobs = 100)
Data_xo <- makeord(Data_xo, ncat = 5)

## introduce missingness
obs    <- rbinom(prod(dim(Data_xo)), 1, .9)
Data_xo <- Data_xo * obs
Data_xo[Data_xo == 0] <- NA

model_xo <- ' f =~ y1 + y2 + y3 + y4 + y5
              f ~ x1 + x2
              y1 ~ x2
              y3 ~ x1 '

try({
fit_xo <- bsem(model_xo, data = Data_xo, fixed.x = TRUE,
               burnin = 100, sample = 200, ordered = TRUE, save.lvs = TRUE)

fit_xob <- sem(model_xo, data = Data_xo, fixed.x = TRUE,
               meanstructure = TRUE, ordered = TRUE, estimator = "pml",
               parameterization = "theta", missing = "pairwise")

expect_inherits(fit_xo, "blavaan",
  info = "One-factor + exo + missing, fixed.x=TRUE: bsem() should return a blavaan object")

expect_true(coef_close(fit_xo, fit_xob),
            info = "One-factor + exo + missing, fixed.x=TRUE: bsem() should return a blavaan object")

expect_true(
  all(diag(vcov(fit_xo)) > 0),
  info = "One-factor + exo + missing: vcov diagonal should be positive"
)

smry_xo <- summary(fit_xo)
expect_true(!is.null(smry_xo),
  info = "One-factor + exo + missing: summary() should return non-null")

fm_xo <- fitMeasures(fit_xo)
expect_true(fm_xo["ppp"] >= 0 & fm_xo["ppp"] <= 1,
  info = "One-factor + exo + missing: PPP should be in [0, 1]")
})


## =============================================================================
## 10. One-factor + exogenous covariates, fixed.x=FALSE
## =============================================================================

try({
fit_xo_fx <- bsem(model_xo, data = Data_xo, fixed.x = FALSE,
                  sample = 100, burnin = 100, target = mytarg, ordered = TRUE)

fit_xob_fx <- sem(model_xo, data = Data_xo, fixed.x = FALSE,
                  conditional.x = FALSE, meanstructure = TRUE, estimator = "pml",
                  ordered = TRUE, parameterization = "theta")

expect_inherits(fit_xo_fx, "blavaan",
  info = "One-factor + exo, fixed.x=FALSE: bsem() should return a blavaan object")

expect_true(coef_close(fit_xo_fx, fit_xob_fx),
            info = "One-factor + exo + missing, fixed.x=FALSE: bsem() should return a blavaan object")

expect_true(
  all(diag(vcov(fit_xo_fx)) > 0),
  info = "One-factor + exo, fixed.x=FALSE: vcov diagonal should be positive"
)

expect_true(fm_xo_fx["ppp"] >= 0 & fm_xo_fx["ppp"] <= 1,
  info = "One-factor + exo + missing: PPP should be in [0, 1]")
})


## =============================================================================
## 11. One-variable ordinal model (intercept + variance only)
## =============================================================================

set.seed(1234)
Data_1v <- makeord(data.frame(y1 = rnorm(100)), ncat = 4)

model_1v <- 'y1 ~ 1; y1 ~~ y1'

try({
fit_1v <- blavaan(model_1v, data = Data_1v, burnin = 100, sample = 100,
                  ordered = TRUE)

fit_1vb <- lavaan(model_1v, data = Data_1v, ordered = TRUE,
                  parameterization = "theta")

expect_inherits(fit_1v, "blavaan",
  info = "One-variable ordinal: blavaan() should return a blavaan object")

expect_equal(
  sort(names(coef(fit_1v))),
  sort(names(coef(fit_1vb))),
  info = "One-variable ordinal: parameter names should match lavaan"
)

expect_true(
  all(diag(vcov(fit_1v)) > 0),
  info = "One-variable ordinal: vcov diagonal should be positive"
)

smry_1v <- summary(fit_1v)
expect_true(!is.null(smry_1v),
  info = "One-variable ordinal: summary() should return non-null")

fm_1v <- fitMeasures(fit_1v)
expect_true("ppp" %in% names(fm_1v),
  info = "One-variable ordinal: fitMeasures() should contain ppp")
})


## =============================================================================
## 12. Simple regression with ordinal outcome, fixed.x=TRUE
## =============================================================================

set.seed(1234)
x1 <- rnorm(100)
y1 <- 0.5 + 2*x1 + rnorm(100)
Data_reg <- makeord(data.frame(y1 = y1, x1 = x1), ncat = 4)

model_reg <- ' y1 ~ x1 '

try({
fit_reg <- bsem(model_reg, data = Data_reg, fixed.x = TRUE,
                sample = 100, burnin = 100, ordered = TRUE)

fit_regb <- sem(model_reg, data = Data_reg, fixed.x = TRUE,
                ordered = TRUE, parameterization = "theta",
                conditional.x = FALSE, estimator = "PML")

expect_inherits(fit_reg, "blavaan",
  info = "Simple ordinal regression fixed.x=TRUE: bsem() should return a blavaan object")

expect_true(coef_close(fit_reg, fit_regb),
            info = "Simple ordinal regression fixed.x=TRUE: bsem() should return a blavaan object")

expect_true(
  all(diag(vcov(fit_reg)) > 0),
  info = "Simple ordinal regression: vcov diagonal should be positive"
)

smry_reg <- summary(fit_reg)
expect_true(!is.null(smry_reg),
  info = "Simple ordinal regression: summary() should return non-null")

fm_reg <- fitMeasures(fit_reg)
expect_true(fm_reg["ppp"] >= 0 & fm_reg["ppp"] <= 1,
  info = "Simple ordinal regression: PPP should be in [0, 1]")
})


## =============================================================================
## 13. Three-variable path analysis with ordinal variables, fixed.x=TRUE
## =============================================================================

set.seed(1234)
x1 <- rnorm(100)
y1 <- 0.5 + 2*x1 + rnorm(100)
y2 <- 0.8 + 0.4*y1 + rnorm(100)
Data_pa <- makeord(data.frame(y1 = y1, y2 = y2, x1 = x1), ncat = c(2, 4, 3))

model_pa <- ' y2 ~ y1; y1 ~ x1 '

try({
fit_pa <- bsem(model_pa, data = Data_pa, fixed.x = TRUE,
               sample = 100, burnin = 100, ordered = TRUE)

fit_pab <- sem(model_pa, data = Data_pa, fixed.x = TRUE,
               ordered = TRUE, parameterization = "theta",
               conditional.x = FALSE, estimator = "PML")

expect_inherits(fit_pa, "blavaan",
  info = "Ordinal path analysis fixed.x=TRUE: bsem() should return a blavaan object")

expect_equal(
  sort(names(coef(fit_pa))),
  sort(names(coef(fit_pab))),
  info = "Ordinal path analysis fixed.x=TRUE: parameter names should match lavaan"
)

expect_true(
  all(diag(vcov(fit_pa)) > 0),
  info = "Ordinal path analysis fixed.x=TRUE: vcov diagonal should be positive"
)

smry_pa <- summary(fit_pa)
expect_true(!is.null(smry_pa),
  info = "Ordinal path analysis fixed.x=TRUE: summary() should return non-null")

fm_pa <- fitMeasures(fit_pa)
expect_true("ppp" %in% names(fm_pa),
  info = "Ordinal path analysis fixed.x=TRUE: fitMeasures() should contain ppp")
})


## =============================================================================
## 14. Three-variable path analysis, fixed.x=FALSE
## =============================================================================

try({
fit_pa_fx <- bsem(model_pa, data = Data_pa, fixed.x = FALSE,
                  sample = 100, burnin = 100, ordered = TRUE)

fit_pab_fx <- sem(model_pa, data = Data_pa, fixed.x = FALSE,
                  ordered = TRUE, parameterization = "theta",
                  conditional.x = FALSE)

expect_inherits(fit_pa_fx, "blavaan",
  info = "Ordinal path analysis fixed.x=FALSE: bsem() should return a blavaan object")

expect_equal(
  sort(names(coef(fit_pa_fx))),
  sort(names(coef(fit_pab_fx))),
  info = "Ordinal path analysis fixed.x=FALSE: parameter names should match lavaan"
)

expect_true(
  all(diag(vcov(fit_pa_fx)) > 0),
  info = "Ordinal path analysis fixed.x=FALSE: vcov diagonal should be positive"
)

smry_pa_fx <- summary(fit_pa_fx)
expect_true(!is.null(smry_pa_fx),
  info = "Ordinal path analysis fixed.x=FALSE: summary() should return non-null")

fm_pa_fx <- fitMeasures(fit_pa_fx)
expect_true("ppp" %in% names(fm_pa_fx),
  info = "Ordinal path analysis fixed.x=FALSE: fitMeasures() should contain ppp")
})

}
