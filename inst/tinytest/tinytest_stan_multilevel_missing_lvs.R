## tinytest file: latent variable (factor score) sampling for two-level
## (multilevel) models with FIML/MAR missing data, Stan backend
##
## Run with: tinytest::run_test_file("tinytest_stan_multilevel_missing_lvs.R")
## or as part of a package: tinytest::test_package("blavaan")
##
## Companion to tinytest_stan_multilevel_missing.R (which covers the FIML/MAR
## *likelihood* only -- parameter estimation, waic/looic, ppp) and
## tinytest_stan_multilevel2.R (which covers complete-data two-level LV
## sampling). This file exercises save.lvs = TRUE /
## blavPredict(type = "lv") / blavInspect(., "lvs") for two-level models
## with real missingness, i.e. samp_lvs_2lev()'s use of lavaan's
## lav_mvn_cl_mi_estep_ranef() (level 2) and genuine within-cluster FIML
## patterns built from lavdata@Mp (level 1). See R/lvgqs.R.

library("tinytest")
library("lavaan")
library("blavaan")

dp_stable <- dpriors(lambda = "normal(1,.4)")

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

data(Demo.twolevel, package = "lavaan")

## y4 is within-only, y1-y3 are "both", w1/w2 are between-only -- same
## variable-class mix as tinytest_stan_multilevel_missing.R's model_wb
model_wb <- '
    level: 1
        fw =~ y1 + y2 + y3 + y4
    level: 2
        fb =~ y1 + y2 + y3
        fb ~~ w1
        fb ~~ w2
'

## =============================================================================
## 1. Column-mapping sanity check (no missingness): lavaan's missing-data
## cluster-posterior function (lav_mvn_cl_mi_estep_ranef(), used when
## real missingness is present) should exactly reproduce the existing
## complete-data cluster-mean function (lav_mvn_cl_em_estep_ranef()) when
## fit with missing = "ml" but no actual NAs -- this isolates the column/
## index mapping (both_idx, y1w space, etc.) from any real missing-data
## statistical questions.
## =============================================================================
try({

fit0 <- lavaan::lavaan(model_wb, data = Demo.twolevel, cluster = "cluster",
                        meanstructure = TRUE, missing = "ml", fixed.x = FALSE,
                        do.fit = FALSE)

lav_implied22l <- getFromNamespace("lav_mvnorm_cluster_implied22l", "lavaan")
lav_estep_old  <- getFromNamespace("lav_mvn_cl_em_estep_ranef", "lavaan")
lav_mi_estep   <- getFromNamespace("lav_mvn_cl_mi_estep_ranef", "lavaan")

modimp <- lav_model_implied(fit0@Model)
implied_g <- lapply(modimp, function(x) x[1:2])

out <- lav_implied22l(fit0@Data@Lp[[1]], implied_g)
clusmns_old <- lav_estep_old(ylp = fit0@SampleStats@YLp[[1]], lp = fit0@Data@Lp[[1]],
                              sigma_w = out$sigma.w, sigma_b = out$sigma.b,
                              sigma_zz = out$sigma.zz, sigma_yz = out$sigma.yz,
                              mu_z = out$mu.z, mu_w = out$mu.w, mu_b = out$mu.b, se = FALSE)

y2g <- rowsum.default(fit0@Data@X[[1]], group = fit0@Data@Lp[[1]]$cluster.idx[[2]],
                       reorder = FALSE, na.rm = FALSE) / fit0@Data@Lp[[1]]$cluster.size[[2]]
mb_j <- lav_mi_estep(y1 = fit0@Data@X[[1]], y2 = y2g,
                      lp = fit0@Data@Lp[[1]], mp = fit0@Data@Mp[[1]],
                      mu_w = implied_g$mean[[1]], sigma_w = implied_g$cov[[1]],
                      mu_b = implied_g$mean[[2]], sigma_b = implied_g$cov[[2]],
                      se = FALSE, impute = TRUE)

expect_equal(dim(clusmns_old), dim(mb_j),
  info = "old (complete-data) and new (missing-data) cluster-mean outputs should have the same shape")
expect_true(max(abs(clusmns_old - mb_j)) < 1e-8,
  info = "lav_mvn_cl_mi_estep_ranef() should exactly reproduce lav_mvn_cl_em_estep_ranef() when nothing is actually missing")

between.idx <- fit0@Data@Lp[[1]]$between.idx[[2]]
raw_between <- fit0@Data@X[[1]][!duplicated(fit0@Data@Lp[[1]]$cluster.idx[[2]]), between.idx]
expect_true(max(abs(attr(mb_j, "z.imputed") - raw_between)) < 1e-8,
  info = "z.imputed should exactly reproduce the raw between-only data when nothing is actually missing")
})

## =============================================================================
## 2. save.lvs = TRUE with real missingness at all three variable classes
## (within-only, both, between-only): blavPredict()/blavInspect() should
## return finite, complete (no NA/NaN) latent variable draws at both levels
## =============================================================================
try({

set.seed(101)
d1 <- Demo.twolevel
d1 <- inject_mcar(d1, c("y1", "y2", "y3", "y4"), rate = 0.15, seed = 101)
d1 <- inject_mcar_between(d1, "w1", rate = 0.15, seed = 102)

bfit1 <- bsem(
  model   = model_wb,
  data    = d1,
  cluster = "cluster",
  fixed.x = FALSE,
  burnin  = 100,
  sample  = 100,
  dp      = dp_stable,
  save.lvs = TRUE
)

expect_inherits(bfit1, "blavaan",
  info = "bsem() with save.lvs=TRUE under two-level FIML/MAR should return a blavaan object")

lv1 <- blavInspect(bfit1, "lvs", level = 1)
lv2 <- blavInspect(bfit1, "lvs", level = 2)

m1 <- as.matrix(do.call(rbind, lv1))
m2 <- as.matrix(do.call(rbind, lv2))

expect_equal(ncol(m1), nrow(d1),
  info = "level-1 lvs should have one column per row of the original data, including any 'empty' (all-level-1-vars-missing) rows")
expect_equal(ncol(m2), length(unique(d1$cluster)),
  info = "level-2 lvs should have one column per cluster")

expect_false(any(is.na(m1)), info = "level-1 lvs should contain no NA")
expect_false(any(is.nan(m1)), info = "level-1 lvs should contain no NaN")
expect_false(any(is.na(m2)), info = "level-2 lvs should contain no NA")
expect_false(any(is.nan(m2)), info = "level-2 lvs should contain no NaN")

pred1 <- blavPredict(bfit1, type = "lv", level = 1)
pred2 <- blavPredict(bfit1, type = "lv", level = 2)
expect_false(any(is.na(pred1)), info = "blavPredict level-1 lvs should contain no NA")
expect_false(any(is.na(pred2)), info = "blavPredict level-2 lvs should contain no NA")

## cross-check against lavaan's own two-level FIML factor scores (EAP vs
## posterior mean, so a loose comparison; both estimate the same quantity)
fscores_lav <- lavPredict(sem(model_wb, data = d1, cluster = "cluster",
                               missing = "ml", fixed.x = FALSE),
                           level = 1)
cor1 <- suppressWarnings(cor(colMeans(m1), fscores_lav[, "fw"], use = "complete.obs"))
expect_true(is.finite(cor1) && abs(cor1) > 0.5,
  info = "posterior mean level-1 factor scores should correlate reasonably with lavaan's FIML factor scores")
})

## =============================================================================
## 3. Complete-data equivalence: save.lvs = TRUE two-level LV sampling with
## NO missing data should be unaffected by the missing-data code path
## (regression check for the "not missing" branch in samp_lvs_2lev())
## =============================================================================
try({

bfit_complete <- bsem(
  model   = model_wb,
  data    = Demo.twolevel,
  cluster = "cluster",
  fixed.x = FALSE,
  burnin  = 100,
  sample  = 100,
  dp      = dp_stable,
  save.lvs = TRUE
)

lv1c <- blavInspect(bfit_complete, "lvs", level = 1)
lv2c <- blavInspect(bfit_complete, "lvs", level = 2)
m1c <- as.matrix(do.call(rbind, lv1c))
m2c <- as.matrix(do.call(rbind, lv2c))

expect_false(any(is.na(m1c)), info = "complete-data level-1 lvs should contain no NA")
expect_false(any(is.na(m2c)), info = "complete-data level-2 lvs should contain no NA")
expect_equal(ncol(m1c), nrow(Demo.twolevel),
  info = "complete-data level-1 lvs should have one column per row")
expect_equal(ncol(m2c), length(unique(Demo.twolevel$cluster)),
  info = "complete-data level-2 lvs should have one column per cluster")
})
