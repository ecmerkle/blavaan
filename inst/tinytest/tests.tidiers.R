# Tests for tidy and glance methods

if (requireNamespace("rstan", quietly = TRUE)) {
  load(system.file("testdata", "sysdata.rda", package = "blavaan"))
  library("lavaan", quietly = TRUE)

  # Test tidy.blavaan with Stan model
  tidy_stan <- tidy.blavaan(fitstan)
  expect_true(inherits(tidy_stan, "data.frame"))
  expected_tidy_cols <- c("term", "estimate", "std.error", "conf.low",
                          "conf.high", "rhat", "ess", "prior")
  expect_true(all(expected_tidy_cols %in% names(tidy_stan)))
  expect_true(nrow(tidy_stan) > 0)

  # Test tidy without confidence intervals
  tidy_noci <- tidy.blavaan(fitstan, conf.int = FALSE)
  expect_false(any(c("conf.low", "conf.high") %in% names(tidy_noci)))

  # Test tidy without rhat and ess
  tidy_nodiag <- tidy.blavaan(fitstan, rhat = FALSE, ess = FALSE)
  expect_false(any(c("rhat", "ess") %in% names(tidy_nodiag)))

  # Test tidy with standardized estimates
  tidy_std <- tidy.blavaan(fitstan, standardized = TRUE)
  expect_true("std.all" %in% names(tidy_std))

  # Test glance.blavaan with Stan model
  glance_stan <- glance.blavaan(fitstan)
  expect_true(inherits(glance_stan, "data.frame"))
  expect_equal(nrow(glance_stan), 1)
  expected_glance_cols <- c("npar", "nobs", "ngroups", "estimator", "nchains")
  expect_true(all(expected_glance_cols %in% names(glance_stan)))
  expect_true(any(c("ppp", "looic") %in% names(glance_stan)))

  # Test glance with fit indices (this is computationally expensive)
  # Only test if blavFitIndices works
  bfi_works <- tryCatch({
    blavFitIndices(fitstan, fit.measures = "BRMSEA")
    TRUE
  }, error = function(e) FALSE)

  if (bfi_works) {
    glance_bfi <- glance.blavaan(fitstan, fit.indices = "BRMSEA")
    expect_true("BRMSEA" %in% names(glance_bfi))
  }

  # Test with classic Stan model
  tidy_stanc <- tidy.blavaan(fitstanc)
  expect_true(inherits(tidy_stanc, "data.frame"))
  expect_true(nrow(tidy_stanc) > 0)

  glance_stanc <- glance.blavaan(fitstanc)
  expect_true(inherits(glance_stanc, "data.frame"))
  expect_equal(nrow(glance_stanc), 1)

  # Test with JAGS model if available
  if (requireNamespace("runjags", quietly = TRUE)) {
    tidy_jags <- tidy.blavaan(fitjags)
    expect_true(inherits(tidy_jags, "data.frame"))
    expect_true(nrow(tidy_jags) > 0)

    glance_jags <- glance.blavaan(fitjags)
    expect_true(inherits(glance_jags, "data.frame"))
    expect_equal(nrow(glance_jags), 1)
  }
}
