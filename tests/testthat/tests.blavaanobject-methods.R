test_that("blavaan object methods work", {
  library(runjags)
  library(rstan)
  library(blavaan)
  fitj <- blavaan:::fitjags
  fits <- blavaan:::fitstan

  expect_equal(fitjags@ParTable$free, fitstan@ParTable$free)

})
