test_that("blavaan object methods work", {
  fitj <- blavaan:::fitjags
  fits <- blavaan:::fitstan

  library("blavaan")
  
  expect_equal(fitjags@ParTable$free, fitstan@ParTable$free)

})
