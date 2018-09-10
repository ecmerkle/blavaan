test_that("blavaan object methods work", {

  if(requireNamespace("rstan", quietly = TRUE) &
     requireNamespace("runjags", quietly = TRUE)){
    load(system.file("testdata", "sysdata.rda", package="blavaan"))
    #fitj <- blavaan:::fitjags
    #fits <- blavaan:::fitstan

    expect_equal(fitjags@ParTable$free, fitstan@ParTable$free)
  }
    
})
