test_that("blavaan object methods work", {

  if(requireNamespace("rstan", quietly = TRUE) &
     requireNamespace("runjags", quietly = TRUE)){
    load(system.file("testdata", "sysdata.rda", package="blavaan"))
    #fitj <- blavaan:::fitjags
    #fits <- blavaan:::fitstan

    # classes
    expect_equal(class(fitjags@external), "list")
    expect_equal(class(fitstan@external), "list")

    ## parameter summaries
    expect_equal(dim(parTable(fitjags)), c(10,20))
    expect_equal(dim(parTable(fitstan)), c(10,25))

    expect_equal(sum(fitjags@ParTable$free > 0, na.rm = TRUE),
                 length(blavInspect(fitjags, 'psrf')))
    expect_equal(sum(fitstan@ParTable$free > 0, na.rm = TRUE),
                 length(blavInspect(fitstan, 'psrf')))
    expect_equal(fitjags@ParTable$free, fitstan@ParTable$free)
    expect_equal(nrow(parTable(fitjags)), nrow(parTable(fitstan)))
    
    expect_error(blavInspect(fitjags, 'blah'))

    ## fitMeasures
    expect_equal(length(fitMeasures(fitjags)),
                 length(fitMeasures(fitstan)))
  }
    
})
