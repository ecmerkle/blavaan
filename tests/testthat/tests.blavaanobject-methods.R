test_that("blavaan object methods work", {

  if(requireNamespace("rstan", quietly = TRUE) &
     requireNamespace("runjags", quietly = TRUE)){
    load(system.file("testdata", "sysdata.rda", package="blavaan"))

    # classes
    expect_equal(class(fitjags@external), "list")
    expect_equal(class(fitstan@external), "list")
    expect_equal(class(fitstanc@external), "list")
    
    ## parameter summaries
    expect_equal(dim(parTable(fitjags)), c(10,20))
    expect_equal(dim(parTable(fitstan)), c(10,24))
    expect_equal(dim(parTable(fitstanc)), c(10,21))
    
    expect_equal(sum(fitjags@ParTable$free > 0, na.rm = TRUE),
                 length(blavInspect(fitjags, 'psrf')))
    expect_equal(sum(fitstan@ParTable$free > 0, na.rm = TRUE),
                 length(blavInspect(fitstan, 'psrf')))
    expect_equal(sum(fitstanc@ParTable$free > 0, na.rm = TRUE),
                 length(blavInspect(fitstanc, 'psrf')))
    expect_equal(fitjags@ParTable$free, fitstan@ParTable$free)
    expect_equal(fitjags@ParTable$free, fitstanc@ParTable$free)
    expect_equal(nrow(parTable(fitjags)), nrow(parTable(fitstan)))
    expect_equal(nrow(parTable(fitjags)), nrow(parTable(fitstanc)))
    
    expect_error(blavInspect(fitjags, 'blah'))

    ## fitMeasures
    expect_equal(length(fitMeasures(fitjags)),
                 length(fitMeasures(fitstan)))
    expect_equal(length(fitMeasures(fitjags)),
                 length(fitMeasures(fitstanc)))

    ## this is how summary() obtains its results, but have not figured out
    ## how to get S4 methods to directly work in testthat
    expect_equal(dim(parameterEstimates(fitjags)), c(10, 6))
    expect_equal(dim(parameterEstimates(fitstan)), c(10, 6))
    expect_equal(dim(parameterEstimates(fitstanc)), c(10, 6))

    ## various blavInspect args
    expect_equal(length(blavInspect(fitjags, 'psrf')),
                 length(blavInspect(fitstan, 'psrf')))

    expect_equal(length(blavInspect(fitjags, 'psrf')),
                 length(blavInspect(fitstanc, 'psrf')))
    
    expect_equal(length(blavInspect(fitjags, 'neff')),
                 length(blavInspect(fitstan, 'neff')))

    expect_equal(length(blavInspect(fitjags, 'neff')),
                 length(blavInspect(fitstanc, 'neff')))

    expect_equal(length(blavInspect(fitjags, 'mcmc')),
                 length(blavInspect(fitstan, 'mcmc')))

    expect_equal(length(blavInspect(fitjags, 'mcmc')),
                 length(blavInspect(fitstanc, 'mcmc')))
    
    expect_equal(length(blavInspect(fitjags, 'start')),
                 length(blavInspect(fitstan, 'start')))

    expect_equal(length(blavInspect(fitjags, 'start')),
                 length(blavInspect(fitstanc, 'start')))
    
    expect_equal(dim(blavInspect(fitjags, 'hpd')),
                 dim(blavInspect(fitstan, 'hpd')))

    expect_equal(dim(blavInspect(fitjags, 'hpd')),
                 dim(blavInspect(fitstanc, 'hpd')))

    expect_equal(dim(standardizedposterior(fitjags)),
                 dim(standardizedposterior(fitstan)))

    expect_equal(dim(standardizedposterior(fitjags)),
                 dim(standardizedposterior(fitstanc)))

    expect_equal(dim(blavInspect(fitstanfs, 'lvmeans')),
                 c(301, 2))

    expect_equal(dim(blavInspect(fitstanfs, 'lvs')[[2]]),
                 c(10, 602))

    ## plots
    expect_silent(p <- plot(fitstan, showplot = FALSE))
    expect_silent(p <- plot(fitstan, 1:4, showplot = FALSE))
    expect_silent(p <- plot(fitstan, plot.type = "hist", showplot = FALSE))
    expect_silent(p <- plot(fitstan, 1:4, plot.type = "dens", showplot = FALSE))
    expect_silent(p <- plot(fitstan, c(2,4), plot.type = "scatter", showplot = FALSE))

    expect_silent(p <- plot(fitstanc, showplot = FALSE))
    expect_silent(p <- plot(fitstanc, 1:4, showplot = FALSE))
    expect_silent(p <- plot(fitstanc, plot.type = "hist", showplot = FALSE))
    expect_silent(p <- plot(fitstanc, 1:4, plot.type = "dens", showplot = FALSE))
    expect_silent(p <- plot(fitstanc, c(2,4), plot.type = "scatter", showplot = FALSE))

    expect_silent(p <- plot(fitjags, showplot = FALSE))
    expect_silent(p <- plot(fitjags, 1:4, showplot = FALSE))
    expect_silent(p <- plot(fitjags, plot.type = "hist", showplot = FALSE))
    expect_silent(p <- plot(fitjags, 1:4, plot.type = "dens", showplot = FALSE))
    expect_silent(p <- plot(fitjags, c(2,4), plot.type = "scatter", showplot = FALSE))
  }
    
})
