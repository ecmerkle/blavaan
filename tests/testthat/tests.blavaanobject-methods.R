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

    HS.model <- ' visual  =~ x1 + x2 + x3
                  textual =~ x4 + x5 + x6 '
    fitlav <- cfa(HS.model, data=HolzingerSwineford1939,
                  meanstructure=TRUE)
    expect_true(cor(blavInspect(fitstanfs, 'lvmeans')[,1],
                    lavPredict(fitlav, type='lv')[,1]) > .95)

    ## plots
    expect_silent(p <- plot(fitstan, showplot = FALSE))
    expect_false(any(p$data$value[grep('~~', p$data$parameter)] < 0)) # check of parameter labels
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

  ## blavFit + ppmc
  discFUN <- list(global = function(fit) {
    fitMeasures(fit, fit.measures = c("cfi","rmsea","srmr","chisq"))
  },
  std.cov.resid = function(fit) lavResiduals(fit, zstat = FALSE,
                                             summary = FALSE)$`1`$cov)
  ppmc_res <- ppmc(fitstan, discFUN = discFUN)
  expect_s4_class(ppmc_res, "blavPPMC")
  ppmc_summ <- summary(ppmc_res, "global", cent = "EAP")
  expect_s3_class(ppmc_summ, "lavaan.data.frame")
  ppmc_summ <- summary(ppmc_res, "std.cov.resid", cent = "MAP",
                       to.data.frame = TRUE, sort.by = "MAP",
                       decreasing = TRUE)
  bf_res <- blavFitIndices(fitstan)
  expect_s4_class(bf_res, "blavFitIndices")
  expect_s3_class(summary(bf_res), "lavaan.data.frame")
})
