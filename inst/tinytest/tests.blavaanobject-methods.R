library("lavaan", quietly = TRUE)
library("blavaan", quietly = TRUE)

## check that jags is installed
hasjags <- !all(Sys.which(c("jags", "jags-terminal")) == "")
  
## short examples to test functionality
set.seed(341)

## seems to be needed for stanclassic:
oopts <- options(future.globals.maxSize = 1.0 * 1e9)
on.exit(options(oopts))

x1 <- rnorm(100)
y1 <- 0.5 + 2*x1 + rnorm(100)
g <- rep(1:2, each=50)
Data <- data.frame(y1 = y1, x1 = x1, g = g)

if (hasjags) {
  model <- ' y1 ~ prior("dnorm(0,1)")*x1 '
  fitjags <- bsem(model, data=Data, fixed.x=TRUE, burnin=20,
                  sample=20, target="jags", group="g", seed=1:3)
}
 
model <- ' y1 ~ prior("normal(0,1)")*x1 '
fitstan <- bsem(model, data=Data, fixed.x=TRUE, burnin=20,
                sample=20, target="stan", group="g", seed=1, meanstructure=TRUE)

fitstanc <- bsem(model, data=Data, fixed.x=TRUE, burnin=20,
                 sample=20, target="stanclassic", group="g", seed=1)

## for checking factor score functionality
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6 '

fitstanfs <- bcfa(HS.model, data=lavaan::HolzingerSwineford1939,
                  burnin=30, sample=10, target="stan",
                  save.lvs=TRUE, n.chains=2, seed=1, meanstructure=TRUE)  
  
options(mc.cores = 1L)

## classes
expect_equal(class(fitstan@external), "list")
expect_equal(class(fitstanc@external), "list")

## parameter summaries
expect_equal(dim(parTable(fitstan)), c(10,24))
expect_equal(dim(parTable(fitstanc)), c(10,21))    
expect_equal(sum(fitstan@ParTable$free > 0, na.rm = TRUE),
             length(blavInspect(fitstan, 'psrf')))
expect_equal(sum(fitstanc@ParTable$free > 0, na.rm = TRUE),
             length(blavInspect(fitstanc, 'psrf')))
    
## fitMeasures
expect_equal(length(fitMeasures(fitstan)),
             length(fitMeasures(fitstanc)))

## this is how summary() obtains its results, but have not figured out
## how to get S4 methods to directly work in testthat
expect_equal(dim(parameterEstimates(fitstan)), c(10, 6))
expect_equal(dim(parameterEstimates(fitstanc)), c(10, 6))

expect_equal(dim(standardizedposterior(fitjags)),
             dim(standardizedposterior(fitstan)))

## various blavInspect args
expect_equal(length(blavInspect(fitstan, 'psrf')),
             length(blavInspect(fitstanc, 'psrf')))

expect_equal(length(blavInspect(fitstan, 'neff')),
             length(blavInspect(fitstanc, 'neff')))

expect_equal(length(blavInspect(fitstan, 'mcmc')),
             length(blavInspect(fitstanc, 'mcmc')))

expect_equal(length(blavInspect(fitstan, 'start')),
             length(blavInspect(fitstanc, 'start')))

expect_equal(dim(blavInspect(fitstan, 'hpd')),
             dim(blavInspect(fitstanc, 'hpd')))

expect_equal(dim(standardizedposterior(fitstan)),
             dim(standardizedposterior(fitstanc)))


if (hasjags) {
  expect_equal(class(fitjags@external), "list")
  expect_equal(dim(parTable(fitjags)), c(10,20))
  expect_equal(sum(fitjags@ParTable$free > 0, na.rm = TRUE),
               length(blavInspect(fitjags, 'psrf')))
  expect_equal(fitjags@ParTable$free, fitstan@ParTable$free)
  expect_equal(fitjags@ParTable$free, fitstanc@ParTable$free)
  expect_equal(nrow(parTable(fitjags)), nrow(parTable(fitstan)))
  expect_equal(nrow(parTable(fitjags)), nrow(parTable(fitstanc)))
  expect_error(blavInspect(fitjags, 'blah'))
  expect_equal(length(fitMeasures(fitjags)),
               length(fitMeasures(fitstan)))
  expect_equal(length(fitMeasures(fitjags)),
               length(fitMeasures(fitstanc)))
  expect_equal(dim(parameterEstimates(fitjags)), c(10, 6))

  expect_equal(length(blavInspect(fitjags, 'psrf')),
               length(blavInspect(fitstan, 'psrf')))
    
  expect_equal(length(blavInspect(fitjags, 'neff')),
               length(blavInspect(fitstan, 'neff')))


  expect_equal(length(blavInspect(fitjags, 'mcmc')),
               length(blavInspect(fitstan, 'mcmc')))

    
  expect_equal(length(blavInspect(fitjags, 'start')),
               length(blavInspect(fitstan, 'start')))

    
  expect_equal(dim(blavInspect(fitjags, 'hpd')),
               dim(blavInspect(fitstan, 'hpd')))
}

  
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

if (hasjags) {
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
expect_equal(class(ppmc_res)[1], "blavPPMC")
ppmc_summ <- summary(ppmc_res, "global", cent = "EAP")
expect_equal(class(ppmc_summ)[1], "lavaan.data.frame")
ppmc_summ <- summary(ppmc_res, "std.cov.resid", cent = "MAP",
                     to.data.frame = TRUE, sort.by = "MAP",
                     decreasing = TRUE)
bf_res <- blavFitIndices(fitstan)
expect_equal(class(bf_res)[1], "blavFitIndices")
expect_equal(class(summary(bf_res))[1], "lavaan.data.frame")

bf_res <- blavFitIndices(fitstan, rescale = "ppmc")
expect_equal(class(bf_res)[1], "blavFitIndices")
expect_equal(class(summary(bf_res))[1], "lavaan.data.frame")  

bf_res <- blavFitIndices(fitstan, rescale = "mcmc")
expect_equal(class(bf_res)[1], "blavFitIndices")
expect_equal(class(summary(bf_res))[1], "lavaan.data.frame")
  
set.seed(341)
  
x1 <- rnorm(100)
y1 <- 0.5 + 2*x1 + rnorm(100)
g <- rep(1:2, each=50)
Data <- data.frame(y1 = y1, x1 = x1, g = g)
  
model <- ' y1 ~ prior("normal(0,1)")*x1 '
fitstanmomentmatch <- bsem(model, data=Data, fixed.x=TRUE, burnin=20, sample=20,
                           mcmcextra=list(data=list(moment_match_k_threshold=0.5)),
                           target="stan", seed=1)
bf_mm_res <- blavFitIndices(fitstanmomentmatch, fit.measures = c("looic"))
expect_equal(class(bf_mm_res)[1], "blavFitIndices")
expect_equal(class(summary(bf_mm_res))[1], "lavaan.data.frame")
expect_true("p_loo" %in% names(bf_mm_res@details$pD))

## blavPredict
expect_error(blavPredict(fitstanc))
if (hasjags) expect_error(blavPredict(fitjags))
  
expect_equal(dim(blavPredict(fitstanfs)[[1]]), c(301,2))
expect_equal(length(blavPredict(fitstanfs)), 20)
expect_equal(dim(blavPredict(fitstanfs, type="lvmeans")), c(301,2))
expect_equal(dim(blavPredict(fitstanfs, type="ov")[[1]]), c(301,6))
expect_equal(dim(blavPredict(fitstanfs, type="ypred")[[1]]), c(301,6))
expect_error(blavPredict(fitstanfs, type="ymis"))
