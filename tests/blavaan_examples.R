## short examples to test functionality
set.seed(341)

x1 <- rnorm(100)
y1 <- 0.5 + 2*x1 + rnorm(100)
g <- rep(1:2, each=50)
Data <- data.frame(y1 = y1, x1 = x1, g = g)

## seemed to help if running this via R CMD check:
## Sys.unsetenv('R_TESTS')
library("blavaan")
model <- ' y1 ~ prior("dnorm(0,1)")*x1 '
fitjags <- bsem(model, data=Data, fixed.x=TRUE, burnin=200,
                sample=200, target="jags", group="g")

model <- ' y1 ~ prior("normal(0,1)")*x1 '
fitstan <- bsem(model, data=Data, fixed.x=TRUE, burnin=200,
                sample=200, target="stan", group="g")

fitstanc <- bsem(model, data=Data, fixed.x=TRUE, burnin=200,
                 sample=200, target="stanclassic", group="g")

## this really blows up file size if kept:
attr(fitstan@external$mcmcout, 'stanmodel') <- NULL
attr(fitstanc@external$mcmcout, 'stanmodel') <- NULL

save(list=c("fitjags", "fitstan", "fitstanc"), file="../inst/testdata/sysdata.rda")

