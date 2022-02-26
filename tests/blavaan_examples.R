## short examples to test functionality
set.seed(341)

x1 <- rnorm(100)
y1 <- 0.5 + 2*x1 + rnorm(100)
g <- rep(1:2, each=50)
Data <- data.frame(y1 = y1, x1 = x1, g = g)

## seemed to help if running this via R CMD check:
## Sys.unsetenv('R_TESTS')
library("blavaan")
## NB! unload DBI package or problems with CRAN!
unloadNamespace("DBI")

## don't care that models are not converged, keeping file size small
model <- ' y1 ~ prior("dnorm(0,1)")*x1 '
fitjags <- bsem(model, data=Data, fixed.x=TRUE, burnin=20,
                sample=20, target="jags", group="g", seed=1:3)

model <- ' y1 ~ prior("normal(0,1)")*x1 '
fitstan <- bsem(model, data=Data, fixed.x=TRUE, burnin=20,
                sample=20, target="stan", group="g", seed=1)

fitstanc <- bsem(model, data=Data, fixed.x=TRUE, burnin=20,
                 sample=20, target="stanclassic", group="g", seed=1)

fitstanmomentmatch <- bsem(
  model, 
  data=Data, 
  fixed.x=TRUE, 
  burnin=20,
  sample=20,
  mcmcextra=list(data=list(moment_match_k_threshold=0.5)),
  target="stan", 
  group="g", 
  seed=1
)

## for checking factor score functionality
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6 '

fitstanfs <- bcfa(HS.model, data=HolzingerSwineford1939,
                  burnin=30, sample=10, target="stan",
                  save.lvs=TRUE, n.chains=2, seed=1)

## this really blows up file size if kept:
attr(fitstan@external$mcmcout, 'stanmodel') <- NULL
attr(fitstanc@external$mcmcout, 'stanmodel') <- NULL
attr(fitstanfs@external$mcmcout, 'stanmodel') <- NULL

save(list=c("fitjags", "fitstan", "fitstanc", "fitstanfs"), file="../inst/testdata/sysdata.rda")

