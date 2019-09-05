# blavaan

blavaan is a free, open source R package for Bayesian latent variable analysis.  It relies on JAGS and Stan to estimate models via MCMC.

The blavaan functions and syntax are similar to lavaan. For example, consider the Political Democracy example from Bollen (1989):

```
library(blavaan)

model <- '
   # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + y2 + y3 + y4
     dem65 =~ y5 + y6 + y7 + y8
   # regressions
     dem60 ~ ind60
     dem65 ~ ind60 + dem60
   # residual covariances
     y1 ~~ y5
     y2 ~~ y4 + y6
     y3 ~~ y7
     y4 ~~ y8
     y6 ~~ y8
'
fit <- bsem(model, data = PoliticalDemocracy)
summary(fit)
```

The development version of blavaan (containing updates not yet on CRAN) can be installed via either of the following commands. Compilation is required in both cases; this may be a problem for users who currently rely on a binary version of blavaan from CRAN.

```
# from github:
remotes::install_github("ecmerkle/blavaan", INSTALL_opts = "--no-multiarch")

# from website:
install.packages("blavaan", repos = "http://faculty.missouri.edu/~merklee", type = "source")
```

For further information, see:

Merkle, E. C., & Rosseel, Y. (2018). [blavaan: Bayesian structural equation models via parameter expansion](https://doi.org/10.18637/jss.v085.i04). Journal of Statistical Software, 85(4), 1–30.
