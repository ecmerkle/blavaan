# Convergence and Efficiency Evaluation

### Introduction

When Bayesian models are estimated with a Markov-Chain Monte Carlo
(MCMC) sampler, the model estimation doesn’t stop when it has achieved
some convergence criteria. It will run as long as desired (determined by
the `burnin` and `sample` arguments), and then you need to evaluate the
convergence and efficiency of the estimated posterior distributions. You
should only analyze the results if convergence has been achieved, as
judged by the metrics described below.

For this example we will use the Industrialization and Political
Democracy example (Bollen 1989).

``` r
model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ a*y1 + b*y2 + c*y3 + d*y4
     dem65 =~ a*y5 + b*y6 + c*y7 + d*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

fit <- bsem(model, data=PoliticalDemocracy,
            std.lv=T, meanstructure=T, n.chains=3,
            burnin=500, sample=1000)
```

### Convergence

The primary convergence diagnostic is $\widehat{R}$, which compares the
between- and within-chain samples of model parameters and other
univariate quantities of interest (Vehtari et al. 2021). If chains have
not mixed well (ie, the between- and within-chain estimates don’t
agree), $\widehat{R}$ is larger than 1. We recommend running at least
three chains by default and only using the posterior samples if
$\widehat{R} < 1.05$ for all the parameters.

`blavaan` presents the $\widehat{R}$ reported by the underlying MCMC
program, either Stan or JAGS (Stan by default). We can obtain the
$\widehat{R}$ from the
[`summary()`](https://rdrr.io/r/base/summary.html) function, and we can
also extract it with the
[`blavInspect()`](http://ecmerkle.github.io/blavaan/reference/blavInspect.md)
function

``` r
blavInspect(fit, "rhat")
```

    ##   ind60=~x1   ind60=~x2   ind60=~x3           a           b           c 
    ##   1.0009401   1.0010621   1.0006748   1.0008258   1.0005812   1.0004656 
    ##           d           a           b           c           d dem60~ind60 
    ##   1.0001249   1.0008258   1.0005812   1.0004656   1.0001249   0.9996787 
    ## dem65~ind60 dem65~dem60      y1~~y5      y2~~y4      y2~~y6      y3~~y7 
    ##   1.0014211   1.0000992   1.0015324   0.9997618   0.9997201   1.0004598 
    ##      y4~~y8      y6~~y8      x1~~x1      x2~~x2      x3~~x3      y1~~y1 
    ##   1.0002699   1.0018175   1.0009606   1.0023951   0.9996752   1.0016669 
    ##      y2~~y2      y3~~y3      y4~~y4      y5~~y5      y6~~y6      y7~~y7 
    ##   0.9992964   0.9991829   0.9994871   1.0003425   0.9997676   0.9996398 
    ##      y8~~y8        x1~1        x2~1        x3~1        y1~1        y2~1 
    ##   1.0011487   1.0004498   1.0012514   1.0021230   1.0034488   1.0017461 
    ##        y3~1        y4~1        y5~1        y6~1        y7~1        y8~1 
    ##   1.0023612   1.0035633   1.0023965   1.0019605   1.0026057   1.0030892

With large models it can be cumbersome to look over all of these
entries. We can instead find the largest $\widehat{R}$ to see if they
are all less than $1.05$

``` r
max(blavInspect(fit, "psrf"))
```

    ## [1] 1.003563

If all $\widehat{R} < 1.05$ then we can establish that the MCMC chains
have converged to a stable solution. If the model has not converged, you
might increase the number of `burnin` iterations

``` r
fit <- bsem(model, data=PoliticalDemocracy,
            std.lv=T, meanstructure=T, n.chains=3,
            burnin=1000, sample=1000)
```

and/or change the model priors with the
[`dpriors()`](http://ecmerkle.github.io/blavaan/reference/dpriors.md)
function. These address issues where the model failed to converge due to
needing more iterations or due to a model misspecification (such as bad
priors). As a rule of thumb, we seldom see a model require more than
1,000 burnin samples in Stan. If your model is not converging after
1,000 burnin samples, it is likely that the default prior distributions
clash with your data. This can happen, e.g., if your variables contain
values in the 100s or 1000s.

### Efficiency

We should also evaluate the efficiency of the posterior samples.
Effective sample size (ESS) is a useful measure for sampling efficiency,
and is well defined even if the chains do not have finite mean or
variance (Vehtari et al. 2021).

In short, the posterior samples produced by MCMC are autocorrelated.
This means that, if you draw 500 posterior samples, you do not have 500
independent pieces of information about the posterior distribution,
because the samples are autocorlated. The ESS metric is [like a currency
conversion,](https://www.johndcook.com/blog/2017/06/27/effective-sample-size-for-mcmc/)
telling you how much your autocorrelated samples are worth if we were to
convert them to independent samples. In `blavaan` we can print it from
the `summary` function with the `neff` argument

``` r
summary(fit, neff=T)
```

We can also extract only those with the
[`blavInspect()`](http://ecmerkle.github.io/blavaan/reference/blavInspect.md)
function

``` r
blavInspect(fit, "neff")
```

    ##   ind60=~x1   ind60=~x2   ind60=~x3           a           b           c 
    ##    1896.038    1785.782    2100.467    1896.494    2062.721    2002.358 
    ##           d           a           b           c           d dem60~ind60 
    ##    2050.703    1896.494    2062.721    2002.358    2050.703    2177.404 
    ## dem65~ind60 dem65~dem60      y1~~y5      y2~~y4      y2~~y6      y3~~y7 
    ##    2773.603    2734.435    2708.021    3068.299    3050.978    2460.702 
    ##      y4~~y8      y6~~y8      x1~~x1      x2~~x2      x3~~x3      y1~~y1 
    ##    3057.625    2002.876    1955.778    1446.659    3184.749    2521.070 
    ##      y2~~y2      y3~~y3      y4~~y4      y5~~y5      y6~~y6      y7~~y7 
    ##    3680.432    3295.853    3165.311    2621.146    3068.907    2215.842 
    ##      y8~~y8        x1~1        x2~1        x3~1        y1~1        y2~1 
    ##    1977.237    1250.536    1119.443    1245.900    1399.031    1972.669 
    ##        y3~1        y4~1        y5~1        y6~1        y7~1        y8~1 
    ##    1636.792    1431.505    1137.751    1254.328    1179.512    1159.257

ESS is a sample size, so it should be at least 100 (optimally, much more
than 100) times the number of chains in order to be reliable and to
indicate that estimates of the posterior quantiles are reliable. In this
example, because we have 3 chains, we would want to see at least
`neff=300` for every parameter.

And we can easily find the lowest ESS with the
[`min()`](https://rdrr.io/r/base/Extremes.html) function:

``` r
min(blavInspect(fit, "neff"))
```

    ## [1] 1119.443

### References

Bollen, Kenneth A. 1989. *Structural Equations with Latent Variables*.
Wiley Series in Probability and Mathematical Statistics. John Wiley &
Sons, Inc.

Vehtari, Aki, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner. 2021. “Rank-Normalization, Folding, and
Localization: An Improved $\widehat{R}$ for Assessing Convergence of
MCMC (with Discussion).” *Bayesian Analysis* 16 (2): 667–718.
<https://doi.org/10.1214/20-BA1221>.
