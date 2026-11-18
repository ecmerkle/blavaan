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

The primary convergence diagnostic is $`\hat{R}`$, which compares the
between- and within-chain samples of model parameters and other
univariate quantities of interest (Vehtari et al. 2021). If chains have
not mixed well (ie, the between- and within-chain estimates don’t
agree), $`\hat{R}`$ is larger than 1. We recommend running at least
three chains by default and only using the posterior samples if
$`\hat{R} < 1.05`$ for all the parameters.

`blavaan` presents the $`\hat{R}`$ reported by the underlying MCMC
program, either Stan or JAGS (Stan by default). We can obtain the
$`\hat{R}`$ from the [`summary()`](https://rdrr.io/r/base/summary.html)
function, and we can also extract it with the
[`blavInspect()`](https://blavaan.org/reference/blavInspect.md) function

``` r

blavInspect(fit, "rhat")
```

    ##   ind60=~x1   ind60=~x2   ind60=~x3           a           b           c 
    ##   1.0013895   1.0036691   1.0028358   0.9997569   1.0000582   1.0002676 
    ##           d           a           b           c           d dem60~ind60 
    ##   1.0008511   0.9997569   1.0000582   1.0002676   1.0008511   1.0011352 
    ## dem65~ind60 dem65~dem60      y1~~y5      y2~~y4      y2~~y6      y3~~y7 
    ##   0.9996525   0.9997125   0.9998416   0.9998988   0.9998335   0.9992729 
    ##      y4~~y8      y6~~y8      x1~~x1      x2~~x2      x3~~x3      y1~~y1 
    ##   0.9993508   0.9999769   1.0013602   1.0025521   0.9995397   0.9991122 
    ##      y2~~y2      y3~~y3      y4~~y4      y5~~y5      y6~~y6      y7~~y7 
    ##   1.0004188   0.9996505   1.0003466   1.0008760   0.9995262   1.0007206 
    ##      y8~~y8        x1~1        x2~1        x3~1        y1~1        y2~1 
    ##   0.9995205   1.0010716   1.0001900   1.0003955   0.9995051   0.9996035 
    ##        y3~1        y4~1        y5~1        y6~1        y7~1        y8~1 
    ##   0.9996716   0.9996339   1.0006826   1.0000459   1.0001920   1.0008223

With large models it can be cumbersome to look over all of these
entries. We can instead find the largest $`\hat{R}`$ to see if they are
all less than $`1.05`$

``` r

max(blavInspect(fit, "psrf"))
```

    ## [1] 1.003669

If all $`\hat{R} < 1.05`$ then we can establish that the MCMC chains
have converged to a stable solution. If the model has not converged, you
might increase the number of `burnin` iterations

``` r

fit <- bsem(model, data=PoliticalDemocracy,
            std.lv=T, meanstructure=T, n.chains=3,
            burnin=1000, sample=1000)
```

and/or change the model priors with the
[`dpriors()`](https://blavaan.org/reference/dpriors.md) function. These
address issues where the model failed to converge due to needing more
iterations or due to a model misspecification (such as bad priors). As a
rule of thumb, we seldom see a model require more than 1,000 burnin
samples in Stan. If your model is not converging after 1,000 burnin
samples, it is likely that the default prior distributions clash with
your data. This can happen, e.g., if your variables contain values in
the 100s or 1000s.

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
[`blavInspect()`](https://blavaan.org/reference/blavInspect.md) function

``` r

blavInspect(fit, "neff")
```

    ##   ind60=~x1   ind60=~x2   ind60=~x3           a           b           c 
    ##    1770.265    1672.734    1820.357    2405.554    2423.187    2554.027 
    ##           d           a           b           c           d dem60~ind60 
    ##    1941.884    2405.554    2423.187    2554.027    1941.884    2748.149 
    ## dem65~ind60 dem65~dem60      y1~~y5      y2~~y4      y2~~y6      y3~~y7 
    ##    3535.890    3151.304    2549.144    2879.880    3530.792    3182.230 
    ##      y4~~y8      y6~~y8      x1~~x1      x2~~x2      x3~~x3      y1~~y1 
    ##    3360.927    1877.001    2439.852    1573.489    4348.182    3277.043 
    ##      y2~~y2      y3~~y3      y4~~y4      y5~~y5      y6~~y6      y7~~y7 
    ##    3724.144    4242.428    3004.572    2076.176    2963.434    2672.514 
    ##      y8~~y8        x1~1        x2~1        x3~1        y1~1        y2~1 
    ##    2002.160    1386.152    1228.174    1233.127    1359.575    1615.111 
    ##        y3~1        y4~1        y5~1        y6~1        y7~1        y8~1 
    ##    1489.897    1172.594    1109.814    1110.313    1187.069    1133.822

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

    ## [1] 1109.814

### References

Bollen, Kenneth A. 1989. *Structural Equations with Latent Variables*.
Wiley Series in Probability and Mathematical Statistics. John Wiley &
Sons, Inc.

Vehtari, Aki, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner. 2021. “Rank-Normalization, Folding, and
Localization: An Improved $`\widehat{R}`$ for Assessing Convergence of
MCMC (with Discussion).” *Bayesian Analysis* 16 (2): 667–718.
<https://doi.org/10.1214/20-BA1221>.
