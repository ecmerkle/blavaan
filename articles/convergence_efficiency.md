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
    ##   1.0000416   1.0004711   1.0000193   0.9994823   1.0002393   0.9999769 
    ##           d           a           b           c           d dem60~ind60 
    ##   1.0005137   0.9994823   1.0002393   0.9999769   1.0005137   1.0005405 
    ## dem65~ind60 dem65~dem60      y1~~y5      y2~~y4      y2~~y6      y3~~y7 
    ##   0.9999811   1.0003824   0.9995924   1.0003955   1.0000458   1.0005992 
    ##      y4~~y8      y6~~y8      x1~~x1      x2~~x2      x3~~x3      y1~~y1 
    ##   0.9997992   1.0003337   1.0015656   1.0004742   0.9997060   0.9997108 
    ##      y2~~y2      y3~~y3      y4~~y4      y5~~y5      y6~~y6      y7~~y7 
    ##   1.0006841   1.0010969   1.0001119   1.0003532   0.9996133   1.0007205 
    ##      y8~~y8        x1~1        x2~1        x3~1        y1~1        y2~1 
    ##   0.9994374   1.0019800   1.0016787   1.0020882   1.0007236   1.0004305 
    ##        y3~1        y4~1        y5~1        y6~1        y7~1        y8~1 
    ##   1.0007000   1.0029648   1.0006838   1.0003314   1.0019794   1.0020937

With large models it can be cumbersome to look over all of these
entries. We can instead find the largest $`\hat{R}`$ to see if they are
all less than $`1.05`$

``` r

max(blavInspect(fit, "psrf"))
```

    ## [1] 1.002965

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
    ##    1589.664    1454.148    1719.240    1904.281    1848.949    2306.692 
    ##           d           a           b           c           d dem60~ind60 
    ##    1785.644    1904.281    1848.949    2306.692    1785.644    2512.656 
    ## dem65~ind60 dem65~dem60      y1~~y5      y2~~y4      y2~~y6      y3~~y7 
    ##    2729.397    3136.284    1781.619    2618.778    2880.128    2193.126 
    ##      y4~~y8      y6~~y8      x1~~x1      x2~~x2      x3~~x3      y1~~y1 
    ##    2033.624    1538.547    2106.719    1461.608    3431.093    2098.000 
    ##      y2~~y2      y3~~y3      y4~~y4      y5~~y5      y6~~y6      y7~~y7 
    ##    3083.955    2887.950    2455.078    2010.992    2270.422    2059.278 
    ##      y8~~y8        x1~1        x2~1        x3~1        y1~1        y2~1 
    ##    1705.189    1309.170    1279.029    1363.984    1479.852    1768.054 
    ##        y3~1        y4~1        y5~1        y6~1        y7~1        y8~1 
    ##    1593.021    1372.880    1381.761    1292.985    1230.600    1310.829

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

    ## [1] 1230.6

### References

Bollen, Kenneth A. 1989. *Structural Equations with Latent Variables*.
Wiley Series in Probability and Mathematical Statistics. John Wiley &
Sons, Inc.

Vehtari, Aki, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner. 2021. “Rank-Normalization, Folding, and
Localization: An Improved $`\widehat{R}`$ for Assessing Convergence of
MCMC (with Discussion).” *Bayesian Analysis* 16 (2): 667–718.
<https://doi.org/10.1214/20-BA1221>.
