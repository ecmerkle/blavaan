# Fit Confirmatory Factor Analysis Models

Fit a Confirmatory Factor Analysis (CFA) model.

## Usage

``` r
bcfa(..., cp = "srs",
     dp = NULL, n.chains = 3, burnin, sample,
     adapt, mcmcfile = FALSE, mcmcextra = list(), inits = "simple",
     convergence = "manual", target = "stan", save.lvs = FALSE,
     wiggle = NULL, wiggle.sd = 0.1, prisamp = FALSE, jags.ic = FALSE,
     seed = NULL, bcontrol = list())
```

## Arguments

- ...:

  Default lavaan arguments. See
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html).

- cp:

  Handling of prior distributions on covariance parameters: possible
  values are `"srs"` (default) or `"fa"`. Option `"fa"` is only
  available for `target="jags"`.

- dp:

  Default prior distributions on different types of parameters,
  typically the result of a call to
  [`dpriors()`](http://ecmerkle.github.io/blavaan/reference/dpriors.md).
  See the
  [`dpriors()`](http://ecmerkle.github.io/blavaan/reference/dpriors.md)
  help file for more information.

- n.chains:

  Number of desired MCMC chains.

- burnin:

  Number of burnin/warmup iterations (not including the adaptive
  iterations, for target="jags"). Defaults to 4000 or target="jags" and
  500 for Stan targets.

- sample:

  The total number of samples to take after burnin. Defaults to 10000
  for target="jags" and 1000 for Stan targets.

- adapt:

  For target="jags", the number of adaptive iterations to use at the
  start of sampling. Defaults to 1000.

- mcmcfile:

  If `TRUE`, the JAGS/Stan model will be written to file (in the
  lavExport directory). Can also supply a character string, which serves
  as the name of the directory to which files will be written.

- mcmcextra:

  A list with potential names `syntax` (unavailable for
  target=`"stan"`), `monitor`, `data`, and `llnsamp`. The `syntax`
  object is a text string containing extra code to insert in the
  JAGS/Stan model syntax. The `data` object is a list of extra data to
  send to the JAGS/Stan model. If `moment_match_k_threshold` is
  specified within `data` the looic of the model will be calculated
  using moment matching. The `monitor` object is a character vector
  containing extra JAGS/Stan parameters to monitor. The `llnsamp` object
  is only relevant to models with ordinal variables, and specifies the
  number of samples that should be drawn to approximate the model
  log-likelihood (larger numbers imply higher accuracy and longer time).
  This log-likelihood is specifically used to compute information
  criteria.

- inits:

  If it is a character string, the options are currently `"simple"`
  (default), `"Mplus"`, `"prior"`, or `"jags"`. In the first two cases,
  parameter values are set as though they will be estimated via ML (see
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html)). The starting
  parameter value for each chain is then perturbed from the original
  values through the addition of random uniform noise. If `"prior"` is
  used, the starting parameter values are obtained based on the prior
  distributions (while also trying to ensure that the starting values
  will not crash the model estimation). If `"jags"`, no starting values
  are specified and JAGS will choose values on its own (and this will
  probably crash Stan targets). You can also supply a list of starting
  values for each chain, where the list format can be obtained from,
  e.g., `blavInspect(fit, "inits")`. Finally, you can specify starting
  values in a similar way to lavaan, using the lavaan `start` argument
  (see the lavaan documentation for all the options there). In this
  case, you should also set `inits="simple"`, and be aware that the same
  starting values will be used for each chain.

- convergence:

  Useful only for `target="jags"`. If `"auto"`, parameters are sampled
  until convergence is achieved (via `autorun.jags()`). In this case,
  the arguments `burnin` and `sample` are passed to `autorun.jags()` as
  `startburnin` and `startsample`, respectively. Otherwise, parameters
  are sampled as specified by the user (or by the `run.jags` defaults).

- target:

  Desired MCMC sampling, with `"stan"` (pre-compiled marginal approach)
  as default. Also available is `"vb"`, which calls the rstan function
  `vb()`. Other options include `"jags"`, `"stancond"`, and
  `"stanclassic"`, which sample latent variables and provide some
  greater functionality (because syntax is written "on the fly"). But
  they are slower and less efficient.

- save.lvs:

  Should sampled latent variables (factor scores) be saved? Logical;
  defaults to FALSE

- wiggle:

  Labels of equality-constrained parameters that should be
  "approximately" equal. Can also be "intercepts", "loadings",
  "regressions", "means".

- wiggle.sd:

  The prior sd (of normal distribution) to be used in approximate
  equality constraints. Can be one value, or (for target="stan") a
  numeric vector of values that is the same length as wiggle.

- prisamp:

  Should samples be drawn from the prior, instead of the posterior
  (`target="stan"` only)? Logical; defaults to FALSE

- jags.ic:

  Should DIC be computed the JAGS way, in addition to the BUGS way?
  Logical; defaults to FALSE

- seed:

  A vector of length `n.chains` (for target `"jags"`) or an integer (for
  target `"stan"`) containing random seeds for the MCMC run. If `NULL`,
  seeds will be chosen randomly.

- bcontrol:

  A list containing additional parameters passed to `run.jags` (or
  `autorun.jags`) or `stan`. See the manpage of those functions for an
  overview of the additional parameters that can be set.

## Details

The `bcfa` function is a wrapper for the more general
[`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan.md)
function, using the following default
[`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html) arguments:
`int.ov.free = TRUE`, `int.lv.free = FALSE`, `auto.fix.first = TRUE`
(unless `std.lv = TRUE`), `auto.fix.single = TRUE`, `auto.var = TRUE`,
`auto.cov.lv.x = TRUE`, `auto.th = TRUE`, `auto.delta = TRUE`, and
`auto.cov.y = TRUE`.

## Value

An object that inherits from class
[lavaan](https://rdrr.io/pkg/lavaan/man/lavaan-class.html), for which
several methods are available, including a `summary` method.

## References

Edgar C. Merkle, Ellen Fitzsimmons, James Uanhoro, & Ben Goodrich
(2021). Efficient Bayesian Structural Equation Modeling in Stan. Journal
of Statistical Software, 100(6), 1-22. URL
http://www.jstatsoft.org/v100/i06/.

Edgar C. Merkle & Yves Rosseel (2018). blavaan: Bayesian Structural
Equation Models via Parameter Expansion. Journal of Statistical
Software, 85(4), 1-30. URL http://www.jstatsoft.org/v85/i04/.

Yves Rosseel (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36. URL
http://www.jstatsoft.org/v48/i02/.

## See also

[`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan.md)

## Examples

``` r
data(HolzingerSwineford1939, package = "lavaan")

# The Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

if (FALSE) { # \dontrun{
fit <- bcfa(HS.model, data = HolzingerSwineford1939)
summary(fit)
} # }

# A short run for rough results
fit <- bcfa(HS.model, data = HolzingerSwineford1939, burnin = 100, sample = 100,
            n.chains = 2)
#> 
#> SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.00027 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.7 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 15
#> Chain 1:            adapt_window = 75
#> Chain 1:            term_buffer = 10
#> Chain 1: 
#> Chain 1: Iteration:   1 / 200 [  0%]  (Warmup)
#> Chain 1: Iteration:  20 / 200 [ 10%]  (Warmup)
#> Chain 1: Iteration:  40 / 200 [ 20%]  (Warmup)
#> Chain 1: Iteration:  60 / 200 [ 30%]  (Warmup)
#> Chain 1: Iteration:  80 / 200 [ 40%]  (Warmup)
#> Chain 1: Iteration: 100 / 200 [ 50%]  (Warmup)
#> Chain 1: Iteration: 101 / 200 [ 50%]  (Sampling)
#> Chain 1: Iteration: 120 / 200 [ 60%]  (Sampling)
#> Chain 1: Iteration: 140 / 200 [ 70%]  (Sampling)
#> Chain 1: Iteration: 160 / 200 [ 80%]  (Sampling)
#> Chain 1: Iteration: 180 / 200 [ 90%]  (Sampling)
#> Chain 1: Iteration: 200 / 200 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.364 seconds (Warm-up)
#> Chain 1:                0.304 seconds (Sampling)
#> Chain 1:                0.668 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000198 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.98 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: WARNING: There aren't enough warmup iterations to fit the
#> Chain 2:          three stages of adaptation as currently configured.
#> Chain 2:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 2:          the given number of warmup iterations:
#> Chain 2:            init_buffer = 15
#> Chain 2:            adapt_window = 75
#> Chain 2:            term_buffer = 10
#> Chain 2: 
#> Chain 2: Iteration:   1 / 200 [  0%]  (Warmup)
#> Chain 2: Iteration:  20 / 200 [ 10%]  (Warmup)
#> Chain 2: Iteration:  40 / 200 [ 20%]  (Warmup)
#> Chain 2: Iteration:  60 / 200 [ 30%]  (Warmup)
#> Chain 2: Iteration:  80 / 200 [ 40%]  (Warmup)
#> Chain 2: Iteration: 100 / 200 [ 50%]  (Warmup)
#> Chain 2: Iteration: 101 / 200 [ 50%]  (Sampling)
#> Chain 2: Iteration: 120 / 200 [ 60%]  (Sampling)
#> Chain 2: Iteration: 140 / 200 [ 70%]  (Sampling)
#> Chain 2: Iteration: 160 / 200 [ 80%]  (Sampling)
#> Chain 2: Iteration: 180 / 200 [ 90%]  (Sampling)
#> Chain 2: Iteration: 200 / 200 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.366 seconds (Warm-up)
#> Chain 2:                0.347 seconds (Sampling)
#> Chain 2:                0.713 seconds (Total)
#> Chain 2: 
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Computing post-estimation metrics (including lvs if requested)...
summary(fit)
#> blavaan 0.5.9.1377 ended normally after 100 iterations
#> 
#>   Estimator                                      BAYES
#>   Optimization method                             MCMC
#>   Number of model parameters                        21
#> 
#>   Number of observations                           301
#> 
#>   Statistic                                 MargLogLik         PPP
#>   Value                                      -3806.293       0.000
#> 
#> Parameter Estimates:
#> 
#> 
#> Latent Variables:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>   visual =~                                                                    
#>     x1                1.000                                                    
#>     x2                0.571    0.115    0.372    0.810    1.017    normal(0,10)
#>     x3                0.752    0.122    0.514    1.026    1.007    normal(0,10)
#>   textual =~                                                                   
#>     x4                1.000                                                    
#>     x5                1.117    0.064    0.996    1.237    1.000    normal(0,10)
#>     x6                0.932    0.061    0.822    1.052    0.995    normal(0,10)
#>   speed =~                                                                     
#>     x7                1.000                                                    
#>     x8                1.260    0.197    0.927    1.682    1.008    normal(0,10)
#>     x9                1.186    0.291    0.775    1.947    1.049    normal(0,10)
#> 
#> Covariances:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>   visual ~~                                                                    
#>     textual           0.390    0.079    0.246    0.544    0.997     lkj_corr(1)
#>     speed             0.240    0.051    0.152    0.337    1.013     lkj_corr(1)
#>   textual ~~                                                                   
#>     speed             0.156    0.046    0.078    0.252    0.990     lkj_corr(1)
#> 
#> Variances:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>    .x1                0.568    0.119    0.339    0.802    1.014 gamma(1,.5)[sd]
#>    .x2                1.141    0.109    0.924    1.377    0.996 gamma(1,.5)[sd]
#>    .x3                0.845    0.102    0.671    1.077    1.000 gamma(1,.5)[sd]
#>    .x4                0.375    0.054    0.278    0.483    1.003 gamma(1,.5)[sd]
#>    .x5                0.454    0.059    0.341    0.573    0.993 gamma(1,.5)[sd]
#>    .x6                0.361    0.044    0.276    0.457    0.991 gamma(1,.5)[sd]
#>    .x7                0.834    0.098    0.661    1.027    1.016 gamma(1,.5)[sd]
#>    .x8                0.495    0.122    0.274    0.731    1.069 gamma(1,.5)[sd]
#>    .x9                0.560    0.114    0.297    0.770    1.054 gamma(1,.5)[sd]
#>     visual            0.800    0.142    0.523    1.106    1.013 gamma(1,.5)[sd]
#>     textual           0.988    0.121    0.762    1.261    0.998 gamma(1,.5)[sd]
#>     speed             0.348    0.093    0.163    0.528    1.006 gamma(1,.5)[sd]
#> 
```
