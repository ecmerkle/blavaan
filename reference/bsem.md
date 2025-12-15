# Fit Structural Equation Models

Fit a Structural Equation Model (SEM).

## Usage

``` r
bsem(..., cp = "srs",
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

The `bsem` function is a wrapper for the more general
[`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan.md)
function, using the following default
[`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html) arguments:
`int.ov.free = TRUE`, `int.lv.free = FALSE`, `auto.fix.first = TRUE`
(unless `std.lv = TRUE`), `auto.fix.single = TRUE`, `auto.var = TRUE`,
`auto.cov.lv.x = TRUE`, `auto.th = TRUE`, `auto.delta = TRUE`, and
`auto.cov.y = TRUE`.

## Value

An object of class
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
# The industrialization and Political Democracy Example
# Bollen (1989), page 332
data(PoliticalDemocracy, package = "lavaan")

model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

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

if (FALSE) { # \dontrun{
# mildly informative priors for mv intercepts and loadings
fit <- bsem(model, data = PoliticalDemocracy,
            dp = dpriors(nu = "normal(5,10)", lambda = "normal(1,.5)"))
summary(fit)
} # }

# A short run for rough results
fit <- bsem(model, data = PoliticalDemocracy, burnin = 100, sample = 100,
            dp = dpriors(nu = "normal(5,10)", lambda = "normal(1,.5)"),
            n.chains = 2)
#> 
#> SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000229 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.29 seconds.
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
#> Chain 1:  Elapsed Time: 0.497 seconds (Warm-up)
#> Chain 1:                0.373 seconds (Sampling)
#> Chain 1:                0.87 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000205 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 2.05 seconds.
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
#> Chain 2:  Elapsed Time: 0.517 seconds (Warm-up)
#> Chain 2:                0.444 seconds (Sampling)
#> Chain 2:                0.961 seconds (Total)
#> Chain 2: 
#> Warning: The largest R-hat is 1.1, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Computing post-estimation metrics (including lvs if requested)...
#> Warning: blavaan WARNING: As specified, the theta covariance matrix is neither diagonal nor unrestricted, so the actual prior might differ from the stated prior. See
#>  https://arxiv.org/abs/2301.08667
summary(fit)
#> blavaan 0.5.9.1376 ended normally after 100 iterations
#> 
#>   Estimator                                      BAYES
#>   Optimization method                             MCMC
#>   Number of model parameters                        31
#>   Number of equality constraints                     3
#> 
#>   Number of observations                            75
#> 
#>   Statistic                                 MargLogLik         PPP
#>   Value                                             NA       0.540
#> 
#> Parameter Estimates:
#> 
#> 
#> Latent Variables:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>   ind60 =~                                                                     
#>     x1                1.000                                                    
#>     x2                2.105    0.143    1.836    2.381    1.005    normal(1,.5)
#>     x3                1.730    0.149    1.452    2.003    0.994    normal(1,.5)
#>   dem60 =~                                                                     
#>     y1                1.000                                                    
#>     y2         (a)    1.175    0.142    0.940    1.436    0.995    normal(1,.5)
#>     y3         (b)    1.159    0.107    0.965    1.352    0.995    normal(1,.5)
#>     y4         (c)    1.236    0.121    1.054    1.536    0.995    normal(1,.5)
#>   dem65 =~                                                                     
#>     y5                1.000                                                    
#>     y6         (a)    1.175    0.142    0.940    1.436    0.995                
#>     y7         (b)    1.159    0.107    0.965    1.352    0.995                
#>     y8         (c)    1.236    0.121    1.054    1.536    0.995                
#> 
#> Regressions:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>   dem60 ~                                                                      
#>     ind60             1.452    0.422    0.705    2.248    1.002    normal(0,10)
#>   dem65 ~                                                                      
#>     ind60             0.596    0.218    0.156    1.035    0.997    normal(0,10)
#>     dem60             0.863    0.072    0.741    1.019    1.020    normal(0,10)
#> 
#> Covariances:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>  .y1 ~~                                                                        
#>    .y5                0.596    0.399   -0.065    1.488    0.998     lkj_corr(1)
#>  .y2 ~~                                                                        
#>    .y4                1.452    0.729    0.181    3.133    0.993       beta(1,1)
#>    .y6                2.180    0.772    0.871    3.806    0.999       beta(1,1)
#>  .y3 ~~                                                                        
#>    .y7                0.827    0.607   -0.270    2.023    0.995     lkj_corr(1)
#>  .y4 ~~                                                                        
#>    .y8                0.427    0.478   -0.488    1.423    0.995       beta(1,1)
#>  .y6 ~~                                                                        
#>    .y8                1.365    0.589    0.192    2.390    0.995       beta(1,1)
#> 
#> Variances:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>    .x1                0.085    0.022    0.047    0.134    1.003 gamma(1,.5)[sd]
#>    .x2                0.143    0.080    0.010    0.327    1.012 gamma(1,.5)[sd]
#>    .x3                0.501    0.109    0.327    0.760    0.991 gamma(1,.5)[sd]
#>    .y1                1.952    0.490    1.076    2.976    0.996 gamma(1,.5)[sd]
#>    .y2                7.966    1.431    5.555   10.806    0.992 gamma(1,.5)[sd]
#>    .y3                5.339    1.180    3.524    8.495    0.998 gamma(1,.5)[sd]
#>    .y4                3.457    0.811    2.013    5.236    0.998 gamma(1,.5)[sd]
#>    .y5                2.397    0.495    1.571    3.695    0.995 gamma(1,.5)[sd]
#>    .y6                5.214    1.003    3.517    7.176    1.020 gamma(1,.5)[sd]
#>    .y7                3.784    0.879    2.406    5.658    1.001 gamma(1,.5)[sd]
#>    .y8                3.497    0.839    2.200    5.148    0.994 gamma(1,.5)[sd]
#>     ind60             0.489    0.093    0.342    0.730    0.995 gamma(1,.5)[sd]
#>    .dem60             4.157    1.114    2.664    6.484    0.994 gamma(1,.5)[sd]
#>    .dem65             0.229    0.196    0.002    0.698    1.063 gamma(1,.5)[sd]
#> 
```
