# Fit Structural Equation Models

Fit a Structural Equation Model (SEM).

## Usage

``` r
bsem(..., cp = "srs",
     dp = NULL, n.chains = 3, burnin, sample,
     adapt, mcmcfile = FALSE, mcmcextra = list(), inits = "simple",
     convergence = "manual", target = getOption("blavaan.target", "stan"),
     save.lvs = FALSE, wiggle = NULL, wiggle.sd = 0.1, prisamp = FALSE,
     jags.ic = FALSE, seed = NULL, bcontrol = list())
```

## Arguments

- ...:

  Default lavaan arguments. See
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html). This also
  accepts sampler-specific arguments that would otherwise need to be
  supplied via `bcontrol` – see `bcontrol` below.

- cp:

  Handling of prior distributions on covariance parameters: possible
  values are `"srs"` (default) or `"fa"`. Option `"fa"` is only
  available for `target="jags"`.

- dp:

  Default prior distributions on different types of parameters,
  typically the result of a call to
  [`dpriors()`](https://blavaan.org/reference/dpriors.md). See the
  [`dpriors()`](https://blavaan.org/reference/dpriors.md) help file for
  more information.

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
  `autorun.jags`), `stan`, `vb`, or (for `target="cmdstan"`) a
  `cmdstanr` model's `$sample()` method. See the manpage of those
  functions for an overview of the additional parameters that can be
  set. Recognized arguments (e.g. `adapt_delta`, `max_treedepth`,
  `cores`, `thin`, `method`, `threads_per_chain`, and other
  target-specific sampler arguments) can instead be supplied directly to
  [`blavaan()`](https://blavaan.org/reference/blavaan.md)/
  [`bcfa()`](https://blavaan.org/reference/bcfa.md)/`bsem()`/[`bgrowth()`](https://blavaan.org/reference/bgrowth.md)
  without wrapping them in `bcontrol`; for
  `target %in% c("stan","stanclassic", "stancond")`, arguments belonging
  to Stan's nested `control=` list (like `adapt_delta`) are detected
  automatically and nested for you. Supplying the same argument both
  directly and inside `bcontrol` is an error. A few names remain
  reserved and must still be supplied via `bcontrol` because they
  collide with another
  [`blavaan()`](https://blavaan.org/reference/blavaan.md)/[`lavaan()`](https://rdrr.io/pkg/lavaan/man/lavaan.html)
  argument: `verbose`, `control` (lavaan's own optimizer control), and
  `cl` (lavaan's bootstrap cluster option).

## Details

The `bsem` function is a wrapper for the more general
[`blavaan`](https://blavaan.org/reference/blavaan.md) function, using
the following default
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

[`blavaan`](https://blavaan.org/reference/blavaan.md)

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
#> Chain 1: Gradient evaluation took 0.000228 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.28 seconds.
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
#> Chain 1:  Elapsed Time: 0.513 seconds (Warm-up)
#> Chain 1:                0.366 seconds (Sampling)
#> Chain 1:                0.879 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'stanmarg' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000199 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.99 seconds.
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
#> Chain 2:  Elapsed Time: 0.349 seconds (Warm-up)
#> Chain 2:                0.315 seconds (Sampling)
#> Chain 2:                0.664 seconds (Total)
#> Chain 2: 
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
#> blavaan 0.6.1.1487 ended normally after 100 iterations
#> 
#>   Estimator                                      BAYES
#>   Optimization method                             MCMC
#>   Number of model parameters                        31
#>   Number of equality constraints                     3
#> 
#>   Number of observations                            75
#> 
#>   Statistic                                 MargLogLik         PPP
#>   Value                                             NA       0.655
#> 
#> Parameter Estimates:
#> 
#> 
#> Latent Variables:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>   ind60 =~                                                                     
#>     x1                1.000                                                    
#>     x2                2.110    0.127    1.877    2.367    0.997    normal(1,.5)
#>     x3                1.721    0.132    1.469    1.962    0.991    normal(1,.5)
#>   dem60 =~                                                                     
#>     y1                1.000                                                    
#>     y2         (a)    1.177    0.136    0.948    1.475    0.994    normal(1,.5)
#>     y3         (b)    1.167    0.108    0.979    1.375    1.001    normal(1,.5)
#>     y4         (c)    1.257    0.127    1.022    1.502    1.002    normal(1,.5)
#>   dem65 =~                                                                     
#>     y5                1.000                                                    
#>     y6         (a)    1.177    0.136    0.948    1.475    0.994                
#>     y7         (b)    1.167    0.108    0.979    1.375    1.001                
#>     y8         (c)    1.257    0.127    1.022    1.502    1.002                
#> 
#> Regressions:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>   dem60 ~                                                                      
#>     ind60             1.362    0.334    0.783    2.035    1.004    normal(0,10)
#>   dem65 ~                                                                      
#>     ind60             0.576    0.238    0.142    1.028    1.020    normal(0,10)
#>     dem60             0.877    0.076    0.738    1.022    0.994    normal(0,10)
#> 
#> Covariances:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>  .y1 ~~                                                                        
#>    .y5                0.622    0.366    0.014    1.446    0.990     lkj_corr(1)
#>  .y2 ~~                                                                        
#>    .y4                1.539    0.721    0.320    2.972    1.002       beta(1,1)
#>    .y6                2.222    0.723    1.096    3.898    1.000       beta(1,1)
#>  .y3 ~~                                                                        
#>    .y7                0.722    0.658   -0.401    1.961    0.999     lkj_corr(1)
#>  .y4 ~~                                                                        
#>    .y8                0.346    0.448   -0.427    1.133    1.024       beta(1,1)
#>  .y6 ~~                                                                        
#>    .y8                1.373    0.614    0.292    2.535    0.996       beta(1,1)
#> 
#> Variances:
#>                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
#>    .x1                0.081    0.019    0.046    0.117    1.000 gamma(1,.5)[sd]
#>    .x2                0.144    0.078    0.025    0.302    1.004 gamma(1,.5)[sd]
#>    .x3                0.500    0.099    0.337    0.710    0.995 gamma(1,.5)[sd]
#>    .y1                2.009    0.511    1.235    2.976    0.996 gamma(1,.5)[sd]
#>    .y2                7.832    1.390    5.468   10.591    0.999 gamma(1,.5)[sd]
#>    .y3                5.160    1.016    3.383    7.594    0.995 gamma(1,.5)[sd]
#>    .y4                3.476    0.771    2.097    5.254    0.992 gamma(1,.5)[sd]
#>    .y5                2.444    0.539    1.555    3.771    0.995 gamma(1,.5)[sd]
#>    .y6                5.250    1.005    3.789    7.639    1.009 gamma(1,.5)[sd]
#>    .y7                3.776    0.795    2.448    5.212    0.994 gamma(1,.5)[sd]
#>    .y8                3.536    0.779    2.350    5.443    1.000 gamma(1,.5)[sd]
#>     ind60             0.491    0.087    0.364    0.656    1.000 gamma(1,.5)[sd]
#>    .dem60             3.903    0.830    2.557    5.440    1.017 gamma(1,.5)[sd]
#>    .dem65             0.165    0.196    0.000    0.763    1.011 gamma(1,.5)[sd]
#> 
```
