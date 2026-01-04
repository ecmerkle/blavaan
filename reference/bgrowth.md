# Fit Growth Curve Models

Fit a Growth Curve model.

## Usage

``` r
bgrowth(..., cp = "srs", dp = NULL, n.chains = 3,
burnin, sample, adapt, mcmcfile = FALSE, mcmcextra = list(), 
inits = "simple", convergence = "manual", target = getOption("blavaan.target", "stan"),
save.lvs = FALSE, wiggle = NULL, wiggle.sd = 0.1, prisamp = FALSE,
jags.ic = FALSE, seed = NULL, bcontrol = list())
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

The `bgrowth` function is a wrapper for the more general
[`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan.md)
function, using the following default
[`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html) arguments:
`meanstructure = TRUE`, `int.ov.free = FALSE`, `int.lv.free = TRUE`,
`auto.fix.first = TRUE` (unless `std.lv = TRUE`),
`auto.fix.single = TRUE`, `auto.var = TRUE`, `auto.cov.lv.x = TRUE`,
`auto.th = TRUE`, `auto.delta = TRUE`, and `auto.cov.y = TRUE`.

## Value

An object of class
[`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan.md), for
which several methods are available, including a `summary` method.

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
if (FALSE) { # \dontrun{
## linear growth model with a time-varying covariate
data(Demo.growth, package = "lavaan")

model.syntax <- '
  # intercept and slope with fixed coefficients
    i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
    s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4

  # regressions
    i ~ x1 + x2
    s ~ x1 + x2

  # time-varying covariates
    t1 ~ c1
    t2 ~ c2
    t3 ~ c3
    t4 ~ c4
'

fit <- bgrowth(model.syntax, data = Demo.growth)
summary(fit)
} # }
```
