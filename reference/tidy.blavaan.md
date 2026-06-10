# Tidy a blavaan object

This function repackages the output from
[`summary()`](https://rdrr.io/r/base/summary.html) into a friendly
`data.frame`.

## Usage

``` r
# S3 method for class 'blavaan'
tidy(
  x,
  estimate.method = c("mean", "median", "mode"),
  conf.int = TRUE,
  rhat = TRUE,
  ess = TRUE,
  priors = TRUE,
  ...
)
```

## Arguments

- x:

  A `blavaan` object.

- estimate.method:

  The posterior summary statistic to compute the point estimates. One of
  `"mean"` (posterior mean, default), `"median"` (posterior median), or
  `"mode"` (posterior mode). Only free parameters are summarized by the
  requested statistic; fixed parameters retain their fixed values.
  Defined/constrained parameters (`:=`) keep the point estimate from the
  standard [`summary()`](https://rdrr.io/r/base/summary.html) output
  regardless of `estimate.method`.

- conf.int:

  Logical indicating whether to include credible intervals. Default is
  `TRUE`.

- rhat:

  Logical indicating whether to include the Rhat convergence diagnostic.
  Default is `TRUE`.

- ess:

  Logical indicating whether to include the effective sample size.
  Default is `TRUE`.

- priors:

  Logical indicating whether to include the prior distribution
  specification. Default is `TRUE`.

- ...:

  Additional arguments (currently ignored).

## Value

A `data.frame` with columns:

- term:

  Parameter name (lhs, op, rhs combined)

- op:

  Operator from the model syntax

- level:

  The level of a parameter estimate

- group:

  Group number (for multigroup models)

- estimate:

  Posterior summary statistic determined by `estimate.method`

- std.error:

  Posterior standard deviation

- conf.low:

  Lower bound of 95% credible interval (if `conf.int = TRUE`)

- conf.high:

  Upper bound of 95% credible interval (if `conf.int = TRUE`)

- std.lv:

  Standardized estimates based on the variances of the (continuous)
  latent variables only

- std.all:

  Standardized estimates based on both the variances of both
  (continuous) observed and latent variables.

- rhat:

  Rhat convergence diagnostic (if `rhat = TRUE`)

- ess:

  Effective sample size (if `ess = TRUE`)

- prior:

  Prior distribution specification (if `priors = TRUE`)

## Examples

``` r
if (FALSE) { # \dontrun{
data(HolzingerSwineford1939, package = "lavaan")

HS.model <- 'visual =~ x1 + x2 + x3'
fit <- bcfa(HS.model, data = HolzingerSwineford1939, seed = 123,
            n.chains = 1, sample = 300)
tidy(fit)
tidy(fit, estimate.method = "median")

data(Demo.twolevel, package = "lavaan")
model <- "
    level: within
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: between
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
"
bfit <- bsem(
  model = model,
  data = Demo.twolevel,
  cluster = "cluster",
  seed = 123,
  n.chains = 1,
  sample = 300,
  target = "stan"
)
tidy(fit)
} # }
```
