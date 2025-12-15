# blavaan Diagnostic Plots

Convenience functions to create plots of blavaan objects, via the
bayesplot package.

## Usage

``` r
# S3 method for class 'blavaan'
plot(x, pars = NULL, plot.type = "trace", showplot = TRUE, ...)
```

## Arguments

- x:

  An object of class `blavaan`.

- pars:

  Parameter numbers to plot, where the numbers correspond to the order
  of parameters as reported by
  [`coef()`](https://rdrr.io/r/stats/coef.html) (also as shown in the
  'free' column of the parTable). If no numbers are provided, all free
  parameters will be plotted.

- plot.type:

  The type of plot desired. This should be the name of a
  [`MCMC`](https://mc-stan.org/bayesplot/reference/MCMC-overview.html)
  function, without the `mcmc_` prefix.

- showplot:

  Should the plot be sent to the graphic device? Defaults to `TRUE`.

- ...:

  Other arguments sent to the bayesplot function.

## Details

In previous versions of blavaan, the plotting functionality was handled
separately for JAGS and for Stan (using plot functionality in packages
runjags and rstan, respectively). For uniformity, all plotting
functionality is now handled by bayesplot. If users desire additional
functionality that is not immediately available, they can extract the
matrix of MCMC draws via `as.matrix(blavInspect(x, 'mcmc'))`.

## Value

An invisible ggplot object that, if desired, can be further customized.

## Examples

``` r
if (FALSE) { # \dontrun{
data(HolzingerSwineford1939, package = "lavaan")

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- bcfa(HS.model, data = HolzingerSwineford1939)

# trace plots of free loadings
plot(fit, pars = 1:6)
} # }
```
