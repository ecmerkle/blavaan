# Glance at a blavaan object

Glance at a blavaan object

## Usage

``` r
# S3 method for class 'blavaan'
glance(x, fit.indices = FALSE, ...)
```

## Arguments

- x:

  A `blavaan` object

- fit.indices:

  Character vector of fit indices to compute from `fitMeasures` and
  `blavFitIndices`. Use `"TRUE"` to also compute `blavFitIndices`
  measures. Default is `"FALSE"`.

- ...:

  Additional arguments (currently ignored).

## Value

A single-row `data.frame` with columns:

- npar:

  Number of estimated parameters

- ngroups:

  Number of groups

- ppp:

  Posterior predictive p-value

- bic:

  Bayesian information criterion

- dic:

  Deviance information criterion

- waic:

  Widely applicable information criterion

- looic:

  Leave-one-out information criterion

- margloglik:

  Marginal log-likelihood

Additional EAP of fit measures from `blavFitIndices` are included if
requested.

## Examples

``` r
if (FALSE) { # \dontrun{
library(blavaan)

HS.model <- 'visual =~ x1 + x2 + x3'
fit <- bcfa(HS.model, data = HolzingerSwineford1939, seed = 123)
glance(fit)
glance(fit, fit.indices = TRUE)  # includes BRMSEA
} # }
```
