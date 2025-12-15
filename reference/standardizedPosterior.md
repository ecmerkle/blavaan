# Standardized Posterior

Standardized posterior distribution of a latent variable model.

## Usage

``` r
standardizedPosterior(object, ...)
```

## Arguments

- object:

  An object of class
  [`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan-class.md).

- ...:

  Additional arguments passed to lavaan's
  [`standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)

## Note

The only allowed
[`standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
arguments are type, cov.std, remove.eq, remove.ineq, and remove.def.
Other arguments are not immediately suited to posterior distributions.

## Value

A matrix containing standardized posterior draws, where rows are draws
and columns are parameters.

## Examples

``` r
if (FALSE) { # \dontrun{
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

fit <- bsem(model, data = PoliticalDemocracy,
            dp = dpriors(nu = "normal(5, 10)"))

standardizedPosterior(fit)
} # }
```
