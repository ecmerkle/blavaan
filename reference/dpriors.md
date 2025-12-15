# Specify Default Prior Distributions

Specify "default" prior distributions for classes of model parameters.

## Usage

``` r
dpriors(..., target = "stan")
```

## Arguments

- ...:

  Parameter names paired with desired priors (see example below).

- target:

  Are the priors for jags, stan (default), or stanclassic?

## Details

The prior distributions always use JAGS/Stan syntax and
parameterizations. For example, the normal distribution in JAGS is
parameterized via the precision, whereas the normal distribution in Stan
is parameterized via the standard deviation.

User-specified prior distributions for specific parameters (using the
`prior()` operator within the model syntax) always override prior
distributions set using `dpriors()`.

The parameter names are:

- nu: Observed variable intercept parameters.

- alpha: Latent variable intercept parameters.

- lambda: Loading parameters.

- beta: Regression parameters.

- itheta: Observed variable precision parameters.

- ipsi: Latent variable precision parameters.

- rho: Correlation parameters (associated with covariance parameters).

- ibpsi: Inverse covariance matrix of blocks of latent variables (used
  for `target="jags"`).

- tau: Threshold parameters (ordinal data only).

- delta: Delta parameters (ordinal data only).

## Value

A character vector containing the prior distribution for each type of
parameter.

## References

Edgar C. Merkle, Ellen Fitzsimmons, James Uanhoro, & Ben Goodrich
(2021). Efficient Bayesian Structural Equation Modeling in Stan. Journal
of Statistical Software, 100(6), 1-22. URL
http://www.jstatsoft.org/v100/i06/.

Edgar C. Merkle & Yves Rosseel (2018). blavaan: Bayesian Structural
Equation Models via Parameter Expansion. Journal of Statistical
Software, 85(4), 1-30. URL http://www.jstatsoft.org/v85/i04/.

## See also

[`bcfa`](http://ecmerkle.github.io/blavaan/reference/bcfa.md),
[`bsem`](http://ecmerkle.github.io/blavaan/reference/bsem.md),
[`bgrowth`](http://ecmerkle.github.io/blavaan/reference/bgrowth.md)

## Examples

``` r
dpriors(nu = "normal(0,10)", lambda = "normal(0,1)", rho = "beta(3,3)")
#>                nu             alpha            lambda              beta 
#>    "normal(0,10)"    "normal(0,10)"     "normal(0,1)"    "normal(0,10)" 
#>             theta               psi               rho             ibpsi 
#> "gamma(1,.5)[sd]" "gamma(1,.5)[sd]"       "beta(3,3)" "wishart(3,iden)" 
#>               tau 
#>   "normal(0,1.5)" 
```
