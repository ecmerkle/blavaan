# Specifying Prior Distributions

There are two ways to specify prior distributions in blavaan. First,
each type of model parameter has a default prior distribution that may
or may not be suitable for your specific situation. You are free to
modify the defaults. Second, the priors for individual model parameters
can be specified in the model syntax. Each is discussed below.

### Defaults

The default priors can be seen via

``` r
dpriors()
```

    ##                nu             alpha            lambda              beta 
    ##    "normal(0,32)"    "normal(0,10)"    "normal(0,10)"    "normal(0,10)" 
    ##             theta               psi               rho             ibpsi 
    ## "gamma(1,.5)[sd]" "gamma(1,.5)[sd]"       "beta(1,1)" "wishart(3,iden)" 
    ##               tau 
    ##   "normal(0,1.5)"

It is important to note that these prior distributions correspond to
Stan parameterizations. These are similar to R parameterizations but not
necessarily exactly the same. The Greek(ish) names above correspond to
the following parameter types (where MV is manifest/observed variable
and LV is latent variable):

    ##                  nu               alpha              lambda                beta 
    ##      "MV intercept"      "LV intercept"           "Loading"        "Regression" 
    ##               theta                 psi                 rho               ibpsi 
    ##      "MV precision"      "LV precision"       "Correlation" "Covariance matrix" 
    ##                 tau 
    ##         "Threshold"

For further information about priors on thresholds, see the [ordinal
modeling
details](http://ecmerkle.github.io/blavaan/articles/ordinal.md).

For `target = "stan"` (the default), priors are currently restricted to
one distribution per parameter type. You can change the prior
distribution parameters (for example, the mean and standard deviation of
a normal), but you cannot change the prior distribution type. The only
exceptions here are the “theta” and “psi” parameters: for those, you can
use the modifiers “\[sd\]”, “\[var\]”, or “\[prec\]” to specify whether
you want the priors to apply to the standard deviation, variance, or
precision. If you require more flexibility in prior specification, you
change the target to either `"stanclassic"` (the old Stan approach) or
`"jags"` (the JAGS approach). Alternatively, you can export the Stan
model via `mcmcfile = TRUE`, edit the file as needed, then fit it via
the rstan package.

To modify prior distributions, we could simply supply a new text string
to [`dpriors()`](http://ecmerkle.github.io/blavaan/reference/dpriors.md)
like this:

``` r
mydp <- dpriors(lambda="normal(1,2)")
mydp
```

    ##                nu             alpha            lambda              beta 
    ##    "normal(0,32)"    "normal(0,10)"     "normal(1,2)"    "normal(0,10)" 
    ##             theta               psi               rho             ibpsi 
    ## "gamma(1,.5)[sd]" "gamma(1,.5)[sd]"       "beta(1,1)" "wishart(3,iden)" 
    ##               tau 
    ##   "normal(0,1.5)"

so that the default prior for loadings is now normal with mean 1 and
standard deviation 2, and the rest of the parameters remain at the
original defaults. The next time we estimate a model (via
[`bsem()`](http://ecmerkle.github.io/blavaan/reference/bsem.md),
[`bcfa()`](http://ecmerkle.github.io/blavaan/reference/bcfa.md),
[`bgrowth()`](http://ecmerkle.github.io/blavaan/reference/bgrowth.md),
or
[`blavaan()`](http://ecmerkle.github.io/blavaan/reference/blavaan.md)),
we would add the argument `dp=mydp` to use this new set of default
priors.

### Individual Parameters

In addition to setting the prior for one type of model parameter, the
user may wish to set the prior of a specific model parameter. This is
accomplished by using the `prior()` modifier within the model
specification. For example, consider the following syntax for the
Holzinger and Swineford (1939) confirmatory factor model:

``` r
HS.model <- ' visual  =~ x1 + prior("normal(1,2)")*x2 + x3
              textual =~ x4 + x5 + prior("normal(3,1.5)")*x6
              speed   =~ x7 + x8 + x9 
              x1 ~~ prior("gamma(3,3)[sd]")*x1 '
```

The loading from `visual` to `x2` now has a normal prior with mean 1 and
standard deviation 2, while the loading from `textual` to `x6` has a
normal prior with mean 3 and standard deviation 1.5. All other loadings
have the default prior distribution.

In the above syntax, we have additionally specified a gamma(3,3) prior
associated with the residual of `x1`. The `[sd]` text at the end of the
distribution says that this prior goes on the residual standard
deviation, as opposed to the residual precision or residual variance.
There exist two more options here: a `[var]` option for the residual
variance, and no brackets for the precision (or you could also use
`[prec]`). This bracketed text can be used for any model
variance/SD/precision parameter and could also be used in default prior
specification if desired.

### Covariance Parameters

One additional note on covariance parameters defined in the model
syntax: the `prior()` syntax specifies a prior on the correlation
associated with the covariance parameter, as opposed to the covariance
itself. The specified distribution should have support on (0,1), and
blavaan automatically translates the prior to an equivalent distribution
with support on (-1,1). It is safest to stick with beta priors here. For
example, the syntax

``` r
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 
              visual ~~ prior("beta(1,1)")*textual '
```

places a Beta(1,1) (uniform) prior on the correlation between the
`visual` and `textual` factors. If desired, we could also specify priors
on the standard deviations (or variances or precisions) of the `visual`
and `textual` factors. Together with the prior on the correlation, these
priors would imply a prior on the covariance between `visual` and
`textual`.

### References

Holzinger, K. J., and F. A. Swineford. 1939. *A Study of Factor
Analysis: The Stability of a Bi-Factor Solution*. Supplementary
Educational Monograph 48. Chicago: University of Chicago Press.
