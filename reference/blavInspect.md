# Inspect or Extract Information from a Fitted blavaan Object

The `blavInspect()` and `blavTech()` functions can be used to
inspect/extract information that is stored inside (or can be computed
from) a fitted blavaan object. This is similar to lavaan's
[`lavInspect()`](https://rdrr.io/pkg/lavaan/man/lavInspect.html)
function.

## Usage

``` r
blavInspect(blavobject, what, ...)

blavTech(blavobject, what, ...)
```

## Arguments

- blavobject:

  An object of class blavaan.

- what:

  Character. What needs to be inspected/extracted? See Details for
  Bayes-specific options, and see
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html)'s
  [`lavInspect()`](https://rdrr.io/pkg/lavaan/man/lavInspect.html) for
  additional options. Note: the `what` argument is not case-sensitive
  (everything is converted to lower case.)

- ...:

  lavaan arguments supplied to
  [`lavInspect()`](https://rdrr.io/pkg/lavaan/man/lavInspect.html); see
  [`lavaan`](https://rdrr.io/pkg/lavaan/man/lavaan.html).

## Details

Below is a list of Bayesian-specific values for the `what` argument;
additional values can be found in the
[`lavInspect()`](https://rdrr.io/pkg/lavaan/man/lavInspect.html)
documentation.

- `"start"`::

  A list of starting values for each chain, unless `inits="jags"` is
  used during model estimation. Aliases: `"starting.values"`, `"inits"`.

- `"rhat"`::

  Each parameter's potential scale reduction factor for convergence
  assessment. Can also use "psrf" instead of "rhat"

- `"ac.10"`::

  Each parameter's estimated lag-10 autocorrelation.

- `"neff"`::

  Each parameters effective sample size, taking into account
  autocorrelation.

- `"mcmc"`::

  An object of class `mcmc` containing the individual parameter draws
  from the MCMC run. Aliases: `"draws"`, `"samples"`.

- `"mcobj"`::

  The underlying run.jags or stan object that resulted from the MCMC
  run.

- `"n.chains"`::

  The number of chains sampled.

- `"cp"`::

  The approach used for estimating covariance parameters (`"srs"` or
  `"fa"`); these are only relevant if using JAGS.

- `"dp"`::

  Default prior distributions used for each type of model parameter.

- `"postmode"`::

  Estimated posterior mode of each free parameter.

- `"postmean"`::

  Estimated posterior mean of each free parameter.

- `"postmedian"`::

  Estimated posterior median of each free parameter.

- `"lvs"`::

  An object of class `mcmc` containing latent variable (factor score)
  draws. In two-level models, use `level = 1` or `level = 2` to specify
  which factor scores you want.

- `"lvmeans"`::

  A matrix of mean factor scores (rows are observations, columns are
  variables). Use the additional `level` argument in the same way.

- `"hpd"`::

  HPD interval of each free parameter. In this case, the `prob` argument
  can be used to specify a number in (0,1) reflecting the desired
  percentage of the interval.

## See also

[`lavInspect`](https://rdrr.io/pkg/lavaan/man/lavInspect.html),
[`bcfa`](http://ecmerkle.github.io/blavaan/reference/bcfa.md),
[`bsem`](http://ecmerkle.github.io/blavaan/reference/bsem.md),
[`bgrowth`](http://ecmerkle.github.io/blavaan/reference/bgrowth.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# The Holzinger and Swineford (1939) example
data(HolzingerSwineford1939, package = "lavaan")

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- bcfa(HS.model, data = HolzingerSwineford1939,
            bcontrol = list(method = "rjparallel"))

# extract information
blavInspect(fit, "psrf")
blavInspect(fit, "hpd", prob = .9)
} # }
```
