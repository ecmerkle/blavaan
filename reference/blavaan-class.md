# Class For Representing A (Fitted) Bayesian Latent Variable Model

The `blavaan` class contains the `lavaan` class, representing a (fitted)
Bayesian latent variable model. It contains a description of the model
as specified by the user, a summary of the data, an internal matrix
representation, and if the model was fitted, the fitting results.

## Objects from the Class

Objects can be created via the
[`bcfa`](http://ecmerkle.github.io/blavaan/reference/bcfa.md),
[`bsem`](http://ecmerkle.github.io/blavaan/reference/bsem.md),
[`bgrowth`](http://ecmerkle.github.io/blavaan/reference/bgrowth.md) or
[`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan.md)
functions.

## Slots

- `version`::

  The lavaan package version used to create this objects

- `call`::

  The function call as returned by
  [`match.call()`](https://rdrr.io/r/base/match.call.html).

- `timing`::

  The elapsed time (user+system) for various parts of the program as a
  list, including the total time.

- `Options`::

  Named list of options that were provided by the user, or filled-in
  automatically.

- `ParTable`::

  Named list describing the model parameters. Can be coerced to a
  data.frame. In the documentation, this is called the \`parameter
  table'.

- `pta`::

  Named list containing parameter table attributes.

- `Data`::

  Object of internal class `"Data"`: information about the data.

- `SampleStats`::

  Object of internal class `"SampleStats"`: sample statistics

- `Model`::

  Object of internal class `"Model"`: the internal (matrix)
  representation of the model

- `Cache`::

  List using objects that we try to compute only once, and reuse many
  times.

- `Fit`::

  Object of internal class `"Fit"`: the results of fitting the model. No
  longer used.

- `boot`::

  List. Unused for Bayesian models.

- `optim`::

  List. Information about the optimization.

- `loglik`::

  List. Information about the loglikelihood of the model (if maximum
  likelihood was used).

- `implied`::

  List. Model implied statistics.

- `vcov`::

  List. Information about the variance matrix (vcov) of the model
  parameters.

- `test`::

  List. Different test statistics.

- `h1`::

  List. Information about the unrestricted h1 model (if available).

- `baseline`::

  List. Information about a baseline model (often the independence
  model) (if available).

- `external`::

  List. Includes Stan or JAGS objects used for MCMC.

## Methods

- coef:

  `signature(object = "blavaan", type = "free")`: Returns the estimates
  of the parameters in the model as a named numeric vector. If
  `type="free"`, only the free parameters are returned. If
  `type="user"`, all parameters listed in the parameter table are
  returned, including constrained and fixed parameters.

- vcov:

  `signature(object = "lavaan")`: returns the covariance matrix of the
  estimated parameters.

- show:

  `signature(object = "blavaan")`: Print a short summary of the model
  fit

- summary:

  `signature(object = "blavaan", header = TRUE, fit.measures = FALSE, estimates = TRUE, ci = TRUE, standardized = FALSE, rsquare = FALSE, std.nox = FALSE, psrf = TRUE, neff = FALSE, postmedian = FALSE, postmode = FALSE, priors = TRUE, bf = FALSE, nd = 3L)`:
  Print a nice summary of the model estimates. If `header = TRUE`, the
  header section (including fit measures) is printed. If
  `fit.measures = TRUE`, additional fit measures are added to the header
  section. If `estimates = TRUE`, print the parameter estimates section.
  If `ci = TRUE`, add confidence intervals to the parameter estimates
  section. If `standardized = TRUE`, the standardized solution is also
  printed. Note that *SE*s and tests are still based on unstandardized
  estimates. Use
  [`standardizedSolution`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
  to obtain *SE*s and test statistics for standardized estimates. If
  `rsquare=TRUE`, the R-Square values for the dependent variables in the
  model are printed. If `std.nox = TRUE`, the `std.all` column contains
  the the `std.nox` column from the parameterEstimates() output. If
  `psrf = TRUE`, potential scale reduction factors (Rhats) are printed.
  If `neff = TRUE`, effective sample sizes are printed. If `postmedian`
  or `postmode` are TRUE, posterior medians or modes are printed instead
  of posterior means. If `priors = TRUE`, parameter prior distributions
  are printed. If `bf = TRUE`, Savage-Dickey approximations of the Bayes
  factor are printed for certain parameters. Nothing is returned (use
  `lavInspect` or another extractor function to extract information from
  a fitted model).

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

[`bcfa`](http://ecmerkle.github.io/blavaan/reference/bcfa.md),
[`bsem`](http://ecmerkle.github.io/blavaan/reference/bsem.md),
[`bgrowth`](http://ecmerkle.github.io/blavaan/reference/bgrowth.md)

## Examples

``` r
if (FALSE) { # \dontrun{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- bcfa(HS.model, data=HolzingerSwineford1939)

summary(fit, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)
coef(fit)
} # }
```
