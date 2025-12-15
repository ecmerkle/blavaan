# SEM Fit Indices for Bayesian SEM

This function provides a posterior distribution of some \\\chi^2\\-based
fit indices to assess the global fit of a latent variable model.

## Usage

``` r
blavFitIndices(object, thin = 1L, pD = c("loo","waic","dic"),
               rescale = c("devM","ppmc","mcmc"),
               fit.measures = "all", baseline.model = NULL)

## S4 method for signature 'blavFitIndices'
# S4 method for class 'blavFitIndices'
summary(object, ...)

# S3 method for class 'bfi'
summary(object, central.tendency = c("mean","median","mode"),
        hpd = TRUE, prob = .90)
```

## Arguments

- object:

  An object of class
  [`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan-class.md).

- thin:

  Optional `integer` indicating how much to thin each chain. Default is
  `1L`, indicating not to thin the chains.

- pD:

  `character` indicating from which information criterion returned by
  `fitMeasures(object)` to use the estimated number of parameters. The
  default is from the leave-one-out information criterion (LOO-IC),
  which is most highly recommended by Vehtari et al. (2017).

- rescale:

  `character` indicating the method used to calculate fit indices. If
  `rescale = "devM"` (default), the Bayesian analog of the \\\chi^2\\
  statistic (the deviance evaluated at the posterior mean of the model
  parameters) is approximated by rescaling the deviance at each
  iteration by subtracting the estimated number of parameters. If
  `rescale = "PPMC"`, the deviance at each iteration is rescaled by
  subtracting the deviance of data simulated from the posterior
  predictive distribution (as in posterior predictive model checking;
  see Hoofs et al., 2017). If `rescale = "MCMC"`, the fit measures are
  simply calculated using
  [`fitMeasures`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html) at
  each iteration of the Markov chain(s), based on the model-implied
  moments at that iteration (NOT advised when the model includes
  informative priors, in which case the model's estimated *pD* will
  deviate from the number of parameters used to calculate *df* in
  [`fitMeasures`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html)).

- fit.measures:

  If `"all"`, all fit measures available will be returned. If only a
  single or a few fit measures are specified by name, only those are
  computed and returned. If `rescale = "devM"` or `"PPMC"`, the
  currently available indices are `"BRMSEA"`, `"BGammaHat"`,
  `"adjBGammaHat"`, `"BMc"`, `"BCFI"`, `"BTLI"`, or `"BNFI"`. If
  `rescale = "MCMC"`, the user may request any indices returned by
  [`fitMeasures`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html) for
  objects of class
  [lavaan](https://rdrr.io/pkg/lavaan/man/lavaan-class.html).

- baseline.model:

  If not `NULL`, an object of class
  [`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan-class.md),
  representing a user-specified baseline model. If a `baseline.model` is
  provided, incremental fit indices (BCFI, BTLI, or BNFI) can be
  requested in `fit.measures`. Ignored if `rescale = "MCMC"`.

- ...:

  Additional arguments to the summary method:

- central.tendency:

  Takes values "mean", "median", "mode", indicating which statistics
  should be used to characterize the location of the posterior
  distribution. By default, all 3 statistics are returned. The posterior
  mean is labeled EAP for *expected a posteriori* estimate, and the mode
  is labeled MAP for *modal a posteriori* estimate.

- hpd:

  A `logical` indicating whether to calculate the highest posterior
  density (HPD) credible interval for each fit index (defaults to TRUE).

- prob:

  The "confidence" level of the credible interval(s) (defaults to 0.9).

## Value

An S4 object of class `blavFitIndices` consisting of 2 slots:

- `@details`:

  A `list` containing the choices made by the user (or defaults; e.g.,
  which values of `pD` and `rescale` were set), as well as the posterior
  distribution of the \\\chi^2\\ (deviance) statistic (rescaled, if
  `rescale = "devM"` or `"PPMC"`).

- `@indices`:

  A `list` containing the posterior distribution of each requested
  `fit.measure`.

The [`summary()`](https://rdrr.io/r/base/summary.html) method returns a
`data.frame` containing one row for each requested `fit.measure`, and
columns containing the specified measure(s) of `central.tendency`, the
posterior *SD*, and (if requested) the HPD credible-interval limits.

## Author

Mauricio Garnier-Villareal (Vrije Universiteit Amsterdam; <mgv@pm.me>)

Terrence D. Jorgensen (University of Amsterdam;
<TJorgensen314@gmail.com>)

## References

`rescale = "PPMC"` based on:

Hoofs, H., van de Schoot, R., Jansen, N. W., & Kant, I. (2017).
Evaluating model fit in Bayesian confirmatory factor analysis with large
samples: Simulation study introducing the BRMSEA. *Educational and
Psychological Measurement*. doi:10.1177/0013164417709314

`rescale = "devM"` based on:

Garnier-Villarreal, M., & Jorgensen, T. D. (2020). Adapting Fit Indices
for Bayesian Structural Equation Modeling: Comparison to Maximum
Likelihood. *Psychological Methods*, 25(1), 46–70.
https://doi.org/dx.doi.org/10.1037/met0000224 (See also
<https://osf.io/afkcw/>)

Other references:

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing, 27*(5), 1413–1432. doi:10.1007/s11222-016-9696-4

## Examples

``` r
 if (FALSE) { # \dontrun{
data(HolzingerSwineford1939, package = "lavaan")

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
## fit target model
fit1 <- bcfa(HS.model, data = HolzingerSwineford1939, 
             n.chains = 2, burnin = 1000, sample = 1000)

## fit null model to calculate CFI, TLI, and NFI
null.model <- c(paste0("x", 1:9, " ~~ x", 1:9), paste0("x", 1:9, " ~ 1"))
fit0 <- bcfa(null.model, data = HolzingerSwineford1939, 
             n.chains = 2, burnin = 1000, sample = 1000)

## calculate posterior distributions of fit indices

## The default method mimics fit indices derived from ML estimation
ML <- blavFitIndices(fit1, baseline.model = fit0)
ML
summary(ML)

## other options:

## - use Hoofs et al.'s (2017) PPMC-based method
## - use the estimated number of parameters from WAIC instead of LOO-IC
PPMC <- blavFitIndices(fit1, baseline.model = fit0,
                       pD = "waic", rescale = "PPMC")
## issues a warning about using rescale="PPMC" with N < 1000 (see Hoofs et al.)

## - specify only the desired measures of central tendency
## - specify a different "confidence" level for the credible intervals
summary(PPMC, central.tendency = c("mean","mode"), prob = .95)



## Access the posterior distributions for further investigation
head(distML <- data.frame(ML@indices))

## For example, diagnostic plots using the bayesplot package:

## distinguish chains
nChains <- blavInspect(fit1, "n.chains")
distML$Chain <- rep(1:nChains, each = nrow(distML) / nChains)

library(bayesplot)
mcmc_pairs(distML, pars = c("BRMSEA","BMc","BGammaHat","BCFI","BTLI"),
           diag_fun = "hist")
## Indices are highly correlated across iterations in both chains

## Compare to PPMC method
distPPMC <- data.frame(PPMC@indices)
distPPMC$Chain <- rep(1:nChains, each = nrow(distPPMC) / nChains)
mcmc_pairs(distPPMC, pars = c("BRMSEA","BMc","BGammaHat","BCFI","BTLI"),
           diag_fun = "dens")
## nonlinear relation between BRMSEA, related to the floor effect of BRMSEA
## that Hoofs et al. found for larger (12-indicator) models

} # }
```
