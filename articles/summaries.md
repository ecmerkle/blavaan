# Model Summaries

Say that we specify a model of the Bollen political democracy data and
draw posterior samples using the following *blavaan* code (where
`save.lvs` saves the latent variable samples for further use):

``` r
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

fit <- bsem(model, data=PoliticalDemocracy, save.lvs = TRUE)
```

We describe here how to summarize the fitted model. The most obvious
functions are [`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html), and
[`vcov()`](https://rdrr.io/r/stats/vcov.html), which all work in a
manner similar to the analogous *lavaan* functions. But instead of
maximum likelihood estimates and standard errors, *blavaan* reports
posterior means and posterior standard deviations. Other summaries that
are unique to Bayesian models include model convergence metrics, model
fit/comparison metrics, and samples of latent variables. These are
discussed below.

### Convergence

Following model estimation, we immediately wish to look at the
“goodness” of the posterior samples, including convergence to a
stationary distribution and autocorrelation. Popular convergence metrics
are available via the
[`blavInspect()`](http://ecmerkle.github.io/blavaan/reference/blavInspect.md)
function:

``` r
blavInspect(fit, 'rhat')
blavInspect(fit, 'neff')
```

where R-hat values near 1.00 indicate convergence, and large effective
sample sizes (hundreds or above) are preferred. For details on these
metrics, see, e.g., the Posterior Analysis section of the [Stan
Reference
Manual](https://mc-stan.org/docs/2_28/reference-manual/index.html).

If the model has definitely not converged (as judged by Rhat), blavaan
will issue multiple warnings. Lack of convergence is sometimes caused by
bad initial values or by a chain that strays to an extreme region of the
posterior space. In these cases, it can be helpful to re-estimate the
model a second time. It is also helpful to specify mildly-informative
priors on loading parameters, so that the chains do not wander to
extreme loading values. For example, if you expect all your variables to
be positively correlated and some loadings are being fixed to 1 for
identification, then Normal(1,.5) would often be a mildly-informative
prior. Otherwise, lack of convergence may imply prior distributions that
severely conflict with the data, or an ill-defined model. It is
sometimes helpful to try to fit the same model in *lavaan*, to observe
whether errors occur there.

### Model Fit & Comparison

Next, we may wish to examine some model fit metrics. While many metrics
are available from the
[`summary()`](https://rdrr.io/r/base/summary.html) output, more are
available from the
[`fitMeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html)
function:

``` r
summary(fit)
fitMeasures(fit)
```

For judging absolute fit, *blavaan* supplies a posterior predictive
p-value that is based on the likelihood ratio statistic. Good-fitting
models have values near 0.5 on this metric. For examining models’
relative fits, *blavaan* supplies the DIC, WAIC, and LOOIC. The latter
two metrics are computed with the help of the *loo* package (Vehtari et
al. 2020). Comparison of multiple models on these criteria is
facilitated via
[`blavCompare()`](http://ecmerkle.github.io/blavaan/reference/blavCompare.md),
which provides standard errors of the difference between two criteria.

Other notable functions include
[`blavFitIndices()`](http://ecmerkle.github.io/blavaan/reference/blavFitIndices.md)
for alternative measures of absolute fit and
[`ppmc()`](http://ecmerkle.github.io/blavaan/reference/ppmc.md) for
general posterior predictive checks.

### Latent Variables & Standardization

An often-discussed advantage of Bayesian models is their abilities to
describe uncertainty in “random” parameters, including random effects
and latent variables. To access this functionality in *blavaan*, users
must set `save.lvs = TRUE` during model estimation, as is done at the
top of this page. After model estimation, uses can access this
information via
[`blavInspect()`](http://ecmerkle.github.io/blavaan/reference/blavInspect.md)
or
[`blavPredict()`](http://ecmerkle.github.io/blavaan/reference/blavPredict.md).
Relevant arguments to
[`blavInspect()`](http://ecmerkle.github.io/blavaan/reference/blavInspect.md)
include `lvmeans` and `lvs`. The former returns posterior means of
latent variables, which are similar to the predictions supplied by
frequentist models. The latter returns posterior samples of latent
variables, so that users could summarize their uncertainties or other
functions of latent variables. These posterior samples are returned as a
list of length `n.chains`, where each list entry has a row per posterior
sample (and number of columns is total number of latent variables in the
model):

``` r
postmns <- blavInspect(fit, what = "lvmeans")
postsamps <- blavInspect(fit, what = "lvs")
```

Some related, but different, information can be obtained by
[`blavPredict()`](http://ecmerkle.github.io/blavaan/reference/blavPredict.md).
This function will also return posterior samples of latent variables,
but in a matrix instead of a list:

``` r
postsamps <- blavPredict(fit, type = "lv")
```

The
[`blavPredict()`](http://ecmerkle.github.io/blavaan/reference/blavPredict.md)
function will also return predictions of observed variables conditioned
on the sampled latent variables. The `type = "yhat"` argument returns
expected values of observed variables conditioned on latent variable
samples; the `type = "ypred"` argument returns posterior predictions of
observed variables including residual noise (essentially `yhat` +
error); and the `type = "ymis"` argument returns posterior predictions
of missing variables conditioned on observed. These expected values and
predictions are returned in list format; for a matrix, see the last line
of code below.

``` r
evpreds <- blavPredict(fit, type = "yhat")
postpreds <- blavPredict(fit, type = "ypred")
mispreds <- blavPredict(fit, type = "ymis")

## convert to matrix from list:
evpreds <- do.call("rbind", evpreds)
```

Finally, not fully related to latent variables: the
[`standardizedPosterior()`](http://ecmerkle.github.io/blavaan/reference/standardizedPosterior.md)
function will return standardized posterior draws. It calls the *lavaan*
function
[`standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
in the background and has some of that function’s flexibility.

### References

Vehtari, Aki, Jonah Gabry, Mans Magnusson, Yuling Yao, Paul-Christian
Bürkner, Topi Paananen, and Andrew Gelman. 2020. “Loo: Efficient
Leave-One-Out Cross-Validation and WAIC for Bayesian Models.”
<https://mc-stan.org/loo/>.
