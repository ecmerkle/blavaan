# Two-level SEM

Starting with version 0.5-1, *blavaan* supports two-level SEM with
random intercepts. The specification and estimation commands are similar
to those of *lavaan*, including use of `level:` in the model
specification and use of the `cluster` argument for estimation.
Consequently, examples involving *lavaan* also generally apply to
*blavaan*, such as the *lavaan* tutorial example below.

``` r
data(Demo.twolevel, package = "lavaan")

model <- '
    level: within
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: between
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
'

bfit <- bsem(model = model, data = Demo.twolevel, cluster = "cluster")
```

Below, we discuss what is currently covered by *blavaan* and some
features that are unique to Bayesian modeling.

### *blavaan* Coverage

As of version 0.5-1, *blavaan* handles two-level, random intercept
models for complete, continuous data. Handling missing data (assuming
missingness at random) will come in a future release. In the meantime,
multiple imputation might be used in combination with the current
*blavaan* functionality (though there is not currently an automatic way
to do it). Alternatively, if there is not much missing data and it
occurs only for lower-level units, listwise deletion could work.

The *blavaan* approach to model estimation mimics the *lavaan* approach,
which uses matrix results (see Rosseel 2021) that enable us to
efficiently evaluate the multilevel SEM likelihood. This will often lead
to more efficient MCMC estimation, as compared to sampling all the level
1 and level 2 latent variables and working with conditional likelihoods
(see Merkle et al. 2021 for discussion of marginal vs conditional
likelihoods).

Similar to single-level models, users can sample latent variables using
the `save.lvs = TRUE` argument in their `bcfa/bsem/bgrowth/blavaan`
commands. Marginal information criteria (marginal over all latent
variables) are also automatically computed, with these information
criteria generally being preferred over those than condition on latent
variables (see Merkle, Furr, and Rabe-Hesketh 2019 for detail in the
context of single-level models).

### Bayes-specific Options

All Bayesian models require prior distributions. The previous *blavaan*
defaults for single-level models are now used for two-level models. You
can continue to use commands like `dpriors(lambda = "normal(1,.5)")` to
specify a Normal(1,.5) prior for all factor loadings and, for two-level
models, that specification will apply to both the level 1 and level 2
loadings. Depending on the model, it may also be useful to specify
priors on individual parameters via the `prior()` argument inside the
model specification syntax. The default prior distributions do not
always work well for observed variables whose values are far from 0. We
continue to encourage users to consider their own prior distributions,
possibly using the `prisamp = TRUE` option to draw samples from the
prior (which could be further used for prior predictive checking).

Model checking also differs between Bayesian and frequentist methods.
Just like it did for one-level models, *blavaan* reports a posterior
predictive p-value for general model assessment. This is computed by
comparing the marginal likelihood of the observed data (marginal over
all latent variables) to the marginal likelihood of artificial data, for
each iteration of MCMC sampling. For finer-grained model assessment, we
encourage users to try
[`ppmc()`](http://ecmerkle.github.io/blavaan/reference/ppmc.md). It
allows you to compute a posterior predictive p-value using your own,
custom model assessment (defined as an R function).

### Concluding Thoughts

We think that the new *blavaan* functionality provides a viable option
for Bayesian two-level SEM, and it should provide a solid base for
future model developments. As always, the underlying Stan files and
supporting data are available via the `mcmcfile = TRUE` argument, and
all the *blavaan* code is available on Github. Bug reports are
appreciated, either at the *blavaan* Google group or as a Github issue.

### References

Merkle, E. C., E. Fitzsimmons, J. Uanhoro, and B. Goodrich. 2021.
“Efficient Bayesian Structural Equation Modeling in Stan.” *Journal of
Statistical Software* 100 (6): 1–22.
<https://www.jstatsoft.org/article/view/v100i06>.

Merkle, E. C., D. Furr, and S. Rabe-Hesketh. 2019. “Bayesian Comparison
of Latent Variable Models: Conditional Versus Marginal Likelihoods.”
*Psychometrika* 84: 802–29. <https://arxiv.org/abs/1802.04452>.

Rosseel, Yves. 2021. “Evaluating the Observed Log-Likelihood Function in
Two-Level Structural Equation Modeling with Missing Data: From Formulas
to R Code.” *Psych* 3 (2): 197–232.
<https://doi.org/10.3390/psych3020017>.
