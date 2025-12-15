# Posterior Predictive Model Checks

This function allows users to conduct a posterior predictive model check
to assess the global or local fit of a latent variable model using any
discrepancy function that can be applied to a
[lavaan](https://rdrr.io/pkg/lavaan/man/lavaan-class.html) model.

## Usage

``` r
ppmc(object, thin = 1, fit.measures = c("srmr","chisq"), discFUN = NULL,
     conditional = FALSE)

# S4 method for class 'blavPPMC'
summary(object, ...)

# S3 method for class 'ppmc'
summary(object, discFUN, dist = c("obs","sim"),
        central.tendency = c("mean","median","mode"),
        hpd = TRUE, prob = .95, to.data.frame = FALSE, diag = TRUE,
        sort.by = NULL, decreasing = FALSE)

# S3 method for class 'blavPPMC'
plot(x, ..., discFUN, element, central.tendency = "",
     hpd = TRUE, prob = .95, nd = 3)

# S3 method for class 'blavPPMC'
hist(x, ..., discFUN, element, hpd = TRUE, prob = .95,
     printLegend = TRUE, legendArgs = list(x = "topleft"),
     densityArgs = list(), nd = 3)

# S3 method for class 'blavPPMC'
pairs(x, discFUN, horInd = 1:DIM, verInd = 1:DIM,
      printLegend = FALSE, ...)
```

## Arguments

- object,x:

  An object of class
  [`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan-class.md).

- thin:

  Optional `integer` indicating how much to thin each chain. Default is
  `1L`, indicating not to thin the chains in `object`.

- fit.measures:

  `character` vector indicating the names of global discrepancy measures
  returned by
  [`fitMeasures`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html).
  Ignored unless `discFUN` is `NULL`, but users may include
  `fitMeasures` in the `list` of discrepancy functions in `discFUN`. For
  ordinal models, the `"logl"` or `"chisq"` computations are done via
  lavaan.

- discFUN:

  `function`, or a `list` of functions, that can be called on an object
  of class [lavaan](https://rdrr.io/pkg/lavaan/man/lavaan-class.html).
  Each function must return an object whose
  [`mode`](https://rdrr.io/r/base/mode.html) is `numeric`, but may be a
  `vector`, `matrix`, or multidimensional `array`. In the `summary` and
  `plot` methods, `discFUN` is a `character` indicating which
  discrepancy function to summarize.

- conditional:

  `logical` indicating whether or not, during artificial data
  generation, we should condition on the estimated latent variables.
  Requires the model to be estimated with `save.lvs = TRUE`.

- element:

  `numeric` or `character` indicating the index (in each `dim`ension of
  the `discFUN` output, if multiple) to plot.

- horInd,verInd:

  Similar to `element`, but a `numeric` or `character` vector indicating
  the indices of a `matrix` to plot in a scatterplot matrix. If
  `horInd==verInd`, histograms will be plotted in the upper triangle.

- dist:

  `character` indicating whether to summarize the distribution of
  `discFUN` on either the `obs`erved or `sim`ulated data.

- central.tendency:

  `character` indicating which statistics should be used to characterize
  the location of the posterior (predictive) distribution. By default,
  all 3 statistics are returned for the `summary` method, but none for
  the `plot` method. The posterior mean is labeled `EAP` for *expected a
  posteriori* estimate, and the mode is labeled `MAP` for *modal a
  posteriori* estimate.

- hpd:

  `logical` indicating whether to calculate the highest posterior
  density (HPD) credible interval for `discFUN`.

- prob:

  The "confidence" level of the credible interval(s).

- nd:

  The number of digits to print in the scatter`plot`.

- to.data.frame:

  `logical` indicating whether the `summary` of a symmetric
  2-dimensional `matrix` returned by `discFUN` should have its unique
  elements stored in rows of a `data.frame` that can be sorted for
  convenience of identifying large discrepancies. If `discFUN` returns
  an asymmetric 2-dimensional `matrix`, the list of matrices returned by
  the `summary` can also be converted to a `data.frame`.

- diag:

  Passed to [`lower.tri`](https://rdrr.io/r/base/lower.tri.html) if
  `to.data.frame=TRUE`.

- sort.by:

  `character`. If `summary` returns a `data.frame`, it can be sorted by
  this column name using [`order`](https://rdrr.io/r/base/order.html).
  Note that if `discFUN` returns an asymmetric 2-dimensional `matrix`,
  each `data.frame` in the returned `list` will be sorted independently,
  so the rows are unlikely to be consistent across summary statistics.

- decreasing:

  Passed to [`order`](https://rdrr.io/r/base/order.html) if
  `!is.null(sort.by)`.

- ...:

  Additional
  `graphical `[`par`](https://rdrr.io/r/graphics/par.html)`ameters` to
  be passed to
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html).

- printLegend:

  `logical`. If `TRUE` (default), a legend will be printed with the
  histogram

- legendArgs:

  `list` of arguments passed to the
  [`legend`](https://rdrr.io/r/graphics/legend.html) function. The
  default argument is a list placing the legend at the top-left of the
  figure.

- densityArgs:

  `list` of arguments passed to the
  [`density`](https://rdrr.io/r/stats/density.html) function, used to
  obtain densities for the `hist` method.

## Value

An S4 object of class `blavPPMC` consisting of 5 `list` slots:

- `@discFUN`:

  The user-supplied `discFUN`, or the `call` to
  [`fitMeasures`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html) that
  returns `fit.measures`.

- `@dims`:

  The dimensions of the object returned by each `discFUN`.

- `@PPP`:

  The posterior predictive *p* value for each `discFUN` element.

- `@obsDist`:

  The posterior distribution of realize values of `discFUN` applied to
  observed data.

- `@simDist`:

  The posterior predictive distribution of values of `discFUN` applied
  to data simulated from the posterior samples.

The [`summary()`](https://rdrr.io/r/base/summary.html) method returns a
`numeric` vector if `discFUN` returns a scalar, a `data.frame` with one
discrepancy function per row if `discFUN` returns a `numeric` vector,
and a `list` with one summary statistic per element if `discFUN` returns
a `matrix` or multidimensional `array`.

The `plot` and `pairs` methods invisibly return `NULL`, printing a plot
(or scatterplot matrix) to the current device.

The `hist` method invisibly returns a `list` or arguments that can be
passed to the function for which the `list` element is named. Users can
edit the arguments in the list to customize their histograms.

## Author

Terrence D. Jorgensen (University of Amsterdam;
<TJorgensen314@gmail.com>)

## References

Levy, R. (2011). Bayesian data–model fit assessment for structural
equation modeling. *Structural Equation Modeling, 18*(4), 663–685.
doi:10.1080/10705511.2011.607723

## Examples

``` r
 if (FALSE) { # \dontrun{
data(HolzingerSwineford1939, package = "lavaan")

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
## fit single-group model
fit <- bcfa(HS.model, data = HolzingerSwineford1939, 
            n.chains = 2, burnin = 1000, sample = 500)
## fit multigroup model
fitg <- bcfa(HS.model, data = HolzingerSwineford1939,
             n.chains = 2, burnin = 1000, sample = 500, group = "school")


## Use fit.measures as a shortcut for global fitMeasures only
## - Note that indices calculated from the "df" are only appropriate under
##   noninformative priors, such that pD approximates the number of estimated
##   parameters counted under ML estimation; incremental fit indices
##   introduce further complications)

AFIs <- ppmc(fit, thin = 10, fit.measures = c("srmr","chisq","rmsea","cfi"))
summary(AFIs)                 # summarize the whole vector in a data.frame
hist(AFIs, element = "rmsea") # only plot one discrepancy function at a time
plot(AFIs, element = "srmr")


## define a list of custom discrepancy functions
## - (global) fit measures
## - (local) standardized residuals

discFUN <- list(global = function(fit) {
                  fitMeasures(fit, fit.measures = c("cfi","rmsea","srmr","chisq"))
                },
                std.cov.resid = function(fit) lavResiduals(fit, zstat = FALSE,
                                                           summary = FALSE)$cov,
                std.mean.resid = function(fit) lavResiduals(fit, zstat = FALSE,
                                                            summary = FALSE)$mean)
out1g <- ppmc(fit, discFUN = discFUN)

## summarize first discrepancy by default (fit indices)
summary(out1g)
## some model-implied correlations look systematically over/underestimated
summary(out1g, discFUN = "std.cov.resid", central.tendency = "EAP")
hist(out1g, discFUN = "std.cov.resid", element = c(1, 7))
plot(out1g, discFUN = "std.cov.resid", element = c("x1","x7"))
## For ease of investigation, optionally export summary as a data.frame,
## sorted by size of average residual
summary(out1g, discFUN = "std.cov.resid", central.tendency = "EAP",
        to.data.frame = TRUE, sort.by = "EAP")
## or sorted by size of PPP
summary(out1g, discFUN = "std.cov.resid", central.tendency = "EAP",
        to.data.frame = TRUE, sort.by = "PPP_sim_LessThan_obs")

## define a list of custom discrepancy functions for multiple groups
## (return each group's numeric output using a different function)

disc2g <- list(global = function(fit) {
                 fitMeasures(fit, fit.measures = c("cfi","rmsea","mfi","srmr","chisq"))
               },
               cor.resid1 = function(fit) lavResiduals(fit, zstat = FALSE,
                                                       type = "cor.bollen",
                                                       summary = FALSE)[[1]]$cov,
               cor.resid2 = function(fit) lavResiduals(fit, zstat = FALSE,
                                                       type = "cor.bollen",
                                                       summary = FALSE)[[2]]$cov)
out2g <- ppmc(fitg, discFUN = disc2g, thin = 2)
## some residuals look like a bigger problem in one group than another
pairs(out2g, discFUN = "cor.resid1", horInd = 1:3, verInd = 7:9) # group 1
pairs(out2g, discFUN = "cor.resid2", horInd = 1:3, verInd = 7:9) # group 2

## print all to file: must be a LARGE picture. First group 1 ...
png("cor.resid1.png", width = 1600, height = 1200)
pairs(out2g, discFUN = "cor.resid1")
dev.off()
## ... then group 2
png("cor.resid2.png", width = 1600, height = 1200)
pairs(out2g, discFUN = "cor.resid2")
dev.off()
} # }
```
