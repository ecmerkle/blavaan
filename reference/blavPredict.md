# Predict the values of latent variables, observed variables, and missing variables.

The purpose of the `blavPredict()` function is to compute various types
of model predictions, conditioned on observed data. This differs
somewhat from
[`lavPredict()`](https://rdrr.io/pkg/lavaan/man/lavPredict.html) in
lavaan.

## Usage

``` r
blavPredict(object, newdata = NULL, type = "lv", level = 1L)
```

## Arguments

- object:

  An object of class
  [`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan-class.md).

- newdata:

  An optional data.frame, containing the same variables as the
  data.frame used when fitting the model in object.

- type:

  A character string. If `"lv"`, estimated values for the latent
  variables in the model are computed. If `"ov"` or `"yhat"`, predicted
  means for the observed variables in the model are computed. If
  `"ypred"` or `"ydist"`, predicted values for the observed variables
  (including residual noise) are computed. If `"ymis"` or `"ovmis"`,
  model predicted values ("imputations") for the missing data are
  computed. See details for further information.

- level:

  For `type = "lv"`, used to specify whether one desires the level 1
  latent variables or level 2 latent variables.

## Details

The [`predict()`](https://rdrr.io/r/stats/predict.html) function calls
the `blavPredict()` function with its default options.

Below, we provide more information about each `type` option. Most
options only work for target="stan", and "number of samples" is defined
as the number of posterior samples across all chains.

`type="lv"`: The posterior distribution of latent variables conditioned
on observed variables. Returns a list with "number of samples" entries,
where each entry is a matrix where rows are observations and columns are
latent variables.

`type="yhat"`: The posterior expected value of observed variables
conditioned on the sampled latent variables. Returns a list with "number
of samples" entries, where each entry is a matrix where rows are
observations and columns are observed variables.

`type="ypred"`: The posterior predictive distribution of observed
variables conditioned on the sampled latent variables (including
residual variability). Returns a list with "number of samples" entries,
where each entry is a data frame where rows are observations and columns
are observed variables.

`type="ymis"`: The posterior predictive distribution of missing values
conditioned on observed variables. Returns a matrix with "number of
samples" rows and "number of missing variables" columns.

## See also

Users may also wish to generate the posterior predictive distribution of
observed data, not conditioned on the latent variables. This would often
be viewed as data from new clusters (people) that were not observed in
the original dataset. For that, see
[`sampleData()`](http://ecmerkle.github.io/blavaan/reference/sampleData.md).

## Examples

``` r
if (FALSE) { # \dontrun{
data(HolzingerSwineford1939, package = "lavaan")

## fit model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- bcfa(HS.model, data = HolzingerSwineford1939, save.lvs = TRUE)
lapply(blavPredict(fit)[1:2], head) # first 6 rows of first 10 posterior samples
head(blavPredict(fit, type = "yhat")[[1]]) # top of first posterior sample

## multigroup models return a list of factor scores (one per group)
mgfit <- bcfa(HS.model, data = HolzingerSwineford1939, group = "school",
              group.equal = c("loadings","intercepts"), save.lvs = TRUE)

lapply(blavPredict(fit)[1:2], head)
head(blavPredict(fit, type = "ypred")[[1]])
} # }
```
