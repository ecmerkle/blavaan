# Sample data from the posterior (or prior) distribution.

The purpose of the `sampleData()` function is to simulate new data from
a model that has already been estimated. This can faciliate posterior
predictive checks, as well as prior predictive checks (setting prisamp =
TRUE during model estimation).

## Usage

``` r
sampleData(object, nrep = NULL, conditional = FALSE, type = "response",
           simplify = FALSE, ...)
```

## Arguments

- object:

  An object of class
  [`blavaan`](http://ecmerkle.github.io/blavaan/reference/blavaan-class.md).

- nrep:

  How many datasets to generate? If not supplied, defaults to the total
  number of posterior samples.

- conditional:

  Logical indicating whether to sample from the distribution that is
  marginal over latent variables (`FALSE`; default) or from the
  distribution that conditions on latent variables (`TRUE`). For `TRUE`,
  you must set `save.lvs = TRUE` during model estimation.

- type:

  The type of data desired (only relevant to ordinal data). The
  `type = "response"` option generates ordinal data. The `type = "link"`
  option generates continuous variables underlying ordinal data (which
  would be cut by thresholds to yield ordinal data).

- simplify:

  For single-group models, should the list structure be simplified? This
  makes each dataset a single list entry, instead of a list within a
  list (which reflects group 1 of dataset 1). Defaults to `FALSE`.

- ...:

  Other arguments, which for now is only `parallel`. Parallelization via
  `future_lapply()` is available by setting `parallel = TRUE`.

## Details

This is a convenience function to generate data for posterior or prior
predictive checking. The underlying code is also used to generate data
for posterior predictive p-value computation.

## See also

This function overlaps with
[`blavPredict()`](http://ecmerkle.github.io/blavaan/reference/blavPredict.md).
The
[`blavPredict()`](http://ecmerkle.github.io/blavaan/reference/blavPredict.md)
function is more focused on generating pieces of data conditioned on
other pieces of observed data (i.e., latent variables conditioned on
observed variables; missing variables conditioned on observed
variables). In contrast, the `sampleData()` function is more focused on
generating new data given the sampled model parameters.

## Examples

``` r
if (FALSE) { # \dontrun{
data(HolzingerSwineford1939, package = "lavaan")

## fit model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- bcfa(HS.model, data = HolzingerSwineford1939)

## 1 dataset generated from the posterior
out <- sampleData(fit, nrep = 1)

## nested lists: 1 list entry per nrep.
## then, within a rep, 1 list entry per group
## so our dataset is here:
dim(out[[1]][[1]])

## 1 posterior dataset per posterior sample:
out <- sampleData(fit)

## obtain the data on x1 across reps and summarize:
x1dat <- sapply(out, function(x) x[[1]][,1])
summary( as.numeric(x1dat) )
} # }
```
