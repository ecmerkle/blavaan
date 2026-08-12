# Posterior Predictive P-value

Obtain the posterior predictive p-value (ppp) for a fitted
[`blavaan`](https://blavaan.org/reference/blavaan-class.md) object,
computing it on demand if it was not already computed during model
estimation.

## Usage

``` r
ppp(object, ...)
```

## Arguments

- object:

  An object of class
  [`blavaan`](https://blavaan.org/reference/blavaan-class.md).

- ...:

  Other arguments passed to the underlying computation for single-level
  models with ordinal variables, or for two-level (multilevel) models
  (currently only `thin` and `parallel` have an effect for either case,
  plus `em_control` for two-level models; see Details).

## Details

For most models, `fitMeasures(object, "ppp")` already returns the ppp
computed during estimation, and `ppp(object)` simply returns that same
value. The exceptions are single-level models (with a single group or
multiple groups) that include ordinal variables, and two-level
(multilevel) models, both fit with `target = "stan"` or `"cmdstan"`: for
these models, `test` defaults to `"none"` and ppp is *not* computed
automatically, because the computation takes noticeably longer than for
other models. Calling `ppp(object)` on such a fit computes it as a
separate step. To have ppp computed automatically during estimation
instead, pass `test = "standard"` (or, for two-level models, any value
other than `"none"`) to
[`bcfa()`](https://blavaan.org/reference/bcfa.md)/[`bsem()`](https://blavaan.org/reference/bsem.md)/[`bgrowth()`](https://blavaan.org/reference/bgrowth.md)/
[`blavaan()`](https://blavaan.org/reference/blavaan.md).

For single-level ordinal/mixed models, the on-demand computation is a
limited-information, pairwise comparison of observed vs model-implied
associations, which better reflects fit for ordinal data than the
previous Stan-computed default.

For two-level (multilevel) models, the on-demand computation is a
saturated-model likelihood-ratio comparison (the same statistic
previously computed inside the compiled Stan program), now computed in R
by calling lavaan's own saturated-model EM algorithm once per retained
posterior draw. Because this per-draw computation is comparatively slow,
increasing `thin` and/or setting `parallel = TRUE` (together with
`future::plan("multicore")` or `future::plan("multisession")`) is
recommended for larger two-level fits; an `em_control` argument (a list
with `tol`, `max_iter`, and `acceleration` elements) is also available
to tune the per-draw EM refit's convergence settings. `fixed.x`
variables are supported at the within level, the between level, or both
(in different variables): they are held fixed at their observed values
in every simulated replicate rather than being resimulated, matching how
fixed.x is already handled for single-level on-demand computations. A
variable that is `fixed.x` at *both* levels simultaneously is not
supported (a pre-existing lavaan limitation, not specific to this
computation) and will raise an error.

## Value

A numeric ppp value (or vector, for multiple imputations/groups,
matching whatever `fitMeasures(object, "ppp")` would otherwise return).

## See also

[`ppmc`](https://blavaan.org/reference/ppmc.md), for more flexible
posterior predictive model checks using arbitrary discrepancy functions.

## Examples

``` r
if (FALSE) { # \dontrun{
data(HolzingerSwineford1939, package = "lavaan")

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

## make some indicators ordinal
HSord <- HolzingerSwineford1939
HSord$x1 <- as.integer(cut(HSord$x1, 3))
HSord$x2 <- as.integer(cut(HSord$x2, 3))

## ppp is not computed by default for this ordinal model
fit <- bcfa(HS.model, data = HSord, ordered = c("x1", "x2"))

## compute it on demand
ppp(fit)
} # }
```
