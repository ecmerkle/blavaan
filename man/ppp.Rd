\name{ppp}
\alias{ppp}
\title{Posterior Predictive P-value}
\description{
  Obtain the posterior predictive p-value (ppp) for a fitted
  \code{\linkS4class{blavaan}} object, computing it on demand if it was not
  already computed during model estimation.
}
\usage{
ppp(object, ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{blavaan}}.}
\item{...}{Other arguments passed to the underlying computation for
  single-level models with ordinal variables, or for two-level (multilevel)
  models (currently only \code{thin} and \code{parallel} have an effect for
  either case, plus \code{em_control} for two-level models; see Details).}
}
\details{
For most models, \code{fitMeasures(object, "ppp")} already returns the ppp
computed during estimation, and \code{ppp(object)} simply returns that same
value. The exceptions are single-level models (with a single group or
multiple groups) that include ordinal variables, and two-level (multilevel)
models, both fit with \code{target = "stan"} or \code{"cmdstan"}: for these
models, \code{test} defaults to \code{"none"} and ppp is \emph{not} computed
automatically, because the computation takes noticeably longer than for
other models. Calling \code{ppp(object)} on such a fit computes it as a
separate step. To have ppp computed automatically during estimation instead,
pass \code{test = "standard"} (or, for two-level models, any value other than
\code{"none"}) to \code{bcfa()}/\code{bsem()}/\code{bgrowth()}/
\code{blavaan()}.

For single-level ordinal/mixed models, the on-demand computation is a
limited-information, pairwise comparison of observed vs model-implied
associations, which better reflects fit for ordinal data than the previous
Stan-computed default.

For two-level (multilevel) models, the on-demand computation is a
saturated-model likelihood-ratio comparison (the same statistic previously
computed inside the compiled Stan program), now computed in R by calling
\pkg{lavaan}'s own saturated-model EM algorithm once per retained posterior
draw. Because this per-draw computation is comparatively slow, increasing
\code{thin} and/or setting \code{parallel = TRUE} (together with
\code{future::plan("multicore")} or \code{future::plan("multisession")}) is
recommended for larger two-level fits; an \code{em_control} argument (a list
with \code{tol}, \code{max_iter}, and \code{acceleration} elements) is also
available to tune the per-draw EM refit's convergence settings. \code{fixed.x}
variables are supported at the within level, the between level, or both (in
different variables): they are held fixed at their observed values in every
simulated replicate rather than being resimulated, matching how fixed.x is
already handled for single-level on-demand computations. A variable that is
\code{fixed.x} at \emph{both} levels simultaneously is not supported (a
pre-existing \pkg{lavaan} limitation, not specific to this computation) and
will raise an error.
}
\value{
A numeric ppp value (or vector, for multiple imputations/groups, matching
whatever \code{fitMeasures(object, "ppp")} would otherwise return).
}
\seealso{
\code{\link{ppmc}}, for more flexible posterior predictive model checks
using arbitrary discrepancy functions.
}
\examples{
\dontrun{
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
}
}
