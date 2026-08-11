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
  single-level models with ordinal variables (currently only \code{thin} and
  \code{parallel} have an effect; see Details).}
}
\details{
For most models, \code{fitMeasures(object, "ppp")} already returns the ppp
computed during estimation, and \code{ppp(object)} simply returns that same
value. The exception is single-level models (with a single group or
multiple groups) that include ordinal variables and were fit with
\code{target = "stan"} or \code{"cmdstan"}: for these models, \code{test}
defaults to \code{"none"} and ppp is \emph{not} computed automatically,
because the computation (a limited-information, pairwise comparison of
observed vs model-implied associations, which better reflects fit for
ordinal data than the previous default) takes noticeably longer than for
other models. Calling \code{ppp(object)} on such a fit computes it as a
separate step. To have ppp computed automatically during estimation instead,
pass \code{test = "standard"} to \code{bcfa()}/\code{bsem()}/\code{bgrowth()}/
\code{blavaan()}.

Posterior predictive checking for two-level (multilevel) models is not
currently supported on demand; for these models, use
\code{test = "standard"} (the default) during estimation to obtain ppp.
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
