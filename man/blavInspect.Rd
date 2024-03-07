\name{blavInspect}
\alias{blavInspect}
\alias{blavTech}
\title{Inspect or Extract Information from a Fitted blavaan Object}
\description{
The \code{blavInspect()} and \code{blavTech()} functions can be used to
inspect/extract information that is stored inside (or can be computed from) a
fitted blavaan object. This is similar to lavaan's \code{lavInspect()} function.}
\usage{
blavInspect(blavobject, what, ...)

blavTech(blavobject, what, ...)
}
\arguments{
\item{blavobject}{An object of class blavaan.}
\item{what}{Character. What needs to be inspected/extracted? See Details for Bayes-specific options, and see \code{\link{lavaan}}'s \code{lavInspect()} for additional options. Note: the \code{what} argument is not case-sensitive (everything is converted to lower case.)}
\item{...}{lavaan arguments supplied to \code{lavInspect()}; see \code{\link{lavaan}}.}
}
\details{
Below is a list of Bayesian-specific values for the \code{what}
argument; additional values can be found in the \code{lavInspect()}
documentation.

\describe{
  \item{\code{"start"}:}{A list of starting values for each chain, unless \code{inits="jags"} is used during model estimation. Aliases: \code{"starting.values"}, \code{"inits"}.}
  \item{\code{"rhat"}:}{Each parameter's potential scale reduction
    factor for convergence assessment. Can also use "psrf" instead of "rhat"}
  \item{\code{"ac.10"}:}{Each parameter's estimated lag-10 autocorrelation.}
  \item{\code{"neff"}:}{Each parameters effective sample size, taking into account autocorrelation.}
  \item{\code{"mcmc"}:}{An object of class \code{mcmc} containing the individual parameter draws from the MCMC run. Aliases: \code{"draws"}, \code{"samples"}.}
  \item{\code{"mcobj"}:}{The underlying run.jags or stan object that resulted from the MCMC run.}
  \item{\code{"n.chains"}:}{The number of chains sampled.}
  \item{\code{"cp"}:}{The approach used for estimating covariance
    parameters (\code{"srs"} or \code{"fa"}); these are only relevant if
    using JAGS.}
  \item{\code{"dp"}:}{Default prior distributions used for each type of model parameter.}
  \item{\code{"postmode"}:}{Estimated posterior mode of each free parameter.}
  \item{\code{"postmean"}:}{Estimated posterior mean of each free parameter.}
  \item{\code{"postmedian"}:}{Estimated posterior median of each free parameter.}
  \item{\code{"lvs"}:}{An object of class \code{mcmc} containing latent variable (factor score) draws. In two-level models, use \code{level = 1} or \code{level = 2} to specify which factor scores you want.}
  \item{\code{"lvmeans"}:}{A matrix of mean factor scores (rows are observations, columns are variables). Use the additional \code{level} argument in the same way.}
  \item{\code{"hpd"}:}{HPD interval of each free parameter. In this case, the \code{prob} argument can be used to specify a number in (0,1) reflecting the desired percentage of the interval.}
}
}
\seealso{
\code{\link{lavInspect}}, \code{\link{bcfa}}, \code{\link{bsem}}, \code{\link{bgrowth}}
}
\examples{
\dontrun{
# The Holzinger and Swineford (1939) example
data(HolzingerSwineford1939, package = "lavaan")

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- bcfa(HS.model, data = HolzingerSwineford1939,
            bcontrol = list(method = "rjparallel"))

# extract information
blavInspect(fit, "psrf")
blavInspect(fit, "hpd", prob = .9)
}
}
