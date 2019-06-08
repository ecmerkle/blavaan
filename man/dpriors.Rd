\name{dpriors}
\alias{dpriors}
\title{Specify default prior distributions}
\description{
Specify "default" prior distributions for classes of model parameters.
}
\usage{
dpriors(..., target = "stan")
}
\arguments{
  \item{...}{Parameter names paired with desired priors (see example
    below).}
  \item{target}{Are the priors for jags, stan (default), or stanclassic?}
}
\details{The prior distributions always use JAGS/Stan syntax and parameterizations. For
  example, the normal distribution in JAGS is parameterized via the
  precision, whereas the normal distribution in Stan is parameterized
  via the standard deviation.

  User-specified prior distributions for specific parameters
  (using the \code{prior()} operator within the model syntax) always
  override prior distributions set using \code{dpriors()}.

  The parameter names are:
  \itemize{
  \item{nu: Observed variable intercept parameters.}
  \item{alpha: Latent variable intercept parameters.}
  \item{lambda: Loading parameters.}
  \item{beta: Regression parameters.}
  \item{itheta: Observed variable precision parameters.}
  \item{ipsi: Latent variable precision parameters.}
  \item{rho: Correlation parameters (associated with covariance parameters).}
  \item{ibpsi: Inverse covariance matrix of
    blocks of latent variables (used for \code{target="jags"}).}
  \item{tau: Threshold parameters (ordinal data only).}
  \item{delta: Delta parameters (ordinal data only).}
  }
}
\value{
A character vector containing the prior distribution for each type of parameter.
}
\references{
Edgar C. Merkle & Yves Rosseel (2018). blavaan: Bayesian Structural
Equation Models via Parameter Expansion. Journal of Statistical
Software, 85(4), 1-30. URL http://www.jstatsoft.org/v85/i04/.
}
\seealso{
\code{\link{bcfa}}, \code{\link{bsem}}, \code{\link{bgrowth}}
}
\examples{
dpriors(nu = "normal(0,10)", lambda = "normal(0,1)", rho = "beta(3,3)")
}
