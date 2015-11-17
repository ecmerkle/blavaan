\name{dpriors}
\alias{dpriors}
\title{Specify default prior distributions}
\description{
Specify "default" prior distributions for classes of model parameters.
}
\usage{
dpriors(nu = "dnorm(0,1e-3)", alpha = "dnorm(0,1e-2)", 
        lambda = "dnorm(0,1e-2)", beta = "dnorm(0,1e-2)", 
        itheta = "dgamma(1,.5)", ipsi = "dgamma(1,.5)", 
        rho = "dbeta(1,1)", ibpsi = "dwish(iden,3)")
}
\arguments{
  \item{nu}{Prior distribution for nu (observed variable intercept) parameters.}
  \item{alpha}{Prior distribution for alpha (latent variable intercept) parameters.}
  \item{lambda}{Prior distribution for lambda (loading) parameters.}
  \item{beta}{Prior distribution for regression parameters.}
  \item{itheta}{Prior distribution for observed variable precision parameters.}
  \item{ipsi}{Prior distribution for latent variable precision parameters.}
  \item{rho}{Prior distribution for correlation parameters (only used under srs approach).}
  \item{ibpsi}{Prior distribution for inverse covariance matrix of
    blocks of latent variables.}
}
\details{User-specified prior distributions for specific parameters
  (using the \code{prior()} operator within the model syntax) always
  override prior distributions set using \code{dpriors()}.}
\value{
A character vector containing the prior distribution for each type of parameter.
}
\references{
Edgar C. Merkle & Yves Rosseel (2015). blavaan: Bayesian Structural
Equation Models via Parameter Expansion.
}
\seealso{
\code{\link{bcfa}}, \code{\link{bsem}}, \code{\link{bgrowth}}
}
\examples{
dpriors(nu = "dunif(0,10)", lambda = "dnorm(0,1e-2) T(0,)", itheta = "dexp(1)")
}
