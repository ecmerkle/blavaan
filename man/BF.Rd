\name{BF}
\alias{BF}
\title{
log-Bayes factor approximation
}
\description{
log-Bayes factor approximation.
}
\usage{
BF(object1, object2, ...)
}
\arguments{
  \item{object1}{An object of class \code{blavaan}.}
  \item{object2}{A second object of class \code{blavaan}.}
  \item{\dots}{Other arguments (unused for now).}
}
\details{
This function approximates the log-Bayes factor of two candidate models using the Laplace approximation to each model's marginal log-likelihood.
}
\value{
The log-Bayes factor approximation, along with each model's approximate marginal log-likelihood.
}
\references{
Raftery, A. E. (1993). Bayesian model selection in structural equation models. In K. A. Bollen & J. S. Long (Eds.), Testing structural equation models (pp. 163-180). Beverly Hills, CA: Sage.
}
\examples{
\dontrun{
hsm1 <- ' visual  =~ x1 + x2 + x3 + x4
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9 '

fit1 <- bcfa(hsm1, data=HolzingerSwineford1939)

hsm2 <- ' visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6 + x7
          speed   =~ x7 + x8 + x9 '

fit2 <- bcfa(hsm2, data=HolzingerSwineford1939)

BF(fit1, fit2)
}
}
