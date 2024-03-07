\name{blavCompare}
\alias{blavCompare}
\alias{BF}
\title{
Bayesian model comparisons
}
\description{
Bayesian model comparisons, including WAIC, LOO, and Bayes factor approximation.
}
\usage{
blavCompare(object1, object2, ...)
}
\arguments{
  \item{object1}{An object of class \code{blavaan}.}
  \item{object2}{A second object of class \code{blavaan}.}
  \item{\dots}{Other arguments to loo().}
}
\details{
This function computes Bayesian model comparison metrics, including a
Bayes factor approximation, WAIC, and LOOIC.

The log-Bayes factor of the two models is based on the Laplace
approximation to each model's marginal log-likelihood.

The WAIC and LOOIC metrics come from the loo package. The ELPD
difference and SE specifically come from loo::loo_compare().
}
\value{
A list containing separate results for log-Bayes factor, WAIC, LOOIC,
and differences between WAIC and LOOIC.
}
\references{
Raftery, A. E. (1993). Bayesian model selection in structural equation models. In K. A. Bollen & J. S. Long (Eds.), Testing structural equation models (pp. 163-180). Beverly Hills, CA: Sage.

Vehtari A., Gelman A., Gabry J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing, 27, 1413-1432.
}
\examples{
\dontrun{
data(HolzingerSwineford1939, package = "lavaan")

hsm1 <- ' visual  =~ x1 + x2 + x3 + x4
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9 '

fit1 <- bcfa(hsm1, data = HolzingerSwineford1939)

hsm2 <- ' visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6 + x7
          speed   =~ x7 + x8 + x9 '

fit2 <- bcfa(hsm2, data = HolzingerSwineford1939)

blavCompare(fit1, fit2)
}
}
