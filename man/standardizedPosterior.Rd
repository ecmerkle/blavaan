\name{standardizedPosterior}
\alias{standardizedPosterior}
\alias{standardizedposterior}
\title{Standardized Posterior}
\description{
Standardized posterior distribution of a latent variable model.}
\usage{
standardizedPosterior(object, ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{blavaan}}.}
\item{...}{Additional arguments passed to lavaan's
  \code{standardizedSolution()}}
}
\note{
The only allowed \code{standardizedSolution()} arguments are type,
cov.std, remove.eq, remove.ineq, and remove.def. Other arguments are not
immediately suited to posterior distributions.
}
\value{
  A matrix containing standardized posterior draws, where rows are draws
  and columns are parameters.
}
\examples{
\dontrun{
model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

fit <- bsem(model, data=PoliticalDemocracy,
            dp=dpriors(nu="dnorm(5,1e-2)"),
            bcontrol=list(method="rjparallel"))

standardizedPosterior(fit)
}
}
