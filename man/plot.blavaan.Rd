\name{plot.blavaan}
\alias{plot.blavaan}
\title{blavaan Diagnostic Plots}
\description{
Convenience functions to create plots of blavaan objects, via the
bayesplot package.
}
\usage{
\method{plot}{blavaan}(x, pars = NULL, plot.type = "trace", showplot = TRUE, ...)
}
\arguments{
  \item{x}{An object of class \code{blavaan}.}
  \item{pars}{Parameter numbers to plot, where the numbers correspond to the order of parameters as reported by \code{coef()} (also as shown in the 'free' column of the parTable). If no numbers are provided, all free parameters will be plotted.}
  \item{plot.type}{The type of plot desired. This should be the name of a \code{\link[bayesplot:MCMC-overview]{MCMC}} function, without the \code{mcmc_} prefix.}
  \item{showplot}{Should the plot be sent to the graphic device? Defaults to \code{TRUE}.}
  \item{...}{Other arguments sent to the bayesplot function.}
}
\details{In previous versions of blavaan, the plotting functionality was
  handled separately for JAGS and for Stan (using plot functionality in
  packages runjags and rstan, respectively). For uniformity, all
  plotting functionality is now handled by bayesplot. If users desire
  additional functionality that is not immediately available, they can extract the matrix of MCMC draws via \code{as.matrix(blavInspect(x, 'mcmc'))}.
}
\value{
An invisible ggplot object that, if desired, can be further customized.
}
\examples{
\dontrun{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- bcfa(HS.model, data=HolzingerSwineford1939)

# trace plots of free loadings
plot(fit, pars = 1:6)
}
}
