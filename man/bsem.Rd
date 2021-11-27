\name{bsem}
\alias{bsem}
\title{Fit Structural Equation Models}
\description{
Fit a Structural Equation Model (SEM).}
\usage{
bsem(..., cp = "srs",
     dp = NULL, n.chains = 3, burnin, sample,
     adapt, mcmcfile = FALSE, mcmcextra = list(), inits = "prior",
     convergence = "manual", target = "stan", save.lvs = FALSE,
     wiggle = NULL, wiggle.sd = 0.1, prisamp = FALSE, jags.ic = FALSE,
     seed = NULL, bcontrol = list())
}
\arguments{
\item{...}{Default lavaan arguments.  See \code{\link{lavaan}}.}
\item{cp}{Handling of prior distributions on covariance parameters:
  possible values are \code{"srs"} (default) or \code{"fa"}. Option \code{"fa"} is only available for \code{target="jags"}.}
\item{dp}{Default prior distributions on different types of
    parameters, typically the result of a call to \code{dpriors()}.
    See the \code{dpriors()} help file for more information.}
\item{n.chains}{Number of desired MCMC chains.}
\item{burnin}{Number of burnin iterations, NOT including the adaptive
  iterations.}
\item{sample}{The total number of samples to take after burnin.}
\item{adapt}{The number of adaptive iterations to use at the start of
  the simulation.}
\item{mcmcfile}{If \code{TRUE}, the JAGS/Stan model will be written to file
  (in the lavExport directory). Can also supply a character
  string, which serves as the name of the directory to which files will be written.}
\item{mcmcextra}{A list with potential names \code{syntax} (unavailable
  for target=\code{"stan"}),
  \code{monitor}, \code{data}, and \code{llnsamp}. The \code{syntax} object is a text string containing extra
  code to insert in the JAGS/Stan model syntax. The \code{data} object
  is a list of extra data to send to the JAGS/Stan model. The \code{monitor} object
  is a character vector containing extra JAGS/Stan parameters to
  monitor. The \code{llnsamp} object is only relevant to models with ordinal
  variables, and specifies the number of samples that should be drawn to approximate
  the model log-likelihood (larger numbers imply higher accuracy and
  longer time). This log-likelihood is specifically used to compute
  information criteria.}
\item{inits}{If it is a character string, the options are currently
    \code{"simple"}, \code{"Mplus"}, \code{"prior"} (default), or \code{"jags"}.  In the first two
    cases, parameter values are set as though they will be estimated via
    ML (see \code{\link{lavaan}}).  The starting parameter value for
    each chain is then perturbed from the original values through the
    addition of random uniform noise.  If \code{"prior"} is used, the starting
    parameter values are obtained based on the prior distributions
    (while also trying to ensure that the starting values will not crash
    the model estimation).  If \code{"jags"}, no starting values are
    specified and JAGS will choose values on its own (and this will probably
    crash Stan targets). You can also supply
    a list of starting values for each chain, where the list format can
    be obtained from, e.g., \code{blavInspect(fit, "inits")}. Finally,
    you can specify starting values in a similar way to lavaan,
    using the lavaan \code{start} argument (see the lavaan
    documentation for all the options there). In this case, you should also set
    \code{inits="simple"}, and be aware that the same starting values
    will be used for each chain.}
\item{convergence}{Useful only for \code{target="jags"}. If \code{"auto"}, parameters are
  sampled until convergence is achieved (via \code{autorun.jags()}). In
  this case, the arguments \code{burnin} and \code{sample} are passed to
  \code{autorun.jags()} as \code{startburnin} and \code{startsample},
  respectively. Otherwise, parameters
  are sampled as specified by the user (or by the \code{run.jags}
  defaults).}
\item{target}{Desired MCMC sampling, with \code{"stan"} (pre-compiled
  marginal approach) as
  default. Other options include \code{"jags"}, \code{"stancond"}, and
  \code{"stanclassic"}, which sample latent variables and provide some
  greater functionality (because syntax is written "on the fly"). But
  they are slower and less efficient.}
\item{save.lvs}{Should sampled latent variables (factor scores) be saved? Logical; defaults to FALSE}
\item{wiggle}{Labels of equality-constrained parameters that should be
  "approximately" equal. Can also be "intercepts", "loadings",
  "regressions", "means".}
\item{wiggle.sd}{The prior sd (of normal distribution) to be used in approximate equality
  constraints. Can be one value, or (for target="stan") a numeric vector
  of values that is the same length as wiggle.}
\item{prisamp}{Should samples be drawn from the prior, instead of the
  posterior (\code{target="stan"} only)? Logical; defaults to FALSE}
\item{jags.ic}{Should DIC be computed the JAGS way, in addition to the BUGS way? Logical; defaults to FALSE}
\item{seed}{A vector of length \code{n.chains} (for target
  \code{"jags"}) or an integer (for target \code{"stan"}) containing random
  seeds for the MCMC run. If \code{NULL}, seeds will be chosen randomly.}
\item{bcontrol}{A list containing additional parameters passed to
    \code{run.jags} (or \code{autorun.jags}) or \code{stan}.  See the manpage of those functions for an
    overview of the additional parameters that can be set.}
}
\details{
    The \code{bsem} function is a wrapper for the more general
    \code{\link{blavaan}} function, using the following default
    \code{\link{lavaan}} arguments:
    \code{int.ov.free = TRUE}, \code{int.lv.free = FALSE},
    \code{auto.fix.first = TRUE} (unless \code{std.lv = TRUE}),
    \code{auto.fix.single = TRUE}, \code{auto.var = TRUE},
    \code{auto.cov.lv.x = TRUE},
    \code{auto.th = TRUE}, \code{auto.delta = TRUE},
    and \code{auto.cov.y = TRUE}.
  }
\value{
  An object of class \code{\linkS4class{lavaan}}, for which several methods
  are available, including a \code{summary} method.
}
\references{
Edgar C. Merkle, Ellen Fitzsimmons, James Uanhoro, & Ben Goodrich (2021). Efficient Bayesian Structural Equation Modeling in Stan. Journal of Statistical
Software, 100(6), 1-22. URL http://www.jstatsoft.org/v100/i06/.

Edgar C. Merkle & Yves Rosseel (2018). blavaan: Bayesian Structural
Equation Models via Parameter Expansion. Journal of Statistical
Software, 85(4), 1-30. URL http://www.jstatsoft.org/v85/i04/.

Yves Rosseel (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36. URL
http://www.jstatsoft.org/v48/i02/.
}
\seealso{
\code{\link{blavaan}}
}
\examples{
\dontrun{
## The industrialization and Political Democracy Example
## Bollen (1989), page 332
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

## unique priors for mv intercepts; parallel chains
fit <- bsem(model, data=PoliticalDemocracy,
            dp=dpriors(nu="normal(5,10)"))
summary(fit)
}
}
