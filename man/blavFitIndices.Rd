\name{blavFitIndices}
\alias{blavFitIndices}
\alias{blavFitIndices-class}
\alias{show,blavFitIndices-method}
\alias{summary,blavFitIndices-method}
\title{SEM Fit Indices for Bayesian SEM}
\description{
This function provides a posterior distribution of some \eqn{\chi^2}-based fit
indices to assess the global fit of a latent variable model.}
\usage{
blavFitIndices(object, pD = c("loo","waic","dic"),
               rescale = c("devM","ppmc","mcmc"),
               fit.measures = "all", baseline.model = NULL)

## S4 method for signature 'blavFitIndices'
\S4method{summary}{blavFitIndices}(object, central.tendency = c("mean","median","mode"),
        hpd = TRUE, prob = .90)
}
\arguments{
  \item{object}{An object of class \code{\linkS4class{blavaan}}.}
  \item{pD}{\code{character} indicating from which information criterion
    returned by \code{fitMeasures(object)} to use the estimated number of
    parameters. The default is from the leave-one-out information criterion
    (LOO-IC), which is most highly recommended by Vehtari et al. (2017).}
  \item{rescale}{\code{character} indicating the method used to calculate
    fit indices. If \code{rescale = "devM"} (default), the Bayesian analog of
    the \eqn{\chi^2} statistic (the deviance evaluated at the posterior mean
    of the model parameters) is approximated by rescaling the deviance at each
    iteration by subtracting the estimated number of parameters. If
    \code{rescale = "PPMC"}, the deviance at each iteration is rescaled by
    subtracting the deviance of data simulated from the posterior predictive
    distribution (as in posterior predictive model checking; see Hoofs et al.,
    2017). If \code{rescale = "MCMC"} (not implemented yet), the fit
    measures are simply calculated at each iteration of the Markov chain(s),
    based on the model-implied moments at that iteration.}
  \item{fit.measures}{If \code{"all"}, all fit measures available will be
    returned. If only a single or a few fit measures are specified by name,
    only those are computed and returned. If \code{rescale = "devM"} or
    \code{"PPMC"}, the currently available indices are \code{"BRMSEA"},
    \code{"BGammaHat"}, \code{"adjBGammaHat"}, \code{"BMc"}, \code{"BCFI"},
    \code{"BTLI"}, or \code{"BNFI"}. If \code{rescale = "MCMC"}, the user may
    request any indices returned by \code{\link[lavaan]{fitMeasures}} for
    objects of class \code{\linkS4class{lavaan}}.}
  \item{baseline.model}{If not \code{NULL}, an object of class 
    \code{\linkS4class{blavaan}}, representing a user-specified baseline model.
    If a \code{baseline.model} is provided, incremental fit indices (BCFI,
    BTLI, or BNFI) can be requested in \code{fit.measures}.}
  \item{central.tendency}{\code{character} indicating which statistics should
    be used to characterize the location of the posterior distribution. By
    default, all 3 statistics are returned. The posterior mean is labeled
    \code{EAP} for \emph{expected a posteriori} estimate, and the mode is
    labeled \code{MAP} for \emph{modal a posteriori} estimate.}
  \item{hpd}{\code{logical} indicating whether to calculate the highest
    posterior density (HPD) credible interval for each fit index.}
  \item{prob}{The "confidence" level of the credible interval(s).}
}
\value{
  An S4 object of class \code{blavFitIndices} consisting of 2 slots:
  \item{\code{@details}}{A \code{list} containing the choices made by the user
    (or defaults; e.g., which values of \code{pD} and \code{rescale} were set),
    as well as the posterior distribution of the \eqn{\chi^2} (deviance)
    statistic (rescaled, if \code{rescale = "devM"} or \code{"PPMC"}).}
  \item{\code{@indices}}{A \code{list} containing the posterior distribution
    of each requested \code{fit.measure}.}
  
  The \code{summary()} method returns a \code{data.frame} containing one row
  for each requested \code{fit.measure}, and columns containing the specified
  measure(s) of \code{central.tendency}, the posterior \emph{SD},
  and (if requested) the HPD credible-interval limits.
}
\author{
Mauricio Garnier-Villareal (Marquette University;
\email{mauricio.garniervillarreal@marquette.edu})

Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\references{
  \code{rescale = "PPMC"} based on:
  
  Hoofs, H., van de Schoot, R., Jansen, N. W., & Kant, I. (2017).
  Evaluating model fit in Bayesian confirmatory factor analysis with large
  samples: Simulation study introducing the BRMSEA. 
  \emph{Educational and Psychological Measurement}. doi:10.1177/0013164417709314
  
  \code{rescale = "devM"} based on:
  
  Garnier-Villarreal, M., & Jorgensen, T. D. (2018).
  \emph{Adapting fit indices for Bayesian SEM: Comparison to maximum likelihood}.
  Unpublished manuscript (see \url{https://osf.io/afkcw/}).
  
  Other references:
  
  Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
  evaluation using leave-one-out cross-validation and WAIC.
  \emph{Statistics and Computing, 27}(5), 1413--1432.
  doi:10.1007/s11222-016-9696-4
  
}
\examples{ \dontrun{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
## fit target model
fit1 <- bcfa(HS.model, data = HolzingerSwineford1939, cp = "fa",
             n.chains = 2, burnin = 1000, sample = 1000)

## fit null model to calculate CFI, TLI, and NFI 
null.model <- c(paste0("x", 1:9, " ~~ x", 1:9), paste0("x", 1:9, " ~ 1"))
fit0 <- bcfa(null.model, data = HolzingerSwineford1939, cp = "fa",
             n.chains = 2, burnin = 1000, sample = 1000)

## calculate posterior distributions of fit indices

## The default method mimics fit indices derived from ML estimation
ML <- blavFitIndices(fit1, baseline.model = fit0)
ML
summary(ML)

## other options:

## - use Hoofs et al.'s (2017) PPMC-based method
## - use the estimated number of parameters from WAIC instead of LOO-IC
PPMC <- blavFitIndices(fit1, baseline.model = fit0,
                       pD = "waic", rescale = "PPMC")
## issues a warning about using rescale="PPMC" with N < 1000 (see Hoofs et al.)

## - specify only the desired measures of central tendency
## - specify a different "confidence" level for the credible intervals
summary(PPMC, central.tendency = c("mean","mode"), prob = .95)



## Access the posterior distributions for further investigation
head(distML <- data.frame(ML@indices))

## For example, diagnostic plots using the bayesplot package:

## distinguish chains
nChains <- blavInspect(fit1, "n.chains")
distML$Chain <- rep(1:nChains, each = nrow(distML) / nChains)

library(bayesplot)
mcmc_pairs(distML, pars = c("BRMSEA","BMc","BGammaHat","BCFI","BTLI"),
           diag_fun = "hist")
## Indices are highly correlated across iterations in both chains

## Compare to PPMC method
distPPMC <- data.frame(PPMC@indices)
distPPMC$Chain <- rep(1:nChains, each = nrow(distPPMC) / nChains)
mcmc_pairs(distPPMC, pars = c("BRMSEA","BMc","BGammaHat","BCFI","BTLI"),
           diag_fun = "dens")
## nonlinear relation between BRMSEA, related to the floor effect of BRMSEA
## that Hoofs et al. found for larger (12-indicator) models

}}

