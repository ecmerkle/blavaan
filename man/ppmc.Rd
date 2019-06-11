\name{ppmc}
\alias{ppmc}
\alias{blavPPMC-class}
\alias{show,blavPPMC-method}
\alias{summary,blavPPMC-method}
\alias{plot,blavPPMC-method}
\alias{hist,blavPPMC-method}
\title{Posterior Predictive Model Checks}
\description{
This function allows users to conduct a posterior predictive model check to
assess the global or local fit of a latent variable model using any discrepancy
function that can be applied to a \code{\linkS4class{lavaan}} model.}
\usage{
ppmc(object, thin = 1, fit.measures = c("srmr","chisq"), discFUN = NULL)

\S4method{summary}{blavPPMC}(object, discFUN, dist = c("obs","sim"),
        central.tendency = c("mean","median","mode"),
        hpd = TRUE, prob = .95)

\S3method{plot}{blavPPMC}(x, ..., discFUN, element, central.tendency = "",
     hpd = TRUE, prob = .95, nd = 3)

\S3method{hist}{blavPPMC}(x, ..., discFUN, element, hpd = TRUE, prob = .95,
     printLegend = TRUE, legendArgs = list(x = "topleft"),
     densityArgs = list(), nd = 3)
}
\arguments{
  \item{object,x}{An object of class \code{\linkS4class{blavaan}}.}
  \item{thin}{Optional \code{integer} indicating how much to thin each chain.
    Default is \code{1L}, indicating not to thin the chains in \code{object}.}
  \item{fit.measures}{\code{character} vector indicating the names of global
    discrepancy measures returned by \code{\link[lavaan]{fitMeasures}}. Ignored
    unless \code{discFUN} is \code{NULL}, but users may include
    \code{fitMeasures} in the \code{list} of discrepancy functions in
    \code{discFUN}. If the first measure is either \code{"logl"} or
    \code{"chisq"}, only the \eqn{\chi^2} fit statistic's posterior (predictive)
    distributions will be returned.}
  \item{discFUN}{\code{function}, or a \code{list} of functions, that can be
    called on an object of class \code{\linkS4class{lavaan}}. Each function
    must return an object whose \code{\link{mode}} is \code{numeric}, but may
    be a \code{vector}, \code{matrix}, or multidimensional \code{array}. In
    the \code{summary} and \code{plot} methods, \code{discFUN} is a
    \code{character} indicating which discrepancy function to summarize.}
  \item{element}{\code{numeric} or \code{character} indicating the index (in
    each \code{dim}ension of the \code{discFUN} output, if multiple) to plot.}
  \item{dist}{\code{character} indicating whether to summarize the distribution
    of \code{discFUN} on either the \code{obs}erved or \code{sim}ulated data.}
  \item{central.tendency}{\code{character} indicating which statistics should
    be used to characterize the location of the posterior (predictive)
    distribution. By default, all 3 statistics are returned for the
    \code{summary} method, but none for the \code{plot} method. The posterior
    mean is labeled \code{EAP} for \emph{expected a posteriori} estimate, and
    the mode is labeled \code{MAP} for \emph{modal a posteriori} estimate.}
  \item{hpd}{\code{logical} indicating whether to calculate the highest
    posterior density (HPD) credible interval for \code{discFUN}.}
  \item{prob}{The "confidence" level of the credible interval(s).}
  \item{nd}{The number of digits to print in the scatter\code{plot}.}
  \item{\dots}{Additional \code{graphical \link[graphics]{par}ameters} to be
    passed to \code{\link[graphics]{plot.default}}.}
  \item{printLegend}{\code{logical}. If \code{TRUE} (default), a legend will
    be printed with the histogram}
  \item{legendArgs}{\code{list} of arguments passed to the
    \code{\link[graphics]{legend}} function.  The default argument is a list
    placing the legend at the top-left of the figure.}
  \item{densityArgs}{\code{list} of arguments passed to the
    \code{\link[stats]{density}} function, used to obtain densities for the
    \code{hist} method.}
}
\value{
  An S4 object of class \code{blavPPMC} consisting of 5 \code{list} slots:
  \item{\code{@discFUN}}{The user-supplied \code{discFUN}, or the \code{call} to
    \code{\link[lavaan]{fitMeasures}} that returns \code{fit.measures}.}
  \item{\code{@dims}}{The dimensions of the object returned by each
    \code{discFUN}.}
  \item{\code{@PPP}}{The posterior predictive \emph{p} value for each
    \code{discFUN} element.}
  \item{\code{@obsDist}}{The posterior distribution of realize values of
    \code{discFUN} applied to observed data.}
  \item{\code{@simDist}}{The posterior predictive distribution of values of
    \code{discFUN} applied to data simulated from the posterior samples.}

  The \code{summary()} method returns a \code{numeric} vector if \code{discFUN}
    returns a scalar, a \code{data.frame} with one discrepancy function per row
    if \code{discFUN} returns a \code{numeric} vector, and a \code{list} with
    one summary statistic per element if \code{discFUN} returns a \code{matrix}
    or multidimensional \code{array}.

  The \code{plot} method invisible returns \code{NULL}, printing a plot to the
    current device.
}
\author{
Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\references{
  Levy, R. (2011). Bayesian data--model fit assessment for structural equation
  modeling. \emph{Structural Equation Modeling, 18}(4), 663--685.
  doi:10.1080/10705511.2011.607723
}
\examples{ \dontrun{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
## fit single-group model
fit <- bcfa(HS.model, data = HolzingerSwineford1939, cp = "fa",
            target = "jags", bcontrol = list(method = "rjparallel"),
            n.chains = 2, burnin = 1000, sample = 500)
## fit multigroup model
fitg <- bcfa(HS.model, data = HolzingerSwineford1939, cp = "fa",
             target = "jags", bcontrol = list(method = "rjparallel"),
             n.chains = 2, burnin = 1000, sample = 500, group = "school")


## use fit.measures as a shortcut for global fitMeasures only

AFIs <- ppmc(object, thin = 10, fit.measures = c("srmr","chisq","rmsea","cfi"))
summary(AFIs)          # summarize all of a discrepancy function's output
plot(AFIs, y = "srmr") # only plot one discrepancy function at a time


## define a list of custom discrepancy functions
## - model-implied correlations
## - (global) fit measures
## - (local) standardized residuals

discFUN <- list(impCor = function(fit) lavInspect(fit, "cor.lv"),
                fitIndices = function(fit) {
                  fitMeasures(fit, fit.measures = c("cfi","rmsea","srmr","chisq"))
                },
                cor.resid = function(fit) lavResiduals(fit)$cov,
                mean.resid = function(fit) lavResiduals(fit)$mean)
out1g <- ppmc(object, discFUN = discFUN)
summary(out1g, discFUN = "impCor")
plot(out1g, y = c("visual","speed"), discFUN = "impCor")


## define a list of custom discrepancy functions for multiple groups
## (must return each group's numeric output is a different function)

disc2g <- list(impCor1 = function(fit) lavInspect(fit, "cor.lv")[[1]],
               impCor2 = function(fit) lavInspect(fit, "cor.lv")[[2]],
               fitIndices = function(fit) {
                 fitMeasures(fit, fit.measures = c("cfi","rmsea","mfi","srmr","chisq"))
               },
               cor.resid1 = function(fit) lavResiduals(fit)[[1]]$cov,
               cor.resid2 = function(fit) lavResiduals(fit)[[2]]$cov)
out2g <- ppmc(fitg, discFUN = disc2g, thin = 2)
summary(out2g, discFUN = "cor.resid1")
plot(out2g, y = c("x3","x9"), discFUN = "cor.resid1") # group 1
plot(out2g, y = c(3, 9), discFUN = "cor.resid2")      # group 2

}}

