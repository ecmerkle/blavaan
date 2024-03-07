\name{ppmc}
\alias{ppmc}
\alias{blavPPMC-class}
\alias{show,blavPPMC-method}
\alias{summary,blavPPMC-method}
\alias{summary.ppmc}
\alias{plot.blavPPMC}
\alias{hist.blavPPMC}
\alias{pairs.blavPPMC}
\title{Posterior Predictive Model Checks}
\description{
This function allows users to conduct a posterior predictive model check to
assess the global or local fit of a latent variable model using any discrepancy
function that can be applied to a \code{\linkS4class{lavaan}} model.}
\usage{
ppmc(object, thin = 1, fit.measures = c("srmr","chisq"), discFUN = NULL,
     conditional = FALSE)

\S4method{summary}{blavPPMC}(object, ...)

\method{summary}{ppmc}(object, discFUN, dist = c("obs","sim"),
        central.tendency = c("mean","median","mode"),
        hpd = TRUE, prob = .95, to.data.frame = FALSE, diag = TRUE,
        sort.by = NULL, decreasing = FALSE)

\S3method{plot}{blavPPMC}(x, ..., discFUN, element, central.tendency = "",
     hpd = TRUE, prob = .95, nd = 3)

\S3method{hist}{blavPPMC}(x, ..., discFUN, element, hpd = TRUE, prob = .95,
     printLegend = TRUE, legendArgs = list(x = "topleft"),
     densityArgs = list(), nd = 3)

\S3method{pairs}{blavPPMC}(x, discFUN, horInd = 1:DIM, verInd = 1:DIM,
      printLegend = FALSE, ...)
}
\arguments{
  \item{object,x}{An object of class \code{\linkS4class{blavaan}}.}
  \item{thin}{Optional \code{integer} indicating how much to thin each chain.
    Default is \code{1L}, indicating not to thin the chains in \code{object}.}
  \item{fit.measures}{\code{character} vector indicating the names of global
    discrepancy measures returned by \code{\link[lavaan]{fitMeasures}}. Ignored
    unless \code{discFUN} is \code{NULL}, but users may include
    \code{fitMeasures} in the \code{list} of discrepancy functions in
    \code{discFUN}. For ordinal models, the \code{"logl"} or \code{"chisq"}
    computations are done via lavaan.}
  \item{discFUN}{\code{function}, or a \code{list} of functions, that can be
    called on an object of class \code{\linkS4class{lavaan}}. Each function
    must return an object whose \code{\link{mode}} is \code{numeric}, but may
    be a \code{vector}, \code{matrix}, or multidimensional \code{array}. In
    the \code{summary} and \code{plot} methods, \code{discFUN} is a
    \code{character} indicating which discrepancy function to
    summarize.}
  \item{conditional}{\code{logical} indicating whether or not, during artificial data
    generation, we should condition on the estimated latent
  variables. Requires the model to be estimated with \code{save.lvs = TRUE}.}
  \item{element}{\code{numeric} or \code{character} indicating the index (in
    each \code{dim}ension of the \code{discFUN} output, if multiple) to plot.}
  \item{horInd,verInd}{Similar to \code{element}, but a \code{numeric} or
    \code{character} vector indicating the indices of a \code{matrix} to plot
    in a scatterplot matrix. If \code{horInd==verInd}, histograms will be
    plotted in the upper triangle.}
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
  \item{to.data.frame}{\code{logical} indicating whether the \code{summary} of a
    symmetric 2-dimensional \code{matrix} returned by \code{discFUN} should have
    its unique elements stored in rows of a \code{data.frame} that can be sorted
    for convenience of identifying large discrepancies. If \code{discFUN}
    returns an asymmetric 2-dimensional \code{matrix}, the list of matrices
    returned by the \code{summary} can also be converted to a \code{data.frame}.}
  \item{diag}{Passed to \code{\link[base]{lower.tri}} if \code{to.data.frame=TRUE}.}
  \item{sort.by}{\code{character}. If \code{summary} returns a \code{data.frame},
    it can be sorted by this column name using \code{\link[base]{order}}. Note
    that if \code{discFUN} returns an asymmetric 2-dimensional \code{matrix},
    each \code{data.frame} in the returned \code{list} will be sorted
    independently, so the rows are unlikely to be consistent across
    summary statistics.}
  \item{decreasing}{Passed to \code{\link[base]{order}} if
    \code{!is.null(sort.by)}.}
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

  The \code{plot} and \code{pairs} methods invisibly return \code{NULL},
  printing a plot (or scatterplot matrix) to the current device.

  The \code{hist} method invisibly returns a \code{list} or arguments that can
  be passed to the function for which the \code{list} element is named.  Users
  can edit the arguments in the list to customize their histograms.
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
data(HolzingerSwineford1939, package = "lavaan")

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
## fit single-group model
fit <- bcfa(HS.model, data = HolzingerSwineford1939, 
            n.chains = 2, burnin = 1000, sample = 500)
## fit multigroup model
fitg <- bcfa(HS.model, data = HolzingerSwineford1939,
             n.chains = 2, burnin = 1000, sample = 500, group = "school")


## Use fit.measures as a shortcut for global fitMeasures only
## - Note that indices calculated from the "df" are only appropriate under
##   noninformative priors, such that pD approximates the number of estimated
##   parameters counted under ML estimation; incremental fit indices
##   introduce further complications)

AFIs <- ppmc(fit, thin = 10, fit.measures = c("srmr","chisq","rmsea","cfi"))
summary(AFIs)                 # summarize the whole vector in a data.frame
hist(AFIs, element = "rmsea") # only plot one discrepancy function at a time
plot(AFIs, element = "srmr")


## define a list of custom discrepancy functions
## - (global) fit measures
## - (local) standardized residuals

discFUN <- list(global = function(fit) {
                  fitMeasures(fit, fit.measures = c("cfi","rmsea","srmr","chisq"))
                },
                std.cov.resid = function(fit) lavResiduals(fit, zstat = FALSE,
                                                           summary = FALSE)$cov,
                std.mean.resid = function(fit) lavResiduals(fit, zstat = FALSE,
                                                            summary = FALSE)$mean)
out1g <- ppmc(fit, discFUN = discFUN)

## summarize first discrepancy by default (fit indices)
summary(out1g)
## some model-implied correlations look systematically over/underestimated
summary(out1g, discFUN = "std.cov.resid", central.tendency = "EAP")
hist(out1g, discFUN = "std.cov.resid", element = c(1, 7))
plot(out1g, discFUN = "std.cov.resid", element = c("x1","x7"))
## For ease of investigation, optionally export summary as a data.frame,
## sorted by size of average residual
summary(out1g, discFUN = "std.cov.resid", central.tendency = "EAP",
        to.data.frame = TRUE, sort.by = "EAP")
## or sorted by size of PPP
summary(out1g, discFUN = "std.cov.resid", central.tendency = "EAP",
        to.data.frame = TRUE, sort.by = "PPP_sim_LessThan_obs")

## define a list of custom discrepancy functions for multiple groups
## (return each group's numeric output using a different function)

disc2g <- list(global = function(fit) {
                 fitMeasures(fit, fit.measures = c("cfi","rmsea","mfi","srmr","chisq"))
               },
               cor.resid1 = function(fit) lavResiduals(fit, zstat = FALSE,
                                                       type = "cor.bollen",
                                                       summary = FALSE)[[1]]$cov,
               cor.resid2 = function(fit) lavResiduals(fit, zstat = FALSE,
                                                       type = "cor.bollen",
                                                       summary = FALSE)[[2]]$cov)
out2g <- ppmc(fitg, discFUN = disc2g, thin = 2)
## some residuals look like a bigger problem in one group than another
pairs(out2g, discFUN = "cor.resid1", horInd = 1:3, verInd = 7:9) # group 1
pairs(out2g, discFUN = "cor.resid2", horInd = 1:3, verInd = 7:9) # group 2

## print all to file: must be a LARGE picture. First group 1 ...
png("cor.resid1.png", width = 1600, height = 1200)
pairs(out2g, discFUN = "cor.resid1")
dev.off()
## ... then group 2
png("cor.resid2.png", width = 1600, height = 1200)
pairs(out2g, discFUN = "cor.resid2")
dev.off()
}}

