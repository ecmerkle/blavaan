\name{blavaan-class}
\docType{class}
\alias{blavaan-class}
\alias{coef,blavaan-method}
\alias{show,blavaan-method}
\alias{predict,blavaan-method}
%\alias{fitted,blavaan-method}
%\alias{fitted.values,blavaan-method}    %% NOTE: no white space here!
\alias{summary,blavaan-method}

\title{Class For Representing A (Fitted) Bayesian Latent Variable Model}
\description{The \code{blavaan} class contains the \code{lavaan}
  class, representing a (fitted) Bayesian latent variable
  model. It contains a description of the model as specified by the user,
  a summary of the data, an internal matrix representation, and if the model
  was fitted, the fitting results.}
\section{Objects from the Class}{
Objects can be created via the
  \code{\link{bcfa}}, \code{\link{bsem}}, \code{\link{bgrowth}} or
  \code{\link{blavaan}} functions.
}
\section{Slots}{
  \describe{
    \item{\code{version}:}{The lavaan package version used to create this objects}
    \item{\code{call}:}{The function call as returned by \code{match.call()}.}
    \item{\code{timing}:}{The elapsed time (user+system) for various parts of
      the program as a list, including the total time.}
    \item{\code{Options}:}{Named list of options that were provided by
      the user, or filled-in automatically.}
    \item{\code{ParTable}:}{Named list describing the model parameters. Can be coerced to a data.frame. In the documentation, this is called the `parameter table'.}
    \item{\code{pta}:}{Named list containing parameter table attributes.}
    \item{\code{Data}:}{Object of internal class \code{"Data"}: information
about the data.}
    \item{\code{SampleStats}:}{Object of internal class \code{"SampleStats"}: sample
      statistics}
    \item{\code{Model}:}{Object of internal class \code{"Model"}: the
      internal (matrix) representation of the model}
    \item{\code{Cache}:}{List using objects that we try to compute only once, and reuse many times.}
    \item{\code{Fit}:}{Object of internal class \code{"Fit"}: the
      results of fitting the model. No longer used.}
    \item{\code{boot}:}{List. Unused for Bayesian models.}
    \item{\code{optim}:}{List. Information about the optimization.}
    \item{\code{loglik}:}{List. Information about the loglikelihood of the model (if maximum likelihood was used).}
    \item{\code{implied}:}{List. Model implied statistics.}
    \item{\code{vcov}:}{List. Information about the variance matrix (vcov) of the model parameters.}
    \item{\code{test}:}{List. Different test statistics.}
    \item{\code{h1}:}{List. Information about the unrestricted h1 model (if available).}
    \item{\code{baseline}:}{List. Information about a baseline model (often
the independence model) (if available).}
    \item{\code{external}:}{List. Includes Stan or JAGS objects used for MCMC.}
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "blavaan", type = "free")}: Returns
      the estimates of the parameters in the model as a named numeric vector.
      If \code{type="free"}, only the free parameters are returned.
      If \code{type="user"}, all parameters listed in the parameter table
      are returned, including constrained and fixed parameters.}
    \item{vcov}{\code{signature(object = "lavaan")}: returns the
      covariance matrix of the estimated parameters.}
    \item{show}{\code{signature(object = "blavaan")}: Print a short summary
      of the model fit}
%    \item{plot}{\code{signature(object = "blavaan")}: S4 method for
%      creating plots. Also see \code{?plot.blavaan}.}
    \item{summary}{\code{signature(object = "blavaan", header = TRUE,
     fit.measures = FALSE, estimates = TRUE, ci = TRUE, 
     standardized = FALSE, rsquare = FALSE, std.nox = FALSE,
     psrf = TRUE, neff = FALSE, postmedian = FALSE, postmode = FALSE,
     priors = TRUE, bf = FALSE, nd = 3L)}:
      Print a nice summary of the model estimates.
      If \code{header = TRUE}, the header section (including fit measures) is
      printed.
      If \code{fit.measures = TRUE}, additional fit measures are added to the
      header section.
      If \code{estimates = TRUE}, print the parameter estimates section.
      If \code{ci = TRUE}, add confidence intervals to the parameter estimates
      section.
      If \code{standardized = TRUE},
      the standardized solution is also printed.  Note that \emph{SE}s and
      tests are still based on unstandardized estimates. Use
      \code{\link[lavaan]{standardizedSolution}} to obtain \emph{SE}s and test
      statistics for standardized estimates.
      If \code{rsquare=TRUE}, the R-Square values for the dependent variables
      in the model are printed.
      If \code{std.nox = TRUE}, the \code{std.all} column contains the
      the \code{std.nox} column from the parameterEstimates() output.
      If \code{psrf = TRUE}, potential scale reduction factors (Rhats)
      are printed.
      If \code{neff = TRUE}, effective sample sizes are printed.
      If \code{postmedian} or \code{postmode} are TRUE, posterior
      medians or modes are printed instead of posterior means.
      If \code{priors = TRUE}, parameter prior distributions are
      printed.
      If \code{bf = TRUE}, Savage-Dickey approximations of the Bayes
      factor are printed for certain parameters.
      Nothing is returned (use
      \code{lavInspect} or another extractor function
      to extract information from a fitted model).}
  }
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
\code{\link{bcfa}}, \code{\link{bsem}}, \code{\link{bgrowth}}
}
\examples{
\dontrun{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- bcfa(HS.model, data=HolzingerSwineford1939)

summary(fit, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)
coef(fit)
}
}
