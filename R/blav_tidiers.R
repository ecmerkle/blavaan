# Tidy and glance methods for blavaan objects
# Following broom ecosystem conventions for Bayesian models

# Import generics from the generics package to ensure compatibility with broom
#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' Tidy a blavaan object
#' 
#' This function repackages the output from \code{summary()} into a friendly \code{data.frame}.
#'
#' @param x A \code{blavaan} object.
#' @param estimate.method The posterior summary statistic to compute the point
#'   estimates. One of \code{"mean"} (posterior mean, default), \code{"median"}
#'   (posterior median), or \code{"mode"} (posterior mode). Only free parameters
#'   are summarized by the requested statistic; fixed parameters retain their
#'   fixed values. Defined/constrained parameters (\code{:=}) keep the point
#'   estimate from the standard \code{summary()} output regardless of
#'   \code{estimate.method}.
#' @param conf.int Logical indicating whether to include credible intervals.
#'   Default is \code{TRUE}.
#' @param standardized Logical indicating whether to include standardized
#'   posterior mean estimates. Default is \code{FALSE}.
#' @param rhat Logical indicating whether to include the Rhat convergence
#'   diagnostic. Default is \code{TRUE}.
#' @param ess Logical indicating whether to include the effective sample size.
#'   Default is \code{TRUE}.
#' @param priors Logical indicating whether to include the prior distribution
#'   specification. Default is \code{TRUE}.
#'
#' @return A \code{data.frame} with columns:
#'   \item{term}{Parameter name (lhs, op, rhs combined)}
#'   \item{op}{Operator from the model syntax}
#'   \item{level} The level of a parameter estimate
#'   \item{group}{Group number (for multigroup models)}
#'   \item{estimate}{Posterior summary statistic determined by \code{estimate.method}}
#'   \item{std.error}{Posterior standard deviation}
#'   \item{conf.low}{Lower bound of 95\% credible interval (if \code{conf.int = TRUE})}
#'   \item{conf.high}{Upper bound of 95\% credible interval (if \code{conf.int = TRUE})}
#'   \item{std.lv}{Standardized estimates based on the variances of the
#'     (continuous) latent variables only}
#'   \item{std.all}{Standardized estimates based on both the variances
#'     of both (continuous) observed and latent variables.}
#'   \item{rhat}{Rhat convergence diagnostic (if \code{rhat = TRUE})}
#'   \item{ess}{Effective sample size (if \code{ess = TRUE})}
#'   \item{prior}{Prior distribution specification (if \code{priors = TRUE})}
#'
#' @importFrom generics tidy
#' @examples
#' \dontrun{
#' data(HolzingerSwineford1939, package = "lavaan")
#'
#' HS.model <- 'visual =~ x1 + x2 + x3'
#' fit <- bcfa(HS.model, data = HolzingerSwineford1939, seed = 123,
#'             n.chains = 1, sample = 300)
#' tidy(fit)
#' tidy(fit, estimate.method = "median")
#' 
#' data(Demo.twolevel, package = "lavaan")
#' model <- "
#'     level: within
#'         fw =~ y1 + y2 + y3
#'         fw ~ x1 + x2 + x3
#'     level: between
#'         fb =~ y1 + y2 + y3
#'         fb ~ w1 + w2
#' "
#' bfit <- bsem(
#'   model = model,
#'   data = Demo.twolevel,
#'   cluster = "cluster",
#'   seed = 123,
#'   n.chains = 1,
#'   sample = 300,
#'   target = "stan"
#' )
#' tidy(fit)
#' }
#'
#' @export
tidy.blavaan <- function(x, estimate.method = c('mean', 'median', 'mode'),
                          conf.int = TRUE, rhat = TRUE, ess = TRUE, priors = TRUE) {

  estimate.method <- match.arg(estimate.method)

  # Get parameter estimates
  PE <- summary(
    x,
    header = FALSE,
    print = FALSE,
    standardized = TRUE,
    neff = ess,
    postmedian = estimate.method == "median",
    postmode = estimate.method == "mode",
    ci = conf.int,
    priors = priors
  )

  # Build the base data frame
  result <- data.frame(
    term = paste0(PE$lhs, PE$op, PE$rhs),
    op = PE$op,
    estimate = PE$est,
    std.error = PE$Post.SD,
    stringsAsFactors = FALSE
  )

  # Add group information for multigroup models
  if (length(unique(PE$group)) > 1) {
    result$group <- as.integer(PE$group)
    result <- result[, c('term', 'op', 'group', 'estimate', 'std.error')]
  }

  # Add block and level information for multilevel models 
  if (length(unique(PE$level)) > 1) {
    result$level <- PE$level
    if ("group" %in% names(result)) {
      result <- result[, c('term', 'op', 'level', 'group', 'estimate', 'std.error')]
    } else 
    {
      result <- result[, c('term', 'op', 'level', 'estimate', 'std.error')]
    }
    
  }

  # Substitute estimate based on estimate.method
  if (estimate.method == "median") {
    result$estimate <- PE$Post.Med
  } else if (estimate.method == "mode") {
    result$estimate <- PE$Post.Mode
  }

  # If conf.int = TRUE, add confidence interval
  if (isTRUE(conf.int)) {
    result$conf.low <- as.numeric(PE$pi.lower)
    result$conf.high <- as.numeric(PE$pi.upper)
  }

  # Add standardized estimates
  result$std.lv <- PE$std.lv
  result$std.all <- PE$std.all

  # Add rhat from PE, if rhat = TRUE
  if (isTRUE(rhat)) {
    result$rhat <- as.numeric(PE$Rhat)
  }

  # Add effective sample size if requested
  if (isTRUE(ess)) {
    result$ess <- PE$neff
  }

  # Add prior information
  if (isTRUE(priors)) {
    result$prior <- PE$prior
  }

  return(result)
}


#' Glance at a blavaan object
#'
#' @param x A \code{blavaan} object
#' @param fit.indices Character vector of fit indices to compute from \code{fitMeasures}
#'    and \code{blavFitIndices}. Use \code{"TRUE"} to also compute \code{blavFitIndices} measures. 
#'    Default is \code{"FALSE"}.
#'   
#' @param ... Additional arguments (currently ignored).
#'
#' @return A single-row \code{data.frame}
#'   with columns:
#'   \item{npar}{Number of estimated parameters}
#'   \item{ngroups}{Number of groups}
#'   \item{ppp}{Posterior predictive p-value}
#'   \item{bic}{Bayesian information criterion}
#'   \item{dic}{Deviance information criterion}
#'   \item{waic}{Widely applicable information criterion}
#'   \item{looic}{Leave-one-out information criterion}
#'   \item{margloglik}{Marginal log-likelihood}
#'   Additional EAP of fit measures from \code{blavFitIndices} are included if requested.
#'
#' @examples
#' \dontrun{
#' library(blavaan)
#'
#' HS.model <- 'visual =~ x1 + x2 + x3'
#' fit <- bcfa(HS.model, data = HolzingerSwineford1939, seed = 123)
#' glance(fit)
#' glance(fit, fit.indices = TRUE)  # includes BRMSEA
#' }
#'
#' @export
glance.blavaan <- function(x, fit.indices = FALSE, ...) {

  if (isTRUE(fit.indices)) {
    bayesFit <- blavFitIndices(x, fit.measures = "all") |> 
      summary()
  }
  
  bayesFitVec <- bayesFit$EAP
  names(bayesFitVec) <- rownames(bayesFit)
  
  result <- c(fitMeasures(x), bayesFitVec)

  return(result)
}