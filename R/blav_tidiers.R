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
#' @param summary.type The posterior summary statistic to use as the point
#'   estimate. One of \code{"mean"} (posterior mean, default), \code{"median"}
#'   (posterior median), or \code{"mode"} (posterior mode). Fixed and
#'   constrained parameters always use the posterior mean regardless.
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
#' @param ... Additional arguments passed to \code{summary()}.
#'
#' @return A \code{data.frame} with columns:
#'   \item{term}{Parameter name (lhs, op, rhs combined)}
#'   \item{op}{Operator from the model syntax}
#'   \item{estimate}{Posterior summary statistic determined by \code{summary.type}}
#'   \item{std.error}{Posterior standard deviation}
#'   \item{conf.low}{Lower bound of 95\% credible interval (if \code{conf.int = TRUE})}
#'   \item{conf.high}{Upper bound of 95\% credible interval (if \code{conf.int = TRUE})}
#'   \item{std.all}{Standardized estimate (if \code{standardized = TRUE})}
#'   \item{rhat}{Rhat convergence diagnostic (if \code{rhat = TRUE})}
#'   \item{ess}{Effective sample size (if \code{ess = TRUE})}
#'   \item{prior}{Prior distribution specification (if \code{priors = TRUE})}
#'   \item{group}{Group number (for multigroup models)}
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
#' tidy(fit, summary.type = "median")
#' tidy(fit, standardized = TRUE)
#' }
#'
#' @export
tidy.blavaan <- function(x, summary.type = c('mean', 'median', 'mode'),
                          conf.int = TRUE, standardized = FALSE,
                         rhat = TRUE, ess = TRUE, priors = TRUE, ...) {

  summary.type <- match.arg(summary.type)

  # Get parameter estimates
  PE <- summary(
    x,
    header = FALSE,
    print = FALSE,
    standardized = standardized,
    neff = ess,
    postmedian = summary.type == "median",
    postmode = summary.type == "mode",
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
  if ("group" %in% names(parameterEstimates(x))) {
    result$group <- PE$group
    result <- result[, c('term', 'op', 'group', 'estimate', 'std.error')]
  }

  # Substitute estimate based on summary.type
  if (summary.type == "median") {
    PE$Post.Med[is.na(PE$Post.Med)] <- PE$est[is.na(PE$Post.Med)]
    result$estimate <- PE$Post.Med
  } else if (summary.type == "mode") {
    PE$Post.Mode[is.na(PE$Post.Mode)] <- PE$est[is.na(PE$Post.Mode)]
    result$estimate <- PE$Post.Mode
  }

  # If conf.int = TRUE, add confidence interval
  if (isTRUE(conf.int)) {
    result$conf.low <- as.numeric(PE$pi.lower)
    result$conf.high <- as.numeric(PE$pi.upper)
  }

  # Add standardized estimates if requested
  if (isTRUE(standardized)) {
    result$std.all <- PE$std.all
  }

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
#' @param fit.indices Character vector of fit indices to compute from
#'   \code{blavFitIndices}. Use \code{"none"} to skip Bayesian fit indices
#'   (faster). Use \code{"default"} for BRMSEA, BGammaHat, and BMc.
#'   Default is \code{"none"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A single-row \code{tibble} (if available) or \code{data.frame}
#'   with columns:
#'   \item{npar}{Number of estimated parameters}
#'   \item{nobs}{Total number of observations}
#'   \item{ngroups}{Number of groups}
#'   \item{estimator}{Estimation method used}
#'   \item{ppp}{Posterior predictive p-value}
#'   \item{looic}{Leave-one-out information criterion}
#'   \item{p_loo}{Effective number of parameters (LOO)}
#'   \item{waic}{Widely applicable information criterion}
#'   \item{p_waic}{Effective number of parameters (WAIC)}
#'   \item{dic}{Deviance information criterion}
#'   \item{p_dic}{Effective number of parameters (DIC)}
#'   \item{bic}{Bayesian information criterion}
#'   \item{margloglik}{Marginal log-likelihood}
#'   \item{converged}{Logical indicating convergence (all Rhat < 1.05)}
#'   \item{nchains}{Number of MCMC chains}
#'   Additional columns from \code{blavFitIndices} are included if requested.
#'
#' @examples
#' \dontrun{
#' library(blavaan)
#'
#' HS.model <- 'visual =~ x1 + x2 + x3'
#' fit <- bcfa(HS.model, data = HolzingerSwineford1939, seed = 123)
#' glance(fit)
#' glance(fit, fit.indices = "default")  # includes BRMSEA etc.
#' }
#'
#' @export
glance.blavaan <- function(x, fit.indices = "none", ...) {

  # Get basic model information
  ngroups <- blavInspect(x, "ngroups")
  nobs <- blavInspect(x, "ntotal")
  nchains <- blavInspect(x, "n.chains")

  # Get fit measures
  bopts <- blavInspect(x, "options")
  test_available <- bopts$test != "none"

  # Initialize result
  result <- data.frame(
    npar = as.integer(lavInspect(x, "npar")),
    nobs = as.integer(nobs),
    ngroups = as.integer(ngroups),
    estimator = "Bayes",
    stringsAsFactors = FALSE
  )

  # Get available fit measures
  if (test_available) {
    fm <- tryCatch(
      blav_fit_measures(x, fit.measures = "all"),
      error = function(e) NULL
    )

    if (!is.null(fm)) {
      # Posterior predictive p-value
      if ("ppp" %in% names(fm)) {
        result$ppp <- unname(fm["ppp"])
      }

      # Information criteria
      if ("looic" %in% names(fm)) {
        result$looic <- unname(fm["looic"])
        result$p_loo <- unname(fm["p_loo"])
      }
      if ("waic" %in% names(fm)) {
        result$waic <- unname(fm["waic"])
        result$p_waic <- unname(fm["p_waic"])
      }
      if ("dic" %in% names(fm)) {
        result$dic <- unname(fm["dic"])
        result$p_dic <- unname(fm["p_dic"])
      }
      if ("bic" %in% names(fm)) {
        result$bic <- unname(fm["bic"])
      }
      if ("margloglik" %in% names(fm)) {
        result$margloglik <- unname(fm["margloglik"])
      }
    }
  }

  # Check convergence (Rhat < 1.05 for all parameters)
  rhat_vals <- tryCatch(
    blavInspect(x, "rhat"),
    error = function(e) NULL
  )
  if (!is.null(rhat_vals)) {
    result$converged <- all(rhat_vals < 1.05, na.rm = TRUE)
  }

  result$nchains <- as.integer(nchains)

  # Add Bayesian fit indices if requested
  if (!identical(fit.indices, "none")) {
    bfi <- tryCatch({
      if (identical(fit.indices, "default")) {
        blavFitIndices(x, fit.measures = c("BRMSEA", "BGammaHat", "BMc"))
      } else {
        blavFitIndices(x, fit.measures = fit.indices)
      }
    }, error = function(e) NULL)

    if (!is.null(bfi)) {
      bfi_summary <- summary(bfi, central.tendency = "mean", hpd = FALSE)
      # Add EAP (posterior mean) for each fit index
      for (i in seq_len(nrow(bfi_summary))) {
        idx_name <- rownames(bfi_summary)[i]
        result[[idx_name]] <- bfi_summary[i, "EAP"]
      }
    }
  }

  return(result)
}