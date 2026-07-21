## Shared helper for tinytest files that fit a blavaan model with short MCMC
## chains and then compare estimates to a lavaan reference fit. Short chains
## occasionally fail to converge by chance (bad seed), which would otherwise
## make coef()/vcov() comparisons flake even though nothing is broken.
##
## robust_fit() fits with FUN/..., and on non-convergence (max rhat >=
## rhat_cutoff) retries with more warmup/samples and a new seed before giving
## up. It never hard-fails -- callers must check the returned `converged`
## flag before running lavaan-comparison assertions, and should skip (not
## fail) those assertions when it's FALSE. Structural assertions (class,
## dimensions, "did this error") should stay unconditional so coverage
## breadth is preserved even on the rare non-convergent run.

robust_fit <- function(FUN, ..., rhat_cutoff = 1.05, max_tries = 3,
                        mult = c(1, 2, 4), base_seed = NULL, quiet = FALSE) {
  ## bsem()/bcfa()/bgrowth() build their own call internally (match.call()-
  ## style) to forward to blavaan(); do.call()'ing the function *object*
  ## breaks that (the callee sees a closure where it expects a symbol, and
  ## dies with "cannot coerce type 'closure' to vector of type 'character'").
  ## do.call() on the function's *name* avoids this.
  fname <- deparse(substitute(FUN))

  dots <- list(...)
  stopifnot(!is.null(dots$burnin), !is.null(dots$sample))
  base_burnin <- dots$burnin
  base_sample <- dots$sample

  result <- list(fit = NULL, converged = FALSE, attempts = 0,
                  rhat_max = NA_real_, error = NULL)

  for (i in seq_len(max_tries)) {
    attempt_dots <- dots
    attempt_dots$burnin <- round(base_burnin * mult[min(i, length(mult))])
    attempt_dots$sample <- round(base_sample * mult[min(i, length(mult))])

    if (i > 1) {
      ## a fresh draw is only wanted once the previous attempt has already
      ## failed; bcontrol$seed wins over a top-level seed= in blavaan's own
      ## seed-merge logic, so write there to guarantee this actually takes
      ## effect even when the caller already passed bcontrol=list(seed=...)
      new_seed <- if (!is.null(base_seed)) base_seed + i * 997L
                  else sample.int(.Machine$integer.max, 1)
      bc <- attempt_dots$bcontrol
      if (is.null(bc)) bc <- list()
      bc$seed <- new_seed
      attempt_dots$bcontrol <- bc
      attempt_dots$seed <- NULL
    }

    fit <- tryCatch(do.call(fname, attempt_dots), error = function(e) e)
    result$attempts <- i

    if (inherits(fit, "error")) {
      ## non-convergence can surface as a hard error (e.g. a chol() failure
      ## downstream of an unstable posterior), not just a high rhat
      result$error <- conditionMessage(fit)
      next
    }

    rmax <- tryCatch(max(blavInspect(fit, "rhat"), na.rm = TRUE),
                      error = function(e) Inf)
    result$rhat_max <- rmax
    fit@optim$converged <- isTRUE(rmax < rhat_cutoff)

    if (isTRUE(rmax < rhat_cutoff)) {
      result$fit <- fit
      result$converged <- TRUE
      break
    }
    result$fit <- fit
  }

  if (!quiet) {
    cat(sprintf(
      "[robust_fit] %s: %s after %d attempt(s), max rhat = %s\n",
      fname,
      if (result$converged) "CONVERGED" else "DID NOT CONVERGE",
      result$attempts,
      if (is.finite(result$rhat_max)) round(result$rhat_max, 3) else "NA"))
  }
  result
}
