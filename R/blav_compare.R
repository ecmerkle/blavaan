blavCompare <- function(object1, object2, ...) {
  ## loo compare code from Mauricio Garnier-Villarreal + old BF() code

  ## Bayes factor approximation based on marginal log-likelihoods
  bf <- object1@test[[1]]$stat - object2@test[[1]]$stat
  res <- c(bf, object1@test[[1]]$stat, object2@test[[1]]$stat)
  names(res) <- c("bf", "mll1", "mll2")
  
  ## FIXME? We already get case_lls in blav_fit_measures and should really
  ## only do it once. But, if we store it in the blavaan object, the size
  ## of that object can get much larger.
  lavopt1 <- object1@Options
  lavopt1$estimator <- "ML"
  ll1 <- case_lls(object1@external$runjags, object1@Model,
                  object1@ParTable, object1@SampleStats,
                  lavopt1, object1@Cache,
                  object1@Data)

  lavopt2 <- object2@Options
  lavopt2$estimator <- "ML"
  ll2 <- case_lls(object2@external$runjags, object2@Model,
                  object2@ParTable, object2@SampleStats,
                  lavopt2, object2@Cache,
                  object2@Data)

  loo1 <- loo(ll1); loo2 <- loo(ll2)
  waic1 <- waic(ll1); waic2 <- waic(ll2)

  diff_loo <- compare(loo1, loo2)
  diff_waic <- compare(waic1, waic2)

  cat("\nWAIC difference & SE (positive values favor object2): \n",
      sprintf("%8.3f", diff_waic[1]),
      sprintf("%8.3f", diff_waic[2]), "\n\n")

  cat("LOO difference & SE (positive values favor object2): \n",
      sprintf("%8.3f", diff_loo[1]),
      sprintf("%8.3f", diff_loo[2]), "\n\n")

  cat("Laplace approximation to the log-Bayes factor\n(experimental; positive values favor object1):",
      sprintf("%8.3f", bf), "\n\n")
  
  res <- list(bf = res, loo = rbind(loo1[1:6], loo2[1:6]),
              diff_loo = diff_loo,
              waic = rbind(waic1[1:6], waic2[1:6]),
              diff_waic = diff_waic)

  invisible(res)
}

BF <- function(object1, object2, ...) {
  cat("BF() is deprecated. Use blavCompare() instead.\n\n")

  blavCompare(object1, object2)
}
