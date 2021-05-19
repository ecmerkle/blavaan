blavCompare <- function(object1, object2, ...) {
  ## loo compare code from Mauricio Garnier-Villarreal + old BF() code
  ## possible TODO: compare using conditional likelihoods, in addition to marginal
  lavopt1 <- blavInspect(object1, "Options")
  lavopt2 <- blavInspect(object2, "Options")
  if((lavopt1$test == "none" & lavopt1$target != "stan") |
     (lavopt2$test == "none" & lavopt2$target != "stan")){
    stop("blavaan ERROR: Models cannot be compared when test='none'")
  }
  targ1 <- lavopt1$target; targ2 <- lavopt2$target
    
  ## Bayes factor approximation based on marginal log-likelihoods
  bf <- object1@test[[1]]$stat - object2@test[[1]]$stat
  res <- c(bf, object1@test[[1]]$stat, object2@test[[1]]$stat)
  names(res) <- c("bf", "mll1", "mll2")
  
  ## FIXME? We already get case_lls in blav_fit_measures and should really
  ## only do it once. But, if we store it in the blavaan object, the size
  ## of that object can get much larger.
  if(targ1 == "stan"){
    ll1 <- loo::extract_log_lik(object1@external$mcmcout)
  } else {
    lavopt1$estimator <- "ML"
    ll1 <- case_lls(object1@external$mcmcout, object1@Model,
                    object1@ParTable, object1@SampleStats,
                    lavopt1, object1@Cache,
                    object1@Data, make_mcmc(object1@external$mcmcout))
  }
  nchain1 <- blavInspect(object1, "n.chains")
  niter1 <- nrow(ll1)/nchain1
  cid1 <- rep(1:nchain1, each=niter1)
  ref1 <- relative_eff(exp(ll1), chain_id = cid1)

  if(targ2 == "stan"){
    ll2 <- loo::extract_log_lik(object2@external$mcmcout)
  } else {
    lavopt2$estimator <- "ML"
    ll2 <- case_lls(object2@external$mcmcout, object2@Model,
                    object2@ParTable, object2@SampleStats,
                    lavopt2, object2@Cache,
                    object2@Data, make_mcmc(object2@external$mcmcout))
  }
  nchain2 <- blavInspect(object1, "n.chains")
  niter2 <- nrow(ll2)/nchain2
  cid2 <- rep(1:nchain2, each=niter2)
  ref2 <- relative_eff(exp(ll2), chain_id = cid2)

  loo1 <- loo(ll1, r_eff=ref1)
  loo2 <- loo(ll2, r_eff=ref2)
  waic1 <- waic(ll1); waic2 <- waic(ll2)

  diff_loo <- loo_compare(loo1, loo2)
  diff_waic <- loo_compare(waic1, waic2)

  cat("\nWAIC estimates: \n",
      paste("object1: ", round( waic1$estimates[3,1], 3) ), "\n",
      paste("object2: ", round( waic2$estimates[3,1], 3) ), "\n" )
  
  cat("\nWAIC difference & SE: \n", 
      sprintf("%8.3f", diff_waic[2, 1]), 
      sprintf("%8.3f", diff_waic[2, 2]), "\n")
  
  cat("\nLOO estimates: \n",
      paste("object1: ", round( loo1$estimates[3,1], 3) ), "\n",
      paste("object2: ", round( loo2$estimates[3,1], 3) ), "\n" )
  
  cat("\nLOO difference & SE: \n", 
      sprintf("%8.3f", diff_loo[2, 1]), 
      sprintf("%8.3f", diff_loo[2, 2]), "\n\n")

  cat("Laplace approximation to the log-Bayes factor\n(experimental; positive values favor object1):",
      sprintf("%8.3f", bf), "\n\n")

  looobj <- list(loo1$estimates, loo2$estimates)
  waicobj <- list(waic1$estimates, waic2$estimates)

  res <- list(bf = res, loo = looobj,
              diff_loo = diff_loo[2,1:2],
              waic = waicobj,
              diff_waic = diff_waic[2,1:2])

  invisible(res)
}

BF <- function(object1, object2, ...) {
  cat("BF() is deprecated. Use blavCompare() instead.\n\n")

  blavCompare(object1, object2)
}
