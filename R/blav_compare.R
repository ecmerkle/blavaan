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
  
  if(targ1 == "stan" && blavInspect(object1, "meanstructure")){
    ll1 <- loo::extract_log_lik(object1@external$mcmcout)
  } else if(blavInspect(object1, "categorical") && lavopt1$test != "none"){
    if("llnsamp" %in% names(lavopt1)){
      cat("blavaan NOTE: These criteria involve likelihood approximations that may be imprecise.\n",
          "You could try running the model again to see how much the criteria fluctuate.\n",
          "You can also manually set llnsamp for greater accuracy (but also greater runtime).\n\n")
    }
    ll1 <- object1@external$casells
  } else {
    lavopt1$estimator <- "ML"
    ll1 <- case_lls(object1@external$mcmcout, make_mcmc(object1@external$mcmcout),
                    object1)
  }
  nchain1 <- blavInspect(object1, "n.chains")
  niter1 <- nrow(ll1)/nchain1
  cid1 <- rep(1:nchain1, each=niter1)
  ref1 <- relative_eff(exp(ll1), chain_id = cid1)

  if(targ2 == "stan" && blavInspect(object2, "meanstructure")){
    ll2 <- loo::extract_log_lik(object2@external$mcmcout)
  } else if(blavInspect(object2, "categorical") && lavopt2$test != "none"){
    if("llnsamp" %in% names(lavopt2)){
      cat("blavaan NOTE: These criteria involve likelihood approximations that may be imprecise.\n",
          "You could try running the model again to see how much the criteria fluctuate.\n",
          "You can also manually set llnsamp for greater accuracy (but also greater runtime).\n\n")
    }
    ll2 <- object2@external$casells
  } else {
    lavopt2$estimator <- "ML"
    ll2 <- case_lls(object2@external$mcmcout, make_mcmc(object2@external$mcmcout),
                    object2)
  }  
  nchain2 <- blavInspect(object1, "n.chains")
  niter2 <- nrow(ll2)/nchain2
  cid2 <- rep(1:nchain2, each=niter2)
  ref2 <- relative_eff(exp(ll2), chain_id = cid2)

  loo1 <- loo(ll1, r_eff=ref1, ...)
  loo2 <- loo(ll2, r_eff=ref2, ...)
  waic1 <- waic(ll1); waic2 <- waic(ll2)

  diff_loo <- loo_compare(loo1, loo2)
  diff_waic <- loo_compare(waic1, waic2)

  cat("\nWAIC estimates: \n",
      paste("object1: ", round( waic1$estimates[3,1], 3) ), "\n",
      paste("object2: ", round( waic2$estimates[3,1], 3) ), "\n" )
  
  cat("\n ELPD difference & SE: \n", 
      sprintf("%8.3f", diff_waic[2, 1]), 
      sprintf("%8.3f", diff_waic[2, 2]), "\n")
  
  cat("\nLOO estimates: \n",
      paste("object1: ", round( loo1$estimates[3,1], 3) ), "\n",
      paste("object2: ", round( loo2$estimates[3,1], 3) ), "\n" )
  
  cat("\n ELPD difference & SE: \n", 
      sprintf("%8.3f", diff_loo[2, 1]), 
      sprintf("%8.3f", diff_loo[2, 2]), "\n\n")

  cat("Laplace approximation to the log-Bayes factor\n(experimental; positive values favor object1):",
      sprintf("%8.3f", bf), "\n\n")

  looobj <- list(loo1, loo2)
  waicobj <- list(waic1, waic2)

  res <- list(bf = res, loo = looobj,
              diff_loo = diff_loo,
              waic = waicobj,
              diff_waic = diff_waic)

  invisible(res)
}

BF <- function(object1, object2, ...) {
  cat("BF() is deprecated. Use blavCompare() instead.\n\n")

  blavCompare(object1, object2)
}
