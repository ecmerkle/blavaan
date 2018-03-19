### Mauricio Garnier-Villareal and Terrence D. Jorgensen
### Last updated: 19 March 2018
### functions applying traditional SEM fit criteria to Bayesian models.
### Inspired by, and adpated from, Rens van de Schoot's idea for BRMSEA, as
### published in http://dx.doi.org/10.1177/0013164417709314


## Public function (eventually, after documenting and peer review)
BayesFIT <- function(fit, rescale = c("devM","ppmc"), fit_null = NULL) {
  rescale <- as.character(rescale[1])
  if (rescale != "devM") stop('Only rescale="devM" is currently available')
  if (rescale == "ppmc" && blavInspect(fit, 'ntotal') < 1000)
    warning("Hoofs et al.'s proposed BRMSEA (and derivative indices based on",
            " the posterior predictive distribution) was only proposed for",
            " evaluating models fit to very large samples (N > 1000).")
  
  chisqs <- as.numeric(apply(fit@external$samplls, 2,
                             function(x) 2*(x[,2] - x[,1])))
  fit_pd <- fitMeasures(fit, 'p_loo')
  # if (rescale = "ppmc") {
  #   reps <- as.numeric(apply(fit@external$replls, 2,
  #                            function(x) 2*(x[,2] - x[,1])))
  # } else reps <- NULL
  
  if (is.null(fit_null)) {
    null_model <- FALSE
    chisq_null <- NULL
    # reps_null <- NULL
    pD_null <- NULL
  } else {
    null_model <- TRUE
    chisq_null <- as.numeric(apply(fit_null@external$samplls, 2,
                                   function(x) 2*(x[,2] - x[,1])))
    if (length(chisqs) != length(chisq_null)) {
      null_model <- FALSE
      chisq_null <- NULL
      # reps_null <- NULL
      pD_null <- NULL
      warning("Incremental fit indices were not calculated.",
              " Save equal number of draws from the posterior of both",
              " the hypothesized and null models.")
    } else {
      # if (rescale = "ppmc") {
      #   reps_null <- as.numeric(apply(fit_null@external$replls, 2,
      #                                 function(x) 2*(x[,2] - x[,1])))
      # }
      pD_null <- fitMeasures(fit_null, 'p_loo')
    }
  }

  ff <- BayesRelFit(obs = chisqs, # reps = reps,
                    nvar = fit@Model@nvar, pD = fit_pd,
                    N = blavInspect(fit, 'ntotal'), Min1 = FALSE,
                    ms = blavInspect(fit, 'meanstructure'),
                    Ngr = blavInspect(fit, 'ngroups'), null_model = null_model,
                    obs_null = chisq_null, # reps_null = reps_null,
                    pD_null = pD_null)
  return(ff)
}


### Hidden function to calculate Bayesian fit indices
### Rens: Bayesian RMSEA adapted from http://dx.doi.org/10.1177/0013164417709314
### MGV:  change the correction method, not longer like Rens
### TDJ:  added argument to select Rens' method ("ppmc") vs. ours ("devM")
### MGV:  modified to also calculate gammahat, adjusted gammahat
### TDJ:  added McDonald's centrality index
### MGV:  If a null model information provided, it calculates CFI, TLI, NFI
BayesRelFit <- function(obs, reps = NULL, nvar, pD, N, ms = TRUE, Min1 = FALSE,
                        Ngr = 1, rescale = c("devM","ppmc"), null_model = TRUE,
                        obs_null = NULL, reps_null = NULL, pD_null = NULL) {
  if (Min1) N <- N - 1
  rescale <- as.character(rescale[1])
  if (rescale == "devM") {
    reps <- pD
    if (!is.null(null_model)) reps_null <- pD_null
  } else stop('Only rescale="devM" is currently available')
  if (rescale == "ppmc" && (is.null(reps) || (null_model && is.null(reps_null))))
    stop('rescale="ppmc" requires non-NULL reps argument (and reps_null, if applicable).')
  
  ## Compute number of modeled moments
  p <- ((nvar * (nvar + 1)) / 2)
  if (ms) p <- p + nvar
  p <- p * Ngr
  ## Substract parameters and estimated parameters
  dif.ppD <- p - pD
  nonc <- obs - reps - dif.ppD # == obs - p when rescale == "devM" because reps = pD
  ## Correct if numerator is smaller than zero
  nonc[nonc < 0] <- 0
  ## Compute BRMSEA
  BRMSEA <- sqrt(nonc / (dif.ppD * N)) * sqrt(Ngr)
  
  ## compute GammaHat and adjusted GammaHat
  gammahat <- nvar / (nvar + 2*nonc/N)
  adjgammahat <- 1 - (p / dif.ppD) * (1 - gammahat)
  
  ## compute McDonald's centrality index
  Mc <- exp(-.5 * nonc/N)
  
  out <- cbind(chi_adj = obs - reps, df_adj = dif.ppD, BRMSEA = BRMSEA,
               BGammaHat = gammahat, adjBGammaHat = adjgammahat, BMc = Mc)
  
  ## calculate fit when null model is provided
  if (null_model) {
    dif.ppD_null <- p - pD_null
    nonc_null <- (obs_null - reps_null) - dif.ppD_null
    
    cfi <- 1 - (nonc / nonc_null)
    tli_null_part <- (obs_null - reps_null) / dif.ppD_null
    tli <- (tli_null_part - (obs - reps) / dif.ppD) / (tli_null_part - 1)
    nfi <- ((obs_null - reps_null) - (obs - reps)) / (obs_null - reps_null)
    
    out <- cbind(out, BCFI = cfi, BTLI = tli, BNFI = nfi)
  }
  
  class(out) <- c("lavaan.matrix","matrix")
  return(out)
}
