### Mauricio Garnier-Villareal and Terrence D. Jorgensen
### Last updated: 12 April 2018
### functions applying traditional SEM fit criteria to Bayesian models.
### Inspired by, and adpated from, Rens van de Schoot's idea for BRMSEA, as
### published in http://dx.doi.org/10.1177/0013164417709314


## -----------------
## Class and Methods
## -----------------

setClass("blavFitIndices",
         representation(details = "list",
                        # chisq = "numeric", df = "numeric", # Bayesian analogs
                        # pD = "character", rescale = "character", # store choices
                        # list of vectors, each storing the posterior of 1 index
                        indices = "list"))

summary.blavFitIndices <- function(object,
                                   central.tendency = c("mean","median","mode"),
                                   hpd = TRUE, prob = .90) {
  if (!is.character(central.tendency)) {
    stop('blavaan ERROR: central.tendency must be a character vector')
  }
  central.tendency <- tolower(central.tendency)
  invalid <- setdiff(central.tendency, c("mean","median","mode","eap","map"))
  if (length(invalid)) {
    stop('blavaan ERROR: Invalid choice(s) of central.tendency:\n"',
         paste(invalid, collapse = '", "'), '"')
  }
  ## collect list of functions to apply to each index
  FUN <- function(x) {
    out <- numeric(0)
    
    if ("mean" %in% central.tendency || "eap" %in% central.tendency) {
      out <- c(out, EAP = mean(x, na.rm = TRUE))
    }
    
    if ("median" %in% central.tendency) {
      out <- c(out, Median = median(x, na.rm = TRUE))
    }
    
    if ("mode" %in% central.tendency || "map" %in% central.tendency) {
      ## can the modeest package be used?
      if (suppressMessages(requireNamespace("modeest", quietly = TRUE))) {
        out <- c(out, MAP =  modeest::mlv(x, method = "kernel", na.rm = TRUE)$M)
      } else {
        ## if not, use the quick-and-dirty way
        dd <- density(x, na.rm = TRUE)
        out <- c(out, MAP = dd$x[which.max(dd$y)])
      }
    }
    
    out <- c(out, SD = sd(x, na.rm = TRUE))
    
    if (hpd) {
      if (!"package:coda" %in% search()) attachNamespace("coda")
      out <- c(out, HPDinterval(as.mcmc(x), prob = prob)[1, ] )
    }
    
    out
  }
  ## apply function to each fit index
  indexSummaries <- lapply(object@indices, FUN)
  out <- as.data.frame(do.call(rbind, indexSummaries))
  class(out) <- c("lavaan.data.frame","data.frame")
  attr(out, "header") <- paste0("Posterior summary statistics and highest",
                                " posterior density (HPD) ", round(prob*100, 2),
                                "% credible intervals for ",
                                object@details$rescale,
                                "-based fit indices:", sep = "")
  out
}
setMethod("summary", "blavFitIndices", summary.blavFitIndices)

setMethod("show", "blavFitIndices", function(object) {
  cat("Posterior mean (EAP) of ", object@details$rescale,
      "-based fit indices:\n\n", sep = "")
  fits <- getMethod("summary","blavFitIndices")(object, "mean", hpd = FALSE)
  out <- fits$EAP
  names(out) <- rownames(fits)
  class(out) <- c("lavaan.vector","numeric")
  print(out)
  invisible(object)
})



## --------------------
## Constructor Function
## --------------------

## Public function (eventually, after documenting and peer review)
blavFitIndices <- function(object, pD = c("loo","waic","dic"),
                           rescale = c("devM","ppmc","mcmc"),
                           fit.measures = "all", baseline.model = NULL) {
  ## check arguments
  pD <- tolower(as.character(pD[1]))
  if (!pD %in% c("loo","waic","dic")) {
    stop('blavaan ERROR: Invalid choice of "pD" argument')
  }

  rescale <- tolower(as.character(rescale[1])) #FIXME: make sure references to "devM" check for "devm"
  if (rescale == "ppmc" && blavInspect(object, 'ntotal') < 1000)
    warning("blavaan WARNING: ",
            "Hoofs et al.'s proposed BRMSEA (and derivative indices based on",
            " the posterior predictive distribution) was only proposed for",
            " evaluating models fit to very large samples (N > 1000).")
  
  chisqs <- as.numeric(apply(object@external$samplls, 2,
                             function(x) 2*(x[,2] - x[,1])))
  fit_pd <- fitMeasures(object, paste0('p_', pD))
  if (rescale == "mcmc") {
    stop('blavaan ERROR: type = "mcmc" not implemented yet.')
    ##FIXME TDJ: Call postpred(), just save different information.
    ## if (!is.null(baseline.model))
    ##      - Pass baseline model (parTable?) as argument? (priors?!)
    ##      - Or calculate using baseline posterior chi-squared?
    ## ff <- postpred(...)$... # extract discFUN() output
    ## out <- new("blavFitIndices", details = list(rescale = rescale,
    ##            customized.baseline = !is.null(baseline.model)), indices = ff)
  } else if (rescale == "ppmc") {
    reps <- postpred(lavpartable = object@ParTable,
                     lavmodel = object@Model,
                     lavoptions = object@Options,
                     lavsamplestats = object@SampleStats,
                     lavdata = object@Data,
                     lavcache = object@Cache,
                     lavjags = object@external$mcmcout,
                     samplls = object@external$samplls)$chisqs[,"reps"]
  } else reps <- NULL
  
  if (is.null(baseline.model)) {
    null_model <- FALSE
    chisq_null <- NULL
    reps_null <- NULL
    pD_null <- NULL
  } else if (rescale != "mcmc") {
    null_model <- TRUE
    chisq_null <- as.numeric(apply(baseline.model@external$samplls, 2,
                                   function(x) 2*(x[,2] - x[,1])))
    if (length(chisqs) != length(chisq_null)) {
      null_model <- FALSE
      chisq_null <- NULL
      reps_null <- NULL
      pD_null <- NULL
      warning("Incremental fit indices were not calculated.",
              " Save equal number of draws from the posterior of both",
              " the hypothesized and null models.")
    } else {
      if (rescale == "ppmc") {
        reps_null <- postpred(lavpartable = baseline.model@ParTable,
                              lavmodel = baseline.model@Model,
                              lavoptions = baseline.model@Options,
                              lavsamplestats = baseline.model@SampleStats,
                              lavdata = baseline.model@Data,
                              lavcache = baseline.model@Cache,
                              lavjags = baseline.model@external$mcmcout,
                              samplls = baseline.model@external$samplls)$chisqs[,"reps"]
      }
      pD_null <- fitMeasures(baseline.model, 'p_loo')
    }
  }
  
  if (rescale != "mcmc") {
    out <- BayesChiFit(obs = chisqs, reps = reps,
                       nvar = object@Model@nvar, pD = fit_pd,
                       N = blavInspect(object, 'ntotal'),
                       Ngr = blavInspect(object, 'ngroups'),
                       Min1 = blavInspect(object, 'options')$mimic == "EQS",
                       ms = blavInspect(object, 'meanstructure'),
                       rescale = rescale, fit.measures = fit.measures,
                       null_model = null_model, obs_null = chisq_null,
                       reps_null = reps_null, pD_null = pD_null)
  }
  
  nChains <- blavInspect(object, 'n.chains')
  out@details <- c(out@details, list(n.chains = nChains))
  out
}



## ----------------
## Hidden functions
## ----------------

### function to calculate fit indices from posterior distribution chi-squared
### Rens: Bayesian RMSEA adapted from http://dx.doi.org/10.1177/0013164417709314
### MGV:  change the correction method, not longer like Rens
### TDJ:  added argument to select Rens' method ("ppmc") vs. ours ("devM")
### MGV:  modified to also calculate gammahat, adjusted gammahat
### TDJ:  added McDonald's centrality index
### MGV:  If a null model information provided, it calculates CFI, TLI, NFI
BayesChiFit <- function(obs, reps = NULL, nvar, pD, N, Ngr = 1,
                        ms = TRUE, Min1 = FALSE,
                        rescale = c("devM","ppmc"), fit.measures = "all",
                        null_model = TRUE, obs_null = NULL,
                        reps_null = NULL, pD_null = NULL) {
  if (!is.character(fit.measures)) {
    stop('blavaan ERROR: fit.measures must be a character vector')
  }
  fit.measures <- tolower(fit.measures)
  if (fit.measures == "all") {
    fit.measures <- c("brmsea","bgammahat","adjbgammahat","bmc")
    if (null_model) fit.measures <- c(fit.measures, "bcfi","btli","bnfi")
  }
  
  if (Min1) N <- N - Ngr
  
  rescale <- tolower(as.character(rescale[1]))
  if (rescale == "devm") {
    reps <- pD
    if (!is.null(null_model)) reps_null <- pD_null
  }
  if (rescale == "ppmc" && (is.null(reps) || (null_model && is.null(reps_null)))) {
    stop('blavaan ERROR: rescale="ppmc" requires non-NULL reps argument (and reps_null, if applicable).')
  }
  ##FIXME TDJ: postpred(discFUN = ) argument would enable type = "mcmc".
  ##           Just calculate any fit indices using each posterior sample.
  ##           Add fit.measures= argument here and above. How to handle CFI...?
  
  ## Compute number of modeled moments
  p <- ((nvar * (nvar + 1)) / 2)
  if (ms) p <- p + nvar
  p <- p * Ngr
  ## Substract parameters and estimated parameters
  dif.ppD <- p - pD
  nonc <- obs - reps - dif.ppD # == obs - p when rescale == "devm" because reps = pD
  ## Correct if numerator is smaller than zero
  nonc[nonc < 0] <- 0
  
  ## assemble results in a vector
  result <- list()
  
  ## Compute BRMSEA
  if ("brmsea" %in% fit.measures) {
    result[["BRMSEA"]] <- sqrt(nonc / (dif.ppD * N)) * sqrt(Ngr)
  }
  
  ## compute GammaHat and adjusted GammaHat
  if ("bgammahat" %in% fit.measures) {
    result[["BGammaHat"]] <- nvar / (nvar + 2*nonc/N)
    if ("adjbgammahat" %in% fit.measures) {
      result[["adjBGammaHat"]] <- 1 - (p / dif.ppD) * (1 - result[["BGammaHat"]])
    }
  } else if ("adjbgammahat" %in% fit.measures) {
    gammahat <- nvar / (nvar + 2*nonc/N)
    result[["adjBGammaHat"]] <- 1 - (p / dif.ppD) * (1 - gammahat)
  }

  ## compute McDonald's centrality index
  if ("bmc" %in% fit.measures) {
    result[["BMc"]] <- exp(-.5 * nonc/N)
  }

  ## calculate incremental fit when null model is provided
  if (null_model) {
    dif.ppD_null <- p - pD_null
    nonc_null <- (obs_null - reps_null) - dif.ppD_null
    
    if ("bcfi" %in% fit.measures) {
      result[["BCFI"]] <- 1 - (nonc / nonc_null)
    }
    if ("btli" %in% fit.measures) {
      tli_null_part <- (obs_null - reps_null) / dif.ppD_null
      result[["BTLI"]] <- (tli_null_part - (obs - reps) / dif.ppD) / (tli_null_part - 1)
    }
    if ("bnfi" %in% fit.measures) {
      result[["BNFI"]] <- ((obs_null - reps_null) - (obs - reps)) / (obs_null - reps_null)
    }
  }
  
  out <- new("blavFitIndices",
             details = list(chisq = obs - reps, df = dif.ppD,
                            pD = pD, rescale = rescale),
             indices = result)
  ## for print methods
  class(out@details$chisq) <- c("lavaan.vector","numeric")
  class(out@details$df) <- c("lavaan.vector","numeric")
  for (i in seq_along(out@indices)) {
    class(out@indices[[i]]) <- c("lavaan.vector","numeric")
  }
  
  out
}




