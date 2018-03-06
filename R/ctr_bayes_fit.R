## functions applying traditional SEM fit criteria to Bayesian models,
## from Mauricio Garnier-Villareal with code contributions from
## Terrence Jorgensen and Rens van de Schoot.

BayesFIT <- function(fit, fit_null=NULL){
  chisqs <- as.numeric(apply(fit@external$samplls, 2,
                             function(x) 2*(x[,2] - x[,1])))
  fit_pd <- fitMeasures(fit, 'p_loo')

  if( is.null(fit_null) ){
    null_model <- FALSE
    chisq_null <- NULL
    pD_null <- NULL
  } else {
    null_model <- TRUE
    chisq_null <- as.numeric(apply(fit_null@external$samplls, 2,
                                   function(x) 2*(x[,2] - x[,1])))
    pD_null <- fitMeasures(fit_null, 'p_loo')
  }

  ff <- BayesRelFit(obs=chisqs,
                    nvar=fit@Model@nvar,
                    pD=fit_pd, N=blavInspect(fit, 'ntotal'),
                    ms=blavInspect(fit, 'meanstructure'),
                    Min1=TRUE, Ngr=blavInspect(fit, 'ngroups'),
                    null_model=null_model, obs_null=chisq_null,
                    pD_null=pD_null)
  
  return(ff)
}


### Bayesian RMSEA from Rens
### MGV: change the correction method, not longer like Rens
### MGV: modified to also calculate gammahat, adjusted gammahat
### TDJ: added McDonald's centrality index
### if a null model information provided, it calculates CFI, TLI, NFI
BayesRelFit <-function(obs, nvar, pD, N, ms = TRUE, Min1 = TRUE,
                       Ngr = 1, null_model=TRUE, obs_null=NULL,
                       pD_null=NULL){

  ## Compute number of modeled moments
  if(ms) p <- (((nvar * (nvar + 1)) / 2) + nvar)
  if(!ms) p <- (((nvar * (nvar + 1))/ 2) + 0)
  p <- p * Ngr
  # # Substract parameters and estimated parameters
  dif.ppD <- p - pD
  nonc <- obs - pD - dif.ppD # == obs - p
  # # Correct if numerator is smaller than zero
  nonc[nonc < 0] <- 0
  # # Compute BRMSEA (with or without the -1 correction)
  if(Min1)
    BRMSEA <- sqrt(nonc / (dif.ppD * (N - 1)))*sqrt(Ngr)
  if(!Min1) BRMSEA <- sqrt(nonc / (dif.ppD * N))*sqrt(Ngr)

  ## compute GammaHat and adjusted GammaHat
  gammahat <- nvar / ( nvar+2* ((nonc)/(N-1))  )
  adjgammahat <- 1 - (((Ngr * nvar * (nvar + 1))/2)/dif.ppD) * (1 - gammahat)
  
  ## compute McDonald's centrality index
  Mc <- exp(-.5 * nonc/(N-1) )
  
  ## calculate fit when null model is provided
  if (null_model) {
    dif.ppD_null <- p - pD_null
    nonc_null <- ( ( obs_null-pD_null ) - dif.ppD_null )
    
    cfi <- (nonc_null - nonc)/nonc_null 
    tli <- ((( obs_null-pD_null )/dif.ppD_null) - (( obs-pD )/dif.ppD)) / (((( obs_null-pD_null )/dif.ppD_null))-1) 
    nfi <- (( obs_null-pD_null ) - ( obs-pD )) / ( obs_null-pD_null )
    
    out <- cbind(chi_adj=obs-pD, df_adj=dif.ppD, 
                 BRMSEA=BRMSEA, BGammaHat=gammahat, adjBGammaHat=adjgammahat, 
                 BMc = Mc, BCFI=cfi, BTLI=tli, BNFI=nfi)
  } else {
    out <- cbind(chi_adj=obs-pD, df_adj=dif.ppD,
                 BRMSEA=BRMSEA, BGammaHat=gammahat, adjBGammaHat=adjgammahat, BMc = Mc)
  }
  
  return(out)
}
