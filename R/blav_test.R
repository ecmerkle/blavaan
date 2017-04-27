blav_model_test <- function(lavmodel       = NULL, 
                            lavpartable    = NULL, 
                            lavsamplestats = NULL, 
                            lavoptions     = NULL, 
                            x              = NULL, 
                            VCOV           = NULL, 
                            lavcache       = NULL,
                            lavdata        = NULL,
                            lavjags        = NULL,
                            samplls        = NULL,
                            jagextra       = NULL,
                            stansumm       = NULL,
                            control        = list()) {


    TEST <- list()

    ## marginal log-likelihood approximation
    ## needs original partable with rhos
    if("syntax" %in% names(jagextra)){
        warning("blavaan WARNING: Marginal log-likelihood cannot be approximated when there is additional JAGS syntax.")
        mll <- NA
    } else {
        mll <- try(margloglik(lavpartable, lavmodel, lavoptions, 
                              lavsamplestats, lavdata, lavcache,
                              lavjags, VCOV, x, stansumm),
                   silent=TRUE)
        if(inherits(mll, "try-error")) mll <- NA
    }

    ppp <- postpred(lavpartable, lavmodel, lavoptions,
                    lavsamplestats, lavdata, lavcache, lavjags,
                    samplls)$ppval

    TEST[[1]] <- list(test="mloglik",
                      stat=as.numeric(mll),
                      stat.group=as.numeric(NA),
                      df=as.integer(NA),
                      refdistr="NA",
                      pvalue=as.numeric(NA))

    TEST[[2]] <- list(test="ppp",
                      ## DIC: 2*ll(theta_hat) - 4*mean(ll(theta_samp))
                      stat=as.numeric(ppp),
                      stat.group=as.numeric(NA),
                      df=as.integer(NA),
                      refdistr="NA",
                      pvalue=as.numeric(NA))

    TEST
}
