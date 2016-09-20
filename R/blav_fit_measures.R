setMethod("fitMeasures", signature(object = "blavaan"),
function(object, fit.measures = "all", baseline.model = NULL) {
    blav_fit_measures(object = object, fit.measures = fit.measures,
                     baseline.model = baseline.model)
})

setMethod("fitmeasures", signature(object = "blavaan"),
function(object, fit.measures = "all", baseline.model = NULL) {
    blav_fit_measures(object = object, fit.measures = fit.measures,
                     baseline.model = baseline.model)
})


#fitMeasures <- fitmeasures <- function(object, fit.measures="all") {
blav_fit_measures <- function(object, fit.measures = "all", 
                              baseline.model = NULL) {

    # has the model converged?
    if(object@Fit@npar > 0L && !object@optim$converged) {
        warning("blavaan WARNING: the chains may not have converged.")
    }

    # do we have a test statistic?
    if(object@Options$test == "none") {
        stop("blavaan ERROR: fit measures cannot be obtained when test=\"none\"")
    }

    if("all" %in% fit.measures) {
       class.flag <- TRUE
    } else {
       class.flag <- FALSE
    }

    # collect info from the lavaan slots
    GLIST <- object@Model@GLIST

    # N versus N-1
    # this affects BIC, RMSEA, cn_01/05, MFI and ECVI
    # Changed 0.5-15: suggestion by Mark Seeto
    if(object@Options$estimator %in% c("ML","PML","FML") && 
       object@Options$likelihood == "normal") {
        N <- object@SampleStats@ntotal
    } else {
        N <- object@SampleStats@ntotal - object@SampleStats@ngroups
    }

    # Change 0.5-13: take into account explicit equality constraints!!
    # reported by Mark L. Taper (affects AIC and BIC)
    npar <- object@Fit@npar
    if(nrow(object@Model@con.jac) > 0L) {
        ceq.idx <- attr(object@Model@con.jac, "ceq.idx")
        if(length(ceq.idx) > 0L) {
            neq <- qr(object@Model@con.jac[ceq.idx,,drop=FALSE])$rank
            npar <- npar - neq
        }
    }

    fx <- object@Fit@fx
    fx.group <- object@Fit@fx.group
    meanstructure <- object@Model@meanstructure
    categorical   <- object@Model@categorical
    multigroup    <- object@Data@ngroups > 1L
    estimator     <- "ML" #object@Options$estimator
    test          <- object@Options$test
    G <- object@Data@ngroups  # number of groups
    ## TODO get these from something like get_ll
    ##X2 <- object@Fit@test[[1]]$stat
    ##df <- object@Fit@test[[1]]$df

    # define 'sets' of fit measures:
    fit.always <- c("npar")

    # logl at posterior means
    fit.logl <- "logl"
    
    # posterior predictive p
    fit.chisq <- "ppp" #, "chisq.pi.lower", "chisq.pi.upper")

    # information criteria
    fit.ic <- c("bic", "dic", "p_dic")

    # check if "loo" is loaded or not
    # pkgs <- names(sessionInfo()[["otherPkgs"]])
    # if( "loo" %in% pkgs ) {
    #if(require("loo")) fit.ic <- c(fit.ic, "waic", "p_waic", "looic", "p_loo")
        fit.ic <- c(fit.ic, "waic", "p_waic", "looic", "p_loo")
    #}

    fit.rmsea <- "rmsea"
    
    # marginal logl
    fit.mll <- "margloglik"

    # lower case
    fit.measures <- tolower(fit.measures)

    # select 'default' fit measures
    if(length(fit.measures) == 1L) {
        if(fit.measures == "default" | fit.measures == "all") {
            if(estimator == "ML") {
                fit.measures <- c(fit.always, fit.logl, fit.chisq, fit.ic, 
                                  fit.rmsea, fit.mll)
            } else {
                stop("blavaan ERROR: Invalid estimator.")

                #if(estimator == "MML") {
                #  fit.measures <- c(fit.always, fit.logl)
            }

        }
    }

    # main container
    indices <- list()

    if("npar" %in% fit.measures) {
        indices["npar"] <- npar
    }

    if("logl" %in% fit.measures) {
        indices["logl"] <- object@Fit@fx
    }
    
    # posterior predictive p
    if("ppp" %in% fit.measures) {
        indices["ppp"] <- object@Fit@test[[2]]$stat
    }
    if(any(c("bic", "dic", "p_dic") %in% fit.measures)) {
        df <- 2*(object@Fit@fx - mean(as.numeric(object@external$runjags$samplls[,,1])))
        indices["bic"] <- -2*object@Fit@fx + npar*log(N)
        indices["dic"] <- -2*object@Fit@fx + 2*df
        indices["p_dic"] <- df
    }
    if(any(c("waic", "p_waic", "looic", "p_loo") %in% fit.measures)) {
        lavopt <- object@Options
        lavopt$estimator <- "ML"
        casells <- case_lls(object@external$runjags, object@Model,
                            object@ParTable, object@SampleStats,
                            lavopt, object@Cache,
                            object@Data)
        fitres <- waic(casells)
        indices["waic"] <- fitres$waic
        indices["p_waic"] <- fitres$p_waic
        fitres <- loo(casells)
        indices["looic"] <- fitres$looic
        indices["p_loo"] <- fitres$p_loo
    }
    if("margloglik" %in% fit.measures) {
        indices["margloglik"] <- object@test[[1]]$stat
    }
    
    out <- unlist(indices[fit.measures])

    if(length(out) > 0L) {
        class(out) <- c("lavaan.vector", "numeric")
    } else {
        return( invisible(numeric(0)) )
    }

    out
}

# print a nice summary of the fit measures
## print.fit.measures <- function(x) {

##    names.x <- names(x)

##    # scaled?
##    scaled <- FALSE #"chisq.scaled" %in% names.x

##    # table fit measures
##    if("C_F" %in% names.x) {
##        cat("\nFull response patterns fit statistics:\n\n")

##        t0.txt <- sprintf("  %-40s", "Observed response patterns (1st group):")
##        t1.txt <- sprintf("  %10i", x["rpat.observed"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        t0.txt <- sprintf("  %-40s", "Total response patterns (1st group):")
##        t1.txt <- sprintf("  %10i", x["rpat.total"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        t0.txt <- sprintf("  %-40s", "Empty response patterns (1st group):")
##        t1.txt <- sprintf("  %10i", x["rpat.empty"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        cat("\n")

##        t0.txt <- sprintf("  %-40s", "C_F Test Statistic")
##        t1.txt <- sprintf("  %10.3f", x["C_F"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        t0.txt <- sprintf("  %-40s", "Degrees of freedom")
##        t1.txt <- sprintf("  %10i", x["C_F.df"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        t0.txt <- sprintf("  %-40s", "P-value")
##        t1.txt <- sprintf("  %10.3f", x["C_F.p.value"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        cat("\n")
##        t0.txt <- sprintf("  %-40s", "C_M Test Statistic")
##        t1.txt <- sprintf("  %10.3f", x["C_M"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        t0.txt <- sprintf("  %-40s", "Degrees of freedom")
##        t1.txt <- sprintf("  %10i", x["C_M.df"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        t0.txt <- sprintf("  %-40s", "P-value")
##        t1.txt <- sprintf("  %10.3f", x["C_M.p.value"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##    }

##    if("C_p" %in% names.x) {
##        cat("\nPairwise tables summary statistic:\n\n")
##        t0.txt <- sprintf("  %-40s", "C_P Test Statistic")
##        t1.txt <- sprintf("  %10.3f", x["C_p"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        t0.txt <- sprintf("  %-40s", "Degrees of freedom")
##        t1.txt <- sprintf("  %10i", x["C_p.df"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        t0.txt <- sprintf("  %-40s", "Bonferroni corrected P-value")
##        t1.txt <- sprintf("  %10.3f", x["C_p.p.value"])
##        t2.txt <- ""
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##    }

##     # cfi/tli
##     if(any(c("cfi","tli","nnfi","rfi","nfi","ifi","rni","pnfi") %in% names.x)) {
##         cat("\nUser model versus baseline model:\n\n")

##         if("cfi" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Comparative Fit Index (CFI)")
##             t1.txt <- sprintf("  %10.3f", x["cfi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["cfi.scaled"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }

##         if("tli" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Tucker-Lewis Index (TLI)")
##             t1.txt <- sprintf("  %10.3f", x["tli"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["tli.scaled"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }

##         if("nnfi" %in% names.x) {
##             t0.txt <- sprintf("  %-42s", "Bentler-Bonett Non-normed Fit Index (NNFI)")
##             t1.txt <- sprintf("  %8.3f", x["nnfi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["nnfi.scaled"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
 
##         if("nfi" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Bentler-Bonett Normed Fit Index (NFI)")
##             t1.txt <- sprintf("  %10.3f", x["nfi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["nfi.scaled"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
  
##         if("nfi" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Parsimony Normed Fit Index (PNFI)")
##             t1.txt <- sprintf("  %10.3f", x["pnfi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["pnfi.scaled"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }

##         if("rfi" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Bollen's Relative Fit Index (RFI)")
##             t1.txt <- sprintf("  %10.3f", x["rfi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["rfi.scaled"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }

##         if("ifi" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Bollen's Incremental Fit Index (IFI)")
##             t1.txt <- sprintf("  %10.3f", x["ifi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["ifi.scaled"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }

##         if("rni" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Relative Noncentrality Index (RNI)")
##             t1.txt <- sprintf("  %10.3f", x["rni"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["rni.scaled"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##     }

##    # likelihood
##    if("logl" %in% names.x) {
##        cat("\nLoglikelihood and Information Criteria:\n\n")
##        t0.txt <- sprintf("  %-40s", "Loglikelihood user model (H0)")
##        t1.txt <- sprintf("  %10.3f", x["logl"])
##        t2.txt <- ifelse(scaled,
##                  sprintf("  %10.3f", x["logl"]), "")
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        #cat(t0.txt, t1.txt, "\n", sep="")
##        if(!is.na(x["scaling.factor.h0"])) {
##            t0.txt <- sprintf("  %-40s", "Scaling correction factor")
##            t1.txt <- sprintf("  %10s", "")
##            t2.txt <- sprintf("  %10.3f", x["scaling.factor.h0"])
##            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##            cat("    for the MLR correction\n")
##        }

##        if("unrestricted.logl" %in% names.x) {
##            t0.txt <- sprintf("  %-40s", "Loglikelihood unrestricted model (H1)")
##            t1.txt <- sprintf("  %10.3f", x["unrestricted.logl"])
##            t2.txt <- ifelse(scaled,
##                      sprintf("  %10.3f", x["unrestricted.logl"]), "")
##            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##            #cat(t0.txt, t1.txt, "\n", sep="")
##            if(!is.na(x["scaling.factor.h1"])) {
##                t0.txt <- sprintf("  %-40s", "Scaling correction factor")
##                t1.txt <- sprintf("  %10s", "")
##                t2.txt <- sprintf("  %10.3f", x["scaling.factor.h1"])
##                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##               cat("    for the MLR correction\n")
##            }
##        }

##        cat("\n")
##        t0.txt <- sprintf("  %-40s", "Number of free parameters")
##        t1.txt <- sprintf("  %10i", x["npar"])
##        t2.txt <- ifelse(scaled,
##                  sprintf("  %10i", x["npar"]), "")
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

##        t0.txt <- sprintf("  %-40s", "Akaike (AIC)")
##        t1.txt <- sprintf("  %10.3f", x["aic"])
##        t2.txt <- ifelse(scaled,
##                  sprintf("  %10.3f", x["aic"]), "")
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        #cat(t0.txt, t1.txt, "\n", sep="")
##        t0.txt <- sprintf("  %-40s", "Bayesian (BIC)")
##        t1.txt <- sprintf("  %10.3f", x["bic"])
##        t2.txt <- ifelse(scaled,
##                  sprintf("  %10.3f", x["bic"]), "")
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        #cat(t0.txt, t1.txt, "\n", sep="")
##        if(!is.na(x["bic2"])) {
##            t0.txt <- sprintf("  %-40s", "Sample-size adjusted Bayesian (BIC)")
##            t1.txt <- sprintf("  %10.3f", x["bic2"])
##            t2.txt <- ifelse(scaled,
##                      sprintf("  %10.3f", x["bic2"]), "")
##            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##    }

##    # RMSEA
##    if("rmsea" %in% names.x) {
##        cat("\nRoot Mean Square Error of Approximation:\n\n")
##        t0.txt <- sprintf("  %-40s", "RMSEA")
##        t1.txt <- sprintf("  %10.3f", x["rmsea"])
##        t2.txt <- ifelse(scaled,
##                  sprintf("  %10.3f", x["rmsea.scaled"]), "")
##        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        if("rmsea.ci.lower" %in% names.x) {
##            t0.txt <- sprintf("  %-38s", "90 Percent Confidence Interval")
##            t1.txt <- sprintf("  %5.3f", x["rmsea.ci.lower"])
##            t2.txt <- sprintf("  %5.3f", x["rmsea.ci.upper"])
##            t3.txt <- ifelse(scaled,
##                      sprintf("       %5.3f  %5.3f", x["rmsea.ci.lower.scaled"],
##                                                x["rmsea.ci.upper.scaled"]), "")
##            cat(t0.txt, t1.txt, t2.txt, t3.txt, "\n", sep="")
##        }
##        if("rmsea.pvalue" %in% names.x) {
##            t0.txt <- sprintf("  %-40s", "P-value RMSEA <= 0.05")
##            t1.txt <- sprintf("  %10.3f", x["rmsea.pvalue"])
##            t2.txt <- ifelse(scaled,
##                  sprintf("  %10.3f", x["rmsea.pvalue.scaled"]), "")
##            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##        }
##    }

##     # SRMR
##     if(any(c("rmr","srmr") %in% names.x)) {
##         cat("\nStandardized Root Mean Square Residual:\n\n")

##         if("rmr" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "RMR")
##             t1.txt <- sprintf("  %10.3f", x["rmr"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["rmr"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##         if("rmr_nomean" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "RMR (No Mean)")
##             t1.txt <- sprintf("  %10.3f", x["rmr_nomean"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["rmr_nomean"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##         if("srmr" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "SRMR")
##             t1.txt <- sprintf("  %10.3f", x["srmr"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["srmr"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##         if("srmr_nomean" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "SRMR (No Mean)")
##             t1.txt <- sprintf("  %10.3f", x["srmr_nomean"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["srmr_nomean"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##     }

##     # WRMR
##     if("wrmr" %in% names.x) {
##         cat("\nWeighted Root Mean Square Residual:\n\n")

##         if("wrmr" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "WRMR")
##             t1.txt <- sprintf("  %10.3f", x["wrmr"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["wrmr"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##     }

##     # Other
##     if(any(c("cn_05","cn_01","gfi","agfi","pgfi","mfi") %in% names.x)) {
##         cat("\nOther Fit Indices:\n\n")

##         if("cn_05" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Hoelter Critical N (CN) alpha=0.05")
##             t1.txt <- sprintf("  %10.3f", x["cn_05"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["cn_05"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##         if("cn_01" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Hoelter Critical N (CN) alpha=0.01")
##             t1.txt <- sprintf("  %10.3f", x["cn_01"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["cn_01"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##         if(any(c("cn_05", "cn_01") %in% names.x)) {
##             cat("\n")
##         }
##         if("gfi" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Goodness of Fit Index (GFI)")
##             t1.txt <- sprintf("  %10.3f", x["gfi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["gfi"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##         if("agfi" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Adjusted Goodness of Fit Index (AGFI)")
##             t1.txt <- sprintf("  %10.3f", x["agfi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["agfi"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##         if("pgfi" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Parsimony Goodness of Fit Index (PGFI)")
##             t1.txt <- sprintf("  %10.3f", x["pgfi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["pgfi"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##         if(any(c("gfi","agfi","pgfi") %in% names.x)) {
##             cat("\n")
##         }
##         if("mfi" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "McDonald Fit Index (MFI)")
##             t1.txt <- sprintf("  %10.3f", x["mfi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["mfi"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }
##         if("mfi" %in% names.x) {
##             cat("\n")
##         }
##         if("ecvi" %in% names.x) {
##             t0.txt <- sprintf("  %-40s", "Expected Cross-Validation Index (ECVI)")
##             t1.txt <- sprintf("  %10.3f", x["ecvi"])
##             t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["ecvi"]), "")
##             cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
##         }

##     }

##     #cat("\n")
## }


