#
# initial version: YR 25/03/2009

short.summary <- function(object) {

    # catch FAKE run
    FAKE <- FALSE
    if(!is.null(object@Model@control$optim.method)) {
        if(tolower(object@Model@control$optim.method) == "none") {
            FAKE <- TRUE
        }
    }

    # Convergence or not?
    if(FAKE) {
        cat(sprintf("blavaan (%s) -- DRY RUN with 0 iterations\n",
                    packageDescription("blavaan", fields="Version")))
    } else if(object@Fit@iterations > 0) {
        if(object@Fit@converged) {
	    cat(sprintf("blavaan (%s) results of %3i samples after %3i adapt+burnin iterations\n",
                    packageDescription("blavaan", fields="Version"),
                    object@external$runjags$sample,
                    object@external$runjags$burnin))
        } else {
            cat(sprintf("** WARNING ** blavaan (%s) did NOT converge after %i adapt+burnin iterations\n", 
                packageDescription("blavaan", fields="Version"),
                object@external$runjags$burnin))
            cat("** WARNING ** Proceed with caution\n")
        }
    } else {
        cat(sprintf("** WARNING ** blavaan (%s) model has NOT been fitted\n",
                    packageDescription("blavaan", fields="Version")))
        #cat("** WARNING ** Estimates below are simply the starting values\n")
    }
    cat("\n")

    # number of free parameters
    #t0.txt <- sprintf("  %-40s", "Number of free parameters")
    #t1.txt <- sprintf("  %10i", object@Fit@npar)
    #t2.txt <- ""
    #cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    #cat("\n")
   
    # listwise deletion?
    listwise <- FALSE
    for(g in 1:object@Data@ngroups) {
       if(object@Data@nobs[[1L]] != object@Data@norig[[1L]]) {
           listwise <- TRUE
           break
       }
    }


    if(object@Data@ngroups == 1L) {
        if(listwise) {
            cat(sprintf("  %-40s", ""), sprintf("  %10s", "Used"), 
                                        sprintf("  %10s", "Total"),
                "\n", sep="")
        }
        t0.txt <- sprintf("  %-40s", "Number of observations")
        t1.txt <- sprintf("  %10i", object@Data@nobs[[1L]])
        t2.txt <- ifelse(listwise,
                  sprintf("  %10i", object@Data@norig[[1L]]), "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    } else {
        if(listwise) {
            cat(sprintf("  %-40s", ""), sprintf("  %10s", "Used"),  
                                        sprintf("  %10s", "Total"),
                "\n", sep="")
        }
        t0.txt <- sprintf("  %-40s", "Number of observations per group")
        cat(t0.txt, "\n")
        for(g in 1:object@Data@ngroups) {
            t.txt <- sprintf("  %-40s  %10i", object@Data@group.label[[g]],
                                              object@Data@nobs[[g]])
            t2.txt <- ifelse(listwise,
                      sprintf("  %10i", object@Data@norig[[g]]), "")
            cat(t.txt, t2.txt, "\n", sep="")
        }
    }
    cat("\n")

    # missing patterns?
    if(object@SampleStats@missing.flag) {
        if(object@Data@ngroups == 1L) {
            t0.txt <- sprintf("  %-40s", "Number of missing patterns")
            t1.txt <- sprintf("  %10i", 
                              object@Data@Mp[[1L]]$npatterns)
            cat(t0.txt, t1.txt, "\n\n", sep="")
        } else {
            t0.txt <- sprintf("  %-40s", "Number of missing patterns per group")
            cat(t0.txt, "\n")
            for(g in 1:object@Data@ngroups) {
                t.txt <- sprintf("  %-40s  %10i", object@Data@group.label[[g]],
                                 object@Data@Mp[[g]]$npatterns)
                cat(t.txt, "\n", sep="")
            }
            cat("\n")
        }
    }

    # Print Chi-square value for the user-specified (full/h0) model

    # other statistics?
    ppp <- TRUE
    if(length(object@Fit@test) == 0) ppp <- FALSE
    shifted <- FALSE

    # 0. heading
    #h.txt <- sprintf("\nChi-square test user model (h0)",
    #                 object@Options$estimator)
    t0.txt <- sprintf("  %-40s", "Statistic")
    t1.txt <- ifelse(ppp,
                     sprintf("  %10s", "MargLogLik"), "") #object@Options$estimator)
    t2.txt <- ifelse(ppp, 
              sprintf("  %10s", "PPP"), "")
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    # default output; change if request other fit stats?

    # 1. test statistics
    t0.txt <- sprintf("  %-40s", "Value")  
    t1.txt <- ifelse(ppp,
                     sprintf("  %10.3f", object@Fit@test[[1]]$stat), "")
    t2.txt <- ifelse(ppp, 
                     sprintf("  %10.3f", object@Fit@test[[2]]$stat), "")
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    # 2. effective number of parameters
    #t0.txt <- sprintf("  %-40s", "Effective # params")
    #t1.txt <- sprintf("  %10i",   object@Fit@test[[1]]$df)
    #t2.txt <- ifelse(dic, 
    #                     ifelse(round(object@Fit@test[[2]]$df) == 
    #                            object@Fit@test[[2]]$df,
    #                            sprintf("  %10i",   object@Fit@test[[2]]$df),
    #                            sprintf("  %10.3f", object@Fit@test[[2]]$df)),
    #                     "")
    #cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    # 3. P-value is skipped


    ## interesting way to add other stuff:
    ## if(object@Options$estimator == "MML") {
    ##     fm <- fitMeasures(object, c("logl", "npar", "aic", "bic", "bic2"))
    ##     print.fit.measures(fm)
    ## }

    #cat("\n")
}

setMethod("show", "blavaan",
function(object) {

    # show only basic information
    short.summary(object)

})

setMethod("summary", "blavaan",
function(object, header       = TRUE,
                 fit.measures = FALSE,
                 estimates    = TRUE,
                 ci           = TRUE,
                 standardized = FALSE,
                 rsquare      = FALSE,
                 std.nox      = FALSE,
                 modindices   = FALSE,
                 psrf         = TRUE,
                 neff         = FALSE,
                 postmedian   = FALSE,
                 postmode     = FALSE,
                 priors       = TRUE,
                 bf           = FALSE,
                 nd = 3L) {

    if(std.nox) standardized <- TRUE

    # print the 'short' summary
    if(header) {
        short.summary(object)
    }

    # only if requested, the fit measures
    if(fit.measures) {
        if(object@Options$test == "none") {
            warning("lavaan WARNING: fit measures not available if test = \"none\"\n\n")
        } else {
            #print.fit.measures( fitMeasures(object, fit.measures="default") )
            # TODO: define proper print function for custom fit measures
        }
    }


    if(estimates) {
        PE <- parameterEstimates(object, zstat = FALSE, ci = TRUE,
                                 standardized = standardized,
                                 rsquare = rsquare,
                                 remove.eq = FALSE, remove.system.eq = TRUE,
                                 remove.ineq = FALSE, remove.def = FALSE,
                                 add.attributes = TRUE)
        if(!("group" %in% names(PE))) PE$group <- 1
        if(standardized && std.nox) {
            PE$std.all <- PE$std.nox
        }

        attributes(PE)$information <- "MCMC"
        attributes(PE)$se <- "MCMC"
        ## 95% HPD; FIXME display blanks for equality-constrained parameters
        ##                (like Std.Err column)

        ## remove rho parameters from PE
        rhos <- grep("rho", object@ParTable$jlabel[object@ParTable$op != "=="])
        if(length(rhos) > 0) PE <- PE[-rhos,]
        ## move rho priors to covariance rows
        rhos <- grep("rho", object@external$runjags$origpt$jlabel)
        covrhos <- grep("@rho", object@external$runjags$origpt$plabel)
        object@ParTable$prior[covrhos] <- object@external$runjags$origpt$prior[rhos]
        ## remove equality constraints from ParTable (rhos removed in blavaan())
        eqc <- which(object@ParTable$op == "==")
        if(length(eqc) > 0){
            newpt <- lapply(object@ParTable, function(x) x[-eqc])
        } else {
            newpt <- object@ParTable
        }

        ## match jags names to partable, then partable to PE
        ptentry <- match(rownames(object@external$runjags$HPD), newpt$jlabel) #object@ParTable$jlabel)
        pte2 <- ptentry[!is.na(ptentry)]
        peentry <- match(with(newpt, paste(lhs[pte2], op[pte2], rhs[pte2], group[pte2], sep="")),
                         paste(PE$lhs, PE$op, PE$rhs, PE$group, sep=""))
        PE$ci.lower[peentry] <- object@external$runjags$HPD[!is.na(ptentry),'Lower95']
        PE$ci.upper[peentry] <- object@external$runjags$HPD[!is.na(ptentry),'Upper95']

        ## NB This is done so that we can remove fixed parameter hpd intervals without
        ##    making changes to lavaan's print.lavaan.parameterEstimates(). But maybe
        ##    this should actually go in the lavaan function.
        char.format <- paste("%", max(8, nd + 5), "s", sep="")
        PE$ci.lower <- round(PE$ci.lower, nd)
        PE$ci.lower[PE$ci.lower == PE$est] <- ""
        PE$ci.lower <- sprintf(char.format, PE$ci.lower)
        PE$ci.upper <- round(PE$ci.upper, nd)
        PE$ci.upper[PE$ci.upper == PE$est] <- ""
        PE$ci.upper <- sprintf(char.format, PE$ci.upper)

        ## FIXME defined parameters never get psrf + others;
        ## see line 200 of lav_print.R
        if(psrf & class(object@external$runjags$psrf) != "character"){
          PE$psrf <- rep(NA, nrow(PE))
          PE$psrf[peentry] <- object@external$runjags$psrf$psrf[!is.na(ptentry),'Point est.']
        }
        if(neff){
          PE$neff <- rep(NA, nrow(PE))
          PE$neff[peentry] <- object@external$runjags$summaries[!is.na(ptentry),'SSeff']
        }
        if(priors){
          PE$prior <- rep(NA, nrow(PE))
          PE$prior[peentry] <- newpt$prior[pte2]
          PE$prior[is.na(PE$prior)] <- ""
        }
        if(postmedian){
          PE$Post.Med <- rep(NA, nrow(PE))
          PE$Post.Med[peentry] <- object@external$runjags$summaries[!is.na(ptentry),'Median']
        }
        if(postmode){
          PE$Post.Mode <- rep(NA, nrow(PE))
          PE$Post.Mode[peentry] <- object@external$runjags$summaries[!is.na(ptentry),'Mode']
          if(all(is.na(PE$Post.Mode))) warning("blavaan WARNING: Posterior modes require installation of the modeest package.")
        }
        if(bf){
          ## we don't know whether priors=TRUE:
          PE2 <- PE
          if(!("prior" %in% names(PE))){
            tmppri <- rep("", nrow(PE))
            tmppri[peentry] <- newpt$prior[pte2]
            PE2$prior <- tmppri
          }
          PE$logBF <- round(SDBF(PE2), nd)
          PE$logBF[is.na(PE$logBF)] <- ""
          PE$logBF <- sprintf(char.format, PE$logBF)
        }
        ## alternative names because this is not ML
        penames <- names(PE)
        ## FIXME we need an est column for print.lavaan if we have constraints
        names(PE)[penames == "est"] <- "Post.Mean"
        names(PE)[penames == "se"] <- "Post.SD"
        names(PE)[penames == "ci.lower"] <- "HPD.025"
        names(PE)[penames == "ci.upper"] <- "HPD.975"
        names(PE)[penames == "psrf"] <- "PSRF"
        print(PE, nd = nd)

    } # parameter estimates

    # modification indices?
    if(modindices) {
        cat("Modification Indices:\n\n")
        object@Options$estimator <- "ML"
        object@Fit@test[[2]] <- NULL
        print( modificationIndices(object, standardized=TRUE) )
    }

})


# setMethod("vcov", "blavaan",
# function(object, labels=TRUE) {
# 
#
#    # check for convergence first!
#    if(object@Fit@npar > 0L && !object@Fit@converged)
#        warning("blavaan WARNING: chains may not have converged, proceed with caution.")
#    
#    VarCov <- object@external$runjags$vcov
#
#    labs <- lav_partable_labels(object@ParTable, type="free")
#    
#    if(labels) rownames(VarCov) <- colnames(VarCov) <- labs
#
#    #if(!redundant) VarCov <- VarCov[!duplicated(labs), !duplicated(labs)]
#
#    class(VarCov) <- c("lavaan.matrix.symmetric", "matrix")
#   
#    VarCov
#})


# logLik (so that we can use the default AIC/BIC functions from stats4(
# setMethod("logLik", "blavaan",
# function(object, ...) {
#    if(object@Options$estimator != "ML") {
#        warning("lavaan WARNING: logLik only available if estimator is ML")
#    }
#    if(object@Fit@npar > 0L && !object@Fit@converged) {
#        warning("lavaan WARNING: model did not converge")
#    }
#    
#    logl.df <- fitMeasures(object, c("logl", "npar", "ntotal"))
#    names(logl.df) <- NULL
#    logl <- logl.df[1]
#    attr(logl, "df") <- logl.df[2]    ### note: must be npar, not df!!
#    attr(logl, "nobs") <- logl.df[3]
#    class(logl) <- "logLik"
#    logl
#})

# nobs
## if(!exists("nobs", envir=asNamespace("stats4"))) {
##     setGeneric("nobs", function(object, ...) standardGeneric("nobs"))
## }
## setMethod("nobs", signature(object = "blavaan"),
## function(object, ...) {
##     object@SampleStats@ntotal
## })

# see: src/library/stats/R/update.R
## setMethod("update", signature(object = "lavaan"),
## function(object, model, ..., evaluate = TRUE) {

##     call <- object@call
##     if(is.null(call))
##         stop("need an object with call slot")

##     extras <- match.call(expand.dots = FALSE)$...

##     if(!missing(model))
##         #call$formula <- update.formula(formula(object), formula.)
##         call$model <- model

##     if(length(extras) > 0) {
##         existing <- !is.na(match(names(extras), names(call)))
##         for(a in names(extras)[existing]) call[[a]] <- extras[[a]]
##         if(any(!existing)) {
##             call <- c(as.list(call), extras[!existing])
##             call <- as.call(call)
##         }
##     }
##     if (evaluate) {
##         eval(call, parent.frame())
##     }
##     else call
## })

plot.blavaan <- function(x, pars, plot.type="trace", ...){
    # NB: arguments go to plot.runjags()
    parnames <- rownames(x@external$runjags$summaries)[pars]
    plot(x@external$runjags, plot.type=plot.type, vars=parnames, ...)
}
    
#setMethod("anova", signature(object = "blavaan"),
#function(object, ...) {
#
#    # NOTE: if we add additional arguments, it is not the same generic
#    # anova() function anymore, and match.call will be screwed up
#
#    # NOTE: we need to extract the names of the models from match.call here,
#    #       otherwise, we loose them in the call stack
#
#    mcall <- match.call(expand.dots = TRUE)
#    dots <- list(...)
#
#    # catch SB.classic and SB.H0
#    SB.classic <- TRUE; SB.H0 <- FALSE
#
#    arg.names <- names(dots)
#    arg.idx <- which(nchar(arg.names) > 0L)
#    if(length(arg.idx) > 0L) {
#        if(!is.null(dots$SB.classic))
#            SB.classic <- dots$SB.classic
#        if(!is.null(dots$SB.H0))
#            SB.H0 <- dots$SB.H0           
#        dots <- dots[-arg.idx]
#    }
#
#    modp <- if(length(dots))
#        sapply(dots, is, "lavaan") else logical(0)
#    mods <- c(list(object), dots[modp])
#    NAMES <- sapply(as.list(mcall)[c(FALSE, TRUE, modp)], deparse)
#
#    # use do.call to handle changed dots
#    ans <- do.call("lavTestLRT", c(list(object = object, 
#                   SB.classic = SB.classic, SB.H0 = SB.H0, 
#                   model.names = NAMES), dots))
#
#    ans
#})


BF <- function(object1, object2, ...) {

    # NB Assumes test slot instead of Fit@test slot
    # Bayes factor approximation based on marginal log-likelihoods
    bf <- object1@test[[1]]$stat - object2@test[[1]]$stat

    cat("Laplace approximation to the log-Bayes factor (experimental):\n",
        sprintf("%8.3f", bf), "\n\n")

    res <- c(bf, object1@test[[1]]$stat, object2@test[[1]]$stat)
    names(res) <- c("bf", "mll1", "mll2")

    invisible(res)
}


SDBF <- function(PE) {
  tmprow <- which(PE$op %in% c("~", "=~"))

  postdens <- dnorm(0, mean=PE$est[tmprow], sd=PE$se[tmprow], log=TRUE)

  pricom <- jagsdist2r(PE$prior[tmprow])
  pridens <- rep(NA, length(tmprow))
  for(i in 1:length(tmprow)){
      tmpdens <- try(eval_prior(pricom[[i]], 0, ""), silent=TRUE)
      
      if(!inherits(tmpdens, "try-error")) pridens[i] <- tmpdens
  }

  bf <- rep(NA, length(PE$op))
  bf[tmprow] <- pridens - postdens
  
  bf
}
