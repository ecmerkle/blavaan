#
# much of this comes from lav_object_methods.R

bl.short.summary <- function(object) {

    # catch FAKE run
    FAKE <- FALSE
    if(!is.null(object@Options$optim.method)) {
        if(tolower(object@Options$optim.method) == "none") {
            FAKE <- TRUE
        }
    }

    class(object) <- "lavaan"
    garb <- capture.output( tmp <- show(object) )
    tmp$test <- NULL
    cat("b")
    print(tmp)
    cat("\n")
  

    # Print margloglik + ppp
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
    bl.short.summary(object)

})


setMethod("summary", signature(object = "blavaan"),
function(object, header       = TRUE,
                 fit.measures = FALSE,
                 estimates    = TRUE,
                 ci           = TRUE,
                 standardized = FALSE,
                 rsquare      = FALSE,
                 std.nox      = FALSE, #TODO: remove deprecated argument in early 2025
                 psrf         = TRUE,
                 neff         = FALSE,
                 postmedian   = FALSE,
                 postmode     = FALSE,
                 priors       = TRUE,
                 bf           = FALSE,
                 nd = 3L) {

    #TODO: remove (deprecated):  if(std.nox) standardized <- TRUE

    # print the 'short' summary
    if(header) {
        bl.short.summary(object)
    }

    # only if requested, the fit measures
    if(fit.measures) {
        bopts <- blavInspect(object, "options")
        if(bopts$test == "none" & bopts$target != "stan") {
            warning("lavaan WARNING: fit measures not available if test = \"none\"", call. = FALSE)
        } else {
            #print.fit.measures( fitMeasures(object, fit.measures="default") )
            # TODO: define proper print function for custom fit measures
        }
    }


  if(estimates) {
        jagtarget <- lavInspect(object, "options")$target == "jags"
        newpt <- object@ParTable
        if(!("group" %in% names(newpt))) newpt$group <- rep(1, length(newpt$lhs))
        if(!("level" %in% names(newpt))) newpt$level <- rep("within", length(newpt$lhs))
        newpt$group[newpt$group == 0] <- 1 # for defined parameters

        if(!jagtarget){
            rhorows <- which(newpt$mat == "rho" | newpt$mat == "lvrho")
            if(length(rhorows) > 0){
                newpt <- lapply(newpt, function(x) x[-rhorows])
                object@ParTable <- lapply(object@ParTable, function(x) x[-rhorows])
            }
        }
        PE <- parameterEstimates(object, se = TRUE, zstat = FALSE,
                                 ci = TRUE,
                                 standardized = standardized,
                                 rsquare = rsquare,
                                 remove.eq = FALSE, remove.system.eq = TRUE,
                                 remove.ineq = FALSE, remove.def = FALSE,
                                 header = TRUE, output = "text")
        
        if(!("group" %in% names(PE))) PE$group <- 1
        if(!("level" %in% names(PE))) PE$level <- "within"
        if(!("psrf" %in% names(PE))) PE$psrf <- NA
        #TODO: remove deprecated argument in early 2025
        # if(standardized && std.nox) {
        #     PE$std.all <- PE$std.nox
        # }
        PE$group[PE$group == 0] <- 1

        if("target" %in% names(object@call)){
            if(object@call$target == "vb"){
                attributes(PE)$information <- "VB"
                attributes(PE)$se <- "VB"
            }
        } else {
            attributes(PE)$information <- "MCMC"
            attributes(PE)$se <- "MCMC"
        }
        ## 95% HPD; FIXME display blanks for equality-constrained parameters
        ##                (like Std.Err column)

        ## match jags names to partable, then partable to PE
        if(jagtarget){
            pte2 <- which(!is.na(newpt$jagpnum))
        } else {
            pte2 <- which(newpt$free > 0)
        }
        peentry <- match(with(newpt, paste(lhs[pte2], op[pte2], rhs[pte2], group[pte2], level[pte2], sep="")),
                         paste(PE$lhs, PE$op, PE$rhs, PE$group, PE$level, sep=""))
        if(!("ci.lower" %in% names(PE))) {
          PE$ci.lower <- PE$ci.upper <- NA
        }
        if(jagtarget){
          if('Lower95' %in% colnames(object@external$mcmcout$HPD)){
            PE$ci.lower[peentry] <- object@external$mcmcout$HPD[newpt$jagpnum[pte2],'Lower95']
            PE$ci.upper[peentry] <- object@external$mcmcout$HPD[newpt$jagpnum[pte2],'Upper95']
          } else {
            PE$ci.lower[peentry] <- rep(NA, length(peentry))
            PE$ci.upper[peentry] <- rep(NA, length(peentry))
          }
        } else {
            parsumm <- rstan::summary(object@external$mcmcout)
            if('2.5%' %in% colnames(parsumm[[1]]) & '97.5%' %in% colnames(parsumm[[1]])){
                PE$ci.lower[peentry] <- parsumm$summary[newpt$stansumnum[pte2],'2.5%']
                PE$ci.upper[peentry] <- parsumm$summary[newpt$stansumnum[pte2],'97.5%']
            } else {
                PE$ci.lower[peentry] <- rep(NA, length(peentry))
                PE$ci.upper[peentry] <- rep(NA, length(peentry))
            }
        }

        ## NB This is done so that we can remove fixed parameter hpd intervals without
        ##    making changes to lavaan's print.lavaan.parameterEstimates(). But maybe
        ##    this should actually go in the lavaan function.
        char.format <- paste("%", max(8, nd + 5), "s", sep="")
        PE$ci.lower <- formatC(as.numeric(PE$ci.lower), digits = nd, format = "f")
        PE$ci.lower[PE$ci.lower == formatC(PE$est, digits = nd, format = "f")] <- ""
        PE$ci.lower <- sprintf(char.format, PE$ci.lower)
        PE$ci.upper <- formatC(as.numeric(PE$ci.upper), digits = nd, format = "f")
        PE$ci.upper[PE$ci.upper == formatC(PE$est, digits = nd, format = "f")] <- ""
        PE$ci.upper <- sprintf(char.format, PE$ci.upper)

        ## FIXME defined parameters never get psrf + others;
        ## see line 200 of lav_print.R
        if(psrf){
          PE$psrf[peentry] <- formatC(as.numeric(newpt$psrf[pte2]), digits = nd, format = "f")
          PE$psrf[is.na(PE$psrf)] <- ""
          PE$psrf <- sprintf(char.format, PE$psrf)
        }
        if(neff){
          PE$neff <- rep(NA, nrow(PE))
          if(jagtarget){
            PE$neff[peentry] <- object@external$mcmcout$summaries[newpt$jagpnum[pte2],'SSeff']
          } else {
            PE$neff[peentry] <- parsumm$summary[newpt$stansumnum[pte2],'n_eff']
          }
        }
        if(priors){
          PE$prior <- rep(NA, nrow(PE))
          PE$prior[peentry] <- newpt$prior[pte2]
          PE$prior[is.na(PE$prior)] <- ""
        }
        if(postmedian){
          PE$Post.Med <- rep(NA, nrow(PE))
          if(jagtarget){
            PE$Post.Med[peentry] <- object@external$mcmcout$summaries[newpt$jagpnum[pte2],'Median']
          } else {
            PE$Post.Med[peentry] <- parsumm$summary[newpt$stansumnum[pte2],'50%']
          }
        }
        if(postmode){
          PE$Post.Mode <- rep(NA, nrow(PE))
          if(jagtarget){
            PE$Post.Mode[peentry] <- object@external$mcmcout$summaries[newpt$jagpnum[pte2],'Mode']
            if(all(is.na(PE$Post.Mode))) warning("blavaan WARNING: Posterior modes require installation of the modeest package.", call. = FALSE)
          } else {
            PE$Post.Mode[peentry] <- blavInspect(object, "postmode")
          }
        }
        if(bf){
          ## we don't know whether priors=TRUE:
          PE2 <- PE
          if(!("prior" %in% names(PE))){
            tmppri <- rep("", nrow(PE))
            tmppri[peentry] <- newpt$prior[pte2]
            PE2$prior <- tmppri
          }
          PE$logBF <- formatC(SDBF(PE2), digits = nd, format = "f")
          PE$logBF[is.na(PE$logBF)] <- ""
          PE$logBF <- sprintf(char.format, PE$logBF)
        }
        
        ## alternative names because this is not ML
        penames <- names(PE)
        ## This could be called "Post.Mean" except constraints
        ## require "est"
        #names(PE)[penames == "est"] <- "Post.Mean"
        #PE$est <- PE$Post.Mean
        if(!('prisamp' %in% names(blavInspect(object, 'options')))){
          ## backwards compatibility before we had prisamp
          names(PE)[penames == "se"] <- "Post.SD"
        } else {
          if(blavInspect(object, 'options')$prisamp){
            names(PE)[penames == "se"] <- "Pri.SD"
          } else {
            names(PE)[penames == "se"] <- "Post.SD"
          }
        }
        names(PE)[penames == "ci.lower"] <- "pi.lower"
        names(PE)[penames == "ci.upper"] <- "pi.upper"
        names(PE)[penames == "psrf"] <- "Rhat"

        print(PE, nd = nd)
    } # parameter estimates
}
)

# NB not absolutely necessary, except for
# bug in lavaan 0.5-23.1097
setMethod("coef", "blavaan",
  function(object, type="free", labels=TRUE) {
    class(object) <- "lavaan"
    callNextMethod(object, type, labels)
  })


plot.blavaan <- function(x, pars=NULL, plot.type="trace", showplot=TRUE, ...){
    # NB: arguments now go to bayesplot functions
    if(length(pars) == 0L){
        pars <- x@ParTable$free
        pars <- pars[pars > 0 & !is.na(pars)]
    }

    if(x@Options$target != "stan"){
        samps <- as.array(blavInspect(x, 'mcmc', add.labels = FALSE), drop = FALSE)
        parnames <- x@ParTable$pxnames[match(pars, x@ParTable$free)]
        samps <- samps[, match(parnames, colnames(samps)), , drop = FALSE]
        ## samps dims must be "iteration, chain, parameter"
        samps <- aperm(samps, c(1, 3, 2))
    } else {
        samps <- as.array(x@external$mcmcout)
        parnums <- x@ParTable$stanpnum[match(pars, x@ParTable$free)]
        samps <- samps[, , parnums, drop = FALSE]
    }
    if(blavInspect(x, 'ngroups') == 1L){
        dimnames(samps)[[3]] <- with(x@ParTable, paste0(lhs,op,rhs)[match(pars, free)])
    } else {
        dimnames(samps)[[3]] <- with(x@ParTable, paste0(lhs,op,rhs,".g",group)[match(pars, free)])
    }

    plfun <- get(paste0("mcmc_", plot.type), asNamespace("bayesplot"))

    pl <- do.call(plfun, c(list(x = samps), list(...)))

    if(showplot) plot(pl)

    invisible(pl)
}

#setMethod("plot", "blavaan", plot.blavaan)

## function/environment for y-axis label of runjags plots
labelfun <- function(var){
  unlist(axis_env$trans[axis_env$trans[,1] == var, 2]) #axis_env$trans[panel.number(),2]
}
  
axis_env <- new.env(parent = emptyenv())


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


## obtain full posterior distribution of standardized parameters
standardizedPosterior <- standardizedposterior <- function(object, ...) {

  dots <- list(...)

  allowargs <- c('type', 'cov.std', 'remove.eq', 'remove.ineq', 'remove.def')
  
  if(!all(names(dots) %in% allowargs)){
    stop(paste0("blavaan ERROR: arguments must be in ",
                paste(allowargs, collapse=" ")))
  }
  
  ## posterior draws in matrix form
  draws <- make_mcmc(object@external$mcmcout)
  draws <- do.call("rbind", draws)

  ## copy of fit, for manipulation
  tf <- object

  ## put a posterior draw in the tf object + compute standardized estimates
  tmp <- fill_params(draws[1,], object@Model, object@ParTable)
  tf@Model <- tmp
  tf@ParTable$est[tf@ParTable$free > 0] <- lav_model_get_parameters(tmp)
  tmp2 <- do.call("standardizedSolution", c(list(object=tf), dots))

  ## use tmp2 object to figure out the right number of columns, then loop
  fullres <- matrix(NA, nrow(draws), nrow(tmp2))
  colnames(fullres) <- with(tmp2, paste0(lhs, op, rhs))
  if("group" %in% colnames(tmp2)) colnames(fullres) <- paste0(colnames(fullres), '.g', tmp2$group)
  fullres[1,] <- tmp2[, 'est.std']

  for(i in 2:nrow(draws)){
    tmp <- fill_params(draws[i,], object@Model, object@ParTable)
    tf@Model <- tmp
    tf@ParTable$est[tf@ParTable$free > 0] <- lav_model_get_parameters(tmp)

    ## compute standardized estimates
    fullres[i,] <- do.call("standardizedSolution", c(list(object=tf), dots))[, 'est.std']
  }

  fullres
}
