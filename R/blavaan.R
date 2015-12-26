# blavaan2
#
# different strategy: call lavaan() first, without fitting

blavaan <- function(...,  # default lavaan arguments
 
                    # bayes-specific stuff
                    ov.cp              = "srs",
                    lv.cp              = "srs",
                    dp                 = dpriors(),
                    n.chains           = 3,
                    burnin             ,
                    sample             ,
                    adapt              ,
                    jagfile            = FALSE,
                    jagextra           = list(),
                    inits              = "prior",
                    jagcontrol         = list()
                   )
{
    ## start timer
    start.time0 <- proc.time()[3]

    # store original call
    mc  <- match.call()

    # catch dot dot dot
    dotdotdot <- list(...); dotNames <- names(dotdotdot)

    # which arguments do we override?
    lavArgsOverride <- c("meanstructure", "missing", "estimator")
    # always warn?
    warn.idx <- which(lavArgsOverride %in% dotNames)
    if(length(warn.idx) > 0L) {
        warning("blavaan WARNING: the following arguments have no effect:\n",
                "                   ", 
                paste(lavArgsOverride[warn.idx], collapse = " "))
    }
    # if do.fit supplied, save it for jags stuff
    jag.do.fit <- TRUE
    if("do.fit" %in% dotNames) jag.do.fit <- dotdotdot$do.fit
    dotdotdot$do.fit <- FALSE        # implies se="none" and test="none"
    dotdotdot$meanstructure <- TRUE
    dotdotdot$missing <- "direct"   # direct/ml creates error? (bug in lavaan?)
    dotdotdot$estimator <- "default" # until 'Bayes' is accepted by lavaan()

    # jags args
    if("debug" %in% dotNames) {
        if(dotdotdot$debug)  {
            ## short burnin/sample
            mc$burnin <- 1000
            mc$sample <- 1000
        }
    }
    jarg <- c("n.chains", "burnin", "sample", "adapt")
    mcj <- match(jarg, names(mc), 0L)
    if(any(mcj > 0)){
        mfj <- as.list(mc[mcj])
        if("sample" %in% names(mc)){
            if(mc$sample*n.chains/5 < 1000) warning("blavaan WARNING: small sample drawn, proceed with caution.\n")
        }
    } else {
        mfj <- list()
    }

    # which argument do we remove/ignore?
    lavArgsRemove <- c("likelihood", "information", "se", "bootstrap",
                       "wls.v", "nacov", "zero.add", "zero.keep.margins",
                       "zero.cell.warn")
    dotdotdot[lavArgsRemove] <- NULL
    warn.idx <- which(lavArgsRemove %in% dotNames)
    if(length(warn.idx) > 0L) {
        warning("blavaan WARNING: the following arguments are ignored:\n",
                "                   ",
                paste(lavArgsRemove[warn.idx], collapse = " "))
    }

    # call lavaan
    LAV <- do.call("lavaan", dotdotdot)

    # check for ordered data
    if(lavInspect(LAV, "categorical")) {
        stop("blavaan ERROR: ordered data is not supported yet.")
    }

    if(LAV@Data@data.type == "moment") {
        stop("blavaan ERROR: full data are required. consider using kd() from package semTools.")
    }

    ineq <- which(LAV@ParTable$op %in% c("<",">"))
    if(length(ineq) > 0) {
        LAV@ParTable <- lapply(LAV@ParTable, function(x) x[-ineq])
        if(jagfile==FALSE) jagfile <- TRUE
        warning("blavaan WARNING: blavaan does not handle inequality constraints.\ntry modifying the exported JAGS code.")
    }

    eqs <- (LAV@Model@ceq.JAC == -1 | LAV@Model@ceq.JAC == 1)
    if(length(dim(eqs)) > 0) {
        compeq <- which(LAV@Model@ceq.rhs != 0 |
                        rowSums(LAV@Model@ceq.JAC != 0) != 2 |
                        rowSums(eqs) != 2)
        if(length(compeq) > 0) {
            eqpars <- which(LAV@ParTable$op == "==")
            LAV@ParTable <- lapply(LAV@ParTable, function(x) x[-eqpars[compeq]])
            warning("blavaan WARNING: blavaan does not currently handle complex equality constraints.\ntry modifying the exported JAGS code.")
        }
    }
    
    # cannot currently use wishart prior with std.lv=TRUE
    if(LAV@Options$auto.cov.lv.x & LAV@Options$std.lv){
        #warning("blavaan WARNING: cannot use Wishart prior with std.lv=TRUE. Reverting to 'srs' priors.")
        LAV@Options$auto.cov.lv.x <- FALSE
    }
    # Check whether there are user-specified priors or equality
    # constraints on lv.x or ov.x. If so, set auto.cov.lv.x = FALSE.
    lv.x <- LAV@pta$vnames$lv.x[[1]]
    ## catch some regressions without fixed x:
    ov.noy <- LAV@pta$vnames$ov.nox[[1]]
    ov.noy <- ov.noy[!(ov.noy %in% LAV@pta$vnames$ov.y)]
    prispec <- "prior" %in% names(LAV@ParTable)
    ndpriors <- rep(FALSE, length(LAV@ParTable$lhs))
    if(prispec) ndpriors <- LAV@ParTable$prior != ""
    con.cov <- any(LAV@ParTable$lhs %in% c(lv.x, ov.noy) &
                   LAV@ParTable$op == "~~" &
                   (LAV@ParTable$free == 0 |
                    ndpriors))
    if(con.cov) LAV@Options$auto.cov.lv.x <- FALSE

    # if std.lv, truncate the prior of each lv's first loading
    if(LAV@Options$std.lv){
        if(ov.cp == "fa" | lv.cp == "fa") stop("blavaan ERROR: ov.cp='fa' and lv.cp='fa' cannot be used with std.lv=TRUE.")
        if(!prispec){
            LAV@ParTable$prior <- rep("", length(LAV@ParTable$id))
        }
        ## first loading for each lv
        loadpt <- LAV@ParTable$op == "=~"
        lvs <- unique(LAV@ParTable$lhs[loadpt])
        fload <- sapply(lvs, function(x) which(LAV@ParTable$lhs[loadpt] == x)[1])

        for(i in 1:length(fload)){
            if(LAV@ParTable$prior[fload[i]] != ""){
                LAV@ParTable$prior[fload[i]] <- paste(LAV@ParTable$prior[fload[i]], "T(0,)", sep="")
            } else {
                LAV@ParTable$prior[fload[i]] <- paste(dp[["lambda"]], "T(0,)", sep="")
            }
        }
    }

    # if jagfile is a directory, vs logical
    if(class(jagfile)=="character"){
        jagdir <- jagfile
        jagfile <- TRUE
    }  else {
        jagdir <- "lavExport"
    }
  
    # extract slots from dummy lavaan object
    lavpartable    <- LAV@ParTable
    lavmodel       <- LAV@Model
    lavdata        <- LAV@Data
    lavoptions     <- LAV@Options
    lavsamplestats <- LAV@SampleStats
    lavcache       <- LAV@Cache
    timing         <- LAV@timing

    # change some 'default' @Options and add some
    lavoptions$estimator <- "Bayes"
    lavoptions$se        <- "standard"
    lavoptions$test <- "standard"
    if("test" %in% dotNames) {
        if(dotdotdot$test == "none") lavoptions$test <- "none"
    } else {
        # if missing data, posterior predictives are way slow
        if(any(is.na(unlist(LAV@Data@X)))) {
            cat("blavaan NOTE: Posterior predictives with missing data are currently very slow.\nConsider setting test=\"none\".\n\n")
        }
    }
    if(!jag.do.fit){
      lavoptions$test <- "none"
      lavoptions$se <- "none"
    }
    lavoptions$missing   <- "ml"
    lavoptions$ov.cp     <- ov.cp
    lavoptions$lv.cp     <- lv.cp

    verbose <- lavoptions$verbose

    # redo estimation + vcov + test
    # 6. estimate free parameters
    start.time <- proc.time()[3]
    x <- NULL
    if(lavmodel@nx.free > 0L) {
        ## convert partable to jags, then run
        jagtrans <- try(lav2jags(model = lavpartable, lavdata = lavdata, 
                                 ov.cp = ov.cp, lv.cp = lv.cp,
                                 lv.x.wish = lavoptions$auto.cov.lv.x,
                                 dp = dp, n.chains = n.chains, jagextra = jagextra, inits = inits),
                        silent = TRUE)

        if(!inherits(jagtrans, "try-error")){
            if(jagfile){
                dir.create(path=jagdir, showWarnings=FALSE)
                cat(jagtrans$model, file = paste(jagdir, "/sem.jag",
                                                 sep=""))
                save(jagtrans, file = paste(jagdir, "/semjags.rda",
                                            sep=""))
            }

            ## merge coefvec with lavpartable
            lavpartable <- mergejag(lavpartable, jagtrans$coefvec)

            ## add extras to monitor, if specified
            sampparms <- jagtrans$coefvec[,1]
            if("monitor" %in% names(jagextra)){
                sampparms <- c(sampparms, jagextra$monitor)
            }
            
            ## let jags set inits; helpful for debugging
            if(inits == "jags") jagtrans$inits <- NA
            rjarg <- with(jagtrans, list(model = paste(model),
                          monitor = sampparms, # "dic", "deviance"),
                          data = data, inits = inits))

            ## user-supplied jags params
            rjarg <- c(rjarg, mfj, jagcontrol)
            ## obtain posterior modes
            if(suppressMessages(requireNamespace("modeest"))) runjags.options(mode.continuous=TRUE)

            if(jag.do.fit){
              res <- try(do.call("run.jags", rjarg))
            } else {
              res <- NULL
            }

            if(inherits(res, "try-error")) {
                dir.create(path=jagdir, showWarnings=FALSE)
                cat(jagtrans$model, file = paste(jagdir, "/sem.jag",
                                                 sep=""))
                save(jagtrans, file = paste(jagdir, "/semjags.rda",
                                            sep=""))
                stop("blavaan ERROR: problem with jags estimation.  The jags model and data have been exported.")
            }
        } else {
            print(jagtrans)
            stop("blavaan ERROR: problem with translation from lavaan to jags.")
        }
        parests <- coeffun(lavpartable, res)
        x <- parests$x

        lavpartable <- parests$lavpartable
        if(jag.do.fit){
          lavmodel <- lav_model_set_parameters(lavmodel, x = x)

          # overwrite lavpartable$est to include def/con values
          lavpartable$est <- lav_model_get_parameters(lavmodel = lavmodel,
                                                      type = "user", extra = TRUE)

          attr(x, "iterations") <- res$sample
          attr(x, "converged") <- TRUE

          attr(x, "control") <- jagcontrol

          attr(x, "fx") <- get_ll(lavmodel = lavmodel, lavpartable = lavpartable,
                                  lavsamplestats = lavsamplestats, lavoptions = lavoptions,
                                  lavcache = lavcache, lavdata = lavdata)[1]

        } else {
          x <- numeric(0L)
          attr(x, "iterations") <- 0L; attr(x, "converged") <- FALSE
          attr(x, "control") <- jagcontrol
          attr(x, "fx") <- get_ll(lavmodel = lavmodel, 
                                  lavpartable = lavpartable,
                                  lavsamplestats = lavsamplestats,
                                  lavoptions = lavoptions,
                                  lavcache = lavcache,
                                  lavdata = lavdata)[1]
          
          lavpartable$est <- lavpartable$start
        }
    }

    ## parameter convergence + implied moments:
    rhorows <- grepl("rho", jagtrans$coefvec[,1])
    parrows <- 1:nrow(jagtrans$coefvec)
    lavimplied <- NULL
    ## compute/store some model-implied statistics
    lavimplied <- lav_model_implied(lavmodel)
    if(jag.do.fit){
      if(any(res$psrf$psrf[parrows[!rhorows],1] > 1.2)) attr(x, "converged") <- FALSE

      ## warn if psrf is large
      if(!attr(x, "converged") && lavoptions$warn) {
        warning("blavaan WARNING: at least one parameter has a psrf > 1.2.")
      }
    }

    ## fx is mean ll, where ll is marginal log-likelihood (integrate out lvs)
    if(lavoptions$test != "none") {
      cat("Computing posterior predictives...\n")
      res$samplls <- samp_lls(res, lavmodel, lavpartable, lavsamplestats,
                              lavoptions, lavcache, lavdata)
    } else {
      res$samplls <- NULL
    }
    
    ## put runjags output in new blavaan slot
    lavjags <- res
    timing$Estimate <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    ## 7. VCOV is now simple
    lavvcov <- list()
    VCOV <- NULL
    if(jag.do.fit){
      VCOV <- diag(parests$sd) %*% parests$vcorr %*% diag(parests$sd)
      lavjags <- c(lavjags, list(vcov = VCOV))

      # store vcov in new @vcov slot
      # strip all attributes but 'dim'
      tmp.attr <- attributes(VCOV)
      VCOV1 <- VCOV
      attributes(VCOV1) <- tmp.attr["dim"]
      lavvcov <- list(vcov = VCOV1)
    }

    timing$VCOV <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    ## 8. "test statistics": marginal log-likelihood, dic
    TEST <- list()
    if(lavoptions$test != "none") { # && attr(x, "converged")) {
        TEST <- blav_model_test(lavmodel            = lavmodel,
                                lavpartable         = lavpartable,
                                lavsamplestats      = lavsamplestats,
                                lavoptions          = lavoptions,
                                x                   = x,
                                VCOV                = VCOV,
                                lavdata             = lavdata,
                                lavcache            = lavcache,
                                lavjags             = lavjags,
                                jagextra            = jagextra)
        if(verbose) cat(" done.\n")
    }
    timing$TEST <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 9. collect information about model fit (S4)
    lavfit <- blav_model_fit(lavpartable = lavpartable,
                             lavmodel    = lavmodel,
                             lavjags     = lavjags,
                             x           = x,
                             VCOV        = VCOV,
                             TEST        = TEST)
    ## add SE and SD-Bayes factor to lavpartable
    ## (code around line 270 of blav_object_methods
    ##  can be removed with this addition near line 923
    ##  of lav_object_methods:
    ## if(object@Options$estimator == "Bayes") {
    ##    LIST$logBF <- PARTABLE$logBF
    ## }
    lavpartable$se <- lavfit@se[lavpartable$id]
    lavpartable$logBF <- SDBF(lavpartable)
    timing$total <- (proc.time()[3] - start.time0)

    ## remove rhos from partable + ses, so lavaan built-ins work
    rhos <- grep("rho", lavpartable$jlabel)
    lavjags <- c(lavjags, list(origpt=lavpartable))
    class(lavjags) <- "runjags"
    if(length(rhos) > 0){
        lavpartable <- lapply(lavpartable, function(x) x[-rhos])
        lavfit@se <- lavfit@se[-rhos]
    }

    # 10. construct lavaan object
    blavaan <- new("blavaan",
                   call         = mc,                  # match.call
                   timing       = timing,              # list
                   Options      = lavoptions,          # list
                   ParTable     = lavpartable,         # list
                   pta          = LAV@pta,             # list
                   Data         = lavdata,             # S4 class
                   SampleStats  = lavsamplestats,      # S4 class
                   Model        = lavmodel,            # S4 class
                   Cache        = lavcache,            # list
                   Fit          = lavfit,              # S4 class
                   boot         = list(),
                   optim        = list(),
                   implied      = lavimplied,          # list
                   vcov         = lavvcov,
                   external     = list(runjags = lavjags), # runjags
                   test         = TEST                 # copied for now
                  )

    # post-fitting check
    if(attr(x, "converged")) {
        lavInspect(blavaan, "post.check")
    }

    blavaan
}

## cfa + sem
bcfa <- bsem <- function(..., ov.cp = "srs", lv.cp = "srs", dp = dpriors(),
    n.chains = 3, burnin, sample, adapt,
    jagfile = FALSE, jagextra = list(), inits = "prior", jagcontrol = list()) {

    dotdotdot <- list(...)
    std.lv <- ifelse(any(names(dotdotdot) == "std.lv"), dotdotdot$std.lv, FALSE)

    mc <- match.call()  
    mc$model.type      = as.character( mc[[1L]] )
    if(length(mc$model.type) == 3L) mc$model.type <- mc$model.type[3L]
    mc$int.ov.free     = TRUE
    mc$int.lv.free     = FALSE
    mc$auto.fix.first  = !std.lv
    mc$auto.fix.single = TRUE
    mc$auto.var        = TRUE
    mc$auto.cov.lv.x   = TRUE
    mc$auto.cov.y      = TRUE
    mc$auto.th         = TRUE
    mc$auto.delta      = TRUE
    mc[[1L]] <- quote(blavaan)

    eval(mc, parent.frame())
}

# simple growth models
bgrowth <- function(..., ov.cp = "srs", lv.cp = "srs", dp = dpriors(),
    n.chains = 3, burnin, sample, adapt,
    jagfile = FALSE, jagextra = list(), inits = "prior", jagcontrol = list()) {

    dotdotdot <- list(...)
    std.lv <- ifelse(any(names(dotdotdot) == "std.lv"), dotdotdot$std.lv, FALSE)

    mc <- match.call()
    mc$model.type      = "growth"
    mc$int.ov.free     = FALSE
    mc$int.lv.free     = TRUE
    mc$auto.fix.first  = !std.lv
    mc$auto.fix.single = TRUE
    mc$auto.var        = TRUE
    mc$auto.cov.lv.x   = TRUE
    mc$auto.cov.y      = TRUE
    mc$auto.th         = TRUE
    mc$auto.delta      = TRUE
    mc[[1L]] <- quote(blavaan)

    eval(mc, parent.frame())
}

