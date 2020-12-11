blavaan <- function(...,  # default lavaan arguments
 
                    # bayes-specific stuff
                    cp                 = "srs",
                    dp                 = NULL,
                    n.chains           = 3,
                    burnin             ,
                    sample             ,
                    adapt              ,
                    mcmcfile            = FALSE,
                    mcmcextra           = list(),
                    inits              = "prior",
                    convergence        = "manual",
                    target             = "stan",
                    save.lvs           = FALSE,
                    wiggle             = NULL,
                    wiggle.sd          = 0.1,
                    prisamp            = FALSE,
                    jags.ic            = FALSE,
                    seed               = NULL,
                    bcontrol         = list()
                   )
{
    ## start timer
    start.time0 <- proc.time()[3]

    # store original call
    mc  <- match.call()
    # catch dot dot dot
    dotdotdot <- list(...); dotNames <- names(dotdotdot)

    # default priors
    if(length(dp) == 0) dp <- dpriors(target = target)

    # burnin/sample/adapt if not supplied (should only occur for direct
    # blavaan call
    sampargs <- c("burnin", "sample", "adapt")
    suppargs <- which(!(sampargs %in% names(mc)))
    if(length(suppargs) > 0){
        if(target == "jags"){
            defiters <- c(4000L, 10000L, 1000L)
        } else {
            defiters <- c(500L, 1000L, 1000L)
        }
        for(i in 1:length(suppargs)){
            assign(sampargs[suppargs[i]], defiters[suppargs[i]])
        }
    }

    # multilevel functionality not available
    if("cluster" %in% dotNames) stop("blavaan ERROR: two-level models are not yet available.")
    if("ordered" %in% dotNames) stop("blavaan ERROR: ordinal models are not yet available.")
  
    # prior predictives only for stan
    if(prisamp) {
      if(target != 'stan') stop("blavaan ERROR: prior predictives currently only work for target='stan'.")
      if(!('test' %in% dotdotdot)) dotdotdot$test <- 'none'
    }
  
    # wiggle sd
    if(any(wiggle.sd <= 0L)) stop("blavaan ERROR: wiggle.sd must be > 0.")

    if(length(wiggle) > 0 & target == 'stancond') stop(paste0("blavaan ERROR: wiggle is currently not available for target ", target))

    if(length(wiggle.sd) > 1){
      if(target != 'stan') stop("blavaan ERROR: parameter-specific wiggle.sd is only available for target='stan'")
      
      if(length(wiggle.sd) != length(wiggle)) stop("blavaan ERROR: length of wiggle.sd must match length of wiggle")
    }
      
    # no current functionality to generate initial values from approximately-equal prior
    if(length(wiggle) > 0 & target %in% c('stanclassic', 'jags') & inherits(inits, 'character')) inits <- "simple"
  
    # ensure rstan/runjags are here. if target is not installed but
    # the other is, then use the other instead.
    if(grepl("stan", target)){
      if(convergence == "auto") stop("blavaan ERROR: auto convergence is unavailable for Stan.")

      if(target == "stan" & length(mcmcextra) > 0) stop("blavaan ERROR: mcmcextra is not available for target='stan'.")
    } else if(target == "jags"){
      if(!pkgcheck("runjags")){
        ## go to rstan if they have it
        if(pkgcheck("rstan")){
          cat("blavaan NOTE: runjags not installed; using rstan instead.\n")
          target <- "stan"
          pkgload("rstan")
        } else {
          stop("blavaan ERROR: runjags package is not installed.")
        }
      } else {
        pkgload("runjags")
      }
    }

    # if seed supplied, check that there is one per chain
    seedlen <- length(seed)
    if(seedlen > 0 & target == "jags" & seedlen != n.chains){
      stop("blavaan ERROR: for JAGS, number of seeds must equal n.chains.")
    }
  
    # capture data augmentation/full information options
    blavmis <- "da"
    if("missing" %in% dotNames) {
      if(dotdotdot$missing %in% c("da", "fi")) {
        blavmis <- dotdotdot$missing
        misloc <- which(dotNames == "missing")
        dotdotdot <- dotdotdot[-misloc]; dotNames <- dotNames[-misloc]
      }
    }

    # covariance priors are now all srs or fa
    cplocs <- match(c("ov.cp", "lv.cp"), dotNames, nomatch = 0)
    if(any(cplocs > 0)){
      cat("blavaan NOTE: Arguments ov.cp and lv.cp are deprecated. Use cp instead. \n")
      if(any(dotdotdot[cplocs] == "fa")){
        cp <- "fa"
        cat("Using 'fa' approach for all covariances.\n")
      } else {
        cat("\n")
      }
      if(!all(dotdotdot[cplocs] %in% c("srs", "fa"))){
        stop("blavaan ERROR: unknown argument to cp.")
      }
      dotdotdot <- dotdotdot[-cplocs]; dotNames <- dotNames[-cplocs]
    }
    ## cannot use lavaan inits with fa priors; FIXME?
    if(cp == "fa" & inherits(inits, 'character')){
      if(inits[1] %in% c("simple", "default")) inits <- "jags"
    }
    if(cp == "fa" & grepl("stan", target)){
      cat("blavaan NOTE: fa priors are not available with stan. srs priors will be used. \n")
    }

    # 'jag' arguments are now mcmcfile, mcmcextra, bcontrol
    jagargs <- c("jagfile", "jagextra")
    barg <- "jagcontrol"
    jaglocs <- match(jagargs, dotNames, nomatch = 0)
    blocs <- match(barg, dotNames, nomatch = 0)
    if(any(jaglocs > 0)){      
      cat(paste0("blavaan NOTE: the following argument(s) are deprecated: ",
                 paste(jagargs[jaglocs > 0], collapse=" "),
                 ".\n        the argument(s) now start with 'mcmc' instead of 'jag'. \n"))
      newargs <- gsub("jag", "mcmc", jagargs[jaglocs > 0])
      for(i in 1:length(newargs)){
        assign(newargs[i], dotdotdot[[jaglocs[jaglocs > 0][i]]])
      }
      dotdotdot <- dotdotdot[-jaglocs]; dotNames <- dotNames[-jaglocs]
    }
    if(any(blocs > 0)){      
      cat(paste0("blavaan NOTE: the following argument is deprecated: ",
                 paste(barg[blocs > 0], collapse=" "),
                 ".\n        the argument now starts with 'b' instead of 'jag'. \n"))
      newargs <- gsub("jag", "b", barg[blocs > 0])
      for(i in 1:length(newargs)){
        assign(newargs[i], dotdotdot[[blocs[blocs > 0][i]]])
      }
      dotdotdot <- dotdotdot[-blocs]; dotNames <- dotNames[-blocs]
    }
  
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
    if("do.fit" %in% dotNames){
        jag.do.fit <- dotdotdot$do.fit
        if(!jag.do.fit){
            burnin <- 0
            sample <- 0
        }
    }
    if("warn" %in% dotNames){
        origwarn <- dotdotdot$warn
    } else {
        origwarn <- TRUE
    }
    dotdotdot$do.fit <- TRUE
    dotdotdot$se <- "none"; dotdotdot$test <- "none"
    # run for 1 iteration to obtain info about equality constraints, for npar
    dotdotdot$control <- list(iter.max = 1); dotdotdot$warn <- FALSE
    dotdotdot$meanstructure <- TRUE
    dotdotdot$missing <- "direct"   # direct/ml creates error? (bug in lavaan?)
    if("ordered" %in% dotNames |
       any(apply(dotdotdot$data, 2, function(x) inherits(x, "ordered")))){
      dotdotdot$missing <- "pairwise" # needed to get missing patterns
    }
    dotdotdot$estimator <- "default"

    # jags args
    if("debug" %in% dotNames) {
        if(dotdotdot$debug)  {
            ## short burnin/sample
            burnin <- 100
            sample <- 100
        }
    }

    mfj <- list(burnin = burnin, sample = sample, adapt = adapt)

    if(target == "stancond") cat("\nblavaan NOTE: target='stancond' is experimental and may not be fully functional\n")
  
    if(convergence == "auto"){
        names(mfj) <- c("startburnin", "startsample", "adapt")
    }
    if(grepl("stan", target)){
        names(mfj) <- c("warmup", "iter", "adapt")
        ## stan iter argument includes warmup:
        mfj$iter <- mfj$warmup + mfj$iter
        mfj <- mfj[-which(names(mfj) == "adapt")]
    }

    if(target == "jags"){
        mfj <- c(mfj, list(n.chains = n.chains))
    } else {
        mfj <- c(mfj, list(chains = n.chains))
    }
    if(convergence == "auto"){
        if(!("startsample" %in% names(mfj))){
            ## bump down default
            mfj$startsample <- 4000
        } else {
            if(mfj$startsample < 4000){
                cat("blavaan NOTE: starting sample was increased to 4000 for auto-convergence\n")
                mfj$startsample <- 4000 # needed for runjags
            }
        }
        sample <- mfj$startsample
        if(!("startburnin" %in% names(mfj))){
            mfj$startburnin <- 4000
            burnin <- 4000
        } else {
            burnin <- mfj$startburnin
        }
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
                paste(lavArgsRemove[warn.idx], collapse = " "), "\n")
    }

    # call lavaan
    mcdebug <- FALSE
    if("debug" %in% dotNames){
      ## only debug mcmc stuff
      mcdebug <- dotdotdot$debug
      dotdotdot <- dotdotdot[-which(dotNames == "debug")]
    }
    LAV <- do.call("lavaan", dotdotdot)

    if(LAV@Data@data.type == "moment") {
        stop("blavaan ERROR: full data are required. consider using kd() from package semTools.")
    }

    # ordinal functionality not available
    if(lavInspect(LAV, 'categorical')) {
        stop("blavaan ERROR: models with ordered variables are not yet available.")
    }

    # save.lvs in a model with no lvs
    if(save.lvs){
        clv <- lavInspect(LAV, 'cov.lv')
        if(is.list(clv)) clv <- clv[[1]]
        if(nrow(clv) == 0) warning("blavaan WARNING: save.lvs=TRUE, but there are no lvs in the model.")
    }
        
    # turn warnings back on by default
    LAV@Options$warn <- origwarn

    # put original do.fit back
    LAV@Options$do.fit <- jag.do.fit
  
    # check for conflicting mv names
    namecheck(LAV@Data@ov.names[[1]])

    # deal with ordinal data, turn some options off
    ordmod <- LAV@Options$categorical
    if(ordmod){
      ## this picks up variables of class ordered that were not
      ## explicitly specified via ordered argument
      dotdotdot$ordered <- LAV@Data@ordered
      if(blavmis == "fi"){
        stop("blavaan ERROR: missing='fi' cannot be used with ordinal data.")
      }
      dotdotdot$test <- "none"
      dotNames <- names(dotdotdot)
    }  
  
    ineq <- which(LAV@ParTable$op %in% c("<",">"))
    if(length(ineq) > 0) {
        LAV@ParTable <- lapply(LAV@ParTable, function(x) x[-ineq])
        if(inherits(mcmcfile, "logical")) mcmcfile <- TRUE
        warning("blavaan WARNING: blavaan does not currently handle inequality constraints.\n try modifying the exported MCMC syntax.")
    }
    eqs <- which(LAV@ParTable$op == "==")
    if(length(eqs) > 0) {
        lhsvars <- rep(NA, length(eqs))
        for(i in 1:length(eqs)){
            lhsvars[i] <- length(all.vars(parse(file="", text=LAV@ParTable$lhs[eqs[i]])))
        }
        if(any(lhsvars > 1)) {
            stop("blavaan ERROR: blavaan does not handle equality constraints with more than 1 variable on the lhs.\n try modifying the constraints.")
        }
    }

    prispec <- "prior" %in% names(LAV@ParTable)
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
    ndpriors <- rep(FALSE, length(LAV@ParTable$lhs))
    if(prispec) ndpriors <- LAV@ParTable$prior != ""
    cov.pars <- (LAV@ParTable$lhs %in% c(lv.x, ov.noy)) & LAV@ParTable$op == "~~"
    con.cov <- any((cov.pars & (LAV@ParTable$free == 0 | ndpriors)) |
                   ((LAV@ParTable$lhs %in% LAV@ParTable$plabel[cov.pars] |
                     LAV@ParTable$rhs %in% LAV@ParTable$plabel[cov.pars]) &
                     LAV@ParTable$op == "=="))
    if(con.cov) LAV@Options$auto.cov.lv.x <- FALSE
    
    # if std.lv, truncate the prior of each lv's first loading
    if(LAV@Options$std.lv){
        if(cp == "fa") stop("blavaan ERROR: 'fa' prior strategy cannot be used with std.lv=TRUE.")
        if(!prispec){
            LAV@ParTable$prior <- rep("", length(LAV@ParTable$id))
        }
        loadpt <- LAV@ParTable$op == "=~"
        lvs <- unique(LAV@ParTable$lhs[loadpt])
        fload <- NULL
        if(length(lvs) > 0){
            for(i in 1:length(lvs)){
                for(k in 1:max(LAV@ParTable$group)){
                    fload <- c(fload, which(LAV@ParTable$lhs == lvs[i] &
                                            LAV@ParTable$op == "=~" &
                                            LAV@ParTable$group == k)[1])
                }
            }

            ## NB truncation doesn't work well in stan. instead
            ##    use generated quantities after the fact.
            trunop <- " T(0,)"
            if(target == "jags"){
                for(i in 1:length(fload)){
                    if(LAV@ParTable$prior[fload[i]] != ""){
                        LAV@ParTable$prior[fload[i]] <- paste(LAV@ParTable$prior[fload[i]], trunop, sep="")
                    }
                }
            }
        }
    }

    # if mcmcfile is a directory, vs list, vs logical
    trans.exists <- FALSE
    if(inherits(mcmcfile, "character")){
        jagdir <- mcmcfile
        mcmcfile <- TRUE
    } else if(inherits(mcmcfile, "list")){
        trans.exists <- TRUE
        ## read syntax file
        jagsyn <- readLines(mcmcfile$syntax)
        ## load jagtrans object
        load(mcmcfile$jagtrans)
        ## add new syntax
        jagtrans$model <- paste(jagsyn, collapse="\n")
        ## make sure we don't rewrite the file:
        mcmcfile <- FALSE
        ## we have no idea what they did, so wipe out the priors
        jagtrans$coefvec$prior <- ""
        jagdir <- "lavExport"
    }  else {
        jagdir <- "lavExport"
    }

    # if inits is list
    initsin <- inits
    if(inherits(inits, "list")) initsin <- "jags"

    # extract slots from dummy lavaan object
    lavpartable    <- LAV@ParTable
    if(!("prior" %in% names(lavpartable))) lavpartable$prior <- rep("", length(lavpartable$lhs))
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
            cat("blavaan NOTE: Posterior predictives with missing data are currently very slow.\n\tConsider setting test=\"none\".\n\n")
        }
    }
    if(!jag.do.fit){
      lavoptions$test <- "none"
      lavoptions$se <- "none"
    }
    lavoptions$missing   <- "ml"
    lavoptions$cp        <- cp
    lavoptions$dp        <- dp
    lavoptions$target    <- target

    verbose <- lavoptions$verbose

    # redo estimation + vcov + test
    # 6. estimate free parameters
    start.time <- proc.time()[3]
    x <- NULL
    if(lavmodel@nx.free > 0L) {
        if(!trans.exists){
            ## convert partable to mcmc syntax, then run
            if(target == "jags"){
                jagtrans <- try(lav2mcmc(model = lavpartable, lavdata = lavdata,
                                         cp = cp, lv.x.wish = lavoptions$auto.cov.lv.x,
                                         dp = dp, n.chains = n.chains,
                                         mcmcextra = mcmcextra, inits = initsin,
                                         blavmis = blavmis, wiggle = wiggle,
                                         wiggle.sd = wiggle.sd, target = "jags"),
                                silent = TRUE)
            } else if(target == "stanclassic"){
                jagtrans <- try(lav2stan(model = LAV,
                                         lavdata = lavdata,
                                         dp = dp, n.chains = n.chains,
                                         mcmcextra = mcmcextra,
                                         inits = initsin, wiggle = wiggle,
                                         wiggle.sd = wiggle.sd, debug = mcdebug),
                                silent = TRUE)
            } else if(target == "stancond"){
                jagtrans <- try(lav2stancond(model = LAV,
                                             lavdata = lavdata,
                                             dp = dp, n.chains = n.chains,
                                             mcmcextra = mcmcextra,
                                             inits = initsin,
                                             debug = mcdebug),
                                silent = TRUE)
            } else {
                l2s <- try(lav2stanmarg(lavobject = LAV, dp = dp,
                                        n.chains = n.chains,
                                        inits = initsin, wiggle = wiggle,
                                        wiggle.sd = wiggle.sd, prisamp = prisamp),
                           silent = TRUE)

                if(!inherits(l2s, "try-error")){
                    lavpartable$prior[as.numeric(rownames(l2s$lavpartable))] <- l2s$lavpartable$prior
                    ldargs <- c(l2s$dat, list(lavpartable = l2s$lavpartable, dumlv = l2s$dumlv,
                                              save_lvs = save.lvs))

                    ## add priors to lavpartable, including wiggle
                    if(length(wiggle) > 0){
                      wigrows <- which(l2s$wigpris != "")
                      lavpartable$prior[as.numeric(rownames(l2s$lavpartable))[wigrows]] <- l2s$wigpris[wigrows]
                    }

                    jagtrans <- try(do.call("stanmarg_data", ldargs), silent = TRUE)

                    if(inherits(jagtrans, "try-error")) stop(jagtrans)

                    ## add lkj for unrestricted psi
                    if(jagtrans$fullpsi == 1L){
                      psirows <- which(l2s$lavpartable$mat == "lvrho")
                      lavpartable$prior[as.numeric(rownames(l2s$lavpartable))[psirows]] <- paste0("lkj_corr(", jagtrans$psi_r_alpha[1], ")")
                    }
                    
                    jagtrans <- list(data = jagtrans,
                                     monitors = c("ly_sign",
                                           #"lx_sign",
                                           "bet_sign", "g_sign",
                                           "Theta_cov", "Theta_var",
                                           "Theta_x_cov", "Theta_x_var",
                                           "Psi_cov", "Psi_var",
                                           #"Ph_cov", "Ph_var",
                                           "Nu_free", "Alpha_free"))
                    if("init" %in% names(l2s)){
                      jagtrans <- c(jagtrans, list(inits = l2s$init))
                    }
                } else {
                    jagtrans <- l2s
                }  
            }
        }

        if(!inherits(jagtrans, "try-error")){
            ## add extras to monitor, if specified
            sampparms <- jagtrans$monitors
            if("monitor" %in% names(mcmcextra)){
                sampparms <- c(sampparms, mcmcextra$monitor)
            }
            if(save.lvs & target != "stan") sampparms <- c(sampparms, "eta")

            if(initsin == "jags"){
                jagtrans$inits <- vector("list", n.chains)
            }
            if(inherits(inits, "list")) jagtrans$inits <- inits

            ## add seed to inits; for stan, just add it to
            ## bcontrol
            if(seedlen > 0 & target == "jags"){
                sdinit <- lapply(seed, function(x) list(.RNG.seed = x,
                                                        .RNG.name = "base::Super-Duper"))
                for(i in 1:n.chains){
                    jagtrans$inits[[i]] <- c(jagtrans$inits[[i]], sdinit[[i]])
                }
            } else if(seedlen > 0 & grepl("stan", target)){
                if(!("seed" %in% names(bcontrol))){
                    bcontrol <- c(bcontrol, list(seed = seed))
                }
            }

            if(mcmcfile){
                dir.create(path=jagdir, showWarnings=FALSE)
                fext <- ifelse(target=="jags", "jag", "stan")
                fnm <- paste0(jagdir, "/sem.", fext)
                if(target=="stan"){
                    cat(stanmodels$stanmarg@model_code, file = fnm)
                } else {
                    cat(jagtrans$model, file = fnm)
                }
                if(target=="jags"){
                    save(jagtrans, file = paste(jagdir, "/semjags.rda",
                                                sep=""))
                } else {
                    stantrans <- jagtrans
                    save(stantrans, file = paste(jagdir, "/semstan.rda",
                                                 sep=""))
                }
            }
          
            if(target == "jags"){
              rjarg <- with(jagtrans, list(model = paste(model),
                                           monitor = sampparms, 
                                           data = data, inits = inits))
            } else if(target %in% c("stanclassic", "stancond")){
              rjarg <- with(jagtrans, list(model_code = model,
                                           pars = sampparms,
                                           data = data,
                                           init = inits))
            } else {
              rjarg <- with(jagtrans, list(object = stanmodels$stanmarg,
                                           data = data,
                                           pars = sampparms,
                                           init = inits))
            }

            ## user-supplied jags params
            rjarg <- c(rjarg, mfj, bcontrol)

            if(target == "jags"){
                ## obtain posterior modes
                ## turned off because runjags needs updated
                if(suppressMessages(requireNamespace("modeest", quietly = TRUE))) runjags::runjags.options(mode.continuous = FALSE)
                runjags::runjags.options(force.summary = TRUE)
            }

            if(jag.do.fit){
                if(target == "jags"){
                    rjcall <- "run.jags"
                } else if(target %in% c("stanclassic", "stancond")){
                    cat("Compiling stan model...")
                    rjcall <- "stan"
                } else {
                    rjcall <- "sampling"
                }
                if(convergence == "auto"){
                    rjcall <- "autorun.jags"
                    if("raftery.options" %in% names(rjarg)){
                        ## autorun defaults
                        if(!("r" %in% names(rjarg$raftery.options))){
                            rjarg$raftery.options$r <- .01
                        }
                        if(!("converge.eps" %in% names(rjarg$raftery.options))){
                            rjarg$raftery.options$converge.eps <- .01
                        }
                    } else {
                        rjarg$raftery.options <- list(r=.01, converge.eps=.01)
                    }
                }
                ## the model is run here:
                res <- try(do.call(rjcall, rjarg))
            } else {
                res <- NULL
            }

            if(inherits(res, "try-error")) {
                if(!trans.exists){
                    dir.create(path=jagdir, showWarnings=FALSE)
                    fext <- ifelse(target=="jags", "jag", "stan")
                    cat(jagtrans$model, file = paste(jagdir, "/sem.",
                                                     fext, sep=""))
                    if(target=="jags"){
                        save(jagtrans, file = paste(jagdir, "/semjags.rda",
                                                    sep=""))
                    } else {
                        stantrans <- jagtrans
                        save(stantrans, file = paste(jagdir, "/semstan.rda",
                                                     sep=""))
                    }
                }
                stop("blavaan ERROR: problem with MCMC estimation. The model syntax and data have been exported.")
            }
        } else {
            print(jagtrans)
            stop("blavaan ERROR: problem with translation from lavaan to MCMC syntax.")
        }

        timing$Estimate <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]

        if(target == "jags"){
          parests <- coeffun(lavpartable, jagtrans$pxpartable, res)
          stansumm <- NA
        } else if(target %in% c("stanclassic", "stancond")){
          parests <- coeffun_stan(lavpartable, jagtrans$pxpartable,
                                  res)
          stansumm <- parests$stansumm
        } else {
          parests <- coeffun_stanmarg(lavpartable, lavInspect(LAV, 'free'), l2s$free2, jagtrans$data, res)
          stansumm <- parests$stansumm
        }
        x <- parests$x
        lavpartable <- parests$lavpartable

        if(jag.do.fit){
            lavmodel <- lav_model_set_parameters(lavmodel, x = x)
            if(target == "jags"){
                attr(x, "iterations") <- res$sample
                sample <- res$sample
                burnin <- res$burnin
            } else {
                wrmup <- ifelse(length(rjarg$warmup) > 0,
                                rjarg$warmup, floor(rjarg$iter/2))
                attr(x, "iterations") <- sample
                if(target == "stan"){
                    ## defined variables come from delta method:
                    lavpartable$est <- lav_model_get_parameters(lavmodel = lavmodel, type = "user", extra = TRUE)
                }
                ## lvs now in R instead of Stan
                if(save.lvs & target == "stan"){
                    stanlvs <- samp_lvs(res, lavmodel, parests$lavpartable, jagtrans$data)
                    if(dim(stanlvs)[3L] > 0){
                        lvsumm <- as.matrix(rstan::monitor(stanlvs, print=FALSE))
                        cmatch <- match(colnames(stansumm), colnames(lvsumm))
                        stansumm <- rbind(stansumm, lvsumm[,cmatch])
                    }
                }
                # burnin + sample already defined, will be saved in
                # @external so summary() can use it:
                #burnin <- wrmup
                #sample <- sample - wrmup
            }
            attr(x, "converged") <- TRUE
        } else {
            x <- numeric(0L)
            attr(x, "iterations") <- 0L
            attr(x, "converged") <- FALSE
            lavpartable$est <- lavpartable$start
            stansumm <- NULL
            stanlvs <- NULL
        }
        attr(x, "control") <- bcontrol

        if(!("ordered" %in% dotNames)) {
            attr(x, "fx") <- get_ll(lavmodel = lavmodel, lavpartable = lavpartable,
                                    lavsamplestats = lavsamplestats, lavoptions = lavoptions,
                                    lavcache = lavcache, lavdata = lavdata)[1]
            if(save.lvs & jag.do.fit) {
                if(target == "jags"){
                    fullpmeans <- summary(make_mcmc(res))[[1]][,"Mean"]
                } else {
                    fullpmeans <- stansumm[,"mean"] #rstan::summary(res)$summary[,"mean"]
                }
                cfx <- get_ll(fullpmeans, lavmodel = lavmodel, lavpartable = lavpartable,
                              lavsamplestats = lavsamplestats, lavoptions = lavoptions,
                              lavcache = lavcache, lavdata = lavdata,
                              lavobject = LAV, conditional = TRUE)[1]
            } else {
                cfx <- NULL
            }
        } else {
            attr(x, "fx") <- as.numeric(NA)
        }
    }

    ## parameter convergence + implied moments:
    lavimplied <- NULL
    ## compute/store some model-implied statistics
    lavimplied <- lav_model_implied(lavmodel)
    if(jag.do.fit & n.chains > 1){
      ## this also checks convergence of monitors from mcmcextra, which may not be optimal
      psrfrows <- which(!is.na(lavpartable$psrf) &
                        !is.na(lavpartable$free) &
                        lavpartable$free > 0)
      if(any(lavpartable$psrf[psrfrows] > 1.2)) attr(x, "converged") <- FALSE

      ## warn if psrf is large
      if(!attr(x, "converged") && lavoptions$warn) {
        warning("blavaan WARNING: at least one parameter has a psrf > 1.2.")
      }
    }

    ## fx is mean ll, where ll is marginal log-likelihood (integrate out lvs)
    if(lavoptions$test != "none") {
      cat("Computing posterior predictives...\n")
      lavmcmc <- make_mcmc(res)
      samplls <- samp_lls(res, lavmodel, lavpartable, lavsamplestats,
                          lavoptions, lavcache, lavdata, lavmcmc)
      if(jags.ic) {
        sampkls <- samp_kls(res, lavmodel, lavpartable,
                            lavsamplestats, lavoptions, lavcache,
                            lavdata, lavmcmc, conditional = FALSE)
      } else {
        sampkls <- NA
      }
      
      if(save.lvs) {
        if(target == "stan"){
          lavmcmc <- make_mcmc(res, stanlvs) ## add on lvs
        }
        csamplls <- samp_lls(res, lavmodel, lavpartable,
                             lavsamplestats, lavoptions, lavcache,
                             lavdata, lavmcmc, lavobject = LAV,
                             conditional = TRUE)
        if(jags.ic) {
          csampkls <- samp_kls(res, lavmodel, lavpartable,
                               lavsamplestats, lavoptions, lavcache,
                               lavdata, lavmcmc, lavobject = LAV,
                               conditional = TRUE)
        }
      }
    } else {
      samplls <- NA
      sampkls <- NA
      csamplls <- NA
      csampkls <- NA
    }

    timing$PostPred <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
  
    ## put runjags output in new blavaan slot
    lavjags <- res

    ## 7. VCOV is now simple
    lavvcov <- list()
    VCOV <- NULL
    if(jag.do.fit){
      dsd <- parests$sd[names(parests$sd) %in% colnames(parests$vcorr)]
      if(length(dsd) > 1) dsd <- diag(dsd)
      VCOV <- dsd %*% parests$vcorr %*% dsd
      rownames(VCOV) <- colnames(VCOV) <- colnames(parests$vcorr)
      #lavjags <- c(lavjags, list(vcov = VCOV))

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
                                samplls             = samplls,
                                jagextra            = mcmcextra,
                                stansumm            = stansumm)
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
    ## TODO fix this for stan
    lavpartable$logBF <- rep(NA, length(lavpartable$se))
    if(target == "jags"){
        lavpartable$logBF <- SDBF(lavpartable)
    }

    ## add monitors in mcmcextra as defined variables (except reserved monitors)
    if(length(mcmcextra$monitor) > 0){
        reservemons <- which(mcmcextra$monitor %in% c('deviance', 'pd', 'popt',
                                                     'dic', 'ped', 'full.pd', 'eta'))

        if(length(reservemons) < length(mcmcextra$monitor)){
            jecopy <- mcmcextra
            if(length(reservemons) > 0) jecopy$monitor <- jecopy$monitor[-reservemons]
            lavpartable <- add_monitors(lavpartable, lavjags, jecopy)
        }
    }

    ## 9b. move some stuff from lavfit to optim, for lavaan 0.5-21
    ##     also create external slot
    optnames <- c('x','npar','iterations','converged','fx','fx.group','logl.group',
                  'logl','control')
    lavoptim <- lapply(optnames, function(x) slot(lavfit, x))
    names(lavoptim) <- optnames

    extslot <- list(mcmcout = lavjags, samplls = samplls,
                    origpt = lavpartable, inits = jagtrans$inits,
                    mcmcdata = jagtrans$data, pxpt = jagtrans$pxpartable,
                    burnin = burnin, sample = sample)
    if(grepl("stan", target)){
      extslot <- c(extslot, list(stansumm = stansumm))
      if(save.lvs & target=="stan") extslot <- c(extslot, list(stanlvs = stanlvs))
    }
    if(jags.ic) extslot <- c(extslot, list(sampkls = sampkls))
    if(save.lvs) {
      extslot <- c(extslot, list(cfx = cfx, csamplls = csamplls))
      if(jags.ic) extslot <- c(extslot, list(csampkls = csampkls))
    }
  
    ## move total to the end
    timing$total <- (proc.time()[3] - start.time0)
    tt <- timing$total
    timing <- timing[names(timing) != "total"]
    timing <- c(timing, list(total = tt))
    
    # 10. construct blavaan object
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
                   optim        = lavoptim,
                   implied      = lavimplied,          # list
                   vcov         = lavvcov,
                   external     = extslot,
                   test         = TEST                 # copied for now
                  )
  
    # post-fitting checks
    if(attr(x, "converged")) {
        lavInspect(blavaan, "post.check")
    }

    if(jag.do.fit & lavoptions$warn){
        if(any(blavInspect(blavaan, 'neff') < 100)){
            warning("blavaan WARNING: Small effective sample sizes (< 100) for some parameters.")
        }
    }
    
    blavaan
}

## cfa + sem
bcfa <- bsem <- function(..., cp = "srs", dp = NULL,
    n.chains = 3, burnin, sample, adapt,
    mcmcfile = FALSE, mcmcextra = list(), inits = "prior",
    convergence = "manual", target = "stan", save.lvs = FALSE, wiggle = NULL,
    wiggle.sd = 0.1, prisamp = FALSE, jags.ic = FALSE, seed = NULL,
    bcontrol = list()) {

    dotdotdot <- list(...)
    std.lv <- ifelse(any(names(dotdotdot) == "std.lv"), dotdotdot$std.lv, FALSE)

    mc <- match.call()  
    mc$model.type      = as.character( mc[[1L]] )
    if(length(mc$model.type) == 3L) mc$model.type <- mc$model.type[3L]
    mc$n.chains        = n.chains
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

    ## change defaults depending on jags vs stan
    sampargs <- c("burnin", "sample", "adapt")
    if(target == "jags"){
        defiters <- c(4000L, 10000L, 1000L)
    } else {
        defiters <- c(500L, 1000L, 1000L)
    }
    suppargs <- which(!(sampargs %in% names(mc)))

    if(length(suppargs) > 0){
        for(i in 1:length(suppargs)){
            mc[[(length(mc)+1)]] <- defiters[suppargs[i]]
            names(mc)[length(mc)] <- sampargs[suppargs[i]]
        }
    }

    eval(mc, parent.frame())
}

# simple growth models
bgrowth <- function(..., cp = "srs", dp = NULL,
    n.chains = 3, burnin, sample, adapt,
    mcmcfile = FALSE, mcmcextra = list(), inits = "prior",
    convergence = "manual", target = "stan", save.lvs = FALSE, wiggle = NULL,
    wiggle.sd = 0.1, prisamp = FALSE, jags.ic = FALSE, seed = NULL,
    bcontrol = list()) {

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

    ## change defaults depending on jags vs stan
    sampargs <- c("burnin", "sample", "adapt")
    if(target == "jags"){
        defiters <- c(4000L, 10000L, 1000L)
    } else {
        defiters <- c(500L, 1000L, 1000L)
    }
    suppargs <- which(!(sampargs %in% names(mc)))

    if(length(suppargs) > 0){
        for(i in 1:length(suppargs)){
            mc[[(length(mc)+1)]] <- defiters[suppargs[i]]
            names(mc)[length(mc)] <- sampargs[suppargs[i]]
        }
    }
    
    eval(mc, parent.frame())
}
