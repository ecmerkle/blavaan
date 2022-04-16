### Main author: Terrence D. Jorgensen
### function to implement a posterior predictor model check using any
### discrepancy function that can be applied to a lavaan object
### Ed Merkle: various additions to keep things working with blavaan


## -----------------
## Class and Methods
## -----------------

setClass("blavPPMC",
         representation(discFUN = "list",  # list of discrepancy functions
                        dims = "list",     # dimensions of each discrepancy function
                        PPP = "list",      # PPP-value for each discFUN
                        obsDist = "list",  # posterior distribution of realized values
                        simDist = "list")) # posterior predictive distribution

summary.ppmc <- function(object, discFUN, dist = c("obs","sim"),
                         central.tendency = c("mean","median","mode"),
                         hpd = TRUE, prob = .95,
                         to.data.frame = FALSE, diag = TRUE,
                         sort.by = NULL, decreasing = FALSE) {
  ## check choices
  if (!is.character(central.tendency)) {
    stop('blavaan ERROR: central.tendency must be a character vector')
  }
  central.tendency <- tolower(central.tendency)
  invalid <- setdiff(central.tendency, c("mean","median","mode","eap","map"))
  if (length(invalid)) {
    stop('blavaan ERROR: Invalid choice(s) of central.tendency:\n"',
         paste(invalid, collapse = '", "'), '"')
  }

  dist <- as.character(dist)[1]

  if (missing(discFUN)) discFUN <- names(object@discFUN)[1]
  if (is.numeric(discFUN)) {
    discFUN <- names(object@discFUN)[ as.integer(discFUN)[1] ]
  } else discFUN <- as.character(discFUN[1])
  if (!discFUN %in% discFUN) stop('Invalid choice of discFUN. Available ',
                                  'choices are:\n\t',
                                  paste0(names(object@discFUN), collapse = "\n\t"))
  XX <- slot(object, paste0(dist, "Dist"))[[discFUN]]

  ## collect list of functions to apply to each index
  getSummaries <- function(x) {
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
        tryMode <- try(modeest::mlv(x, method = "kernel", na.rm = TRUE),
                       silent = TRUE)
        if (!inherits(tryMode, "try-error") && is.numeric(tryMode)) {
          out <- c(out, MAP = tryMode[1])
        }
      }
      ## if the mode was not obtained above, use the quick-and-dirty way
      if (is.na(out["MAP"])) {
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

  ## apply function to selected discFUN, then append PPP
  if (length(XX[[1]]) == 1L) {
    ## scalar, store in vector
    out <- c(getSummaries(do.call(c, XX)),
             "PPP_sim>obs" = object@PPP[[discFUN]],
             "PPP_sim<obs" = 1 - object@PPP[[discFUN]])
    class(out) <- c("lavaan.vector","numeric")

  } else if (is.null(dim(XX[[1]]))) {
    ## vector, store in data.frame (1 row per discrepancy, 1 column per summary)
    out <- data.frame(t(apply(do.call(rbind, XX), 2, getSummaries)))
    out <- cbind(out, "PPP_sim_GreaterThan_obs" = object@PPP[[discFUN]],
                 "PPP_sim_LessThan_obs" = 1 - object@PPP[[discFUN]])
    ## sort?
    if (length(sort.by)) {
      sort.by <- as.character(sort.by)[1]
      if (sort.by %in% colnames(out)) {
        ORD <- order(out[ , sort.by], decreasing = decreasing)
        out <- out[ORD, ]
      }
    }
    class(out) <- c("lavaan.data.frame","data.frame")
    attr(out, "header") <- paste0("Posterior summary statistics",
                                  ifelse(hpd, {
                                    paste0(" and highest posterior density (HPD) ",
                                           round(prob*100, 2),
                                           "% credible intervals")
                                  }, ""),
                                  " for the posterior ",
                                  if (dist == "sim") "predictive ",
                                  "distribution of ",
                                  if (dist == "obs") "realized ",
                                  "discrepancy-function values based on ",
                                  if (dist == "obs") "observed ",
                                  if (dist == "sim") {
                                    "data simulated from the posterior, "
                                  } else "data, ",
                                  "along with posterior predictive p values ",
                                  "to test hypotheses in either direction:\n",
                                  sep = "")
  } else {
    ## multidimensional array (including matrices)
    ## vectorize to apply same calculation as above, reapply attributes to each
    XXvecs <- lapply(XX, as.numeric)
    XXmat <- apply(do.call(rbind, XXvecs), 2, getSummaries)
    out <- sapply(rownames(XXmat), function(n) XXmat[n,], simplify = FALSE)
    for (nn in seq_along(out)) attributes(out[[nn]]) <- attributes(XX[[1]])

    ## 2-dimensional matrix provides further options for customizing output
    if (length(dim(XX[[1]])) == 2L) {

      if (isSymmetric(out[[1]]) && to.data.frame) {
        ## store unique elements of symmetric matrix in data.frame
        XXidx <- which(lower.tri(out[[1]], diag = diag), arr.ind = TRUE)
        XXnames <- rownames(out[[1]])
        DF <- data.frame(moment = paste0(XXnames[XXidx[,2]], "~~", XXnames[XXidx[,1]]))
        for (nn in names(out)) DF[[nn]] <- out[[nn]][XXidx]
        DF$PPP_sim_GreaterThan_obs <- object@PPP[[discFUN]][XXidx]
        DF$PPP_sim_LessThan_obs <- 1 - object@PPP[[discFUN]][XXidx]
        out <- DF
        ## sort?
        if (length(sort.by)) {
          sort.by <- as.character(sort.by)[1]
          if (sort.by %in% colnames(out)) {
            ORD <- order(out[ , sort.by], decreasing = decreasing)
            out <- out[ORD, ]
          }
        }
        class(out) <- c("lavaan.data.frame","data.frame")
        attr(out, "header") <- paste0("Posterior summary statistics",
                                      ifelse(hpd, {
                                        paste0(" and highest posterior density (HPD) ",
                                               round(prob*100, 2),
                                               "% credible intervals")
                                      }, ""),
                                      " for the posterior ",
                                      if (dist == "sim") "predictive ",
                                      "distribution of ",
                                      if (dist == "obs") "realized ",
                                      "discrepancy-function values based on ",
                                      if (dist == "obs") "observed ",
                                      if (dist == "sim") {
                                        "data simulated from the posterior, "
                                      } else "data, ",
                                      "along with posterior predictive p values ",
                                      "to test hypotheses in either direction:\n",
                                      sep = "")

      } else if (to.data.frame) {
        ## add PPP matrices to list of output
        out[["PPP_sim>obs"]] <- object@PPP[[discFUN]]
        out[["PPP_sim<obs"]] <- 1 - object@PPP[[discFUN]]

        ## asymmetric, but each matrix could be a (sortable) data.frame
        for (mm in seq_along(out)) {
          out[[mm]] <- as.data.frame(out[[mm]])
          class(out[[mm]]) <- c("lavaan.data.frame","data.frame")
          if (names(out)[[mm]] %in% c("EAP","Median","MAP","SD")) {
            attr(out[[mm]], "header") <- paste0("The ",
                                                ifelse(names(out)[[mm]] == "EAP", "mean",
                                                ifelse(names(out)[[mm]] == "MAP", "mode",
                                                ifelse(names(out)[[mm]] == "SD", "SD","median"))),
                                          " of the posterior ",
                                          if (dist == "sim") "predictive ",
                                          "distribution of ",
                                          if (dist == "obs") "realized ",
                                          "discrepancy-function values based on ",
                                          if (dist == "obs") "observed ",
                                          if (dist == "sim") {
                                            "data simulated from the posterior, "
                                          } else "data:\n",
                                          sep = "")
          } else if (names(out)[[mm]] %in% c("lower","upper")) {
            attr(out[[mm]], "header") <- paste0("Highest posterior density (HPD) ",
                                                round(prob*100, 2), "% ",
                                                toupper(names(out)[[mm]]),
                                                " credible-interval limits for the posterior ",
                                                if (dist == "sim") "predictive ",
                                                "distribution of ",
                                                if (dist == "obs") "realized ",
                                                "discrepancy-function values based on ",
                                                if (dist == "obs") "observed ",
                                                if (dist == "sim") {
                                                  "data simulated from the posterior, "
                                                } else "data:\n",
                                                sep = "")
          } else {
            attr(out[[mm]], "header") <- paste0("Posterior predictive p values for ",
                                                "testing directional hypotheses:\n\n",
                                                sep = "")
          }

          if (length(sort.by)) {
            sort.by <- as.character(sort.by)[1]
            if (sort.by %in% colnames(out[[mm]])) {
              ORD <- order(out[[mm]][ , sort.by], decreasing = decreasing)
              out[[mm]] <- out[[mm]][ORD, ]
            }
          }

        }

      } else {
        ## add PPP matrices to list of output
        out[["PPP_sim>obs"]] <- object@PPP[[discFUN]]
        out[["PPP_sim<obs"]] <- 1 - object@PPP[[discFUN]]

        ## convert all matrices to lavaan for printing options
        for (mm in seq_along(out)) {
          class(out[[mm]]) <- c(ifelse(isSymmetric(out[[mm]]),
                                       "lavaan.matrix.symmetric",
                                       "lavaan.matrix"), "matrix")
        }

      }

    } else {
      ## must be a multidimensional array
      ## add PPP to list of output
      out[["PPP_sim>obs"]] <- object@PPP[[discFUN]]
      out[["PPP_sim<obs"]] <- 1 - object@PPP[[discFUN]]
    }

  }

  out
}

summary.blavPPMC <- function(object, ...) {
  summary.ppmc(object, ...)
}

setMethod("summary", "blavPPMC", summary.blavPPMC)


setMethod("show", "blavPPMC", function(object) {
  ## header states what discrepancy functions are available
  cat('Posterior predictive model checks available, using the following ',
      'discrepancy functions:\n\t',
      paste0(names(object@discFUN), '\t(',
             ifelse(sapply(object@dims, length) == 1, "length == ", "dim: ["),
             lapply(object@dims, paste, collapse = ","),
             ifelse(sapply(object@dims, length) == 1, "", "]"), ')',
             collapse = "\n\t"),
      "\n\n", sep = "")
  ## Instruct how to obtain summaries
  cat('To request numerical summaries of the posterior distribution of realized ',
      'values of a discrepancy function (such as measures of central tendency, ',
      'credible intervals, and the posterior predictive p value in reference ',
      'to the posterior predictive distribution), use the summary() method. ',
      'The "discFUN" argument can be used to name a particular discrepancy ',
      'function of interest, or else the first one will be returned.\n\n',
      sep = "")
  ## Instruct how to obtain plots
  cat('To request graphical summaries of the posterior distribution of realized ',
      'values of a discrepancy function (i.e., plot realized values by expected ',
      'values, as represented by the posterior predictive distribution), use ',
      'the plot() method. The "discFUN" argument must be used to name a ',
      'particular discrepancy function of interest (or else the first one will ',
      'be returned), and the "element" argument is used to select a particular ',
      'element returned by that function (or else the first element will be ',
      'returned). \n\nSee the ?ppmc help page for examples.\n\n',
      sep = "")

  invisible(object)
})


## S3 plot() method
plot.blavPPMC <- function(x, ..., discFUN, element, central.tendency = "",
                          hpd = TRUE, prob = .95, nd = 3) {
  ## check discFUN
  if (missing(discFUN)) discFUN <- names(x@discFUN)[1]
  if (is.numeric(discFUN)) {
    discFUN <- names(x@discFUN)[ as.integer(discFUN)[1] ]
  } else discFUN <- as.character(discFUN[1])
  if (!discFUN %in% discFUN) stop('Invalid choice of discFUN. Available ',
                                  'choices are:\n\t',
                                  paste0(names(x@discFUN), collapse = "\n\t"))

  ## extract element from discFUN
  if (missing(element)) element <- 1L
  element <- as.vector(element)
  if (length(element) > 1L) {
    if (length(element) != length(dim(x@obsDist[[discFUN]][[1]])))
      stop('element must have length 1 or the number of dimensions of the ',
           'discFUN object.')
    element <- matrix(element, nrow = 1) # make sure it is a row matrix
  }
  OBS <- sapply(x@obsDist[[discFUN]], function(i) i[element])
  SIM <- sapply(x@simDist[[discFUN]], function(i) i[element])
  PPP <- x@PPP[[discFUN]][element]
  if (is.numeric(element)) {
    if (is.null(dim(x@PPP[[discFUN]]))) {
      ## vector, try names()
      NAMES <- names(x@PPP[[discFUN]])[element]
    } else {
      ## multidimensional array, try dimnames()
      NAMES <- mapply(function(i, j) i[j],
                      i = dimnames(x@PPP[[discFUN]]), j = element)
    }
    if (length(NAMES) < length(element)) NAMES <- element
  } else NAMES <- as.character(element)
  idx <- paste0(discFUN, "[",
                if (is.character(NAMES)) '"',
                paste(NAMES,
                      collapse = ifelse(is.character(NAMES), '","', ", ")),
                if (is.character(NAMES)) '"',
                "]")

  ## specify arguments for plot.default (scatterplot)
  dots <- list(...)
  if (is.null(dots$main)) dots$main <- bquote(PPP[sim>obs]~"="~.(round(PPP, nd))~(PPP[obs>sim]~"="~.(round(1 - PPP, nd))))
  if (is.null(dots$xlab)) dots$xlab <- paste('Realized Values of', idx)
  if (is.null(dots$ylab)) dots$ylab <- 'Posterior Predicted Values'
  if (is.null(dots$xlim)) dots$xlim <- range(c(OBS, SIM))
  if (is.null(dots$ylim)) dots$ylim <- range(c(OBS, SIM))
  dots$x <- OBS
  dots$y <- SIM
  ## call plot
  do.call(graphics::plot.default, dots)
  dots$x <- dots$y <- dots$lty <- NULL

  ## plot diagonal line
  dots$a <- 0
  dots$b <- 1
  do.call(graphics::abline, dots)
  dots$a <- dots$b <- NULL

  ## collect list of functions to apply to each index
  getSummaries <- function(x) {
    out <- numeric(0)

    if (any(central.tendency[1] %in% c("mean","eap"))) {
      out <- c(out, EAP = mean(x, na.rm = TRUE))
    }

    if (any(central.tendency[1] %in% "median")) {
      out <- c(out, Median = median(x, na.rm = TRUE))
    }

    if (any(central.tendency[1] %in% c("mode","map"))) {
      ## can the modeest package be used?
      if (suppressMessages(requireNamespace("modeest", quietly = TRUE))) {
        tryMode <- try(modeest::mlv(x, method = "kernel", na.rm = TRUE),
                       silent = TRUE)
        if (!inherits(tryMode, "try-error") && is.numeric(tryMode)) {
          out <- c(out, MAP = tryMode[1])
        }
      }
      if (is.na(out["MAP"])) {
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
  obsSumm <- getSummaries(OBS)
  simSumm <- getSummaries(SIM)
  ## plot central.tendency
  if (length(central.tendency)) if (central.tendency[1] != "") {
    dots$h <- simSumm[1]
    dots$v <- obsSumm[1]
    dots$lty <- "dashed"
    do.call(graphics::abline, dots)
  }
  ## plot HPD-CI limits
  if (hpd) {
    #dots$h <- simSumm[c("lower","upper")]
    dots$v <- simSumm[c("lower","upper")]
    dots$lty <- "dotted"
    do.call(graphics::abline, dots)
    polygon(x = c(simSumm["lower"], simSumm["lower"], simSumm["upper"], simSumm["upper"]),
            y = c(dots$ylim[1], dots$ylim[2], dots$ylim[2], dots$ylim[1]),
            col = adjustcolor(if (is.null(dots$col)) "grey" else dots$col,
                              alpha.f = 0.25), border = NA)
  }

  invisible(NULL)
}


## S3 hist() method
hist.blavPPMC <- function(x, ..., discFUN, element, hpd = TRUE, prob = .95,
                          printLegend = TRUE, legendArgs = list(x = "topleft"),
                          densityArgs = list(), nd = 3) {
  ## check discFUN
  if (missing(discFUN)) discFUN <- names(x@discFUN)[1]
  if (is.numeric(discFUN)) {
    discFUN <- names(x@discFUN)[ as.integer(discFUN)[1] ]
  } else discFUN <- as.character(discFUN[1])
  if (!discFUN %in% discFUN) stop('Invalid choice of discFUN. Available ',
                                  'choices are:\n\t',
                                  paste0(names(x@discFUN), collapse = "\n\t"))

  ## extract element from discFUN
  if (missing(element)) element <- 1L
  element <- as.vector(element)
  if (length(element) > 1L) {
    if (length(element) != length(dim(x@obsDist[[discFUN]][[1]])))
      stop('element must have length 1 or the number of dimensions of the ',
           'discFUN object.')
    element <- matrix(element, nrow = 1) # make sure it is a row matrix
  }
  OBS <- sapply(x@obsDist[[discFUN]], function(i) i[element])
  SIM <- sapply(x@simDist[[discFUN]], function(i) i[element])
  PPP <- x@PPP[[discFUN]][element]
  if (is.numeric(element)) {
    if (is.null(dim(x@PPP[[discFUN]]))) {
      ## vector, try names()
      NAMES <- names(x@PPP[[discFUN]])[element]
    } else {
      ## multidimensional array, try dimnames()
      NAMES <- mapply(function(i, j) i[j],
                      i = dimnames(x@PPP[[discFUN]]), j = element)
    }
    if (length(NAMES) < length(element)) NAMES <- element
  } else NAMES <- as.character(element)
  idx <- paste0('Discrepancy Function: ', discFUN, "[",
                if (is.character(NAMES)) '"',
                paste(NAMES,
                      collapse = ifelse(is.character(NAMES), '","', ", ")),
                if (is.character(NAMES)) '"',
                "]")

  ## grab density limits from 2 densities to set ylim
  densityArgs$x <- SIM
  H0 <- do.call(stats::density.default, densityArgs)
  densityArgs$x <- OBS
  H1 <- do.call(stats::density.default, densityArgs)

  ## specify arguments for hist()
  dots <- list(...)
  dots$x <- H1$x # plot OBS first
  dots$y <- H1$y
  dots$type = "l"
  dots$lty <- 'solid'
  if (is.null(dots$main)) dots$main <- bquote(PPP[sim>obs]~"="~.(round(PPP, nd))~(PPP[obs>sim]~"="~.(round(1 - PPP, nd))))
  if (is.null(dots$xlab)) dots$xlab <- idx
  if (is.null(dots$ylab)) dots$ylab <- 'Density'
  if (is.null(dots$yaxt)) dots$yaxt <- 'n' # suppress y axis
  if (is.null(dots$xlim)) dots$xlim <- range(c(H1$x, H0$x))
  if (is.null(dots$ylim)) dots$ylim <- range(c(H1$y, H0$y))
  if (is.null(dots$lwd)) dots$lwd <- 1.5
  if (is.null(dots$bty)) dots$bty <- 'n'

  ## specify arguments for lines
  lineArgs <- dots
  lineArgs$type <- NULL
  lineArgs$x <- H0$x
  lineArgs$y <- H0$y
  lineArgs$lty <- 'dashed'

  ## specify arguments for abline, if HPD-CI limits requested
  if (hpd) {
    if (!"package:coda" %in% search()) attachNamespace("coda")
    ## add HPD-CI limits
    ablineArgs <- lineArgs
    ablineArgs$x <- NULL
    ablineArgs$y <- NULL
    ablineArgs$v <- HPDinterval(as.mcmc(SIM), prob = prob)[1, ]
    ablineArgs$lty <- 'dotted'
  } else ablineArgs <- NULL

  ## specify arguments for legend, if requested
  if (printLegend) {
    if (is.null(legendArgs$bty)) legendArgs$bty <- "n"
    if (is.null(legendArgs$box.lty)) legendArgs$box.lty <- 0
    legendArgs$lty <- c("solid", "blank", "dashed")
    if (hpd) legendArgs$lty <- c(legendArgs$lty, "blank", "dotted")
    legendArgs$lwd <- c(dots$lwd, 0, dots$lwd)
    if (hpd) legendArgs$lwd <- c(legendArgs$lwd, 0, dots$lwd)
    if (is.null(legendArgs$col) && !is.null(dots$col)) {
      legendArgs$col <- c(dots$col, "", dots$col)
      if (hpd) legendArgs$col <- c(legendArgs$col, 0, dots$col)
    }
    legendArgs$legend <- c("Realized\nValues", "", "Posterior\nPredictions")
    if (hpd) legendArgs$legend <- c(legendArgs$legend, "",
                                    paste0("Expected ", round(100*prob),
                                           "%\nCredible Limits"))
  } else legendArgs <- NULL

  ## plot
  do.call(graphics::plot.default, dots)      # OBS
  do.call(graphics::lines.default, lineArgs) # SIM
  if (hpd) do.call(graphics::abline, ablineArgs)
  if (printLegend) do.call(graphics::legend, legendArgs)

  ## return arguments to create plot (and optionally, legend)
  invisible(list(plot = dots, lines = lineArgs, abline = ablineArgs,
                 legend = legendArgs, density = densityArgs))
}

## S3 hist() method
pairs.blavPPMC <- function(x, discFUN, horInd = 1:DIM, verInd = 1:DIM,
                           printLegend = FALSE, ...) {
  ## check discFUN
  if (missing(discFUN)) discFUN <- names(x@discFUN)[1]
  if (is.numeric(discFUN)) {
    discFUN <- names(x@discFUN)[ as.integer(discFUN)[1] ]
  } else discFUN <- as.character(discFUN[1])
  if (!discFUN %in% discFUN) stop('Invalid choice of discFUN. Available ',
                                  'choices are:\n\t',
                                  paste0(names(x@discFUN), collapse = "\n\t"))

  ## check that pairs() can be applied to a 2-dim matrix
  DIM <- dim(x@PPP[[discFUN]])
  if (!is.matrix(x@PPP[[discFUN]]))
    stop('The pairs() method can only be applied when discFUN returns a',
         '2-dimensional matrix/array. "', discFUN, '" ',
         if (is.null(DIM)) 'is a vector.' else {
           paste('has', length(DIM), 'dimensions.')
         })

  ## set up grid
  DIM # evaluate promise
  opar <- par(mfrow = c(length(horInd), length(verInd))); on.exit(par(opar))

  ## symmetric?
  SYM <- length(horInd) == length(verInd)
  if (SYM) SYM <- all(sort(horInd) == sort(verInd))
  if (SYM) verInd <- horInd

  ## loop over indices
  for (RR in seq_along(horInd)) for (CC in seq_along(verInd)) {
    if (RR > CC || !SYM) {
      ## scatterplots in lower triangle (or whole thing if !SYM)
      plot(x, discFUN = discFUN,
           element = paste0("x", c(horInd[RR], verInd[CC])))

    } else if (RR == CC) {
      ## names on diagonal
      plot.new()
      legend("center", paste0("x", RR), bty = "n", cex = 5)

    } else {
      ## histograms above the diagonal
      hist(x, discFUN = discFUN, printLegend = printLegend,
           element = paste0("x", c(horInd[RR], verInd[CC])))
    }

  }

  invisible(NULL)
}



## --------------------
## Constructor Function
## --------------------

## Public function, wrapper around hidden function postpred()
ppmc <- function(object, thin = 1, fit.measures = c("srmr","chisq"),
                 #baseline.model = NULL,
                 discFUN = NULL, conditional = FALSE) {

  ## check custom discrepancy function(s)
  if (!is.null(discFUN)) {
    if (!is.list(discFUN)) {
      discFUN <- list(discFUN)
      funcName <- substitute(discFUN) # if the function is an object, save name
      if (is.name(funcName)) names(discFUN) <- as.character(funcName)
    }
    if (!all(sapply(discFUN, is.function)))
      stop('blavaan ERROR: The "discFUN" argument must be a (list of) function(s).')
  }

  ## differentiate between multiple possible chisq stats
  if (length(fit.measures) == 1L & blavInspect(object, "categorical")) {
    if (fit.measures == "chisq") {
      ## add another so that the dwls chisq from lavaan is used
      ## (due to the way blav_model_loglik is structured)
      fit.measures <- c("chisq", "chisq.scaled")
    } else if (fit.measures == "marglogl") {
      ## approximate the marginal lrt
      fit.measures <- "chisq"
    }
  }

  jagtarget <- lavInspect(object, "options")$target == "jags"

  if(jagtarget){
    etas <- any(object@external$mcmcout$monitor == "eta")
  } else {
    etas <- any(grepl("^eta", rownames(object@external$stansumm)))
  }
  if (etas & conditional) {
    if (!length(discFUN)) warning("blavaan WARNING: conditional=TRUE has no effect if you do not supply a discFUN.")
  }
  
  out <- postpred(samplls = object@external$samplls, lavobject = object, 
                  measure = fit.measures, thin = thin, discFUN = discFUN)

  ## "out" is a list:
  ## - $ppval contains PPP-values, nested in a list per discFUN if length(discFUN) > 1
  ## - $ppdist is a list containing:
  ##    - $obs list of realized values (one per sample), nested within a
  ##           list per discFUN if length(discFUN) > 1
  ##    - $reps same as $obs, but the posterior predictive distribution

  ## make "out" conform to blavPPMC class representation
  if (length(discFUN)) {
    DIMS <- PPP <- OBS <- SIM <- vector("list", length(discFUN))
    ## loop over discrepancies
    for (d in seq_along(discFUN)) {
      PPP[[d]] <- out$ppval[[d]]
      DIMS[[d]] <- if (is.null(dim(PPP[[d]]))) length(PPP[[d]]) else dim(PPP[[d]])
      OBS[[d]] <- out$ppdist$obs[[d]]
      SIM[[d]] <- out$ppdist$reps[[d]]
    }

  } else {
    ## scalar or vector from fitMeasures(), store as only element in a list
    PPP <- list(out$ppval)
    DIMS <- if (is.null(dim(PPP[[1]]))) {
      list(length(PPP[[1]]))
    } else list(dim(PPP[[1]]))
    OBS <- list(out$ppdist$obs)
    SIM <- list(out$ppdist$reps)
  }

  ## to store fit.measures choices in the discFUN slot, manufacture a call
  if (is.null(discFUN)) {
    discFUN <- eval(parse(text = paste0('list(fit.indices = function(fit) ',
                                        'lavaan::fitMeasures(fit, fit.measures = c("',
                                        paste(fit.measures, collapse = '","'),
                                        '")))')))
  }

  ## names for each discFUN
  if (is.null(names(discFUN))) names(discFUN) <- paste0("discFUN", seq_along(discFUN))

  ## pass names to other lists
  for (d in seq_along(discFUN)) {
    names(DIMS) <- names(PPP) <- names(OBS) <- names(SIM) <- names(discFUN)
  }

  obj <- new("blavPPMC",
             discFUN = discFUN, # "list", # list of discrepancy functions
             dims = DIMS,   # "list", # dimensions of each discrepancy function
             PPP = PPP,     # "list", # PPP-value for each discFUN
             obsDist = OBS, # "list", # posterior distribution of realized values
             simDist = SIM) # "list"  # posterior predictive distribution
  obj
}






