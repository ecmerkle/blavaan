set_phantoms <- function(partable, parvec, ov.names, lv.names, ov.names.x, lv.names.x, ov.cp, lv.cp, lv.x.wish, ngroups) {
  ## Add phantom lvs for covariance parameters

  ## exclude lv.x if we are using dmnorm/dwish:
  if(lv.x.wish){
    covpars <- which(partable$op == "~~" &
                     partable$lhs != partable$rhs &
                     partable$group == 1 &
                     !(partable$lhs %in% ov.names.x &
                       partable$free == 0) &
                     !(partable$lhs %in% lv.names.x))
  } else {
    if(lv.cp == "srs"){
      covpars <- which(partable$op == "~~" &
                       partable$lhs != partable$rhs &
                       partable$group == 1 &
                       !(partable$lhs %in% ov.names.x &
                         partable$free == 0))
    }
    if(lv.cp == "fa"){
      covpars <- which(partable$op == "~~" &
                       partable$lhs != partable$rhs &
                       partable$group == 1 &
                       !(partable$lhs %in% ov.names.x &
                         partable$free == 0) &
                       !(partable$lhs %in% lv.names &
                         partable$free == 0))
    }
  }

  ## remove covariances fixed to zero under fa
  ## if(lv.cp == "fa"){
  ##   fcovs <- which(partable$free[covpars] == 0)
  ##   zcovs <- which(partable$ustart[covpars][fcovs] == 0)
  ##   covpars <- covpars[-fcovs[zcovs]]
  ## }

  blkrow <- rep(NA, length(partable$id))

  ## Hold jags code to convert phantom parameters back to
  ## the original error variances & covariances
  TXT <- ""

  ## Only do this if covpars exist
  if(length(covpars) > 0){
    ## add to model matrices
    ## added entries in lambda vs in beta
    nmvcovs <- sum(partable$lhs[covpars] %in% ov.names)
    nlvcovs <- length(covpars) - nmvcovs
    patts <- attributes(parvec)
    for(k in 1:ngroups){
      if(!("lambda" %in% patts$mmNames[[k]]) & nmvcovs > 0){
        lcolstart <- 0
        attributes(parvec)$mmNames[[k]] <- c(patts$mmNames[[k]],
                                             "lambda")
        attributes(parvec)$mmRows[[k]] <- c(patts$mmRows[[k]], lambda=length(ov.names))
        attributes(parvec)$mmCols[[k]] <- c(patts$mmCols[[k]], lambda=nmvcovs)
      } else {
        lcolstart <- patts$mmCols[[k]]["lambda"]
        attributes(parvec)$mmCols[[k]]["lambda"] <- patts$mmCols[[k]]["lambda"] + nmvcovs
      }
      if(!("beta" %in% patts$mmNames[[k]]) & nlvcovs > 0){
        bcolstart <- 0
        attributes(parvec)$mmNames[[k]] <- c(patts$mmNames[[k]], "beta")
        attributes(parvec)$mmRows[[k]] <- c(patts$mmRows[[k]], beta=
nlvcovs)
        attributes(parvec)$mmCols[[k]] <- c(patts$mmCols[[k]], beta=nlvcovs)
      } else {
        bcolstart <- patts$mmCols[[k]]["beta"]
        attributes(parvec)$mmRows[[k]]["beta"] <- patts$mmRows[[k]]["beta"] + nlvcovs
        attributes(parvec)$mmCols[[k]]["beta"] <- patts$mmCols[[k]]["beta"] + nlvcovs
      }
      if(!("psi" %in% patts$mmNames[[k]])){
        thetcolstart <- 0
        attributes(parvec)$mmNames[[k]] <- c(patts$mmNames[[k]], "psi")
        attributes(parvec)$mmRows[[k]] <- c(patts$mmRows[[k]], psi=length(covpars))
        attributes(parvec)$mmCols[[k]] <- c(patts$mmCols[[k]], psi=length(covpars))
      } else {
        thetcolstart <- patts$mmCols[[k]]["psi"]
        attributes(parvec)$mmRows[[k]]["psi"] <- patts$mmRows[[k]]["psi"] + length(covpars)
        attributes(parvec)$mmCols[[k]]["psi"] <- patts$mmCols[[k]]["psi"] + length(covpars)
      }
    }
    
    cprm <- NULL

    ## which covariances are under srs?
    ridx <- 1:length(covpars)
        
    for(k in 1:ngroups){
      rhoind <- 1
      tlcs <- lcolstart
      tbcs <- bcolstart
      ttcs <- thetcolstart
      for(i in 1:length(covpars)){
        ## Find the row for group k (needed for multiple groups)
        covparg <- which(partable$op == "~~" &
                         partable$lhs == partable$lhs[covpars[i]] &
                         partable$rhs == partable$rhs[covpars[i]] &
                         partable$group == k)

        ## Is this constrained equal to a previous parameter?
        eq.const <- FALSE
        grp.idx <- k
        eq.idx <- which(partable$op == "==" & partable$rhs == partable$plabel[covparg])
        if(length(eq.idx) > 0){
          eq.const <- TRUE
          ## TODO? assumes it is equal to another covariance; do any models
          ## restrict covariances to be equal to other types of parameters?
          full.idx <- which(partable$plabel == partable$lhs[eq.idx])
          old.idx <- which(partable$lhs[covpars] == partable$lhs[full.idx[1]] &
                           partable$rhs[covpars] == partable$rhs[full.idx[1]])
          old.ridx <- ridx[old.idx]
          grp.idx <- partable$group[full.idx[1]]
        }
          
        tmprows <- nrow(partable) + 1:3
        phname <- paste(".phant", i, sep="")
        partable <- rbind(partable, blkrow, blkrow, blkrow)

        partable$group[tmprows] <- k

        partable$rhs[tmprows[1]] <- partable$lhs[covpars[i]]
        partable$rhs[tmprows[2]] <- partable$rhs[covpars[i]]
        ## Decide on =~ (ov) vs ~ (lv)
        if(partable$lhs[covpars[i]] %in% ov.names){
          partable$lhs[tmprows[1]] <- phname
          partable$op[tmprows[1]] <- "=~"
          partable$rhs[tmprows[1]] <- partable$lhs[covpars[i]]
          partable$mat[tmprows[1]] <- "lambda"
          partable$row[tmprows[1]] <- match(partable$lhs[covpars[i]],
                                            patts$mmDimNames[[k]]$lambda[[1]])
          tlcs <- tlcs + 1
          partable$col[tmprows[1]] <- tlcs
          tmpv1 <- paste("theta[", match(partable$lhs[covpars[i]], patts$mmDimNames[[k]]$lambda[[1]]), ",", k, "]", sep="")
          if(eq.const){
            oldv1 <- paste("theta[", match(partable$lhs[full.idx], patts$mmDimNames[[k]]$lambda[[1]]), ",", grp.idx, "]", sep="")
          }
          ctype <- "ov"
        } else {
          partable$lhs[tmprows[1]] <- partable$lhs[covpars[i]]
          partable$op[tmprows[1]] <- "~"
          partable$rhs[tmprows[1]] <- phname
          partable$mat[tmprows[1]] <- "beta"
          tbcs <- tbcs + 1
          partable$row[tmprows[1]] <- tbcs
          partable$col[tmprows[1]] <- match(partable$lhs[covpars[i]],
                                            patts$mmDimNames[[k]]$psi[[1]])
          tmpv1 <- paste("psi[", match(partable$lhs[covpars[i]], patts$mmDimNames[[k]]$psi[[1]]), ",", k, "]", sep="")
          if(eq.const){
            oldv1 <- paste("psi[", match(partable$lhs[full.idx], patts$mmDimNames[[k]]$psi[[1]]), ",", grp.idx, "]", sep="")
          }
          ctype <- "lv"
        }
        if(partable$rhs[covpars[i]] %in% ov.names){
          partable$lhs[tmprows[2]] <- phname
          partable$op[tmprows[2]] <- "=~"
          partable$rhs[tmprows[2]] <- partable$rhs[covpars[i]]
          partable$mat[tmprows[2]] <- "lambda"
          partable$row[tmprows[2]] <- match(partable$rhs[covpars[i]],
                                            patts$mmDimNames[[k]]$lambda[[1]])
          tlcs <- tlcs + 1
          partable$col[tmprows[2]] <- tlcs
          tmpv2 <- paste("theta[", match(partable$rhs[covpars[i]], ov.names), ",", k, "]", sep="")
          if(eq.const){
            oldv2 <- paste("theta[", match(partable$rhs[full.idx], ov.names), ",", grp.idx, "]", sep="")
          }
        } else {
          partable$lhs[tmprows[2]] <- partable$rhs[covpars[i]]
          partable$op[tmprows[2]] <- "~"
          partable$rhs[tmprows[2]] <- phname
          partable$mat[tmprows[2]] <- "beta"
          tbcs <- tbcs + 1
          partable$row[tmprows[2]] <- tbcs
          partable$col[tmprows[2]] <- match(partable$rhs[covpars[i]],
                                            patts$mmDimNames[[k]]$psi[[1]])
          tmpv2 <- paste("psi[", match(partable$rhs[covpars[i]], patts$mmDimNames[[k]]$psi[[1]]), ",", k, "]", sep="")
          if(eq.const){
            oldv2 <- paste("psi[", match(partable$rhs[full.idx], patts$mmDimNames[[k]]$psi[[1]]), ",", grp.idx, "]", sep="")
          }
        }

        ## Decide what priors to use
        ttcs <- ttcs + 1
        if((ctype == "ov" & ov.cp == "srs") | (ctype == "lv" & lv.cp == "srs")){
          ## srs priors
          partable$free[tmprows[1:3]] <- 0
          partable$exo[tmprows[1:3]] <- 0
          partable$ustart[tmprows[1]] <- paste("sqrt(abs(rho[", rhoind, ",", k, "])*", tmpv1, ")", sep="")
          partable$ustart[tmprows[2]] <- paste("(-1 + 2*step(rho[", rhoind, ",", k,
                                             "]))*sqrt(abs(rho[", rhoind, ",", k, "])*", tmpv2, ")", sep="")
          partable$ustart[tmprows[3]] <- 1
          partable$mat[tmprows[3]] <- "theta"
          partable$row[tmprows[3]] <- partable$col[tmprows[3]] <- ttcs

          partable$id[covparg] <- paste("rho[", rhoind, ",", k, "]", sep="")
          partable$plabel[tmprows] <- paste(".p", tmprows,
                                                   ".", sep="")

          partable$lhs[tmprows[3]] <- partable$rhs[tmprows[3]] <- phname
          partable$op[tmprows[3]] <- "~~"

          ## add to rho index
          rhoind <- rhoind + 1
        } else {
          ## factor analysis priors
          partable$mat[tmprows[3]] <- "theta"
          partable$row[tmprows[3]] <- partable$col[tmprows[3]] <- ttcs
          if(partable$free[covparg] == 0){
            if(partable$ustart[covparg] != 0) stop("blavaan ERROR: Cannot fix covariances to nonzero values under fa priors.\n")
            partable$free[tmprows[1:3]] <- 0
            partable$ustart[tmprows[3]] <- 1
            partable$ustart[tmprows[1:2]] <- 0
          } else {
            partable$plabel[tmprows[3]] <- partable$plabel[covparg]
            partable$free[tmprows[3]] <- partable$free[covparg]
            partable$exo[tmprows[1:3]] <- 0
            partable$plabel[tmprows[1:2]] <- paste(".p", tmprows[1:2], ".", sep="")
            partable$free[tmprows[1:2]] <- tmprows[1:2]
          }
          ## TODO tmprows[3] should get some prior associated with covariance param.
          partable$prior[tmprows[1:3]] <- ""
          cprm <- c(cprm, covparg)

          partable$lhs[tmprows[3]] <- partable$rhs[tmprows[3]] <- phname
          partable$op[tmprows[3]] <- "~~"
        }
        ## equality constraints
        if(eq.const){
          ## set loading/regression parameters equal to others
          partable$free[tmprows[1:2]] <- tmprows[1:2]
          old.labels <- which(partable$op %in% c("=~", "~") &
                              partable$group == grp.idx &
                              (grepl(paste(".phant", old.ridx, sep=""),
                               partable$lhs) |
                               grepl(paste(".phant", old.ridx, sep=""),
                               partable$rhs)))

          partable <- rbind(partable, blkrow, blkrow)
          nr <- nrow(partable)
          partable$lhs[(nr-1):nr] <- partable$plabel[old.labels[1:2]]
          partable$op[(nr-1):nr] <- "=="
          partable$rhs[(nr-1):nr] <- partable$plabel[tmprows[1:2]]
          partable$user[(nr-1):nr] <- 2
          partable$free[(nr-1):nr] <- 0
          partable$group[(nr-1):nr] <- 0

          if((ctype == "ov" & ov.cp == "fa") | (ctype == "lv" & lv.cp == "fa")){
            old.label <- which(partable$op == "~~" &
                               partable$group == grp.idx &
                               grepl(paste(".phant", old.ridx, sep=""), partable$lhs))

            partable <- rbind(partable, blkrow)
            nr <- nrow(partable)
            partable$lhs[nr] <- partable$plabel[old.label]
            partable$op[nr] <- "=="
            partable$rhs[nr] <- partable$plabel[tmprows[3]]
            partable$user[nr] <- 2
            partable$free[nr] <- 0
            partable$group[nr] <- 0
          }
        }
      }
    }
    ## Now remove rows of fa covariance parameters
    if(!is.null(cprm)) partable <- partable[-cprm,]
  }

  ## FIXME?
  ## Remove covariances associated with fixed x
  covpars <- which(partable$op == "~~" &
                   partable$lhs != partable$rhs &
                   partable$group == 1 &
                   partable$lhs %in% ov.names.x &
                   partable$free == 0)

  if(length(covpars) > 0) partable <- partable[-covpars,]

  list(partable = partable, parvec = parvec)
}

set_mv0 <- function(partable, ov.names, ngroups) {
    ## If any mvs have fixed 0 variance (single indicator lv),
    ## move the fixed 0 variance to the lv
    mv0 <- which(partable$op == "~~" &
                 partable$lhs %in% ov.names &
                 partable$rhs == partable$lhs &
                 partable$group == 1 &
                 partable$free == 0 &
                 partable$ustart == 0)

    if(length(mv0) > 0){
        ovn <- partable$lhs[mv0]
    
        for(i in 1:length(ovn)){
            for(j in 1:ngroups){
                mvloc <- which(partable$op == "~~" &
                               partable$lhs == ovn[i] &
                               partable$rhs == partable$lhs &
                               partable$group == j &
                               partable$free == 0 &
                               partable$ustart == 0)
                
                lvloc <- which(partable$op == "=~" &
                               partable$rhs == ovn[i] &
                               partable$group == j)

                ## This should always be length 1 because, if it is also
                ## an indicator of another lv (with other mvs), then
                ## we can estimate the variance
                if(length(lvloc) > 1){
                    stop(paste("blavaan ERROR: Problem with", ovn[i],
                               "as a single indicator of a latent variable.\n"))
                }

                lvname <- partable$lhs[lvloc]

                ## TODO? check for covariances with the lv?
                lvvar <- which(partable$lhs == lvname &
                               partable$rhs == lvname &
                               partable$op == "~~" &
                               partable$group == j)

                tmpfree <- partable$free[lvvar]
                tmpustart <- partable$ustart[lvvar]
                tmpplabel <- partable$plabel[lvvar]
                tmpstart <- partable$start[lvvar]

                partable$free[lvvar] <- partable$free[mvloc]
                partable$ustart[lvvar] <- partable$ustart[mvloc]
                partable$plabel[lvvar] <- partable$plabel[mvloc]
                partable$start[lvvar] <- partable$plabel[mvloc]

                partable$free[mvloc] <- tmpfree
                partable$ustart[mvloc] <- tmpustart
                partable$plabel[mvloc] <- tmpplabel
                partable$start[mvloc] <- tmpstart
            }
        }
    }

    partable
}
