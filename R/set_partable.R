set_phantoms <- function(partable, ov.names, lv.names, ov.names.x, lv.names.x, ov.cp, lv.cp, lv.x.wish, ngroups) {
  ## Add phantom lvs for covariance parameters

  ## first: parameter matrices + indexing
  partable <- lavMatrixRepresentation(partable, add.attributes = TRUE)

  ## add prior column if it doesn't exist
  if(is.na(match("prior", names(partable)))) partable$prior <- rep("", length(partable$id))
    
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
  facovs <- NULL

  ## Only do this if covpars exist
  if(length(covpars) > 0){
    ## add to model matrices
    ## added entries in lambda vs in beta
    nmvcovs <- sum(partable$lhs[covpars] %in% ov.names)
    nlvcovs <- length(covpars) - nmvcovs
    patts <- attributes(partable)
    for(k in 1:ngroups){
      if(!("lambda" %in% patts$mmNames[[k]]) & nmvcovs > 0){
        lcolstart <- 0
        attributes(partable)$mmNames[[k]] <- c(patts$mmNames[[k]],
                                             "lambda")
        attributes(partable)$mmRows[[k]] <- c(patts$mmRows[[k]], lambda=length(ov.names))
        attributes(partable)$mmCols[[k]] <- c(patts$mmCols[[k]], lambda=nmvcovs)
      } else {
        lcolstart <- patts$mmCols[[k]]["lambda"]
        attributes(partable)$mmCols[[k]]["lambda"] <- patts$mmCols[[k]]["lambda"] + nmvcovs
      }
      if(!("beta" %in% patts$mmNames[[k]]) & nlvcovs > 0){
        bcolstart <- 0
        attributes(partable)$mmNames[[k]] <- c(patts$mmNames[[k]], "beta")
        attributes(partable)$mmRows[[k]] <- c(patts$mmRows[[k]], beta=
nlvcovs)
        attributes(partable)$mmCols[[k]] <- c(patts$mmCols[[k]], beta=nlvcovs)
      } else {
        bcolstart <- patts$mmCols[[k]]["beta"]
        attributes(partable)$mmRows[[k]]["beta"] <- patts$mmRows[[k]]["beta"] + nlvcovs
        attributes(partable)$mmCols[[k]]["beta"] <- patts$mmCols[[k]]["beta"] + nlvcovs
      }

      if(!("psi" %in% patts$mmNames[[k]])){
        psicolstart <- 0
        attributes(partable)$mmNames[[k]] <- c(patts$mmNames[[k]], "psi")
        attributes(partable)$mmRows[[k]] <- c(patts$mmRows[[k]], psi=length(covpars))
        attributes(partable)$mmCols[[k]] <- c(patts$mmCols[[k]], psi=length(covpars))
      } else {
        psicolstart <- patts$mmCols[[k]]["psi"]
        attributes(partable)$mmRows[[k]]["psi"] <- patts$mmRows[[k]]["psi"] + length(covpars)
        attributes(partable)$mmCols[[k]]["psi"] <- patts$mmCols[[k]]["psi"] + length(covpars)
      }
    }
    
    cprm <- NULL

    ## which covariances are under srs?
    ridx <- 1:length(covpars)

    for(k in 1:ngroups){
      tlcs <- lcolstart
      tbcs <- bcolstart
      tpcs <- psicolstart
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
          v1var <- which(partable$lhs == partable$lhs[covpars[i]] &
                         partable$rhs == partable$lhs[covpars[i]] &
                         partable$group == k &
                         partable$op == "~~")
          tmpv1 <- paste(partable$mat[v1var], "[", partable$row[v1var], ",", partable$col[v1var], ",", k,
                         "]", sep="")
          if(eq.const){
            oldr <- match(partable$lhs[full.idx], patts$mmDimNames[[k]]$lambda[[1]])
            oldv1 <- paste(partable$mat[v1var], "[", oldr , ",", oldr, ",", grp.idx, "]", sep="")
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
          tmpv1 <- paste("psi[", partable$col[tmprows[1]], ",", partable$col[tmprows[1]],
                                       ",", k, "]", sep="")
          if(eq.const){
            oldr <- match(partable$lhs[full.idx], patts$mmDimNames[[k]]$psi[[1]])
            oldv1 <- paste("psi[", oldr, ",", oldr, ",", grp.idx, "]", sep="")
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
          v2var <- which(partable$lhs == partable$rhs[covpars[i]] &
                         partable$rhs == partable$rhs[covpars[i]] &
                         partable$group == k &
                         partable$op == "~~")
          tmpv2 <- paste(partable$mat[v2var], "[", partable$row[v2var], ",", partable$col[v2var], ",", k, "]", sep="")
          if(eq.const){
            oldr <- match(partable$rhs[full.idx], patts$mmDimNames[[k]]$lambda[[1]])
            oldv2 <- paste(partable$mat[v2var], "[", oldr, ",", oldr, ",", grp.idx, "]", sep="")
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
          tmpv2 <- paste("psi[", partable$col[tmprows[2]], ",", partable$col[tmprows[2]],
                                       ",", k, "]", sep="")
          if(eq.const){
            oldr <- match(partable$lhs[full.idx], patts$mmDimNames[[k]]$psi[[1]])
            oldv2 <- paste("psi[", oldr, ",", oldr, ",", grp.idx, "]", sep="")
          }
        }

        ## TODO tmprows[3] should get some prior associated with covariance param for fa priors?
        partable$prior[tmprows[1:3]] <- ""
        
        ## Decide what priors to use
        tpcs <- tpcs + 1
        if((ctype == "ov" & ov.cp == "srs") | (ctype == "lv" & lv.cp == "srs")){
          rhomat <- "rho"
          if(partable$mat[covpars[i]] == "psi") rhomat <- "lvrho"
          rhoind <- paste(partable$row[covpars[i]], ",", partable$col[covpars[i]], sep="")
          ## srs priors
          partable$free[tmprows[1:3]] <- 0
          partable$exo[tmprows[1:3]] <- 0
          partable$ustart[tmprows[1]] <- paste("sqrt(abs(", rhomat, "[", rhoind, ",", k, "])*", tmpv1, ")", sep="")
          partable$ustart[tmprows[2]] <- paste("(-1 + 2*step(", rhomat, "[", rhoind, ",", k,
                                             "]))*sqrt(abs(", rhomat, "[", rhoind, ",", k, "])*", tmpv2, ")", sep="")
          partable$ustart[tmprows[3]] <- 1
          partable$mat[tmprows[3]] <- "psi"
          partable$row[tmprows[3]] <- partable$col[tmprows[3]] <- tpcs

          partable$id[covparg] <- paste(rhomat, "[", rhoind, ",", k, "]", sep="")
          partable$plabel[tmprows] <- paste(".p", tmprows,
                                                   ".", sep="")

          partable$lhs[tmprows[3]] <- partable$rhs[tmprows[3]] <- phname
          partable$op[tmprows[3]] <- "~~"
        } else {
          ## factor analysis priors
          partable$mat[tmprows[3]] <- "psi"
          partable$row[tmprows[3]] <- partable$col[tmprows[3]] <- tpcs
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
          partable$mat[(nr-1):nr] <- ""
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
            partable$mat[(nr-1):nr] <- ""
            partable$user[nr] <- 2
            partable$free[nr] <- 0
            partable$group[nr] <- 0
          }
        }
      }
    }
    ## Now remove rows of fa covariance parameters
    if(!is.null(cprm)){
      facovs <- partable[cprm,]
      partable <- partable[-cprm,]
    }
  }

  ## FIXME?
  ## Remove covariances associated with fixed x
  covpars <- which(partable$op == "~~" &
                   partable$lhs != partable$rhs &
                   partable$group == 1 &
                   partable$lhs %in% ov.names.x &
                   partable$free == 0)

  if(length(covpars) > 0) partable <- partable[-covpars,]

  ## Add parameter numbers now that we have phantoms
  parnums <- rep(NA, nrow(partable))
  parrows <- which(partable$op != "==")
  parnums[parrows] <- 1:length(parrows)
  partable$parnums <- parnums    
    
  list(partable = partable, facovs = facovs)
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

set_phanvars <- function(partable, ov.names, lv.names, ov.cp, lv.cp, ngroups){
    ## once we have defined all phantoms, go back to set equality
    ## constraints on phantom variances

    vnames <- c(ov.names, lv.names)
    for(i in 1:length(vnames)){
        for(k in 1:ngroups){
            if(vnames[i] %in% ov.names){
                phanlvs <- which(partable$rhs == vnames[i] &
                                 grepl(".phant", partable$lhs) &
                                 partable$op == "=~" &
                                 partable$group == k)
            } else {
                phanlvs <- which(partable$lhs == vnames[i] &
                                 grepl(".phant", partable$rhs) &
                                 partable$op == "~" &
                                 partable$group == k)
            }
            if(length(phanlvs) > 0){
                vvar <- which(partable$lhs == vnames[i] &
                              partable$lhs == partable$rhs &
                              partable$op == "~~" &
                              partable$group == k)

                if(ov.cp == "srs"){
                    eqconst <- paste(partable$mat[vvar], "[", partable$row[vvar],
                                     ",", partable$col[vvar], ",", k, "]", sep="")
                    for(j in 1:length(phanlvs)){
                        eqconst <- paste(eqconst, " - (", partable$ustart[phanlvs[j]],
                                         ")^2", sep="")
                    }
                    partable <- rbind(partable, partable[vvar,])
                    partable$parnums[nrow(partable)] <- max(partable$parnums, na.rm=TRUE) + 1
                    partable$mat[vvar] <- paste(partable$mat[vvar], "star", sep="")
                } else {
                    ## fa priors TODO change around so it is faster
                    eqconst <- paste(partable$mat[vvar], "star[", partable$row[vvar],
                                     ",", partable$col[vvar], ",", k, "]", sep="")
                    for(j in 1:length(phanlvs)){
                        if(vnames[i] %in% ov.names){
                          phname <- partable$lhs[phanlvs[j]]
                        } else {
                          phname <- partable$rhs[phanlvs[j]]
                        }
                        phanvar <- which(partable$lhs == phname &
                                     partable$lhs == partable$rhs &
                                     partable$op == "~~" &
                                     partable$group == k)
                        eqconst <- paste(eqconst, " + (", partable$mat[phanlvs[j]], "[",
                                         partable$row[phanlvs[j]], ",",
                                         partable$col[phanlvs[j]], ",", k, "]^2*",
                                         partable$mat[phanvar], "[",
                                         partable$row[phanvar], ",", partable$col[phanvar],
                                         ",", k, "])", sep="")
                    }
                    partable <- rbind(partable, partable[vvar,])
                    partable$mat[nrow(partable)] <- paste(partable$mat[vvar], "star",
                                                          sep="")
                    partable$parnums[nrow(partable)] <- max(partable$parnums, na.rm=TRUE) + 1
                }
                partable$free[vvar] <- 0
                partable$ustart[vvar] <- eqconst
            }
        }
    }
    partable
}
