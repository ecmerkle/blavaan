set_stancovs <- function(partable, ov.names, ov.names.x, dp) {
  ## Add phantom lvs for covariance parameters

  ## first: parameter matrices + indexing
  partable <- lavMatrixRepresentation(partable, add.attributes = TRUE)
  ## for defined parameters
  defpar <- which(partable$op == ":=")
  if(length(defpar) > 0){
    partable$mat[defpar] <- "def"
    partable$row[defpar] <- 1:length(defpar)
    partable$col[defpar] <- 1
    partable$group[defpar] <- 1
  }

  ## add prior column if it doesn't exist
  if(is.na(match("prior", names(partable)))) partable$prior <- rep("", length(partable$id))
  
  covpars <- which(partable$op == "~~" &
                   partable$lhs != partable$rhs &
                   !(partable$lhs %in% ov.names.x &
                     partable$free == 0))

  blkrow <- rep(NA, length(partable$id))
  partable$rhoidx <- rep(NA, length(partable$id))

  ## Only do this if covpars exist
  if(length(covpars) > 0){
    ## add to model matrices
    nmvcovs <- sum(partable$lhs[covpars] %in% ov.names)
    nlvcovs <- length(covpars) - nmvcovs

    mvcov <- 0
    lvcov <- 0
    
    for(i in 1:length(covpars)){      
      ## Is this constrained equal to a previous parameter?
      eq.const <- FALSE
      eq.idx <- which(partable$op == "==" & partable$rhs == partable$plabel[covpars[i]])
      if(length(eq.idx) > 0){
        eq.const <- TRUE
        ## TODO? assumes it is equal to another covariance; do any models
        ## restrict covariances to be equal to other types of parameters?
        full.idx <- which(partable$plabel == partable$lhs[eq.idx])
        old.idx <- partable$rhoidx[full.idx]
      }
      
      tmprows <- nrow(partable) + 1
      partable <- rbind(partable, blkrow)

      ## TODO? should 'block' ever differ from 'group'?
      partable$group[tmprows] <- partable$block[tmprows] <-
        partable$group[covpars[i]]

      partable$lhs[tmprows] <- partable$lhs[covpars[i]]
      partable$rhs[tmprows] <- partable$rhs[covpars[i]]
      
      ## Decide on =~ (ov) vs ~ (lv)
      if(partable$lhs[covpars[i]] %in% ov.names){
        if(!eq.const){
          mvcov <- mvcov + 1
          covidx <- mvcov
        }
        partable$mat[tmprows] <- "rho"
      } else {
        if(!eq.const){
          lvcov <- lvcov + 1
          covidx <- lvcov
        }
        partable$mat[tmprows] <- "lvrho"
      }
      partable$op[tmprows] <- "~~"
      partable$row[tmprows] <- partable$row[covpars[i]]
      partable$col[tmprows] <- partable$col[covpars[i]]
      partable$group[tmprows] <- partable$group[covpars[i]]

      v1var <- which(partable$lhs == partable$lhs[covpars[i]] &
                     partable$rhs == partable$lhs[covpars[i]] &
                     partable$group == partable$group[covpars[i]] &
                     partable$op == "~~")
      tmpv1 <- paste(partable$mat[v1var], "[", partable$row[v1var], ",", partable$col[v1var], ",", partable$group[v1var],
                         "]", sep="")
      
      v2var <- which(partable$lhs == partable$rhs[covpars[i]] &
                     partable$rhs == partable$rhs[covpars[i]] &
                     partable$group == partable$group[covpars[i]] &
                     partable$op == "~~")
      tmpv2 <- paste(partable$mat[v2var], "[", partable$row[v2var], ",", partable$col[v2var], ",", partable$group[v2var], "]", sep="")

      if(partable$prior[covpars[i]] != ""){
        partable$prior[tmprows] <- partable$prior[covpars[i]]
      } else {
        partable$prior[tmprows] <- ""
      }

      if(eq.const){
        partable$ustart[covpars[i]] <- paste0(partable$mat[full.idx],
                                              "[",
                                              partable$row[full.idx],
                                              ",", partable$col[full.idx],
                                              ",", partable$group[full.idx],
                                              "]")
        partable$ustart[tmprows] <- paste0(partable$ustart[covpars[i]], "/sqrt(", tmpv1, "*", tmpv2, ")")
      } else {
        partable$rhoidx[tmprows] <- partable$rhoidx[covpars[i]] <- covidx
        partable$ustart[covpars[i]] <- paste0(partable$mat[tmprows],
                                              "[", partable$row[tmprows],
                                              ",", partable$col[tmprows],
                                              ",", partable$group[tmprows],
                                              "] * sqrt(", tmpv1,
                                              " * ", tmpv2, ")")
      }
      partable$free[tmprows] <- partable$free[covpars[i]]
      partable$free[covpars[i]] <- 0
      partable$plabel[tmprows] <- paste(".p", tmprows, ".", sep="")
      partable$label[tmprows] <- ""
    }

    ## put covariances last, so that they appear last in
    ## the defined parameter block (they are functions of
    ## other parameters)
    ptcov <- partable[covpars,]
    partable <- partable[-covpars,]
    partable <- rbind(partable, ptcov)
  }

  ## FIXME?
  ## Remove covariances associated with fixed x
  covpars <- which(partable$op == "~~" &
                   partable$lhs != partable$rhs &
                   partable$group == 1 &
                   partable$lhs %in% ov.names.x &
                   partable$free == 0)

  if(length(covpars) > 0) partable <- partable[-covpars,]

  partable
}
