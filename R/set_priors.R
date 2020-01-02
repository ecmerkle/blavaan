set_parvec <- function(TXT2, partable, dp, cp, lv.x.wish, lv.names.x, target="jags"){
    ## tabs
    t1 <- paste(rep(" ", 2L), collapse="")
    t2 <- paste(rep(" ", 4L), collapse="")
    t3 <- paste(rep(" ", 6L), collapse="")

    eqop <- ifelse(target == "stan", "=", "<-")
    commop <- ifelse(target == "stan", "// ", "# ")
    eolop <- ifelse(target == "stan", ";", "")
  
    ## parameter assignments separate from priors
    TXT3 <- paste("\n", t1, commop, "Priors", sep="")

    ## find parameters with wishart priors
    wishpars <- NULL
    if(lv.x.wish & length(lv.names.x) > 1){
      wishpars <- which(partable$lhs %in% lv.names.x &
                        partable$rhs %in% lv.names.x &
                        partable$op == "~~")
    }

    ## parameter numbers that need priors
    partable$freeparnums <- rep(0, length(partable$parnums))
    parvecnum <- 0
    for(i in 1:nrow(partable)){
        miscignore <- (partable$mat[i] == "") | (i %in% wishpars)

        eqpar <- which(partable$rhs == partable$plabel[i] &
                       partable$op == "==")
        defeq <- partable$op[i] %in% c("==", ":=") &
                 grepl("\\+|-|/|\\*|\\(|\\)|\\^", partable$rhs[i])
        compeq <- which(partable$lhs == partable$plabel[i] &
                        partable$op %in% c("==", ":=") &
                        grepl("\\+|-|/|\\*|\\(|\\)|\\^", partable$rhs))
        fixed <- partable$free[i] == 0 & partable$op[i] != ":="
        if(length(eqpar) > 0 | defeq | length(compeq) > 0 | fixed |
           miscignore){
            next
        } else {
            parvecnum <- parvecnum + 1
            partable$freeparnums[i] <- parvecnum
        }
    }

    for(i in 1:nrow(partable)){
        if((partable$mat[i] != "" & !(i %in% wishpars)) | partable$op[i] == ":="){            
            ## to find equality constraints
            eqpar <- which(partable$rhs == partable$plabel[i] &
                           partable$op == "==")

            ## only complex equality constraints and defined parameters;
            ## rhs needs math expression
            defeq <- partable$op[i] %in% c("==", ":=") &
                     grepl("\\+|-|/|\\*|\\(|\\)|\\^", partable$rhs[i])
            compeq <- which(partable$lhs == partable$plabel[i] &
                            partable$op %in% c("==", ":=") &
                            grepl("\\+|-|/|\\*|\\(|\\)|\\^", partable$rhs))
            ## TODO block prior associated with lv.x.wish
            ##      put entries of parvec in matrix for dwish?
            ## TODO check for inequality constraints here?

            ## correlation parameter under srs
            if(grepl("rho", partable$id[i])){
                rhoinf <- strsplit(partable$id[i], "[, \\[^\\]]+", perl=TRUE)
                partable$mat[i] <- rhoinf[[1]][1]
                partable$row[i] <- rhoinf[[1]][2]
                partable$col[i] <- rhoinf[[1]][3]
            }
          
            ## start parameter assignment
            TXT2 <- paste(TXT2, "\n", t1, partable$mat[i], "[",
                          partable$row[i], ",", partable$col[i],
                          ",", partable$group[i], "] ", eqop,
                          " ", sep="")
            if(grepl("rho", partable$id[i]) & partable$free[i] > 0){
              TXT2 <- paste(TXT2, "-1 + 2*", sep="")
            }
          
            if(partable$free[i] == 0 & partable$op[i] != ":="){
                if(is.na(partable$ustart[i])){
                    ## exo
                    TXT2 <- paste(TXT2, partable$start[i], eolop,
                                  sep="")
                } else {
                    TXT2 <- paste(TXT2, partable$ustart[i], eolop,
                                  sep="")
                }
            } else if(length(eqpar) > 0){
                eqpar <- which(partable$plabel == partable$lhs[eqpar])
                ## in case it is an "expanded" variance
                if(length(eqpar) > 1){
                    if(length(eqpar) > 2) stop("blavaan ERROR: problem with parameter equality constraints")
                    eqpar <- eqpar[which(partable$freeparnums[eqpar] > 0)]
                }
                if(partable$freeparnums[eqpar] == 0){
                    eqtxt <- paste(partable$mat[eqpar], "[",
                                   partable$row[eqpar], ",",
                                   partable$col[eqpar], ",",
                                   partable$group[eqpar], "]",
                                   eolop, sep="")
                } else {
                    eqtxt <- paste("parvec[",
                                   partable$freeparnums[eqpar],
                                   "]", eolop, sep="")
                }

                vpri <- grepl("\\[var\\]", partable$prior[eqpar])
                spri <- grepl("\\[sd\\]", partable$prior[eqpar])
                if(!vpri & (grepl("theta", partable$mat[i]) | grepl("psi", partable$mat[i]))){
                    sq <- ifelse(spri, "2", "-1")
                    TXT2 <- paste(TXT2, "pow(", eqtxt, ",", sq,
                                  ")", eolop, sep="")
                } else {
                    TXT2 <- paste(TXT2, eqtxt, sep="")
                }
            } else if(defeq | length(compeq) > 0){
                if(length(compeq) == 0) compeq <- i
                ## constraints with one parameter label on lhs
                ## FIXME? cannot handle, e.g., b1 + b2 == 2
                ## see lav_partable_constraints.R
                rhsvars <- all.vars(parse(file="",
                                          text=partable$rhs[compeq]))
                if(compeq == i){
                    pvnum <- match(rhsvars, partable$label)
                } else {
                    pvnum <- match(rhsvars, partable$plabel)
                }

                rhstrans <- paste(partable$mat[pvnum], "[",
                                  partable$row[pvnum], ",",
                                  partable$col[pvnum], ",",
                                  partable$group[pvnum], "]",
                                  sep="")

                oldjageq <- partable$rhs[compeq]
                transtab <- as.list(rhstrans)
                names(transtab) <- rhsvars
                jagexpr <- parse(text=oldjageq)[[1]]
                jageq <- do.call("substitute", list(jagexpr,
                                                    transtab))
                jageq <- paste(deparse(jageq, width.cutoff = 500), collapse="")

                jageq <- gsub('\"', '', jageq)

                TXT2 <- paste(TXT2, jageq, eolop, sep="")
            } else {
                ## needs a prior
                TXT3 <- paste(TXT3, "\n", t1, "parvec[",
                              partable$freeparnums[i], "]", sep="")
                if(partable$prior[i] == ""){
                    if(partable$mat[i] == "lvrho"){
                        partype <- grep("rho", names(dp))
                    } else if(grepl("star", partable$mat[i])){
                        pname <- paste("i", strsplit(partable$mat[i], "star")[[1]][1], sep="")
                        partype <- grep(pname, names(dp))
                    } else {
                        partype <- grep(partable$mat[i], names(dp))
                    }
                    if(length(partype) > 1) partype <- partype[1] # due to psi and ibpsi

                    partable$prior[i] <- dp[partype]
                }
                jagpri <- strsplit(partable$prior[i], "\\[")[[1]][1]
                vpri <- grepl("\\[var\\]", partable$prior[i])
                spri <- grepl("\\[sd\\]", partable$prior[i])
                if(!vpri & (grepl("theta", partable$mat[i]) | grepl("psi", partable$mat[i]))){
                    sq <- ifelse(spri, "2", "-1")
                    TXT2 <- paste(TXT2, "pow(parvec[",
                                  partable$freeparnums[i], "],", sq,
                                  ")", eolop, sep="")
                } else {
                    TXT2 <- paste(TXT2, "parvec[",
                                  partable$freeparnums[i],
                                  "]", eolop, sep="")
                }
                TXT3 <- paste(TXT3, " ~ ", jagpri, eolop, sep="")
            }
        }
    }

    ## deal with wishart priors
    if(lv.x.wish & length(lv.names.x) > 1){
      nlvx <- length(lv.names.x)
      ngroups <- max(partable$group, na.rm = TRUE)

      TXT3 <- paste(TXT3, "\n", t1, "for(k in 1:", ngroups,
                    ") {\n", t2, "ibpsi[1:", nlvx, ",1:", nlvx,
                    ",k] ~ dwish(iden,", nlvx+1, ")\n", sep="")
      
      TXT3 <- paste(TXT3, t2, "bpsi[1:", nlvx, ",1:", nlvx, ",k] <- inverse(ibpsi[1:",
                    nlvx, ",1:", nlvx, ",k])\n", t1, "}\n", sep="")

      for(i in 1:length(wishpars)){
        tmppar <- wishpars[i]
        wishrow <- which(lv.names.x == partable$lhs[tmppar])
        wishcol <- which(lv.names.x == partable$rhs[tmppar])

        partable$prior[tmppar] <- dp[["ibpsi"]]
        
        TXT2 <- paste(TXT2, "\n", t1, partable$mat[tmppar], "[",
                      partable$row[tmppar], ",", partable$col[tmppar],
                      ",", partable$group[tmppar], "] ", eqop,
                      " bpsi[", wishrow, ",", wishcol, ",",
                      partable$group[tmppar], "]", eolop, sep="")
      }
    }
      
    ## now define inferential covariances and priors for inferential
    ## variances, if needed
    covs <- unique(partable$lhs[grep(".phant", partable$lhs)])
    
    if(length(covs) > 0){
      TXT2 <- paste(TXT2, "\n\n", t1, commop, "Inferential covariances", sep="")
        for(i in 1:length(covs)){
            for(k in 1:max(partable$group)){
                varlocs <- which(((partable$lhs == covs[i] &
                                   partable$op == "=~") |
                                  (partable$rhs == covs[i] &
                                   partable$op == "~")) &
                                 partable$group == k)
                vartxt <- "star"
                vars <- partable$rhs[varlocs]
                ## catch where we need lhs instead of rhs
                lhsvars <- grepl(".phant", vars)
                if(any(lhsvars)){
                    vars[lhsvars] <- partable$lhs[varlocs[lhsvars]]
                }

                if(length(varlocs) == 0){
                    ## lv
                    varlocs <- which(partable$rhs == covs[i] &
                                     partable$op == "~" &
                                     partable$group == k)
                    vars <- partable$lhs[varlocs]
                }

                var1 <- which(partable$lhs == vars[1] &
                              partable$lhs == partable$rhs &
                              partable$group == partable$group[varlocs[1]] &
                              grepl(vartxt, partable$mat))
                var2 <- which(partable$lhs == vars[2] &
                              partable$lhs == partable$rhs &
                              partable$group == partable$group[varlocs[1]] &
                              grepl(vartxt, partable$mat))

                matname <- ifelse(grepl("theta", partable$mat[var1]), "theta", "psi")
                phpars <- which((partable$lhs == covs[i] |
                                 partable$rhs == covs[i]) &
                                partable$group == k)
                if(length(phpars) == 1){
                    phpars <- which(partable$rhs == covs[i] &
                                    partable$group == k)
                }

                ## covariances
                TXT2 <- paste(TXT2, "\n", t1, matname, "[", partable$row[var1],
                              ",", partable$row[var2], ",", partable$group[varlocs[1]], "] ", eqop, " ",
                              partable$mat[phpars[1]], "[", partable$row[phpars[1]], ",",
                              partable$col[phpars[1]], ",", partable$group[phpars[1]], "]*",
                              partable$mat[phpars[2]], "[", partable$row[phpars[2]], ",",
                              partable$col[phpars[2]], ",", partable$group[phpars[2]], "]*",
                              partable$mat[phpars[3]], "[", partable$row[phpars[3]], ",",
                              partable$col[phpars[3]], ",", partable$group[phpars[3]], "]",
                              eolop, sep="")
            }

        }
        TXT2 <- paste(TXT2, "\n", sep="")
    }
    list(TXT2 = TXT2, TXT3 = TXT3, partable = partable)
}
