set_parvec <- function(TXT2, partable, dp, cp, lv.x.wish, lv.names.x){
    ## tabs
    t1 <- paste(rep(" ", 2L), collapse="")
    t2 <- paste(rep(" ", 4L), collapse="")
    t3 <- paste(rep(" ", 6L), collapse="")

    ## parameter assignments separate from priors
    TXT3 <- paste("\n\n", t1, "# Priors/constraints", sep="")

    ## find parameters with wishart priors
    wishpars <- NULL
    if(lv.x.wish & length(lv.names.x) > 1){
      wishpars <- which(partable$lhs %in% lv.names.x &
                        partable$rhs %in% lv.names.x &
                        partable$op == "~~")
    }

    for(i in 1:nrow(partable)){
        if(partable$mat[i] != "" & !(i %in% wishpars)){
            ## to find equality constraints
            eqpar <- which(partable$rhs == partable$plabel[i] &
                           partable$op == "==")
            ## only complex equality constraints; rhs needs math expression
            compeq <- which(partable$lhs == partable$label[i] &
                            partable$op == "==" &
                            grepl("\\+|-|/|\\*|\\(|\\)|\\^", partable$rhs))
            ## TODO block prior associated with lv.x.wish
            ##      put entries of parvec in matrix for dwish?
            ## TODO check for inequality constraints here?

            TXT3 <- paste(TXT3, "\n", t1, "parvec[",
                          partable$parnums[i], "]", sep="")

            ## correlation parameter under srs
            if(grepl("rho", partable$id[i])){
                rhoinf <- strsplit(partable$id[i], "[, \\[^\\]]+", perl=TRUE)
                partable$mat[i] <- rhoinf[[1]][1]
                partable$row[i] <- rhoinf[[1]][2]
                partable$col[i] <- rhoinf[[1]][3]
                if(partable$free[i] == 0){
                    partable$ustart[i] <- (as.numeric(partable$ustart[i]) + 1)/2
                }
            }
          
            if(partable$free[i] == 0){
                TXT3 <- paste(TXT3, " <- ", sep="")
                if(is.na(partable$ustart[i])){
                    ## exo
                    TXT3 <- paste(TXT3, partable$start[i], sep="")
                } else {
                    TXT3 <- paste(TXT3, partable$ustart[i], sep="")
                }
            } else if(length(eqpar) > 0){
                eqpar <- which(partable$plabel == partable$lhs[eqpar])
                TXT3 <- paste(TXT3, " <- parvec[", partable$parnums[eqpar],
                              "]", sep="")
            } else if(length(compeq) > 0){
                ## constraints with one parameter label on lhs
                ## FIXME? cannot handle, e.g., b1 + b2 == 2
                ## see lav_partable_constraints.R
                rhsvars <- all.vars(parse(file="",
                                          text=partable$rhs[compeq]))
                pvnum <- match(rhsvars, partable$label)

                rhstrans <- paste("parvec[", partable$parnums[pvnum], "]",
                                  sep="")

                jageq <- partable$rhs[compeq]
                for(j in 1:length(rhsvars)){
                    jageq <- gsub(rhsvars[j], rhstrans[j], jageq)
                }

                TXT3 <- paste(TXT3, " <- ", jageq, sep="")
            } else {
                ## needs a prior
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

                vpri <- grepl("\\[var\\]", partable$prior[i])
                spri <- grepl("\\[sd\\]", partable$prior[i])
                jagpri <- strsplit(partable$prior[i], "\\[")[[1]][1]
                if(!vpri & (grepl("theta", partable$mat[i]) | grepl("psi", partable$mat[i]))){
                    sq <- ifelse(spri, "2", "-1")
                    TXT3 <- paste(TXT3, " <- pow(pvec", partable$parnums[i], ",", sq,
                                  ")\n", sep="")
                    TXT3 <- paste(TXT3, t1, "pvec", partable$parnums[i], " ~ ",
                                  jagpri, sep="")
                } else {
                    TXT3 <- paste(TXT3, " ~ ", jagpri, sep="")
                }
            }

            TXT2 <- paste(TXT2, "\n", t1, partable$mat[i], "[",
                          partable$row[i], ",", partable$col[i],
                          ",", partable$group[i], "] <- ", sep="")
            if(grepl("rho", partable$id[i])) TXT2 <- paste(TXT2, "-1 + 2*", sep="")
            TXT2 <- paste(TXT2, "parvec[", partable$parnums[i], "]", sep="")
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
                      ",", partable$group[tmppar], "] <- bpsi[",
                      wishrow, ",", wishcol, ",",
                      partable$group[tmppar], "]", sep="")
      }
    }
  
    ## add priors/constraints after model parameter declarations
    TXT2 <- paste(TXT2, TXT3, sep="")
    
    ## now define inferential covariances and priors for inferential
    ## variances, if needed
    covs <- unique(partable$lhs[grep(".phant", partable$lhs)])
    
    if(length(covs) > 0){
        TXT2 <- paste(TXT2, "\n\n", t1, "# Inferential covariances", sep="")
        for(i in 1:length(covs)){
            for(k in 1:max(partable$group)){
                varlocs <- which(partable$lhs == covs[i] &
                                 partable$op == "=~" &
                                 partable$group == k)
                vartxt <- "star"
                vars <- partable$rhs[varlocs]
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
                phpars <- which(partable$lhs == covs[i] &
                                partable$group == k)
                if(length(phpars) == 1){
                    phpars <- which(partable$rhs == covs[i] &
                                    partable$group == k)
                }

                ## covariances
                TXT2 <- paste(TXT2, "\n", t1, matname, "[", partable$row[var1],
                              ",", partable$row[var2], ",", partable$group[varlocs[1]], "] <- ",
                              partable$mat[phpars[1]], "[", partable$row[phpars[1]], ",",
                              partable$col[phpars[1]], ",", partable$group[phpars[1]], "]*",
                              partable$mat[phpars[2]], "[", partable$row[phpars[2]], ",",
                              partable$col[phpars[2]], ",", partable$group[phpars[2]], "]*",
                              partable$mat[phpars[3]], "[", partable$row[phpars[3]], ",",
                              partable$col[phpars[3]], ",", partable$group[phpars[3]], "]",
                              sep="")
            }

        }
    }
    list(TXT2 = TXT2, partable = partable)
}
