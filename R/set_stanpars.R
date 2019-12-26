set_stanpars <- function(TXT2, partable, nfree, dp, lv.names.x){
    ## tabs
    t1 <- paste(rep(" ", 2L), collapse="")
    t2 <- paste(rep(" ", 4L), collapse="")
    t3 <- paste(rep(" ", 6L), collapse="")
  
    eqop <- "="
    commop <- "// "
    eolop <- ";"
  
    ## parameter assignments separate from priors
    TXT3 <- paste("\n", t1, commop, "Priors", sep="")

    ## parameter numbers that need priors
    partable$freeparnums <- rep(0, length(partable$id))
    matparnums <- rep(0, length(nfree))
    parvecnum <- 0

    ## get free parameter numbers separately for each parameter type
    for(i in 1:nrow(partable)){
        miscignore <- partable$mat[i] == ""

        eqpar <- which((partable$rhs == partable$plabel[i] &
                       partable$op == "==") |
                       (grepl("rho", partable$mat[i]) &
                        is.na(partable$rhoidx[i])))
        compeq <- which(partable$lhs == partable$plabel[i] &
                        partable$op %in% c("==", ":=") &
                        grepl("\\+|-|/|\\*|\\(|\\)|\\^", partable$rhs))
        fixed <- partable$free[i] == 0 & partable$op[i] != ":="
        if(length(eqpar) > 0 | length(compeq) > 0 | fixed |
           miscignore){
            next
        } else {
            partype <- match(partable$mat[i], names(nfree))
            matparnums[partype] <- matparnums[partype] + 1
            partable$freeparnums[i] <- matparnums[partype]
        }
    }

    for(i in 1:nrow(partable)){
        if(partable$mat[i] != "" | partable$op[i] == ":="){
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

            ## TODO check for inequality constraints here?
          
            ## start parameter assignment
            TXT2 <- paste(TXT2, "\n", t1, partable$mat[i], "[",
                          partable$row[i], ",", partable$col[i],
                          ",", partable$group[i], "] ", eqop,
                          " ", sep="")

            if(grepl("rho", partable$mat[i]) & is.na(partable$ustart[i]) & partable$free[i] > 0){
                TXT2 <- paste(TXT2, "-1 + 2*", sep="")
            }
          
            if((partable$free[i] == 0 & partable$op[i] != ":=") |
               (grepl("rho", partable$mat[i]) & !is.na(partable$ustart[i]))){
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
                if(partable$freeparnums[eqpar] == 0){
                    eqtxt <- paste(partable$mat[eqpar], "[",
                                   partable$row[eqpar], ",",
                                   partable$col[eqpar], ",",
                                   partable$group[eqpar], "]",
                                   sep="")
                } else {
                    eqtxt <- paste(partable$mat[eqpar], "free[",
                                   partable$freeparnums[eqpar],
                                   "]", sep="")
                }

                vpri <- grepl("\\[var\\]", partable$prior[eqpar])
                spri <- grepl("\\[sd\\]", partable$prior[eqpar])
                if(!vpri & (grepl("theta", partable$mat[i]) | grepl("psi", partable$mat[i]))){
                    sq <- ifelse(spri, "2", "-1")
                    TXT2 <- paste(TXT2, "pow(", eqtxt, ",", sq,
                                  ")", eolop, sep="")
                } else {
                    TXT2 <- paste(TXT2, eqtxt, eolop, sep="")
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
                    
                rhstrans <- paste(partable$mat[pvnum], "[", partable$row[pvnum],
                                  ",", partable$col[pvnum], ",", partable$group[pvnum],
                                  "]", sep="")
                ## defined variables involved in another equality
                defvars <- which(partable$mat[pvnum] == "def")
                if(length(defvars) > 0){
                    defpt <- pvnum[defvars]
                    rhstrans[defvars] <- paste0(partable$mat[defpt],
                                                "[", partable$row[defpt],
                                                ",", partable$col[defpt],
                                                ",", partable$group[defpt],
                                                "]")
                }

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
                TXT3 <- paste(TXT3, "\n", t1, "target += ", sep="")
                if(partable$prior[i] == ""){
                    if(partable$mat[i] == "lvrho"){
                        partype <- grep("rho", names(dp))
                    } else if(grepl("star", partable$mat[i])){
                        pname <- paste("i", strsplit(partable$mat[i], "star")[[1]][1], sep="")
                        partype <- grep(pname, names(dp))
                    } else if(grepl("UNC", partable$mat[i])){
                        pname <- strsplit(partable$mat[i], "UNC")[[1]][1]
                        partype <- grep(pname, names(dp))
                    } else {
                        partype <- grep(partable$mat[i], names(dp))
                    }
                    if(length(partype) > 1) partype <- partype[1] # due to psi and ibpsi
                    partable$prior[i] <- dp[partype]

                    ## if rho, re-add prior to cov row
                    if(partable$mat[i] %in% c("rho", "lvrho")){
                        covr <- grep(paste0(partable$mat[i], "[",
                                            partable$row[i], ",",
                                            partable$col[i], ",",
                                            partable$group[i], "]"),
                                     partable$ustart, fixed = TRUE)
                        partable$prior[covr] <- dp[partype]
                    }
                }
                jagpri <- strsplit(partable$prior[i], "\\[")[[1]][1]
                vpri <- grepl("\\[var\\]", partable$prior[i])
                spri <- grepl("\\[sd\\]", partable$prior[i])
                if(vpri){
                    jagpri <- strsplit(partable$prior[i], "\\[var")[[1]][1]
                } else if(spri){
                    jagpri <- strsplit(partable$prior[i], "\\[sd")[[1]][1]
                } else {
                    jagpri <- partable$prior[i]
                }
                splpri <- unlist(strsplit(jagpri, "\\("))
                jagpdist <- paste0(splpri[1], "_lpdf(")
                jagpparm <- paste(splpri[-1], collapse = "(")
                if(!vpri & (grepl("theta", partable$mat[i]) | grepl("psi", partable$mat[i]))){
                    sq <- ifelse(spri, "2", "-1")
                    TXT2 <- paste(TXT2, "pow(", partable$mat[i], "free[",
                                  partable$freeparnums[i], "],", sq,
                                  ")", eolop, sep="")
                } else {
                    TXT2 <- paste(TXT2, partable$mat[i], "free[",
                                  partable$freeparnums[i],
                                  "]", eolop, sep="")
                }
              
              TXT3 <- paste0(TXT3, jagpdist, partable$mat[i], "free[",
                            partable$freeparnums[i], "] | ", jagpparm, eolop)
            }
        }
    }

    list(TXT2 = TXT2, TXT3 = TXT3, partable = partable)
}
