set_priors <- function(priorres, partable, i, varnames, ngroups, type="int", dp, blk, j=NULL,
                       p=NULL, nov=NULL, lv.names=NULL, lv.names.x=NULL, ov.cp=NULL,
                       lv.cp=NULL, lv.x.wish=NULL, mvcovs=NULL) {
  ## Write prior distributions/constraints to TXT2
  ## ... includes j for loadings, regressions,
  ##       lv.nox.int + lv.nox.reg (is lv.nox[j])
  ##     includes p for lv.nox.reg, which is rhs.idx[p]
  ##     includes nov, lv.names, lv.names.x, ov.cp, lv.cp, lv.x.wish, mvcovs for vars
  ##     includes nov, lv.names, lv.names.x, ov.cp, lv.cp, lv.x.wish for covs
  TXT2 <- priorres$TXT2
  coefvec <- priorres$coefvec
  
  ## tabs
  t1 <- paste(rep(" ", 2L), collapse="")
  t2 <- paste(rep(" ", 4L), collapse="")
  t3 <- paste(rep(" ", 6L), collapse="")
  
  if(type=="int"){
    for(k in 1:ngroups){
      int.idx2 <- which(partable$op == "~1" &
                        partable$lhs == varnames[i] &
                        partable$group == k)
      ## fixed?
      if(partable$free[int.idx2]==0){
         TXT2 <- paste(TXT2, t1, "nu[", i, ",", k, "] <- ",
                       partable$ustart[int.idx2], "\n", sep="")
      }
      ## equality constraint?
      else if(any(partable$op == "==" &
                  partable$rhs == partable$plabel[int.idx2])){
        oth.eq <- which(partable$op == "==" &
                        partable$rhs == partable$plabel[int.idx2])
        old.idx <- which(partable$plabel == partable$lhs[oth.eq])
        TXT2 <- paste(TXT2, t1, "nu[", i, ",", k, "] <- nu[",
                      old.idx, ",", partable$group[old.idx], "]\n",
                      sep="")
        coefvec[partable$free[int.idx2],] <- c(paste("nu[", i, ",", k, "]", sep=""),
                                               partable$plabel[int.idx2], "")
      } else {
        if(partable$prior[int.idx2] == ""){
          tmppri <- dp[["nu"]]
        } else {
          tmppri <- partable$prior[int.idx2]
        }
        if(blk){
          TXT2 <- paste(TXT2, t1, "nu[", i, ",", k, "] ~ ", tmppri, "\n", sep="")
        }
        coefvec[partable$free[int.idx2],] <- c(paste("nu[", i, ",", k, "]", sep=""),
                                               partable$plabel[int.idx2],
                                               tmppri)
      }
    } # end k loop
  } # end int

  if(type=="loadings"){
    for(k in 1:ngroups){
      ## For getting lambda's row number:
      lam.idx <- which(partable$op == "=~" &
                       partable$rhs == varnames[i] &
                       partable$group == 1)[j]
      ## For looking up fixed values/constraints:
      lam.idx2 <- which(partable$op == "=~" &
                        partable$rhs == varnames[i] &
                        partable$group == k)[j]
      ## is this associated with a phantom lv?
      phvar <- grep(".phant", partable$lhs[lam.idx2])
      ## fixed?
      if(partable$free[lam.idx2]==0){
        TXT2 <- paste(TXT2, t1, "lambda[", lam.idx, ",", k, "] <- ",
                      partable$ustart[lam.idx2], "\n", sep="")
      }
      ## equality constraint?
      else if(any(partable$op == "==" &
                  partable$rhs == partable$plabel[lam.idx2])){
        oth.eq <- which(partable$op == "==" &
                        partable$rhs == partable$plabel[lam.idx2])
        old.idx <- which(partable$plabel == partable$lhs[oth.eq])
        g1.idx <- which(partable$lhs == partable$lhs[old.idx] &
                        partable$rhs == partable$rhs[old.idx] &
                        partable$op == partable$op[old.idx] &
                        partable$group == 1)
        TXT2 <- paste(TXT2, t1, "lambda[", lam.idx, ",", k,
                      "] <- lambda[", g1.idx, ",",
                      partable$group[old.idx], "]\n",
                      sep="")
        ## avoid phantoms:
        if(length(grep(".phant", partable$lhs[lam.idx2])) == 0){
           coefvec[partable$free[lam.idx2],] <- c(paste("lambda[", lam.idx, ",", k, "]",
                                                        sep=""),
                                                  partable$plabel[lam.idx2], "")
         }
      } else {
        if(partable$prior[lam.idx2] == ""){
          tmppri <- dp[["lambda"]]
        } else {
          tmppri <- partable$prior[lam.idx2]
        }
        if(blk){
          TXT2 <- paste(TXT2, t1, "lambda[", lam.idx, ",", k, "] ~ ", tmppri, "\n",
                        sep="")
        }
        if(length(phvar) == 0){
          coefvec[partable$free[lam.idx2],] <- c(paste("lambda[", lam.idx, ",", k, "]",
                                                       sep=""),
                                                 partable$plabel[lam.idx2],
                                                 tmppri)
        }
      }
    } # end k loop
  } # end loadings

  if(type=="regressions"){
    for(k in 1:ngroups){
      r.idx2 <- which(partable$lhs == varnames[i] &
                      partable$rhs == partable$rhs[j] &
                      partable$group == k)
      ## fixed?
      if(partable$free[r.idx2]==0){
        TXT2 <- paste(TXT2, t1, "beta[", j, ",", k, "] <- ",
                      partable$ustart[r.idx2], "\n", sep="")
      }
      ## equality constraint?
      else if(any(partable$op == "==" &
                  partable$rhs == partable$plabel[r.idx2])){
        oth.eq <- which(partable$op == "==" &
                        partable$rhs == partable$plabel[r.idx2])
        old.idx <- which(partable$plabel == partable$lhs[oth.eq])
        TXT2 <- paste(TXT2, t1, "beta[", j, ",", k,
                      "] <- beta[", old.idx, ",",
                      partable$group[old.idx], "]\n",
                      sep="")
        coefvec[partable$free[r.idx2],] <- c(paste("beta[", j, ",", k, "]", sep=""),
                                             partable$plabel[r.idx2])
      } else {
        if(partable$prior[r.idx2] == ""){
          tmppri <- dp[["beta"]]
        } else {
          tmppri <- partable$prior[r.idx2]
        }
        if(blk){
          TXT2 <- paste(TXT2, t1, "beta[", j, ",", k, "] ~ ", tmppri, "\n", sep="")
        }
        coefvec[partable$free[r.idx2],] <- c(paste("beta[", j, ",", k, "]", sep=""),
                                             partable$plabel[r.idx2],
                                             tmppri)
      }
    } # end k loop
  } # end regressions

  if(type=="lv.nox.int"){
    for(k in 1:ngroups){
      int.idx2 <- which(partable$op == "~1" &
                        partable$lhs == varnames[j] &
                        partable$group == k)
      ## fixed?
      if(partable$free[int.idx2]==0){
        TXT2 <- paste(TXT2, t1, "alpha[", j, ",", k, "] <- ",
                      partable$ustart[int.idx2], "\n", sep="")
      }
      ## equality constraint?
      else if(any(partable$op == "==" &
                  partable$rhs == partable$plabel[int.idx2])){
        ## equal to a previous parameter?  If so, this parameter
        ## gets <-; otherwise it gets a prior
        oth.eq <- which(partable$op == "==" &
                        partable$rhs == partable$plabel[int.idx2])
        old.idx <- which(partable$plabel == partable$lhs[oth.eq])
        TXT2 <- paste(TXT2, t1, "alpha[", j, ",", k, "] <- alpha[",
                      old.idx, ",", partable$group[old.idx], "]\n",
                      sep="")
        coefvec[partable$free[int.idx2],] <- c(paste("alpha[", j, ",", k, "]", sep=""),
                                               partable$plabel[int.idx2], "")
      } else {
        if(partable$prior[int.idx2] == ""){
          tmppri <- dp[["alpha"]]
        } else {
          tmppri <- partable$prior[int.idx2]
        }
        if(blk){
          TXT2 <- paste(TXT2, t1, "alpha[", j, ",", k, "] ~ ", tmppri, "\n", sep="")
        }
        coefvec[partable$free[int.idx2],] <- c(paste("alpha[", j, ",", k, "]", sep=""),
                                               partable$plabel[int.idx2],
                                               tmppri)
      }
    } # end k loop
  } # end lv.nox.int

  if(type=="lv.nox.reg"){
    for(k in 1:ngroups){
      r.idx2 <- which(partable$lhs == j &
                      partable$rhs == partable$rhs[p] &
                      partable$group == k)
      phvar <- grep(".phant", partable$rhs[r.idx2])

      ## fixed?
      if(partable$free[r.idx2]==0){
        TXT2 <- paste(TXT2, t1, "beta[", p, ",", k, "] <- ",
                      partable$ustart[r.idx2], "\n", sep="")
      }
      ## equality constraint?
      else if(any(partable$op == "==" &
                  partable$rhs == partable$plabel[r.idx2])){
        ## equal to a previous parameter?  If so, this parameter
        ## gets <-; otherwise it gets a prior
        oth.eq <- which(partable$op == "==" &
                        partable$rhs == partable$plabel[r.idx2])
        old.idx <- which(partable$plabel == partable$lhs[oth.eq])
        TXT2 <- paste(TXT2, t1, "beta[", p, ",", k, "] <- beta[", old.idx, ",",
                      partable$group[old.idx], "]\n", sep="")
        ## avoid phantoms:
        if(length(phvar) == 0){
          coefvec[partable$free[r.idx2],] <- c(paste("beta[", p, ",", k, "]", sep=""),
                                               partable$plabel[r.idx2], "")
        }
      } else {
        if(partable$prior[r.idx2] == "" & length(phvar) > 0){
          ## get here if cp=="fa"
          TXT2 <- paste(TXT2, t1, "beta[", p, ",", k, "] ~ dnorm(0, 1e-4)\n", sep="")
        } else {
          if(partable$prior[r.idx2] == ""){
            tmppri <- dp[["beta"]]
            if(length(phvar) == 0 & blk){
              TXT2 <- paste(TXT2, t1, "beta[", p, ",", k, "] ~ ", tmppri, "\n", sep="")
            }
          } else {
            tmppri <- partable$prior[r.idx2]
            if(blk){
              TXT2 <- paste(TXT2, t1, "beta[", p, ",", k, "] ~ ", tmppri, "\n", sep="")
            }
          }
          if(length(phvar) == 0){
            coefvec[partable$free[r.idx2],] <- c(paste("beta[", p, ",", k, "]", sep=""),
                                                 partable$plabel[r.idx2],
                                                 tmppri)
          }
        }
      }
    } # end k loop
  } # end lv.nox.reg

  if(type=="vars"){
    ## non-exogenous manifest residuals
    tvname <- "invtheta"
    if(ov.cp == "fa" & mvcovs > 0) tvname <- "invthetstar"
    for(i in 1:nov){
      for(k in 1:ngroups){
        mv.idx <- which(partable$lhs == varnames[i] &
                        partable$lhs == partable$rhs &
                        partable$group == k)
        if(partable$free[mv.idx] == 0){
          TXT2 <- paste(TXT2, t1, tvname, "[", i, ",", k, "] <- ",
                        partable$ustart[mv.idx], "\n", sep="")
        } else if(any(partable$op == "==" &
                      partable$rhs == partable$plabel[mv.idx])){
          oth.eq <- which(partable$op == "==" &
                          partable$rhs == partable$plabel[mv.idx])
          old.idx <- which(partable$plabel == partable$lhs[oth.eq])
          row.idx <- match(partable$lhs[old.idx], varnames)
          TXT2 <- paste(TXT2, t1, tvname, "[", i, ",", k, "] <- ",
                        tvname, "[", row.idx, ",",
                        partable$group[old.idx], "]\n", sep="")
          coefvec[partable$free[mv.idx],] <- c(paste("theta[", i, ",", k, "]", sep=""),
                                               partable$plabel[mv.idx], "")
        } else {
          if(partable$prior[mv.idx] == "" & !grepl("\\[", dp[["itheta"]])){
            TXT2 <- paste(TXT2, t1, tvname, "[", i, ",", k, "] ~ ", dp[["itheta"]], "\n", sep="")
            coefvec[partable$free[mv.idx],] <- c(paste("theta[", i, ",", k, "]", sep=""),
                                                 partable$plabel[mv.idx],
                                                 dp[["itheta"]])
          } else {
            parpri <- tolower(partable$prior[mv.idx])
            vpri <- grepl("\\[var\\]", parpri)
            spri <- grepl("\\[sd\\]", parpri)
            vdpri <- grepl("\\[var\\]", dp[["itheta"]])
            sdpri <- grepl("\\[sd\\]", dp[["itheta"]])
            if(vpri | spri | vdpri | sdpri){
              parname <- ifelse((vpri | vdpri), "thetvar", "thetsd")
              if(!(vpri | spri) & (vdpri | sdpri)) parpri <- dp[["itheta"]]
              sq <- ifelse((vpri | vdpri), "", "^2")
              TXT2 <- paste(TXT2, t1, tvname, "[", i, ",", k,
                            "] <- 1/", parname, i, k, sq, "\n", sep="")
              TXT2 <- paste(TXT2, t1, parname, i, k, " ~ ",
                            strsplit(parpri, "\\[")[[1]][1], "\n", sep="")
            } else {
              TXT2 <- paste(TXT2, t1, tvname, "[", i, ",", k, "] ~ ",
                            parpri, "\n", sep="")
            }
            coefvec[partable$free[mv.idx],] <- c(paste("theta[", i, ",", k, "]", sep=""),
                                                 partable$plabel[mv.idx],
                                                 parpri)
          }
        }
      } # k
    } # i

    TXT2 <- paste(TXT2, "\n", t1, "for(j in 1:", nov, ") {\n", t2,
                  "for(k in 1:", ngroups, ") {\n",
                  t3, "theta[j,k] <- 1/invtheta[j,k]\n",
                  t2, "}\n", t1, "}\n", sep="")

    ## lvs
    nlv <- length(varnames) - nov
    lvname <- "invpsi"
    if(lv.cp == "fa") lvname <- "invpsistar"
    lvstart <- 1
    if(nlv > 0L){
      TXT2 <- paste(TXT2, "\n", sep="")

      nlvx <- length(lv.names.x)
      if(lv.x.wish & nlvx > 1){
        lvstart <- nlvx + 1
        
        TXT2 <- paste(TXT2, t1, "for(k in 1:", ngroups, ") {\n", t2,
                      "ibpsi[1:", nlvx, ",1:", nlvx, ",k] ~ dwish(iden,", nlvx+1, ")\n",
                      sep="")
        TXT2 <- paste(TXT2, t2, "bpsi[1:", nlvx, ",1:", nlvx, ",k] <- inverse(ibpsi[1:",
                      nlvx, ",1:", nlvx, ",k])\n", sep="")
        
        TXT2 <- paste(TXT2, t2, "for(j in 1:", nlvx, ") {\n", t3,
                      lvname, "[j,k] <- 1/bpsi[j,j,k]\n", t2, "}\n",
                      t1, "}\n", sep="")

        ## now add to coefvec
        for(i in 1:nlvx){
          for(k in 1:ngroups){
            lv.idx <- which(partable$lhs == varnames[(nov+i)] &
                            partable$lhs == partable$rhs &
                            partable$group == k)

            coefvec[partable$free[lv.idx],] <- c(paste("psi[", i, ",", k, "]", sep=""),
                                                 partable$plabel[lv.idx],
                                                 dp[["ibpsi"]])
          }
        }
      }

      if(lvstart <= nlv){
        mu.ind.start <- ifelse(lv.x.wish & nlvx > 1, nlvx, 0)
        ptrows <- which(partable$lhs %in% varnames[(nov + lvstart):(nov + nlv)] &
                        !grepl(":", partable$rhs) &
                        partable$lhs == partable$rhs &
                        partable$group == 1)
        lvpt <- partable[ptrows,]
        for(i in lvstart:nlv){
          tmp.eq <- which(partable$op == "==" &
                          partable$lhs == varnames[(nov+i)] &
                          grepl(":", partable$rhs))
          if(length(tmp.eq) > 0) next
          mu.ind <- which(lvpt$lhs == varnames[(nov+i)]) + mu.ind.start
          for(k in 1:ngroups){
            lv.idx <- which(partable$lhs == varnames[(nov+i)] &
                            partable$lhs == partable$rhs &
                            partable$group == k)

            if(partable$free[lv.idx] == 0){
              TXT2 <- paste(TXT2, t1, lvname, "[", mu.ind, ",", k, "] <- ",
                            partable$ustart[lv.idx], "\n", sep="")
            }
            else if(any(partable$op == "==" &
                        partable$rhs == partable$plabel[lv.idx])){
              oth.eq <- which(partable$op == "==" &
                              partable$rhs == partable$plabel[lv.idx])
              pt.idx <- which(partable$plabel == partable$lhs[oth.eq])
              old.idx <- which(lvpt$lhs == partable$lhs[pt.idx]) + mu.ind.start
              TXT2 <- paste(TXT2, t1, lvname, "[", mu.ind, ",", k, "] <- ",
                            lvname, "[", old.idx, ",",
                            partable$group[pt.idx], "]\n", sep="")
              coefvec[partable$free[lv.idx],] <- c(paste("psi[", mu.ind, ",", k, "]", sep=""),
                                                   partable$plabel[lv.idx], "")
            }
            else{
              if(partable$prior[lv.idx] == "" & !grepl("\\[", dp[["ipsi"]])){
                TXT2 <- paste(TXT2, t1, lvname, "[", mu.ind, ",", k, "] ~ ", dp[["ipsi"]], "\n", sep="")
                if(length(grep(".phant", partable$lhs[lv.idx])) == 0){
                  coefvec[partable$free[lv.idx],] <- c(paste("psi[", mu.ind, ",", k, "]", sep=""),
                                                       partable$plabel[lv.idx],
                                                       dp[["ipsi"]])
                }
              } else {
                parpri <- tolower(partable$prior[lv.idx])
                vpri <- grepl("\\[var\\]", parpri)
                spri <- grepl("\\[sd\\]", parpri)
                vdpri <- grepl("\\[var\\]", dp[["ipsi"]])
                sdpri <- grepl("\\[sd\\]", dp[["ipsi"]])
                if(vpri | spri | vdpri | sdpri){
                  parname <- ifelse((vpri | vdpri), "psivar", "psisd")
                  if(!(vpri | spri) & (vdpri | sdpri)) parpri <- dp[["ipsi"]]
                  sq <- ifelse((vpri | vdpri), "", "^2")
                  TXT2 <- paste(TXT2, t1, lvname, "[", mu.ind, ",", k, "] <- 1/",
                                parname, mu.ind, k, sq, "\n", sep="")
                  TXT2 <- paste(TXT2, t1, parname, mu.ind, k, " ~ ",
                                strsplit(parpri, "\\[")[[1]][1], "\n", sep="")
                } else {
                  TXT2 <- paste(TXT2, t1, lvname, "[", mu.ind, ",", k, "] ~ ",
                                partable$prior[lv.idx], "\n", sep="")
                }
                if(length(grep(".phant", partable$lhs[lv.idx])) == 0){
                  coefvec[partable$free[lv.idx],] <- c(paste("psi[", mu.ind, ",", k, "]", sep=""),
                                                       partable$plabel[lv.idx],
                                                       parpri)
                }
              }
            }
          } # end k
        } # end i
      } # end if lvstart >= nlv      

      TXT2 <- paste(TXT2, "\n", t1, "for(j in 1:", nlv, ") {\n", t2,
                    "for(k in 1:", ngroups, ") {\n",
                    t3, "psi[j,k] <- 1/invpsi[j,k]\n",
                    t2, "}\n", t1, "}\n", sep="")
    } # nlv > 0
  } # end vars

  if(type=="covs"){
    phnames <- unique(partable$lhs[grep(".phant", partable$lhs)])
    varnames <- varnames[!(varnames %in% phnames)]
    nlvx <- length(lv.names.x)

    if(ov.cp == "srs"){
      covpar1 <- which(partable$op == "~~" &
                       partable$lhs != partable$rhs &
                       partable$lhs %in% varnames[1:nov] &
                       partable$group == 1)
    } else {
      tmpid <- which(partable$rhs %in% varnames[1:nov] &
                     partable$lhs %in% phnames)
      tmpph <- unique(partable$lhs[tmpid])
      covpar1 <- which(partable$op == "~~" &
                       partable$lhs %in% tmpph &
                       partable$group == 1)
    }

    if(lv.cp == "srs"){
      covpar2 <- which(partable$op == "~~" &
                       partable$lhs != partable$rhs &
                       partable$lhs %in% varnames[(nov+1):length(varnames)] &
                       partable$group == 1)
    } else {
      tmpid <- which(partable$rhs %in% phnames &
                     partable$lhs %in% varnames[(nov+1):length(varnames)])
      tmpph <- unique(partable$rhs[tmpid])
      covpar2 <- which(partable$op == "~~" &
                       partable$lhs %in% tmpph &
                       partable$group == 1)
      ## add lv.x covariances if they come from dmnorm/dwish
      if(lv.x.wish & nlvx > 1){
        covparw <- which(partable$op == "~~" &
                         partable$lhs != partable$rhs &
                         partable$lhs %in% lv.names.x &
                         partable$rhs %in% lv.names.x)
        covpar2 <- c(covpar2, covparw)
      }
    }
    covpars <- c(covpar1, covpar2)

    if(length(covpars) > 0){
      if(!all(partable$lhs[covpars] == partable$rhs[covpars])){
        TXT2 <- paste(TXT2, "\n", t1, "# correlations/covariances \n", sep="")
      }
      for(k in 1:ngroups){
        for(i in 1:length(covpars)){
          c.idx <- which(partable$op == "~~" &
                         partable$lhs == partable$lhs[covpars[i]] &
                         partable$rhs == partable$rhs[covpars[i]] &
                         partable$group == k)
          ## if these are in dmnorm/dwish, take from bpsi
          if(lv.x.wish & nlvx > 1 & partable$lhs[covpars[i]] %in% lv.names.x &
             partable$rhs[covpars[i]] %in% lv.names.x){
            row.idx <- match(partable$lhs[covpars[i]], lv.names.x)
            col.idx <- match(partable$rhs[covpars[i]], lv.names.x)
            TXT2 <- paste(TXT2, t1, "cov[", i, ",", k, "] <- bpsi[",
                          row.idx, ",", col.idx, ",", k, "]\n", sep="")
            coefvec[partable$free[c.idx],] <- c(paste("cov[", i, ",", k, "]", sep=""),
                                                partable$plabel[c.idx], dp[["ibpsi"]])
          }
          else if(partable$free[c.idx] == 0){
            TXT2 <- paste(TXT2, t1, partable$id[c.idx],
                          " <- ", partable$ustart[c.idx], "\n",
                          sep="")
          }
          else if(any(partable$op == "==" &
                      partable$rhs == partable$rhs[c.idx])){
            ## FIXME probably needs "cov[]" instead of plabel below
            oth.eq <- which(partable$op == "==" &
                            partable$rhs == partable$rhs[c.idx])
            old.idx <- which(partable$plabel == partable$lhs[oth.eq])
            TXT2 <- paste(TXT2, t1, "cov[", i, ",", k, "] <- ",
                          partable$plabel[old.idx], "]\n", sep="")
            coefvec[partable$free[c.idx],] <- c(paste("cov[", i, ",", k, "]", sep=""),
                                                partable$plabel[c.idx], "")
          }
          else if(partable$lhs[c.idx] == partable$rhs[c.idx]){ # fa priors
            tmppri <- partable$prior[c.idx]
            coefvec[partable$free[c.idx],] <- c(paste("cov[", i, ",", k, "]", sep=""),
                                                partable$plabel[c.idx], tmppri)
          }
          else{ # srs priors
            coefvec[partable$free[c.idx],1:2] <- c(paste("cov[", i, ",", k, "]", sep=""),
                                                   paste(partable$plabel[c.idx], "@", partable$id[c.idx], sep=""))
            if(grepl("sqrt", partable$id[c.idx]) & grepl("/", partable$id[c.idx])){
              ## equality constraint
              TXT2 <- paste(TXT2, " ", partable$id[c.idx], "\n")
              coefvec[partable$free[c.idx],3] <- ""
            } else {
              tmprc <- strsplit(partable$id[c.idx], "\\[")[[1]][2]
              if(partable$prior[c.idx] == ""){
                TXT2 <- paste(TXT2, t1, partable$id[c.idx],
                              " <- -1 + 2*rstar[", tmprc, "\n", sep="")
                TXT2 <- paste(TXT2, t1, "rstar[", tmprc, " ~ ", dp[["rho"]], "\n", sep="")
                coefvec[partable$free[c.idx],3] <- dp[["rho"]]
              } else {
                TXT2 <- paste(TXT2, t1, partable$id[c.idx],
                              " <- -1 + 2*rstar[", tmprc, "\n", sep="")
                TXT2 <- paste(TXT2, t1, "rstar[", tmprc,
                              " ~ ", partable$prior[c.idx], "\n",
                              sep="")
                coefvec[partable$free[c.idx],3] <- partable$prior[c.idx]
              }
            }
          }
        } # i
      } # k
      TXT2 <- paste(TXT2, "\n", sep="")

      ## add rho parameters to coefvec, for margloglik()
      rps <- grep("rho", coefvec[,2])
      if(length(rps) > 0){
        for(j in rps){
          pnames <- strsplit(coefvec[j,2], "@")
          jname <- strsplit(pnames[[1]][2], " <-")[[1]][1]
          coefvec <- rbind(coefvec, c(jname,
                                      pnames[[1]][1], coefvec[j,3]))
        }
      }
    }
  } # covs
  
  list(TXT2=TXT2, coefvec=coefvec)
}

block_priors <- function(priorres, partable) {
    coefvec <- priorres$coefvec
    TXT2 <- priorres$TXT2
    
    ## tabs
    t1 <- paste(rep(" ", 2L), collapse="")
    t2 <- paste(rep(" ", 4L), collapse="")
    t3 <- paste(rep(" ", 6L), collapse="")
    
    TXT2 <- paste(TXT2, "\n\n", t1, "# Blocked mvn priors", sep="")

    nblk <- max(partable$blk, na.rm=TRUE)

    for(i in 1:nblk){
        blkpars <- partable$plabel[partable$blk == i]

        ## jags names/priors
        jnames <- coefvec[coefvec[,2] %in% blkpars,1]
        jpris <- coefvec[coefvec[,2] %in% blkpars,3]
        jpris <- strsplit(jpris, "[, ()]+")
        jmns <- sapply(jpris, function(x) x[2])
        jprecs <- sapply(jpris, function(x) x[3])

        npars <- length(jnames)

        TXT2 <- paste(TXT2, "\n", t1, "blk", i, " ~ dmnorm(blkmu", i, ",blkprec", i, ")", sep="")
        for(j in 1:npars){
            TXT2 <- paste(TXT2, "\n", t1, jnames[j], " <- blk", i, "[", j, "]", sep="")
            TXT2 <- paste(TXT2, "\n", t1, "blkmu", i, "[", j, "] <- ", jmns[j], sep="")
            TXT2 <- paste(TXT2, "\n", t1, "blkprec", i, "[", j, ",", j, "] <- ", jprecs[j], sep="")
        }
        TXT2 <- paste(TXT2, "\n", t1, "for(i in 1:", npars-1, "){\n", t2,
                      "for(j in (i+1):", npars, "){\n", t3, "blkprec", i, "[i,j] <- 0\n",
                      t3, "blkprec", i, "[j,i] <- 0\n", t2, "}\n", t1, "}\n", sep="")
    }

    list(TXT2=TXT2, coefvec=coefvec)
}
