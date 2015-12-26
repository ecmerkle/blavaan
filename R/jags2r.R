jagsdist2r <- function(priors, direction = 'jags2r'){
    ## Convert univariate JAGS distributions to R, relying
    ## on extra packages if needed.  Partially inspired by
    ## LeBauer et al (2013), R Journal, 207-209.

    tabs <- transtables()
    disttrans <- tabs$disttrans
    jag2rfuns <- tabs$jag2rfuns
  
    priargs <- strsplit(priors, "[, ()]+")

    newargs <- lapply(priargs, function(x){
        rnum <- match(x[1], disttrans[,1])
        if(length(x) == 0){
            res <- ""
        } else if(grepl("dwish", x[1])){
            res <- x
        } else {
            trun <- which(x == "T")
            sdvar <- grep("\\[", x)
            if(length(trun) > 0 | length(sdvar) > 0){
                trun <- min(c(trun, sdvar))
                trunargs <- x[trun[1]:length(x)]
                x <- x[1:(trun[1]-1)]
            }
            ## distribution name
            if(is.na(rnum)){
                ## Table 6.4 of jags manual
                aliases <- c("dbin","dbinom","dchisqr","dchisq",
                             "dnegbin","dnbinom","dweib", "dweibull",
                             "ddirch", "ddirich")
                aliases <- t(matrix(aliases, 2, 5))
                rnum <- match(x[1], aliases[,2])
                if(is.na(rnum)){
                    stop("blavaan ERROR: Bad prior specification.")
                } else {
                    rname <- aliases[rnum]
                }
            } else {
                rname <- disttrans[rnum,2]
            }

            ## parameter changes
            rpars <- jag2rfuns[[rnum]](as.numeric(x[2:length(x)]))

            res <- c(rname, rpars)

            ## add truncation info
            if(length(trun) > 0){
                res <- c(res, trunargs)
            }
        }        
        res
    })

    newargs
}

## Might be interesting to also convert R distributions
## to JAGS, but I don't see an immediate application.
##rdist2jags <- function(priors){
##    res <- jagsdist2r(priors, direction = 'r2jags')
##    res
##}

transtables <- function(){
    ## maintain name/parameter translations between r and jags

    ## distribution names for jags and R:
    disttrans <- c("dbeta","dbeta","dchisqr","dchisq","ddexp","ddexp","dexp","dexp","df","df","dgamma","dgamma","dgen.gamma",NA #"rmutil::dggamma",
                   ,"dlogis","dlogis","dlnorm","dlnorm","dnchisqr","dchisq","dnorm","dnorm","dpar",NA,"dt","dt","dunif","dunif","dweib","dweibull",#"dbetabin","VGAM::dbetabinom",
                   "dbern","dbinom","dbin","dbinom","dcat",NA #TODO "dmultinom"
                  ,"dhyper","dhyper","dnegbin","dnbinom","dpois","dpois")#,"ddirch","MCMCpack::ddirichlet","dmnorm","dmnorm","dwish","dwish","dmt",NA,"dmulti","dmultinom")
    ## TODO: How to input non-scalar parameters
    ##       for multivariate distributions? 
    disttrans <- data.frame(t(matrix(disttrans,2,length(disttrans)/2)), stringsAsFactors = FALSE)

    names(disttrans) <- c("jags","r")

    ## functions translating jags parameters to R
    ## element 11 is dnorm, still need others
    jag2rfuns <- vector("list", length=nrow(disttrans))
    for(i in 1:length(jag2rfuns)){
        jag2rfuns[[i]] <- identity
    }

    ## for dnorm + dlnorm:
    dnloc <- which(disttrans$jags == "dnorm")
    jag2rfuns[[dnloc]] <- function(x){
        x[2] <- 1/sqrt(x[2])
        x
    }
    jag2rfuns[[which(disttrans$jags == "dlnorm")]] <- jag2rfuns[[dnloc]]

    ## dlogis + ddexp:
    dlogloc <- which(disttrans$jags == "dlogis")
    jag2rfuns[[dlogloc]] <- function(x){
        x[2] <- 1/x[2]
        x
    }
    jag2rfuns[[which(disttrans$jags == "ddexp")]] <- jag2rfuns[[dlogloc]]

    ## dbin/dbern
    binloc <- which(disttrans$jags == "dbin")
    jag2rfuns[[binloc]] <- function(x) x[2:1]
    bernloc <- which(disttrans$jags == "dbern")
    jag2rfuns[[bernloc]] <- function(x) c(1,x)

    ## others with reversed parameters
    nbloc <- which(disttrans$jags == "dnegbin")
    jag2rfuns[[nbloc]] <- jag2rfuns[[binloc]]
    ## dgamma does not have reversed parameters,
    ## despite LeBauer:
    ##gloc <- which(disttrans$jags == "dgamma")
    ##jag2rfuns[[gloc]] <- jag2rfuns[[binloc]]

    ## weibull
    wbloc <- which(disttrans$jags == "dweib")
    jag2rfuns[[wbloc]] <- function(x){
        ## LeBauer, p. 208
        par2 <- x[2]^(-1/x[1])
        c(x[1], par2)
    }
    
    ## hypergeometric: what to do with noncentral?
    hyloc <- which(disttrans$jags == "dhyper")
    jag2rfuns[[hyloc]] <- function(x){
        if(as.numeric(x[4]) != 1){
            warning("blavaan WARNING: Fit measures with noncentral hypergeometric priors are inaccurate.")
        }
        x[1:3]
    }

    ## t distribution
    tloc <- which(disttrans$jags == "dt")
    jag2rfuns[[tloc]] <- function(x){
        x[c(3,1,2)]
    }
    
    ## beta, poisson, exponential, uniform identical
  
    list(disttrans = disttrans, jag2rfuns = jag2rfuns)
}

## define ddexp(), rdexp() here
ddexp <- function(x, mu = 0, scale = 1, log = FALSE){
    if(scale <= 0) stop("blavaan ERROR: Negative scale parameter to ddexp().")

    dens <- -log(2*scale) - abs(x - mu)/scale
    if(!log) dens <- exp(dens)

    dens
}

rdexp <- function(n, mu = 0, scale = 1){
    if(scale <= 0) stop("blavaan ERROR: Negative scale parameter to rdexp().")

    U <- runif(n, -.5, .5)

    X <- mu - scale * sign(U) * log(1 - 2 * abs(U))

    X
}

## Empirical comparison of univariate distributions in jags vs r
if(FALSE){
    library(runjags)
    source("jags2r.R")
    tt <- transtables()

    compres <- vector("list", 21)
    for(i in 1:21){
        tmppri <- paste(tt$disttrans$jags[i], "(.75", sep="")
        if(tt$disttrans$jags[i] %in% c("dchisqr", "dexp", "dbern", "dpois", "dcat")){
            tmppri <- paste(tmppri, ")", sep="")
        } else if (tt$disttrans$jags[i] %in% c("dgen.gamma", "dbetabin", "dt")){
            tmppri <- paste(tmppri, ",1,5)", sep="")
        } else if (tt$disttrans$jags[i] == "dhyper"){
            tmppri <- paste(tt$disttrans$jags[i], "(3,5,4,1)", sep="")
        } else {
            tmppri <- paste(tmppri, ",5)", sep="")
        }
        tmpmod <- paste("model{\n  y ~ ", tmppri, "\n}", sep="")

        jagres <- run.jags(tmpmod, monitor="y", data=list(x=rep(1,10)), adapt=10,
                           burnin=10, sample=4000, n.chains=1)

        rparam <- jagsdist2r(tmppri)
        if(is.na(rparam[[1]][1])) next
        rfun <- gsub("^d", "r", rparam[[1]][1])
        rargs <- list(10000, as.numeric(rparam[[1]][2]))
        if(length(rparam[[1]]) > 2){
            for(j in 3:length(rparam[[1]])){
                rargs <- c(rargs, list(as.numeric(rparam[[1]][j])))
            }
        }
        ## special handling of t distribution because R
        ## doesn't allow us to set mean/variance
        ## Could also use rt.scaled() from metRology package.
        if(rfun == "rt"){
            mnprec <- c(rargs[[3]][1], rargs[[4]][1])
            rargs[[4]] <- NULL
            rargs[[3]] <- NULL
        }
        rres <- do.call(rfun, rargs)
        if(rfun == "rt"){
            rres <- rres/sqrt(as.numeric(mnprec[2])) + as.numeric(mnprec[1])
        }
        tmpres <- t(matrix(c(summary(as.numeric(jagres$mcmc[[1]])),
                             summary(rres)), 6, 2))
        
        compres[[i]] <- tmpres
    }

    ## check results
    for(i in 1:21){
        if(is.null(compres[[i]])) next
        tmpd <- diff(compres[[i]][,3])
        if(tmpd > .05) print(i)
    }
}
