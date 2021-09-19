dpriors <- function(..., target="stan"){
  userspec <- list(...)

  jagpres <- pkgcheck("runjags")
  stanpres <- pkgcheck("rstan")

  if(jagpres & !stanpres){
    dp <- do.call("jagpriors", userspec)
  } else if(stanpres & !jagpres){
    dp <- do.call("stanpriors", userspec)
  } else if(length(userspec) > 0){
    ## check whether they are supplying jags or stan distributions
    jagdists <- transtables()$disttrans[,'jags']
    ## add other jags dists not in the translation table
    jagdists <- c(jagdists, 'dbetabin', 'ddirch', 'dmnorm',
                  'dwish', 'dmt', 'dmulti',
                  'dbinom', 'dchisq', 'dggamma', # aliases
                  'dnbinom', 'dweibull', 'ddirich')

    userjags <- sapply(jagdists, function(x) grep(x, userspec))

    ## > 1 match can occur for things like ddexp:
    if(length(unlist(userjags)) >= length(userspec)){
      if(target == "jags"){
        dp <- do.call("jagpriors", userspec)
      } else {
        stop("blavaan ERROR: JAGS distributions sent to dpriors(), but target != 'jags'")
      }
    } else if(length(unlist(userjags)) == 0){
      if(target == "jags") stop("blavaan ERROR: target='jags', but no jags distributions were found")
      ## assume they wanted stan
      if(target %in% c("stanclassic", "stancond")){
        dp <- do.call("stanclassicpriors", userspec)
      } else {
        dp <- do.call("stanpriors", userspec)
      }
    } else {
      stop("blavaan ERROR: Distributions sent to dpriors() do not match target.")
    }
  } else {
    ## nothing is user specified, just use target
    if(target == "jags"){
      dp <- do.call("jagpriors", userspec)
    } else if(target %in% c("stanclassic", "stancond")){
      dp <- do.call("stanclassicpriors", userspec)
    } else {
      dp <- do.call("stanpriors", userspec)
    }
  }
  
  dp
}

jagpriors <- function(nu="dnorm(0,1e-3)", alpha="dnorm(0,1e-2)",
                    lambda="dnorm(0,1e-2)", beta="dnorm(0,1e-2)",
                    itheta="dgamma(1,.5)[prec]", ipsi="dgamma(1,.5)[prec]",
                    rho="dbeta(1,1)", ibpsi="dwish(iden,3)",
                    tau="dnorm(0,.44)", delta="dgamma(1,.5)[prec]"){

  dp <- c(nu=nu, alpha=alpha, lambda=lambda, beta=beta,
          itheta=itheta, ipsi=ipsi, rho=rho, ibpsi=ibpsi,
          tau=tau, delta=delta)
  
  dp
}

## see ?stan::expose_stan_functions for obtaining margloglik info
stanpriors <- function(nu="normal(0,32)",
                       alpha="normal(0,10)", lambda="normal(0,10)",
                       beta="normal(0,10)", theta="gamma(1,.5)[sd]",
                       psi="gamma(1,.5)[sd]", rho="beta(1,1)",
                       ibpsi="wishart(3,iden)",
                       tau="normal(0,1.5)",
                       delta="gamma(1,.5)[sd]"){

  dp <- c(nu=nu, alpha=alpha, lambda=lambda, beta=beta,
          theta=theta, psi=psi, rho=rho, ibpsi=ibpsi,
          tau=tau, delta=delta)

  dp
}

stanclassicpriors <- function(nu="normal(0,1000^.5)",
                              alpha="normal(0,10)", lambda="normal(0,10)",
                              beta="normal(0,10)", itheta="gamma(1,.5)[prec]",
                              ipsi="gamma(1,.5)[prec]", rho="beta(1,1)",
                              ibpsi="wishart(3,iden)",
                              tau="normal(0,1.5)",
                              delta="gamma(1,.5)[prec]"){

  dp <- c(nu=nu, alpha=alpha, lambda=lambda, beta=beta,
          itheta=itheta, ipsi=ipsi, rho=rho, ibpsi=ibpsi,
          tau=tau, delta=delta)

  dp
}
