dpriors <- function(..., target="jags"){
  userspec <- list(...)

  if(target == "jags"){
    dp <- do.call("jagpriors", userspec)
  } else {
    dp <- do.call("stanpriors", userspec)
  }
  
  dp
}

jagpriors <- function(nu="dnorm(0,1e-3)", alpha="dnorm(0,1e-2)",
                    lambda="dnorm(0,1e-2)", beta="dnorm(0,1e-2)",
                    itheta="dgamma(1,.5)", ipsi="dgamma(1,.5)",
                    rho="dbeta(1,1)", ibpsi="dwish(iden,3)",
                    tau="dnorm(0,.1)", delta="dgamma(1,.5)"){

  dp <- c(nu=nu, alpha=alpha, lambda=lambda, beta=beta,
          itheta=itheta, ipsi=ipsi, rho=rho, ibpsi=ibpsi,
          tau=tau, delta=delta)
  
  dp
}

## see ?stan::expose_stan_functions for obtaining margloglik info
stanpriors <- function(nu="normal(0,1000^.5)",
                       alpha="normal(0,10)", lambda="normal(0,10)",
                       beta="normal(0,10)", itheta="gamma(1,.5)",
                       ipsi="gamma(1,.5)", rho="beta(1,1)",
                       ibpsi="wishart(3,iden)",
                       tau="normal(0,10^.5)",
                       delta="gamma(1,.5)"){

  dp <- c(nu=nu, alpha=alpha, lambda=lambda, beta=beta,
          itheta=itheta, ipsi=ipsi, rho=rho, ibpsi=ibpsi,
          tau=tau, delta=delta)

  dp
}
