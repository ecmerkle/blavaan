dpriors <- function(nu="dnorm(0,1e-3)", alpha="dnorm(0,1e-2)",
                    lambda="dnorm(0,1e-2)", beta="dnorm(0,1e-2)",
                    itheta="dgamma(1,.5)", ipsi="dgamma(1,.5)",
                    rho="dbeta(1,1)", ibpsi="dwish(iden,3)",
                    tau="dnorm(0,.1)", delta="dgamma(1,.5)"){

  dp <- c(nu=nu, alpha=alpha, lambda=lambda, beta=beta,
          itheta=itheta, ipsi=ipsi, rho=rho, ibpsi=ibpsi,
          tau=tau, delta=delta)
  
  dp
}
