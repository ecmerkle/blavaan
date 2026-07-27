## tinytest file: reduce_sum() (use_wcp) correctness, target="stan" only
##
## Run with: tinytest::run_test_file("test_stan_reduce_sum.R")

library("tinytest")
library("lavaan")
library("blavaan")

## only meaningful for target="stan": the compiled rstan stanmodel is reused
## directly with two different data lists (use_wcp on/off) at a fixed
## parameter vector, with no sampling needed
if(Sys.getenv("test_target", unset = "stan") == "stan") {

model <- ' visual  =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed   =~ x7 + x8 + x9 '

fit0 <- bcfa(model, data = HolzingerSwineford1939, target = "stan",
             n.chains = 1, burnin = 5, sample = 5)

sf0  <- fit0@external$mcmcout
dat1 <- fit0@external$mcmcdata
dat1$use_wcp   <- 1L
dat1$grainsize <- 1L

sf1 <- rstan::sampling(getFromNamespace("stanmodels", "blavaan")$stanmarg,
                        data = dat1, chains = 0)

upars <- rstan::unconstrain_pars(sf0, rstan::get_inits(sf0)[[1]])

lp0 <- rstan::log_prob(sf0, upars, adjust_transform = TRUE)
lp1 <- rstan::log_prob(sf1, upars, adjust_transform = TRUE)
expect_equal(lp0, lp1, tolerance = 1e-8)

g0 <- rstan::grad_log_prob(sf0, upars, adjust_transform = TRUE)
g1 <- rstan::grad_log_prob(sf1, upars, adjust_transform = TRUE)
expect_equal(g0, g1, tolerance = 1e-8)

}
