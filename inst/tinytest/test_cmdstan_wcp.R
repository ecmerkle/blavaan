## tinytest file: use_wcp / threads_per_chain defaulting, target="cmdstan" only
##
## Compiling with cpp_options=list(stan_threads=TRUE) is a fresh cmdstan
## recompile, too slow/heavyweight for routine CI; only run with cmdstanr
## installed and Sys.setenv(blavaan_slow_tests = "true")

library("tinytest")
library("lavaan")
library("blavaan")

if (Sys.getenv("blavaan_slow_tests") == "true" &&
    requireNamespace("cmdstanr", quietly = TRUE)) {

model <- ' visual  =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed   =~ x7 + x8 + x9 '

## use_wcp set without bcontrol$threads_per_chain: should warn (reduce_sum()
## would otherwise silently run single-threaded) but still fit successfully
expect_warning(
  fit_wcp <- bcfa(model, data = HolzingerSwineford1939, target = "cmdstan",
                  n.chains = 1, burnin = 20, sample = 20,
                  mcmcextra = list(data = list(use_wcp = 1L, grainsize = 1L))),
  "threads_per_chain"
)
expect_true(inherits(fit_wcp, "blavaan"))

## explicitly supplying threads_per_chain should silence the warning
fit_wcp2 <- bcfa(model, data = HolzingerSwineford1939, target = "cmdstan",
                 n.chains = 1, burnin = 20, sample = 20,
                 mcmcextra = list(data = list(use_wcp = 1L, grainsize = 1L)),
                 bcontrol = list(threads_per_chain = 1L))
expect_true(inherits(fit_wcp2, "blavaan"))

## direct (non-bcontrol) threads_per_chain (issue #57) should also silence
## the warning; cmdstanr's $sample() takes it flat, so no control= nesting
## is expected here
fit_wcp3 <- bcfa(model, data = HolzingerSwineford1939, target = "cmdstan",
                 n.chains = 1, burnin = 20, sample = 20,
                 mcmcextra = list(data = list(use_wcp = 1L, grainsize = 1L)),
                 threads_per_chain = 1L)
expect_true(inherits(fit_wcp3, "blavaan"))
expect_equal(fit_wcp3@optim$control$threads_per_chain, 1L)
expect_true(is.null(fit_wcp3@optim$control$control))

}
