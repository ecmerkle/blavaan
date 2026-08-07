## tinytest file: direct sampler-argument passthrough for bcontrol (issue #57)
##
## Fast checks using target="stan" (precompiled stanmarg model, tiny
## burnin/sample) that recognized sampler arguments can be supplied
## directly to bcfa()/bsem()/bgrowth() instead of nested inside
## bcontrol=list(...), and that the old nested style keeps working.

library("tinytest")
library("lavaan")
library("blavaan")

model <- ' visual  =~ x1 + x2 + x3 '

## 1. direct passthrough of a Stan control=list(...) sub-field
fit1 <- bcfa(model, data = HolzingerSwineford1939, target = "stan",
             adapt_delta = 0.8, n.chains = 1, burnin = 20, sample = 20,
             seed = 1)
expect_true(inherits(fit1, "blavaan"))
expect_equal(fit1@optim$control$control$adapt_delta, 0.8)

## 2. old fully-nested bcontrol style still works
fit2 <- bcfa(model, data = HolzingerSwineford1939, target = "stan",
             bcontrol = list(control = list(adapt_delta = 0.85)),
             n.chains = 1, burnin = 20, sample = 20, seed = 1)
expect_true(inherits(fit2, "blavaan"))
expect_equal(fit2@optim$control$control$adapt_delta, 0.85)

## 3. mixed: direct arg + bcontrol-nested arg on different keys merge cleanly
fit3 <- bcfa(model, data = HolzingerSwineford1939, target = "stan",
             adapt_delta = 0.7,
             bcontrol = list(control = list(max_treedepth = 12)),
             n.chains = 1, burnin = 20, sample = 20, seed = 1)
expect_true(inherits(fit3, "blavaan"))
expect_equal(fit3@optim$control$control$adapt_delta, 0.7)
expect_equal(fit3@optim$control$control$max_treedepth, 12)

## 4. error: same argument supplied both directly and inside bcontrol
expect_error(
  bcfa(model, data = HolzingerSwineford1939, target = "stan",
       adapt_delta = 0.7, bcontrol = list(adapt_delta = 0.9),
       n.chains = 1, burnin = 20, sample = 20),
  "supplied both directly"
)

## 5. error: same argument supplied flat and nested within one bcontrol= call
expect_error(
  bcfa(model, data = HolzingerSwineford1939, target = "stan",
       bcontrol = list(adapt_delta = 0.7, control = list(adapt_delta = 0.9)),
       n.chains = 1, burnin = 20, sample = 20),
  "supplied both flat and inside"
)

## 6. reserved names (verbose, cl) are not intercepted -- still governed by
## lavaan's own option / still require bcontrol=list(...)
fit6 <- bcfa(model, data = HolzingerSwineford1939, target = "stan",
             verbose = FALSE, n.chains = 1, burnin = 20, sample = 20,
             seed = 1)
expect_true(inherits(fit6, "blavaan"))

## 7. target="vb": recognized vb argument stays flat (no control= nesting)
fitvb <- bcfa(model, data = HolzingerSwineford1939, target = "vb",
              adapt_engaged = TRUE, seed = 1)
expect_true(inherits(fitvb, "blavaan"))
expect_true(is.null(fitvb@optim$control$control))
expect_equal(fitvb@optim$control$adapt_engaged, TRUE)
