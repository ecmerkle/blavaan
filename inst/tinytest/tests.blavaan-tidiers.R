## Tests for the broom-style tidy() / glance() methods (R/blav_tidiers.R).
## Models are fit with very short chains: these tests check the *structure*
## of the returned objects, not the numeric estimates.

library("lavaan", quietly = TRUE)
library("blavaan", quietly = TRUE)

set.seed(341)
options(future.globals.maxSize = 1.0 * 1e9)
options(mc.cores = 1L)

## ---------------------------------------------------------------------------
## single-group, single-level model
## ---------------------------------------------------------------------------
data(HolzingerSwineford1939, package = "lavaan")
HS.model <- 'visual =~ x1 + x2 + x3'
fit <- bcfa(HS.model, data = HolzingerSwineford1939, target = "stan",
            seed = 123, n.chains = 2, burnin = 20, sample = 20)

td <- tidy(fit)
expect_true(is.data.frame(td))

## default call exposes all of the optional columns
expect_true(all(c("term", "op", "estimate", "std.error",
                  "conf.low", "conf.high", "std.lv", "std.all",
                  "rhat", "ess", "prior") %in% names(td)))

## one row per summarized parameter
expect_equal(nrow(td), nrow(parameterEstimates(fit)))

expect_true(is.numeric(td$estimate))
expect_true(is.numeric(td$std.error))

## median / mode point estimates must also work. These are regression tests
## for summary(neff = TRUE) / summary(postmedian = TRUE), which read the
## object@external$stansumm matrix.
td_med <- tidy(fit, estimate.method = "median")
expect_equal(nrow(td_med), nrow(td))
expect_true(is.numeric(td_med$estimate))

td_mode <- tidy(fit, estimate.method = "mode")
expect_equal(nrow(td_mode), nrow(td))

## fixed parameters (e.g. the visual=~x1 reference loading, fixed to 1) retain
## their fixed value under every estimate.method rather than becoming NA, and
## no estimate is dropped to NA by the median/mode substitution.
fixed_row <- which(td$term == "visual=~x1")
expect_equal(td$estimate[fixed_row], 1)
expect_equal(td_med$estimate[fixed_row], 1)
expect_equal(td_mode$estimate[fixed_row], 1)
expect_false(anyNA(td_med$estimate))
expect_false(anyNA(td_mode$estimate))

## toggling the optional pieces off drops the corresponding columns
td_min <- tidy(fit, conf.int = FALSE, rhat = FALSE, ess = FALSE, priors = FALSE)
expect_false(any(c("conf.low", "conf.high", "rhat", "ess", "prior") %in% names(td_min)))

## summary(print = FALSE) returns the parameter table invisibly and prints
## nothing (regression test for the restored 'print' argument)
out <- capture.output(s <- summary(fit, header = FALSE, print = FALSE))
expect_equal(length(out), 0L)
expect_true(is.data.frame(s))

## glance() returns the fit measures
g <- glance(fit)
expect_true(inherits(g, "data.frame"))
expect_true(nrow(g) > 1L)

## ---------------------------------------------------------------------------
## multigroup model -> tidy() gains a 'group' column
## ---------------------------------------------------------------------------
fit_g <- bcfa(HS.model, data = HolzingerSwineford1939, group = "school",
              target = "stan", seed = 123, n.chains = 2,
              burnin = 20, sample = 20)
td_g <- tidy(fit_g)
expect_true("group" %in% names(td_g))

## ---------------------------------------------------------------------------
## two-level model -> tidy() gains a 'level' column
## ---------------------------------------------------------------------------
data(Demo.twolevel, package = "lavaan")
model2 <- '
    level: within
        fw =~ y1 + y2 + y3
    level: between
        fb =~ y1 + y2 + y3
'
bfit <- bsem(model2, data = Demo.twolevel, cluster = "cluster",
             target = "stan", seed = 123, n.chains = 2,
             burnin = 20, sample = 20)
td_2l <- tidy(bfit)
expect_true("level" %in% names(td_2l))
