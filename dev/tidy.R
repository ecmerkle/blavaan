# Model ------------------------------------------------------------------
data(HolzingerSwineford1939, package = "lavaan")
HS.model <- 'visual =~ x1 + x2 + x3'
fit <- bcfa(
  HS.model,
  data = HolzingerSwineford1939,
  seed = 123,
  n.chains = 1,
  sample = 300
)

# Standard ---------------------------------------------------------------
debugonce(getMethod("summary", "blavaan"))

summary(fit)


# Desired ----------------------------------------------------------------
# TODO: Compare against how summary() does it
debugonce(tidy.blavaan)
tidy(fit)
tidy(fit, conf.int = TRUE, conf.level = 0.95)
tidy(fit, standardized = TRUE)

# debugonce(compute_blavaan_ci)
# debugonce(get_blavaan_rhat)
# debugonce(get_blavaan_ess)
# debugonce(get_blavaan_priors)


# Two Level --------------------------------------------------------------
data(Demo.twolevel, package = "lavaan")
model <- '
    level: within
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: between
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
'
bfit <- bsem(
  model = model,
  data = Demo.twolevel,
  cluster = "cluster",
  seed = 123,
  n.chains = 1,
  sample = 300,
  target = "stan"
)

debugonce(tidy.blavaan)
tidy(fit)
