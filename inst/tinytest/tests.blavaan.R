set.seed(8675309)
x1 <- rnorm(100)
x2 <- rnorm(100)
y1 <- 0.5 + 2*x1 + rnorm(100)
Data <- data.frame(y1 = y1, x1 = x1, x2 = x2)

model <- ' y1 ~ x1 '

## auto convergence in stan
expect_error(bsem(model, data=Data, fixed.x=TRUE, target="stan", convergence="auto"))

## seed length != # chains for jags
expect_error(bsem(model, data=Data, fixed.x=TRUE, seed=1, target="jags"))

## unknown cp
expect_error(bsem(model, data=Data, ov.cp="blah", fixed.x=TRUE))

## cp/std.lv clash
expect_error(bsem(model, data=Data, fixed.x=TRUE, std.lv=TRUE, cp="fa"))

model2 <- ' y1 ~ b1*x1 + b2*x2
            b1 + b2 == 0 '

## equality constraint with multiple variables on lhs
expect_error(bsem(model2, data=Data, fixed.x=TRUE))

model2 <- ' y1 ~ b1*x1 + b2*x2
            b1 == -b2/2 '
fit <- bsem(model2, data=Data, target='jags', adapt=1,
            burnin=1, sample=3)
## ensure that == constraints are being respected
expect_true(round(2*fit@Fit@x[1] + fit@Fit@x[2], 5) == 0L)

## do.fit=FALSE
fit <- bsem(model, data=Data, fixed.x=TRUE, adapt=2,
            burnin=2, sample=2, do.fit=FALSE)
expect_equal(class(fit)[1], "blavaan")

## mcmcextra
fit <- bsem(model, data=Data, save.lvs=TRUE, do.fit=FALSE,
            mcmcextra=list(data=list(emiter=101, llnsamp=78)))
expect_equal(class(fit)[1], "blavaan")
expect_equal(fit@external$mcmcdata$emiter, 101L)
expect_equal(fit@Options$llnsamp, 78L)

## vb
fit <- bsem(model, data=Data, target="vb")
expect_equal(class(fit)[1], "blavaan")

## named variable that clashes
names(Data)[1] <- "lambda"
model2 <- ' lambda ~ b1*x1 + b2*x2 '
expect_error(bsem(model2, data=Data))

## one prior on variance, one on sd (problem for target="stan" only)
## and check that defined parameters translate
names(Data)[1] <- "y1"
model3 <- ' y1 ~ a*x1
            x2 ~ b*x1
            y1 ~~ prior("gamma(1,.5)[sd]")*y1
            x2 ~~ prior("gamma(1,.5)[var]")*x2
            pprod := a/b '
expect_error(bsem(model3, data=Data, target="stan"))

## priors are wrong form but will not throw error until estimation
fit <- bsem(model3, data=Data, target="jags", do.fit=FALSE)
expect_equal(class(fit)[1], "blavaan")

fit <- bsem(model3, data=Data, target="stanclassic", do.fit=FALSE)
expect_equal(class(fit)[1], "blavaan")
  
## unknown prior
expect_error(bsem(model, data=Data, dp=dpriors(psi="mydist(1,.5)")))

## wiggle argument
expect_error(bsem(model3, data=Data, wiggle='a', wiggle.sd=0))  ## sd=0 not allowed
expect_error(bsem(model3, data=Data, wiggle='sponge'))          ## sd is string
expect_error(bsem(model3, data=Data, wiggle='b', wiggle.sd=c(1,2))) ## 2 sds, but 1 wiggle
expect_error(bsem(model3, data=Data, wiggle=c('a','b'), wiggle.sd=c(.2,.3), target='jags'))
expect_error(bsem(model3, data=Data, wiggle=c('a','b'), wiggle.sd=c(.2,.3), target='stanclassic')) ## wiggle.sd of length > 1 not allowed for these targets
  
HS.model <- ' visual  =~ x1 + x2 + x3 '

expect_equal(class(bcfa(HS.model, data=HolzingerSwineford1939, target="stan", do.fit=FALSE, group="school", group.equal=c("intercepts","loadings"), wiggle=c("intercepts"), wiggle.sd=.1))[1], "blavaan")
expect_equal(class(bcfa(HS.model, data=HolzingerSwineford1939, target="stanclassic", do.fit=FALSE, group="school", group.equal=c("intercepts","loadings"), wiggle=c("intercepts"), wiggle.sd=.1))[1], "blavaan")
expect_equal(class(bcfa(HS.model, data=HolzingerSwineford1939, target="jags", do.fit=FALSE, group="school", group.equal=c("intercepts","loadings"), wiggle=c("intercepts"), wiggle.sd=.1))[1], "blavaan")

## moment match mcmcextra
set.seed(341)

x1 <- rnorm(100)
y1 <- 0.5 + 2*x1 + rnorm(100)
g <- rep(1:2, each=50)
Data <- data.frame(y1 = y1, x1 = x1, g = g)

model <- ' y1 ~ prior("normal(0,1)")*x1 '
fitstanmomentmatch <- bsem(
  model, 
  data=Data, 
  fixed.x=TRUE, 
  burnin=20,
  sample=20,
  mcmcextra=list(data=list(moment_match_k_threshold=0.5)),
  target="stan",
  seed=1
)
momentmatch_mcobj <- blavInspect(fitstanmomentmatch, "mcobj")
expect_true("Lambda_y_free" %in% names(momentmatch_mcobj@par_dims))
expect_equal(
  fitstanmomentmatch@external$mcmcdata$moment_match_k_threshold,
  0.5
)
