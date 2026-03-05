## tinytest file: Testing blavaan models with Stan backend
##
## Run with: tinytest::run_test_file("tinytest_stan2.R")
## or as part of a package: tinytest::test_package("blavaan")

library("tinytest")
library("lavaan")
library("blavaan")

set.seed(12345)

mytarg <- "stan"

## ── Helpers ──────────────────────────────────────────────────────────────────

## compare logliks in stan vs lavaan
compll <- function(blavobject) {
  lls <- loo::extract_log_lik(blavobject@external$mcmcout)
  lavmcmc <- blavaan:::make_mcmc(blavobject@external$mcmcout)
  blavobject@Options$target <- "jags"
  blavobject@Options$estimator <- "ML"
  lll <- blavaan:::get_ll(lavmcmc[[1]][1,], blavobject)
  out <- FALSE

  if(round(sum(lls[1,]),4) == round(lll[1],4)) out <- TRUE

  out
}


## ensure lavaan function names have not changed
lavmvh1 <- getFromNamespace("lav_mvnorm_missing_h1_estimate_moments", "lavaan")
expect_true(inherits(lavmvh1, "function"))
lavmvll <- getFromNamespace("lav_mvnorm_missing_loglik_data", "lavaan")
expect_true(inherits(lavmvll, "function"))
lav2ll <- getFromNamespace("lav_mvnorm_cluster_loglik_samplestats_2l", "lavaan")
expect_true(inherits(lav2ll, "function"))
lavd <- getFromNamespace("lav_lavdata", "lavaan")
expect_true(inherits(lavd, "function"))
lav_eeta <- getFromNamespace("lav_model_eeta", "lavaan")
expect_true(inherits(lav_eeta, "function"))
lav_implied22l <- getFromNamespace("lav_mvnorm_cluster_implied22l", "lavaan")
expect_true(inherits(lav_implied22l, "function"))
lav_estep <- getFromNamespace("lav_mvnorm_cluster_em_estep_ranef", "lavaan")
expect_true(inherits(lav_estep, "function"))
lavigh <- getFromNamespace("lav_integration_gauss_hermite", "lavaan")
expect_true(inherits(lavigh, "function"))
lavdmvnorm <- getFromNamespace("lav_mvnorm_dmvnorm", "lavaan")
expect_true(inherits(lavdmvnorm, "function"))



## =============================================================================
## 1. Ex 1
## =============================================================================
try({

## Ex 1
model <- ' 
       # latent variable definitions
         ind60 =~ x1 + x2 + x3
         dem60 =~ y1 + a*y2 + b*y3 + c*y4
         dem65 =~ y5 + a*y6 + b*y7 + c*y8
     
       # regressions
         dem60 ~ ind60
         dem65 ~ ind60 + dem60
     
       # residual correlations
         y1 ~~ y5
         y2 ~~ y4 + y6
         y3 ~~ y7
         y4 ~~ y8
         y6 ~~ y8
     '
data(PoliticalDemocracy, package='lavaan')

expect_true(inherits(bsem(model, data=PoliticalDemocracy, target=mytarg, do.fit=FALSE, jagfile=TRUE), "blavaan"))
## mcmcextra arguments
expect_true(inherits(bsem(model, data=PoliticalDemocracy, target=mytarg, do.fit=FALSE, jagfile=TRUE,
                          mcmcextra = list(data=list(emiter=30), monitor="log_lik_rep_sat")), "blavaan"))

fit <- bsem(model, data=PoliticalDemocracy, target=mytarg, burnin=20, sample=20, prisamp=TRUE, wiggle="a", wiggle.sd=.1)
expect_true(inherits(fit, "blavaan"))

## sample.cov
fit <- bsem(model, sample.cov=cov(PoliticalDemocracy), sample.nobs=nrow(PoliticalDemocracy), burnin=200, sample=200, dp=dpriors(lambda="normal(1,.3)"))
expect_true(inherits(fit, "blavaan"))
expect_true(fitMeasures(fit, 'p_dic') > 24 && fitMeasures(fit, 'p_dic') < 28)

## meanstructure/sample.mean arguments without data
expect_error(fit <- bsem(model, sample.cov=cov(PoliticalDemocracy), sample.nobs=nrow(PoliticalDemocracy), sample.mean=colMeans(PoliticalDemocracy), burnin=10, sample=10, dp=dpriors(lambda="normal(1,.3)")))

## sample.mean/sample.cov arguments with data
fit <- bsem(model, data=PoliticalDemocracy, sample.cov=cov(PoliticalDemocracy), sample.nobs=nrow(PoliticalDemocracy), sample.mean=colMeans(PoliticalDemocracy), burnin=200, sample=200, dp=dpriors(lambda="normal(1,.3)"))
fitb <- sem(model, data=PoliticalDemocracy, meanstructure=TRUE)

expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
## for target="stan", ensure starting values of residual correlations are close to 0
expect_true(all(abs(fit@external$inits[[1]]$Theta_r_free) < .1))

## meanstructure
fit <- bsem(model, data=PoliticalDemocracy, target=mytarg, burnin=200, sample=200, bcontrol=list(refresh=-1, cores=3), seed=1163, dp=dpriors(lambda="normal(1,1)"), meanstructure=TRUE)
fit@optim$converged <- TRUE

expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))

expect_silent(summary(fit, neff=TRUE, prior=FALSE))
expect_true(fitMeasures(fit, 'ppp') < .8 & fitMeasures(fit, 'ppp') > .2)
expect_true(compll(fit))

## send in burnin/sample + save lvs
bin <- 200
samp <- 200
pd60 <- PoliticalDemocracy[1:60,]
newd <- PoliticalDemocracy[61:75,]
fit <- bsem(model, data=pd60, burnin=bin, sample=samp, save.lvs=TRUE, target='stan', jags.ic=TRUE)
expect_silent(tmp <- blavInspect(fit, 'lvmeans'))
expect_true(!is.na(tmp[1,1]))
expect_silent(tmp <- blavPredict(fit, type = 'yhat'))
expect_silent(tmp <- blavPredict(fit, type = 'ypred'))
expect_error(blavPredict(fit, type = 'ymis')) # no missing data
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_identical(class(blavPredict(fit))[1], "list")
expect_silent(tmp <- blavInspect(fit, "postmode"))
expect_silent(tmp <- blavInspect(fit, "postmedian"))
expect_silent(tmp <- blavPredict(fit, newdata = newd, type = "ov"))
expect_silent(tmp <- blavPredict(fit, newdata = newd, type = "ypred"))
expect_silent(tmp <- blavPredict(fit, newdata = newd))
lvmn <- Reduce("+", tmp)/length(tmp)
mns <- rowMeans(newd[,9:11])
expect_true(cor(mns, lvmn[,1]) > .95)
mns <- rowMeans(newd[,5:8])
expect_true(cor(mns, lvmn[,3]) > .95)


## same but with meanstructure
fit <- bsem(model, data=PoliticalDemocracy, burnin=bin, sample=samp, save.lvs=TRUE, target='stan', jags.ic=TRUE, meanstructure=TRUE, dp=dpriors(lambda="normal(1,.3)"))
expect_silent(tmp <- blavInspect(fit, 'lvmeans'))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_identical(class(blavPredict(fit))[1], "list")

## std.lv=TRUE now works
fit <- bsem(model, data=PoliticalDemocracy, dp=dpriors(nu="normal(5,10)", target=mytarg), std.lv=TRUE, burnin=200, sample=200, target=mytarg, bcontrol=list(seed=1163))
fitb <- sem(model, data=PoliticalDemocracy, std.lv=TRUE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))

## same but with meanstructure
fit <- bsem(model, data=PoliticalDemocracy, dp=dpriors(nu="normal(5,10)", target=mytarg), std.lv=TRUE, burnin=200, sample=200, target=mytarg, bcontrol=list(seed=1163), meanstructure=TRUE)
fitb <- sem(model, data=PoliticalDemocracy, std.lv=TRUE, meanstructure=TRUE)

expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_true(compll(fit))


## missing data
set.seed(9619)
mis <- matrix(rbinom(prod(dim(PoliticalDemocracy)), 1, .9), nrow(PoliticalDemocracy), ncol(PoliticalDemocracy))
pd <- PoliticalDemocracy*mis
pd[pd==0] <- NA

## lvs lead to chol error, but only if model did not converge
fitm <- bsem(model, data=pd, dp=dpriors(lambda="normal(1,.3)", target=mytarg), burnin=200, sample=200, target=mytarg, save.lvs=TRUE)
fitmb <- sem(model, data=pd, missing='ml', meanstructure=TRUE)

fitm@optim$converged <- TRUE
expect_true(all(abs(coef(fitm) - coef(fitmb)) < 1))
expect_true(all(sqrt(diag(vcov(fitm))) - sqrt(diag(vcov(fitmb))) < 1))
expect_true(fitMeasures(fitm, 'ppp') < .8 & fitMeasures(fitm, 'ppp') > .2)
expect_true(compll(fitm))


## check that 'lvs' and 'lvmeans' match with missingness
tmp <- blavInspect(fitm, 'lvmeans')
tmp2 <- lavPredict(fitmb)
expect_true(cor(as.numeric(tmp), as.numeric(tmp2)) > .99)
tmp3 <- blavInspect(fitm, 'lvs')
tmp3 <- do.call("rbind", tmp3)
expect_true(cor(as.numeric(tmp), as.numeric(colMeans(tmp3))) > .99)


## variation where MLL should be computed even though theta is not diagonal/full
model <- "
  # latent variable definitions
  dem60 =~ y1 + a*y2
  dem65 =~ y5 + a*y6

  # regressions
  dem65 ~ dem60

  # residual correlations
  y1 ~~ y5
"

bfit <- bsem(model, data = PoliticalDemocracy, n.chains = 1, burnin = 50, sample = 100)
expect_true(!is.na(fitMeasures(bfit, "margloglik")[[1]]))

## variation where MLL should *not* be computed when theta is not diagonal/full
model <- "
  # latent variable definitions
  dem60 =~ y1 + a*y2
  dem65 =~ y5 + a*y6

  # regressions
  dem65 ~ dem60

  # residual correlations
  y1 ~~ y5
  y1 ~~ y6
"

bfit <- bsem(model, data = PoliticalDemocracy, n.chains = 1, burnin = 50, sample = 100)
expect_true(is.na(fitMeasures(bfit, "margloglik")[[1]]))





}) ## end try

## =============================================================================
## 2. CFA
## =============================================================================
try({

## CFA
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=20, burnin=10, prisamp=TRUE)
expect_true(inherits(fit2, "blavaan"))
fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=200, burnin=100)
fit2b <- cfa(HS.model, data=HolzingerSwineford1939)

fit2@optim$converged <- TRUE
expect_true(all(abs(coef(fit2) - coef(fit2b)) < 1))
expect_true(all(sqrt(diag(vcov(fit2))) - sqrt(diag(vcov(fit2b))) < 1))
expect_true(fitMeasures(fit2, 'ppp') < .8)


## multi group missing with exo fixed.x
hs39 <- HolzingerSwineford1939
obs <- rbinom(nrow(hs39)*9, 1, .9)
hs39[,7:15] <- hs39[,7:15] * obs
hs39[hs39 == 0] <- NA

exomod <- ' visual  =~ x1 + x2 + x3
            textual =~ x4 + x5 + x6
            speed   =~ x7 + x8 + x9
            visual ~ agemo
            textual ~ agemo
            speed ~ agemo '

fit3 <- bcfa(exomod, data=hs39, target=mytarg, sample=100, burnin=100, group="school", group.equal=c("loadings","intercepts"))
fit3b <- cfa(exomod, data=hs39, group="school", group.equal=c("loadings","intercepts"), missing="ml")

expect_true(inherits(fit3, "blavaan"))
fit3@optim$converged <- TRUE
expect_true(all(abs(coef(fit3) - coef(fit3b)) < 1))
expect_true(all(sqrt(diag(vcov(fit3))) - sqrt(diag(vcov(fit3b))) < 1))
expect_true(fitMeasures(fit3, 'ppp') < .8)

newd <- hs39[c(1:6, 157:162),]
tmp <- blavPredict(fit3, newdata = newd)
lvmn <- Reduce("+", tmp)/length(tmp)
lavfs <- lavPredict(fit3b, newdata = newd)
expect_true(cor(as.numeric(lvmn), as.numeric(do.call("rbind", lavfs))) > .98)



}) ## end try

## =============================================================================
## 3. variation for blavCompare()
## =============================================================================
try({

## variation for blavCompare()
hsm2 <- ' visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6 + x7
          speed   =~ x7 + x8 + x9 '
     
fit2b <- bcfa(hsm2, data=HolzingerSwineford1939, target=mytarg, sample=200, burnin=100)

tmp <- blavCompare(fit2, fit2b)
expect_true(inherits(tmp, "list"))


## same but meanstructure
fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=200, burnin=100, meanstructure=TRUE)
fit2b <- cfa(HS.model, data=HolzingerSwineford1939, meanstructure=TRUE)
expect_true(compll(fit2))

fit2b <- bcfa(hsm2, data=HolzingerSwineford1939, target=mytarg, sample=200, burnin=100, meanstructure=TRUE)

tmp <- blavCompare(fit2, fit2b)
expect_true(inherits(tmp, "list"))



}) ## end try

## =============================================================================
## 4. sample.cov
## =============================================================================
try({

## sample.cov
expect_error(bcfa(HS.model, sample.cov=cov(HolzingerSwineford1939[,paste0('x', 1:9)]), sample.nobs=nrow(HolzingerSwineford1939), meanstructure=TRUE))
fit2 <- bcfa(HS.model, sample.cov=cov(HolzingerSwineford1939[,paste0('x', 1:9)]), sample.nobs=nrow(HolzingerSwineford1939), target=mytarg, sample=200, burnin=100, dp=dpriors(lambda="normal(1,.3)"))
expect_true(inherits(fit2, "blavaan"))
expect_true(fitMeasures(fit2, 'p_dic') > 19 && fitMeasures(fit2, 'p_dic') < 23)



}) ## end try

## =============================================================================
## 5. std.lv=TRUE, where need to change sign of lv covs
## =============================================================================
try({

## std.lv=TRUE, where need to change sign of lv covs
fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=200, burnin=100, std.lv=TRUE, dp=dpriors(psi="gamma(1,1)", target='stan')) # for stanclassic: dp=dpriors(ipsi="gamma(.1,.1)", target='stan'))
fit2b <- cfa(HS.model, data=HolzingerSwineford1939, std.lv=TRUE)

fit2@optim$converged <- TRUE
expect_true(all(abs(coef(fit2) - coef(fit2b)) < 1))
expect_true(all(sqrt(diag(vcov(fit2))) - sqrt(diag(vcov(fit2b))) < 1))
expect_true(fitMeasures(fit2, 'ppp') < .8)


## meanstructure=TRUE
fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=200, burnin=100, std.lv=TRUE, dp=dpriors(psi="gamma(1,1)", target='stan'), meanstructure=TRUE)
fit2b <- cfa(HS.model, data=HolzingerSwineford1939, std.lv=TRUE, meanstructure=TRUE)
expect_true(compll(fit2))



}) ## end try

## =============================================================================
## 6. variation with equality constraints + priors
## =============================================================================
try({

## variation with equality constraints + priors
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9
              x2 ~~ v1*x2 + prior("gamma(2,2)[sd]") * x2
              x3 ~~ v1*x3 + prior("gamma(2,2)[sd]") * x3 '

fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=100, burnin=100)
expect_true(inherits(fit2, "blavaan"))
fit2b <- cfa(HS.model, data=HolzingerSwineford1939, parser='old')

fit2@optim$converged <- TRUE
expect_true(all(abs(coef(fit2) - coef(fit2b)) < 1))
expect_true(all(sqrt(diag(vcov(fit2))) - sqrt(diag(vcov(fit2b))) < 1))
expect_true(fitMeasures(fit2, 'ppp') < .8)




}) ## end try

## =============================================================================
## 7. variation with equality constraints using ==
## =============================================================================
try({

## variation with equality constraints using ==
HS.model <- ' visual  =~ x1 + p1*x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + p2*x8 + x9
              x2 ~~ v1*x2 + prior("gamma(2,2)[sd]") * x2
              x3 ~~ v1*x3 + prior("gamma(2,2)[sd]") * x3
              p1 == p2 '

fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=100, burnin=100)
expect_true(inherits(fit2, "blavaan"))
fit2b <- cfa(HS.model, data=HolzingerSwineford1939, parser='old')

fit2@optim$converged <- TRUE
expect_true(all(abs(coef(fit2) - coef(fit2b)) < 1))
expect_true(all(sqrt(diag(vcov(fit2))) - sqrt(diag(vcov(fit2b))) < 1))
expect_true(fitMeasures(fit2, 'ppp') < .8)




}) ## end try

## =============================================================================
## 8. Multi-group cfa with std.lv and group.equal
## =============================================================================
try({

## Multi-group cfa with std.lv and group.equal
HS.model <- ' visual  =~ x1 + x2 + x3 '

fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=200, burnin=100, group="school", group.equal="loadings", dp=dpriors(lambda="normal(1,.3)"))
fit2b <- cfa(HS.model, data=HolzingerSwineford1939, group="school", group.equal="loadings")

fit2@optim$converged <- TRUE
expect_true(all(abs(coef(fit2) - coef(fit2b)) < 1))
expect_true(all(sqrt(diag(vcov(fit2))) - sqrt(diag(vcov(fit2b))) < 1))
expect_identical(class(fitMeasures(fit2))[1], "lavaan.vector")
expect_true(compll(fit2))

## meanstructure
fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=200, burnin=100, group="school", group.equal="loadings", dp=dpriors(lambda="normal(1,.3)"), meanstructure=TRUE)
fit2b <- cfa(HS.model, data=HolzingerSwineford1939, group="school", group.equal="loadings", meanstructure=TRUE)
fit2@optim$converged <- TRUE
expect_true(all(abs(coef(fit2) - coef(fit2b)) < 1))
expect_true(all(sqrt(diag(vcov(fit2))) - sqrt(diag(vcov(fit2b))) < 1))
expect_identical(class(fitMeasures(fit2))[1], "lavaan.vector")
expect_true(compll(fit2))




}) ## end try

## =============================================================================
## 9. Multi-group sample.cov
## =============================================================================
try({

## Multi-group sample.cov
hs39 <- HolzingerSwineford1939
sc <- list(cov(hs39[hs39$school=='Pasteur',paste0('x',1:9)]), cov(hs39[hs39$school=='Grant-White',paste0('x',1:9)]))
sn <- with(hs39, tapply(x1, school, length))
fit2c <- bcfa(HS.model, sample.cov=sc, sample.nobs=sn, target=mytarg, sample=200, burnin=100, dp=dpriors(lambda="normal(1,.3)"))
expect_true(inherits(fit2c, 'blavaan'))



}) ## end try

## =============================================================================
## 10. Multi-group cfa + missing + save.lvs
## =============================================================================
try({

## Multi-group cfa + missing + save.lvs
set.seed(1050)
obs <- matrix(rbinom(nrow(HolzingerSwineford1939)*3, 1, prob=.9), nrow(HolzingerSwineford1939), 3)
hs39 <- HolzingerSwineford1939[,c(paste0("x",1:3))] * obs
hs39[hs39==0L] <- NA
hs39$school <- HolzingerSwineford1939$school

fit2m <- bcfa(HS.model, data=hs39, target=mytarg, sample=100, burnin=100, group="school", group.equal="loadings", save.lvs=TRUE)
fit2bm <- cfa(HS.model, data=hs39, group="school", group.equal="loadings", missing="ml")

fit2m@optim$converged <- TRUE
expect_true(all(abs(coef(fit2m) - coef(fit2bm)) < 1))
expect_true(all(sqrt(diag(vcov(fit2m))) - sqrt(diag(vcov(fit2bm))) < 1))
expect_identical(class(fitMeasures(fit2m))[1], "lavaan.vector")
expect_true(compll(fit2m))
expect_identical(class(blavPredict(fit2m))[1], "list")
expect_true(inherits(blavCompare(fit2m, fit2m), "list"))

lvout <- blavInspect(fit2m, 'lvs')
expect_true(inherits(lvout, "mcmc.list"))
expect_true(sum(duplicated(lvout[[1]][1,])) == 0)

lvmns <- blavInspect(fit2m, 'lvmeans')
lavlv <- lavPredict(fit2bm)
expect_true(cor(lvmns[1:150,1], lavlv[[1]][1:150]) > .98)



}) ## end try

## =============================================================================
## 11. Similar model, but missing data only in one group
## =============================================================================
try({

## Similar model, but missing data only in one group
hs39 <- HolzingerSwineford1939
hs39$x1[1] <- hs39$x3[1] <- NA
fit2m <- bcfa(HS.model, data=hs39, target=mytarg, sample=100, burnin=100, group="school", group.equal="loadings")

fit2bm <- cfa(HS.model, data=hs39, group="school", group.equal="loadings", missing="ml")

fit2m@optim$converged <- TRUE
expect_true(all(abs(coef(fit2m) - coef(fit2bm)) < 1))
expect_true(all(sqrt(diag(vcov(fit2m))) - sqrt(diag(vcov(fit2bm))) < 1))
expect_identical(class(fitMeasures(fit2m))[1], "lavaan.vector")
expect_true(compll(fit2m))




}) ## end try

## =============================================================================
## 12. multi-group wiggle
## =============================================================================
try({

## multi-group wiggle
fit3 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=200, burnin=100, group="school", group.equal=c("intercepts","loadings"), wiggle=c("intercepts"), wiggle.sd=.1, save.lvs=TRUE)
fit3b <- cfa(HS.model, data=HolzingerSwineford1939, group="school", group.equal=c("intercepts","loadings"), meanstructure=TRUE)
tmp <- blavInspect(fit3, 'lvmeans')
tmp2 <- lavPredict(fit3b)
expect_true(cor(as.numeric(tmp[HolzingerSwineford1939$school=='Pasteur',]), as.numeric(tmp2[[1]])) > .99)
expect_true(fitMeasures(fit3, 'ppp') < .8)
expect_true(compll(fit3))



}) ## end try

## =============================================================================
## 13. Multi-group + missing + save.lvs, with > 1 lv
## =============================================================================
try({

## Multi-group + missing + save.lvs, with > 1 lv
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6 '

set.seed(1051)
obs <- matrix(rbinom(nrow(HolzingerSwineford1939)*6, 1, prob=.9), nrow(HolzingerSwineford1939), 6)
hs39 <- HolzingerSwineford1939[,c(paste0("x",1:6))] * obs
hs39[hs39==0L] <- NA
hs39$school <- HolzingerSwineford1939$school

fit2m <- bcfa(HS.model, data=hs39, target=mytarg, sample=100, burnin=100, group="school", group.equal="loadings", save.lvs=TRUE, dp=dpriors(lambda="normal(1,.3)"))
fit2bm <- cfa(HS.model, data=hs39, group="school", group.equal="loadings", missing="ml")

fit2m@optim$converged <- TRUE
expect_true(all(abs(coef(fit2m) - coef(fit2bm)) < 1))
expect_true(all(sqrt(diag(vcov(fit2m))) - sqrt(diag(vcov(fit2bm))) < 1))
expect_identical(class(fitMeasures(fit2m))[1], "lavaan.vector")
expect_true(compll(fit2m))

lvout <- blavInspect(fit2m, 'lvs')
expect_true(inherits(lvout, "mcmc.list"))
expect_true(sum(duplicated(lvout[[1]][1,])) == 0)

lvmns <- blavInspect(fit2m, 'lvmeans')
lavlv <- lavPredict(fit2bm)
expect_true(cor(lvmns[1:150,1], lavlv[[1]][1:150]) > .98)
expect_true(cor(lvmns[1:150,1], lvout[[1]][1,1:150]) > .75)



}) ## end try

## =============================================================================
## 14. full monty
## =============================================================================
try({

## full monty
data(StereotypeThreat, package="psychotools")
StereotypeThreat <- transform(StereotypeThreat, 
                              group = interaction(ethnicity, condition))
StereotypeThreat <- StereotypeThreat[order(StereotypeThreat$group),]

model <- ' f1 =~ abstract + verbal + c(l1,l1,l1,l4)*numerical
           f1 ~  c(maj,min1,maj,min2)*1 + c(NA,0,NA,0)*1
           abstract ~ c(ar1,ar2,ar3,ar3)*1
           numerical  ~ c(na1,na1,na1,na4)*1
           numerical ~~ c(e1,e1,e1,e4)*numerical
           f1 ~~ c(v1.maj,v1.min,v1.maj,v1.min)*f1
         '

fit <- bcfa(model, data=StereotypeThreat, group="group",
            group.equal=c("loadings", "residuals", "intercepts"),
            burnin=200, sample=200, target=mytarg)

fitb <- cfa(model, data=StereotypeThreat, group="group",
            start="simple",
            group.equal=c("loadings", "residuals", "intercepts"))

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(fitMeasures(fit, 'ppp') < .8)
expect_true(compll(fit))

## make sure wiggle does not lead to strict equality constraints
fit <- bcfa(model, data=StereotypeThreat, group="group",
            group.equal=c("loadings", "residuals", "intercepts"),
            wiggle="intercepts", burnin=200, sample=200, target=mytarg)
expect_true(length(unique(coef(fit)[grepl('ar3',names(coef(fit)))])) == 2) # look at ar3

fit <- bcfa(model, data=StereotypeThreat, group="group",
            group.equal=c("loadings", "residuals", "intercepts"),
            wiggle="loadings", burnin=200, sample=200, target=mytarg)
expect_true(length(unique(coef(fit)[grepl('.p2.',names(coef(fit)))])) == 4)
expect_true(length(unique(coef(fit)[grepl('l1',names(coef(fit)))])) == 3)
expect_true(fitMeasures(fit, 'ppp') < .8)
expect_true(compll(fit))

fit <- bcfa(model, data=StereotypeThreat, group="group",
            group.equal=c("loadings", "residuals", "intercepts"),
            wiggle=c("intercepts","loadings"), burnin=200, sample=200, target=mytarg)
expect_true(inherits(fit, 'blavaan'))

fit <- bcfa(model, data=StereotypeThreat, group="group",
            group.equal=c("loadings", "residuals", "intercepts"),
            wiggle=c("intercepts","loadings"), wiggle.sd=c(.2,.4), burnin=200, sample=200,
            target=mytarg)
expect_true(inherits(fit, 'blavaan'))

## should error:
expect_error(fit <- bcfa(model, data=StereotypeThreat, group="group",
            group.equal=c("loadings", "residuals", "intercepts"),
            wiggle="residuals", burnin=200, sample=200, target=mytarg))

expect_error(fit <- bcfa(model, data=StereotypeThreat, group="group",
            group.equal=c("loadings", "residuals", "intercepts"),
            wiggle=c("intercepts","loadings"), wiggle.sd=c(.2,.4,.6), burnin=200, sample=200,
            target=mytarg))



}) ## end try

## =============================================================================
## 15. Regression model with 0 constraint
## =============================================================================
try({

## Regression model with 0 constraint
data(StereotypeThreat, package="psychotools")
stt <- StereotypeThreat
## use these lines if want to remove missing:
mis <- apply(is.na(stt[,c(5,8,11)]),1,sum)
stt <- stt[mis == 0,]

model <- ' verbal ~ vintelligence + gpa
           vintelligence ~~ 0*gpa '

fit <- bsem(model, data=stt, fixed.x=FALSE, sample=200, target=mytarg, burnin=200)
fitc <- sem(model, data=stt, fixed.x=FALSE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitc)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitc))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))

## meanstructure
fit <- bsem(model, data=stt, fixed.x=FALSE, sample=200, target=mytarg, burnin=200, meanstructure=TRUE)
fitc <- sem(model, data=stt, fixed.x=FALSE, meanstructure=TRUE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitc)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitc))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(compll(fit))


## std.lv=TRUE with dummy lvs only
fit <- bsem(model, data=stt, fixed.x=FALSE, sample=200, target=mytarg, burnin=200, std.lv=TRUE)
fitb <- sem(model, data=stt, fixed.x=FALSE, std.lv=TRUE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))




}) ## end try

## =============================================================================
## 16. Ex 2c
## =============================================================================
try({

## Ex 2c
## lv covariance and residual covariance in same model + multi-group
## equality constraints on residual covariances
## (can constrain across groups or within groups)
HS.modelc <- ' visual  =~ x1 + x2 + x3
               textual =~ x4 + x5 + x6
               speed   =~ x7 + x8 + x9
               x1 ~~ c(a,b)*x4
               x2 ~~ c(a,b)*x5 '

fit2c <- bcfa(HS.modelc, data=HolzingerSwineford1939, group="school", burnin=200, sample=200, target=mytarg, dp=dpriors(lambda='normal(1,.5)'))
fit2d <- cfa(HS.modelc, data=HolzingerSwineford1939, group="school")

fit2c@optim$converged <- TRUE
expect_true(all(abs(coef(fit2c) - coef(fit2d)) < 1))
expect_true(all(sqrt(diag(vcov(fit2c))) - sqrt(diag(vcov(fit2d))) < 1))
expect_identical(class(fitMeasures(fit2c))[1], "lavaan.vector")
expect_true(all(blavInspect(fit2c, 'rhat') < 1.05))
## expect_true(compll(fit2c)) ## FIXME: correlations are constrained equal but not SDs


## equality constraints on loadings and std.lv
HS.modeld <- ' visual  =~ c(a,b)*x1 + x2 + c(e,f)*x3
               textual =~ c(a,b)*x4 + x5 + c(e,f)*x6
               speed   =~ x7 + x8 + x9 '

fit2c <- bcfa(HS.modeld, data=HolzingerSwineford1939, group="school", burnin=200, sample=200, target=mytarg, dp=dpriors(lambda='normal(1,.5)'), std.lv = TRUE)
fit2d <- cfa(HS.modeld, data=HolzingerSwineford1939, group="school", std.lv = TRUE)

fit2c@optim$converged <- TRUE
expect_true(all(abs(coef(fit2c) - coef(fit2d)) < 1))
expect_true(all(sqrt(diag(vcov(fit2c))) - sqrt(diag(vcov(fit2d))) < 1))
expect_identical(class(fitMeasures(fit2c))[1], "lavaan.vector")
expect_true(all(blavInspect(fit2c, 'rhat') < 1.05))
expect_true(compll(fit2c)) ## FIXME: correlations are constrained equal but not SDs




}) ## end try

## =============================================================================
## 17. Ex 5b
## =============================================================================
try({

## Ex 5b, manually fixing covariances to zero in cfa
hsm <- ' visual  =~ x1 + x2 + x3
         textual =~ x4 + x5 + x6
         speed   =~ x7 + x8 + x9
         visual ~~ 0*textual
         textual ~~ 0*speed '

fit5c <- bcfa(hsm, data=HolzingerSwineford1939, sample=200, burnin=200, target=mytarg, dp=dpriors(lambda='normal(1,.5)'))
fit5e <- cfa(hsm, data=HolzingerSwineford1939)

fit5c@optim$converged <- TRUE
expect_true(all(abs(coef(fit5c) - coef(fit5e)) < 1))
expect_true(all(sqrt(diag(vcov(fit5c))) - sqrt(diag(vcov(fit5e))) < 1))
expect_identical(class(fitMeasures(fit5c))[1], "lavaan.vector")
expect_true(all(blavInspect(fit5c, 'rhat') < 1.05))
#expect_true(compll(fit5c)) # comment out because meanstructure=F, so compares wishart with mvn



}) ## end try

## =============================================================================
## 18. std.lv=TRUE with equality constraints
## =============================================================================
try({

## std.lv=TRUE with equality constraints, handling sign indeterminacy
mod2 <- ' visual  =~ c("l1","l1")*x1 + c("l2","l2")*x2 + c("l3","l3")*x3
          textual =~ c("l4","l4")*x4 + c("l5","l5")*x5 + c("l6","l6")*x6
          visual ~~ c(1, NA)*visual
          textual ~~ c(1, NA)*textual '

fit2 <- bcfa(mod2, data = HolzingerSwineford1939, burnin=200, sample=200, group = "school", std.lv = TRUE, meanstructure = TRUE)
fit <- cfa(mod2, data = HolzingerSwineford1939, group = "school", std.lv = TRUE, meanstructure = TRUE)

fit2@optim$converged <- TRUE
expect_true(all(abs(coef(fit2) - coef(fit)) < 1))
expect_true(all(sqrt(diag(vcov(fit2))) - sqrt(diag(vcov(fit))) < 1))
expect_identical(class(fitMeasures(fit2))[1], "lavaan.vector")
expect_true(all(blavInspect(fit2, 'rhat') < 1.05))
expect_true(compll(fit2))



}) ## end try

## =============================================================================
## 19. Ex 6
## =============================================================================
try({

## Ex 6
model <- ' i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
           s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4 '
fit6 <- bgrowth(model, data=Demo.growth, burnin=200, sample=200, target=mytarg)
fit6b <- growth(model, data=Demo.growth)

fit6@optim$converged <- TRUE
expect_true(all(abs(coef(fit6) - coef(fit6b)) < 1))
expect_true(all(sqrt(diag(vcov(fit6))) - sqrt(diag(vcov(fit6b))) < 1))
expect_identical(class(fitMeasures(fit6))[1], "lavaan.vector")
expect_true(all(blavInspect(fit6, 'rhat') < 1.05))
expect_true(compll(fit6))



}) ## end try

## =============================================================================
## 20. finding subblocks of psi matrix for lkj
## =============================================================================
try({

## finding subblocks of psi matrix for lkj
set.seed(1234)
pop.model <- '
  eta =~ y1 + y2
  xi1 =~ y3 + y4 + y5
  xi2 =~ y10 + y11 + y12
  xi3 =~ y16 + y17 + y18
  xi4 =~ y19 + y20 + y21

  eta ~ .8*xi1 + .8*xi2 + .8*xi3 + .8*xi4 '

Data <- simulateData(pop.model, sample.nobs = 150)

model <- '
  eta =~ y1 + y2
  xi1 =~ y3 + y4 + y5
  xi2 =~ y10 + y11 + y12
  xi3 =~ y16 + y17 + y18
  xi4 =~ y19 + y20 + y21

  eta ~ xi1 + xi2 + xi3 + xi4
  xi1 ~~ prior("lkj_corr(4)") * xi2
'

fit <- bsem(model, data = Data, burnin = 150, sample = 150, dp = dpriors(lambda = "normal(1,.5)"), meanstructure = TRUE)
fitb <- sem(model, data = Data, meanstructure = TRUE, parser='old')

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(compll(fit))



}) ## end try

## =============================================================================
## 21. variation with two lv blocks
## =============================================================================
try({

## variation with two lv blocks
model <- '
  eta =~ y1 + y2 + y3
  xi1 =~ y4 + y5
  xi2 =~ y10 + y11 + y12
  xi3 =~ y16 + y17 + y18
  xi4 =~ y19 + y20 + y21

  eta ~ xi2 + xi3 + xi4
  xi1 ~ xi2 + xi3 + xi4
  eta ~~ xi1
  xi2 ~~ prior("lkj_corr(4)") * xi3
'

fit <- bsem(model, data = Data, burnin = 150, sample = 150, dp = dpriors(lambda = "normal(1,.5)"), meanstructure = TRUE)
fitb <- sem(model, data = Data, meanstructure = TRUE, parser='old')

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(compll(fit))




}) ## end try

## =============================================================================
## 22. variation with three lv blocks
## =============================================================================
try({

## variation with three lv blocks of different sizes
## FIXME margloglik() warning tmpmat[lower.tri(tmpmat)] <- lavpartable$est[covpars]
set.seed(1234)
pop.model <- '
  eta =~ y1 + y2
  xi1 =~ y3 + y4 + y5
  xi2 =~ y10 + y11 + y12
  xi3 =~ y16 + y17 + y18
  xi4 =~ y19 + y20 + y21
  xi5 =~ y22 + y23 + y24
  xi6 =~ y25 + y26 + y27
  xi7 =~ y28 + y29 + y30
  xi8 =~ y31 + y32 + y33

  eta ~ .8*xi2 + .8*xi3 + .8*xi4
  xi1 ~ .8*xi2 + .8*xi3 + .8*xi4
  eta ~~ xi1

  xi5 ~ .8*eta + .8*xi1
  xi6 ~ .8*eta + .8*xi1
  xi7 ~ .8*eta + .8*xi1
  xi8 ~ .8*eta + .8*xi1

  xi5 ~~ xi6 + xi7 + xi8
  xi6 ~~ xi7 + xi8
  xi7 ~~ xi8
  '

Data <- simulateData(pop.model, sample.nobs = 150)

model <- '
  eta =~ y1 + y2
  xi1 =~ y3 + y4 + y5
  xi2 =~ y10 + y11 + y12
  xi3 =~ y16 + y17 + y18
  xi4 =~ y19 + y20 + y21
  xi5 =~ y22 + y23 + y24
  xi6 =~ y25 + y26 + y27
  xi7 =~ y28 + y29 + y30
  xi8 =~ y31 + y32 + y33

  eta ~ xi2 + xi3 + xi4
  xi1 ~ xi2 + xi3 + xi4
  eta ~~ xi1

  xi5 ~ eta + xi1
  xi6 ~ eta + xi1
  xi7 ~ eta + xi1
  xi8 ~ eta + xi1

  xi5 ~~ xi6 + xi7 + xi8
  xi6 ~~ xi7 + xi8
  xi7 ~~ xi8
  '

fit <- bsem(model, data = Data, burnin = 150, sample = 150, dp = dpriors(lambda = "normal(1,.3)"), meanstructure = TRUE)
fitb <- sem(model, data = Data, meanstructure = TRUE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(compll(fit))



}) ## end try

## =============================================================================
## 23. write the model so that blocks of psi are not evident from ordering
## =============================================================================
try({

## write the model so that blocks of psi are not evident from ordering
model <- '
  eta =~ y1 + y2
  xi1 =~ y3 + y4 + y5
  xi2 =~ y10 + y11 + y12
  xi5 =~ y22 + y23 + y24
  xi6 =~ y25 + y26 + y27
  xi3 =~ y16 + y17 + y18
  xi7 =~ y28 + y29 + y30
  xi4 =~ y19 + y20 + y21
  xi8 =~ y31 + y32 + y33

  eta ~ xi2 + xi3 + xi4
  xi1 ~ xi2 + xi3 + xi4
  eta ~~ xi1

  xi5 ~ eta + xi1
  xi6 ~ eta + xi1
  xi7 ~ eta + xi1
  xi8 ~ eta + xi1

  xi5 ~~ xi6 + xi7 + xi8
  xi6 ~~ xi7 + xi8
  xi7 ~~ xi8
  '

fit <- bsem(model, data = Data, burnin = 150, sample = 150, dp = dpriors(lambda = "normal(1,.3)"), meanstructure = TRUE)
fitb <- sem(model, data = Data, meanstructure = TRUE)




}) ## end try

## =============================================================================
## 24. Start of Yves's examples / one-factor (no eXo)
## =============================================================================
try({

## Start of Yves's examples
## one-factor (no eXo)
set.seed(1234)
pop.model <- ' f =~ 0.7*y1 + 0.7*y2 + 0.7*y3 + 0.7*y4 + 0.7*y5 '
Data <- simulateData(pop.model, sample.nobs=100)

model <- ' f =~ y1 + y2 + y3 + y4 + y5
           f ~~ prior("gamma(1, 2.5)")*f '
fit <- bsem(model, data=Data, fixed.x=TRUE, burnin=200, sample=200, target=mytarg)
fitb <- sem(model, data=Data, fixed.x=TRUE, parser='old')

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)

## same but meanstructure
fit <- bsem(model, data=Data, fixed.x=TRUE, burnin=200, sample=200, target=mytarg, meanstructure=TRUE)
fitb <- sem(model, data=Data, fixed.x=TRUE, meanstructure=TRUE, parser='old')
expect_true(compll(fit))



}) ## end try

## =============================================================================
## 25. one-factor + eXo
## =============================================================================
try({

## one-factor + eXo
set.seed(1234) ## this seed previously failed on the fixed.x=FALSE model
pop.model <- ' f =~ 0.7*y1 + 0.7*y2 + 0.7*y3 + 0.7*y4 + 0.7*y5
               f ~ (-2.3)*x1 + 0.8*x2
               y1 ~ 0.2*x2 
               y3 ~ 0.7*x1 '
Data <- simulateData(pop.model, sample.nobs=100)

model <- ' f =~ y1 + y2 + y3 + y4 + y5
           f ~ x1 + x2
           y1 ~ x2
           y3 ~ x1 '

## (missing for testing)
datfull <- Data
obs <- rbinom(prod(dim(Data)), 1, .9)
Data <- Data*obs
Data[Data==0] <- NA

fit <- bsem(model, data=Data, fixed.x=TRUE, burnin=200, sample=200, target=mytarg)
fitb <- sem(model, data=Data, fixed.x=TRUE, missing='ml')

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit, postmedian = TRUE, postmode = TRUE), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)
expect_true(compll(fit))

## sample.cov
fit <- bsem(model, sample.cov=cov(datfull), sample.nobs=nrow(datfull), fixed.x=TRUE, burnin=200, sample=200)
expect_true(inherits(fit, 'blavaan'))
expect_true(fitMeasures(fit, 'p_dic') > 12 && fitMeasures(fit, 'p_dic') < 16)


## one-factor + eXo, fixed.x=FALSE (using previous model)
fit <- bsem(model, data=Data, fixed.x=FALSE, sample=200, burnin=200, target=mytarg)
fitb <- sem(model, data=Data, fixed.x=FALSE, missing='ml')

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)
expect_true(compll(fit))


## sample.cov
fit <- bsem(model, sample.cov=cov(datfull), sample.nobs=nrow(datfull), fixed.x=FALSE, burnin=200, sample=200)
expect_true(inherits(fit, 'blavaan'))
expect_true(fitMeasures(fit, 'p_dic') > 14 && fitMeasures(fit, 'p_dic') < 18)



}) ## end try

## =============================================================================
## 26. two-factor (no eXo)
## =============================================================================
try({

## two-factor (no eXo)
## convergence problems for sample.nobs=100; this is 
## due to priors on factor variances.  so custom
## priors are specified.
set.seed(1234)
pop.model <- ' f1 =~ 0.7*y1 + 0.7*y2 + 0.5*y3 
               f2 =~ 0.6*y4 + 0.6*y5 + 0.5*y6 '
Data <- simulateData(pop.model, sample.nobs=500)

model <- ' f1 =~ y1 + y2 + y3 
           f2 =~ y4 + y5 + y6
           f1 ~~ prior("gamma(1,1)[sd]")*f1
           f2 ~~ prior("gamma(1,1)[sd]")*f2 '
fit <- bsem(model, data=Data, fixed.x=TRUE, burnin=200, sample=200, target=mytarg, dp=dpriors(lambda="normal(1,.5)"))
fitb <- sem(model, data=Data, fixed.x=TRUE, parser='old')

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.1))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)



}) ## end try

## =============================================================================
## 27. two-factor + eXo (fixed x)
## =============================================================================
try({

## two-factor + eXo (fixed x)
set.seed(1234)
pop.model <- ' f1 =~ 0.7*y1 + 0.7*y2 + 0.5*y3 
               f2 =~ 0.6*y4 + 0.6*y5 + 0.5*y6
               f1 ~ (-0.3)*x1 + 0.8*x2
               f2 ~ 1*x2 
               y1 ~ 0.2*x2 
               y3 ~ 0.7*x1 '
Data <- simulateData(pop.model, sample.nobs=500)

model <- ' f1 =~ y1 + y2 + y3 
           f2 =~ y4 + y5 + y6
           f1 ~ x1 + x2
           f2 ~ x2
           y1 ~ x2
           y3 ~ x1 '
fit <- bsem(model, data=Data, fixed.x=TRUE, sample=400, target=mytarg, burnin=200)
fitb <- sem(model, data=Data, fixed.x=TRUE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)



}) ## end try

## =============================================================================
## 28. two-factor + eXo (non-fixed x)
## =============================================================================
try({

## two-factor + eXo (non-fixed x)
fit <- bsem(model, data=Data, fixed.x=FALSE, sample=200, target=mytarg, burnin=200)
fitb <- sem(model, data=Data, fixed.x=FALSE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)

## meanstructure
fit <- bsem(model, data=Data, fixed.x=FALSE, sample=200, target=mytarg, burnin=200, meanstructure=TRUE)
fitb <- sem(model, data=Data, fixed.x=FALSE, meanstructure=TRUE)
expect_true(compll(fit))



}) ## end try

## =============================================================================
## 29. 1-variable model
## =============================================================================
try({

## 1-variable model (test 219)
set.seed(1234)
Data <- data.frame(y1 = rnorm(100))

## intercept + variance only
model <- 'y1 ~ 1; y1 ~~ y1' 
fit <- blavaan(model, data=Data, target=mytarg)
fitb <- lavaan(model, data=Data)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)
expect_true(compll(fit))



}) ## end try

## =============================================================================
## 30. 2-variable model: simple regression
## =============================================================================
try({

## 2-variable model: simple regression (fixed x)
set.seed(1234)
x1 <- rnorm(100)
y1 <- 0.5 + 2*x1 + rnorm(100)
Data <- data.frame(y1 = y1, x1 = x1)

model <- ' y1 ~ x1 '
fit <- bsem(model, data=Data, fixed.x=TRUE, sample=200, target=mytarg, burnin=200)
fitb <- sem(model, data=Data, fixed.x=TRUE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)


fit <- bsem(model, data=Data, fixed.x=TRUE, sample=200, target=mytarg, burnin=200, meanstructure=TRUE)
fitb <- sem(model, data=Data, fixed.x=TRUE, meanstructure=TRUE)
expect_true(compll(fit))


## regression with non-fixed x
fit <- bsem(model, data=Data, fixed.x=FALSE, sample=200, target=mytarg, burnin=200)
fitb <- sem(model, data=Data, fixed.x=FALSE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)

fit <- bsem(model, data=Data, fixed.x=FALSE, sample=200, target=mytarg, burnin=200, meanstructure=TRUE)
fitb <- sem(model, data=Data, fixed.x=FALSE, meanstructure=TRUE)
expect_true(compll(fit))


# 3-variable model: simple path analysis
set.seed(1234)
x1 <- rnorm(100)
y1 <- 0.5 + 2*x1 + rnorm(100)
y2 <- 0.8 + 0.4*y1 + rnorm(100)
Data <- data.frame(y1 = y1, y2 = y2, x1 = x1)

model <- ' y2 ~ y1; y1 ~ x1 '
fit <- bsem(model, data=Data, fixed.x=TRUE, sample=600, target=mytarg, burnin=200)
fitb <- sem(model, data=Data, fixed.x=TRUE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)


## path analysis with non-fixed x
fit <- bsem(model, data=Data, fixed.x=FALSE, sample=200, target=mytarg, burnin=200)
fitb <- sem(model, data=Data, fixed.x=FALSE)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)




}) ## end try

## =============================================================================
## 31. LCS model
## =============================================================================
try({

data("Demo.growth")
colnames(Demo.growth)[1:4] = c("X1","X2","X3","X4")

#LCS specification (constant change)
constantX.syntax <-'

#X

#latent variables
lX1 =~ 1*X1; lX2 =~ 1*X2; lX3 =~ 1*X3; lX4 =~ 1*X4;

#autoregressions
lX2 ~ 1*lX1; lX3 ~ 1*lX2; lX4 ~ 1*lX3;

#change - delta; d
dX1 =~ 1*lX2; dX2 =~ 1*lX3; dX3 =~ 1*lX4;

#intercept and slope
intX =~ 1*lX1;
slopeX =~ 1*dX1 + 1*dX2 + 1*dX3;

#residuals equal
X1 ~~ residX*X1; X2 ~~ residX*X2; X3 ~~ residX*X3; X4 ~~ residX*X4;


#manifest means @0
X1 ~ 0*1; X2 ~0*1; X3 ~ 0*1; X4 ~ 0*1;

#auto-proportions
dX1 ~ 0*lX1; dX2 ~ 0*lX2; dX3 ~ 0*lX3;

#slope and intercept means
slopeX ~ 1;
intX ~ 1;

#Latent variances and covariance
slopeX ~~ slopeX;
intX ~~ intX;
slopeX ~~ intX;

#means and vars @0
lX1 ~ 0*1; lX2 ~0*1; lX3 ~ 0*1; lX4 ~ 0*1;
dX1 ~ 0*1; dX2 ~0*1; dX3 ~ 0*1;

lX1 ~~ 0*lX1; lX2 ~~ 0*lX2; lX3 ~~ 0*lX3; lX4 ~~ 0*lX4;
dX1 ~~ 0*dX1; dX2 ~~ 0*dX2; dX3 ~~ 0*dX3;

'
fitb <- lavaan(constantX.syntax, data = Demo.growth[,1:4])

fit <- blavaan(constantX.syntax, data = Demo.growth, sample=200, burnin=200, target=mytarg)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)
expect_true(compll(fit))


## growth model above with missing data
set.seed(560) ## to make sure no fully missing cases
dgmis <- Demo.growth[,1:4] * matrix(rbinom(prod(dim(Demo.growth[,1:4])),
                                           1, .9), nrow(Demo.growth),
                                    4)
dgmis[dgmis==0] <- NA

fit <- blavaan(constantX.syntax, data = dgmis, burnin=200, sample=300, target=mytarg, save.lvs=TRUE)

fitb <- lavaan(constantX.syntax, data = dgmis, missing='ml')

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(fitMeasures(fit,'ppp') > .2 & fitMeasures(fit,'ppp') < .8)
expect_true(compll(fit))

tmp <- blavInspect(fit, 'lvmeans')
tmp2 <- lavPredict(fitb)
expect_true(cor(as.numeric(tmp), as.numeric(tmp2)) > .99)



}) ## end try

## =============================================================================
## 32. defined parameters
## =============================================================================
try({

## defined parameters
set.seed(1234)
X <- rnorm(100)
M <- 0.5*X + rnorm(100)
Y <- 0.7*M + rnorm(100)
Data <- data.frame(X = X, Y = Y, M = M)
Data$g <- rep(1:2, each = 50)
model <- ' # direct effect
             Y ~ c(c1,c2)*X
           # mediator
             M ~ c(a1,a2)*X
             Y ~ c(b1,b2)*M
             M ~~ M
             Y ~~ Y

             Y ~ 1
             M ~ 1
           # indirect effect (a*b)
             ab := a1*b1
           # total effect
             total := c1 + (ab)
           # diff
             d2 := a2 - a1
         '
fit <- bsem(model, data = Data, burnin=200, sample=200, target=mytarg, group = "g")
fitb <- sem(model, data = Data, group = "g")

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit, type='user') - coef(fitb, type='user')) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))
expect_true(inherits(blavInspect(fit, 'mcmc'), 'mcmc.list'))
expect_true(inherits(blavInspect(fit, 'hpd', prob = .9), 'matrix'))
expect_true(nrow(blavInspect(fit, 'hpd')) == 17) ## are we picking up the defined parameters

## now with missing data
set.seed(1002)
obs <- rbinom(prod(dim(Data)), 1, .8)
datmis <- Data*obs
datmis[datmis==0] <- NA

fit <- bsem(model, data = datmis, burnin=200, sample=200, target='stan', group = "g")
fitb <- sem(model, data = datmis, missing='ml', group = "g")

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit, type='user') - coef(fitb, type='user')) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))



}) ## end try

## =============================================================================
## 33. non-recursive (cyclic) model
## =============================================================================
try({

## non-recursive (cyclic) model
set.seed(1234)
pop.model <- ' dv1 ~ 0.3*dv3
               dv2 ~ 0.5*dv1
               dv3 ~ 0.7*dv2 '
Data <- simulateData(pop.model, sample.nobs=500)

model <- ' dv1 ~ dv3
           dv2 ~ dv1
           dv3 ~ dv2 '

fit <- bsem(model, data=Data, burnin=200, sample=200, target='stan', dp=dpriors(beta="normal(.5,.5)"))
fitb <- sem(model, data = Data)

fit@optim$converged <- TRUE
expect_true(all(abs(coef(fit) - coef(fitb)) < 1))
expect_true(all(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(fitb))) < 1))
expect_identical(class(fitMeasures(fit))[1], "lavaan.vector")
expect_true(all(blavInspect(fit, 'rhat') < 1.05))
expect_true(inherits(summary(fit), 'matrix'))




}) ## end try

## =============================================================================
## 34. analytic results from roy levy
## =============================================================================
try({

## analytic results from roy levy
test.data <- as.data.frame(c(70, 80, 96))
colnames(test.data) <- "x"

bcfa.model.expressions.stan <- '

  # Conditional probability of the data
  theta =~ 1*x

  # Prior for the error variance
  # In this case fixed to specified value
  x ~~ 16*x

  # Prior for the true score variance
  # In this case fixed to specified value
  theta ~~ 36*theta

  # Prior for the true score mean
  # In this case fixed to specified value
  theta ~ 80*1

  #theta ~ prior("normal(80,.001)")*1

  # Fix observable intercept to 0
  x ~ prior("normal(0,.001)")*1

'

# For x=70, posterior is N(73.08, 3.33)
# For x=80, posterior is N(80.00, 3.33)
# For x=96, posterior is N(91.08, 3.33)

fit <- bcfa(bcfa.model.expressions.stan, burnin=300, sample=3000, target="stanclassic",  save.lvs=TRUE, test="none", data=test.data)
tmp <- blavInspect(fit, 'lvmeans')
expect_true(max(abs(tmp - c(73.08, 80, 91.08))) < .25)

fit <- bcfa(bcfa.model.expressions.stan, burnin=300, sample=2000, save.lvs=TRUE, test="none", data=test.data)
tmp <- blavInspect(fit, 'lvmeans')
expect_true(max(abs(tmp - c(73.08, 80, 91.08))) < .25)


test.data <- as.data.frame(c(91, 85, 72, 87, 71, 77, 88, 94, 84, 92))
colnames(test.data) <- "x"

blavaan.model.expressions.stan <- '
  # "Equation" for the model
  # This specifies a mean (intercept) for the variable

  x ~ 80*1

  # Variance for the  variable
  x ~~ prior("gamma(5,150)[prec]")*x
'

## analytic posterior of variance is inv-gamma(10,534.5) with post mean 59.39 and sd 21
fit <- blavaan(blavaan.model.expressions.stan, burnin=300, sample=4000, target="stanclassic",  test="none", data=test.data)
expect_true(abs(59.39 - fit@external$stansumm['theta[1,1,1]', 'mean']) < 2*(21/sqrt(10000))) ## posterior sd is 21, so se is about 21/sqrt(nsamps), and we want to be within 2 ses
expect_true(abs(21 - fit@external$stansumm['theta[1,1,1]', 'sd']) < 1)

fit <- blavaan(blavaan.model.expressions.stan, burnin=1000, sample=10000, test="none", data=test.data)
expect_true(abs(59.39 - fit@external$stansumm['Theta_var[1]', 'mean']) < 2*(21/sqrt(30000)))
expect_true(abs(21 - fit@external$stansumm['Theta_var[1]', 'sd']) < 1)



}) ## end try

## =============================================================================
## 35. cross-loadings, sign constraints, std.lv
## =============================================================================
try({

## cross-loadings, sign constraints, std.lv
HS.model <- ' visual  =~ x1 + x2 + x3 + prior("normal(0,.1)")*x4 + prior("normal(0,.1)")*x5 + prior("normal(0,.1)")*x6
              textual =~ x4 + x5 + x6 + prior("normal(0,.1)")*x1 + prior("normal(0,.1)")*x2 + prior("normal(0,.1)")*x3 '

fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, target=mytarg, sample=200, burnin=100, std.lv=TRUE)
expect_true(all(blavInspect(fit2, 'rhat') < 1.05))





}) ## end try

## =============================================================================
## 36. fit indices
## =============================================================================
try({

## fit indices
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
## fit target model
fit1 <- bcfa(HS.model, data = HolzingerSwineford1939, 
             n.chains = 2, burnin = 100, sample = 100)
expect_true(inherits(fit1, "blavaan"))

## fit null model to calculate CFI, TLI, and NFI
null.model <- c(paste0("x", 1:9, " ~~ x", 1:9), paste0("x", 1:9, " ~ 1"))
fit0 <- bcfa(null.model, data = HolzingerSwineford1939, 
             n.chains = 2, burnin = 100, sample = 100)
expect_true(inherits(fit0, "blavaan"))

ML <- blavFitIndices(fit1, baseline.model = fit0)
expect_true(inherits(ML, "blavFitIndices"))

PPMC <- blavFitIndices(fit1, baseline.model = fit0, pD = "waic", rescale = "PPMC")
expect_true(inherits(PPMC, "blavFitIndices"))

expect_true(inherits(summary(PPMC, central.tendency = c("mean","mode"), prob = .95), "data.frame"))


## Access the posterior distributions for further investigation
distML <- data.frame(ML@indices)

nChains <- blavInspect(fit1, "n.chains")
distML$Chain <- rep(1:nChains, each = nrow(distML) / nChains)

p <- do.call(bayesplot::mcmc_pairs, list(x = distML, pars = c("BRMSEA","BMc","BGammaHat","BCFI","BTLI"),
                                         diag_fun = "hist"))
expect_true(inherits(p, "bayesplot_grid"))

## Compare to PPMC method
distPPMC <- data.frame(PPMC@indices)
distPPMC$Chain <- rep(1:nChains, each = nrow(distPPMC) / nChains)
p <- do.call(bayesplot::mcmc_pairs, list(x = distPPMC, pars = c("BRMSEA","BMc","BGammaHat","BCFI","BTLI"), diag_fun = "dens"))
expect_true(inherits(p, "bayesplot_grid"))

}) ## end try
