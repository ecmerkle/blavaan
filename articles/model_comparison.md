# Model Comparison

### Introduction

The traditional method for model comparison in frequentist SEM (fSEM) is
the $\chi^{2}$ (Likelihood Ratio Test) and its variations. But for BSEM,
we would take the Bayesian model comparison methods, and apply them to
SEM.

Specifically, we will focus on two information criteria, (1) Widely
Applicable Information Criterion (WAIC), and (2) Leave-One-Out
cross-validation (LOO).

These methods intend to evaluate the out-of-sample predictive accuracy
of the models, and compare that performance. This is the ability to
predict a datapoint that hasn’t been used in the **training** model
(McElreath 2020)

For this example we will use the Industrialization and Political
Democracy example (Bollen 1989).

``` r
model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ a*y1 + b*y2 + c*y3 + d*y4
     dem65 =~ a*y5 + b*y6 + c*y7 + d*y8

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

fit1 <- bsem(model, data=PoliticalDemocracy,
            std.lv=T, meanstructure=T, n.chains=3,
            burnin=500, sample=1000)
```

### Widely Applicable Information Criterion

WAIC (Watanabe 2010) can be seen as a fully Bayesian generalization of
the Akaike Information Criteria (AIC), where we have a measure of
uncertainty/information of the model prediction for each row in the data
across all posterior draws. This is the Log-Pointwise-Predictive-Density
(lppd). The WAIC is defined as

$$\text{WAIC} = - 2\text{lppd} + 2\text{efp}_{\text{WAIC}},$$ The first
term involves the log-likelihoods of observed data (marginal over latent
variables) and the second term is the effective number of parameters.
The first term, $\text{lppd}$, is estimated as:

$$\widehat{\text{lppd}} = \sum\limits_{i = 1}^{n}\log(\frac{1}{S}\sum\limits_{S = 1}^{S}f\left( y_{i}|\theta^{S} \right))$$

where $S$ is the number of posterior draws and
$f\left( y_{i}|\theta^{S} \right)$ is the density of observation $i$
with respect to the parameter sampled at iteration $s$.

The effective number of parameter ($\text{efp}_{\text{WAIC}}$) is
calculated as:

$$\text{efp}_{\text{WAIC}} = \sum\limits_{i = 1}^{n}\text{var}_{s}\left( \log f\left( y_{i}|\theta \right) \right)$$

A separate variance is estimated for each observation $i$ across the $S$
posterior draws.

### Leave-One-Out cross-validation

The $\text{LOO}$ measures the predictive density of each observation
holding out one observation at the time and use the rest of the
observations to update the prior. This estimation is calculated via
(Vehtari, Gelman, and Gabry 2017):

$$\text{LOO} = - 2\sum\limits_{i = 1}^{n}\log(\frac{\sum\limits_{s = 1}^{S}w_{i}^{s}f\left( y_{i}|\theta^{s} \right)}{\sum\limits_{s = 1}^{s}w_{i}^{s}})$$

Where the $w_{i}^{s}$ are Pareto-smoothed sampling weights based on the
relative magnitude of individual $i$ density function across the $S$
posterior samples.

The $\text{LOO}$ effective number of parameters involves the
$\text{lppd}$ term from $\text{WAIC}$:

$$\text{efp}_{\text{LOO}} = \text{lppd} + \text{LOO}/2$$

### Model comparison

As both WAIC and LOO approximate the models’ performance across
posterior draws, we are able to calculate a standard error for them and
for model comparisons involving them.

The model differences estimate the differences across the Expected
Log-Pointwise-Predictive-Density (elpd), and the standard error of the
respective difference.

There are no clear cutoff rules on how to interpret and present these
comparisons, and the researchers need to use their expert knowledge as
part of the decision process. The best recommendation is to present the
differences in elpd $\Delta\text{elpd}$, the standard error, and the
ratio between them. If the ratio is at least $2$ can be consider
evidence of differences between the models, and a ratio of $4$ would be
considered stronger evidence.

For the first example, we will compare the standard political democracy
model, with a model where all factor regressions are fixed to $0$.

``` r
model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ a*y1 + b*y2 + c*y3 + d*y4
     dem65 =~ a*y5 + b*y6 + c*y7 + d*y8

  # regressions
    dem60 ~ 0*ind60
    dem65 ~ 0*ind60 + 0*dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

fit2 <- bsem(model, data=PoliticalDemocracy,
            std.lv=T, meanstructure=T, n.chains=3,
            burnin=500, sample=1000)
```

Once we have the $2$ models, we can compare them with the `blavCompare`

``` r
bc12 <- blavCompare(fit1, fit2)
```

By looking into this comparison object, you can see the WAIC, LOO,
estimates, and the respective differences between them. As these are
information criteria, the **best** model is the one with the lowest
value

``` r
bc12
```

    ## $bf
    ##   bf mll1 mll2 
    ##   NA   NA   NA 
    ## 
    ## $loo
    ## $loo[[1]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##          Estimate   SE
    ## elpd_loo  -1606.0 19.5
    ## p_loo        37.4  2.9
    ## looic      3211.9 39.0
    ## ------
    ## MCSE of elpd_loo is 0.2.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.7, 1.3]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## $loo[[2]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##          Estimate   SE
    ## elpd_loo  -1646.7 18.9
    ## p_loo        34.4  2.7
    ## looic      3293.4 37.8
    ## ------
    ## MCSE of elpd_loo is 0.2.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.6, 1.2]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## 
    ## $diff_loo
    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -40.8       7.9  
    ## 
    ## $waic
    ## $waic[[1]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##           Estimate   SE
    ## elpd_waic  -1605.7 19.5
    ## p_waic        37.1  2.8
    ## waic        3211.3 39.0
    ## 
    ## 38 (50.7%) p_waic estimates greater than 0.4. We recommend trying loo instead. 
    ## 
    ## $waic[[2]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##           Estimate   SE
    ## elpd_waic  -1646.4 18.8
    ## p_waic        34.1  2.7
    ## waic        3292.8 37.7
    ## 
    ## 32 (42.7%) p_waic estimates greater than 0.4. We recommend trying loo instead. 
    ## 
    ## 
    ## $diff_waic
    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -40.7       7.9

In this case we can see that model 1 has lower LOOIC, and the ratio
shows that the LOO differences is $5$ SEs of magnitude. This indicates
that the model with the estimated regressions is better

``` r
abs(bc12$diff_loo[,"elpd_diff"] / bc12$diff_loo[,"se_diff"])
```

    ##   model1   model2 
    ##      NaN 5.156885

Now, lets look at an example with a smaller difference between models,
where only the smallest regression (`dem65~ind60`) is fixed to $0$.

``` r
model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ a*y1 + b*y2 + c*y3 + d*y4
     dem65 =~ a*y5 + b*y6 + c*y7 + d*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ 0*ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

fit3 <- bsem(model, data=PoliticalDemocracy,
            std.lv=T, meanstructure=T, n.chains=3,
            burnin=500, sample=1000)
bc13 <- blavCompare(fit1, fit3)
```

When we see the LOOIC, we see that the difference between the two models
is minimal, and the ratio is $0.21$. This indicates that the models are
functionally equivalent. In a case like this, it is up to the
researchers to decide which model is a **better** representation, and
theoretically stronger.

``` r
bc13
```

    ## $bf
    ##   bf mll1 mll2 
    ##   NA   NA   NA 
    ## 
    ## $loo
    ## $loo[[1]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##          Estimate   SE
    ## elpd_loo  -1606.0 19.5
    ## p_loo        37.4  2.9
    ## looic      3211.9 39.0
    ## ------
    ## MCSE of elpd_loo is 0.2.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.7, 1.3]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## $loo[[2]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##          Estimate   SE
    ## elpd_loo  -1605.9 19.4
    ## p_loo        36.6  2.8
    ## looic      3211.8 38.8
    ## ------
    ## MCSE of elpd_loo is 0.2.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.5, 1.2]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## 
    ## $diff_loo
    ##        elpd_diff se_diff
    ## model2  0.0       0.0   
    ## model1 -0.1       1.0   
    ## 
    ## $waic
    ## $waic[[1]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##           Estimate   SE
    ## elpd_waic  -1605.7 19.5
    ## p_waic        37.1  2.8
    ## waic        3211.3 39.0
    ## 
    ## 38 (50.7%) p_waic estimates greater than 0.4. We recommend trying loo instead. 
    ## 
    ## $waic[[2]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##           Estimate   SE
    ## elpd_waic  -1605.6 19.3
    ## p_waic        36.3  2.8
    ## waic        3211.2 38.7
    ## 
    ## 36 (48.0%) p_waic estimates greater than 0.4. We recommend trying loo instead. 
    ## 
    ## 
    ## $diff_waic
    ##        elpd_diff se_diff
    ## model2  0.0       0.0   
    ## model1 -0.1       1.0

``` r
abs(bc13$diff_loo[,"elpd_diff"] / bc13$diff_loo[,"se_diff"])
```

    ##     model2     model1 
    ##        NaN 0.05527664

Lets do one last model, where only the largest regression
(`dem65~dem60`) is fixed to $0$.

``` r
model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ a*y1 + b*y2 + c*y3 + d*y4
     dem65 =~ a*y5 + b*y6 + c*y7 + d*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + 0*dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

fit4 <- bsem(model, data=PoliticalDemocracy,
            std.lv=T, meanstructure=T, n.chains=3,
            burnin=500, sample=1000)
bc14 <- blavCompare(fit1, fit4)
```

In this case, by looking at the LOOIC, we see that model one is better
(lower value), and the ratio of the difference shows that the model is
$5$ SEs in magnitude. Indicating that there is evidence of model
predictive differences

``` r
bc14
```

    ## $bf
    ##   bf mll1 mll2 
    ##   NA   NA   NA 
    ## 
    ## $loo
    ## $loo[[1]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##          Estimate   SE
    ## elpd_loo  -1606.0 19.5
    ## p_loo        37.4  2.9
    ## looic      3211.9 39.0
    ## ------
    ## MCSE of elpd_loo is 0.2.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.7, 1.3]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## $loo[[2]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##          Estimate   SE
    ## elpd_loo  -1629.5 19.9
    ## p_loo        37.9  3.0
    ## looic      3259.0 39.7
    ## ------
    ## MCSE of elpd_loo is 0.2.
    ## MCSE and ESS estimates assume MCMC draws (r_eff in [0.6, 1.2]).
    ## 
    ## All Pareto k estimates are good (k < 0.7).
    ## See help('pareto-k-diagnostic') for details.
    ## 
    ## 
    ## $diff_loo
    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -23.5       4.0  
    ## 
    ## $waic
    ## $waic[[1]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##           Estimate   SE
    ## elpd_waic  -1605.7 19.5
    ## p_waic        37.1  2.8
    ## waic        3211.3 39.0
    ## 
    ## 38 (50.7%) p_waic estimates greater than 0.4. We recommend trying loo instead. 
    ## 
    ## $waic[[2]]
    ## 
    ## Computed from 3000 by 75 log-likelihood matrix.
    ## 
    ##           Estimate   SE
    ## elpd_waic  -1629.3 19.8
    ## p_waic        37.6  3.0
    ## waic        3258.6 39.7
    ## 
    ## 37 (49.3%) p_waic estimates greater than 0.4. We recommend trying loo instead. 
    ## 
    ## 
    ## $diff_waic
    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -23.6       4.0

``` r
abs(bc14$diff_loo[,"elpd_diff"] / bc14$diff_loo[,"se_diff"])
```

    ##   model1   model2 
    ##      NaN 5.943907

### Bayes factor

In the Bayesian literature you will make use of the Bayes factor (BF) to
compare models. There are a number of criticisms related to the use of
the BF in BSEM, including (1) the BF is unstable for large models (like
most SEMs), (2) it is highly sensitive to model priors, (3) it requires
strong priors to have stable estimation of it, (4) it can require large
number of posterior draws, (5) the estimation using the marginal
likelihood ignores a lot of information from the posterior
distributions. For more details on this discussion please see Tendeiro
and Kiers (2019) and Schad et al. (2022). These criticisms lead us to
recommend against use of the BF in everyday BSEM estimation. For
researchers who commit to their prior distributions and who commit to
exploring the noise in their computations, the BF can used to describe
the relative odds of one model over another, which is more intuitive
than some other model comparison metrics.

### Summary

We recommend the use of LOO or WAIC as general model comparison metrics
for BSEM. They allow us to estimate the models’ out-of-sample predictive
accuracy, and the respective differences across posterior draws. They
also provide us uncertainty estimates in the comparison.

In most cases LOO and WAIC will lead to similar results, and LOO is
recommended as the most stable metric (Vehtari, Gelman, and Gabry 2017).
In general, a $\Delta\text{elpd}$ of at least $2$ standard errors and
preferably $4$ standard errors can be interpreted as evidence of
differential predictive accuracy.

### References

Bollen, Kenneth A. 1989. *Structural Equations with Latent Variables*.
Wiley Series in Probability and Mathematical Statistics. John Wiley &
Sons, Inc.

McElreath, Richard. 2020. *Statistical Rethinking: A Bayesian Course
with Examples in R and Stan*. 2nd ed. CRC Texts in Statistical Science.
Boca Raton: Taylor; Francis, CRC Press.

Schad, Daniel J., Bruno Nicenboim, Paul-Christian Bürkner, Michael
Betancourt, and Shravan Vasishth. 2022. “Workflow Techniques for the
Robust Use of Bayes Factors.” *Psychological Methods*, March.
<https://doi.org/10.1037/met0000472>.

Tendeiro, Jorge N., and Henk A. L. Kiers. 2019. “A Review of Issues
about Null Hypothesis Bayesian Testing.” *Psychological Methods* 24 (6):
774–95. <https://doi.org/10.1037/met0000221>.

Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2017. “Practical Bayesian
Model Evaluation Using Leave-One-Out Cross-Validation and WAIC.”
*Statistics and Computing* 27 (5): 1413–32.
<https://doi.org/10.1007/s11222-016-9696-4>.

Watanabe, Sumio. 2010. “Asymptotic Equivalence of Bayes Cross Validation
and Widely Applicable Information Criterion in Singular Learning
Theory.” *Journal of Machine Learning Research* 11: 3571–94.
