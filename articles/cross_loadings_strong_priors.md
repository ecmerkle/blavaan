# Cross-loadings with strong priors

### Introduction

An advantage of BSEM is that we can use priors to set up **soft**
constraints in the model, by estimating a parameter with a strong prior.
This way the parameter is estimated, but the prior will restrict the
possible values.

This was suggested by Muthén and Asparouhov (2012), as a way to estimate
all possible cross-loadings in a CFA. This way, if the posterior
distribution of the restricted parameters includes values outside of the
strong prior, it can be interpreted as a model modification. This means
that the parameters should be less restricted, or that the prior
distribution should be relaxed.

In this tutorial we present how to estimate a CFA where all possible
cross-loadings are restricted by strong priors.

### Cross-loadings

We will show an example with the Holzinger and Swineford (1939) data.
First we will estimate the regular model with no cross-loadings and
default priors.

``` r

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit_df <- bcfa(HS.model, data=HolzingerSwineford1939, 
            std.lv=TRUE, meanstructure=T)
```

We can see the overall model results with the
[`summary()`](https://rdrr.io/r/base/summary.html) function, looking at
the posterior distribution for the factor loadings, correlations,
intercepts and variances.

``` r

summary(fit_df)
```

    ## blavaan 0.6.1.1487 ended normally after 1000 iterations
    ## 
    ##   Estimator                                      BAYES
    ##   Optimization method                             MCMC
    ##   Number of model parameters                        30
    ## 
    ##   Number of observations                           301
    ## 
    ##   Statistic                                 MargLogLik         PPP
    ##   Value                                      -3870.955       0.000
    ## 
    ## Parameter Estimates:
    ## 
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   visual =~                                                                    
    ##     x1                0.911    0.090    0.738    1.096    1.000    normal(0,10)
    ##     x2                0.500    0.083    0.337    0.667    1.000    normal(0,10)
    ##     x3                0.660    0.079    0.507    0.817    0.999    normal(0,10)
    ##   textual =~                                                                   
    ##     x4                1.003    0.059    0.891    1.120    0.999    normal(0,10)
    ##     x5                1.115    0.063    0.999    1.240    1.000    normal(0,10)
    ##     x6                0.928    0.055    0.819    1.038    0.999    normal(0,10)
    ##   speed =~                                                                     
    ##     x7                0.619    0.078    0.462    0.766    1.000    normal(0,10)
    ##     x8                0.737    0.079    0.581    0.884    1.001    normal(0,10)
    ##     x9                0.680    0.079    0.535    0.843    1.003    normal(0,10)
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   visual ~~                                                                    
    ##     textual           0.449    0.066    0.314    0.574    1.000     lkj_corr(1)
    ##     speed             0.459    0.084    0.291    0.620    1.003     lkj_corr(1)
    ##   textual ~~                                                                   
    ##     speed             0.278    0.070    0.140    0.412    1.001     lkj_corr(1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .x1                4.936    0.069    4.800    5.070    0.999    normal(0,32)
    ##    .x2                6.088    0.068    5.954    6.220    0.999    normal(0,32)
    ##    .x3                2.251    0.065    2.124    2.377    1.000    normal(0,32)
    ##    .x4                3.059    0.069    2.927    3.198    1.000    normal(0,32)
    ##    .x5                4.339    0.076    4.185    4.490    0.999    normal(0,32)
    ##    .x6                2.184    0.065    2.061    2.314    1.001    normal(0,32)
    ##    .x7                4.186    0.064    4.067    4.312    1.000    normal(0,32)
    ##    .x8                5.526    0.060    5.409    5.643    1.000    normal(0,32)
    ##    .x9                5.373    0.060    5.255    5.492    1.001    normal(0,32)
    ##     visual            0.000                                                    
    ##     textual           0.000                                                    
    ##     speed             0.000                                                    
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .x1                0.551    0.132    0.275    0.794    1.000 gamma(1,.5)[sd]
    ##    .x2                1.153    0.108    0.956    1.376    1.000 gamma(1,.5)[sd]
    ##    .x3                0.859    0.098    0.674    1.053    0.999 gamma(1,.5)[sd]
    ##    .x4                0.379    0.050    0.284    0.480    1.000 gamma(1,.5)[sd]
    ##    .x5                0.456    0.059    0.344    0.575    1.000 gamma(1,.5)[sd]
    ##    .x6                0.364    0.046    0.278    0.458    1.000 gamma(1,.5)[sd]
    ##    .x7                0.818    0.091    0.649    1.011    1.001 gamma(1,.5)[sd]
    ##    .x8                0.501    0.095    0.327    0.695    1.003 gamma(1,.5)[sd]
    ##    .x9                0.569    0.096    0.375    0.749    1.005 gamma(1,.5)[sd]
    ##     visual            1.000                                                    
    ##     textual           1.000                                                    
    ##     speed             1.000

Next, we will add all possible cross-loadings with a strong prior of
$`N(0, \sigma = 0.08)`$. The prior centers the loadings around 0 and
allows them little space to move.

``` r

HS.model.cl<-' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 
    
              ## Cross-loadings
              visual =~  prior("normal(0,.08)")*x4 + prior("normal(0,.08)")*x5 + prior("normal(0,.08)")*x6 + prior("normal(0,.08)")*x7 + prior("normal(0,.08)")*x8 + prior("normal(0,.08)")*x9
              textual =~ prior("normal(0,.08)")*x1 + prior("normal(0,.08)")*x2 + prior("normal(0,.08)")*x3 + prior("normal(0,.08)")*x7 + prior("normal(0,.08)")*x8 + prior("normal(0,.08)")*x9 
              speed =~ prior("normal(0,.08)")*x1 + prior("normal(0,.08)")*x2 + prior("normal(0,.08)")*x3 + prior("normal(0,.08)")*x4 + prior("normal(0,.08)")*x5 + prior("normal(0,.08)")*x6'

fit_cl <- bcfa(HS.model.cl, data=HolzingerSwineford1939, 
            std.lv=TRUE, meanstructure=T)
```

It is important that, for each factor, the first variable after `=~` is
one whose loading we expect to be far from 0. So, in the above model, we
specified the regular cfa first (whose loadings we expect to be larger),
then the loadings with small-variance priors on a separate line. This is
important because, in blavaan, the first loading is either constrained
to be positive or fixed to 1 (depending on `std.lv`). If the posterior
distribution of that constrained loading is centered near 0, we may
experience identification problems. Reverse-coded variables can also be
problematic here, because a positive constraint on a reverse-coded
loading can lead other loadings to assume negative values. If you use
informative priors in this situation, then you should verify that the
prior density is on the correct side of 0.

After estimation, you can look at the
[`summary()`](https://rdrr.io/r/base/summary.html) of this model and
evaluate the cross-loadings. You can specifically see whether any of the
cross-loadings seem large enough to suggest that they should be kept in
the model, by looking at the posterior mean (`Estimate`) and credible
interval.

``` r

summary(fit_cl)
```

    ## blavaan 0.6.1.1487 ended normally after 1000 iterations
    ## 
    ##   Estimator                                      BAYES
    ##   Optimization method                             MCMC
    ##   Number of model parameters                        48
    ## 
    ##   Number of observations                           301
    ## 
    ##   Statistic                                 MargLogLik         PPP
    ##   Value                                      -3858.848       0.122
    ## 
    ## Parameter Estimates:
    ## 
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   visual =~                                                                    
    ##     x1                0.763    0.100    0.577    0.969    1.000    normal(0,10)
    ##     x2                0.566    0.089    0.394    0.740    1.000    normal(0,10)
    ##     x3                0.770    0.095    0.589    0.964    1.000    normal(0,10)
    ##   textual =~                                                                   
    ##     x4                0.985    0.062    0.868    1.108    1.000    normal(0,10)
    ##     x5                1.155    0.073    1.017    1.303    1.000    normal(0,10)
    ##     x6                0.894    0.060    0.782    1.015    1.000    normal(0,10)
    ##   speed =~                                                                     
    ##     x7                0.728    0.087    0.563    0.903    1.000    normal(0,10)
    ##     x8                0.793    0.086    0.632    0.967    1.000    normal(0,10)
    ##     x9                0.543    0.073    0.399    0.685    1.001    normal(0,10)
    ##   visual =~                                                                    
    ##     x4                0.031    0.058   -0.082    0.142    0.999   normal(0,.08)
    ##     x5               -0.073    0.064   -0.197    0.052    1.000   normal(0,.08)
    ##     x6                0.062    0.056   -0.050    0.168    1.000   normal(0,.08)
    ##     x7               -0.131    0.067   -0.262    0.000    1.000   normal(0,.08)
    ##     x8               -0.007    0.066   -0.141    0.119    0.999   normal(0,.08)
    ##     x9                0.191    0.060    0.075    0.309    1.000   normal(0,.08)
    ##   textual =~                                                                   
    ##     x1                0.109    0.066   -0.021    0.239    1.000   normal(0,.08)
    ##     x2                0.005    0.060   -0.115    0.122    1.000   normal(0,.08)
    ##     x3               -0.087    0.063   -0.212    0.037    0.999   normal(0,.08)
    ##     x7                0.016    0.062   -0.109    0.137    1.001   normal(0,.08)
    ##     x8               -0.039    0.063   -0.168    0.081    1.000   normal(0,.08)
    ##     x9                0.033    0.054   -0.073    0.141    1.000   normal(0,.08)
    ##   speed =~                                                                     
    ##     x1                0.041    0.065   -0.085    0.167    1.000   normal(0,.08)
    ##     x2               -0.049    0.063   -0.176    0.070    0.999   normal(0,.08)
    ##     x3                0.027    0.064   -0.098    0.151    1.001   normal(0,.08)
    ##     x4               -0.006    0.056   -0.116    0.102    1.000   normal(0,.08)
    ##     x5                0.006    0.060   -0.110    0.126    1.001   normal(0,.08)
    ##     x6               -0.003    0.054   -0.107    0.102    1.000   normal(0,.08)
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   visual ~~                                                                    
    ##     textual           0.381    0.091    0.196    0.553    1.000     lkj_corr(1)
    ##     speed             0.358    0.111    0.132    0.553    1.000     lkj_corr(1)
    ##   textual ~~                                                                   
    ##     speed             0.259    0.105    0.041    0.463    1.001     lkj_corr(1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .x1                4.936    0.068    4.801    5.066    0.999    normal(0,32)
    ##    .x2                6.089    0.068    5.957    6.222    0.999    normal(0,32)
    ##    .x3                2.250    0.068    2.115    2.381    0.999    normal(0,32)
    ##    .x4                3.063    0.066    2.935    3.192    1.000    normal(0,32)
    ##    .x5                4.342    0.074    4.199    4.490    0.999    normal(0,32)
    ##    .x6                2.186    0.062    2.064    2.312    1.000    normal(0,32)
    ##    .x7                4.187    0.065    4.059    4.313    1.000    normal(0,32)
    ##    .x8                5.529    0.060    5.409    5.645    0.999    normal(0,32)
    ##    .x9                5.374    0.059    5.258    5.490    1.000    normal(0,32)
    ##     visual            0.000                                                    
    ##     textual           0.000                                                    
    ##     speed             0.000                                                    
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .x1                0.675    0.107    0.455    0.881    1.000 gamma(1,.5)[sd]
    ##    .x2                1.092    0.107    0.897    1.313    1.000 gamma(1,.5)[sd]
    ##    .x3                0.720    0.112    0.497    0.938    0.999 gamma(1,.5)[sd]
    ##    .x4                0.386    0.051    0.290    0.492    0.999 gamma(1,.5)[sd]
    ##    .x5                0.413    0.063    0.292    0.539    0.999 gamma(1,.5)[sd]
    ##    .x6                0.374    0.044    0.293    0.464    1.000 gamma(1,.5)[sd]
    ##    .x7                0.714    0.097    0.520    0.903    1.000 gamma(1,.5)[sd]
    ##    .x8                0.436    0.097    0.238    0.604    1.001 gamma(1,.5)[sd]
    ##    .x9                0.588    0.066    0.466    0.728    1.000 gamma(1,.5)[sd]
    ##     visual            1.000                                                    
    ##     textual           1.000                                                    
    ##     speed             1.000

We suggest to not simply look at whether the CI excludes 0 (similar to
the null hypothesis), but to evaluate whether the minimum value of the
CI (the value closer to 0) is far enough away from 0 to be relavant
instead of just **different** from 0.

### Caveats

The model with all possible cross-loadings should not be kept as the
final analysis model, but should be used as a step to make decisions
about model changes. This for two main reasons, (1) this model is
overfitted and would present *good* overall fit just due to the
inclusion of a lot of nuisance parameters. In this example the posterior
predictive p-value goes from ppp = 0 to ppp = 0.122, and is not that the
model is better theoretically but that we are inflating the model fit.
And (2), the addition of small-variance priors can prevent detection of
important misspecifications in Bayesian confirmatory factor analysis, as
it can obscure underlying problems in the model by diluting it through a
large number of nuisance parameters (Jorgensen et al. 2019).

### References

Holzinger, K. J., and F. A. Swineford. 1939. *A Study of Factor
Analysis: The Stability of a Bi-Factor Solution*. Supplementary
Educational Monograph 48. Chicago: University of Chicago Press.

Jorgensen, Terrence D, Mauricio Garnier-Villarreal, Sunthud
Pornprasertmanit, and Jaehoon Lee. 2019. “Small-Variance Priors Can
Prevent Detecting Important Misspecifications in Bayesian Confirmatory
Factor Analysis.” In *Quantitative Psychology: The 83rd Annual Meeting
of the Psychometric Society, New York, NY, 2018*, edited by Marie
Wiberg, Steven Culpepper, Rianne Janssen, Jorge González, and Dylan
Molenaar, vol. 265. Springer Proceedings in Mathematics & Statistics.
Springer. <https://doi.org/10.1007/978-3-030-01310-3_23>.

Muthén, Bengt, and Tihomir Asparouhov. 2012. “Bayesian Structural
Equation Modeling: A More Flexible Representation of Substantive
Theory.” *Psychological Methods* 17 (3): 313–35.
<https://doi.org/10.1037/a0026802>.
