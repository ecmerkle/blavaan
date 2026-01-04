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

    ## blavaan 0.5.9.1382 ended normally after 1000 iterations
    ## 
    ##   Estimator                                      BAYES
    ##   Optimization method                             MCMC
    ##   Number of model parameters                        30
    ## 
    ##   Number of observations                           301
    ## 
    ##   Statistic                                 MargLogLik         PPP
    ##   Value                                      -3871.000       0.000
    ## 
    ## Parameter Estimates:
    ## 
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   visual =~                                                                    
    ##     x1                0.910    0.087    0.746    1.083    1.000    normal(0,10)
    ##     x2                0.500    0.082    0.343    0.661    0.999    normal(0,10)
    ##     x3                0.662    0.079    0.508    0.819    1.000    normal(0,10)
    ##   textual =~                                                                   
    ##     x4                1.001    0.058    0.891    1.118    1.000    normal(0,10)
    ##     x5                1.114    0.063    0.993    1.241    1.000    normal(0,10)
    ##     x6                0.927    0.054    0.822    1.035    1.001    normal(0,10)
    ##   speed =~                                                                     
    ##     x7                0.616    0.077    0.460    0.761    1.001    normal(0,10)
    ##     x8                0.733    0.079    0.579    0.889    1.001    normal(0,10)
    ##     x9                0.680    0.081    0.524    0.843    1.001    normal(0,10)
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   visual ~~                                                                    
    ##     textual           0.450    0.065    0.315    0.572    1.000     lkj_corr(1)
    ##     speed             0.463    0.087    0.290    0.627    1.000     lkj_corr(1)
    ##   textual ~~                                                                   
    ##     speed             0.279    0.072    0.132    0.416    1.000     lkj_corr(1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .x1                4.936    0.067    4.802    5.064    1.000    normal(0,32)
    ##    .x2                6.088    0.067    5.958    6.221    0.999    normal(0,32)
    ##    .x3                2.250    0.065    2.125    2.377    0.999    normal(0,32)
    ##    .x4                3.059    0.067    2.931    3.183    1.001    normal(0,32)
    ##    .x5                4.341    0.074    4.195    4.488    1.000    normal(0,32)
    ##    .x6                2.185    0.063    2.062    2.307    0.999    normal(0,32)
    ##    .x7                4.185    0.062    4.063    4.303    1.000    normal(0,32)
    ##    .x8                5.526    0.059    5.412    5.641    1.000    normal(0,32)
    ##    .x9                5.374    0.058    5.260    5.490    1.000    normal(0,32)
    ##     visual            0.000                                                    
    ##     textual           0.000                                                    
    ##     speed             0.000                                                    
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .x1                0.555    0.127    0.292    0.792    1.000 gamma(1,.5)[sd]
    ##    .x2                1.149    0.105    0.956    1.361    1.000 gamma(1,.5)[sd]
    ##    .x3                0.859    0.100    0.673    1.058    1.000 gamma(1,.5)[sd]
    ##    .x4                0.379    0.049    0.287    0.479    1.000 gamma(1,.5)[sd]
    ##    .x5                0.455    0.060    0.347    0.580    1.000 gamma(1,.5)[sd]
    ##    .x6                0.363    0.046    0.281    0.458    1.000 gamma(1,.5)[sd]
    ##    .x7                0.821    0.092    0.654    1.010    1.000 gamma(1,.5)[sd]
    ##    .x8                0.502    0.096    0.313    0.694    1.001 gamma(1,.5)[sd]
    ##    .x9                0.567    0.096    0.369    0.743    1.001 gamma(1,.5)[sd]
    ##     visual            1.000                                                    
    ##     textual           1.000                                                    
    ##     speed             1.000

Next, we will add all possible cross-loadings with a strong prior of
$N(0,\sigma = 0.08)$. The prior centers the loadings around 0 and allows
them little space to move.

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

    ## blavaan 0.5.9.1382 ended normally after 1000 iterations
    ## 
    ##   Estimator                                      BAYES
    ##   Optimization method                             MCMC
    ##   Number of model parameters                        48
    ## 
    ##   Number of observations                           301
    ## 
    ##   Statistic                                 MargLogLik         PPP
    ##   Value                                      -3858.976       0.124
    ## 
    ## Parameter Estimates:
    ## 
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   visual =~                                                                    
    ##     x1                0.766    0.098    0.581    0.967    1.001    normal(0,10)
    ##     x2                0.565    0.093    0.387    0.749    1.000    normal(0,10)
    ##     x3                0.768    0.096    0.585    0.961    1.000    normal(0,10)
    ##   textual =~                                                                   
    ##     x4                0.985    0.064    0.862    1.112    1.000    normal(0,10)
    ##     x5                1.156    0.071    1.021    1.300    1.000    normal(0,10)
    ##     x6                0.893    0.060    0.778    1.013    1.001    normal(0,10)
    ##   speed =~                                                                     
    ##     x7                0.724    0.085    0.561    0.894    1.000    normal(0,10)
    ##     x8                0.795    0.083    0.635    0.963    1.003    normal(0,10)
    ##     x9                0.542    0.072    0.404    0.690    1.000    normal(0,10)
    ##   visual =~                                                                    
    ##     x4                0.030    0.059   -0.083    0.145    1.000   normal(0,.08)
    ##     x5               -0.073    0.062   -0.196    0.044    1.000   normal(0,.08)
    ##     x6                0.061    0.055   -0.047    0.171    1.000   normal(0,.08)
    ##     x7               -0.128    0.064   -0.256   -0.000    1.000   normal(0,.08)
    ##     x8               -0.007    0.067   -0.139    0.122    1.000   normal(0,.08)
    ##     x9                0.192    0.059    0.075    0.309    0.999   normal(0,.08)
    ##   textual =~                                                                   
    ##     x1                0.107    0.067   -0.028    0.235    1.000   normal(0,.08)
    ##     x2                0.005    0.060   -0.118    0.123    1.000   normal(0,.08)
    ##     x3               -0.088    0.062   -0.211    0.028    1.000   normal(0,.08)
    ##     x7                0.017    0.062   -0.105    0.138    1.001   normal(0,.08)
    ##     x8               -0.038    0.061   -0.158    0.084    1.001   normal(0,.08)
    ##     x9                0.033    0.054   -0.072    0.137    1.000   normal(0,.08)
    ##   speed =~                                                                     
    ##     x1                0.041    0.065   -0.084    0.163    0.999   normal(0,.08)
    ##     x2               -0.048    0.064   -0.172    0.076    0.999   normal(0,.08)
    ##     x3                0.027    0.065   -0.102    0.154    0.999   normal(0,.08)
    ##     x4               -0.008    0.055   -0.117    0.102    1.000   normal(0,.08)
    ##     x5                0.006    0.060   -0.115    0.120    1.000   normal(0,.08)
    ##     x6               -0.002    0.052   -0.106    0.099    1.001   normal(0,.08)
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##   visual ~~                                                                    
    ##     textual           0.380    0.093    0.186    0.553    1.001     lkj_corr(1)
    ##     speed             0.354    0.109    0.132    0.555    1.000     lkj_corr(1)
    ##   textual ~~                                                                   
    ##     speed             0.256    0.106    0.039    0.455    1.001     lkj_corr(1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .x1                4.936    0.068    4.808    5.071    0.999    normal(0,32)
    ##    .x2                6.088    0.068    5.952    6.224    0.999    normal(0,32)
    ##    .x3                2.250    0.066    2.122    2.379    1.000    normal(0,32)
    ##    .x4                3.062    0.067    2.925    3.192    1.000    normal(0,32)
    ##    .x5                4.341    0.075    4.191    4.485    1.000    normal(0,32)
    ##    .x6                2.186    0.063    2.062    2.313    0.999    normal(0,32)
    ##    .x7                4.184    0.065    4.060    4.310    1.000    normal(0,32)
    ##    .x8                5.527    0.060    5.411    5.637    1.000    normal(0,32)
    ##    .x9                5.374    0.059    5.260    5.491    0.999    normal(0,32)
    ##     visual            0.000                                                    
    ##     textual           0.000                                                    
    ##     speed             0.000                                                    
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper     Rhat    Prior       
    ##    .x1                0.673    0.106    0.458    0.879    1.000 gamma(1,.5)[sd]
    ##    .x2                1.090    0.107    0.895    1.313    1.000 gamma(1,.5)[sd]
    ##    .x3                0.719    0.111    0.498    0.935    1.001 gamma(1,.5)[sd]
    ##    .x4                0.387    0.050    0.292    0.492    0.999 gamma(1,.5)[sd]
    ##    .x5                0.411    0.062    0.293    0.538    0.999 gamma(1,.5)[sd]
    ##    .x6                0.373    0.044    0.290    0.463    0.999 gamma(1,.5)[sd]
    ##    .x7                0.715    0.095    0.533    0.899    1.001 gamma(1,.5)[sd]
    ##    .x8                0.432    0.094    0.239    0.608    1.002 gamma(1,.5)[sd]
    ##    .x9                0.589    0.068    0.466    0.735    1.000 gamma(1,.5)[sd]
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
predictive p-value goes from ppp = 0 to ppp = 0.124, and is not that the
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
Molenaar, 265:255–63. Springer Proceedings in Mathematics & Statistics.
New York, NY, US: Springer.
<https://doi.org/10.1007/978-3-030-01310-3_23>.

Muthén, Bengt, and Tihomir Asparouhov. 2012. “Bayesian Structural
Equation Modeling: A More Flexible Representation of Substantive
Theory.” *Psychological Methods* 17 (3): 313–35.
<https://doi.org/10.1037/a0026802>.
