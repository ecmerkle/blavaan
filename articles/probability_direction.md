# Probability of Direction

### Introduction

The Probability of Direction (pd) is an index of effect existence,
ranging from 0% to 100%, representing the certainty with which an effect
goes in a particular direction (i.e., is positive or negative) (Makowski
et al. 2019). Beyond its simplicity of interpretation, understanding and
computation, this index also presents other interesting properties: *It
is independent from the model: It is solely based on the posterior
distributions and does not require any additional information from the
data or the model.* It is robust to the scale of both the response
variable and the predictors. \*It is strongly correlated with the
frequentist p-value, and can thus be used to draw parallels and give
some reference to readers non-familiar with Bayesian statistics.

Can be interpreted as the probability that a parameter (described by its
posterior distribution) is above or below a chosen cutoff, an explicit
hypothesis. It is mathematically defined as the proportion of the
posterior distribution that satisfies the specified hypothesis. Although
differently expressed, this index is fairly similar (i.e., is strongly
correlated) to the frequentist p-value.

### Probability of Direction (pd)

For this example we will use the Industrialization and Political
Democracy example (Bollen 1989). We will first estimate the latent
regression model

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

fit <- bsem(model, data=PoliticalDemocracy,
            std.lv=T, meanstructure=T)
```

We can then look at the overall model results with the `summary`
function, in this case we are also asking for the standardized
estimates, and $`R^2`$

``` r

summary(fit, standardize=T, rsquare=T)
```

    ## blavaan 0.6.1.1484 ended normally after 1000 iterations
    ## 
    ##   Estimator                                      BAYES
    ##   Optimization method                             MCMC
    ##   Number of model parameters                        42
    ##   Number of equality constraints                     4
    ## 
    ##   Number of observations                            75
    ## 
    ##   Statistic                                 MargLogLik         PPP
    ##   Value                                             NA       0.025
    ## 
    ## Parameter Estimates:
    ## 
    ## 
    ## Latent Variables:
    ##                    Estimate  Post.SD pi.lower pi.upper   Std.lv  Std.all
    ##   ind60 =~                                                              
    ##     x1                0.707    0.073    0.576    0.860    0.707    0.922
    ##     x2                1.540    0.147    1.265    1.850    1.540    0.973
    ##     x3                1.281    0.144    1.021    1.589    1.281    0.872
    ##   dem60 =~                                                              
    ##     y1         (a)    1.459    0.174    1.128    1.806    1.788    0.761
    ##     y2         (b)    1.740    0.222    1.325    2.177    2.133    0.588
    ##     y3         (c)    1.815    0.203    1.423    2.224    2.225    0.700
    ##     y4         (d)    1.950    0.194    1.587    2.347    2.391    0.792
    ##   dem65 =~                                                              
    ##     y5         (a)    1.459    0.174    1.128    1.806    2.303    0.812
    ##     y6         (b)    1.740    0.222    1.325    2.177    2.747    0.772
    ##     y7         (c)    1.815    0.203    1.423    2.224    2.866    0.836
    ##     y8         (d)    1.950    0.194    1.587    2.347    3.080    0.871
    ##      Rhat    Prior       
    ##                          
    ##     1.000    normal(0,10)
    ##     1.000    normal(0,10)
    ##     1.000    normal(0,10)
    ##                          
    ##     0.999    normal(0,10)
    ##     1.000    normal(0,10)
    ##     1.000    normal(0,10)
    ##     1.001    normal(0,10)
    ##                          
    ##     0.999                
    ##     1.000                
    ##     1.000                
    ##     1.001                
    ## 
    ## Regressions:
    ##                    Estimate  Post.SD pi.lower pi.upper   Std.lv  Std.all
    ##   dem60 ~                                                               
    ##     ind60             0.709    0.172    0.396    1.069    0.578    0.578
    ##   dem65 ~                                                               
    ##     ind60             0.246    0.181   -0.115    0.608    0.156    0.156
    ##     dem60             0.867    0.131    0.614    1.127    0.673    0.673
    ##      Rhat    Prior       
    ##                          
    ##     1.001    normal(0,10)
    ##                          
    ##     1.000    normal(0,10)
    ##     1.000    normal(0,10)
    ## 
    ## Covariances:
    ##                    Estimate  Post.SD pi.lower pi.upper   Std.lv  Std.all
    ##  .y1 ~~                                                                 
    ##    .y5                0.748    0.444   -0.040    1.694    0.748    0.296
    ##  .y2 ~~                                                                 
    ##    .y4                1.733    0.837    0.265    3.494    1.733    0.320
    ##    .y6                2.202    0.805    0.717    3.967    2.202    0.332
    ##  .y3 ~~                                                                 
    ##    .y7                1.332    0.689    0.137    2.878    1.332    0.312
    ##  .y4 ~~                                                                 
    ##    .y8                0.388    0.489   -0.566    1.391    0.388    0.121
    ##  .y6 ~~                                                                 
    ##    .y8                1.074    0.737   -0.259    2.672    1.074    0.273
    ##      Rhat    Prior       
    ##                          
    ##     1.000     lkj_corr(1)
    ##                          
    ##     1.000       beta(1,1)
    ##     1.000       beta(1,1)
    ##                          
    ##     1.001     lkj_corr(1)
    ##                          
    ##     1.000       beta(1,1)
    ##                          
    ##     1.000       beta(1,1)
    ## 
    ## Intercepts:
    ##                    Estimate  Post.SD pi.lower pi.upper   Std.lv  Std.all
    ##    .x1                5.054    0.090    4.881    5.234    5.054    6.597
    ##    .x2                4.791    0.186    4.421    5.159    4.791    3.025
    ##    .x3                3.555    0.171    3.229    3.900    3.555    2.421
    ##    .y1                5.452    0.273    4.921    5.993    5.452    2.321
    ##    .y2                4.241    0.424    3.410    5.053    4.241    1.170
    ##    .y3                6.552    0.376    5.813    7.291    6.552    2.061
    ##    .y4                4.443    0.359    3.738    5.158    4.443    1.472
    ##    .y5                5.128    0.337    4.478    5.803    5.128    1.807
    ##    .y6                2.960    0.422    2.187    3.794    2.960    0.832
    ##    .y7                6.182    0.404    5.373    6.976    6.182    1.803
    ##    .y8                4.032    0.419    3.207    4.864    4.032    1.141
    ##     ind60             0.000                               0.000    0.000
    ##    .dem60             0.000                               0.000    0.000
    ##    .dem65             0.000                               0.000    0.000
    ##      Rhat    Prior       
    ##     1.002    normal(0,32)
    ##     1.002    normal(0,32)
    ##     1.002    normal(0,32)
    ##     1.001    normal(0,32)
    ##     1.000    normal(0,32)
    ##     1.001    normal(0,32)
    ##     1.003    normal(0,32)
    ##     1.001    normal(0,32)
    ##     1.000    normal(0,32)
    ##     1.002    normal(0,32)
    ##     1.002    normal(0,32)
    ##                          
    ##                          
    ##                          
    ## 
    ## Variances:
    ##                    Estimate  Post.SD pi.lower pi.upper   Std.lv  Std.all
    ##    .x1                0.088    0.023    0.046    0.138    0.088    0.149
    ##    .x2                0.135    0.082    0.003    0.312    0.135    0.054
    ##    .x3                0.516    0.106    0.343    0.751    0.516    0.239
    ##    .y1                2.320    0.590    1.260    3.611    2.320    0.420
    ##    .y2                8.598    1.556    5.997   12.014    8.598    0.654
    ##    .y3                5.153    1.058    3.389    7.558    5.153    0.510
    ##    .y4                3.401    0.902    1.833    5.402    3.401    0.373
    ##    .y5                2.746    0.638    1.673    4.190    2.746    0.341
    ##    .y6                5.112    1.127    3.189    7.767    5.112    0.404
    ##    .y7                3.538    0.891    2.080    5.521    3.538    0.301
    ##    .y8                3.015    0.927    1.220    4.990    3.015    0.241
    ##     ind60             1.000                               1.000    1.000
    ##    .dem60             1.000                               0.665    0.665
    ##    .dem65             1.000                               0.401    0.401
    ##      Rhat    Prior       
    ##     1.002 gamma(1,.5)[sd]
    ##     1.000 gamma(1,.5)[sd]
    ##     1.000 gamma(1,.5)[sd]
    ##     1.000 gamma(1,.5)[sd]
    ##     1.001 gamma(1,.5)[sd]
    ##     1.001 gamma(1,.5)[sd]
    ##     1.000 gamma(1,.5)[sd]
    ##     1.000 gamma(1,.5)[sd]
    ##     1.000 gamma(1,.5)[sd]
    ##     1.001 gamma(1,.5)[sd]
    ##     0.999 gamma(1,.5)[sd]
    ##                          
    ##                          
    ##                          
    ## 
    ## R-Square:
    ##                    Estimate
    ##     x1                0.851
    ##     x2                0.946
    ##     x3                0.761
    ##     y1                0.580
    ##     y2                0.346
    ##     y3                0.490
    ##     y4                0.627
    ##     y5                0.659
    ##     y6                0.596
    ##     y7                0.699
    ##     y8                0.759
    ##     dem60             0.335
    ##     dem65             0.599

To calculate the probability of direction we will use a function from
the package `brms` (Bürkner 2017)

``` r

library(brms)
```

And we will need to extract the posterior draws as a matrix,

``` r

mc_out <- as.matrix(blavInspect(fit, "mcmc", add.labels = FALSE))
dim(mc_out)
```

    ## [1] 3000   42

``` r

colnames(mc_out)
```

    ##  [1] "ly_sign[1]"    "ly_sign[2]"    "ly_sign[3]"    "ly_sign[4]"   
    ##  [5] "ly_sign[5]"    "ly_sign[6]"    "ly_sign[7]"    "ly_sign[4]"   
    ##  [9] "ly_sign[5]"    "ly_sign[6]"    "ly_sign[7]"    "bet_sign[1]"  
    ## [13] "bet_sign[2]"   "bet_sign[3]"   "Theta_cov[1]"  "Theta_cov[2]" 
    ## [17] "Theta_cov[3]"  "Theta_cov[4]"  "Theta_cov[5]"  "Theta_cov[6]" 
    ## [21] "Theta_var[1]"  "Theta_var[2]"  "Theta_var[3]"  "Theta_var[4]" 
    ## [25] "Theta_var[5]"  "Theta_var[6]"  "Theta_var[7]"  "Theta_var[8]" 
    ## [29] "Theta_var[9]"  "Theta_var[10]" "Theta_var[11]" "Nu_free[1]"   
    ## [33] "Nu_free[2]"    "Nu_free[3]"    "Nu_free[4]"    "Nu_free[5]"   
    ## [37] "Nu_free[6]"    "Nu_free[7]"    "Nu_free[8]"    "Nu_free[9]"   
    ## [41] "Nu_free[10]"   "Nu_free[11]"

It is also important to note that the parameters in the posterior draws
are named after the `Stan` underlying object names, instead of the
`(b)lavaan` parameter names. This is due to the argument
`add.labels = FALSE` and is used here to avoid trouble with parameter
names that have tildes or equal signs in them. You can see what each
parameter name equates to with the
[`partable()`](https://rdrr.io/pkg/lavaan/man/parTable.html) function,
as follows

``` r

pt <- partable(fit)[,c("lhs","op","rhs","pxnames")]
pt
```

    ##      lhs op   rhs       pxnames
    ## 1  ind60 =~    x1    ly_sign[1]
    ## 2  ind60 =~    x2    ly_sign[2]
    ## 3  ind60 =~    x3    ly_sign[3]
    ## 4  dem60 =~    y1    ly_sign[4]
    ## 5  dem60 =~    y2    ly_sign[5]
    ## 6  dem60 =~    y3    ly_sign[6]
    ## 7  dem60 =~    y4    ly_sign[7]
    ## 8  dem65 =~    y5    ly_sign[4]
    ## 9  dem65 =~    y6    ly_sign[5]
    ## 10 dem65 =~    y7    ly_sign[6]
    ## 11 dem65 =~    y8    ly_sign[7]
    ## 12 dem60  ~ ind60   bet_sign[1]
    ## 13 dem65  ~ ind60   bet_sign[2]
    ## 14 dem65  ~ dem60   bet_sign[3]
    ## 15    y1 ~~    y5  Theta_cov[1]
    ## 16    y2 ~~    y4  Theta_cov[2]
    ## 17    y2 ~~    y6  Theta_cov[3]
    ## 18    y3 ~~    y7  Theta_cov[4]
    ## 19    y4 ~~    y8  Theta_cov[5]
    ## 20    y6 ~~    y8  Theta_cov[6]
    ## 21    x1 ~~    x1  Theta_var[1]
    ## 22    x2 ~~    x2  Theta_var[2]
    ## 23    x3 ~~    x3  Theta_var[3]
    ## 24    y1 ~~    y1  Theta_var[4]
    ## 25    y2 ~~    y2  Theta_var[5]
    ## 26    y3 ~~    y3  Theta_var[6]
    ## 27    y4 ~~    y4  Theta_var[7]
    ## 28    y5 ~~    y5  Theta_var[8]
    ## 29    y6 ~~    y6  Theta_var[9]
    ## 30    y7 ~~    y7 Theta_var[10]
    ## 31    y8 ~~    y8 Theta_var[11]
    ## 32 ind60 ~~ ind60          <NA>
    ## 33 dem60 ~~ dem60          <NA>
    ## 34 dem65 ~~ dem65          <NA>
    ## 35    x1 ~1          Nu_free[1]
    ## 36    x2 ~1          Nu_free[2]
    ## 37    x3 ~1          Nu_free[3]
    ## 38    y1 ~1          Nu_free[4]
    ## 39    y2 ~1          Nu_free[5]
    ## 40    y3 ~1          Nu_free[6]
    ## 41    y4 ~1          Nu_free[7]
    ## 42    y5 ~1          Nu_free[8]
    ## 43    y6 ~1          Nu_free[9]
    ## 44    y7 ~1         Nu_free[10]
    ## 45    y8 ~1         Nu_free[11]
    ## 46 ind60 ~1                <NA>
    ## 47 dem60 ~1                <NA>
    ## 48 dem65 ~1                <NA>
    ## 49  .p4. ==  .p8.          <NA>
    ## 50  .p5. ==  .p9.          <NA>
    ## 51  .p6. == .p10.          <NA>
    ## 52  .p7. == .p11.          <NA>

For this example we will focus on the regressions between factors

``` r

pt[pt$op=="~",]
```

    ##      lhs op   rhs     pxnames
    ## 12 dem60  ~ ind60 bet_sign[1]
    ## 13 dem65  ~ ind60 bet_sign[2]
    ## 14 dem65  ~ dem60 bet_sign[3]

Now, we can calculate pd, with the
[`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)
function from `brms`. We can ask specific questions of the posterior
distributions, for example if we want to know what proportion of the
regression `dem65~ind60` is higher than 0. The function requires 2
arguments, the posterior draws (`mc_out`) and a hypothesis
(`bet_sign[2] > 0`). We are also adding the \`\`alpha\`\`\` argument
that specifies the size for the credible intervals

``` r

hypothesis(mc_out, "bet_sign[2] > 0", alpha = 0.05)
```

    ## New names:
    ## • `ly_sign[4]` -> `ly_sign[4]...4`
    ## • `ly_sign[5]` -> `ly_sign[5]...5`
    ## • `ly_sign[6]` -> `ly_sign[6]...6`
    ## • `ly_sign[7]` -> `ly_sign[7]...7`
    ## • `ly_sign[4]` -> `ly_sign[4]...8`
    ## • `ly_sign[5]` -> `ly_sign[5]...9`
    ## • `ly_sign[6]` -> `ly_sign[6]...10`
    ## • `ly_sign[7]` -> `ly_sign[7]...11`

    ## Hypothesis Tests for class :
    ##          Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (bet_sign[2]) > 0     0.25      0.18    -0.04     0.55      10.49      0.91
    ##   Star
    ## 1     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

The estimate presents the mean of the posterior distribution, and the
respective measures of variability (deviation and credible interval).
`Post.Prob` is the pd under the stated hypothesis, so in this example we
can say that 91% of the posterior distribution of `dem65~ind60` is lower
than 0. This is equivalent to the one-tail test. And `Evid.Ratio` is the
evidence ratio for the hypothesis, when the hypothesis is of the form
$`a > b`$, the evidence ratio is the ratio of the posterior probability
of $`a > b`$ and the posterior probability of $`a < b`$

In another example, we want to know what proportion of the regression
`dem60~ind60` is higher than 0. Here we can see that 100% of the
posterior probability is higher than 0, in such a case
`Evid.Ratio = Inf`, this will happens when the whole distribution
fulfills the hypothesis.

``` r

hypothesis(mc_out, "bet_sign[1] > 0", alpha = 0.05)
```

    ## New names:
    ## • `ly_sign[4]` -> `ly_sign[4]...4`
    ## • `ly_sign[5]` -> `ly_sign[5]...5`
    ## • `ly_sign[6]` -> `ly_sign[6]...6`
    ## • `ly_sign[7]` -> `ly_sign[7]...7`
    ## • `ly_sign[4]` -> `ly_sign[4]...8`
    ## • `ly_sign[5]` -> `ly_sign[5]...9`
    ## • `ly_sign[6]` -> `ly_sign[6]...10`
    ## • `ly_sign[7]` -> `ly_sign[7]...11`

    ## Hypothesis Tests for class :
    ##          Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (bet_sign[1]) > 0     0.71      0.17     0.44     1.01        Inf         1
    ##   Star
    ## 1    *
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

In another possible case of interest, you could use this to test
equalities between parameters, for example we can test if `dem60~ind60`
is higher than `dem65~ind60`. Here we see 97% of the posteriors state
that `dem60~ind60` is higher than `dem65~ind60`, and the mean of the
difference (`dem60~ind60 - dem65~ind60`) is `Estimate=0.46`

``` r

hypothesis(mc_out, "bet_sign[1] - bet_sign[2] > 0", alpha = 0.05)
```

    ## New names:
    ## • `ly_sign[4]` -> `ly_sign[4]...4`
    ## • `ly_sign[5]` -> `ly_sign[5]...5`
    ## • `ly_sign[6]` -> `ly_sign[6]...6`
    ## • `ly_sign[7]` -> `ly_sign[7]...7`
    ## • `ly_sign[4]` -> `ly_sign[4]...8`
    ## • `ly_sign[5]` -> `ly_sign[5]...9`
    ## • `ly_sign[6]` -> `ly_sign[6]...10`
    ## • `ly_sign[7]` -> `ly_sign[7]...11`

    ## Hypothesis Tests for class :
    ##                 Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
    ## 1 (bet_sign[1]-bet_... > 0     0.46      0.26     0.04     0.91         29
    ##   Post.Prob Star
    ## 1      0.97    *
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

### Region of Practical Equivalence (ROPE)

Note that so far we have only tested the hypothesis against 0, which
would be equivalent to the frequentist null hypothesis tests. But we can
test against any other. Bayesian inference is not based on statistical
significance, where effects are tested against “zero”. Indeed, the
Bayesian framework offers a probabilistic view of the parameters,
allowing assessment of the uncertainty related to them. Thus, rather
than concluding that an effect is present when it simply differs from
zero, we would conclude that the probability of being outside a specific
range that can be considered as “practically no effect” (i.e., a
negligible magnitude) is sufficient. This range is called the region of
practical equivalence (ROPE).

Indeed, statistically, the probability of a posterior distribution being
different from 0 does not make much sense (the probability of it being
different from a single point being infinite). Therefore, the idea
underlining ROPE is to let the user define an area around the null value
enclosing values that are equivalent to the null value for practical
purposes (Kruschke and Liddell 2018)

For these examples, we would change the value tested, a common
recommendations is to use `|0.1|` as the minimally relevant value for
standardized regressions, in this case we find that `0.79` proportion of
the posterior is above `0.1`

``` r

hypothesis(mc_out, "bet_sign[2] > .1", alpha = 0.05)
```

    ## New names:
    ## • `ly_sign[4]` -> `ly_sign[4]...4`
    ## • `ly_sign[5]` -> `ly_sign[5]...5`
    ## • `ly_sign[6]` -> `ly_sign[6]...6`
    ## • `ly_sign[7]` -> `ly_sign[7]...7`
    ## • `ly_sign[4]` -> `ly_sign[4]...8`
    ## • `ly_sign[5]` -> `ly_sign[5]...9`
    ## • `ly_sign[6]` -> `ly_sign[6]...10`
    ## • `ly_sign[7]` -> `ly_sign[7]...11`

    ## Hypothesis Tests for class :
    ##               Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
    ## 1 (bet_sign[2])-(.1) > 0     0.15      0.18    -0.14     0.45       3.82
    ##   Post.Prob Star
    ## 1      0.79     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

### 89% vs. 95% CI

Most commonly and from the frequentist tradition you will see the use of
the 95% Credible interval. Using 89% is another popular choice, and used
to be the default for a long time. How did it start?

Naturally, when it came about choosing the CI level to report by
default, people started using 95%, the arbitrary convention used in the
frequentist world. However, some authors suggested that 95% might not be
the most appropriate for Bayesian posterior distributions, potentially
lacking stability if not enough posterior samples are drawn (McElreath
2020).

The proposition was to use 90% instead of 95%. However, recently,
McElreath (2020) suggested that if we were to use arbitrary thresholds
in the first place, why not use 89%? Moreover, 89 is the highest prime
number that does not exceed the already unstable 95% threshold. What
does it have to do with anything? Nothing, but it reminds us of the
total arbitrariness of these conventions (McElreath 2020).

You can use this as the argument `alpha` argument in the `hypothesis`
function, or as the interpretation values for `Post.Prob`

### Caveats

Although this allows testing of hypotheses in a similar manner as in the
frequentist null-hypothesis testing framework, we strongly argue against
using arbitrary cutoffs (e.g., p \< .05) to determine the ‘existence’ of
an effect.

ROPE is sensitive to scale, so be aware that the value of interest is
representative in the respective scale. For this, standardized
parameters are useful to have in a commonly used scale.

### References

Bollen, Kenneth A. 1989. *Structural Equations with Latent Variables*.
Wiley Series in Probability and Mathematical Statistics. John Wiley &
Sons, Inc.

Bürkner, Paul-Christian. 2017. “brms: An R Package for Bayesian
Multilevel Models Using Stan.” *Journal of Statistical Software* 80 (1):
1–28. <https://doi.org/10.18637/jss.v080.i01>.

Kruschke, John K., and Torrin M. Liddell. 2018. “The Bayesian New
Statistics: Hypothesis Testing, Estimation, Meta-Analysis, and Power
Analysis from a Bayesian Perspective.” *Psychonomic Bulletin & Review*
25 (1): 178–206. <https://doi.org/10.3758/s13423-016-1221-4>.

Makowski, Dominique, Mattan S. Ben-Shachar, S. H. Annabel Chen, and
Daniel Lüdecke. 2019. “Indices of Effect Existence and Significance in
the Bayesian Framework.” *Frontiers in Psychology* 10: 2767.
<https://doi.org/10.3389/fpsyg.2019.02767>.

McElreath, Richard. 2020. *Statistical Rethinking: A Bayesian Course
with Examples in R and Stan*. 2nd ed. CRC Texts in Statistical Science.
Taylor; Francis, CRC Press.
