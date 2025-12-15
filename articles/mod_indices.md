# Modification indices

### Introduction

In SEM, one of the first steps is to evaluate the model’s global fit.
After global fit, we need to evaluate the local fit of a model, meaning
how the model reproduces specific correlations between observed
variables.

There are a couple of common methods for this, (a) testing for high
residual correlations, or (b) modification indices. This tutorial
focuses on the second. Modification indices test the **likely** change
in the model fit if a single parameter is added to the model that was
not originally included. This test can be carried out for every possible
parameter that was not included (Bentler 1990).

### Modification Indices

Modification indices present different **indices** to quantify the
effect of each parameter, and we will focus on two here. These are (a)
the modification index (MI) or Lagrange multiplier, which estimates the
extent to which the model’s chi-square ($\chi^{2}$) test statistic would
decrease if a parameter were added to the model and freely estimated,
and (b) standardized expected parameter change (SEPC), which is the
approximated standardized value of the parameter if it were to be
estimated in the model (Whittaker 2012; Garnier-Villarreal and Jorgensen
2024).

MI presents the possible effect on the overall model, and SEPC presents
the effect size for the missed parameter.

We will show an example with the Holzinger and Swineford (1939) model.
You first estimate your SEM/CFA model as usual

``` r
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- bcfa(HS.model, data=HolzingerSwineford1939, std.lv=TRUE)
```

Then we would need to write a **discrepancy** function to collect the
modification indices. The list below contains two functions that
estimate and save the MI and SEPC.

``` r
discFUN <- list(
  mod.ind_mi = function(object){
    temp <- modificationindices(object, free.remove = F)
    mods <- temp$mi
    names(mods) <- paste0(temp$lhs, temp$op, temp$rhs)
    return(mods)
  },
  mod.ind_sepc.all = function(object){
    temp <- modificationindices(object, free.remove = F)
    sepc.all <- temp$sepc.all
    names(sepc.all) <- paste0(temp$lhs, temp$op, temp$rhs)
    return(sepc.all)
  }
)
```

Then we will pass this function to the
[`ppmc()`](http://ecmerkle.github.io/blavaan/reference/ppmc.md) function
of *blavaan*. With this function, the MI and SEPC are computed for each
posterior sample, leading to posterior distributions for each of them.

``` r
out <- ppmc(fit, discFUN = discFUN)
```

Then we view the top 5 parameters arrange by the posterior mean (EAP)
MI, which in this case shows that the parameter having the highest
impact in overall model fit (according to EAP) is **visual=~x9**, the
cross-loading from the Visual factor to item **x9**.

``` r
summary(out, prob=.9, discFUN = "mod.ind_mi", sort.by="EAP", decreasing=T)[1:5,]
```

    ## 
    ## Posterior summary statistics and highest posterior density (HPD) 90% credible intervals
    ##  for the posterior distribution of realized discrepancy-function values based on observed data, 
    ##  along with posterior predictive p values to test hypotheses in either direction:
    ## 
    ## 
    ##               EAP Median    MAP     SD  lower  upper PPP_sim_GreaterThan_obs
    ## visual=~x9 35.215 35.295 36.207 11.263 16.921 53.858                   0.017
    ## x7~~x8     33.485 35.939 39.163 14.513  6.690 52.589                   0.069
    ## x8~~x9     26.685 12.448  2.617 40.677  0.000 70.149                   0.319
    ## x4~~x6     20.637  7.156  2.110 34.437  0.000 54.681                   0.456
    ## visual=~x7 18.742 16.911 11.939 10.682  3.880 32.609                   0.014
    ##            PPP_sim_LessThan_obs
    ## visual=~x9                0.983
    ## x7~~x8                    0.931
    ## x8~~x9                    0.681
    ## x4~~x6                    0.544
    ## visual=~x7                0.986

But according to the posterior median, the parameter that would have the
highest impact would be the residual correlation between indicators
**x7** and **x8**

``` r
summary(out, prob=.9, discFUN = "mod.ind_mi", sort.by="Median", decreasing=T)[1:5,]
```

    ## 
    ## Posterior summary statistics and highest posterior density (HPD) 90% credible intervals
    ##  for the posterior distribution of realized discrepancy-function values based on observed data, 
    ##  along with posterior predictive p values to test hypotheses in either direction:
    ## 
    ## 
    ##                EAP Median    MAP     SD  lower  upper PPP_sim_GreaterThan_obs
    ## x7~~x8      33.485 35.939 39.163 14.513  6.690 52.589                   0.069
    ## visual=~x9  35.215 35.295 36.207 11.263 16.921 53.858                   0.017
    ## visual=~x7  18.742 16.911 11.939 10.682  3.880 32.609                   0.014
    ## x8~~x9      26.685 12.448  2.617 40.677  0.000 70.149                   0.319
    ## textual=~x1 10.991  9.971  9.913  8.125  0.000 22.101                   0.219
    ##             PPP_sim_LessThan_obs
    ## x7~~x8                     0.931
    ## visual=~x9                 0.983
    ## visual=~x7                 0.986
    ## x8~~x9                     0.681
    ## textual=~x1                0.781

The MI is still recommended as the best metric to indicate which
parameter is best to include next, and we can use the SEPC to evaluate
the **likely** effect size for the respective parameters.

``` r
summary(out, prob=.9, discFUN = "mod.ind_sepc.all", sort.by="EAP", decreasing=T)[1:5,]
```

    ## 
    ## Posterior summary statistics and highest posterior density (HPD) 90% credible intervals
    ##  for the posterior distribution of realized discrepancy-function values based on observed data, 
    ##  along with posterior predictive p values to test hypotheses in either direction:
    ## 
    ## 
    ##               EAP Median   MAP    SD lower upper PPP_sim_GreaterThan_obs
    ## x7~~x8      0.807  0.786 0.746 0.369 0.481 1.226                   0.043
    ## visual=~x9  0.521  0.498 0.483 0.133 0.335 0.703                   0.009
    ## textual=~x1 0.271  0.297 0.323 0.165 0.044 0.508                   0.129
    ## x1~~x9      0.247  0.246 0.246 0.040 0.198 0.299                   0.022
    ## x2~~x3      0.222  0.222 0.218 0.040 0.174 0.280                   0.023
    ##             PPP_sim_LessThan_obs
    ## x7~~x8                     0.957
    ## visual=~x9                 0.991
    ## textual=~x1                0.871
    ## x1~~x9                     0.978
    ## x2~~x3                     0.977

Here we see that for the 2 highest parameters, the likely SEPC is
x7\~~x8 = 0.807250019181851 and visual=~x9 = 0.520720661163884. With
this information we can decide to include one of these new parameters in
the model (one at the time). For this example, because factor loadings
have a larger impact on the model-implied covariance matrix, I would
choose **visual=~x9**

``` r
HS.model <- ' visual  =~ x1 + x2 + x3 + x9
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, std.lv=TRUE)
```

And you can check if the added parameter has the expected impact on
overall fit with the
[`blavFitIndices()`](http://ecmerkle.github.io/blavaan/reference/blavFitIndices.md)
and the [`summary()`](https://rdrr.io/r/base/summary.html) functions.

It is important to consider also the theoretical relevance of the
suggested parameters, and to ensure that they make sense, instead of
just adding parameters until having **good** fit.

### Summary

You can see more details about the application an test of these indices
in Bayesian SEM in Garnier-Villarreal and Jorgensen (2024).

In this tutorial we show how to calculate the MI and SEPC across
posterior distributions, and evaluate which parameters can be added.

With the [`ppmc()`](http://ecmerkle.github.io/blavaan/reference/ppmc.md)
function we are able to calculate relevant information after model
estimation, and build posterior distributions of them.

The general recommendations are to use MI to identify the most likely
parameter to add, and SEPC as the effect size of the new parameter
(Garnier-Villarreal and Jorgensen 2024).

### References

Bentler, P. M. 1990. “Fit Indexes, Lagrange Multipliers, Constraint
Changes and Incomplete Data in Structural Models.” *Multivariate
Behavioral Research* 25 (2): 163–72.
<https://doi.org/10.1207/s15327906mbr2502_3>.

Garnier-Villarreal, Mauricio, and Terrence D Jorgensen. 2024.
“Evaluating Local Model Misspecification with Modification Indices in
Bayesian Structural Equation Mo.” *Structural Equation Modeling: A
Multidisciplinary Journal*, 1–15.
<https://doi.org/10.1080/10705511.2024.2413128>.

Holzinger, K. J., and F. A. Swineford. 1939. *A Study of Factor
Analysis: The Stability of a Bi-Factor Solution*. Supplementary
Educational Monograph 48. Chicago: University of Chicago Press.

Whittaker, Tiffany A. 2012. “Using the Modification Index and
Standardized Expected Parameter Change for Model Modification.” *The
Journal of Experimental Education* 80 (1): 26–44.
<https://doi.org/10.1080/00220973.2010.531299>.
