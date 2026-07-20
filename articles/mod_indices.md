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
extent to which the model’s chi-square ($`\chi^2`$) test statistic would
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
[`ppmc()`](https://blavaan.org/reference/ppmc.md) function of *blavaan*.
With this function, the MI and SEPC are computed for each posterior
sample, leading to posterior distributions for each of them.

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
    ## visual=~x9 35.042 35.026 35.775 10.815 17.976 53.079                   0.012
    ## x7~~x8     33.204 35.462 40.307 14.081  7.722 54.138                   0.071
    ## x8~~x9     27.097 12.622  2.792 40.934  0.000 69.737                   0.300
    ## x4~~x6     21.598  6.976  2.497 40.440  0.000 57.579                   0.449
    ## visual=~x7 18.375 16.671 12.465  9.844  3.505 32.370                   0.014
    ##            PPP_sim_LessThan_obs
    ## visual=~x9                0.988
    ## x7~~x8                    0.929
    ## x8~~x9                    0.700
    ## x4~~x6                    0.551
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
    ## x7~~x8      33.204 35.462 40.307 14.081  7.722 54.138                   0.071
    ## visual=~x9  35.042 35.026 35.775 10.815 17.976 53.079                   0.012
    ## visual=~x7  18.375 16.671 12.465  9.844  3.505 32.370                   0.014
    ## x8~~x9      27.097 12.622  2.792 40.934  0.000 69.737                   0.300
    ## textual=~x1 11.256 10.215  6.573  8.260  0.000 22.667                   0.222
    ##             PPP_sim_LessThan_obs
    ## x7~~x8                     0.929
    ## visual=~x9                 0.988
    ## visual=~x7                 0.986
    ## x8~~x9                     0.700
    ## textual=~x1                0.778

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
    ##               EAP Median   MAP    SD  lower upper PPP_sim_GreaterThan_obs
    ## x7~~x8      0.801  0.784 0.724 0.336  0.488 1.230                   0.036
    ## visual=~x9  0.517  0.498 0.468 0.132  0.336 0.691                   0.008
    ## textual=~x1 0.267  0.297 0.319 0.186 -0.006 0.505                   0.131
    ## x1~~x9      0.246  0.245 0.240 0.035  0.199 0.299                   0.020
    ## x2~~x3      0.220  0.221 0.217 0.037  0.163 0.274                   0.026
    ##             PPP_sim_LessThan_obs
    ## x7~~x8                     0.964
    ## visual=~x9                 0.992
    ## textual=~x1                0.869
    ## x1~~x9                     0.980
    ## x2~~x3                     0.974

Here we see that for the 2 highest parameters, the likely SEPC is x7x8 =
0.800788445388977 and visual=~x9 = 0.516773687948915. With this
information we can decide to include one of these new parameters in the
model (one at the time). For this example, because factor loadings have
a larger impact on the model-implied covariance matrix, I would choose
**visual=~x9**

``` r

HS.model <- ' visual  =~ x1 + x2 + x3 + x9
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit2 <- bcfa(HS.model, data=HolzingerSwineford1939, std.lv=TRUE)
```

And you can check if the added parameter has the expected impact on
overall fit with the
[`blavFitIndices()`](https://blavaan.org/reference/blavFitIndices.md)
and the [`summary()`](https://rdrr.io/r/base/summary.html) functions.

It is important to consider also the theoretical relevance of the
suggested parameters, and to ensure that they make sense, instead of
just adding parameters until having **good** fit.

### Summary

You can see more details about the application an test of these indices
in Bayesian SEM in Garnier-Villarreal and Jorgensen (2024).

In this tutorial we show how to calculate the MI and SEPC across
posterior distributions, and evaluate which parameters can be added.

With the [`ppmc()`](https://blavaan.org/reference/ppmc.md) function we
are able to calculate relevant information after model estimation, and
build posterior distributions of them.

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
