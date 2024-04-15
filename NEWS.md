
# Version 0.5-4
## New features
* New function sampleData() for generating data from a blavaan object.

* Functionality for the newdata argument in blavPredict(), which can generate lv (and other) predictions for new data from a model that has already been estimated (continuous data only, target = "stan" only).

* Refinements to two-level models (handle missingness via listwise deletion; lvs only at one level; improved messaging)

* Bugs from 0.5-3 are fixed.

# Version 0.5-3
## New features
* Functionality to find unrestricted blocks of the model's psi matrix (lv covariance matrix). lkj priors are assigned to these unrestricted blocks, improving the positive definite issue described in the "Opaque priors" paper.

* Improved functionality for obtaining posterior modes via, e.g., summary(., postmode = TRUE)

* blavCompare() messaging is improved to clarify ELPD differences, and the function returns more output.

* Bug fix in two-level models with within-only observed variables, messaging added for unstable ppp.

* When extracting posterior draws via blavInspect(., "mcmc"), column names now match lavaan parameter names. For old behavior involving Stan parameter names, use argument add.labels = FALSE

* Bugs from 0.5-2 are fixed.

## Bugs/glitches discovered after the release:
* Some models with exogenous covariates, fixed.x=TRUE, and missing data fail to converge and yield implausible parameter values (reported by DeAnne Hunter).

# Version 0.5-2
## New features
* This is a maintenance release, primarily adding the new array declaration syntax in Stan models (syntax that became available in the new version of rstan).

## Bugs/glitches discovered after the release:
* blavCompare() does not work with models that have meanstructure = FALSE (reported by Pedro Ribeiro).

* for target="jags", posterior modes cannot be obtained via postmode = TRUE (reported by Giada Venaruzzo).

* models with both continuous and ordinal variables fail for cases where all ordinal variables are missing (reported by Sonja Winter).

* certain equality constraints involving named parameters fail for target="stan" (reported by Niels Skovgaard-Olsen)

# Version 0.5-1
## New features
* Two-level models are now supported (for complete, continuous data) via the cluster argument.

## Bugs/glitches discovered after the release:
* For two-level model specification, the levels have to be labeled "within" and "between". This is more restrictive than lavaan specification.

* For target="jags", latent variable extraction via blavInspect(, "lvs") fails (reported by Joseph Saraceno).

# Version 0.4-8
## New features
* This is a maintenance release with bug fixes and some changes in compiler settings

## Bugs/glitches discovered after the release:
* For certain models with residual correlations and/or correlated factors, the initial values under target='stan' lead to non-positive definite matrices (reported by Yuanyuan Hu).

* For models where a latent variable is regressed on an observed variable (lv ~ ov), the latent variable samples do not account for the mean of the observed variable (they are centered around 0 and off by a constant).


# Version 0.4-7
## New features
* This is primarily an update to address a C++14 vs C++17 compilation issue identified by CRAN

* But bugs from 0.4-6 have also been fixed

## Bugs/glitches discovered after the release:
* Sampling from the priors (prisamp = TRUE) fails for models with meanstructure = FALSE; the posterior is still estimated (reported by Armel Brizuela Rodríguez).

* For target = "jags", models with a single-indicator latent variable, where the latent variable is regressed on other variables, return incorrect parameter estimates (reported by Brad Cosentino).


# Version 0.4-6
## New features
* For target = "stan", meanstructure=FALSE is allowed, along with use of sample.cov and sample.nobs instead of raw data

* Users are warned about priors on covariance matrices that are neither diagonal nor unrestricted

* For models where observed variable intercepts appear in the latent intercept vector (alpha), default priors come from the observed intercept vector nu (as the user would expect)

* inits = "simple" is now default (instead of "prior"), to address some convergence problems

* For stan targets, ":=" can now be used as an identity function

* For target = "stan", fix the missing data issue from 0.4-3 (complete data in one group but not the other)

* Column names are added to blavPredict(, type="lv")

## Bugs/glitches discovered after the release:
* blavFitIndices() and save.lvs = TRUE do not work correctly for models without meanstructure. Workaround is to use meanstructure = TRUE in the model estimation command (reported by Charles Hofacker).

* The lavaan summary() method is sometimes called instead of the blavaan summary() method (reported by multiple users, with Shu Fai Cheung providing helpful examples).



# Version 0.4-3
## New features
* For target = "stan", most models should run faster than they did in earlier versions (use of sufficient statistics)

* Posterior summaries are faster for ordinal models (using mnormt::sadmvn() by default)

* Variational Bayes option added: target="vb", which uses rstan::vb()

* cmdstanr functionality added: target="cmdstanr", which uses the model from target="stan"

* Fix blavInspect(., "lvs"/"lvmeans") for multiple groups + missing data

* Fixes to ppmc() for ordinal models; blavFitIndices() turned off for ordinal models (more research needed)

* loo() moment matching available by passing mcmcextra = list(data = list(moment_match_k_threshold))

## Bugs/glitches discovered after the release:
* target = "stan" fails when there are complete data in one group and missing data in another group (reported by Ronja Runge).

* blavPredict(, type="ymis") still not available for models with ordinal variables



# Version 0.4-1
## New features
* Functionality for ordinal observed variables is now available.

* For models with missing data, posterior summaries have been sped up (log-likelihood computations now done in Stan).

##  Bugs/glitches discovered after the release:
* blavPredict(, type="ymis") is not working for models with ordinal variables

* blavInspect(, 'lvs') or (, 'lvmeans') can fail for models with a combination of multiple groups, missing values, and excluded cases

* blavFitIndices() and ppmc() are not working for models with ordinal variables, or may indicate excessively bad fit

* blavFitIndices(, rescale="mcmc") fails



# Version 0.3-18
## New features
* This version adds a reference to the new JSS paper, including DOI, and corrects an inconsistent version dependency. There are no other changes as compared to 0.3-17.


# Version 0.3-17
## New features
* This is a maintenance release to correct major bugs from the previous version.


# Version 0.3-16
## New features
* blavPredict() function added for predicting latent variables and missing data.

* Some posterior summaries are sped up. (and fitMeasures are available when test="none")

* bug fixes from the previous version.

##  Bugs/glitches discovered after the release:
* For certain models with missing data, ppp-values are incorrect (sometimes equaling 1.0).

* For target="stan", some multiple group models fail when some cases are missing all observed variables (reported by DeAnne Hunter).


# Version 0.3-15
## New features
* Added an S3 summary() method for ppmc

* Posterior intervals summary() bug is fixed

##  Bugs/glitches discovered after the release:
* The summary() method for ppmc() and fitIndices() does not always work correctly.

* A Jacobian was incorrect for target="stan", when (non-default) priors were placed on precisions or variances instead of on standard deviations. This could impact estimates of posterior variability (reported by Roy Levy).


# Version 0.3-14
* (version 0.3-13 violated a CRAN policy)

## New features
* This is a maintenance release in response to a change in package Matrix.

##  Bugs/glitches discovered after the release:
* Posterior intervals are NA in summary(). Workarounds are to use parameterEstimates() (intervals assuming posterior normality) or to compute them yourself using the posterior samples (`blavInspect(fit, "mcmc")')
	

# Version 0.3-12
* (version 0.3-11 failed Windows CRAN checks)

## New features
* vector values of wiggle.sd are allowed for different priors on approximate
      equality constraints

* logical argument "prisamp" added, for sampling from a model's prior

* for target="stan", lkj prior is used for unrestricted lv correlation matrices

* default priors for conditional approaches (targets jags and stanclassic) revert to being placed on precisions (as opposed to SDs), for improvement in sampling efficiency
	


# Version 0.3-10
## New features
* save.lvs=TRUE works for missing data under target="stan"

* new arguments "wiggle" and "wiggle.sd" for approximate equality constraints under target="stan"
   
## Bugs/glitches discovered after the release:
* plot labels for target="stan" are sometimes incorrect (displaying a parameter different from the panel label).

* complex equality constraints are sometimes ignored (target="jags" or "stanclassic")

* equality constraints with std.lv=TRUE sometimes fail (target="stan")

* placing priors on variances or precisions yields incorrect results (target="stan"; reported by Roy Levy)


# Version 0.3-9
## New features
* improvements to save.lvs=TRUE for target="stan".
   
* target="stancond" is added, which is an experimental, noncentered Stan approach.
   
* bug fixes for prior settings and std.lv in target="stan", and defined parameters.
   
## Bugs/glitches discovered after the release:
* For target="stan", problems with sampling lvs when there are multiple groups or missing data.

* Errors for blavCompare() and blavFitIndices() due to version updates of other packages.

* For target="stan", some models with std.lv=TRUE would not converge.


# Version 0.3-8
## New features
* post-estimation, posterior predictive computations are sped up considerably.

* 0.3-7 bugs fixed.

## Bugs/glitches discovered after the release:
* For target="stan" and std.lv=TRUE, estimation fails for certain (growth) models (reported by Mauricio Garnier-Villareal).

* Some defined variables fail for target="jags" and "stanclassic" (reported by Mariëlle Zondervan-Zwijnenburg).

* User-specified priors sometimes are placed on the wrong parameter, related to the 0.3-7 bug (reported by Mauricio Garnier-Villareal).

* The dpriors() issue from 0.3-3 remains.

# Version 0.3-7
## New features
* for target="stan", gamma priors can now be placed on user's choice of
      variances, standard deviations, or precisions.

* plot() now works uniformly across Stan and JAGS, relying on bayesplot.

* post-MCMC parallelization is now handled via future.apply package
      (requires an extra "plan" command from user, but works on windows).

* 0.3-6 bugs fixed.

## Bugs/glitches discovered after the release:
* blavInspect(, 'lvmeans') returns rows in the wrong order for target="stan" (reported by Mehdi Momen).

* User-specified priors sometimes are placed on the wrong parameter, for target="stan" (reported by Enrico Toffalini).

* The dpriors() issue from 0.3-3 remains.

# Version 0.3-6
## New features
* this fixes the stan plot bug from 0.3-5.

## Bugs/glitches discovered after the release:
* user-specified priors on correlation parameters are silently ignored for target="stan" (reported by James Uanhoro).

* save.lvs=TRUE does not work for target="stan" (reported by Mauricio Garnier-Villareal).

* The dpriors() issue from 0.3-3 remains.

# Version 0.3-5
## New features
* target="stan" is now the default, using a pre-compiled Stan model instead of "on the fly" code.

* ppmc() function added by Terrence Jorgensen, facilitating posterior predictive checks.

* default priors are changed from gamma on precisions to gamma on standard deviations.

## Bugs/glitches discovered after the release:
* The Stan plot method silently fails (reported by Matt Yalch).

* The dpriors() issue from 0.3-3 remains.

# Version 0.3-4
## New features
* Add function standardizedPosterior() for standardizing posterior draws.

* Turn off posterior modes for target="jags", due to conflict between current versions of runjags and modeest.

* Rearrange posterior predictive internals.

## Bugs/glitches discovered after the release:
* The dpriors() issue from 0.3-3 remains.

* For target="jags", lv means obtained from blavInspect() (via argument 'lvmeans') are incorrect. (reported by Mauricio Garnier-Villareal)

* Use of plot() with target="stan" causes problems for future blavInspect() calls.

# Version 0.3-3
## New features
* For convergence="auto", max time was previously 5 min (undocumented). It is now Inf.
	
* Axis labels (parameter names) are now more sensible on convergence plots.

* Relative effective sample size now used to compute loo/waic SEs, and some SEs are now returned via fitMeasures().

* Added unit testing via package testthat.

* Fixed bugs from 0.3-2 (with exception of identity assignments using ':=')

## Bugs/glitches discovered after the release:
* Use of 'dpriors()': some observed variable precisions assigned latent precision (ipsi) prior; some latent means assigned observed mean (nu) prior.

# Version 0.3-2
## New features
* Conditional (on latent variables) information criteria available when save.lvs = TRUE.

* Experimental function 'blavFitIndices()' added for Bayesian versions of SEM metrics, contributed by Terrence Jorgensen.

* blavaan "intelligently" chooses target, if either runjags or rstan (but not both) is installed.

* Fixed bugs from 0.3-1, especially related to missing data in Stan.

## Bugs/glitches discovered after the release:
* Errors for Stan models with std.lv=TRUE, and an observed variable regressed on a latent variable (reported by Bo Zhang).

* Error for identity assignments using ':=' (reported by Marco Tullio Liuzza).
	
* Explicitly adding the argument 'do.fit=TRUE' fails (reported by Esteban Montenegro).
	

# Version 0.3-1
## New features
* Stan export now supported; use target="stan".

* Improved handling of complex models, including growth/change models.

* Sampling of factor scores (lvs) available via 'save.lvs=TRUE'. Samples/means can be obtained by supplying arguments 'lvs' and 'lvmeans' to 'blavInspect()'.

* Fixed bugs from 0.2-4.

## Bugs/glitches discovered after the release:
* Errors for Stan models with missing data, when there are exogenous ("x") variables.

* Errors for multi-group Stan models with std.lv=TRUE.
