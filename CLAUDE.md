# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

blavaan is an R package for Bayesian latent variable analysis (CFA, SEM, growth curve models). It mirrors the
`lavaan` syntax/workflow but fits models via MCMC using Stan (default) or JAGS. It's a compiled package (Rcpp +
StanHeaders), so most workflows require a working C++ toolchain.

## Development workflow

There is no test runner script beyond standard R package tooling.

- **Load the package for interactive development** (see `dev/autorun_for_devs.R`):
  ```r
  pkgload::load_all(".", export_all = TRUE, compile = TRUE, debug = FALSE)   # first load, or after touching src/ or inst/stan
  pkgload::load_all(".", export_all = TRUE, debug = TRUE, compile = FALSE)   # subsequent reloads, R-only changes
  ```
  Recompiling the Stan model (`inst/stan/stanmarg.stan`) is slow; only pass `compile = TRUE` when C++/Stan sources
  changed.
- **Full build/check** (requires compilation):
  ```sh
  R CMD build .
  R CMD INSTALL --no-multiarch .
  R CMD check <tarball>
  ```
- **Run the full test suite** (tinytest, only runs for dev-version builds — see `tests/tinytest.R`):
  ```r
  tinytest::test_package("blavaan")
  ```
- **Run a single test file** (after `pkgload::load_all(".")`):
  ```r
  tinytest::run_test_file("inst/tinytest/tests.blavaan.R")
  ```
  Test files live in `inst/tinytest/` (not `tests/testthat/`). Several tests check `Sys.which("jags")` and skip
  JAGS-dependent assertions if JAGS isn't installed locally; `rjags`/`runjags` are `Suggests`, not hard
  dependencies. `cmdstanr` (also `Suggests`) is needed for `target = "cmdstan"` tests.
- CI (`.github/workflows/R-CMD-check.yaml`) runs `rcmdcheck` on macOS/Ubuntu with `rjags`/`runjags` installed as
  extra packages.
- Stray `*~` backup files and `.Rhistory` files scattered through the tree are editor artifacts, not part of the
  package — ignore them, don't try to reconcile them with the real source.
- **Versioning convention**: the last segment of `Version:` in `DESCRIPTION` tracks the commit count (e.g. the
  `1443` in `0.5-10.1443`). When making a commit, update that segment to the number the new commit will become
  (`git rev-list --count HEAD` computed before committing, plus 1, since `HEAD` doesn't yet include the commit
  being made) and include the `DESCRIPTION` change in that same commit. The next commit after that bumps it
  again by one, and so on.
- **Commit messages**: keep the `Co-Authored-By: Claude ...` trailer, but drop the `<noreply@anthropic.com>`
  email address — it's useless here.
- **When preparing a new release for CRAN**: check that the whitelisted `bcontrol` passthrough argument
  names (constants in `R/blav_utils.R`) still match the current `rstan::sampling()`/`stan()`/`vb()`,
  `cmdstanr`'s `$sample()`, and `run.jags()`/`autorun.jags()` signatures. These names are hardcoded rather
  than introspected at runtime, so they can silently drift out of date as those packages add, rename, or
  deprecate arguments.

## Architecture

### Entry points and dispatch

`blavaan()` in `R/blavaan.R` (~1300 lines) is the core function; `bcfa()`/`bsem()` (aliases of each other) and
`bgrowth()` at the bottom of that file are thin convenience wrappers with CFA/SEM/growth-specific defaults, all
delegating to `blavaan()`. The `target` argument selects the estimation backend and is threaded through nearly the
whole call: `"stan"` (default), `"cmdstan"`, `"jags"`, `"stancond"` (experimental), `"stanclassic"` (legacy).
`blavaan()` validates the target/argument combination up front (e.g. `convergence = "auto"` is JAGS-only,
`wiggle` isn't supported for `"stancond"`, multilevel `cluster` models need `"stan"`/`"cmdstan"`) before handing
off to model generation.

### Model translation pipeline

A fitted-syntax `lavaan` parameter table is turned into backend-specific model code/data. Each backend has its own
generator, all named `lav2<backend>`:

- `R/lav_export_stanmarg.R` → `lav2stanmarg()` — the **default** path for `target = "stan"`/`"cmdstan"`. Builds
  data/inits for the precompiled *marginal-likelihood* Stan program (latent variables integrated out), rather than
  generating fresh Stan code per model.
- `R/lav_export_stancond.R` → `lav2stancond()` — conditional (latent variables sampled explicitly) Stan approach;
  experimental.
- `R/lav_export_stanclassic.R` → `lav2stan()` — older approach that generates a full Stan program's source text per
  model.
- `R/lav_export_mcmc.R` → `lav2mcmc()` — JAGS/BUGS code generator for `target = "jags"`. Note: this file's Stan
  support is dead code — it predates rstan's 2023 array-syntax change and blavaan no longer uses it to produce
  Stan models.

Shared helpers used by the generators to build up the parameter table (adding phantom latent variables for
covariance parameters, assigning priors, generating init values) live in `R/set_partable.R`, `R/set_priors.R`,
`R/set_stanpars.R`, `R/set_stancovs.R`, and `R/set_inits.R`. Default priors per target come from `R/dpriors.R`.

### The precompiled Stan model

`inst/stan/stanmarg.stan` is the single Stan program backing the default `target = "stan"` path; it `#include`s
helper functions from `inst/stanfuns/*.stan` (e.g. `sem_lv.stan`, `sem_mean.stan`, `fill_lower.stan`). It's
compiled at package build time via `rstantools`/`Rcpp` into `src/stanExports_stanmarg.*`; `R/stanmodels.R` (marked
"Generated by rstantools. Do not edit by hand.") loads the compiled `stanmodel` object at package load. Changing
`.stan` sources requires a recompile (`compile = TRUE` in `load_all()`, or a fresh `R CMD INSTALL`) — R-level
changes alone won't pick them up. `R/stanmarg_data.R` builds the data list passed to this model.

### Object model and post-fit functionality

`blavaan` (S4, `R/00class.R`) extends `lavaan`'s S4 class, so most `lavaan` accessors/methods work on blavaan
fits for free; blavaan-specific behavior is layered on top:

- `R/blav_fit.R` — assembles the lavaan-compatible `Fit` slot from raw MCMC output.
- `R/blav_object_inspect.R`, `R/blav_object_methods.R` — `blavInspect()`/`blavTech()` and S4 methods
  (`summary`, `coef`, `show`, `predict`) specific to Bayesian fits.
- `R/blav_tidiers.R` — `tidy()`/`glance()` (broom-style) S3 methods.
- `R/blav_predict.R`, `R/postpred.R` — `blavPredict()` and posterior/prior predictive sampling
  (`sampleData()`).
- `R/ctr_ppmc.R` — `ppmc()`, posterior predictive model checking, plus the `blavPPMC` class and its
  summary/plot/hist/pairs methods.
- `R/ctr_bayes_fit.R` — despite the name, defines the `blavFitIndices` class and `blavFitIndices()` (Bayesian
  analogs of traditional SEM fit indices like RMSEA/CFI).
- `R/blav_fit_measures.R`, `R/margloglik.R` — fit measures and marginal log-likelihood computation, used by
  `R/blav_compare.R` (`blavCompare()`, `BF()`) for Bayesian model comparison.
- `R/blav_adapt_quad.R` — adaptive quadrature (used for ordinal-data likelihoods).
- `R/lvgqs.R` — latent-variable generated quantities.
- `R/jags2r.R` — converts JAGS/`runjags` MCMC output into blavaan's internal representation.

### Multi-target abstraction pattern

Because the same conceptual model must be expressed for up to five different backends, most non-trivial changes to
model-building logic (new syntax feature, new prior type, new parameter matrix) touch several files in parallel:
the relevant `lav_export_*.R` generator(s) for the affected target(s), plus the shared `set_*.R` helpers if the
parameter table structure changes. When modifying behavior, check whether it needs to be mirrored across
`stanmarg`, `stancond`, `stanclassic`, and `jags` (`lav2mcmc`) paths, or is genuinely specific to one target.
