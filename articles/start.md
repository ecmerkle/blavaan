# Getting Started with blavaan

This page contains some tips for getting started with the *blavaan*
package.

### Installation

*blavaan* can be installed from CRAN in the usual way:

``` r
install.packages("blavaan")
```

In some situations, you may wish to install *blavaan* from GitHub. The
GitHub version sometimes contains bug fixes that are not yet on CRAN,
though it can also be less stable. To install from GitHub, use the
following command.

``` r
remotes::install_github("ecmerkle/blavaan", INSTALL_opts = "--no-multiarch")
```

This command requires that your system can compile Stan models, which is
not guaranteed if you usually install *blavaan* from CRAN. If you are
having trouble, it may help to look at the [RStan Getting Started
page.](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)

### Commands and Syntax

The *blavaan* package depends on the *lavaan* package for model
specification and for some computations. This means that, if you already
know *lavaan*, then you should already be able to do many things in
*blavaan*. In particular, many *blavaan* commands add the letter “b” to
the start of the *lavaan* command. For example,
[`sem()`](https://rdrr.io/pkg/lavaan/man/sem.html) becomes
[`bsem()`](http://ecmerkle.github.io/blavaan/reference/bsem.md), and
[`lavInspect()`](https://rdrr.io/pkg/lavaan/man/lavInspect.html) becomes
[`blavInspect()`](http://ecmerkle.github.io/blavaan/reference/blavInspect.md).
It is also sometimes possible to use a *lavaan* command on a *blavaan*
object, though the results may not always be what you expect.

With these details in mind, look at the [lavaan
tutorial](https://lavaan.ugent.be/tutorial/index.html) for many examples
of models. You can translate many of those examples to *blavaan* by
adding a “b” to the start of the commands. And look at the other pages
here, to learn about the additional *blavaan* arguments that are
specific to Bayesian methods.
