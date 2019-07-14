.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("This is ",paste(pkgname, version))
    packageStartupMessage("Major changes from blavaan 0.3-4:\n  - The default target is now 'stan' instead of 'jags'.\n  - To use JAGS, you now need to explicitly set target = 'jags'.\n  - To use the old Stan approach, set target = 'stanclassic'.\n  - Priors on variance parameters now default to gamma(1,.5) on the SD.")
}

.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}
