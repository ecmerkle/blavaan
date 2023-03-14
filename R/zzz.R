.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("This is ",paste(pkgname, version))
    packageStartupMessage('On multicore systems, we suggest use of future::plan("multicore") or\n', '  future::plan("multisession") for faster post-MCMC computations.')
}

