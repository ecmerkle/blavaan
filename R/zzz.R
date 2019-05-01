.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("This is ",paste(pkgname, version))
    packageStartupMessage("Note: The default target is now 'stan', which may break code using old blavaan versions.\n (To use JAGS, you now need to explicitly state target = 'jags'.)")
}

