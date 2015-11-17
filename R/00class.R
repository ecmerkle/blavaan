#setOldClass("runjags")
setClass("blavaan",
    #slots = c(
    #  runjags     = "runjags"             # output from run.jags()
    #),
    contains = "lavaan"
)
