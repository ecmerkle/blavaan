library("testthat")
library("blavaan")

## test only if it is a development version
if (length(strsplit(packageDescription("blavaan")$Version, "\\.")[[1]]) > 2) {
  test_check("blavaan")
}
