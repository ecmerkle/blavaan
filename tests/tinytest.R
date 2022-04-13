## test only if it is a development version
if (length(strsplit(packageDescription("blavaan")$Version, "\\.")[[1]]) > 2 &
    requireNamespace("tinytest", quietly=TRUE)) {
  tinytest::test_package("blavaan")
}