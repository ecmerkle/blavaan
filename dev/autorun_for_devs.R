# Pretty error happening
library(rlang)
rlang::global_handle()

# rlang::last_trace()
# Could always reread https://adv-r.hadley.nz/debugging.html#debugging and https://rstudio.github.io/r-manuals/r-exts/Debugging.html#debugging-r-code
library(lobstr)

# lobstr::tree()


# library(cli)
cli::pretty_print_code()

# To make things go fast
# options(blavaan.target = "cmdstan")

# First time loading
#' pkgload::load_all(".", export_all = FALSE, compile = TRUE, debug = FALSE)

# Already loaded
pkgload::load_all(".", export_all = TRUE, debug = TRUE, compile = FALSE)


