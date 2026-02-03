## We do not use "useDynLib" in NAMESPACE but control it in
## .onLoad() and .onUnload() manually.

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("Rgraphviz")

# .onLoad <- function(libname, pkgname) {
#   invisible(suppressPackageStartupMessages(
#     sapply(c("acopula","copula","doParallel","prodlim",
#              "quantreg","stats","survival"),
#            requireNamespace, quietly = TRUE)
#   ))
#
#
#
#   library.dynam("CopulaSCR", pkgname, libname)
# }
#
# .onUnload <- function(libpath) {
#   library.dynam.unload("CopulaSCR", libpath)
# }


