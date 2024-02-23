# Specify package imports
#' @import collapse
#' @import dplyr
#' @rawNamespace import(stats, except = c(filter, lag, D))
#' @rawNamespace import(data.table, except = c(first, last, between, fdroplevels))
NULL

#-----

.onLoad <- function (libname, pkgname) {

  # Create default option value for number of cores
  options(fusionModel.cores = max(1L, parallel::detectCores() - 1L))

  # Package startup message
  packageStartupMessage("fusionModel v", utils::packageVersion("fusionModel"), " | https://github.com/ummel/fusionModel")

}
