# Specify package imports
#' @import dplyr
#' @rawNamespace import(stats, except = c(filter, lag))
#' @rawNamespace import(data.table, except = c(first, last, between))
NULL

#-----

.onLoad <- function (libname, pkgname) {

  packageStartupMessage("fusionModel v", utils::packageVersion("fusionModel"), " | https://github.com/ummel/fusionModel")

}
