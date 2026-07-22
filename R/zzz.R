.onLoad <- function(libname, pkgname) {

  # Set default option for available processing cores safely
  cores <- parallel::detectCores()
  if (is.na(cores) || cores < 1L) {
    cores <- 1L
  } else {
    cores <- max(1L, cores - 1L)
  }

  options(fusionModel.cores = cores)

  invisible()
}

.onAttach <- function(libname, pkgname) {

  # Package startup message displayed only when attached via library()
  packageStartupMessage(
    "fusionModel v", utils::packageVersion("fusionModel"),
    " | https://github.com/ummel/fusionModel"
  )

  invisible()
}
