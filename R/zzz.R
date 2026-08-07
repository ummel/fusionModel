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
  current_ver <- utils::packageVersion("fusionModel")

  # Always print package info first
  packageStartupMessage("fusionModel v", current_ver, " | https://github.com/ummel/fusionModel")

  # Direct raw URL check (bypasses GitHub API rate limits)
  remote_ver <- tryCatch({
    raw_url <- "https://raw.githubusercontent.com/ummel/fusionModel/master/DESCRIPTION"

    # Set a strict 2-second timeout to prevent slowing down package load on slow connections
    opts <- options(timeout = 2)
    on.exit(options(opts), add = TRUE)

    desc_lines <- readLines(raw_url, warn = FALSE)
    ver_line <- grep("^Version:", desc_lines, value = TRUE)

    if (length(ver_line) > 0) {
      package_version(trimws(sub("^Version:", "", ver_line[1])))
    } else {
      NULL
    }
  }, error = function(e) NULL)

  # Alert based on remote version check
  if (!is.null(remote_ver)) {
    if (remote_ver > current_ver) {
      packageStartupMessage("A newer version (v", remote_ver, ") is available. Run devtools::install_github('ummel/fusionModel') to upgrade.")
    } else {
      packageStartupMessage("You are using the latest version.")
    }
  }
}
