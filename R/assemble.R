#' Assemble microdata across surveys
#'
#' @description
#' For fusionACS usage only. Provides a safe and efficient way to assemble (merge) fused microdata across surveys to return a single data table with the requested variables. The requested variables can come from any fused (donor) survey and/or the American Community Survey (ACS). The necessary variables are automatically and efficiently read from the appropriate local file and safely merged on household and/or person ID variables, optionally collapsing or expanding records as necessary depending on the \code{respondent} argument. Assumes (and checks for) a local \code{/fusionData} directory with appropriate file structure and conventions.
#'
#' @param year Integer. One or more years for which to return results (i.e. the ACS recipient year).
#' @param var Character. Name of one or more variables to return. May contain household- and/or person-level variables. See Details.
#' @param respondent Character. Should household or person-level microdata be returned?
#' @param M Integer. The first \code{M} implicates are returned for fused variables. Set \code{M = Inf} to return all implicates. Ignored if \code{var} contains only ACS variables (i.e. no implicates)
#' @param cores Integer. Number of cores used by the \code{\link[fst]{fst-package}} when reading from disk.
#'
#' @details The \code{var} argument can contain a mix of household- and/or person-level variables. When \code{respondent = "household"}, the reference person (i.e. head of household) value is returned for any person-level variables. When \code{respondent = "person"}, the values of any household-level variables are replicated for each person in the household.
#' @details A warning is issued if any \code{var} cannot be located in the available local files.
#'
#' @return A keyed data table containing the following columns, in addition to the variables named in \code{var}:
#' \describe{
#'  \item{M}{Implicate number. See \code{\link{fuse}}.}
#'  \item{year}{Year of the ACS recipient microdata.}
#'  \item{hid}{ACS household ID using fusionACS convention.}
#'  \item{pid}{ACS person ID using fusionACS convention, if \code{respondent = "person"}.}
#'  \item{weight}{ACS microdata primary sample weight.}
#'  }
#'
#' @examples
#' # NOTE: Requires local /fusionData directory containing the necessary ACS and .fsd files
#' test <- assemble(year = 2018:2019,
#'                  var = c("dollarel", "hincp", "agep", "state"),
#'                  respondent = "household",
#'                  M = 1)
#' @export

#------------

# library(tidyverse)
# library(data.table)
# source("R/utils.R")
# setwd( "/home/kevin/Documents/Projects/fusionData")
#
# test <- assemble(year = 2015,
#                  #var = c("dollarel" ,"dollarng", "dollarfo" ,"dollarlp" ,"hincp", "pov_ratio", "agep", "ref_race5", "rac1p", "hisp", "state", "puma10"),
#                  var = c("hincp", "pov_ratio", "state", "puma10", "agep"),
#                  respondent = "household",
#                  M = 1,
#                  cores = 2)

#-----

assemble <- function(year,
                     var,
                     respondent = "household",
                     M = 1,
                     cores = 1) {

  # Check validity of the working directory path
  # Checks if "/fusionData" is part of the path, as this is required
  b <- strsplit(full.path(getwd()), .Platform$file.sep, fixed = TRUE)[[1]]
  i <- which(b == "fusionData")
  if (length(i) == 0) stop("'/fusionData' is not part of the working directory path; this is required.")
  dir <- paste(b[1:i], collapse = .Platform$file.sep)

  # respondent identifier ("H" or "P")
  rtype <- substring(toupper(respondent), 1, 1)

  # Initial argument check
  stopifnot({
    year >= 2000 & year %% 1 == 0
    is.character(var) & length(var)
    rtype %in% c("H", "P")
    M > 0 & M %% 1 == 0
    cores > 0 & cores %% 1 == 0
  })

  t0 <- Sys.time()
  fst::threads_fst(nr_of_threads = cores)
  setDTthreads(threads = cores)

  # Universal variables always returned in output, if possible
  uvar <- c('M', 'year', 'hid', 'pid', 'weight')

  hh <- rtype == "H"
  v <- setdiff(var, uvar)
  Mimp <- M

  result <- lapply(year, function(yr) {

    # Get file paths to both ACS and fused microdata for 'yr'
    pa <- list.files(file.path(dir, "survey-processed/ACS"), pattern = paste0(yr, "_._processed.fst"), recursive = TRUE, full.names = TRUE)
    pc <- list.files(file.path(dir, "survey-processed/ACS"), pattern = paste0(yr, "_._custom.fst"), recursive = TRUE, full.names = TRUE)
    pf <- rev(list.files(file.path(dir, "fusion"), pattern = paste0(yr, "_._fused.fsd"), recursive = TRUE, full.names = TRUE))
    fpaths <- c(pa, pc, pf)

    d <- data.table()
    for (x in fpaths) {

      a <- substring(basename(x), 1, 4) == "ACS_"
      r <- rev(strsplit(x, "_")[[1]])[[2]]
      xn <- names(fst::fst(x))
      keep <- setdiff(intersect(xn, v), names(d))  # Excludes any variables already in 'd'
      if (a & rtype == r) keep <- c('weight', keep)
      if (length(keep)) {
        # NOTE: This could probably be made faster in the case of implicates, since not all needed to be loaded initially
        dt <- fst::read_fst(x, columns = intersect(xn, c('M', 'year', 'hid', 'pid', keep)), as.data.table = TRUE)
        if ("M" %in% names(dt) & is.finite(Mimp)) dt <- dt[M <= Mimp]
      } else {
        dt <- data.table()
      }

      # If household-level data requested, merge the reference person value for any person-level variables
      # TO DO: Allow different summary metric besides just the reference person value
      if (hh & rtype != r & nrow(dt) > 0) {
        i <- !duplicated(dt$hid)  # Retains the reference person data; first entry within each 'hid'
        dt <- select(dt[i], -pid)
      }

      # NOTE: If person-level data is requested, the merge() below causes any household-level variables to be replicated for each person
      # Since it is an inner merge, persons not in the household data (i.e. group quarter population) are discarded
      # If no household variables need to be added, then the returned person microdata will include the group quarter individuals
      if (nrow(dt) > 0) {
        dt <- setkeyv(dt, cols = intersect(c('M', 'year', 'hid', 'pid'), names(dt)))         # Set keys for merge
        if (nrow(d) == 0) {
          d <- dt
        } else {
          d <- merge(d, dt, by = intersect(key(d), key(dt)), allow.cartesian = TRUE)
        }
      }
      rm(dt)

    }

    return(d)

  }) %>%
    rbindlist() %>%
    setcolorder(neworder = intersect(c(uvar, var), names(.))) %>%  # Order columns
    setkeyv(cols = intersect(c('M', 'year', 'hid', 'pid'), names(.)))  # Set keys

  # Any 'var' missing in the output? If so, report as warning.
  miss <- setdiff(v, names(result))
  if (length(miss)) warning("Could not locate the following variable(s): ", paste(miss, collapse = ", "))

  return(result)

}
