#' Assemble fusionACS microdata across surveys
#'
#' @description
#' For fusionACS usage only. Provides a safe and efficient way to assemble (merge) fused microdata across surveys to return a single data table with the requested variables. The requested variables can come from any fused (donor) survey and/or the American Community Survey (ACS). The necessary variables are automatically and efficiently read from the appropriate local file and safely merged on household and/or person ID variables, optionally collapsing or expanding records as necessary depending on the \code{respondent} argument. Assumes (and checks for) a local \code{/fusionData} directory with appropriate file structure and conventions.
#'
#' @param year Integer. One or more years for which to return results (i.e. the ACS recipient year).
#' @param var Character. Name of one or more variables to return. May contain household- and/or person-level variables. See Details.
#' @param respondent Character. Should \code{"household"}- or \code{"person"}-level microdata be returned?
#' @param M Integer. The first \code{M} implicates are returned for fused variables. Set \code{M = Inf} to return all implicates. Ignored if \code{var} contains only ACS variables (i.e. no implicates)
#' @param df Data frame. Data frame used to identify a subset of rows to return. Default is to return all rows.
#' @param cores Integer. Number of cores used by the \code{\link[fst]{fst-package}} when reading from disk.
#' @param source Character Specifies where to look for \code{var}: all available microdata (\code{source = "all"}); only ACS microdata (\code{source = "ACS"}); or only fused microdata (\code{source = "fused"}). Note that no observation weights are returned if \code{source = "fused"}, since weights are stored in the ACS microdata.
#' @param silent Logical. If \code{FALSE}, a warning is issued if any \code{var} cannot be located in available local files.
#'
#' @details The \code{var} argument can contain a mix of household- and/or person-level variables. When \code{respondent = "household"}, the reference person (i.e. head of household) value is returned for any person-level variables. When \code{respondent = "person"}, the values of any household-level variables are replicated for each person in the household.
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

# year = 2019
# var = c("dollarel")
# respondent = "household"
# M = 2
# df = NULL
# cores = 2
# source = "fused"
# silent= TRUE

#-----

assemble <- function(year,
                     var,
                     respondent = "household",
                     M = 1,
                     df = NULL,
                     cores = 1,
                     source = "all",
                     silent = FALSE) {

  # Check validity of the working directory path
  # Checks if "/fusionData" is part of the path, as this is required
  b <- strsplit(full.path(getwd()), .Platform$file.sep, fixed = TRUE)[[1]]
  i <- which(b == "fusionData")
  if (length(i) == 0) stop("'/fusionData' is not part of the working directory path; this is required.")
  dir <- paste(b[1:i], collapse = .Platform$file.sep)

  # Respondent identifier ("H" or "P")
  rtype <- substring(toupper(respondent), 1, 1)

  # Initial argument check
  stopifnot({
    year >= 2005 & year %% 1 == 0
    is.character(var) & length(var) > 0
    rtype %in% c("H", "P")
    is.null(df) | is.data.frame(df)
    M > 0
    cores > 0 & cores %% 1 == 0
    tolower(source) %in% c('all', 'acs', 'fused')
    is.logical(silent)
  })

  # Universal variables always returned in output, if possible
  uvar <- c('M', 'year', 'hid', 'pid', 'weight')

  fst::threads_fst(nr_of_threads = cores)
  setDTthreads(threads = cores)
  #if (is.null(df)) df <- data.table()
  hh <- rtype == "H"
  v <- setdiff(var, uvar)

  #-----

  result <- lapply(year, function(yr) {

    # Get file paths to both ACS and fused microdata for 'yr'
    pa <- list.files(file.path(dir, "survey-processed/ACS"), pattern = paste0(yr, "_._processed.fst"), recursive = TRUE, full.names = TRUE)
    pa <- sort(pa, decreasing = rtype == "P") # Sort the processed ACS paths, to place either Household or Person microdata first
    pc <- list.files(file.path(dir, "survey-processed/ACS"), pattern = paste0(yr, "_._custom.fst"), recursive = TRUE, full.names = TRUE)
    pf <- rev(list.files(file.path(dir, "fusion"), pattern = paste0(yr, "_._fused.fsd"), recursive = TRUE, full.names = TRUE))
    fpaths <- switch(tolower(source),
                     all = c(pa, pc, pf),
                     acs = c(pa, pc),
                     fused = pf)

    d <- data.table()
    for (x in fpaths) {

      a <- substring(basename(x), 1, 4) == "ACS_"
      r <- rev(strsplit(x, "_")[[1]])[[2]]
      xn <- fst::metadata_fst(x)$columnNames
      keep <- if (!all(c('year', 'hid') %in% xn)) {
        warning("Skipping file ",  basename(x), " due to irregular file structure")
        NULL
      } else {
        # Excludes any variables already in 'd'
        temp <- intersect(xn, c(v, names(df)))
        if (a & rtype == r) temp <- c(temp, 'weight')
        setdiff(temp, names(d))
      }

      # Load requested variables from disk
      if (length(keep)) {
        dt <- fusionModel::read_fsd(path = x,
                                    columns = intersect(xn, c('M', 'year', 'hid', 'pid', keep)),
                                    M = M,
                                    df = if (is.null(df)) NULL else select(df, any_of(xn)),
                                    cores = cores)
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
      # Since it is an left merge, persons not in the household data (i.e. group quarter population) will contain NA values
      # If no household variables need to be added, then the returned person microdata will include the group quarter individuals
      if (nrow(dt) > 0) {
        dt <- setkeyv(dt, cols = intersect(c('M', 'year', 'hid', 'pid'), names(dt)))
        if (nrow(d) == 0) {
          d <- dt
        } else {
          #d <- merge(d, dt, by = intersect(key(d), key(dt)), allow.cartesian = TRUE)
          d <- collapse::join(x = d,
                              y = dt,
                              on = intersect(key(d), key(dt)),
                              how = "left",  # See NOTE above about left join and household vs. person data
                              multiple = TRUE,
                              verbose = FALSE)
        }
      }
      rm(dt)

    }

    return(d)

  }) %>%
    rbindlist() %>%
    setcolorder(neworder = intersect(c(uvar, var), names(.)))

  # Set keys
  if (nrow(result) > 0) setkeyv(result, cols = intersect(c('M', 'year', 'hid', 'pid'), names(result)))

  # Any 'var' missing in the output? If so, report as warning.
  miss <- setdiff(v, names(result))
  if (length(miss) & !silent) warning("Could not locate the following variable(s): ", paste(miss, collapse = ", "))

  return(result)

}
