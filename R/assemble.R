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

  # Capture the function call; added as attribute in the final output
  mcall <- match.call()

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
  hh <- rtype == "H"

  # Specified source for each 'var', if provided
  vsrc <- ifelse(str_detect(var, ":"), str_extract(var, "^[^:]+"), "")
  v <- ifelse(str_detect(var, ":"), str_extract(var, "(?<=:).*"), var)
  i <- !v %in% uvar
  vsrc <- vsrc[i]
  v <- v[i]

  #-----

  result <- lapply(year, function(yr) {

    # Get file paths to ACS microdata for 'yr'
    pa <- list.files(file.path(dir, "survey-processed/ACS"), pattern = paste0(yr, "_._processed.fst"), recursive = TRUE, full.names = TRUE)
    pa <- sort(pa, decreasing = rtype == "P") # Sort the processed ACS paths, to place either Household or Person microdata first
    pc <- list.files(file.path(dir, "survey-processed/ACS"), pattern = paste0(yr, "_._custom.fst"), recursive = TRUE, full.names = TRUE)

    # Get file path(s) to fused microdata for 'yr'
    pf <- rev(list.files(file.path(dir, "fusion"), pattern = paste0(yr, "_._fused.fsd"), recursive = TRUE, full.names = TRUE))

    # Select paths based on 'source' argument
    fpaths <- switch(tolower(source),
                     all = c(pa, pc, pf),
                     acs = c(pa, pc),
                     fused = pf)

    # The survey name (e.g. RECS_2020) associated with each path
    survey <- str_extract(basename(fpaths), "^[^_]*_[^_]*")

    #---

    # Check the donor survey 'fpaths' for duplicate occurrences of requested variables
    # Stop with error if conflict(s) detected
    # This can occur if there are identically named variables across surveys OR multiple vintages of a donor survey are fused to the same ACS vintage (e.g. RECS 2015 and RECS 2015 fused to ACS 2015-2019)
    i <- substring(survey, 1, 4) != "ACS_"
    vlist <- lapply(fpaths[i], function(x) intersect(v[vsrc == ""], fst::metadata_fst(x)$columnNames))
    u <- table(unlist(vlist))
    check <- lapply(vlist, function(x) intersect(x, names(u)[u > 1]))
    names(vlist) <- names(check) <- survey[i]
    if (any(lengths(check) > 0)) {
      check <- check[lengths(check) > 0]
      error.msg <- paste(capture.output(str(check)), collapse = "\n")
      stop("Conflicting 'var' names. The following variables are present in more than one source file:\n", error.msg, "\nRevise 'var' to use a colon to specify the source; e.g. ", paste0("'", names(check)[1], ":", check[[1]][1], "'"))
    }

    # Update the 'vsrc' object to assign the source survey for 'v' that are unassigned
    temp <- lapply(vlist, function(x) intersect(x, names(u)[u == 1]))
    for (i in seq_along(temp)) {
      for (j in temp[[i]]) {
        vsrc[v == j] <- names(vlist)[i]
      }
    }

    # Data frame with variable sources
    dv <- data.frame(year = yr, var = v, source = vsrc)

    #---

    d <- data.table()
    for (x in fpaths) {

      a <- substring(basename(x), 1, 4) == "ACS_"
      r <- rev(strsplit(x, "_")[[1]])[[2]]
      xn <- fst::metadata_fst(x)$columnNames

      keep <- if (!all(c('year', 'hid') %in% xn)) {
        warning("Skipping file ",  basename(x), " due to irregular file structure")
        NULL
      } else {
        i <- vsrc == str_extract(basename(x), "^[^_]*_[^_]*") | vsrc == ""
        temp <- intersect(xn, c(v[i], names(df)))
        if (a & rtype == r) temp <- c(temp, 'weight')
        setdiff(temp, names(d))  # Excludes any variables already in 'd'
      }

      # Load requested variables from disk
      if (length(keep)) {
        dt <- fusionModel::read_fsd(path = x,
                                    columns = intersect(xn, c('M', 'year', 'hid', 'pid', keep)),
                                    M = M,
                                    df = if (any(xn %in% names(df))) select(df, any_of(xn)) else NULL,
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
        setkeyv(dt, cols = intersect(c('M', 'year', 'hid', 'pid'), names(dt)))
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
          setkeyv(d, cols = key(dt))
        }
      }
      rm(dt)

    }

    #return(d)
    return(list(d, dv))

  })

  # Extract and rbind the attribute data.frames
  attr.df <- result %>%
    purrr::map(2) %>%
    rbindlist()

  # Extract and rbind the microdata output
  result <- result %>%
    purrr::map(1) %>%
    rbindlist() %>%
    setcolorder(neworder = intersect(c(uvar, var), names(.))) %>%
    setattr(name = "origin", value = attr.df) %>%
    setattr(name = "assemble", value = mcall)

  # Set keys
  if (nrow(result) > 0) setkeyv(result, cols = intersect(c('M', 'year', 'hid', 'pid'), names(result)))

  # Any 'var' missing in the output? If so, report as warning.
  miss <- setdiff(v, names(result))
  if (length(miss) & !silent) warning("Could not locate the following variable(s): ", paste(miss, collapse = ", "))

  return(result)

}
