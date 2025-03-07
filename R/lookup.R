#' Lookup details about available variables
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
# year = 2015
# var <- c("hincp", "pov_ratio", "state", "puma10", "agep", "btuel", "cooltype", "dollarel")
# year = NULL
# var = NULL

lookup <- function(var = NULL, year = NULL) {

  # Loading fusionData dictionary (not necessarily available)
  # data(dictionary)
  # names(dictionary) <- tolower(names(dictionary))

  # Check validity of the working directory path
  # Checks if "/fusionData" is part of the path, as this is required
  b <- strsplit(full.path(getwd()), .Platform$file.sep, fixed = TRUE)[[1]]
  i <- which(b == "fusionData")
  if (length(i) == 0) stop("'/fusionData' is not part of the working directory path; this is required.")
  dir <- paste(b[1:i], collapse = .Platform$file.sep)

  # Get file path(s) to ACS and fused microdata, based on values in 'year'
  fpaths <- lapply(if (is.null(year)) "" else year, function(yr) {
    pa <- list.files(file.path(dir, "survey-processed/ACS"), pattern = paste0(yr, "_._processed.fst"), recursive = TRUE, full.names = TRUE)
    pc <- list.files(file.path(dir, "survey-processed/ACS"), pattern = paste0(yr, "_._custom.fst"), recursive = TRUE, full.names = TRUE)
    pf <- rev(list.files(file.path(dir, "fusion"), pattern = paste0(yr, "_._fused.fsd"), recursive = TRUE, full.names = TRUE))
    c(pa, pc, pf)
  }) %>%
    unlist()

  # Extract the 'var' present in each file in 'fpaths'
  vlist <- if (is.null(var)) {
    lapply(fpaths, function(x) fst::metadata_fst(x)$columnNames)
  } else {
    lapply(fpaths, function(x) intersect(var, fst::metadata_fst(x)$columnNames))
  }

  acs.year <- str_extract(basename(fpaths) , "(\\d{4})(?=_[^_]*_[^_]*$)")
  rtype <- str_extract(basename(fpaths), ".(?=_[^_]*$)")
  rtype <- case_when(rtype == "H" ~ "Household", rtype == "P" ~ "Person")
  source <- basename(fpaths) %>%
    str_extract("^[^_]*_[^_]*") %>%
    str_split(pattern = "_", n = 2, simplify = TRUE) %>%
    as.data.frame() %>%
    setNames(c('survey', 'vintage')) %>%
    mutate(source = 1:n())

  out <- vlist %>%
    tibble::enframe(name = "source", value = "var") %>%
    left_join(source, by = join_by(source)) %>%
    mutate(respondent = rtype,
           acs_year = acs.year,
           path = gsub(dir, "", fpaths, fixed = TRUE),
           source = NULL) %>%
    tidyr::unnest(var) %>%
    filter(!var %in% c('M', 'year', 'hid', 'weight'), !str_detect(var, "^rep_\\d*$")) %>%
    #left_join(dictionary, by = join_by(var == variable, survey, vintage, respondent)) %>%
    #select(var, acs_year, respondent, survey, vintage, description, type, values, path)
    select(var, acs_year, respondent, survey, vintage, path) %>%
    arrange(across(everything()))

  return(out)

}

