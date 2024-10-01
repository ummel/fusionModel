#' Read ORNL UrbanPop data from disk
#'
#' @description NOTE: For fusionACS internal use only! Efficiently read pre-processed ORNL UrbanPop data from disk. Function arguments only make sense if you are familiar with the structure of the processed UrbanPop data.
#' @param path Character. File path to \code{\link[fst]{fst}} file containing UrbanPop data.
#' @param year Integer. Year(s) to select.
#' @param state Integer. State FIPS code(s) to select.
#' @param county Integer. County FIPS code(s) to select.
#' @param tract_bg Integer. Tract and block group FIPS code(s) to select.
#' @param hid Integer. ACS-PUMS household ID(s) to select.
#' @param df Data frame containing at least 'year' and/or 'state' columns. Provides unique combinations of the above argument values to return. \code{df} is used to perform an inner merge on an initial subset of data based on 'year' or 'state'.
#' @param cores Integer. Number of cores used by \code{\link[fst]{read_fst}}.
#' @details Provides an efficient and fast way to load a subset of UrbanPop data into memory. An initial subset is read using \code{state} or \code{year} to restrict rows. Then \code{\link[data.table]{data.table}} operations (subset or merge) are used to efficiently reduce to the final returned subset.
#' @return A keyed \code{\link[data.table]{data.table}}.
#' @examples
#'up.path <- "~/Documents/Projects/fusionData/urbanpop/Processed national UrbanPop.fst"
#'
#'out <- read_up(path = up.path, year = 2015, state = 4)
#'unique(dplyr::select(out, year, state))
#'
#'out <- read_up(path = up.path, year = 2015:2016, state = c(2, 12, 15))
#'unique(dplyr::select(out, year, state))
#'
#'up.df <- data.frame(year = c(2015, 2018), state = c(8, 5), county = c(1, 3))
#'out <- read_up(path = up.path, df = up.df)
#'unique(dplyr::select(out, year, state, county))
#' @export

read_up <- function(path,
                  year = NULL,
                  state = NULL,
                  county = NULL,
                  tract_bg = NULL,
                  hid = NULL,
                  df = NULL,
                  cores = max(1, parallel::detectCores(logical = FALSE) - 1)) {

  stopifnot({
    endsWith(path, ".fst")
    file.exists(path)
  })

  n <- fst::threads_fst()
  fst::threads_fst(nr_of_threads = cores)
  d <- fst::fst(path)
  N <- nrow(d)

  if (!is.null(df)) {
    stopifnot(is.data.frame(df))
    stopifnot(nrow(df) > 0)
    if (!any(c("state", "year") %in% names(df))) stop("You must provide columns for 'year' and/or 'state' in 'df'")
    year <- unique(df$year)
    state <- unique(df$state)
    df <- unique(df)
    cat("Valid 'df' provided for merge; ignorning any other arguments supplied\n")
  }

  # The initial (expensive) load of data uses either 'state' or 'year', whichever appears likely to result in smaller initial load to RAM
  if (!is.null(state) & uniqueN(state) / 50 <= uniqueN(year) / 4) {
    d <- d[d[["state"]] %in% as.integer(state), ]
    state <- NULL
  } else {
    if (is.null(year)) {
      stop("You must provide argument 'year' and/or 'state'")
    } else {
      d <- d[d[["year"]] %in% as.integer(year), ]
      year <- NULL
    }
  }

  d <- setDT(d)
  if (is.null(df)) {
    # https://stackoverflow.com/questions/11612235/select-rows-from-a-data-frame-based-on-values-in-a-vector
    if (!is.null(year)) {X <- year; d <- d[J(X), on = .(year), nomatch = 0]}
    if (!is.null(state)) {X <- state; d <- d[J(X), on = .(state), nomatch = 0]}
    if (!is.null(county)) {X <- county; d <- d[J(X), on = .(county), nomatch = 0]}
    if (!is.null(tract_bg)) {X <- tract_bg; d <- d[J(X), on = .(tract_bg), nomatch = 0]}
    if (!is.null(hid)) {X <- hid; d <- d[J(X), on = .(hid), nomatch = 0]}
  } else {
    df <- df[, intersect(names(df), names(d))]
    cat("Using the following 'df' variables for merge:", paste(names(df), collapse = ", "), "\n")
    d <- merge.data.table(d, df, by = names(df), sort = FALSE)
  }

  setkey(d, state, county, tract_bg, year, hid)
  dpct <- signif(100 * (nrow(d) / N), 3)
  cat("Returning ", ifelse(dpct > 1, dpct, "< 1"), "% of the total UrbanPop observations\n", sep = "")
  fst::threads_fst(nr_of_threads = n)  # Reset number of threads
  return(d)

}
