#' Read fusion output from disk
#'
#' @description
#' Read fusion output that was written directly to disk via \code{\link{fuse}}.
#'
#' @param fsd Character. File path ending in \code{.fsd} as produced by call to \code{\link{fuse}}.
#' @param cores Integer. Number of cores used by \code{\link[data.table]{fread}}.
#'
#' @details The file written by \code{\link{fuse}} includes specially formatted metadata in the column names that \code{read_fsd} uses to construct the appropriate output values and classes. Reading \code{fsd} with a different file reader function will not give correct results.
#'
#' @return A \code{\link[data.table]{data.table}} with integer column "M" indicating the implicate assignment of each observation. Note that the ordering of recipient observations is consistent within implicates, so do not change the row order if using with \code{\link{analyze}} or \code{\link{validate}}.
#'
#' @examples
#' # Build a fusion model using RECS microdata
#' # Note that "fusion_model.fsn" will be written to working directory
#' ?recs
#' fusion.vars <- c("electricity", "natural_gas", "aircon")
#' predictor.vars <- names(recs)[2:12]
#' fsn.path <- train(data = recs, y = fusion.vars, x = predictor.vars)
#'
#' # Write fusion output directly to disk
#' # Note that "results.fsd" will be written to working directory
#' recipient <- recs[predictor.vars]
#' sim <- fuse(data = recipient, fsn = fsn.path, M = 5, csv = "results.fsd")
#'
#' # Read the fusion output saved to disk
#' sim <- read_fsd(sim)
#' head(sim)
#'
#' @export

read_fsd <- function(fsd, cores = data.table::getDTthreads()) {

  # Input 'd' can be a data frame for calling within fuse()
  # Otherwise, load the .fsd file from disk using fread()
  if (is.data.frame(fsd)) {
    d <- fsd
  } else {
    # Have to rename .fsd to .csv.gz (temporarily) for fread() to recognize the file
    stopifnot(endsWith(fsd, ".fsd"))
    csv <- sub("\\.fsd$", ".csv.gz", fsd)
    file.rename(from = fsd, to = csv)
    d <- data.table::fread(file = csv, data.table = TRUE, nThread = cores)
    file.rename(from = csv, to = fsd)
  }

  # Parse the column name metadata
  dclass <- sapply(d, class)
  temp <- lapply(names(d), strsplit, split = "|!|", fixed = TRUE)
  temp <- unlist(temp, recursive = FALSE)
  ynames <- sapply(temp, function(x) x[1])
  yclass <- sapply(temp, function(x) x[2])
  yordered <- sapply(temp, function(x) as.logical(as.integer(x[3])))
  ylevels <- lapply(temp, function(x) {
    if (length(x) < 4) NA else unlist(strsplit(x[4], "|%|", fixed = TRUE))
  })

  # Ensure simulated variables are correct data type with appropriate labels/levels
  for (i in seq_along(ynames)) {
    if (yclass[i] == "factor") {
      lev <- ylevels[[i]]
      set(d, i = NULL, j = i, value = factor(lev[d[[i]] + 1], levels = lev, ordered = yordered[i]))
    }
    if (yclass[i] == "logical" & dclass[i] != "logical") set(d, i = NULL, j = i, value = as.logical(d[[i]]))
    if (yclass[i] == "integer" & dclass[i] != "integer") set(d, i = NULL, j = i, value = as.integer(d[[i]]))
  }

  setnames(d, ynames)
  return(d)

}
