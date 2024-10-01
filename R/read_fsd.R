#' Read fusion output from disk
#'
#' @description Read fusion output that was written directly to disk via \code{\link{fuse}}.
#' @param fsd Character. File path ending in \code{.fsd} as produced by call to \code{\link{fuse}}.
#' @param columns Character. Column names to read. The default is to read all columns.
#' @param cores Integer. Number of cores used by \code{\link[fst]{read_fst}}.
#' @details As of version 2.3.0, this is simply a convenient wrapper around \code{\link[fst]{read_fst}}, since fusion output data files (.fsd) are actually native \code{\link[fst]{fst}} files under the hood.
#' @return A \code{\link[data.table]{data.table}} with integer column "M" indicating the implicate assignment of each observation. Note that the ordering of recipient observations is consistent within implicates, so do not change the row order if using with \code{\link{analyze}} or \code{\link{validate}}.
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

read_fsd <- function(fsd,
                     columns = NULL,
                     cores = max(1, parallel::detectCores(logical = FALSE) - 1)) {

  stopifnot({
    endsWith(fsd, ".fsd")
    file.exists(fsd)
  })

  n <- fst::threads_fst()
  fst::threads_fst(nr_of_threads = cores)
  d <- fst::read_fst(path = fsd,
                     columns = columns,
                     as.data.table = TRUE)
  fst::threads_fst(nr_of_threads = n)  # Reset number of threads

  return(d)

}
