#' Read fusion output from disk
#'
#' @description Efficiently read fusion output that was written to disk, optionally returning a subset of rows and/or columns. Since a \code{.fsd} file is simply a \code{\link[fst]{fst}} file under the hood, this function also works on any \code{.fst} file.
#' @param path Character. Path to a \code{.fsd} (or \code{.fst}) file, typically produced by \code{\link{fuse}}.
#' @param columns Character. Column names to read. The default is to return all columns.
#' @param df Data frame. Data frame used to identify a subset of rows to return. Default is to return all rows.
#' @param cores Integer. Number of cores used by \code{\link[fst]{fst}}.
#' @details If \code{df} is provided and the file size on disk is less than 250 MB, then a full read and inner \code{\link[collapse]{join}} is performed. For larger files, a manual read of the required rows is performed, using \code{\link[collapse]{fmatch}} for the matching operation.
#' @return A \code{\link[data.table]{data.table}}; keys are preserved if present in the on-disk data. When \code{path} points to a \code{.fsd} file, it includes an integer column "M" indicating the implicate assignment of each observation (unless explicitly ignored by \code{columns}).
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
#' sim <- fuse(data = recipient, fsn = fsn.path, M = 5, fsd = "results.fsd")
#'
#' # Read the fusion output saved to disk
#' sim <- read_fsd(sim)
#' head(sim)
#'
#' @export

read_fsd <- function(path,
                     columns = NULL,
                     df = NULL,
                     cores = max(1, parallel::detectCores(logical = FALSE) - 1)) {

  stopifnot({
    endsWith(path, ".fsd") | endsWith(path, ".fst")
    file.exists(path)
    is.null(columns) | is.character(columns)
    is.null(df) | is.data.frame(df)
    cores > 0 & cores %% 1 == 0
  })

  require(collapse, quietly = TRUE)
  n <- fst::threads_fst()
  fst::threads_fst(nr_of_threads = cores)

  if (is.null(df)) {
    d <- fst::read_fst(path, columns = columns, as.data.table = TRUE)
  } else {
    meta <- fst::metadata_fst(path)
    stopifnot(all(names(df) %in% meta$columnNames))
    df <- unique(qDT(df))
    # If the file size is less than 250 MB, simply read the full data and subset via collapse::join()
    if (file.size(path) / 1e6 < 250) {
      d <- fst::read_fst(path, columns = columns, as.data.table = TRUE)
      d <- collapse::join(d, df, how = "inner", verbose = FALSE)
    } else {
      if (is.null(columns)) columns <- meta$columnNames else stopifnot(all(columns %in% meta$columnNames))
      d <- fst::fst(path)
      m <- qDT(d[names(df)])
      i <- m %iin% df  # Collapse equivalent of which(x %in% table) using fmatch()
      d <- qDT(d[i, setdiff(columns, names(m))])
      d <- cbind(d, m[i, ])
      setcolorder(d, columns)
      setkeyv(d, cols = meta$keys)
    }
  }

  # Reset number of threads
  fst::threads_fst(nr_of_threads = n)

  return(d)

}
