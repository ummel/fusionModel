#' Read fusion output from disk
#'
#' @description Efficiently read fusion output that was written to disk, optionally returning a subset of rows and/or columns. Since a \code{.fsd} file is simply a \code{\link[fst]{fst}} file under the hood, this function also works on any \code{.fst} file.
#' @param path Character. Path to a \code{.fsd} (or \code{.fst}) file, typically produced by \code{\link{fuse}}.
#' @param columns Character. Column names to read. The default is to return all columns.
#' @param M Integer. The first \code{M} implicates are returned. Set \code{M = Inf} to return all implicates. Ignored if \code{M} column not present in data.
#' @param df Data frame. Data frame used to identify a subset of rows to return. Default is to return all rows.
#' @param cores Integer. Number of cores used by \code{\link[fst]{fst}}.
#' @details If \code{df} is provided and the file size on disk is less than 100 MB, then a full read and inner \code{\link[collapse]{join}} is performed. For larger files, a manual read of the required rows is performed, using \code{\link[collapse]{fmatch}} for the matching operation.
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

#-----------

# TEST
# library(tidyverse)
# library(collapse)
# library(data.table)
# setwd("/home/kevin/Documents/Projects/fusionData")
#
# test <- read_fsd("fusion/RECS/2020/2019/RECS_2020_2019_H_fused.fsd",
#                  columns = c("M", "dollarel"),
#                  M = 5,
#                  df = df,
#                  cores = 3)

#-----------

read_fsd <- function(path,
                     columns = NULL,
                     M = 1,
                     df = NULL,
                     cores = max(1, parallel::detectCores(logical = FALSE) - 1)) {

  stopifnot({
    endsWith(path, ".fsd") | endsWith(path, ".fst")
    file.exists(path)
    is.null(columns) | is.character(columns)
    is.null(df) | is.data.frame(df)
    M >= 1
    cores > 0 & cores %% 1 == 0
  })

  #require(collapse, quietly = TRUE)
  n <- fst::threads_fst()
  fst::threads_fst(nr_of_threads = cores)

  meta <- fst::metadata_fst(path)
  if (is.null(columns)) {
    columns <- meta$columnNames
  } else {
    columns <- unique(columns)
    stopifnot(all(columns %in% meta$columnNames))
  }

  # TURNED OFF FOR TESTING
  # Add 'M' column to 'df' to subset on the number of implicates
  # if (is.finite(M) & "M" %in% v) {
  #   df <- if (nrow(df) == 0) {
  #     data.table(M = 1:M)
  #   } else {
  #     df %>%
  #       select(any_of(v)) %>%
  #       select(-any_of("M")) %>%
  #       unique() %>%
  #       slice(rep(1:nrow(.), M)) %>%
  #       mutate(M = rep(1:M, each = nrow(.) / M))
  #   }
  # }

  #-----

  # Determine which rows have the requested implicates (M)
  # Since the data are assumed to be sorted by M, this should yield consecutive integers

  if ('M' %in% meta$columnNames) {
    d <- fst::fst(path)
    i <- which(d$M <= M)
    stopifnot(!is.unsorted(i))
  } else {
    i <- 1L:meta$nrOfRows  # Return all rows in 'd', if no implicates column 'M'
  }

  #-----

  if (is.null(df)) {

    d <- fst::read_fst(path, columns = columns, from = i[1], to = i[length(i)], as.data.table = TRUE)

  } else {

    # Check for validity of 'df' column names
    stopifnot(all(names(df) %in% meta$columnNames))
    df <- collapse::funique(df)

    # If the file size is less than 100 MB, simply read the full data and subset via collapse::join()
    if (file.size(path) / 1e6 < 100) {

      d <- fst::read_fst(path, columns = unique(c(columns, names(df))), from = i[1], to = i[length(i)], as.data.table = TRUE)
      d <- collapse::join(d, df, how = "inner", verbose = FALSE)
      d <- d[, ..columns]

    } else {

      d <- fst::fst(path)
      m <- qDT(d[i, names(df)])
      i <- i[m %iin% df] # Uses 'collapse' package equivalent of which(x %in% table) using fmatch()
      d <- qDT(d[i, setdiff(columns, names(m)), drop = FALSE])
      j <- intersect(names(m), columns)
      if (length(j)) d <- cbind(d, m[i, ..j])

    }
  }

  # Set column order and data.table keys
  setcolorder(d, neworder = columns)
  suppressWarnings(setkeyv(d, cols = intersect(meta$keys, columns)))

  # Reset number of threads
  fst::threads_fst(nr_of_threads = n)

  return(d)

}
