#' Read Fusion Output Files From Disk
#'
#' @description
#' `read_fsd()` efficiently loads microdata fusion outputs (`.fsd` files) or standard Fast
#' Storage (`.fst`) datasets into memory. Because fusion files can become extremely large
#' when generating multiple implicates, `read_fsd()` leverages partial reading capabilities
#' to load only the specific rows, columns, or implicates needed, drastically reducing memory
#' footprint and processing time.
#'
#' @param path Character. Path to a valid `.fsd` or `.fst` file on disk, typically created
#'   by calling \code{\link{fuse}} with the `fsd` argument specified.
#' @param columns Character vector. Specific column names to read from disk. The default
#'   (`NULL`) loads all available columns in the dataset.
#' @param M Numeric or Integer. The maximum number of implicates to read (e.g., `M = 5` loads
#'   implicates `1` through `5`). Set `M = Inf` to load all implicates available in the file.
#'   Default is `1`. Ignored if the target dataset does not contain an `M` column.
#' @param df Data frame or data.table. Optional filtering subset. If provided, only rows in
#'   the file matching the combination of key values in `df` will be returned. Default is `NULL`.
#' @param cores Integer. Number of CPU threads to assign to `fst` for multi-threaded decompression.
#'   Defaults to available logical cores minus one (minimum 1).
#'
#' @details
#' The `.fsd` format (Fusion Data file) is functionally identical to the high-performance binary
#' format supplied by the \pkg{fst} package, optimized for fast random access and compression.
#'
#' **Subsetting Mechanics:**
#' \itemize{
#'   \item **Implicate Filtering:** When the dataset contains multiple implicates (identified by column `M`),
#'     `read_fsd()` calculates contiguous row indices corresponding to `M <= max_M` before reading,
#'     preventing unnecessary disk I/O for unused implicates.
#'   \item **Row Subsetting (`df`):** If a matching dataset `df` is provided:
#'     \itemize{
#'       \item For files **< 100 MB**, the relevant implicate range is fully loaded into memory and filtered
#'         using an efficient inner join via \code{\link[collapse]{join}}.
#'       \item For files **>= 100 MB**, low-memory index matching is performed out-of-core using high-speed
#'         C-level matching via \code{\link[collapse]{fmatch}} (`%iin%`), loading only matching record indices.
#'     }
#' }
#'
#' @return A \code{\link[data.table]{data.table}} containing the requested columns and rows.
#'   Original `data.table` key attributes are automatically preserved if present in the source file.
#'
#' @examples
#' \dontrun{
#' # Build a fusion model using RECS microdata
#' fusion.vars <- c("electricity", "natural_gas", "aircon")
#' predictor.vars <- names(recs)[2:12]
#' fsn.path <- train(data = recs, y = fusion.vars, x = predictor.vars)
#'
#' # Write fusion implicates directly to disk during fusion
#' recipient <- recs[predictor.vars]
#' fsd_file <- file.path(tempdir(), "results.fsd")
#' fuse(data = recipient, fsn = fsn.path, M = 3, fsd = fsd_file)
#'
#' # Read only the first implicate (M = 1) for all variables
#' sim_m1 <- read_fsd(path = fsd_file, M = 1)
#' head(sim_m1)
#'
#' # Read specific columns across all 3 implicates
#' sim_sub <- read_fsd(
#'   path = fsd_file,
#'   columns = c("M", "electricity"),
#'   M = 3
#' )
#' head(sim_sub)
#' }
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

  # Temporarily configure parallel decompression threads in 'fst'
  n <- fst::threads_fst()
  fst::threads_fst(nr_of_threads = cores)

  meta <- fst::metadata_fst(path)
  if (is.null(columns)) {
    columns <- meta$columnNames
  } else {
    columns <- unique(columns)
    stopifnot(all(columns %in% meta$columnNames))
  }

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

      # For larger files, perform low-memory out-of-core index matching using collapse::fmatch via %iin%
      d <- fst::fst(path)
      m <- qDT(d[i, names(df)])
      i <- i[m %iin% df] # Uses 'collapse' package equivalent of which(x %in% table) using fmatch()
      d <- qDT(d[i, setdiff(columns, names(m)), drop = FALSE])
      j <- intersect(names(m), columns)
      if (length(j)) d <- cbind(d, m[i, ..j])

    }
  }

  # Set column order and restore data.table keys
  setcolorder(d, neworder = columns)
  suppressWarnings(setkeyv(d, cols = intersect(meta$keys, columns)))

  # Reset number of 'fst' decompression threads to initial state
  fst::threads_fst(nr_of_threads = n)

  return(d)

}
