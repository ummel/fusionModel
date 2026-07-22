#' Impute Missing Data via Microdata Fusion
#'
#' @description
#' Imputes missing values (\code{NA}) in a data frame using sequential microdata fusion.
#' Wraps iterative calls to \code{\link{train}} and \code{\link{fuse}} under the hood
#' to provide a fast, non-parametric imputation workflow powered by LightGBM gradient boosting.
#'
#' @param data Data frame (or \code{data.table}) containing missing values (\code{NA}) to be imputed.
#' @param weight Character string, optional. Name of the observation sampling weight
#'   column in \code{data}. If \code{NULL} (default), uniform weights equal to 1 are assumed.
#'   Weight column must not contain missing values.
#' @param ignore Character vector, optional. Names of columns in \code{data} to ignore.
#'   Ignored variables are neither imputed nor utilized as predictors in the underlying fusion models.
#' @param cores Integer. Number of physical CPU cores used for parallel computation
#'   by \code{fst}, \code{data.table}, and \code{LightGBM}. Defaults to physical cores minus 1
#'   (\code{parallel::detectCores(logical = FALSE) - 1L}).
#'
#' @details
#' \code{impute()} automates missing data imputation across complex microdata datasets
#' through the following steps:
#'
#' \describe{
#'   \item{1. Sequential Ordering}{Variables containing missing values are sorted by total
#'     missingness and imputed sequentially, starting with the variable with the fewest
#'     \code{NA} values. Once a variable is imputed, its complete values immediately become
#'     available as predictors for subsequent target variables.}
#'   \item{2. Predictor Screening}{For larger datasets (>10,000 rows or >20 columns),
#'     a rank-based Spearman correlation matrix is computed to pre-screen predictors.
#'     To maintain efficiency while preserving strong predictor signals, candidate
#'     predictors for each target are filtered to those meeting an absolute correlation threshold
#'     (> 0.025), capped at a maximum of 30 predictors (and a minimum of 5).}
#'   \item{3. Stratified Training Subsampling}{For variables with large non-missing sample sizes,
#'     a stratified sample of non-missing training observations (ranging between 5,000 and 50,000
#'     rows) is selected to speed up model fitting.}
#'   \item{4. Model Fitting & Fusion}{A LightGBM fusion model is trained using 80% validation
#'     splitting (\code{nfolds = 0.8}) to determine optimal tree depth and early stopping.
#'     Missing values for that target variable are then filled via \code{\link{fuse}} with
#'     a single implicate (\code{M = 1}).}
#' }
#'
#' Because \code{LightGBM} natively accommodates \code{NA} values within predictor columns,
#' partially complete predictor variables can be used effectively during model training.
#'
#' \strong{Parallel Processing Note:} Because underlying computations in \code{LightGBM},
#' \code{fst}, and \code{data.table} rely on OpenMP multithreading, it is recommended to manage
#' parallelism via the \code{cores} parameter rather than wrapping \code{impute()} inside
#' outer forked parallel loops.
#'
#' @return A data frame (or \code{data.table}, matching the input class of \code{data})
#'   with all missing values imputed. Original column ordering and variable data types
#'   are preserved.
#'
#' @examples
#' \dontrun{
#' library(fusionModel)

#' # Load sample microdata and introduce random missing values (NA)
#' data(recs)
#' sample_data <- recs[, 2:7]
#' set.seed(123)
#' miss_mask <- replicate(ncol(sample_data), runif(nrow(sample_data)) < runif(1, 0.05, 0.25))
#' sample_data[miss_mask] <- NA
#'
#' # Check missing value counts per column
#' colSums(is.na(sample_data))
#'
#' # Impute missing values
#' imputed_data <- impute(sample_data)
#'
#' # Verify that no NA values remain
#' anyNA(imputed_data)
#' head(imputed_data)
#' }
#'
#' @export

#---

# library(tidyverse)
# library(data.table)
# source("R/utils.R")
# data <- fst::read_fst("~/Documents/Projects/fusionData/test_hus.fst")
# weight = "WGTP"
# cores = 1
#
# #ignore <- names(select(data[1, ], SERIALNO, PUMA, WGTP1:SERIALNO_original))
# # Specify this in ASC processing script
# ignore <- setdiff(names(data), c(y, "WGTP", "puma_rent", "puma_value", "puma_income", "puma_mortgage", "ST", "NP", "ACR", "BLD", "FS", "HFL", "HHL", "BDSP", "BDS", "RMSP", "RMS", "TEN", "VEH", "YBL", "YRBLT", "HINCP", "FES", "WIF", "R18", "R65"))
#
# test <- impute(data, "WGTP", ignore = ignore)

#---

impute <- function(data,
                   weight = NULL,
                   ignore = NULL,
                   cores = parallel::detectCores(logical = FALSE) - 1L) {

  t0 <- Sys.time()

  fst::threads_fst(nr_of_threads = cores)
  setDTthreads(threads = cores)

  stopifnot(exprs = {
    is.data.frame(data)
    ncol(data) > 1
    anyNA(data)
    any(is.null(weight), length(weight) == 1 & weight %in% names(data))
    any(is.null(ignore), all(ignore %in% names(data)))
    cores > 0 & cores %% 1 == 0 & cores <= parallel::detectCores(logical = FALSE)
  })

  dnames <- names(data)
  data.dt <- is.data.table(data)
  d <- as.data.table(data)
  rm(data)

  if (is.null(weight)) {
    weight = "W_.._"
    d[, W_.._ := 1L]
  } else {
    if (anyNA(d[[weight]])) stop("NA values are not allowed in 'weight'")
  }

  miss <- sort(colSums(is.na(d)))
  y <- names(miss)[miss > 0 & miss < nrow(d)]
  y <- setdiff(y, ignore)
  if (!length(y)) stop("No un-ignored columns with NA values to impute")
  temp.fsn <- paste0(tempfile(), ".fsn")

  #---

  # Predictor prescreen step if 'd' is sufficiently large
  #if (length(x) & (nrow(d) > 10e3 | ncol(d) > 10)) {
  if ((nrow(d) > 10e3 | ncol(d) > 20)) {

    d2 <- copy(d)
    d2[, (weight) := NULL]

    # Convert 'd2' to plausible ranks for correlation screening
    # All output columns should be NA, integer, or logical
    # NA's in input are preserved in output
    for (i in 1:ncol(d2)) {
      z <- d2[[i]]
      if (is.numeric(z)) {
        # Ties method 'dense' ensures integer output with minimum of 1 and maximum of length(na.omit(z))
        z <- frank(z, ties.method = "dense", na.last = "keep")
      } else {
        # Converts ordered factors to integer
        if (is.ordered(z)) {
          z <- as.integer(z)
        } else {
          if (!is.logical(z)) {
            # Converts character and un-ordered factors to TRUE for the most-common (non-NA) value and FALSE otherwise
            # zt <- collapse::qtable(z, na.exclude = TRUE)
            # z <- z == names(which.max(zt))
            z <- z == collapse::fmode(z, na.rm = TRUE)
          }
        }
      }
      # Replace any NA's with the median value
      # This makes the cor() operation below much faster
      if (anyNA(z)) z[is.na(z)] <- median(z, na.rm = TRUE)
      # Update values in 'd2'
      set(d2, j = i, value = z)
    }

    # Randomly down-sample when 'd2' is large
    if (nrow(d2) > 100e3) d2 <- d2[sample.int(nrow(d2), 100e3), ]

    # Correlation matrix
    # Note that Spearman (rank) correlations are used (data pre-ranked) to reduce effect of outliers
    ok <- setdiff(names(d2), ignore)
    cmat <- suppressWarnings(cor(d2[, ..y], d2[, ..ok]))
    cmat[is.na(cmat)] <- 0

    # Initial correlation screening, based on absolute correlation value
    xlist <- lapply(y, function(v) {

      p <- cmat[v, ]
      names(p) <- colnames(cmat)
      p <- p[names(p) != v]
      p <- sort(abs(p), decreasing = TRUE)

      # Restrict to predictors that meet arbitrary absolute correlation threshold (> 0.025)
      # Limit to maximum 30 potential predictors for each model (hard cap)
      # Attempt to retain minimum 5 potential predictors regardless of correlation
      xv <- names(p[p > 0.025])
      xv <- xv[1:min(30, length(xv))]
      if (length(xv) < 5) xv <- names(p)[1:min(5, length(p))]

      return(xv)

    })

    rm(d2)

  } else {

    xlist <- lapply(y, function(v) setdiff(names(d), c(v, ignore, weight)))

  }

  xlist <- setNames(xlist, y)

  #---

  # Coerce any character variables in 'y' or 'xlist' to unordered factor
  # This is necessary to avoid errors in train() and fuse(), which expect factors
  # The 'cconv' columns are converted to character prior to returning final function output
  ccols <- names(which(sapply(d, is.character)))
  cconv <- intersect(ccols, unique(c(y, unlist(xlist))))
  if (length(cconv) > 0)  d[, (cconv) := lapply(.SD, factor), .SDcols = cconv]

  # Check number of factor levels in 'y' variables
  nlev <- sapply(d[, ..y], nlevels)
  bad <- names(nlev)[nlev >= 200]
  if (length(bad)) {
    cli::cli_alert_warning("Detected categorical imputation variable(s) with more than 200 levels (can be slow):\n{paste(bad, collapse = '\n')}")
  }

  #---

  pb <- txtProgressBar(max = length(y), style = 3)

  for (i in 1:length(y)) {

    # Response and predictor variables
    v <- y[i]
    xv <- xlist[[i]]
    vtrain <- c(weight, v, xv)

    # Observations to impute
    imp <- is.na(d[[v]])

    # Training dataset
    dtrain <- d[!imp, ..vtrain]

    # If sample size is large, use a stratified sample of the training data
    # If there are few missing values, then fewer training observations are possible (10x the number missing)
    # Restricts training sample to no less than 5k and no more than 50k observations
    maxn <- min(max(5e3, 10 * sum(imp)), 50e3)
    if (nrow(dtrain) > maxn) {
      keep <- stratify(y = dtrain[[v]],
                       ycont = is.numeric(dtrain[[v]]),
                       tfrac = maxn / nrow(dtrain),
                       ntiles = 20)
      dtrain <- dtrain[keep, ]
    }

    # Train fusion model (saved to 'temp.fsn' temporary file)
    # Suppress output to the console
    invisible(capture.output(
      suppressMessages(
        fusionModel::train(data = dtrain,
                           y = v,
                           x = xv,
                           fsn = temp.fsn,
                           weight = weight,
                           nfolds = 0.8,
                           nquantiles = 2,
                           cores = cores)
      )
    ))

    # Perform fusion/imputation
    # Suppress output to the console
    invisible(capture.output(
      suppressMessages(
        p <- fusionModel::fuse(data = d[imp, ],
                               fsn = temp.fsn,
                               M = 1,
                               cores = cores)
      )
    ))

    # Update 'd' with imputed values
    set(d, i = which(imp), j = v, value = p[[2]])

    # Update progress bar
    setTxtProgressBar(pb, value = i)

  }

  # Close progress bar
  close(pb)

  # Remove temporary fusion model
  unlink(temp.fsn)

  # Check for NA's in output
  stopifnot(!anyNA(d[, ..y]))

  # If any character columns were converted to factor, convert back to character
  if (length(cconv) > 0)  d[, (cconv) := lapply(.SD, as.character), .SDcols = cconv]

  # Ensure output column order and class matches input 'data'
  suppressWarnings(set(d, j = "W_.._", value = NULL))  # Removes the placeholder weight variable, if present
  setcolorder(d, dnames)
  if (!data.dt) d <- as.data.frame(d)

  # Report processing time
  tout <- difftime(Sys.time(), t0)
  cli::cli_alert_success("Total processing time: {signif(as.numeric(tout), 3)} {attr(tout, 'units')}")

  return(d)

}
