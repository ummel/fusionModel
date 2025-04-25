#' Impute missing data via fusion
#'
#' @description
#' A universal missing data imputation tool that wraps successive calls to \code{\link{train}} and \code{\link{fuse}} under the hood. Designed for simplicity and ease of use.
#' @param data A data frame with missing values.
#' @param weight Optional name of observation weights column in \code{data}.
#' @param ignore Optional names of columns in \code{data} to ignore. These variables are neither imputed nor used as predictors.
#' @param cores Number of physical CPU cores used by \code{\link[lightgbm]{lightgbm}}. LightGBM is parallel-enabled on all platforms if OpenMP is available.
#' @details Variables with missing values are imputed sequentially, beginning with the variable with the fewest missing values. Since LightGBM models accommodate NA values in the predictor set, all available variables are used as potential predictors (excluding \code{ignore} variables). For each call to \code{\link{train}}, 80% of observations are randomly selected for training and the remaining 20% are used as a validation set to determine an appropriate number of tree learners. All LightGBM model parameters are kept at the sensible default values in \code{\link{train}}. Since \code{\link[lightgbm]{lightgbm}} uses OpenMP multithreading, it is not advisable to use \code{\link{impute}} inside a forked/parallel process when \code{cores > 1}.
#' @return A data frame with all missing values imputed.
#' @examples
#' # Create data frame with random NA values
#' ?recs
#' data <- recs[, 2:7]
#' miss <- replicate(ncol(data), runif(nrow(data)) < runif(1, 0.01, 0.3))
#' data[miss] <- NA
#' colSums(is.na(data))
#'
#' # Impute the missing values
#' result <- impute(data)
#' anyNA(result)
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
            zt <- table2(z, na.rm = TRUE)
            z <- z == names(which.max(zt))
          }
        }
      }
      set(d2, j = i, value = z)
    }

    # Randomly down-sample when 'd2' is large
    if (nrow(d2) > 100e3) d2 <- d2[sample.int(nrow(d2), 100e3), ]

    # Correlation matrix
    # Note that Spearman (rank) correlations are used (data pre-ranked) to reduce effect of outliers
    ok <- setdiff(names(d2), ignore)
    cmat <- suppressWarnings(cor(d2[, ..y], d2[, ..ok], use = "pairwise.complete.obs"))
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

  #---

  pb <- txtProgressBar(max = length(y), style = 3)

  for (i in 1:length(y)) {

    # Response and predictor variables
    v <- y[i]
    #xv <- c(xlist[[i]], setdiff(y, v))
    #xv <- setdiff(xv, ignore)
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
      fusionModel::train(data = dtrain,
                         y = v,
                         x = xv,
                         fsn = temp.fsn,
                         weight = weight,
                         nfolds = 0.8,
                         cores = cores)
    ))

    # Perform fusion/imputation
    # Suppress output to the console
    invisible(capture.output(
      p <- fusionModel::fuse(data = d[imp, ],
                             fsn = temp.fsn,
                             M = 1,
                             cores = cores)
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
  cat("Total processing time:", signif(as.numeric(tout), 3), attr(tout, "units"), "\n", sep = " ")

  return(d)

}
