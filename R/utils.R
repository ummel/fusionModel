# Function to "clean" a numeric vector by reducing to significant digits and converting to integer, if possible
cleanNumeric <- function(x, convert = TRUE, ...) {
  x <- signifDigits(x, ...)
  if (convert) x <- convertInteger(x)
  return(x)
}

#------------------

# Function to return numeric vector rounded to reasonable significant digits
# Returns a significant digit-ized result that is within 'tol' (percent) of the original value for all observations
# If minimize = TRUE, function will try converting x to Z-scores first and 'tol' assesed relative to the Z-scores, then return result that minimizes number of unique values
signifDigits <- function(x, tol = 0.001, minimize = FALSE) {

  intFUN <- function(x, orig = x) {
    out <- rep(NA, length(x))
    out[x == 0 | is.na(x) | is.infinite(x)] <- 0
    i <- 1
    while (any(is.na(out))) {
      ind <- which(is.na(out))
      y <- x[ind]
      z <- abs(signif(y, i) - y) / abs(y)
      ok <- ind[z <= tol]
      out[ok] <- i
      i <- i + 1
    }
    return(signif(orig, out))
  }

  x1 <- intFUN(x)

  if (!minimize) {
    return(x1)
  } else {
    x2 <- intFUN(scale(x), x)
    if (length(unique(x1)) <= length(unique(x2))) return(x1) else return(x2)
  }

}

#------------------

# Function to convert a numeric vector to integer, if possible
convertInteger <- function(x) {
  if (all(x[!is.na(x)] %% 1 == 0) & max(x, na.rm = TRUE) < 2*10^9) {
    return(as.integer(round(x)))
  } else {
    return(x)
  }
}

#------------------

# Function returns TRUE if 'x' has only one non-NA value
novary <- function(x) length(unique(na.omit(x))) == 1

#------------------

# Function to detect and impute any missing values in 'data'
# Performs median imputation of continuous variables and frequency-weighted sampling of categorical variables
imputationValue <- function(x, na.ind) {
  if (is.numeric(x)) {
    m <- median(x, na.rm = TRUE)
    m <- ifelse(is.integer(x), as.integer(round(m)), m)
  } else {
    tab <- table(x) / sum(!na.ind)
    m <- sample(names(tab), size = sum(na.ind), replace = TRUE, prob = tab)
  }
  return(m)
}

#------------------

# Function to treat integer and numeric as equal when checking for identical classes in fuse()
sameClass <- function(x, y) {
  if (x[1] == "integer") x <- "numeric"
  if (y[1] == "integer") y <- "numeric"
  identical(x, y)
}

#------------------

# Weighted mean; slightly faster than weighted.mean()
wmean <- function(x, w) {
  w <- w / sum(w)
  sum(w * x) / sum(w)
}

#------------------

# Weighted standard deviation
# Equivalent to Hmisc::wtd.var() with normwt = TRUE and taking sqrt() of result
wsd <- function(x, w) {
  w <- (w / sum(w)) * length(x)
  sw <- sum(w)
  xbar <- sum(w * x) / sw
  sqrt(sum(w * ((x - xbar) ^ 2)) / (sw - 1))
}

#------------------

# Detect if a numeric variable is likely to be zero-inflated
# Returns TRUE or FALSE
inflated <- function(x, threshold = 0.9) {
  if (sum(x == 0) >= 0.01 * length(x)) {
    d1 <- density(x)
    d2 <- density(x[x != 0], bw = d1$bw, from = min(d1$x), to = max(d1$x))
    z <- which.min(abs(d1$x))
    d2$y[z] / d1$y[z] < threshold  # Arbitrary threshold for detecting zero-inflated distribution
  } else {
    FALSE
  }
}

#------------------

# Function to integerize real (non-integer) positive weights
# 'mincor' refers to the minimum allowable Pearson correlation between 'x' and the integerized version of 'x'
# Function will also handle 'x' that is constant or already integer
# integerize <- function(x, mincor = 0.999) {
#   stopifnot(all(x > 0))
#   if (sd(x) == 0) {
#     return(rep(1L, length(x)))
#   } else {
#     p <- 0
#     i <- 0
#     r <- max(x) / min(x)
#     while (p < mincor) {
#       i <- i + 1
#       mx <- ifelse(is.integer(x), r, max(r, 10 ^ i))
#       z <- 1 + mx * ((x - min(x)) / r)
#       z <- as.integer(round(z))
#       p <- cor(x, z)
#     }
#     return(z)
#   }
# }

# Examples
# x <- rlnorm(1e3)
# xint <- integerize(x)
# cor(x, xint)
#
# x <- 1:10
# xint <- integerize(x)
# all.equal(x, xint)
#
# x <- rep(0.1, 10)
# xint <- integerize(x)
# unique(xint)

#------------------

# Convert data frame to 'dgCMatrix' sparse matrix for use by lightgbm
# Convert factors to numeric, by reference (efficient)
# Sets minimal categorical integer value to zero
# Converts to Matrix class 'dgCMatrix'
# See here: https://www.gormanalysis.com/blog/sparse-matrix-construction-and-use-in-r/
tomat <- function(data) {
  dmat <- as.data.table(data)
  for (v in names(dmat)) {
    if (is.factor(dmat[[v]])) set(dmat, i = NULL, j = v, value = as.integer(dmat[[v]]) - 1L)
    if (is.logical(dmat[[v]])) set(dmat, i = NULL, j = v, value = as.integer(dmat[[v]]))
  }
  dmat <- as(as.matrix(dmat), "dgCMatrix")
  return(dmat)
}

#------------------

# Normalize continuous variables prior to KNN
# Assumes that there are no NA's in 'x'
normalize <- function(x, center, scale, eps = 0.001) {
  y <- (x - center) / scale
  z <- (y - qnorm(eps)) / (2 * qnorm(1 - eps))
  return(z)
}

#------------------

# One-hot encoding of data.table for use in fitting LASSO models
# Based on: https://github.com/ben519/mltools/blob/master/R/one_hot.R
# Note that some arguments have been dropped and their original default values are assumed

one_hot <- function(dt, sparse_matrix = TRUE, sparsifyNAs = FALSE, naCols = FALSE) {

  if (!is.data.table(dt)) dt <- as.data.table(dt)

  # Convert ordered factors and logicals to integer
  dt <- mutate_if(dt, is.ordered, as.integer)
  dt <- mutate_if(dt, is.logical, as.integer)

  OHEID <- NULL
  cols <- colnames(dt)[which(sapply(dt, function(x) is.factor(x) & !is.ordered(x)))]
  if (length(cols) == 0) return(dt)
  tempDT <- dt[, cols, with = FALSE]
  tempDT[, `:=`(OHEID, .I)]
  for (col in cols) set(tempDT, j = col, value = factor(paste(col, tempDT[[col]], sep = ".."), levels = paste(col, levels(tempDT[[col]]), sep = "..")))
  melted <- data.table::melt(tempDT, id = "OHEID", value.factor = T, na.rm = TRUE)
  newCols <- data.table::dcast(melted, OHEID ~ value, drop = T, fun.aggregate = length)
  newCols <- newCols[tempDT[, list(OHEID)]]
  newCols[is.na(newCols[[2]]), `:=`(setdiff(paste(colnames(newCols)), "OHEID"), 0L)]
  if (!sparsifyNAs | naCols) {
    na_cols <- character(0)
    for (col in cols) if (any(is.na(tempDT[[col]])))
      na_cols <- c(na_cols, col)
    if (!sparsifyNAs)
      for (col in na_cols) newCols[is.na(tempDT[[col]]), `:=`(intersect(levels(tempDT[[col]]), colnames(newCols)), NA_integer_)]
    if (naCols)
      for (col in na_cols) newCols[, `:=`(eval(paste0(col, "_NA")), is.na(tempDT[[col]]) * 1L)]
  }
  result <- cbind(dt, newCols[, !"OHEID"])
  possible_colnames <- character(0)
  for (col in colnames(dt)) {
    possible_colnames <- c(possible_colnames, col)
    if (col %in% cols) {
      possible_colnames <- c(possible_colnames, paste0(col, "_NA"))
      possible_colnames <- c(possible_colnames, paste(levels(tempDT[[col]])))
    }
  }
  sorted_colnames <- intersect(possible_colnames, colnames(result))
  setcolorder(result, sorted_colnames)
  result <- result[, !cols, with = FALSE]

  if (sparse_matrix) result <- as(as.matrix(result), "dgCMatrix")

  return(result)
}

#------------------

# Create stratified training set or cross-validation fold assignment
# ycont: logical. Should 'y' be treated as continuous?
# tfrac: Either a fraction of training data to retain (if less than 1) or the number of CV folds if greater than 1
# ntiles: Number of buckets to break continuous 'y' into for stratified sampling
# cv_list: If TRUE, a list of length 'tfrac' is returned with indices of fold assignment; otherwise, a single vector of length(y) with values 1:tfrac giving the fold assignment

stratify <- function(y, ycont, tfrac, ntiles, cv_list = FALSE) {
  stopifnot((tfrac > 0 & tfrac <= 1) | (tfrac > 1 & tfrac %% 1 == 0))
  if (ycont) y <- dplyr::ntile(y, ntiles)
  if (tfrac <= 1) { # Training set indicator (logical)
    out <- vector(mode = "logical", length = length(y))
    for (i in unique(y)) {
      ind <- y == i
      N <- sum(ind)
      out[ind] <- data.table::frank(runif(N), ties.method = "random") <= round(tfrac * N)
    }
  }
  if (tfrac > 1) {  # Cross-validation folds (list of integer row indices)
    out <- vector(mode = "integer", length = length(y))
    for (i in unique(y)) {
      ind <- y == i
      out[ind] <- dplyr::ntile(runif(sum(ind)), tfrac)
    }
    if (cv_list) {
      out <- data.table::as.data.table(out)[, list(list(.I)), by = out]
      out <- out$V1
    }
  }
  return(out)
}
