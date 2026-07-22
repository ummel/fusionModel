#' Clean a Numeric Vector by Digit Reduction and Integer Coercion
#'
#' @description
#' `cleanNumeric()` prepares a numeric vector for memory efficiency and optimal modeling
#' compatibility. It attempts to safely coerce numeric values to standard R integers if
#' doing so preserves precision within specified tolerance thresholds.
#'
#' @param x A numeric vector to clean.
#' @param tol Numeric scalar specifying the allowable relative percentage difference when
#'   reducing precision via \code{\link{signifDigits}}. Default is `0.001`.
#' @param minimize Logical. If `TRUE`, attempts Z-score transformation first to minimize
#'   the count of distinct values. Default is `FALSE`.
#' @param threshold Numeric scalar between 0 and 1 indicating the required fraction of
#'   integer-coercible values needed to trigger integer conversion. Default is `0.999`.
#'
#' @return A cleaned numeric vector, converted to integer where possible or rounded to
#'   significant digits.
#'
#' @keywords internal
#' @noRd
cleanNumeric <- function(x, tol = 0.001, minimize = FALSE, threshold = 0.999) {
  x <- convertInteger(x, threshold = 1)  # Coerces to integer, if possible
  if (is.double(x)) {
    x <- signifDigits(x, tol = tol, minimize = minimize)
    x <- convertInteger(x, threshold = threshold)
  }
  return(x)
}

#------------------

#' Round Numeric Vector to Relative Significant Digits
#'
#' @description
#' Rounds a numeric vector to a dynamic number of significant digits such that all returned
#' values fall within a relative error tolerance of the original inputs.
#'
#' @param x A numeric vector.
#' @param tol Numeric scalar defining the maximum allowable relative deviation
#'   `abs(approx - orig) / abs(orig)`. Default is `0.001`.
#' @param minimize Logical. If `TRUE`, evaluates scale transformations (Z-scores) to find
#'   a representation that minimizes the total count of unique values. Default is `FALSE`.
#'
#' @return A numeric vector rounded to allowable significant digits.
#'
#' @keywords internal
#' @noRd
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
    if (data.table::uniqueN(x1) <= data.table::uniqueN(x2)) return(x1) else return(x2)
  }

}

#------------------

#' Conditionally Convert Numeric Vector to Integer
#'
#' @description
#' Evaluates whether values in a numeric vector fall within 32-bit signed integer bounds
#' (`.Machine$integer.max`) and satisfy an integer-coercion threshold fraction.
#'
#' @param x A numeric or logical vector.
#' @param threshold Numeric scalar between 0 and 1 representing the minimum proportion of
#'   values that must be whole numbers to trigger conversion. Default is `0.99`.
#'
#' @return An integer vector if conditions are satisfied; otherwise the original vector or
#'   logical vector (if all values are `NA`).
#'
#' @keywords internal
#' @noRd
convertInteger <- function(x, threshold = 0.99) {
  if (collapse::allNA(x)) {
    x <- as.logical(x)
  } else {
    ok32 <- max(x, na.rm = TRUE) <= .Machine$integer.max
    if (ok32) {
      chk <- collapse::na_rm(x) %% 1 == 0
      if (sum(chk) / length(chk) >= threshold) {
        x <- as.integer(round(x))
      }
    }
  }
  return(x)
}

#-------------------

#' Check if a Vector Lacks Variation
#'
#' @description
#' Identifies whether a vector contains zero non-trivial variance (i.e., has at most one
#' unique non-`NA` value or is entirely `NA`).
#'
#' @param x An atomic vector.
#'
#' @return Logical scalar (`TRUE` if the vector has <= 1 unique non-`NA` value).
#'
#' @keywords internal
#' @noRd
novary <- function(x) data.table::uniqueN(x, na.rm = TRUE) <= 1

#------------------

#' Calculate Imputation Replacement Values
#'
#' @description
#' Derives replacement values for missing data points (`NA`). Uses median value replacement
#' for numeric vectors and weighted random sampling based on non-missing frequency
#' distributions for categorical vectors.
#'
#' @param x Vector containing missing values.
#' @param na.ind Logical vector indicating location of `NA` values (`TRUE` for missing).
#'
#' @return Imputed values matching the count of `NA` occurrences (`sum(na.ind)`).
#'
#' @keywords internal
#' @noRd
imputationValue <- function(x, na.ind) {
  if (is.numeric(x)) {
    m <- median(x, na.rm = TRUE)
    m <- ifelse(is.integer(x), as.integer(round(m)), m)
  } else {
    tab <- table2(x) / sum(!na.ind)
    m <- sample(names(tab), size = sum(na.ind), replace = TRUE, prob = tab)
  }
  return(m)
}

#------------------

#' Compare Equality of Variable Classes
#'
#' @description
#' Evaluates whether two class descriptors are identical, treating `"integer"` and `"numeric"`
#' as equivalent numeric types for data fusion operations.
#'
#' @param x Character vector representing class attributes of first object.
#' @param y Character vector representing class attributes of second object.
#'
#' @return Logical scalar indicating structural class equivalence.
#'
#' @keywords internal
#' @noRd
sameClass <- function(x, y) {
  if (x[1] == "integer") x <- "numeric"
  if (y[1] == "integer") y <- "numeric"
  identical(x, y)
}

#------------------

#' Fast Weighted Mean Computation
#'
#' @description
#' Efficient wrapper around \code{\link[matrixStats]{weightedMean}} for rapid weighted average
#' calculations.
#'
#' @param x Numeric vector of values.
#' @param w Numeric vector of non-negative weights matching length of `x`.
#'
#' @return Numeric scalar representing weighted arithmetic mean.
#'
#' @keywords internal
#' @noRd
wmean <- function(x, w) {
  matrixStats::weightedMean(x, w)
}

#------------------

#' Weighted Standard Deviation Calculation
#'
#' @description
#' Computes sample standard deviation incorporating observation weights. Results match
#' `Hmisc::wtd.var()` with normalized weighting enabled (`normwt = TRUE`).
#'
#' @param x Numeric vector of sample observations.
#' @param w Numeric vector of weights.
#'
#' @return Numeric scalar standard deviation.
#'
#' @keywords internal
#' @noRd
wsd <- function(x, w) {
  if (anyNA(x)) {
    i <- !is.na(x)
    x <- x[i]
    w <- w[i]
  }
  w <- (w / sum(w)) * length(x)
  sw <- sum(w)
  xbar <- sum(w * x) / sw
  sqrt(sum(w * ((x - xbar) ^ 2)) / (sw - 1))
}

#------------------

#' Detect Zero-Inflated Numeric Distribution
#'
#' @description
#' Evaluates density profile comparison to classify if a numeric variable exhibits
#' structural zero inflation.
#'
#' @param x Vector to evaluate.
#' @param threshold Numeric scalar threshold for density ratio test. Default is `0.9`.
#'
#' @return Logical scalar (`TRUE` if distribution is zero-inflated).
#'
#' @keywords internal
#' @noRd
inflated <- function(x, threshold = 0.9) {
  if (is.numeric(x)) {
    if (sum(x == 0) >= 0.01 * length(x)) {
      d1 <- density(x)
      d2 <- density(x[x != 0], bw = d1$bw, from = min(d1$x), to = max(d1$x))
      z <- which.min(abs(d1$x))
      d2$y[z] / d1$y[z] < threshold  # Threshold for detecting zero-inflated distribution
    } else {
      FALSE
    }
  } else {
    FALSE
  }
}

#-------------------

#' Convert Real Positive Weights to Integers
#'
#' @description
#' Iteratively scales and rounds continuous positive sampling weights into integers while
#' preserving a target Pearson correlation threshold relative to original continuous values.
#'
#' @param x Numeric vector of positive weights.
#' @param mincor Numeric scalar specifying target minimum allowable correlation between original
#'   and integerized weights. Default is `0.999`.
#'
#' @return Integer vector of integerized weights.
#'
#' @keywords internal
#' @noRd
integerize <- function(x, mincor = 0.999) {
  stopifnot(all(x > 0))
  if (sd(x) == 0) {
    return(rep(1L, length(x)))
  } else {
    p <- 0
    i <- 0
    r <- max(x) / min(x)
    while (p < mincor) {
      i <- i + 1
      mx <- ifelse(is.integer(x), r, max(r, 10 ^ i))
      z <- 1 + mx * ((x - min(x)) / r)
      z <- as.integer(round(z))
      p <- cor(x, z)
    }
    return(z)
  }
}

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

#' Convert Data Frame to Matrix Representation
#'
#' @description
#' Efficiently converts factors and logical columns in a data frame or `data.table` to integer
#' values before casting into a numeric matrix.
#'
#' @param data Data frame or `data.table` to transform.
#'
#' @return Standard R numeric matrix object.
#'
#' @keywords internal
#' @noRd
to_mat <- function(data) {
  dmat <- as.data.table(data)
  for (v in names(dmat)) {
    if (is.factor(dmat[[v]])) set(dmat, i = NULL, j = v, value = as.integer(dmat[[v]]) - 1L)
    if (is.logical(dmat[[v]])) set(dmat, i = NULL, j = v, value = as.integer(dmat[[v]]))
  }
  dmat <- as.matrix(dmat)
  return(dmat)
}

#------------------

#' Normalize Continuous Variables for Distance Computation
#'
#' @description
#' Centers and scales input vectors before standard normal quantile compression into a
#' bounded unit metric space.
#'
#' @param x Numeric vector to normalize.
#' @param center Location shift parameter (typically mean or median).
#' @param scale Scale parameter (typically standard deviation or MAD).
#' @param eps Small offset numeric scalar bound preventing infinite quantile estimates.
#'   Default is `0.001`.
#'
#' @return Bounded, scaled numeric vector.
#'
#' @keywords internal
#' @noRd
normalize <- function(x, center, scale, eps = 0.001) {
  y <- (x - center) / scale
  z <- (y - qnorm(eps)) / (2 * qnorm(1 - eps))
  return(z)
}

#' Denormalize Bounded Values Back to Original Scale
#'
#' @description
#' Reverses continuous transformation applied by \code{\link{normalize}}.
#'
#' @param z Normalized numeric vector.
#' @param center Original centering location shift parameter used in normalization.
#' @param scale Original scaling parameter used in normalization.
#' @param eps Epsilon boundary offset used during initial normalization. Default is `0.001`.
#'
#' @return Numeric vector in original units and scale.
#'
#' @keywords internal
#' @noRd
denormalize <- function(z, center, scale, eps = 0.001) {
  y <- 2 * z * qnorm(1 - eps) + qnorm(eps)
  y * scale + center
}

#------------------

#' Efficient One-Hot Encoding of Factor Variables
#'
#' @description
#' Encodes factor columns in a data frame into binary dummy indicator columns using optimized
#' sparse matrix representations.
#'
#' @param data A `data.frame` or `data.table` containing columns to encode.
#' @param dropOriginal Logical. If `TRUE`, drops source factor columns after encoding. Default
#'   is `TRUE`.
#' @param dropUnusedLevels Logical. If `TRUE`, removes factor levels with zero occurrence in the
#'   data. Default is `FALSE`.
#'
#' @return Object of same class as `data` containing additional dummy columns, with a link
#'   data frame attached under the `"one_hot_link"` attribute.
#'
#' @keywords internal
#' @noRd
one_hot <- function(data, dropOriginal = TRUE, dropUnusedLevels = FALSE) {

  stopifnot(is.data.frame(data))
  dt <- inherits(data, "data.table")
  data <- as.data.table(data)

  # Identify factor variables in 'data' to hot encode
  y <- names(which(sapply(data, is.factor)))

  # Proceed only if there are factors to one-hot encode
  if (length(y) > 0) {

    dy <- data[, ..y]
    if (dropOriginal) data[, c(y) := NULL]

    # Construct levels linkage
    dlink <- lapply(y, function(v) {
      lev <- levels(dy[[v]])
      data.frame(original = v, dummy = paste(v, lev, sep = ""), level = lev)
    })
    dlink <- do.call(rbind, dlink)

    # Matrix package implementation (sparse matrix construction)
    dy <- Matrix::sparse.model.matrix(~ 0 + .,
                                      data = dy,
                                      contrasts.arg = lapply(dy, contrasts, contrasts = FALSE),
                                      drop.unused.levels = FALSE)

    dy <- as.matrix(dy)

    # Drop dummy columns with no 1's
    if (dropUnusedLevels) {
      keep <- colSums(dy) > 0
      dy <- dy[, keep]
    }

    data <- cbind(data, dy)
    if (!dt) data <- as.data.frame(data)
    attr(data, "one_hot_link") <- dlink
  }
  return(data)

}

#------------------

#' Stratified Partitioning and Cross-Validation Fold Assignment
#'
#' @description
#' Generates stratified subsamples or K-fold cross-validation fold assignments while preserving
#' outcome variable empirical distribution profiles.
#'
#' @param y Outcome vector used for distribution stratification.
#' @param ycont Logical indicating whether `y` is continuous (`TRUE`) or discrete/categorical
#'   (`FALSE`).
#' @param tfrac Sampling fraction (if <= 1) or integer count of cross-validation folds (if > 1).
#' @param ntiles Number of quantile bins to bucket continuous `y` into prior to stratified sampling.
#' @param cv_list Logical. If `TRUE` and `tfrac > 1`, returns list of fold row indices; otherwise
#'   returns vector of integer fold tags. Default is `FALSE`.
#'
#' @return Logical selection vector (if `tfrac <= 1`), integer fold vector, or list of index
#'   vectors (if `cv_list = TRUE`).
#'
#' @keywords internal
#' @noRd
stratify <- function(y, ycont, tfrac, ntiles, cv_list = FALSE) {
  stopifnot((tfrac > 0 & tfrac <= 1) | (tfrac > 1 & tfrac %% 1 == 0))
  if (ycont) y <- dplyr::ntile(y, ntiles)
  if (tfrac <= 1) { # Create training set (logical vector indicating training set observations)
    out <- vector(mode = "logical", length = length(y))
    for (i in unique(y)) {
      ind <- y == i
      N <- sum(ind)
      out[ind] <- data.table::frank(runif(N), ties.method = "random") <= round(tfrac * N)
    }
  }
  if (tfrac > 1) {  # Assign cross-validation folds (list of integer row indices in each fold)
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

#------------------

#' Validate and Clean Modeling Data Frame Inputs
#'
#' @description
#' Checks specified target and predictor variables within a dataset. Ensures character columns
#' are converted, drops zero-variance columns, removes sparse zero-inflated targets relative
#' to CV fold settings, and optionally imputes missing values.
#'
#' @param data Input data frame containing features and targets.
#' @param y Character vector of target variable names.
#' @param x Character vector of predictor variable names.
#' @param nfolds Optional integer number of CV folds used for checking minimal non-zero count
#'   thresholds. Default is `NULL`.
#' @param impute Logical. If `TRUE`, imputes missing predictor values using \code{\link{imputationValue}}.
#'   Default is `FALSE`.
#'
#' @return Cleaned data frame with invalid variables dropped and missing values optionally imputed.
#'
#' @keywords internal
#' @noRd
checkData <- function(data, y, x, nfolds = NULL, impute = FALSE) {

  # Check for character-type variables; stop with error if any detected
  xc <- sapply(data[c(x, y)], is.character)
  if (any(xc)) stop("Coerce character variables to factor:\n", paste(names(which(xc)), collapse = ", "))

  # Remove (with warning) any zero-variance 'y' variables
  constant <- names(which(sapply(data[y], novary)))
  if (length(constant)) {
    y <- setdiff(y, constant)
    data <- select(data, -all_of(constant))
    warning("Removed zero-variance 'y' variable(s):\n", paste(constant, collapse = ", "), immediate. = TRUE)
  }

  # Remove (with warning) any zero-inflated 'y' variables with fewer than 'nfolds' non-zero values
  if (all(!is.null(nfolds), nfolds > 0)) {
    toofew <- names(which(sapply(data[y], function(x) if (inflated(x)) sum(x != 0) < nfolds else FALSE)))
    if (length(toofew)) {
      data <- select(data, -all_of(toofew))
      warning("Removed 'y' variable(s) with too few non-zero values:\n", paste(toofew, collapse = ", "), immediate. = TRUE)
    }
  }

  # Remove (with message) any zero-variance 'x' variables
  constant <- names(which(sapply(data[x], novary)))
  if (length(constant)) {
    x <- setdiff(x, constant)
    data <- select(data, -all_of(constant))
    cli::cli_alert_info("Removed zero-variance 'x' variable(s): {paste(constant, collapse = ', ')}")[cite: 1]
  }

  # Detect and impute any missing values in 'x' variables
  if (impute) {
    na.cols <- names(which(sapply(data[x], anyNA)))
    if (length(na.cols) > 0) {
      for (j in na.cols) {
        xj <- data[[j]]
        ind <- is.na(xj)
        data[ind, j] <-  imputationValue(xj, ind)
      }
    }
  }

  return(data)

}

#------------------

#' Compute Memory Size of R Object in Megabytes
#'
#' @description
#' Evaluates system memory overhead assigned to an object.
#'
#' @param x Any R object.
#'
#' @return Numeric size in megabytes (MB).
#'
#' @keywords internal
#' @noRd
objectSize <- function(x) as.numeric(utils::object.size(x)) / 1048576

#------------------

#' Query System Free Physical Memory
#'
#' @description
#' Queries OS system utilities across Windows, macOS, Linux, and HPC SLURM environments to estimate
#' available system RAM.
#'
#' @return Numeric scalar estimating available memory in megabytes (MB).
#'
#' @keywords internal
#' @noRd
freeMemory <- function() {
  gc()
  sys <- Sys.info()["sysname"]
  # Windows
  if (sys == "Windows") {
    x <- system2("wmic", args = "OS get FreePhysicalMemory /Value", stdout = TRUE)
    x <- x[grepl("FreePhysicalMemory", x)]
    x <- gsub("FreePhysicalMemory=", "", x, fixed = TRUE)
    x <- gsub("\r", "", x, fixed = TRUE)
    as.numeric(x) / 1e3
  } else {
    # Mac OS
    if (sys == "Darwin") {
      x <- system("vm_stat", intern = TRUE)
      pagesize <- x[grepl("Mach Virtual Memory Statistics",x)]
      pagesize <- gsub("Mach Virtual Memory Statistics: (page size of", "", pagesize, fixed = TRUE)
      pagesize <- gsub("bytes)", "", pagesize, fixed = TRUE)
      x <- x[grepl("Pages free: ", x)]
      x <- gsub("Pages free: ", "", x, fixed = TRUE)
      x <- gsub(".", "", x, fixed = TRUE)
      as.numeric(x) * as.numeric(pagesize) / (1024 ^ 2)
    } else {
      ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
      if (is.na(ncores)) {
        # Linux system assumed as backstop
        x <- system('grep MemAvailable /proc/meminfo', intern = TRUE)
        x <- strsplit(x, "\\s+")[[1]][2]
        as.numeric(x) / 1024
      } else {
        # Yale HPC setting
        ncores * as.integer(Sys.getenv("SLURM_MEM_PER_CPU"))
      }
    }
  }
}

#------------------

# # NOT USED
# # Structured downsampling of 2-dimensional (x, y) inputs
# # Compared to random sampling, guarantees a more uniform sampling of the space and oversampling of unusual observations
# # Returns adjusted weights such that total cluster weights are respected in the downsample
# # Will return no less than N observations, but can return significantly more due to oversampling of small/unusual clusters
#
# # Example usage
# # library(fusionModel)
# # x <- recs$electricity
# # y <- recs$square_feet
# # test <- downsample(x = recs$electricity, y = recs$square_feet)
# # dim(test)
# # plot(recs$electricity, recs$square_feet)
# # points(test[1:2], col = 2)
# # summary(test$w)
#
# downsample <- function(x,
#                        y,
#                        w = rep(1, length(x)),
#                        N = 0.1 * length(x),
#                        K = 30,
#                        min_samp = 30) {
#
#   stopifnot({
#     length(x) == length(y)
#     length(x) == length(w)
#     !anyNA(c(x, y, w))
#     N >= 1
#     K >= 1
#     min_samp >= 1
#   })
#
#   N0 <- length(x)
#
#   # Keep copies of inputs
#   x0 <- x
#   y0 <- y
#
#   # Scale x and y prior to clustering
#   x <- scale(x)
#   y <- scale(y)
#
#   # Average required adjustment factor
#   adj <- N / N0
#
#   # k-means
#   km <- kmeans(cbind(x, y), centers = K, iter.max = K * 5)
#
#   cnt <- apply(km$centers, 2, weighted.mean, w = km$size)
#
#   # Distance of each cluster center from the dataset center
#   dc <- sqrt((km$centers[, 1] - cnt[1]) ^ 2 + (km$centers[, 2] - cnt[2]) ^ 2)
#   ofct <- dc / weighted.mean(dc, km$size)
#   cadj <- adj * ofct  # Cluster-specific adjustment factor
#
#   # Clustered downsample
#   samp <- sapply(1:nrow(km$centers), function(cl) {
#     i <- which(km$cluster == cl)
#     cnt <- km$centers[cl, ]  # Center of the cluster
#     d <- sqrt((x[i] - cnt[1]) ^ 2 + (y[i] - cnt[2]) ^ 2)
#     j <- seq(from = 1, to = length(d), length.out = max(min_samp, length(d) * cadj[cl]))
#     j <- unique(round(j))
#     k <- i[order(d)][j]
#     wout <- w[k] * (sum(w[i]) / sum(w[k]))
#     cbind(x = x0[k], y = y0[k], w = wout)
#   })
#
#   samp <- do.call(rbind, samp)
#   return(as.data.frame(samp))
#
# }

#------------------

#' Univariate K-Means Clustering for Continuous Variables
#'
#' @description
#' Partitions a continuous vector into ordered factor bins via deterministic 1D k-means clustering.
#'
#' @param x Continuous numeric vector or ordered factor.
#' @param k Integer number of target clusters.
#'
#' @return An ordered factor vector of cluster assignments ranging from 1 to `k`.
#'
#' @keywords internal
#' @noRd
uniCluster <- function(x, k) {
  stopifnot(!is.character(x))
  if (is.ordered(x)) {
    x <- as.integer(x)
  } else {
    if (is.factor(x)) stop("Unordered factor not allowed in uniCluster()")
  }
  k <- min(k, data.table::uniqueN(x))
  iter <- 0
  halt <- FALSE

  # Ensures the kmeans() results are deterministic by specifying initial centers
  ind <- round(seq(from = 1, to = data.table::uniqueN(x), length.out = k))
  initial <- sort(unique(x))[ind]
  km <- suppressWarnings(kmeans(x, centers = initial, iter.max = 50))

  cl <- data.table::frank(km$centers, ties.method = "dense")[km$cluster]
  cl <- factor(cl, levels = 1:k, ordered = TRUE)
  return(cl)
}

#------------------

#' Compute Weighted Mean and Sampling Variance
#'
#' @description
#' Estimates weighted mean/proportion and corresponding variance following Cochran (1977)
#' formulation for continuous metrics or survey proportion formulations for binary vectors.
#'
#' @param x Continuous or binary numeric vector.
#' @param w Non-negative sampling weights vector matching length of `x`.
#'
#' @return Numeric vector of length 2 containing `c(weighted_mean, variance)`.
#'
#' @keywords internal
#' @noRd
weightedVAR <- function(x, w) {

  # Continuous case
  if (any(!x %in% 0:1)) {
    n <- length(x)
    wbar <- mean(w)
    xWbar <- wmean(x, w)
    wvar <- n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
  } else {

    # Binary case
    w <- w / sum(w)
    sw2 <- sum(w ^ 2)
    xWbar <- sum(w * x)
    wvar <- xWbar * (1 - xWbar) * sw2
  }
  return(c(xWbar, wvar))
}

#------------------

#' Bivariate Conditional Quantile Smoother
#'
#' @description
#' Estimates conditional target quantiles across continuous predictors using decision tree
#' partition binning and generalized additive model (GAM) smoothing.
#'
#' @param x Predictor numeric vector.
#' @param y Target numeric outcome vector.
#' @param qu Numeric vector of target quantiles (e.g. `0.5`) or string `"mad"` to derive median
#'   and median absolute deviations. Default is `0.5`.
#' @param xout Optional numeric vector specifying prediction target locations along `x`. Default
#'   is `NULL`.
#' @param a Quantile bucket scaling exponent governing recursive tree node minimum bin counts. Default is `0.3`.
#'
#' @return Data frame containing columns `x` (evaluation values), `y` (smoothed quantile), and `q` (quantile tag).
#'
#' @keywords internal
#' @noRd
smoothQuantile <- function(x,
                           y,
                           qu = 0.5,
                           xout = NULL,
                           a = 0.3) {


  stopifnot({
    length(x) == length(y)
    qu == "mad" | all(qu > 1, qu < 1)
    is.null(xout) | length(xout) > 1
    a > 0
  })

  # Remove NA observations
  na <- is.na(x) | is.na(y)
  x <- x[!na]
  y <- y[!na]

  # Order 'x' and create other necessary inputs
  n <- length(x)
  i <- order(x)
  x <- x[i]
  y <- y[i]
  d <- data.table(x, y)
  qmad <- qu[1] == "mad"
  qu <- if (qmad) c(0.25, 0.5, 0.75) else sort(unique(qu))

  # Check if 'xout' is outside natural range of 'x'
  if (!is.null(xout)) {
    if (any(xout < min(x) | xout > max(x))) warning("'xout' contains values outside range of 'x' (be careful)")
  }

  # x-values at which to predict smoothed values
  mout <- if (is.null(xout) | length(xout) == 2) {
    rng <- if (is.null(xout)) range(x) else range(xout)
    xx <- if (is.null(xout)) x else x[x >= rng[1] & x <= rng[2]]
    mn <- diff(rng) / (0.1 * mad(xx, constant = 1))
    mn <- min(mn, 200)
    seq(from = rng[1], to = rng[2], length.out = ceiling(mn))
  } else {
    xout
  }

  f <- min(1, max(0.05, 30e3 / n))
  df <- if (f < 1) {
    ind <- round(seq(from = 1, to = n, by = (n - 1) / (round(f * n) - 1)))
    ind[length(ind)] <- n
    d[ind, ]
  } else {
    d
  }

  # Minimum number of observations per rpart bin
  mb <- max(5, ceiling(n ^ a))
  mbr <- max(5, ceiling(mb * f))

  # Winsorize y-values in 'df' with obviously extreme values prior to the rpart() call
  dr <- copy(df)
  zy <- (dr$y - median(dr$y)) / mad(dr$y)
  zy[is.na(zy)] <- 0
  ymax <- 3.5 * mad(dr$y) + median(dr$y)
  dr$y[zy > 3.5] <- ymax
  dr$y[zy < -3.5] <- -ymax

  # Fit rpart to determine bin breakpoints along 'x'
  fit <- rpart::rpart(formula = y ~ x, data = dr,
                      xval = 0, cp = 0.001,
                      maxcompete = 0, maxsurrogate = 0,
                      minbucket = mbr, minsplit = 2 * mbr)

  r <- rle(as.integer(as.factor(fit$where)))
  rsum <- c(1, cumsum(r$lengths))

  # Bin assignments and bin metrics calculation
  for (i in 1:3) {
    p <- c(0, 0.5, -0.5)[i]
    br <- pmax(1, round(p * diff(rsum) + rsum[-(length(rsum))]))
    br <- c(br, max(rsum))
    if (p > 0) br <- c(min(rsum), br)
    bin <- cut(x, breaks = unique(df$x[br]), labels = FALSE, include.lowest = TRUE)
    set(d, j = paste0("bin", i), value = bin)
  }

  d <- lapply(1:3, function(i) {
    d[, list(s = .N, m = sum(x), q = quantile(y, probs = qu, type = 1)), by = eval(paste0("bin", i))]
  })
  d <- rbindlist(d, use.names = FALSE, idcol = "window")
  d[, m := m / s]  # Mean x-value of each bin
  d[, qu := rep(qu, times = nrow(d) / length(qu))]

  # Retain only those bins with sufficient number of observations
  d <- d[s >= mb, ]

  if (qmad) {
    f <- function(x) c(x[2] - x[1], x[2], x[3] - x[2])
    d[, q := f(q), by = c("window", "bin1")]
  }

  d[m == min(m), m := min(x)]
  d[m == max(m), m := max(x)]

  # Fit GAM(s) and return the smoothed quantile values for each 'mout' x-value
  out <- sapply(qu, function(i) {
    df <- d[qu == i, ]
    fit <- mgcv::gam(formula = q ~ s(m), data = df, weights = d$s)
    predict(fit, newdata = data.table(m = mout))
  })

  # Ensure no quantile crossing
  if (ncol(out) > 1 & !qmad) out <- t(apply(out, 1, sort))

  # Final output data frame
  out <- data.frame(x = mout, y = as.vector(out), q = rep(qu, each = length(mout)))

  return(out)

}

#------------------

# Robust mean smoother
# TO DO -- create separate outliers() function that embeds in smoothMean with options to make robust on the fly

# acs <- fst::read_fst("~/Documents/Projects/fusionData/survey-processed/ACS/2019/ACS_2019_H_processed.fst")
# x <- acs$hincp
# y <- acs$elep
#
# test <- smoothQuantile(x, y, a= 0.001)
# test <- smoothMean(x, y, xout = c(0, 200e3))

# smoothMean <- function(x, y, xout = NULL,
#                        zmax = 3.5, a = 0.3, b = 0.1,
#                        se = FALSE, verbose = TRUE) {
#
#   #n <- length(x)
#   i <- order(x)
#   x <- x[i]
#   y <- y[i]
#
#   # Eliminate obviously extreme y values?
#   zy <- abs(y - median(y)) / mad(y)
#   zy[is.na(zy)] <- 0
#   iy <- zy > zmax ^ 2
#
#   ix <- rep(FALSE, length(iy))
#
#   # Initial outliers
#   oi <- iy | ix
#   x <- x[!oi]
#   y <- y[!oi]
#
#   # If 'y' is already monotonic or has no variance, simply return a linear interpolation
#   simple <- sorted(y) | var(y) == 0
#   xout <- if (simple & is.null(xout)) range(x) else xout
#
#   # Data table used below
#   df <- data.table(x, y)
#
#   # Detect and remove outliers
#   if (!simple) {
#
#     # Return smoothed median and lower and upper MAD
#     d <- smoothQuantile(x, y, qu = "mad", xout = NULL, a = a, b = b)
#     d <- split(d, d$q)
#
#     # Smoothed median for each x value
#     m <- approx(x = d[[2]]$x, y = d[[2]]$y, xout = x)$y
#
#     # Upper MAD
#     i <- d[[3]]$y > 0
#     d1 <- approx(x = d[[3]]$x[i], y = d[[3]]$y[i], xout = x)$y
#     d1[d1 == 0] <- NA
#
#     # Lower MAD
#     i <- d[[1]]$y > 0
#     d2 <- approx(x = d[[1]]$x[i], y = d[[1]]$y[i], xout = x)$y
#     d2[d2 == 0] <- NA
#
#     # Detect the y outliers
#     out1 <- y > m + zmax * d1
#     out2 <- y < m - zmax * d2
#     oy <- out1 | out2
#
#     # Report total number of outliers
#     nout <- sum(oy, na.rm = TRUE)
#     if (verbose) cat("Removed ", nout, " outliers (", paste0(signif(100 * nout / n, 3), "%"), ")\n", sep = "")
#
#     # Remove the outlier observations
#     df <- df[!oy]
#
#   }
#
#   # x-values at which to predict smoothed values
#   xout <- if (is.null(xout)) {
#     d[[1]]$x  # Use the x-values returned by quantile smoother
#   } else {
#     if (length(xout) == 2) {
#       xx <- x[x >= min(xout) & x <= max(xout)]
#       rng <- range(xx)
#       mn <- diff(rng) / (b * mad(xx, constant = 1))
#       mn <- min(mn, 200)  # Hard max on number of output x-values
#       seq(from = rng[1], to = rng[2], length.out = ceiling(mn))
#     } else {
#       unique(sort(xout))
#     }
#   }
#
#   # Final predicted smooth output
#   result <- if (simple) {
#
#     out <- data.frame(x = xout,
#                       y = suppressWarnings(approx(x = df$x, y = df$y, xout = xout)$y),
#                       se = 0)
#     if (!se) out$se <- NULL
#     out
#
#   } else {
#
#     # Fit smoothed mean (gam or bam)
#     big <- nrow(df) > 30e3
#     fun <- ifelse(big, mgcv::bam,  mgcv::gam)
#     if (verbose) cat("Fitting smoothed mean using", ifelse(big, "mgcv::bam()", "mgcv::gam()"), "\n")
#     fit <- fun(y ~ s(x), data = df, discrete = big)
#     data.frame(x = xout, predict(fit, data.table(x = xout), se.fit = se))
#
#   }
#
#   # Clean up and return result
#   names(result) <- c("x", "y", if (se) "se" else NULL)
#   row.names(result) <- NULL
#
#   return(result)
#
# }

#------------------

# Create 'maxr_fun' function used inside validate() to derive the cut-off point for ratio 'r' (ubar / b)
# CI's using 'r' larger than that returned by maxr_fun() are invalid
# The cutoff point increases with the number of implicates

# mseq <- seq(2, 100, by = 2)
# out <- sapply(mseq, function(m, p = 0.95) {
#
#   # Ratio of ubar to b
#   R <- seq(0, 1, length.out = 1e4)
#
#   # Degrees of freedom
#   df <- (m - 1) * (1 - R * m  / (m + 1)) ^ 2
#
#   # Calculate SE and associated CI
#   ubar <- 1
#   se <- sqrt(ubar / R * (1 + 1 / m) - ubar)
#   ci <- se * qt(p, df)
#
#   R[which.min(ci)]
#
# })
#
# maxr_fun <- approxfun(mseq, out, rule = 2)
# rm(mseq, out)

#------------------

#' Check Vector Monotonicity
#'
#' @description
#' Tests whether vector values are sorted monotonously (either strictly non-decreasing or
#' non-increasing within floating point limits).
#'
#' @param x Numeric or atomic vector.
#'
#' @return Logical scalar (`TRUE` if sorted monotonically).
#'
#' @keywords internal
#' @noRd
sorted <- function(x) {
  xd <- diff(x)
  all(xd >= -sqrt(.Machine$double.eps)) | all(xd <= sqrt(.Machine$double.eps))
}

#------------------

#' Retrieve Function Call with Preserved Formals Defaults
#'
#' @description
#' Extends \code{\link[base]{match.call}} to preserve default function arguments that were not
#' explicitly declared in the call stack.
#'
#' @param ... Passed to internal call expansion.
#' @param exclude Character vector of parameter argument names to exclude from returned call.
#'   Default is `NULL`.
#'
#' @return Language call object.
#'
#' @keywords internal
#' @noRd
match.call.defaults <- function(..., exclude = NULL) {
  call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
  formals <- evalq(formals(), parent.frame(1))
  for (i in setdiff(names(formals), c(names(call), exclude)))
    call[i] <- list(formals[[i]])
  match.call(sys.function(sys.parent()), call)
}

#-------------------

#' Standardize Absolute File Path Expressions
#'
#' @description
#' Normalizes system file paths using current platform-specific directory separators
#' (`.Platform$file.sep`).
#'
#' @param path Character vector of path locations.
#' @param mustWork Logical indicating whether to throw an error if path does not exist. Default is `NA`.
#'
#' @return Character vector of normalized file path targets.
#'
#' @keywords internal
#' @noRd
full.path <- function(path, mustWork = NA) {
  normalizePath(path = path, winslash = .Platform$file.sep, mustWork = mustWork)
}

#-------------------

#' Fast Weighted Frequency Counting
#'
#' @description
#' High-performance weighted tabulation using `data.table` grouping structures.
#'
#' @param x Atomic vector of values to aggregate.
#' @param w Optional numeric weights vector matching length of `x`. Default is `NULL`.
#' @param na.rm Logical. If `TRUE`, drops `NA` categories from returned counts. Default is `FALSE`.
#'
#' @return Named numeric vector of frequency totals.
#'
#' @keywords internal
#' @noRd
table2 <- function(x, w = NULL, na.rm = FALSE) {
  require(data.table)
  stopifnot(is.atomic(x))
  if (is.null(w)) {
    ds <- setDT(list(x = x), key = "x")
    ds <- ds[, .N, by = "x"]
  } else {
    stopifnot(is.numeric(w) & length(w) == length(x))
    ds <- setDT(list(x = x, w = w), key = "x")
    ds <- ds[, .(N = sum(w)), by = "x"]
  }
  if (na.rm) ds <- na.omit(ds)
  return(setNames(ds$N, ds$x))
}

#-------------------

#' Consolidate Low-Frequency Factor Levels
#'
#' @description
#' Groups infrequent factor levels into a consolidated level label (e.g., `"_Other_"`), respecting
#' optional observation weighting.
#'
#' @param x Factor vector.
#' @param nmax Maximum integer count of distinct factor levels to retain. Default is `5`.
#' @param w Optional numeric weight vector matching length of `x`. Default is `NULL`.
#' @param other_level String label assigned to consolidated sparse levels. Default is `"_Other_"`.
#'
#' @return Unordered factor vector containing at most `nmax` distinct levels.
#'
#' @keywords internal
#' @noRd
lumpFactor <- function(x, nmax = 5, w = NULL, other_level = "_Other_") {
  stopifnot(is.factor(x))
  if (is.null(w)) w <- rep.int(1L, length(x))
  stopifnot(is.numeric(w))
  if (data.table::uniqueN(x, na.rm = TRUE) >= nmax) {
    p <- table2(x, w = w, na.rm = TRUE)
    p <- sort(p, decreasing = TRUE) / sum(p)
    keep <- intersect(levels(x), names(p)[1:(nmax - 1)])  # Retains original factor ordering
    x <- as.character(x)
    x[!x %in% keep] <- other_level
    x <- factor(x, levels = c(keep, other_level))
  }
  return(x)
}

#-------------------

#' Check for Infinite Numeric Values
#'
#' @description
#' Efficiently determines whether any vector elements contain `Inf` or `-Inf` values.
#'
#' @param x Atomic vector to evaluate.
#'
#' @return Logical scalar (`TRUE` if infinite values are detected).
#'
#' @keywords internal
#' @noRd
anyInf <- function(x) {
  if (is.double(x)) collapse::anyv(x, Inf) else FALSE
}
