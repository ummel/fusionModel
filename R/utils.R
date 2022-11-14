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
  if (all(x[!is.na(x)] %% 1 == 0) & max(x, na.rm = TRUE) <= .Machine$integer.max) {
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
  matrixStats::weightedMean(x, w)
}

#------------------

# Weighted standard deviation
# Equivalent to Hmisc::wtd.var() with normwt = TRUE and taking sqrt() of result
# Equivalent to na.rm = TRUE
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

# Detect if a numeric variable is likely to be zero-inflated
# Returns TRUE or FALSE
inflated <- function(x, threshold = 0.9) {
  if (is.numeric(x)) {
    if (sum(x == 0) >= 0.01 * length(x)) {
      d1 <- density(x)
      d2 <- density(x[x != 0], bw = d1$bw, from = min(d1$x), to = max(d1$x))
      z <- which.min(abs(d1$x))
      d2$y[z] / d1$y[z] < threshold  # Arbitrary threshold for detecting zero-inflated distribution
    } else {
      FALSE
    }
  } else {
    FALSE
  }
}

#-------------------

# Function to integerize real (non-integer) positive weights
# 'mincor' refers to the minimum allowable Pearson correlation between 'x' and the integerized version of 'x'
# Function will also handle 'x' that is constant or already integer

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

# Convert data frame to 'dgCMatrix' sparse matrix for use by lightgbm
# Convert factors to numeric, by reference (efficient)
# Sets minimal categorical integer value to zero
# Converts to Matrix class 'dgCMatrix'
# See here: https://www.gormanalysis.com/blog/sparse-matrix-construction-and-use-in-r/
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

# Normalize continuous variables prior to KNN
# Assumes that there are no NA's in 'x'
normalize <- function(x, center, scale, eps = 0.001) {
  y <- (x - center) / scale
  z <- (y - qnorm(eps)) / (2 * qnorm(1 - eps))
  return(z)
}

# Back-transform normalized values
denormalize <- function(z, center, scale, eps = 0.001) {
  y <- 2 * z * qnorm(1 - eps) + qnorm(eps)
  y * scale + center
}

#------------------

# Function to efficiently one-hot encode factor variables in a data frame
# Based on suggestion here: https://stackoverflow.com/questions/39905820/how-to-one-hot-encode-factor-variables-with-data-table
# If dt = TRUE, a data.table is returned; otherwise, data.frame
# If dropOriginal = TRUE, the original factor columns are dropped
# If dropUnusedLevels = TRUE, unused factor levels are dropped
# Adds attribute "one_hot_link" to output data.frame/table

one_hot <- function(data, y = NULL, dt = TRUE, dropOriginal = TRUE, dropUnusedLevels = FALSE) {

  stopifnot(is.data.frame(data))
  if (dt) data <- as.data.table(data)

  # If no 'y' provided, one-hot encode all factors in 'data'
  if (is.null(y)) y <- names(which(sapply(data, is.factor)))

  # Proceed only if there are factors to one-hot encode
  if (length(y) > 0) {

    # Subset variables depending on input type
    if (class(data)[1] == "data.frame") {
      dy <- data[y]
      if (dropOriginal) data <- data[setdiff(names(data), y)]
    } else {
      dy <- data[, ..y]
      if (dropOriginal) data[, c(y) := NULL]
    }

    # Construct levels linkage
    # Must be consistent with how Matrix::sparse.model.matrix names the dummy variables
    dlink <- lapply(y, function(v) {
      lev <- levels(dy[[v]])
      data.frame(original = v, dummy = paste(v, lev, sep = ""), level = lev)
    })
    dlink <- do.call(rbind, dlink)

    # One-hot encode to sparse matrix (requires Matrix package for speed; faster than pure data.table implementation)
    # https://stackoverflow.com/questions/4560459/all-levels-of-a-factor-in-a-model-matrix-in-r
    dy <- Matrix::sparse.model.matrix(~ 0 + .,
                                      data = dy,
                                      contrasts.arg = lapply(dy, contrasts, contrasts = FALSE),
                                      drop.unused.levels = FALSE)
                                      #drop.unused.levels = dropUnusedLevels)

    # Could return sparse matrix, but appending to original data frame by default
    dy <- as.matrix(dy)

    # Drop dummy columns with no 1's
    if (dropUnusedLevels) {
      keep <- colSums(dy != 0) > 0
      dy <- dy[, keep]
    }

    data <- cbind(data, dy)
    attr(data, "one_hot_link") <- dlink
  }
  return(data)

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

#------------------

# Function to check a 'data' object against inputs 'y' and 'x'
# Detects character columns, no-variance columns, and missing data columns (imputes as needed)
# This simply wraps a code chunk that was present in train(), blockchain(), and prescreen()

checkData <- function(data, y, x) {

  # Check for character-type variables; stop with error if any detected
  xc <- sapply(data[c(x, y)], is.character)
  if (any(xc)) stop("Coerce character variables to factor:\n", paste(names(which(xc)), collapse = ", "))

  # Check for no-variance (constant) variables
  # Stop with error if any 'y' are constant; remove constant 'x' with message
  constant <- names(which(sapply(data[y], novary)))
  if (length(constant)) stop("Zero-variance 'y' variable(s) detected (remove them):\n", paste(constant, collapse = ", "))
  constant <- names(which(sapply(data[x], novary)))
  if (length(constant)) {
    x <- setdiff(x, constant)
    data <- select(data, -all_of(constant))
    cat("Removed zero-variance 'x' variable(s):\n", paste(constant, collapse = ", "), "\n")
  }

  # Detect and impute any missing values in 'x' variables
  na.cols <- names(which(sapply(data[x], anyNA)))
  if (length(na.cols) > 0) {
    cat("Missing values imputed for the following 'x' variable(s):\n", paste(na.cols, collapse = ", "), "\n")
    for (j in na.cols) {
      xj <- data[[j]]
      ind <- is.na(xj)
      data[ind, j] <-  imputationValue(xj, ind)
    }
  }

  return(data)

}

#------------------

# Function to return in-memory size of object in Mb
objectSize <- function(x) as.numeric(utils::object.size(x)) / 1048576

#------------------

# Function to return free (actually, "available") system memory in Mb
# Different commands used on Linux vs. MacOS vs. Windows
# https://stackoverflow.com/questions/27788968/how-would-one-check-the-system-memory-available-using-r-on-a-windows-machine
# https://stackoverflow.com/questions/6457290/how-to-check-the-amount-of-ram
# On free vs available: https://haydenjames.io/free-vs-available-memory-in-linux/
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
      # Linux system assumed as backstop
      # See output of 'cat /proc/meminfo'
      x <- system('grep MemAvailable /proc/meminfo', intern = TRUE)
      x <- strsplit(x, "\\s+")[[1]][2]
      as.numeric(x) / 1024
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
#   # NOT USED: Initial brute downsample step, if necessary
#   # Could be useful if the inputs are really big, making kmeans too slow
#   # s0 <- if (N_init < N) sample.int(N, size = N_init) else 1:N
#   # x <- x[s0]
#   # y <- y[s0]
#   # if (!is.null(w)) w <- w[s0] * (sum(w) / sum(w[s0]))  # Preserve the original total sample weight
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
#   # Measure of the unusual-ness of each cluster (relative distance from all other cluster centers)
#   # Used to calculate 'cadj'; cluster-specific downsampling factor
#
#   # Center of the (scaled) dataset
#   #cnt <- c(weighted.mean(x, w), weighted.mean(y, w))  # Only necessary b/c scale() above is not weighted
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
#     #k <- s0[k]
#     wout <- w[k] * (sum(w[i]) / sum(w[k]))
#     cbind(x = x0[k], y = y0[k], w = wout)
#   })
#
#   samp <- do.call(rbind, samp)
#   return(as.data.frame(samp))
#
# }

#------------------

# Function to cluster a univariate continuous variable into 'k' clusters
# Returns an ordered factor with k ordered levels
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
  while (!halt) {
    iter <- iter + 1
    km <- suppressWarnings(kmeans(x, centers = k, nstart = 3 * iter, iter.max = 20))
    halt <- km$ifault == 0
  }
  cl <- data.table::frank(km$centers, ties.method = "dense")[km$cluster]
  cl <- factor(cl, levels = 1:k, ordered = TRUE)
  return(cl)
}

#------------------

# This function returns weighted mean and approximate variance about the weighted mean
# 'x' must be numeric but can be continuous or binary
# In binary case, it return the (weighted) proportion and standard error
# Note the difference in standard error calculation for continuous vs. binary case
weightedVAR <- function(x, w) {

  # Continuous case
  # Code: https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
  # Computes the variance of a weighted mean following Cochran 1977: https://www.sciencedirect.com/science/article/abs/pii/135223109400210C
  if (any(!x %in% 0:1)) {
    n <- length(x)
    wbar <- mean(w)
    xWbar <- wmean(x, w)
    wvar <- n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
  } else {

    # Binary case
    # https://stats.stackexchange.com/questions/159204/how-to-calculate-the-standard-error-of-a-proportion-using-weighted-data
    w <- w / sum(w)
    sw2 <- sum(w ^ 2)
    xWbar <- sum(w * x)
    wvar <- xWbar * (1 - xWbar) * sw2
  }
  return(c(xWbar, wvar))
}

#------------------

# Function to quickly estimate the conditional quantile smoother for a bivariate (x, y) scatterplot
# The results are not deterministic, but they are pretty stable if 'n' is higher enough
# The returned curve is a LOESS fit to k * n observations of x and median(y)|x derived via repeated kmeans clustering of the 'x' values
# The LOESS fit weights cluster observations by sqrt(cluster size); this down-weights small-sample medians

# Test example from ?qgam
# library(MASS)
# x <- mcycle$times
# y <- mcycle$accel
# xout <- sort(unique(x))
# fit.qgam <- qgam::qgam(accel~s(times), data = mcycle, qu = 0.5)
# pred.qgam <- predict(fit.qgam, data.frame(times = xout))
# plot(x, y)
# lines(xout, pred.qgam, col = 2)
#
# test <- smoothQuantile(x, y)
# lines(test$x, test$y, col = 3)
#
# data(recs)
# x <- recs$square_feet
# y <- recs$electricity
# xout <- xout0 <-  unique(sort(x)[seq(from = 1, to = length(x), length.out = min(length(x), 200))])
# fit.qgam <- qgam::qgam(y~s(x), data = data.frame(x, y), qu = 0.5)
# pred.qgam <- predict(fit.qgam, data.frame(x = xout))
# plot(x, y)
# lines(xout0, pred.qgam, col = 2)
#
# test <- smoothQuantile(x, y)
# lines(test$x, test$y, col = 3)

smoothQuantile <- function(x,
                           y,
                           qu = 0.5,
                           xout = NULL,
                           a = 0.3,
                           b = 0.1) {

  stopifnot({
    !anyNA(x)
    !anyNA(x)
    length(x) == length(y)
    qu == "mad" | all(qu > 1, qu < 1)
    is.null(xout) | length(xout) > 1
    a > 0
    b > 0
  })

  qmad <- qu[1] == "mad"
  qu <- if (qmad) c(0.25, 0.5, 0.75) else sort(unique(qu))
  n <- length(x)
  i <- order(x)
  x <- x[i]
  y <- y[i]
  d <- data.table(x, y)
  mb <- max(5, ceiling(n ^ a))

  # Check if 'xout' is outside natural range of 'x'
  if (!is.null(xout)) {
    if (any(xout < min(x) | xout > max(x))) warning("'xout' contains values outside range of 'x' (be careful)")
  }

  #---

  # We can initially sample 'x' for the rpart() call
  #  since we only care about estimating plausible bin breakpoints
  f <- min(1, max(0.05, 30e3 / n))
  df <- if (f < 1) {
    ind <- round(seq(from = 1, to = n, by = (n - 1) / (round(f * n) - 1)))
    ind[length(ind)] <- n
    d[ind, ]
  } else {
    d
  }

  #---

  # Fit rpart to determine bin breakpoints along 'x'
  mbr <- max(5, ceiling(mb * f))
  fit <- rpart::rpart(formula = y ~ x, data = df,
                      xval = 0, cp = 0.001,
                      maxcompete = 0, maxsurrogate = 0,
                      minbucket = mbr, minsplit = 2 * mbr)

  r <- rle(as.integer(as.factor(fit$where)))
  rsum <- c(1, cumsum(r$lengths))

  # Initial check if there is sufficient number of bins (require at least 10)
  #k <- 3 * (length(rsum) - 1)
  #if (k < 10) stop("Too few bins (", k, "); try using smaller 'a' value\n", sep = "")

  # Check plot -- useful diagnostic
  # fit.mean <- mgcv::bam(y ~ s(x), data = df, discrete = TRUE)
  # temp <- df[sample.int(n = nrow(df), size = min(n, 10e3), replace = TRUE), ]
  # plot(temp)
  # lines(df$x, fitted(fit.mean), col = 2)
  # abline(v = df$x[rsum], col = 3)

  #---

  # Get the bin assignments and compute the bin metrics
  # Uses data.table operations for speed

  # This creates 3 binning vectors in 'd'
  # Each row/observation is assigned to 3 different bins/windows
  for (i in 1:3) {
    p <- c(0, 0.5, -0.5)[i]
    br <- pmax(1, round(p * diff(rsum) + rsum[-(length(rsum))]))
    br <- c(br, max(rsum))
    if (p > 0) br <- c(min(rsum), br)
    bin <- cut(x, breaks = unique(df$x[br]), labels = FALSE, include.lowest = TRUE)
    set(d, j = paste0("bin", i), value = bin)
  }

  # This computes the bin-specific outcomes (count, x-mean, and y-quantiles)
  d <- lapply(1:3, function(i) {
    d[, list(s = .N, m = sum(x), q = quantile(y, probs = qu)), by = eval(paste0("bin", i))]
  })
  d <- rbindlist(d, use.names = FALSE, idcol = "window")
  d[, m := m / s]  # Mean x-value of each bin
  d[, qu := rep(qu, times = nrow(d) / length(qu))]

  # Retain only those bins with sufficient number of observations
  d <- d[s >= mb, ]

  # Final check if there is sufficient number of bins (require at least 10)
  #k <- nrow(d) / length(qu)
  #if (k < 10) stop("Too few bins (", k, "); try using smaller 'a' value\n", sep = "")

  # If qmad, replace 25th and 75th percentiles with each bin's observed lower and upper MAD, respectively
  # This ensures the GAM smooth output reflects the lower and upper MAD rather than the quantiles themselves
  if (qmad) {
    f <- function(x) c(x[2] - x[1], x[2], x[3] - x[2])
    d[, q := f(q), by = c("window", "bin1")]
  }

  # Order bins by mean x-value
  #setorder(d, m)

  # Force the minimum and maximum x-values into the bin results
  # This ensures that the GAM fit attempts to smooth out to the extreme x-values
  #if (!qmad) {
  d[m == min(m), m := min(x)]
  d[m == max(m), m := max(x)]
  #}

  #---

  # x-values at which to predict smoothed values
  mout <- if (is.null(xout) | length(xout) == 2) {
    rng <- if (is.null(xout)) range(d$m) else range(xout)
    xx <- if (is.null(xout)) x else x[x >= rng[1] & x <= rng[2]]
    mn <- diff(rng) / (b * mad(xx, constant = 1))
    mn <- min(mn, 200)  # Hard max on number of output x-values
    seq(from = rng[1], to = rng[2], length.out = ceiling(mn))
  } else {
    xout
  }

  #---

  # Fit GAM(s) and return the smoothed quantile values for each 'mout' x-value
  out <- sapply(qu, function(i) {
    df <- d[qu == i, ]
    fit <- mgcv::gam(formula = q ~ s(m), data = df, weights = sqrt(df$s))  # Weighted by square root of bin sample size
    predict(fit, newdata = data.table(m = mout))
  })

  # Ensure there is no quantile crossing in the results
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

smoothMean <- function(x, y, xout = NULL,
                       zmax = 3.5, a = 0.3, b = 0.1,
                       se = FALSE, verbose = TRUE) {

  n <- length(x)
  i <- order(x)
  x <- x[i]
  y <- y[i]

  #---

  # Eliminate obviously extreme y values?
  zy <- abs(y - median(y)) / mad(y)
  iy <- zy > zmax ^ 2

  # Eliminate obviously extreme x values?
  zx <- abs(x - median(x)) / mad(x)
  ix <- zx > zmax ^ 2

  # Initial outliers
  oi <- iy | ix
  x <- x[!oi]
  y <- y[!oi]

  #---

  # If 'y' is already monotonic or has no variance, simply return a linear interpolation
  simple <- sorted(y) | var(y) == 0
  xout <- if (simple & is.null(xout)) range(x) else xout

  # Data table used below
  df <- data.table(x, y)

  #---

  # Detect and remove outliers
  if (!simple) {

    # Return smoothed median and lower and upper MAD
    d <- smoothQuantile(x, y, qu = "mad", xout = NULL, a = a, b = b)
    d <- split(d, d$q)

    # Smoothed median for each x value
    m <- approx(x = d[[2]]$x, y = d[[2]]$y, xout = x)$y

    # Upper MAD
    i <- d[[3]]$y > 0
    d1 <- approx(x = d[[3]]$x[i], y = d[[3]]$y[i], xout = x)$y
    d1[d1 == 0] <- NA

    # Lower MAD
    i <- d[[1]]$y > 0
    d2 <- approx(x = d[[1]]$x[i], y = d[[1]]$y[i], xout = x)$y
    d2[d2 == 0] <- NA

    # See the median and the cutoff lines for upper and lower outliers
    # Points falling outside the envelope are considered outliers
    # plot(x, y)
    # lines(x, m, col = 2)
    # lines(x, m + zmax * d1, col = 3)
    # lines(x, m - zmax * d2, col = 4)

    # Detect the y outliers
    out1 <- y > m + zmax * d1
    out2 <- y < m - zmax * d2
    oy <- out1 | out2

    # Report total number of outliers
    nout <- sum(oi) + sum(oy)
    #nout <- sum(oy)
    if (verbose) cat("Removed ", nout, " outliers (", paste0(signif(100 * nout / n, 3), "%"), ")\n", sep = "")

    # Remove the outlier observations
    df <- df[!oy]

  }

  #---

  # x-values at which to predict smoothed values
  xout <- if (is.null(xout)) {
    d[[1]]$x  # Use the x-values returned by quantile smoother
  } else {
    if (length(xout) == 2) {
      xx <- x[x >= min(xout) & x <= max(xout)]
      rng <- range(xx)
      mn <- diff(rng) / (b * mad(xx, constant = 1))
      mn <- min(mn, 200)  # Hard max on number of output x-values
      seq(from = rng[1], to = rng[2], length.out = ceiling(mn))
    } else {
      unique(sort(xout))
    }
  }

  #---

  # Final predicted smooth output
  result <- if (simple) {

    out <- data.frame(x = xout,
                      y = suppressWarnings(approx(x = df$x, y = df$y, xout = xout)$y),
                      se = 0)
    if (!se) out$se <- NULL
    out

  } else {

    # Fit smoothed mean (gam or bam)
    big <- nrow(df) > 30e3
    fun <- ifelse(big, mgcv::bam,  mgcv::gam)
    if (verbose) cat("Fitting smoothed mean using", ifelse(big, "mgcv::bam()", "mgcv::gam()"), "\n")
    fit <- fun(y ~ s(x), data = df, discrete = big)
    data.frame(x = xout, predict(fit, data.table(x = xout), se.fit = se))

  }

  # Clean up and return result
  names(result) <- c("x", "y", if (se) "se" else NULL)
  row.names(result) <- NULL
  return(result)

}

#------------------

# Create 'maxr_fun' function used inside validate() to derive the cut-off point for ratio 'r' (ubar / b)
# CI's using 'r' larger than that returned by maxr_fun() are invalid
# The cutoff point increases with the number of implicates

mseq <- seq(2, 100, by = 2)
out <- sapply(mseq, function(m, p = 0.95) {

  # Ratio of ubar to b
  R <- seq(0, 1, length.out = 1e4)

  # Degrees of freedom (checked; this is correct)
  df <- (m - 1) * (1 - R * m  / (m + 1)) ^ 2

  # Calculate SE and associated CI
  ubar <- 1  # Doesn't matter; changes magnitude but shape of subsequent curves
  se <- sqrt(ubar / R * (1 + 1 / m) - ubar)  # Estimated standard error
  ci <- se * qt(p, df)  # Estimated CI half-width

  # For a given ubar, the CI width declines with R up to some point and then begins to increase
  # We expect CI to decline with R, since increasing R means b is getting smaller relative to ubar (greater confidence)
  # But the eventual increase doesn't make sense, so we return the 'R' value at which the CI is minimized
  #plot(R, ci, type = "l")

  # The ratio (ubar / b) at which CI width is minimized
  # Ratios beyond this point return rapidily increasing CI widths
  R[which.min(ci)]

})

# plot(mseq, out, type = "l")
maxr_fun <- approxfun(mseq, out, rule = 2)
rm(mseq, out)

#------------------

# Function to check if values are sorted; is.unsorted() checks only for increasing
sorted <- function(x) {
  xd <- diff(x)
  all(xd >= -sqrt(.Machine$double.eps)) | all(xd <= sqrt(.Machine$double.eps))
}



