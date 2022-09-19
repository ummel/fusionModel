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

#------------------

# One-hot encoding of data.table for use in fitting LASSO models
# Based on: https://github.com/ben519/mltools/blob/master/R/one_hot.R
# Note that some arguments have been dropped and their original default values are assumed
# Default assumptions:
# dropCols = TRUE -- Remove original columns
# dropUnusedLevels = TRUE -- Remove columns of all zeros
# naCols = FALSE

one_hot <- function(dt, sparse_matrix = TRUE, sparsifyNAs = TRUE, ord_to_int = FALSE, sep = "..") {

  if (!is.data.table(dt)) dt <- as.data.table(dt)

  # Convert ordered factors and logicals to integer
  if (ord_to_int) dt <- mutate_if(dt, is.ordered, as.integer)
  dt <- mutate_if(dt, is.logical, as.integer)

  OHEID <- NULL
  cols <- colnames(dt)[which(sapply(dt, function(x) is.factor(x)))]

  out <- if (length(cols) > 0) {
    tempDT <- dt[, cols, with = FALSE]
    tempDT[, `:=`(OHEID, .I)]
    for (col in cols) set(tempDT, j = col, value = factor(paste(col, tempDT[[col]], sep = sep), levels = paste(col, levels(tempDT[[col]]), sep = sep)))
    melted <- data.table::melt(tempDT, id = "OHEID", value.factor = T, na.rm = TRUE)
    newCols <- data.table::dcast(melted, OHEID ~ value, drop = T, fun.aggregate = length)
    newCols <- newCols[tempDT[, list(OHEID)]]
    newCols[is.na(newCols[[2]]), `:=`(setdiff(paste(colnames(newCols)), "OHEID"), 0L)]
    if (!sparsifyNAs) {
      na_cols <- character(0)
      for (col in cols) if (any(is.na(tempDT[[col]]))) na_cols <- c(na_cols, col)
      if (!sparsifyNAs) for (col in na_cols) newCols[is.na(tempDT[[col]]), `:=`(intersect(levels(tempDT[[col]]), colnames(newCols)), NA_integer_)]
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
    result[, !cols, with = FALSE]
  } else{
    dt
  }
  if (sparse_matrix) out <- Matrix::Matrix(as.matrix(out), sparse = TRUE)
  return(out)
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

# Structured downsampling of 2-dimensional (x, y) inputs
# Compared to random sampling, guarantees a more uniform sampling of the space and oversampling of unusual observations
# Returns adjusted weights such that total cluster weights are respected in the downsample
# Will return no less than N observations, but can return significantly more due to oversampling of small/unusual clusters

# Example usage
# library(fusionModel)
# x <- recs$electricity
# y <- recs$square_feet
# test <- downsample(x = recs$electricity, y = recs$square_feet)
# dim(test)
# plot(recs$electricity, recs$square_feet)
# points(test[1:2], col = 2)
# summary(test$w)

downsample <- function(x,
                       y,
                       w = rep(1, length(x)),
                       N = 0.1 * length(x),
                       K = 30,
                       min_samp = 30) {

  stopifnot({
    length(x) == length(y)
    length(x) == length(w)
    !anyNA(c(x, y, w))
    N >= 1
    K >= 1
    min_samp >= 1
  })

  N0 <- length(x)

  # NOT USED: Initial brute downsample step, if necessary
  # Could be useful if the inputs are really big, making kmeans too slow
  # s0 <- if (N_init < N) sample.int(N, size = N_init) else 1:N
  # x <- x[s0]
  # y <- y[s0]
  # if (!is.null(w)) w <- w[s0] * (sum(w) / sum(w[s0]))  # Preserve the original total sample weight

  # Keep copies of inputs
  x0 <- x
  y0 <- y

  # Scale x and y prior to clustering
  x <- scale(x)
  y <- scale(y)

  # Average required adjustment factor
  adj <- N / N0

  # k-means
  km <- kmeans(cbind(x, y), centers = K, iter.max = K * 5)

  # Measure of the unusual-ness of each cluster (relative distance from all other cluster centers)
  # Used to calculate 'cadj'; cluster-specific downsampling factor

  # Center of the (scaled) dataset
  #cnt <- c(weighted.mean(x, w), weighted.mean(y, w))  # Only necessary b/c scale() above is not weighted
  cnt <- apply(km$centers, 2, weighted.mean, w = km$size)

  # Distance of each cluster center from the dataset center
  dc <- sqrt((km$centers[, 1] - cnt[1]) ^ 2 + (km$centers[, 2] - cnt[2]) ^ 2)
  ofct <- dc / weighted.mean(dc, km$size)
  cadj <- adj * ofct  # Cluster-specific adjustment factor

  # Clustered downsample
  samp <- sapply(1:nrow(km$centers), function(cl) {
    i <- which(km$cluster == cl)
    cnt <- km$centers[cl, ]  # Center of the cluster
    d <- sqrt((x[i] - cnt[1]) ^ 2 + (y[i] - cnt[2]) ^ 2)
    j <- seq(from = 1, to = length(d), length.out = max(min_samp, length(d) * cadj[cl]))
    j <- unique(round(j))
    k <- i[order(d)][j]
    #k <- s0[k]
    wout <- w[k] * (sum(w[i]) / sum(w[k]))
    cbind(x = x0[k], y = y0[k], w = wout)
  })

  samp <- do.call(rbind, samp)
  return(as.data.frame(samp))

}
