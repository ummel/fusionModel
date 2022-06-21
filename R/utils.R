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

# Create nice ggplot2 barplot of variable importance in rpart() model
# varImp <- function(m) {
#   imp <- m$variable.importance / sum(m$variable.importance)
#   enframe(imp) %>%
#     filter(value > 0) %>%
#     ggplot(aes(x = reorder(name, value), y = value)) +
#     geom_bar(stat = "identity") +
#     coord_flip() +
#     theme_bw(11)
# }

#------------------

# Remove elements of rpart() model object to reduce size while still allowing prediction
slimRpart <- function(m) {
  m$call <- NULL
  m$numresp <- NULL
  m$parms <- NULL
  m$functions <- NULL
  m$ordered <- NULL
  m$y <- NULL
  #attr(m, "xlevels") <- NULL
  if ("splits" %in% names(m)) m$splits[, c('improve', 'adj')] <- 0
  if ("frame" %in% names(m)) m$frame[, c('n', 'wt', 'dev', 'complexity')] <- 0
  if ("terms" %in% names(m)) attributes(m$terms)['.Environment'] <- NULL
  return(m)
}

#------------------

# Remove elements of biglm() model object to reduce size while still allowing prediction
# slimBiglm <- function(m) {
#   m$call <- NULL
#   m$assign <- NULL
#   m$df.resid <- NULL
#   m$weights <- NULL
#   m$n <- NULL
#   m$names <- NULL
#   return(m)
# }

#------------------

# Validate input variable vector against a given set of names
# 'x' can include regular expressions
validNames <- function(x, nms, exclude = FALSE) {
  rgx <- setdiff(x, nms)
  v <- c(intersect(x, nms), unlist(lapply(rgx, function(x) grep(x, nms, value = TRUE))))
  if (exclude) v <- setdiff(nms, v)
  out <- v[na.omit(match(nms, v))]
  return(out)
}

#------------------

# Return rpart() model prediction node for each observation in 'newdata'
# Based on: https://github.com/cran/rpart.plot/blob/master/R/rpart.predict.R
predictNode <- function(object, newdata) {

  where <-
    if (missing(newdata)) {
      object$where
    } else {
      if(is.null(attr(newdata, "terms"))) {
        Terms <- delete.response(object$terms)
        newdata <- model.frame(Terms, newdata, na.action = na.pass,
                               xlev = attr(object, "xlevels"))
        if(!is.null(cl <- attr(Terms, "dataClasses")))
          .checkMFClasses(cl, newdata, TRUE)
      }
      newdata <- getFromNamespace("rpart.matrix", ns="rpart")(newdata)
      getFromNamespace("pred.rpart", ns="rpart")(object, newdata)
    }

  # nn <- as.numeric(rownames(object$frame)[where])
  #
  # # Index assigning each row in 'pred' to a node number
  # node <- match(nn, as.integer(row.names(object$frame)))

  #return(node)
  return(where)

}

#------------------

# TO DO: REPLACE use in train() with faster application
# Function to create appropriate matrix object (perhaps dummy variables) for variable 'v' in data frame 'data'
# This used within both train() and fuse()
# matFun <- function(v, data) {
#   stopifnot(length(v) == 1L)
#   x <- data[[v]]
#   if (is.numeric(x) | is.ordered(x)) {
#     m <- data.table::frank(x)
#     dim(m) <- c(length(m), 1L)  # See "Note" in ?matrix about converting vector to matrix
#     dimnames(m) <- list(NULL, v)
#   } else {
#     u <- levels(x)
#     m <- matrix(data = 0L, nrow = length(x), ncol = length(u), dimnames = list(NULL, paste0(v, u)))
#     for (i in 1:ncol(m)) m[x == u[i], i] <- 1L
#   }
#   return(m)
# }

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
wmean <- function(x, w) sum(w * x) / sum(w)

#------------------

# Weighted standard deviation
# Equivalent to Hmisc::wtd.var() with normwt = TRUE and taking sqrt() of result
# wsd <- function(x, w) {
#   w <- w * length(x) / sum(w)
#   sw <- sum(w)
#   xbar <- sum(w * x) / sw
#   sqrt(sum(w * ((x - xbar) ^ 2)) / (sw - 1))
# }

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

#------------------

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
