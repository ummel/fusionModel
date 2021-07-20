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

# Remove elements of rpart() mode object to reduce size while still allowing prediction
#sapply(m, object.size)
slimRpart <- function(m) {
  m$where <- NULL
  m$call <- NULL
  #m$cptable <- NULL
  m$numresp <- NULL
  m$parms <- NULL
  m$functions <- NULL
  m$ordered <- NULL
  m$y <- NULL
  attr(m, "xlevels") <- NULL
  if ("splits" %in% names(m)) m$splits[, c('improve', 'adj')] <- 0
  if ("frame" %in% names(m)) m$frame[, c('n', 'wt', 'dev', 'complexity')] <- 0
  if ("terms" %in% names(m)) attributes(m$terms)['.Environment'] <- NULL
  return(m)
}

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

  nn <- as.numeric(rownames(object$frame)[where])

  # Index assigning each row in 'pred' to a node number
  node <- match(nn, as.integer(row.names(object$frame)))

  return(node)

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
    if (is.integer(x)) as.integer(round(m)) else m
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
