#' Train a fusion model on donor data
#'
#' @description
#' Train a fusion model on donor data.
#'
#' @param data Data frame. Donor dataset. All categorical variables should be factors and ordered whenever possible.
#' @param y Character vector. The variables to fused to recipient dataset. Can include \link[base:regex]{regular expressions}.
#' @param x Character vector. Predictor variables common to donor and recipient. Can include regular expressions. If NULL, all variables other than those in \code{y} and \code{weight} are used. Only one of \code{x} and \code{ignore} can be non-NULL.
#' @param ignore Character vector. Alternative way to specify predictor variables. Can include regular expressions. If non-NULL, all variables other than those in \code{y}, \code{weight}, and \code{ignore} are used. Only one of \code{x} and \code{ignore} can be non-NULL.
#' @param weight Character vector. Name of the observation weights column. If NULL (default), uniform weights are assumed.
#' @param mc Logical. Should multicore processing be used? Implementation uses \code{\link[parallel]{mclapply}}, so \code{mc = TRUE} will fail on Windows.
#' @param maxcats Positive integer or NULL. Maximum number of levels allowed in an unordered factor predictor variable when the response (fusion) variable is also an unordered factor. Prevents excessive \code{\link[rpart]{rpart}} computation time. A K-means clustering strategy is used to cluster the predictor to no more than \code{maxcats} levels.
#' @param lasso Numeric (0-1) or NULL. Controls extent of (fast) predictor variable pre-screening via LASSO regression. If NULL (default), no pre-screening is performed. See Details.
#' @param ... Optional arguments passed to \code{\link[rpart]{rpart}} to control tree-building. By default \code{cp = 0}, \code{xval = 0}, \code{minbucket = 50} (\code{minbucket = 10} for discrete models), and all other arguments are left at default values.
#'
#' @details When \code{lasso} is non-NULL, predictor variables are "pre-screened" using LASSO regression via \code{\link[glmnet]{glmnet}} prior to fitting a \code{\link[rpart]{rpart}} model. Predictors with a LASSO coefficient of zero are excluded from consideration. This can speed up tree-fitting considerably when \code{data} is large. Lower values of \code{lasso} are more aggressive at excluding predictors; the LASSO \emph{lambda} is chosen such that the deviance explained is at least \code{lasso}-% of the maximum. To ensure the LASSO step is fast, pre-screening is only used for numeric and ordered factor response variables. No predictor pre-screening occurs when the response is an unordered factor.
#'
#' @return A list containing trained model information to be passed to \link{fuse}.
#'
#' @examples
#' donor <- recs
#' recipient <- subset(recs, select = c(division, urban_rural, climate, income, age, race))
#' fusion.vars <- setdiff(names(donor), names(recipient))
#' fit <- train(data = donor, y = fusion.vars)
#' @export

#---------------------

# Manual testing

# library(fusionModel)
# source("R/utils.R")
# source("R/fitRpart.R")
# source("R/fitDensity.R")
# source("R/fusionOrder.R")
#
# data <- recs
# #recipient <- subset(recs, select = c(division, urban_rural, climate, income, age, race, education, employment))
# recipient <- subset(recs, select = c(division, urban_rural, climate, income))
# y = setdiff(names(data), names(recipient))
# weight = NULL
# x <- NULL
# ignore = NULL
# maxcats = 10
# mc = TRUE
# lasso = 0.9

#---------------------

train <- function(data,
                  y,
                  x = NULL,
                  ignore = NULL,
                  weight = NULL,
                  mc = FALSE,
                  maxcats = NULL,
                  lasso = NULL,
                  ...) {

  stopifnot(exprs = {
    is.data.frame(data)
    !missing(y)
    is.logical(mc)
    is.null(maxcats) | (maxcats > 1 & maxcats %% 1 == 0)
    is.null(lasso) | (lasso >= 0 & lasso <= 1)
  })

  # Check that no more than one of 'x' or 'ignore' is specified
  if (!is.null(x) & !is.null(ignore)) stop("Only one of 'x' or 'ignore' may be non-NULL")

  # Set default rpart() arguments and modify with user-supplied arguments, if necessary
  rpart.args <- list(cp = 0, xval = 0, minbucket = 50, maxcompete = 0)  # Default rpart() arguments
  rpart.default <- TRUE
  if (!missing(...)) {
    rpart.default <- FALSE
    allowed <- c(names(formals(rpart::rpart)), names(formals(rpart::rpart.control)))
    allowed <- setdiff(allowed, c("formula", "data", "weights", "subset", "na.action", "method"))  # Allowable arguments for '...'
    usr.args <- list(...)
    miss <- setdiff(names(usr.args), allowed)
    if (length(miss) > 0) stop("Improper arguments passed to '...':\n", paste(miss, collapse = ", "))
    rpart.args[names(usr.args)] <- usr.args
  }

  # Number of cores to use in parallel
  # If NULL, defaults to lapply() behavior via pblapply()
  mc.cores <- if (mc) max(1L, parallel::detectCores() - 1L) else NULL

  #-----

  # Detect fusion variables
  nms <- names(data)
  yvars <- validNames(y, nms, exclude = FALSE)

  # Detect predictor variables
  xvars <- setdiff(nms, c(yvars, weight))
  if (!is.null(ignore)) xvars <- validNames(ignore, xvars, exclude = TRUE)

  # Check validity of variables
  stopifnot(exprs = {
    length(yvars) > 0
    length(xvars) > 0
    !anyNA(data[yvars])  # Require non-NA values in response values
  })

  #-----

  # Create and/or check observation weights
  # Weights are normalized (mean = 1) to avoid integer overflow issues
  w <- "_(wgt)_"
  if (is.null(weight)) {
    data[[w]] <- rep(1L, nrow(data))
  } else {
    stopifnot(exprs = {
      weight %in% names(data)
      !weight %in% names(c(xvars, yvars))
      !anyNA(data[[weight]])
      all(data[[weight]] >= 0)
    })
    data[[w]] <- data[[weight]] / mean(data[[weight]])  # Scaled weight to avoid numerical issues with larger integer weights
    data[[weight]] <- NULL
  }

  #-----

  # Check for no-variance variables; stop with error if any detected
  x <- sapply(data[c(xvars, yvars)], novary)
  if (any(x)) stop("The following variables have no variance:\n", paste(names(which(x)), collapse = ", "))

  # Check for character-type variables; stop with error if any detected
  x <- sapply(data[c(xvars, yvars)], is.character)
  if (any(x)) stop("Please coerce character variables to factor:\n", paste(names(which(x)), collapse = ", "))

  #-----

  # Check that the 'xvars' and 'yvars' contain only syntactically valid names
  stopifnot(all(make.names(yvars, unique = TRUE) == yvars))
  stopifnot(all(make.names(xvars, unique = TRUE) == xvars))

  # Print to console
  cat("Identified ", length(yvars), " fusion variables\n")
  cat("Identified ", length(xvars), " predictor variables\n")

  # Limit 'data' to the necessary variables
  data <- data[c(w, xvars, yvars)]

  #-----

  # Identify which of the 'yvars' are continuous
  ycont <- names(which(sapply(data[yvars], is.numeric)))

  # Observed "outer" range" (min/max) for continuous 'yvars'
  youter <- lapply(data[ycont], range, na.rm = TRUE)

  # The "inner range" for continuous 'yvars'; i.e. largest negative and smallest positive value
  # If there are no negative or positive values, 0 is returned
  yinner <- lapply(data[ycont], function(x) c(ifelse(any(x < 0), max(x[x < 0], na.rm = TRUE), 0), ifelse(any(x > 0), min(x[x > 0], na.rm = TRUE), 0)))

  # Extract data classes and levels for the 'yvars'
  x <- data[yvars]
  yclass <- lapply(x, class)
  ylevels <- lapply(x[grepl("factor", yclass)], levels)

  # Extract data classes and levels for the 'xvars'
  x <- data[xvars]
  xclass <- lapply(x, class)
  xlevels <- lapply(x[grepl("factor", xclass)], levels)

  #----------------------------------
  #----------------------------------

  # DETERMINE FUSION ORDER

  cat("Determining order of fusion variables...\n")

  # Fit full models, including yvars as predictors
  # Only 'variable.importance' slot is returned
  full.varimp <- pbapply::pblapply(X = yvars, FUN = function(y) {
    m <- fitRpart(y = y,
                  x = c(xvars, setdiff(yvars, y)),
                  w = w,
                  data = data,
                  maxcats = maxcats,
                  args = if (rpart.default & !y %in% ycont) replace(rpart.args, 3, 10) else rpart.args, # Sets default minbucket = 10 in the discrete case
                  lasso.threshold = lasso)
    gc()  # Perhaps useful when run in parallel
    m$variable.importance
  }, cl = mc.cores) %>%
    setNames(yvars)

  # Find preferred order of yvars
  yord <- fusionOrder(varimp = full.varimp)

  # Cleanup
  rm(full.varimp)

  #----------------------------------
  #----------------------------------

  # FIT MODELS

  # Assemble dummy matrix for 'xvars' and 'yvars'
  ranks <- lapply(c(xvars, yord), matFun, data = data)
  ranks.var <- unlist(mapply(rep, x = c(xvars, yord), each = sapply(ranks, ncol)), use.names = FALSE)  # Variable associated with each column of 'ranks' matrix
  ranks <- do.call(cbind, ranks)

  # Placeholder list for 'ycor' with slots for continuous and ordered factor fusion variables (but not unordered factors)
  ycor.vars <- names(which(sapply(yclass, function(x) x[1] != "factor")))

  #-----

  cat("Building fusion models...\n")

  # Function to fit rpart() model and calculate rank correlations for i-th response variable in 'yord'
  fitModel <- function(i) {

    y <- yord[i]

    # Is 'y' continuous?
    cont <- y %in% ycont

    # Prior response variables in the fusion sequence
    yprior <- if (i == 1) NULL else yord[1:(i - 1)]

    #--------

    # Row index giving location of non-NA response values
    ind <- which(!is.na(data[[y]]))
    d <- data[ind, ]

    # Fit model
    m <- fitRpart(y = y,
                  x = c(xvars, yprior),
                  w = w,
                  data = d,
                  maxcats = maxcats,
                  args = if (rpart.default & !cont) replace(rpart.args, 3, 10) else rpart.args, # Sets default minbucket = 10 in the discrete case
                  lasso.threshold = lasso)

    #--------

    # If 'y' is continuous...
    if (cont) {

      # Create list in 'm' to hold the quantile values associated with a set of 'N' known percentiles; see fitDensity()
      # The 'Q' values describe the shape of the quantile function (i.e. inverse CDF)
      nodes <- sort(unique(m$where))

      # Placeholder matrix for the smoothed quantile values
      Qn <- 500  # Ideally, would be a train() control argument
      m$Q <- matrix(0L, nrow = Qn, ncol = length(nodes))
      colnames(m$Q) <- nodes

      # Fit density to observations in each node
      for (n in nodes) {

        # Index identifying observations in node 'n'
        nind <- m$where == n

        # Node density result
        fd <- fitDensity(x = d[nind, y], w = d[nind, w], inner.range = yinner[[y]], outer.range = youter[[y]], N = Qn)

        # Add quantile function values to 'm'
        m$Q[, as.character(n)] <- fd

      }

    }

    #--------

    # Calculate (rank) correlation between 'y' and existing variables in 'ranks' (if necessary; only continuous and ordered factors)
    # NOTE: Correlation calculation restricted to observations where 'y' is non-zero
    if (y %in% ycor.vars) {
      ind <- which(as.numeric(data[[y]]) > 0)
      j <- which(ranks.var == y)  # Column in 'ranks' associated with 'y'
      p <- suppressWarnings(cor(ranks[ind, 1:(j - 1)], ranks[ind, j]))[, 1]  # Will silently return NA if no variance in a particular column
      p <- cleanNumeric(p[!is.na(p)], tol = 0.001)
    } else {
      p <- NULL
    }

    # Cleanup
    gc()  # Perhaps useful when run in parallel

    return(list(m = slimRpart(m), p = p))

  }

  #-----

  # Call fitModel() for each 'yord', possibly in parallel
  mout <- pbapply::pblapply(
    X = 1:length(yord),
    FUN = fitModel,
    cl = mc.cores)
  names(mout) <- yord

  #-----

  # Assemble final result
  result <- list(
    models = map(mout, "m"),
    xclass = xclass,
    xlevels = xlevels,
    yclass = yclass,
    ylevels = ylevels,
    yinner = yinner,
    ycor = compact(map(mout, "p"))
  )

  return(result)

}
