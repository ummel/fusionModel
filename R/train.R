#' Train a fusion model on donor data
#'
#' @description
#' Train a fusion model on donor data.
#'
#' @param data Data frame. Donor dataset. All categorical variables must be factors and ordered whenever possible.
#' @param y Character vector. The variables to fuse to a recipient dataset. Can include \link[base:regex]{regular expressions}.
#' @param x Character vector. Predictor variables common to donor and eventual recipient. Can include regular expressions. If NULL, all variables other than those in \code{y} and \code{weight} are used. Only one of \code{x} and \code{ignore} can be non-NULL.
#' @param ignore Character vector. Alternative way to specify predictor variables. Can include regular expressions. If non-NULL, all variables other than those in \code{y}, \code{weight}, and \code{ignore} are used. Only one of \code{x} and \code{ignore} can be non-NULL.
#' @param weight Character vector. Name of the observation weights column. If NULL (default), uniform weights are assumed.
#' @param cores Integer. Number of cores used for potential parallel operations. Passed to \code{cl} argument of \code{\link[pbapply]{pblapply}}. Ignored on Windows systems.
#' @param lasso Numeric (0-1) or NULL. Controls extent of predictor variable pre-screening via LASSO regression. If NULL, no pre-screening is performed. Default value \code{lasso = 1} invokes the least-restrictive LASSO. See Details.
#' @param maxcats Positive integer or NULL. Maximum number of levels allowed in an unordered factor predictor variable when the response (fusion) variable is also an unordered factor. Prevents excessive \code{\link[rpart]{rpart}} computation time. A K-means clustering strategy is used to cluster the predictor to no more than \code{maxcats} levels.
#' @param complexity Numeric. Passed to \code{cp} argument of \code{\link[rpart]{rpart.control}} to control complexity of decision trees.
#' @param node.obs Numeric vector of length 2. Minimum number of observations in tree nodes. First number is for numeric response variables; second for categorical. Each is passed to \code{minbucket} argument of \code{\link[rpart]{rpart.control}}.
#' @param initial Numeric vector of length 2. Controls speed/performance of the initial model-fitting routune used to determine fusion order. First number is proportion of \code{data} observations to randomly sample. Second number is the \code{cp} argument passed to \code{\link[rpart]{rpart.control}}. See Details.
#' @param ... Optional arguments passed to \code{\link[rpart]{rpart}} to control tree-building. By default \code{cp = 0}, \code{xval = 0}, \code{minbucket = 50} (\code{minbucket = 10} for discrete models), and all other arguments are left at default values.
#'
#' @details When \code{lasso} is non-NULL, predictor variables are "pre-screened" using LASSO regression via \code{\link[glmnet]{glmnet}} prior to fitting a \code{\link[rpart]{rpart}} model. Predictors with a LASSO coefficient of zero are excluded from consideration. This can speed up tree-fitting considerably when \code{data} is large. Lower values of \code{lasso} are more aggressive at excluding predictors; the LASSO \emph{lambda} is chosen such that the deviance explained is at least \code{lasso}-% of the maximum. To ensure the LASSO step is fast, pre-screening is only used for numeric, logical, and ordered factor response variables (the latter integerized).
#' @details Since determination of the fusion order only makes use of variable importance results from the initial (fully-specified) models, employing a random subset of observations and less complex models (controlled via \code{initial}) can yield a competitive fusion order at less expense.
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
#
# library(fusionModel)
# source("R/utils.R")
# source("R/fitRpart.R")
# source("R/fitDensity.R")
# source("R/fusionOrder.R")
# source("R/detectDependence.R")
#
# Inputs for ?train example - with some modification for harder test cases
# data <- recs
#
# i <- sample(which(data$propane > 0), size = 100, replace = FALSE)
# data$propane[i] <- -data$propane[i]
# data$propane_btu[i] <- -data$propane_btu[i]
# data$propane_expend[i] <- -data$propane_expend[i]
#
# data <- data %>%
#   mutate(cenac_age = ifelse(centralac_age == "Less than 2 years old", 1, 0),
#          cenac_age = ifelse(centralac_age %in% c("2 to 4 years old", "5 to 9 years old"), 4.5, cenac_age),
#          cenac_age = ifelse(centralac_age %in% c("10 to 14 years old", "15 to 19 years old", "20 years or older"), 17, cenac_age))
#
# recipient <- subset(recs, select = c(division, urban_rural, climate, income, age, race))
# y = setdiff(names(data), names(recipient))
# weight = NULL
# x <- NULL
# ignore = NULL
# maxcats = 10
# cores = 1
# lasso = 1
# complexity = 0.0001
# node.obs = c(50, 10)
# initial <- c(0.5, 0.001)

#---------------------

train <- function(data,
                  y,
                  x = NULL,
                  ignore = NULL,
                  weight = NULL,
                  cores = 1,
                  lasso = 1,
                  maxcats = 10,
                  complexity = 0,
                  node.obs = c(50, 10),
                  initial = c(0.5, 0.001),
                  ...) {

  stopifnot(exprs = {
    is.data.frame(data)
    !missing(y)
    is.null(maxcats) | (maxcats > 1 & maxcats %% 1 == 0)
    cores > 0 & cores %% 1 == 0
    is.null(lasso) | (lasso > 0 & lasso <= 1)
    complexity >= 0
    length(node.obs) == 2 & all(node.obs > 0)
    is.numeric(initial) & length(initial) == 2
  })

  # Check that no more than one of 'x' or 'ignore' is specified
  if (!is.null(x) & !is.null(ignore)) stop("Only one of 'x' or 'ignore' may be non-NULL")

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
    #!anyNA(data[yvars])  # Require non-NA values in response values
  })

  #-----

  # Create and/or check observation weights
  # Weights are normalized (mean = 1) to avoid integer overflow issues
  w <- "._wgt_."
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

  # Check for character-type variables; stop with error if any detected
  x <- sapply(data[c(xvars, yvars)], is.character)
  if (any(x)) stop("Please coerce character variables to factor:\n", paste(names(which(x)), collapse = ", "))

  #-----

  # Check that the 'xvars' and 'yvars' contain only syntactically valid names
  stopifnot(all(make.names(yvars, unique = TRUE) == yvars))
  stopifnot(all(make.names(xvars, unique = TRUE) == xvars))

  # Print variable information to console
  cat(length(yvars), "fusion variables\n")
  cat(length(xvars), "initial predictor variables\n")
  cat(nrow(data), "observations\n")

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

  # Placeholder list for 'ycor' with slots for continuous and ordered factor fusion variables (but not unordered factors)
  ycor.vars <- names(which(sapply(yclass, function(x) x[1] != "factor")))

  #-----

  # Extract data classes and levels for the 'xvars'
  x <- data[xvars]
  xclass <- lapply(x, class)
  xlevels <- lapply(x[grepl("factor", xclass)], levels)

  #-----

  # Check for no-variance (constant) variables; record them in 'sim.constant' and remove from 'data'
  constant <- names(which(sapply(data[c(xvars, yvars)], novary)))
  if (length(constant)) {
    sim.constant <- map(data[intersect(constant, yvars)], ~ na.omit(unique(.x)))
    yvars <- setdiff(yvars, constant)
    xvars <- setdiff(xvars, constant)
    cat("Detected", length(constant), "constant variable(s):\n", paste(constant, collapse = ", "), "\n")
  } else {
    sim.constant <- NULL
  }

  #-----

  # Detect if unordered factors with many levels could cause problems and issue immediate warning to console
  if (is.null(maxcats) & any(map_chr(yclass, 1L) == "factor") & any(lengths(xlevels) > 15)) {
    warning("\nBe careful: Unordered factors are present that could cause long compute times.\nSee 'maxcats' argument in ?train.\n", immediate. = TRUE)
    # NOT USED: Prompt user to confirm at console that they want to continue
    #   cat("Hold up! Unordered factors are present that could cause long compute times. See 'maxcats' argument in ?train.\n")
    #   continue <- readline(prompt = "Do you want to continue with 'maxcats = NULL'? (Y/N):\n")
    #   if (tolower(continue) != "y") stop(call. = FALSE)
  }

  #-----

  # Detect and impute any missing values in 'data'
  na.cols <- names(which(sapply(data, anyNA)))
  if (length(na.cols) > 0) {
    #cat("Imputing missing values in predictor variables...\n")
    cat("Missing values imputed for the following variable(s):\n", paste(na.cols, collapse = ", "), "\n")
    for (j in na.cols) {
      x <- data[[j]]
      ind <- is.na(x)
      set(data, i = which(ind), j = j, value = imputationValue(x, ind))
    }
  }

  #-----

  # Detect 1-to-1 linear dependencies among the variables
  # This allows certain variables to be excluded from the modeling process
  # The information needed to simulate these variables in subsequent fuse() call is retained in the 'sim.lm' and 'sim.merge' lists
  cat("Searching for derivative relationships...\n")

  # Detect near-perfect bivariate relationships where the effective R-squared is above some (high) 'threshold'
  # List of slimmed 'biglm' models that can be used for simple simulation of variables identified in 'names(sim.lm)'
  sim.lm <- detectCorrelationDependence(data,
                                        fvars = yvars,
                                        exclude = grep("..", names(data), value = TRUE, fixed = TRUE),
                                        threshold = 0.99)

  # Detect perfect bivariate categorical relationships
  # List of data frames that can be used to add/merge derivative variables when variables in 'names(sim.merge)' are known
  # NOTE: 'exclude' argument below excludes fusionACS spatial predictors ("..") by default! Not particularly safe.
  #sim.merge <- detectCategoricalDependence(data, fvars = yvars, exclude = names(sim.lm))  # This is universally appropriate but slow with many fusionACS spatial predictors
  sim.merge <- detectCategoricalDependence(data,
                                           fvars = yvars,
                                           exclude = c(names(sim.lm), grep("..", names(data), value = TRUE, fixed = TRUE)),
                                           cores = cores)

  # Vector of variables identified as derivative of other variables
  derivative <- unique(c(names(sim.lm), unlist(map(sim.merge, ~ names(.x)[-1]))))

  # Variables identified as derivative can be removed from 'data'
  if (length(derivative) > 0) {
    data <- select(data, -any_of(derivative))
    cat("Detected", length(derivative), "derivative variable(s) that can be omitted from modeling\n")
  }

  #----------------------------------
  #----------------------------------

  # RANKS MATRIX

  #cat("Building ranks matrix...\n")

  # Unordered factor variables among all retained variables
  unordered <- sapply(c(xclass, yclass), function(x) x[1] == "factor")
  unordered <- unordered[names(unordered) %in% names(data)]

  # Build 'ranks' data.table for 'xvars' that are NOT unordered factors
  ranks <- subset(data, select = names(which(!unordered)))
  for (v in names(ranks)) data.table::set(ranks, j = v, value = data.table::frank(ranks[[v]], ties.method = "average"))

  # Variable associated with each column of 'ranks' matrix
  # Updated below when dummy variable columns are created
  ranks.var <- names(ranks)

  # Create dummy variable columns in 'ranks' for the 'xvars' that ARE unordered factors
  for (v in names(which(unordered))) {
    dt <- subset(data, select = v)
    u <- xlevels[[v]]
    newv <- paste0(v, u)
    ranks.var <- c(ranks.var, rep(v, length(u)))
    ranks[newv] <- lapply(u, function(x) as.integer(dt == x))  # Works if 'ranks' is data.frame or data.table
  }

  # Safety check
  stopifnot(length(ranks.var) == ncol(ranks))

  #----------------------------------
  #----------------------------------

  # data.table with zero/nonzero binary version of numeric variables
  dzero <- as.data.table(data)
  vnum <- setdiff(names(which(sapply(data, is.numeric))), w)
  for (v in vnum) set(dzero, j = v, value = dzero[[v]] == 0)

  # Make dzero unique? Would this be faster?

  #----------------------------------
  #----------------------------------

  # DETERMINE FUSION ORDER

  cat("Determining order of fusion variables...\n")

  # Adjusts 'yvars' for removal of variables that can be deduced by 1-to-1 linear relationships
  yvars <- intersect(yvars, names(data))

  # Fit full models, including yvars as predictors
  # Only 'variable.importance' slot is returned
  full.varimp <- pbapply::pblapply(X = yvars, FUN = function(y) {
    m <- fitRpart(y = y,
                  x = c(xvars, setdiff(yvars, y)),
                  w = w,
                  data = data,
                  n = initial[1] * nrow(data),
                  maxcats = maxcats,
                  linear = FALSE,
                  lasso.threshold = lasso,
                  args = list(cp = initial[2],
                              minbucket = ifelse(y %in% ycont, node.obs[1], node.obs[2]),
                              xval = 0,
                              maxcompete = 0)
    )
  }, cl = cores) %>%
    setNames(yvars)

  # Find preferred order of yvars
  yord <- fusionOrder(varimp = map(full.varimp, "variable.importance"))

  # Cleanup
  rm(full.varimp)
  gc()

  #----------------------------------
  #----------------------------------

  # FIT MODELS

  cat("Building fusion models...\n")

  # Function to fit rpart() model and calculate rank correlations for i-th response variable in 'yord'
  fitModel <- function(i) {

    y <- yord[i]

    # Is 'y' continuous?
    cont <- y %in% ycont

    # Prior response variables in the fusion sequence
    yprior <- if (i == 1) NULL else yord[1:(i - 1)]

    #--------

    # Restrict data to rows with non-NA response values
    ind <- which(!is.na(data[[y]]))
    d <- data[ind, c(w, y, yprior, xvars)]

    #--------

    # Detect if there is a numeric predictor with zero/nonzero structure that is derivative of 'y'
    # If a zero/nonzero binary version of a predictor can explain 'y', then we want to include that binary version and force its inclusion
    # Example: one of the predictors is natural gas consumption (numeric) and 'y' is categorical indicating if household uses natural gas
    xbin <- catDependence(data = dzero,
                          targets = y,
                          candidates = intersect(vnum, c(xvars, yprior)),
                          cores = 1L)

    # Update 'd' with binary version of 'xbin' derivative variables
    # In such cases, the binary variable generally
    if (is.null(xbin)) {
      xbinary <- NULL
    } else {
      xbin <- xbin[[1]][1]  # Retains only a single explanator variable, even if multiple exist
      xbinary <- paste0(xbin, "_BINARY_")
      d <- cbind(d, setnames(x = dzero[, ..xbin], new = xbinary))

      # Extract any unique 'y' values associated with binary predictor
      xmerge <- d %>%
        select(all_of(c(xbinary, y))) %>%
        distinct() %>%
        group_by_at(xbinary) %>%
        filter(n() == 1) %>%
        ungroup()

      # Subset 'd' so that decision tree only fit to subset of values not explained by 'xmerge'
      d <- subset(d, !d[[y]] %in% xmerge[[y]])

    }

    #--------

    # Fit prediction model
    m <- fitRpart(y = y,
                  x = c(xvars, yprior, xbinary),
                  w = w,
                  data = d,
                  maxcats = maxcats,
                  n = nrow(d),
                  linear = TRUE,
                  lasso.threshold = lasso,
                  args = list(cp = complexity,
                              minbucket = ifelse(y %in% ycont, node.obs[1], node.obs[2]),
                              xval = 0,
                              maxcompete = 0)
    )

    #--------

    # Add 'xbin' slot indicating which numeric variables where binarized prior to model-fitting
    # Add 'xmerge' slot indicating which binary numeric predictors perfectly explain an outcome value
    # If the response variable is categorical, remove the "known" level from the 'ylevels' attribute of the rpart object (causes prediction errors otherwise)
    if (!is.null(xbinary)) {
      m$xbin <- xbin
      m$xmerge <- xmerge
      if (class(m) == "rpart" & (is.factor(xmerge[[2]]) | is.logical(xmerge[[2]]))) {
        attr(m, "ylevels") <- setdiff(attr(m, "ylevels"), xmerge[[2]])
      }
    }

    #--------

    # If 'y' is continuous...
    if (cont) {

      # Add additional slots to rpart model output
      if (class(m) == "rpart") {

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

    }

    #--------

    # Calculate (rank) correlation between 'y' and existing variables in 'ranks' (if necessary; only continuous and ordered factors)
    if (y %in% ycor.vars) {
      #ind <- which(as.numeric(data[[y]]) != 0)  # Restrict correlation calculation to non-zero observations
      ind <- 1:nrow(data)  # No restriction on correlation calculation
      j <- ranks.var == y  # Column in 'ranks' associated with 'y'
      k <- ranks.var %in% c(xvars, yprior)  # Columns in 'ranks' associated with predictors of 'y'
      p <- suppressWarnings(cor(ranks[ind, k], ranks[ind, j]))[, 1]  # Will silently return NA if no variance in a particular column
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
    cl = cores) %>%
    setNames(yord)

  # test <- sapply(mout, class)
  # which(test == "try-error")
  # check <- lapply(which(test == "try-error"), fitModel)

  # Troubleshoot
  #for (i in 1:length(yord)) fitModel(i)

  # Check certain outputs
  # map(mout, ~ .x$m$xbin)
  # table(map_chr(mout, ~class(.x$m)))

  #-----

  # Assemble final result
  result <- list(
    models = map(mout, "m"),
    xclass = xclass,
    xlevels = xlevels,
    yclass = yclass,
    ylevels = ylevels,
    yinner = yinner,
    youter = youter,
    ycor = compact(map(mout, "p")),
    derivative = list(constant = sim.constant,
                      model = sim.lm,
                      merge = sim.merge)
  )

  return(result)

}
