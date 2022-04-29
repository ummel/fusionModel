#' Train a fusion model using Classification and Regression Trees (CART)
#'
#' @description
#' Train a fusion model on donor data using CART implemented via \code{\link[rpart]{rpart}}.
#'
#' @param data Data frame. Donor dataset. Categorical variables must be factors and ordered whenever possible.
#' @param y Character vector. The variables to fuse to a recipient dataset.
#' @param x Character vector. Predictor variables common to donor and eventual recipient.
#' @param weight Character vector. Name of the observation weights column. If NULL (default), uniform weights are assumed.
#' @param order Character vector. Order in which to fuse \code{y}. If NULL (default), a pseudo-optimal order is determined internally.
#' @param deriv Logical. Should algorithm check for derivative relationships prior to building fusion models?
#' @param smoothed Logical. Should the synthetic values be smoothed via KDE? Default is \code{FALSE}.
#' @param cores Integer. Number of cores used for parallel operations. Passed to \code{cl} argument of \code{\link[pbapply]{pblapply}}. Ignored on Windows systems.
#' @param lasso Numeric (0-1) or NULL. Controls extent of predictor variable pre-screening via LASSO regression. If NULL (default), no screening is performed. \code{lasso = 1} invokes the least-restrictive screening. See Details.
#' @param maxcats Positive integer or NULL. Maximum number of levels allowed in an unordered factor predictor variable when the response (fusion) variable is also an unordered factor. Prevents excessive \code{\link[rpart]{rpart}} computation time. A K-means clustering strategy is used to cluster the predictor to no more than \code{maxcats} levels.
#' @param complexity Numeric. Passed to \code{cp} argument of \code{\link[rpart]{rpart.control}} to control complexity of decision trees.
#' @param cvfolds Integer. Number of cross-validation folds used by \code{\link[rpart]{rpart}} to determine optimal tree complexity. Default is no cross-validation (\code{cvfolds = 0}).
#' @param cvfactor Numeric. Controls how decision trees are pruned when \code{cvfolds > 0}. \code{cvfactor = 0} (the default) selects the tree complexity that minimizes the cross-validation error. \code{cvfactor = 1} is equivalent to Breiman's "1-SE" rule.
#' @param node.obs Numeric vector of length 2. Minimum number of observations in tree nodes. First number is for numeric response variables; second for categorical. Each is passed to \code{minbucket} argument of \code{\link[rpart]{rpart.control}}.
#' @param initial Numeric vector of length 2. Controls speed/performance of the initial model-fitting routune used to determine fusion order. First number is proportion of \code{data} observations to randomly sample. Second number is the \code{cp} argument passed to \code{\link[rpart]{rpart.control}}. See Details.
#'
#' @details When \code{lasso} is non-NULL, predictor variables are "pre-screened" using LASSO regression via \code{\link[glmnet]{glmnet}} prior to fitting a \code{\link[rpart]{rpart}} model. Predictors with a LASSO coefficient of zero are excluded from consideration. This can speed up tree-fitting considerably when \code{data} is large. Lower values of \code{lasso} are more aggressive at excluding predictors; the LASSO \emph{lambda} is chosen such that the deviance explained is at least \code{lasso}-% of the maximum. To ensure the LASSO step is fast, pre-screening is only used for numeric, logical, and ordered factor response variables (the latter integerized).
#' @details Since determination of the fusion order only makes use of variable importance results from the initial (fully-specified) models, employing a random subset of observations and less complex models (controlled via \code{initial}) can yield a competitive fusion order at less expense.
#'
#' @return A list containing trained model information to be passed to \link{fuse}.
#'
#' @examples
#' # Build a fusion model using RECS microdata
#' ?recs
#' fusion.vars <- c("electricity", "natural_gas", "aircon")
#' predictor.vars <- names(recs)[2:12]
#' fit <- train(data = recs, y = fusion.vars, x = predictor.vars)
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
# data <- recs[1:28]
#
# recipient <- subset(recs, select = c(weight, division, urban_rural, climate, income, age, race, hh_size, televisions))
# y = setdiff(names(data), names(recipient))
# weight <- "weight"
# x <- setdiff(names(recipient), weight)
# maxcats = 10
# cores = 1
# lasso = 1
# complexity = 0.0001
# node.obs = c(50, 10)
# initial <- c(0.5, 0.001)

# fusion.cart <- train(data = data, y = y, x = x, weight = weight)
# saveRDS(fusion.cart, "fusion_cart.rds")

#---------------------

train <- function(data,
                  y,
                  x,
                  weight = NULL,
                  order = NULL,
                  deriv = TRUE,
                  smoothed = FALSE,
                  cores = 1,
                  lasso = NULL,
                  maxcats = 12,
                  complexity = 0,
                  cvfolds = 0,
                  cvfactor = 0,
                  node.obs = c(20, 10),
                  initial = c(1, 0)) {

  stopifnot(exprs = {
    is.data.frame(data)
    is.character(y)
    is.character(x)
    is.logical(smoothed)
    cores > 0 & cores %% 1 == 0
    is.null(lasso) | (lasso > 0 & lasso <= 1)
    is.null(maxcats) | (maxcats > 1 & maxcats %% 1 == 0)
    complexity >= 0
    cvfolds >= 0 & cvfolds %% 1 == 0
    cvfactor >= 0 & length(cvfactor) == 1
    length(node.obs) == 2 & all(node.obs > 0)
    is.numeric(initial) & length(initial) == 2
  })

  # Check that no more than one of 'x' or 'ignore' is specified
  #if (!is.null(x) & !is.null(ignore)) stop("Only one of 'x' or 'ignore' may be non-NULL")

  #-----

  # Detect fusion variables
  nms <- names(data)
  yvars <- validNames(y, nms, exclude = FALSE)

  # Detect predictor variables, dependent on how 'x' and 'ignore' arguments are specified
  xvars <- setdiff(nms, c(yvars, weight))
  #xvars <- validNames(ignore, xvars, exclude = TRUE)
  #if (!is.null(x)) xvars <- validNames(x, xvars, exclude = FALSE)
  xvars <- validNames(x, xvars, exclude = FALSE)

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

  if (length(yvars) > 1 & deriv) {

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

  } else {
    sim.lm <- NULL
    sim.merge <- NULL
  }

  #----------------------------------
  #----------------------------------

  # # RANKS MATRIX
  #
  # #cat("Building ranks matrix...\n")
  #
  # # Unordered factor variables among all retained variables
  # unordered <- sapply(c(xclass, yclass), function(x) x[1] == "factor")
  # unordered <- unordered[names(unordered) %in% names(data)]
  #
  # # Build 'ranks' data.table for 'xvars' that are NOT unordered factors
  # ranks <- subset(data, select = names(which(!unordered)))
  # for (v in names(ranks)) data.table::set(ranks, j = v, value = data.table::frank(ranks[[v]], ties.method = "average"))
  #
  # # Variable associated with each column of 'ranks' matrix
  # # Updated below when dummy variable columns are created
  # ranks.var <- names(ranks)
  #
  # # Create dummy variable columns in 'ranks' for the 'xvars' that ARE unordered factors
  # for (v in names(which(unordered))) {
  #   dt <- subset(data, select = v)
  #   u <- xlevels[[v]]
  #   newv <- paste0(v, u)
  #   ranks.var <- c(ranks.var, rep(v, length(u)))
  #   ranks[newv] <- lapply(u, function(x) as.integer(dt == x))  # Works if 'ranks' is data.frame or data.table
  # }
  #
  # # Safety check
  # stopifnot(length(ranks.var) == ncol(ranks))

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

  if (length(yvars) > 1 & is.null(order)) {

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
                    cvfactor = cvfactor,
                    args = list(cp = initial[2],
                                minbucket = ifelse(y %in% ycont, node.obs[1], node.obs[2]),
                                xval = cvfolds,
                                maxcompete = 0)
      )
    }, cl = cores) %>%
      setNames(yvars)

    # Find preferred order of yvars
    yord <- fusionOrder(varimp = map(full.varimp, "variable.importance"))

    # Cleanup
    rm(full.varimp)
    gc()

  } else {

    if (length(yvars) == 1) {
      yord <- yvars
    } else {
      stopifnot(all(yvars %in% order))
      yord <- order
    }

  }

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
                  cvfactor = cvfactor,
                  args = list(cp = complexity,
                              minbucket = ifelse(y %in% ycont, node.obs[1], node.obs[2]),
                              xval = cvfolds,
                              maxcompete = 0)
    )

    #--------

    # Add 'xbin' slot indicating which numeric variables where binarized prior to model-fitting
    # Add 'xmerge' slot indicating which binary numeric predictors perfectly explain an outcome value
    # If the response variable is categorical, remove the "known" level from the 'ylevels' attribute of the rpart object (causes prediction errors otherwise)
    if (!is.null(xbinary)) {
      m$xbin <- xbin
      m$xmerge <- xmerge
      #if (class(m) == "rpart" & ((is.factor(xmerge[[2]] & !is.ordered(xmerge[[2]])) | is.logical(xmerge[[2]]))) {  # Not clear if is.logical condition is correct here! Not thoroughly tested
      if (class(m) == "rpart" & (is.factor(xmerge[[2]]) & !is.ordered(xmerge[[2]]))) {  # Omitting is.logical condition until clear it is necessary
        attr(m, "ylevels") <- setdiff(attr(m, "ylevels"), xmerge[[2]])
      }
    }

    #--------

    # If 'y' is continuous...
    if (cont) {

      # Add additional slots to rpart model output
      if (class(m) == "rpart") {

        # Unique node identifiers
        nodes <- sort(unique(m$where))

        #---

        if (smoothed) {

          # Placeholder matrix for the smoothed quantile values
          # The 'Q' values describe the shape of the quantile function (i.e. inverse CDF)
          Qn <- 500  # Ideally, would be a train() control argument
          m$Q <- matrix(0L, nrow = Qn, ncol = length(nodes))
          colnames(m$Q) <- nodes

          # Fit density to observations in each node
          for (n in nodes) {

            # Index identifying observations in node 'n'
            nind <- m$where == n

            # Node density result
            fd <- fitDensity(x = d[nind, y][[1]],
                             w = d[nind, w][[1]],
                             inner.range = yinner[[y]],
                             outer.range = youter[[y]],
                             N = Qn)

            # Add quantile function values to 'm'
            m$Q[, as.character(n)] <- fd

          }

        } else {

          # Return 'Q' as named list, each slot containing a 2-column matrix with y-values and observation weights found in each leaf node
          m$Q <- vector(mode = "list", length = length(nodes))
          names(m$Q) <- as.character(nodes)
          for (n in nodes) {
            nind <- m$where == n
            m$Q[[as.character(n)]] <- cbind(d[nind, y], d[nind, w])
          }

        }

      }

    }

    #--------

    # # Calculate (rank) correlation between 'y' and existing variables in 'ranks' (if necessary; only continuous and ordered factors)
    # if (y %in% ycor.vars) {
    #   #ind <- which(as.numeric(data[[y]]) != 0)  # Restrict correlation calculation to non-zero observations
    #   ind <- 1:nrow(data)  # No restriction on correlation calculation
    #   j <- ranks.var == y  # Column in 'ranks' associated with 'y'
    #   k <- ranks.var %in% c(xvars, yprior)  # Columns in 'ranks' associated with predictors of 'y'
    #   p <- suppressWarnings(cor(ranks[ind, k], ranks[ind, j]))[, 1]  # Will silently return NA if no variance in a particular column
    #   p <- cleanNumeric(p[!is.na(p)], tol = 0.001)
    # } else {
    #   p <- NULL
    # }

    # Cleanup
    gc()  # Perhaps useful when run in parallel

    #return(list(m = slimRpart(m), p = p))
    return(list(m = slimRpart(m)))

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
    #ycor = compact(map(mout, "p")),
    derivative = list(constant = sim.constant,
                      model = sim.lm,
                      merge = sim.merge)
  )

  return(result)

}
