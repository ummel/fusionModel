#' Train a fusion model
#'
#' @description
#' Train a fusion model on "donor" data using sequential \href{https://lightgbm.readthedocs.io/en/latest/}{LightGBM} models to model the conditional distributions. The resulting fusion model (.fsn file) can be used with \code{\link{fuse}} to simulate outcomes for a "recipient" dataset.
#'
#' @param data Data frame. Donor dataset. Categorical variables must be factors and ordered whenever possible.
#' @param y Character or list. Variables in \code{data} to eventually fuse to a recipient dataset. Variables are fused in the order provided. If \code{y} is a list, each entry is a character vector possibly indicating multiple variables to fuse as a block.
#' @param x Character or list. Predictor variables in \code{data} common to donor and eventual recipient. If a list, each slot specifies the \code{x} predictors to use for each \code{y}.
#' @param fsn Character. File path where fusion model will be saved. Must use \code{.fsn} suffix.
#' @param weight Character. Name of the observation weights column in \code{data}. If NULL (default), uniform weights are assumed.
#' @param nfolds Numeric. Number of cross-validation folds used for LightGBM model training. Or, if \code{nfolds < 1}, the fraction of observations to use for training set; remainder used for validation (faster than cross-validation).
#' @param nquantiles Numeric. Number of quantile models to train for continuous \code{y} variables (along with the conditional mean). \code{nquantiles} evenly-distributed percentiles are used. The default \code{nquantiles = 3} is usually sufficient and yields quantile models for the 16.6th, 50th, and 83.3rd percentiles.
#' @param nclusters Numeric. Maximum number of k-means clusters to use. Higher is better but at computational cost. \code{nclusters = 0} or \code{nclusters = Inf} turn off clustering.
#' @param krange Numeric. Minimum and maximum number of nearest neighbors to use for construction of continuous conditional distributions. Higher \code{max(krange)} is better but at computational cost.
#' @param hyper List. LightGBM hyperparameters to be used during model training. If \code{NULL}, default values are used. See Details and Examples.
#' @param fork Logical. Should parallel processing via forking be used, if possible? See Details.
#' @param cores Integer. Number of physical CPU cores used for parallel computation. When \code{fork = FALSE} or on Windows platform (since forking is not possible), the fusion variables/blocks are processed serially but LightGBM uses \code{cores} for internal multithreading via OpenMP. On a Unix system, if \code{fork = TRUE}, \code{cores > 1}, and \code{cores <= length(y)} then the fusion variables/blocks are processed in parallel via \code{\link[parallel]{mclapply}}.
#'
#' @details When \code{y} is a list, each slot indicates either a single variable or, alternatively, multiple variables to fuse as a block. Variables within a block are sampled jointly from the original donor data during fusion. See Examples.
#' @details The fusion model written to \code{fsn} is a zipped archive created by \code{\link[zip]{zip}} containing models and data required by \code{\link{fuse}}.
#' @details The \code{hyper} argument can be used to specify the LightGBM hyperparameter values over which to perform a "grid search" during model training. \href{https://lightgbm.readthedocs.io/en/latest/Parameters.html}{See here} for the full list of parameters. For each combination of hyperparameters, \code{nfolds} cross-validation is performed using \code{\link[lightgbm]{lgb.cv}} with an early stopping condition. The parameter combination with the lowest loss function value is used to fit the final model via \code{\link[lightgbm]{lgb.train}}. The more candidate parameter values specified in \code{hyper}, the longer the processing time. If \code{hyper = NULL}, a single set of parameters is used LightGBM default values. Typically, users will only have reason to specify the following parameters via \code{hyper}:
#' @details \itemize{
#'   \item boosting
#'   \item num_leaves
#'   \item bagging_fraction
#'   \item feature_fraction
#'   \item max_depth
#'   \item min_data_in_leaf
#'   \item num_iterations
#'   \item learning_rate
#'   \item max_bin
#'   \item min_data_in_bin
#'   \item max_cat_threshold
#'  }
#' @details Testing with small-to-medium size datasets suggests that forking is typically faster than OpenMP multithreading (the default). However, forking will sometimes "hang" (continue to run with no CPU usage or error message) if an OpenMP process has been previously used in the same session. The issue appears to be related to Intel's OpenMP implementation (\href{https://github.com/Rdatatable/data.table/issues/2418}{see here}). This can be triggered when other operations are called before \code{train()} that use \code{\link[data.table]{data.table}} or \code{\link[fst]{fst}} in multithread mode. If you experience hanged forking, try calling \code{data.table::setDTthreads(1)} and \code{fst::threads_fst(1)} immediately after \code{library(fusionModel)} in a new session.
#'
#' @return A fusion model object (.fsn) is saved to \code{fsn}.
#'
#' @examples
#' # Build a fusion model using RECS microdata
#' # Note that "fusion_model.fsn" will be written to working directory
#' ?recs
#' fusion.vars <- c("electricity", "natural_gas", "aircon")
#' predictor.vars <- names(recs)[2:12]
#' fsn.path <- train(data = recs, y = fusion.vars, x = predictor.vars)
#'
#' # When 'y' is a list, it can specify variables to fuse as a block
#' fusion.vars <- list("electricity", "natural_gas", c("heating_share", "cooling_share", "other_share"))
#' fusion.vars
#' train(data = recs, y = fusion.vars, x = predictor.vars)
#'
#' # Specify a single set of LightGBM hyperparameters
#' train(data = recs, y = fusion.vars, x = predictor.vars,
#'       hyper = list(boosting = "goss",
#'                    feature_fraction = 0.8,
#'                    num_iterations = 300
#'       ))
#'
#' # Specify a range of LightGBM hyperparameters to search over
#' # This takes longer, because there are more models to test
#' train(data = recs, y = fusion.vars, x = predictor.vars,
#'       hyper = list(num_leaves = c(10, 30),
#'                    feature_fraction = c(0.7, 0.9),
#'                    num_iterations = 50
#'       ))
#' @export

#---------------------

# Manual testing
# library(fusionModel)
# library(data.table)
# library(dplyr)
# source("R/utils.R")
# source("R/fitLGB.R")

# # Inputs for testing - with some modification for harder test cases
# data <- recs[1:28]
# recipient <- subset(recs, select = c(weight, income:heat_type))
# y = setdiff(names(data), names(recipient))[1:8]
# weight <- "weight"
# x <- setdiff(names(recipient), weight)
# fsn = "fusion_model.fsn"
# nfolds = 5
# nquantiles = 3
# nclusters = 2000
# hyper = NULL
# krange = c(10, 300)
# fork = FALSE
# cores = 2

# Try with variable block
#y <- c(list(y[1:2]), y[-c(1:2)])

#---------------------

train <- function(data,
                  y,
                  x,
                  fsn = "fusion_model.fsn",
                  weight = NULL,
                  nfolds = 5,
                  nquantiles = 3,
                  nclusters = 2000,
                  krange = c(10, 500),
                  hyper = NULL,
                  fork = FALSE,
                  cores = 1) {

  t0 <- Sys.time()

  # TO DO: Make data.table operations throughout (just applies to pre-loop checks)
  if (is.data.table(data)) data <- as.data.frame(data)

  # Create correct 'y' and 'ylist', depending on input type
  if (is.list(y)) {
    ylist <- y
    y <- unlist(ylist)
  } else {
    ylist <- as.list(y)
  }

  # Create correct 'x' and 'xlist', depending on input type
  if (is.list(x)) {
    xlist <- x
    x <- unique(unlist(x))
  } else {
    xlist <- rep(list(x), length(ylist))
  }

  #---

  # Check validity of other inputs
  stopifnot(exprs = {
    is.data.frame(data)
    all(y %in% names(data))
    !any(c("M", "W..", "R..") %in% y)  # Reserved variable names
    all(x %in% names(data))
    length(ylist) == length(xlist)
    length(intersect(y, x)) == 0
    is.character(fsn) & endsWith(fsn, ".fsn")
    is.null(weight) | (length(weight) == 1 & weight %in% names(data) & !weight %in% c(y, x))
    nfolds > 0
    nquantiles > 0
    nclusters >= 0
    all(krange >= 5) & length(krange) == 2
    is.null(hyper) | is.list(hyper)
    is.logical(fork)
    cores > 0 & cores %% 1 == 0 & cores <= parallel::detectCores(logical = FALSE)
  })

  #---

  # Calculate the percentiles to use for quantile models
  ptiles <- seq(from = 1 / nquantiles / 2, length.out = nquantiles, by = 2 * 1 / nquantiles / 2)

  # Determine if parallel forking will be used
  # This forces use of OpenMP if there are more cores than fusion steps
  fork <- fork & .Platform$OS.type == "unix" & cores > 1 & length(ylist) > 1 & cores <= length(ylist)

  # Check 'fsn' path and create parent directories, if necessary
  dir <- normalizePath(dirname(fsn), mustWork = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  fsn <- file.path(dir, basename(fsn))

  # Temporary directory to save outputs to
  td <- tempfile()
  dir.create(td, showWarnings = FALSE)

  # Determine fusion variable directory prefixes for saving to disk
  pfixes <- formatC(seq_along(ylist), width = nchar(length(ylist)), format = "d", flag = "0")

  #-----

  # Create and/or check observation weights
  # LightGBM can handel numeric weights; ctree() requires integer weights
  # A different version is created for each
  W.lgb <- if (is.null(weight)) {
    rep(1L, nrow(data))
  } else {
    if (anyNA(data[[weight]])) stop("Missing (NA) values are not allowed in 'weight'")
    data[[weight]] / mean(data[[weight]]) # Scaled weights to avoid numerical issues
  }

  # Integerized version of the observation weights
  W.int <- integerize(W.lgb, mincor = 0.999)

  #-----

  # Check data and variable name validity
  if (anyNA(data[y])) stop("Missing (NA) values are not allowed in 'y'")

  # Check that the 'xvars' and 'yvars' contain only syntactically valid names
  bad <- setdiff(c(x, y), make.names(c(x, y)))
  if (length(bad)) stop("Fix invalid column names (see ?make.names):\n", paste(bad, collapse = ", "))

  # Check for character-type variables; stop with error if any detected
  xc <- sapply(data[c(x, y)], is.character)
  if (any(xc)) stop("Coerce character variables to factor:\n", paste(names(which(xc)), collapse = ", "))

  # Check for character-type variables; stop with error if any detected
  # Check for no-variance (constant) variables
  # Detect and impute any missing values in 'x' variables
  data <- checkData(data, y, x)
  x <- intersect(x, names(data))

  #-----

  # Extract data classes and levels for the 'yvars'
  d <- data[y]
  yclass <- lapply(d, class)
  ylevels <- lapply(d[grepl("factor", yclass)], levels)

  # Extract data classes and levels for the 'xvars'
  d <- data[x]
  xclass <- lapply(d, class)
  xlevels <- lapply(d[grepl("factor", xclass)], levels)

  # Determine LightGBM response variable types
  ytype <- ifelse(sapply(data[y], is.factor), "multiclass", "continuous")
  ytype <- ifelse(sapply(data[y], function(x) is.logical(x) | length(levels(x)) == 2), "binary", ytype)

  # Nominal/unordered categorical variables (both x and y variables)
  # Passed to LightGBM so it knows which variables to treat as unordered factors
  nominal <- names(select(data, where(is.factor) & !where(is.ordered)))

  rm(d)

  #-----

  # Print variable information to console
  xvars <- x
  yvars <- y
  cat(length(yvars), "fusion variables\n")
  cat(length(xvars), "initial predictor variables\n")
  cat(nrow(data), "observations\n")

  # Report if using different sets of predictor variables across the fusion variables
  if (all(lengths(xlist) == length(x))) {
    cat("Using all available predictors for each fusion variable\n")
  } else {
    cat("Using specified set of predictors for each fusion variable\n")
  }

  # Limit 'data' to the necessary variables
  data <- data[c(xvars, yvars)]

  # Coerce 'data' to sparse numeric matrix for use with LightGBM
  dmat <- to_mat(data)
  rm(data)

  # Write to disk (donor.fst)
  # Save the necessary response/fusion variables and observation weights to disk
  # This information is used be fuse() to select simulated values
  # This should include all fusion variables other than solo categorical (multiclass or logical) variables
  ysave <- y[ytype == "continuous" | y %in% unlist(ylist[lengths(ylist) > 1])]
  dtemp <- as.data.frame(dmat[, ysave, drop = FALSE])
  dtemp$W <- W.int
  fst::write_fst(x = dtemp, path = file.path(td, "donor.fst"), compress = 100)
  rm(dtemp)

  #-----

  # Set the 'hyper.default' object
  if (is.null(hyper)) hyper <- list()

  # Default hyperparameter values, per LightGBM documentation
  # Parameters: https://lightgbm.readthedocs.io/en/latest/Parameters.html
  # Parameter tuning: https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html
  # More details: https://sites.google.com/view/lauraepp/parameters
  hyper.default <- list(
    boosting = "gbdt",
    num_leaves = 31,
    bagging_fraction = 1,
    feature_fraction = 1,
    min_data_in_leaf = 20,
    num_iterations = 100,
    learning_rate = 0.1,
    max_depth = -1,
    max_bin = 255,
    min_data_in_bin = 3,
    max_cat_threshold = 32
  )

  # Use default hyperparameters, if not specified by user
  for (v in names(hyper.default)) {
    if (!v %in% names(hyper)) {
      hyper[[v]] <- hyper.default[[v]]
    }
  }

  # Set the number of LightGBM threads
  # If forking, LightGBM uses single core internally
  hyper$num_threads <- ifelse(fork, 1L, cores)

  # The 'dataset' parameters 'max_bin', 'min_data_in_bin', and 'max_cat_threshold' can only have a single value (they are not eligible to be varied within fitLGB)
  # Note that 'feature_pre_filter' is forced to FALSE if there are multiple 'min_data_in_leaf' values within 'hyper'
  # https://lightgbm.readthedocs.io/en/latest/Parameters.html#dataset-parameters
  for (v in c("max_bin", "min_data_in_bin", "max_cat_threshold")) {
    if (length(hyper[[v]]) > 1) {
      hyper[[v]] <- hyper[[v]][1]
      cat("Only one", v, "value allowed. Using:", hyper[[v]], "\n")
    }
  }
  dparams <- list(max_bin = hyper$max_bin,
                  min_data_in_bin = hyper$min_data_in_bin,
                  max_cat_threshold = hyper$max_cat_threshold,
                  feature_pre_filter = length(hyper$min_data_in_leaf) == 1)
  hyper$max_bin <- NULL
  hyper$min_data_in_bin <- NULL
  hyper$max_cat_threshold <- NULL

  # Create hyperparameter grid to search
  hyper.grid <- hyper %>%
    expand.grid() %>%
    mutate(
      bagging_fraction = ifelse(boosting == "goss", 1, bagging_fraction),  # Bagging not possible with 'goss'
      bagging_freq = ifelse(bagging_fraction == 1, 0, 1)  # Turn on bagging if indicated by 'bagging_fraction'
    ) %>%
    distinct() %>%
    split(seq(nrow(.)))

  #-----

  # Function to build LightGBM prediction model for step 'i' in 'ylist'
  buildFun <- function(i, verbose = FALSE) {

    v <- ylist[[i]]
    block <- length(v) > 1

    # Print message to console
    if (verbose) cat("Training step ", i, " of ", length(pfixes), ": ", paste(v, collapse = ", "), "\n", sep = "")

    # 'y' variables from prior clusters to be included as predictors
    yv <- if (i == 1) NULL else unlist(ylist[1:(i - 1)])

    # Full set of predictor variables, including 'y' from clusters earlier in sequence
    xv <- c(xlist[[i]], yv)

    path <- file.path(td, pfixes[i])
    dir.create(path)

    # Placeholder objects for necessary outputs
    cd <- yi <- wi <- ycenter <- yscale <- kcenters <- nn <- r2 <- NULL

    #---

    for (y in v) {

      # Response variable type and values
      Y <- dmat[, y]
      type <- ytype[[y]]

      # Placeholder objects
      zc <- mc <- qc <- dtrain <- dvalid <- NULL
      ti <- pi <- rep(TRUE, length(Y))

      #-----

      # Build zero model, if necessary
      if (type == "continuous" & inflated(Y)) {

        # Indices to identify which observations to use for subsequent training (ti) and prediction (pi)
        ti <- Y != 0
        if (!block) pi <- ti

        #---

        # Build LightGBM datasets for zero-inflated model

        # List indicating assignment of folds OR vector indicating training observations when nfolds <= 1
        # List indicating random assignment of folds
        cv.folds <- stratify(y = (Y == 0), ycont = FALSE, tfrac = nfolds, cv_list = TRUE)

        # Create full LGB training dataset with all available observations
        dfull <- lightgbm::lgb.Dataset(data = dmat[, xv],
                                       label = as.integer(Y == 0),
                                       weight = W.lgb,
                                       categorical_feature = intersect(xv, nominal),
                                       params = dparams) %>%
          lightgbm::lgb.Dataset.construct()

        # Create 'dtrain' and 'dvalid' sets, if requested
        if (nfolds <= 1) {
          ind <- which(cv.folds)
          dtrain <- lightgbm::lgb.Dataset(data = dmat[ind, xv],
                                          label = as.integer(Y == 0)[ind],
                                          weight = W.lgb[ind],
                                          categorical_feature = intersect(xv, nominal),
                                          params = dparams) %>%
            lightgbm::lgb.Dataset.construct()

          dvalid <- lightgbm::lgb.Dataset(data = dmat[-ind, xv],
                                          label = as.integer(Y == 0)[-ind],
                                          weight = W.lgb[-ind],
                                          categorical_feature = intersect(xv, nominal),
                                          params = dparams,
                                          reference = dtrain) %>%
            lightgbm::lgb.Dataset.construct()
        }

        #---

        # Set the loss function and performance metric used with lightGBM
        params.obj <- list(objective = "binary", metric = "binary_logloss")

        # Fit model
        zmod <- fitLGB(dfull = dfull,
                       dtrain = dtrain,
                       dvalid = dvalid,
                       hyper.grid = hyper.grid,
                       params.obj = params.obj,
                       cv.folds = cv.folds)

        # Save LightGBM mean model (m.txt) to disk
        lightgbm::lgb.save(booster = zmod, filename = file.path(path, paste0(y, "_z.txt")))

        # Conditional probability of zero for the non-zero training observations
        # Note that 'zc' will be NULL in case of single continuous variables (zeros are simulated directly by 'zmod')
        if (block) {
          zc <- matrix(predict(object = zmod, data = dmat[pi, xv], reshape = TRUE))
          colnames(zc) <- paste0(y, "_z")
        }

        rm(zmod)

      }

      #---

      # Build LightGBM datasets for mean and quantile models

      # List indicating assignment of folds OR vector indicating training observations when nfolds <= 1
      cv.folds <- stratify(y = Y[ti], ycont = (type == "continuous"), tfrac = nfolds, ntiles = 10, cv_list = TRUE)

      # Create full LGB training dataset with all available observations
      dfull <- lightgbm::lgb.Dataset(data = dmat[ti, xv],
                                     label = Y[ti],
                                     weight = W.lgb[ti],
                                     categorical_feature = intersect(xv, nominal),
                                     params = dparams) %>%
        lightgbm::lgb.Dataset.construct()

      # Create 'dtrain' and 'dvalid' sets, if requested
      if (nfolds <= 1) {
        ind <- which(ti)[cv.folds]
        dtrain <- lightgbm::lgb.Dataset(data = dmat[ind, xv],
                                        label = Y[ind],
                                        weight = W.lgb[ind],
                                        categorical_feature = intersect(xv, nominal),
                                        params = dparams) %>%
          lightgbm::lgb.Dataset.construct()

        dvalid <- lightgbm::lgb.Dataset(data = dmat[-ind, xv],
                                        label = Y[-ind],
                                        weight = W.lgb[-ind],
                                        categorical_feature = intersect(xv, nominal),
                                        params = dparams,
                                        reference = dtrain) %>%
          lightgbm::lgb.Dataset.construct()
      }

      #---

      # Set the loss function and performance metric used with lightGBM
      # This is here so that the Huber loss 'alpha' is set based on observed MAD
      # See full list of built-in metrics: https://lightgbm.readthedocs.io/en/latest/Parameters.html#metric-parameters
      # And here about Huber: https://stats.stackexchange.com/questions/465937/how-to-choose-delta-parameter-in-huber-loss-function
      params.obj <- switch(type,
                           #continuous = list(objective = "regression", metric = "huber", alpha = 1.35 * mad(Y)),  # Will throw error if mad(Y) is zero
                           continuous = list(objective = "regression", metric = "l2"),
                           multiclass = list(objective = "multiclass", metric = "multi_logloss", num_class = 1L + max(Y)),
                           binary = list(objective = "binary", metric = "binary_logloss"))

      #-----

      # Fit mean response model
      mmod <- fitLGB(dfull = dfull,
                     dtrain = dtrain,
                     dvalid = dvalid,
                     hyper.grid = hyper.grid,
                     params.obj = params.obj,
                     cv.folds = cv.folds)

      # Save LightGBM mean model (m.txt) to disk
      # NOTE: mmod$best_iter and mmod$best_score custom attributes are not retained in save
      lightgbm::lgb.save(booster = mmod, filename = file.path(path, paste0(y, "_m.txt")))

      # Record the conditional mean/expectation results in metadata
      # This is necessary because lgb.save does not retain them in .txt model object
      # yvalid <- c(yvalid, mmod$best_score)
      # yiters <- c(yiters, mmod$best_iter)

      #-----

      # Predict conditional mean/probabilities for the training observations
      if (block | type == "continuous") {
        mc <- predict(object = mmod, data = dmat[pi, xv], reshape = TRUE)
        if (!is.matrix(mc)) mc <- matrix(mc)
        colnames(mc) <- paste0(y, "_m", 1:ncol(mc))
      }

      rm(mmod)

      #----

      # Fit quantile models
      if (type == "continuous") {

        # Fit quantile models for each value in 'ptiles'
        qmods <- vector(mode = "list", length = length(ptiles))
        names(qmods) <- paste0("q", 1:length(ptiles))

        for (k in seq_along(ptiles)) {
          qmods[[k]] <- fitLGB(dfull = dfull,
                               dtrain = dtrain,
                               dvalid = dvalid,
                               hyper.grid = hyper.grid,
                               params.obj = list(objective = "quantile", metric = "quantile", alpha = ptiles[k]),
                               cv.folds = cv.folds)
        }

        # Predict conditional quantiles for full dataset
        qc <- matrix(data = NA, nrow = sum(pi), ncol = length(ptiles))
        for (k in 1:length(ptiles)) qc[, k] <- predict(object = qmods[[k]], data = dmat[pi, xv])
        colnames(qc) <- paste0(y, "_", names(qmods))

        # Save LightGBM quantile models (q**.txt) to disk
        for (k in seq_along(ptiles)) lightgbm::lgb.save(booster = qmods[[k]], filename = file.path(path, paste0(colnames(qc)[k], ".txt")))
        rm(qmods)

      }

      #----

      # Construct objects needed for subsequent nearest-neighbor and clustering operations
      # Conditional expectations for each training observation
      # Note that 'zc' only exists if 'y' is part of a block; NULL otherwise
      if (block | type == "continuous") {
        cd <- cbind(cd, data.table(zc, mc, qc)) # Add conditonal expectations to retained 'cd' object
        yi <- cbind(yi, Y[pi])  # Observed values associated with 'cd'
        wi <- cbind(wi, W.int[pi]) # Observed weights associated with 'cd'
      }

    } # Done processing 'y'; loop repeats with remaining variables in 'v' (if any)

    #-----

    if (!is.null(cd)) {

      # Center and scale the conditional expectations
      # Scaling is only applied to conditional expectations that are not probabilities
      # The simple check below for values between 0 and 1 is not entirely safe but unlikely to cause problems in practice
      for (j in 1:ncol(cd)) {
        x <- cd[[j]]
        if (all(x >= 0 & x <= 1)) {  # Checks if the column looks to be conditional probabilities (no scaling applied)
          ncenter <- NA
          nscale <- NA
        } else {
          ncenter <- median(x)
          nscale <- mad(x)
          if (nscale == 0) nscale <- sd(x) # Use sd() if mad() returns zero
          if (nscale == 0) nscale <- 1  # If scale is still zero, set to 1 so that normalized values will all be 0.5 (prevents errors)
          eps <- 0.001
          x <- (x - ncenter) / nscale
          x <- (x - qnorm(eps)) / (2 * qnorm(1 - eps))
          set(cd, i = NULL, j = j, value = x)
        }
        ycenter <- c(ycenter, ncenter)
        yscale <- c(yscale, nscale)
      }
      names(ycenter) <- names(yscale) <- names(cd)

      #-----

      # Potentially partition the training observations into 'nclusters' clusters via kmeans
      if (nclusters == 0) nclusters <- Inf
      ncd <- data.table::uniqueN(cd)
      nclusters <- min(nclusters, ncd)

      if (ncd > nclusters) {

        # Distance of each row in 'cd' from the approximate median point
        temp <- as.matrix(cd)
        for (j in 1:ncol(temp)) temp[, j] <- (temp[, j] - median(temp[, j])) ^ 2
        temp <- sqrt(rowSums(temp))

        # Initial cluster centers for k-means
        # This allows cluster selection to be deterministic and initial centers evenly spread through the distribution
        qt <- quantile(temp, probs = seq(from = 0, to = 1, length.out = nclusters), type = 1)
        ind <- unique(match(qt, temp))

        # Perform k-means to identify optimal cluster centers
        km <- stats::kmeans(x = cd, centers = cd[ind, ], iter.max = 30)
        kcenters <- kc <- km$centers

      } else {
        kcenters <- kc <- as.matrix(cd)
      }

      #-----

      # Find the indices of nearest neighbors in

      if (block) {

        # If block = TRUE, then we simply retaining a fixed number (30) of nearest neighbors
        # For each cluster center, find the 30 nearest neighbors in 'd' (approximate match; eps = 0.1)
        nn <- RANN::nn2(data = cd, query = kcenters, k = 30, eps = 0.1)
        nn <- nn$nn.idx

      } else {

        # If block = FALSE, then we constructing conditional expectation neighborhood with up to max(krange) neighbors
        # We are dealing with a single continuous fusion variable in this case

        # For each cluster center, find the max(krange) nearest neighbors in 'd' (approximate match; eps = 0.1)
        nn <- RANN::nn2(data = cd, query = kcenters, k = min(max(krange), nrow(cd)), eps = 0.1)
        nn <- nn$nn.idx

        # Convert the 'kc' cluster centers back to original response units
        # Original units are necessary for computing the objective function below
        for (j in 1:ncol(kc)) kc[, j] <- denormalize(z = kc[, j], center = ycenter[j], scale = yscale[j], eps = 0.001)

        # Create 'm' and 'w' matrices with the nearest-neighbor values and weights
        m <- nn; m[] <- yi[m]  # Neighbor values matrix
        w <- nn; w[] <- wi[w]  # Neighbor weights matrix

        # Cumulative weighted means
        denom <- matrixStats::rowCumsums(w)
        m1 <- matrixStats::rowCumsums(m * w) / denom

        # Approximate conditional standard deviation, based on conditional quantiles
        sdc <- abs(as.vector(matrixStats::rowDiffs(kc[, c(2, ncol(kc))]))) / diff(qnorm(range(ptiles)))  # Uses the max ptiles as the moment (need to test with larger number of ptiles)
        z <- (m1 - kc[, 1]) / sdc  # Z-score using the assumed SD of the conditional distribution under assumption of normality

        # Mean expectation error
        merr <- 1 - dnorm(z) / dnorm(0)

        # Quantile expectation error
        qerr <- lapply(2:ncol(kc), function(i) {
          tau <- ptiles[i - 1]  # Target percentile
          emax <- ifelse(tau > 0.5, tau, 1 - tau)  # Maximum possible error
          p <- matrixStats::rowCumsums((m <= kc[, i]) * w) / denom
          abs(p - tau) / emax
        })

        # Visualize the maximum error as function of 'tau' (percentile)
        # tau <- seq(0, 1, length.out = 100)
        # plot(tau, ifelse(tau > 0.5, tau, 1 - tau))
        # plot(tau, pmax(tau, 1 - tau))  # Identical

        #---

        # Average error (mean and quantile) across the conditional values
        # Note that this is effectively a weighted mean using the weights applied previous via 'colweight' vector
        err <- merr + Reduce("+", qerr)

        # Moving average of the 'err' values
        # This smooths the noise when the number of neighbors is low
        # Lagged, 5-neighbor window
        err <- sapply(5:ncol(err), function(j) matrixStats::rowMeans2(err, cols = (j - 4):j))

        # Set full columns of 'err' to Inf to ensure min(krange) number of selected neighbors is respected
        if (min(krange) > 5) err[, 1:(min(krange) - 5)] <- Inf

        # Each row's minimum error column, prior to enforcing 'tol'
        b0 <- max.col(-err, ties.method = "first")

        # Enforce absolute and relative error tolerance
        # Any error within 'tol' of the minimum is considering identical to the minimum and all are set to zero
        # This reduces the number of neighbor indices that need to be retained without affecting results materially
        kerror <- matrixStats::rowMins(err, na.rm = TRUE)
        err[err <= pmax(kerror, 0.01)] <- 0  # Absolute tolerance
        err[err <= pmax(kerror * 1.05, kerror + 0.005)] <- 0  # Relative tolerance

        # Error-minimizing column index after absolute tolerance fix
        # Can compare to 'b0' to see how 'tol' affects number of neighbors
        kbest <- max.col(-err, ties.method = "first")

        # Final "best k"
        # Add 4 because the first moving average window (i.e column 1 in 'err') includes the 5 nearest neighbors
        kbest <- kbest + 4L

        # Manual check of results
        # r <- 1724  # Row number
        # vj <- m[r, 1:kbest[r]]
        # c(mean(vj), quantile(vj, probs = ptiles))
        # kc[r, ]
        # plot(density(vj, from = 0))
        # abline(v = kc[r,], col = 2)

        #---

        # Set neighbors beyond the error-minimizing k to NA
        # Reduces size on disk when 'nn' matrix is saved
        for (j in 1:ncol(err)) nn[, j + 4] <- replace(nn[, j + 4], kbest < (j + 4), NA)

        # This operation can be used in fuse() to recreate 'kbest' quickly from 'nn'
        # test <- matrixStats::rowCounts(is.na(nn), value = FALSE)
        # all.equal(test, kbest)

        # Drop any columns in 'nn' that are only NA's
        # Not necessary to retain on disk
        keep <- matrixStats::colSums2(nn, na.rm = TRUE) > 0
        if (any(!keep)) nn <- nn[, keep]

        # Convert the index values in 'nn' to "complete data" indices
        # This allows fuse() to lookup response values using the full dataset, rather than worry about non-zero subsets
        # Note that this is unnecessary in blocked case, b/c 'pi' necessarily contains all indices (no zero model)
        nn[] <- which(pi)[nn]

        #---

        # Calculate the R-squared value between each cluster's mean expectation and the mean value of the 'kbest' neighbors
        # These values should be pretty close overall; report result to console
        mb <- matrixStats::rowCollapse(m1, kbest)  # Extract the neighborhood mean associated with 'best' k
        ssr <- sum((kc[, 1] - mb) ^ 2)
        sst <- sum((kc[, 1] - mean(kc[, 1])) ^ 2)
        r2 <- 1 - ssr / sst  # r-squared
        if (verbose) cat("-- R-squared of cluster means:", round(r2, 3), "\n")
        #plot(kc[, 1], mb); abline(0, 1, col = 2)

        # Report to console info about the preferred number of neighbors across clusters
        if (verbose) {
          cat("-- Number of neighbors in each cluster:\n")
          print(summary(kbest))
        }

      }

      #-----

      # If one wanted to smooth the neighborhood values...
      # Play with kernel density
      # x <- na.omit(m[2196, ])  # Neighborhood values
      # xw <- na.omit(w[2196, ])
      # ds <- density(x, weights = xw / sum(xw), from = 0, to = 1.2 * max(x))
      # plot(ds)
      # weighted.mean(x, xw)
      # weighted.mean(ds$x, ds$y)
      # quantile(x, probs = ptiles)
      # Use weighted quantile function to see how density matches 'kc'
      # Then could optimize bandwith to minimize quantile error
      # The issue is probably how to estimate it quickly; maybe a manual calculation
      # See here: https://rpubs.com/mcocam12/KDF_byHand


    }

    #---

    # Assemble output list
    out <- list(xpreds = xv,
                ycenter = ycenter,
                yscale = yscale,
                kcenters = kcenters,
                kneighbors = nn,
                cluster_mean_r2 = r2)

    return(out)

  }

  #-----

  # Apply buildFun() to each index in 'ylist', using forked parallel processing or serial (depending on 'fork' variable)
  # NOTE: pblapply() was imposed significant overhead, so using straight mclapply for the time being
  if (fork) {
    cat("Processing ", length(pfixes), " training steps in parallel via forking (", cores, " cores)", "\n", sep = "")
    out <- parallel::mclapply(X = 1:length(ylist),
                              FUN = buildFun,
                              mc.cores = cores,
                              mc.preschedule = FALSE,
                              verbose = FALSE)
  } else {
    if (cores > 1) cat("Using OpenMP multithreading within LightGBM (", cores, " cores)", "\n", sep = "")
    out <- lapply(X = 1:length(ylist),
                  FUN = buildFun,
                  verbose = TRUE)
  }

  #-----

  # Assemble metadata
  metadata <- list(
    xclass = xclass,
    xlevels = xlevels,
    yclass = yclass,
    ylevels = ylevels,
    ylist = ylist,
    ytype = ytype,
    ptiles = ptiles,
    dnames = colnames(dmat),
    nobs = nrow(dmat),
    timing = difftime(Sys.time(), t0),
    version = list(fusionModel = utils::packageVersion("fusionModel"), R = getRversion()),
    call = match.call.defaults()
  )

  # Add the metadata information returned by buildFun()
  metadata <- c(metadata, purrr::transpose(out))

  # Save metadata to disk
  saveRDS(metadata, file = file.path(td, "metadata.rds"))

  #-----

  # Save final model object to disk
  # Requires 'zip' package: https://cran.r-project.org/web/packages/zip/index.html
  # Zip up all of the model directories
  zip::zip(zipfile = fsn, files = list.files(td, recursive = TRUE), root = td, mode = "mirror", include_directories = TRUE)  # Zips to the temporary directory
  file.copy(from = list.files(td, "\\.fsn$", full.names = TRUE), to = fsn, overwrite = TRUE)  # Copy .zip/.fsn file to desired location
  cat("Fusion model saved to:\n", fsn, "\n")
  unlink(td)

  # Report processing time
  tout <- difftime(Sys.time(), t0)
  cat("Total processing time:", signif(as.numeric(tout), 3), attr(tout, "units"), "\n", sep = " ")

  # Return .fsn file path invisibly
  invisible(fsn)

}
