#' Train a fusion model
#'
#' @description
#' Train a fusion model on "donor" data using sequential \href{https://lightgbm.readthedocs.io/en/latest/}{LightGBM} models to model the characteristics of conditional distributions. The resulting fusion model (.fsn file) can used with \link{fuse} to simulate outcomes for a "recipient" dataset.
#'
#' @param data Data frame. Donor dataset. Categorical variables must be factors and ordered whenever possible.
#' @param y Character or list. Variables in \code{data} to eventually fuse to a recipient dataset. Variables are fused in the order provided. If \code{y} is a list, each entry is a character vector possibly indicating multiple variables to fuse as a block.
#' @param x Character. Predictor variables in \code{data} common to donor and eventual recipient.
#' @param file Character. File where fusion model will be saved. Must use \code{.fsn} suffix.
#' @param weight Character. Name of the observation weights column in \code{data}. If NULL (default), uniform weights are assumed.
#' @param nfolds Integer. Number of cross-validation folds used for LightGBM model training.
#' @param ptiles Numeric. One or more percentiles for which quantile models are trained for continuous \code{y} variables (along with the conditional mean).
#' @param hyper List. LightGBM hyperparameters to be used during model training. If \code{NULL}, default values are used. See Details and Examples.
#' @param threads Integer. Number of threads used for LightGBM parallel operations. \code{threads = 0} will use all threads detected by OpenMP. NOTE: This may change in future.
#'
#' @details When \code{y} is a list, each slot indicates either a single variable or, alternatively, multiple variables to fuse as a block. Variables within a block are sampled jointly from the original donor data during fusion. See Examples.
#' @details The fusion model written to \code{file} is a zipped archive created by \code{\link[zip]{zip}} containing models and data required by \link{fuse}.
#'
#' @details The \code{hyper} argument can be used to specify the LightGBM hyperparameter values over which to perform a "grid search" during model training. \href{https://lightgbm.readthedocs.io/en/latest/Parameters.html}{See here} for the full list of parameters. For each combination of hyperparameters, \code{nfolds} cross-validation is performed using \code{\link[lightgbm]{lgb.cv}} with an early stopping condition. The parameter combination with the lowest loss function value is used to fit the final model via \code{\link[lightgbm]{lgb.train}}. The more candidate parameter values specified in \code{hyper}, the longer the processing time. If \code{hyper = NULL}, a single set of parameters is used LightGBM default values. Typically, users will only have reason to specify the following parameters via \code{hyper}:
#' @details \itemize{
#'   \item boosting
#'   \item num_leaves
#'   \item bagging_fraction
#'   \item feature_fraction
#'   \item min_data_in_leaf
#'   \item num_iterations
#'   \item learning_rate
#'   }
#'
#' @return A fusion model object (.fsn) is saved to \code{file}.
#'
#' @examples
#' # Build a fusion model using RECS microdata
#' # Note that "test_model.fsn" will be written to working directory
#' ?recs
#' fusion.vars <- c("electricity", "natural_gas", "aircon")
#' predictor.vars <- names(recs)[2:12]
#' train(data = recs, y = fusion.vars, x = predictor.vars, file = "test_model.fsn")
#'
#' # When 'y' is a list, it can specify variables to fuse as a block
#' fusion.vars <- list("electricity", "natural_gas", c("heating_share", "cooling_share", "other_share"))
#' fusion.vars
#' train(data = recs, y = fusion.vars, x = predictor.vars, file = "test_model.fsn")
#'
#' # Specify a single set of LightGBM hyperparameters
#' train(data = recs, y = fusion.vars, x = predictor.vars, file = "test_model.fsn",
#'       hyper = list(boosting = "goss",
#'                    feature_fraction = 0.8,
#'                    num_iterations = 300
#'       ))
#'
#' # Specify a range of LightGBM hyperparameters to search over
#' # This takes longer, because there are more models to test
#' train(data = recs, y = fusion.vars, x = predictor.vars, file = "test_model.fsn",
#'       hyper = list(num_leaves = c(10, 30),
#'                    feature_fraction = c(0.8, 0.9, 1),
#'                    num_iterations = 50
#'       ))
#' @export

#---------------------

# Manual testing
# library(fusionModel)
#
# source("R/utils.R")
# source("R/stratify.R")
# source("R/fitLGB.R")

# # Inputs for testing - with some modification for harder test cases
# data <- recs[1:28]
# recipient <- subset(recs, select = c(weight, division, urban_rural, climate, income, age, race, hh_size, televisions))
# y = setdiff(names(data), names(recipient))
# weight <- "weight"
# x <- setdiff(names(recipient), weight)
# nfolds <- 5L
# ptiles <- c(0.25, 0.75)
# file = "fusion_model_test.fsn"
# threads = 1
# hyper <- list()
#
# # Create clustering of 'y' variables
# y0 <- y
# y <- list("education", "square_feet", "natural_gas", "renter")
# y <- c(y, list(setdiff(y0, unlist(y))))
#
# # Full call
# train(data, y, x, file)
#
# # From 'fuse5.r'
# test <- fuse(data = data, fsn_file = file)

#---------------------

train <- function(data,
                  y,
                  x,
                  file = "fusion_model.fsn",
                  weight = NULL,
                  nfolds = 5,
                  ptiles = c(0.25, 0.75),
                  hyper = NULL,
                  threads = 1) {

  stopifnot(exprs = {
    is.data.frame(data)
    all(unlist(y) %in% names(data))
    all(x %in% names(data))
    is.character(file) & endsWith(file, ".fsn")
    is.null(weight) | weight %in% names(data)
    length(unique(c(x, y, weight))) == length(c(x, y, weight))
    nfolds > 0 & nfolds %% 1 == 0
    is.numeric(ptiles) & all(ptiles > 0 & ptiles < 1)
    is.null(hyper) | is.list(hyper)
    threads >= 0 & threads %% 1 == 0
  })

  if (is.null(hyper)) hyper <- list()
  if (is.data.table(data)) data <- as.data.frame(data)

  # Check 'file' path and create parent directories, if necessary
  dir <- normalizePath(dirname(file), mustWork = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  file <- file.path(dir, basename(file))

  # Create correct 'y' and 'yord', depending on input type
  if (is.list(y)) {
    yord <- y
    y <- unlist(yord)
  } else {
    yord <- as.list(y)
  }

  # Cluster number assignment for each 'y' variable
  cn <- unlist(mapply(rep, seq_along(yord), each = lengths(yord)))

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

  #-----

  # Check data and variable name validity
  if (anyNA(data[y])) stop("Missing (NA) values are not allowed in 'y'")

  # Check for character-type variables; stop with error if any detected
  xc <- sapply(data[c(x, y)], is.character)
  if (any(xc)) stop("Coerce character variables to factor:\n", paste(names(which(xc)), collapse = ", "))

  # Check that the 'xvars' and 'yvars' contain only syntactically valid names
  bad <- setdiff(c(x, y), make.names(c(x, y)))
  if (length(bad)) cat("Fix invalid column names (see ?make.names):\n", paste(bad, collapse = ", "))

  # Check for no-variance (constant) variables
  # Stop with error if any 'y' are constant; remove constant 'x' with message
  constant <- names(which(sapply(data[y], novary)))
  if (length(constant)) stop("Zero-variance 'y' variable(s) detected (remove them):\n", paste(constant, collapse = ", "), "\n")
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
      x <- data[[j]]
      ind <- is.na(x)
      data[ind, j] <-  imputationValue(x, ind)
    }
  }

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

  # Limit 'data' to the necessary variables
  data <- data[c(yvars, xvars)]

  # Coerce 'data' to sparse numeric matrix for use with LightGBM
  dmat <- tomat(data)

  # Determine response variable prefixes for saving to disk
  pfixes <- formatC(seq_along(yord), width = nchar(length(yord)), format = "d", flag = "0")

  #-----

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
    num_threads = threads  # 0 means default number of threads in OpenMP
  )

  # Use default hyperparameters, if not specified by user
  for (v in names(hyper.default)) {
    if (!v %in% names(hyper)) {
      hyper[[v]] <- hyper.default[[v]]
    }
  }

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

  # Placeholder lists for metadata outputs
  xlist.lgb <- ycenter <- yscale <- colweight <- vector(mode = "list", length = length(yord))
  meddist <- vector(mode = "numeric", length = length(yord))

  # Temporary directory to save lightGBM models to
  td <- tempfile()
  dir.create(td, showWarnings = FALSE)

  #-----

  # Sequentially build LightGBM prediction models
  cat("Building LightGBM models...\n")

  for (i in 1:length(yord)) {

    v <- yord[[i]]

    block <- length(v) > 1

    # 'y' variables from prior clusters to be included as predictors
    yv <- if (i == 1) NULL else unlist(yord[1:(i - 1)])

    # Full set of predictor variables, including 'y' from clusters earlier in sequence
    # Assign the x predictors to 'xlist'
    xv <- as.vector(c(xvars, yv))
    xlist.lgb[[i]] <- xv

    path <- file.path(td, pfixes[i])
    dir.create(path)

    cdata <- NULL

    for (y in v) {

      # Response variable type and values
      Y <- dmat[, y]
      type <- ytype[[y]]

      zc <- NULL
      mc <- NULL
      qc <- NULL
      ti <- rep(TRUE, length(Y))
      pi <- ti

      #-----

      if (type == "continuous" & inflated(Y)) {

        ti <- Y != 0
        if (!block) pi <- ti

        # Create the LGB training dataset
        dfull <- lightgbm::lgb.Dataset(data = dmat[, xv],
                                       label = as.integer(Y == 0),
                                       weight = W.lgb,
                                       categorical_feature = intersect(xv, nominal))
        lightgbm::lgb.Dataset.construct(dfull)

        # List indicating random assignment of folds
        cv.folds <- stratify(y = (Y == 0), ycont = FALSE, tfrac = nfolds, cv_list = TRUE)

        # Set the loss function and performance metric used with lightGBM
        params.obj <- list(objective = "binary", metric = "binary_logloss")

        # Fit model
        zmod <- fitLGB(data.lgb = dfull, hyper.grid = hyper.grid, params.obj = params.obj, cv.folds = cv.folds)

        # Save LightGBM mean model (m.txt) to disk
        lightgbm::lgb.save(booster = zmod, filename = file.path(path, paste0(y, "_z.txt")))

        # Conditional probability of zero for the non-zero training observations
        if (block) {
          zc <- matrix(predict(object = zmod, data = dmat[pi, xv], reshape = TRUE))
          colnames(zc) <- paste0(y, "_z")
        }

        rm(zmod)

      }

      #-----

      # Create the LGB training dataset
      dfull <- lightgbm::lgb.Dataset(data = dmat[ti, xv],
                                     label = Y[ti],
                                     weight = W.lgb[ti],
                                     categorical_feature = intersect(xv, nominal))
      lightgbm::lgb.Dataset.construct(dfull)

      # List indicating random assignment of folds
      cv.folds <- stratify(y = Y[ti], ycont = (type == "continuous"), tfrac = nfolds, ntiles = 10, cv_list = TRUE)

      # Set the loss function and performance metric used with lightGBM
      # This is here so that the Huber loss 'alpha' is set based on observed MAD
      # See full list of built-in metrics: https://lightgbm.readthedocs.io/en/latest/Parameters.html#metric-parameters
      # And here about Huber: https://stats.stackexchange.com/questions/465937/how-to-choose-delta-parameter-in-huber-loss-function
      params.obj <- switch(type,
                           #continuous = list(objective = "regression", metric = "huber", alpha = 1.35 * mad(Y)),  # Will throw error if mad(Y) is zero...
                           continuous = list(objective = "regression", metric = "l2"),
                           multiclass = list(objective = "multiclass", metric = "multi_logloss", num_class = 1L + max(Y)),
                           binary = list(objective = "binary", metric = "binary_logloss"))

      # Fit model
      mmod <- fitLGB(data.lgb = dfull, hyper.grid = hyper.grid, params.obj = params.obj, cv.folds = cv.folds)

      # Save LightGBM mean model (m.txt) to disk
      lightgbm::lgb.save(booster = mmod, filename = file.path(path, paste0(y, "_m.txt")))

      # Conditional mean/probabilities for the training observations
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
          qmods[[k]] <- fitLGB(data.lgb = dfull,
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

      # Assemble the conditional predictions data.table
      # and normalize columns if 'y' is continuous
      if (block | type == "continuous") {
        d <- data.table(zc, mc, qc)
        if (type == "continuous") {
          cnorm <- c(colnames(mc), colnames(qc))
          ncenter <- sapply(d[, cnorm, with = FALSE], median)
          nscale <- sapply(d[, cnorm, with = FALSE], mad)
          nscale[nscale == 0] <- sapply(d[, cnorm[nscale == 0], with = FALSE], sd) # Use sd() if mad() returns zero
          for (j in cnorm) {
            if (nscale[[j]] == 0) {
              set(d, i = NULL, j = j, value = NULL)  # If the scale parameter is still zero, remove the column
            } else {
              val <- normalize(d[[j]], center = ncenter[[j]], scale = nscale[[j]])
              set(d, i = NULL, j = j, value = val)
            }
          }
          ycenter[[i]] <- c(ycenter[[i]], unlist(ncenter))
          yscale[[i]] <- c(yscale[[i]], unlist(nscale))
        }

        # Add column weights to be kept in metadata
        wtemp <- rep((1 / length(v)) / ncol(d), ncol(d))
        names(wtemp) <- names(d)
        colweight[[i]] <- c(colweight[[i]], wtemp)

        cdata <- cbind(cdata, d)
      }

    }

    # Once 'y' is processed, write the 'donor' data frame to disk, if necessary
    if (!is.null(cdata)) {

      # Apply the column weights to 'cdata'
      for (j in names(cdata)) set(cdata, i = NULL, j = j, value = cdata[[j]] * colweight[[i]][[j]])

      # Calculate and retain the median distance between all donor observations
      # This is used to scale the search distance used by nn2() within fuse()
      meddist[i] <- median(as.numeric(dist(cdata)))

      # Clean final donor data, add the observed response variable values, and write to disk
      # TO DO: Make data.table operation for efficiency?
      cdata <- cdata %>%
        mutate_all(cleanNumeric, tol = 0.001) %>%
        cbind(data[pi, v, drop = FALSE])
      fst::write_fst(x = cdata, path = file.path(path, "donor.fst"), compress = 100)

    }

  }

  #-----

  # Assemble metadata and save to disk
  metadata <- list(
    xclass = xclass,
    xlevels = xlevels,
    yclass = yclass,
    ylevels = ylevels,
    yorder = yord,
    ytype = ytype,
    xlist.lgb = xlist.lgb,
    ycenter = ycenter,
    yscale = yscale,
    colweight = colweight,
    meddist = meddist,
    ptiles = ptiles,
    call = match.call()
  )
  saveRDS(metadata, file = file.path(td, "metadata.rds"))

  #-----

  # Save final model object to disk
  # !!!NOTE Requires 'zip' package: https://cran.r-project.org/web/packages/zip/index.html
  # Zip up all of the models directories
  zip::zip(zipfile = file, files = list.files(td, recursive = TRUE), root = td, mode = "mirror", include_directories = TRUE)  # Zips to the temporary directory
  file.copy(from = list.files(td, "\\.fsn$", full.names = TRUE), to = file, overwrite = TRUE)  # Copy .zip/.fsn file to desired location
  cat("Fusion model saved to:", file)
  unlink(td)
  invisible(file)

}
