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
#' @param nfolds Numeric. Number of cross-validation folds used for LightGBM model training. Or, if \code{nfolds < 1}, the fraction of observations to use for training set; remainder used for validation (faster than cross-validation).
#' @param ptiles Numeric. One or more percentiles for which quantile models are trained for continuous \code{y} variables (along with the conditional mean).
#' @param hyper List. LightGBM hyperparameters to be used during model training. If \code{NULL}, default values are used. See Details and Examples.
#' @param cores Integer. Number of physical CPU cores used for parallel computation. If \code{cores > 1} on a Unix system, the fusion variables/blocks are processed in parallel via \code{\link[parallel]{mclapply}}. On Windows (since forking is not possible), the fusion variables/blocks are processed serially but LightGBM uses \code{cores} for internal multithreading via OpenMP.
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
#'   \item max_bin
#'   \item min_data_in_bin
#'  }
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
# source("R/utils.R")
# source("R/fitLGB.R")

# # Inputs for testing - with some modification for harder test cases
# data <- recs[1:28]
# data <- bind_rows(data, data, data, data)
# recipient <- subset(recs, select = c(weight, division, urban_rural, climate, income, age, race, hh_size, televisions))
# y = setdiff(names(data), names(recipient))
# weight <- "weight"
# x <- setdiff(names(recipient), weight)
# nfolds <- 0.75
# ptiles <- c(0.25, 0.75)
# file = "fusion_model_test.fsn"
# cores = 1
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
                  ptiles = c(0.165, 0.835),
                  hyper = NULL,
                  cores = 1) {

  t0 <- Sys.time()

  stopifnot(exprs = {
    is.data.frame(data)
    all(unlist(y) %in% names(data))
    !any(c("M", "W..", "R..") %in% unlist(y))  # Reserved variable names
    all(x %in% names(data))
    length(intersect(y, x)) == 0
    is.character(file) & endsWith(file, ".fsn")
    is.null(weight) | (length(weight) == 1 & weight %in% names(data) & !weight %in% c(y, x))
    nfolds > 0  # Not entirely safe
    is.numeric(ptiles) & all(ptiles > 0 & ptiles < 1)
    is.null(hyper) | is.list(hyper)
    cores > 0 & cores %% 1 == 0 & cores <= parallel::detectCores(logical = FALSE)
  })

  # TO DO: Make data.table operations throughout (just applies to pre-loop checks)
  if (is.data.table(data)) data <- as.data.frame(data)

  # Create correct 'y' and 'yord', depending on input type
  if (is.list(y)) {
    yord <- y
    y <- unlist(yord)
  } else {
    yord <- as.list(y)
  }

  # Determine if parallel forking will be used
  # This forces use of OpenMP if there are more cores than fusion steps
  unix <- .Platform$OS.type == "unix"
  fork <- unix & cores > 1 & length(yord) > 1 & cores <= length(yord)

  # Check 'file' path and create parent directories, if necessary
  dir <- normalizePath(dirname(file), mustWork = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  file <- file.path(dir, basename(file))

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

  # Limit 'data' to the necessary variables
  data <- data[c(xvars, yvars)]

  # Coerce 'data' to sparse numeric matrix for use with LightGBM
  dmat <- tomat(data)
  rm(data)

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
    min_data_in_bin = 3
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

  # The 'dataset' parameters 'max_bin' and 'min_data_in_bin' can only have a single value (the are not eligible to be varied within the CV routine)
  # Note that 'feature_pre_filter' is set to FALSE to allow multiple 'min_data_in_leaf' values within 'hyper'
  # https://lightgbm.readthedocs.io/en/latest/Parameters.html#dataset-parameters
  for (v in c("max_bin", "min_data_in_bin")) {
    if (length(hyper[[v]]) > 1) {
      hyper[[v]] <- hyper[[v]][1]
      cat("Only one", v, "value allowed. Using:", hyper[[v]], "\n")
    }
  }
  dparams <- list(max_bin = hyper$max_bin,
                  min_data_in_bin = hyper$min_data_in_bin,
                  feature_pre_filter = FALSE)
  hyper$max_bin <- NULL
  hyper$min_data_in_bin <- NULL

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

  # Determine response variable prefixes for saving to disk
  pfixes <- formatC(seq_along(yord), width = nchar(length(yord)), format = "d", flag = "0")

  # Placeholder lists for metadata outputs
  # lgbpred <- ycenter <- yscale <- colweight <- vector(mode = "list", length = length(yord))
  # meddist <- vector(mode = "numeric", length = length(yord))
  lgbpred <- ycenter <- yscale <- colweight <- meddist <- NULL
  #meddist <- vector(mode = "numeric", length = length(yord))

  # Temporary directory to save lightGBM models to
  td <- tempfile()
  dir.create(td, showWarnings = FALSE)

  #-----

  # Function to build LightGBM prediction model for step 'i' in 'yord'
  buildFun <- function(i, verbose = FALSE) {

    # Ensure OpenMP is not multithreading if forking
    # This should be done automatically by the packages, so this is a safety check
    if (fork) {
      threads_fst(1L)
      setDTthreads(1L)
    }

    v <- yord[[i]]
    block <- length(v) > 1

    # Print message to console
    if (verbose) cat("Training step ", i, " of ", length(pfixes), ": ", paste(v, collapse = ", "), "\n", sep = "")

    # 'y' variables from prior clusters to be included as predictors
    yv <- if (i == 1) NULL else unlist(yord[1:(i - 1)])

    # Full set of predictor variables, including 'y' from clusters earlier in sequence
    # Assign the x predictors to 'lgbpred'
    xv <- as.vector(c(xvars, yv))
    #lgbpred[[i]] <- xv
    lgbpred <- xv

    path <- file.path(td, pfixes[i])
    dir.create(path)

    cdata <- NULL

    for (y in v) {

      # Response variable type and values
      Y <- dmat[, y]
      type <- ytype[[y]]

      # Placeholders
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
                           #continuous = list(objective = "regression", metric = "huber", alpha = 1.35 * mad(Y)),  # Will throw error if mad(Y) is zero...
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
      lightgbm::lgb.save(booster = mmod, filename = file.path(path, paste0(y, "_m.txt")))

      # Predict conditional mean/probabilities for the training observations
      if (block | type == "continuous") {

        mc <- predict(object = mmod, data = dmat[pi, xv], reshape = TRUE)
        if (!is.matrix(mc)) mc <- matrix(mc)
        colnames(mc) <- paste0(y, "_m", 1:ncol(mc))

        # If 'y' is multiclass (not binary) drop the least common response value (not needed for kNN)
        if (type == "multiclass") mc <- mc[, -which.min(Matrix::colMeans(mc))]

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

      # Assemble the conditional predictions data.table ('cdata')
      # and normalize columns if 'y' is continuous
      if (block | type == "continuous") {

        d <- data.table(zc, mc, qc)

        # #---TEST----
        # # How well do different 'k' capture the implied conditional distribution?
        # p0 <- copy(d)
        # #---END TEST----

        if (type == "continuous") {
          cnorm <- c(colnames(mc), colnames(qc))
          ncenter <- sapply(d[, cnorm, with = FALSE], median)
          nscale <- sapply(d[, cnorm, with = FALSE], mad)
          nscale[nscale == 0] <- sapply(d[, cnorm[nscale == 0], with = FALSE], sd) # Use sd() if mad() returns zero
          nscale[nscale == 0] <- 1  # If scale is still zero, set to 1 so that normalized values will all be 0.5 (prevents errors)
          for (j in cnorm) {
            val <- normalize(d[[j]], center = ncenter[[j]], scale = nscale[[j]])
            set(d, i = NULL, j = j, value = val)
          }
          #ycenter[[i]] <- c(ycenter[[i]], unlist(ncenter))
          #yscale[[i]] <- c(yscale[[i]], unlist(nscale))
          ycenter <- c(ycenter, unlist(ncenter))
          yscale <- c(yscale, unlist(nscale))
        }

        # Add column weights to be kept in metadata
        # Use square root of the raw weights
        wtemp <- rep((1 / length(v)) / ncol(d), ncol(d))  # Raw weights
        wtemp <- sqrt(wtemp)
        names(wtemp) <- names(d)
        #colweight[[i]] <- c(colweight[[i]], wtemp)
        colweight <- c(colweight, wtemp)

        # Apply the column weights to 'd'
        for (j in names(d)) set(d, i = NULL, j = j, value = d[[j]] * wtemp[[j]])

        # #---TEST----
        # # How well do different 'k' capture the implied conditional distribution?
        # dd <- unique(d)
        # test <- RANN::nn2(dd, dd, k = 200 + 1)
        # test <- test$nn.idx[, -1]  # Remove the exact match NN
        #
        # res <- 20  # Number of k-values to check
        # kseq <- unique(ceiling(seq(from = ncol(test) / res, to = ncol(test), length.out = res)))
        # #o <- 2  # Obs to check (manual)
        # out <- lapply(1:500, function(o) {
        #   check <- sapply(kseq, function(K) {
        #     # Check first 'K' neighbors
        #     j <- test[o, 1:K]
        #     jv <- dmat[j, y]  # The neighbors' response value (numeric, in this case)
        #
        #     # Scaled
        #     cnt <- as.numeric(p0[o, 1]) # The mean
        #     scl <- diff(range(p0[o, -1]))  # Inter-ptile range (?)
        #     #scl <- scl / (qnorm(max(ptiles)) - qnorm(min(ptiles)))  # NO DIFFERENT that above for opt K; Estimated sd assumed N(), using min/max quantile values
        #     pv <- (p0[o, ] - cnt) / scl
        #     jv <- (jv - cnt) / scl
        #     kv <- c(weighted.mean(jv, W.lgb[j]), weightedQuantile(jv, w = W.lgb[j], p = ptiles))
        #     score <- rowMeans(abs(kv - pv))
        #
        #     # Not scaled
        #     # pv <- p0[o, ]  # NOT scaled
        #     # kv <- c(weighted.mean(jv, W.lgb[j]), weightedQuantile(jv, w = W.lgb[j], p = ptiles))
        #     # score <- rowMeans(abs((kv - pv) / pv))
        #
        #     pv; kv # Manual compare
        #     score
        #   })
        #   #plot(kseq, check, type = "b")
        #   check
        # })
        #
        # # Extract the optimal 'k' for each observation
        # out <- sapply(out, function(x) kseq[which.min(x)])
        # #---END TEST----

        cdata <- cbind(cdata, d)
        rm(d)
      }

    } # Done processing 'y'; loop repeats with remaining variables in 'v' (if any)

    #-----

    # Once all of the 'y' are processed, write the 'donor' data frame to disk, if necessary
    if (!is.null(cdata)) {

      # Calculate and retain the median distance between all donor observations
      # This is used to scale the search distance used by nn2() within fuse()
      # Take a random sample to prevent excessive compute time with large datasets
      ind <- sample.int(nrow(cdata), min(nrow(cdata), 5e3))
      cdist <- as.numeric(dist(cdata[ind, ]))
      #meddist[i] <- median(cdist)
      meddist <- median(cdist)

      # Reduce precision of donor conditional values for better on-disk compression
      for (j in names(cdata)) set(cdata, i = NULL, j = j, value = signif(round(cdata[[j]], 3), 3))

      # Add the original/observed response values; converted to integer whenever possible
      for (j in v) set(cdata, i = NULL, j = j, value = if ("numeric" %in% yclass[[j]]) dmat[pi, j] else as.integer(dmat[pi, j]))

      # Add integerized sample weight for each observation
      #setkey(cdata)  # Only relevant for weight collapse below (not used)
      set(cdata, i = NULL, j = "W..", value = W.int[pi])

      # Add original row number for each observation, if necessary
      # This is only necessary if 'pi' has any FALSE values; only applicable for single continuous variable that is zero-inflated
      # Required within 'fuse' when 'ignore_self = TRUE'
      if (any(!pi)) set(cdata, i = NULL, j = "R..", value = which(pi))

      # Collapse 'cdata' to unique observations, summing the weight ("W..") column
      # NOT USED: not compatible with 'ignore_self' option in fuse()
      #cdata <- cdata[, .(W.. = sum(W..)), by = key(cdata)]

      # Write to disk (donor.fst)
      fst::write_fst(x = cdata, path = file.path(path, "donor.fst"), compress = 100)
      rm(cdata)

    }

    out <- list(lgbpred = lgbpred, ycenter = ycenter, yscale = yscale, colweight = colweight, meddist = meddist)
    return(out)

  }

  #-----

  # Apply buildFun() to each index in 'yord', using forked parallel processing or serial (depending on 'fork' variable)
  # NOTE: pblapply() was imposed significant overhead, so using straight mclapply for the time being
  if (fork) {
    cat("Processing ", length(pfixes), " training steps in parallel (", cores, " cores)...", "\n", sep = "")
    out <- parallel::mclapply(X = 1:length(yord),
                              FUN = buildFun,
                              mc.cores = cores,
                              mc.preschedule = FALSE,
                              verbose = FALSE)
  } else {
    out <- lapply(X = 1:length(yord), buildFun, verbose = TRUE)
  }

  #-----

  # Assemble metadata
  metadata <- list(
    xclass = xclass,
    xlevels = xlevels,
    yclass = yclass,
    ylevels = ylevels,
    yorder = yord,
    ytype = ytype,
    ptiles = ptiles,
    dnames = colnames(dmat),
    timing = difftime(Sys.time(), t0),
    version = list(fusionModel = packageVersion("fusionModel"), R = getRversion()),
    call = match.call()
  )

  # Add the metadata information returned by buildFun()
  metadata <- c(metadata, purrr::transpose(out))

  # Save metadata to disk
  saveRDS(metadata, file = file.path(td, "metadata.rds"))

  #-----

  # Save final model object to disk
  # !!!NOTE Requires 'zip' package: https://cran.r-project.org/web/packages/zip/index.html
  # Zip up all of the models directories
  zip::zip(zipfile = file, files = list.files(td, recursive = TRUE), root = td, mode = "mirror", include_directories = TRUE)  # Zips to the temporary directory
  file.copy(from = list.files(td, "\\.fsn$", full.names = TRUE), to = file, overwrite = TRUE)  # Copy .zip/.fsn file to desired location
  cat("Fusion model saved to:", file, "\n")
  unlink(td)

  # Report processing time
  tout <- difftime(Sys.time(), t0)
  cat("Total processing time:", signif(as.numeric(tout), 3), attr(tout, "units"), "\n", sep = " ")

  # Return .fsn file path invisibly
  invisible(file)

}
