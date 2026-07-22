#' Train a Data Fusion Model
#'
#' @description
#' Trains a statistical fusion model on a "donor" dataset using sequential
#' \href{https://lightgbm.readthedocs.io/en/latest/}{LightGBM} gradient boosting
#' models to capture conditional distributions. The resulting fitted model archive
#' (\code{.fsn} file) contains the conditional expectations, candidate donor pool
#' indices, and metadata required by \code{\link{fuse}} to simulate synthetic outcomes
#' onto a "recipient" dataset.
#'
#' @param data A data frame (or \code{data.table}) containing the donor microdata.
#'   Categorical variables should be formatted as factors (ordered whenever applicable).
#' @param y Character vector or list. Variable(s) in \code{data} to be fused to a
#'   recipient dataset. Variables are modeled and fused sequentially in the order
#'   provided. If \code{y} is a list, each element can be a character vector representing
#'   a "block" of variables to be sampled jointly during fusion to preserve
#'   multivariate dependence.
#' @param x Character vector or list. Predictor variable name(s) common to both donor
#'   and recipient datasets. If a list, each element specifies the predictor set to use
#'   for the corresponding element in \code{y}. If a character vector, earlier \code{y}
#'   variables in the sequence are automatically appended as predictors for subsequent
#'   \code{y} models.
#' @param fsn Character string. File path where the trained fusion model archive will
#'   be saved. Must end with the \code{.fsn} extension. Default is \code{"fusion_model.fsn"}.
#' @param weight Character string. Name of the column in \code{data} containing survey
#'   or sampling weights. If \code{NULL} (default), uniform observation weights are assumed.
#' @param nfolds Numeric. Number of cross-validation folds used during LightGBM model
#'   tuning. If \code{0 < nfolds < 1}, it represents the proportion of observations allocated
#'   to training, with the remainder used for validation (substantially faster than
#'   full cross-validation). Default is \code{5}.
#' @param nquantiles Numeric. Number of quantile models to fit for continuous \code{y}
#'   variables in addition to the conditional mean model. Specified quantiles are
#'   evenly spaced. For example, \code{nquantiles = 2} (default) models the 25th and 75th
#'   percentiles. Even values are recommended, as the conditional mean already captures
#'   central tendency.
#' @param nclusters Numeric. Maximum number of \eqn{k}-means clusters used to group donor
#'   observations in conditional expectation space. Higher values increase donor selection
#'   precision at the expense of memory and processing time. Set to \code{0} or \code{Inf}
#'   to skip clustering (i.e., treat every donor row as a cluster center). Default is \code{2000}.
#' @param krange Numeric vector of length 2. Specifies the minimum and maximum number of
#'   nearest neighbors (\eqn{k}) evaluated when selecting optimal candidate pools for
#'   continuous conditional distributions. Default is \code{c(10, 500)}.
#' @param hyper List. LightGBM hyperparameter grid or custom values to evaluate during
#'   training. If \code{NULL} (default), a standardized baseline configuration optimized for
#'   survey fusion is used. See Details.
#' @param fork Logical. If \code{TRUE}, uses parallel processing via process forking
#'   (\code{\link[parallel]{mclapply}}) across fusion steps. Only supported on Unix/Linux/macOS platforms.
#'   Default is \code{FALSE}.
#' @param cores Integer. Number of physical CPU cores allocated to computation. When
#'   \code{fork = FALSE} or on Windows, fusion steps are processed serially, while
#'   LightGBM utilizes \code{cores} for internal OpenMP multithreading. When \code{fork = TRUE}
#'   on Unix, independent fusion steps are executed concurrently across \code{cores}.
#'   Default is \code{1}.
#'
#' @details
#' \subsection{Sequential and Block Modeling}{
#'   Data fusion proceeds sequentially through the variables specified in \code{y}. To
#'   preserve complex dependencies among tightly coupled outcomes (e.g., fuel expenditure
#'   shares), supply those variables as a character vector within a list element of \code{y}.
#'   Block variables are predicted jointly, and donor observations within blocks are sampled
#'   en masse during fusion.
#' }
#'
#' \subsection{Hyperparameter Optimization}{
#'   If a list of vectors is supplied to \code{hyper}, \code{train()} performs grid search
#'   across all parameter combinations using \eqn{V}-fold cross-validation (\code{\link[lightgbm]{lgb.cv}})
#'   with early stopping. The optimal parameter combination (minimizing loss) is selected to fit
#'   the final booster.
#'
#'   When \code{hyper = NULL}, default hyperparameters applied include:
#'   \itemize{
#'     \item \code{boosting = "gbdt"}
#'     \item \code{data_sample_strategy = "goss"}
#'     \item \code{num_leaves = 31}
#'     \item \code{feature_fraction = 0.8}
#'     \item \code{max_depth = 5}
#'     \item \code{min_data_in_leaf = max(10, round(0.001 * nrow(data)))}
#'     \item \code{num_iterations = 2500}
#'     \item \code{learning_rate = 0.1}
#'     \item \code{max_bin = 255}
#'     \item \code{min_data_in_bin = 3}
#'     \item \code{max_cat_threshold = 32}
#'   }
#' }
#'
#' \subsection{Parallel Execution & OpenMP Considerations}{
#'   On Unix-like operating systems, process forking via \code{fork = TRUE} can yield faster
#'   execution times than OpenMP multithreading. However, if OpenMP threads have already been
#'   initialized in the current R session (e.g., by \code{\link[data.table]{data.table}} or
#'   \code{\link[fst]{fst}}), forking may hang. If this occurs, execute \code{data.table::setDTthreads(1)}
#'   and \code{fst::threads_fst(1)} immediately after launching R before calling \code{train()}.
#' }
#'
#' @return Returns the file path to the saved \code{.fsn} archive invisibly.
#'
#' @references
#' Ummel, K., et al. (2024). Multidimensional well-being of US households at a fine spatial scale using fused household surveys.
#' \emph{Scientific Data}, 11(142). \doi{10.1038/s41597-023-02788-7}
#'
#' @seealso \code{\link{fuse}}, \code{\link{validate}}, \code{\link{fitLGB}}
#'
#' @examples
#' \dontrun{
#' # Load sample RECS survey dataset supplied with fusionModel
#' data(recs)

#' # Define fusion targets and common predictors
#' fusion.vars <- c("electricity", "natural_gas", "aircon")
#' predictor.vars <- names(recs)[2:12]

#' # 1. Basic model training (saves output to working directory)
#' fsn.path <- train(data = recs, y = fusion.vars, x = predictor.vars)

#' # 2. Block fusion: Preserving joint dependency across component shares
#' fusion.vars.block <- list(
#'   "electricity",
#'   "natural_gas",
#'   c("heating_share", "cooling_share", "other_share")
#' )
#' train(data = recs, y = fusion.vars.block, x = predictor.vars, fsn = "model_block.fsn")

#' # 3. Custom predictor specification per fusion step
#' xlist <- list(predictor.vars[1:4], predictor.vars[2:8], predictor.vars)
#' train(data = recs, y = fusion.vars.block, x = xlist, fsn = "model_xlist.fsn")

#' # 4. Hyperparameter override (Random Forest boosting alternative)
#' train(data = recs, y = fusion.vars, x = predictor.vars,
#'       hyper = list(boosting = "rf", feature_fraction = 0.6, max_depth = 10))

#' # 5. Grid search hyperparameter tuning
#' train(data = recs, y = fusion.vars, x = predictor.vars,
#'       hyper = list(max_depth = c(5, 10), feature_fraction = c(0.7, 0.9)))
#' }
#' @export

#---------------------

# Manual testing
# library(fusionModel)
# library(tidyverse)
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
# nfolds = 0.75
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
                  nquantiles = 2,
                  nclusters = 2000,
                  krange = c(10, 500),
                  hyper = NULL,
                  fork = FALSE,
                  cores = 1) {

  t0 <- Sys.time()

  # TO DO: Make data.table operations throughout (just applies to pre-loop checks)
  # TO DO: Make collapse operations where possible
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
    cli::cli_inform("Using specified set of predictors for each fusion variable")
  } else {
    # Sequence chaining: Append previously predicted 'y' variables as predictors for subsequent 'y' variables
    xlist <- lapply(1:length(ylist), function(i) c({if (i == 1) NULL else unlist(ylist[1:(i - 1)])}, x))
    cli::cli_inform("Using all available predictors for each fusion variable")
  }
  x <- setdiff(x, y)

  #---

  # Validate inputs against required logical constraints
  stopifnot(exprs = {
    is.data.frame(data)
    all(y %in% names(data))
    !any(c("M", "W..", "R..") %in% y)  # Reserved internal column names
    all(lengths(ylist) > 0)
    all(x %in% names(data))
    all(lengths(xlist) > 0)
    length(ylist) == length(xlist)
    length(intersect(y, x)) == 0
    is.character(fsn) & endsWith(fsn, ".fsn")
    is.null(weight) | (length(weight) == 1 & weight %in% names(data) & !weight %in% c(y, x))
    nfolds > 0
    nquantiles > 0
    nclusters >= 0
    all(krange >= 5) & length(krange) == 2
    is.null(hyper) | (is.list(hyper) & !anyDuplicated(names(hyper)))
    is.logical(fork)
    cores %in% 1:parallel::detectCores(logical = FALSE)
  })

  if (!any(x %in% xlist[[1]])) stop("There must be at least one 'x' predictor variable assigned to the first 'y' variable")

  #---

  # Calculate the percentiles to use for quantile models
  ptiles <- seq(from = 1 / nquantiles / 2, length.out = nquantiles, by = 2 * 1 / nquantiles / 2)

  # Determine if parallel forking will be used
  # Forces use of OpenMP if there are more cores than fusion steps
  fork <- fork & .Platform$OS.type == "unix" & cores > 1 & length(ylist) > 1 & cores <= length(ylist)

  # Check 'fsn' path and create parent directories, if necessary
  dir <- full.path(dirname(fsn), mustWork = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  fsn <- file.path(dir, basename(fsn))

  # Temporary directory to save intermediate outputs to
  td <- tempfile()
  dir.create(td, showWarnings = FALSE)

  #-----

  # Create and/or check observation weights
  W.lgb <- if (is.null(weight)) {
    rep(1L, nrow(data))
  } else {
    if (anyNA(data[[weight]])) stop("Missing (NA) values are not allowed in 'weight'")
    if (any(data[[weight]] < 0)) cli::cli_alert_warning("Setting negative observation weights to zero")
    set(data, j = weight, value = pmax(0, data[[weight]]))
    data[[weight]] / mean(data[[weight]]) # Scaled weights to avoid numerical issues
  }

  # Integerized version of the observation weights (requires less space when saved to .fsn object)
  W.int <- integerize(W.lgb, mincor = 0.999)

  #-----

  # Check data and variable name validity
  if (anyNA(data[y])) stop("Missing (NA) values are not allowed in 'y'")

  # Check that 'x' and 'y' contain only syntactically valid names
  bad <- setdiff(c(x, y), make.names(c(x, y)))
  if (length(bad)) stop("Fix invalid column names (see ?make.names):\n", paste(bad, collapse = ", "))

  # Check for character-type variables; stop with error if any detected
  xc <- sapply(data[c(x, y)], is.character)
  if (any(xc)) stop("Coerce character variables to factor:\n", paste(names(which(xc)), collapse = ", "))

  # 4/28/25: Removing use of checkData(); zero-variance allowed to pass through
  # Initial testing suggests LightGBM handles zero-variance cases without error

  # Determine fusion variable directory prefixes for saving to disk
  pfixes <- formatC(seq_along(ylist), width = nchar(length(ylist)), format = "d", flag = "0")

  #-----

  # Extract data classes and levels for the 'y' variables
  d <- data[y]
  yclass <- lapply(d, class)
  ylevels <- lapply(d[grepl("factor", yclass)], levels)

  # Extract data classes and levels for the 'x' variables
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

  # Output summary info to console
  cli::cli_inform(c(
    "i" = "{length(y)} fusion variable{?s}",
    "i" = "{length(x)} predictor variable{?s}",
    "i" = "{nrow(data)} observation{?s}"
  ))

  # Limit 'data' to the necessary variables
  data <- data[c(x, y)]

  # Coerce 'data' to sparse numeric matrix for use with LightGBM
  dmat <- to_mat(data)
  rm(data)

  # Write donor response/fusion variables and observation weights to disk (donor.fst)
  # Used by fuse() to select simulated donor values (excludes solo categorical response variables)
  ysave <- y[ytype == "continuous" | y %in% unlist(ylist[lengths(ylist) > 1])]
  dtemp <- as.data.frame(dmat[, ysave, drop = FALSE])
  dtemp$W <- W.int
  fst::write_fst(x = dtemp, path = file.path(td, "donor.fst"), compress = 95)
  rm(dtemp)

  #-----

  # Set the 'hyper.default' object
  if (is.null(hyper)) hyper <- list()

  # Default hyperparameter values, per LightGBM documentation
  hyper.default <- list(
    boosting = "gbdt",
    data_sample_strategy = "goss",
    num_leaves = 31,
    feature_fraction = 0.8,
    min_data_in_leaf = max(10, round(0.001 * nrow(dmat))),
    num_iterations = 2500,
    learning_rate = 0.1,
    max_depth = 5,
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
  # If forking, LightGBM uses single core internally per thread
  hyper$num_threads <- ifelse(fork, 1L, cores)

  # The 'dataset' parameters can only have a single value (not eligible for grid search)
  for (v in c("min_data_in_leaf", "max_bin", "min_data_in_bin", "max_cat_threshold")) {
    if (length(hyper[[v]]) > 1) {
      hyper[[v]] <- hyper[[v]][1]
      cli::cli_alert_warning("Only one {.val {v}} value allowed. Using: {.val {hyper[[v]]}}")
    }
  }

  # Construct LightGBM dataset parameter set
  dparams <- list(max_bin = hyper$max_bin,
                  min_data_in_bin = hyper$min_data_in_bin,
                  max_cat_threshold = hyper$max_cat_threshold,
                  min_data_in_leaf = hyper$min_data_in_leaf,
                  feature_pre_filter = TRUE)
  # Remove 'dparams' hyperparameters from 'hyper' object
  hyper[names(dparams)] <- NULL

  # Create hyperparameter grid to search
  hyper.grid <- hyper %>%
    expand.grid() %>%
    distinct() %>%
    split(seq(nrow(.)))

  #-----

  # Internal function to build LightGBM prediction model for step 'i' in 'ylist'
  buildFun <- function(i, verbose = FALSE) {

    v <- ylist[[i]]
    block <- length(v) > 1

    # Print message to console
    if (verbose) {
      cli::cli_inform("-- Training step {i} of {length(pfixes)}: {paste(v, collapse = ', ')}")
    }

    # 'y' variables from prior clusters to be included as predictors
    yv <- if (i == 1) NULL else unlist(ylist[1:(i - 1)])

    # Full set of predictor variables
    xv <- xlist[[i]]

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
      hyper.results <- list()

      #-----

      # Build zero-inflation model if continuous variable exhibits structural zeros
      if (type == "continuous" & inflated(Y)) {

        # Indices to identify which observations to use for subsequent training (ti) and prediction (pi)
        ti <- Y != 0
        if (!block) pi <- ti

        #---

        # Build LightGBM datasets for zero-inflated model
        cv.folds <- stratify(y = (Y == 0), ycont = FALSE, tfrac = nfolds, cv_list = TRUE)

        # Create full LGB training dataset with all available observations
        dfull <- lightgbm::lgb.Dataset(data = dmat[, xv, drop = FALSE],
                                       label = as.integer(Y == 0),
                                       weight = W.lgb,
                                       categorical_feature = intersect(xv, nominal),
                                       params = dparams) %>%
          lightgbm::lgb.Dataset.construct()

        # Create 'dtrain' and 'dvalid' sets if validation split requested (nfolds < 1)
        if (is.logical(cv.folds)) {
          ind <- which(cv.folds)
          dtrain <- lightgbm::lgb.Dataset(data = dmat[ind, xv, drop = FALSE],
                                          label = as.integer(Y == 0)[ind],
                                          weight = W.lgb[ind],
                                          categorical_feature = intersect(xv, nominal),
                                          params = dparams) %>%
            lightgbm::lgb.Dataset.construct()

          dvalid <- lightgbm::lgb.Dataset(data = dmat[-ind, xv, drop = FALSE],
                                          label = as.integer(Y == 0)[-ind],
                                          weight = W.lgb[-ind],
                                          categorical_feature = intersect(xv, nominal),
                                          params = dparams,
                                          reference = dtrain) %>%
            lightgbm::lgb.Dataset.construct()
        }

        #---

        # Set loss function and performance metric for binary zero-inflation model
        params.obj <- list(objective = "binary", metric = "binary_logloss")

        # Fit model
        zmod <- fitLGB(dfull = dfull,
                       dtrain = dtrain,
                       dvalid = dvalid,
                       hyper.grid = hyper.grid,
                       params.obj = params.obj,
                       cv.folds = cv.folds)

        hyper.results$z <- zmod$record_evals

        # Save LightGBM zero model (_z.txt) to disk
        lightgbm::lgb.save(booster = zmod, filename = file.path(path, paste0(y, "_z.txt")))

        # Predict conditional probability of zero for block items
        if (block) {
          zc <- matrix(predict(object = zmod, newdata = dmat[pi, xv]))
          colnames(zc) <- paste0(y, "_z")
        }

        rm(zmod)

      }

      #---

      # Build LightGBM datasets for mean and quantile models
      cv.folds <- stratify(y = Y[ti], ycont = (type == "continuous"), tfrac = nfolds, ntiles = 10, cv_list = TRUE)

      # Create full LGB training dataset with all available observations
      dfull <- lightgbm::lgb.Dataset(data = dmat[ti, xv, drop = FALSE],
                                     label = Y[ti],
                                     weight = W.lgb[ti],
                                     categorical_feature = intersect(xv, nominal),
                                     params = dparams) %>%
        lightgbm::lgb.Dataset.construct()

      # Create 'dtrain' and 'dvalid' sets if validation split requested
      if (is.logical(cv.folds)) {
        ind <- which(ti)[cv.folds]
        dtrain <- lightgbm::lgb.Dataset(data = dmat[ind, xv, drop = FALSE],
                                        label = Y[ind],
                                        weight = W.lgb[ind],
                                        categorical_feature = intersect(xv, nominal),
                                        params = dparams) %>%
          lightgbm::lgb.Dataset.construct()

        dvalid <- lightgbm::lgb.Dataset(data = dmat[-ind, xv, drop = FALSE],
                                        label = Y[-ind],
                                        weight = W.lgb[-ind],
                                        categorical_feature = intersect(xv, nominal),
                                        params = dparams,
                                        reference = dtrain) %>%
          lightgbm::lgb.Dataset.construct()
      }

      #---

      # Set loss function and metric appropriate for variable response type
      params.obj <- switch(type,
                           continuous = list(objective = "regression", metric = "l2"),
                           multiclass = list(objective = "multiclass", metric = "multi_logloss", num_class = 1L + max(Y)),
                           binary = list(objective = "binary", metric = "binary_logloss"))

      #-----

      # Fit conditional mean/probability model
      mmod <- fitLGB(dfull = dfull,
                     dtrain = dtrain,
                     dvalid = dvalid,
                     hyper.grid = hyper.grid,
                     params.obj = params.obj,
                     cv.folds = cv.folds)

      hyper.results$m <- mmod$record_evals

      # Save LightGBM mean model (_m.txt) to disk
      lightgbm::lgb.save(booster = mmod, filename = file.path(path, paste0(y, "_m.txt")))

      #-----

      # Predict conditional mean/probabilities for training observations
      if (block | type == "continuous") {
        mc <- predict(object = mmod, newdata = dmat[pi, xv, drop = FALSE])
        if (!is.matrix(mc)) mc <- matrix(mc)
        colnames(mc) <- paste0(y, "_m", 1:ncol(mc))
      }

      rm(mmod)

      #----

      # Fit quantile regression models for continuous variables
      if (type == "continuous") {

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

        for (k in seq_along(ptiles)) hyper.results[[names(qmods)[k]]] <- qmods[[k]]$record_evals

        # Predict conditional quantiles across training set
        qc <- matrix(data = NA, nrow = sum(pi), ncol = length(ptiles))
        for (k in seq_along(ptiles)) qc[, k] <- predict(object = qmods[[k]], newdata = dmat[pi, xv, drop = FALSE])
        colnames(qc) <- paste0(y, "_", names(qmods))

        # Save LightGBM quantile models (_q**.txt) to disk
        for (k in seq_along(ptiles)) lightgbm::lgb.save(booster = qmods[[k]], filename = file.path(path, paste0(colnames(qc)[k], ".txt")))
        rm(qmods)

      }

      #----

      # Append predicted conditional expectations, observed values, and weights
      if (block | type == "continuous") {
        cd <- cbind(cd, data.table(zc, mc, qc))
        yi <- cbind(yi, Y[pi])
        wi <- cbind(wi, W.int[pi])
      }

    } # Repeat for remaining variables in block 'v'

    #-----

    if (!is.null(cd)) {

      # Standardize conditional expectations (median centering + MAD scaling)
      # Skip columns representing probabilities (bounded in [0, 1])
      for (j in 1:ncol(cd)) {
        x <- cd[[j]]
        if (all(x >= 0 & x <= 1)) {
          ncenter <- NA
          nscale <- NA
        } else {
          ncenter <- median(x)
          nscale <- mad(x)
          if (nscale == 0) nscale <- sd(x)
          if (nscale == 0) nscale <- 1
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

      # Partition training observations into 'nclusters' centers via k-means in expectation space
      if (nclusters == 0) nclusters <- Inf
      ncd <- data.table::uniqueN(cd)
      nclusters <- min(nclusters, ncd)

      if (ncd > nclusters) {

        # Distance of each row in 'cd' from median point
        temp <- as.matrix(cd)
        for (j in 1:ncol(temp)) temp[, j] <- (temp[, j] - median(temp[, j])) ^ 2
        temp <- sqrt(rowSums(temp))

        # Deterministic initial cluster centers evenly spread across expectation distribution
        qt <- quantile(temp, probs = seq(from = 0, to = 1, length.out = nclusters), type = 1)
        ind <- unique(match(qt, temp))

        # Perform k-means to compute cluster centers
        km <- stats::kmeans(x = cd, centers = cd[ind, ], iter.max = 30)
        kcenters <- kc <- km$centers

      } else {
        kcenters <- kc <- as.matrix(cd)
      }

      #-----

      # Identify candidate nearest neighbors in donor pool surrounding cluster centers
      if (block) {

        # For variable blocks: retain fixed candidate pool of 30 nearest neighbors
        nn <- RANN::nn2(data = cd, query = kcenters, k = 30, eps = 0.1)
        nn <- nn$nn.idx

      } else {

        # For single continuous variables: construct optimal candidate pools up to max(krange)
        nn <- RANN::nn2(data = cd, query = kcenters, k = min(max(krange), nrow(cd)), eps = 0.1)
        nn <- nn$nn.idx

        # Convert 'kc' cluster centers back to unscaled response units for error evaluation
        for (j in 1:ncol(kc)) kc[, j] <- denormalize(z = kc[, j], center = ycenter[j], scale = yscale[j], eps = 0.001)

        # Construct neighbor value and weight matrices
        m <- nn; m[] <- yi[m]
        w <- nn; w[] <- wi[w]
        w <- w / mean(w)

        # Compute cumulative weighted candidate means
        denom <- matrixStats::rowCumsums(w)
        m1 <- matrixStats::rowCumsums(m * w) / denom

        # Estimate conditional standard deviation from outer quantiles
        sdc <- abs(as.vector(matrixStats::rowDiffs(kc[, c(2, ncol(kc))]))) / diff(qnorm(range(ptiles)))
        z <- (m1 - kc[, 1]) / sdc

        # Compute conditional expectation loss terms (mean + quantile errors)
        merr <- 1 - dnorm(z) / dnorm(0)

        qerr <- lapply(2:ncol(kc), function(i) {
          tau <- ptiles[i - 1]
          emax <- ifelse(tau > 0.5, tau, 1 - tau)
          p <- matrixStats::rowCumsums((m <= kc[, i]) * w) / denom
          abs(p - tau) / emax
        })

        # Aggregate error across mean and quantile predictions
        err <- merr + Reduce("+", qerr)

        if (ncol(err) > 5) {
          # Apply 5-neighbor moving average to smooth local neighbor noise
          err <- sapply(5:ncol(err), function(j) matrixStats::rowMeans2(err, cols = (j - 4):j))
          if (min(krange) > 5) err[, 1:(min(krange) - 5)] <- Inf
          delta <- 4L
        } else {
          delta <- 0L
        }

        # Enforce tolerance bounds around minimal error to select optimal k
        kerror <- matrixStats::rowMins(err, na.rm = TRUE)
        err[err <= pmax(kerror, 0.01)] <- 0
        err[err <= pmax(kerror * 1.05, kerror + 0.005)] <- 0

        kbest <- max.col(-err, ties.method = "first")
        kbest <- kbest + delta

        # Trim neighbor indices beyond optimal k to optimize disk space
        for (j in 1:ncol(err)) nn[, j + delta] <- replace(nn[, j + delta], kbest < (j + delta), NA)

        keep <- matrixStats::colSums2(nn, na.rm = TRUE) > 0
        if (any(!keep)) nn <- nn[, keep]

        # Convert relative neighborhood indices to absolute row indices in donor dataset
        nn[] <- which(pi)[nn]

        # Compute R-squared between cluster expectations and optimal neighborhood means
        mb <- matrixStats::rowCollapse(m1, kbest)
        ssr <- sum((kc[, 1] - mb) ^ 2)
        sst <- sum((kc[, 1] - mean(kc[, 1])) ^ 2)
        r2 <- 1 - ssr / sst
        if (verbose) cli::cli_inform("  -- R-squared of cluster means: {round(r2, 3)}")

        if (verbose) {
          cli::cli_inform("  -- Number of neighbors in each cluster:")
          print(summary(kbest))
        }

      }

    }

    #---

    # Assemble output object for current fusion step
    out <- list(xpreds = xv,
                ycenter = ycenter,
                yscale = yscale,
                kcenters = kcenters,
                kneighbors = nn,
                cluster_mean_r2 = r2,
                hyper_results = hyper.results)

    return(out)

  }

  #-----

  # Execute step-wise model building across fusion sequence (forked or serial execution)
  if (fork) {
    cli::cli_inform("Processing {length(pfixes)} training steps in parallel via forking ({cores} cores)")
    out <- parallel::mclapply(X = 1:length(ylist),
                              FUN = buildFun,
                              mc.cores = cores,
                              mc.preschedule = FALSE,
                              verbose = FALSE)
  } else {
    if (cores > 1) cli::cli_inform("Using OpenMP multithreading within LightGBM ({cores} cores)")
    out <- lapply(X = 1:length(ylist),
                  FUN = buildFun,
                  verbose = TRUE)
  }

  #-----

  # Collect package and run metadata
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
    call = match.call()
  )

  # Merge step outputs into metadata structure
  metadata <- c(metadata, purrr::transpose(out))

  # Save metadata object to temporary archive path
  saveRDS(metadata, file = file.path(td, "metadata.rds"))

  #-----

  # Compress model files into output .fsn archive
  zip::zip(zipfile = fsn, files = list.files(td, recursive = TRUE), root = td, mode = "mirror", include_directories = TRUE)
  file.copy(from = list.files(td, "\\.fsn$", full.names = TRUE), to = fsn, overwrite = TRUE)

  cli::cli_alert_success("Fusion model saved to: {.file {fsn}}")
  unlink(td)

  # Report total execution time
  tout <- difftime(Sys.time(), t0)
  cli::cli_inform("Total processing time: {signif(as.numeric(tout), 3)} {attr(tout, 'units')}")

  # Return .fsn file path invisibly
  invisible(fsn)

}
