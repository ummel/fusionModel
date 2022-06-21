#' Fuse variables to a recipient dataset using CART fusion model
#'
#' @description
#' Fuse variables to a recipient dataset using a fusion model object produced by \code{train()}. \code{fuseM()} provides a convenience wrapper for generating multiple implicates.
#'
#' @param data Data frame. Recipient dataset. All categorical variables should be factors and ordered whenever possible. Data types and levels are strictly validated against predictor variables defined in \code{train.object}.
#' @param train.object Output from a successful call to \link{train}.
#' @param induce Logical. Experimental. Should simulated values be adjusted to induce better agreement with observed rank correlations in donor? Warning: \code{induce = TRUE} can be slow for large datasets.
#' @param induce.ignore Character. If \code{induce = TRUE}, an optional vector of fusion and/or predictor variables for which correlation should NOT be induced. Can include \link[base:regex]{regular expressions}. The default value (\code{induce.ignore = NULL}) induces correlation across all variables.
#' @param verbose Logical. Should updates be printed to console?
#' @param ... Arguments passed to \code{fuse()}.

#' @param M Integer. Number of implicates to simulate.
#' @param cores Integer. Number of cores used for parallel operations.
#'
#' @return For \code{fuse()}, a data frame with same number of rows as \code{data} and one column for each synthetic fusion variable. The order of the columns reflects the order in which they where fused.
#' @return For \code{fuseM()}, a data frame with number of rows equal to \code{M * nrow(data)}. Integer column ".M" indicates implicate assignment of each observation. Note that the ordering of recipient observations is consistent within implicates, so do not change the row order if using with \code{analyze()}.
#'
#' @examples
#' # Build a fusion model using RECS microdata
#' ?recs
#' fusion.vars <- c("electricity", "natural_gas", "aircon")
#' predictor.vars <- names(recs)[2:12]
#' fit <- train(data = recs, y = fusion.vars, x = predictor.vars)
#'
#' # Generate single implicate of synthetic 'fusion.vars',
#' #  using original RECS data as the recipient
#' recipient <- recs[predictor.vars]
#' sim <- fuse(data = recipient, train.object = fit)
#' head(sim)
#'
#' # Generate multiple implicates
#' sim <- fuseM(data = recipient, train.object = fit, M = 5)
#' head(sim)
#' table(sim$.M)
#' @export

#---------------------

# Manual testing

# library(fusionModel)
# source("R/utils.R")

# Example inputs
# donor <- recs
# data <- subset(recs, select = c(division, urban_rural, climate, income, age, race))
# induce = FALSE
# induce.ignore = NULL
# fusion.vars <- setdiff(names(donor), names(recipient))
# train.object <- train(data = donor, y = fusion.vars)

# data = readRDS("~/Documents/Projects/fusionData/recs_recipient.rds")
# train.object <- readRDS("~/Documents/Projects/fusionData/recs_fit.rds")

#---------------------

fuseCART <- function(data,
                     train.object,
                     induce = FALSE,
                     induce.ignore = NULL) {

  stopifnot(exprs = {
    is.data.frame(data)
    is.logical(induce)
    !(!induce & !is.null(induce.ignore))  # Nonsensical input
  })

  verbose <- TRUE

  #-----

  # Coerce 'data' to data.table, if necessary
  data <- data.table::as.data.table(data)

  # Check that predictor variables are present
  xclass <- train.object$xclass
  xvars <- names(xclass)
  xlevels <- train.object$xlevels
  miss <- setdiff(xvars, names(data))
  if (length(miss) > 0) stop("The following predictor variables are missing from 'data':\n", paste(miss, collapse = ", "))

  # Restrict 'data' to the xvars and ensure correct ordering of columns consistent with names(xclass)
  data <- subset(data, select = xvars)

  #-----

  # Check for appropriate class/type of predictor variables

  xtest <- lapply(data, class)
  miss <- !map2_lgl(xclass, xtest, sameClass)
  if (any(miss)) stop("Incompatible data type for the following predictor variables:\n", paste(names(miss)[miss], collapse = ", "))

  # Check for appropriate levels of factor predictor variables
  xtest <- lapply(subset(data, select = names(xlevels)), levels)
  miss <- !map2_lgl(xlevels, xtest, identical)
  if (any(miss)) stop("Incompatible levels for the following predictor variables\n", paste(names(miss)[miss], collapse = ", "))

  #-----

  # Names and order of modeled variables to be fused
  yord <- c(names(train.object$models), names(train.object$derivative$model))

  # Identify continuous yvars
  ycont <- names(which(sapply(train.object$yclass, function(x) x[1] %in% c("integer", "numeric"))))

  # Create the percentile values associated with the quantile function values
  # Note: Only really necessary if using smoothing
  if (length(ycont) > 0) {
    Qx <- map(train.object$models[ycont], ~ nrow(.x$Q)) %>% compact()
    Qx <- if (length(Qx) > 0) {
      Qx <- dnorm(seq(-3, 3, length.out = Qx[[1]] - 1))
      Qx <- c(0, cumsum(Qx / sum(Qx)))
    }
  }

  #-----

  # Set 'induce.vars'; these are the variables for which correlation will be induced
  if (induce) {
    induce.vars <- if (is.null(induce.ignore)) {
      c(xvars, yord)
    } else {
      validNames(induce.ignore, c(xvars, yord), exclude = TRUE)
    }
    stopifnot(length(induce.vars) > 0)
    if (verbose) cat("Will attempt to induce correlation for a total of", length(induce.vars), "variables\n")
  }

  #-----

  # Detect and impute any missing values in 'data'
  na.cols <- names(which(sapply(data, anyNA)))
  if (length(na.cols) > 0) {
    if (verbose) cat("Missing values imputed for the following predictors:\n", paste(na.cols, collapse = ", "), "\n")
    for (j in na.cols) {
      x <- data[[j]]
      ind <- is.na(x)
      data.table::set(data, i = which(ind), j = j, value = imputationValue(x, na.ind = ind))
    }
  }

  #-----

  # Build 'ranks' data.table for the 'xvars', if 'induce = TRUE'
  # TO DO: Restrict the variables for which this is applicable?
  # This is a memory-efficient implementation using data.table

  if (induce) {

    if (verbose) cat("Building ranks matrix (induce = TRUE)...\n")

    # Correlation variables to retain in initial 'ranks' data.table, based on 'induce.vars' argument
    retain <- intersect(induce.vars, xvars)

    # Unordered factor variables among retained 'xvars'
    xunordered <- sapply(xclass[retain], function(x) x[1] == "factor")

    # Build 'ranks' data.table for 'xvars' that are NOT unordered factors
    ranks <- subset(data, select = names(which(!xunordered)))
    for (v in names(ranks)) data.table::set(ranks, j = v, value = data.table::frank(ranks[[v]], ties.method = "average"))

    # Create dummy variable columns in 'ranks' for the 'xvars' that ARE unordered factors
    for (v in names(which(xunordered))) {
      dt <- subset(data, select = v)
      u <- xlevels[[v]]
      newv <- paste0(v, u)
      data.table::set(ranks, j = newv, value = lapply(u, function(x) as.integer(dt == x)))
    }

    # Scale all variable ranks for unit variance and zero mean
    # Makes computations simpler in induceCor()
    for (v in names(ranks)) data.table::set(ranks, j = v, value = as.vector(scale(ranks[[v]])))

    # Clean up
    rm(dt)
    gc()

  }

  #-----

  # Extract names of the numeric predictor variables for which binary versions are required

  # Variables for which binary versions are used in models
  xbin.model <- train.object$models %>%
    map("xbin") %>%
    unlist(use.names = FALSE)

  # Variables for which binary versions are used in merges
  temp <- train.object$derivative$merge %>%
    map(names) %>%
    unlist(use.names = FALSE)
  xbin.merge <- grep("_BINARY_", temp, value = TRUE)
  xbin.merge <- gsub("_BINARY_", "", xbin.merge)

  # Full vector of variables for which binary versions are required
  xbin <- unique(c(xbin.model, xbin.merge))

  # Create binary versions of any 'xbin' among the initial predictor variables
  vbin <- intersect(xbin, names(data))
  for (v in vbin) data.table::set(data, j = paste0(v, "_BINARY_"), value = data[[v]] == 0)

  #-------------------------
  #-------------------------

  if (verbose) cat("Fusing donor variables to recipient...\n")

  # Progress bar printed to console
  if (verbose) pb <- pbapply::timerProgressBar(min = 0, max = length(yord), char = "+", width = 50, style = 3)

  for (y in yord) {

    cont <- y %in% ycont

    yclass <- train.object$yclass[[y]]

    m <- if (y %in% names(train.object$models)) train.object$models[[y]] else train.object$derivative$model[[y]]

    #-----

    # Extract the 'xmerge' object, if available
    # This is used to merge known 'y' values up front, prior to probabilistic simulation
    # Object 'ind' gives the row numbers for which values are not already known and must be simulated
    xmerge <- m$xmerge
    if (!is.null(xmerge)) {
      data <- merge(data, xmerge, by = names(xmerge)[1], all.x = TRUE, sort = FALSE)
      ind <- which(is.na(data[[y]]))
    } else {
      ind <- 1L:nrow(data)
    }

    #-----

    # If 'y' is continuous...
    if (cont) {

      if (class(m) == "rpart") {

        smoothed <- is.null(names(m$Q))

        # Vector of nodes in model 'm'
        nodes <- as.integer(c(names(m$Q), colnames(m$Q)))

        # Predicted node for rows in 'data'
        pnode <- predictNode(object = m, newdata = data[ind, ])
        gc()

        # Catch and fix/kludge rare case of missing node (unclear why this occurs)
        # Note that the 'popsynth' package has a fix for the same issue: https://github.com/cran/synthpop/blob/master/R/functions.syn.r
        # Erroneous nodes are re-assigned to the valid node with the closest predicted 'yval'
        miss <- setdiff(pnode, nodes)
        for (n in miss) pnode[pnode == n] <- nodes[which.min(abs(m$frame[n, "yval"] - m$frame[nodes, "yval"]))]
        stopifnot(all(pnode %in% nodes))

        #---

        # Placeholder vector for simulated values
        S <- vector(mode = "numeric", length = length(pnode))

        # Fit density to observations in each node
        for (n in nodes) {

          # Index identifying observations in node 'n'
          i <- pnode == n

          if (any(i)) {

            if (smoothed) {

              # Extract inputs needed for quantile function and proportion of zeros
              Q <- m$Q[, as.character(n)]

              # Randomly simulate values from the conditional distribution
              # Note that this is repeated as necessary to ensure 's' does not contain any values already assigned (exclusively) via 'xmerge'
              f <- approxfun(x = Qx, y = Q)
              s <- f(runif(n = sum(i)))
              while (any(s %in% xmerge[[y]])) {
                j <- s %in% xmerge[[y]]
                s[j] <- f(runif(n = sum(j)))
              }

            } else {

              # Randomly sample the observed values within the node
              # Note that this is repeated as necessary to ensure 's' does not contain any values already assigned (exclusively) via 'xmerge'
              Q <- m$Q[[as.character(n)]]
              s <- sample(x = Q[, 1], size = sum(i), replace = TRUE, prob = Q[, 2])
              while (any(s %in% xmerge[[y]])) {
                j <- s %in% xmerge[[y]]
                s[j] <- sample(x = Q[, 1], size = sum(j), replace = TRUE, prob = Q[, 2])
              }

            }

            # Assign simulated values to 'S'
            S[i] <- s

          }

        }

      } else {

        # Make predictions using linear (biglm) model in 'm'
        fobj <- formula(paste("~", as.character(m$terms)[3L]))
        newmf <- model.frame(formula = fobj, data[ind, ])
        newmm <- model.matrix(fobj, newmf)
        S <- drop(newmm %*% replace_na(coef(m), 0))

        #---

        # Adjust simulated values to enforce the "outer.range" constraint
        # This isn't strictly necessary with rpart() simulated values, because the constraint is enforced in the quantile values themselves
        # It might be relevant, however, when a linear model is used for simulation (? - unclear)
        # outer.range <- train.object$youter[[y]]
        # S <- pmin(pmax(S, outer.range[1]), outer.range[2])

        # Adjust simulated values to enforce the "inner.range" constraint
        inner.range <- train.object$yinner[[y]]
        S[S > inner.range[1] & S < 0] <- inner.range[1]
        S[S > 0 & S < inner.range[2]] <- inner.range[2]

      }

      # Ensure simulated column is correct data type
      if (yclass == "integer") S <- as.integer(round(S))

    }

    #-----

    if (!cont) {

      if (class(m) == "rpart") {

        # Add the clustered predictors, if necessary
        km <- m$kmeans.xwalk
        if (!is.null(km)) {
          for (d in km) {
            k <- match(data[[names(d)[1]]], d[[1]])
            data.table::set(data, j = names(d)[2], value = d[k, 2])
          }
        }

        # Class probabilities
        p <- predict(object = m, newdata = data[ind, ])
        gc()

        # Simulated value
        ptile <- runif(n = length(ind))
        for (i in 2:ncol(p)) p[, i] <- p[, i - 1] + p[, i]
        S <- rowSums(ptile > p) + 1L
        S <- colnames(p)[S]

      } else {

        # Make predictions using linear (biglm) model in 'm'
        # Note that 'm' predicts the integerized factor/logical values, so it must be converted to the correct label
        fobj <- formula(paste(as.character(m$terms)[-2], collapse = ""))
        newmf <- model.frame(formula = fobj, data[ind, ])
        newmm <- model.matrix(fobj, newmf)
        S <- drop(newmm %*% replace_na(coef(m), 0))
        S <- round(S)  # Round prediction to integer

        lev <- if ("factor" %in% yclass) train.object$ylevels[[y]] else c(FALSE, TRUE)  # Factor/logical labels
        S <- pmin(pmax(S, 1), length(lev))  # Ensure plausible integer values
        S <- lev[S]  # Convert to factor level

      }

      # Ensure simulated vector is correct data type
      # 'S' is a character vector by default; must be coerced to factor or logical
      if ("factor" %in% yclass) S <- factor(S, levels = train.object$ylevels[[y]], ordered = "ordered" %in% yclass)
      if ("logical" %in% yclass) S <- as.logical(S)

    }

    #-----

    # Assign simulated vector to 'data'
    data.table::set(data, i = ind, j = y, value = S)

    #-----

    # Proceed to induce correlation and/or update 'ranks' matrix, if requested
    if (induce) {

      # Create appropriate column(s) in 'ranks' for simulated 'y' values
      # NOTE that this is only done if 'y' is in 'induce.vars'; otherwise, its rank can be ignored
      if (y %in% induce.vars) {
        if (yclass[1] == "factor") {
          u <- levels(S)
          newv <- paste0(y, u)
          data.table::set(ranks, j = newv, value = lapply(u, function(x) as.integer(S == x)))
        } else {
          data.table::set(ranks, j = y, value = as.vector(scale(data.table::frank(S, ties.method = "average"))))
        }
      }

      #-----

      # For continuous and ordered factors, adjust initial simulated values to better match known rank correlation with other variables
      if (y %in% names(train.object$ycor)) {

        # Target rank correlations
        rho <- train.object$ycor[[y]]

        # Restrict target correlation to variables present in 'ranks'
        rho <- rho[names(rho) %in% names(ranks)]

        # Restrict correlation correction to non-zero observations in 'y'
        # NOTE: Setting of this option should be consistent with analogous line in train()
        #ind <- which(as.numeric(data[[y]]) != 0)

        # No restriction on correlation correction
        ind <- 1:nrow(data)

        # Attempt to induce target rank correlations
        Yout <- induceCor(data = data.table::copy(ranks[ind, ]), rho = rho, y = y, scale.data = FALSE)

        # Only updated y-values if the correlation adjustment was successful (sigma2 >= 0)
        if (Yout$sigma2 >= 0) {

          # Re-order original y data to match ranks in Y (this preserve the original distribution)
          Y <- sort(data[ind, ][[y]])[data.table::frank(Yout$Y, ties.method = "random")]

          # NOTE: The confirmation code below has not been updated in some time (probably broken)
          # Before and after rank correlations compared to 'rho' target correlation
          # plot(rho, cor(ranks[, -..y], ranks[[y]])[, 1])  # Before
          # plot(rho, cor(ranks[, -..y], Yout$Y)[, 1])  # After
          # abline(0, 1)

          # Confirm that univariate distribution is unchanged
          # hist(data[i, y])
          # hist(Y)

          # Comparing before and after y values, original scale
          # plot(data[i, y], Y)
          # abline(0, 1, col = 2)
          # cor(data[i, y], Y)

          # Update original 'y' in 'data' with adjusted simulated values
          data.table::set(data, i = ind, j = y, value = Y)

          # Update the 'ranks' matrix with ranks derived from adjusted 'Y'
          if (y %in% names(ranks)) {
            if (yclass[1] == "factor") {
              dt <- subset(data, select = y)
              u <- levels(dt)
              newv <- paste0(y, u)
              data.table::set(ranks, j = newv, value = lapply(u, function(x) as.integer(dt == x)))
              #data.table::set(ranks, j = newv, value = lapply(u, function(x) as.integer(Y == x)))  # Safe only when 'ind' is 1:nrow(data)
            } else {
              data.table::set(ranks, j = y, value = as.vector(scale(data.table::frank(data[[y]], ties.method = "average"))))
              #data.table::set(ranks, j = y, value = as.vector(scale(data.table::frank(Y, ties.method = "average"))))  # Safe only when 'ind' is 1:nrow(data)
            }
          }

        }

      }

    }

    #-----

    # If 'y' is one of the 'xbin', add a binary version of simulated values to 'data'
    if (y %in% xbin) {
      data.table::set(data, j = paste0(y, "_BINARY_"), value = data[[y]] == 0)
      xbin <- setdiff(xbin, y)
    }

    #-----

    # Update for() loop progress bar
    if (verbose) pbapply::setTimerProgressBar(pb, match(y, yord))

  }

  # Close progress bar
  if (verbose) pbapply::closepb(pb)

  #-------------------------
  #-------------------------

  # Add any constant fusion variables
  for (v in names(train.object$derivative$constant)) {
    set(data, j = v, value = train.object$derivative$constant[[v]])
  }

  # Merge categorical fusion variables that are identified by a 1-to-1 linear relationship with another categorical variable
  for (m in train.object$derivative$merge) {
    data <- merge(data, m, by = names(m)[1], all.x = TRUE, sort = FALSE)
  }

  #-----

  # Simulation complete
  # Return only the fusion variables, in the order in which they were added to 'data'
  return(data %>%
           subset(select = c(yord, setdiff(names(train.object$yclass), yord))) %>%
           as.data.frame()
  )

}
