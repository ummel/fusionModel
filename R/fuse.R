#' Fuse variables to a recipient dataset
#'
#' @description
#' Fuse variables to a recipient dataset.
#'
#' @param data Data frame. Recipient dataset. All categorical variables should be factors and ordered whenever possible. Data types and levels are strictly validated against predictor variables defined in \code{train.object}.
#' @param train.object Output from successfull call to \link{train}.
#' @param induce Logical. Experimental. Should simulated values be adjusted to induce better agreement with observed rank correlations in donor? \code{induce = TRUE} can be slow for large datasets.
#' @param induce.vars Character. If \code{induce = TRUE}, an optional vector of fusion and/or predictor variables for which correlation should be induced. The default value (\code{induce.vars = NULL}) induces correlation across all variables.
#' @param use.biglm Logical. If \code{induce = TRUE}, \code{use.biglm = TRUE} will use \code{\link[biglm]{biglm}} from the \href{https://cran.r-project.org/web/packages/biglm/index.html}{biglm package} for the necessary OLS regressions. This can be faster and more memory efficiency for large datasets. The default (\code{use.biglm = FALSE}) uses \code{\link[stats]{.lm.fit}}, which is still quite fast in most cases.
#'
#' @return A data frame with same number of rows as \code{data} and one column for each synthetic fusion variable defined in \code{train.object}. The order of the columns reflects the order in which they where fused.
#' @examples
#' donor <- recs
#' recipient <- subset(recs, select = c(division, urban_rural, climate, income, age, race))
#' fusion.vars <- setdiff(names(donor), names(recipient))
#' fit <- train(data = donor, y = fusion.vars)
#' sim <- fuse(data = recipient, train.object = fit)
#' @export

#---------------------

# Manual testing

# library(fusionModel)
# source("R/utils.R")
# donor <- recs
# data <- recipient <- subset(recs, select = c(urban_rural, income, age, race, education))
# induce = TRUE
# train.object <- train(data = donor, y = setdiff(names(donor), names(recipient)))

# Big data
# data <- readRDS("~/Documents/Projects/fusionData/rec_data.rds")
# train.object <- readRDS("~/Documents/Projects/fusionData/fit.rds")
# induce <- TRUE
#
# # EXAMPLE usage
# #induce.vars <- NULL
# induce.vars <- c(yord, grep("__", xvars, fixed = TRUE, value = TRUE))

#---------------------

fuse <- function(data,
                 train.object,
                 induce = TRUE,
                 induce.vars = NULL,
                 use.biglm = FALSE) {

  stopifnot(exprs = {
    is.data.frame(data)
    #class(train.object) == ...
    is.logical(induce)
    !(!induce & !is.null(induce.vars))  # Nonsensical input
    !(!induce & use.biglm)  # Nonsensical input

  })

  # Check if 'biglm' package is required/installed
  if (use.biglm & !"biglm" %in% installed.packages()[, "Package"]) {
    stop("The 'biglm' package must be installed when 'use.biglm = TRUE'")
  }

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
  miss <- !map2_lgl(xclass, xtest, identical)
  if (any(miss)) stop("Incompatible data type for the following predictor variables:\n", paste(names(miss)[miss], collapse = ", "))

  # Check for appropriate levels of factor predictor variables
  xtest <- lapply(subset(data, select = names(xlevels)), levels)
  miss <- !map2_lgl(xlevels, xtest, identical)
  if (any(miss)) stop("Incompatible levels for the following predictor variables\n", paste(names(miss)[miss], collapse = ", "))

  #-----

  # Names and order of variables to be fused
  yord <- names(train.object$models)

  # Identify continuous yvars
  ycont <- names(which(sapply(train.object$yclass, function(x) x[1] %in% c("integer", "numeric"))))

  # Create the percentile values associated with the quantile function values
  if (length(ycont) > 0) {
    Qx <- dnorm(seq(-3, 3, length.out = nrow(train.object$models[[ycont[1]]]$Q) - 1))
    Qx <- c(0, cumsum(Qx / sum(Qx)))
  }

  #-----

  # Set default 'induce.vars' if initially NULL
  if (induce & is.null(induce.vars)) induce.vars <- c(xvars, yord)

  # Check that all 'induce.vars' are valid
  miss <- setdiff(induce.vars, c(xvars, yord))
  if (any(miss)) stop("The following 'induce.vars' are not valid:\n", paste(names(miss)[miss], collapse = ", "))

  #-----

  # Detect and impute any missing values in 'data'
  na.cols <- names(which(sapply(data, anyNA)))
  if (length(na.cols) > 0) {
    cat("Imputing missing values in predictor variables...\n")
    warning("Missing values were imputed for the following variables: ", paste(na.cols, collapse = ", "))
    for (j in na.cols) {
      x <- data[[j]]
      ind <- is.na(x)
      data.table::set(data, i = which(ind), j = j, value = imputationValue(x, ind))
    }
  }

  #-----

  # Build 'ranks' data.table for the 'xvars', if 'induce = TRUE'
  # TO DO: Restrict the variables for which this is applicable?
  # This is a memory-efficient implementation using data.table

  if (induce) {

    cat("Building ranks matrix (induce = TRUE)...\n")

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

  cat("Fusing donor variables to recipient...\n")

  # Progress bar printed to console
  pb <- txtProgressBar(min = 0, max = length(yord), style = 3)

  for (y in yord) {

    cont <- y %in% ycont

    yclass <- train.object$yclass[[y]]

    m <- train.object$models[[y]]

    #-----

    # If 'y' is continuous...
    if (cont) {

      # Vector of nodes in model 'm'
      nodes <- as.integer(colnames(m$Q))

      # Predicted node for rows in 'data'
      pnode <- predictNode(object = m, newdata = data)
      gc()

      # Catch and fix rare case of missing node (unclear why this might occur)
      miss <- setdiff(pnode, nodes)
      for (n in miss) pnode[pnode == n] <- nodes[which.min(abs(n - nodes))]
      stopifnot(all(pnode %in% nodes))

      # Placeholder vector for simulated values
      S <- vector(mode = "numeric", length = nrow(data))

      # Fit density to observations in each node
      for (n in nodes) {

        # Index identifying observations in node 'n'
        ind <- pnode == n

        if (any(ind)) {

          # Extract inputs needed for quantile function and proportion of zeros
          Q <- m$Q[, as.character(n)]

          # Simulated value
          S[ind] <- approx(x = Qx, y = Q, xout = runif(n = sum(ind)))$y

        }

      }

      # Adjust simulated values to enforce the "inner.range" constraint
      inner.range <- train.object$yinner[[y]]
      S[S > inner.range[1] & S < 0] <- inner.range[1]
      S[S > 0 & S < inner.range[2]] <- inner.range[2]

      # Ensure simulated column is correct data type
      if (yclass == "integer") S <- round(as.integer(S))

    }

    #-----

    if (!cont) {

      # Add the clustered predictors, if necessary
      km <- m$kmeans.xwalk
      if (!is.null(km)) {
        for (d in km) {
          ind <- match(data[[names(d)[1]]], d[[1]])
          data.table::set(data, j = names(d)[2], value = d[ind, 2])
        }
      }

      # Class probabilities
      p <- predict(object = m, newdata = data)
      gc()

      # Simulated value
      ptile <- runif(n = nrow(data))
      for (i in 2:ncol(p)) p[, i] <- p[, i - 1] + p[, i]
      for (i in 1:ncol(p)) p[, i] <- ptile > p[, i]
      S <- rowSums(p) + 1L
      S <- colnames(p)[S]

      # Ensure simulated vector is correct data type
      # 'S' is a character vector by default; must be coerced to factor
      if ("factor" %in% yclass) S <- factor(S, levels = train.object$ylevels[[y]], ordered = "ordered" %in% yclass)

    }

    #-----

    # Assign final simulated vector
    data.table::set(data, j = y, value = S)

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

        # Attempt to induce target rank correlations
        Yout <- induceCor(data = data.table::copy(ranks), rho = rho, y = y, scale.data = FALSE, use.biglm = use.biglm)

        # Only updated y-values if the correlation adjustment was successful (sigma2 >= 0)
        if (Yout$sigma2 >= 0) {

          # Before and after rank correlations compared to 'rho' target correlation
          # plot(rho, cor(ranks[, -..y], ranks[[y]])[, 1])  # Before
          # plot(rho, cor(ranks[, -..y], Yout$Y)[, 1])  # After
          # abline(0, 1)

          # Re-order original y data to match ranks in Y (this preserve the original distribution)
          Y <- sort(data[[y]])[data.table::frank(Yout$Y, ties.method = "random")]

          # Confirm that univariate distribution is unchanged
          # hist(data[i, y])
          # hist(Y)

          # Comparing before and after y values, original scale
          # plot(data[i, y], Y)
          # abline(0, 1, col = 2)
          # cor(data[i, y], Y)

          # Update original 'y' data with adjusted simulated values
          data.table::set(data, j = y, value = Y)

          # Update the 'ranks' matrix with ranks derived from adjusted 'Y'
          if (y %in% names(ranks)) {
            if (yclass[1] == "factor") {
              u <- levels(Y)
              newv <- paste0(y, u)
              data.table::set(ranks, j = newv, value = lapply(u, function(x) as.integer(Y == x)))
            } else {
              data.table::set(ranks, j = y, value = as.vector(scale(data.table::frank(Y, ties.method = "average"))))
            }
          }

        }

      }

    }

    #-----

    # Update for() loop progress bar
    setTxtProgressBar(pb, match(y, yord))

  }

  # Close progress bar
  close(pb)

  #-----

  # Simulation complete
  # Return only the fusion variables
  return(as.data.frame(subset(data, select = yord)))

}
