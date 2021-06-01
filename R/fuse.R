#' Fuse variables to a recipient dataset
#'
#' @description
#' Fuse variables to a recipient dataset.
#'
#' @param data Data frame. Recipient dataset. All categorical variables should be factors and ordered whenever possible. Data types and levels are strictly validated against predictor variables defined in \code{train.object}.
#' @param train.object Output from successfull call to \link{train}.
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

#---------------------

fuse <- function(data, train.object) {

  stopifnot(exprs = {
    is.data.frame(data)
    #class(train.object) == ...
  })

  # Check that predictor variables are present
  xclass <- train.object$xclass
  xvars <- names(xclass)
  xlevels <- train.object$xlevels
  miss <- setdiff(xvars, names(data))
  if (length(miss) > 0) stop("The following predictor variables are missing from 'data':\n", paste(miss, collapse = ", "))

  # Check for appropriate class/type of predictor variables
  xtest <- lapply(data[xvars], class)
  miss <- !map2_lgl(xclass, xtest, identical)
  if (any(miss)) stop("Incompatible data type for the following predictor variables:\n", paste(names(miss)[miss], collapse = ", "))

  # Check for appropriate levels of factor predictor variables
  xtest <- lapply(data[names(xlevels)], levels)
  miss <- !map2_lgl(xlevels, xtest, identical)
  if (any(miss)) stop("Incompatible levels for the following predictor variables\n", paste(names(miss)[miss], collapse = ", "))

  # Restrict 'data' to the xvars
  data <- data[xvars]

  #-----

  # Names and order of variables to be fused
  yord <- names(train.object$models)

  # Placeholder columns for the fusion output
  data[yord] <- NA

  # Identify continuous yvars
  ycont <- names(which(sapply(train.object$yclass, function(x) x[1] %in% c("integer", "numeric"))))

  # Create the percentile values associated with the quantile function values
  if (length(ycont) > 0) {
    Qx <- dnorm(seq(-3, 3, length.out = nrow(train.object$models[[ycont[1]]]$Q) - 1))
    Qx <- c(0, cumsum(Qx / sum(Qx)))
  }

  # Assemble 'ranks' matrix for the xvars
  ranks <- lapply(xvars, matFun, data = data)
  ranks <- do.call(cbind, ranks)

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
      for (d in km) {
        ind <- match(data[[names(d)[1]]], d[[1]])
        data[names(d)[2]] <- d[ind, 2]
      }

      # Class probabilities
      p <- predict(object = m, newdata = data)

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
    data[[y]] <- S

    # Create 'ranks' matrix addition for 'y'
    yrank <- matFun(y, data)

    #-----

    # For continuous and ordered factors, adjust initial simulated values to better match known rank correlation with other variables
    if (y %in% names(train.object$ycor)) {

      # Target (rank) correlations
      rho <- train.object$ycor[[y]]

      # Restrict adjusted to non-zero observations
      i <- as.numeric(data[[y]]) != 0

      # Identify and remove any no-variance columns in 'X'
      X <- ranks[i, names(rho)]
      nv <- apply(X, MARGIN = 2, FUN = novary)
      if (any(nv)) {
        X <- X[, -which(nv)]
        rho <- rho[-which(nv)]
      }

      # Attempt to induce desired correlation
      Y <- induceCor(x = X, rho = rho, y = yrank[i, ])

      # Before and after rank correlations compared to 'rho' target correlation
      # plot(rho, cor(X, yrank[i, ])[, 1])  # Before
      # plot(rho, cor(X, Y)[, 1])  # After
      # abline(0, 1)

      # Re-order original y data to match ranks in Y (this preserve the original distribution)
      Y <- sort(data[i, y])[rank(Y, ties.method = "random")]

      # Confirm that univariate distribution is unchanged
      # hist(data[i, y])
      # hist(Y)

      # Comparing before and after y values, original scale
      # plot(data[i, y], Y)
      # abline(0, 1, col = 2)
      # cor(data[i, y], Y)

      # Update original 'y' data with adjusted simulated values
      data[i, y] <- Y

      # Update 'yrank' object to reflect adjustment
      yrank <- matFun(y, data)

    }

    #-----

    # Update the 'ranks' matrix
    ranks <- cbind(ranks, yrank)

    # Update progress bar
    setTxtProgressBar(pb, match(y, yord))

  }

  #-----

  # Simulation complete
  # Return only the fusion variables
  result <- data[yord]
  return(result)

}
