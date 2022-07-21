#' Block and chain fusion variables
#'
#' @description
#' Given a set of prospective fusion variables, \code{blockchain()} derives a pseudo-optimal "chain" (fusion sequence/order) for the variables and (optionally) groups them into "blocks" for joint fusion of multiple variables at once. The output can be passed directly to the \code{y} argument of \code{\link{train}}.
#'
#' @param data Data frame. Donor dataset. All categorical variables should be factors and ordered whenever possible.
#' @param y Character or list. Variables in \code{data} to eventually fuse to a recipient dataset. If \code{y} is a list, any pre-specified blocks are preserved in the output.
#' @param x Character. Predictor variables in \code{data} common to donor and eventual recipient.
#' @param delta Numeric. Controls how aggressively variables are grouped into blocks. \code{delta = 0} results in no new blocks. See Details.
#' @param maxsize Integer. Maximum number of variables allowed in a block (excluding pre-specified blocks).
#' @param weight Character. Name of the observation weights column in \code{data}. If NULL (default), uniform weights are assumed.
#' @param nfolds Integer > 3. Number of cross-validation folds used to fit LASSO models.
#' @param criterion Character. Either \code{"min"} or \code{"1se"}. Determines what level of cross-validated model complexity is returned by \code{\link[glmnet]{cv.glmnet}}. All else equal, \code{criterion = "1se"} will lead to more aggressive blocking.
#' @param cores Integer. Number of cores used. Only applicable on Unix systems.
#'
#' @details The algorithm uses cross-validated LASSO models fit via [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html). It first builds a "complete" model for each *y* (fusion) variable using *all* other variables as predictors; the cross-validated model skill is the maximum possible. Next, a "minimal" model is fit using only the *x* predictors; the model skill is divided by the maximum to create a "score" metric. The *y* with the maximum score is assigned to the first position in the fusion chain and included as a predictor in all subsequent models. The remaining *y* are assigned to the chain in the same, greedy fashion.
#' @details When a fusion variable (*y1*) is selected for inclusion in the chain, its score is compared to that of the previous iteration; i.e. its score prior to including the preceding fusion variable (*y0*) as a predictor. If the score does not improve by at least \code{delta}, then *y1* is grouped into a block with *y0*. The general logic here is that chaining makes sense if/when it adds substantial explanatory power (i.e. when *y0* helps predict *y1*). If chaining does not appear to do this, then the default preference is to fuse the variables jointly as a block.
#' @examples
#' ?recs
#' fusion.vars <- names(recs)[13:18]
#' predictor.vars <- names(recs)[2:12]
#' yorder <- blockchain(data = recs, y = fusion.vars, x = predictor.vars)
#' yorder
#'
#' # 'y' can be a list with pre-specified blocks that are preserved in the output
#' fusion.vars <- list("electricity", "natural_gas", c("heating_share", "cooling_share", "other_share"))
#' yorder <- blockchain(data = recs, y = fusion.vars, x = predictor.vars)
#' yorder
#' @export

#---------------------

# library(fusionModel)
# source("R/utils.R")
#
# # Example data
# data <- recs[1:26]
# recipient <- subset(data, select = c(division, urban_rural, climate, income, age, race, hh_size, renter))
# x <- names(recipient)
# y <- setdiff(names(data), c(x, "weight"))
# weight <- NULL
# maxsize <- 4
# nfolds <- 10
# delta <- 0.01
# criterion <- "min"
# cores <- 2
#
# y <- as.list(y)
# y[[15]] <- unlist(y[c(8, 15:17)])
# y <- y[c(1:7,9:15)]

#-----

blockchain <- function(data,
                       y,
                       x,
                       delta = 0.01,
                       maxsize = 4,
                       weight = NULL,
                       nfolds = 10,
                       criterion = c("1se", "min"),
                       cores = 1) {

  stopifnot(exprs = {
    is.data.frame(data)
    all(unlist(y) %in% names(data))
    all(x %in% names(data))
    length(x) > 1  # glmnet requires at least two predictor variables
    delta >= 0
    maxsize >= 1 & maxsize %% 1 == 0
    is.null(weight) | weight %in% names(data)
    nfolds > 3 & nfolds %% 1 == 0  # glmnet requires nfolds by > 3 (recommends 10)
    is.character(criterion) & all(criterion %in% c("min", "1se"))
    cores >= 1 & cores %% 1 == 0
  })

  # TO DO: Make data.table operations throughout (just applies to pre-checks)
  if (is.data.table(data)) data <- as.data.frame(data)
  if (length(criterion) > 1) criterion <- criterion[1]

  if (is.list(y)) {
    input <- y
    y <- unlist(input)
  } else {
    input <- as.list(y)
  }

  W <- if (is.null(weight)) {
    rep(1L, nrow(data))
  } else {
    data[[weight]] / mean(data[[weight]])
  }

  # Check for character-type variables; stop with error if any detected
  xc <- sapply(data[c(x, y)], is.character)
  if (any(xc)) stop("Coerce character variables to factor:\n", paste(names(which(xc)), collapse = ", "))

  # Detect and impute any missing values in 'x' variables
  na.cols <- names(which(sapply(data[x], anyNA)))
  if (length(na.cols) > 0) {
    cat("Missing values imputed for the following 'x' variable(s):\n", paste(na.cols, collapse = ", "), "\n")
    for (j in na.cols) {
      xj <- data[[j]]
      ind <- is.na(xj)
      data[ind, j] <-  imputationValue(xj, ind)
    }
  }

  # Which 'y' variable are continuous?
  ycont <- names(which(sapply(data[y], is.numeric)))

  cli::cli_progress_step("Preparing data")
  d <- data[c(x, y)]
  d <- mutate_if(d, is.ordered, as.integer)
  d <- mutate_if(d, is.logical, as.integer)
  d <- mutate_if(d, is.factor, ~ factor(.x, levels = intersect(levels(.x), unique(.x))))  # This is necessary to ensure that the levels are actually present in the data
  lev <- lapply(d, levels)
  rm(data)

  # Observed class proportions/probabilities for the 'y' variables (if unordered factor; 1 otherwise)
  yweight <- lapply(y, function(v) {
    if (!is.numeric(d[[v]])) {
      xt <- xtabs(formula(paste("W~", v)), data = d)
      xt / sum(xt)
    } else {
      1
    }
  })
  names(yweight) <- y

  #-----

  # Get one-hot expanded variable names
  # Note that the 'sep' argument must match that used in one_hot()
  V <- lapply(names(lev), function(v) if (is.null(lev[[v]])) v else paste(v, lev[[v]], sep = ".."))
  names(V) <- names(lev)

  # One-hot encode unordered factors and convert to sparse matrix
  d <- one_hot(d)
  stopifnot(all(unlist(V) == colnames(d)))

  #-----

  # Names of the stable 'x' predictor variables
  X <- unlist(V[x])

  # Names of the individual 'y' columns to be modeled
  Y <- unlist(V[y])

  # Create stratified cross-validation fold assignments
  cli::cli_progress_step("Constructing cross-validation folds")
  cv.foldid <- lapply(Y, function(y) {
    stratify(d[, y], ycont = y %in% ycont, tfrac = nfolds, ntiles = 10)
  })

  #-----

  # Fit Gaussian lasso model function
  # lasso <- function(y, x) {
  #   m <- glmnet::glmnet(
  #     x = d[, x],
  #     y = d[, y],
  #     family = "gaussian",
  #     weights = w,
  #     alpha = 1)
  #   #i <- which(m$dev.ratio / max(m$dev.ratio) >= 0.95)[1] # Index of preferred lambda
  #   # plot(m$dev.ratio, type = "l")
  #   # abline(v = i, h = m$dev.ratio[i], lty = 2)
  #   # m$lambda[i]
  #   #m$dev.ratio[i] # Akin to R-squared
  #   max(m$dev.ratio)
  # }

  # Cross-validated LASSO
  lassoCV <- function(y, x, return_x = FALSE) {

    mcv <- glmnet::cv.glmnet(
      x = d[, x],
      y = d[, y],
      weights = W,
      type.measure = "mae",  # Use mean absolute error for robustness
      foldid = cv.foldid[[y]],
      family = "gaussian",
      alpha = 1)

    # Get the coefficients...


    # DEPRECATED FOR NOW -- pain to set up
    # xkeep <- if (return_x) {
    #   m <- glmnet::glmnet(
    #     x = d[, x],
    #     y = d[, y],
    #     weights = w,
    #     family = "gaussian",
    #     alpha = 1,
    #     lambda = mcv$lambda)
    #   cf <- coef(m, s = mcv$lambda.1se)
    #   nonzero <- sapply(V, function(x) any(x %in% X[cf[-1, ] != 0]))
    #   names(which(nonzero))
    # } else {
    #   NULL
    # }

    # Coefficients of preferred model
    #cf <- coef(mcv$glmnet.fit, s = mcv$lambda.1se)
    #cf <- coef(mcv$glmnet.fit, s = mcv$lambda.min)

    # weighted.r2 <- function(obs, pred, w) {
    #   res <- obs - pred
    #   SSe <- sum(w * res ^ 2)
    #   SSt <- sum(w * (obs - weighted.mean(obs, w)) ^ 2)
    #   1 - SSe / SSt
    # }
    #
    # r2 <- weighted.r2(obs = d[, y], pred = predict(object = mcv, newx = d[, x], s = "lambda.1se"), w = w)
    # return(r2)

    # Mean cross-validated error using lambda with error with 1SE of the minimum
    cvm <- mcv$cvm[which(mcv$lambda == mcv[[paste0("lambda.", criterion)]])]
    return(cvm)

  }

  # test <- lassoCV("education", X)
  # test <- lassoCV("electricity", X)

  #-----

  # # Testing -- compare gaussian to multinomial
  #
  # v <- "heat_type"
  #
  # m <- glmnet::glmnet(
  #   x = d[, X],
  #   y = recs[[v]],
  #   family = "multinomial",
  #   weights = w,
  #   alpha = 1)
  # #i <- which(m$dev.ratio / max(m$dev.ratio) >= 0.95)[1] # Index of preferred lambda
  # #m$dev.ratio[i]
  #
  # f <- function(v) {
  #   m <- sapply(V[[v]], lasso, x = X)
  #   w <- yweight[[v]]
  #   sum(m * w)
  # }
  # f(v)

  #-----

  # Placeholder for results sequence
  N <- length(input)
  ord <- vector(mode = "list", length = N)
  ratio <- vector(mode = "list", length = N)

  # Extract R2 for the "full" model that includes all possible predictors
  cli::cli_progress_step("Fitting complete models")
  Y <- unlist(V[y])
  mfull <- parallel::mclapply(input, function(g) {
    sapply(g, function(v) {
      yv <- unlist(V[v])
      w <- yweight[[v]]
      m <- sapply(yv, lassoCV, x = c(X, setdiff(Y, unlist(V[g]))))
      sum(m * w)
    })
  }, mc.cores = cores)

  #-----

  # Selection sequence
  cli::cli_progress_bar("Determining order and blocks", total = N, clear = TRUE)
  i <- 0
  while (length(input) > 0) {

    i <- i + 1
    Y <- unlist(V[unlist(ord)])

    m0 <- parallel::mclapply(input, function(g) {
      sapply(g, function(v) {
        yv <- unlist(V[v])
        w <- yweight[[v]]
        m <- sapply(yv, lassoCV, x = c(X, Y))
        sum(m * w)
      })
    }, mc.cores = cores)

    # How does mfull error compare to m0 error?
    # Variable with high 'rel' indicates current iteration error is similar to the minimum possible (full model)
    stopifnot(all(lengths(m0) == lengths(mfull)))
    rel <- sapply(seq_along(m0), function(i) mean(mfull[[i]] / m0[[i]]))

    # Select the input for highest 'rel'
    b <- which.max(rel)  # best

    check <-

      # Check if variable should be added to block
      if (i > 1) {
        bsize <- length(ord[[i - 1]]) + length(input[[b]])  # Subsequent block size if best input variables are added to current block (cannot exceed 'maxsize')
        if (bsize <= maxsize) {
          if (pmax(0, rel[b] - rel0[b]) < delta) {
            i <- i - 1
            rel[-b] <- rel0[-b]
          }
        }
      }

    ratio[[i]] <- c(ratio[[i]], rel[b])  # Could print or return this result?
    ord[[i]] <- c(ord[[i]], input[[b]])
    input <- input[-b]
    mfull <- mfull[-b]
    rel0 <- rel[-b]

    cli::cli_progress_update()

  }

  cli::cli_progress_done()
  cli::cli_alert_success("Determining order and blocks")

  #-----

  # Return preferred order and blocks
  ord <- if (all(lengths(ord) == 1)) unlist(ord) else purrr::compact(ord)
  return(ord)

}
