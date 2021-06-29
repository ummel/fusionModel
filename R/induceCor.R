# Returns a vector with specified linear correlations with matrix of existing variables 'x'
# Note: Target correlation may be impossible, in which case an error is returned

# Source:
# https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
# https://stats.stackexchange.com/questions/444039/whuber-s-generation-of-a-random-variable-with-fixed-covariance-structure

# # Example using discrete x's
# x <- recs[c("income")]
# y <- y0 <- rank(recs$electricity)
# threshold = 1e-12
#
# # Convert 'x' to dummies...
# x <- model.matrix(~ 0 + ., data = x)[, -1]
#
# # Target correlations
# rho <- cor(x, y)[, 1] + runif(ncol(x), -0.2, 0.2)

# TEST CASE
# x <- matrix(c(rep(1, 1000), rep(0, 2000)), ncol = 1)
# x <- cbind(x, matrix(c(rep(1, 2000), rep(0, 1000)), ncol = 1))
# y <- y0 <- c(runif(1000, 0.5, 1), runif(2000, 0.2, 0.6))
# rho <- c(0.7, 0.5)
# threshold <- 1e-12

#---------------------------

# data: a data.table
# rho: named vector of target correlations; all names must be in 'data'

induceCor <- function(data, rho, y = NULL, scale.data = FALSE, use.biglm = FALSE, threshold = 1e-12) {

  # if (is.data.frame(x)) x <- as.matrix(x)
  # if (is.vector(x)) x <- matrix(x, ncol = 1)

  stopifnot(exprs = {
    is.data.frame(data)
    all(names(rho) %in% names(data))
    is.null(y) | y %in% names(data)
    class(rho) == "numeric"
    all(rho <= 1)
    all(rho >= -1)
  })

  d <- length(rho)
  n <- nrow(data)

  # Random initial Gaussian for 'y', if none provided
  if (is.null(y)) data.table::set(data, j = y, value = rnorm(n))

  # Scale 'y'; retain mean and standard deviation to undo scaling at very end
  y.mu <- mean(data[[y]])
  y.sd <- sd(data[[y]])

  # Reorder 'data' columns so the first columns match names in 'rho'
  # This moves any unspecified columns (e.g. 'y') to after the 'rho'-name columns
  cols <- c(y, names(rho))
  suppressWarnings(data.table::set(data, j = setdiff(names(data), cols), value = NULL))  # Remove unnecessary columns
  data.table::setcolorder(data, cols)  # Order columns
  data.table::setnames(data, new = make.names(names(data)))  # Ensure syntactically valid names

  # Scale all variables, if requested
  # Makes computations simpler according to original code
  # NOTE: When called within fuse(), the 'ranks' input data are already scaled
  if (scale.data) for (v in names(data)) data.table::set(data, j = v, value = as.vector(scale(data[[v]])))

  #-----

  # Remove the effects of `y` on `x`

  if (use.biglm) {

    # Use biglm() to compute residuals
    # Is faster and more memory efficient for especially large 'data'
    mod <- biglm::biglm(formula = as.formula(paste(names(data)[1], "~", paste(names(data)[-1], collapse = " + "))), data = data)
    x <- as.matrix(data[, -1L])
    y <- data[[y]]
    rm(data)
    gc()

    beta <- coef(mod)
    beta[is.na(beta)] <- 0
    e <- y - (beta[1L] + drop(x %*% beta[-1L]))  # Model residuals; drop intercept

  } else {

    # Use .lm.fit() to compute residuals
    # On speed of various conventional lm() implementations: https://stackoverflow.com/questions/25416413/is-there-a-faster-lm-function
    # NOTE: .lm.fit() is faster than either lm.fit() or lm()
    x <- as.matrix(data[, -1L])
    e <- .lm.fit(x = x, y = data[[1L]])$residuals
    y <- data[[y]]
    rm(data)
    gc()

  }

  gc()

  #-----

  # Calculate the coefficient `sigma` of `e` so that the correlation of
  # `x` with the linear combination x.dual %*% rho + sigma*e is the desired vector.
  if (d == 1) {
    x.dual <- with(svd(x), (n - 1) * u %*% matrix(ifelse(d > threshold, 1 / d, 0)) %*% t(v))  # This prevents non-conformable array error when d = 1
  } else {
    x.dual <- with(svd(x), (n - 1) * u %*% diag(ifelse(d > threshold, 1 / d, 0)) %*% t(v))
  }
  sigma2 <- c((1 - rho %*% cov(x.dual) %*% rho) / var(e))

  # Return the linear combination
  if (sigma2 >= 0) {
    sigma <- sqrt(sigma2)
    z <- x.dual %*% rho + sigma * e
  } else {
    z <- y  # Return original values
    #stop("Correlations are impossible.")
  }

  # Transform 'z' to original units of 'y'
  result <- list(Y = as.numeric(z * y.sd + y.mu), sigma2 = sigma2)

  return(result)

}
