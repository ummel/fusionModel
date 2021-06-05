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

induceCor <- function(x, rho, y = NULL, threshold = 1e-12) {

  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.vector(x)) x <- matrix(x, ncol = 1)

  stopifnot(exprs = {
    class(rho) == "numeric"
    length(rho) == ncol(x)
    all(rho <= 1)
    all(rho >= -1)
  })

  d <- ncol(x)
  n <- nrow(x)

  # Makes computations simpler
  x <- scale(x)

  # Random initial Gaussian, if none provided
  if (is.null(y)) y <- rnorm(n)
  y.mu <- mean(y)
  y.sd <- sd(y)
  y <- scale(y)

  # Remove the effects of `y` on `x`
  # On speed of various lm() implementations: https://stackoverflow.com/questions/25416413/is-there-a-faster-lm-function
  #e <- residuals(lm(y ~ x))
  #e <- lm.fit(x = x, y = y)$residuals  # Considerably faster than lm()
  e <- .lm.fit(x = x, y = y)$residuals  # Even faster than lm.fit

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
  result <- as.numeric(z * y.sd + y.mu)

  return(result)

}
