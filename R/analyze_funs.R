weighted_mean <- function(x, w1, w2) {

  if (!(is.numeric(x) | is.logical(x)) | length(x) != length(w1) | length(x) != length(w2)) stop("Input vectors are invalid")

  out <- if (length(x) == 0) {
    rep(NA, 3)
  } else {

    # If we need mean and variance for BOTH sets of weights:

    # Combine the weights
    # w <- cbind(w1 / sum(w1), w2 / sum(w2))
    #
    # # Continuous case (ONLY)
    # # Code: https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
    # # Computes the variance of a weighted mean following Cochran 1977: https://www.sciencedirect.com/science/article/abs/pii/135223109400210C
    # n <- nrow(w)
    # xWbar <- colSums(x * w)  # Weighted mean, since the weights are already scaled
    # wbar <- colMeans(w)
    #
    # # Single weights case - confirms vectorized version
    # # W <- w[, 1]
    # # n/((n-1)*sum(W)^2)*(sum((W*x-wbar[1]*xWbar[1])^2)-2*xWbar[1]*sum((W-wbar[1])*(W*x-wbar[1]*xWbar[1]))+xWbar[1]^2*sum((W-wbar[1])^2))
    # # W <- w[, 2]
    # # n/((n-1)*sum(W)^2)*(sum((W*x-wbar[2]*xWbar[2])^2)-2*xWbar[2]*sum((W-wbar[2])*(W*x-wbar[2]*xWbar[2]))+xWbar[2]^2*sum((W-wbar[2])^2))
    #
    # a <- sweep(x * w, 2, wbar * xWbar)  # Equiv. to: w*x-wbar*xWbar
    # b <- sweep(w, 2, wbar) # Equiv. to: w-wbar
    #
    # wvar <- n/((n-1)*colSums(w)^2) * (colSums(a ^ 2) - 2*xWbar*colSums(b*a) + xWbar^2*colSums(b^2))

    #---

    W <- w1

    if (any(!x %in% 0:1)) {

      # Continuous case
      # Code: https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
      # Computes the variance of a weighted mean following Cochran 1977: https://www.sciencedirect.com/science/article/abs/pii/135223109400210C
      n <- length(x)
      xWbar <- wmean(x, W)
      wbar <- mean(W)
      wvar <- n/((n-1)*sum(W)^2)*(sum((W*x-wbar*xWbar)^2)-2*xWbar*sum((W-wbar)*(W*x-wbar*xWbar))+xWbar^2*sum((W-wbar)^2))

    } else {

      # Binary case
      # https://stats.stackexchange.com/questions/159204/how-to-calculate-the-standard-error-of-a-proportion-using-weighted-data
      w <- W / sum(W)
      xWbar <- sum(w * x)
      sw2 <- sum(w ^ 2)
      wvar <- xWbar * (1 - xWbar) * sw2

    }

    # Secondary (replicate) weighted mean
    xWbar2 <- wmean(x, w2)

    # Return the primary weighted mean, secondary (replicate) weighted mean, and primary weights within-sample variance
    c(xWbar, xWbar2, wvar)

  }

  return(out)

}

#------------------------

weighted_sum <- function(x, w1, w2) {
  if (!(is.numeric(x) | is.logical(x)) | length(x) != length(w1) | length(x) != length(w2)) stop("Input vectors are invalid")
  out <- weighted_mean(x, w1, w2)
  s <- sum(w1)
  out[1] <- out[1] * s
  out[2] <- out[2] * sum(w2)
  out[3] <- out[3] * s ^ 2
  return(out)
}

#------------------------

# nboot: number of bootstrap samples for small N variance estimation
weighted_median <- function(x, w1, w2, nboot = 200) {

  n <- length(x)

  out <- if (var(x) > 0) {

    # Median estimate using 'w1'
    med1 <- matrixStats::weightedMedian(x, w1)

    # Median estimate using 'w2'
    med2 <- matrixStats::weightedMedian(x, w2)

    # If 'x' is large, use approximation; otherwise, do bootstrap estimate of the sampling variance
    wvar <- if (n > 1000 & med1 != 0) {

      # Approximate variance of sample median for large N
      # https://stats.stackexchange.com/questions/45124/central-limit-theorem-for-sample-medians
      pdf <- density(x, weights = w1 / sum(w1))
      fu <- approx(pdf$x, pdf$y, xout = med1)$y
      1 / (4 * n * fu ^ 2)

    } else {

      # Draw bootstrap samples
      bdraw <- replicate(n = nboot, sample(x, size = n, replace = TRUE))

      # Bootstrap weighted median estimates using 'w1'
      bq <- matrixStats::colWeightedMedians(bdraw, w = w1)

      # Bootstrapped variance estimate
      var(bq)

    }

    c(med1, med2, wvar)

  } else {
    c(x[1], x[1], 0)
  }

  return(out)

}

# REGRESSION - TO DO
# f_glm <- function(data, formula, w1, w2) {
#
#   fit1 <- stats::glm(formula = formula, data = data, weights = w1)
#   cf1 <- coef(summary(fit1))
#
#   fit2 <- stats::glm(formula = formula, data = data, weights = w2)
#   cf2 <- coef(summary(fit2))
#
# }

# Other cases to do....
# regression (y ~ x1 + x2) -- inputs a formula
# bivariate smooth (y ~ x) -- inputs a formula
# density -- inputs a length 1 character
# Custom function -- inputs a character of any length or NULL (?)
