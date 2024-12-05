#' Create a monotonic relationship between two variables
#'
#' @description
#' \code{monotonic()} returns modified values of input vector \code{y} that are smoothed, monotonic, and consistent across all values of input \code{x}. It was designed to be used post-fusion when one wants to ensure a plausible relationship between consumption (\code{x}) and expenditure (\code{y}), under the assumption that all consumers face an identical, monotonic pricing structure. By default, the mean of the returned values is forced to equal the original mean of \code{y} (\code{preserve = TRUE}). The direction of monotonicity (increasing or decreasing) is detected automatically, so use cases are not limited to consumption and expenditure variables.
#' @param x Numeric.
#' @param y Numeric.
#' @param w Numeric. Optional observation weights.
#' @param preserve Logical. Preserve the original mean of the \code{y} values in the returned values?
#' @param expend Logical. Assume \code{y} is an expenditure variable? If \code{TRUE}, a safety check is implemented to ensure \code{y > 0} when \code{x > 0}.
#' @param fast Logical. If \code{TRUE}, only \code{\link[scam]{supsmu}} is used with coercion of result to monotone.
#' @param nmax Integer. Maximum number of observations to use for smoothing. Set lower for faster computation. \code{nmax = Inf} eliminates sampling.
#' @param plot Logical. Plot the (sampled) data points and derived monotonic relationship?
#' @details The initial smoothing is accomplished via \code{\link[scam]{supsmu}} with the result coerced to monotone. If \code{fast = FALSE} and the coercion step modifies the values too much, a second smooth is attempted via a \code{\link[scam]{scam}} model with either a monotone increasing or decreasing constraint. If the SCAM fails to fit, the function falls back to \code{\link[stats]{lm}} with simple linear predictions. If \code{y = 0} when \code{x = 0} (as typical for consumption-expenditure variables), then that outcome is enforced in the result. The input data are randomly sampled to no more than \code{nmax} observations, if necessary, for speed.
#' @return A numeric vector of modified \code{y} values. Optionally, a plot showing the returned monotonic relationship.
#' @examples
#' y <- monotonic(x = recs$propane_btu, y = recs$propane_expend, plot = TRUE)
#' mean(recs$propane_expend)
#' mean(y)
#' @export

#---------

# TEST
# library(tidyverse)
# library(data.table())
#
# d <- fusionModel::read_fsd("~/Downloads/RECS_2020_2019_fused_UP.fsd")
# acs <- fst::read_fst("~/Documents/Projects/fusionData/survey-processed/ACS/2019/ACS_2019_H_processed.fst", columns = c('weight', 'state', 'puma10'))
# d <- cbind(d, acs)
# system.time(
#   d[, `:=`(dollarel_z =  monotonic(x = btuel, y = dollarel, w = weight),
#            dollarng_z =  monotonic(x = btung, y = dollarng, w = weight),
#            dollarlp_z =  monotonic(x = btulp, y = dollarlp, w = weight),
#            dollarfo_z =  monotonic(x = btufo, y = dollarfo, w = weight)),
#     by = .(state, puma10)]
# )

#---------

monotonic <- function(x,
                      y,
                      w = NULL,
                      preserve = TRUE,
                      expend = TRUE,
                      fast = TRUE,
                      nmax = 5000,
                      plot = FALSE) {

  stopifnot(exprs = {
    length(x) == length(y)
    is.numeric(x) & !anyNA(x)
    is.numeric(y) & !anyNA(y)
    is.null(w) | length(w) == length(x)
    is.logical(preserve)
    is.logical(expend)
    is.logical(fast)
    nmax > 1
    is.logical(plot)
  })

  if (is.null(w)) w <- rep.int(1L, length(x))
  ymean <- weighted.mean(y, w)
  yint <- is.integer(y)
  ymin <- min(y[y != 0])
  x0 <- x
  w0 <- w
  if (expend & any(x < 0 | y < 0)) warning("'expend = TRUE' but detected negative values in 'x' and/or 'y'")

  # If zeros in 'x' (almost) always produce zeros in 'y', restrict to non-zero observations in 'x'
  force.zero <- FALSE
  if (any(x == 0) & sum(y[x == 0] == 0) / sum(x == 0) > 0.995) {
    force.zero <- TRUE
    i <- c(match(0, x), which(x != 0))  # Retains first instance of zero in 'x'
    x <- x[i]
    y <- y[i]
    w <- w[i]
  }

  # Remove outliers in 'x' and 'y'
  # ok <- abs(x - median(x)) / mad(x) < 3.5 & abs(y - median(y)) / mad(y) < 3.5
  # ok[is.na(ok)] <- TRUE
  # if (!all(ok)) {
  #   i <- match(range(x), x)  # Retains first instance of min and max 'x'
  #   i <- unique(c(i, which(ok)))
  #   x <- x[i]
  #   y <- y[i]
  #   w <- w[i]
  # }

  # If necessary, sample the data for speed
  n <- length(x)
  if (n > nmax) {
    i <- match(range(x), x)  # Retains first instance of min and max 'x'
    i <- c(i, sample.int(n = n, size = nmax - 2))  # Downsample to 'nmax' observations
    x <- x[i]
    y <- y[i]
    w <- w[i]
  }

  # Initial smooth via stats::supsmu()
  # 'span' set based on recommendation in ?supsmu
  m <- stats::supsmu(x, y, wt = w, span = ifelse(length(x) < 40, 0.3, "cv"))
  xu <- m$x
  inc <- suppressWarnings(cor(m$x, m$y) >= 0)
  if (is.na(inc)) inc <- TRUE
  p <- sort(m$y, decreasing = !inc) # Force monotonic predictions

  # Check if supsmu() smoothed and monotonic output is sufficient
  # Ideally, coercion to monotonic via sort() does not cause significant difference between 'p' and m$y
  fail <- sum(abs((p - m$y) / m$y) > 0.05) / length(p)  # Percent of observations with more than 5% absolute error
  if (is.na(fail)) fail <- Inf
  if (!fast & fail > 0.05 & length(p) >= 100) {
    # Attempt to fit SCAM model with monotonic constraint
    m <- try(scam::scam(y ~ s(x, bs = ifelse(inc, "mpi", "mpd")), data = data.frame(x, y), weights = w), silent = TRUE)
    if (inherits(m, "scam")) {
      p <- as.vector(predict(m, newdata = data.frame(x = xu), type = "response", newdata.guaranteed = TRUE))
    } else {
      # If SCAM model fails, fall back to OLS model and make simple linear predictions
      m <- stats::lm(y ~ x, weights = w)
      p <- as.vector(suppressWarnings(predict(m, newdata = data.frame(x = xu))))
    }
  }

  # If 'y' is assumed to be expenditure, ensure that 'p' values meet some minimum positive value
  # This is to prevent possibility of zero values for 'y' where x > 0
  if (expend) p <- pmax(p, ymin)

  # Make 'y' predictions for all original 'x'
  yout <- if (length(xu) == 1) rep(p, length(x0)) else approx(xu, p, xout = x0, rule = 2)$y

  # If necessary, set values to zero when 'x' is zero
  if (force.zero) yout[x0 == 0] = 0

  # Safety check
  if (anyNA(yout)) stop("NA values in result vector")

  # If requested, adjustment factor to ensure mean of transformed 'y' matches original mean value
  yadj <- 1  # Defined for use in plotting code, below, if 'preserve = FALSE'
  if (preserve) {
    yadj <- ymean / weighted.mean(yout, w0)
    if (is.na(yadj)) yadj <- 1  # Catch divide by zero case
    yout <- yout * yadj
  }

  # If 'y' input is integer, force output to integer
  if (yint) yout <- as.integer(round(yout))

  # Optional plot of transformation
  if (plot) {
    plot(x, y, col = "#00000033", ylim = range(c(y, p)), xlab = "x", ylab = "y")
    lines(xu, p * yadj, col = "red")
  }

  return(yout)

}
