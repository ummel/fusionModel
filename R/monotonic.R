#' Ensure a monotonic relationship between two variables
#'
#' @description
#' \code{monotonic()} returns modified values of input vector \code{y} that are smoothed, monotonic, and consistent across all values of input \code{x}. It was designed to be used post-fusion when one wants to ensure a plausible relationship between consumption (\code{x}) and expenditure (\code{y}), under the assumption that all consumers face an identical, monotonic pricing structure. By default, the mean of the returned values is forced to equal the original mean of \code{y} (\code{preserve_mean = TRUE}). The direction of monotonicity (increasing or decreasing) is detected automatically, so use cases are not limited to consumption and expenditure variables.
#' @param x Numeric.
#' @param y Numeric.
#' @param w Numeric. Optional observation weights.
#' @param preserve_mean Logical. Preserve the original mean of the \code{y} values in the returned values?
#' @param preserve_type Logical. Preserve the original data type of the \code{y} values in the returned values?
#' @param plot Logical. Plot the (sampled) data points and derived monotonic relationship?
#' @details The initial smoothing is accomplished via \code{\link[scam]{supsmu}} with the result coerced to monotone. If the coercion step modifies the values too much, a second smooth is attempted via a \code{\link[scam]{scam}} model with either a monotone increasing or decreasing constraint. If the SCAM fails to fit, the function falls back to \code{\link[stats]{lm}} with simple linear predictions. If \code{y = 0} when \code{x = 0} (as typical for consumption-expenditure variables), then that outcome is enforced in the result. The input data are randomly sampled to no more than 10,000 observations, if necessary, for speed.
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
                      preserve_mean = TRUE,
                      preserve_type = TRUE,
                      plot = FALSE) {

  stopifnot(exprs = {
    length(x) == length(y)
    is.numeric(x) & !anyNA(x)
    is.numeric(y) & !anyNA(y)
    is.null(w) | length(w) == length(x)
    is.logical(preserve_mean)
    is.logical(preserve_type)
    is.logical(plot)
  })

  if (is.null(w)) w <- rep.int(1L, length(x))
  ymean <- weighted.mean(y, w)
  ytype <- storage.mode(y)
  x0 <- x
  w0 <- w

  # If zeros in 'x' (almost) always produce zeros in 'y', restrict to non-zero observations in 'x'
  force.zero <- FALSE
  if (any(x == 0) & sum(y[x == 0] == 0) / sum(x == 0) > 0.999) {
    force.zero <- TRUE
    i <- c(which(x != 0), match(0, x))  # Retains first instance of zero in 'x'
    x <- x[i]
    y <- y[i]
    w <- w[i]
  }

  # Sample the data for speed, if necessary
  n <- length(x)
  if (n > 10e3) {
    i <- match(range(x), x)  # Retains first instance of min and max 'x'
    i <- c(i, sample.int(n = n, size = n - 2))
    x <- x[i]
    y <- y[i]
    w <- w[i]
  }

  m <- stats::supsmu(x, y, wt = w)
  xu <- m$x
  inc <- suppressWarnings(cor(m$x, m$y) >= 0)
  if (is.na(inc)) inc <- TRUE
  p <- sort(m$y, decreasing = !inc) # Force monotonic predictions
  delta <- mean(abs((p - m$y) / m$y))
  if (is.na(delta)) delta <- 0
  if (delta > 0.005) {
    m <- try(scam::scam(y ~ s(x, bs = ifelse(inc, "mpi", "mpd")), data = data.frame(x, y), weights = w), silent = TRUE)
    if (inherits(m, "scam")) {
      p <- as.vector(predict(m, newdata = data.frame(x = xu), type = "response", newdata.guaranteed = TRUE))
    } else {
      # If SCAM model fails, fall back to OLS model and make simple linear predictions
      m <- stats::lm(y ~ x, weights = w)
      p <- as.vector(suppressWarnings(predict(m, newdata = data.frame(x = xu))))
    }
  }

  # Make 'y' predictions for all original 'x'
  yout <- if (length(xu) == 1) rep(p, length(x0)) else approx(xu, p, xout = x0, rule = 2)$y

  # If necessary, set values to zero when 'x' is zero
  if (force.zero) yout[x0 == 0] = 0

  # Safety check
  if (anyNA(yout)) stop("NA values in result vector")

  # If requested, adjustment factor to ensure mean of transformed 'y' matches original mean value
  yadj <- 1  # Defined for use in plotting code, below, if 'preserve_mean = FALSE'
  if (preserve_mean) {
    yadj <- ymean / weighted.mean(yout, w0)
    if (is.na(yadj)) yadj <- 1  # Catch divide by zero case
    yout <- yout * yadj
  }

  # If requested, force output type to match input type
  if (preserve_type) storage.mode(yout) <- ytype

  # Optional plot of transformation
  if (plot) {
    plot(x, y, col = "#00000033", ylim = range(c(y, p)), xlab = "x", ylab = "y")
    lines(xu, p * yadj, col = "red")
  }

  return(yout)

}
