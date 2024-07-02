#' Ensure a monotonic relationship between two variables
#'
#' @description
#' \code{makeMonotonic()} returns modified values of input vector \code{y} that are smoothed, monotonic, and consistent across all values of input \code{x}. It was designed to be used post-fusion when one wants to ensure a plausible relationship between consumption (\code{x}) and expenditure (\code{y}), under the assumption that all consumers face an identical, monotonic pricing structure. By default, the mean of the returned (\code{y}) values is forced to equal the original mean value (\code{preserve = TRUE}). The direction of monotonicity (increasing or decreasing) is detected automatically, so use cases are not necessarily limited to consumption and expenditure variables.
#' @param x Numeric.
#' @param y Numeric.
#' @param w Numeric. Optional observation weights.
#' @param N Integer. Size of random sample to use (if necessary) to reduce computation time.
#' @param preserve Logical. Preserve the original mean of the \code{y} values in the returned values?
#' @param plot Logical. Plot the (sampled) data points and derived monotonic relationship?
#' @details The smoothing is accomplished via a \code{\link[scam]{scam}} model with either a monotone increasing or decreasing constraint, depending on the correlation of the input data. An additional step checks the derivative of the smoothed predictions, eliminates any observations with outlier derivative values (i.e. unusually small or large "jumps"), and then fits a monotonic spline to derive the final relationship. If the SCAM model fails to fit, the function falls back to \code{\link[stats]{loess}} with predictions coerced to monotonic. If LOESS fails, the function falls back to \code{\link[stats]{lm}} with simple linear predictions. If \code{y = 0} when \code{x = 0} (as typical for consumption-expenditure variables), then that outcome is enforced in the result.
#' @return A numeric vector of modified \code{y} values. Optionally, a plot showing the returned monotonic relationship.
#' @examples
#' y <- makeMonotonic(x = recs$propane_btu, y = recs$propane_expend, plot = TRUE)
#' mean(recs$propane_expend)
#' mean(y)
#' @export

#---------

# TEST
# library(dplyr)
# library(data.table())
#
# d <- fusionModel::read_fsd("~/Downloads/RECS_2020_2019_fused_UP.fsd")
# acs <- fst::read_fst("~/Documents/Projects/fusionData/survey-processed/ACS/2019/ACS_2019_H_processed.fst")
#
# test <- d %>%
#   filter(M == 1) %>%
#   cbind(acs) %>%
#   filter(state == 25) %>%
#   as.data.table()
#
# test[, dollarlp_z :=  makeMonotonic(x = btulp, y = dollarlp, w = weight), by = .(state, puma10)]

#---------

makeMonotonic <- function(x,
                          y,
                          w = NULL,
                          N = 5000,
                          preserve = TRUE,
                          plot = FALSE) {

  stopifnot({
    length(x) == length(y)
    is.numeric(x) & !anyNA(x)
    is.numeric(y) & !anyNA(y)
    is.null(w) | length(w) == length(x)
    N >= 100
    is.logical(preserve)
    is.logical(plot)
  })

  if (is.null(w)) w <- rep.int(1L, length(x))
  x0 <- x; w0 <- w
  xu <- sort(unique(x))
  #yu <- sort(unique(y))
  yrng <- range(y)
  ymean <- weighted.mean(y, w)

  # If zeros in 'x' (almost) always produce zeros in 'y', restrict to non-zero observations in 'x'
  if (any(x == 0) & sum(y[x == 0] == 0) / sum(x == 0) > 0.999) {
    force.zero <- TRUE
    i <- x != 0
    x <- x[i]
    y <- y[i]
    w <- w[i]
  } else {
    force.zero <- FALSE
  }

  # Sample the data for fitting SCAM model
  n <- length(x)
  if (n > N) {
    i <- sample.int(n = n, size = N)
    x <- x[i]
    y <- y[i]
    w <- w[i]
  }

  # Wrapper around tryCatch() that traps warnings as well as errors
  try2 <- function(...) tryCatch(..., error = function(e) e, warning = function(w) w)

  # Is y decreasing or increasing with x?
  inc <- suppressWarnings(cor(x, y) > 0)

  # Fit monotonic SCAM model
  m <- try2(scam::scam(y ~ s(x, bs = ifelse(inc, "mpi", "mpd")), data = data.frame(x, y), weights = w))
  # If the SCAM model is successful, use the model to predict 'y' for unique 'x' values
  if (inherits(m, "scam")) {
    p <- as.vector(predict(m, newdata = data.frame(x = xu), type = "response"))
  } else {
    # If SCAM model fails, fall back to LOESS model and coerce the predictions to monotonic
    m <- try2(stats::loess(y ~ x, weights = w, control = loess.control(surface = "direct")))
    if (inherits(m, "loess")) {
      p <- predict(m, newdata = data.frame(x = xu))
      p <- sort(p) # Force monotonic predictions
      if (!inc) p <- rev(p)
    } else {
      # If LOESS model fails, fall back to OLS model and make simple linear predictions
      m <- stats::lm(y ~ x, weights = w)
      p <- as.vector(suppressWarnings(predict(m, newdata = data.frame(x = xu))))
    }
  }

  if (length(p) > 1) {

    # Fit monotonic spline to the predictions
    spf <- splinefun(xu, p, method = "monoH.FC")

    # Set observations with zero-value or outlier derivatives to NA
    spd <- spf(xu, deriv = 1)
    spd[spd == 0] <- NA
    z <- (spd - median(spd, na.rm = TRUE)) / mad(spd, na.rm = TRUE)
    spd[abs(z) > 3.5] <- NA

    # Re-fit monotonic spline, interpolating over any initially "flat" (zero derivative) or unusually steep locales
    xu2 <- xu[!is.na(spd)]
    p2 <- p[!is.na(spd)]
    spf <- splinefun(xu2, p2, method = "monoH.FC")

    # Make 'y' predictions for all 'x'
    yout <- if (length(xu) < 1000) {
      spf(x0)  # Exact but comparatively slow
    } else {
      # Linear approximation is faster for large 'xu'
      approx(xu, spf(xu), xout = x0)$y
    }

  } else {
    yout <- rep(p, length(x0))
  }

  # If necessary, set values to zero when 'x' is zero
  if (force.zero) yout[x0 == 0] = 0

  # Adjustment factor to ensure mean of transformed 'y' matches original mean value
  if (preserve) {
    yadj <- ymean / weighted.mean(yout, w0)
    yout <- yout * yadj
  }

  # NOT IMPLEMENTED
  # Optional step: Force the predicted values to nearest observed/original y-value
  # ind <- findInterval(x = yout, vec = c(-Inf, yu[-1] - diff(yu) / 2, Inf), rightmost.closed = TRUE)
  # yout <- yout[ind]
  # mean(yout)  # !!! Can change the mean considerably

  # Optional plot of transformation
  if (plot) {
    plot(x, y, col = "#00000033", xlim = range(xu), ylim = yrng, xlab = "x", ylab = "y")
    if (length(x0) > 100) {
      yp <- spf(xu)
      if (force.zero) yp[xu == 0] <- 0
      lines(xu, yp * yadj, col = "red")
    } else {
      if (length(x0) == 1) {
        points(x0, yout, col = "red")
      } else {
        lines(x0, yout, col = "red")
      }
    }
  }

  # Return transformed y-values
  return(yout)

}
