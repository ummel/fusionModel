#' Enforce Monotonic Relationships Between Numerical Variables
#'
#' @description
#' Smoothes and transforms a numeric target vector (\code{y}) relative to an ordering
#' vector (\code{x}) to enforce a strict, monotonic relationship (either consistently
#' increasing or decreasing) across all observed values of \code{x}. Designed primarily
#' for post-fusion adjustment of energy consumption (\code{x}) and expenditure (\code{y})
#' variables to guarantee a plausible pricing structure (e.g., higher fuel volume never
#' yields lower total cost). By default, the weighted mean of \code{y} is preserved in
#' the returned vector (\code{preserve = TRUE}).
#'
#' @param x Numeric vector. The independent ordering variable (e.g., energy consumption in BTU).
#'   Must not contain missing values (\code{NA}).
#' @param y Numeric vector of the same length as \code{x}. The dependent response
#'   variable to be transformed (e.g., fuel expenditure in dollars). Must not contain
#'   missing values (\code{NA}).
#' @param w Numeric vector of the same length as \code{x}, optional. Observation sampling
#'   weights. If \code{NULL} (default), uniform weights equal to 1 are assumed.
#' @param preserve Logical. Should the original weighted mean of \code{y} be preserved in the
#'   returned numeric vector? Defaults to \code{TRUE}.
#' @param expend Logical. Treat \code{y} as an expenditure variable tied to consumption \code{x}?
#'   If \code{TRUE} (default), safety checks enforce physical consistency: negative values in
#'   \code{x} or \code{y} throw an error, non-zero expenditures when \code{x == 0} are zeroed,
#'   and zero expenditures when \code{x > 0} are raised to the minimum observed positive expenditure.
#' @param fast Logical. If \code{TRUE} (default), rapid smoothing via Friedman's super-smoother
#'   (\code{\link[stats]{supsmu}}) is performed and directly coerced to monotonicity via
#'   sorting. If \code{FALSE}, a Shape Constrained Additive Model (\code{\link[scam]{scam}}) is
#'   attempted if the initial super-smoother coercion results in excessive error.
#' @param nmax Integer. Maximum number of observations sampled for model fitting to optimize
#'   computational speed. Defaults to \code{5000}. Set to \code{Inf} to disable sampling.
#' @param plot Logical. Should a diagnostic scatterplot of the sampled input points and the
#'   fitted monotonic relationship (in red) be rendered to the active graphics device? Defaults to \code{FALSE}.
#'
#' @details
#' \code{monotonic()} provides non-parametric post-processing to rectify logical inconsistencies
#' in microdata fusion output (such as negative marginal prices or non-monotonic tariff curves).
#'
#' \strong{Algorithmic Workflow:}
#' \describe{
#'   \item{1. Physical Sanity Checks}{If \code{expend = TRUE}, boundary conditions are enforced
#'     to ensure zero consumption yields zero expenditure and positive consumption yields positive
#'     expenditure.}
#'   \item{2. High-Speed Subsampling}{If \code{length(x) > nmax}, extreme boundaries (minimum and
#'     maximum) are preserved while intermediate points are randomly down-sampled to \code{nmax}
#'     for efficient smoothing.}
#'   \item{3. Super-Smoother Fit & Direction Detection}{An initial non-parametric smooth is fit
#'     via \code{\link[stats]{supsmu}}. The overall correlation between \code{x} and predicted
#'     \code{y} determines whether the relationship should be monotonic increasing or decreasing.}
#'   \item{4. Monotonic Coercion & SCAM Fallback}{Predicted values are sorted to enforce monotonicity.
#'     If \code{fast = FALSE} and sorting introduces substantial error (>5% relative divergence on
#'     over 5% of points), a constrained spline is fit via \code{\link[scam]{scam}}. If SCAM fitting
#'     fails to converge, linear regression (\code{\link[stats]{lm}}) serves as the ultimate fallback.}
#'   \item{5. Interpolation & Mean Preservation}{Fitted values are mapped back to the complete
#'     original \code{x} vector using linear interpolation (\code{\link[stats]{approx}}). If
#'     \code{preserve = TRUE}, the resulting vector is scaled so its weighted mean matches
#'     \code{weighted.mean(y, w)}.}
#' }
#'
#' @return A numeric vector of modified, monotonic \code{y} values matching the length and order
#'   of input \code{x}. If input \code{y} was integer-typed, the returned vector is rounded to integer.
#'
#' @examples
#' \dontrun{
#' library(fusionModel)
#' data(recs)
#'
#' # Enforce a monotonic pricing curve between propane consumption and expenditure
#' adjusted_expend <- monotonic(
#'   x = recs$propane_btu,
#'   y = recs$propane_expend,
#'   w = recs$weight,
#'   plot = TRUE
#' )
#'
#' # Compare original and adjusted weighted means
#' weighted.mean(recs$propane_expend, recs$weight)
#' weighted.mean(adjusted_expend, recs$weight)
#' }
#'
#' @export

#---------

# TEST
# library(tidyverse)
# library(data.table)
#
# d <- fusionModel::read_fsd("~/Downloads/RECS_2020_2019_fused_UP.fsd")
# acs <- fst::read_fst("~/Documents/Projects/fusionData/survey-processed/ACS/2019/ACS_2019_H_processed.fst", columns = c('weight', 'state', 'puma10'))
# d <- cbind(d, acs)
# system.time(
#    d[, `:=`(dollarel_z =  monotonic(x = btuel, y = dollarel, w = weight),
#             dollarng_z =  monotonic(x = btung, y = dollarng, w = weight),
#             dollarlp_z =  monotonic(x = btulp, y = dollarlp, w = weight),
#             dollarfo_z =  monotonic(x = btufo, y = dollarfo, w = weight)),
#      by = .(state, puma10)]
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
  ymean <- weighted.mean(y, w, na.rm = TRUE)
  yint <- is.integer(y)
  ymin <- if (all(y == 0)) 0 else min(y[y != 0])
  x0 <- x
  w0 <- w

  # If 'expend = TRUE', check for violations
  # If any issues are detected, a helpful warning is issued
  if (expend) {
    if (any(x < 0)) stop("'expend = TRUE' but detected negative values in 'x'")
    if (any(y < 0)) stop("'expend = TRUE' but detected negative values in 'y'")
    i <- x == 0 & y != 0
    if (any(i)) {
      y[i] <- 0L
      cli::cli_warn("Set {sum(i)} non-zero y-value(s) ({round(100 * sum(i) / length(y), 2)}%) to zero where x == 0 because 'expend = TRUE'")
    }
    i <- x > 0 & y == 0
    if (any(i)) {
      y[i] <- ymin
      cli::cli_warn("Set {sum(i)} zero y-value(s) ({round(100 * sum(i) / length(y), 2)}%) to observed non-zero minimum where x > 0 because 'expend = TRUE'")
    }
  }

  # If 'expend = TRUE' OR zeros in 'x' (almost) always produce zeros in 'y', restrict to non-zero observations in 'x'
  force.zero <- FALSE
  if (expend | (any(x == 0) & sum(y[x == 0] == 0) / sum(x == 0) > 0.995)) {
    force.zero <- TRUE
    i <- c(match(0, x), which(x != 0))  # Retains first instance of zero in 'x'
    x <- x[i]
    y <- y[i]
    w <- w[i]
  }

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

  # First time: If expend = TRUE, ensure that the output values meet some minimum positive value when x > 0
  if (expend) p[xu > 0] <- pmax(p[xu > 0], ymin)

  # If necessary, set values to zero when 'x' is zero
  if (force.zero) p[xu == 0] = 0

  # Make 'y' predictions for all original 'x'
  yout <- if (length(xu) == 1) rep(p, length(x0)) else approx(xu, p, xout = x0, rule = 2)$y

  # If requested, adjustment factor to ensure mean of transformed 'y' matches original mean value
  yadj <- 1  # Defined for use in plotting code, below, if 'preserve = FALSE'
  if (preserve) {
    yadj <- ymean / weighted.mean(yout, w0)
    if (!is.finite(yadj)) yadj <- 1  # Catch divide by zero case (mean not strictly preserved in this case)
    yout <- yout * yadj
  }

  # If 'y' input is integer, force output to integer
  # May cause input mean of 'y' to not be strictly preserved in output
  if (yint) yout <- as.integer(round(yout))

  # May cause input mean of 'y' to NOT be strictly preserved in output
  if (expend) yout[x0 > 0] <- pmax(yout[x0 > 0], ymin)

  # Optional plot of transformation
  if (plot) {
    plot(x, y, col = "#00000033", ylim = range(c(y, p)), xlab = "x", ylab = "y")
    lines(xu, p * yadj, col = "red")
  }

  # Safety check
  if (anyNA(yout)) stop("NA values in result vector")

  return(yout)

}
