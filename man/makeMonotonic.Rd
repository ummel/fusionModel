% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeMonotonic.R
\name{makeMonotonic}
\alias{makeMonotonic}
\title{Ensure a monotonic relationship between two variables}
\usage{
makeMonotonic(x, y, w = NULL, N = 5000, preserve = TRUE, plot = FALSE)
}
\arguments{
\item{x}{Numeric.}

\item{y}{Numeric.}

\item{w}{Numeric. Optional observation weights.}

\item{N}{Integer. Size of random sample to use (if necessary) to reduce computation time.}

\item{preserve}{Logical. Preserve the original mean of the \code{y} values in the returned values?}

\item{plot}{Logical. Plot the original data points and returned monotonic relationship?}
}
\value{
A numeric vector of modified \code{y} values. Optionally, a plot showing the returned monotonic relationship.
}
\description{
\code{makeMonotonic()} returns modified values of input vector \code{y} that are smoothed, monotonic, and consistent across all values of input \code{x}. It was designed to be used post-fusion when one wants to ensure a plausible relationship between consumption (\code{x}) and expenditure (\code{y}), under the assumption that all consumers face an identical, monotonic pricing structure. By default, the mean of the returned (\code{y}) values is forced to equal the original mean value (\code{preserve = TRUE}). The direction of monotonicity (increasing or decreasing) is detected automatically, so use cases are not necessarily limited to consumption and expenditure variables.
}
\details{
The smoothing is accomplished via a \code{\link[scam]{scam}} model with either a monotone increasing or decreasing constraint, depending on the correlation of the input data. An additional step checks the derivative of the smoothed predictions, eliminates any observations with outlier derivative values (i.e. unusually small or large "jumps"), and then fits a monotonic spline to derive the final relationship. If \code{y = 0} when \code{x = 0} (as typical for consumption-expenditure variables), then that outcome is simply forced in the result.
}
\examples{
y <- makeMonotonic(x = recs$propane_btu, y = recs$propane_expend, plot = TRUE)
mean(recs$propane_expend)
mean(y)
}
