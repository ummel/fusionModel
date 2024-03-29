% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepXY.R
\name{prepXY}
\alias{prepXY}
\title{Prepare the 'x' and 'y' inputs}
\usage{
prepXY(
  data,
  y,
  x,
  weight = NULL,
  cor_thresh = 0.05,
  lasso_thresh = 0.95,
  xmax = 100,
  xforce = NULL,
  fraction = 1,
  cores = 1
)
}
\arguments{
\item{data}{Data frame. Training dataset. All categorical variables should be factors and ordered whenever possible.}

\item{y}{Character or list. Variables in \code{data} to eventually fuse to a recipient dataset. If \code{y} is a list, each entry is a character vector possibly indicating multiple variables to fuse as a block.}

\item{x}{Character. Predictor variables in \code{data} common to donor and eventual recipient.}

\item{weight}{Character. Name of the observation weights column in \code{data}. If NULL (default), uniform weights are assumed.}

\item{cor_thresh}{Numeric. Predictors that exhibit less than \code{cor_thresh} absolute Spearman (rank) correlation with a \code{y} variable are screened out prior to the LASSO step. Fast exclusion of predictors that the LASSO step probably doesn't need to consider.}

\item{lasso_thresh}{Numeric. Controls how aggressively the LASSO step screens out predictors. Lower value is more aggressive. \code{lasso_thresh = 0.95}, for example, retains predictors that collectively explain at least 95\% of the deviance explained by a "full" model.}

\item{xmax}{Integer. Maximum number of predictors returned by LASSO step. Does not strictly control the number of final predictors returned (especially for categorical \code{y} variables), but useful for setting a (very) soft upper bound. Lower \code{xmax} can help control computation time if a large number of \code{x} pass the correlation screen. \code{xmax = Inf} imposes no restriction.}

\item{xforce}{Character. Subset of \code{x} variables to "force" as included predictors in the results.}

\item{fraction}{Numeric. Fraction of observations in \code{data} to randomly sample. For larger datasets, sampling often has minimal effect on results but speeds up computation.}

\item{cores}{Integer. Number of cores used. Only applicable on Unix systems.}
}
\value{
List with named slots "y" and "x". Each is a list of the same length. Former gives the preferred fusion order. Latter gives the preferred sets of predictor variables.
}
\description{
Optional-but-useful function to: 1) provide a plausible ordering of the 'y' (fusion) variables and 2) identify the subset of 'x' (predictor) variables likely to be consequential during subsequent model training. Output can be passed directly to \code{\link{train}}. Most useful for large datasets with many and/or highly-correlated predictors. Employs an absolute Spearman rank correlation screen and then LASSO models (via \code{\link[glmnet]{glmnet}}) to return a plausible ordering of 'y' and the preferred subset of 'x' variables associated with each.
}
\examples{
y <- names(recs)[c(14:16, 20:22)]
x <- names(recs)[2:13]

# Fusion variable "blocks" are respected by prepXY()
y <- c(list(y[1:2]), y[-c(1:2)])

# Do the prep work...
prep <- prepXY(data = recs, y = y, x = x)

# The result can be passed to train()
train(data = recs, y = prep$y, x = prep$x)

}
