#' Fast initial screening of predictor variables
#'
#' @description
#' Quickly detect which predictor variables in a proposed analysis are likely to be consequential in more detailed models. Useful when there are many predictor variables and not all are likely to be material to subsequent model-fitting. \code{prescreen()} employs an absolute correlation screen and then a LASSO model to identify a plausible predictor subset for each response variable.
#'
#' @param data Data frame. Training dataset. All categorical variables should be factors and ordered whenever possible.
#' @param y Character. Response variable(s).
#' @param x Character. Predictor variables.
#' @param weight Character. Name of the observation weights column in \code{data}. If NULL (default), uniform weights are assumed.
#' @param fraction Numeric. Fraction of observations in \code{data} to randomly sample. For larger datasets, sampling often has minimal effect on results but speeds up computation.
#' @param cor_thresh Numeric. Predictors that exhibit less than \code{cor_thresh} absolute correlation with a \code{y} variable are screened out prior to the LASSO step.
#' @param lasso_thresh Numeric. Controls how aggressively the LASSO step screens out predictors. Lower value is more aggressive. \code{lasso_thresh = 0.95}, for example, retains predictors that collectively explain at least 95% of the deviance explained by a "full" model.
#' @param cores Integer. Number of cores used. Only applicable on Unix systems.
#'
#' @return Character vector giving the \code{x} variables that passed both the correlation and LASSO screens when predicting at least one of the \code{y} variables.
#'
#' @examples
#' ?recs
#' fusion.vars <- names(recs)[13:18]
#' predictor.vars <- names(recs)[2:12]
#' xretain <- prescreen(data = recs, y = fusion.vars, x = predictor.vars)
#' xretain
#'
#' @export

prescreen <- function(data,
                      y,
                      x,
                      weight = NULL,
                      fraction = 1,
                      cor_thresh = 0.02,
                      lasso_thresh = 0.95,
                      cores = 1) {

  stopifnot(exprs = {
    is.data.frame(data)
    all(y %in% names(data))
    all(x %in% names(data))
    is.null(weight) | weight %in% names(data)
    fraction > 0 & fraction <= 1
    cor_thresh > 0 & cor_thresh <= 1
    lasso_thresh > 0 & lasso_thresh <= 1
    cores >= 1 & cores %% 1 == 0
  })

  if (is.data.table(data)) data <- as.data.frame(data)

  # Check for character-type variables; stop with error if any detected
  # Check for no-variance (constant) variables
  # Detect and impute any missing values in 'x' variables
  data <- checkData(data, y, x)
  x <- intersect(x, names(data))

  # Observation weights vector
  W <- if (is.null(weight)) {
    rep(1L, nrow(data))
  } else {
    data[[weight]] / mean(data[[weight]])
  }

  #-----

  # Sample 'data', if requested
  if (fraction < 1) {
    samp <- sample.int(n = nrow(data), size = round(nrow(data) * fraction))
    data <- data[samp, ]
    W <- W[samp]
  }

  # Prepare for one-hot encoding
  data <- data[c(y, x)] %>%
    mutate_if(is.ordered, as.integer) %>%
    mutate_if(is.numeric, data.table::frank) %>%   # Rank numeric variables
    mutate_if(is.logical, as.integer)

  # One-hot encode the response and predictors (separately)
  Y <- one_hot(data[y], sep = "||")
  X <- one_hot(data[x], sep = "||")
  rm(data)

  #-----

  # Loop using each column in 'Y' as a response variable
  # Output are the 'x' variables used by any one of the LASSO models
  retain <- parallel::mclapply(colnames(Y), function(v) {

    # Initial correlation screening, based on absolute correlation value
    p <- abs(cor(Y[, v], as.matrix(X)))
    vx <- which(p > cor_thresh)  # Arbitrary correlation threshold
    if (length(vx) < 20) vx <- order(p, decreasing = TRUE)[1:min(20, length(p))]  # Ensure minimum number of predictors are passed to glmnet()

    m <- glmnet::glmnet(
      x = X[, vx],
      y = Y[, v],
      weights = W,
      family = "gaussian",
      alpha = 1)

    i <- which(m$dev.ratio / max(m$dev.ratio) >= lasso_thresh)[1]  # Index of preferred lambda, based on arbitrary lasso threshold
    cf <- coef(m, s = m$lambda[i])
    keep <- names(which(cf[-1, 1] != 0))  # Ignore intercept term
    keep <- purrr::map_chr(strsplit(keep, "||", fixed = TRUE), 1)  # Extract the original variable name
    keep <- unique(keep)
    return(keep)

  }, mc.cores = cores) %>%
    unlist() %>%
    unique()

  # Order the retained variables to ordering of original 'x'
  ord <- order(match(retain, x))
  retain <- retain[ord]

  cat("Retained", length(retain), "of", length(x), "predictor variables.")
  return(retain)

}
