#' Prepare the 'x' and 'y' inputs
#'
#' @description
#' Optional-but-useful function to: 1) provide a plausible ordering of the 'y' (fusion) variables and 2) identify the subset of 'x' (predictor) variables likely to be consequential during subsequent model training. Output can be passed directly to \code{\link{train}}. Most useful for large datasets with many and/or highly-correlated predictors. Employs an absolute Spearman rank correlation screen and then LASSO models (via \code{\link[glmnet]{glmnet}}) to return a plausible ordering of 'y' and the preferred subset of 'x' variables associated with each.
#'
#' @param data Data frame. Training dataset. All categorical variables should be factors and ordered whenever possible.
#' @param y Character or list. Variables in \code{data} to eventually fuse to a recipient dataset. If \code{y} is a list, each entry is a character vector possibly indicating multiple variables to fuse as a block.
#' @param x Character. Predictor variables in \code{data} common to donor and eventual recipient.
#' @param weight Character. Name of the observation weights column in \code{data}. If NULL (default), uniform weights are assumed.
#' @param cor_thresh Numeric. Predictors that exhibit less than \code{cor_thresh} absolute Spearman (rank) correlation with a \code{y} variable are screened out prior to the LASSO step. Fast exclusion of predictors that the LASSO step probably doesn't need to consider.
#' @param lasso_thresh Numeric. Controls how aggressively the LASSO step screens out predictors. Lower value is more aggressive. \code{lasso_thresh = 0.95}, for example, retains predictors that collectively explain at least 95% of the deviance explained by a "full" model.
#' @param xmax Integer. Maximum number of predictors returned by LASSO step. Does not strictly control the number of final predictors returned (especially for categorical \code{y} variables), but useful for setting a (very) soft upper bound. Lower \code{xmax} can help control computation time if a large number of \code{x} pass the correlation screen. \code{xmax = Inf} imposes no restriction.
#' @param xforce Character. Subset of \code{x} variables to "force" as included predictors in the results.
#' @param fraction Numeric. Fraction of observations in \code{data} to randomly sample. For larger datasets, sampling often has minimal effect on results but speeds up computation.
#' @param cores Integer. Number of cores used. Only applicable on Unix systems.
#'
#' @return List with named slots "y" and "x". Each is a list of the same length. Former gives the preferred fusion order. Latter gives the preferred sets of predictor variables.
#'
#' @examples
#' y <- names(recs)[c(14:16, 20:22)]
#' x <- names(recs)[2:13]
#'
#' # Fusion variable "blocks" are respected by prepXY()
#' y <- c(list(y[1:2]), y[-c(1:2)])
#'
#' # Do the prep work...
#' prep <- prepXY(data = recs, y = y, x = x)
#'
#' # The result can be passed to train()
#' train(data = recs, y = prep$y, x = prep$x)
#'
#' @export

#-----

# library(fusionModel)
# library(data.table)
# library(tidyverse)
# source("R/utils.R")
#
# data <- select(recs, -starts_with("rep_"))
# y <- names(recs)[c(13:16, 20:22)]
# x <- setdiff(names(data), c(y, "weight"))
# weight = "weight"
# fraction = 1
# cor_thresh = 0.05
# lasso_thresh = 0.95
# cores = 1
# xmax = 5
# xforce = NULL
#
# # Let y have a block
# y <- c(list(y[1:2]), y[-c(1:2)])
#
# test <- prepXY(data, y, x, weight = "weight")

#-----

prepXY <- function(data,
                   y,
                   x,
                   weight = NULL,
                   cor_thresh = 0.05,
                   lasso_thresh = 0.95,
                   xmax = 100,
                   xforce = NULL,
                   fraction = 1,
                   cores = 1) {

  t0 <- Sys.time()

  # Create correct 'y' and 'ylist', depending on input type
  if (is.list(y)) {
    ylist <- y
    y <- unlist(ylist)
  } else {
    ylist <- as.list(y)
  }

  stopifnot(exprs = {
    is.data.frame(data)
    all(y %in% names(data))
    all(x %in% names(data))
    length(intersect(y, x)) == 0
    all(lengths(ylist) > 0)
    is.null(weight) | weight %in% names(data)
    is.null(xforce) | all(xforce %in% x)
    xmax >= 1
    fraction > 0 & fraction <= 1
    cor_thresh > 0 & cor_thresh <= 1
    lasso_thresh > 0 & lasso_thresh <= 1
    cores >= 1 & cores %% 1 == 0
  })

  # TO DO: Make operations data.table for efficiency
  if (is.data.table(data)) data <- as.data.frame(data)

  # Check for character-type variables; stop with error if any detected
  # Check for no-variance (constant) variables
  # Detect and impute any missing values in 'x' variables
  data <- checkData(data = data, y = y, x = x, nfolds = NULL, impute = TRUE)

  # In case any variables are removed by checkData(), update 'x', 'y', and 'ylist'
  x <- intersect(x, names(data))
  for (i in 1:length(ylist)) ylist[[i]] <- intersect(ylist[[i]], names(data))
  ylist <- purrr::compact(ylist)
  y <- unlist(y)

  #---

  # Observation weights vector
  W <- if (is.null(weight)) {
    rep(1L, nrow(data))
  } else {
    data[[weight]] / mean(data[[weight]])
  }
  set(data, j = weight, value = NULL)

  #-----

  # Which 'y' variables are zero-inflated?
  yinf <- names(which(sapply(data[y], inflated)))

  if (length(yinf)) {

    # Create '*_zero' versions of the zero-inflated variables
    dinf <- data[yinf] %>%
      mutate_all(~ .x != 0) %>%
      setNames(paste0(yinf, "_zero"))  # Variable has "_zero" suffix, but it actually indicates when the original 'y' variable is NON-zero.

    # Update the zero-inflated variables in data to have NA instead of zero values
    data <- data %>%
      cbind(dinf)

    # Update 'ylist' to include '*_zero' versions in block with original zero-inflated variables
    ylist <- lapply(ylist, function(v) if (any(v %in% yinf)) c(v, paste0(intersect(v, yinf), "_zero")) else v)
    y <- unlist(ylist)

  }

  #-----

  # Sample 'data', if requested
  if (fraction < 1) {
    samp <- sample.int(n = nrow(data), size = round(nrow(data) * fraction))
    data <- data[samp, ]
    W <- W[samp]
  }

  #-----

  # Assemble the 'Z' matrix with all required variables

  X <- data[x] %>%
    mutate_if(is.ordered, as.integer) %>%
    mutate_if(is.logical, as.integer) %>%
    mutate_if(is.factor, lumpFactor, nmax = 5) %>%
    one_hot(dropUnusedLevels = TRUE)
  xlink <- attr(X, "one_hot_link")
  xcols <- names(X)

  Y <- data[y] %>%
    mutate_if(is.logical, as.integer) %>%
    mutate_if(is.factor, lumpFactor, nmax = 5) %>%
    one_hot(dropOriginal = TRUE, dropUnusedLevels = TRUE)
  ylink <- attr(Y, "one_hot_link")
  yfactor <- names(which(sapply(data[y], is.factor)))
  ycols <- names(Y)

  xylink <- rbind(xlink, ylink)
  rm(data)
  Z <- as.matrix(cbind(Y, X))  # Could make sparse?
  rm(X, Y)

  #-----

  # TO DO: Switch to using this code instead of 'X' and 'Y' blocks above
  # 8/9/24:TEST ALT
  # MOVE this could into separate function? It is also in impute()
  # d2 <- copy(data)
  #
  # # Convert 'd2' to plausible ranks for correlation screening
  # # All output columns should be NA, integer, or logical
  # # NA's in input are preserved in output
  # for (i in 1:ncol(d2)) {
  #   z <- d2[[i]]
  #   if (is.numeric(z)) {
  #     # Ties method 'dense' ensures integer output with minimum of 1 and maximum of length(na.omit(z))
  #     z <- frank(z, ties.method = "dense", na.last = "keep")
  #   } else {
  #     if (is.ordered(z)) {
  #       z <- as.integer(z)
  #     } else {
  #       if (!is.logical(z)) {
  #         # Converts character and un-ordered factors to TRUE for the most-common (non-NA) value and FALSE otherwise
  #         zt <- table2(z, na.rm = TRUE)
  #         z <- z == names(which.max(zt))
  #       }
  #     }
  #   }
  #   set(d2, j = i, value = z)
  # }

  #-----

  # intersect() call restricts to factor levels present in 'Z' (some levels can be dropped when 'data' is randomly subsampled)
  vc <- lapply(y, function(v) {
    out <- if (v %in% yfactor) {
      vl <- filter(ylink, original == v)$dummy
      intersect(vl, colnames(Z))
    } else {
      v
    }
    return(out)
  }) %>%
    setNames(y)

  #-----

  # Determine the x-predictors that pass absolute correlation threshold for each y
  # The 'Zr' matrix contains the ranks, so the correlation threshold refers to Spearman (rank) correlation
  cat("Identifying 'x' that pass absolute Spearman correlation threshold\n")
  Zr <- matrixStats::colRanks(Z, ties.method = "average", preserveShape = TRUE, useNames = TRUE)
  xok <- parallel::mclapply(unlist(vc), function(v) {

    # Initial correlation screening, based on absolute correlation value
    p <- abs(suppressWarnings(cor(Zr[, v], Zr[, xcols], use = "pairwise.complete.obs")))
    p[is.na(p)] <- 0
    vx <- which(p > cor_thresh)  # Arbitrary correlation threshold

    # Ensure some minimum number of predictors are passed to glmnet()
    # If there are too few predictors, glmnet() may fail
    if (length(vx) < 20) vx <- order(p, decreasing = TRUE)[1:min(20, length(p))]

    return(xcols[vx])

  }, mc.cores = cores) %>%
    setNames(unlist(vc))
  rm(Zr)

  #-----

  # Wrapper function for fitting a glmnet LASSO model
  # Used repeatedly in looped calls below
  gfit <- function(y, x) {
    i <- if (y %in% yinf) Z[, y] != 0 else rep(TRUE, nrow(Z))
    suppressWarnings({
      glmnet::glmnet(
        x = Z[i, x],
        y = Z[i, y],
        weights = W[i],
        family = "gaussian",
        pmax = min(xmax, length(x)),
        alpha = 1)
    })
  }

  # Testing with 'pmax' argument enabled
  # y = "square_feet"
  # x = setdiff(xcols, y)
  # m <- gfit(y, x)
  # i <- which(m$dev.ratio / max(m$dev.ratio) >= ifelse(m$jerr == 0, lasso_thresh, 1))[1] # Preferred lambda index value for each model fit, based on supplied lasso threshold
  # cf <- coef(m, s = m$lambda[i])
  # sum(cf[, 1] != 0)

  #-----

  cat("Fitting full models for each 'y'\n")

  # Weights
  ywgt <- matrixStats::colWeightedMeans(Z, W, cols = ycols)

  # Fit the "full" models for each fusion variable/block
  rmax <- parallel::mclapply(ylist, function(yvar) {
    sapply(yvar, function(v) {
      V <- vc[[v]]  # Column names of dummies in case where 'v' is a factor
      fits <- lapply(V, function(yv) gfit(y = yv, x = c(xok[[yv]], setdiff(ycols, c(yvar, V)))))
      r2 <- sapply(fits, function(m) max(m$dev.ratio))
      if (length(r2) > 1) r2 <- sum(r2 * ywgt[V]) # Only applies weighting when 'r2' contains multiple values (i.e. unique value for each factor level within a variable)
      r2
    }) %>%
      mean()  # Returns mean of R2 in case of multiple 'yvar'
  }, mc.cores = cores) %>%
    simplify2array() %>%
    setNames(ylist)

  #-----

  cat("Iteratively constructing preferred fusion order\n")

  # Start building the preferred fusion order...
  ord <- NULL  # Vector with preferred fusion variable sequence
  xpred <- NULL

  # Print loop progress to console?
  for (i in 1:length(ylist)) {

    # Candidate y variables remaining to add to 'ord'
    ycand <- setdiff(ylist, ord)

    out <- parallel::mclapply(ycand, function(yvar) {

      # This wrapper is necessary to handle cases of blocked 'yvar' with 2+ fusion variables OR case of zero-inflated fusion variable with a "_zero" version included.
      out2 <- lapply(yvar, function(v) {
        V <- vc[[v]]
        fits <- lapply(V, function(yv) {
          xx <- xok[[yv]]  # x-predictors that pass minimum correlation threshold (or best-n to meet minimum number of predictors)
          m <-  gfit(y = yv, x = c(xx, unlist(vc[unlist(ord)])))
          i <- which(m$dev.ratio / max(m$dev.ratio) >= ifelse(m$jerr == 0, lasso_thresh, 1))[1] # Preferred lambda index value for each model fit, based on supplied lasso threshold
          r2 <- m$dev.ratio[i]
          cf <- coef(m, s = m$lambda[i])
          xk <- names(which(Matrix::rowSums(cf != 0) > 0)[-1])  # Predictors with non-zero coefficients
          xx <- setdiff(xx, xk)  # Remaining zero-coefficient x-predictors, in order of correlation preference
          if (length(xk) < 20 & length(xx) > 0) xk <- c(xk, xx[1:min(20 - length(xk), length(xx))])  # Adds zero-coefficient predictors to achieve some minimum number
          list(r2 = r2, xk = xk)
        })

        # Extract R-squared for each model and calculate weighted mean across factor levels, if necessary
        r2 <- sapply(fits, function(m) m$r2)
        if (length(r2) > 1) r2 <- sum(r2 * ywgt[V])

        # Extract full set of useful x-predictors
        xk <- lapply(fits, function(m) m$xk) %>%
          unlist() %>%
          unique()

        # Return R2 and x-predictor results
        list(r2 = r2, xk = xk)

      })

      #---

      # Handle blocked 'v' case (combine results across individual variables)
      out2 <- if (length(out2) > 1) {
        list(r2 = mean(sapply(out2, function(x) x$r2)),
             xk = unique(unlist(lapply(out2, function(x) x$xk))))
      } else {
        out2[[1]]
      }

      return(out2)

    }, mc.cores = cores) %>%
      setNames(ycand)

    r2 <- sapply(out, function(m) m$r2)
    score <- r2 / rmax[names(r2)]
    best <- which.max(score)

    # Update 'ord' with best next fusion variable(s) in the chain
    ord <- c(ord, ycand[best])

    # Extract the predictor variables to be used for 'best' fusion variable
    keep <- out[[best]]$xk
    # if (!is.null(xlink)) {
    #   i <- keep %in% xlink$dummy
    #   keep[i] <- filter(xlink, dummy %in% keep[i])$original
    if (!is.null(xylink)) {
      # i <- keep %in% xylink$dummy
      # keep[i] <- filter(xylink, dummy %in% keep[i])$original
      i <- xylink$original[match(keep, xylink$dummy)]
      i[is.na(i)] <- keep[is.na(i)]
      keep <- i
    }
    keep <- unique(keep)
    #keep <- setdiff(keep, ycols)  # Remove any fusion variables from the preferred predictor set
    xpred <- c(xpred, list(keep))

  }

  #---

  # Remove any "*_zero" variables from the 'ord' result or 'xpred' results
  ord <- lapply(ord, function(v) setdiff(v, paste0(yinf, "_zero")))
  xpred <- lapply(xpred, function(v) setdiff(v, paste0(yinf, "_zero")))

  # Force inclusion of 'xforce' predictor variables
  xpred <- lapply(xpred, function(v) unique(c(v, xforce)))

  # Nicely name and order the x-predictors list in order of original 'x'
  # names(xpred) <- sapply(ord, paste, collapse = " | ")
  # xpred <- lapply(xpred, function(v) v[order(match(v, x))])

  # Results list
  result <- list(y = ord, x = xpred)

  # The full set of variables being retained, stored as attribute
  pvars <- unique(unlist(xpred))
  #stopifnot(all(pvars %in% x))
  attr(result, "xpredictors") <- intersect(x, pvars)
  attr(result, "xforce") <- xforce
  attr(result, "xoriginal") <- x  # Original, full set of potential predictor variables
  #cat("Retained", length(pvars), "of", length(x), "potential predictor variables\n")
  cat("Retained", length(intersect(x, pvars)), "of", length(x), "potential predictor variables\n")

  # Report processing time
  tout <- difftime(Sys.time(), t0)
  cat("Total processing time:", signif(as.numeric(tout), 3), attr(tout, "units"), "\n", sep = " ")

  return(result)

}
