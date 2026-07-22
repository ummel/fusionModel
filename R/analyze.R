#' Analyze Fused Microdata Implicates (Legacy)
#'
#' @description
#' Calculates point estimates and associated margins of error (MOE) for analyses
#' performed on fused or synthetic microdata. Supports calculation of means,
#' proportions, sums, counts, and medians, with optional breakdown across population subgroups.
#'
#' \strong{Legacy Notice:} This function documents and implements the point estimate
#' and variance estimation methodology used in early microdata fusion publications (e.g., 2024).
#' A newer, improved statistical approach has since been introduced in the \pkg{fusionACS}
#' package (\href{https://ummel.github.io/fusionACS/articles/methods.html}{fusionACS Methods}).
#' While \code{analyze()} is retained in \code{fusionModel} for legacy support and
#' replication purposes, users are strongly encouraged to review the updated \pkg{fusionACS}
#' workflow for new analysis pipelines.
#'
#' @param x List. A named list specifying the desired analysis type(s) and associated
#'   target variable(s). Supported analysis types include \code{"mean"}, \code{"sum"},
#'   and \code{"median"}. Example: \code{x = list(mean = c("v1", "v2"), median = "v3")}.
#'   Target variables that are factors automatically return proportions (for \code{"mean"})
#'   or counts (for \code{"sum"}). Target variables must exist in \code{implicates},
#'   \code{static}, or be created by a custom \code{fun}.
#' @param implicates Data frame (or \code{data.table}). Synthetic/fused microdata
#'   containing implicates, typically produced by \code{\link{fuse}}. Implicates must be
#'   row-stacked and identified by an integer column named \code{"M"}.
#' @param static Data frame (or \code{data.table}), optional. Static (non-synthetic)
#'   variables that do not vary across implicates. Must satisfy
#'   \code{nrow(static) == nrow(implicates) / max(implicates$M)} and match the row ordering of \code{implicates}.
#' @param weight Character, optional. Name of the primary sampling weight column in
#'   \code{static}. If \code{NULL} (default), uniform weights equal to 1 are assumed.
#' @param rep_weights Character vector, optional. Vector of replicate weight column names
#'   in \code{static}. If provided, standard errors reflect additional sampling weight
#'   uncertainty across replicates.
#' @param by Character vector, optional. Column name(s) present in \code{implicates}
#'   or \code{static} defining population subgroups for stratified estimation. If \code{NULL},
#'   analysis is performed over the full sample.
#' @param fun Function, optional. A custom transformation function applied to input data
#'   prior to analysis. Must return a \code{data.frame} containing custom derived variables.
#' @param var_scale Numeric. Scaling factor applied to unadjusted replicate weight
#'   variance, determined by the survey design. Default is \code{4} (appropriate for ACS and RECS).
#' @param cores Integer. Number of CPU cores for parallel processing. (Applicable to Unix-based systems; defaults to \code{1}).
#'
#' @details
#' Inputs are checked for consistent row dimensions and implicate structures.
#' Estimates and standard errors are computed independently for each implicate.
#' The final point estimate represents the simple mean of estimates across implicates (\eqn{M}).
#' Standard errors and degrees of freedom are pooled across implicates using Rubin's (1987) rules.
#'
#' When replicate weights (\code{rep_weights}) are provided, standard errors for each implicate
#' account for sampling design variance. Within-implicate variance is evaluated around
#' the point estimate (equivalent to setting \code{mse = TRUE} in \code{survey::svrepdesign}).
#'
#' When replicate weights are absent, within-implicate variance for means relies on Cochran's (1977)
#' ratio variance approximation (Gatz & Smith, 1995). Proportions use a weighted
#' variance formula, while medians use an asymptotic density approximation (for large \eqn{N})
#' or bootstrap sampling (for small \eqn{N}).
#'
#' @return A \code{\link[data.table]{data.table}} containing summary results grouped
#'   by any \code{by} variables. Columns include:
#'   \describe{
#'     \item{\code{by...}}{Subgroup groupings, if \code{by} was specified.}
#'     \item{\code{N}}{Number of observations per implicate in the analyzed subgroup.}
#'     \item{\code{y}}{Name of the analyzed target variable.}
#'     \item{\code{level}}{Factor level label (if target variable \code{y} was a factor/logical).}
#'     \item{\code{type}}{Metric type: \code{"mean"}, \code{"proportion"}, \code{"sum"}, \code{"count"}, or \code{"median"}.}
#'     \item{\code{est}}{Pooled point estimate across implicates.}
#'     \item{\code{moe}}{Margin of error corresponding to a 90% confidence interval.}
#'   }
#'
#' @references
#' Cochran, W. G. (1977). \emph{Sampling Techniques} (3rd ed.). John Wiley & Sons.
#'
#' Gatz, D. F., & Smith, L. (1995). The Standard Error of a Weighted Mean Concentration — I. Bootstrapping vs Other Methods. \emph{Atmospheric Environment}, 29(11), 1185–1193.
#'
#' Rubin, D. B. (1987). \emph{Multiple Imputation for Nonresponse in Surveys}. John Wiley & Sons.
#'
#' @examples
#' \dontrun{
#' library(fusionModel)
#'
#' # Build a fusion model using RECS microdata
#' fusion.vars <- c("electricity", "natural_gas", "aircon")
#' predictor.vars <- names(recs)[2:12]
#' fsn.path <- train(data = recs, y = fusion.vars, x = predictor.vars)
#'
#' # Generate 30 implicates of the 'fusion.vars' using original RECS as the recipient
#' sim <- fuse(data = recs, fsn = fsn.path, M = 30)
#' head(sim)
#'
#' # Full-sample analysis across multiple targets and metrics
#' result <- analyze(
#'   x = list(
#'     mean = c("natural_gas", "aircon"),
#'     median = "electricity",
#'     sum = c("electricity", "aircon")
#'   ),
#'   implicates = sim,
#'   weight = "weight"
#' )
#' head(result)
#'
#' # Mean electricity consumption by climate zone and urban/rural status
#' result1 <- analyze(
#'   x = list(mean = "electricity"),
#'   implicates = sim,
#'   static = recs,
#'   weight = "weight",
#'   by = c("climate", "urban_rural")
#' )
#'
#' # Subgroup analysis incorporating sample weight uncertainty via replicate weights
#' result2 <- analyze(
#'   x = list(mean = "electricity"),
#'   implicates = sim,
#'   static = recs,
#'   weight = "weight",
#'   rep_weights = paste0("rep_", 1:96),
#'   by = c("climate", "urban_rural")
#' )
#'
#' # Custom derivation function prior to analysis
#' my_fun <- function(data) {
#'   kwh_per_ft2 <- data$electricity / data$square_feet
#'   use_natural_gas <- data$natural_gas > 0
#'   data.frame(kwh_per_ft2, use_natural_gas)
#' }
#'
#' result_custom <- analyze(
#'   x = list(mean = c("kwh_per_ft2", "use_natural_gas", "electricity")),
#'   implicates = sim,
#'   static = recs,
#'   weight = "weight",
#'   fun = my_fun
#' )
#' }
#'
#' @export

analyze <- function(x,
                    implicates,
                    static = NULL,
                    weight = NULL,
                    rep_weights = NULL,
                    by = NULL,
                    fun = NULL,
                    var_scale = 4,
                    cores = 1) {

  t0 <- Sys.time()

  stopifnot({
    is.data.frame(implicates)
    is.null(static) | is.data.frame(static)
    is.null(weight) | weight %in% names(static)
    is.null(rep_weights) | all(rep_weights %in% names(static))
    is.null(by) | (all(by %in% c(names(implicates), names(static))))
    var_scale > 0
    cores > 0 & cores %% 1 == 0
  })


  sim <- data.table(implicates)
  rm(implicates)

  # Check validity of 'x' argument
  stopifnot(all(names(x) %in% c("mean", "sum", "median")))

  # Variables for which estimates will be computed
  fx <- sapply(x, class) == "formula"
  fvars <- if (is.null(fun)) {
    unique(unlist(c(x[!fx], lapply(x, all.vars))))
  } else {
    setdiff(c(names(sim), names(static)), c(by, weight, rep_weights))  # If custom function supplied, retain all variables for the moment
  }
  stopifnot(all(fvars %in% c(names(sim), names(static))))
  if (any(by %in% fvars)) stop("Analysis variables and 'by' variables must be distinct")

  #---

  # Check that input dimensions are consistent with one another
  Mimp <- max(sim$M)
  N <- nrow(sim)  / Mimp
  nM <- sim[, .N, by = M]
  stopifnot(all(nM$M %in% seq_len(Mimp)))
  stopifnot(all(nM$N == N))
  if (!is.null(static)) {
    if (nrow(static) != N) stop("The number of rows in 'static' should equal the number of rows in 'sim' per implicate (", N, ")")
  }
  if (Mimp > 1) {
    cli::cli_alert_info("Using {Mimp} implicates")
  } else {
    cli::cli_alert_warning("Only 1 implicate; returning NA for 'moe'")
  }

  #---

  # Prepare 'static' data.table
  if (is.null(static)) {
    static <- data.table(W = rep(1, N))
    cli::cli_alert_info("Assuming uniform sample weights")
  } else {
    static <- as.data.table(static)
    if (is.null(weight) & is.null(rep_weights)) {
      static[, W := 1]
      cli::cli_alert_info("Assuming uniform sample weights")
    } else {
      setnames(static, weight, "W")
      static[, W := as.double(W)]
      rm(weight)
    }
  }

  #---

  # Extract replicate weights if supplied
  # Assigns a unique replicate weight set per implicate index (M)
  WR <- NULL
  if (!is.null(rep_weights)) {
    if (length(rep_weights) < Mimp) stop("Number of replicate weights must be >= number of implicates\n")
    wkeep <- rep_weights[1:Mimp]
    cli::cli_alert_info("Using {Mimp} of {length(rep_weights)} replicate weights")
    WR <- static[, ..wkeep]
    WR <- as.double(unlist(WR))  # Make double to avoid integer overflow
  }

  #---

  # Subset and (optionally) one-hot encode 'sim'
  vsim <- c("M", intersect(names(sim), c(by, fvars)))  # Can include 'by' variables? Why not?
  sim <- sim[, ..vsim]

  # Subset and (optionally) one-hot encode 'static'
  vstatic <- c("W", intersect(names(static), c(by, fvars)))
  vstatic <- setdiff(vstatic, vsim)  # Prioritizes variables in 'sim' if names are duplicated
  static <- static[, ..vstatic]

  # Combine static and simulated data
  sim <- cbind(sim, static[rep(1:N, Mimp)])

  # If replicate weight exist, add them to 'sim'; different set of weights for each implicate
  if (is.null(WR)) sim[, WR := W] else sim[, WR := WR]

  rm(static, WR)

  #---

  # If a custom function is supplied, apply it to 'sim' before proceeding
  # 'fun' returns a data frame that is then cbind'd to original 'sim'
  if (!is.null(fun)) {
    cli::cli_alert_info("Applying custom 'fun' to data")
    fout <- fun(sim)
    if (!is.data.frame(fout)) stop("Custom 'fun' must return a data.frame")
    drop <- intersect(names(fout), names(sim))
    if (length(drop)) sim[, c(drop) := NULL]
    sim <- cbind(sim, fout)
  }

  #---

  # Variables to retain prior to potential one-hot encoding
  # This ensures that only the variables specified in 'x' are retained
  fvars <- unique(unlist(c(x[!fx], lapply(x, all.vars))))
  stopifnot(all(fvars %in% names(sim)))
  keep <- c("M", "W", "WR", by, fvars)
  sim <- sim[, ..keep]

  # Convert logical targets to factors so they are processed as categorical indicators (proportions/counts)
  flgl <- names(which(sapply(sim[, ..fvars], is.logical)))
  for (v in flgl) set(sim, j = v, value = as.factor(sim[[v]]))

  #---

  # Check that 'fvars' for median analysis are all numeric
  check <- x$median
  if (length(check)) {
    if (!all(sapply(sim[, ..check], is.numeric))) stop("Variables for median analyses must be numeric")
  }

  #---

  # Create one-hot encoded dummy variables for factor targets requested under "mean" or "sum"
  hotf <- names(x) %in% c("mean", "sum")
  othx <- unique(unlist(lapply(x[!hotf], function(x) if (is.vector(x)) x else all.vars(x))))
  hotx <- unique(unlist(x[hotf]))
  hotx <- names(which(sapply(sim[, ..hotx], is.factor)))
  if (length(hotx) > 0) {
    temp <- one_hot(sim[, ..hotx], dropOriginal = TRUE)
    sim <- cbind(sim, temp)
    attr(sim, "one_hot_link") <- attr(temp, "one_hot_link")
    rm(temp)
  }

  # Look up table for variables that are one-hot encoded
  flink <- attr(sim, "one_hot_link")

  # Which 'hotx' original columns can be dropped?
  # This occurs if a 'hotx' variable is not requested by any other analysis
  drop <- setdiff(hotx, othx)
  if (length(drop) > 0) sim[, c(drop) := NULL]

  # TO DO: Print status update to console about which variables being used...

  # Key the data.tables for speed in 'by' operations below
  setkeyv(sim, c("M", by))

  #---

  # Function to compute confidence interval bounds
  # p = 0.95 entails a 90% confidence interval
  # calcCI <- function(d, p = 0.95) {
  #   d %>%
  #     mutate(
  #       lwr = est - qt(p, df) * se, # CI lower bound
  #       upr = est + qt(p, df) * se  # CI upper bound
  #     )
  # }

  # Function to compute margin of error (MOE)
  # p = 0.95 entails a 90% confidence interval
  calcMOE <- function(d, p = 0.95) {
    d %>%
      mutate(
        moe = se * qt(p, df)
      )
  }

  #-----

  # Process a particular combination of subsetting variables

  #processSubset <- function(iset) {

  # svar <- scomb[[iset]]
  #
  # # Calculate the share of observations within each subset, using the observed data
  # subshr <- obs[, .(share = .N / N), by = svar]
  # subshr <- subshr[share >= min_size / N]
  # subshr$id <- paste(iset, 1:nrow(subshr), sep = "_")
  #
  # # Restrict 'obs' and 'sim' to subsets with at least 'min_size' observations
  # obs <- obs[subshr, on = svar]
  # sim <- sim[subshr, on = svar]
  #
  #---

  # MEAN
  y <- x$mean
  out.mean <- if (length(y)) {
    yf <- intersect(y, hotx)
    if (length(yf)) {
      y <- setdiff(y, yf)
      y <- c(y, filter(flink, original %in% yf)$dummy)
    }
    out <- sim[, lapply(.SD, FUN = weighted_mean, w1 = W, w2 = WR), by = c("M", by), .SDcols = y]
    out[, type := "mean"]
    out
  } else {
    NULL
  }

  # SUM
  y <- x$sum
  out.sum <- if (length(y)) {
    yf <- intersect(y, hotx)
    if (length(yf)) {
      y <- setdiff(y, yf)
      y <- c(y, filter(flink, original %in% yf)$dummy)
    }
    out <- sim[, lapply(.SD, FUN = weighted_sum, w1 = W, w2 = WR), by = c("M", by), .SDcols = y]
    out[, type := "sum"]
    out
  } else {
    NULL
  }

  # MEDIAN
  y <- x$median
  out.median <- if (length(y)) {
    out <- sim[, lapply(.SD, FUN = weighted_median, w1 = W, w2 = WR), by = c("M", by), .SDcols = y]
    out[, type := "median"]
    out
  } else {
    NULL
  }

  # GLM (TO DO)

  # Combine analysis results
  d <- rbindlist(list(out.mean, out.sum, out.median), fill = TRUE)

  #----

  d$metric <- rep(c("estimate1", "estimate2", "variance"), times = nrow(d) / 3)
  d <- melt(d, id.vars = c("M", "type", "metric", by), variable.name = "y", na.rm = TRUE)
  d <- dcast(d, ... ~ metric, value.var = "value")

  # Pool point estimates and within/between-implicate variances across implicates
  # est: Point estimates (mean of estimates across implicates)
  # ubar: Mean of the within-implicate variances
  # b1: Variance of estimates across implicates (implicate variance)
  # b2: Additional variance component capturing sample weight uncertainty from replicate weights
  d <- d[, .(est = mean(estimate1),
             ubar = mean(variance),
             b1 = var(estimate1),
             #b1 = (sum(estimate1 ^ 2 + variance) / Mimp) - mean(estimate1) ^ 2,  # Conservative approximation of var(estimate1); does not risk b1 < ubar
             b2 = pmax(0, var(estimate2) - var(estimate1)) * var_scale * (Mimp - 1) / Mimp),
         by = c(by, "type", "y")]

  # Calculate standard error, degrees of freedom, and margin of error using Rubin's (1987) rules
  # Note that if ubar = 0, then b = 0 necessarily; the maxr adjustment results in b = NA; ifelse() below forces b = 0 in this case
  # 'df' is NA when se = 0, so it is forced to 1 so CI calculation does not produce error
  #maxr <- maxr_fun(Mimp)
  d <- mutate(d,
              b = b1 + b2,
              # b = ifelse(ubar / b > maxr, ubar / maxr, b),
              # b = ifelse(ubar == 0, 0, b),

              # Rubin 1987
              se = sqrt(ubar + (1 + Mimp^(-1)) * b),
              r = (1+Mimp^(-1))*b/ubar,
              df = (Mimp-1)*(1+r^(-1))^2,
              r = NULL

              # Reiter and Raghunathan (2007)
              # se = sqrt(b * (1 + 1 / Mimp) - ubar),
              # df = (Mimp - 1) * (1 - Mimp * ubar / ((Mimp + 1) * b)) ^ 2,
              # df = ifelse(se == 0, 1, df)
  ) %>%
    calcMOE()

  #---

  # Process each subset combination in 'scomb'
  # cat("Processing validation analyses for", yN, "fusion variables\n")
  # comp <- parallel::mclapply(seq_along(scomb), processSubset, mc.cores = cores) %>%
  #   rbindlist()

  # Number of observations in each group
  # The total number of observations must sum to nrow(static)
  out.N <- sim[, list(N = as.integer(.N / Mimp)), by = by]

  # Add the number of observations in each group
  if (is.null(by)) {
    d$N <- as.integer(out.N)
  } else {
    d <- d[out.N, on = by]
  }

  #---

  # Map one-hot encoded dummy outputs back to original factor column names and levels
  if (is.data.frame(flink)) {
    d <- merge(x = d, y = flink, by.x = "y", by.y = "dummy", all.x = TRUE, sort = FALSE)
    d <- d %>%
      mutate(y = ifelse(is.na(level), y, original),
             type = ifelse(is.na(level), type, ifelse(type == "mean", "proportion", "count"))) %>%
      select(-original)
  } else {
    d$level <- NA
  }

  # # Enforce lower and upper CI bounds on nominal variables that report proportions
  # if (fun == "mean") {
  #   d <- d %>%
  #     mutate(lwr = ifelse(is.na(level), lwr, pmax(0, lwr)),
  #            upr = ifelse(is.na(level), upr, pmin(1, upr)))
  # }

  # Enforce lower and upper CI bounds on nominal variables that report sums
  # if (fun == "sum") {
  #   d <- d %>%
  #     mutate(lwr = ifelse(is.na(level), lwr, pmax(0, lwr)),
  #            upr = ifelse(is.na(level), upr, pmin(wtotal, upr)))
  # }

  # Format final output structure
  keep <- c(by, "N", "y", "level", "type", "est", "moe")
  result <- d[, ..keep]
  result$y <- as.character(result$y)
  setorderv(result, c(by, "y", "type"))
  # class(result) <- c("validate", class(result))

  # Report processing time
  tout <- difftime(Sys.time(), t0)
  cli::cli_alert_success("Total processing time: {signif(as.numeric(tout), 3)} {attr(tout, 'units')}")

  return(result)

}

#------------------
#------------------

# Below are helper functions for weighted statistics

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

      # Continuous case variance calculation via ratio approximation (Cochran, 1977)
      n <- length(x)
      xWbar <- wmean(x, W)
      wbar <- mean(W)
      wvar <- n/((n-1)*sum(W)^2)*(sum((W*x-wbar*xWbar)^2)-2*xWbar*sum((W-wbar)*(W*x-wbar*xWbar))+xWbar^2*sum((W-wbar)^2))

    } else {

      # Binary case variance calculation for weighted proportions
      w <- W / sum(W)
      xWbar <- sum(w * x)
      sw2 <- sum(w ^ 2)
      wvar <- xWbar * (1 - xWbar) * sw2

    }

    # Secondary (replicate) weighted mean
    xWbar2 <- wmean(x, w2)

    # Return primary weighted mean, secondary (replicate) weighted mean, and primary within-sample variance
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

  out <- if (n > 1 & var(x) > 0) {

    # Median estimate using primary weights 'w1'
    med1 <- matrixStats::weightedMedian(x, w1)

    # Median estimate using replicate weights 'w2'
    med2 <- matrixStats::weightedMedian(x, w2)

    # Calculate variance: density approximation for large sample sizes, bootstrapping for small sample sizes
    wvar <- if (n > 1000 & med1 != 0) {

      # Approximate variance of sample median for large N
      pdf <- density(x, weights = w1 / sum(w1))
      fu <- approx(pdf$x, pdf$y, xout = med1)$y
      1 / (4 * n * fu ^ 2)

    } else {

      # Bootstrap sampling for small sample sizes
      bdraw <- replicate(n = nboot, sample(x, size = n, replace = TRUE))

      # Bootstrap weighted median estimates using primary weights 'w1'
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
