#' Plot validation results
#'
#' @description
#' Creates and optionally saves to disk representative plots of validation results returned by \code{\link{validate}}. Requires the suggested \code{\link{ggplot2}} package. This function is (by default) called within \code{\link{validate}}. Can be useful on its own to save graphics to disk or generate plots for a subset of fusion variables.
#'
#' @param valid Object returned by \code{\link{validate}}.
#' @param y Character. Fusion variables to use for validation graphics. Useful for plotting partial validation results. Default is to use all fusion variables present in \code{valid}.
#' @param path Character. Path to directory where .png graphics are to be saved. Directory is created if necessary. If NULL (default), no files are saved to disk.
#' @param cores Integer. Number of cores used. Only applicable on Unix systems.
#' @param ... Arguments passed to \code{\link[ggplot2]{ggsave}} to control .png graphics saved to disk.
#'
#' @details Validation results are visualized to convey expected, typical (median) performance of the fusion variables. That is, how well do the simulated data match the observed data with respect to point estimates and confidence intervals for population subsets of various size?
#'
#' Plausible error metrics are derived from the input validation data for plotting. For comparison of point estimates, the error metric is absolute percent error for continuous variables; in the categorical case it is absolute error scaled such that the maximum possible error is 1. Since these metrics are not strictly comparable, the all-variable plots denote categorical fusion variables with dotted lines.
#'
#' For a given fusion variable, the error metric will exhibit variation (often quite skewed) even for subsets of comparable size, due to the fact that each subset looks at a unique partition of the data. In order to convey how expected, typical performance varies with subset size, the smoothed median error conditional on subset size is approximated and plotted.
#'
#' @return A list with "plots", "smooth", and "data" slots. The "plots" slot contains the following \code{\link[ggplot2]{ggplot}} objects:
#' \itemize{
#'   \item est: Discrepancy (absolute percent error) in point estimates.
#'   \item moe: Discrepancy (absolute percent error) in 90% margin of error.
#'   \item bias: Bias (ratio of simulated-to-observed) in 90% margin of error.
#'   \item Additional named slots (one for each of the fusion variables) contain the plots described above with scatterplot results.
#'  }
#' "smooth" is a data frame with the plotting values used to produce the smoothed median plots.
#' "data" is a data frame with the complete validation results as returned by the original call to \code{\link{validate}}.
#'
#' @examples
#' # Build a fusion model using RECS microdata
#' # Note that "fusion_model.fsn" will be written to working directory
#' fusion.vars <- c("electricity", "natural_gas", "aircon")
#' predictor.vars <- names(recs)[2:12]
#' fsn.path <- train(data = recs,
#'                   y = fusion.vars,
#'                   x = predictor.vars,
#'                   weight = "weight")
#'
#' # Fuse back onto the donor data (multiple implicates)
#' sim <- fuse(data = recs,
#'             file = fsn.path,
#'             M = 30)
#'
#' # Calculate validation results but do not generate plots
#' valid <- validate(observed = recs,
#'                   implicates = sim,
#'                   subset_vars = c("income", "education", "race", "urban_rural"),
#'                   weight = "weight",
#'                   plot = FALSE)
#'
#' # Create validation plots
#' valid <- plot_valid(valid)
#'
#' # View some of the plots
#' valid$plots$est
#' valid$plots$moe
#' valid$plots$electricity$bias
#'
#' # Can also save the plots to disk at creation
#' # Will save .png files to 'valid_plots' folder in working directory
#' # Note that it is fine to pass a 'valid' object with existing $plots slot
#' # In that case, the plots are simply re-generated
#' vplots <- plot_valid(valid,
#'                      path = file.path(getwd(), "valid_plots"),
#'                      width = 8, height = 6)
#'
#' @export

#-----

# From ?validate example
# library(fusionModel)
# library(dplyr)
# library(data.table)
# source("R/utils.R")
# fusion.vars <- c("electricity", "natural_gas", "aircon")
# predictor.vars <- names(recs)[2:12]
# sim <- fuse(data = recs,
#             file = "test_model.fsn",
#             M = 40,
#             ignore_self = TRUE)
# valid <- validate(observed = recs,
#                   simulated = sim,
#                   subset_vars = c("income", "education", "race", "urban_rural"),
#                   weight = "weight",
#                   plot = FALSE)
#
# test <- plot_valid(valid)

#------------------

plot_valid <- function(valid,
                       y = NULL,
                       path = NULL,
                       cores = 1,
                       ...) {

  # Check if ggplot is installed
  suppressMessages(ok <- require(ggplot2, quietly = TRUE, warn.conflicts = FALSE))
  if (!ok) stop("package 'ggplot2' must be installed for plot_valid() to work")

  if (inherits(valid, "validate")) {
    if (!is.data.frame(valid)) valid <- valid$data  # This forces 'valid' to just the data frame with validation results
  } else {
    stop("'valid' must be an object generated by validate()")
  }

  # Check inputs
  stopifnot({
    all(y %in% valid$y)
    is.null(path) | is.character(path)
  })

  # Create 'directory', if necessary
  if (!is.null(path)) {
    path <- normalizePath(path, mustWork = FALSE)
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  }

  #---

  # Restrict fusion (y) variables, if requested
  if (!is.null(y)) valid <- filter(valid, y %in% !!y)

  #---

  # Confidence interval overlap
  # See Compare.CI() here: https://github.com/cran/synthpop/blob/master/R/compare.syn.r
  # CIoverlap(10, 20, 5, 25)
  # CIoverlap <- function(lwr_obs, upr_obs, lwr_sim, upr_sim) {
  #   L <- pmax(lwr_obs, lwr_sim)
  #   U <- pmin(upr_obs, upr_sim)
  #   0.5 * (((U - L) / (upr_obs - lwr_obs)) + ((U - L) / (upr_sim - lwr_sim)))
  # }

  #---

  errorFun <- function(obs, sim) {
    out <- suppressWarnings(log(sim / obs))  # ln(Q); will be NA if 'est.obs' is negative and Inf if it is zero
    mape <- abs((sim - obs) / obs)  # Standard MAPE can handle case when 'est.obs' is negative (still gives Inf when zero, though)
    i <- !is.finite(out)
    out[i] <- mape[i]
    out[!is.finite(out)] <- NA
    return(out)
  }

  # Calculate the error metrics for plotting
  vest <- valid %>%
    mutate(
      cont = is.na(level),
      est = abs(errorFun(obs = est.obs, sim = est.sim)),
      moe = abs(errorFun(obs = moe.obs, sim = moe.sim)),
      bias = abs(moe.sim / est.sim) / abs(moe.obs / est.obs)
    ) %>%
    group_by(id, share, y, cont) %>%
    summarize_at(c("est", "moe", "bias"), .funs = mean) %>%
    ungroup() %>%
    mutate_if(is.numeric, ~ ifelse(is.finite(.x), .x, NA)) %>%
    na.omit() # Retain only observations with complete, non-NA error metrics

  #---

  # Conditional median of the error metrics
  cat("Smoothing validation metrics\n")
  y <- unique(vest$y)
  qest <- parallel::mclapply(y, function(v) {

    d <- filter(vest, y == v)

    est <- smoothMean(x = d$share, y = d$est, verbose = FALSE)
    moe  <- smoothMean(x = d$share, y = d$moe, verbose = FALSE)
    bias <- smoothMean(x = d$share, y = d$bias, verbose = FALSE)

    out <- rbindlist(list(est = est, moe = moe, bias = bias), idcol = "metric")
    out <- out[y >= 0, ]  # Drop rows where the smoothed value goes negative (not valid)
    out[, CAT := !d$cont[1]]
    out[, VAR := v]
    setnames(out, "x", "SHR")
    return(out)

  }, mc.cores = cores) %>%
    data.table::rbindlist()

  #---

  # Restrict the smooth results to x-range available for all of the variables
  # This ensures that the plots have a same/consistent x-axis for each of the y variables

  #xrng <- range(qest$SHR)
  xrng <- qest %>%
    group_by(VAR) %>%
    summarize(min = min(SHR), max = max(SHR), .groups = "drop") %>%
    summarize(min = max(min), max = min(max)) %>%
    unlist()

  qest <- qest %>%
    filter(SHR >= xrng[1], SHR <= xrng[2])

  #---

  # Create plot objects
  cat("Creating ggplot2 graphics", ifelse(is.null(path), "", "and saving .png files to disk"), "\n")

  # X-axis definition
  b <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1) # Nice break marks on square root scale
  b <- b[max(1, which(b <= xrng[1])[1], na.rm = TRUE) : which(b >= xrng[2])[1]]
  pct <- function(x) paste0(as.character(x * 100), "%")
  xaxis <- scale_x_continuous(name = "Subset size (percent of total population)",
                              limits = xrng, breaks = b, trans = "sqrt", labels = pct)

  # ggplot elements to add to all plots
  plot.all <- list(xaxis,
                   theme_bw(),
                   theme(plot.title = element_text(face = "bold")),
                   guides(color = guide_legend(title = "Fusion variable"),
                          linetype = guide_legend(title = "Categorical")))

  #---

  # ggplot elements to add to each specific plot
  p1.add <- list(scale_y_continuous(name = "Mean absolute percent error", n.breaks = 8, limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05)), labels = pct),
                 labs(subtitle = "Discrepancy in point estimates"))

  p2.add <- list(scale_y_continuous(name = "Mean absolute percent error", n.breaks = 8, limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05)), labels = pct),
                 labs(subtitle = "Discrepancy in 90% margin of error (MOE)"))

  p3.add <- list(scale_y_continuous(name = "Mean ratio of simulated-to-observed MOE", n.breaks = 8, limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05))),
                 geom_hline(yintercept = 1, linetype = 2),
                 labs(subtitle = "Bias in 90% margin of error (MOE)"))

  #---

  # Multiple variable plots; no standard error shading
  p1 <- ggplot(filter(qest, metric == "est"), aes(x = SHR, y = y, color = VAR, linetype = CAT)) + p1.add + plot.all + geom_line()
  p2 <- ggplot(filter(qest, metric == "moe"), aes(x = SHR, y = y, color = VAR, linetype = CAT)) + p2.add + plot.all + geom_line()
  p3 <- ggplot(filter(qest, metric == "bias"), aes(x = SHR, y = y, color = VAR, linetype = CAT)) + p3.add + plot.all + geom_line()

  out1 <- list(p1, p2, p3)
  names(out1) <- c("est", "moe", "bias")

  if (!is.null(path)) {
    for (i in names(out1)) {
      suppressMessages(ggsave(filename = paste0("allvars_", i, ".png"), plot = out1[[i]], path = path, ...))
    }
  }

  #---

  # Single variable scatterplots

  out2 <- lapply(y, function(v) {

    fct <- 10  # Controls extent of extreme value inclusion in scatterplots
    pdata <- filter(qest, VAR == v)
    ed <- filter(vest, y == v, share >= xrng[1], share <= xrng[2])

    r <- median(ed$est) + c(-1, 1) * fct * mad(ed$est)
    p1 <- ggplot(filter(pdata, metric == "est"), aes(x = SHR, y = y)) +
      geom_point(data = filter(ed, est >= r[1], est <= r[2]),
                 aes(x = share, y = est), shape = 1) +
      #geom_point(data = filter(ed, error_est < max(pdata$est + ofct * pdata$peSE)), aes(x = share, y = error_pe), shape = 1) +
      # geom_ribbon(aes(ymin = pmax(pe - peSE * qt(0.975, DFpe), 0),
      #                 ymax = pmin(pe + peSE * qt(0.975, DFpe), ifelse(categorical, 1, Inf))),
      #             fill = "gray", alpha = 0.25, color = "red", linetype = "dotted") +
      p1.add + plot.all + labs(title = v) + geom_line(color = "red")

    r <- median(ed$moe) + c(-1, 1) * fct * mad(ed$moe)
    p2 <- ggplot(filter(pdata, metric == "moe"), aes(x = SHR, y = y)) +
      geom_point(data = filter(ed, moe >= r[1], moe <= r[2]),
                 aes(x = share, y = moe), shape = 1) +
      p2.add + plot.all + labs(title = v) + geom_line(color = "red")

    r <- median(ed$bias) + c(-1, 1) * fct * mad(ed$bias)
    p3 <- ggplot(filter(pdata, metric == "bias"), aes(x = SHR, y = y)) +
      geom_point(data = filter(ed, bias >= r[1], bias <= r[2]),
                 aes(x = share, y = bias), shape = 1) +
      p3.add + plot.all + labs(title = v) + geom_line(color = "red")

    out <- list(p1, p2, p3)
    names(out) <- c("est", "moe", "bias")

    # Save plots to disk
    if (!is.null(path)) {
      for (i in names(out)) {
        suppressMessages(ggsave(filename = paste0(v, "_", i, ".png"), plot = out[[i]], path = path, ...))
      }
    }

    return(out)

  }) %>%
    setNames(y)

  #-----

  # Report location of output plots, if requested
  if (!is.null(path)) cat("Plots saved to:", path, "\n")

  # Assemble final result
  result <- list(plots = c(out1, out2), smooth = qest, data = valid)
  if (!inherits(valid, "validate")) class(result) <- c("validate", class(result))
  return(result)

}
