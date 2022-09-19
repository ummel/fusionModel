#' Plot validation results
#'
#' @description
#' Creates and optionally saves to disk representative plots of validation results returned by \code{\link{validate}}. Requires the \code{\link[qgam]{qgam}}, \code{\link[ggplot2]{ggplot2}}, and \code{\link[scales]{scales}} packages.
#'
#' @param valid List returned by \code{\link{validate}}.
#' @param y Character. Fusion variables to use for validation graphics. Useful for plotting partial validation results. Default is to use all fusion variables present in \code{valid}.
#' @param path Character. Path to directory where .png graphics are to be saved. Directory is created if necessary. If NULL (default), no files are saved to disk.
#' @param cores Integer. Number of cores used. Only applicable on Unix systems.
#' @param ... Arguments passed to \code{\link[ggplot2]{ggsave}} to control .png graphics saved to disk.
#'
#' @details The validation results are visualized with the goal of conveying expected, typical performance across the fusion variables. That is, how well do the simulated data match the observed data with respect to point estimates, pairwise correlations, and uncertainty/standard errors for both the full sample and population subsets of various size? See \code{\link{validate}} for details and rationale of subset creation.
#' @details The input validation results are converted to plausible error metrics for plotting. For comparison of point estimates, the error metric is absolute percent error for continuous variables and classification error for categorical (i.e. the proportion of the population that is misclassified in the simulated data compared to the observed). Since these metrics are similar in interpretation, results for continuous and categorical variables are plotted together.
#' @details The error metric is absolute error when comparing correlation coefficients. Uncertainty estimates are compared using a t-test of equal point estimates and the ratio of simulated-to-observed standard errors. The latter allows the user to detect systematic bias in standard errors.
#' @details For a given fusion variable, the error metric will exhibit variation even for subsets of comparable size, due to the fact that each subset looks at a unique partition of the data. In order to convey how expected performance varies with subset size, an additive quantile regression model (see \code{\link[qgam]{qgam}}) is used to calculate a smoothed, robust relationship between subset size and the validation metric. The resulting curve gives the expected, typical (median) performance, conditional on subset size.
#' @details A console message indicating "outer Newton did not converge fully" van be produced (via \code{\link[qgam]{qgam}}; it cannot be suppressed). This can generally be ignored, though it is always good practice to check the variable-specific plots to ensure the smoothed fits are plausible.
#'
#' @return A list with the following named slots, each containing a \code{\link[ggplot2]{ggplot}} object visualizing validation results for all of the fusion variables.
#' @return \itemize{
#'   \item plot1: Validity of point estimates (absolute percent error).
#'   \item plot2: Validity of pairwise correlation coefficients (absolute error).
#'   \item plot3: Validity of \emph{t}-test of equality of observed and simulated point estimates (probability of p-value < 0.05).
#'   \item plot4: Validity of standard errors (ratio of simulated SE to observed SE).
#'   \item Additional named slots (one for each of the fusion variables) contain the same four plots described above with scatterplot of subset-specific results and a shaded 95% confidence interval of the smoothed median.
#'   \item est_smooth: data frame with raw smoothing results for estimates (expert use only)
#'   \item cor_smooth: data frame with raw smoothing results for correlations (expert use only)
#'  }
#'
#' @examples
#' # Build a fusion model using RECS microdata
#' # Note that "test_model.fsn" will be written to working directory
#' fusion.vars <- c("electricity", "natural_gas", "aircon")
#' predictor.vars <- names(recs)[2:12]
#' train(data = recs,
#'       y = fusion.vars,
#'       x = predictor.vars,
#'       file = "test_model.fsn",
#'       weight = "weight")
#'
#' # Generate multiple implicates
#' sim <- fuse(data = recs,
#'             file = "test_model.fsn",
#'             M = 30,
#'             ignore_self = TRUE)
#'
#' # Calculate validation results
#' valid <- validate(observed = recs[, c(fusion.vars, predictor.vars)],
#'                   simulated = sim,
#'                   sample_weights = recs$weight,
#'                   replicate_weights = subset(recs, select = rep_1:rep_96))
#'
#' # Create validation plots
#' vplots <- plot_valid(valid)
#' names(vplots)
#'
#' # View some of the plots
#' vplots$plot1
#' vplots$electricity$plot1
#' vplots$aircon$plot2
#'
#' # Can also save the plots to disk at creation
#' # Will save .png files to 'valid_plots' folder in working directory
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
#             M = 20,
#             ignore_self = TRUE)
# valid <- validate(observed = recs[, c(fusion.vars, predictor.vars)],
#                   simulated = sim,
#                   sample_weights = recs$weight,
#                   replicate_weights = subset(recs, select = rep_1:rep_96))

#------------------

# PLOTTING FUNCTION
# Generate plots of validate() output

plot_valid <- function(valid,
                       y = NULL,
                       path = NULL,
                       cores = 1,
                       ...) {

  # Check if necessary suggested packages are present
  suppressMessages(ok <- require(qgam, quietly = TRUE, warn.conflicts = FALSE))
  if (!ok) stop("package 'qgam' must be installed")
  suppressMessages(ok <- require(ggplot2, quietly = TRUE, warn.conflicts = FALSE))
  if (!ok) stop("package 'ggplot2' must be installed")
  suppressMessages(ok <- require(scales, quietly = TRUE, warn.conflicts = FALSE))
  if (!ok) stop("package 'scales' must be installed")

  t0 <- Sys.time()

  # Check inputs
  stopifnot({
    all.equal(names(valid), c("estimates", "correlation"))
    all(y %in% valid$estimates$y)
    is.null(path) | is.character(path)
    cores > 0 & cores %% 1 == 0
  })

  # Create 'directory', if necessary
  if (!is.null(path)) {
    path <- normalizePath(path, mustWork = FALSE)
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  }

  #---

  fvars <- levels(valid$estimates$y)
  se <- "se0" %in% names(valid$estimates)

  # Restrict fusion (y) variables, if requested
  if (!is.null(y)) {
    fvars <- y
    valid$estimates <- filter(valid$estimates, y %in% fvars)
    valid$correlation <- filter(valid$correlation, y %in% fvars)
  }

  #---

  # Point estimate discrepancy
  vest <- valid$estimates %>%
    mutate(
      error_pe = ifelse(is.na(level), abs((m1 - m0) / m0), abs(m1 - m0)),
    )

  # If observed standard errors are available...
  # Standard error discrepancy
  # Ratio of simulated SE to observed SE
  # Statistical test p-value
  vest <- if (se) {
    vest %>% mutate(
      se0 = ifelse(round(se0, 4) == 0, 0, se0),  # Ensures that 'se0' values very close to zero are treated as such
      #error_se = ifelse(se0 == 0, NA, abs((se1 - se0) / m0)),
      #error_se = abs((se1 - se0) / se0),
      #ratio = (se1 / m1) / (se0 / m0),
      ratio = se1 / se0,
      fail = ifelse(se0 == 0, NA, pvalue < 0.05)
    )
  } else {
    vest %>% mutate(
      ratio = NA, fail = NA
    )
  }

  # Summarize point estimate validation metrics by fusion variable and subset
  vest <- vest %>%
    mutate_if(is.numeric, ~ ifelse(!is.finite(.x), NA, .x)) %>%
    mutate(cont = is.na(level)) %>%
    group_by(iset, share, y) %>%
    mutate(wgt = 1 / n(),
           error_pe = sum(error_pe) / ifelse(cont, 1, 2)) %>%
    ungroup()
  #group_by(iset, share, y, level, cont, wgt, error_pe)
  # summarize(
  #   #cont = is.na(level[1]),
  #   #error_pe = sum(error_pe) / ifelse(cont, 1, 2),
  #   #error_se = mean(error_se),
  #   #pvalue = mean(pvalue),
  #   fail = pvalue < 0.05,
  #   ratio = mean(ratio),
  #   .groups = "drop")

  #---

  # Smooth with median via qgam
  #xseq <- data.frame(share = seq(min(est$share), 1, length.out = 100))
  xseq <- seq(min(vest$share), 1, length.out = 100)

  # QGAM FOR ESTIMATES
  cat("Smoothing validation metrics for point estimates\n")
  qest <- parallel::mclapply(fvars, function(v) {

    d <- filter(vest, y == v)

    pred.pe <- qgam2(x = d$share, y = d$error_pe, w = d$wgt, xout = xseq)[-1] %>% setNames(c("pe", "peSE"))
    pred.tt <- if (se) qgam2(x = d$share, y = d$fail, w = d$wgt, xout = xseq)[-1] %>% setNames(c("fa", "faSE")) else NULL
    pred.ra <- if (se) qgam2(x = d$share, y = d$ratio, w = d$wgt, xout = xseq)[-1] %>% setNames(c("ra", "raSE")) else NULL

    #pred.se <- if (se) qgam2(x = d$share, y = d$error_se, xout = xseq)[-1] %>% setNames(c("se", "seSE")) else NULL
    #pred.pv <- if (se) qgam2(x = d$share, y = d$pvalue, xout = xseq)[-1] %>% setNames(c("pv", "pvSE")) else NULL

    # Compile smoothed output for all validation metrics
    c(pred.pe, pred.tt, pred.ra) %>%
      as.data.table() %>%
      mutate_all(~ pmax(0, .x)) %>%  # Enforce minimum value of zero
      mutate(categorical = !d$cont[1],
             VAR = v,
             SHR = xseq,
             DFpe = sum(is.finite(d$error_pe)) - 1L,
             DFse = ifelse(se, sum(is.finite(d$ratio)) - 1L, NA))

  }, mc.cores = cores) %>%
    data.table::rbindlist()

  #---

  # Summarize correlation validation metrics by fusion variable and subset
  vcor <- valid$correlation %>%
    group_by(iset, share, y) %>%
    mutate(wgt = 1 / n()) %>%
    ungroup() %>%
    filter(is.finite(error)) %>%
    left_join(select(vest, y, cont) %>% distinct(), by = "y")

  #---

  # QGAM FOR CORRELATION
  cat("Smoothing validation metrics for correlation coefficients\n")
  qcor <- parallel::mclapply(fvars, function(v) {
    d <- filter(vcor, y == v)
    qgam2(x = d$share, y = d$error, xout = xseq)[-1] %>%
      setNames(c("ce", "ceSE")) %>%
      mutate_all(~ pmax(0, .x)) %>%  # Enforce minimum value of zero
      mutate(categorical = !d$cont[1],
             VAR = v,
             SHR = xseq,
             DF = nrow(d) - 1L)
  }, mc.cores = cores) %>%
    data.table::rbindlist()

  #----

  # Create plot objects
  cat("Creating ggplot objects", ifelse(is.null(path), "", "and saving .png files to disk"), "\n")

  # X-axis definition
  b <- c(0.01, 0.03, 0.05, 0.1, 0.2, 0.4, 0.7, 1) # Nice break marks on square root scale
  xaxis <- scale_x_continuous(name = "Subset size (percent of total population)",
                              limits = c(min(qest$SHR), 1), breaks = b, trans = "sqrt", labels = scales::percent)

  # ggplot elements to add to all plots
  plot.all <- list(xaxis,
                   theme_bw(),
                   theme(plot.title = element_text(face = "bold")),
                   guides(color = guide_legend(title = "Fusion variable")))

  # ggplot elements to add to each specific plot
  p1.add <- list(scale_y_continuous(name = "Median absolute percent error", n.breaks = 8, limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05)), labels = scales::percent),
                 labs(subtitle = "Validity of point estimates"))
  p2.add <- list(scale_y_continuous(name = "Median absolute correlation error", n.breaks = 8, limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05))),
                 labs(subtitle = "Validity of pairwise correlations"))
  # p3.add <- list(scale_y_continuous(name = "Median p-value", n.breaks = 8, limits = c(0, 1), expand = expansion(mult = c(0.01, 0.05))),
  #                labs(subtitle = bquote('Validity of t-test with '*H[0]*': observed and simulated point estimates are equal')))
  p3.add <- list(scale_y_continuous(name = "Probability of p-value < 0.05", n.breaks = 8, limits = c(0, 1), expand = expansion(mult = c(0.01, 0.05))),
                 labs(subtitle = bquote('Validity of t-test with '*H[0]*': observed and simulated point estimates are equal')))
  p4.add <- list(scale_y_continuous(name = "Median ratio of simulated SE to observed SE", n.breaks = 8, limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05))),
                 geom_hline(yintercept = 1, linetype = 2),
                 labs(subtitle = "Validity of standard errors"))
  # p5.add <- list(scale_y_continuous(name = "Median absolute percent error (relative to obs. estimate)", n.breaks = 8, limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05)), labels = scales::percent),
  #                labs(subtitle = "Validity of standard errors (comparison #2)"))

  #----

  # Multiple variable plots; no standard error shading
  p1 <- ggplot(qest, aes(x = SHR, y = pe, color = VAR, linetype = categorical)) + p1.add + plot.all + geom_line()
  p2 <- ggplot(qcor, aes(x = SHR, y = ce, color = VAR, linetype = categorical)) + p2.add + plot.all + geom_line()
  #p3 <- ggplot(qest, aes(x = SHR, y = pmin(pv, 1), color = VAR, linetype = categorical)) + p3.add + plot.all + geom_line()
  p3 <- ggplot(qest, aes(x = SHR, y = fa, color = VAR, linetype = categorical)) + p3.add + plot.all + geom_line()
  p4 <- ggplot(qest, aes(x = SHR, y = ra, color = VAR, linetype = categorical)) + p4.add + plot.all + geom_line()
  #p5 <- ggplot(qest, aes(x = SHR, y = se, color = VAR, linetype = categorical)) + p5.add + plot.all

  out1 <- list(p1, p2)
  if (se) out1 <- c(out1, list(p3, p4))
  names(out1) <- paste0("plot", 1:length(out1))

  if (!is.null(path)) {
    for (i in 1:length(out1)) {
      suppressMessages(ggsave(filename = paste0("allvars", "_plot", i, ".png"), plot = out1[[i]], path = path, ...))
    }
  }

  #---

  # Single variable plots; with standard error shading

  out2 <- lapply(fvars, function(v) {

    # Outlier factor; controls extent of extreme values to exclude from plot
    ofct <- 15

    ### Plot of typical point estimate error for a single fusion variable
    pdata <- filter(qest, VAR == v)
    ed <- filter(vest, y == v)
    p1 <- ggplot(pdata, aes(x = SHR, y = pe)) +
      geom_point(data = filter(ed, error_pe < max(pdata$pe + ofct * pdata$peSE)), aes(x = share, y = error_pe), shape = 1) +
      geom_ribbon(aes(ymin = pmax(pe - peSE * qt(0.975, DFpe), 0),
                      ymax = pmin(pe + peSE * qt(0.975, DFpe), ifelse(categorical, 1, Inf))),
                  fill = "gray", alpha = 0.25, color = "red", linetype = "dotted") +
      p1.add + plot.all + labs(title = v) + geom_line(color = "red")

    ### Plot of typical correlation error for a single fusion variable
    cdata <- filter(qcor, VAR == v)
    cd <- filter(vcor, y == v)
    p2 <- ggplot(cdata, aes(x = SHR, y = ce)) +
      geom_point(data = filter(cd, error < max(cdata$ce + ofct * cdata$ceSE)), aes(x = share, y = error), shape = 1) +
      geom_ribbon(aes(ymin = pmax(ce - ceSE * qt(0.975, DF), 0),
                      ymax = pmin(ce + ceSE * qt(0.975, DF), 2)),
                  fill = "gray", alpha = 0.25, color = "red", linetype = "dotted") +
      p2.add + plot.all + labs(title = v) + geom_line(color = "red")

    ### Plot of probability of failing t-test
    p3 <- ggplot(pdata, aes(x = SHR, y = fa)) +
      geom_point(data = filter(ed, is.finite(fail)), aes(x = share, y = as.integer(fail)), shape = 1) +
      geom_ribbon(aes(ymin = pmax(fa - faSE * qt(0.975, DFse), 0),
                      ymax = pmin(fa + faSE * qt(0.975, DFse), 1)),
                  fill = "gray", alpha = 0.25, color = "red", linetype = "dotted") +
      p3.add + plot.all + labs(title = v) + geom_line(color = "red")

    ### Plot of ratio of simulated SE to observed SE
    p4 <- ggplot(pdata, aes(x = SHR, y = ra)) +
      geom_point(data = filter(ed, ratio < max(pdata$ra + ofct * pdata$raSE)), aes(x = share, y = ratio), shape = 1) +
      geom_ribbon(aes(ymin = ra - raSE * qt(0.975, DFse),
                      ymax = ra + raSE * qt(0.975, DFse)),
                  fill = "gray", alpha = 0.25, color = "red", linetype = "dotted") +
      p4.add + plot.all + labs(title = v) + geom_line(color = "red")

    ### Plot of typical uncertainty error for a single fusion variable
    # p5 <- ggplot(pdata, aes(x = SHR, y = se)) +
    #   geom_point(data = filter(ed, error_se < max(pdata$se + ofct * pdata$seSE)), aes(x = share, y = error_se), shape = 1) +
    #   geom_ribbon(aes(ymin = pmax(se - seSE * qt(0.975, DFse), 0),
    #                   ymax = pmin(se + seSE * qt(0.975, DFse), ifelse(categorical, 1, Inf))),
    #               fill = "gray", alpha = 0.25, color = "black", linetype = "dotted") +
    #   p5.add + plot.all + labs(title = v)

    out <- list(p1, p2)
    if (se) out <- c(out, list(p3, p4))

    # Save plots to disk
    if (!is.null(path)) {
      for (i in 1:length(out)) {
        suppressMessages(ggsave(filename = paste0(v, "_plot", i, ".png"), plot = out[[i]], path = path, ...))
      }
    }

    names(out) <- paste0("plot", 1:length(out))
    return(out)

  }) %>%
    setNames(fvars)

  #-----

  # Report location of output plots, if requested
  if (!is.null(path)) cat("Plots saved to:", path, "\n")

  # Report processing time
  tout <- difftime(Sys.time(), t0)
  cat("Total processing time:", signif(as.numeric(tout), 3), attr(tout, "units"), "\n", sep = " ")

  # Assemble final result
  result <- c(out1, out2)
  result$est_smooth <- qest
  result$cor_smooth <- qcor
  return(result)

}

#-----

# Fit conditional median smooth via qgam()

# Example usage
# x <- recs$electricity
# y <- recs$square_feet
# w <- recs$weight
# N <- 1000
# xout <- seq(1, max(x), length.out = 200)
# test <- qgam2(x, y, w, N = 1000, xout = seq(min(x), max(x), length.out = 200))

qgam2 <- function(x,
                  y,
                  xout,
                  w = rep(1, length(x)),
                  N = 1000,
                  plot = FALSE) {

  # Remove non-finite values from inputs
  i <- is.finite(x) & is.finite(y) & is.finite(w)
  x <- x[i]; y <- y[i]; w <- w[i]

  d <- if (length(x) > N) {
    downsample(x, y, w, N = N)
  } else {
    data.frame(x, y, w)
  }

  # Collapse duplicates
  d <- d %>%
    group_by(x, y) %>%
    summarize(w = sum(w), .groups = "drop")

  # Fit smoothed conditional quantile(s)
  # qgam: https://mfasiolo.github.io/qgam/articles/qgam.html
  # See here for guidance: https://stats.stackexchange.com/questions/243367/smoothing-methods-for-gam-in-mgcv-package
  # And on choosing k: https://stats.stackexchange.com/questions/359568/choosing-k-in-mgcvs-gam
  # Binomial case: https://stats.stackexchange.com/questions/487135/how-to-use-a-gam-to-predict-the-probability-in-binomial-data-as-a-function-of-pr
  xterm <- ifelse(all(d$x > 0), "log(x)", "x")  # log 'x' if all values are positive
  fobj <- as.formula(paste0("y~s(", xterm, ", k = k, bs = 'tp')"))
  ok <- FALSE
  k <- -1
  while (!ok & k < min(nrow(d), 200)) {
    suppressWarnings({
      sink <- utils::capture.output({
        fit <- if (is.logical(d$y)) {
          mgcv::gam(form = fobj, data = d, weights = w, method = "REML", family = binomial("logit"))
        } else {
          qgam::qgam(form = fobj, data = d, qu = 0.5, argGam = list(weights = d$w, method = "REML"))
        }
        chk <- mgcv::k.check(fit)
        bad <- chk[2] / chk[1] > 0.8  # Check if 'k' should be increased
        if (bad) k <- 2 * (1 + chk[1]) else ok <- TRUE
      })
    })
  }

  # Generate predictions for 'xout'
  pred <- predict(object = fit, newdata = data.frame(x = xout), type = "response", se.fit = TRUE)
  pred <- as.data.frame(pred)
  pred <- setNames(pred, c("y", "se"))
  pred$x <- xout
  pred <- subset(pred, select = c(x, y, se))

  # Plot smooth, if requested
  if (plot) {
    plot(d[-3])
    lines(pred$x, pred$y, type = "l", col = 2)
  }

  return(pred)

}
