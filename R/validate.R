#' Validate fusion output
#'
#' @description
#' Performs specific (non-general) internal validation exercises on fused microdata to estimate how well the simulated data reflect statistical patterns in an analogous observed dataset. This provides a standard approach to validating a fusion model created by \code{\link{train}}. See Examples for recommended usage.
#'
#' @param observed Data frame. Observed data against which to validate the \code{simulated} variables. Typically the same dataset used to \code{\link{train}} the fusion model that one seeks to validate. Need not be limited to fusion variables (see Details).
#' @param simulated Data frame. Simulated fusion variables. Typically the output from \code{\link{fuse}} containing multiple implicates. The implicates should be row-stacked and identified by integer column "M".
#' @param sample_weights Numeric. Vector of primary sampling weights with length equal to \code{nrow(observed)}.
#' @param replicate_weights Data frame. Each column is a vector of survey replicate weights. \code{nrow(replicate_weights)} must equal \code{length(sample_weights)}.
#' @param nsets Integer. Number of population subsets to analyze. A higher number generally produces more reliable validation results but takes longer to compute. See Details.
#' @param var_scale Scalar. Factor by which to scale the standard variance in the presence of replicate weights. This is determined by the survey design associated with \code{observed}. The default (\code{var.scale = 4}) is appropriate for both RECS and ACS.
#' @param min_size Integer. Minimum number of observations allowed in a subset.
#' @param cores Integer. Number of cores used. Only applicable on Unix systems.
#'
#' @details The objective of \code{validate()} is to reveal the utility of the synthetic data across myriad analyses. Utility here is based on comparison of point estimates and standard errors derived using multiple-implicate synthetic data with those derived using the original donor data and replicate weights. The specific analyses tested include variable levels (means and proportions) and bivariate relationships (Pearson correlation coefficients) for both the full sample and relevant population subsets of varying size. This allows us to estimate how each of the synthetic variables perform in analyses with real-world relevance, at varying levels of complexity.
#' @details In effect, \code{validate()} performs a large number of analyses of the kind that the \code{\link{analyze}} function is designed to do on a one-by-one basis. Being purpose-built for particular kinds of analysis, \code{validate()} does this much more efficiently than repeated calls to \code{\link{analyze}}. It carries out hundreds or thousands of potential analyses separately for both the synthetic and observed data.
#' @details We assume that users are most likely to analyze a variable across subsets of the population that produce meaningful differences in outcomes. No (sane) analyst is interested in the mean of Y for a random subset, since the expected value within and without the subset is identical. But it does make sense to subset the population into (for example) homeowners and renters, if there is reason to believe that Y varies between the two groups.
#' @details Correlations among the provided variables are used to prioritize subsetting variables with larger absolute correlations with the fusion variables (i.e. subsets of likely real-world relevance). In addition, subsetting variables are selected so as to over-represent smaller subsets. The universe of subsets containing 95% of observations exhibit comparatively little variation compared to those containing just 5%; the latter subsets result in noisier validation results. Consequently, the returned analyses are biased toward smaller subsets in order to derive more reliable estimates of small-sample utility.
#' @details \code{validate()} attempts to create \code{nsets} total subsets. The mean/proportions of each fusion variable are returned for each of the \code{nsets} subsets, as is the pairwise (bivariate) Pearson correlation coefficient error between each fusion variable and all other variables provided in \code{observed}. See Value for details.
#' @details Note that less than \code{nsets} subsets may be returned if the provided variables do not allow for it. For purposes of subset creation, each continuous variable is converted to ten cumulative binary variables on the basis of a univariate \code{\link[stats]{kmeans}} clustering.
#'
#' @return A list with slots named "estimates" and "correlation". This data does not usually need to be interrogated, since users are more likely to pass the result to \code{\link{plot_valid}} for visualization.
#' @return \itemize{
#'   \item \emph{estimates}: Data frame giving observed (m0, se0) and synthetic (m1, se1) point estimates and standard errors for each fusion variable (y) and population subset (iset). Column 'pvalue' gives the p-value associated with a statistical test where the null hypothesis is that 'm0' and 'm1' are identical.
#'   \item \emph{correlation}: Data frame summarizing the difference (error) between synthetic and observed data for pairwise correlation coefficients of each fusion variable (y) and all other variables (correlate). Column 'error' gives the mean absolute error; 'error_w' gives the mean absolute error weighted by the square of the observed correlation (i.e. weighted by variance explained). In the case of categorical 'y' or 'correlate', the underlying correlations are calculated using binary vectors and the returned results reflect the mean absolute error across all levels.
#'}
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
#' # Fuse back onto the donor data (multiple implicates)
#' # Since the result is intended for validation, we set 'ignore_self = TRUE'
#' sim <- fuse(data = recs,
#'             file = "test_model.fsn",
#'             M = 20,
#'             ignore_self = TRUE)
#'
#' # Calculate validation results
#' valid <- validate(observed = recs[, c(fusion.vars, predictor.vars)],
#'                   simulated = sim,
#'                   sample_weights = recs$weight,
#'                   replicate_weights = subset(recs, select = rep_1:rep_96))
#'
#' @export

#-----

# # RECS EXAMPLE FOR DOCS
# library(fusionModel)
# fusion.vars <- c("electricity", "natural_gas", "aircon")
# predictor.vars <- names(recs)[2:12]
# train(data = recs, y = fusion.vars, x = predictor.vars, file = "test_model.fsn", weight = "weight")
#
# observed <- recs[, c(fusion.vars, predictor.vars)]
# simulated <- fuse(recs, file = "test_model.fsn", M = 20, cores = 2)
#
# test <- validate(observed = recs[, c(fusion.vars, predictor.vars)],
#                  simulated = simulated,
#                  sample_weights = recs$weight,
#                  replicate_weights = subset(recs, select = rep_1:rep_96),
#                  nsets = 100,
#                  var_scale = 4,
#                  cores = 2)

#-----

# library(dplyr)
# library(data.table)
# library(purrr)
# library(ggplot2)
# source("R/utils.R")
#
# observed <- recs[, c(fusion.vars, predictor.vars)]
# simulated <- fuse(recs, file = "test_model.fsn", M = 20, cores = 2)
#
# sample_weights = recs$weight
# replicate_weights = dplyr::select(recs, starts_with("rep_"))
# nsets = 100
# var_scale = 4
# cores = 2
#
# sample_weights = NULL
# replicate_weights = NULL

#-----

validate <- function(observed,
                     simulated,
                     sample_weights = NULL,
                     replicate_weights = NULL,
                     nsets = 200,
                     var_scale = 4,
                     min_size = 20,
                     cores = 1) {

  t0 <- Sys.time()

  # Check inputs
  stopifnot({
    is.data.frame(observed)
    is.data.frame(simulated)
    is.data.frame(replicate_weights) | is.matrix(replicate_weights)
    all(setdiff(names(simulated), "M") %in% names(observed))
    is.null(sample_weights) | any(is.numeric(sample_weights))
    is.null(replicate_weights) | any(ncol(replicate_weights) >= 2)
    var_scale > 0
    cores > 0 & cores %% 1 == 0
  })

  # Check that input dimensions are consistent with one another
  N <- nrow(observed)
  mcnt <- max(simulated$M)
  tb <- table(simulated$M)
  stopifnot(N == nrow(simulated) / mcnt)
  stopifnot(all(names(tb) %in% seq_len(mcnt)))
  stopifnot(all(tb == N))
  if (!is.null(observed) & any(nrow(observed) != N)) stop("'observed' has incorrect number of rows")
  w <- if (is.null(sample_weights)) {
    cat("Assuming uniform 'sample_weights`\n")
    rep(1L, N)
  } else {
    sample_weights / mean(sample_weights)
  }
  if (length(w) != N) stop("'sample_weights' is incorrect length")
  if (is.null(replicate_weights)) {
    cat("No 'replicate_weights' provided; observed standard errors omitted\n")
  } else {
    if (nrow(replicate_weights) != N) stop("'replicate_weights' has incorrect number of rows")
    if(!is.matrix(replicate_weights)) replicate_weights <- as.matrix(replicate_weights)
  }
  wr <- replicate_weights
  rm(sample_weights, replicate_weights)

  # Remove any ".." in input variable names to avoid conflict with one_hot() dummy output names
  # NOTE: This will cause variable names in output to (potentially) not match input (TO DO: fix)
  names(observed) <- gsub("..", "", names(observed), fixed = TRUE)
  names(simulated) <- gsub("..", "", names(simulated), fixed = TRUE)

  # Indices of simulation implicates to iterate over
  # TO DO: Make smaller/more efficient
  Mind <- sapply(unique(simulated$M), function(i) which(simulated$M == i))

  # Identify factor/logical variables for which to drop smallest level
  # Function to drop smallest factor level for binary variables
  #obs.fct <- names(which(sapply(observed, is.factor) | sapply(observed, is.logical)))
  # ffun <- function(x) {
  #   lev <- if (is.factor(x)) levels(x) else unique(x)
  #   if (length(lev) == 2) {
  #     factor(x, levels = setdiff(lev, names(sort(table(x)))[1]), ordered = FALSE)
  #   } else {
  #     x
  #   }
  # }

  # NECESSARY????
  # Function to drop smallest factor level for binary variables
  # obs.bin <- names(which(purrr::map_lgl(observed, ~ length(unique(.x)) == 2)))
  # binFun <- function(x) {
  #   lev <- if (is.factor(x)) levels(x) else unique(x)
  #   factor(x, levels = setdiff(lev, names(sort(table(x)))[1]), ordered = FALSE)
  # }

  #---

  # Create the observed data matrix
  x <- observed %>%
    select(-any_of(colnames(wr))) %>%
    #mutate_at(obs.bin, binFun) %>%
    select(any_of(names(simulated)), everything()) %>%
    one_hot(sparse_matrix = FALSE) %>%
    as.matrix()

  # Simulation data matrix
  y <- simulated %>%
    one_hot(sparse_matrix = FALSE) %>%
    select(all_of(intersect(colnames(x), names(.)))) %>%
    as.matrix()

  rm(simulated)

  #-----

  fvars <- colnames(y)

  # The continuous (non-binary) fusion variables
  fcont <- names(which(apply(x[, fvars], 2, function(x) any(!x %in% c(0, 1)))))

  # The binary (non-continuous) fusion variables
  fbin <- setdiff(fvars, fcont)

  # The non-fusion (other) variables
  ovars <- setdiff(colnames(x), fvars)

  # Safety check
  stopifnot(all(colnames(y) %in% colnames(x)))

  # Number of validation observations
  #N <- nrow(x)

  #-----

  # Confirm
  #tapply(X, cl, mean)
  #plot(X, cl)
  uniCluster <- function(x, k) {
    k <- min(k, length(unique(x)))
    iter <- 0
    halt <- FALSE
    while (!halt) {
      iter <- iter + 1
      km <- suppressWarnings(kmeans(x, centers = k, algorithm = "Hartigan-Wong", nstart = iter))
      halt <- km$ifault == 0
    }
    cl <- data.table::frank(km$centers, ties.method = "dense")[km$cluster]
    cl[cl == k] <- NA
    cl <- factor(cl, levels = 1:(k - 1), ordered = TRUE)  # Drops the final level (set to NA)
    return(cl)
  }

  #---

  # CLEAN THIS UP...redundant with generation in 'x' above...
  # Create 'x2' with continuous observed variables converted to discretixed dummies
  obs.num <- names(which(sapply(observed, is.numeric)))
  obs.ord <- c(obs.num, names(which(sapply(observed, is.ordered))))
  x2 <- observed %>%
    select(-any_of(colnames(wr))) %>%
    #mutate_at(obs.bin, binFun) %>%  # This is repetitive with creation of 'x'; make cleaner
    mutate_at(obs.num, uniCluster, k = 10) %>%  # TO DO: Set 'k' automatically
    #select(any_of(names(simulated)), everything()) %>%
    one_hot(sparse_matrix = FALSE) %>%
    as.matrix()

  # Convert 'obs.ord' variables to cumulative binary variables
  for (v in obs.ord) {
    j <- grep(paste0("^", v, "\\.\\."), colnames(x2))
    if (length(j) > 1) x2[, j] <- matrixStats::rowCummaxs(x2[, j])
  }

  # Restrict 'x2' to subsetting variables that meet minimum subset size requirement
  ssize <- matrixStats::colSums2(x2)
  ok <- which(ssize >= min_size & ssize <= (N - min_size))
  x2 <- x2[, ok]
  ssize <- matrixStats::colSums2(x2 * w) / sum(w)

  # What is the correlation between fusion variables and the subset variables?
  x2v0 <- purrr::map_chr(strsplit(colnames(x2), "..", fixed = TRUE), 1)
  xcor <- matrixStats::allocMatrix(nrow = length(fvars), ncol = ncol(x2), value = NA)
  for (i in seq_along(fvars)) {
    v <- fvars[[i]]
    v0 <- strsplit(v, "..", fixed = TRUE)[[1]][1]
    j <- x2v0 != v0
    xcor[i, j] <- suppressWarnings(cor(x[, v], x2[, j]))
  }

  # How to summarize correlation of each subsetting variable?
  xcor <- matrixStats::colMeans2(abs(xcor), na.rm = TRUE)
  #xcor <- matrixStats::colMaxs(abs(xcor), na.rm = TRUE)  # Could use mean or max or median
  #xcor <- matrixStats::colMedians(abs(xcor), na.rm = TRUE)  # Could use mean or max or median
  xcor[is.na(xcor)] <- 0
  stopifnot(length(xcor) == ncol(x2))

  vint <- nsets / 20
  temp <- data.frame(xv = rep(colnames(x2), 2),
                     cor = rep(xcor, 2),
                     size = c(ssize, 1 - ssize),
                     id = rep(c(1, 0), each = length(ssize))) %>%
    mutate(ventile = findInterval(1 - size ^ 0.30, seq(0, 1, length.out = 21)),
           ventile = factor(ventile, levels = 1:20)) %>%
    arrange(-cor) %>%
    group_by(ventile) %>%
    slice(1:min(n(), vint))

  # Create logical matrix indicating the subsets to analyze (smat)
  #smat <- matrixStats::allocMatrix(nrow = N, ncol = nsets, value = FALSE)

  # Create logical matrix indicating the subsets to analyze
  #for (i in 1:nrow(temp)) smat[, i] <- x2[, temp$xv[i]] == temp$id[i]]
  smat <- sapply(1:nrow(temp), function(i) x2[, temp$xv[i]] == temp$id[i])
  #cat("Selected", ncol(smat), "empirical subsets (\n")
  cat("Selected ", ncol(smat), " empirical subsets (", nsets, " requested)\n", sep = "")

  # # Add random subsets if there are too few in 'smat'
  # if (ncol(smat) < nsets) {
  #   rfill <- vint - table(temp$ventile)
  #   rfill <- rfill[rfill > 0]
  #   srnd <- lapply(1:length(rfill), function(i) {
  #     j <- rfill[i]
  #     k <- as.integer(names(rfill)[i])
  #     r <- qexp(seq(0.95, 0, length.out = 20)[k] + runif(j, 0, 0.05), rate = 5)
  #     r <- pmax(r, min_size / N)
  #     sapply(r, function(x) sample(c(T, F), size = N, prob = c(x, 1 - x), replace = TRUE))
  #   })
  #   srnd <- do.call(cbind, srnd)
  #   smat <- cbind(smat, srnd)
  #   cat("Generated", ncol(srnd), "additional random subsets\n")
  #   rm(srnd)
  # }

  # Probability of selection
  # sub_size = 0.2  # INPUT
  # rate <- 1 / sub_size
  # ep <- dexp(ssize, rate) / rate
  # P <- rep(xcor, 2) * c(ep, 1 - ep)
  # #plot(c(ssize, 1 - ssize), P)
  # #plot(ssize, ep)
  # cat("Identified", length(P), "candidate subsets\n")

  # Sample a subsetting variable in 'x2' -- NOTE: Don't love the restriction of nsets...
  # nsets <- min(nsets, length(P))
  # z <- sample.int(length(P), size = nsets, prob = P, replace = FALSE)

  # Order the columns in 'smat' from smallest subset to largest
  smat <- smat[, order(matrixStats::colSums2(smat * w), decreasing = TRUE)]

  # Force largest subset to include full sample
  smat[, 1] <- TRUE

  # Distribution of subset sizes
  #summary(colSums(smat) / nrow(smat))
  #hist(colSums(smat) / nrow(smat))

  # Update 'nsets'
  nsets <- ncol(smat)

  #stopifnot(ncol(smat) == nsets)
  rm(observed, x2)

  #-----

  # Correlation template

  # Template data.frame used for assigning variable names and other auxiliary info to correlation results
  template <- expand.grid(v1 = fvars, v2 = c(fvars, ovars), stringsAsFactors = FALSE) %>%
    mutate(n = 1:n(),
           uni = purrr::map2_chr(v1, v2, ~ paste(sort(c(.x, .y)), collapse = "--"))) %>%
    distinct(uni, .keep_all = TRUE) %>%
    filter(v1 != v2) %>%  # Remove self-correlation (i.e. r = 1)
    mutate(v1b = purrr::map_chr(strsplit(v1, "..", fixed = TRUE), 1),
           v2b = purrr::map_chr(strsplit(v2, "..", fixed = TRUE), 1)) %>%
    group_by(v1b, v2b) %>%
    mutate(vwgt = 1 / n()) %>%  # This enables equal weighting of each fusion variable, regardless of number of factor levels
    group_by(v1b) %>%
    mutate(vwgt = vwgt / sum(vwgt)) %>%  # This enables equal weighting of each fusion variable, regardless of number of factor levels
    ungroup()

  # Index for subsetting correlation results with replicate expression (below)
  keep <- template$n

  # Check that total weight is same across all fusion variables (v1b)
  # template %>%
  #   group_by(v1b) %>%
  #   summarize(total_vwgt = sum(vwgt))

  #-----

  # Standard error of a weighted mean for a continuous variable (SE's for proportions are calculated differently)
  # Code: https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
  # Computes the variance of a weighted mean following Cochran 1977 definition
  # Original paper: https://www.sciencedirect.com/science/article/abs/pii/135223109400210C
  weightedSE <- function(x, w) {
    #w <- w / mean(w)  # To prevent possible overflow
    n <- length(x)
    xWbar <- matrixStats::weightedMean(x, w)
    wbar <- mean(w)
    var.se <- n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
    return(sqrt(var.se))
  }

  pooledSE <- function(est, se) {
    stopifnot(identical(dim(est), dim(se)))
    m <- ncol(est)
    ubar <- matrixStats::rowMeans2(se ^ 2)  # Within-implicate variance of estimate ('Ubar')
    b <- matrixStats::rowVars(est)
    t <- ubar + (1 + m ^ (-1)) * b  # Total variance, of estimate; the square root of 't' is the standard error (see below)
    r <- (1 + m ^ (-1)) * b / ubar
    degf <- (m - 1) * (1 + r ^ (-1)) ^ 2  # Rubin's formula assuming infinite sample size (TO DO: Reiter 2007 correction below for finite sample)
    se <- sqrt(t)
    return(cbind(se, degf))
  }

  #-----

  # Perform calculations over 'nsets' iterations
  cat("Calculating validation results for", nsets, "subsets\n")
  comp <- parallel::mclapply(1:nsets, function(iset) {

    # Calculate point estimates for training and simulated data for a random subset
    # TO DO: Move this out of function?
    #ss <- sample.int(N, S, replace = TRUE)
    ss <- which(smat[, iset])

    # Pre-compute subsets
    ws <- w[ss]
    wrs <- wr[ss, ]
    xs <- x[ss, ]
    xso <- xs[, ovars]
    Yind <- Mind[ss, ]

    #----------

    # Observed weighted mean, ignoring replicate weights
    # This is the point estimate as typically calculated for surveys with replicate weights; i.e. equivalent to mse = TRUE in ?svrepdesign
    # Suitable for both continuous and binary cases
    m0 <- matrixStats::colWeightedMeans(x = xs[, fvars], w = ws)
    se0 <- NULL
    df0 <- NULL

    # Calculate SE of observed fusion variables using replicate weights
    if (!is.null(wrs)) {
      wrs.tot <- matrixStats::colSums2(wrs)
      se0 <- sapply(fvars, function(v) {
        mu <- matrixStats::colSums2(wrs * xs[, v]) / wrs.tot  # Mean values across replicates
        sqrt(var_scale * sum((mu - m0[v]) ^ 2) / length(mu)) # Equivalent to mse = TRUE in ?svrepdesign
        #sqrt(var_scale * sum((mu - mean(mu)) ^ 2) / length(mu))  # Equivalent to mse = FALSE in ?svrepdesign
      })
      df0 <- length(ss) - 1
    }

    #---

    # Simulated weighted mean of the fusion variables
    # Suitable for both continuous and binary cases
    sws2 <- sum((ws / sum(ws)) ^ 2)
    temp <- 1:ncol(Mind) %>%
      sapply(function(i) {

        ind <- Yind[, i]
        mu <- matrixStats::colWeightedMeans(x = y[ind, fvars], w = ws)
        se <- mu  # Placeholder for SE results

        # Standard error of simulated weighted mean for NUMERIC fusion variables
        # https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
        se[fcont] <- apply(y[ind, fcont], 2, weightedSE, w = ws)

        # Standard error of weighted mean for BINARY fusion variables
        # https://stats.stackexchange.com/questions/159204/how-to-calculate-the-standard-error-of-a-proportion-using-weighted-data
        se[fbin] <- sqrt(mu[fbin] * (1 - mu[fbin]) * sws2)

        c(mu, se)

      })

    j <- 1:length(fvars)
    m1 <- matrixStats::rowMeans2(temp, rows = j)
    pse <- pooledSE(est = temp[j, ], se = temp[-j, ])
    se1 <- pse[, "se"]
    df1 <- pse[, "degf"]

    #----------

    # Correlation analysis
    # NOTE: weights::wtd.cors() will return all NA's if the supplied weights are too small in magnitude

    # Observed correlation using only the primary sample weights
    r0 <- weights::wtd.cors(x = xs[, fvars], y = xs, weight = ws)
    #r0[!is.finite(r0)] <- 0   # Replace zero-variance cases with zero correlation?
    r0 <- as.vector(r0)[keep]

    # Not clear if I should calculate average over replicates (slower) or just use point estimate above
    # ALT: Or compute average correlation across replicate weights
    # temp <- 1:ncol(wrs) %>%
    #   sapply(function(i) {
    #     weights::wtd.cors(x = xs[, fvars], y = xs, weight = wrs[, i])
    #   })
    # temp[!is.finite(temp)] <- 0   # Replace zero-variance cases with zero correlation?
    # r0 <- matrixStats::rowMeans2(temp, rows = keep)
    #se0 <- matrixStats::rowSds(temp, rows = keep)

    #---

    # Simulated correlations (mean across implicates)
    temp <- 1:ncol(Mind) %>%
      sapply(function(i) {
        ind <- Yind[, i]
        ytemp <- y[ind, fvars]
        weights::wtd.cors(x = ytemp, y = cbind(ytemp, xso), weight = ws)
      })
    #temp[!is.finite(temp)] <- 0   # Replace zero-variance cases with zero correlation?
    r1 <- matrixStats::rowMeans2(temp, rows = keep)
    #se1 <- matrixStats::rowSds(temp, rows = keep)

    #---

    out <- list(share = sum(ws) / sum(w),
                mu = cbind(m0, se0, df0, m1, se1, df1),
                cor = cbind(r0, r1)
    )

  }, mc.cores = cores)

  #-----

  # Save 'comp' to disk
  # saveRDS(comp, "CEI_validate_comp.rds", compress = TRUE)
  # saveRDS(template, "CEI_validate_template.rds", compress = TRUE)

  #-----

  # Confidence interval overlap
  # See Compare.CI() here: https://github.com/cran/synthpop/blob/master/R/compare.syn.r
  # CIoverlap(10, 20, 5, 25)
  # CIoverlap <- function(lwr_obs, upr_obs, lwr_sim, upr_sim) {
  #   L <- pmax(lwr_obs, lwr_sim)
  #   U <- pmin(upr_obs, upr_sim)
  #   0.5 * (((U - L) / (upr_obs - lwr_obs)) + ((U - L) / (upr_sim - lwr_sim)))
  # }

  #-----

  # Two-sample t-test
  # https://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
  # https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha
  Ttest <- function(m1, m2, se1, se2, df1, df2) {
    se <- sqrt(se1 ^ 2 + se2 ^ 2)
    df <- ((se1 ^ 2 + se2 ^ 2) ^ 2) / ((se1 ^ 2) ^ 2 / df1 + (se2 ^ 2) ^ 2 / df2)
    tstat <- (m1 - m2) / se
    pval <- 2 * pt(-abs(tstat), df)
    return(pval)
  }

  # # Confirm accuracy
  # x1 <- rnorm(100)
  # x2 <- rnorm(200)
  # t.test(x1, x2)
  # Ttest(m1 = mean(x1),
  #        m2 = mean(x2),
  #        se1 = sd(x1) / sqrt(length(x1)),
  #        se2 = sd(x2) / sqrt(length(x2)),
  #        df1 = length(x1) - 1,
  #        df2 = length(x2) - 1)


  #-----

  # Compile output from mean analysis
  est.out <- do.call(rbind, purrr::map(comp, "mu")) %>%
    as.data.table(keep.rownames = "v1") %>%
    mutate(
      iset = rep(1:nsets, each = nrow(.) / nsets),
      share = rep(purrr::map_dbl(comp, "share"), each = nrow(.) / nsets),

      y = purrr::map_chr(strsplit(v1, "..", fixed = TRUE), 1),
      level = sapply(v1, function(x) sub("^.*\\.\\.", "", x)),
      level = ifelse(y == level, NA, level))

  # P-value of test for difference in means
  if (!is.null(wr)) est.out <- mutate(est.out, pvalue = round(Ttest(m0, m1, se0, se1, df0, df1), 4))

  # deprecated...
  # ftype = ifelse(y %in% fcont, "cont", "cat"),
  # ftype = ifelse(y %in% fbin, "bin", ftype),

  # fmax = rep(matrixStats::colMaxs(y[, -1]), times = nrow(.) / length(fvars)),
  # fmin = rep(matrixStats::colMins(y[, -1]), times = nrow(.) / length(fvars)),

  # Calculate the confidence interval
  # lci0 = m0 + qt(0.025, df = df0) * se0,
  # uci0 = m0 - qt(0.025, df = df0) * se0,
  # lci1 = m1 + qt(0.025, df = df1) * se1,
  # uci1 = m1 - qt(0.025, df = df1) * se1,
  #
  # # Enforce zero-bound on SE, if plausible
  # lci0 = ifelse(fmin == 0, pmax(0, lci0)),
  # lci1 = ifelse(fmin == 0, pmax(0, lci1)),
  # uci0 = ifelse(fmax == 0, pmin(0, uci0)),
  # uci1 = ifelse(fmax == 0, pmin(0, uci1)),

  # Absolute percent error in point estimates
  # ape = abs((m1 - m0) / m0),
  # ape = ifelse(!is.finite(ape), NA, ape),

  # MOVE TO PLOTTING FUNCTION....
  # Point estimate discrepancy; percent for continuous variables and absolute otherwise
  # error_pe = ifelse(is.na(level), abs((m1 - m0) / m0), abs(m1 - m0)),
  # error_pe = ifelse(!is.finite(error_pe), NA, error_pe),
  #
  # # Standard error discrepancy; percent relative to the observed mean
  # error_se = abs((se1 - se0) / m0),  # Denominator is observed mean
  # error_se = ifelse(!is.finite(error_se), NA, error_se),

  est.out <- est.out %>%
    arrange(iset, y) %>%
    select(iset, share, y, level, m0, m1, any_of(c("se0", "se1", "pvalue"))) %>%
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.numeric, cleanNumeric, tol = 0.001)

  #-----

  # Compile output from correlation analysis

  cor.out <- do.call(rbind, purrr::map(comp, "cor")) %>%
    as.data.table() %>%
    mutate(y = rep(template$v1b, nsets),
           correlate = rep(template$v2b, nsets),
           share = rep(purrr::map_dbl(comp, "share"), each = nrow(template)),
           iset = rep(1:nsets, each = nrow(template)),
           r0 = ifelse(!is.finite(r0), NA, r0),
           r1 = ifelse(!is.finite(r1), NA, r1)) %>%
    group_by(iset, share, y, correlate) %>%
    summarize(error = mean(abs(r1 - r0), na.rm = TRUE),
              error_w = weighted.mean(abs(r1 - r0), w = r0 ^ 2, na.rm = TRUE),
              .groups = "drop") %>%
    mutate_all(~ ifelse(is.na(.x), NA, .x)) %>%
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.numeric, cleanNumeric, tol = 0.001)

  #-----

  # Report processing time
  tout <- difftime(Sys.time(), t0)
  cat("Total processing time:", signif(as.numeric(tout), 3), attr(tout, "units"), "\n", sep = " ")

  # Prepare final result
  #rm(comp)
  out <- list(estimates = est.out,
              correlation = cor.out)
  return(out)

}

