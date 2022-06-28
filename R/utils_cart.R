# Functions used internally by trainCART() and/or fuseCART()

#-----------------------------------------
#-----------------------------------------

# TEST CASES
# w <- NULL
# N <- 1000
#
# x <- as.integer(c(rnorm(10, -5, 2), rep(0, 10), rnorm(20, 20, 50)))
# inner.range <- c(0.8, 0.8) * c(max(x[x < 0]), min(x[x > 0]))
# outer.range <- c(1.2, 1.2) * range(x)
#
# x <- as.integer(c(rnorm(20, 20, 5), rep(0, 10)))
# inner.range <- c(0, 10)
# outer.range <- c(0, 1.2 * max(x))

fitDensity <- function(x,
                       w = NULL,
                       N = 500,
                       inner.range = c(0, 0),
                       outer.range = c(-Inf, Inf)) {

  # Check the arguments
  stopifnot(exprs = {
    is.numeric(x)
    !anyNA(x)
    is.null(w) | length(w) == length(x)
    N > 0 & N %% 1 == 0
    is.numeric(inner.range)
    length(inner.range) == 2
    is.numeric(outer.range)
    length(outer.range) == 2
  })

  # Set default for 'w', if necessary
  if (is.null(w)) {
    w <- rep(1, length(x)) / length(x)
  } else {
    w <- w / sum(w)
  }

  # Copy of original 'x' vector
  x0 <- x

  #-----

  if (novary(x)) {

    #P <- rep(0.5, length(x))
    Q <- rep(x[1], N)

  } else {

    # Fit density
    # NOTE: bw = "nrd" throws error if there are too few 'x' values, but bw = "nrd" does not
    d <- density(x, weights = w, bw = "nrd0")

    x <- d$x
    y <- d$y
    dx <- diff(x)
    y <- y[-1L] - diff(y) / 2
    y <- y / sum(y * dx)  # This adjusts the densities so the cumulative probability sums to 1
    pcum <- c(0, cumsum(y * dx))  # Cumulative probability at each value of 'x'

    #plot(pcum, x)

    # Quantile values associated with each of 'N' percentiles given by 'p'
    # See plot(diff(p)) to confirm that the tails are over-sampled by design
    # This is because the tails are typically more spread out than the center of the distribution (i.e. we want more resolution in the tails)
    p <- dnorm(seq(-3, 3, length.out = N - 1))
    p <- c(0, cumsum(p / sum(p)))
    Q <- suppressWarnings(spline(pcum, x, xout = p, method = "hyman")$y)

    #-----

    # Restrict the min and max of 'Q' (outer range)
    Q[Q < outer.range[1]] <- outer.range[1]
    Q[Q > outer.range[2]] <- outer.range[2]

    # Set the correct number of zeros in 'Q' (if any)
    if (any(x0 == 0)) {
      pneg <- weighted.mean(x0 < 0, w)
      ppos <- weighted.mean(x0 > 0, w)
      ind <- which.min(abs(p - pneg)):which.min(abs(p - (1 - ppos)))
      Q[ind] <- 0L
    }

    # Restrict the "inner range" of 'Q'
    Q[Q > inner.range[1] & Q < 0] <- inner.range[1]
    Q[Q > 0 & Q < inner.range[2]] <- inner.range[2]

    # Reduce precision of Q
    Q <- if (is.integer(x0)) {
      as.integer(round(Q))
    } else {
      cleanNumeric(Q, tol = 0.001)
    }

    #plot(p, Q)

    #----

    # Conditional percentile for each of the original 'x' values
    #P <- suppressWarnings(approx(Q, p, xout = x0)$y)

    #----

    # TESTING: Draw random sample from the inverse CDF
    # ptile <- runif(1e3)
    # z <- approx(x = p, y = Q, xout = ptile)$y
    # table(z) / length(z)
    # weighted.mean(x0 == 0, w)
    # mean(z == 0)
    # summary(x0)
    # summary(z)
    # plot(p, Q, type = "l")
    # lines(sort(ptile), sort(z), col = 2)


  }

  # Return P and Q along with the (weighted) proportion of zeros in 'x'
  # return(list(P = P,
  #             Q = Q))

  return(Q)

}

#-----------------------------------------
#-----------------------------------------

fitRpart <- function(y, x, w, data, n = NULL, maxcats = NULL, linear = TRUE, lasso.threshold = 1, cvfactor = 0, args) {

  # Turned off since not necessary when used within train()
  # stopifnot(exprs = {
  #   length(y) == 1
  #   y %in% names(data)
  #   all(x %in% names(data))
  #   !y %in% x
  #   length(w) == 1
  #   w %in% names(data)
  #   is.numeric(maxcats)
  #   maxcats > 1 & maxcats %% 1 == 0
  #   is.null(n) | (n > length(x) & n <= nrow(data))
  # })

  # Is 'y' continuous?
  ycon <- is.numeric(data[[y]])
  unordered <- !ycon & !is.ordered(data[[y]])

  # Restrict data to non-NA response (y) values
  data <- data %>%
    select(all_of(c(w, y, x))) %>%
    filter(!is.na(data[[y]]))

  # Downsample 'data', if requested
  if (is.null(n)) n <- nrow(data)
  if (n < nrow(data)) {
    data <- data %>%
      slice_sample(n = n, replace = TRUE)
  }

  #-----

  # Output from the LASSO step can be:
  # NULL or character vector: 'x' variables to ignore in subsequent rpart() call
  # biglm() model that predicts response variable using only linear relationships

  lasso.out <- if (is.null(lasso.threshold)) {
    NULL
  } else {
    fitLASSO(y = y,
             x = x,
             w = w,
             data = data,
             lasso.threshold = lasso.threshold,
             linear = linear)
  }

  #-----

  if (class(lasso.out) == "biglm") {

    return(lasso.out)

  } else {

    # Collapse categorical predictor categories, if necessary
    # Only necessary if the response variable (y) is categorical
    kmeans.xwalk <- NULL
    if (unordered  & !is.null(maxcats)) {

      # Number of categorical levels for each allowable, non-numeric predictor
      cats.count <- sapply(data[x], function(x) length(levels(x)))

      # Categorical predictor variables that need to be collapsed prior to fitting rpart() model
      xlarge <- names(cats.count[cats.count > maxcats])

      if (length(xlarge) > 0) {

        for (i in xlarge) collapseCategorical(x = i, y = y, w = w, data = data, n = maxcats)

        kmeans.xwalk <- lapply(xlarge, collapseCategorical, y = y, w = w, data = data, n = maxcats)
        names(kmeans.xwalk) <- xlarge
        for (d in kmeans.xwalk) {
          v <- names(d)[1]
          data <- left_join(data, d, by = v)
          x[x == v] <- names(d)[2]
        }
        data[xlarge] <- NULL
      } else {
        kmeans.xwalk <- NULL
      }

    }

    #-----

    # Formula object
    # NOTE that 'x' predictors are potentially excluded via 'lasso.out'
    xpred <- setdiff(x, lasso.out)
    fobj <- as.formula(paste0(y, "~", paste(xpred, collapse = "+")))

    # NOT USED: Experimenting with forcing aggressive splitting of certain predictors via 'cost' argument
    # cost <- rep(1, length(xpred))
    # cost <- replace(cost, force %in% xpred, 1e-14)  # Arbitrary low cost value to force early inclusion of 'force' variables

    # Arguments list for rpart, passed via do.call()
    args.list <- list(formula = fobj,
                      data = data,
                      weights = data[[w]],
                      method = ifelse(ycon, "anova", "class"),
                      minsplit = ceiling(args$minbucket * 2))

    # Adds any custom arguments specified in 'args'
    args.list <- c(args.list, args)

    # Call rpart() with specified arguments
    m <- do.call(rpart::rpart, args = args.list)

    # If cross-validation used, select the pruned tree that is within cvfactor-SE of the minimum cross-validation error
    # Note that this technically forces at least one split in the tree
    # See here: https://stats.stackexchange.com/questions/17904/one-standard-error-rule-for-variable-selection
    if ("xerror" %in% colnames(m$cptable)) {
      imin <- which.min(m$cptable[, "xerror"])
      xse <- m$cptable[imin, "xstd"] / sqrt(args$xval - 1)  # Approximate SE derived from cross-validated SD and number of folds
      ind <- max(which(m$cptable[1:imin, "xerror"] >= m$cptable[imin, "xerror"] + cvfactor * xse))
      m <- rpart::prune(m, cp = m$cptable[max(2, ind), "CP"])
    }

    #-----

    if ("variable.importance" %in% names(m)) {
      m$variable.importance <- m$variable.importance / sum(m$variable.importance)
    } else {
      m$variable.importance <- rep(0L, length(x))
      names(m$variable.importance) <- x
    }

    # Add the 'y' variable names to rpart object
    m$yvar <- y

    # Add the collapsed predictor crosswalk, if present
    if (unordered) m$kmeans.xwalk <- kmeans.xwalk

    m <- slimRpart(m)

    return(m)

  }

}

#-----------------------------------------
#-----------------------------------------

collapseCategorical <- function(x, y, w, data, n) {

  mm <- model.matrix.lm(object = as.formula(paste("~ 0 +", y)), data = data, na.action = "na.pass")

  d <- data[c(w, x)] %>%
    cbind(mm) %>%
    na.omit() %>%
    group_by_at(x) %>%
    summarize_at(colnames(mm), ~ weighted.mean(.x, !!sym(w), na.rm = TRUE))

  # Derive kmeans() clusters
  maxn <- nrow(distinct(d[-1L]))
  if (maxn > n) {
    k <- kmeans(x = d[-1L], centers = n, nstart = 30)
  } else {
    k <- group_by_at(d, -1L) %>% mutate(cluster = cur_group_id())
  }

  # Return crosswalk between original 'x' and the cluster assignment
  d[paste0(x, "__clus")] <- factor(paste("cluster", k$cluster, sep = "_"))
  d <- d[c(1, ncol(d))]
  return(d)

}

#-----------------------------------------
#-----------------------------------------

fitLASSO <- function(y, x, w, data, lasso.threshold, linear) {

  Y <- data[[y]]

  # Prepare 'Y' and 'type' inputs to glmnet() based on the response variable
  if (is.numeric(Y)) {
    type <- "gaussian"
  } else {
    if (length(levels(Y)) == 2 | is.ordered(Y)) {
      type <- "gaussian"
      Y <- as.integer(Y)
    } else {
      type <- "skip"  # Skip unordered factor response
    }
  }

  # Assign the transformed y-value to 'data'
  data[[y]] <- Y

  #----

  out <- NULL

  if (type != "skip") {

    # Fit LASSO glmnet model
    m <- suppressWarnings(
      glmnet::glmnet(
        x = data[x],
        y = data[[y]],
        family = type,
        weights = data[[w]],
        alpha = 1,
        lambda.min.ratio = 1e-5,
        pmax = length(x) - 1  # Ensures there is at least one non-zero variable along with the intercept
      )
    )

    #-----

    # If the LASSO is highly-predictive, attempt to return a slimmed linear model with estimated variable importance
    if (linear & any(m$dev.ratio > 0.99)) {  # Should specify threshold

      # Determine which predictors have non-zero coefficients in the LASSO results
      cf <- glmnet::coef.glmnet(m, s = min(m$lambda))
      lm.x <- setdiff(rownames(cf)[as.vector(cf != 0)], "(Intercept)")

      # Fit linear model using only the non-zero predictors
      fobj <- formula(paste(y, "~", paste(lm.x, collapse = "+")))
      m.biglm <- biglm::biglm(formula = fobj,
                              data = data,
                              weights = ~ ._wgt_.)

      # Check that the R2 of the final linear model still exceeds threshold; if not, proceed with the glmnet LASSO results (below)
      # Variable importance is approximated by absolute value of the t-statistic (see here: https://www.rdocumentation.org/packages/caret/versions/6.0-88/topics/varImp)
      sum.biglm <- try(summary(m.biglm), silent = TRUE)
      if (class(sum.biglm) != "try-error"){
        if (sum.biglm$rsq > 0.99) {
          cf <- sum.biglm$mat
          tval <- replace_na(abs(cf[, "Coef"] / cf[, "SE"]), 0)
          #tval[sum.biglm$mat[, "p"] > 0.05] <- 0
          tval.mean <- setNames(tapply(tval, m.biglm$assign, FUN = max), c("Intercept", lm.x))
          tval.mean <- tval.mean[-1]  # Drop the intercept
          vimp <- tval.mean / sum(tval.mean)
          vimp <- vimp[vimp > 0]
          m.biglm$variable.importance <- vimp
          out <- m.biglm
        }
      }
    }

    #-----

    # If no highly-predictive linear model is possible (or desired) use LASSO results to determine which predictors to ignore when fitting rpart()
    if (is.null(out)) {

      # Otherwise, determine which predictor variables can be dropped prior to fitting decision tree

      # Find preferred lambda
      ind <- which(m$dev.ratio / max(m$dev.ratio) >= lasso.threshold)[1]
      # plot(m$dev.ratio, type = "l")
      # abline(v = ind, h = m$dev.ratio[ind], lty = 2)
      # m$lambda[ind]

      # Get model coefficients for preferred lambda
      mcoef <- glmnet::coef.glmnet(m, s = m$lambda[ind])

      # In categorical response case, sum coefficients across levels
      if (is.list(mcoef)) mcoef <- Reduce("+", mcoef)

      # Predictor variables to ignore (LASSO coefficient = 0)
      out <- rownames(mcoef)[-1L][as.vector(mcoef == 0)[-1L]]

    }

    return(out)

  }

}

#-----------------------------------------
#-----------------------------------------

# Determine 'yord' of models
fusionOrder <- function(varimp) {

  yvars <- names(varimp)

  var.imp <- map2_dfr(.x = varimp, .y = yvars, .f = getVarImp, yvars = yvars)

  # How important are the xvars for each y?
  # When the xvars have high collective importance, y can be moved forward...
  ximp <- var.imp %>%
    filter(predictor == "_xvars_") %>%
    select(response, importance) %>%
    rename(ximp = importance,
           yvar = response)

  yord <- vector(mode = "character", length = length(yvars))
  for (i in 1:(length(yord) - 1)) {

    yimp <- var.imp %>%
      filter(predictor != "_xvars_") %>%
      group_by(predictor) %>%
      summarize(yimp = mean(importance)) %>%
      rename(yvar = predictor)

    ranking <- left_join(yimp, ximp, by = "yvar") %>%
      arrange(-(ximp + yimp))

    vmax <- ranking$yvar[1]
    yord[i] <- vmax
    var.imp <- filter(var.imp, predictor != vmax, response != vmax)

  }

  # Add the final imputation variable
  yord[length(yord)] <- setdiff(yvars, yord)

  return(yord)

}

#-----------------------------------------
#-----------------------------------------

getVarImp <- function(vimp, y, yvars) {
  vimp %>%
    tibble::enframe(name = "predictor", value = "importance") %>%
    mutate(predictor = ifelse(predictor %in% yvars, predictor, "_xvars_"),
           predictor = factor(predictor, levels = c(yvars, "_xvars_"))) %>%
    group_by(predictor) %>%
    summarize(importance = sum(importance), .groups = 'drop') %>%
    complete(predictor, fill = list(importance = 0)) %>%
    mutate(predictor = as.character(predictor),
           response = y) %>%
    filter(predictor != response)
}

#-----------------------------------------
#-----------------------------------------

# Detect bivariate relationships where the correlation is above some (high) threshold
# For these cases, fit a linear model that explains variable A as a function of B
# Variable A does not need to be simulated, but the linear model must be retained to add A at end of simulation in subsequent call to fuse()
# Specified threshold (R-squared value) for functional equivalence of two variables

detectCorrelationDependence <- function(data, fvars, exclude, threshold = 0.99) {

  out <- NULL

  # Coerce factor to integer for correlation calculation
  dtemp <- data %>%
    select(-any_of(exclude)) %>%
    mutate_if(is.factor, as.integer) %>%
    mutate_if(is.logical, as.integer)

  # Detect (near-)perfect bivariate correlations
  # Compares absolute Pearson correlation to 'threshold'
  # Results data frame 'correlated' gives the variables in each detected dyad
  #z <- suppressWarnings(cor(dtemp, use = ifelse(anyNA(dtemp), "pairwise.complete.obs", "everything")))
  cmat <- suppressWarnings(cor(dtemp, method = "pearson"))

  #cmat[lower.tri(cmat, diag = TRUE)] <- NA
  #z <- apply(cmat, MARGIN = 1, FUN = function(x) names(which(abs(x) > sqrt(threshold)))) %>% compact()

  z <- list()
  for (y in fvars) {
    x <- cmat[rownames(cmat) == y, colnames(cmat) != y]
    x <- abs(x) - sqrt(threshold)
    x <- x[x >= 0]
    m <- names(which.max(x))  # Retains only 1 potential match
    z[[y]] <- m
    if (length(m)) cmat[, y] <- 0
  }

  if (length(z) > 0) {

    correlated <- z %>%
      enframe("retained", "derivative") %>%
      unnest(derivative)

    # Function to fit biglm model for row in "correlated"
    mfun <- function(i) {
      y <- correlated$derivative[i]
      x <- correlated$retained[i]
      fobj <- formula(paste(y, "~", x))
      m <- biglm::biglm(formula = fobj,
                        data = dtemp,
                        weights = ~ ._wgt_.)
      #m <- slimBiglm(m)  # Not really necessary
      return(m)
    }

    # List of slimmed 'biglm' models that can be used for simple simulation of variables identified in 'names(sim.lm)'
    out <- 1:nrow(correlated) %>%
      lapply(mfun) %>%
      setNames(correlated$derivative)

  }
  return(out)
}

#-----------------------------------------
#-----------------------------------------

# Workhorse function
# targets: variables for which we want to know if derivatives exist
# candidates: non-numeric variables that should be considered as potential derivatives of the 'targets'

catDependence <- function(data, targets = NULL, candidates = NULL, cores = 1L) {

  if (is.null(targets)) {
    targets <- names(data)
  } else {
    stopifnot(all(targets %in% names(data)))
  }

  # Restrict 'targets' to eliminate any columns with all unique values (e.g. household ID)
  targets <- data %>%
    subset(select = targets) %>%
    select_if(~ length(unique(na.omit(.x))) < length(.x)) %>%
    names()

  # Non-numeric variables we test to see if they are potential derivatives of 'targets'. This
  # excludes numeric variables, which are assumed to never be derivative of
  # another variable.
  if (is.null(candidates)) {
    candidates <- data %>%
      select_if(~ !is.numeric(.x) & length(unique(na.omit(.x))) > 1) %>%
      names()
  } else {
    stopifnot(all(candidates %in% names(data)))
    stopifnot(all(!sapply(subset(data, select = candidates), is.numeric)))
  }

  #-----

  # TO DO: Write as data.table compatible so input data.table object is preserved and does not require re-conversion
  # Convert factors to integer and coerce to data.table for speed
  DT <- data %>%
    select(all_of(unique(c(targets, candidates)))) %>%
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.factor, as.integer) %>%
    data.table()

  #---

  # Function to iterate over each 'x'...
  #x <- targets[3]

  # For variable 'x', detect any 'candidates' that can be deduced from 'x'
  checkX <- function(x) {

    y <- setdiff(candidates, x)

    # Data subset to work with
    dt <- DT %>%
      subset(!is.na(DT[[x]]), select = c(x, y)) %>%
      setkeyv(x)

    # If 'x' is numeric, first identify which 'y' variables are sorted after sorting on 'x'
    if (is.numeric(data[[x]])) {
      temp <- is.na(dt)
      for (v in y) dt[ , (v) := .GRP, by = v]
      dt[temp] <- NA  # Reassign NA values after group ID's inserted
      check <- unlist(!dt[, lapply(.SD, is.unsorted, na.rm = TRUE)])[1, -1]
      y <- names(which(check))
    }

    # Is a given 'y' deducible from 'x'? 'deduce' returns a logical. 'y' is
    # deducible if each 'x' value is associated with a single 'y' value
    if (length(y) > 0) {
      if (length(y) < ncol(dt) - 1) dt <- dt[, c(x, y), with = FALSE]  # Restrict to remaining 'y' candidates, if necessary
      f <- function(x) {length(unique(na.omit(x))) == 1}
      check1 <- as.vector(!dt[ , lapply(.SD, f), .SDcols = y])  # This check returns FALSE if a 'y' variable has only 1 non-NA value (this is considered incomplete information to determine derivative status)
      check2 <- colSums(!dt[, lapply(.SD, f), by = x][, -1]) == 0  # This check returns TRUE if each 'x' value is associated with a single (non-NA) 'y' value
      y <- names(which(check1 & check2))
    }

    out <- NULL
    if (length(y) > 0) out <- setNames(list(y), x)

    return(out)

  }

  #-----

  # Prevent progress bar output to console when cores = 1
  # This is useful when catDependence() is called inside a parallel process, like in train()
  # See Examples for ?pboptions
  if (cores == 1) {
    opb <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(opb))
  }

  result <- pbapply::pblapply(targets, FUN = checkX, cl = cores) %>%
    compact() %>%
    unlist(recursive = FALSE)

  return(result)

}

#-----------------------------------------
#-----------------------------------------

# Detect dependencies among categorical variables
# Returns a list of data frames that can be used to add/merge the derivative variables when the parent is known

detectCategoricalDependence <- function(data, fvars, exclude, suffix = "_BINARY_", cores = 1L) {

  out <- NULL

  # Resulting list indicates bivariate dependence among categorical variables
  # Following example means: Variable 'have_ac' is derivative of 'aircon'; i.e. 'have_ac' is precisely known if 'aircon' is known
  #   $aircon
  #   [1] "have_ac"

  dtemp <- select(data, -any_of(exclude))
  vcat <- names(select(dtemp, !where(is.numeric)))
  vnum <- setdiff(names(dtemp), vcat)
  dtemp <- dtemp %>%
    mutate_if(is.numeric, ~ .x == 0) %>%
    distinct()

  dd <- catDependence(dtemp,
                      targets = NULL, # The 'predictor' variables in this context
                      candidates = intersect(fvars, vcat),
                      cores = cores)

  if (length(dd) > 0) {

    # NECESSARY?
    # Remove dependencies where the "predictor" is itself explained by another predictor (recursive cases)
    #dd <- dd[setdiff(names(dd), unlist(dd))]

    # Check that the 'fvars' are only "explained" by one relationship in 'dd' (no need for multiple)
    explained <- unlist(dd)
    dd <- dd[unique(match(explained, unlist(dd)))]

    # List of data frames that can be used to add/merge derivative variables when variables in 'names(sim.merge)' are known
    out <- seq_along(dd) %>%
      map(~ dtemp %>%
            select(names(dd)[.x], dd[[.x]]) %>%
            distinct()) %>%
      setNames(names(dd))

    # Rename output data frame using (for example) "_BINARY_" suffix for numeric variables
    for (v in names(out)) {
      if (v %in% vnum) {
        names(out[[v]])[1] <- paste0(names(out[[v]])[1], suffix)
      }
    }

  }

  return(out)

}

#-----------------------------------------
#-----------------------------------------

# Remove elements of rpart() model object to reduce size while still allowing prediction
slimRpart <- function(m) {
  m$call <- NULL
  m$numresp <- NULL
  m$parms <- NULL
  m$functions <- NULL
  m$ordered <- NULL
  m$y <- NULL
  #attr(m, "xlevels") <- NULL
  if ("splits" %in% names(m)) m$splits[, c('improve', 'adj')] <- 0
  if ("frame" %in% names(m)) m$frame[, c('n', 'wt', 'dev', 'complexity')] <- 0
  if ("terms" %in% names(m)) attributes(m$terms)['.Environment'] <- NULL
  return(m)
}

#-----------------------------------------
#-----------------------------------------

# Remove elements of biglm() model object to reduce size while still allowing prediction
# slimBiglm <- function(m) {
#   m$call <- NULL
#   m$assign <- NULL
#   m$df.resid <- NULL
#   m$weights <- NULL
#   m$n <- NULL
#   m$names <- NULL
#   return(m)
# }

#-----------------------------------------
#-----------------------------------------

# Validate input variable vector against a given set of names
# 'x' can include regular expressions
validNames <- function(x, nms, exclude = FALSE) {
  rgx <- setdiff(x, nms)
  v <- c(intersect(x, nms), unlist(lapply(rgx, function(x) grep(x, nms, value = TRUE))))
  if (exclude) v <- setdiff(nms, v)
  out <- v[na.omit(match(nms, v))]
  return(out)
}

#-----------------------------------------
#-----------------------------------------

# Return rpart() model prediction node for each observation in 'newdata'
# Based on: https://github.com/cran/rpart.plot/blob/master/R/rpart.predict.R
predictNode <- function(object, newdata) {

  where <-
    if (missing(newdata)) {
      object$where
    } else {
      if(is.null(attr(newdata, "terms"))) {
        Terms <- delete.response(object$terms)
        newdata <- model.frame(Terms, newdata, na.action = na.pass,
                               xlev = attr(object, "xlevels"))
        if(!is.null(cl <- attr(Terms, "dataClasses")))
          .checkMFClasses(cl, newdata, TRUE)
      }
      newdata <- getFromNamespace("rpart.matrix", ns="rpart")(newdata)
      getFromNamespace("pred.rpart", ns="rpart")(object, newdata)
    }

  # nn <- as.numeric(rownames(object$frame)[where])
  #
  # # Index assigning each row in 'pred' to a node number
  # node <- match(nn, as.integer(row.names(object$frame)))

  #return(node)
  return(where)

}
