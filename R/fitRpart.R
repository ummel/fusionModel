fitRpart <- function(y, x, w, data, n = NULL, maxcats = NULL, linear = TRUE, lasso.threshold = 1, args) {

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
                      method = ifelse(ycon, "anova", "class"))

    # Adds any custom arguments specified in 'args'
    args.list <- c(args.list, args)

    # Call rpart() with specified arguments
    m <- do.call(rpart::rpart, args = args.list)

    # If cross-validation used, select the pruned tree that minimized cross-validation error
    if ("xerror" %in% colnames(m$cptable)) {
      m <- rpart::prune(m, cp = m$cptable[which.min(m$cptable[, "xerror"]), "CP"])
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

#--------------------
#--------------------

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

#--------------------
#--------------------

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
