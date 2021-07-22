fitRpart <- function(y, x, w, data, maxcats = NULL, lasso.threshold = NULL, args) {

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
  # })

  # Is 'y' continuous?
  ycon <- is.numeric(data[[y]])
  unordered <- !ycon & !is.ordered(data[[y]])

  # Restrict data to non-NA response (y) values
  data <- data %>%
    select(all_of(c(w, y, x))) %>%
    filter(!is.na(data[[y]]))

  #-----

  # Collapse categorical predictor categories, if necessary
  # Only necessary if the response variable (y) is categorical
  kmeans.xwalk <- NULL
  if (unordered  & !is.null(maxcats)) {

    # Number of categorical levels for each allowable, non-numeric predictor
    cats.count <- sapply(data[x], function(x) length(levels(x)))

    # Categorical predictor variables that need to be collapsed prior to fitting rpart() model
    xlarge <- names(cats.count[cats.count > maxcats])

    if (length(xlarge) > 0) {
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

  # If requested, screen 'x' predictors via LASSO regression
  lasso.ignore <- if (is.null(lasso.threshold)) {
    NULL
  } else {
    LASSOignore(y = y,
                x = x,
                w = w,
                data = data,
                lasso.threshold = lasso.threshold)
  }

  #-----

  # Formula object
  # NOTE that 'x' predictors are potentially excluded via 'lasso.ignore'
  fobj <- as.formula(paste0(y, "~", paste(setdiff(x, lasso.ignore), collapse = "+")))

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

  return(m)

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

LASSOignore <- function(y, x, w, data, lasso.threshold) {

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

  # Fit LASSO model
  if (type == "skip") {
    xdrop <- NULL
  } else {
    m <- suppressWarnings(
      glmnet::glmnet(
        x = data[x],
        y = Y,
        family = type,
        weights = data[[w]],
        alpha = 1,
        pmax = length(x) - 1,  # Ensures there is at least one non-zero variable along with the intercept
        type.multinomial = "grouped"  # Only relevant when 'type' = "multinomial"
      )
    )

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
    xdrop <- rownames(mcoef)[-1L][as.vector(mcoef == 0)[-1L]]

  }

  return(xdrop)

}
