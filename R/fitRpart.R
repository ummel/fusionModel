fitRpart <- function(y, x, w, data, maxcats = 10, args = NULL) {

  stopifnot(exprs = {
    length(y) == 1
    y %in% names(data)
    all(x %in% names(data))
    !y %in% x
    length(w) == 1
    w %in% names(data)
    is.numeric(maxcats)
    maxcats > 1 & maxcats %% 1 == 0
  })

  # Is 'y' continuous?
  ycon <- is.numeric(data[[y]])

  # Restrict data to non-NA response (y) values
  data <- data %>%
    select(all_of(c(w, y, x))) %>%
    filter(!is.na(data[[y]]))

  #-----

  # Collapse categorical predictor categories, if necessary
  # Only necessary if the response variable (y) is categorical
  kmeans.xwalk <- NULL
  if (!ycon) {

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

  # Formual object
  fobj <- as.formula(paste0(y, "~", paste(x, collapse = "+")))

  # Fit rpart() model
  args.list <- c(list(formula = fobj,
                      data = data,
                      weights = data[[w]],
                      method = ifelse(ycon, "anova", "class")),
                 args)

  m <- do.call(rpart::rpart, args = args.list)

  # m <- rpart::rpart(formula = fobj,
  #                   data = data,
  #                   weights = data[[w]],
  #                   method = ifelse(ycon, "anova", "class"),
  #                   ...)
  # minbucket = max(50, ceiling(0.0001 * nrow(data))),
  # cp = 0,  # Set sufficiently low without imposing computing cost
  # xval = 0)

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
  if (!ycon) m$kmeans.xwalk <- kmeans.xwalk

  return(m)

}

#---------

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
    k <- kmeans(x = d[-1L], centers = n)
  } else {
    k <- group_by_at(d, -1L) %>% mutate(cluster = cur_group_id())
  }

  # Return crosswalk between original 'x' and the cluster assignment
  d[paste0(x, "__clus")] <- factor(paste("cluster", k$cluster, sep = "_"))
  d <- d[c(1, ncol(d))]
  return(d)

}
