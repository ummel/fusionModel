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

#-------------------------------

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

#-------------------------------

# Detect dependencies among categorical variables
# Returns a list of data frames that can be used to add/merge the derivative variables when the parent is known

detectCategoricalDependence <- function(data, fvars, exclude, cores = 1L) {

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

    # Rename output data frame using "_BINARY_" suffix for numeric variables
    for (v in names(out)) {
      if (v %in% vnum) {
        names(out[[v]])[1] <- paste0(names(out[[v]])[1], "_BINARY_")
      }
    }

  }

  return(out)

}
