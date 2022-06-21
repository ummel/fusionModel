# # Example data
# data <- recs[1:26]
# recipient <- subset(data, select = c(division, urban_rural, climate, income, age, race, hh_size, renter))
# x <- names(recipient)
# y <- setdiff(names(data), c(x, "weight"))
#
# cl <- cluster(data[y], w = data$weight)
#
# # TEST BIG DATA
# #data <- data[sample.int(nrow(data), size = 100e3, replace = TRUE),]
#
# test <- chain(data, cl, x, w = "weight")

#-----

chain <- function(data, y, x, w = NULL) {

  stopifnot(exprs = {
    is.data.frame(data)
    !anyNA(data)
    all(unlist(y) %in% names(data))
    all(x %in% names(data))
    is.null(w) | w %in% names(data)
  })

  if (is.list(y)) {
    input <- y
    y <- unlist(input)
  } else {
    input <- as.list(y)
  }

  w <- if (is.null(w)) {
    rep(1L, nrow(data))
  } else {
    data[[w]] / mean(data[[w]])
  }

  d <- data[c(x, y)]
  d <- mutate_if(d, is.numeric, rank)
  d <- mutate_if(d, is.ordered, as.integer)
  d <- mutate_if(d, is.logical, as.integer)
  lev <- lapply(d, levels)

  # Get one-hot expanded variable names
  V <- lapply(names(lev), function(v) if (is.null(lev[[v]])) v else paste(v, lev[[v]], sep = "_"))
  names(V) <- names(lev)

  # Observed class proportions/probabilities for the 'y' variables (if nominal)
  yweight <- lapply(y, function(v) {
    yv <- V[[v]]
    w <- if (length(yv) == 1) {
      1
    } else {
      xt <- xtabs(paste0("~", v), data = d)
      xt / sum(xt)
    }
  })
  names(yweight) <- y
  rm(data)

  # One-hot encode unordered factors
  d <- mltools::one_hot(as.data.table(d))
  stopifnot(all(unlist(V) == colnames(d)))

  # Create sparse matrix for input to glmnet()
  d <- as(as.matrix(d), "dgCMatrix")

  # Fit Lasso
  lasso <- function(y, x) {
    m <- glmnet::glmnet(
      x = d[, x],
      y = d[, y],
      family = "gaussian",
      weights = w,
      alpha = 1)
    i <- which(m$dev.ratio / max(m$dev.ratio) >= 0.95)[1] # Index of preferred lambda
    # plot(m$dev.ratio, type = "l")
    # abline(v = i, h = m$dev.ratio[i], lty = 2)
    # m$lambda[i]
    m$dev.ratio[i] # Akin to R-squared
  }

  #-----

  # Stable 'x' predictor variables
  X <- unlist(V[x])

  # Group of y-variables
  #g <- input[[9]]

  # Placeholder for results sequence
  N <- length(input)
  ord <- vector(mode = "list", length = N)
  ratio <- vector(mode = "numeric", length = N)

  # Extract R2 for the "full" model that includes all possible regressors
  Y <- unlist(V[y])  # All y's are in play for full model
  mfull <- lapply(input, function(g) {
    sapply(g, function(v) {
      yv <- unlist(V[v])
      m <- sapply(yv, lasso, x = c(X, setdiff(Y, unlist(V[g]))))
      w <- yweight[[v]]
      sum(m * w)
    })
  })

  # Selection sequence
  for (i in 1:(N - 1)) {
    Y <- unlist(V[unlist(ord)])
    m0 <- lapply(input, function(g) {
      sapply(g, function(v) {
        yv <- unlist(V[v])
        m <- sapply(yv, lasso, x = c(X, Y))
        w <- yweight[[v]]
        sum(m * w)
      })
    })

    # How strong is m0 relative to mfull?
    stopifnot(all(lengths(m0) == lengths(mfull)))
    rel <- sapply(m0, mean) / sapply(mfull, mean)

    # Select the input for highest 'rel'
    b <- which.max(rel)  # best
    ratio[i] <- rel[b]
    ord[[i]] <- input[[b]]
    input <- input[-b]
    mfull <- mfull[-b]
  }

  # Assign final entry
  ratio[N] <- 1
  ord[[N]] <- input[[1]]

  # Return results
  #return(list(order = ord, relperf = ratio))
  #print(signif(setNames(ratio, seq_along(ratio)), 4))
  if (all(lengths(ord) == 1)) ord <- unlist(ord)
  return(ord)

}
