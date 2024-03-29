
fitLGB <- function(dfull, dtrain = NULL, dvalid = NULL, cv.folds = NULL, hyper.grid, params.obj) {

  perf <- if (is.null(dvalid)) {

    # If full cross-validation is requested...
    lapply(hyper.grid, FUN = function(x) {
      sink <- utils::capture.output({
        mod <- lightgbm::lgb.cv(
          params = c(as.list(x), params.obj),
          data = dfull,
          folds = cv.folds,
          early_stopping_rounds = 1L,
          verbose = -1L
        )
      })
      c(mod$best_score, mod$best_iter)
    })

  } else {

    # If training/test-set validation is requested...
    lapply(hyper.grid, FUN = function(x) {
      p <- c(as.list(x), params.obj)
      sink <- utils::capture.output({
        mod <- lightgbm::lgb.train(
          params = p,
          data = dtrain,
          valids = list(valid = dvalid),
          early_stopping_rounds = 1L,
          verbose = -1L
        )
      })
      c(mod$best_score, mod$best_iter)
    })

  }

  # Compare validation set performance across hyper-parameter sets
  comp <- do.call(rbind, perf)
  opt <- which.min(comp[, 1])
  params.opt <- hyper.grid[[opt]]
  params.opt$num_iterations <- as.integer(comp[opt, 2])

  # Fit final model using full dataset and optimal parameter values
  mod <- lightgbm::lgb.train(
    params = c(params.opt, params.obj),
    data = dfull,
    verbose = -1L
  )

  # Add the optimal validation score and number of iterations to the 'mod' object
  # This are NA and -1, by default, which doesn't provide any useful information
  mod$best_score <- as.numeric(comp[opt, 1])
  mod$best_iter <- params.opt$num_iterations

  # Plot the evolution of the loss function
  #plot(unlist(mod$record_evals[[2]]))

  return(mod)

}

#---
#!!! NOTE: It appears callbacks are not (yet) exported but that could be coming soon.
# See here: https://github.com/microsoft/LightGBM/pull/5018
# https://github.com/Microsoft/LightGBM/blob/master/R-package/R/callback.R
# https://stackoverflow.com/questions/54027734/adding-callbacks-to-lightgbm-in-r
# IMPORTANT -- this is presumably the way to set the 'min_delta' argument: https://github.com/ummel/fusionModel/issues/24
# IN PROGRESS: https://github.com/microsoft/LightGBM/pull/5018#pullrequestreview-892006092
# Add here: https://github.com/microsoft/LightGBM/pull/5123  ( appears to be merged?!?!?)
# mod.cv <- lightgbm::lgb.cv(
#   params = c(params.global, obj, min_data_in_leaf = min.leaf),
#   data = dfull,
#   nrounds = 500L,  # See note at top about min_delta. Could set higher once min_delta is allowed...
#   folds = cv.folds,
#   early_stopping_rounds = 1L,
#   verbose = -1L,
#   callbacks = list(cb_early_stop(verbose = FALSE))  # !!! EXAMPLE -- may suppress console output
# )
#---

