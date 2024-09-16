fitLGB <- function(dfull, dtrain = NULL, dvalid = NULL, cv.folds = NULL, hyper.grid, params.obj) {

  # If full cross-validation is requested...
  if (is.list(cv.folds)) {
    perf <-  lapply(hyper.grid, FUN = function(x) {
      sink <- utils::capture.output({
        mod <- lightgbm::lgb.cv(
          params = c(as.list(x), params.obj),
          data = dfull,
          folds = cv.folds,
          early_stopping_rounds = 2L,
          verbose = -1L
        )
      })
      data.frame(best_score = mod$best_score, best_iter = mod$best_iter)
    })
  }

  # If training/test-set validation is requested...
  if (is.logical(cv.folds)) {
    perf <- lapply(hyper.grid, FUN = function(x) {
      p <- c(as.list(x), params.obj)
      sink <- utils::capture.output({
        mod <- lightgbm::lgb.train(
          params = p,
          data = dtrain,
          valids = list(valid = dvalid),
          early_stopping_rounds = 2L,
          verbose = -1L
        )
      })
      data.frame(best_score = mod$best_score, best_iter = mod$best_iter)
    })
  }

  # # If no validation is requested; i.e. over-fitting scenario
  # if (is.null(dvalid) & is.null(cv.folds)) {
  #   perf <- lapply(hyper.grid, FUN = function(x) {
  #     p <- c(as.list(x), params.obj)
  #     sink <- utils::capture.output({
  #       mod <- lightgbm::lgb.train(
  #         params = p,
  #         data = dfull,
  #         verbose = -1L
  #       )
  #     })
  #     # Can't get it to return training loss
  #     # Return the training log-loss at the maximum number of iterations
  #     #train.evals <- unlist(mod$record_evals$train)
  #     #c(min(train.evals), which.min(train.evals))
  #     c(NA, NA)
  #   })
  # }

  # Compare validation set performance across hyper-parameter sets
  comp <- bind_rows(perf)
  opt <- which.min(comp$best_score)
  # This is not ideal -- overfitting process should return training loss and work with multiple hypergrid options (instead of just choosing #1)
  if (length(opt) == 0) {
    params.opt <- hyper.grid[[1]]
  } else {
    params.opt <- hyper.grid[[opt]]
    params.opt$num_iterations <- as.integer(comp$best_iter[[opt]])
  }

  # Fit final model using full dataset and optimal parameter values
  mod <- lightgbm::lgb.train(
    params = c(params.opt, params.obj),
    data = dfull,
    verbose = -1L
  )

  # Add the optimal validation score and number of iterations to the 'mod' object
  # These are NA and -1, by default, which doesn't provide any useful information
  mod$best_score <- as.numeric(comp$best_score[opt])
  mod$best_iter <- params.opt$num_iterations

  # Storing hyper results in 'record_evals' slot, since adding a custom slot is not allowed
  hyper.results <- cbind(bind_rows(hyper.grid), comp, final_model = FALSE)
  hyper.results$final_model[opt] <- TRUE
  mod$record_evals <- hyper.results

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

