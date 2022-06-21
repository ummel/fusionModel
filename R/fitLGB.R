fitLGB <- function(data.lgb, hyper.grid, params.obj, cv.folds) {

  cv.perf <- hyper.grid %>%
    lapply(FUN = function(x) {
      sink <- capture.output({
        mod.cv <- lightgbm::lgb.cv(
          params = c(as.list(x), params.obj),
          data = data.lgb,
          folds = cv.folds,
          early_stopping_rounds = 1L,
          verbose = -1L
        )
      })
      c(mod.cv$best_score, mod.cv$best_iter)
    })

  # Compare cross-validated performance across hyper-parameter sets
  comp <- do.call(rbind, cv.perf)
  opt <- which.min(comp[, 1])
  params.opt <- hyper.grid[[opt]]
  params.opt$num_iterations <- comp[opt, 2]

  mod <- lightgbm::lgb.train(
    params = c(params.opt, params.obj),
    data = data.lgb,
    verbose = -1L
  )

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

