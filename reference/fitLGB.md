# Internal LightGBM Hyperparameter Tuning and Model Fitting

`fitLGB()` is an internal helper function invoked during
[`train`](https://ummel.github.io/fusionModel/reference/train.md) to
perform hyperparameter grid searches and train final gradient boosting
models using the lightgbm package. It evaluates candidate hyperparameter
sets via K-fold cross-validation or a single validation split, selects
the combination minimizing the loss metric, and returns a final model
trained on the complete dataset.

## Usage

``` r
fitLGB(
  dfull,
  dtrain = NULL,
  dvalid = NULL,
  cv.folds = NULL,
  hyper.grid,
  params.obj
)
```

## Arguments

- dfull:

  An [`lgb.Dataset`](https://rdrr.io/pkg/lightgbm/man/lgb.Dataset.html)
  object containing the full dataset (predictor matrix and target
  outcome) used for the final model fit.

- dtrain:

  An optional
  [`lgb.Dataset`](https://rdrr.io/pkg/lightgbm/man/lgb.Dataset.html)
  object containing the training subset when validation-set evaluation
  is performed. Default is `NULL`.

- dvalid:

  An optional
  [`lgb.Dataset`](https://rdrr.io/pkg/lightgbm/man/lgb.Dataset.html)
  object containing the validation subset used for early stopping and
  performance comparison. Default is `NULL`.

- cv.folds:

  A list of integer vectors specifying predefined fold indices for
  K-fold cross-validation (via
  [`lgb.cv`](https://rdrr.io/pkg/lightgbm/man/lgb.cv.html)), or a
  logical flag (`TRUE`) when a single validation split is used.

- hyper.grid:

  A list of lists or data frame/grid where each element represents a
  specific hyperparameter combination (e.g., learning rate, tree depth,
  feature fraction) to evaluate.

- params.obj:

  A named list of global LightGBM parameters, including objective
  setting, metric definitions, and thread counts.

## Value

An object of class `lgb.Booster` containing the final trained model,
augmented with `$best_score`, `$best_iter`, and `$record_evals` metadata
summarizing the tuning process.

## Details

Internal Fitting Routine for LightGBM Models

**Workflow:**

- **Cross-Validation / Validation Evaluation:** Loops through each
  candidate row in `hyper.grid`. Depending on whether `cv.folds` is a
  list of fold indices or a logical flag, it executes
  [`lgb.cv`](https://rdrr.io/pkg/lightgbm/man/lgb.cv.html) or
  [`lgb.train`](https://rdrr.io/pkg/lightgbm/man/lgb.train.html) with
  early stopping enabled (`early_stopping_rounds = 2L`).

- **Optimal Selection:** Compares `best_score` across all candidate
  hyperparameter sets, identifying the set that achieved the lowest
  evaluation loss and its associated optimal number of boosting
  iterations (`best_iter`).

- **Final Model Training:** Refits the model on the full dataset
  (`dfull`) using the optimal hyperparameters and fixed iteration count
  (`num_iterations = best_iter`).

- **Result Augmentation:** Attaches a metadata data.frame containing the
  grid search results to the returned `lgb.Booster` object's
  `$record_evals` slot.

## See also

[`train`](https://ummel.github.io/fusionModel/reference/train.md)
