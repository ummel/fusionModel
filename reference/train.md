# Train a Data Fusion Model

Trains a statistical fusion model on a "donor" dataset using sequential
[LightGBM](https://lightgbm.readthedocs.io/en/latest/) gradient boosting
models to capture conditional distributions. The resulting fitted model
archive (`.fsn` file) contains the conditional expectations, candidate
donor pool indices, and metadata required by
[`fuse`](https://ummel.github.io/fusionModel/reference/fuse.md) to
simulate synthetic outcomes onto a "recipient" dataset.

## Usage

``` r
train(
  data,
  y,
  x,
  fsn = "fusion_model.fsn",
  weight = NULL,
  nfolds = 5,
  nquantiles = 2,
  nclusters = 2000,
  krange = c(10, 500),
  hyper = NULL,
  fork = FALSE,
  cores = 1
)
```

## Arguments

- data:

  A data frame (or `data.table`) containing the donor microdata.
  Categorical variables should be formatted as factors (ordered whenever
  applicable).

- y:

  Character vector or list. Variable(s) in `data` to be fused to a
  recipient dataset. Variables are modeled and fused sequentially in the
  order provided. If `y` is a list, each element can be a character
  vector representing a "block" of variables to be sampled jointly
  during fusion to preserve multivariate dependence.

- x:

  Character vector or list. Predictor variable name(s) common to both
  donor and recipient datasets. If a list, each element specifies the
  predictor set to use for the corresponding element in `y`. If a
  character vector, earlier `y` variables in the sequence are
  automatically appended as predictors for subsequent `y` models.

- fsn:

  Character string. File path where the trained fusion model archive
  will be saved. Must end with the `.fsn` extension. Default is
  `"fusion_model.fsn"`.

- weight:

  Character string. Name of the column in `data` containing survey or
  sampling weights. If `NULL` (default), uniform observation weights are
  assumed.

- nfolds:

  Numeric. Number of cross-validation folds used during LightGBM model
  tuning. If `0 < nfolds < 1`, it represents the proportion of
  observations allocated to training, with the remainder used for
  validation (substantially faster than full cross-validation). Default
  is `5`.

- nquantiles:

  Numeric. Number of quantile models to fit for continuous `y` variables
  in addition to the conditional mean model. Specified quantiles are
  evenly spaced. For example, `nquantiles = 2` (default) models the 25th
  and 75th percentiles. Even values are recommended, as the conditional
  mean already captures central tendency.

- nclusters:

  Numeric. Maximum number of \\k\\-means clusters used to group donor
  observations in conditional expectation space. Higher values increase
  donor selection precision at the expense of memory and processing
  time. Set to `0` or `Inf` to skip clustering (i.e., treat every donor
  row as a cluster center). Default is `2000`.

- krange:

  Numeric vector of length 2. Specifies the minimum and maximum number
  of nearest neighbors (\\k\\) evaluated when selecting optimal
  candidate pools for continuous conditional distributions. Default is
  `c(10, 500)`.

- hyper:

  List. LightGBM hyperparameter grid or custom values to evaluate during
  training. If `NULL` (default), a standardized baseline configuration
  optimized for survey fusion is used. See Details.

- fork:

  Logical. If `TRUE`, uses parallel processing via process forking
  ([`mclapply`](https://rdrr.io/r/parallel/mclapply.html)) across fusion
  steps. Only supported on Unix/Linux/macOS platforms. Default is
  `FALSE`.

- cores:

  Integer. Number of physical CPU cores allocated to computation. When
  `fork = FALSE` or on Windows, fusion steps are processed serially,
  while LightGBM utilizes `cores` for internal OpenMP multithreading.
  When `fork = TRUE` on Unix, independent fusion steps are executed
  concurrently across `cores`. Default is `1`.

## Value

Returns the file path to the saved `.fsn` archive invisibly.

## Details

### Sequential and Block Modeling

Data fusion proceeds sequentially through the variables specified in
`y`. To preserve complex dependencies among tightly coupled outcomes
(e.g., fuel expenditure shares), supply those variables as a character
vector within a list element of `y`. Block variables are predicted
jointly, and donor observations within blocks are sampled en masse
during fusion.

### Hyperparameter Optimization

If a list of vectors is supplied to `hyper`, `train()` performs grid
search across all parameter combinations using \\V\\-fold
cross-validation
([`lgb.cv`](https://rdrr.io/pkg/lightgbm/man/lgb.cv.html)) with early
stopping. The optimal parameter combination (minimizing loss) is
selected to fit the final booster.

When `hyper = NULL`, default hyperparameters applied include:

- `boosting = "gbdt"`

- `data_sample_strategy = "goss"`

- `num_leaves = 31`

- `feature_fraction = 0.8`

- `max_depth = 5`

- `min_data_in_leaf = max(10, round(0.001 * nrow(data)))`

- `num_iterations = 2500`

- `learning_rate = 0.1`

- `max_bin = 255`

- `min_data_in_bin = 3`

- `max_cat_threshold = 32`

### Parallel Execution & OpenMP Considerations

On Unix-like operating systems, process forking via `fork = TRUE` can
yield faster execution times than OpenMP multithreading. However, if
OpenMP threads have already been initialized in the current R session
(e.g., by
[`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html) or
[`fst`](http://www.fstpackage.org/reference/fst.md)), forking may hang.
If this occurs, execute `data.table::setDTthreads(1)` and
`fst::threads_fst(1)` immediately after launching R before calling
`train()`.

## References

Ummel, K., et al. (2024). Multidimensional well-being of US households
at a fine spatial scale using fused household surveys. *Scientific
Data*, 11(142).
[doi:10.1038/s41597-023-02788-7](https://doi.org/10.1038/s41597-023-02788-7)

## See also

[`fuse`](https://ummel.github.io/fusionModel/reference/fuse.md),
[`validate`](https://ummel.github.io/fusionModel/reference/validate.md),
[`fitLGB`](https://ummel.github.io/fusionModel/reference/fitLGB.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load sample RECS survey dataset supplied with fusionModel
data(recs)
# Define fusion targets and common predictors
fusion.vars <- c("electricity", "natural_gas", "aircon")
predictor.vars <- names(recs)[2:12]
# 1. Basic model training (saves output to working directory)
fsn.path <- train(data = recs, y = fusion.vars, x = predictor.vars)
# 2. Block fusion: Preserving joint dependency across component shares
fusion.vars.block <- list(
  "electricity",
  "natural_gas",
  c("heating_share", "cooling_share", "other_share")
)
train(data = recs, y = fusion.vars.block, x = predictor.vars, fsn = "model_block.fsn")
# 3. Custom predictor specification per fusion step
xlist <- list(predictor.vars[1:4], predictor.vars[2:8], predictor.vars)
train(data = recs, y = fusion.vars.block, x = xlist, fsn = "model_xlist.fsn")
# 4. Hyperparameter override (Random Forest boosting alternative)
train(data = recs, y = fusion.vars, x = predictor.vars,
      hyper = list(boosting = "rf", feature_fraction = 0.6, max_depth = 10))
# 5. Grid search hyperparameter tuning
train(data = recs, y = fusion.vars, x = predictor.vars,
      hyper = list(max_depth = c(5, 10), feature_fraction = c(0.7, 0.9)))
} # }
```
