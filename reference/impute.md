# Impute Missing Data via Microdata Fusion

Imputes missing values (`NA`) in a data frame using sequential microdata
fusion. Wraps iterative calls to
[`train`](https://ummel.github.io/fusionModel/reference/train.md) and
[`fuse`](https://ummel.github.io/fusionModel/reference/fuse.md) under
the hood to provide a fast, non-parametric imputation workflow powered
by LightGBM gradient boosting.

## Usage

``` r
impute(
  data,
  weight = NULL,
  ignore = NULL,
  cores = parallel::detectCores(logical = FALSE) - 1L
)
```

## Arguments

- data:

  Data frame (or `data.table`) containing missing values (`NA`) to be
  imputed.

- weight:

  Character string, optional. Name of the observation sampling weight
  column in `data`. If `NULL` (default), uniform weights equal to 1 are
  assumed. Weight column must not contain missing values.

- ignore:

  Character vector, optional. Names of columns in `data` to ignore.
  Ignored variables are neither imputed nor utilized as predictors in
  the underlying fusion models.

- cores:

  Integer. Number of physical CPU cores used for parallel computation by
  `fst`, `data.table`, and `LightGBM`. Defaults to physical cores minus
  1 (`parallel::detectCores(logical = FALSE) - 1L`).

## Value

A data frame (or `data.table`, matching the input class of `data`) with
all missing values imputed. Original column ordering and variable data
types are preserved.

## Details

`impute()` automates missing data imputation across complex microdata
datasets through the following steps:

- 1\. Sequential Ordering:

  Variables containing missing values are sorted by total missingness
  and imputed sequentially, starting with the variable with the fewest
  `NA` values. Once a variable is imputed, its complete values
  immediately become available as predictors for subsequent target
  variables.

- 2\. Predictor Screening:

  For larger datasets (\>10,000 rows or \>20 columns), a rank-based
  Spearman correlation matrix is computed to pre-screen predictors. To
  maintain efficiency while preserving strong predictor signals,
  candidate predictors for each target are filtered to those meeting an
  absolute correlation threshold (\> 0.025), capped at a maximum of 30
  predictors (and a minimum of 5).

- 3\. Stratified Training Subsampling:

  For variables with large non-missing sample sizes, a stratified sample
  of non-missing training observations (ranging between 5,000 and 50,000
  rows) is selected to speed up model fitting.

- 4\. Model Fitting & Fusion:

  A LightGBM fusion model is trained using 80% validation splitting
  (`nfolds = 0.8`) to determine optimal tree depth and early stopping.
  Missing values for that target variable are then filled via
  [`fuse`](https://ummel.github.io/fusionModel/reference/fuse.md) with a
  single implicate (`M = 1`).

Because `LightGBM` natively accommodates `NA` values within predictor
columns, partially complete predictor variables can be used effectively
during model training.

**Parallel Processing Note:** Because underlying computations in
`LightGBM`, `fst`, and `data.table` rely on OpenMP multithreading, it is
recommended to manage parallelism via the `cores` parameter rather than
wrapping `impute()` inside outer forked parallel loops.

## Examples

``` r
if (FALSE) { # \dontrun{
library(fusionModel)
# Load sample microdata and introduce random missing values (NA)
data(recs)
sample_data <- recs[, 2:7]
set.seed(123)
miss_mask <- replicate(ncol(sample_data), runif(nrow(sample_data)) < runif(1, 0.05, 0.25))
sample_data[miss_mask] <- NA

# Check missing value counts per column
colSums(is.na(sample_data))

# Impute missing values
imputed_data <- impute(sample_data)

# Verify that no NA values remain
anyNA(imputed_data)
head(imputed_data)
} # }
```
