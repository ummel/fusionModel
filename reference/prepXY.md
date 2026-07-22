# Prepare Optimal Fusion Order and Screen Predictor Variables

Determines a data-driven, sequential fusion order for recipient target
variables (`y`) and screens out uninformative predictor variables (`x`)
prior to model fitting. Designed primarily for large-scale donor
datasets containing many and/or collinear variables, `prepXY()` pairs an
initial rank-correlation screening step with LASSO regularization via
[`glmnet`](https://rdrr.io/pkg/glmnet/man/glmnet.html). The generated
output can be passed directly to
[`train`](https://ummel.github.io/fusionModel/reference/train.md).

## Usage

``` r
prepXY(
  data,
  y,
  x,
  weight = NULL,
  cor_thresh = 0.05,
  lasso_thresh = 0.95,
  xmax = 100,
  xforce = NULL,
  fraction = 1,
  cores = 1
)
```

## Arguments

- data:

  Data frame (or `data.table`) containing donor training data.
  Categorical variables should be represented as factors (ordered
  whenever applicable) to ensure optimal dummy encoding and LASSO
  evaluation.

- y:

  Character vector or list. Variable names in `data` intended for fusion
  to a recipient dataset. If passed as a list, individual elements can
  contain character vectors of multiple variables to force their fusion
  together as a unified block.

- x:

  Character vector. Predictor variable names present in `data` that are
  common to both the donor and recipient datasets.

- weight:

  Character string, optional. Name of the observation sampling weight
  column in `data`. If `NULL` (default), uniform weights equal to 1 are
  assumed.

- cor_thresh:

  Numeric value between 0 and 1. Predictors exhibiting an absolute
  Spearman rank correlation below `cor_thresh` relative to a target `y`
  variable are filtered out prior to LASSO optimization. Defaults to
  `0.05`.

- lasso_thresh:

  Numeric value between 0 and 1. Controls predictor screening
  aggressiveness during LASSO regularization. Lower values screen more
  aggressively. For example, `lasso_thresh = 0.95` (default) retains the
  subset of candidate predictors that collectively account for at least
  95% of the deviance explained by a full LASSO model.

- xmax:

  Integer. Soft ceiling on the maximum number of predictors returned by
  the LASSO step. Serves as a performance guardrail when candidate `x`
  pools are very large. Set to `Inf` to disable upper-bound constraints.
  Defaults to `100`.

- xforce:

  Character vector, optional. Subset of `x` predictor variable names to
  unconditionally retain across all target variables, bypassing
  correlation and LASSO screens.

- fraction:

  Numeric value strictly greater than 0 and less than or equal to 1.
  Fraction of observations in `data` to randomly sample during
  screening. Sampling significantly decreases computation times on large
  microdata files with minimal impact on variable selection. Defaults to
  `1` (full dataset).

- cores:

  Integer. Number of physical CPU cores used for parallel execution via
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html). Applicable on
  Unix-like operating systems (Linux/macOS). Defaults to `1`.

## Value

A named list containing two primary slots:

- y:

  A list of character vectors indicating the recommended, sequential
  order for fusing target variables.

- x:

  A list of character vectors of equal length to `y`, specifying the
  preferred subset of `x` predictors associated with each target
  variable step.

Additional diagnostic attributes are attached to the output list:

- xpredictors:

  Character vector of all unique common predictors retained across any
  of the target steps.

- xforce:

  The character vector of forced predictor variables provided by the
  user.

- xoriginal:

  The original vector of candidate `x` predictor variables passed into
  the function.

## Details

`prepXY()` establishes a disciplined, empirical sequence for microdata
fusion while reducing dimensionality before full model training in
[`train`](https://ummel.github.io/fusionModel/reference/train.md).

**Methodological Overview:**

- 1\. Zero-Inflation & Factor Handling:

  Zero-inflated numeric target variables are automatically split into a
  binary indicator (`*_zero`) and a non-zero sub-model to handle
  spike-at-zero distributions. High-cardinality factors are lumped to
  manage dummy expansion.

- 2\. Spearman Rank Correlation Screen:

  A fast rank-based correlation matrix is calculated between all `y`
  target levels and candidate `x` predictors. Variables falling below
  `cor_thresh` are screened out early.

- 3\. Full Model Baseline Fitting:

  LASSO models (`alpha = 1`) are fitted for all candidate target
  variables against the remaining predictor pool to establish maximum
  achievable deviance explained (\\R\_{max}^2\\).

- 4\. Iterative Chain Construction:

  Target variables are greedily ordered by identifying which variable
  achieves the highest fraction of its total potential deviance
  explained using only common `x` predictors and previously selected `y`
  targets in the chain. Predictors meeting the `lasso_thresh` deviance
  ratio are retained for that step.

The resulting list matches the structural expectation of
[`train`](https://ummel.github.io/fusionModel/reference/train.md),
enabling direct down-stream pipeline integration.

## Examples

``` r
if (FALSE) { # \dontrun{
library(fusionModel)
data(recs)

# Select candidate target (y) and predictor (x) variables
y <- names(recs)[c(14:16, 20:22)]
x <- names(recs)[2:13]

# Group first two y variables into a joint fusion block
y_blocked <- c(list(y[1:2]), y[-c(1:2)])

# Run prepXY to determine preferred fusion ordering and predictor screening
prep <- prepXY(data = recs, y = y_blocked, x = x, cor_thresh = 0.05, lasso_thresh = 0.95)

# Pass the prepared lists directly to train()
trained_model <- train(data = recs, y = prep$y, x = prep$x)
} # }
```
