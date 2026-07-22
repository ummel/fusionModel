# Analyze Fused Microdata Implicates (Legacy)

Calculates point estimates and associated margins of error (MOE) for
analyses performed on fused or synthetic microdata. Supports calculation
of means, proportions, sums, counts, and medians, with optional
breakdown across population subgroups.

**Legacy Notice:** This function documents and implements the point
estimate and variance estimation methodology used in early microdata
fusion publications (e.g., 2024). A newer, improved statistical approach
has since been introduced in the fusionACS package ([fusionACS
Methods](https://ummel.github.io/fusionACS/articles/methods.html)).
While `analyze()` is retained in `fusionModel` for legacy support and
replication purposes, users are strongly encouraged to review the
updated fusionACS workflow for new analysis pipelines.

## Usage

``` r
analyze(
  x,
  implicates,
  static = NULL,
  weight = NULL,
  rep_weights = NULL,
  by = NULL,
  fun = NULL,
  var_scale = 4,
  cores = 1
)
```

## Arguments

- x:

  List. A named list specifying the desired analysis type(s) and
  associated target variable(s). Supported analysis types include
  `"mean"`, `"sum"`, and `"median"`. Example:
  `x = list(mean = c("v1", "v2"), median = "v3")`. Target variables that
  are factors automatically return proportions (for `"mean"`) or counts
  (for `"sum"`). Target variables must exist in `implicates`, `static`,
  or be created by a custom `fun`.

- implicates:

  Data frame (or `data.table`). Synthetic/fused microdata containing
  implicates, typically produced by
  [`fuse`](https://ummel.github.io/fusionModel/reference/fuse.md).
  Implicates must be row-stacked and identified by an integer column
  named `"M"`.

- static:

  Data frame (or `data.table`), optional. Static (non-synthetic)
  variables that do not vary across implicates. Must satisfy
  `nrow(static) == nrow(implicates) / max(implicates$M)` and match the
  row ordering of `implicates`.

- weight:

  Character, optional. Name of the primary sampling weight column in
  `static`. If `NULL` (default), uniform weights equal to 1 are assumed.

- rep_weights:

  Character vector, optional. Vector of replicate weight column names in
  `static`. If provided, standard errors reflect additional sampling
  weight uncertainty across replicates.

- by:

  Character vector, optional. Column name(s) present in `implicates` or
  `static` defining population subgroups for stratified estimation. If
  `NULL`, analysis is performed over the full sample.

- fun:

  Function, optional. A custom transformation function applied to input
  data prior to analysis. Must return a `data.frame` containing custom
  derived variables.

- var_scale:

  Numeric. Scaling factor applied to unadjusted replicate weight
  variance, determined by the survey design. Default is `4` (appropriate
  for ACS and RECS).

- cores:

  Integer. Number of CPU cores for parallel processing. (Applicable to
  Unix-based systems; defaults to `1`).

## Value

A [`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html)
containing summary results grouped by any `by` variables. Columns
include:

- `by...`:

  Subgroup groupings, if `by` was specified.

- `N`:

  Number of observations per implicate in the analyzed subgroup.

- `y`:

  Name of the analyzed target variable.

- `level`:

  Factor level label (if target variable `y` was a factor/logical).

- `type`:

  Metric type: `"mean"`, `"proportion"`, `"sum"`, `"count"`, or
  `"median"`.

- `est`:

  Pooled point estimate across implicates.

- `moe`:

  Margin of error corresponding to a 90% confidence interval.

## Details

Inputs are checked for consistent row dimensions and implicate
structures. Estimates and standard errors are computed independently for
each implicate. The final point estimate represents the simple mean of
estimates across implicates (\\M\\). Standard errors and degrees of
freedom are pooled across implicates using Rubin's (1987) rules.

When replicate weights (`rep_weights`) are provided, standard errors for
each implicate account for sampling design variance. Within-implicate
variance is evaluated around the point estimate (equivalent to setting
`mse = TRUE` in `survey::svrepdesign`).

When replicate weights are absent, within-implicate variance for means
relies on Cochran's (1977) ratio variance approximation (Gatz & Smith,
1995). Proportions use a weighted variance formula, while medians use an
asymptotic density approximation (for large \\N\\) or bootstrap sampling
(for small \\N\\).

## References

Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.). John Wiley &
Sons.

Gatz, D. F., & Smith, L. (1995). The Standard Error of a Weighted Mean
Concentration — I. Bootstrapping vs Other Methods. *Atmospheric
Environment*, 29(11), 1185–1193.

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons.

## Examples

``` r
if (FALSE) { # \dontrun{
library(fusionModel)

# Build a fusion model using RECS microdata
fusion.vars <- c("electricity", "natural_gas", "aircon")
predictor.vars <- names(recs)[2:12]
fsn.path <- train(data = recs, y = fusion.vars, x = predictor.vars)

# Generate 30 implicates of the 'fusion.vars' using original RECS as the recipient
sim <- fuse(data = recs, fsn = fsn.path, M = 30)
head(sim)

# Full-sample analysis across multiple targets and metrics
result <- analyze(
  x = list(
    mean = c("natural_gas", "aircon"),
    median = "electricity",
    sum = c("electricity", "aircon")
  ),
  implicates = sim,
  weight = "weight"
)
head(result)

# Mean electricity consumption by climate zone and urban/rural status
result1 <- analyze(
  x = list(mean = "electricity"),
  implicates = sim,
  static = recs,
  weight = "weight",
  by = c("climate", "urban_rural")
)

# Subgroup analysis incorporating sample weight uncertainty via replicate weights
result2 <- analyze(
  x = list(mean = "electricity"),
  implicates = sim,
  static = recs,
  weight = "weight",
  rep_weights = paste0("rep_", 1:96),
  by = c("climate", "urban_rural")
)

# Custom derivation function prior to analysis
my_fun <- function(data) {
  kwh_per_ft2 <- data$electricity / data$square_feet
  use_natural_gas <- data$natural_gas > 0
  data.frame(kwh_per_ft2, use_natural_gas)
}

result_custom <- analyze(
  x = list(mean = c("kwh_per_ft2", "use_natural_gas", "electricity")),
  implicates = sim,
  static = recs,
  weight = "weight",
  fun = my_fun
)
} # }
```
