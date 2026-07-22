# Enforce Monotonic Relationships Between Numerical Variables

Smoothes and transforms a numeric target vector (`y`) relative to an
ordering vector (`x`) to enforce a strict, monotonic relationship
(either consistently increasing or decreasing) across all observed
values of `x`. Designed primarily for post-fusion adjustment of energy
consumption (`x`) and expenditure (`y`) variables to guarantee a
plausible pricing structure (e.g., higher fuel volume never yields lower
total cost). By default, the weighted mean of `y` is preserved in the
returned vector (`preserve = TRUE`).

## Usage

``` r
monotonic(
  x,
  y,
  w = NULL,
  preserve = TRUE,
  expend = TRUE,
  fast = TRUE,
  nmax = 5000,
  plot = FALSE
)
```

## Arguments

- x:

  Numeric vector. The independent ordering variable (e.g., energy
  consumption in BTU). Must not contain missing values (`NA`).

- y:

  Numeric vector of the same length as `x`. The dependent response
  variable to be transformed (e.g., fuel expenditure in dollars). Must
  not contain missing values (`NA`).

- w:

  Numeric vector of the same length as `x`, optional. Observation
  sampling weights. If `NULL` (default), uniform weights equal to 1 are
  assumed.

- preserve:

  Logical. Should the original weighted mean of `y` be preserved in the
  returned numeric vector? Defaults to `TRUE`.

- expend:

  Logical. Treat `y` as an expenditure variable tied to consumption `x`?
  If `TRUE` (default), safety checks enforce physical consistency:
  negative values in `x` or `y` throw an error, non-zero expenditures
  when `x == 0` are zeroed, and zero expenditures when `x > 0` are
  raised to the minimum observed positive expenditure.

- fast:

  Logical. If `TRUE` (default), rapid smoothing via Friedman's
  super-smoother ([`supsmu`](https://rdrr.io/r/stats/supsmu.html)) is
  performed and directly coerced to monotonicity via sorting. If
  `FALSE`, a Shape Constrained Additive Model
  ([`scam`](https://rdrr.io/pkg/scam/man/scam.html)) is attempted if the
  initial super-smoother coercion results in excessive error.

- nmax:

  Integer. Maximum number of observations sampled for model fitting to
  optimize computational speed. Defaults to `5000`. Set to `Inf` to
  disable sampling.

- plot:

  Logical. Should a diagnostic scatterplot of the sampled input points
  and the fitted monotonic relationship (in red) be rendered to the
  active graphics device? Defaults to `FALSE`.

## Value

A numeric vector of modified, monotonic `y` values matching the length
and order of input `x`. If input `y` was integer-typed, the returned
vector is rounded to integer.

## Details

`monotonic()` provides non-parametric post-processing to rectify logical
inconsistencies in microdata fusion output (such as negative marginal
prices or non-monotonic tariff curves).

**Algorithmic Workflow:**

- 1\. Physical Sanity Checks:

  If `expend = TRUE`, boundary conditions are enforced to ensure zero
  consumption yields zero expenditure and positive consumption yields
  positive expenditure.

- 2\. High-Speed Subsampling:

  If `length(x) > nmax`, extreme boundaries (minimum and maximum) are
  preserved while intermediate points are randomly down-sampled to
  `nmax` for efficient smoothing.

- 3\. Super-Smoother Fit & Direction Detection:

  An initial non-parametric smooth is fit via
  [`supsmu`](https://rdrr.io/r/stats/supsmu.html). The overall
  correlation between `x` and predicted `y` determines whether the
  relationship should be monotonic increasing or decreasing.

- 4\. Monotonic Coercion & SCAM Fallback:

  Predicted values are sorted to enforce monotonicity. If `fast = FALSE`
  and sorting introduces substantial error (\>5% relative divergence on
  over 5% of points), a constrained spline is fit via
  [`scam`](https://rdrr.io/pkg/scam/man/scam.html). If SCAM fitting
  fails to converge, linear regression
  ([`lm`](https://rdrr.io/r/stats/lm.html)) serves as the ultimate
  fallback.

- 5\. Interpolation & Mean Preservation:

  Fitted values are mapped back to the complete original `x` vector
  using linear interpolation
  ([`approx`](https://rdrr.io/r/stats/approxfun.html)). If
  `preserve = TRUE`, the resulting vector is scaled so its weighted mean
  matches `weighted.mean(y, w)`.

## Examples

``` r
if (FALSE) { # \dontrun{
library(fusionModel)
data(recs)

# Enforce a monotonic pricing curve between propane consumption and expenditure
adjusted_expend <- monotonic(
  x = recs$propane_btu,
  y = recs$propane_expend,
  w = recs$weight,
  plot = TRUE
)

# Compare original and adjusted weighted means
weighted.mean(recs$propane_expend, recs$weight)
weighted.mean(adjusted_expend, recs$weight)
} # }
```
