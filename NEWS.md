# fusionModel 0.3.0

* Training and fusion process now utilizes either `rpart()` or `biglm()` models, the latter selected only when highly-predictive (R-squared > 0.99). Relevant changes made to internal function `fitRpart()`.

* `train()` includes detection of pairwise "derivative" variables; i.e. cases where one variable is fully dependent upon another. This has two flavors: either a linear model when Pearson correlation is > 0.99^2 or explicit merge for categorical derivatives. `fuse()` handles simulation accordingly. See `detectDependence.R`.

* Related to above, `train()` attempts to detect "structural zeros" in donor data and `fuse()` preserves that structure in the fusion output.

* `train()` has new arguments to control speed of initial model-fitting and to explicitly pass `cp` and `minbucket` arguments to `rpart()`. Ellipsis argument removed. See `?train`.

* `train()` allows user to set number of cores via new `cores` argument. Defaults to one. Old argument `mc` removed.

* `use.biglm` argument removed from `fuse()`; it is now default when `induce = TRUE`.

* `recs` example data expanded to include more variables.

* README updated.

* Various bug fixes.

# fusionModel 0.2.0

* `train()` is now fully parallel when `mc = TRUE` on UNIX-like systems.

* `lasso` argument added to `train()` enabling fast pre-screening of predictors. See `?train`.

* Fixed bug preventing disabling of `maxcats` argument in `train()`. Now disabled by default (`maxcats = NULL`).

* `induce` argument added to `fuse()` to (optionally) turn off induced rank correlation step.

* [`pbapply`](https://cran.r-project.org/web/packages/pbapply/index.html) package now required for multicore and progress bars.

* [`glmnet`](https://cran.r-project.org/web/packages/glmnet/index.html) package now suggested for LASSO pre-screening.

# fusionModel 0.1.0

* Initial beta.
