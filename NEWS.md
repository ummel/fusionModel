# fusionModel 0.2.0

* `train()` is now fully parallel when `mc = TRUE` on UNIX-like systems.

* `lasso` argument added to `train()` enabling fast pre-screening of predictors. See `?train`.

* Fixed bug preventing disabling of `maxcats` argument in `train()`. Now disabled by default (`maxcats = NULL`).

* `induce` argument added to `fuse()` to (optionally) turn off induced rank correlation step.

* [`pbapply`](https://cran.r-project.org/web/packages/pbapply/index.html) package now required for multicore and progress bars.

* [`glmnet`](https://cran.r-project.org/web/packages/glmnet/index.html) package now suggested for LASSO pre-screening.

# fusionModel 0.1.0

* Initial beta.
