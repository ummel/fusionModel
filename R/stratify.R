stratify <- function(y, ycont, tfrac, ntiles) {
  stopifnot((tfrac > 0 & tfrac <= 1) | (tfrac > 1 & tfrac %% 1 == 0))
  if (ycont) y <- dplyr::ntile(y, ntiles)
  if (tfrac <= 1) { # Training set indicator (logical)
    out <- vector(mode = "logical", length = length(y))
    for (i in unique(y)) {
      ind <- y == i
      N <- sum(ind)
      out[ind] <- data.table::frank(runif(N), ties.method = "random") <= round(tfrac * N)
    }
  }
  if (tfrac > 1) {  # Cross-validation folds (list of integer row indices)
    out <- vector(mode = "integer", length = length(y))
    for (i in unique(y)) {
      ind <- y == i
      out[ind] <- dplyr::ntile(runif(sum(ind)), tfrac)
    }
    out <- data.table::as.data.table(out)[, list(list(.I)), by = out]
    out <- out$V1
  }
  return(out)
}

#-----

# Examples...
# y <- sample.int(5, 1e3, replace = TRUE)
# test <- stratify(y, ycont = FALSE, tfrac = 0.75)
# table(test)
#
# test <- stratify(y, ycont = FALSE, tfrac = 5)
# lengths(test)
#
# y <- runif(1e3)
# test <- stratify(y, ycont = TRUE, tfrac = 0.75, ntiles = 10)
# table(test)
#
# test <- stratify(y, ycont = TRUE, tfrac = 5, ntiles = 10)
# lengths(test)
