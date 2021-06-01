# TEST CASES
# w <- NULL
# N <- 1000
#
# x <- as.integer(c(rnorm(10, -5, 2), rep(0, 10), rnorm(20, 20, 50)))
# inner.range <- c(0.8, 0.8) * c(max(x[x < 0]), min(x[x > 0]))
# outer.range <- c(1.2, 1.2) * range(x)
#
# x <- as.integer(c(rnorm(20, 20, 5), rep(0, 10)))
# inner.range <- c(0, 10)
# outer.range <- c(0, 1.2 * max(x))

#-----

fitDensity <- function(x,
                       w = NULL,
                       N = 500,
                       inner.range = c(0, 0),
                       outer.range = c(-Inf, Inf)) {

  # Check the arguments
  # TURNED OFF for use in train(), since arguments are guaranteed to be valid (slightly faster)
  # stopifnot(exprs = {
  #   is.numeric(x)
  #   !anyNA(x)
  #   is.null(w) | length(w) == length(x)
  #   N > 0
  #   N %% 1 == 0
  #   is.numeric(inner.range)
  #   length(inner.range) == 2
  #   is.numeric(outer.range)
  #   length(outer.range) == 2
  # })

  # Set default for 'w', if necessary
  if (is.null(w)) {
    w <- rep(1, length(x)) / length(x)
  } else {
    w <- w / sum(w)
  }

  # Copy of original 'x' vector
  x0 <- x

  #-----

  if (novary(x)) {

    #P <- rep(0.5, length(x))
    Q <- rep(x[1], N)

  } else {

    # Fit density
    # NOTE: bw = "nrd" throws error if there are too few 'x' values, but bw = "nrd" does not
    d <- density(x, weights = w, bw = "nrd0")

    x <- d$x
    y <- d$y
    dx <- diff(x)
    y <- y[-1L] - diff(y) / 2
    y <- y / sum(y * dx)  # This adjusts the densities so the cumulative probability sums to 1
    pcum <- c(0, cumsum(y * dx))  # Cumulative probability at each value of 'x'

    # Quantile values associated with each of 'N' percentiles given by 'p'
    # See plot(diff(p)) to confirm that the tails are over-sampled by design
    # This is because the tails are typically more spread out than the center of the distribution (i.e. we want more resolution in the tails)
    p <- dnorm(seq(-3, 3, length.out = N - 1))
    p <- c(0, cumsum(p / sum(p)))
    Q <- suppressWarnings(spline(pcum, x, xout = p, method = "hyman")$y)

    #-----

    # Restrict the min and max of 'Q'
    Q[Q < outer.range[1]] <- outer.range[1]
    Q[Q > outer.range[2]] <- outer.range[2]

    # Set the correct number of zeros in 'Q' (if any)
    if (any(x0 == 0)) {
      pneg <- weighted.mean(x0 < 0, w)
      ppos <- weighted.mean(x0 > 0, w)
      ind <- which.min(abs(p - pneg)):which.min(abs(p - (1 - ppos)))
      Q[ind] <- 0L
    }

    # Restrict the "inner range" of 'Q'
    Q[Q > inner.range[1] & Q < 0] <- inner.range[1]
    Q[Q > 0 & Q < inner.range[2]] <- inner.range[2]

    # Reduce precision of Q
    Q <- if (is.integer(x0)) {
      as.integer(round(Q))
    } else {
      cleanNumeric(Q, tol = 0.001)
    }

    #plot(p, Q)

    #----

    # Conditional percentile for each of the original 'x' values
    #P <- suppressWarnings(approx(Q, p, xout = x0)$y)

    #----

    # TESTING: Draw random sample from the inverse CDF
    # ptile <- runif(1e3)
    # z <- round(approx(x = p, y = Q, xout = ptile)$y)
    # table(z) / length(z)
    # weighted.mean(x0 == 0, w)
    # mean(z == 0)
    # summary(x0)
    # summary(z)
    # plot(p, Q, type = "l")
    # lines(sort(ptile), sort(z), col = 2)


  }

  # Return P and Q along with the (weighted) proportion of zeros in 'x'
  # return(list(P = P,
  #             Q = Q))

  return(Q)

}
