#' Compare two analyses of synthetic data
#'
#' @description
#' Calculates confidence interval overlap coefficients for comparing two results from \link{analyze()}.
#'
#' @param analysis1 Results tibble from \link{analyze()}.
#' @param analysis2 Results tibble from \link{analyze()}.
#'
#' @details Using the 95% confidence intervals for each statistic common to both \code{analysis1} and \code{analysis2}, the overlap coefficient is calculated as proposed by Karr et al. (2006). Overlap equals one if the two confidence intervals are identical, zero if the intervals barely touch, and is increasingly negative as two confidence intervals move further apart.
#'
#' @return A tibble reporting the \code{overlap} coefficient and associated confidence intervals of each statistic.
#'
#' @references  A. F Karr, C. N Kohnen, A Oganian, J. P Reiter & A. P Sanil. (2006). A Framework for Evaluating the Utility of Data Altered to Protect Confidentiality, The American Statistician, 60:3, 224-232.
#'
#' @examples
#' # See usage in Examples of ?analyze().
#' @export

compare <- function(analysis1, analysis2) {
  join.by <- names(analysis1)[1:(which(names(analysis1) == "estimate") - 1)]
  stopifnot(all(join.by %in% names(analysis2)))
  temp <- full_join(analysis1, analysis2, by = join.by, suffix = c(".1", ".2"))
  stopifnot(nrow(temp) < nrow(analysis1) + nrow(analysis2))
  comp <- temp %>%
    mutate(overlap = CIoverlap(lwr_obs = lower_ci.1, upr_obs = upper_ci.1,  lwr_sim = lower_ci.2, upr_sim = upper_ci.2)) %>%
    select(any_of(join.by), overlap, starts_with("lower_ci"), starts_with("upper_ci")) %>%
    arrange_at(join.by)
  return(comp)
}

# Confidence interval overlap
# See Compare.CI() here: https://github.com/cran/synthpop/blob/master/R/compare.syn.r
CIoverlap <- function(lwr_obs, upr_obs, lwr_sim, upr_sim) {
  overlap.lower <- pmax(lwr_obs, lwr_sim)
  overlap.upper <- pmin(upr_obs, upr_sim)
  0.5 * (((overlap.upper - overlap.lower) / (upr_obs - lwr_obs)) + ((overlap.upper - overlap.lower) / (upr_sim - lwr_sim)))
}

# Check
#CIoverlap(1, 2, 1, 2)
#CIoverlap(1, 2, 2, 3)
#CIoverlap(1, 2, 3, 4)
