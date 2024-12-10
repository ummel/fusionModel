#' Analyze fusionACS microdata
#'
#' @description
#' For fusionACS internal use only. Calculation of point estimates and associated uncertainty (margin of error) for analyses using ACS and/or fused donor survey variables.
#' Efficiently computes means, medians, sums, proportions, and counts, optionally across population subgroups.
#' The use of native ACS weights or ORNL UrbanPop synthetic population weights is automatically determined given the requested geographic resolution.
#' Requires a local \code{/fusionData} directory in the working directory path with assumed file structure and conventions.
#'
#' @param analyses List. Specifies the desired analyses. Each analysis is a formula. See Details and Examples.
#' @param year Integer. One or more years for which microdata are pooled to compute \code{analyses} (i.e. ACS recipient year). Currently defaults to \code{year = 2015:2019}, if the \code{by} variables indicate a sub-PUMA analysis requiring UrbanPop weights.
#' @param respondent Character. Should the \code{analyses} be computed using \code{"household"}- or \code{"person"}-level microdata?
#' @param by Character. Optional variable(s) that collectively define the set of population subgroups for which each analysis is computed. Can be a mix of geographic (e.g. census tract) and/or socio-demographic microdata variables (e.g. poverty status); the latter may be existing variables on disk or custom variables created on-the-fly via \code{fun()}. If \code{NULL}, analysis is done for the whole (national) sample.
#' @param area Call. Optional unquoted call specifying a geographic area within which to compute the \code{analyses}. Useful for restricting the study area to a manageable size.
#' @param fun Function. Optional function for creating custom microdata variables that cannot be accommodated in \code{analyses}. Must take \code{data} and (optionally) \code{weight} as the only function arguments and must return a \code{data.frame} with number of rows equal to \code{nrow(data)}. See Details and Examples.
#' @param M Integer. The first \code{M} implicates are used. Set \code{M = Inf} to use all available implicates.
#' @param R Integer. The first \code{R} replicate weights are used. Set \code{R = Inf} to use all available replicate weights.
#' @param cores Integer. Number of cores used for multithreading in \code{\link[collapse]{collapse-package}} functions.
#' @param version_up Integer. Use \code{version_up = 1} to access national, single-implicate weights. Use \code{version_up = 2} to access 10-replicate weights for 17 metro areas.
#' @param force_up Logical. If \code{TRUE}, force use of UrbanPop weights even if the requested analysis can be done using native ACS weights.
#'
#' @details Allowable geographic units of analysis specified in \code{by} are currently limited to: region, division, state, cbsa10, puma10, county10, cousubfp10 (county subdivision), zcta10 (zip code), tract10 (census tract), and bg10 (block group).
#'
#' @details The final point estimates are the mean estimates across implicates. The final margin of error is derived from the pooled standard error across implicates, calculated using Rubin's pooling rules (1987). The within-implicate standard error's are calculated using the replicate weights.
#'
#' @details Each entry in the `analyses` list is a \code{\link[stats]{formula}} of the format `Z ~ F(E)`, where `Z` is an optional, user-friendly name for the analysis, `F` is an allowable “outer function”, and `E` is an “inner expression” containing one or more microdata variables. For example:
#'
#' @details `mysum ~ mean(Var1 + Var2)`
#' @details In this case, the outer function is mean(). Allowable outer functions are: mean(), sum(), median(), sd(), and var(). When the inner expression contains more than one variable, it is first evaluated and then `F()` is applied to the result. In this case, an internal variable `X = Var1 + Var2` is generated across all observations, and then `mean(X)` is computed.
#' @details If no inner expression is desired, the `analyses` list can use the following convenient syntax to apply a single outer function to multiple variables:
#' @details `mean = c("Var1", "Var2")`
#' @details The inner expression can also utilize any function that takes variable names as arguments and returns a vector with the same length as the inputs. This is useful for defining complex operations in a separate function (e.g. microsimulation). For example:
#' @details `myfun = function(Var1, Var2) {Var1 + Var2}`
#' @details `mysum ~ mean(myfun(Var1, Var2))`
#' @details The use of sum() or mean() with an inner expression that returns a categorical vector automatically results in category-wise weighted counts and proportions, respectively. For example, the following analysis would fail if evaluated literally, since mean() expects numeric input but the inner expression returns character. But this is interpreted as a request to return weighted proportions for each categorical outcome.
#' @details `myprop ~ mean(ifelse(Var1 > 10 , 'Yes', 'No'))`
#' @details `analyze_fusionACS()` uses "fast" versions of the allowable outer functions, as provided by \code{\link[collapse]{fast-statistical-functions}} in the `collapse` package. These functions are highly optimized for weighted, grouped calculations. In addition, outer functions mean(), sum(), and median() enjoy the use of platform-independent multithreading across columns when `cores > 1`. Analyses with numerical inner expressions are processed using a series of calls to \code{\link[collapse]{collap}} with unique observation weights. Analyses with categorical inner expressions utilize a series of calls to \code{\link[collapse]{fsum}}.
#'
#' @return A tibble reporting analysis results, possibly across subgroups defined in \code{by}. The returned quantities include:
#' @return \describe{
#'  \item{lhs}{Optional analysis name; the "left hand side" of the analysis formula.}
#'  \item{rhs}{The "right hand side" of the analysis formula.}
#'  \item{type}{Type of analysis: sum, mean, median, prop(ortion) or count.}
#'  \item{level}{Factor levels for categorical analyses; NA otherwise.}
#'  \item{N}{Mean number of valid microdata observations across all implicates and replicates; i.e. the sample size used to construct the estimate.}
#'  \item{est}{Point estimate; mean estimate across all implicates and replicates.}
#'  \item{moe}{Margin of error associated with the 90% confidence interval.}
#'  \item{se}{Standard error of the estimate.}
#'  \item{df}{Degrees of freedom used to calculate the margin of error.}
#'  \item{cv}{Coefficient of variation; conventional scale-independent measure of estimate reliability. Calculated as: \code{100 * moe / 1.645 / est}}
#'  \item{rshare}{Share of \code{moe} attributable to replicate weight uncertainty (as opposed to uncertainty across implicates).}
#'  }
#'
#' @references Rubin, D.B. (1987). \emph{Multiple imputation for nonresponse in surveys}. Hoboken, NJ: Wiley.
#'
#' @examples
#' # Analysis using ACS native weights for year 2017, by PUMA, in South Atlantic Census Division
#' # Uses all available implicates and replicate weights
#' test <- analyze_fusionACS(analyses = list(high_burden ~ mean(dollarel / hincp > 0.05)),
#'                           year = 2017,
#'                           by = "puma10",
#'                           area = division == "South Atlantic")
#'
#' # Analysis using UrbanPop 2015-2019 weights, by tract, in Utah (actually Salt Lake City metro given current UrbanPop data)
#' # Uses 5 (of possible 20) fusion implicates for RECS "dollarel" variable
#' # Uses 5 (of possible 10) UrbanPop replicate weights
#' test <- analyze_fusionACS(analyses = list(median_burden ~ median(dollarel / hincp)),
#'                           year = 2015:2019,
#'                           by = "tract10",
#'                           area = state_name == "Utah",
#'                           M = 5,
#'                           R = 5)
#'
#' # User function to create custom variables from microdata
#' # Variables explicitly referenced in my_fun() are automatically loaded into 'data' within analyze_fusionACS()
#' # Variables returned by my_fun() may be used in 'by' or inner expressions of 'analyses'
#' my_fun <- function(data) {
#'   require(tidyverse, quietly = TRUE)
#'   data %>%
#'     mutate(elderly = agep >= 65,
#'            energy_expend = dollarel + dollarfo + dollarlp + dollarng,
#'            energy_burden = energy_expend / hincp,
#'            energy_burden = ifelse(hincp < 5000, NA, energy_burden)) %>%
#'     select(elderly, energy_burden, energy_expend)
#' }
#'
#' # Analysis using UrbanPop 2015-2019 weights, by zip code and elderly head of household, in Atlanta CBSA
#' test <- analyze_fusionACS(analyses = list(energy_burden ~ mean(energy_burden),
#'                                           at_risk ~ mean(energy_burden > 0.075 | acequipm_pub == "No air conditioning")),
#'                           year = 2015:2019,
#'                           by = c("zcta10", "elderly"),
#'                           area = cbsa10 == "12060",
#'                           fun = my_fun,
#'                           M = 5,
#'                           R = 5)
#'
#' @export

#---------------
#
# library(collapse)
# library(tidyverse)
# library(data.table)
# source("R/utils.R")

# library(fusionModel)
# setwd("/home/kevin/Documents/Projects/fusionData")
#
# test <- analyze_fusionACS(analyses = list(elec_expend ~ mean(dollarel)),
#                           respondent = "household",
#                           year = 2019,
#                           by = "state",
#                           area = region == "West",
#                           M = 5,
#                           R = 5,
#                           cores = 3)
#
#
# analyses = list(elec_expend ~ mean(dollarel), myvar ~ mean(scalee), svar ~ mean(scaleb))
# analyses = list(myvar1 ~ mean(dollarel > 500), myvar2 ~ sum(dollarel > 500), m1 ~ mean(dollarel), m2~median(dollarel))
# analyses = list(myvar1 ~ mean(dollarel > 500), myvar2 ~ sum(dollarel > 500), m1 ~ mean(dollarel), m2 ~ median(dollarel), m3 ~ sum(dollarel))
# respondent = "household"
# year = 2019
# by = "tract10"
# area = substitute(state == 49)
# M = 5
# R = 5
# cores = 3
# fun = NULL
# version_up = 2
# force_up = FALSE

#---------------
#
# library(collapse)
# library(tidyverse)
# library(data.table)
# source("R/utils.R")
#
# setwd("/home/kevin/Documents/Projects/fusionData")

# # Function inputs
# year = 2015:2019 # Way to automate this?
#
# #var <- c("dollarel" ,"dollarng", "dollarfo" ,"dollarlp" ,"hincp", "agep", "pov_ratio", "weight")  # 'weight' is specified here, because it is needed in fun()
#
# # Allow user to use custom function to generate the desired analysis variables
# fun <- function(data, weight) {
#   data %>%
#     mutate(recs_util = dollarel + dollarng + dollarfo + dollarlp,
#            recs_burden = ifelse(hincp > 0, (dollarel + dollarng + dollarfo + dollarlp) / hincp, NA),
#            #myvar = pov_ratio < 2 & recs_burden > matrixStats::weightedMedian(recs_burden, w = weight, na.rm = TRUE),
#            white = ref_race5 == "White") %>%
#     select(recs_util, recs_burden, white)
# }
#
# #by = "county10"
# #by = c("ref_race5", "tract10")
# by = c("ref_race5", "state")
# #area = substitute(state_name == "Louisiana")
# #area = substitute(cbsa10 == "12580")
# area = substitute(region == "West")
# cores = 3
# M <- 2
# R <- 2
# respondent = "household"
#
# # TEST: Get variables used in 'analyses'
# analyses <- list(
#   ~ mean(dollarel + dollarng),
#   ~ median(recs_util),
#   ~ mean(recs_burden > 0.075)
#   #~ mean(myvar)
# )
#
# #---
#
# # Tract level?
#
# test <- analyze_fusionACS(
#   analyses = list(
#     ~ mean(dollarel + dollarng),
#     ~ median(recs_util),
#     ~ mean(recs_burden > 0.075)),
#   respondent = "household",
#   year = 2015:2019,
#   by = c("tract10"),
#   #area = cbsa10 == 14460,
#   area = state %in% 24:26,
#   fun = fun,
#   M = Inf,  #
#   R = Inf,
#   cores = 3
# )
#
# # Reliability?
# x <- findInterval(test$cv, c(-Inf, 12, 30, Inf))
# table(x) / length(x)
#
#
# #----
#
# # Using native ACS replicate weights
# test1 <- analyze_fusionACS(
#   analyses = list(
#     ~ mean(dollarel + dollarng),
#     ~ median(recs_util),
#     ~ mean(recs_burden > 0.075)),
#   respondent = "household",
#   year = 2015:2019,
#   by = c("ref_race5", "puma10"),
#   area = state %in% 24:26,
#   fun = fun,
#   M = 1,  # Assuming single implicate; i.e. like calculating ACS MOE using just replicate weights
#   R = 10,  # 10 to match urbanpop?
#   cores = 3
# )
#
# # Using UrbanPop weights
# test2 <- analyze_fusionACS(
#   analyses = list(
#     ~ mean(dollarel + dollarng),
#     ~ median(recs_util),
#     ~ mean(recs_burden > 0.075)),
#   respondent = "household",
#   year = 2015:2019,
#   by = c("ref_race5", "puma10"),
#   area = state %in% 24:26,
#   fun = fun,
#   M = 1,
#   R = 10,  # Only ten available in UrbanPop
#   cores = 3
# )
#
# # Temp: Restrict to PUMA's in urbanpop data
# test2 <- filter(test2, puma10 %in% test1$puma10)
# test1 <- filter(test1, puma10 %in% test2$puma10)
#
# # TEMP: restrict to similar N?
# i <- abs(test1$N - test2$N) / test1$N < 0.1
# test1 <- test1[i, ]
# test2 <- test2[i, ]
#
# # How to results compare?
# plot(test1$N, test2$N)
# abline(0,1)
#
# plot(test1$est, test2$est)
# abline(0,1)
#
# # Relative MOE; There does not appear to be any bias (on average/mean), but the median value using UrbanPop may be a bit lower (~10%)
# summary(test1$cv)
# summary(test2$cv)
# plot(test1$cv, test2$cv)
# abline(0, 1)
#
# # Reliability categories
# cv.breaks <- c(-Inf, 12, 30, Inf)
# table(findInterval(test1$cv, cv.breaks))
# table(findInterval(test2$cv, cv.breaks))

#-------------------------------
#-------------------------------

# TO DO: Option to turn off MOE calculation

analyze_fusionACS <- function(analyses,
                              year,
                              respondent = "household",
                              by = NULL,
                              area = NULL,
                              fun = NULL,
                              M = Inf,
                              R = Inf,
                              cores = 1,
                              version_up = 2,
                              force_up = FALSE) {

  # Check validity of the working directory path
  # Checks if "/fusionData" is part of the path, as this is required
  b <- strsplit(full.path(getwd()), .Platform$file.sep, fixed = TRUE)[[1]]
  i <- which(b == "fusionData")
  if (length(i) == 0) stop("'/fusionData' is not part of the working directory path; this is required.")
  dir <- paste(b[1:i], collapse = .Platform$file.sep)

  t0 <- Sys.time()
  fst::threads_fst(nr_of_threads = cores)
  setDTthreads(threads = cores)

  # Set collapse package option to silently ignore inapplicable 'nthreads' argument
  # See '...' entry in ?collap; applicable to fsd() and fvar()
  options(collapse_unused_arg_action = "none")

  # Check the 'area' expression to determine how to parse it
  # If 'area' is passed a 'call', it is not modified
  # This is useful for special user input like: area = str2lang(paste0("state == '", st.obj, "'"))
  # The more common case is for 'area' to follow usage like in subset()
  # See here: http://adv-r.had.co.nz/Computing-on-the-language.html
  check <- try(is.call(area), silent = TRUE)
  if (inherits(check, "try-error") | check == FALSE) {
    area <- substitute(area)
    if (is.character(area)) area <- str2lang(area)
  }
  area.vars <- all.vars(area)

  # Respondent identifier ("H" or "P")
  rtype <- substring(toupper(respondent), 1, 1)

  # Initial argument check
  stopifnot({
    rlang::is_formula(analyses) | is.list(analyses)
    rtype %in% c("H", "P")
    all(year >= 2005) & all(diff(year) == 1)
    is.null(by) | is.character(by)
    is.null(area) | is.call(area)
    is.null(fun) | is.function(fun)
    M >= 1 & M %% 1 == 0
    R >= 0 & R %% 1 == 0
    cores > 0 & cores %% 1 == 0
    is.logical(force_up)
    version_up %in% c(1, 2)
  })

  #---

  # Prepare and check 'analyses' input
  if (!is.list(analyses)) analyses <- list(analyses)

  # Attempt to convert any non-formula entries in 'analyses' into a plausible formula
  # This applies to legacy analysis formulation of the kind:
  #  analyses <- list(mean = c("natural_gas", "aircon"), median = "electricity")
  # The code below converts these to an equivalent formula and assigns LHS as concatenation of analysis variable name and outer function
  analyses <- lapply(seq_along(analyses), function(i) {
    x <- analyses[[i]]
    if (!rlang::is_formula(x)) {
      f <- names(analyses)[i]  # The requested outer function
      fobj <- paste0(gsub(" ", "_", str_squish(x)), "_", f, "~", f, "(`", x, "`)")
      lapply(fobj, as.formula)
    } else {
      x
    }
  }) %>%
    unlist()

  #---

  # Check that all 'analyses' are formulas
  if (!all(sapply(analyses, rlang::is_formula))) stop("All 'analyses' must be formulas")

  # Remove any duplicate analyses
  dupe <- duplicated(sapply(analyses, rlang::f_rhs))
  if (any(dupe)) {
    cat("Removed duplicated analyses:\n ", paste0(analyses[dupe], collapse = "\n "), "\n", sep = "")
    analyses <- analyses[!dupe]
  }

  # Check for duplicate analysis names (LHS of formula)
  if (anyDuplicated(purrr::compact(sapply(analyses, rlang::f_lhs)))) stop("Detected duplicate LHS analysis names (must be unique)")

  # Get all potential variables required to perform analyses
  ylist <- lapply(analyses, function(x) all.vars(rlang::f_rhs(x)))
  avars <- unique(unlist(ylist, use.names = FALSE))

  # An analysis variable cannot be in 'by'
  err <- intersect(avars, by)
  if (length(err)) stop("Analysis variables cannot also be in the 'by' argument: ", paste(err, collapse = ", "))

  #----

  # Parse 'analyses' to get the outer function, inner expression, and LHS analysis name
  # https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
  alist <- lapply(seq_along(analyses), function(i) {
    lhs <- as.character(rlang::f_lhs(analyses[[i]]))  # LHS analysis name
    rhs <- rlang::f_text(analyses[[i]])
    rhs <- gsub('\"', "'", rhs, fixed = TRUE)  # Remove any escaped quotes
    f <- sub("\\((.+)$", "", rhs)  # Outer function
    iexp <- sub(paste0(f, "("), "", rhs, fixed = TRUE)
    iexp <- substr(iexp, 1, nchar(iexp) - 1) # Inner expression
    list(f, iexp, lhs, rhs)  # Return required elements
  })

  # Inner expression of each analysis
  aexp <- purrr::map_chr(alist, 2)

  # Internal names of the "analysis variables"; "A..1", "A..2", etc.
  anames <- paste0("A..", match(aexp, unique(aexp)))
  names(alist) <- anames

  # Outer function of each analysis
  afun <- purrr::map_chr(alist, 1)

  # Abbreviation of the inner expression with function appended
  # Used below to assign LHS when none is provided
  afun <- gsub("proportion", "prop", afun)
  lhs.abb <- paste(abbreviate(gsub('`', '', purrr::map_chr(alist, 2))), afun, sep = "_")

  # Convert outer functions, if necessary, and check that the requested function is allowed
  afun <- gsub("count", "sum", afun)  # Alternative way of requesting a sum
  afun <- gsub("prop", "mean", afun)
  invalid <- !afun %in% c('sum', 'mean', 'median')  # Valid outer functions
  if (any(invalid)) stop("Outer functions must be sum(), mean(), or median()")

  # Transform outer functions to "fast" versions in 'collapse' package
  afun[] <- paste0("f", afun)

  # 'afun' MUST include the analysis variable ("A..1", etc) as names attribute
  stopifnot(!is.null(names(afun)))

  # NOT USED: Arguments to pass to the outer function in collap()
  # This is only applicable to "prob" argument for "quantile" (NOT USED)
  #aarg <- sapply(alist, function(x) ifelse(length(x[[2]]), x[[2]], NA))

  # The "ANALYSIS" label used to identify each analysis (combination of function and analysis variable, separated by single dot)
  alabel <- paste(afun, names(afun), sep = ".")

  # LHS and RHS of each analysis; assigned to the final output
  alhs <- sapply(alist, function(x) ifelse(length(x[[3]]), x[[3]], NA))
  arhs <- purrr::map_chr(alist, 4)

  # If no LHS provided, assign an abbreviation based on the inner expression
  alhs[is.na(alhs)] <- lhs.abb[is.na(alhs)]

  #-----

  # Extract input variables required by 'fun' user function

  fvars <- if (!is.null(fun)) {
    fargs <- names(formals(fun))
    if (!all(fargs %in% c('data', 'weight'))) stop("Supplied fun() must have 'data' and (optionally) 'weight' as arguments)")
    setdiff(all.vars(body(fun)), c("data", avars))
  } else {
    NULL
  }

  #-----

  # Create 'geocon' data frame containing crosswalk between block group (geoid), state, PUMA, and any geographic 'by' variables
  # This is used to subset/filter the ACS microdata for the requested geographic extent AND (possibly) merge UrbanPop weights at the block group level

  # Report target area
  cat("Geographic area:", if(is.null(area)) "national" else deparse(area), "\n")

  # Crosswalk linking block groups to the target geography ('gtarget'), possibly subsetted by 'area' filter argument
  geocon <- fst::read_fst(path = file.path(dir, "geo-processed/concordance/geo_concordance.fst"),
                          columns = unique(c(area.vars, 'region', 'division', 'state_name', 'state', 'cbsa10', 'puma10', 'county10', 'cousubfp10', 'tract10', 'bg10', 'zcta10'))) %>%
    mutate(keep = eval(area)) %>%  # Adds 'keep' column only if area is non-NULL
    filter(bg10 != "None") %>%  # Not sure why there are some entries with "None" for 'bg10' variable
    distinct()

  # Identify which of the 'by' variable(s) is geographic
  gtarget <- intersect(by, names(geocon))
  stopifnot(length(gtarget) %in% c(0, 1))
  if (length(gtarget) > 0) {

    # Create crosswalk between block groups (geoid) and the target geography (gtarget)
    # What are the allowable geographic variables for purposes of downscaling?
    # Not a comprehensive list: https://www.census.gov/programs-surveys/geography/guidance/geo-identifiers.html
    cat("Geographic unit of analysis:", gtarget, "\n")
    if (is.null(area)) geocon$keep <- TRUE
    geocon <- geocon %>%
      fmutate(geoid = paste0(state, county10, tract10, bg10),
              gtarget = switch(gtarget,
                               region = region,
                               division = division,
                               state_name = state_name,
                               state = state,
                               cbsa10 = cbsa10,
                               puma10 = paste0(state, puma10),
                               county10 = substring(geoid, 1, 5),
                               cousubfp10 = paste0(state, county10, cousubfp10),
                               tract10 = substring(geoid, 1, 11),
                               bg10 = geoid,
                               zcta10 = zcta10)
      ) %>%
      group_by(gtarget) %>%
      mutate(keep = all(keep)) %>%  # 'keep' = TRUE only if ALL of a target geography's block groups are within the requested 'area' (prevents returning a spatial subset of a geography)
      ungroup() %>%
      filter(keep) %>%
      distinct(geoid, state, puma10, gtarget) %>%
      as.data.table()

    if (nrow(geocon) == 0) stop("No geographic intersections identified. Check if your 'by' and 'area' arguments are compatible.")

    # Determine how many PUMA's are associated with each value of the target geography
    # Used to determine if ACS replicate weights can be used instead of UrbanPop
    # If the target geography is perfectly delineated by PUMA's, there is no need to use UrbanPop; can simply use native ACS replicate weights
    acs.check <- geocon %>%
      distinct(puma10, state, gtarget) %>%
      add_count(puma10, state)

    # Use urbanpop weights?
    use.up <- any(acs.check$n != 1) | force_up
    rm(acs.check)

  } else {

    cat("No geographic unit of analysis in 'by'", gtarget, "\n")
    geocon <- geocon %>%
      distinct(state, puma10) %>%
      mutate(gtarget = NA) %>%
      as.data.table()

    use.up <- force_up

  }

  #-----

  if (use.up) {

    # Use UrbanPop dataset to get assignment of household ID (hid) to the target geography (gtarget), possibly for multiple UrbanPop replicate weights

    # Adjust variables in 'geocon' to match variables names in UrbanPop dataset
    geocon <- geocon %>%
      select(-state, -puma10) %>%  # Can remove 'state' and 'puma10' since using UrbanPop
      tidyr::separate_wider_position(geoid, widths = c(state = 2, county10 = 3, tract10 = 6, bg10 = 1)) %>%  # This function is experimental but works for now...
      mutate_at(vars(state, county10, tract10, bg10), as.integer)

    # Path to UrbanPop data
    up.path <- ifelse(version_up == 1, "urbanpop/v1/2015-2019", "urbanpop/v2/2015-2019")

    # If using urbanpop, any 'year' input is overridden
    cat("Using UrbanPop weights for 2015-2019 period (version ", version_up,")\n", sep = "")
    year <- 2015:2019

    # Get paths to required .fst UrbanPop files, identified by state and replicate
    upaths <- lapply(unique(geocon$state), function(s) {
      p <- list.files(path = file.path(up.path, str_pad(s, width = 2, pad = 0)), full.names = TRUE)
      if (length(p)) p[1:ifelse(is.finite(R), min(R + 1, length(p)), length(p))] else NULL
    }) %>%
      unlist()

    if (length(upaths) == 0) {
      stop("Could not identify any UrbanPop .fst files to load within: ", up.path)
    } else {
      cat("Loading UrbanPop data from disk\n")
    }

    # Load UrbanPop data from disk
    # Note that read_fsd() uses gtemp to subset the data as it is loaded (reduces memory consumption)
    gtemp <- select(geocon, -gtarget)
    up <- upaths %>%
      lapply(fusionModel::read_fsd, df = gtemp, cores = cores) %>%
      rbindlist() %>%
      #fcountv(cols = names(gtemp), name = "N0", add = TRUE) %>%
      collapse::join(geocon, how = "inner", verbose = FALSE)  # This simply adds the 'gtarget' variable

    # Temp safety check
    stopifnot(all(year %in% up$year))

    # Check geographic coverage in UrbanPop data
    # At least 95% of all block groups for a given 'gtarget' must be present in UrbanPop in order for it to be included in results
    check <- fcountv(geocon, cols = "gtarget") %>%
      left_join(fcountv(unique(up, by = names(geocon)), cols = "gtarget"), by = "gtarget") %>%
      mutate(coverage = N.y / N.x) %>%
      filter(coverage >= 0.95)

    if (nrow(check) > 0) {
      cat("UrbanPop provides coverage for ", nrow(check), " of ", uniqueN(geocon$gtarget), " (", 100 * signif(nrow(check) / uniqueN(geocon$gtarget), 3), "%) of ", gtarget, " within the requested 'area'\n", sep = "")
    } else {
      stop("UrbanPop does not provide any coverage for ", gtarget, " within the requested 'area'")
    }

    # Collapse UrbanPop 'weight' variable by household, gtarget, and rep
    # TO DO: Skip summarize() if 'gtarget' is block group (unnecessary computation)?
    up <- up[gtarget %iin% check$gtarget, list(weight = sum(weight)), by = list(year, hid, gtarget, rep)]

    # Make weights wide
    # Fill missing with NA values here for more efficient creation of 'ncount' below
    up <- data.table::dcast(up, ... ~ rep, value.var = "weight", fill = NA)

    # Set the replicate weight variable names in 'static'
    i <- which(names(up) == "0")
    wvars <- paste0("REP__", seq(0, ncol(up) - i))
    setnames(up, old = i:ncol(up), new = wvars)

  } else {

    # Report using ACS weights
    cat("Using ACS weights for", paste(unique(range(year)), collapse = "-"), "period\n")

  }

  #---

  # TO DO: Rename
  # Unique year-hid identifiers in the urbanpop or state-puma identifiers in the geocon data
  up.hid <- if (use.up) {
    up %>%
      select(year, hid) %>%
      unique() %>%
      setkey()
  } else {
    geocon %>%
      select(state, puma10, gtarget) %>%
      unique() %>%
      setkey()
  }

  rm(geocon)

  #-----

  # Use assemble() to load ACS variables, possibly including replicate weights
  cat("Loading variables from ACS microdata\n")
  w <- if (use.up) NULL else c('weight', if (R > 0) paste0("rep_", 1:min(R, 80)) else NULL)
  static <- fusionModel::assemble(year = year,  # Better way to automate this?
                                  var = c(avars, fvars, by, w),  # Loads any ACS variables it can
                                  respondent = rtype,
                                  df = up.hid,
                                  cores = cores,
                                  source = "ACS",
                                  silent = TRUE)

  #-----

  # If Urban pop, merge static with 'up' weight and transform weights into wide columns
  if (use.up) {

    # Merge 'static' with UrbanPop weight 'up'
    static <- static %>%
      select(-weight) %>%  # Remove static 'weight', since UrbanPop data has its own primary weight variable
      collapse::join(up,
                     on = intersect(names(static), c('year', 'hid')),
                     how = "inner",
                     multiple = TRUE,
                     verbose = FALSE) %>%
      select(-any_of(gtarget)) %>% # Remove existing column for gtarget, if present (e.g. "puma10" in case of force_up = TRUE)
      setnames(old = "gtarget", new = gtarget) # Rename 'gtarget' to the actual target geography (e.g. 'tract10')

    rm(up)

  } else {

    # Add the geographic target in ACS-only case
    if (!allNA(up.hid$gtarget)) {
      static <- static %>%
        collapse::join(up.hid, how = "inner", on = c('state', 'puma10'), verbose = FALSE) %>%   # This simply adds the 'gtarget' variable
        select(-any_of(gtarget)) %>%  # Remove 'puma10', for example, if it is the target geographic variable
        setnames(old = "gtarget", new = gtarget) # Rename 'gtarget' to the actual target geography (e.g. 'puma10')
    }

    # If using ACS replicate weights, simply rename them
    # Set weight variables names for 'static' in ACS replicate weights case
    wvars <- paste0("REP__", seq(0, length(w) - 1))
    setnames(static, old = w, new = wvars)

  }

  #-----

  # Load any required fusion variables, possibly with multiple implicates
  vsim <- setdiff(c(avars, fvars, by), names(static))
  if (length(vsim)) {
    cat("Loading variables from fused microdata\n")
    sim <- fusionModel::assemble(year = unique(static$year),
                                 var = vsim,
                                 respondent = rtype,
                                 M = M,
                                 df = select(static, year, hid) %>% unique(),
                                 cores = cores,
                                 source = "fused",
                                 silent = TRUE)
  } else {
    # In case where ONLY ACS variables are used, this creates a placeholder 'sim' to avoid code breakage
    sim <- static %>%
      select(year, hid) %>%
      mutate(M = 1L)
  }

  #-----

  # Check that input dimensions are consistent with one another
  # Add here: Something about 'year' being processed/included
  Mv <- sim$M
  Mimp <- max(Mv)
  N <- nrow(sim) / Mimp
  nM <- collapse::fcountv(Mv)
  stopifnot(!is.unsorted(sim$M))
  stopifnot(all(nM$M %in% seq_len(Mimp)))
  stopifnot(all(nM$N == N))

  #-------

  # For 'sim' and 'static' data inputs to be compatible, they must have the same number of rows per implicate
  # This section ensures that the row-order of 'sim' and 'static' are in agreement w.r.t 'hid' and 'year'

  # Variables required to perform subsequent analyses
  req <- c(avars, fvars, by, wvars)

  if (use.up) {

    # The issue here is that 'static' (in the UrbanPop case) can assign a given household (year, hid) to multiple target geographies.
    # However, 'sim' contains just one record for each household-implicate (M) combination
    # We need to expand 'sim' so that the ordering of households within each implicate matches the ordering within 'static'

    # Much faster than match() with manual ID for large data
    ind <- collapse::fmatch(fselect(static, year, hid), fsubset(sim, M == 1, year, hid))

    # Create row indices for 'sim'
    isim <- lapply(0:(Mimp - 1), function(x) ind + N * x)  # Expand indices to include all 'Mimp' implicates

    # IS THIS REALLY NECESSARY?
    # Create row indices for 'static'
    istatic <- (1:nrow(static))[!is.na(ind)]

    # Safety check
    stopifnot(sum(lengths(isim)) / length(istatic) == Mimp)

    # Expand 'sim' data
    # The use of 'isim' ensures 'sim' has row ordering that matches 'static'
    v <- intersect(names(sim), c("M", req))
    sim <- ss(sim, i = unlist(isim), j = v, check = FALSE)

    # Expand 'static' data
    # If same variable name appears in both 'sim' and 'static', use version in 'sim'
    # The use of 'istatic' ensures 'static' has row ordering that matches 'sim'
    v <- intersect(names(static), req)
    v <- setdiff(v, names(sim))
    static <- ss(static, i = istatic, j = v, check = FALSE)

    # Safety check on dimensions
    stopifnot(nrow(sim) / nrow(static) == Mimp)
    rm(ind, isim, istatic)

  } else {

    # In ACS case, simply restrict to minimum set of required variables
    v <- intersect(names(sim), c("M", req))
    sim <- get_vars(sim, v)

    v <- intersect(names(static), req)
    v <- setdiff(v, names(sim))
    static <- get_vars(static, v)

  }

  #---

  # Combine 'sim' and any non-weight variables from 'static'
  # After this operation, 'sim' contains the necessary analysis variables and 'static' contains only the weight variables
  v <- setdiff(names(static), wvars)
  if (length(v)) {
    add_vars(sim) <- alloc(get_vars(static, v), Mimp) %>% rbindlist()
    get_vars(static, v) <- NULL
  }

  # Convert non-positive weights to NA and convert double to integer for memory efficiency
  for (v in names(static)) {
    x <- static[[v]]
    x[x <= 0] <- NA
    if (is.double(x)) x <- as.integer(round(x))
    set(static, j = v, value = x)
  }

  #-------

  # TO DO: Better reporting to console
  n <- length(analyses)
  cat(n, ifelse(n == 1, "analysis", "analyses"), "to perform\n")
  cat("across", Mimp, "implicates\n")
  cat("using", length(wvars) - 1, "replicate weights\n")

  #-------

  # Apply the custom user 'fun'

  if (is.function(fun)) {

    cat("Applying user fun() to microdata\n")

    # Apply function, optionally providing observation weights
    if ("weight" %in% fargs) {
      temp <- fun(data = sim, weight = rep(static$REP__0, Mimp))
    } else {
      temp <- fun(data = sim)
    }

    # Current safety check: fun() must return same number of rows
    stopifnot(nrow(temp) == nrow(sim))

    # If there are identically-named variables in 'temp' and 'sim', remove from latter
    drop <- intersect(names(temp), names(sim))
    if (length(drop)) get_vars(sim, drop) <- NULL

    # Add the variables creates by fun() to 'sim'
    # Retain only the variables needed to conduct analyses
    add_vars(sim) <- temp
    rm(temp)
    drop <- setdiff(names(sim), c('M', avars, by, gtarget))
    if (length(drop)) get_vars(sim, drop) <- NULL

    # Safety check on dimensions
    stopifnot(nrow(sim) / nrow(static) == Mimp)

    cat("Successfully applied user fun() to microdata\n")

  }

  # Check if all required analysis variables are present in 'sim' prior to evaluating inner expressions
  miss <- setdiff(avars, names(sim))
  if (length(miss)) stop("The following analysis variables are not present in 'sim': ", paste(miss, collapse = ", "))

  #-------

  # 'solo' analyses are those with no inner expression modification (can simply rename the target variable)
  solo <- sapply(aexp, function(x) x %in% names(sim))

  # Evaluation of any "inner expressions" in 'analyses' that create custom variables
  ind <- !solo & !duplicated(anames)
  if (any(ind)) cat("Evaluating inner expressions:\n")
  for (i in which(ind)) {

    # Report to console
    cat(" ~", aexp[[i]], "\n")

    # Create a function for each inner expression
    # https://stackoverflow.com/questions/29243037/conversion-of-expression-to-function-in-r
    f <- function(x) x
    body(f) <- parse(text = aexp[[i]])
    x <- ylist[[i]]
    fmls <- vector(mode = "list", length = length(x))
    names(fmls) <- x
    rlang::fn_fmls(f) <- fmls

    # Add the evaluation result to right-side of 'sim' and name the column
    add_vars(sim) <- do.call(f, args = get_vars(sim, x))
    setnames(sim, old = ncol(sim), new = names(alist)[i])

    # Remove any variables in 'x' that are not needed for subsequent analyses
    # This is done for optimal memory management
    j <- setdiff(c(which(solo), pmax(i, which(!solo))), i)
    drop <- setdiff(x, unlist(ylist[j]))
    get_vars(sim, drop) <- NULL

  }

  # Rename any solo-variable inner expression variables
  # This simply changes the original name to the "A__X" convention
  ind <- solo & !duplicated(anames)
  if (any(ind)) setnames(sim, old = unlist(ylist[ind]), new = anames[ind])

  # Convert any character type analysis variables to factors for memory efficiency
  # TO DO: Move earlier for better memory efficiency? Would typically just affect gtarget variable?
  if (length(char_vars(sim, return = "names"))) char_vars(sim) <- dapply(char_vars(sim), as.factor)

  #---------------

  # Determine which analyses are categorical and have outer function "mean" or "sum"
  acat <- which(anames %in% cat_vars(sim, return = "names"))
  invalid <- !afun[acat] %in% c("fmean", "fsum")
  if (any(invalid)) stop("User-created categorical variables must use sum() or mean():\n", paste(analyses[invalid], collapse = "\n"))

  # Numerical analyses to be performed
  anum <- setdiff(seq_along(analyses), acat)

  # Compute once for use below
  ind <- rep.int(seq_row(static), times = nrow(sim) / nrow(static))

  #---------------

  # Do the categorical variable calculations
  if (length(acat)) {

    cat("Computing estimates for categorical analyses:\n ~", paste(sapply(analyses[acat], rlang::f_text), collapse = "\n ~ "), "\n")

    # Main calculation
    cout <- lapply(unique(anames[acat]), function(v) {

      # Requested functions (mean or sum) for analysis variable 'v'
      f <- afun[names(afun) %iin% v]

      # Is a count requested for 'v'?
      vcount <- "fsum" %in% f

      # Identify which summary functions should be used by collap()
      # NOTE: For unknown reason, using a list of function -- list(fsum) -- causes collap() to be much slower than using a quoted character vector
      flist <- c('fsum', 'fnobs')
      if (use.up & R > 0 & vcount) flist <- c(flist, 'fvar')

      # Subset for missing values in 'v', if necessary
      if (anyNA(sim[[v]])) {
        ok <- whichNA(sim[[v]], invert = TRUE)
        W <- ss(static, i = ind[ok], j = names(static), check = FALSE)
        grp <- GRP(sim[ok], by = c("M", by, v))
      } else {
        W <- ss(static, i = ind, j = names(static), check = FALSE)
        grp <- GRP(sim, by = c("M", by, v))
      }

      # Do grouped calculations
      out <- suppressWarnings(
        collapse::collap(X = W,
                         by = grp,
                         FUN = flist,
                         nthreads = cores, # This is passed to underlying fast function (e.g. fmean), if applicable
                         stable.algo = FALSE, # This is passed to fvar() to use faster but numerically unstable variance computation (see ?fvar Details)
                         give.names = TRUE)
      )

      # List of measure variables passed to melt()
      mvars <- list(EST = grep("^fsum", names(out), value = TRUE),
                    N = grep("^fnobs", names(out), value = TRUE),
                    VAR = grep("^fvar", names(out), value = TRUE))
      mvars <- mvars[lengths(mvars) > 0]

      # NOTE: Warning message from data.table about 'measure.vars is a list with length=1' is OK -- should double-check code once expected package changes go into effect
      out %>%
        melt(measure.vars = mvars, variable.name = "REP") %>%
        fmutate(REP = as.integer(REP) - 1L,
                #REP = if (is.factor(REP)) {as.integer(gsub("_", "", str_sub(REP, start = -2)))} else {REP - 1L},  # Check this when data.table package is updated!
                #REP = as.integer(gsub("^REP__", "", REP)) - 1L, # NOTE: In future release of data.table, the gsub() call here will become unnecessary
                ANALYSIS = rep(list(paste(f, v, sep = ".")), nrow(.)),
                aname = factor(v)) %>%
        setnames(old = v, new = "level")

    }) %>%
      rbindlist()

    #---

    # Expand the 'ANALYSIS' list column to account for possible presence of both mean/proportion and sum/count of same column in 'sim'
    # This expands 'cout' to include a 2nd copy of the summed weights, if necessary
    # Identical to: unnest(cout, cols = ANALYSIS) -- but faster
    u <- unlist(cout$ANALYSIS)
    l <- lengths(cout$ANALYSIS)
    if (any(l > 1)) cout <- cout[rep(seq_row(cout), times = l)]
    cout[, ANALYSIS := factor(u)]
    rm(u, l)

    #---

    # Complete the data by including zero estimates for unobserved combinations
    # Correct tidyr result (perhaps sorted)
    #cout <- tidyr::complete(cout, tidyr::nesting(!!!rlang::syms(by), M, REP), tidyr::nesting(ANALYSIS, level), fill = list(EST = 0))

    # Alternative approach (faster)
    # COULD wrap this in its own generic function
    # https://stackoverflow.com/questions/25888706/r-data-table-cross-join-not-working
    CJ.dt = function(X, Y) {
      k = NULL
      X = X[, c(k=1, .SD)]
      setkey(X, k)
      Y = Y[, c(k=1, .SD)]
      setkey(Y, NULL)
      X[Y, allow.cartesian=TRUE][, k := NULL][]
    }
    u1 <- cout %>%
      get_vars(c(by, "M", "REP")) %>%
      funique()
    u2 <- cout %>%
      get_vars(c("ANALYSIS", "level")) %>%
      funique()
    temp <- CJ.dt(u1, u2)
    temp <- fsetdiff(temp, get_vars(cout, names(temp)))
    temp <- add_vars(temp, EST = alloc(0L, nrow(temp)), VAR = alloc(0, nrow(temp)), N = alloc(0L, nrow(temp)))
    cout <- rbind(cout, temp, fill = TRUE)

    #---

    # Calculate estimate and standard variance, based on the type of analysis (proportion or count) requested

    # Which rows in 'cout' refer to mean/proportion analyses?
    i <- grepl("^fmean", cout$ANALYSIS)

    # For proportions...
    if (any(i)) {

      # Coerce all estimates to double to avoid type conflicts when mixing integer and double estimates
      cout[, EST := as.double(EST)]

      # Calculate proportions and total number of observations (summed across levels), by group
      cout[i, `:=`(EST = EST / sum(EST), n = sum(N)), by = c('M', 'REP', 'ANALYSIS', by)]

      # For proportions, divide by 'n' (number of observations in each group)
      cout[i, VAR := EST * (1 - EST) / n]

    }

    # For counts...
    # The estimates are unaffected
    # The standard variance is the raw variance multiplied by 'N' (number of observations in each group-level)
    if (!all(i)) cout[!i, VAR := VAR * N]

    # Retain only necessary variables
    keep <- c('M', 'REP', by, 'ANALYSIS', 'level', 'N', 'EST', 'VAR')
    cout <- cout[, ..keep]

  } else {
    cout <- data.table()
  }

  #----------

  # Remove unnecessary columns in 'sim' before proceeding
  drop <- unique(anames[setdiff(acat, anum)])
  get_vars(sim, drop) <- NULL

  #----------

  # Process the 'anum' analyses in a collap() call

  if (length(anum)) {

    cat("Computing estimates for numerical analyses:\n ~", paste(sapply(analyses[anum], rlang::f_text), collapse = "\n ~ "), "\n")


    # Set any 'Inf' values in analysis variables to NA
    for (v in unique(anames[anum])) {
      if (anyv(sim[[v]], Inf)) {
        i <- whichv(sim[[v]], Inf)
        set(sim, i = i, j = v, value = NA)
        warning("Set ", length(i), " infinite values (", signif(100 * length(i) / nrow(sim), 3), "%) to NA for analysis variable ", v, " in 'sim'")
      }
    }

    #---

    # Compute number of non-NA observations
    ncount <- lapply(unique(anames[anum]), function(v) {

      # Subset for missing values in 'v', if necessary
      if (anyNA(sim[[v]])) {
        ok <- whichNA(sim[[v]], invert = TRUE)
        W <- ss(static, i = ind[ok], j = names(static), check = FALSE)
        grp <- GRP(sim[ok], by = c("M", by))
      } else {
        W <- ss(static, i = ind, j = names(static), check = FALSE)
        grp <- GRP(sim, by = c("M", by))
      }

      # Calculate number of non-NA observations
      fnobs(W, g = grp) %>%
        cbind(grp$groups) %>%
        data.table::melt(id.vars = c('M', by), variable.name = "REP", value.name = "N") %>%
        mutate(aname = factor(v),
               REP = as.integer(gsub("^REP__", "", REP)))

    }) %>%
      rbindlist()

    #---

    # Create the grouping object
    grp <- GRP(sim, by = c("M", by), sort = FALSE)

    # Retain only the analysis variables; 'M' and 'by' variables not necessary since this information is stored in 'grp'
    #sim <- get_vars(sim, vars = setdiff(names(sim), c('M', by)))

    # Create list to pass to 'custom' argument of collap()
    fn <- afun[anum]
    fnames <- unique(fn)
    flist <- lapply(fnames, function(f) {
      v <- names(which(fn == f))
      match(v, names(sim))
    })
    names(flist) <- fnames

    # Add calculation of variance via collapse::fvar()
    if (use.up & R > 0) flist$fvar <- unique(unlist(flist, use.names = FALSE))

    # Once 'ncount' is created, it may be necessary to set NA in 'static' to zeros
    # This is only necessary if the outer functions include 'fmedian' (will throw error); NA's are allowed for the other outer functions
    if (any(afun == "fmedian")) setv(static, NA, 0L)

    #---

    # Major group-wise computations via collap()
    nout <- lapply(wvars, function(w) {
      collapse::collap(X = sim,
                       by = grp,
                       w = ss(static, i = ind, j = w, check = FALSE)[[1]],
                       custom = flist,
                       keep.w = FALSE,
                       nthreads = cores,  # This is passed to underlying fast function (e.g. fmean), if applicable
                       stable.algo = FALSE, # This is passed to fvar() to use faster but numerically unstable variance computation (see ?fvar Details)
                       give.names = TRUE) %>%
        fmutate(REP = as.integer(gsub("^REP__", "", w)))  # Add 'REP' weights identifier and return result
    }) %>%
      rbindlist()  # Combine results and convert to data.table

    #---

    # Integer output columns (converted below)
    intv <- setdiff(names(which(vtypes(nout) == "integer")), c("M", "REP", by))

    # List of measure variables passed to melt()
    temp <- grep("^fvar", names(nout), value = TRUE)
    mvars <- list(EST = setdiff(names(nout), c('M', 'REP', by, temp)))
    if (length(temp)) mvars$VAR = str_replace_all(mvars$EST, "^fmean.", "fvar.") %>% str_replace_all("^fmedian.", "fvar.") %>% str_replace_all("^fsum.", "fvar.")

    # Reshape the results to long format and assemble output
    # NOTE: Warning message from data.table about 'measure.vars is a list with length=1' is OK -- should double-check code once expected package changes go into effect
    nout <- nout %>%
      ftransformv(vars = intv, FUN = as.double) %>%  # Convert integer columns to double to avoid class conflicts in melt()
      melt(measure.vars = mvars, variable.name = "ANALYSIS") %>%
      fmutate(ANALYSIS = if (length(mvars) == 1) {ANALYSIS} else {mvars$EST[ANALYSIS]},
              aname = str_remove(ANALYSIS, "^fmean.") %>% str_remove("^fmedian.") %>% str_remove("^fsum."),
              level = NA) %>%    # Placeholder for compatibility with categorical analyses
      join(ncount, on = c('M', 'REP', 'aname', by), verbose = FALSE)

    rm(grp, intv, mvars, ncount)

    #---

    # Calculate the standard variance of each estimate

    # Placeholder NA column if no variance calculation required
    if (!"VAR" %in% names(nout)) nout[, VAR := NA_integer_]

    # For means, the standard variance is var(x) / N
    nout[grepl("^fmean", nout$ANALYSIS), VAR := VAR / N]

    # For medians, we assume the standard error is sqrt(pi / 2) times the standard error of the mean
    # This means with adjust the standard variance by factor of pi/2
    # See here: https://stats.stackexchange.com/questions/59838/standard-error-of-the-median/61759#61759
    # NOTE: There is a technique to approximate the confidence interval for a percentile of any distribution, but it is much slower than the adjusted variance approximation.
    # see here: https://stats.stackexchange.com/questions/99829/how-to-obtain-a-confidence-interval-for-a-percentile/284970#284970
    nout[grepl("^fmedian", nout$ANALYSIS), VAR := (pi / 2) * VAR / N]

    # For sums, the standard variance is var(x) * N
    nout[grepl("^fsum", nout$ANALYSIS), VAR := VAR * N]

    #---

    # Retain only necessary variables
    keep <- c('M', 'REP', by, 'ANALYSIS', 'level', 'N', 'EST', 'VAR')
    nout <- nout[, ..keep]

  } else {
    nout <- data.table()
  }

  #----------

  # Remove the principle data objects
  rm(sim, static, ind)

  # Combine analysis output data frames
  # TO DO -- IMPROVE NAMES!
  result <- rbind(nout, cout, fill = TRUE)
  rm(nout, cout)

  #---

  # Combine everything and compute final results

  if (R == 0) {
    cat("Computing final point estimates (margin of error not possible when R = 0)\n")
  } else {
    cat("Computing final point estimates and margin of error\n")
  }

  #---

  # ubar: Mean of the observation-based (i.e. classical) standard variances
  # b: Variance of the across-replicate point estimates (zero in case of only one replicate)

  # Compute final estimates
  result <- result %>%
    filter(!is.na(EST)) %>%  # Necessary in case an estimate was impossible for a particular sub-group
    group_by_at(c(by, "level", "ANALYSIS", "M")) %>%

    # Summarize estimate and variance across weights (REP), for each subgroup-implicate
    summarize(

      # Number of estimates, which is equivalent to the available weights (i.e. 1 + the number of replicate weights)
      # The number of estimates can vary across groups depending on the structure of missing data and zero weights
      Rn = n(),

      # Mean number of microdata observations
      N = mean(N),

      # Compute the central estimate, across weights (REP)
      est = if (use.up) {
        # When using UrbanPop, due to presence of zero weights, there is no guarantee that the primary weight (REP__0) returns a valid estimate
        # Consequently, prefer to calculate the estimate as the mean across all weights for which an estimate is available
        mean(EST)
      } else  {
        # Conventional ACS approach to estimate; use the estimate associated with the primary weight
        EST[REP == 0]
      },

      # Compute the variance of the estimates, across weights (REP)
      # Note that Rn can differ from R + 1 due to possible NA's in 'sim' and the structure of zeros in weights for small population subgroups
      var = if (use.up) {
        # Pooled standard variance (Rubin)
        # ubar = mean(VAR)
        # b = var(EST)
        # ubar + (1 + Rn^(-1)) * b
        mean(VAR) + (1 + Rn^(-1)) * var(EST)
      } else {
        # Conventional ACS approach to variance; compute variances around the primary estimate
        # Equivalent to 'mse = TRUE' in survey::svrepdesign()
        # The denominator of the first term is supposed to be R (i.e. the number of replicate weights); we subtract 1 because Rn includes the primary weight as well
        (4 / (Rn - 1)) * sum((EST[REP != 0] - est) ^ 2)
      },

      # "M" must be the final grouping variable so that it is dropped here
      .groups = "drop_last"

    ) %>%

    summarize(

      Rn = Rn[1],
      Mn = n(),

      # Mean number of microdata observations
      N = round(mean(N)),

      # Pooled standard variance (Rubin)
      ubar = mean(var),
      b = ifelse(Mn == 1, 0, var(est)),
      var = ubar + (1 + Mn^(-1)) * b,

      # Final pooled estimate
      est = mean(est),

      .groups = "drop"

    ) %>%

    mutate(

      # Final standard error (Rubin)
      se = sqrt(var),

      # Uncertainty estimation suppressed (NA) if multiple weights are not present
      se = ifelse(Rn == 1, NA, se),

      # Rubin's 'r' used in calculation of degrees of freedom
      # Set to Inf if ubar = 0 (would return NA otherwise)
      r = ifelse(ubar == 0, Inf, (1 + Mn^(-1)) * b / ubar),

      # Degrees of freedom assuming infinite population (original Rubin 1987 degrees of freedom)
      # Confirmation of 'df' in case of replicate weights only: https://stats.stackexchange.com/questions/380467/degrees-of-freedom-for-successive-differences-weights
      #df = if (Mimp == 1) {pmax(1, length(wvars) - 2L)} else  {(Mimp - 1) * (1 + r^(-1)) ^ 2},
      #df = if (Mimp == 1 & !use.up) R - 1L else (Mn - 1) * (1 + r^(-1)) ^ 2,
      #df = ifelse(r == 0, Rn - 1, df),  # If 'r' is zero, there is no variance across the implicates, so we derive 'df' from the number of weights only
      df = ifelse(Mn == 1, Inf, (Mn - 1) * (1 + r^(-1)) ^ 2),

      # DEPRECATED
      # Alternative code if attempting to implement Barnard and Rubin (1999) finite population correction
      # vm = (Mimp - 1) * (1 + r^(-1)) ^ 2,  # Degrees of freedom assuming infinite sample (original Rubin 1987 degrees of freedom)
      # vobs = (1 - r / (r + 1)) * ((dfcom + 1) / (dfcom + 3)) * dfcom,
      # df = ifelse(Mimp == 1, length(wvars) - 2, (vm^(-1) + vobs^(-1))^(-1)),  # Final degrees of freedom using Barnard and Rubin (1999) finite population correction

      # Final 90% confidence interval (p = 0.95 means a 90% confidence interval)
      # Based on the PUMS "User Verification" estimates provided by Census Bureau, it appears they use a t-score of 1.645 when using replicate weights,
      #  but using the 'df' calculation (above) is more correct as explained here: https://stats.stackexchange.com/questions/380467/degrees-of-freedom-for-successive-differences-weights
      moe = se * suppressWarnings(qt(p = 0.95, df)),

      # Coefficient of variation (https://sites.tufts.edu/gis/files/2013/11/Amercian-Community-Survey_Margin-of-error-tutorial.pdf)
      cv = 100 * (moe / 1.645) / est,

      # Calculate share of uncertainty attributable to replicate weights (rshare)
      # Identical to rshare below in most cases
      rshare = sqrt(ubar / se ^ 2),

      # Calculate share of MOE attributable to replicate weights (rshare)
      # The replicate-only MOE is equivalent to using just 'ubar' with 'df' equal to the number of ACS replicate weights minus 1
      # i.e. we assume no variance in estimates across implicates (no modeling uncertainty)
      # The use of pmax() simply ensures that the 'rshare' value cannot exceed 1, due to difference in degrees of freedom
      #rshare = if (Mimp == 1) {1L} else {(sqrt(ubar) * qt(p = 0.95, df = pmax(df, length(wvars)- 2L))) / moe}

    ) %>%

    mutate(i = fmatch(ANALYSIS, alabel), # Returns the analysis number (index); possibly multiple, associated with each ANALYSIS
           lhs = alhs[i],
           rhs = arhs[i],
           aname = anames[i],
           type = sub("f", "", afun[i]),            # Determine the type of analytical result returned
           type = ifelse(!is.na(level) & type == "mean",  "prop", type),
           type = ifelse(!is.na(level) & type == "sum",  "count", type)) %>%

    arrange(i, !!!rlang::syms(by), level) %>%  # Arranges rows according to original 'analyses' order
    select(lhs, rhs, type, all_of(by), level, N, est, moe, se, df, cv, rshare) %>%
    #select(lhs, rhs, type, all_of(by), level, N, ubar, b, r, est, moe, se, df, cv, rshare) %>%
    mutate_all(tidyr::replace_na, replace = NA) %>%  # Replaces NaN from zero division with normal NA
    mutate_if(is.double, convertInteger, threshold = 1)

  #-----

  # Report processing time
  tout <- difftime(Sys.time(), t0)
  cat("Total processing time:", signif(as.numeric(tout), 3), attr(tout, "units"), "\n", sep = " ")

  return(result)

}
