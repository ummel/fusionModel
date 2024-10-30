#' Analyze fusion output to produce small area estimates
#'
#' @description
#' Calculation of point estimates and associated margin of error for small area analyses using variables fused to ACS microdata.
#' ORNL UrbanPop synthetic population is used to re-weight the ACS microdata for individual small areas.
#' Allowable small area geographic units currently limited to: block groups, Census tracts, zip codes, and counties.
#' Efficiently computes means, proportions, sums, counts, medians, standard deviations, and variances, optionally across population subgroups.
#'
#' @param analyses List. Specifies the desired analyses. See Details and Examples. Variables referenced in `analyses` must be in `implicates` or associated ACS microdata.
#' @param implicates Data frame or file path. Implicates of synthetic (fused) variables; typically the output from \link{fuse}. Can be more than one file path to use multiple implicates.
#' @param urbanpop File path to UrbanPop synthetic population data.
#' @param by Character. Optional column name(s) in \code{implicates} or \code{static} (typically factors) that collectively define the set of population subgroups for which each analysis is executed. If \code{NULL}, analysis is done for the whole sample.
#' @param area Specify a geographic area within which to perform the \code{analyses}. Useful for restricting the study area to a manageable size given local computing resources.
#' @param var_scale Scalar. Factor by which to scale the unadjusted replicate weight variance. This is determined by the survey design. The default (\code{var_scale = 4}) is appropriate for ACS and RECS.
#' @param cores Integer. Number of cores used for multithreading in \code{\link[collapse]{collapse-package}} functions.
#'
#' @details The final point estimates are the mean estimates across implicates. The final margin of error is derived from the pooled standard error across implicates, calculated using Rubin's pooling rules (1987). The within-implicate standard error's are calculated using the replicate weights and `var_scale`.
#' @details Each entry in the `analyses` list is a \code{\link[stats]{formula}} of the format `Z ~ F(E)`, where `Z` is an optional, user-friendly name for the analysis, `F` is an allowable “outer function”, and `E` is an “inner expression” containing one or more microdata variables. For example:
#' @details `mysum ~ mean(Var1 + Var2)`
#' @details In this case, the outer function is mean(). Allowable outer functions are: mean(), sum(), median(), sd(), and var(). When the inner expression contains more than one variable, it is first evaluated and then `F()` is applied to the result. In this case, an internal variable `X = Var1 + Var2` is generated across all observations, and then `mean(X)` is computed.
#' @details If no inner expression is desired, the `analyses` list can use the following convenient syntax to apply a single outer function to multiple variables:
#' @details `mean = c("Var1", "Var2")`
#' @details The inner expression can also utilize any function that takes variable names as arguments and returns a vector with the same length as the inputs. This is useful for defining complex operations in a separate function (e.g. microsimulation). For example:
#' @details `myfun = function(Var1, Var2) {Var1 + Var2}`
#' @details `mysum ~ mean(myfun(Var1, Var2))`
#' @details The use of sum() or mean() with an inner expression that returns a categorical vector automatically results in category-wise weighted counts and proportions, respectively. For example, the following analysis would fail if evaluated literally, since mean() expects numeric input but the inner expression returns character. But this is interpreted as a request to return weighted proportions for each categorical outcome.
#' @details `myprop ~ mean(ifelse(Var1 > 10 , 'Yes', 'No'))`
#' @details `analyze2()` uses "fast" versions of the allowable outer functions, as provided by \code{\link[collapse]{fast-statistical-functions}} in the `collapse` package. These functions are highly optimized for weighted, grouped calculations. In addition, outer functions mean(), sum(), and median() enjoy the use of platform-independent multithreading across columns when `cores > 1`. Analyses with numerical inner expressions are processed using a series of calls to \code{\link[collapse]{collap}} with unique observation weights. Analyses with categorical inner expressions utilize a series of calls to \code{\link[collapse]{fsum}}.
#'
#' @return A tibble reporting analysis results, possibly across subgroups defined in \code{by}. The returned quantities include:
#' @return \describe{
#'  \item{lhs}{Optional analysis name; the "left hand side" of the analysis formula.}
#'  \item{rhs}{The "right hand side" of the analysis formula.}
#'  \item{type}{Type of analysis: sum, mean, median, prop(ortion) or count.}
#'  \item{level}{Factor levels for categorical analyses; NA or omitted otherwise.}
#'  \item{est}{Point estimate; mean estimate across implicates.}
#'  \item{moe}{Margin of error associated with the 90% confidence interval.}
#'  \item{rshare}{Share of MOE attributable to replicate weights (as opposed to variance across implicates).}
#'  }
#'
#' @references Rubin, D.B. (1987). \emph{Multiple imputation for nonresponse in surveys}. Hoboken, NJ: Wiley.
#'
#' @examples
#' # Build a fusion model using RECS microdata
#' fusion.vars <- c("electricity", "natural_gas", "aircon", "insulation")
#' predictor.vars <- names(recs)[2:12]
#' fsn.path <- train(data = recs, y = fusion.vars, x = predictor.vars)
#'
#' # Generate 30 implicates of the 'fusion.vars' using original RECS as the recipient
#' recipient <- recs[c(predictor.vars, "weight", paste0("rep_", 1:96))]
#' sim <- fuse(data = recipient, fsn = fsn.path, M = 30)
#' head(sim)
#'
#'#-----
#'
#'# Example of custom pre-processing function
#'myfun <- function(v1, v2, v3) v1 + v2 + v3
#'
#'# Various ways to specify analyses...
#'my.analyses <- list(
#'   # Return means for 'electricity' and proportions for 'aircon'
#'   mean = c("electricity", "aircon"),
#'   # Identical to mean = "electricity"; duplicate analyses automatically removed
#'   electricity ~ mean(electricity),
#'   # Simple addition in the inner expression
#'   mysum ~ sum(electricity + natural_gas),
#'   # Standard deviation of electricity
#'   sd = "electricity",
#'   # Unnamed analyses (no left-hand side in formula)
#'   ~ var(electricity + natural_gas),
#'   ~ mean(insulation),  # Proportions
#'   ~ sum(insulation),  # Counts
#'   # Proportions involving manipulation of >1 variable
#'   myprop ~ mean(aircon != "No air conditioning" & insulation < "Adequately insulated"),
#'   # Custom inner function
#'   mycustom ~ median(myfun(electricity, natural_gas, v3 = 100))
#' )
#'
#'# Do the requested analyses, by "division"
#'result <- analyze2(
#'  analyses = my.analyses,
#'  implicates = sim,
#'  static = recipient,
#'  weight = "weight",
#'  rep_weights = paste0("rep_", 1:96),
#'  by = "division"
#')
#'head(result)
#'
#'#-----
#'
#'# To calculate a conditional estimate, set unused/ignored observations to NA
#'# All outer functions execute with 'na.rm = TRUE'
#'# Example: mean natural_gas conditional on natural_gas > 0
#'# data.table::fifelse() is much faster than base::ifelse() for large data
#'result <- analyze2(
#'  analyses = ~mean(data.table::fifelse(natural_gas > 0, natural_gas, NA_real_)),
#'  implicates = sim,
#'  static = recipient,
#'  weight = "weight",
#'  rep_weights = paste0("rep_", 1:96),
#'  by = "division"
#')
#' @export

#-----

# TEST VALUE
# library(collapse)
# library(tidyverse)
# library(data.table)
# source("R/utils.R")

# result <- analyze2(~sum(electricity > 10e3),
#                    implicates = sim,
#                    static = recs,
#                    weight = "weight",
#                    rep_weights = paste0("rep_", 1:96),
#                    by = "division")

# Check raw fusion data
# m1 <- read_parquet("urbanpop/fusion/RECS_2020/year=2015/M=1/part-0.parquet")
# m2 <- read_parquet("urbanpop/fusion/RECS_2020/year=2015/M=2/part-0.parquet")

# urban pop testing
#analyses = ~mean(fsstatus)
#analyses = ~mean(btuel)

#implicates <- list.files("~/Documents/Projects/fusionModel/ORNL Concept Paper/fusion", pattern = ".parquet$", full.names = TRUE)
#implicates <- "urbanpop/fusion/RECS_2020"
#fvars <- open_dataset(implicates) %>% names()

# TESTING 'area' argument
#area <- state_name == "Alabama" | state %in% c("02", "03")  # Actual input
#area <- substitute(state_name == "Alabama" | state %in% c("02", "04"))  # Input for testing
#area <- substitute(state_name == "Colorado")  # Input for testing
#area <- substitute(state == 36)
#area <- substitute(cbsa10 == 10420)  # Akron, Ohio
#area <- NULL

#by = "tract10"
#by = "zcta10"
#by = c("tract")
#by = "bg"
#by = c("np", "puma")
#by = "county10"
#by = "state"
#by = "puma10"

# Can var_scale be eliminated as argument, since we know we are using ACS weights?
#var_scale = 4  # NOT sure if this is accurate in urbanpop case (need to do some tests)
#cores = 2

# NEW TESTING: 10/24/24
# analyses = list(~mean(btuel), myvar ~ mean(hotma))
# implicates = c("~/Documents/Projects/fusionData/fusion/RECS/2020/2015/RECS_2020_2015_H_fused.fsd",
#                "~/Documents/Projects/fusionData/fusion/RECS/2020/2015/RECS_2020_2015_H_fused.fsd")
# by = "puma10"
# area = substitute(state_name == "Colorado")
# cores = 2

#-------------------------------
#-------------------------------

# TO DO: Option to turn off MOE calculation

analyze3 <- function(analyses,
                     implicates,
                     by = NULL,
                     area = NULL,
                     cores = 1) {

  t0 <- Sys.time()
  fst::threads_fst(nr_of_threads = cores)
  setDTthreads(threads = cores)

  # Check validity of the working directory path
  # Checks if "/fusionData" is part of the path, as this is required
  b <- strsplit(full.path(getwd()), .Platform$file.sep, fixed = TRUE)[[1]]
  i <- which(b == "fusionData")
  if (length(i) == 0) stop("'/fusionData' is not part of the working directory path; this is required.")
  dir <- paste(b[1:i], collapse = .Platform$file.sep)

  # Check the 'area' expression to determine how to parse it
  # If 'area' is passed a 'call', it is not modified
  # This is useful for special user input like: area = str2lang(paste0("state == '", st.obj, "'"))
  # The more common case is for 'area' to follow usage like in subset()
  # See here: http://adv-r.had.co.nz/Computing-on-the-language.html
  check <- try(is.call(area), silent = TRUE)
  if (inherits(check, "try-error") | check == FALSE) {
    area <- substitute(area)
    if (is.character(area)) {
      area <- str2lang(area)
    }
  }

  # Initial argument check
  stopifnot({
    rlang::is_formula(analyses) | is.list(analyses)
    is.data.frame(implicates) | is.character(implicates)  # TO DO: Better check
    is.null(by) | is.character(by)
    is.null(area) | is.call(area)
    cores > 0 & cores %% 1 == 0
  })

  # I believe this can be hard-coded, since replicate weight calculations are always ACS-based
  var_scale <- 4

  # REMOVE?
  # Handle NULL "area" argument
  # See NSE example for subset() function: http://adv-r.had.co.nz/Computing-on-the-language.html
  #if (!is.null(area)) area <- substitute(area)

  # TEMPORARY: Hard-coded path to UrbanPop synthetic population dataset'
  # Parquet dataset directory
  #urbanpop <- "~/Documents/Projects/fusionModel/ORNL Concept Paper/urbanpop/2015-2019/"  # For original Concept paper
  #urbanpop <- "~/Documents/Projects/fusionData/urbanpop/weights"  # For RECS 2020 testing
  urbanpop <- file.path(dir, "urbanpop/weights")

  #---

  # Prepare and check 'analyses' input
  if (!is.list(analyses)) analyses <- list(analyses)

  # Attempt to convert any non-formula entries in 'analyses' into a plausible formula
  # This applies to legacy analysis formulation of the kind:
  #  analyses <- list(mean = c("natural_gas", "aircon"), median = "electricity")
  # The code below converts these to an equivalent formula provided that the function referenced is in .FAST_STAT_FUN
  analyses <- lapply(seq_along(analyses), function(i) {
    x <- analyses[[i]]
    if (!rlang::is_formula(x)) {
      f <- names(analyses)[i]  # The requested outer function
      fobj <- paste("~", f, "(", x, ")")  # No LHS name in this case
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
  y <- unique(unlist(ylist, use.names = FALSE))

  # An analysis variable cannot be in 'by'
  err <- intersect(y, by)
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
    # # Get any argument to outer fun; currently only 'prob' for quantile()
    # temp <- as.list(str2lang(rhs))
    # #temp <- temp[intersect(names(temp), rlang::fn_fmls_names(get(f)))]
    # temp <- temp[intersect(names(temp), "prob")]  # Hard-coded to accept only 'prob' argument
    # arg <- paste(names(temp), temp, sep = " = ")
    # if (length(arg)) {
    #   iexp <- sub(paste0(", ", arg), "", iexp)  # Remove 'arg' from the inner expression string
    #   arg <- sub("prob = ", "n = ", arg) # Switch 'prob' argument to 'n' for compatibility with collapse::fnth()
    # }
    list(f, iexp, lhs, rhs)  # Return required elements
  })

  # Inner expression of each analysis
  aexp <- purrr::map_chr(alist, 2)

  # Internal names of the "analysis variables"; "A..1", "A..2", etc.
  anames <- paste0("A..", match(aexp, unique(aexp)))
  names(alist) <- anames

  # Outer function of each analysis; check that the requested function is allowed
  afun <- purrr::map_chr(alist, 1)
  invalid <- !afun %in% c('sum', 'mean', 'median', 'sd', 'var')  # Valid outer functions
  if (any(invalid)) stop("Outer functions must be sum(), mean(), median(), sd(), or var()")

  # Transform outer functions to "fast" versions in 'collapse' package
  afun[] <- paste0("f", afun)

  # 'afun' MUST include the analysis variable ("A..1", etc) as names attribute
  stopifnot(!is.null(names(afun)))

  # NOT USED: Arguments to pass to the outer function in collap()
  # This is only applicable to "prob" argument for "quantile"
  #aarg <- sapply(alist, function(x) ifelse(length(x[[2]]), x[[2]], NA))

  # The "ANALYSIS" label used to identify each analysis (combination of function and analysis variable, separated by single dot)
  alabel <- paste(afun, names(afun), sep = ".")

  # LHS and RHS of each analysis; assigned to the final results
  alhs <- sapply(alist, function(x) ifelse(length(x[[3]]), x[[3]], NA))
  arhs <- purrr::map_chr(alist, 4)

  #-----

  # Create observation weights data frame by aggregating UrbanPop block group weights to the target geography

  # Report target area
  cat("Geographic area:", if(is.null(area)) "national" else deparse(area), "\n")

  # Crosswalk linking block groups to the target geography ('gtarget'), possibly subsetted by 'area' filter argument
  area.vars <- if(is.null(area)) NULL else all.vars(area)
  geocon <- fst::read_fst(path = file.path(dir, "geo-processed/concordance/geo_concordance.fst"),
                          columns = unique(c(area.vars, 'region', 'division', 'state_name', 'state', 'puma10', 'county10', 'tract10', 'bg10', 'zcta10'))) %>%
    mutate(keep = eval(area)) %>%  # Adds 'keep' column only if area is non-NULL
    filter(bg10 != "None") %>%  # Not sure why there are some entries with "None" for 'bg10' variable
    distinct()

  # Identify which of the 'by' variable(s) is geographic
  # TO DO: Catch case of NO geographic 'by' variables!!
  gtarget <- intersect(by, names(geocon))
  stopifnot(length(gtarget) %in% c(0, 1))
  if (length(gtarget) == 0) {
    stop("Argument 'by' must include one of:\n", paste(names(geocon), collapse = "\n"))
  } else {
    cat("Geographic unit of analysis:", gtarget, "\n")
  }

  # Create crosswalk between block groups (geoid) and the target geography (gtarget)
  # What are the allowable geographic variables for purposes of downscaling?
  # Not a comprehensive list: https://www.census.gov/programs-surveys/geography/guidance/geo-identifiers.html
  if (is.null(area)) geocon$keep <- TRUE
  geocon <- geocon %>%
    fmutate(geoid = paste0(state, county10, tract10, bg10),
            gtarget = switch(gtarget,
                             region = region,
                             division = division,
                             state_name = state_name,
                             state = state,
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
    distinct(geoid, state, puma10, gtarget)

  # TO DO: Confirm check if the 'area' filter violates the integrity of any of the 'gtarget' value
  # For example, if area = state == 36 and by = "division", it should fail because 'area' does not contain a complete division

  #-----

  # Determine how many PUMA's are associated with each value of the target geography
  # Used to determine if ACS replicate weights can be used instead of UrbanPop (i.e. if target geography is PUMA's or larger)
  acs.check <- geocon %>%
    distinct(puma10, state, gtarget) %>%
    add_count(puma10, state)

  # Useful?
  # Extract the years contained in the supplied 'implicates' dataset
  # These are the year values that will also be loaded for ACS and/or UrbanPop data
  # bname <- basename(implicates)
  # sim.years <- if (substring(bname, 1, 4) == "year") {
  #   as.integer(substring(bname, 6))
  # } else {
  # sim.years <- implicates %>%
  #   arrow::open_dataset() %>%
  #   distinct(year) %>%
  #   collect() %>%
  #   pull(year)

  sim.years <- implicates %>%
    strsplit(split = "_", fixed = TRUE) %>%
    purrr::map_chr(~rev(.)[3]) %>%
    as.integer() %>%
    unique()

  #---

  # Create 'static' object

  # TO DO: Auto-detect ACS years and H or P based on 'implicates'
  #acs.paths <- list.files("ORNL Concept Paper/acs", pattern = ".parquet", recursive = TRUE, full.names = TRUE)
  #acs.paths <- "~/Documents/Projects/fusionData/survey-processed/ACS"
  acs.paths <- file.path(dir, "survey-processed/ACS")

  # If using native ACS weights...
  # Determine if analysis can be done using PUMA-level ACS replicate weights only (UrbanPop NOT required)
  if (all(acs.check$n == 1)) {

    # TO DO: Apply 0.2 (i.e. 1/5) adjustment factor (at end) to estimates that are totals, when ACS native weights are used
    # Since we are using five sets of 1-year replicate weights, the total weighted population exceeds true total by factor of 5
    # ALT: Could just use the actual 5-year sample replicate weights? Does it matter? That would require additional data to be processed/loaded.

    # Necessary???
    # acs.geo <- switch(gtarget,
    #                   #region = "region",  # Turned off for now
    #                   #division = "division", # Turned off for now
    #                   state = "state",
    #                   puma10 = c("state", "puma10"))

    # Add required ACS variables to 'static'
    # This code is a bit of a kluge, since it works with the ACS data in .fst format (see below for code if they are in .parquet format)
    static <- lapply(sim.years, function(year) {
      x <- list.files(acs.paths, pattern = paste0(year, "_H_processed.fst"), recursive = TRUE, full.names = TRUE)
      d <- fst::fst(x)
      d <- d[, intersect(names(d), c(y, by, "hid", "year", "state", "puma10", "np", "weight", paste0("rep_", 1:80))), drop = FALSE]
      #d$year <- as.integer(year)
      #names(d)[1] <- "hid"
      #if (year < 2018) d$hid <- as.integer(substring(d$hid, first = 2)) # Remove the leading "1" from the ACS integer ID, since urbanPop does not distinguish between HU and GW prior to 2018
      return(d)
    }) %>%
      bind_rows()

    # ALT: Using Arrow database for ACS microdata...
    # static <- acs.paths %>%
    #   arrow::open_dataset() %>%
    #   select(year, hid, any_of(c(y, by, "state", "puma10", "np")), weight, starts_with("rep_")) %>%
    #   inner_join(acs.check %>% select(-n)) %>%
    #   select(-all_of(acs.geo)) %>%
    #   collect() %>%
    #   setnames(old = c("weight", paste0("rep_", 1:80)), new = paste0("REP__", 0:80)) %>%
    #   setnames(old = "gtarget", new = gtarget)

    # Add the 'gtarget' variable via merge on state and PUMA
    static <- static %>%
      inner_join(select(acs.check, -n), by = join_by(puma10, state)) %>%
      select(-state, -puma10) %>%
      setnames(old = "gtarget", new = gtarget)

    # Update names of weight variables
    # This is because the processed ACS microdata has a 'weight' variable for central weight and then 'rep_1', 'rep_2', etc. for the replicates
    static <- setnames(static, old = c("weight", paste0("rep_", 1:80)), new = paste0("REP__", 0:80))

    # TO DO: Document...
    # Divide by number of 'sim.years' to account for aggregation of annual PUMS files
    pop.check <- static %>%
      group_by_at(by) %>%
      summarize(N = n(),  # Number of sampled households
                hh = as.integer(sum(REP__0) / length(sim.years)),  # Households in the population
                pop = as.integer(sum(REP__0 * np) / length(sim.years)),  # Approximate number of persons in the population
                pshare = NA, # NA because we don't have any GQ data to add here (just using occupied households; if calculated, pshare would be very high for the likely target geographies)
                .groups = "drop")

    # NOTE: When 'area = NULL' and the data are national, the total 'pop' returned from PUMS HH data seems a bit low (maybe 5%). Why?
    # https://www.census.gov/quickfacts/fact/table/US
    #sum(pop.check$pop) / 1e6

    # Drops 'np' if it is not one of the 'y' or 'by' variables
    if (!"np" %in% c(y, by)) static$np <- NULL
    #select(year, hid, any_of(by), starts_with("REP__"))

    # If uncertainty turned off....
    # !!!! TEMPORARY !!!! Restricting to single ACS weight (REP__0), because we aren't doing MOE estimation with current test runs (since UrbanPop has only one replicate, anyway)
    #static <- static %>% select(-any_of(paste0("REP__", 1:80)))

    # Weight variable names
    wvars <- grep("^REP__", names(static), value = TRUE)

  } else {

    # TEMP
    #bup <- geocon

    # Adjust variables in 'geocon' to match format in current national UrbanPop dataset
    geocon <- geocon %>%
      select(-state, -puma10) %>%  # Can remove 'state' and 'puma10' since using UrbanPop
      tidyr::separate_wider_position(geoid, widths = c(state = 2, county = 3, tract_bg = 7)) %>%  # This function is experimental but works for now...
      mutate_at(vars(state, county, tract_bg), as.integer)

    # Use UrbanPop dataset to get assignment of household ID (hid) to the target geography (gtarget), possibly for multiple UrbanPop replicate weights
    # TO DO: Skip summarize() if 'gtarget' is block group (unnecessary computation)?
    static <- urbanpop %>%
      arrow::open_dataset() %>%
      #inner_join(geocon %>% select(geoid, gtarget), by = join_by(geoid)) %>%
      #group_by(year, hid, gtarget, rep) %>%
      inner_join(geocon %>% select(gtarget, state, county, tract_bg), by = join_by(state, county, tract_bg)) %>%
      mutate(rep = 0L) %>% # !!! TEMPORARY: Set 'rep' variable set to zero (central weight) because current national UrbanPop only has single replicate but code below expects a 'rep' column
      group_by(year, hid, gtarget, rep) %>%
      summarize(weight = sum(weight), .groups = "drop") %>%
      collect()

    # TO DO: Could this ACS-related code chunk be pulled out of if/then. Looks same as code chunk used in ACS-only case above...
    # Add required ACS variables to 'static'
    # This code is a bit of a kluge, since it works with the ACS data in .fst format (see below for code if they are in .parquet format)
    static <- lapply(sim.years, function(year) {
      x <- list.files(file.path("survey-processed/ACS", year), pattern = "H_processed.fst", recursive = TRUE, full.names = TRUE)
      d <- fst::fst(x)
      d <- d[, intersect(names(d), c(names(d)[1], "np", y, by)), drop = FALSE]
      d$year <- as.integer(year)
      names(d)[1] <- "hid"
      if (year < 2018) d$hid <- as.integer(substring(d$hid, first = 2)) # Remove the leading "1" from the ACS integer ID, since urbanPop does not distinguish between HU and GW prior to 2018
      return(d)
    }) %>%
      bind_rows() %>%
      #right_join(static, by = join_by(year, hid)) %>%
      collapse::join(static, on = c('year', 'hid'), how = "right", multiple = TRUE, verbose = FALSE) %>%  # Identical to right_join() but a bit faster
      setnames(old = "gtarget", new = gtarget)

    # Easier to have the ACS data stored as as parquet files...
    # static <- acs.paths %>%
    #   arrow::open_dataset() %>%
    #   select(year, hid, any_of(c(y, by, "np"))) %>%
    #   right_join(static, by = join_by(year, hid)) %>%
    #   collect() %>%
    #   setnames(old = "gtarget", new = gtarget)

    # TO DO: Document...
    # Note that 'np' will be NA for group quarter units in UrbanPop that don't have a match in the fusionACS processed ACS microdata
    # The 'pshare' value should be close to 1 for most geographic areas
    pop.check <- static %>%
      filter(rep == 0) %>%  # Using the "central" UrbanPop weights
      group_by_at(by) %>%
      summarize(N = sum(is.finite(np)), # Number of sampled households (non-GQ) for the central weight (rep == 0)
                hh = sum(weight[is.finite(np)]), # Number of households (non-GQ) in the population
                pop = sum(weight * np, na.rm = TRUE),  # Number of individuals (non-GQ) in the population
                pop0 = sum(weight * replace_na(np, 1L))) %>%   # Number of individuals (including GQ) in the population
      mutate(pshare = pop / pop0)

    # Make replicate weights wide
    static <- static %>%
      filter(is.finite(np)) %>%  # Removes group-quarter observations
      select(year, hid, any_of(c(y, by)), rep, weight) %>%  # Drops 'np' unless it is not one of the 'by' or 'y' variables
      qDT() %>%
      data.table::dcast(... ~ rep, value.var = "weight", fill = 0L)

    # Set the replicate weight variable names in 'static'
    i <- which(names(static) == "0")
    wvars <- paste0("REP__", seq(0, ncol(static) - i))
    setnames(static, old = i:ncol(static), new = wvars)

  }

  rm(acs.check)
  rm(geocon)

  #---

  # Load 'implicates' data, restricted to households in 'static' and required variables
  #if (!inherits(sim, "ArrowObject")) stop("Unable to open 'implicates' as an Arrow dataset object")
  # NOTE: This is fairly slow at national scale...faster to use an upfront filter() call on 'hid'?
  # When the geographic extent is large (e.g. national), this is quite slow (if using all implicates) because it effectively has to load everything

  # Arrow implementation
  # sim <- implicates %>%
  #   arrow::open_dataset() %>%
  #   #filter(M == 1) %>%   # !!!!! TEMP -- much faster for testing by restricting to single implicate !!!!!!
  #   select(M, year, hid, any_of(c(y, by))) %>%
  #   inner_join(static %>% select(year, hid) %>% distinct(), by = join_by(year, hid)) %>%
  #   arrange(M, year, hid) %>%
  #   collect()

  # fst/.fsd implementation
  sim <- lapply(implicates, function(x) {
    fusionModel::read_fsd(fsd = x,
                          columns = intersect(c('M', 'year', 'hid', y, by), names(fst::fst(x))),
                          cores = cores) %>%
      #filter(M == 1) %>%   # !!!!! TEMPORARY -- much faster for testing by restricting to single implicate !!!!!!
      inner_join(static %>% select(year, hid) %>% distinct(), by = join_by(year, hid))
  }) %>%
    rbindlist() %>%
    arrange(M, year, hid)

  #---

  # Required variables to perform analysis
  req <- c(y, by, wvars)

  # Ensure all required variables are present in either 'sim' or 'static'
  miss <- setdiff(req, c(names(sim), names(static)))
  if (length(miss)) stop("Could not find required variables in 'implicates' or 'static': ", paste(miss, collapse = ", "))

  #---

  # Check that input dimensions are consistent with one another
  # Add here: Something about the 'years' being processed/included
  Mv <- sim$M
  Mimp <- max(Mv)
  N <- nrow(sim) / Mimp
  nM <- collapse::fcountv(Mv)
  stopifnot(!is.unsorted(sim$M))
  stopifnot(all(nM$M %in% seq_len(Mimp)))
  stopifnot(all(nM$N == N))
  n <- length(analyses)
  cat(n, ifelse(n == 1, "analysis", "analyses"), "to perform\n")
  cat("across", Mimp, "implicates\n")
  cat("using", length(wvars) - 1, "replicate weights\n")  # TO DO: Better text

  #---

  # For 'sim' and 'static' data inputs to be compatible, they must have the same number of rows per implicate
  # This section is critical, since it ensures that the row-order of 'sim' and 'static' are in agreement w.r.t household ID

  # TEMP
  # bup.static <- static
  # bup.sim <- sim

  # The problem here is that 'static' (in the UrbanPop case) can assign a given household (year, hid) to multiple target geographies.
  # However, 'sim' contains just one record for each household-implicate (M) combination
  # We need to expand 'sim' so that the ordering of households within each implicate matches the ordering within 'static'

  # Indicates matching rows in 'static' to first occurrence in 'sim'
  # static <- fmutate(static, ID = paste0(year, hid), .keep = "unused")
  # sim <- fmutate(sim, ID = paste0(year, hid), .keep = "unused") # (!!!) This is quite slow for bigger data -- better way to do the match() call?
  # ind <- match(static$ID, sim$ID[1:N])

  # Much faster than match() with manual ID for large data
  ind <- collapse::fmatch(fselect(static, year, hid), fsubset(sim, M == 1, year, hid))

  # Create row indices for 'sim'
  isim <- lapply(0:(Mimp - 1), function(x) ind + N * x)  # Expand indices to include all 'Mimp' implicates

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

  #---

  # Combine 'sim' and any non-weight variables from 'static'
  # After this operation, 'sim' contains the necessary analysis variables and 'static' contains only the weight variables
  v <- setdiff(names(static), wvars)
  if (length(v)) {
    add_vars(sim) <- alloc(get_vars(static, v), Mimp) %>% rbindlist()
    get_vars(static, v) <- NULL
  }

  # Convert any double weight variables to integer for memory efficiency
  # Not a critical issue if using ACS replicate weights, since they should already be integer
  dbl <- names(which(sapply(static, is.double)))
  if (length(dbl)) get_vars(static, dbl) <- dapply(get_vars(static, dbl), as.integer)

  #----

  # 'solo' analyses are those with no inner expression modification (can simply rename the target variable)
  solo <- sapply(aexp, function(x) x %in% y)

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
  if (length(char_vars(sim, return = "names"))) char_vars(sim) <- dapply(char_vars(sim), as.factor)

  #---

  # Determine which analyses are categorical and have outer function "mean" or "sum"
  acat <- which(anames %in% cat_vars(sim, return = "names"))
  invalid <- !afun[acat] %in% c("fmean", "fsum")
  if (any(invalid)) stop("Inner expressions that return categorical results must use sum() or mean():\n", paste(analyses[invalid], collapse = "\n"))

  # Numerical analyses to be performed
  anum <- setdiff(seq_along(analyses), acat)

  #---

  # Process the 'anum' analyses in a collap() call

  if (length(anum)) {

    cat("Computing estimates for numerical analyses:\n ~", paste(sapply(analyses[anum], rlang::f_text), collapse = "\n ~ "), "\n")

    # Create list to pass to 'custom' argument of collap()
    fun <- afun[anum]
    fnames <- unique(fun)
    flist <- lapply(fnames, function(f) {
      v <- names(which(fun == f))
      match(v, names(sim))
    })
    names(flist) <- fnames

    # Create the grouping object
    grp <- GRP(sim, by = c("M", by), sort = FALSE)

    # Set collapse package option to silently ignore unapplicable 'nthreads' argument
    # See '...' entry in ?collap; applicable to fsd() and fvar()
    options(collapse_unused_arg_action = "none")

    # Do computation via collap()
    nout <- lapply(wvars, function(w) {

      # This is preferable but fails if length(flist) == 1
      collapse::collap(X = sim,
                       by = grp,
                       #w = ss(static, i = sim$ID, j = w, check = FALSE)[[1]],
                       w = rep(ss(static, j = w, check = FALSE)[[1]], times = nrow(sim) / nrow(static)),  # !! TO DO: Check if this is correct!
                       custom = flist,
                       keep.w = FALSE,
                       nthreads = cores,
                       give.names = TRUE) %>%
        fmutate(REP = as.integer(gsub("^REP__", "", w)))  # Add 'REP' weights identifier and return result

    }) %>%
      rbindlist()  # Converts results to data.table

    # Is this necessary? - don't think so
    # TEMPORARY - TO DO: check collapse package github for status?
    # If there is only 1 analysis in 'flist', collap() does not automatically append the function to the output column
    # Have to do it manually -- might be fixed in development version; waiting for it to become official
    # if (length(flist) == 1) {
    #   v <- names(sim)[flist[[1]]]
    #   nout <- setnames(nout, old = v, new = paste(names(flist)[1], v, sep = "."))
    # }

    # Reshape the results to long format
    # fct.to.int <- setdiff(cat_vars(out, return = "names"), by) # Identify factor columns to be converted to integers prior to the melt() call so that the "EST" column will be strictly numeric (no class conflict)
    int.to.num <- setdiff(names(which(sapply(nout, is.integer))), c("M", "REP", by))
    nout <- nout %>%
      ftransformv(vars = int.to.num, FUN = as.numeric) %>%  # Convert integer columns to numeric to avoid class conflicts in melt()
      melt(id.vars = c("M", "REP", by), variable.name = "ANALYSIS", value.name = "EST") %>%
      fmutate(level = NA)  # Placeholder for compatibility with categorical analyses

  } else {
    nout <- data.frame()
  }

  #---

  # Remove unnecessary columns in 'sim' before proceeding
  drop <- unique(anames[anum])
  get_vars(sim, drop) <- NULL

  #---

  # Do the categorical variable calculations
  if (length(acat)) {

    cat("Computing estimates for categorical analyses:\n ~", paste(sapply(analyses[acat], rlang::f_text), collapse = "\n ~ "), "\n")

    # Faster alternative to original dcast() code
    cvars <- unique(anames[acat])
    svars <- setdiff(names(sim), c("M", cvars))
    temp <- ss(sim, i = sim$M == 1, j = svars, check = FALSE)
    #sim <- qDF(rsplit(sim, ~M, use.names = FALSE, cols = cvars))  # Not correct with multiple 'cvars'
    sim <- bind_cols(rsplit(sim, ~M, use.names = FALSE, cols = cvars, simplify = FALSE), .name_repair = "minimal")  # Not sure if this is fastest, but it works.
    #names(sim) <- paste(cvars, 1:Mimp, sep = "_")
    names(sim) <- unlist(lapply(1:Mimp, function(x) paste(cvars, x, sep = "_")))
    add_vars(sim) <- temp
    rm(temp)

    # # OLD code block...
    # Respondent/record ID for each row in 'sim'
    # N <- sum(sim$M == 1)
    # sim <- add_vars(sim, ID = rep(1L:N, Mimp))
    #
    # # I think dcast is somwhat memory inefficient for this particular operation if there are many variables?
    # cvars <- unique(anames[acat])
    # bsim <- intersect(names(sim), by)  # 'by' variables in 'sim'
    # fsim <- as.formula(paste("ID", ifelse(length(bsim), "+", ""), paste(bsim, collapse = "+"), "~M"))
    # sim <- qDT(get_vars(sim, c("ID", "M", by, cvars)))
    # sim <- data.table::dcast(data = sim, formula = fsim, value.var = cvars)
    # if (length(cvars) == 1) setnames(sim, old = as.character(1:Mimp), new = paste(cvars, 1:Mimp, sep = "_"))

    # IS THIS NECESSARY ANY MORE?
    # This internal ss() call is only necessary if dcast() results in more rows in 'sim' than in 'static'
    # When does this occur???
    # sim <- if (nrow(sim) == N) {
    #   add_vars(sim, static)
    # } else {
    #   add_vars(sim, collapse::ss(static, sim$ID, check = FALSE))
    # }
    # rm(static)
    # get_vars(sim, "ID") <- NULL

    # Add 'static' to 'sim'
    add_vars(sim) <- static
    rm(static)

    #---

    #ccols <- unlist(lapply(cvars, function(x) paste(x, 1:Mimp, sep = "_")))
    ccols <- setdiff(names(sim), c(svars, wvars))  # Maybe safe?
    cout <- lapply(ccols, function(v) {
      #for (v in ccols) {
      s <- strsplit(v, "_", fixed = TRUE)[[1]]
      sim %>%
        fgroup_by(c(by, v), sort = FALSE) %>%
        get_vars(wvars) %>%
        fsum(nthreads = cores) %>%
        qDT() %>%
        data.table::melt(id.vars = c(by, v), variable.name = "REP", value.name = "EST") %>%
        na_omit(cols = v) %>% # Remove rows where there is no observable outcome (NA) for the analysis variable
        rename(level = !!v) %>%
        fmutate(ANALYSIS = paste("fsum", s[1], sep = "."),
                M = as.integer(s[2]),
                REP = as.integer(gsub("^REP__", "", REP)))
    }) %>%
      rbindlist()

    # Correct tidyr result (perhaps sorted)
    # Complete the data by including zero estimates for unobserved combinations
    #cout <- tidyr::complete(cout, tidyr::nesting(!!!rlang::syms(by), M, REP), tidyr::nesting(ANALYSIS, level), fill = list(EST = 0))

    # Alternative approach (faster)
    # COULD wrap this in its own generic function
    # TO DO: Move to utils?
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
    temp <- add_vars(temp, EST = alloc(0, nrow(temp)))
    cout <- rbind(cout, temp)

  } else {
    cout <- data.frame()
  }

  #---

  # Remove the principle data objects
  suppressWarnings(rm(sim, static))

  # Combine analysis output data frames
  # TO DO -- IMPROVE NAMES!
  result <- rbind(nout, cout)
  rm(nout, cout)

  #---

  # Combine everything and compute final results
  # TO DO: Adjust wording for absence of MOE calculation!
  cat("Computing final point estimates and margin-of-error\n")

  # 'i0' gives the ANALYSIS values in 'result' that were requested
  # If an analysis wasn't explicitly requested, we know it must be a proportion (i.e. mean requested but sum returned in 'result')
  i0 <- match(unique(result$ANALYSIS), alabel)

  # Determine degrees of freedom and share of population, by each separate analysis
  temp <- result %>%
    filter(REP == 0, M == 1) %>%
    group_by_at(c(by, "ANALYSIS")) %>%
    tally(name = "nlevels") %>%
    ungroup() %>%
    left_join(pop.check, by = by) %>%
    mutate(dfcom = pmax(1, N - nlevels)) %>%  # The 'complete data' degrees of freedom used for Barnard and Rubin finite population correction
    select(all_of(by), ANALYSIS, N, hh, pop, pshare, dfcom)

  # Compute final estimates
  result <- result %>%
    group_by_at(c(by, "level", "ANALYSIS", "M")) %>%
    # Compute metrics across replicates (REP)
    summarize(
      est = if (any(REP == 0)) EST[REP == 0] else mean(EST), # Central estimate is associated with the primary weight, if present; otherwise, mean of the replicate estimates
      var = (var_scale / sum(REP != 0)) * sum((EST[REP != 0] - est) ^ 2),
      .groups = "drop_last"  # "M" must be the final grouping variable so that it is dropped here
    ) %>%
    # Compute metrics across implicates (M)
    summarize(
      ubar = mean(var), # Mean of the within-implicate (replicate-based) variances
      b = if (Mimp == 1) 0 else var(est), # Variance of the across-implicate point estimates (zero if there is only 1 implicate)
      est = mean(est), # Final point estimate; average across the implicates
      .groups = "drop"
    ) %>%
    left_join(temp, by = c(by, "ANALYSIS")) %>%  # Add 'N', 'pop', and 'pshare' columns
    mutate(
      se = sqrt(ubar + (1 + Mimp^(-1)) * b), # Final standard error (Rubin)
      r = ifelse(ubar == 0, Inf, (1 + Mimp^(-1)) * b / ubar),  # Set to Inf if ubar = 0 (would return NA otherwise)
      df = if (Mimp == 1) {length(wvars) - 2L} else  {(Mimp - 1) * (1 + r^(-1)) ^ 2},  # Degrees of freedom assuming infinite sample (original Rubin 1987 degrees of freedom)

      # Alternative code if attempting to implement Barnard and Rubin (1999) finite population correction
      # vm = (Mimp - 1) * (1 + r^(-1)) ^ 2,  # Degrees of freedom assuming infinite sample (original Rubin 1987 degrees of freedom)
      # vobs = (1 - r / (r + 1)) * ((dfcom + 1) / (dfcom + 3)) * dfcom,
      # df = ifelse(Mimp == 1, length(wvars) - 2, (vm^(-1) + vobs^(-1))^(-1)),  # Final degrees of freedom using Barnard and Rubin (1999) finite population correction

      # Final 90% confidence interval (p = 0.95 means a 90% confidence interval)
      moe = se * qt(p = 0.95, df),

      # Calculate share of MOE attributable to replicate weights (rshare)
      # The replicate-only MOE is equivalent to using just 'ubar' with 'df' equal to the number of ACS replicate weights minus 1
      # i.e. we assume no variance in estimates across implicates (no modeling uncertainty)
      # Confirmation of 'df' value: https://stats.stackexchange.com/questions/380467/degrees-of-freedom-for-successive-differences-weights
      # The use of pmax() simply ensures that the 'rshare' value cannot exceed 1, due to difference in degrees of freedom
      rshare = if (Mimp == 1) {1L} else {(sqrt(ubar) * qt(p = 0.95, df = pmax(df, length(wvars) - 2))) / moe}

    ) %>%

    # This returns the analysis number (index); possibly multiple, associated with each ANALYSIS
    # Categorical analyses are all "fsum" at this point, and can be matched against requested analyses (sum, mean, or both)
    # If mean or both requested, the unnest() below creates the correct number of entries for each 'i', replicating as necessary
    mutate(i = lapply(ANALYSIS, function(x) unique(na.omit(c(match(x, alabel), match(sub("fsum", "fmean", x), alabel)))))) %>%
    tidyr::unnest(i) %>%
    group_by_at(c(by, "i")) %>%

    # Adjustment factor for rows to be converted to proportions
    # If 'i' is in 'i0' then no conversion to proportion is necessary; otherwise divide each estimate by the total
    mutate(adj = ifelse(i %in% i0, 1, 1 / sum(est))) %>%
    ungroup() %>%
    mutate(lhs = alhs[i],
           rhs = arhs[i],
           est = est * adj,  # Adjustment for proportions applied here
           moe = moe * adj,
           # Determine the type of analytical result returned
           type = sub("f", "", afun[i]),
           type = ifelse(!is.na(level) & type == "mean",  "prop", type),
           type = ifelse(!is.na(level) & type == "sum",  "count", type)) %>%
    arrange(i, !!!rlang::syms(by), level) %>%  # Arranges rows according to original 'analyses' order
    select(lhs, rhs, type, all_of(by), N, hh, pop, pshare, level, est, moe, se, df, rshare) %>%
    mutate_all(tidyr::replace_na, replace = NA) %>%   # Replaces NaN from zero division with normal NA
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.numeric, cleanNumeric, tol = 0.001) %>%
    mutate(level = if (all(is.na(level))) NULL else as.character(level)) # Remove 'level' if it contains no information

  # Report processing time
  tout <- difftime(Sys.time(), t0)
  cat("Total processing time:", signif(as.numeric(tout), 3), attr(tout, "units"), "\n", sep = " ")

  return(result)

}
