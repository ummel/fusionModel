#' Analyze fusion output
#'
#' @description
#' Calculation of point estimates and associated margin of error for analyses using fused/synthetic microdata with replicate weights. Efficiently computes means, proportions, sums, counts, medians, standard deviations, and variances, optionally across population subgroups.
#' This differs from \link{analyze} in that it requires replicate weights and calculates uncertainty using full replicate weight variance (no approximation).
#'
#' @param analyses List. Specifies the desired analyses. See Details and Examples. Variables referenced in `analyses` must be in `implicates` or `static`.
#' @param implicates Data frame or file path. Implicates of synthetic (fused) variables; typically the output from \link{fuse}. The implicates should be row-stacked and identified by integer column "M". If a file path to a ".fst" file, only the necessary columns are read into memory.
#' @param static Data frame or file path. Static variables that do not vary across implicates; typically the "recipient" microdata passed to \link{fuse}. At a minimum, `static` must contain `weight` and `rep_weights`. If a file path to a ".fst" file, only the necessary columns are read into memory. Note that \code{nrow(static) = nrow(implicates) / max(implicates$M)} and the row-ordering is assumed to be consistent between \code{static} and \code{implicates}.
#' @param weight Character. Name of the primary observation weights column in \code{static}.
#' @param rep_weights Character. Vector of replicate weight columns in \code{static}.
#' @param by Character. Optional column name(s) in \code{implicates} or \code{static} (typically factors) that collectively define the set of population subgroups for which each analysis is executed. If \code{NULL}, analysis is done for the whole sample.
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
#'# Do the requeted analyses, by "division"
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

#-------------------------------
#-------------------------------

analyze2 <- function(analyses,
                     implicates,
                     static,
                     weight,
                     rep_weights,
                     by = NULL,
                     var_scale = 4,
                     cores = 1) {

  t0 <- Sys.time()

  # Initial argument check
  stopifnot({
    rlang::is_formula(analyses) | is.list(analyses)
    is.data.frame(implicates) | is.character(implicates)
    is.data.frame(static) | is.character(static)
    is.character(weight) & length(weight) == 1
    is.null(rep_weights) | is.character(rep_weights)
    is.null(by) | is.character(by)
    var_scale > 0
    cores > 0 & cores %% 1 == 0
  })

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
  if (length(err)) stop("Analysis variables cannot also be in 'by': ", paste(err, collapse = ", "))

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

  # Outer function of each analysis; check that it is allowed
  afun <- purrr::map_chr(alist, 1)
  invalid <- !afun %in% c('sum', 'mean', 'median', 'sd', 'var')  # Valid outer functions
  if (any(invalid)) stop("Outer functions must be sum(), mean(), median(), sd(), or var()")

  # Transform outer functions to "fast" versions, if possible
  # This includes changing "quantile" to "fnth"
  afun[] <- paste0("f", afun)
  #afun[afun == "fquantile"] <- "fnth"  # NOT USED

  # NOTE: 'afun' MUST include the analysis variable ("A..1", etc) as names attribute
  stopifnot(!is.null(names(afun)))

  # Arguments to pass to the outer function in collap()
  # This is only applicable to "prob" argument for "quantile"
  #aarg <- sapply(alist, function(x) ifelse(length(x[[2]]), x[[2]], NA))

  # The "ANALYSIS" label used for each analysis (combination of function and analysis variable, separated by single dot)
  alabel <- paste(afun, names(afun), sep = ".")

  # LHS and RHS of each analysis; assigned to the final results
  alhs <- sapply(alist, function(x) ifelse(length(x[[3]]), x[[3]], NA))
  arhs <- purrr::map_chr(alist, 4)

  #----

  # Get only the variables in 'implicates' that need to be retained/loaded
  sim <- if (is.data.frame(implicates)) {
    qDF(implicates)
  } else {
    if (!endsWith(implicates, ".fsd") & !endsWith(implicates, ".fst")) stop("'implicates' file path must end with .fsd or .fst")
    cat("Reading required 'implicates' data from disk\n")
    if (endsWith(implicates, ".fsd")) fusionModel::read_fsd(implicates, )  # This is slow...
    if (endsWith(implicates, ".fst")) fst::fst(implicates) # This returns just a fst_table (FAST)
  }
  rm(implicates)
  v <- intersect(names(sim), c("M", y, by))
  sim <- if (is.data.frame(sim)) get_vars(sim, v) else qDF(sim[v])

  #----

  # Get only the variables in 'static' that need to be retained/loaded
  static <- if (is.data.frame(static)) {
    qDF(static)
  } else {
    if (!endsWith(static, ".fst")) stop("'static' file path must end with .fst")
    cat("Reading required 'static' data from disk\n")
    fst::fst(static)
  }
  v <- intersect(names(static), c(y, by, weight, rep_weights))
  v <- setdiff(v, names(sim))  # If same variable name appears in both 'implicates' and 'static', use version in 'implicates'
  if (!all(c(weight, rep_weights) %in% v)) stop("All 'weight' and 'rep_weights' must be in 'static'")
  static <- if (is.data.frame(static)) get_vars(static, v) else qDF(static[v])

  #----

  # Ensure all required variables are present in either 'sim' or 'static'
  req <- c(y, by, weight, rep_weights)
  miss <- setdiff(req, c(names(sim), names(static)))
  if (length(miss)) stop("Could not find required variables in 'implicates' or 'static': ", paste(miss, collapse = ", "))

  # Standardize name of weight variable in 'static'
  wvars <- paste0("REP__", 0:length(rep_weights))
  setnames(static, c(weight, rep_weights), wvars)

  #----

  # Check that input dimensions are consistent with one another
  Mimp <- max(sim$M)
  N <- nrow(sim)  / Mimp
  nM <- collapse::fcount(sim, M)
  stopifnot(all(nM$M %in% seq_len(Mimp)))
  stopifnot(all(nM$N == N))
  if (!is.null(static)) stopifnot(nrow(static) == N)
  n <- length(analyses)
  cat(n, ifelse(n == 1, "analysis", "analyses"), "to perform\n")
  cat(Mimp, "implicates\n")
  cat(length(rep_weights), "replicate weights\n")

  # Respondent/record ID for each row in 'sim'
  sim <- add_vars(sim, ID = rep(1L:N, Mimp))

  #----

  # Combine 'sim' and any non-weight variables from 'static'
  # After this operation, 'sim' contains the necessary analysis variables and 'static' contains only the weight variables
  v <- setdiff(names(static), wvars)
  if (length(v)) {
    sim <- add_vars(sim, alloc(get_vars(static, v), Mimp) %>% rbindlist())
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

  # Ensure all categorical analysis variables are factors for memory efficiency
  if (length(cat_vars(sim, return = "names"))) cat_vars(sim) <- dapply(cat_vars(sim), as.factor)

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

    # Set collapse package option to silently ignore inapplicable 'nthreads' argument
    # See '...' entry in ?collap; applicable to fsd() and fvar()
    options(collapse_unused_arg_action = "none")

    nout <- lapply(wvars, function(w) {

      # This is preferable but fails if length(flist) == 1
      collapse::collap(X = sim,
                       by = grp,
                       w = ss(static, i = sim$ID, j = w, check = FALSE)[[1]],
                       custom = flist,
                       keep.w = FALSE,
                       nthreads = cores) %>%
        fmutate(REP = as.integer(gsub("^REP__", "", w)))  # Add 'REP' weights identifier and return result

    }) %>%
      rbindlist()  # Converts results to data.table

    # TEMPORARY
    # If there is only 1 analysis in 'flist', collap() does not automatically append the function to the output column
    # Have to do it manually -- might be fixed in development version; waiting for it to become official
    if (length(flist) == 1) {
      v <- names(sim)[flist[[1]]]
      nout <- setnames(nout, old = v, new = paste(names(flist)[1], v, sep = "."))
    }

    # Reshape the results to long format
    # fct.to.int <- setdiff(cat_vars(out, return = "names"), by) # Identify factor columns to be converted to integers prior to the melt() call so that the "EST" column will be strictly numeric (no class conflict)
    int.to.num <- setdiff(names(which(sapply(nout, is.integer))), c("M", "REP", by))
    nout <- nout %>%
      #ftransformv(vars = fct.to.int, FUN = as.numeric) %>%
      ftransformv(vars = int.to.num, FUN = as.numeric) %>%  # Convert integer columns to numeric to avoid class conflicts in melt()
      melt(id.vars = c("M", "REP", by), variable.name = "ANALYSIS", value.name = "EST") %>%
      fmutate(level = NA)  # Placeholder for compatibility with categorical analyses

  } else {
    nout <- data.frame()
  }

  #---

  # # NO LONGER NECESSARY
  # # If a categorical result comes out of collap(), how to process it?
  # # This occurs whenever a FUN like fmode, fmin, fax, works on a categorical analysis variable
  # other.cat <- lapply(fct.to.int, function(v) {
  #   i <- match(v, alabel)  # Analysis number
  #   lev <- levels(sim[[anames[i]]])
  #   other %>%
  #     ss(i = other$ANALYSIS %in% v, check = FALSE) %>%
  #     fmutate(level = factor(lev[EST], levels = lev), EST = 1L) %>%
  #     tidyr::complete(tidyr::nesting(!!!rlang::syms(by), M, REP, ANALYSIS), level, fill = list(EST = 0L))
  # }) %>%
  #   rbindlist()

  # Non-categorical variables returned by collap()
  # other <- other %>%
  #   ss(i = !other$ANALYSIS %in% fct.to.int, check = FALSE)

  # Remove unnecessary columns in 'sim' before proceeding
  drop <- unique(anames[anum])
  get_vars(sim, drop) <- NULL

  #---

  # MOVE THIS To end (after other calculations) for better memory handling
  # Do the categorical variable calculations
  if (length(acat)) {

    cat("Computing estimates for categorical analyses:\n ~", paste(sapply(analyses[acat], rlang::f_text), collapse = "\n ~ "), "\n")

    cvars <- unique(anames[acat])
    bsim <- intersect(names(sim), by)  # 'by' variables in 'sim'
    fsim <- as.formula(paste("ID", ifelse(length(bsim), "+", ""), paste(bsim, collapse = "+"), "~M"))
    sim <- qDT(get_vars(sim, c("ID", "M", by, cvars)))
    sim <- data.table::dcast(data = sim, formula = fsim, value.var = cvars)
    if (length(cvars) == 1) setnames(sim, old = as.character(1:Mimp), new = paste(cvars, 1:Mimp, sep = "_"))

    # This internal ss() call is only necessary if melt() results in more rows in 'sim' than in 'static'
    sim <- if (nrow(sim) == N) {
      add_vars(sim, static)
    } else {
      add_vars(sim, collapse::ss(static, sim$ID, check = FALSE))
    }
    rm(static)
    get_vars(sim, "ID") <- NULL

    #---

    ccols <- unlist(lapply(cvars, function(x) paste(x, 1:Mimp, sep = "_")))
    cout <- lapply(ccols, function(v) {
      s <- strsplit(v, "_", fixed = TRUE)[[1]]
      test <- sim %>%
        fgroup_by(c(by, v), sort = FALSE) %>%
        get_vars(wvars) %>%
        fsum(nthreads = cores) %>%
        melt(id.vars = c(by, v), variable.name = "REP", value.name = "EST") %>%
        na_omit(cols = v) %>% # Remove rows where there is no observable outcome (NA) for the analysis variable
        rename(level = !!v) %>%
        fmutate(ANALYSIS = paste("fsum", s[[1]], sep = "."),
                M = as.integer(s[2]),
                REP = as.integer(gsub("^REP__", "", REP)))
    }) %>%
      rbindlist() %>%
      # Complete the data by including zero estimates for unobserved combinations
      # NOTE: The call to complete() is generally slow compared to the fsum() above
      tidyr::complete(tidyr::nesting(!!!rlang::syms(by), M, REP), tidyr::nesting(ANALYSIS, level), fill = list(EST = 0))

  } else {
    cout <- data.frame()
  }

  #---

  # Combine everything and compute final results
  cat("Computing final point estimates and margin-of-error\n")

  # Remove the principle data objects
  suppressWarnings(rm(sim, static))

  # Combine analysis output data frames
  # TO DO -- IMPROVE NAMES!
  result <- rbind(nout, cout)
  rm(nout, cout)

  # 'i0' gives the ANALYSIS values in 'result' that were requested
  # If it wasn't requested, we know it must be a proportion (i.e. mean requested but suym returned in 'cout')
  i0 <- match(unique(result$ANALYSIS), alabel)

  # Compute final estimates
  result <- result %>%
    group_by_at(c(by, "level", "ANALYSIS", "M")) %>%
    summarize(est = EST[REP == 0], # Estimate using the primary weight
              var = (var_scale / sum(REP != 0)) * sum((EST[REP != 0] - est) ^ 2),
              .groups = "drop_last") %>%  # "M" should be the final grouping variable; automatically dropped by summarize()
    summarize(ubar = mean(var), # Mean of the within-implicate (replicate-based) variances
              b = var(est), # Variance of the across-implicate point estimates
              est = mean(est), # Final point estimate; average across the implicates
              .groups = "drop") %>%
    mutate(se = sqrt(ubar + (1 + Mimp^(-1)) * b), # Final standard error (Rubin)
           r = ifelse(ubar == 0, Inf, (1+Mimp^(-1))*b/ubar),  # Set to Inf if ubar = 0 (NA otherwise)
           df = (Mimp-1)*(1+r^(-1))^2,  # Final degrees of freedom
           moe = se * qt(p = 0.95, df),  # Final 90% confidence interval (p = 0.95 means a 90% confidence interval)
           # Calculate share of MOE attributable to replicate weights
           # The replicate-weight-only MOE is equivalent to use just 'ubar' with ; no variance in estimates across implicates (no modeling uncertainty)
           rshare = (sqrt(ubar) * qt(p = 0.95, df = Inf)) / moe) %>%
    # This returns the analysis number (index); possibly multiple, associated with each ANALYSIS
    # Categorical analyses are all "fsum" at this point, and can be matched against requested analyses (sum, mean, or both)
    # If mean or both requested, the unnest() below creates the correct number of entries for each 'i', replicating as necessary
    mutate(i = lapply(ANALYSIS, function(x) unique(na.omit(c(match(x, alabel), match(sub("fsum", "fmean", x), alabel)))))) %>%
    tidyr::unnest(i) %>%
    group_by_at(c(by, "i")) %>%
    # Adjustment factor for rows to be converted to proportions
    # If 'i' is in 'i0' then no conversion to proportion is necessary
    mutate(adj = ifelse(i %in% i0, 1, 1 / sum(est))) %>%
    ungroup() %>%
    mutate(lhs = alhs[i],
           rhs = arhs[i],
           est = est * adj,  # Adjustment for proportions applied here
           moe = moe * adj,
           # Assign the type of analytical result returned
           type = sub("f", "", afun[i]),
           type = ifelse(!is.na(level) & type == "mean",  "prop", type),
           type = ifelse(!is.na(level) & type == "sum",  "count", type)) %>%
    arrange(i, !!!rlang::syms(by), level) %>%  # Arranges rows according to original 'analyses' order
    select(lhs, rhs, type, all_of(by), level, est, moe, rshare) %>%
    mutate_all(tidyr::replace_na, replace = NA) %>%   # Replaces NaN from zero division with normal NA
    mutate_if(is.character, as.factor) %>%
    mutate(level = if (all(is.na(level))) NULL else as.character(level)) # Remove 'level' if it contains no information

  # Report processing time
  tout <- difftime(Sys.time(), t0)
  cat("Total processing time:", signif(as.numeric(tout), 3), attr(tout, "units"), "\n", sep = " ")

  return(result)

}
