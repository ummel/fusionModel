# Lookup details about available variables

lookup <- function(var = NULL, year = NULL) {

  # Loading fusionData dictionary (not necessarily available)
  # data(dictionary)
  # names(dictionary) <- tolower(names(dictionary))

  # Check validity of the working directory path
  # Checks if "/fusionData" is part of the path, as this is required
  b <- strsplit(full.path(getwd()), .Platform$file.sep, fixed = TRUE)[[1]]
  i <- which(b == "fusionData")
  if (length(i) == 0) stop("'/fusionData' is not part of the working directory path; this is required.")
  dir <- paste(b[1:i], collapse = .Platform$file.sep)

  # Get file path(s) to ACS and fused microdata, based on values in 'year'
  fpaths <- lapply(if (is.null(year)) "" else year, function(yr) {
    pa <- list.files(file.path(dir, "survey-processed/ACS"), pattern = paste0(yr, "_._processed.fst"), recursive = TRUE, full.names = TRUE)
    pc <- list.files(file.path(dir, "survey-processed/ACS"), pattern = paste0(yr, "_._custom.fst"), recursive = TRUE, full.names = TRUE)
    pf <- rev(list.files(file.path(dir, "fusion"), pattern = paste0(yr, "_._fused.fsd"), recursive = TRUE, full.names = TRUE))
    c(pa, pc, pf)
  }) %>%
    unlist()

  # Extract the 'var' present in each file in 'fpaths'
  vlist <- if (is.null(var)) {
    lapply(fpaths, function(x) fst::metadata_fst(x)$columnNames)
  } else {
    lapply(fpaths, function(x) intersect(var, fst::metadata_fst(x)$columnNames))
  }

  acs.year <- str_extract(basename(fpaths) , "(\\d{4})(?=_[^_]*_[^_]*$)")
  rtype <- str_extract(basename(fpaths), ".(?=_[^_]*$)")
  rtype <- case_when(rtype == "H" ~ "Household", rtype == "P" ~ "Person")
  source <- basename(fpaths) %>%
    str_extract("^[^_]*_[^_]*") %>%
    str_split(pattern = "_", n = 2, simplify = TRUE) %>%
    as.data.frame() %>%
    setNames(c('survey', 'vintage')) %>%
    mutate(source = 1:n())

  out <- vlist %>%
    tibble::enframe(name = "source", value = "var") %>%
    left_join(source, by = join_by(source)) %>%
    mutate(respondent = rtype,
           acs_year = acs.year,
           path = gsub(dir, "", fpaths, fixed = TRUE),
           source = NULL) %>%
    tidyr::unnest(var) %>%
    filter(!var %in% c('M', 'year', 'hid', 'weight'), !str_detect(var, "^rep_\\d*$")) %>%
    #left_join(dictionary, by = join_by(var == variable, survey, vintage, respondent)) %>%
    #select(var, acs_year, respondent, survey, vintage, description, type, values, path)
    select(var, acs_year, respondent, survey, vintage, path) %>%
    arrange(across(everything()))

  return(out)

}

