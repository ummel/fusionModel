# Determine 'yord' of models
fusionOrder <- function(varimp) {

  yvars <- names(varimp)
  #yvars <- map_chr(models, "yvar")

  var.imp <- map2_dfr(.x = varimp, .y = yvars, .f = getVarImp, yvars = yvars)
  #var.imp <- map_dfr(models, getVarImp, yvars = yvars)

  # How important are the xvars for each y?
  # When the xvars have high collective importance, y can be moved forward...
  ximp <- var.imp %>%
    filter(predictor == "_xvars_") %>%
    select(response, importance) %>%
    rename(ximp = importance,
           yvar = response)

  yord <- vector(mode = "character", length = length(yvars))
  for (i in 1:(length(yord) - 1)) {

    yimp <- var.imp %>%
      filter(predictor != "_xvars_") %>%
      group_by(predictor) %>%
      summarize(yimp = mean(importance)) %>%
      rename(yvar = predictor)

    ranking <- left_join(yimp, ximp, by = "yvar") %>%
      arrange(-(ximp + yimp))

    vmax <- ranking$yvar[1]
    yord[i] <- vmax
    var.imp <- filter(var.imp, predictor != vmax, response != vmax)

  }

  # Add the final imputation variable
  yord[length(yord)] <- setdiff(yvars, yord)

  return(yord)

}

#----------

getVarImp <- function(vimp, y, yvars) {
  #m$variable.importance %>%
  vimp %>%
    tibble::enframe(name = "predictor", value = "importance") %>%
    mutate(predictor = ifelse(predictor %in% yvars, predictor, "_xvars_"),
           predictor = factor(predictor, levels = c(yvars, "_xvars_"))) %>%
    group_by(predictor) %>%
    summarize(importance = sum(importance), .groups = 'drop') %>%
    complete(predictor, fill = list(importance = 0)) %>%
    mutate(predictor = as.character(predictor),
           response = y) %>%
    filter(predictor != response)
}
