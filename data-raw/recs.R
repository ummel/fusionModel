recs <- readRDS("~/Documents/Projects/fusionData/survey-processed/RECS/2015/RECS_2015_H_processed.rds") |>
  dplyr::select(
    #weight,
    moneypy,
    hhage,
    householder_race,
    education,
    employhh,
    kownrent,
    nhsldmem,
    division,
    uatyp10,
    iecc_climate_pub,
    typehuq,
    yearmaderange,
    totsqft_en,
    adqinsul,
    tvcolor,
    agerfri1,
    equipm,
    cooltype,
    scalee,
    scaleg,
    kwh,
    cufeetng,
    gallonlp,
    btufo,
    sdescent
  ) |>
  dplyr::rename(
    #sample_weight = weight,
    income = moneypy,
    age = hhage,
    race = householder_race,
    tenure = kownrent,
    employment = employhh,
    hh_size = nhsldmem,
    urban_rural = uatyp10,
    climate = iecc_climate_pub,
    home_type = typehuq,
    year_built = yearmaderange,
    square_feet = totsqft_en,
    insulation = adqinsul,
    televisions = tvcolor,
    refrigerator_age = agerfri1,
    heating = equipm,
    aircon = cooltype,
    disconnect = scalee,
    unhealthy = scaleg,
    electricity = kwh,
    natural_gas = cufeetng,
    propane = gallonlp,
    fuel_oil = btufo
  ) |>
  dplyr::mutate(
    race = ifelse(sdescent == "Yes", "Latino Alone", as.character(race)),
    race = factor(trimws(gsub("Alone", "", race))),
    sdescent = NULL
  ) |>
  as.data.frame()

usethis::use_data(recs, overwrite = TRUE)
