recs <- fst::read_fst("~/Documents/Projects/fusionData/survey-processed/RECS/2015/RECS_2015_H_processed.fst") |>
  dplyr::select(
    weight,
    moneypy,
    hhage,
    householder_race,
    education,
    employhh,
    nhsldmem,
    recs_division,
    ur12,
    recs_iecc_zone,
    kownrent,
    typehuq,
    yearmaderange,
    totsqft_en,
    adqinsul,
    fuelheat,
    cooltype,
    agecenac,
    tvcolor,
    scalee,
    kwh,
    cufeetng,
    gallonfo,
    gallonlp,
    btulp,
    dollarlp,
    sdescent,
    rep_1:rep_96
  ) |>
  dplyr::rename(
    income = moneypy,
    age = hhage,
    race = householder_race,
    employment = employhh,
    hh_size = nhsldmem,
    division = recs_division,
    urban_rural = ur12,
    climate = recs_iecc_zone,
    renter = kownrent,
    home_type = typehuq,
    year_built = yearmaderange,
    square_feet = totsqft_en,
    insulation = adqinsul,
    heating = fuelheat,
    aircon = cooltype,
    centralac_age = agecenac,
    televisions = tvcolor,
    disconnect = scalee,
    electricity = kwh,
    natural_gas = cufeetng,
    fuel_oil = gallonfo,
    propane = gallonlp,
    propane_btu = btulp,
    propane_expend = dollarlp
  ) |>
  dplyr::mutate(
    urban_rural = factor(ifelse(urban_rural == "U", "Urban", "Rural")),
    renter = renter != "Owned or being bought by someone in your household",
    use_ng = natural_gas > 0,
    have_ac = aircon != "No air conditioning",
    race = ifelse(sdescent == "Yes", "Latino Alone", as.character(race)),
    race = factor(trimws(gsub("Alone", "", race))),
    sdescent = NULL
  ) |>
  select(-starts_with("rep_"), starts_with("rep_"))
  #as.data.frame()

usethis::use_data(recs, overwrite = TRUE)
