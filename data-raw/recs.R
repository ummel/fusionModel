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
    fuelheat,
    adqinsul,
    cooltype,
    agecenac,
    tvcolor,
    scalee,
    totsqft_en,
    kwh,
    cufeetng,
    gallonfo,
    gallonlp,
    btulp,
    dollarlp,
    sdescent,
    totalbtu,
    totalbtusph,
    btuelahucol,
    btuelcol,
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
    heat_type = fuelheat,
    insulation = adqinsul,
    aircon = cooltype,
    centralac_age = agecenac,
    televisions = tvcolor,
    disconnect = scalee,
    square_feet = totsqft_en,
    electricity = kwh,
    natural_gas = cufeetng,
    fuel_oil = gallonfo,
    propane = gallonlp,
    propane_btu = btulp,
    propane_expend = dollarlp
  ) |>
  dplyr::mutate(
    electricity = as.integer(round(electricity)),  # Make integer for testing purposes
    heating_share = round(totalbtusph / totalbtu, 3),
    cooling_share = round((btuelahucol + btuelcol) / totalbtu, 3),
    other_share = 1 - heating_share - cooling_share,
    urban_rural = factor(ifelse(urban_rural == "U", "Urban", "Rural")),
    renter = renter != "Owned or being bought by someone in your household",
    use_ng = natural_gas > 0,
    have_ac = aircon != "No air conditioning",
    race = ifelse(sdescent == "Yes", "Latino Alone", as.character(race)),
    race = factor(trimws(gsub("Alone", "", race)))
  ) |>
  select(-sdescent, -totalbtu, -totalbtusph, -btuelahucol, -btuelcol) |>
  select(-starts_with("rep_"), starts_with("rep_"))

usethis::use_data(recs, overwrite = TRUE)
