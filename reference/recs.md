# Data from the 2015 Residential Energy Consumption Survey (RECS)

Pre-processed, household-level microdata containing a selection of 31
variables derived from the 2015 RECS, plus 96 survey replicate weights.
A variety of data types are included (numeric, factor, logical). There
are no missing values (`NA`). Variable names have been standardized from
the original EIA microdata format.

## Usage

``` r
recs
```

## Format

A tibble (`tbl_df`) with 5,686 rows and 124 variables:

- weight:

  Primary sampling weight.

- income:

  Annual gross household income for the last year.

- age:

  Respondent age.

- race:

  Respondent race.

- education:

  Highest education completed by respondent.

- employment:

  Respondent employment status.

- hh_size:

  Number of household members.

- division:

  Census Division.

- urban_rural:

  Census 2010 Urban Type.

- climate:

  IECC Climate Code.

- renter:

  Is household renting the home?

- home_type:

  Type of housing unit.

- year_built:

  Range when housing unit was built.

- square_feet:

  Total square footage.

- insulation:

  Level of insulation.

- heating:

  Main space heating fuel.

- aircon:

  Type of air conditioning equipment used.

- centralac_age:

  Age of central air conditioner.

- televisions:

  Number of televisions used.

- disconnect:

  Frequency of receiving disconnect notice.

- electricity:

  Total annual electricity usage, in kilowatthours.

- natural_gas:

  Total annual natural gas usage, in hundred cubic feet.

- fuel_oil:

  Total annual fuel oil/kerosene usage, in gallons.

- propane:

  Total annual propane usage, in gallons.

- propane_btu:

  Total annual propane usage, in thousand Btu.

- propane_expend:

  Total annual propane expenditure, in dollars.

- heating_share:

  Share of total energy used for space heating.

- cooling_share:

  Share of total energy used for cooling (AC and fans).

- other_share:

  Share of total energy used for other end-uses.

- use_ng:

  Logical indicating if household uses natural gas.

- have_ac:

  Logical indicating if household has air conditioning.

- rep_1, rep_2, ..., rep_96:

  Replicate weights (rep_1 through rep_96) for uncertainty estimation.

## References

U.S. Energy Information Administration (EIA). (2015). *Residential
Energy Consumption Survey (RECS)*.
<https://www.eia.gov/consumption/residential/data/2015/>

## Examples

``` r
data(recs)
head(recs)
```
