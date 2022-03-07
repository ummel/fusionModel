#' Data from the 2015 Residential Energy Consumption Survey (RECS)
#'
#' Pre-processed, household-level microdata containing a selection of 28
#' variables derived from the 2015 RECS. A variety of data types are included.
#' There are no missing values. Variable names have been altered from the original.
#'
#' @format A tibble with 5,686 rows and 124 variables:
#' \describe{
#'   \item{weight}{Primary sampling weight}
#'   \item{income}{Annual gross household income for the last year}
#'   \item{age}{Respondent age}
#'   \item{race}{Respondent race}
#'   \item{education}{Highest education completed by respondent}
#'   \item{employment}{Respondent employment status}
#'   \item{hh_size}{Number of household members}
#'   \item{division}{Census Division}
#'   \item{urban_rural}{Census 2010 Urban Type}
#'   \item{climate}{IECC Climate Code}
#'   \item{renter}{Is household renting the home?}
#'   \item{home_type}{Type of housing unit}
#'   \item{year_built}{Range when housing unit was built}
#'   \item{square_feet}{Total square footage}
#'   \item{insulation}{Level of insulation}
#'   \item{heating}{Main space heating fuel}
#'   \item{aircon}{Type of air conditioning equipment used}
#'   \item{centralac_age}{Age of central air conditioner}
#'   \item{televisions}{Number of televisions used}
#'   \item{disconnect}{Frequency of receiving disconnect notice}
#'   \item{electricity}{Total annual electricity usage, in kilowatthours}
#'   \item{natural_gas}{Total annual natural gas usage, in hundred cubic feet}
#'   \item{fuel_oil}{Total annual fuel oil/kerosene usage, in gallons}
#'   \item{propane}{Total annual propane usage, in gallons}
#'   \item{propane_btu}{Total annual propane usage, in thousand Btu}
#'   \item{propane_expend}{Total annual propane expenditure, in dollars}
#'   \item{use_ng}{Logical indicating if household uses natural gas}
#'   \item{have_ac}{Logical indicating if household has air conditioning}
#'   \item{rep_1:rep_96}{Replicate weights for uncertainty estimation}
#'   }
#' @source \url{https://www.eia.gov/consumption/residential/data/2015/}
"recs"
