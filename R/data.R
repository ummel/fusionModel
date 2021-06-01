#' Data from the 2015 Residential Energy Consumption Survey (RECS)
#'
#' Pre-processed, household-level microdata containing a selection of 24
#' variables derived from the 2015 RECS. A variety of data types are included.
#' There are no missing values.
#'
#' @format A data frame with 5,686 rows and 24 variables:
#' \describe{
#'   \item{income}{Annual gross household income for the last year}
#'   \item{age}{Respondent age}
#'   \item{race}{Respondent race}
#'   \item{education}{Highest education completed by respondent}
#'   \item{employment}{Respondent employment status}
#'   \item{tenure}{Does respondent own or rent?}
#'   \item{hh_size}{Number of household members}
#'   \item{division}{Census Division}
#'   \item{urban_rural}{Census 2010 Urban Type}
#'   \item{climate}{IECC Climate Code}
#'   \item{home_type}{Type of housing unit}
#'   \item{year_built}{Range when housing unit was built}
#'   \item{square_feet}{Total square footage}
#'   \item{insulation}{Level of insulation}
#'   \item{televisions}{Number of televisions used}
#'   \item{refrigerator_age}{Age of most-used refrigerator}
#'   \item{heating}{Main space heating equipment type}
#'   \item{aircon}{Type of air conditioning equipment used}
#'   \item{disconnect}{Frequency of receiving disconnect notice}
#'   \item{unhealthy}{Frequency of keeping home at unhealthy temperature}
#'   \item{electricity}{Total site electricity usage, in kilowatthours}
#'   \item{natural_gas}{Total natural gas usage, in hundred cubic feet}
#'   \item{propane}{Total propane usage, in gallons}
#'   \item{fuel_oil}{Total fuel oil/kerosene usage, in thousand Btu}
#'   }
#' @source \url{https://www.eia.gov/consumption/residential/data/2015/}
"recs"
