% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fusionOutput.R
\name{fusionOutput}
\alias{fusionOutput}
\title{Generate output files resulting from fusion}
\usage{
fusionOutput(
  donor,
  respondent,
  acs_year,
  fusion_vars,
  M = 1,
  fsn = NULL,
  rlocation = NULL,
  note = NULL,
  test_mode = TRUE,
  validation = TRUE,
  ncores = 1,
  margin = 4,
  ...
)
}
\arguments{
\item{M}{Integer. Desired number of fusion implicates. If \code{M = NULL} (default) it is internally set to 40 or, if \code{test_mode = TRUE}, 2 implicates.}

\item{note}{Character. Optional note supplied by user. Inserted in the log file for reference.}

\item{test_mode}{Logical. If \code{test_mode = TRUE} (default), the result files are always saved within a "/fusion_" directory in \code{output} (possibly created); faster hyperparameters are used for \code{\link[fusionModel]{train}}; and the internal validation step is skipped by default.}

\item{validation}{Logical or integer. Controls execution of internal validation (Steps 3 and 4). If \code{validation = 0} or \code{FALSE}, neither step is performed (default when \code{test_mode = TRUE}). If \code{1}, only Step 3. If \code{2} or \code{TRUE}, both Steps 3 and 4.}

\item{ncores}{Integer. Number of physical CPU cores used for parallel computation.}

\item{margin}{Numeric. Passed to same argument in \code{\link[fusionModel]{fuse}}.}

\item{...}{Optional, non-default arguments passed to \code{\link[fusionModel]{train}}. For example, \code{fork = TRUE} to enable forked parallel processing.}

\item{input}{Character. Path to directory containing files created by \code{fusionInput}.}

\item{output}{Character. Optional path to directory where output files will be saved. If \code{output = NULL} (default), the output directory is automatically constructed from \code{input}.}
}
\value{
Saves resulting \code{output} data files to appropriate local directory. Also saves a .txt log file alongside data files that records console output from \code{fusionOutput}.
}
\description{
Handles all operations needed to "do fusion" using input files generated by a successful call to \code{fusionInput}. Trains a fusion model, generates internal validation results, and then simulates multiple implicates for recipient microdata.
}
\details{
The function checks arguments and determines the file path to the appropriate \code{output} directory (creating it if necessary). The output files are always placed within the appropriate directory hierarchy, based on the donor and recipient information detected in the \code{input} file names. In practice, \code{output} need only be specified if working in an environment where the output files need to located somewhere different from the input files.

The function executes the following steps:
\enumerate{
\item \strong{Load training data inputs}. Loads donor training microdata and results of \code{\link[fusionModel]{prepXY}}.
\item \strong{Run fusionModel::train()}. Calls \code{\link[fusionModel]{train}} using sensible defaults and hyperparameters. If \code{test_mode = TRUE}, the hyperparameters are designed to do a fast/rough-and-ready model training.
\item \strong{Fuse onto training data for internal validation}. Optional step (see \code{validation} argument). Fuses multiple implicates to original donor training data using \code{\link[fusionModel]{fuse}}. Results saved to disk.
\item \strong{Run fusionModel::validate()}. Optional step (see \code{validation} argument). Passes previous step's results to \code{\link[fusionModel]{validate}}. Results saved to disk.
\item \strong{Fuse onto prediction data}. Fuses multiple implicates to supplied input prediction data using \code{\link[fusionModel]{fuse}}. Results saved to disk.
\item \strong{fusionOutput() is finished!} Upon completion, a log file named \code{"outputlog.txt"} is written to \code{output} for reference.
}
}
\examples{
# Since 'test_mode = TRUE' by default, this will affect files in local '/fusion_' directory
dir <- fusionInput(donor = "RECS_2015",
                   recipient = "ACS_2015",
                   respondent = "household",
                   fuse = c("btung", "btuel", "cooltype"),
                   force = c("moneypy", "householder_race", "education", "nhsldmem", "kownrent", "recs_division"),
                   note = "Hello world. Reminder: running in test mode by default.")

# List files in the /input directory
list.files(dir)

# Using default settings
out <- fusionOutput(input = dir)
list.files(out)

}
