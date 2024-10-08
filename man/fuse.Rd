% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuse.R
\name{fuse}
\alias{fuse}
\title{Fuse variables to a recipient dataset}
\usage{
fuse(
  data,
  fsn,
  fsd = NULL,
  M = 1,
  retain = NULL,
  kblock = 10,
  margin = 2,
  cores = 1
)
}
\arguments{
\item{data}{Data frame. Recipient dataset. All categorical variables should be factors and ordered whenever possible. Data types and levels are strictly validated against predictor variables defined in \code{fsn}.}

\item{fsn}{Character. Path to fusion model file (.fsn) generated by \code{\link{train}}.}

\item{fsd}{Character. Optional fusion output file to be created ending in \code{.fsd} (i.e. "fused data"). This is a compressed binary file that can be read using the \code{\link[fst]{fst}} package. If \code{fsd = NULL} (the default), the fusion results are returned as a \code{\link[data.table]{data.table}}.}

\item{M}{Integer. Number of implicates to simulate.}

\item{retain}{Character. Names of columns in \code{data} that should be retained in the output; i.e. repeated across implicates. Useful for retaining ID or weight variables for use in subsequent analysis of fusion output.}

\item{kblock}{Integer. Fixed number of nearest neighbors to use when fusing variables in a block. Must be >= 5 and <= 30. Not applicable for variables fused on their own (i.e. no block).}

\item{margin}{Numeric. Safety margin used when estimating how many implicates can be processed in memory at once. Set higher if \code{fuse()} experiences a memory shortfall. Alternatively, can be set to a negative value to manually specify the number of chunks to use. For example, \code{margin = -3} splits \code{M} implicates into three chunks of approximately equal size.}

\item{cores}{Integer. Number of cores used. LightGBM prediction is parallel-enabled on all systems if OpenMP is available.}
}
\value{
If \code{fsd = NULL}, a \code{\link[data.table]{data.table}} with number of rows equal to \code{M * nrow(data)}. Integer column "M" indicates implicate assignment of each observation. Note that the ordering of recipient observations is consistent within implicates, so do not change the row order if using with \code{\link{analyze}}.

If \code{fsd} is specified, the path to .fsd file where results were written. Metadata for column classes and factor levels are stored in the column names. \code{\link{read_fsd}} should be used to load files saved via the \code{fsd} argument.
}
\description{
Fuse variables to a recipient dataset using a .fsn model produced by \code{\link{train}}. Output can be passed to \code{\link{analyze}} and \code{\link{validate}}.
}
\details{
TO UPDATE.
}
\examples{
# Build a fusion model using RECS microdata
# Note that "fusion_model.fsn" will be written to working directory
?recs
fusion.vars <- c("electricity", "natural_gas", "aircon")
predictor.vars <- names(recs)[2:12]
fsn.path <- train(data = recs, y = fusion.vars, x = predictor.vars)

# Generate single implicate of synthetic 'fusion.vars',
#  using original RECS data as the recipient
recipient <- recs[predictor.vars]
sim <- fuse(data = recipient, fsn = fsn.path)
head(sim)

# Calling fuse() again produces different results
sim <- fuse(data = recipient, fsn = fsn.path)
head(sim)

# Generate multiple implicates
sim <- fuse(data = recipient, fsn = fsn.path, M = 5)
head(sim)
table(sim$M)

# Optionally, write results directly to disk
# Note that "results.fsd" will be written to working directory
sim <- fuse(data = recipient, fsn = fsn.path, M = 5, fsd = "results.fsd")
sim <- read_fsd(sim)
head(sim)

}
