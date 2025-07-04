% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fillDefaultsFromFishBase.R
\name{fillDefaultsFromFishBase}
\alias{fillDefaultsFromFishBase}
\title{Fill Species Traits Using FishBase}
\usage{
fillDefaultsFromFishBase(
  species_params,
  scientific_name_col = "Scientific_name",
  overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{species_params}{A data frame containing at least a column with scientific names for each species. Other
columns (e.g. common names or user-defined identifiers) may be included and will be preserved. This is
typically a minimal table created at the start of model construction (e.g., from a species list), not yet a
full \code{species_params} object.}

\item{scientific_name_col}{A string giving the column name in \code{species_params} that holds the scientific names.
Defaults to \code{"Scientific_name"}.}

\item{overwrite}{Logical (default = \code{FALSE}). If \code{TRUE}, existing values in the data frame are replaced by FishBase values.
If \code{FALSE}, only missing values are filled.}

\item{verbose}{Logical (default = \code{TRUE}). If \code{TRUE}, prints a summary of filled species.}
}
\value{
The same data frame, augmented with the following columns where data are available:
\itemize{
\item \code{a}, \code{b}: Length–weight allometric coefficients
\item \code{Length}: Maximum length (cm)
\item \code{w_max}: Maximum weight (calculated from Length, a, and b)
\item \code{l_mat}: Length at maturity (median)
\item \code{age_mat}: Age at maturity (median)
\item \code{w_mat}: Weight at maturity (from l_mat, a, and b)
}
}
\description{
Populates key life-history traits in a data frame of species using FishBase, via the \code{rfishbase} package.
}
\details{
This function is typically used early in model setup, to populate a table of species names with biological parameters needed for Mizer models. It supports cases where only species names and scientific names are initially known.

This is a convenience function that queries FishBase through the \code{rfishbase} package.
Traits are pulled from the \code{estimate()}, \code{maturity()}, and \code{species()} tables. Maturity traits
are summarised by median within FishBase \code{SpecCode} groups. Derived weights are calculated from
length–weight relationships.

Existing columns in the input will only be overwritten if \code{overwrite = TRUE}.
}
\note{
Requires the \code{rfishbase} package. Install it with \code{install.packages("rfishbase")}.
}
\examples{
\dontrun{
# Minimal example with only scientific names
species_df <- data.frame(
  Scientific_name = c("Merluccius merluccius", "Scomber scombrus")
)
enriched <- fillDefaultsFromFishBase(species_df)

# Optional: include user-defined species labels
species_df <- data.frame(
  species = c("Hake", "Mackerel"),
  Scientific_name = c("Merluccius merluccius", "Scomber scombrus")
)
enriched <- fillDefaultsFromFishBase(species_df)
}

}
\seealso{
\code{\link[rfishbase:estimate]{rfishbase::estimate()}}, \code{\link[rfishbase:maturity]{rfishbase::maturity()}}, \code{\link[rfishbase:species]{rfishbase::species()}}
}
