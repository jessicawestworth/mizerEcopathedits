% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchExtMort.R
\name{matchExtMortOnce}
\alias{matchExtMortOnce}
\title{Match the external mortality to that required by the ecopath production rate.}
\usage{
matchExtMortOnce(params, steady = TRUE)
}
\arguments{
\item{params}{A MizerParams object}

\item{steady}{Whether to return the model to a steady state after adjusting
the external mortality. Default is TRUE.}
}
\value{
A MizerParams object with the external mortality matched
}
\description{
Adjust the external mortality so that the biomass loss due to mortality
(including the fisheries yield) matches the somatic production.
}
\seealso{
Other match functions: 
\code{\link{matchCatch}()},
\code{\link{matchConsumption}()},
\code{\link{matchDiet}()},
\code{\link{matchGonadicProportionOnce}()},
\code{\link{matchProduction}()},
\code{\link{matchRespirationOnce}()}
}
\concept{match functions}
