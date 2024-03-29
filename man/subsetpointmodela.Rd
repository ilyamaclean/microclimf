% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataprep.R
\name{subsetpointmodela}
\alias{subsetpointmodela}
\title{Subsets outputs from point microclimate model run over arrays}
\usage{
subsetpointmodela(pointmodela, tstep = "month", what = "tmax", days = NA)
}
\arguments{
\item{pointmodela}{a list of model outputs from the point microclimate model applied over arrays as returned by [runpointmodela()].}

\item{tstep}{one of `year` or `month` (see details)}

\item{what}{one of `tmax`, `tmin` or `tmedian` (maximum, minimum or median temperature repsectively - see details)}

\item{days}{optionally a vector of the days in the time sequence  to return data for (if provided tstep is ignored)}
}
\value{
A list with the same format as `pointmodela` but with specified values
'only selected.
}
\description{
The function `subsetpointmodela` is the equivalent of [subsetpointmodel()] used
for outputs from [`runpointmodela()]
}
\seealso{
[subsetpointmodel()], [runpointmodela()]
}
