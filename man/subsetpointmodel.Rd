% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataprep.R
\name{subsetpointmodel}
\alias{subsetpointmodel}
\title{Subsets outputs from point microclimate model}
\usage{
subsetpointmodel(pointmodel, tstep = "month", what = "tmax", days = NA)
}
\arguments{
\item{pointmodel}{a list of model outputs from the point microclimate model as returned by [runpointmodel()].}

\item{tstep}{one of `year` or `month` (see details)}

\item{what}{one of `tmax`, `tmin` or `tmedian` (maximum, minimum or median temperature respectively - see details)}

\item{days}{optionally a vector of the days in the time sequence  to return data for (if provided `tstep` is ignored)}
}
\value{
A list with the same format as `pointmodel` but with specified values
'only selected.
}
\description{
The function `subsetpointmodel` provides a means of selecting monthly or
yearly values from the ouputs of the point microclimate model
}
\details{
setting 'tstep' to `year` identifies the day in each year with the e.g.
the hottest or coldest hourly temperature, and 'tstep' to `month` ideas the day
in each month in each year with e.g. the hottest or coldest hourly temperature.
Values for all hours of that day are returned, to ensure that the ground heat flux
in the grid microclimate model can be estimated. If `what` is set to `tmax` or `tmin`
the hottest or coldest hour within each month or year are identified. if If `what`
is set to `tmedian` hourly temperatures within the month or year are ranked and
the median hour identified.
}
\examples{
# Extract all hourly values for day in which hottest hour in each month occurs
sub_micropoint_hr <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")
Tcanopy <- sub_micropoint_hr$microp$Tc
tme <- as.POSIXct(sub_micropoint_hr$weather$obs_time, tz = "UTC")
# Plot
plot(Tcanopy ~ tme, type = "l")
}
\seealso{
[runpointmodel()]
}
