% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonhourly.R
\name{runmicro_dy}
\alias{runmicro_dy}
\title{Run microclimate model (daily)}
\usage{
runmicro_dy(
  micro_dy,
  reqhgt,
  expand = TRUE,
  pai_a = NA,
  tfact = 1.5,
  surfwet = NA,
  slr = NA,
  apr = NA,
  hor = NA,
  twi = NA,
  wsa = NA
)
}
\arguments{
\item{reqhgt}{height above ground at which model outputs are needed (m).}

\item{expand}{optional logical indicating whether to expand daily values to hourly (see details).}

\item{pai_a}{an optional array of plant area index values above `reqhgt` (see details)}

\item{tfact}{coefficient determining sensitivity of soil moisture to variation
in topographic wetness (see [soilmdistribute()])}

\item{surfwet}{an optional single numeric value of array of values specifying the proportion
of the canopy surface that should be treated as wet surface (modelled if not supplied)}

\item{slr}{an optional SpatRaster object of slope values (Radians). Calculated from
dtm if not supplied, but outer cells will be NA.}

\item{apr}{an optional SpatRaster object of aspect values (Radians). Calculated from
dtm if not supplied, but outer cells will be NA.}

\item{hor}{an optional array of the tangent of the angle to the horizon in
24 directions. Calculated from dtm if not supplied, but outer cells will be NA.}

\item{twi}{optional SpatRaster object of topographic wetness index values.
Calculated from `dtm` of not supplied, but outer cells will be NA.}

\item{wsa}{an optional array of wind shelter coefficients in 8 directions.
Calculated from dtm if not supplied, but outer cells will be NA.}

\item{microd}{object of class microindaily as returned by [modelin_dy()]}
}
\value{
if expand = TRUE, an object of class microout with the following components:

`Tz` Array of air temperatures at height `reqhgt` (deg C). Identical to `T0`
if `reqhgt` = 0.

`tleaf` Array of leaf temperatures at height `reqhgt` (deg C).
NA if `reqhgt` greater than canopy height or `reqhgt` <= 0.

`T0` Array of ground surface temperatures (deg C)

`relhum` Array of relative humidities at height `reqhgt` (percentage).
NA if `reqhgt` <= 0.

`windspeed` Array of wind speeds at height `reqhgt` (m/s).
NA if `reqhgt` <= 0.

`Rdirdown` Array of downward direct shortwave radiation incident on
horizontal surface (W/m^2)

`Rdifdown` Array of downward diffuse shortwave radiation incident on
horizontal surface (W/m^2)

`Rlwdown` Array of downward longwave radiation incident on horizontal
surface (W/m^2)

`Rswup` Array of upward shortwave radiation (assumed diffuse) incident
on underside of horizontal surface (W/m^2)

`Rlwup` Array of upward longwave radiation incident on underside of
horizontal surface (W/m^2)

if expand = FALSE, an object of class microutdaily, list with the following components:

mout_mn an object of class microut for minimum daily temperatures

mout_mx an object of class microut for maximum daily temperatures
}
\description{
The function `runmicro_dy` runs the microclimate model in daily
time increments, with the option to expand to hourly
}
\details{
If expand = TRUE, daily minima and maxima are expanded to hourly using values in the hourly
weather dataset. See also details for [runmicro_hr()].
}
\seealso{
[runmicro_hr()] for running microclimate model in hourly time-steps and
[runmicro_big()] for running the microclimate model over large areas as tiles.
}
