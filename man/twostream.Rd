% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/componentfunctions.R
\name{twostream}
\alias{twostream}
\title{Applies two-stream radiation model}
\usage{
twostream(
  micro,
  reqhgt = 0.05,
  pai_a = NA,
  tfact = 1.5,
  slr = NA,
  apr = NA,
  hor = NA,
  twi = NA
)
}
\arguments{
\item{micro}{an object of class `micro` as returned by [modelin()]}

\item{reqhgt}{height at which temperatures are required (m). Negative if below ground}

\item{pai_a}{plant area index above `reqhgt`. Determined from total `pai` if not supplied}

\item{tfact}{coefficient determining sensitivity of soil moisture to variation
in topographic wetness}

\item{slr}{an optional SpatRaster object of slope values (Radians). Calculated from
the dtm in micro if not supplied, but outer cells will be NA.}

\item{apr}{an optional SpatRaster object of aspect values (Radians). Calculated from
the dtm in micro if not supplied, but outer cells will be NA.}

\item{hor}{an optional array of the tangent of the angle to the horizon in
24 directions. Calculated from the dtm in micro if not supplied, but outer cells will be NA.}

\item{twi}{optional SpatRaster object of topographic wetness index values.
Calculated from he dtm in micro if not supplied, but outer cells will be NA.}
}
\value{
an object of class micro with additional terms added for subsequent modelling
}
\description{
The function `twostream` applies a varient of the Sellers two-stream
radiation model to calaculate long and showtwave radiation absorbed by the ground and
canopy and shortwave radiation absorbed by leafs.
}
