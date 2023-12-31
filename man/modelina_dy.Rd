% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runfunctions.R
\name{modelina_dy}
\alias{modelina_dy}
\title{Create object of class microindaily with weather data as an array}
\usage{
modelina_dy(
  climarray,
  rainarray,
  tme,
  r,
  altcorrect = 0,
  vegp,
  soilc,
  dtm,
  merid = 0,
  dst = 0,
  runchecks = FALSE
)
}
\arguments{
\item{climarray}{a list of arrays of weather variables (see details). See also [nctoarray()]}

\item{rainarray}{an array of daily rainfall (see details)}

\item{tme}{an object of class POSIXlt giving the dates and times for each weather variable stroed in the array}

\item{r}{a raster object giving with the resolution, spatial extent, and projection of the weather data (see details)}

\item{altcorrect}{a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction, see details)}

\item{vegp}{an object of class vegparams as returned by [vegpfromhab()] (see details)}

\item{soilc}{an object of class soilcharac as returned by [soilcfromtype()]}

\item{dtm}{a RasterLayer onject of elevations (see details)}

\item{merid}{optionally, longitude of local time zone meridian (decimal degrees)}

\item{dst}{optionally, numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).}

\item{runchecks}{optional logical indicating whether to call [checkinputs()] to run
checks on format and units of input data.}
}
\description{
The function `modelin` creates an object of class microin
which unpacks various component inputs and reformats as required
for running the model in hourly timesteps. Here it is assumed that the input
weather data are as arrays - i.e. variable in space
}
\details{
The units of `climarray` must follow those in the dataset `climdata`.
It must be a list with each component of the list an array, named using the same
names as the column headers in weather (e.g. temp for temperature), excluding `obs_time`.
Dimensions 1 and 2 of the array must be the same as `r` and dimension 3 must have
the same length as `tme`. If `r` has a different resolution to `dtm` the climate
data are resampled to match the resolution of `dtm`. The array of Plant Area index values in `vegp` must
of the same x and y dims as `dtm` but can contain any number of repeated
measures up to the number of entries in `tme`. Data are interpolated to the
time increment of `tme`. Other vegetation paramaters, including vegetation
height are assumed time-invarient. The RasterLayer datasets in `soilc` must have
the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
metres and the coordinate reference system must be defined. If `altcorrect`>0,
and the dimenions of `r` are not identical to those of `dtm`, the elevation
difference between each pixel of the dtm and the dtm coarsed to the resolution of
`r` is calaculated and an elevational lapse rate correction is applied to the
temperature data to accoutn for these elevation differences. If `altcorrect`=1,
a fixed lapse rate of 5 degrees per 100m is applied. If `altcorrect`=2, humidity-dependent
lapse rates are calaculate and applied.
}
\seealso{
[modelin_dy()], [modelina()], [nctoarray()]
}
