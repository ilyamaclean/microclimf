% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runfunctions.R
\name{modelin_dy}
\alias{modelin_dy}
\title{Create object of class microindaily}
\usage{
modelin_dy(
  weather,
  rainfall,
  vegp,
  soilc,
  dtm,
  windhgt = 2,
  merid = 0,
  dst = 0,
  runchecks = TRUE
)
}
\arguments{
\item{weather}{a data.frame of hourly weather variables (see details)}

\item{rainfall}{a vector of daily rainfall}

\item{vegp}{an object of class vegparams as returned by [vegpfromhab()] (see details)}

\item{soilc}{an object of class soilcharac as returned by [soilcfromtype()]}

\item{dtm}{a RasterLayer onject of elevations (see details)}

\item{merid}{optionally, longitude of local time zone meridian (decimal degrees)}

\item{dst}{optionally, numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).}
}
\description{
The function `modelin` creates an object of class microindaily
which unpacks various component inputs and reformats as required
for running the model in daily timesteps
}
\details{
The format and and units of `weather` must follow that in the example
dataset `climdata`. The array of Plant Area index values in `vegp` must
of the same x and y dims as `dtm` but can contain any number of repeated
measures up to the number of entries in `weather`. Data are interpolated to the
time increment of `weather`. Other vegetation paramaters, including vegetation
height are assumed time-invarient. The RasterLayer datasets in `soilc` must have
the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
metres and the coordinate reference system must be defined.
}
\seealso{
[inputchecks()], [modelin()]
}
