% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rimplementation.R
\name{modelin}
\alias{modelin}
\title{Create object of class microin}
\usage{
modelin(
  micropoint,
  vegp,
  soilc,
  dtm,
  dtmc = NA,
  altcorrect = 0,
  runchecks = TRUE
)
}
\arguments{
\item{micropoint}{either an object of class `micropoint` as returned
by [runpointmodel()] or [subsetpointmodel()] or a list of micropoint
objects as returned by [runpointmodela()] or [subsetpointmodela()],
which assumes input climate data are spatially variable.}

\item{vegp}{an object of class vegparams as returned by [vegpfromhab()] (see details)}

\item{soilc}{an object of class soilcharac as returned by [soilcfromtype()]}

\item{dtm}{a SpatRaster object of elevations in metres (see details)}

\item{dtmc}{a SpatRaster object giving the resolution, spatial extent, and projection
of the weather data used in `micropoint`. Ignored if climate data provided
to `micropoint` as a data.frame. Must give elevations in metres if
climate data were arrays and `altcorrect` > 0 or if setting runchecks to TRUE.}

\item{altcorrect}{a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction, see details)}

\item{runchecks}{optional logical indicating whether to call [checkinputs()] to run}
}
\description{
The function `modelin` creates an object of class microin
which unpacks various component inputs and reformats as required
for running the R version of the model.
}
\details{
The array of Plant Area index values in `vegp` must
of the same x and y dims as `dtm` but can contain any number of repeated
measures up to the number of entries in `weather`. Data are interpolated to the
time increment of `weather`. Other vegetation paramaters, including vegetation
height are assumed time-invarient. If these are time-varient, use the c++ version
of the model. The SpatRaster datasets in `soilc` must have
the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
metres and the coordinate reference system must be defined.
If not NA, the units of dtmc must match dtm and must be elevation data if an
altitude correction is applied. If `altcorrect`>0, the elevation
difference between each pixel of dtm and dtmc is calculated and
an elevation lapse rate correction is applied to the temperature and
pressure data to account for these elevation differences. If `altcorrect`= 1, a fixed
lapse rate of 5 degrees per 100m is applied to the temperature data. If
`altcorrect`= 2, humidity-dependent lapse rates are calculated and applied.
See also details for  modelin.
}
\examples{
# ================================================================= #
# Preparing model inputs using climate data provided as a data.frame
# ================================================================= #
# run and subset point model
micropoint<-runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc)
micropoint <- subsetpointmodel(micropoint)
micro <- modelin(micropoint, vegp, soilc, dtmcaerth)
# ================================================================= #
# Preparing model inputs using climate data provided as arrays
# ================================================================= #
# Create dummy array data
.ta<-function(x,dtm,xdim=5,ydim=5) {
  a<-array(rep(x,each=ydim*xdim),dim=c(ydim,xdim,length(x)))
  .rast(a,dtm)
}
dtm <- rast(dtmcaerth)
climarrayr<-list(temp = .ta(climdata$temp, dtm),
  relhum = .ta(climdata$relhum, dtm),
  pres = .ta(climdata$pres, dtm),
  swdown = .ta(climdata$swdown, dtm),
  difrad = .ta(climdata$difrad, dtm),
  lwdown = .ta(climdata$lwdown, dtm),
  windspeed = .ta(climdata$windspeed, dtm),
  winddir = .ta(climdata$winddir, dtm),
  precip = .ta(climdata$precip, dtm))
tme <- as.POSIXlt(climdata$obs_time, tz="UTC")
# Run and subset point model array
micropointa <- runpointmodela(climarrayr, tme, reqhgt = 0.05, dtmcaerth, vegp, soilc)
micropointa <- subsetpointmodela(micropointa)
# Prepare grid model input
dtmc <- resample(dtm, climarrayr$temp[[1]])
micro <- modelin(micropointa, vegp, soilc, dtm, dtmc)
}
