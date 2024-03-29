% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataprep.R
\name{leafrfromalb}
\alias{leafrfromalb}
\title{Derives leaf reflectance from albedo}
\usage{
leafrfromalb(pai, gref, x, alb, ltrr = 1, out = "lref")
}
\arguments{
\item{pai}{a SpatRaster of plant area index values}

\item{gref}{a SpatRaster of ground reflectance values}

\item{x}{a SpatRaster of the ratio of vertical to horizontal projections of leaf foliage}

\item{alb}{a SpatRaster of white-sky albedo}

\item{ltrr}{an optional numeric value giving an approximate estimate of the ratio of leaf transmittance to leaf reflectance (e.g. value of 1 makes leaf transmittance equal to reflectance). See details}
}
\value{
leaf reflectance in the range 0 to 1,
}
\description{
Derives leaf reflectance from albedo
}
\details{
the microclimate model is not unduly sensitive to `lttr` so if unknown, an apprxoimate
value or the default can be used.
}
\examples{
pai <- .rast(vegp$pai[,,9],rast(dtmcaerth)) # Plant Area Index in Sep (month in whihc albedo image was flow)
gref <- rast(soilc$groundr)
x <- rast(vegp$x)
alb <- rast(albedo)
leafr <- leafrfromalb(pai, gref, x, alb)
plot(leafr)
}
