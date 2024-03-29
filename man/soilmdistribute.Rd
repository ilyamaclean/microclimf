% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/componentfunctions.R
\name{soilmdistribute}
\alias{soilmdistribute}
\title{Spatially distribute soil moisture by topographic wetness index}
\usage{
soilmdistribute(micro, tfact = 1.5, twi = NA)
}
\arguments{
\item{micro}{an object of class microin as returned by e.g. modelin}

\item{tfact}{coefficient determining sensitivity of soil moisture to variation
in topographic wetness}

\item{twi}{optional SpatRaster object of topographic wetness index values.
Calculated from he dtm in micro of not supplied, but outer cells will be NA.}
}
\value{
a 3D array of soil moistures with the same x and y dims as `dtm` and z
equivelent to length(soilm)
}
\description{
The function `soilmdistribute` spatially distrubutes soil moisture by the Bevan and Kirkby
topographic wetness index
}
