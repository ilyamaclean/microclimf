#' A one m resolution dataset of white sky albedos
#'
#' A spatial dataset of white sky albedos for Caerthillean Cove: an area bounded by 169475, 169525, 12475, 12525
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system (CRS: 27700).
#'
#' @format A PackedSpatRaster object with 50 rows and 50 columns
#' @source derived using aerial imagery from \url{https://digimap.edina.ac.uk/}
"albedo"
#' A data frame of hourly weather
#'
#' A data frame of hourly weather in 2017 at Caerthillian Cove, Lizard, Cornwall (49.96807N, 5.215668W)
#'
#' @format a data frame with the following elements:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times}
#'  \item{temp}{temperature (degrees C)}
#'  \item{relhum}{relative humidity (percentage)}
#'  \item{pres}{atmospheric press (kPa)}
#'  \item{swdown}{Total downward shortwave radiation (W / m^2)}
#'  \item{difrad}{Diffuse radiation (W / m^2)}
#'  \item{lwdown}{Total downward longwave radiation (W / m^2)}
#'  \item{windspeed}{Wind speed (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#'  \item{precip}{Precipitation (mm)}
#' }
"climdata"
#'
#' A one m resolution digital terrain dataset.
#'
#' A spatial dataset of elevations (m) for Caerthillean Cove: an area bounded by 169475, 169525, 12475, 12525
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system (CRS: 27700).
#'
#' @format A PackedSpatRaster object with 50 rows and 50 columns
#' @source \url{http://www.tellusgb.ac.uk/}
"dtmcaerth"
#' A dataset of global climate variables
#'
#' A global dataset containing containing the following climate variables averaged
#' over the period 2008 to 2017. Used by [vegpfromhab()] to generate plant arrea index values.
#'
#' @format An array with 94 rows, 192 columns and the following five climate variables:
#' \describe{
#'   \item{1}{mean annual temperature (ºC)}
#'   \item{2}{coefficient of variation in temperature (K)}
#'   \item{3}{mean annual temperature (mm per year)}
#'   \item{4}{coefficient of variation in annual rainfall (mm per 0.25 days)}
#'   \item{5}{numeric month with the most rainfall (1-12)}
#' }
#' @source \url{http://www.ncep.noaa.gov/}
"globclim"
#' A one m resolution dataset of habitat types
#'
#' A spatial dataset of habitat types expressed as integers for Caerthillean Cove: an area bounded by 169475, 169525, 12475, 12525
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system (CRS: 27700). Integer values
#' correspond to habitat types listed in [vegpfromhab()].
#'
#' @format A PackedSpatRaster object with 50 rows and 50 columns
"habitats"
#' A dataset of soil characteristics.
#'
#' An object of class soilcharac for the area bounded by 169475, 169525, 12475, 12525
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system (CRS: 27700).
#'
#' @format A list of  the following objects:
#' \describe{
#'   \item{soiltype}{a PackedSpatRaster object of numeric integers corresponding to soil types given in soilparameters}
#'   \item{groundr}{a PackedSpatRaster object of soil reflectance of shortwave radiation}
#'}
"soilc"
#' A table of soil parameters
#'
#' A table of soil parameters for different soil types
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Soil.type}{description of soil type}
#'   \item{Number}{an integer value corresponding to each soil type}
#'   \item{Smax}{Volumetric water content at saturation (m^3 / m^3)}
#'   \item{Smin}{Residual water content (m^3 / m^3)}
#'   \item{Ksat}{Saturated hydraulic conductivity (kg / m^3 / day)}
#'   \item{b}{Shape parameter for Campbell soil moisture model (dimensionless, > 1)}
#'   \item{psi_e}{Matric potential (J / m^3)}
#'   \item{Vq}{Volumetric quartz content of soil}
#'   \item{Vm}{Volumetric mineral content of soil}
#'   \item{Vo}{Volumetric organic content of soil}
#'   \item{Mc}{Mass fraction of clay}
#'   \item{rho}{Soil bulk density (Mg / m^3)}
#'   \item{mult}{soil moisture model radiation coefficient}
#'   \item{rmu}{soil moisture model rainfall coefficient}
#'   \item{a}{soil moisture model deeper layer multiplier coefficient}
#'   \item{pwr}{soil moisture model deeper layer power coefficient}
#' }
#' @source: \url{https://onlinelibrary.wiley.com/doi/full/10.1002/ird.1751}
"soilparameters"
#' A dataset of vegetation parameters.
#'
#' An object of class vegparams for the area for the area bounded by 169475, 169525, 12475, 12525
#' (xmin, xmax, ymin, ymax) using the Ordnance Survey GB Grid Reference system (CRS: 27700).
#'
#' @format A list of  the following objects:
#' \describe{
#'   \item{pai}{a PackedSpatRaster of monthly plant area index values, 12 layers}
#'   \item{hgt}{a PackedSpatRaster object of vegetation heights (m)}
#'   \item{x}{a PackedSpatRaster object of the ratio of vertical to horizontal projections of leaf foliage}
#'   \item{gsmax}{a PackedSpatRaster object maximum stomatal conductances (mol / m^2 / s)}
#'   \item{leafr}{a PackedSpatRaster object of leaf reflectance values}
#'   \item{clump}{a PackedSpatRaster of monthly values between 0 and 1 indicating the fraction of radiation passing through larger gaps in the canopy, here estimated using [clumpestimate()], 12 layers}
#'   \item{leafd}{a PackedSpatRaster object of leaf diameters (m)}
#'   \item{leaft}{a PackedSpatRaster object of leaf transmittance}
#' }
#' "vegp"
#'
#' A dataset of soilparameters parameters for the point model.
#'
#' An object of class groundparams
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Soil.type}{description of soil type}
#'   \item{Smax}{Volumetric water content at saturation (m^3 / m^3)}
#'   \item{Smin}{Residual water content (m^3 / m^3)}
#'   \item{alpha}{Shape parameter of the van Genuchten model (cm^-1)}
#'   \item{n}{Pore size distribution parameter (dimensionless, > 1)}
#'   \item{Ksat}{Saturated hydraulic conductivity (cm / day)}
#'   \item{Vq}{Volumetric quartz fraction of soil}
#'   \item{Vm}{Volumetric mineral fraction of soil}
#'   \item{Vo}{Volumetric organic fraction of soil}
#'   \item{Mc}{Mass fraction of clay}
#'   \item{rho}{Soil bulk density (Mg / m^3)}
#'   \item{b}{Shape parameter for Campbell model (dimensionless, > 1)}
#'   \item{psi_e}{Matric potential (J / m^3)}
#'   \item{mult}{soil moisture model radiation coefficient}
#'   \item{rmu}{soil moisture model rainfall coefficient}
#'   \item{a}{soil moisture model deeper layer multiplier coefficient}
#'   \item{pwr}{soil moisture model deeper layer power coefficient}
#' }
#' @source: \url{https://onlinelibrary.wiley.com/doi/full/10.1002/ird.1751}
"soilparamsp"
#' A 2.5 degree resolution global dataset indicating snow environments
#'
#' A spatial dataset of indicating snow environments (0 = Alpine, 1 = maritime, 2 = Prairie, 3 = Taiga, 4 = Tundra)
#'
#' @format A PackedSpatRaster object with 72 rows and 144 columns
#' @source derived from \url{https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world}
"snowenv"
