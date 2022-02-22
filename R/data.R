#' A data frame of hourly weather
#'
#' A data frame of hourly weather in 2017 at Caerthillean Cove, Lizard, Cornwall (49.96807N, 5.215668W)
#'
#' @format a data frame with the following elements:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times}
#'  \item{temp}{temperature (degrees C)}
#'  \item{relhum}{relative humidity (percentage)}
#'  \item{pres}{atmospheric press (kPa)}
#'  \item{swrad}{Total incoming shortwave radiation (W / m^2)}
#'  \item{difrad}{Diffuse radiation (W / m^2)}
#'  \item{skyem}{Sky emissivity (0-1)}
#'  \item{windspeed}{Wind speed (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#' }
"climdata"
#'
#' A one m resolution digital terrain dataset.
#'
#' A dataset of elevations (m) for Caerthillean Cove: an area bounded by 169475, 169525, 12475, 12525
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system (CRS: 27700).
#'
#' @format A RasterLayer object with 50 rows and 50 columns
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
#' A raster layer of habitat types expressed as integers for Caerthillean Cove: an area bounded by 169475, 169525, 12475, 12525
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system (CRS: 27700). Integer values
#' correspond to habitat types listed in [vegpfromhab()].
#'
#' @format A RasterLayer object with 50 rows and 50 columns
"habitats"
#' Daily rainfall
#'
#' A vector of daily rainfall for 169500, 12500 (x, y) using the Ordance Survey GB Grid Reference system (CRS: 27700).
#'
#' @format A vector of daily rainfall (mm/day)
"rainfall"
#' A dataset of soil characteristics.
#'
#' An object of class soilcharac for the area bounded by 169475, 169525, 12475, 12525
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system (CRS: 27700).
#'
#' @format A list of  the following objects:
#' \describe{
#'   \item{soiltype}{a RasterLayer object of numeric integers corresponding to soil types given in soilparameters}
#'   \item{groundr}{a RasterLayer object of soil reflectance of shortwave radiation}
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
#'   \item{int}{temperature model intercept}
#'   \item{t1}{temperature model radiation coefficient}
#'   \item{t2}{temperature model soil moisture coefficient}
#'   \item{t3}{temperature model wind coefficient}
#'   \item{t4}{temperature model soil moisture * radiation interaction coefficient}
#'   \item{t5}{temperature model soil moisture * wind interaction coefficient}
#'   \item{t6}{temperature model radiation * wind interaction coefficient}
#'   \item{t7}{temperature model three-way interaction coefficient}
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
#' (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference system (CRS: 27700).
#'
#' @format A list of  the following objects:
#' \describe{
#'   \item{pai}{an array of monthly plant area index values}
#'   \item{hgt}{a RasterLayer object of vegetation heights (m)}
#'   \item{x}{ratio of vertical to horizontal projections of leaf foliage}
#'   \item{gsmax}{a Rasterlayer object maximum stomatal conductances (mol / m^2 / s)}
#'   \item{leafr}{a Rasterlayer object of leaf reflectance values}
#'   \item{clump}{a Rasterlayer object of values between 0 and 1 indicating the degree of canopy clumpiness}
#'   \item{leafd}{a Rasterlayer object of leaf diameters (m)}
#'}
"vegp"

