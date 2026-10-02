#' White-sky albedo at Caerthillian Cove
#'
#' One-metre-resolution white-sky albedo for a 50 x 50 m area of Caerthillian
#' Cove, Cornwall, on the British National Grid (EPSG:27700; extent 169475--169525
#' E, 12475--12525 N).
#'
#' @format A 50 x 50 cell \code{PackedSpatRaster}.
#' @source Derived from aerial imagery available through \url{https://digimap.edina.ac.uk/}.
"albedo"

#' Example hourly weather from Caerthillian Cove
#'
#' Hourly weather observations for 2017 at Caerthillian Cove, Lizard, Cornwall
#' (49.96807 N, 5.215668 W).
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{obs_time}{Date and time.}
#'   \item{temp}{Air temperature (deg C).}
#'   \item{relhum}{Relative humidity (percent).}
#'   \item{pres}{Atmospheric pressure (kPa).}
#'   \item{swdown}{Downward shortwave radiation (W/m2).}
#'   \item{difrad}{Diffuse shortwave radiation (W/m2).}
#'   \item{lwdown}{Downward longwave radiation (W/m2).}
#'   \item{windspeed}{Wind speed (m/s).}
#'   \item{winddir}{Wind direction (degrees).}
#'   \item{precip}{Precipitation (mm).}
#' }
"climdata"

#' Digital terrain model for Caerthillian Cove
#'
#' One-metre-resolution elevation (m) for a 50 x 50 m area of Caerthillian
#' Cove, Cornwall, on the British National Grid (EPSG:27700; extent 169475--169525
#' E, 12475--12525 N).
#'
#' @format A 50 x 50 cell \code{PackedSpatRaster}.
#' @source \url{http://www.tellusgb.ac.uk/}
"dtmcaerth"

#' Global climate summary data
#'
#' Global climate summaries averaged over 2008--2017 on a 94 x 192 grid.
#'
#' @format An array with 94 rows, 192 columns and five climate variables:
#' \describe{
#'   \item{1}{Mean annual temperature (deg C).}
#'   \item{2}{Coefficient of variation in temperature (K).}
#'   \item{3}{Mean annual rainfall (mm per year).}
#'   \item{4}{Coefficient of variation in annual rainfall (mm per 0.25 days).}
#'   \item{5}{Calendar month with the greatest rainfall (1--12).}
#' }
#' @source \url{http://www.ncep.noaa.gov/}
"globclim"

#' Soil properties at Caerthillian Cove
#'
#' Spatial soil class and ground shortwave reflectance for the Caerthillian Cove
#' example landscape.
#'
#' @format A list with:
#' \describe{
#'   \item{soiltype}{Integer soil-class \code{PackedSpatRaster}; codes correspond to \code{soilparamstable$Number}.}
#'   \item{groundr}{Ground/litter shortwave reflectance \code{PackedSpatRaster}.}
#' }
"soilc"

#' Vegetation properties at Caerthillian Cove
#'
#' Spatial vegetation structure and optical properties for the Caerthillian Cove
#' example landscape. Bare-ground cells have zero height and plant area index,
#' habitat class 16, and no value for the other vegetation properties.
#'
#' @format A list containing:
#' \describe{
#'   \item{pai}{Monthly plant area index, 12-layer \code{PackedSpatRaster}.}
#'   \item{hgt}{Vegetation height (m), \code{PackedSpatRaster}.}
#'   \item{x}{Leaf-angle distribution parameter, \code{PackedSpatRaster}.}
#'   \item{gsmax}{Maximum stomatal conductance (mol/m2/s), \code{PackedSpatRaster}.}
#'   \item{leafr}{Leaf shortwave reflectance, \code{PackedSpatRaster}.}
#'   \item{clump}{Monthly canopy clumping factor, 12-layer \code{PackedSpatRaster}.}
#'   \item{leaft}{Leaf shortwave transmittance, \code{PackedSpatRaster}.}
#'   \item{habitat}{Integer habitat class (1--16) from which plant functional
#'     types are assigned, \code{PackedSpatRaster}: here 6 closed shrubland,
#'     7 open shrubland, 10 short grassland, 11 tall grassland and 16 barren or
#'     sparsely vegetated.}
#' }
"vegp"

#' Campbell soil physical parameters by soil type
#'
#' Soil physical parameters used by the Campbell soil-water and soil-heat
#' calculations. This is a Campbell-only table: it has no van Genuchten
#' \code{n} column, since the model derives the Campbell hydraulic
#' conductivity exponent internally as \code{2*b+3} rather than storing it
#' separately. Classes 1-11 follow Campbell's texture ordering from sand to
#' clay; \code{Silt} is class 12.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{Soil.type}{Soil type.}
#'   \item{Number}{Integer soil-class code, as used by \code{soilc$soiltype}.}
#'   \item{Smax}{Volumetric water content at saturation (m3/m3).}
#'   \item{Smin}{Residual volumetric water content (m3/m3).}
#'   \item{Ksat}{Saturated hydraulic conductivity (kg s/m3).}
#'   \item{Vq}{Volumetric quartz content.}
#'   \item{Vm}{Volumetric mineral content.}
#'   \item{Vo}{Volumetric organic content.}
#'   \item{Mc}{Clay mass fraction.}
#'   \item{rho}{Soil bulk density (Mg/m3).}
#'   \item{b}{Campbell soil-water retention parameter (dimensionless).}
#'   \item{psi_e}{Air-entry matric-potential magnitude (J/kg), stored as a positive value.}
#' }
#' @source \url{https://onlinelibrary.wiley.com/doi/full/10.1002/ird.1751}
"soilparamstable"

#' Vegetation parameters by plant functional type
#'
#' Structural, radiative, photosynthetic and hydraulic parameters for the plant
#' functional types represented by the model, mirroring the convention used in
#' JULES. Each row is one parameter; PFT columns contain its stored value and
#' \code{multiplier} converts that value to the working units reported in
#' \code{units}. \code{rpmin} (minimum whole-plant hydraulic resistance) is
#' stored as a precomputed constant per type rather than derived at run time,
#' since it depends only on static structural parameters. \code{h} (canopy
#' height) is kept as its own row because it is also used elsewhere in the
#' model (canopy structure, wind profile).
#' \code{root50} and \code{root95} are the depths (m) above which 50\% and 95\% of
#' roots lie, from Schenk & Jackson (2002, Ecological Monographs 72:311-328,
#' Table 4), and set the depth profile of root water uptake.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{varname}{Model parameter name.}
#'   \item{Description}{Parameter description.}
#'   \item{units}{Working units.}
#'   \item{BET.Tr}{Broadleaf evergreen tree, tropical.}
#'   \item{BET.Te}{Broadleaf evergreen tree, temperate.}
#'   \item{BDT}{Broadleaf deciduous tree.}
#'   \item{NET}{Needleleaf evergreen tree.}
#'   \item{NDT}{Needleleaf deciduous tree.}
#'   \item{C3}{C3 grass.}
#'   \item{C4}{C4 grass.}
#'   \item{ESh}{Evergreen shrub.}
#'   \item{DSh}{Deciduous shrub.}
#'   \item{multiplier}{Conversion factor from stored values to working units.}
#' }
"PFTparams"

#' Global snow-environment classes
#'
#' A 2.5-degree global classification of snow environments: 0 = Alpine,
#' 1 = maritime, 2 = Prairie, 3 = Taiga, and 4 = Tundra.
#'
#' @format A 72 x 144 cell \code{PackedSpatRaster}.
#' @source Derived from \url{https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world}.
"snowenv"
