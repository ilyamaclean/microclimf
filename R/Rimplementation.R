#' Create object of class microin
#'
#' @description The function `modelin` creates an object of class microin
#' which unpacks various component inputs and reformats as required
#' for running the R version of the model.
#' @param micropoint either an object of class `micropoint` as returned
#' by [runpointmodel()] or [subsetpointmodel()] or a list of micropoint
#' objects as returned by [runpointmodela()] or [subsetpointmodela()],
#' which assumes input climate data are spatially variable.
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a SpatRaster object of elevations in metres (see details)
#' @param dtmc a SpatRaster object giving the resolution, spatial extent, and projection
#' of the weather data used in `micropoint`. Ignored if climate data provided
#' to `micropoint` as a data.frame. Must give elevations in metres if
#' climate data were arrays and `altcorrect` > 0 or if setting runchecks to TRUE.
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction, see details)
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' @details The array of Plant Area index values in `vegp` must
#' of the same x and y dims as `dtm` but can contain any number of repeated
#' measures up to the number of entries in `weather`. Data are interpolated to the
#' time increment of `weather`. Other vegetation paramaters, including vegetation
#' height are assumed time-invarient. If these are time-varient, use the c++ version
#' of the model. The SpatRaster datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined.
#' If not NA, the units of dtmc must match dtm and must be elevation data if an
#' altitude correction is applied. If `altcorrect`>0, the elevation
#' difference between each pixel of dtm and dtmc is calculated and
#' an elevation lapse rate correction is applied to the temperature and
#' pressure data to account for these elevation differences. If `altcorrect`= 1, a fixed
#' lapse rate of 5 degrees per 100m is applied to the temperature data. If
#' `altcorrect`= 2, humidity-dependent lapse rates are calculated and applied.
#' See also details for  modelin.
#' @export
#' @rdname modelin
#' @examples
#' # ================================================================= #
#' # Preparing model inputs using climate data provided as a data.frame
#' # ================================================================= #
#' # run and subset point model
#' micropoint<-runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc)
#' micropoint <- subsetpointmodel(micropoint)
#' micro <- modelin(micropoint, vegp, soilc, dtmcaerth)
#' # ================================================================= #
#' # Preparing model inputs using climate data provided as arrays
#' # ================================================================= #
#' # Create dummy array data
#' .ta<-function(x,dtm,xdim=5,ydim=5) {
#'   a<-array(rep(x,each=ydim*xdim),dim=c(ydim,xdim,length(x)))
#'   .rast(a,dtm)
#' }
#' dtm <- rast(dtmcaerth)
#' climarrayr<-list(temp = .ta(climdata$temp, dtm),
#'   relhum = .ta(climdata$relhum, dtm),
#'   pres = .ta(climdata$pres, dtm),
#'   swdown = .ta(climdata$swdown, dtm),
#'   difrad = .ta(climdata$difrad, dtm),
#'   lwdown = .ta(climdata$lwdown, dtm),
#'   windspeed = .ta(climdata$windspeed, dtm),
#'   winddir = .ta(climdata$winddir, dtm),
#'   precip = .ta(climdata$precip, dtm))
#' tme <- as.POSIXlt(climdata$obs_time, tz="UTC")
#' # Run and subset point model array
#' micropointa <- runpointmodela(climarrayr, tme, reqhgt = 0.05, dtmcaerth, vegp, soilc)
#' micropointa <- subsetpointmodela(micropointa)
#' # Prepare grid model input
#' dtmc <- resample(dtm, climarrayr$temp[[1]])
#' micro <- modelin(micropointa, vegp, soilc, dtm, dtmc)
modelin <- function(micropoint, vegp, soilc, dtm, dtmc = NA, altcorrect = 0, runchecks = TRUE) {
  if (class(micropoint) == "micropoint") {
    micro<-.modelin(micropoint,vegp,soilc,dtm,runchecks)
  } else {
    micro<-.modelina(micropoint,vegp,soilc,dtm,dtmc,altcorrect,runchecks)
  }
  return(micro)
}
#' Spatially distribute soil moisture by topographic wetness index
#'
#' The function `soilmdistribute` spatially distrubutes soil moisture by the Bevan and Kirkby
#' topographic wetness index
#'
#' @param micro an object of class microin as returned by e.g. modelin
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness
#' @return a 3D array of soil moistures with the same x and y dims as `dtm` and z
#' equivelent to length(soilm)
#' @import terra
#' @export
#' @rdname soilmdistribute
#' @examples
#' micropoint <- runpointmodel(climdata, 0.05, dtmcaerth, vegp, soilc)
#' micropoint <- subsetpointmodel(micropoint)
#' micro <- modelin(micropoint, vegp, soilc, dtmcaerth)
#' micro <- soilmdistribute(micro) # default tfact
#' plot(rast(micro$soilm[,,1]))
soilmdistribute <- function(micro, tfact = 1.5) {
  if (class(micro$dtm)[1] == "PackedSpatRaster") micro$dtm<-rast(dtm)
  twi<-.topidx(micro$dtm)
  Smin<-.rta(micro$Smin,dim(micro$soilm)[3])
  Smax<-.rta(micro$Smax,dim(micro$soilm)[3])
  s<-which(micro$soilm>=Smax)
  micro$soilm[s]<-Smax[s]-0.001
  s<-which(micro$soilm<=Smin)
  micro$soilm[s]<-Smin[s]+0.001
  theta<-(micro$soilm-Smin)/(Smax-Smin)
  lt<-log(theta/(1-theta))
  ltwi<-log(.is(twi))/tfact
  me<-mean(ltwi,na.rm=T)
  smout<-lt+.rta(ltwi-me,dim(micro$soilm)[3])
  smout<-1/(1+exp(-smout))
  smout<-smout*(Smax-Smin)+Smin
  smout<-.rast(smout, micro$dtm)
  smout<-mask(smout,micro$dtm)
  micro$soilm<-as.array(smout)
  micro$progress<-1
  return(micro)
}
#' Applies two-stream radiation model
#'
#' @description The function `twostream` applies a varient of the Sellers two-stream
#' radiation model to calaculate long and showtwave radiation absorbed by the ground and
#' canopy and shortwave radiation absorbed by leafs.
#'
#' @param micro an object of class `micro` as returned by [modelin()]
#' @param reqhgt height at which temperatures are required (m). Negative if below ground
#' @param pai_a a SpatRaster of plant area index values above `reqhgt`. Determined from
#' total `pai` if not supplied. Must patch the dimensions of `vegp$pai` if supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' @return an object of class micro with additional terms added for subsequent modelling
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @rdname twostream
twostream<-function(micro, reqhgt = 0.05, pai_a = NA, tfact=1.5) {
  if (micro$progress<1) micro<-soilmdistribute(micro,tfact)
  # Calculate slope and aspect
  slope<-terra::terrain(micro$dtm,'slope')
  aspect<-terra::terrain(micro$dtm,'aspect')
  slope[is.na(slope)]<-0
  aspect[is.na(aspect)]<-0
  slope<-.is(mask(slope,micro$dtm))
  aspect<-.is(mask(aspect,micro$dtm))
  # Calculate obstime data.frame
  obstime<-data.frame(year=micro$tme$year+1900,month=micro$tme$mon+1,day=micro$tme$mday,
                      hour=micro$tme$hour+micro$tme$min/60+micro$tme$sec/3600)
  # Calculate si etc
  micro<-solargrid(slope,aspect,obstime,micro)
  # Calculate shadowmask
  hor<-array(NA,dim=c(dim(micro$dtm)[1:2],24))
  for (i in 1:24) hor[,,i]<-.horizon(micro$dtm,(i-1)*15)
  i<-round(micro$sazi/15,0)+1; i[i==25]<-1
  hora<-hor[,,i]
  # Calculate terrain shading
  shadowmask<-hora*0+1
  shadowmask[hora>tan(pi/2-micro$zen)]<-0
  micro$si<-micro$si*shadowmask
  # Calculate sky view
  msl<-tan(apply(atan(hor),c(1,2),mean))
  svf<-0.5*cos(2*msl)+0.5
  micro$svfa<-.rta(svf,length(obstime$year))
  # === (1e) Calculate leaf area above
  pai_a<-.expandpaia(pai_a, micropoint)
  fd<-.foliageden(reqhgt,micro$veghgt,micro$pai,pai_a)
  micro$leafden<-fd$leafden
  micro$paia<-fd$pai_a
  micro<-twostreamgrid(reqhgt,micro)
  micro$progress<-2
  # Clean micro
  micro$si<-NULL
  micro$svfa<-NULL
  micro$lat<-NULL
  micro$long<-NULL
  micro$vegx<-NULL
  micro$lref<-NULL
  micro$ltra<-NULL
  micro$clump<-NULL
  return(micro)
}
#' Downscales wind
#'
#' @description The function `wind` downscales wind speed accounting for
#' vegetation and terrain
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above ground for which wind speeds are wanted. If negative (below ground) wind friction velocity only is returned
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @return an object of class microin with the following components added:
#' \describe{
#'   \item{uf}{wind friction velocity (m/s)}
#'   \item{uz}{if `reqhgt > 0` wind speed at height `reqhgt` (m/s)}
#'   \item{gHa}{convective conductance (mol/m^2/s)}
#' }
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @rdname wind
#' @details In downscaling wind, two processes are accounted for. Firstly the drag effects
#' of vegetation on wind , which ultimately dictate the wind height profile. Secondly,
#' the effects of local vegetation and terrain on wind speed. Terrain effects are calculated
#' by applying the topographic shelter coefficient described in Maclean et al (2019) Methods
#' in Ecology and Evolution, 10:280-290.
wind <- function(micro, reqhgt = 0.05, pai_a = NA, tfact = 1.5) {
  if (micro$progress < 2) micro<-twostream(micro,reqhgt,pai_a,tfact)
  # Calculate wind shelter coefficient
  s<-1
  if (res(micro$dtm)[1]<=100) s<-10
  wsa<-.windsheltera(micro$dtm,micro$zref,s)
  dsm<-.is(micro$dtm)+micro$veghgt[,,1]
  micro$ws<-.windshelter(micro$winddir,dsm,2,s,wsa)
  micro<-windgrid(reqhgt, micro)
  micro$progress<-3
  # Clean micro
  micro$d<-NULL
  micro$zm<-NULL
  micro$u2<-NULL
  micro$ws<-NULL
  return(micro)
}
#' Calculates ground surface temperature (hourly)
#'
#' @description The function `soiltemp` estimates ground surface temperature
#'
#' @param micro an object of class `microin` as returned by [modelin()] (see details)
#' @param reqhgt height above ground at which model outputs are needed (m).
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @return an object of class micro with the following components added:
#' \describe{
#'   \item{T0}{ground surface temperature (deg C)}
#'   \item{G}{ground heat flux (W/m^2)}
#'   \item{...}{additional terms needed for subsequent modelling}
#' }
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @rdname soiltemp
#' @seealso [wind()]
#' @details if [wind()] has not been run to add additional elements to `micro`
#' it is automatically called.
soiltemp  <- function(micro, reqhgt = 0.05, pai_a = NA, tfact = 1.5) {
  # run two-stream and wind functions if not run
  if (micro$progress<3) {
    micro<-wind(micro,reqhgt,pai_a,tfact)
  }
  # generate soilparamsp
  micro<-soiltempgrid(micro)
  micro$progress<-4
  # Clean micro
  micro$rho<-NULL
  micro$Vm<-NULL
  micro$Vq<-NULL
  micro$Mc<-NULL
  micro$soilb<-NULL
  micro$psi_e<-NULL
  return(micro)
}
#' Estimate temperature and humidity at specified height above ground
#'
#' @description The function `aboveground` runs the above ground component of
#' the microclimate model.
#'
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above ground for which wind speeds are wanted. If negative (below ground) wind friction velocity only is returned
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @return a list of the following:
#' \describe{
#'   \item{Tz}{Air temperatures at height `reqhgt` (deg C)}
#'   \item{tleaf}{Leaf temperatures at height `reqhgt` (deg C)}
#'   \item{T0}{Ground surface temperatures (deg C)}
#'   \item{relhum}{Relative humidities at height `reqhgt` (percentage)}
#'   \item{windspeed}{Wind speeds at height `reqhgt` (m/s)}
#'   \item{Rdirdown}{Flux density of downward direct radiation at `reqhgt` (W/m^2 - on the horizontal}
#'   \item{Rdirdown}{Flux density of downward diffuse radiation at `reqhgt` (W/m^2)}
#'   \item{Rlwdown}{Flux density of downward longwave radiation at `reqhgt` (W/m^2)}
#'   \item{Rswup}{Flux density of upward shorwtave radiation at `reqhgt` (W/m^2), assumed diffuse}
#'   \item{Rlwup}{Flux density of downward longwave radiation at `reqhgt` (W/m^2)}
#' }
#' @seealso [belowground()]
#' @details `pai_a` is used to calculate the radiation intercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calculated from total plant area index by
#' assuming leaf density within the canopy is uniformly vertically distributed. If supplied
#' it must have the same dimensions as micro$pai. I.e. with the same x and y dims as the
#' the supplied dtm and values for each hour as the z dimension.
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @rdname aboveground
aboveground<- function(micro, reqhgt = 0.05, pai_a = NA, tfact = 1.5) {
  if (micro$progress<4) {
    micro<-soiltemp(micro,reqhgt,pai_a,tfact)
  }
  h<-dim(micro$tc)[3]
  micro$Smin<-.rta(micro$Smin,h)
  micro$Smax<-.rta(micro$Smax,h)
  mout<-abovegrid(reqhgt, micro)
  return(mout)
}
#' Estimate temperature at specified height below ground
#'
#' @description The function `belowground` runs the below ground component of
#' the microclimate model.
#'
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above ground for which wind speeds are wanted. If negative (below ground) wind friction velocity only is returned
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @return a list of the following:
#' \describe{
#'   \item{Tz}{Soil temperatures at height `reqhgt` below ground (deg C)}
#'   \item{T0}{Ground surface temperatures (deg C)}
#'   \item{soilm}{Soil moisture in top 10 cm of soil (volumetric water fraction)}
#' }
#' @seealso [aboveground()]
#' @details `pai_a` is used to calculate the radiation intercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calculated from total plant area index by
#' assuming leaf density within the canopy is uniformly vertically distributed. If supplied
#' it must have the same dimensions as micro$pai. I.e. with the same x and y dims as the
#' the supplied dtm and values for each hour as the z dimension.
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @rdname belowground
belowground<- function(micro, reqhgt = -0.05, pai_a = NA, tfact = 1.5) {
  if (micro$progress<4) {
    micro<-soiltemp(micro,reqhgt,pai_a,tfact)
  }
  year<-micro$tme$year[1]+1900
  hiy<-ifelse(year%%4==0,366*24,365*24)
  dif<-round((as.numeric(micro$tme[25])-as.numeric(micro$tme[24]))/3600,0)
  complete<-FALSE
  if (dif == 1) complete<-TRUE
  mout<-belowgrid(reqhgt, micro, hiy, complete)
  return(mout)
}
