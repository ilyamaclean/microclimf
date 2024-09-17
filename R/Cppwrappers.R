#' Runs point microclimate model
#'
#' The function `runpointmodel` runs the point microclimate model
#'
#' @param weather a data.frame of weather variables (see details)
#' @param reqhgt height for which temperatures are needed (used only when reqhgt < 0 to calculate tmeperature below ground)
#' @param dtm a SpatRaster of elevations for the study area
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param zref height above ground (m) of temperature measurements in weather
#' @param windhgt height above ground (m) of wind speed data in weather
#' @param soilm optional vector of soil moisture values in upper 10 cm of the soil (calculated if not supplied)
#' @param matemp optionally mean annual temperature. Only used for refining below-ground temperature estimates and
#' calculate as the mean of weather data temperatures if not provided.
#' @param dTmx optional maximum amount by which canopy or ground surface temperatures can exceed air temperatures.
#' Included to ensure model convergence:
#' @param maxiter optional integer indicating the maximum number of iterations (see details)
#' @param yearG optional logical indicating whether or not to include annual ground heat flux cycle
#' Reduces tile effects:
#' @param lat optional central latitude of study area (removes tile effects when running in tiles)
#' @param long optional central longitude of study area (removes tile effects when running in tiles)
#' Used by runpointmodela:
#' @param vegp_p optional vector of point model vegetation parameters
#' @param groundp_p = optional vector of point model ground parameters
#' @param soiltype = optional modal soil type
#' @param mxhgt = optional maximum height of vegetation
#' @return a list of the following:
#' (1) weather - a data.frame of weather variables, but with temperature and
#' wind speed height-adjusted to be above canopy if necessary (see details)
#' (2) dfo - a data.frame of microclimate point model outputs required for running grid model
#' (3) if reqhgt < 0, a vector iof point model temperatures below ground
#' (4) lat - the latitude of the centroid of the study area for which the point model was run (decimal degrees)
#' (5) long - the longitude of the centroid of the study area for which the point model was run (decimal degrees)
#' (6) zref - the height to which weather data have been height adjusted (see details)
#' (7) subs & tmeorig used subsequently by grid model to handle time-variant vegetation inputs
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. As not all wind measurements are at reference height, the height of the wind speed measurement
#' must be specified if not 2 m. To enable calculation of below-canopy wind and temperature profiles
#' in tall canopy, the wind speed and temperature data are  adjusted to give values for a height at
#' the maximum vegetation height if the tallest vegetation exceeds two metres.
#' For doing so a stand vegetation surface typical of that in which a weather station would be located is assumed.
#' The parameter `maxiter` sets the maximum number of times the model is iterated to achieve
#' convergence. Increasing this value improves accuracy at the expense of computation time.
#' The array of Plant Area index values and clumping factors in `vegp` must
#' have the same x and y dims as `dtm` but can contain any number of repeated
#' measures up to the number of entries in `weather`. Data are interpolated to the
#' time increment of `weather`. Other vegetation paramaters are assumed time-invarient.
#' @import stats
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @rdname runpointmodel
#' @examples
#' # Run model:
#' micropoint<-runpointmodel(climdata,0.05,dtmcaerth,vegp,soilc)
#' # Plot canopy heat exchange surface temperature
#' microp<-micropoint$dfo
#' plot(microp$Tc,type="l") # temperature of canopy surface
runpointmodel<-function(weather, reqhgt = 0.05, dtm, vegp, soilc, runchecks = TRUE,
                        zref = 2, windhgt = zref, soilm = NA, matemp = NA, dTmx = 25,
                        maxiter = 20, yearG = TRUE, lat = NA, long = NA, vegp_p = NA,
                        groundp_p = NA, soiltype = NA, mxhgt = NA) {
  # Unpack and check variables
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  if (runchecks) {
    rc<-checkinputs(weather,vegp,soilc,dtm)
    weather<-rc$weather
    vegp<-rc$vegp
    soilc<-rc$soilc
  }
  if (is.na(matemp)) matemp<-mean(weather$temp)
  # do wind height adjustment if necessary
  if (zref != windhgt) {
    weather$windspeed<-weather$windspeed*log(67.8*zref-5.42)/log(67.8*windhgt-5.42)
  }
  # Create date data.frame of obstime
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,hour=tme$hour+tme$min/60+tme$sec/3600)
  weather$obs_time<-NULL
  # Create point model-type vegetation properites input
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  if (class(vegp_p) == "logical") vegp_p<-.sortvegp(vegp,method="P")
  if (class(groundp_p) == "logical") groundp_p<-.sortsoilc(soilc,method="P")
  # Peform weather height adjustment
  if (class(mxhgt) == "logical") mxhgt<-.mfr(vegp$hgt,max)
  zout<-ifelse(mxhgt>2,mxhgt,2)
  if (class(vegp$x)[1] == "PackedSpatRaster") vegp$x<-rast(vegp$x)
  if (class(lat) == "logical") {
    ll<-.latlongfromraster(vegp$x)
  } else ll<-data.frame(lat=lat,long=long)
  if (zout > zref) {
    tst<-1
    ctr<-0
    while (tst > 0) {
      weather2<-weatherhgtCpp(obstime,weather,zref,zout,zout,ll$lat,ll$long)
      tt<-mean(weather2$temp)
      if (is.na(tt) == FALSE) tst<-0
      ctr<-ctr+1
      if (ctr > 5) tst<-0
    }
    tt<-mean(weather2$temp)
    if (is.na(tt)) {
      weather<-weather
    } else {
      weather<-weather2
    }
    zref<-zout
  }
  # Set minimum wind speed to avoid convergence issues
  weather$windspeed[weather$windspeed<0.5]<-0.5
  # Run soil moisture model if not provided
  if (class(soilm)=="logical") {
    if (is.na(soiltype)) {
      ii<-.getmode(.is(soilc$soiltype))
    } else ii<-soiltype
    soilm<-soilmCpp(weather,soilparamsp$rmu[ii],soilparamsp$mult[ii],soilparamsp$pwr[ii],
                    soilparamsp$Smax[ii],soilparamsp$Smin[ii],soilparamsp$Ksat[ii],
                    soilparamsp$a[ii])
    soilm<-stats::spline(soilm,n=length(weather$temp))$y
  }
  # Run big Leaf model
  microp<-BigLeafCpp(obstime,weather,vegp_p,groundp_p,soilm,ll$lat,ll$long,dTmx,zref,maxiter,0.5,0.5,0.1,yearG)
  # Create data.frame for running post processing
  dfp<-data.frame(windspeed=weather$windspeed,tc=weather$temp,rh=weather$relhum,
                  pk=weather$pres,uf=microp$uf,soilm=soilm,RabsG=microp$RabsG)
  dfo<-as.data.frame(pointmprocess(dfp,zref,vegp_p[1],vegp_p[2],groundp_p[5],groundp_p[6],groundp_p[7],groundp_p[8]))
  dfo$G<-microp$G
  dfo$soilm<-soilm
  dfo$Tg<-microp$Tg
  dfo$Tc<-microp$Tc
  if (reqhgt < 0) {
    Tbz<-.soilbelowT(dfo,reqhgt)
  } else Tbz<-NA
  # return outputs
  weather$obs_time<-tme
  subs<-c(1:length(dfo$Tg))
  out<-list(weather=weather,dfo=dfo,Tbz=Tbz,lat=ll$lat,long=ll$long,zref=zref,subs=subs,tmeorig=tme,matemp=matemp)
  class(out)<-"micropoint"
  return(out)
}
#' Runs point microclimate model over each grid cell of an array
#'
#' The function `runpointmodela` runs the point microclimate model for
#' each grid cell of a climarray provided as an array
#'
#' @param climarrayr a list of weather variables provided as SpatRasters (see details)
#' @param reqhgt height for which temperatures are needed (used only when reqhgt < 0 to calculate tmeperature below ground)
#' @param tme POSIXlt object giving the dates and times for each weather variable stored in the array
#' @param dtm a SpatRaster of elevations for the study area
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param zref height above ground (m) of temperature measurements in weather
#' @param windhgt height above ground (m) of wind speed data in weather
#' @param soilm optional vector of soil moisture values in upper 10 cm of the soil (calculated if not supplied)
#' @param matemp optionally a single numeric value of approximate mean annual temperatures averaged across study region.
#' Only used for refining below-ground temperature estimates and calculate as the mean of
#' weather data temperatures if not provided.
#' @param dTmx optional maximum amount by which canopy or ground surface temperatures can exceed air temperatures.
#' Included to ensure model convergence
#' @param maxiter optional integer indicating the maximum number of iterations (see details)
#' @param yearG optional logical indicating whether or not to include annual ground heat flux cycle
#' @param lat optional central latitude of study area (removes tile effects when running in tiles)
#' @param long optional central longitude of study area (removes tile effects when running in tiles)
#' @return a list of point microclimate model outputs (as returned by runmicropoint) but for each grid cell
#' @details The units of `climarrayr` must follow those in the dataset `climdata`.
#' It must be a list with each component of the list a SpatRaster, named using the same
#' names as the column headers in climdata (e.g. temp for temperature), excluding `obs_time`.
#' As not all wind measurements are at reference height, the height of the wind speed measurement
#' must be specified if not 2 m. To enable calculation of below-canopy wind and temperature profiles
#' in tall canopy, the wind speed and temperature data are  adjusted to give values for a height at
#' the maximum vegetation height if the tallest vegetation exceeds two metres.
#' For doing so a stand vegetation surface typical of that in which a weather station would be located is assumed.
#' The parameter `maxiter` sets the maximum number of times the model is iterated to achieve
#' convergence. Increasing this value improves accuracy at the expense of computation time.
#' @import stats
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @rdname runpointmodela
#' @examples
#' # ======== Create dummy array datasets ========= #
#' .ta<-function(x,dtm,xdim=5,ydim=5) {
#'    a<-array(rep(x,each=ydim*xdim),dim=c(ydim,xdim,length(x)))
#'    .rast(a,dtm)
#' }
#' dtm<-rast(dtmcaerth)
#' climarrayr<-list(temp=.ta(climdata$temp,dtm),
#'                  relhum=.ta(climdata$relhum,dtm),
#'                  pres=.ta(climdata$pres,dtm),
#'                  swdown=.ta(climdata$swdown,dtm),
#'                  difrad=.ta(climdata$difrad,dtm),
#'                  lwdown=.ta(climdata$lwdown,dtm),
#'                  windspeed=.ta(climdata$windspeed,dtm),
#'                  winddir=.ta(climdata$winddir,dtm),
#'                  precip=.ta(climdata$precip,dtm))
#' # Run model
#' tme<-as.POSIXlt(climdata$obs_time,tz="UTC")
#' # Takes ~15 seconds to run
#' pointmodela <- runpointmodela(climarrayr, tme, reqhgt = 0.05, dtm, vegp, soilc)
runpointmodela<-function(climarrayr, tme, reqhgt = 0.05, dtm, vegp, soilc, matemp = NA, zref = 2, windhgt = 2, soilm = NA, dTmx = 25, maxiter = 20, yearG = TRUE)  {
  # Unpack fine-res rasters
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  mxhgt<-max(.is(vegp$hgt),na.rm=T)
  # Sort out all the rasters
  # Calculate modal wind direction
  wdir<-apply(.is(climarrayr$winddir),3,.getmode)
  if (crs(climarrayr$temp) != crs(dtm))  climarrayr$temp<-project(climarrayr$temp,crs(dtm))
  if (crs(climarrayr$relhum) != crs(dtm))  climarrayr$relhum<-project(climarrayr$relhum,crs(dtm))
  if (crs(climarrayr$pres) != crs(dtm))  climarrayr$pres<-project(climarrayr$pres,crs(dtm))
  if (crs(climarrayr$swdown) != crs(dtm))  climarrayr$swdown<-project(climarrayr$swdown,crs(dtm))
  if (crs(climarrayr$difrad) != crs(dtm))  climarrayr$difrad<-project(climarrayr$difrad,crs(dtm))
  if (crs(climarrayr$lwdown) != crs(dtm))  climarrayr$lwdown<-project(climarrayr$lwdown,crs(dtm))
  if (crs(climarrayr$windspeed) != crs(dtm))  climarrayr$windspeed<-project(climarrayr$windspeed,crs(dtm))
  if (crs(climarrayr$precip) != crs(dtm))  climarrayr$precip<-project(climarrayr$precip,crs(dtm))
  r<-climarrayr$temp
  pb <- utils::txtProgressBar(min = 0, max = dim(r)[1]*dim(r)[2], style = 3)
  # Resample vegetation and soil characteristics
  vegp<-.resamplev(vegp,r)
  af<-res(r)[1]/res(soilc$soiltype)[1]
  soilc$soiltype<-aggregate(soilc$soiltype,af,"modal",na.rm=T)
  soilc$groundr<-aggregate(soilc$groundr,af,na.rm=T)
  soilc$soiltype<-resample(soilc$soiltype,r,method="mode")
  soilc$groundr<-resample(soilc$groundr,r)
  ll<-.latslonsfromr(r)
  lats<-ll$lats
  lons<-ll$lons
  k<-1
  pointo<-list()
  soiltype<-.getmode(.is(soilc$soiltype))
  for (i in 1:dim(r)[1]) {
    for (j in 1:dim(r)[2]) {
      climdf<-.todf(climarrayr,i,j,tme,wdir)
      if (class(soilm) == "logical") {
        soilmo<-NA
      } else soilmo<-soilm[i,j,]
      if (is.na(climdf$temp[1]) == FALSE & is.na(.is(vegp$hgt)[i,j]) == FALSE) {
        vegp_p<-.tovp(vegp,i,j)
        gp<-.togp(soilc,i,j)
        groundp_p<-gp$groundp_p
        pointo[[k]]<-runpointmodel(climdf,reqhgt,dtm,vegp,soilc,runchecks = FALSE,zref,windhgt,soilmo,matemp,dTmx,maxiter,yearG,
                                   lats[i,j],lons[i,j],vegp_p,groundp_p,soiltype,mxhgt)
      } else pointo[[k]]<-NA
      k<-k+1
      utils::setTxtProgressBar(pb,(i-1)*dim(r)[2]+j)
    }
  }
  return(pointo)
}
#' Runs the grid microclimate model
#'
#' @description The function `runmicro` runs the grid version of the microclimate model
#' @param micropoint an object of class micropoint or a list of objects of class
#' micropoint as returned by [runpointmodel()], [runpointmodela()], [subsetpointmodel()]
#' or [subsetpointmodela()]
#' @param reqhgt  height above (postive) or below (negative) ground for which microclimate variables are required (m)
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a SpatRaster object of elevations in metres (see details)
#' @param dtmc a SpatRaster object giving the resolution, spatial extent, and projection
#' of the climate data used when running [micropointa()]. Ignored if climate data
#' used for running the point model are provided as a data.frame. Must
#' give elevations in metres if `altcorrect` > 0 or if setting runchecks to TRUE.
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction, see details)
#' @param snow optional logical indicating whether to account for snow (TRUE = yes)
#' @param snowmod optional list of snow model outputs as returned by [runsnowmodel()]. Only required if `snow` set to TRUE. Must match entries in `micropoint`.
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' @param pai_a an array of plant area index values above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()]).
#' @param out optional vector of logicals indicating which variables to
#' return ordered as for the listed outputs when `rehgt > 0'` (e.g. `out[1] = TRUE` indicates that
#' `Tz` is returned, `out[2]` that `tleaf` is returned etc). By default all variables
#' are returned.
#' @param slr an optional SpatRaster object of slope values (degrees). Calculated from the
#' dtm in micro if not supplied, but outer cells assumed flat.
#' @param apr an optional SpatRaster object of aspect values (degrees). Calculated from
#' the dtm ion micro if not supplied, but outer cells assumed flat.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from the dtm in micro if not supplied, but edge effects not accounted for.
#' @param twi optional SpatRast object of topographic wetness index values.
#' Calculated from the dtm in micro if not supplied, , but edge effects not accounted for.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from the dtm in micro if not supplied, , but edge effects not accounted for.
#' @param sva optional SpatRaster object of the skyview factor. Calculated from the
#' dtm in micro if not supplied, but cannot account for edge effects.
#' @param method on of `R` or `Cpp` (see details)
#' @return If `reqhgt > 0`:
#' \describe{
#'   \item{Tz}{Air temperatures at height `reqhgt` (deg C)}
#'   \item{tleaf}{Leaf temperatures at height `reqhgt` (deg C)}
#'   \item{relhum}{Relative humidities at height `reqhgt` (percentage)}
#'   \item{soilm}{Volumtric water fraction in 10 cm of soil (m^3/m^3)}
#'   \item{windspeed}{Wind speeds at height `reqhgt` (m/s)}
#'   \item{Rdirdown}{Flux density of downward direct radiation at `reqhgt` (W/m^2 - on the horizontal}
#'   \item{Rdifdown}{Flux density of downward diffuse radiation at `reqhgt` (W/m^2)}
#'   \item{Rlwdown}{Flux density of downward longwave radiation at `reqhgt` (W/m^2)}
#'   \item{Rswup}{Flux density of upward shorwtave radiation at `reqhgt` (W/m^2), assumed diffuse}
#'   \item{Rlwup}{Flux density of downward longwave radiation at `reqhgt` (W/m^2)}
#' }
#' If `reqhgt == 0`:
#' \describe{
#'   \item{Tz}{Soil surface temperatures (deg C)}
#'   \item{soilm}{Volumtric water fraction in 10 cm of soil (m^3/m^3)}
#'   \item{Rdirdown}{Flux density of downward direct radiation at `reqhgt` (W/m^2 - on the horizontal}
#'   \item{Rdirdown}{Flux density of downward diffuse radiation at `reqhgt` (W/m^2)}
#'   \item{Rlwdown}{Flux density of downward longwave radiation at `reqhgt` (W/m^2)}
#'   \item{Rswup}{Flux density of upward shorwtave radiation at `reqhgt` (W/m^2), assumed diffuse}
#'   \item{Rlwup}{Flux density of downward longwave radiation at `reqhgt` (W/m^2)}
#' }
#' If `reqhgt < 0`:
#' \describe{
#'   \item{Tz}{Soil temperatures at depth `-reqhgt` (deg C)}
#'   \item{soilm}{Volumtric water fraction in 10 cm of soil (m^3/m^3)}
#' }
#' Returned variables are also contingent on `out`.
#' @details
#' `pai_a` is used to calculate the radiation intercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calculated from total plant area index by
#' assuming a realistic shape to the vertical profile foliage within the canopy. If supplied,
#' `pai_a` must have the same dimensions as micro$pai. I.e. with the same x and y
#' dims as the the supplied dtm and values for each hour as the z dimension. The parameter `surfwet`
#' determines how much of the canopy should be treated as wet surface when calculating
#' latent heat fluxes. The units of dtmc must match dtm and must be elevation data if an
#' altitude correction is applied. If `altcorrect`>0, the elevation
#' difference between each pixel of dtm and dtmc is calculated and
#' an elevation lapse rate correction is applied to the temperature and
#' pressure data to account for these elevation differences. If `altcorrect`= 1, a fixed
#' lapse rate of 5 degrees per 100m is applied to the temperature data. If
#' `altcorrect`= 2, humidity-dependent lapse rates are calculated and applied.
#' @import terra
#' @import stats
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @rdname runmicro
#' @examples
#' library(terra)
#' # ** First run point model setting reqhgt to 5 cm above ground using inbuilt datasets
#' micropoint <- runpointmodel(climdata, 0.05, dtmcaerth, vegp, soilc)
#' # ** Subset inbuilt point model to get hottest and coldest days in each month
#' micropoint_mx <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")
#' micropoint_mn <- subsetpointmodel(micropoint, tstep = "month", what = "tmin")
#' # ** Run grid model for hottest and coldest days (takes ~20 seconds to run)
#' mout_mx <- runmicro(micropoint_mx, 0.05, vegp, soilc, dtmcaerth)
#' mout_mn <- runmicro(micropoint_mn, 0.05, vegp, soilc, dtmcaerth)
#' # Plot air temperatures on hottest hour
#' mypal <- colorRampPalette(c("darkblue", "blue", "green", "yellow",
#' "orange", "red"))(255)
#' plot(rast(mout_mx$Tz[,,134]), col = mypal)
#' # Plot mean of monthly max and min
#' mairt<-apply((mout_mn$Tz + mout_mx$Tz) / 2, c(1,2),mean)
#' plot(rast(mairt), col = mypal)
#' # Remove vegetation effects and run again for one cm above ground
#' vegp2 <- vegp
#' vegp2$pai <- rast(vegp2$pai) * 0
#' vegp2$hgt <- rast(vegp2$hgt) * 0
#' mout_mx <- runmicro(micropoint_mx, 0.005, vegp2, soilc, dtmcaerth)
#' plot(rast(mout_mx$Tz[,,134]), col = mypal)
runmicro <- function(micropoint, reqhgt, vegp, soilc, dtm, dtmc = NA, altcorrect = 0,
                     snow = FALSE, snowmod = NA, runchecks = TRUE, pai_a = NA, tfact = 1.5,
                     out = rep(TRUE, 10), slr = NA, apr = NA, hor = NA, twi = NA,
                     wsa = NA, svf = NA, method = "Cpp") {
  if (snow) {
    if (class(micropoint) == "micropoint") { # data.frame climate input
      mout<-.runmicrosnow1(micropoint,reqhgt,vegp,soilc,dtm,snowmod,runchecks,pai_a,
                           tfact,out,slr,apr,hor,twi,wsa,svf)
    } else {  # array climate input
      mout<-.runmicrosnow2(micropoint,reqhgt,vegp,soilc,dtm,dtmc,snowmod,altcorrect,
                           runchecks,pai_a,tfact,out,slr,apr,hor,twi,wsa,svf)
    }
  } else {
    mout<-.runmicronosnow(micropoint,reqhgt,vegp,soilc,dtm,dtmc,altcorrect,runchecks,
                    pai_a,tfact,out,slr,apr,hor,twi,wsa,svf)
  }
  return(mout)
}
#' runmicro on big areas
#'
#' The function `runmicro_big` tiles larger studies and saves outputs for each tile
#' @description The function `runmicro` runs the grid version of the microclimate model
#' @param micropoint an object of class micropoint or a list of objects of class
#' micropoint as returned by [runpointmodel()], [runpointmodela()], [subsetpointmodel()]
#' or [subsetpointmodela()]
#' @param reqhgt  height above (postive) or below (negative) ground for which microclimate variables are required (m)
#' @param pathout a file directory to which to save data. Data saved to a subdirectory
#' called `microut` in this directory. Default is working directory.
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a SpatRaster object of elevations in metres (see details)
#' @param dtmc a SpatRaster object giving the resolution, spatial extent, and projection
#' of the climate data used when running [micropointa()]. Ignored if climate data
#' used for running the point model are provided as a data.frame. Must
#' give elevations in metres if `altcorrect` > 0 or if setting runchecks to TRUE.
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction, see details)
#' @param tilesize optional integer of the number of pixels in x and y of each tile (returned tiles are square).
#' Calculated automatically based on data size of not provided.
#' @param toverlap optional integer specifying the number of pixels of overlap between adjacent tiles. Default 0,
#' Set to > 0 and use [mosaicblend()] if output appears to have tiling effects.
#' @param writeasnc optional logical indicating whether to write output data as netCDF4 files (default TRUE). Can only be used
#' if all entries of `out` are true. Alternatively data are saved as RDS files.
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' @param pai_a an array of plant area index values above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param out optional vector of logicals indicating which variables to
#' return ordered as for the listed outputs when `reqhgt > 0'` (e.g. `out[1] = TRUE` indicates that
#' `Tz` is returned, `out[2]` that `tleaf` is returned etc). By default all variables
#' are returned.
#' @param silent optional logical indicating whether to report on progress (default FALSE - progress reported).
#' @import terra
#' @import ncdf4
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @details if `writeasnc = TRUE` and all entries in `out` are `TRUE`, a ncdf4
#' files for each tile are written out, unless the corresponding tile comprises
#' all NAs. To save memory, data are written out as integers: temperatures and wind
#' speeds are multiplied by 100 prior to doing so. If `writeasnc = FALSE` or
#' not all variables are required, then an RDS file is written out. a dtm (as a
#' rwapped SpatRaster) of the corresponding tile is attached to the model outputs prior to doing so to
#' aid with georeferencing.
#' see [runmicro()]
runmicro_big <- function(micropoint, reqhgt, pathout = getwd(), vegp, soilc, dtm, dtmc = NA, altcorrect = 0,
                         tilesize = NA, toverlap = 0, writeasnc = FALSE, runchecks = TRUE, pai_a = NA, tfact = 1.5,
                         out = rep(TRUE,10)) {
  # ============== Unpack and check data ===================================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  .checkbiginputs(dtm,vegp,soilc) # check no odd NAs
  # Clean variables

  # ========= Determine appropriate tile size if not specified =============== #
  # Calculate tile size
  if (is.na(tilesize) == TRUE) {
    if (class(micropoint) == "micropoint") {
      nt<-length(micropoint$weather$temp)
      tme<-as.POSIXlt(micropoint$weather$obs_time,tz="UTC")
    } else {
      for (i in 1:length(micropoint)) {
        onepoint<-micropoint[[i]]
        if (class(onepoint) != "logical") {
          nt<-length(onepoint$weather$temp)
          tme<-as.POSIXlt(onepoint$weather$obs_time,tz="UTC")
        }
      }
    }
    osize<-sqrt(20000000/nt)-2*toverlap
    sizeo<-c(10,20,50,100,200,500,1000,2000)
    tilesize<-sizeo[which.min(abs(osize-sizeo))]
  }
  dms<-dim(dtm)[1:2]
  rws<-ceiling(dms[1]/tilesize)
  cls<-ceiling(dms[2]/tilesize)
  nn<-rws*cls+5
  pb <- utils::txtProgressBar(min = 0, max = nn, style = 3)
  # ============== Create directory for storing data ========================= #
  # Create directory for storing data
  path2<-paste0(pathout,"microut/")
  dir.create(path2,showWarnings = FALSE)
  # Calculate universal variables
  slr<-terrain(dtm,v="slope")
  apr<-terrain(dtm,v="aspect")
  slr[is.na(slr)]<-0
  apr[is.na(apr)]<-0
  slr<-mask(slr,dtm)
  apr<-mask(apr,dtm)
  utils::setTxtProgressBar(pb, 1)
  twi<-.topidx(dtm)
  utils::setTxtProgressBar(pb, 2)
  hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) hor[,,i]<-.horizon(dtm,(i-1)*15)
  utils::setTxtProgressBar(pb, 3)
  dsm<-dtm+vegp$hgt
  wsa<-.windsheltera(dsm,8,NA)
  utils::setTxtProgressBar(pb, 4)
  msl<-tan(apply(atan(hor),c(1,2),mean))
  svfa<-.rast(0.5*cos(2*msl)+0.5,dtm)
  utils::setTxtProgressBar(pb, 5)
  ctr<-5
  for (rw in 1:rws) {
    for (cl in 1:cls) {
      dtmi<-.croprast(dtm,rw,cl,tilesize,toverlap)
      v<-as.vector(dtmi)
      v<-v[is.na(v)==F]
      if (length(v)>1) {
        slri<-crop(slr,ext(dtmi))
        apri<-crop(apr,ext(dtmi))
        twii<-crop(twi,ext(dtmi))
        sfvi<-crop(svfa,ext(dtmi))
        hori<-as.array(crop(.rast(hor,dtm),ext(dtmi)))
        wsai<-as.array(crop(.rast(wsa,dtm),ext(dtmi)))
        vegpi<-.listcrop(vegp,ext(dtmi))
        soilci<-.listcrop(soilc,ext(dtmi))
        mout<-runmicro(micropoint,reqhgt,vegpi,soilci,dtmi,dtmc,altcorrect,snow=FALSE,snowmod=NA,
                       runchecks,pai_a,tfact,out,slri,apri,hori,twii,wsai,svf=svfi)
        mout$tme<-tme
        # Write output
        rwt<-ifelse(rw<10,paste0("0",rw),paste0("",rw))
        clt<-ifelse(cl<10,paste0("0",cl),paste0("",cl))
        if (writeasnc) {
          sel<-which(out==FALSE)
          if (length(sel) > 0) {
            warning("Can only write as nc with all variables in out set to TRUE. Writing as RDS\n")
            mout$dtm<-wrap(dtmi)
            fo<-paste0(path2,"area_",rwt,"_",clt,".RDS")
            saveRDS(mout,fo)
          }  else {
            fo<-paste0(path2,"area_",rwt,"_",clt,".nc")
            writetonc(mout,fo,dtmi,reqhgt)
          }
        } else {
          mout$dtm<-wrap(dtmi)
          fo<-paste0(path2,"area_",rwt,"_",clt,".RDS")
          saveRDS(mout,fo)
        }
      } # end NA check
      ctr<-ctr+1
      utils::setTxtProgressBar(pb, ctr)
    } # end column
  } # end row
}
#' Generates microclimate equivalent of bioclim variables
#'
#' The function `runbioclim` runs the microclimate model to produce
#' a list of SpatRasters each equivalent to the 19 bioclimate variables
#' produced by www.worldclim.org.
#' @param climdata a data.frame or list of arrays of weather variables as for
#' [runpointmodel()] or [runpointmodela()]. Can be for more than one year.
#' @param reqhgt height (m) for which microclimate variables are needed (<0 if below ground)
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a SpatRaster of elevations for the study area
#' @param dtmc a SpatRaster object giving the resolution, spatial extent, and
#' projection of the climate data used when running [micropointa()]. Ignored if
#' climate data are provided as a data.frame. Must give elevations in metres if
#' `altcorrect` > 0 or if setting runchecks to TRUE.
#' @param tme  POSIXlt object giving the dates and times for each weather variable.
#' Used only if `climdata` provided as arrays.
#' @param temp one of `air` or `leaf` indicating whether outputs represent leaf or
#' air temperatures.
#' @param zref height above ground (m) of temperature measurements in climdata
#' @param windhgt height above ground (m) of wind speed data inclimdata.
#' @param soilm optional vector of soil moisture values in upper 10 cm of the soil (calculated if not supplied)
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' @param altcorrect a single numeric value indicating whether to apply an
#' elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed
#' lapse rate correction, 2 = humidity-dependent variable lapse rate correction,
#' see [runmicro()]).
#' @param pai_a an array of plant area index values above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param out optional vector of logicals indicating which bioclim variables to return.
#' Default all 19.
#' @param vegpisannual optional logical used when vegetation varies temporally to
#' indicate whether vegetation data correspond to values for a full annual cycle.
#' Relevent if `climdata` represent data for more than one year.
#' @return a multilayer SpatRast of the following:
#' \describe{
#'  \item{BIO1}{Mean of monthy median temperatures (degrees C)}
#'  \item{BIO2}{Mean diurnal temperature range (degrees C)}
#'  \item{BIO3}{Isothermality (BIO2/BIO7) (×100)}
#'  \item{BIO4}{Temperature Seasonality (standard deviation of monthly median temperatures × 100)}
#'  \item{BIO5}{Maximum temperature (degrees C)}
#'  \item{BIO6}{Minimum temperature (degrees C)}
#'  \item{BIO7}{Temperature Annual Range (BIO5-BIO6) (degrees C)}
#'  \item{BIO8}{Mean of monthly median temperatures in wettest three months (degrees C)}
#'  \item{BIO9}{Mean of monthly median temperatures in driest three months (degrees C)}
#'  \item{BIO10}{Mean of monthly median temperatures in warmest three months (degrees C)}
#'  \item{BIO11}{Mean of monthly median temperatures in coldest three months (degrees C)}
#'  \item{BIO12}{Mean of monthly soil moistures (m^3 / m^3)}
#'  \item{BIO13}{Wettest soil moisture of days with median temperature (m^3 / m^3)}
#'  \item{BIO14}{Driest soil moisture of days with median temperature(m^3 / m^3)}
#'  \item{BIO15}{Soil moisture seasonality (Coefficient of Variation) on day in each month with median temperature}
#'  \item{BIO16}{Mean soil moisture of wettest three months (m^3 / m^3)}
#'  \item{BIO17}{Mean soil moisture of driest three months (m^3 / m^3)}
#'  \item{BIO18}{Mean soil moisture of warmest three months (m^3 / m^3)}
#'  \item{BIO19}{Mean soil moisture of coldest three months (m^3 / m^3)}
#' }
#' @details
#' To enhance computational efficiency the microclimate model is run for selected days only. Thus, to compute
#' mean annual temperature, the mean ambient temperature of each day in the input weather data is calculated,
#' the day with median temperatures in each month selected and the mean across months calculated. This is
#' not, strictly speaking, the same as the mean temperature, but differences are likely to be minor, and for each year of
#' data supplied, there is an approximately 30-fold gain in computational efficiency by calculating
#' BIO1 in this way. Similarly, to calculate maximum temperature (BIO5), the day of the year with the hottest ambient
#' temperature is selected, and microclimate temperatures calculated on this day only. This ignores
#' the possibility that on a slightly cooler, but sunnier day, microclimate temperatures may be hotter
#' hotter at certain locations. If `hourly = TRUE` all hours within a given day are selected and calculations
#' performed on hourly data. If If `hourly = FALSE` only the hours corresponding to times when hourly temperatures
#' are at their daily maximum and minimum and selected. This results in a c. 10-fold increase in computational
#' efficiency, but cannot pick out areas where terrain results in near-ground temperatures reaching
#' a maximum later in the afternoon than the peak in ambient temperature. If weather data for more than one year are supplied,
#' only one set of median, maximum and minimum monthly temperature data are selected representing an average
#' across years. Resultant, there is little computational penalty if providing data for multiple years
#' in comparison to one year of data.
#' @import stats
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @rdname runbioclim
#' @examples
#' # Run bioclim model with default inputs
#' bioclim <- runbioclim(climdata, 0.05, vegp, soilc, dtm)
#' plot(bioclim[[1]]) # mean annual temperature
runbioclim <- function(climdata, reqhgt, vegp, soilc, dtm, dtmc = NA, tme = NA, temp="air",
           zref = 2, windhgt = zref, soilm = NA, runchecks = TRUE, altcorrect=0,
           pai_a = NA, tfact = 1.5, out = rep(TRUE,19), vegpisannual = TRUE) {
  up<-.unpack(dtm,vegp,soilc)
  dmx<-.vegpdmx(up$vegp)
  if (dmx == 1) {  # temporally invarient vegetation
    if (class(climdata) == "data.frame") { # data.frame climate input
      mout<-.runbioclim1(climdata,reqhgt,vegp,soilc,dtm,temp,zref,windhgt,soilm,runchecks,
                         pai_a,tfact,out)
    } else {  # array climate input
      mout<-.runbioclim2(climdata,tme,reqhgt,vegp,soilc,dtm,dtmc,temp,zref,windhgt,
                         soilm,runchecks,altcorrect,pai_a,tfact,out)
    }
  } else { # time variant vegetation
    if (class(climdata) == "data.frame") { # data.frame climate input
      mout<-.runbioclim3(climdata,reqhgt,vegp,soilc,dtm,temp,zref,windhgt,soilm,runchecks,
                         pai_a,tfact,out,vegpisannual)
    } else { # array climate input
      mout<-.runbioclim4(climdata,tme,reqhgt,vegp,soilc,dtm,dtmc,temp,zrefwindhgt,
                         soilm,runchecks,altcorrect,pai_a,tfact,out,vegpisannual)
    } # end if array
  }  # end if time variant
  return(mout)
}
#' Runs the grid snow model
#'
#' @description The function `runsnowmodel` runs the inbuilt snow model returning
#' snow depth at each time increment
#' @param weather a data.frame or list of arrays of weather variables as provided to [runpointmodel()]
#' or [runpointmodela()] (see details)
#' @param micropoint an object of class micropoint or a list of objects of class
#' micropoint as returned by [runpointmodel()], [runpointmodela()], [subsetpointmodel()]
#' or [subsetpointmodela()]
#' @param reqhgt  height above (postive) or below (negative) ground for which microclimate variables are required (m)
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a SpatRaster object of elevations in metres (see details)
#' @param dtmc a SpatRaster object giving the resolution, spatial extent, and projection
#' of the climate data used when running [micropointa()]. Ignored if climate data
#' used for running the point model are provided as a data.frame. Must
#' give elevations in metres if `altcorrect` > 0 or if setting runchecks to TRUE.
#' @param tme POSIXlt object giving the dates and times for each weather variable stored in the array. Only rquired if `weather` is a list of arrays.
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction)
#' @param snowenv one of `Maritime`, `Prairie`, `Tundra`, `Taiga` or `Alpine`.
#' Used to compute snow density as a function of snow age following Sturm et al (2010) J Hydrometeorol, 11: 1380–94.
#' @param method one of `fast` or `slow` and used when `micropoint` is a subset of values
#' returned by e.g. [subsetpointmodel()]. If `fast` the full snow model is only run over days over which microclimate
#' estimates are required and an approximation method is used to estimate evolving snow depth during days
#' in which the model is not run. If `slow` the full model is run for every hour and the resulting
#' output subset to match `micropoint`.
#' @param snowinitd a single numeric value or matrix of values indicating initial snow depths (m)
#' at the start of the model run
#' @param snowinita a single numeric value or matrix of values indicating initial snow age (hours)
#' at the start of the model run
#' @param zref height above ground (m) of temperature measurements in weather
#' @param windhgt height above ground (m) of wind speed data in weather
#' @param stfact optional parameter indicating sensitivity of spatial snow re-distribution to terrain (0 = no terrain effect).
#' @return
#' \describe{
#'   \item{Tc}{Average temperature of snowpack (deg C)}
#'   \item{Tg}{Temperature of ground snow surface (deg C)}
#'   \item{groundsnowdepth}{Depth of ground-lying snow (m)}
#'   \item{totalSWE}{Total snow water equivelent of snow pack (mm)}
#'   \item{snowden}{Snow density (Kg / m^3)}
#' }
#' @details weather data provided to `runsnowmodel` must comprise hourly data for every hour,
#' even if the input `micropoint` is a subset version of the point model as returned by e.g.
#' [subsetpointmodel()].
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @rdname runsnowmodel
#' @examples
#' climdata$temp <- climdata$temp - 8 # Make it colder so there is snow
#' # Run full snow model for every hour with default settings (takes ~90 seconds)
#' micropoint <- runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc) # Make it colder so there is snow
#' smod <- runsnowmodel(climdata, micropoint, vegp, soilc, dtmcaerth)
#' # Plot mean snow depth through time
#' sdepth <- apply(smod$groundsnowdepth, 3, mean, na.rm = TRUE)
#' plot(sdepth, type = "l")
#' # Plot spatial variation in snow depth when snow is at its deepest
#' n<-which.max(sdepth)
#' plot(rast(smod$groundsnowdepth[,,n]))  # spatial snow depth at maximum
#' # subset point model and run again
#' micropoint <- subsetpointmodel(micropoint)
#' # Run model using method = slow (takes ~90 seconds again as full model run and then subset)
#' smod <- runsnowmodel(climdata, micropoint, vegp, soilc, dtmcaerth, method = "slow")
#' # Run model using method = fast (takes ~4 seconds)
#' smod <- runsnowmodel(climdata, micropoint, vegp, soilc, dtmcaerth, method = "fast")
runsnowmodel<-function(weather, micropoint, vegp, soilc, dtm, dtmc = NA, tme = NA, altcorrect = 0,
                       snowenv="Taiga", method="fast", snowinitd = 0,  snowinita = 0,
                       zref = 2, windhgt = zref, stfact = 0.01) {
  if (class(micropoint) == "micropoint") { # data.frame weather input
    if (length(micropoint$subs) == length(micropoint$tmeorig)) { # no subset required
      smod<-.snowmodel1(micropoint$weather,dtm,vegp,soilc,snowenv,snowinitd,snowinita,micropoint$zref,micropoint$zref,stfact)
    } else { # subset of model required
      if (method == "fast") {
        smod<-.snowmodelq1(weather,dtm,vegp,soilc,micropoint$subs,snowenv,snowinitd,snowinita,zref,windhgt,stfact)
      }  else {  # method = slow
        smod<-.snowmodel1(weather,dtm,vegp,soilc,snowenv,snowinitd,snowinita,zref,windhgt,stfact)
        smod<-subsetsnowmodel(smod,micropoint$subs)
      }  # end fast / slow test
    } # end subset required test
  } else {  # climate as array
    for (i in 1:length(micropoint)) {
      onepoint<-micropoint[[i]]
      if (class(onepoint) != "logical") {
        subs<-onepoint$subs
        tmeorig<-onepoint$tmeorig
      }
    }
    if (length(subs) == length(tmeorig)) { # no subset required
      smod<-.snowmodel2(micropoint$weather,tme,dtm,dtmc,vegp,soilc,altcorrect,snowenv,
                        snowinitd,snowinita,micropoint$zref,micropoint$zref,stfact)
    } else { # subset of model required
      if (method == "fast") {
        smod<-.snowmodelq2(weather,tme,dtm,dtmc,vegp,soilc,subs,altcorrect,snowenv,
                           snowinitd,snowinita,zref,windhgt,stfact)
      } else {# method = slow
        smod<-.snowmodel2(weather,tme,dtm,dtmc,vegp,soilc,altcorrect,snowenv,
                          snowinitd,snowinita,zref,windhgt,stfact)
        smod<-subsetsnowmodel(smod,subs)
      } # end fast / slow test
    } # end subset test
  } # end array / data.frame test
  return(smod)
}

