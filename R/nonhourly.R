#' Runs point microclimate model over arrays
#'
#' The function `runpointmodela` runs the point microclimate model over every cell of an array of climate data
#'
#' @param climarray a list of arrays of weather variables (see details and [nctoclimarray()])
#' @param precarray an array of daily precipitation (see details)
#' @param tme an object of class POSIXlt giving the dates and times for each weather variable stored in the array
#' @param reqhgt height for which temperatures are needed (used only when reqhgt < 0 to calculate tmeperature below ground)
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param windhgt height above ground of wind speed data in weather
#' @param soilm optional array of soil moisture values in upper 10 cm of the soil (calculated if not supplied)
#' @param dTmx optional maximum amount by which canopy or ground surface temperatures can exceed air temperatures.
#' Included to ensure model convergence
#' @param maxiter optional integer indicating the maximum number of iterations (see details)
#' @return a list of point model outputs (as returned by [runpointmodel()] for each grid cell
#' @details The units of `climarray` must follow those in the dataset `climdata`.
#' It must be a list with each component of the list an array, named using the same
#' names as the column headers in weather (e.g. temp for temperature), excluding `obs_time`.
#' As not all wind measurements are at reference height, the height of the wind speed measurement
#' must be specified if not 2 m. To enable calculation of below-canopy wind and temperature profiles
#' in tall canopy, the wind speed and temperature data are  adjusted to give values for a height at
#' the maximum vegetation height if the tallest vegetation exceeds two metres.
#' For doing so a stand vegetation surface typical of that in which a weather station would be located is assumed.
#' The parameter `maxiter` sets the maximum number of times the model is iterated to achieve
#' convergence. Increasing this value improves accuracy at the expense of computation time.
#' @export
runpointmodela<-function(climarray, precarray, tme, reqhgt = 0.05, vegp, soilc, windhgt = 2, soilm = NA, dTmx = 25, maxiter = 10)  {
  .todf<-function(climarray,i,j,tme) {
    climdf<-with(climarray,data.frame(obs_time=tme,
                                      temp=temp[i,j,],
                                      relhum=relhum[i,j,],
                                      pres=pres[i,j,],
                                      swrad=swrad[i,j,],
                                      difrad=difrad[i,j,],
                                      skyem=skyem[i,j,],
                                      windspeed=windspeed[i,j,],
                                      winddir=winddir[i,j,]))
    climdf
  }
  k<-1
  pointo<-list()
  for (i in 1:dim(precarray)[1]) {
    for (j in 1:dim(precarray)[2]) {
      climdf<-.todf(climarray,i,j,tme)
      prec<-precarray[i,j,]
      if (class(soilm) == "logical") {
        soilmo<-NA
      } else soilmo<-soilm[i,j,]
      if (is.na(climdf$temp[1]) == FALSE) {
        pointo[[k]]<-runpointmodel(climdf,prec,reqhgt,vegp,soilc,windhgt,soilmo,dTmx,maxiter)
      } else pointo[[k]]<-NA
      k<-k+1
    }
  }
  return(pointo)
}
#' Create object of class microindaily
#'
#' @description The function `modelin_dy` creates an object of class microindaily
#' which unpacks various component inputs and reformats as required
#' for running the model in daily timesteps
#'
#' @param micropoint an object of class micropoint as returned by [runpointmodel()] or [subsetpointmodel()]
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a PackedSpatRaster or SpatRaster object of elevations (see details)
#' @param windhgt height above ground of wind speed measurement (m) in weather dataset (see details).
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' checks on format and units of input data.
#' @param xyf optional spatial smoothing factor applied in calculation of surface
#' roughness and zero-plane displacement heights (see details)
#' @details see details for modelin.
#'
#' @seealso [checkinputs()], [modelin()], [modelina_dy()]
#'
#' @import terra
#' @export
modelin_dy<-function(micropoint, vegp, soilc, dtm, runchecks = TRUE, xyf = 1) {
  # Extract daily min and max values from micropoint
  ptd<-.subsetpointtoday(micropoint)
  micropoint_mn<-ptd$micropoint_mn
  micropoint_mx<-ptd$micropoint_mx
  # Create model in objects
  micro_mn<-modelin(micropoint_mn,vegp,soilc,dtm,runchecks,xyf)
  micro_mx<-modelin(micropoint_mx,vegp,soilc,dtm,runchecks,xyf)
  out<-list(micro_mn=micro_mn,micro_mx=micro_mx,micropoint=micropoint)
  class(out)<-"microindaily"
  return(out)
}
#' Create object of class microin with weather data as an array
#'
#' @description The function `modelina` creates an object of class microin
#' which unpacks various component inputs and reformats as required
#' for running the model in hourly timesteps. Here it is assumed that the input
#' weather data are as arrays - i.e. variable in space
#' @param micropointa a list of objects of class micropoint as returned by [runpointmodela()] or [subsetpointmodela()]
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a SpatRaster object of elevations in metres (see details)
#' @param dtmc a SpatRaster object giving the resolution, spatial extent, and projection
#' of the weather data used when running [micropointa()]. Must give elevations in metres if
#' `altcorrect` > 0 or if setting runchecks to TRUE.
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction, see details)
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' @param xyf optional spatial smoothing factor applied in calculation of surface
#' roughness and zero-plane displacement heights (see [modelin()])
#' @param precarray an array of daily precipitation (see details)
#' @details The units of dtmc must match dtm and must be elevation data if an
#' altitude correction is applied. If `altcorrect`>0, the elevation
#' difference between each pixel of dtm and dtmc is calculated and
#' an elevation lapse rate correction is applied to the temperature and
#' pressure data to account for these elevation differences. If `altcorrect`= 1, a fixed
#' lapse rate of 5 degrees per 100m is applied to the temperature data. If
#' `altcorrect`= 2, humidity-dependent lapse rates are calculated and applied.
#' See also details for  modelin.
#' @seealso [modelin()], [modelina_dy()], [nctoclimarray()]
#'
#' @import terra
#' @export
#' @examples
#' library(terra)
#' # ====================== Create dummy array datasets ===================== #
#' ta <- function(x, xdim = 5, ydim = 5) {
#'   array(rep(x, each = ydim * xdim), dim = c(ydim, xdim, length(x)))
#' }
#' # ~~ Create list of climate variable arrays ~~ #
#' climarray <- list(temp = ta(climdata$temp),
#'                   relhum = ta(climdata$relhum),
#'                   pres = ta(climdata$pres),
#'                   swrad = ta(climdata$swrad),
#'                   difrad = ta(climdata$difrad),
#'                   skyem = ta(climdata$skyem),
#'                   windspeed = ta(climdata$windspeed),
#'                   winddir = ta(climdata$winddir))
#' # ~~ Create precipitation array ~~ #
#' precarray <- ta(rainfall)
#' # ============= Run point microclimate model for each grid cell ========== #
#' tme <- as.POSIXlt(climdata$obs_time, tz="UTC")
#' micropoinaltcorrectta <- runpointmodela(climarray, precarray, tme, reqhgt = 0.05, vegp, soilc)
#' micropointa<-subsetpointmodela(micropointa, tstep = "month", what = "tmax")
#' # =================== Run model input function  ========================== #
#' dtmc <- aggregate(rast(dtmcaerth), 10) # Coarse resolution dtm
#' micro <- modelina(micropointa, vegp, soilc, dtmcaerth, dtmc)
modelina <- function(micropointa, vegp, soilc, dtm, dtmc, altcorrect = 0, runchecks = TRUE, xyf = 1) {
  .cca<-function(weather,varn,h,r,rfi) {
    a<-array(NA,dim=c(dim(r)[1:2],h))
    k<-1
    for (i in 1:dim(a)[1]) {
      for (j in 1:dim(a)[2]) {
        onew<-weather[[k]]
        s<-which(names(onew)==varn)
        a[i,j,]<-onew[,s]
        k<-k+1
      }
    }
    ro<-.rast(a,r)
    ro<-resample(ro,rfi)
    as.array(ro)
  }
  weather<-list()
  precip<-list()
  for (i in 1:length(micropointa)) {
    onepoint<-micropointa[[i]]
    weather[[i]]<-onepoint$weather
    precip[[i]]<-onepoint$precip
    if (runchecks) {
      rc<-checkinputs(weather[[i]],precip[[i]],vegp,soilc,dtmc,FALSE,tstep=onepoint$tstep)
      weather[[i]]<-rc$weather
      precip[[i]]<-rc$precip
      vegp<-rc$vegp
      soilc<-rc$soilc
    }
  }
  if (class(dtmc)[1] == "PackedSpatRaster") {
    r<-rast(dtmc)
  } else r<-dtmc
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  # Turn input climate variables into arrays
  rfi<-dtm
  tme<-as.POSIXlt(weather[[1]]$obs_time,tz="UTC")
  h<-length(tme)
  tc<-.cca(weather,"temp",h,r,rfi)
  pk<-.cca(weather,"pres",h,r,r)
  relhum<-.cca(weather,"relhum",h,r,rfi)
  estl<-.satvap(tc)
  ea<-estl*relhum/100
  tdew<-.dewpoint(ea,tc)
  # Apply altitudinal correction
  if (altcorrect == 0) {  # No altitudinal correction
    pk<-as.array(resample(.rast(pk,r),rfi))
  } else { # Altitudinal correction applied
    dtmc[is.na(dtmc)]<-0
    psl<-pk/(((293-0.0065*.rta(dtmc,h))/293)^5.26)
    psl<-as.array(resample(.rast(psl,r),rfi))
    pk<-psl*(((293-0.0065*.rta(rfi,h))/293)^5.26)
    dc<-resample(dtmc,rfi)
    elevd<-.rta(dc-dtm,h)
    if (altcorrect==1) {  # Fixed lapse rate
      tcdif<-elevd*(5/1000)
    } else { # Humidity-dependent lapse rate
      lr<-.lapserate(tc,ea,pk)
      tcdif<-lr*elevd
    }
    tc<-tcdif+tc
  }
  difr<-.cca(weather,"difrad",h,r,rfi)
  swrad<-.cca(weather,"swrad",h,r,rfi)
  di<-swrad-difr
  di[di<0]<-0
  jd<-.jday(tme)
  lt<-with(tme,hour+min/60+sec/3600)
  ll<-.latslonsfromr(rfi)
  lats<-ll$lats
  lons<-ll$lons
  sa<-.solalt(.vta(lt,rfi),.rta(lats,h),.rta(lons,h),.vta(jd,rfi))
  ze<-90-sa
  si<-cos(ze*(pi/180))
  si[si<0]<-0
  dirr<-di/si
  dirr[is.na(dirr)]<-0
  dirr[dirr>1352]<-1352
  dp<-difr/swrad
  dp[is.na(dp)]<-0.5; dp[dp<0]<-0; dp[dp>1]<-1
  skyem<-.cca(weather,"skyem",h,r,rfi)
  # Wind speed and direction
  u2<-.cca(weather,"windspeed",h,r,r)
  wd<-.cca(weather,"winddir",h,r,r)
  wu<-u2*cos(wd*pi/180)
  wv<-u2*sin(wd*pi/180)
  wu<-as.array(resample(.rast(wu,r),rfi))
  wv<-as.array(resample(.rast(wv,r),rfi))
  u2<-sqrt(wu^2+wv^2)
  wdir<-(atan2(wv,wu)*180/pi)%%360
  # Turn point model variables into arrays
  ufp<-array(NA,dim=c(dim(r)[1:2],h))
  soilmp<-ufp
  Gp<-ufp
  k<-1
  for (i in 1:dim(Gp)[1]) {
    for (j in 1:dim(Gp)[2]) {
      ufp[i,j,]<-micropointa[[k]]$microp$uf
      soilmp[i,j,]<-micropointa[[k]]$soilm
      Gp[i,j,]<-micropointa[[k]]$microp$G
      k<-k+1
    }
  }
  ufp<-as.array(resample(.rast(ufp,r),rfi))
  soilmp<-as.array(resample(.rast(soilmp,r),rfi))
  Gp<-as.array(resample(.rast(Gp,r),rfi))
  if (class(micropointa[[1]]$Tbz) != "logical")  {
    Tbp<-array(NA,dim=c(dim(r)[1:2],h))
    Tg<-Tbp
    k<-1
    for (i in 1:dim(Tg)[1]) {
      for (j in 1:dim(Tg)[2]) {
        Tg[i,j,]<-micropointa[[k]]$microp$Tg
        Tbp[i,j,]<-micropointa[[k]]$Tbz
        k<-k+1
      }
    }
    Tg<-as.array(resample(.rast(Tg,r),rfi))
    Tbp<-as.array(resample(.rast(Tbp,r),rfi))
  } else {
    Tbp<-NA
    Tg<-NA
  }
  T0p<-array(NA,dim=c(dim(r)[1:2],h))
  k<-1
  for (i in 1:dim(T0p)[1]) {
    for (j in 1:dim(T0p)[2]) {
      T0p[i,j,]<-.point0(micropointa[[k]])
      k<-k+1
    }
  }
  T0p<-as.array(resample(.rast(T0p,r),rfi))
  # Obtain derived variables
  vegx<-.rta(vegp$x,h)
  lref<-.rta(vegp$leafr,h)
  ltra<-.rta(vegp$leaft,h)
  gref<-.rta(soilc$groundr,h)
  # Calculate zero plane displacement height and roughness length
  dzm<-.sortrough(vegp,xyf)
  if (micropointa[[1]]$tstep=="hour") {
    pai<-.unpackpai(vegp$pai,h/24)
    d<-.unpackpai(dzm$d,h/24)
    zm<-.unpackpai(dzm$zm,h/24)
    xx<-rep(c(1:(h/24)),each=24)
    pai<-pai[,,xx]
    d<-d[,,xx]
    zm<-zm[,,xx]
    if (length(vegp$clump) > 1) {
      clump<- .unpackpai(vegp$clump,h/24)
      clump<-clump[,,xx]
    }
  } else {
    pai<-.unpackpai(vegp$pai,h)
    d<-.unpackpai(dzm$d,h)
    zm<-.unpackpai(dzm$zm,h)
    if (length(vegp$clump) > 1) {
      clump<- .unpackpai(vegp$clump,h)
    }
  }
  soilp<-.soilinit(soilc)
  # Set limit
  pai[pai<0.001]<-0.001
  lref[lref<0.001]<-0.001
  lref[lref>0.99]<-0.99
  ltra[ltra<0.001]<-0.001
  ltra[ltra>0.99]<-0.99
  # Calculate wind direction array
  out<-list(# Misc variables:
    tme=tme,zref=micropointa[[1]]$zref,dtm=dtm,lat=mean(lats),long=mean(lons),
    # Weather variables:
    tc=tc,difr=difr,dirr=dirr,skyem=skyem,estl=estl,ea=ea,tdew=tdew,pk=pk,u2=u2,
    # Point model variables:
    ufp=ufp,soilmp=soilmp,Gp=Gp,Tbp=Tbp,Tg=Tg,T0p=T0p,
    vegp_p=micropointa[[1]]$vegp_p,groundp_p=micropointa[[1]]$groundp_p,DDp=micropointa[[1]]$DD,
    # Vegetation parameters:
    pai=pai,vegx=vegx,lref=lref,ltra=ltra,veghgt=vegp$hgt,gsmax=vegp$gsmax,clump=clump,leafd=vegp$leafd,
    # Soil thermal parameters:
    gref=gref,rho=soilp$rho,Vm=soilp$Vm,Vq=soilp$Vq,Mc=soilp$Mc,
    # Soil hydraulic parameters:
    soilb=soilp$soilb,psi_e=soilp$psi_e,Smax=soilp$Smax,Smin=soilp$Smin,
    # Derived variables
    d=d,zm=zm,
    # Climate and timestep variables variables
    winddir=wdir,tstep="hour",progress=0)
  class(out) <-"microin"
  return(out)
}
#' Create object of class microindaily with weather data as an array
#'
#' @description The function `modelina_dy` creates an object of class microindaily
#' which unpacks various component inputs and reformats as required
#' for running the model in hourly timesteps. Here it is assumed that the input
#' weather data are as arrays - i.e. variable in space
#' @param micropointa a list of objects of class micropoint as returned by [runpointmodela()] or [subsetpointmodela()]
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a SpatRaster object of elevations in metres (see details)
#' @param dtmc a SpatRaster object giving the resolution, spatial extent, and projection
#' of the weather data used when running [micropointa()]. Must give elevations in metres if
#' `altcorrect` > 0.
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction, see details)
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' @param xyf optional spatial smoothing factor applied in calculation of surface
#' roughness and zero-plane displacement heights (see [modelin()])
#' @details see details of modelin and modelina.
#' @seealso [modelin()], [modelina()], [nctoclimarray()]
#'
#' @import terra
#' @export
modelina_dy<-function(micropointa, vegp, soilc, dtm, dtmc, altcorrect = 0, runchecks = TRUE, xyf = 1) {
  micropointa_mn<-list()
  micropointa_mx<-list()
  n<-length(micropointa)
  for (i in 1:n) {
    ptd<-.subsetpointtoday(micropointa[[i]])
    micropointa_mn[[i]]<-ptd$micropoint_mn
    micropointa_mx[[i]]<-ptd$micropoint_mx
  }
  # Create model in objects
  micro_mn<-modelina(micropointa_mn,vegp,soilc,dtm,dtmc,altcorrect,runchecks,xyf)
  micro_mx<-modelina(micropointa_mx,vegp,soilc,dtm,dtmc,altcorrect,runchecks,xyf)
  out<-list(micro_mn=micro_mn,micro_mx=micro_mx,micropoint=micropointa[[trunc(n/2)]])
  class(out)<-"microindaily"
  return(out)
}
#' Run microclimate model (daily)
#'
#' @description The function `runmicro_dy` runs the microclimate model in daily
#' time increments, with the option to expand to hourly
#'
#' @param microd object of class microindaily as returned by [modelin_dy()]
#' @param reqhgt height above ground at which model outputs are needed (m).
#' @param expand optional logical indicating whether to expand daily values to hourly (see details).
#' @param pai_a an optional array of plant area index values above `reqhgt` (see details)
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param surfwet an optional single numeric value of array of values specifying the proportion
#' of the canopy surface that should be treated as wet surface (modelled if not supplied)
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied, but outer cells will be NA.
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from dtm if not supplied, but outer cells will be NA.
#' @return if expand = TRUE, an object of class microout with the following components:
#' @return `Tz` Array of air temperatures at height `reqhgt` (deg C). Identical to `T0`
#' if `reqhgt` = 0.
#' @return `tleaf` Array of leaf temperatures at height `reqhgt` (deg C).
#' NA if `reqhgt` greater than canopy height or `reqhgt` <= 0.
#' @return `T0` Array of ground surface temperatures (deg C)
#' @return `relhum` Array of relative humidities at height `reqhgt` (percentage).
#' NA if `reqhgt` <= 0.
#' @return `windspeed` Array of wind speeds at height `reqhgt` (m/s).
#' NA if `reqhgt` <= 0.
#' @return `Rdirdown` Array of downward direct shortwave radiation incident on
#' horizontal surface (W/m^2)
#' @return `Rdifdown` Array of downward diffuse shortwave radiation incident on
#' horizontal surface (W/m^2)
#' @return `Rlwdown` Array of downward longwave radiation incident on horizontal
#' surface (W/m^2)
#' @return `Rswup` Array of upward shortwave radiation (assumed diffuse) incident
#' on underside of horizontal surface (W/m^2)
#' @return `Rlwup` Array of upward longwave radiation incident on underside of
#' horizontal surface (W/m^2)
#' @return if expand = FALSE, an object of class microutdaily, list with the following components:
#' @return mout_mn an object of class microut for minimum daily temperatures
#' @return mout_mx an object of class microut for maximum daily temperatures
#' @seealso [runmicro_hr()] for running microclimate model in hourly time-steps and
#' [runmicro_big()] for running the microclimate model over large areas as tiles.
#'
#' @details
#' If expand = TRUE, daily minima and maxima are expanded to hourly using values in the hourly
#' weather dataset. See also details for [runmicro_hr()].
#'
#' @import terra
#' @export
runmicro_dy <- function(micro_dy, reqhgt, expand = TRUE, pai_a = NA, tfact = 1.5, surfwet = NA,
                        slr = NA, apr = NA, hor = NA, twi = NA, wsa = NA) {
  # Calculate soil surface temperature and soil moisture
  microd<-soiltemp_dy(micro_dy,reqhgt,pai_a,tfact,slr,apr,hor,twi,wsa)
  # Run above ground
  if (reqhgt > 0) {
    mout_mn<-temphumE(microd$micro_mn,reqhgt,pai_a,tfact,surfwet,slr,apr,hor,twi,wsa)
    mout_mx<-temphumE(microd$micro_mx,reqhgt,pai_a,tfact,surfwet,slr,apr,hor,twi,wsa)
  }
  if (reqhgt == 0) {
    mout_mn<-with(microd$micro_mn,list(Tz=T0,tleaf=NA,T0=T0,soilm=soild,
                                       relhum=NA,windspeed=NA,Rdirdown=radGdir,Rdifdown=radGdif,
                                       Rlwdown=radGlw,Rswup=(1-gref)*radGsw,Rlwup=0.97*5.67*10^-8*(T0+273.15)^4))
    mout_mx<-with(microd$micro_mx,list(Tz=T0,tleaf=NA,T0=T0,soilm=soild,
                                       relhum=NA,windspeed=NA,Rdirdown=radGdir,Rdifdown=radGdif,
                                       Rlwdown=radGlw,Rswup=(1-gref)*radGsw,Rlwup=0.97*5.67*10^-8*(T0+273.15)^4))
  }
  # Run below ground
  if (reqhgt < 0) {
    microd<-below_dy(microd,reqhgt,pai_a,tfact,slr,apr,hor,twi,wsa)
    mout_mn<-with(microd$micro_mn,list(Tz=Tz,tleaf=NA,T0=T0,soilm=soild,
                                         relhum=NA,windspeed=NA,Rdirdown=NA,Rdifdown=NA,Rlwdown=NA,
                                         Rswup=NA,Rlwup=NA))
    mout_mx<-with(microd$micro_mx,list(Tz=Tz,tleaf=NA,T0=T0,soilm=soild,
                                         relhum=NA,windspeed=NA,Rdirdown=NA,Rdifdown=NA,Rlwdown=NA,
                                         Rswup=NA,Rlwup=NA))
  }
  if (expand) {
    mout<-.expandclim(mout_mn,mout_mx,microd$micropoint,reqhgt)
    mout$tme<-as.POSIXlt(microd$micropoint$weather$obs_time,tz="UTC")
    class(mout)<-"microout"
  } else {
    mout_mn$tme<-microd$micro_mn$tme
    mout_mx$tme<-microd$micro_mx$tme
    mout<-list(mout_mn=mout_mn,mout_mx=mout_mx)
    class(mout)<-"microoutdaily"
  }
  return(mout)
}
#' runmicro on big areas
#'
#' The function `runmicro_big` tiles larger studies and saves outputs for each tile
#'
#'
#' @param weather a data.frame of weather variables (see details)
#' @param precip a vector of daily precipitation
#' @param reqhgt height for which temperatures are needed (m, negative if below ground surface)
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a PackedSpatRaster or SpatRaster object of elevations (see details)
#' @param pathout a file directory to which to save data. Default saves to subfolder 'microut' in wording directory.
#' @param hourly optional logical indicating whether to run the model in hourly or daily mode (see details).
#' @param expand optional logical indicating whether to expand model outputs to hourly if run in daily mode
#' (ignored if `hourly = TRUE`).
#' @param tilesize optional intiger of the number of pixels in x and y of each tile (returned tiles are square).
#' Calculated automatically based on data size of not provided.
#' @param writeasnc optional logical indicating whether to write output data as netCDF4 files (default TRUE)
#' @param subsetmodel optional logical indicating whether to subset the model to return e.g. monthly values.
#' See also [subsetpointmodel()]
#' @param tstep one of `year` or `month` or `bioclim`. Used only if `subsetmodel == TRUE` (see [subsetpointmodel()])
#' @param what one of `tmax`, `tmin` or `tmedian`. Used only if `subsetmodel == TRUE` (see [subsetpointmodel()])
#' @param days optionally a vector of the days in the time sequence to return data for.
#' Used only if `subsetmodel == TRUE`. If provided `tstep` is ignored. See [subsetpointmodel()]).
#' @param windhgt height above ground of wind speed data in weather (m) (see details)
#' @param soilm optional vector of soil moisture values in upper 10 cm of the soil (calculated if not supplied)
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run checks on format and units of input data.
#' @param xyf optional spatial smoothing factor applied in calculation of surface roughness and zero-plane displacement
#' heights (see details).
#' @param pai_a an optional array of plant area index values above `reqhgt` (see details)
#' @param tfact coefficient determining sensitivity of soil moisture to variation in topographic wetness (see [soilmdistribute()])
#' @param surfwet an optional single numeric value of array of values specifying the proportion
#' of the canopy surface that should be treated as wet surface (see details)
#' @param silent optional logical indicating whether to report on progress (default FALSE - progress reported).
#' @seealso [runmicro_biga()] for running the microclimate model over large areas
#' as tiles with climate input data provided as arrays.
#' @return  if `hourly = TRUE` or `expand = TRUE` and `writeasnc = TRUE` the function writes writes a seperate netcdf4 file to disk
#' for each tile containing georeferenced microclimate data as returned by [runmicro_hr]. If `hourly = FALSE` and `expand = FALSE`
#' seperate netcdf4 files are written for the hours correspondiing to maximum and minimum temperature within each day. If
#' `writeasnc = FALSE` objects of class 'microout' or `microoutdaily` are written to file as RDS files for each tile with
#' an additional wrapped SpatRast of 'dtm' (see [terra::wrap()]) cropped to corresponding tile added to the output to enable subsequent
#' georeferencing of data. Objects of class 'microout' or `microoutdaily` are both lists.
#' @import terra ncdf4
#' @export
#'
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
#' The SpatRaster datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined. In subsequent
#' downscaling of wind, the drag effects of vegetation, determined by vegetation height
#' and foliage area are accounted for can calculated at this stage. In so doing, it is
#' necessary to accommodate the possibility that the wind speed is not just affected
#' by the surface roughness in each pixel,  but also by vegetation surrounding the
#' location. This is accommodated for by applying `xyf` which effectively smooths
#' the surface roughness coefficients using `terra::aggregate` where `xyf` is the
#' aggregation factor. If `xyf` is set to `NA`, the roughness coefficients are averaged
#' across the study area.  `pai_a` is used to calculate the radiation intercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calculated from total plant area index by
#' assuming a realistic shape to the vertical profile foliage within the canopy. If supplied, `pai_a` must have the
#' same dimensions as micro$pai. I.e. with the same x and y dims as the the
#' supplied dtm and values for each hour as the z dimension. The parameter `surfwet`
#' determines how much of the canopy should be treated as wet surface when calculating
#' latent heat fluxes. However, except when extremely droughted, the matric potential of leaves
#' is such that `surfwet` ~ 1. If set to NA, surface wetness is modelled.
runmicro_big <- function(weather, precip, reqhgt, vegp, soilc, dtm, pathout,
                         hourly = TRUE, expand = FALSE, tilesize = NA, writeasnc = TRUE,
                         subsetmodel = TRUE, tstep = "month", what = "tmax", days = NA,
                         windhgt = 2, soilm = NA, runchecks = TRUE,
                         xyf = 1, pai_a = NA, tfact = 1.5, surfwet = NA, silent = FALSE) {
  # Unpack data
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  .checkbiginputs(dtm,vegp,soilc) # check no odd NAs
  # Run point model
  if (silent == FALSE) cat("Running point microclimate model\n")
  ll<-.latlongfromraster(dtm)
  micropoint<-runpointmodel(weather,precip,reqhgt,vegp,soilc,windhgt,soilm,lat=ll$lat,long=ll$long)
  if (subsetmodel)   micropoint<-subsetpointmodel(micropoint,tstep,what,days)
  # Create directory for storing data
  path2<-paste0(pathout,"microut/")
  dir.create(path2,showWarnings = FALSE)
  # Calculate tile size
  if (is.na(tilesize) == TRUE) {
    nt<-length(micropoint$weather$temp)
    n<-nt*dim(dtm)[1]*dim(dtm)[2]
    reqt<-n/2000000   # Number of tiles needed
    if (hourly==FALSE & expand==FALSE) reqt<-reqt/12
    osize<-sqrt((dim(dtm)[1]*dim(dtm)[2])/reqt) # Optimal size
    sizeo<-c(10,20,50,100,200,500,1000,2000)
    tilesize<-sizeo[which.min(abs(osize-sizeo))]
  }
  # Calculate universal variables
  if (silent == FALSE) cat("Computing terrain variables over whole area\n")
  slr<-terrain(dtm,v="slope",unit="radians")
  apr<-terrain(dtm,v="aspect",unit="radians")
  twi<-.topidx(dtm)
  hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) hor[,,i]<-.horizon(dtm,(i-1)*15)
  dsm<-dtm+vegp$hgt
  wsa<-.windsheltera(dsm,8,NA)
  dms<-dim(dtm)[1:2]
  rws<-ceiling(dms[1]/tilesize)
  cls<-ceiling(dms[2]/tilesize)
  if (silent == FALSE) cat(paste0("Running model over ",rws*cls," tiles\n"))
  for (rw in 1:rws) {
    for (cl in 1:cls) {
      dtmi<-.croprast(dtm,rw,cl,tilesize)
      v<-as.vector(dtmi)
      v<-v[is.na(v)==F]
      if (length(v)>1) {
        slri<-crop(slr,ext(dtmi))
        apri<-crop(apr,ext(dtmi))
        twii<-crop(twi,ext(dtmi))
        hori<-as.array(crop(.rast(hor,dtm),ext(dtmi)))
        wsai<-as.array(crop(.rast(wsa,dtm),ext(dtmi)))
        vegpi<-.vegcrop(vegp,dtmi)
        soilci<-.soilcrop(soilc,dtmi)
        if (silent==FALSE) cat(paste0("Running model for tile ",rw," ",cl,"\n"))
        if (hourly) {
          micro<-modelin(micropoint,vegpi,soilci,dtmi,runchecks,xyf)
          mout<-runmicro_hr(micro,reqhgt,pai_a,tfact,surfwet,slri,apri,hori,twii,wsai)
        } else {
          microd<-modelin_dy(micropoint,vegpi,soilci,dtmi,runchecks,xyf)
          mout<-runmicro_dy(microd,reqhgt,expand,pai_a,tfact,surfwet,slri,apri,hori,twii,wsai)
        }
        # Write output
        rwt<-ifelse(rw<10,paste0("0",rw),paste0("",rw))
        clt<-ifelse(cl<10,paste0("0",cl),paste0("",cl))
        if (silent == FALSE) cat(paste0("Writing model outputs for tile ",rw," ",cl,"\n"))
        if (writeasnc) {
          if (class(mout) == "microoutdaily") {
            fo1<-paste0(path2,"area_",rwt,"_",clt,"_min.nc")
            fo2<-paste0(path2,"area_",rwt,"_",clt,"_max.nc")
            writetonc(mout$mout_mn,fo1,dtmi,reqhgt)
            writetonc(mout$mout_mx,fo2,dtmi,reqhgt)
          } else {
            fo<-paste0(path2,"area_",rwt,"_",clt,".nc")
            writetonc(mout,fo,dtmi,reqhgt)
          }
        } else {
          mout$dtm<-wrap(dtmi)
          fo<-paste0(path2,"area_",rwt,"_",clt,".RDS")
          saveRDS(mout,fo)
        }
      } # end NA check
    } # end column
  } # end row
}
#' runmicro on big areas with climate data provided as arrays
#'
#' The function `runmicro_biga` tiles larger studies and saves outputs for each tile
#'
#' @param climarray a data.frame of weather variables (see details)
#' @param precarray an array of daily precipitation (see details)
#' @param tme an object of class POSIXlt giving the dates and times for each weather variable stored
#' in the array. Must be in UTC timezone.
#' @param reqhgt height for which temperatures are needed (m, negative if below ground surface)
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a PackedSpatRaster or SpatRaster object of elevations (see details)
#' @param dtmc a SpatRaster object giving the resolution, spatial extent, and projection
#' of the weather data used when running [micropointa()]. Must give elevations in metres if
#' `altcorrect` > 0.
#' @param pathout a file directory to which to save data. Default saves to subfolder 'microut' in wording directory.
#' @param hourly optional logical indicating whether to run the model in hourly or daily mode (see details).
#' @param expand optional logical indicating whether to expand model outputs to hourly if run in daily mode
#' (ignored if `hourly = TRUE`).
#' @param tilesize optional intiger of the number of pixels in x and y of each tile (returned tiles are square).
#' Calculated automatically based on data size of not provided.
#' @param writeasnc optional logical indicating whether to write output data as netCDF4 files (default TRUE)
#' @param subsetmodel optional logical indicating whether to subset the model to return e.g. monthly values.
#' See also [subsetpointmodel()]
#' @param tstep one of `year` or `month` or `bioclim`. Used only if `subsetmodel == TRUE` (see [subsetpointmodel()])
#' @param what one of `tmax`, `tmin` or `tmedian`. Used only if `subsetmodel == TRUE` (see [subsetpointmodel()])
#' @param days optionally a vector of the days in the time sequence to return data for.
#' Used only if `subsetmodel == TRUE`. If provided `tstep` is ignored. See [subsetpointmodel()]).
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to
#' temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable
#' lapse rate correction, see details)
#' @param windhgt height above ground of wind speed data in weather (m) (see details)
#' @param soilm optional vector of soil moisture values in upper 10 cm of the soil (calculated if not supplied)
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run checks on format and units of input data.
#' @param xyf optional spatial smoothing factor applied in calculation of surface roughness and zero-plane displacement
#' heights (see details)
#' @param pai_a an optional array of plant area index values above `reqhgt` (see details)
#' @param tfact coefficient determining sensitivity of soil moisture to variation in topographic wetness (see [soilmdistribute()])
#' @param surfwet an optional single numeric value of array of values specifying the proportion
#' of the canopy surface that should be treated as wet surface (see details)
#' @param silent optional logical indicating whether to report on progress (default FALSE - progress reported).
#' @seealso [runmicro_biga()] for running the microclimate model over large areas
#' as tiles with climate input data provided as a data.frame.
#' @return  if `hourly = TRUE` or `expand = TRUE` and `writeasnc = TRUE` the function writes writes a seperate netcdf4 file to disk
#' for each tile containing georeferenced microclimate data as returned by [runmicro_hr]. If `hourly = FALSE` and `expand = FALSE`
#' seperate netcdf4 files are written for the hours correspondiing to maximum and minimum temperature within each day. If
#' `writeasnc = FALSE` objects of class 'microout' or `microoutdaily` are written to file as RDS files for each tile with
#' an additional wrapped SpatRast of 'dtm' (see [terra::wrap()]) cropped to corresponding tile added to the output to enable subsequent
#' georeferencing of data. Objects of class 'microout' or `microoutdaily` are both lists.
#' @import terra ncdf4
#' @export
#'
#' @details
#' The units of `climarray` must follow those in the dataset `climdata`.
#' It must be a list with each component of the list an array, named using the same
#' names as the column headers in weather (e.g. temp for temperature), excluding `obs_time`.
#' As not all wind measurements are at reference height, the height of the wind speed measurement
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
#' The SpatRaster datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined. In subsequent
#' downscaling of wind, the drag effects of vegetation, determined by vegetation height
#' and foliage area are accounted for can calculated at this stage. In so doing, it is
#' necessary to accommodate the possibility that the wind speed is not just affected
#' by the surface roughness in each pixel,  but also by vegetation surrounding the
#' location. This is accommodated for by applying `xyf` which effectively smooths
#' the surface roughness coefficients using `terra::aggregate` where `xyf` is the
#' aggregation factor. If `xyf` is set to `NA`, the roughness coefficients are averaged
#' across the study area.  `pai_a` is used to calculate the radiation intercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calculated from total plant area index by
#' assuming a realistic shape to the vertical profile foliage within the canopy. If supplied, `pai_a` must have the
#' same dimensions as micro$pai. I.e. with the same x and y dims as the the
#' supplied dtm and values for each hour as the z dimension. The parameter `surfwet`
#' determines how much of the canopy should be treated as wet surface when calculating
#' latent heat fluxes. However, except when extremely droughted, the matric potential of leaves
#' is such that `surfwet` ~ 1. If set to NA, surfess wetness is modelled.
#' The units of dtmc must match dtm and must be elevation data if an
#' altitude correction is applied. If `altcorrect`>0, the elevation
#' difference between each pixel of dtm and dtmc is calculated and
#' an elevation lapse rate correction is applied to the temperature and
#' pressure data to account for these elevation differences. If `altcorrect`= 1, a fixed
#' lapse rate of 5 degrees per 100m is applied to the temperature data. If
#' `altcorrect`= 2, humidity-dependent lapse rates are calculated and applied.
runmicro_biga <- function(climarray, precarray, tme, reqhgt, vegp, soilc, dtm, dtmc,
                          pathout = getwd(), hourly = FALSE, expand = FALSE, tilesize = NA, writeasnc = TRUE,
                          subsetmodel = TRUE, tstep = "month", what = "tmax", days = NA,
                          altcorrect = 1, windhgt = 2, soilm = NA, runchecks = TRUE,
                          xyf = 1, pai_a = NA, tfact = 1.5, surfwet = NA, silent = FALSE) {
  # Run point model
  if (silent == FALSE) cat("Running point microclimate model over all grid cells of climate array\n")
  micropointa<-runpointmodela(climarray,precarray,tme,reqhgt,vegp,soilc,windhgt,soilm,maxiter=10)
  if (subsetmodel) micropointa<-subsetpointmodela(micropointa,tstep,what,days)
  # Create directory for storing data
  path2<-paste0(pathout,"microut/")
  dir.create(path2,showWarnings = FALSE)
  # Unpack data
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  # Calculate tile size
  if (is.na(tilesize) == TRUE) {
    nt<-length(micropointa[[1]]$weather$temp)
    n<-nt*dim(dtm)[1]*dim(dtm)[2]
    reqt<-n/2000000   # Number of tiles needed
    if (hourly==FALSE & expand==FALSE) reqt<-reqt/12
    osize<-sqrt((dim(dtm)[1]*dim(dtm)[2])/reqt) # Optimal size
    sizeo<-c(10,20,50,100,200,500,1000,2000)
    tilesize<-sizeo[which.min(abs(osize-sizeo))]
  }
  # Calculate universal variables
  if (silent == FALSE) cat("Computing terrain variables over whole area\n")
  slr<-terrain(dtm,v="slope",unit="radians")
  apr<-terrain(dtm,v="aspect",unit="radians")
  twi<-.topidx(dtm)
  hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) hor[,,i]<-.horizon(dtm,(i-1)*15)
  dsm<-dtm+vegp$hgt
  wsa<-.windsheltera(dsm,8,NA)
  dms<-dim(dtm)[1:2]
  rws<-ceiling(dms[1]/tilesize)
  cls<-ceiling(dms[2]/tilesize)
  if (silent == FALSE) cat(paste0("Running model over ",rws*cls," tiles\n"))
  for (rw in 1:rws) {
    for (cl in 1:cls) {
      dtmi<-.croprast(dtm,rw,cl,tilesize)
      v<-as.vector(dtmi)
      v<-v[is.na(v)==F]
      if (length(v)>1) {
        slri<-crop(slr,ext(dtmi))
        apri<-crop(apr,ext(dtmi))
        twii<-crop(twi,ext(dtmi))
        hori<-as.array(crop(.rast(hor,dtm),ext(dtmi)))
        wsai<-as.array(crop(.rast(wsa,dtm),ext(dtmi)))
        vegpi<-.vegcrop(vegp,dtmi)
        soilci<-.soilcrop(soilc,dtmi)
        if (silent==FALSE) cat(paste0("Running model for tile ",rw," ",cl,"\n"))
        if (hourly) {
          micro<-modelina(micropointa,vegpi,soilci,dtmi,dtmc,altcorrect,runchecks,xyf)
          mout<-runmicro_hr(micro,reqhgt,pai_a,tfact,surfwet,slri,apri,hori,twii,wsai)
        } else {
          microd<-modelina_dy(micropointa,vegpi,soilci,dtmi,dtmc,altcorrect,runchecks,xyf)
          mout<-runmicro_dy(microd,reqhgt,expand,pai_a,tfact,surfwet,slri,apri,hori,twii,wsai)
        }
        # Write output
        rwt<-ifelse(rw<10,paste0("0",rw),paste0("",rw))
        clt<-ifelse(cl<10,paste0("0",cl),paste0("",cl))
        if (silent == FALSE) cat(paste0("Writing model outputs for tile ",rw," ",cl,"\n"))
        if (writeasnc) {
          if (class(mout) == "microoutdaily") {
            fo1<-paste0(path2,"area_",rwt,"_",clt,"_min.nc")
            fo2<-paste0(path2,"area_",rwt,"_",clt,"_max.nc")
            writetonc(mout$mout_mn,fo1,dtmi,reqhgt)
            writetonc(mout$mout_mx,fo2,dtmi,reqhgt)
          } else {
            fo<-paste0(path2,"area_",rwt,"_",clt,".nc")
            writetonc(mout,fo,dtmi,reqhgt)
          }
        } else {
          mout$dtm<-wrap(dtmi)
          fo<-paste0(path2,"area_",rwt,"_",clt,".RDS")
          saveRDS(mout,fo)
        }
      } # end NA check
    } # end column
  } # end row
}
