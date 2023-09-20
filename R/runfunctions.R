#' Runs point microclimate model
#'
#' The function `runpointmodel` runs the point microclimate model
#'
#' @param weather a data.frame of weather variables (see details)
#' @param precip a vector of daily precipitation
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param windhgt height above ground of wind speed data in weather
#' @param reqhgt height for which temperatures are needed (used only when reqhgt < 0 to calculate tmeperature below ground)
#' @param soilm optional vector of soil moisture values in upper 10 cm of the soil (calculated if not supplied)
#' @param dTmx optional maximum amount by which canopy or ground surface temperatures can exceed air temperatures.
#' Included to ensure model convergence
#' @param maxiter optional integer indicating the maximum number of iterations (see details)
#' @return a list of the following:
#' (1) weather - a data.frame of weather variables, but with temperature and
#' wind speed height-adjusted to be above canopy if necessary (see details)
#' (2) precip - a vector of daily precipitation (same as input)
#' (3) microp - an object of class pointmicro as returned bgy [microiter::RunBigLeaf()].
#' (4) soilm - a vector of hourly soil moisture fractions in the upper 10 cm of the soil
#' (5) Tbz - a vector of below ground temperatures. NA if reqhgt > 0
#' (5) vegp - a list of vegetation paremeters used by the point model
#' (6) groundp - a list of ground paremeters used by the point model
#' (7) soiltype - soil type assumed when running the point model
#' (7) lat - latitude (decimal degrees) of study area centre. Determined from input raster
#' (8) long - longitude (decimal degrees) of study area centre. Determined from input raster
#  (9) zref - height to which input weather adjusted in metres (see details)
#' (10) tstep - here set to `hour` (see [subsetpointmodel()])
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
#' @export
#' @examples
#' # example code
#' # Run model:
#' micropoint<-runpointmodel(climdata,rainfall,vegp,soilc)
#' # Plot canopy heat exchange surface temperature
#' microp<-micropoint$microp
#' plot(microp$Tc,type="l") # temperature of canopy surface
runpointmodel<-function(weather, precip, vegp, soilc, windhgt = 2, reqhgt = 0.1, soilm = NA, dTmx = 25, maxiter = 100)  {
  # Convert weather data format
  w2<-weather
  w2$lwdown<-weather$skyem*5.67*10^-8*(weather$temp+273.15)^4
  w2$precip<-rep(precip,each=24)/24
  w2$swdown<-weather$swrad
  w2$swrad<-NULL
  w2$skyem<-NULL
  # Extract vegetation properites
  vegp_p<-list(h=.mfr(vegp$hgt),
               pai=mean(vegp$pai,na.rm=T),
               x=.mfr(vegp$x),
               clump=mean(vegp$clump,na.rm=T),
               lref=.mfr(vegp$leafr),
               ltra=.mfr(vegp$leaft),
               leafd=.mfr(vegp$leafd),
               em=0.97,
               gsmax=.mfr(vegp$gsmax),
               q50=100)
  # Extract modal soil type
  if (class(soilc$soiltype)[1] == "PackedSpatRaster") soilc$soiltype<-rast(soilc$soiltype)
  sn<-.getmode(as.vector(soilc$soiltype))
  soiltype<-soilparameters$Soil.type[sn]
  # Create list of soil parameters
  groundp_p<-list(gref=.mfr(soilc$groundr),
                  slope=0,aspect=180,em=0.97,
                  rho=soilparameters$rho[sn],
                  Vm=soilparameters$Vm[sn],
                  Vq=soilparameters$Vq[sn],
                  Mc=soilparameters$Mc[sn],
                  b=soilparameters$b[sn],
                  Psie=soilparameters$psi_e[sn],
                  Smax=soilparameters$Smax[sn],
                  Smin=soilparameters$Smin[sn],
                  alpha=microiter::soilparams$Smin[sn],
                  n=microiter::soilparams$n[sn],
                  Ksat=soilparameters$Ksat[sn])
  # Peform weather height adjustment
  mxhgt<-.mfr(vegp$hgt,max)
  zref<-ifelse(mxhgt>2,mxhgt,2)
  zref<-max(zref,windhgt)
  if (class(vegp$x)[1] == "PackedSpatRaster") vegp$x<-rast(vegp$x)
  ll<-.latlongfromraster(vegp$x)
  if (zref > 2) {
    # Extract latitude and longitude
    w2<-microiter::weather_hgt(w2,zin=2,uzin=windhgt,zout=zref,ll$lat,ll$long)
  }
  # Set minimum wind spped to avoid convergence issues
  w2$windspeed[w2$windspeed<0.5]<-0.5
  # Run soil moisture model if not provided
  if (class(soilm)=="logical") soilm<-microiter::soilmmodel(w2, soiltype)
  # Run big Leaf model
  microp <- microiter::RunBigLeaf(w2,vegp_p,groundp_p,soilm,ll$lat,ll$long,dTmx,zref,maxiter=maxiter)
  # Sort out weather
  w2$skyem<-weather$skyem
  w2$swrad<-w2$swdown
  w2$swdown<-NULL
  w2$lwdown<-NULL
  w2$precip<-NULL
  if (reqhgt < 0) {
    # Calculate soil diffusivity
    cs<-with(groundp_p,2400*rho/2.64+4180*soilm) # specific heat of soil in J/kg/K
    ph<-with(groundp_p,(rho*(1-soilm)+soilm)*1000)   # bulk density in kg/m3
    frs<-with(groundp_p,Vm+Vq)
    c1<-with(groundp_p,(0.57+1.73*Vq+0.93*Vm)/(1-0.74*Vq-0.49*Vm)-2.8*frs*(1-frs))
    c2<-1.06*groundp_p$rho*soilm
    c3<-1+2.6*groundp_p$Mc^-0.5
    c4<-0.03+0.7*frs^2
    k<-c1+c2*soilm-(c1-c4)*exp(-(c3*soilm)^4) # Thermal conductivity   W/m/K
    ka<-k/(cs*ph)
    # Calculate damping depth
    omdy<-(2*pi)/(24*3600)
    DD<-sqrt(2*ka/omdy)
    # Calculate n
    n<- -118.35*reqhgt/DD
    nmn<-floor(min(n))
    nmx<-ceiling(max(n))
    # Calculate temperatures
    Tnmn<-.ma(microp$Tg,nmn)
    Tnmx<-.ma(microp$Tg,nmx)
    # Calculate wgt (from Tnmx)
    wgt<-(n-nmn)/(nmx-nmn)
    Tb<-wgt*Tnmx+(1-wgt)*Tnmn
    # Calculate tmean weight
    wgt<- 0.041596*(reqhgt/mean(DD))+0.87142
    Tbz<-wgt*Tb+(1-wgt)*mean(microp$Tg)
  } else {
    Tbz<-NA
    DD<-NA
  }
  # return outputs
  return(list(weather=w2,precip=precip,microp=microp,soilm=soilm,Tbz=Tbz,DD=DD,
              vegp_p=vegp_p,groundp_p=groundp_p,soiltype=soiltype,
              lat=ll$lat,long=ll$long,zref=zref,tstep="hour"))
}
#' Create object of class microin with weather data as data.frame
#'
#' @description The function `modelin` creates an object of class microin
#' which unpacks various component inputs and reformats as required
#' for running the model in hourly timesteps. Here it is assumed that the input
#' weather data are a data.frame - i.e. not spatially variable.
#'
#' @param micropoint an object of class micropoint as returned by [runpointmodel()] or [subsetpointmodel()]
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a PackedSpatRaster or SpatRaster object of elevations (see details)
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' checks on format and units of input data.
#' @param xyf optional spatial smoothing factor applied in calculation of surface
#' roughness and zero-plane displacement heights (see details)
#' @details The array of Plant Area index values in `vegp` must
#' of the same x and y dims as `dtm` but can contain any number of repeated
#' measures up to the number of entries in `weather`. Data are interpolated to the
#' time increment of `weather`. Other vegetation paramaters, including vegetation
#' height are assumed time-invarient. The SpatRaster datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined. In subsequent
#' downscaling of wind, the drag effects of vegetation, determined by vegetation height
#' and foliage area are accounted for can calculated at this stage. In so doing, it is
#' necessary to accommodate the possibility that the wind speed is not just affected
#' by the surface roughness in each pixel,  but also by vegetation surrounding the
#' location. This is accommodated for by applying `xyf` which effectively smooths
#' the surface roughness coefficients using `terra::aggregate` where `xyf` is the
#' aggregation factor. If `xyf` is set to `NA`, the roughness coefficients are averaged
#' across the study area.
#'
#' @seealso [checkinputs()], [modelin_dy()], [modelina()]
#'
#' @import terra
#' @export
modelin <- function(micropoint, vegp, soilc, dtm, runchecks = TRUE, xyf = 1) {
  weather<-micropoint$weather
  precip<-micropoint$precip
  if (runchecks) {
    rc<-checkinputs(weather,precip,vegp,soilc,dtm,FALSE,tstep=micropoint$tstep)
    weather<-rc$weather
    precip<-rc$precip
    vegp<-rc$vegp
    soilc<-rc$soilc
  }
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  # Turn input climate variables into arrays
  r<-dtm
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  h<-length(tme)
  tc<-.vta(weather$temp,r)
  difr<-.vta(weather$difrad,r)
  di<-weather$swrad-weather$difrad
  di[di<0]<-0
  jd<-.jday(tme)
  lt<-with(tme,hour+min/60+sec/3600)
  sa<-.solalt(lt,micropoint$lat,micropoint$long,jd)
  ze<-90-sa
  si<-cos(ze*(pi/180))
  si[si<0]<-0
  dni<-di/si
  dni[is.na(dni)]<-0
  dni[dni>1352]<-1352
  dirr<-.vta(dni,r)
  dp<-weather$difrad/weather$swrad
  dp[is.na(dp)]<-0.5; dp[dp<0]<-0; dp[dp>1]<-1
  dp<-.vta(dp,dtm)
  skyem<-.vta(weather$skyem,r)
  u2<-.vta(weather$windspeed,r)
  # Turn point model variables into arrays
  ufp<-.vta(micropoint$microp$uf,r)
  soilmp<-.vta(micropoint$soilm,r)
  Gp<-.vta(micropoint$microp$G,r)
  if (class(micropoint$Tbz) != "logical")  {
    Tbp<-.vta(micropoint$Tbz,r)
    Tg<-.vta(micropoint$microp$Tg,r)
  } else {
    Tbp<-NA
    Tg<-NA
  }
  xx<-.point0(micropoint)
  T0p<-.vta(.point0(micropoint),r)
  # Obtain derived variables
  estl<-.satvap(weather$temp)
  ea<-(weather$relhum/100)*estl
  tdew<-.dewpoint(ea,weather$temp)
  estl<-.vta(estl,r)
  ea<-.vta(ea,r)
  tdew<-.vta(tdew,r)
  vegx<-.rta(vegp$x,h)
  lref<-.rta(vegp$leafr,h)
  ltra<-.rta(vegp$leaft,h)
  pk<-.vta(weather$pres,r)
  gref<-.rta(soilc$groundr,h)
  # Calculate zero plane displacement height and roughness length
  dzm<-.sortrough(vegp,xyf)
  if (micropoint$tstep=="hour") {
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
  out<-list(# Misc variables:
    tme=tme,zref=micropoint$zref,dtm=dtm,lat=micropoint$lat,long=micropoint$long,
    # Weather variables:
    tc=tc,difr=difr,dirr=dirr,skyem=skyem,estl=estl,ea=ea,tdew=tdew,pk=pk,u2=u2,
    # Point model variables:
    ufp=ufp,soilmp=soilmp,Gp=Gp,Tbp=Tbp,Tg=Tg,T0p=T0p,vegp_p=micropoint$vegp_p,groundp_p=micropoint$groundp_p,DDp=micropoint$DD,
    # Vegetation parameters:
    pai=pai,vegx=vegx,lref=lref,ltra=ltra,veghgt=vegp$hgt,gsmax=vegp$gsmax,clump=clump,leafd=vegp$leafd,
    # Soil thermal parameters:
    gref=gref,rho=soilp$rho,Vm=soilp$Vm,Vq=soilp$Vq,Mc=soilp$Mc,
    # Soil hydraulic parameters:
    soilb=soilp$soilb,psi_e=soilp$psi_e,Smax=soilp$Smax,Smin=soilp$Smin,
    # Derived variables
    d=d,zm=zm,
    # Climate and timestep variables variables
    winddir=weather$winddir,tstep="hour",progress=0)
  class(out) <-"microin"
  return(out)
}
#' Run microclimate model (hourly)
#'
#' @description The function `runmicro_hr` runs the microclimate model in hourly
#' time increments
#'
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above (or below) ground at which model outputs are needed (m). Negative values indicate below ground surface.
#' @param pai_a an optional array of plant area index values above `reqhgt` (see details)
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from dtm if not supplied, but outer cells will be NA.
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
#' @param surfwet an optional single numeric value of array of values specifying the proportion
#' of the canopy surface that should be treated as wet surface (see details)
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
#'
#' @seealso [runmicro_dy()] for faster running microclimate model in daily time-steps,
#' including the option to expand daily outputs to hourly using the input diurnal
#' temperature cycle and [runmicro_big()] for running the microclimate model over large areas
#' as tiles.
#'
#' @details `pai_a` is used to calculate the radiation intercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calculated from total plant area index by
#' assuming a realistic shape to the vertical profile foliage within the canopy. If supplied,
#' `pai_a` must have the same dimensions as micro$pai. I.e. with the same x and y
#' dims as the the supplied dtm and values for each hour as the z dimension. The parameter `surfwet`
#' determines how much of the canopy should be treated as wet surface when calculating
#' latent heat fluxes. However, except when extremely droughted, the matric potential of leaves
#' is such that `surfwet` ~ 1. If set to NA, surfess wetness is modelled.
#' @import terra
#' @export
#'
#' @examples
#' library(terra)
#' library(zoo)
#' library(abind)
#' # ** First run point model (NB not run as uses inbuilt dataset)
#' # micropoint<-runpointmodel(climdata,rainfall,vegp,soilc,reqhgt = 0.05)
#' # ** Subset inbuilt point model to get hottest days in each month
#' micropoint_mx<-subsetpointmodel(micropoint, tstep = "month", what = "tmax")
#' micropoint_mn<-subsetpointmodel(micropoint, tstep = "month", what = "tmin")
#' # ** Create model inputs with inbuilt datasets
#' micro_mx<-modelin(micropoint_mx,vegp,soilc,dtmcaerth)
#' micro_mn<-modelin(micropoint_mn,vegp,soilc,dtmcaerth)
#' # Run model 5 cm above ground
#' mout_mx<-runmicro_hr(micro_mx, 0.05)
#' mout_mn<-runmicro_hr(micro_mn, 0.05)
#' # Plot air temperatures on hottest hour
#' plot(rast(mout_mx$Tz[,,134]))
#' # Plot mean of monthly max and min
#' mairt<-apply((mout_mn$Tz + mout_mx$Tz) / 2, c(1,2),mean)
#' plot(rast(mairt))
runmicro_hr <- function(micro, reqhgt, pai_a = NA, tfact = 1.5, slr = NA, apr = NA,
                        hor = NA, wsa = NA, twi = NA, surfwet = NA) {
  # Calculate soil surface temperature and soil moisture
  micro<-soiltemp_hr(micro,reqhgt,pai_a,tfact,slr,apr,hor,wsa,twi)
  # Run above ground
  if (reqhgt > 0) {
    mout<-temphumE(micro,reqhgt,pai_a,tfact,slr,apr,hor,wsa,twi,surfwet)
  }
  # Run at ground level
  if (reqhgt == 0) {
    mout<-list(Tz=micro$T0,tleaf=NA,T0=micro$T0,soilm=micro$soild,
               relhum=NA,windspeed=NA,Rdirdown=micro$radGdir,Rdifdown=micro$radGdif,
               Rlwdown=micro$radGlw,Rswup=(1-micro$gref)*micro$radGsw,
               Rlwup=0.97*5.67*10^-8*(micro$T0+273.15)^4)
  }
  # Run below ground
  if (reqhgt < 0) {
    micro<-below_hr(micro,reqhgt,pai_a,tfact,slr,apr,hor,wsa,twi)
    mout<-list(Tz=micro$Tz,tleaf=NA,T0=micro$T0,soilm=micro$soild,
               relhum=NA,windspeed=NA,Rdirdown=NA,Rdifdown=NA,Rlwdown=NA,
               Rswup=NA,Rlwup=NA)
  }
  mout$tme<-micro$tme
  class(mout)<-"microout"
  return(mout)
}
