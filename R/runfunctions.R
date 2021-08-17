#' Check format and values in model inputs
#'
#' The function `checkinputs` checks all the inputs into the model and returns errors
#' and warnings
#'
#' @param weather a data.frame of weather variables (see details)
#' @param rainfall a vector of daily rainfall
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a RasterLayer onject of elevations (see details)
#' @param merid optionally, longitude of local time zone meridian (decimal degrees)
#' @param dst optionally, numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#' @param daily optional logical indicating whether input weather data are daily
#' @return a list of checked inputs
#'
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. The array of Plant Area index values in `vegp` must
#' of the same x and y dims as `dtm` but can contain any number of repeated
#' measures up to the number of entries in `weather`. Data are interpolated to the
#' time increment of `weather`. Other vegetation paramaters, including vegetation
#' height are assumed time-invarient. The RasterLayer datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined.
#'
#' @export
#' @import raster
#'
#' @examples
#' # No warnings or errors given:
#' checks<-checkinputs(climdata, rainfall, vegp, soilc, dtmcaerth)
#' # Warning given (not run)
#' # weather<-climdata
#' # weather$relhum[1]<-101
#' # checks<-checkinputs(weather, rainfall, vegp, soilc, dtmcaerth)
#' # Error given (NB not run)
#' # weather<-climdata
#' # weather$pres<-weather$pres*1000
#' # checks<-checkinputs(weather, rainfall, vegp, soilc, dtmcaerth)
#' # Error given for vegp (not run)
#' # vegp2<-vegp
#' # vegp2$clump<-1
#' # checks<-checkinputs(climdata, rainfall, vegp2, soilc, dtmcaerth)
#' # Warning given for vegp (not run)
#' # vegp2<-vegp
#' # vegp2$pai<-vegp$pai*10
#' # checks<-checkinputs(climdata, rainfall, vegp2, soilc, dtmcaerth)
checkinputs <- function(weather, rainfall, vegp, soilc, dtm, merid = 0, dst = 0, daily = FALSE) {
  check.names<-function(nms,char) {
    sel<-which(nms==char)
    if (length(sel) == 0) stop(paste0("Cannot find ",char," in weather"))
  }
  check.vals<-function(x,mn,mx,char,unit) {
    sel<-which(is.na(x))
    if (length(sel)>0) stop(paste0("Missing values in weather$",char))
    sel<-which(x<mn)
    if (length(sel)>0) {
      mx<-min(x)
      stop(paste0(mx," outside range of typical ",char," values. Units should be ",unit))
    }
    sel<-which(x>mx)
    if (length(sel)>0) {
      mx<-max(x)
      stop(paste0(mx," outside range of typical ",char," values. Units should be ",unit))
    }
  }
  check.mean<-function(x,mn,mx,char,unit) {
    me<-mean(x,na.rm=T)
    if (me<mn | me>mx) stop(paste0("Mean ",char," of ",me," implausible. Units should be ",unit))
  }
  up.lim<-function(x,mx,char) {
    mxx<-max(x,na.rm=T)
    x[x>mx]<-mx
    if (mxx>mx) warning(paste0(char," values capped at ",mx))
    x
  }
  if (class(dtm)[1] != "RasterLayer") stop("dtm must be a raster")
  # check names
  nms<-names(weather)
  check.names(nms,"obs_time")
  check.names(nms,"temp")
  check.names(nms,"relhum")
  check.names(nms,"pres")
  check.names(nms,"swrad")
  check.names(nms,"difrad")
  check.names(nms,"skyem")
  check.names(nms,"windspeed")
  check.names(nms,"winddir")
  # get si
  if (is.na(crs(dtm))) stop("dtm must have a coordinate reference system specified")
  ll<-.latlongfromraster(dtm)
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  # check tme
  sel<-which(is.na(tme))
  if (length(sel)>0) stop("Cannot recognise all obs_time in weather")
  jd<-.jday(tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  sa<-.solalt(lt,ll$lat,ll$long,jd,merid,dst)
  ze<-90-sa
  si<-cos((90-sa)*(pi/180))
  si[si<0]<-0
  # Check weather  values
  check.vals(weather$temp,-50,65,"temperature","deg C")
  weather$relhum<-up.lim(weather$relhum,100,"relative humidity")
  check.vals(weather$relhum,0,100,"relative humidity","percentage (0-100)")
  check.mean(weather$relhum,5,100,"relative humidity","percentage (0-100)")
  check.vals(weather$pres,84,109,"pressure","kPa ~101.3")
  check.vals(weather$swrad,0,1350,"shortwave radiation","W / m^2")
  check.vals(weather$difrad,0,1350,"diffuse radiation","W / m^2")
  weather$skyem<-up.lim(weather$skyem,1,"sky emissivity")
  check.mean(weather$skyem,0.3,0.99,"sky emissivity","in range 0 - 1")
  check.vals(weather$skyem,0.3,1,"sky emissivity","in range 0 - 1")
  check.vals(weather$windspeed,0,100,"wind speed","m/s")
  if (max(weather$windspeed)>30) warning("Maximum wind speed seems quite high. Check units are m/s and for 2 m above ground")
  # Check direct radiation at low solar angles
  dirr<-weather$swrad-weather$difrad
  sel<-which(dirr<0)
  if (length(sel)>0) {
    weather$difrad[sel]<-weather$swrad[sel]
    warning("Diffuse radiation values higher than shortwave radiation, and so was set to shortwave radiation values")
  }
  dni<-dirr/si
  dni[is.na(dni)]<-0
  dni2<-ifelse(dni>1352,1352,dni)
  dirr2<-dni2*si
  dif<-dirr-dirr2
  if (max(dif) > 0.01) {
    warning("Computed direct radiation values too high near dawn / dusk. Capped, and excess assigned as diffuse radiation")
    weather$difrad<-weather$difrad+dif
  }
  dirr<-weather$swrad-weather$difrad
  sel<-which(dirr<0)
  weather$swrad[sel]<-weather$difrad[sel]
  weather$swrad<-ifelse(weather$swrad>1352,1352,weather$swrad)
  weather$difrad<-ifelse(weather$difrad>1352,1352,weather$difrad)
  # Wind direction
  mn<-min(weather$winddir)
  mx<-max(weather$winddir)
  if (mn<0 | mx>360) {
    weather$winddir<-weather$winddir%%360
    warning("wind direction adjusted to range 0-360 using modulo operation")
  }
  # Check other variables
  hrs<-dim(weather)[1]
  if (daily) {
    dys<-hrs
  } else dys<-hrs/24
  if(dys != floor(dys)) stop ("weather needs to include data for entire days (24 hours)")
  if(length(rainfall) != dys) stop("duration of rainfall sequence doesn't match weather'")
  # Check vegp data
  xy<-dim(dtm)[1:2]
  if (dim(vegp$pai)[1] != xy[1]) stop("y dimension of vegp$pai does not match dtm")
  if (dim(vegp$pai)[2] != xy[2]) stop("x dimension of vegp$pai does not match dtm")
  if (class(vegp$x)[1] != "RasterLayer") stop("vegp$x must be a raster")
  if (class(vegp$gsmax)[1] != "RasterLayer") stop("vegp$gsmax must be a raster")
  if (class(vegp$leafr)[1] != "RasterLayer") stop("vegp$leafr must be a raster")
  if (class(vegp$clump)[1] != "RasterLayer") stop("vegp$clump must be a raster")
  if (class(vegp$leafd)[1] != "RasterLayer") stop("vegp$leafd must be a raster")
  if (dim(vegp$x)[1] != xy[1]) stop("y dimension of vegp$x does not match dtm")
  if (dim(vegp$x)[2] != xy[2]) stop("x dimension of vegp$x does not match dtm")
  if (dim(vegp$gsmax)[1] != xy[1]) stop("y dimension of vegp$gsmax does not match dtm")
  if (dim(vegp$gsmax)[2] != xy[2]) stop("x dimension of vegp$gsmax does not match dtm")
  if (dim(vegp$leafr)[1] != xy[1]) stop("y dimension of vegp$leafr does not match dtm")
  if (dim(vegp$leafr)[2] != xy[2]) stop("x dimension of vegp$leafr does not match dtm")
  if (dim(vegp$clump)[1] != xy[1]) stop("y dimension of vegp$clump does not match dtm")
  if (dim(vegp$clump)[2] != xy[2]) stop("x dimension of vegp$clump does not match dtm")
  if (dim(vegp$leafd)[1] != xy[1]) stop("y dimension of vegp$leafd does not match dtm")
  if (dim(vegp$leafd)[2] != xy[2]) stop("x dimension of vegp$leafd does not match dtm")
  if (dim(vegp$x)[3]>1) stop("time variant vegp$x not supported")
  if (dim(vegp$gsmax)[3]>1) stop("time variant vegp$x not supported")
  if (dim(vegp$leafr)[3]>1) stop("time variant vegp$leafr not supported")
  if (dim(vegp$clump)[3]>1) stop("time variant vegp$clump not supported")
  if (dim(vegp$leafd)[3]>1) stop("time variant vegp$leafd not supported")
  # Check soil data
  xx<-unique(getValues(soilc$soiltype))
  if (max(xx,na.rm=T) > 11) stop("Unrecognised soil type")
  if (min(xx,na.rm=T) < 1) stop("Unrecognised soil type")
  if (class(soilc$soiltype)[1] != "RasterLayer") stop("soilc$soiltype must be a raster")
  if (class(soilc$groundr)[1] != "RasterLayer") stop("soilc$groundr must be a raster")
  if (dim(soilc$soiltype)[1] != xy[1]) stop("y dimension of soilc$soiltype does not match dtm")
  if (dim(soilc$soiltype)[2] != xy[2]) stop("x dimension of soilc$soiltype does not match dtm")
  if (dim(soilc$groundr)[1] != xy[1]) stop("y dimension of soilc$groundr does not match dtm")
  if (dim(soilc$groundr)[2] != xy[2]) stop("x dimension of soilc$groundr does not match dtm")
  if (dim(soilc$soiltype)[3]>1) stop("time variant soilc$soiltype not supported")
  if (dim(soilc$groundr)[3]>1) stop("time variant soilc$groundr not supported")
  xx<-getValues(soilc$groundr)
  xx<-xx[is.na(xx)==F]
  check.vals(xx,0,1,"soil reflectivity","range 0 to 1")
  xx<-getValues(vegp$leafr)
  xx<-xx[is.na(xx)==F]
  check.vals(xx,0,1,"leaf reflectivity","range 0 to 1")
  xx<-getValues(vegp$clump)
  xx<-xx[is.na(xx)==F]
  check.vals(xx,0,1,"vegetation clumping factor","range 0 to 1")
  xx<-getValues(vegp$leafd)
  xx<-xx[is.na(xx)==F]
  check.vals(xx,0,5,"leaf diamater","metres")
  xx<-getValues(vegp$gsmax)
  xx<-xx[is.na(xx)==F]
  check.vals(xx,0,2,"maximum stomatal conductance","mol / m^2 /s")
  if (mean(xx)>1) warning(paste0("Mean leaf diameter of ",mean(xx)," seems large. Check units are in metres"))
  # PAI
  if (class(vegp$pai)[1] != "array") {
    if (class(vegp$pai)[1] == "RasterLayer") {
      vegp$pai<-array(getValues(vegp$pai,format="matrix"))
      warning("vegp$pai assumed time-invariant and converted to an array")
    } else stop("vegp$pai must be an array or a raster")
  }
  xx<-vegp$pai
  xx<-xx[is.na(xx)==F]
  if (min(xx)<0) stop("Minimum vegp$pai must be greater than or equal to zero")
  if (max(xx)>15) warning(paste0("Maximum vegp$pai of ",max(xx)," seems high"))
  return(list(weather=weather,rainfall=rainfall,vegp=vegp,soilc=soilc))
}

#' Create object of class microin with weather data as data.frame
#'
#' @description The function `modelin` creates an object of class microin
#' which unpacks various component inputs and reformats as required
#' for running the model in hourly timesteps. Here it is assumed that the input
#' weather data are a data.frame - i.e. not spatially variable.
#'
#' @param weather a data.frame of weather variables (see details)
#' @param rainfall a vector of daily rainfall
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a RasterLayer object of elevations (see details)
#' @param windhgt height of wind speed measurement (m) in weather dataset (see details).
#' @param merid optionally, longitude of local time zone meridian (decimal degrees)
#' @param dst optionally, numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#' @param runchecks optional logical indicating whteher to call [checkinputs()] to run
#' checks on format and units of input data.
#' @param daily optional logical indicating whether `weather` is daily or hourly
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. The array of Plant Area index values in `vegp` must
#' of the same x and y dims as `dtm` but can contain any number of repeated
#' measures up to the number of entries in `weather`. Data are interpolated to the
#' time increment of `weather`. Other vegetation paramaters, including vegetation
#' height are assumed time-invarient. The RasterLayer datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined. As not all wind
#' measurements are at reference height, the height of the wind speed measurement
#' must be specified if not 2 m. To enable calculation of below-canopy wind profiles
#' in tall canopy, the wind speed is adjusted to give values for a height 2 m above
#' the maximum vegetation height, using the wind-height profile for a reference grass
#' surface.
#'
#' @seealso [checkinputs()], [modelin_dy()], [modelina()]
#'
#' @import raster sp
#' @export
modelin <- function(weather, rainfall, vegp, soilc, dtm, windhgt = 2, merid = 0, dst = 0, runchecks = TRUE,daily = FALSE) {
  if (runchecks) {
    rc<-checkinputs(weather,rainfall,vegp,soilc,dtm,merid,dst,daily)
    weather<-rc$weather
    rainfall<-rc$rainfall
    vegp<-rc$vegp
    soilc<-rc$soilc
  }
  # correct wind height
  mxhgt<-max(.is(vegp$hgt),na.rm=T)
  weather$windspeed<-.windcorrect(weather$windspeed,windhgt,mxhgt)
  r<-dtm
  ll<-.latlongfromraster(r)
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  h<-length(tme)
  tc<-.vta(weather$temp,r)
  difr<-.vta(weather$difrad,r)
  di<-weather$swrad-weather$difrad
  di[di<0]<-0
  jd<-.jday(tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  sa<-.solalt(lt,ll$lat,ll$long,jd,merid,dst)
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
  # == (0b) Derived variables
  estl<-.satvap(weather$temp)
  ea<-(weather$relhum/100)*estl
  tdew<-.dewpoint(ea,weather$temp)
  estl<-.vta(estl,r)
  ea<-.vta(ea,r)
  tdew<-.vta(tdew,r)
  vegx<-.rta(vegp$x,h)
  lref<-.rta(vegp$leafr,h)
  pk<-.vta(weather$pres,r)
  gref<-.rta(soilc$groundr,h)
  pai=.unpackpai(vegp$pai,h)
  clump<-.rta(vegp$clump,h)
  soilp<-.soilinit(soilc)
  out<-list(tme=tme,tc=tc,difr=difr,dirr=dirr,dp=dp,skyem=skyem,
            estl=estl,ea=ea,tdew=tdew,pk=pk,pai=pai,vegx=vegx,lref=lref,veghgt=vegp$hgt,
            gsmax=vegp$gsmax,clump=clump,gref=gref,rho=soilp$rho,Vm=soilp$Vm,leafd=vegp$leafd,
            Vq=soilp$Vq,Mc=soilp$Mc,soilb=soilp$soilb,psi_e=soilp$psi_e,Smax=soilp$Smax,
            dtm=dtm,lat=ll$lat,long=ll$long,merid=merid, dst=dst,
            climdata=weather,prec=rainfall,soilc=soilc)
  class(out) <-"microin"
  out
}
#' Create object of class microin with weather data as an array
#'
#' @description The function `modelin` creates an object of class microin
#' which unpacks various component inputs and reformats as required
#' for running the model in hourly timesteps. Here it is assumed that the input
#' weather data are as arrays - i.e. variable in space
#'
#' @param climarray a list of arrays of weather variables (see details). See also [nctoarray()]
#' @param rainarray an array of daily rainfall (see details)
#' @param tme an object of class POSIXlt giving the dates and times for each weather variable stroed in the array
#' @param r a raster object giving with the resolution, spatial extent, and projection of the weather data (see details)
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction, see details)
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a RasterLayer onject of elevations (see details)
#' @param merid optionally, longitude of local time zone meridian (decimal degrees)
#' @param dst optionally, numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' checks on format and units of input data.
#' @param daily optional logical indicating whether `climarray` is daily or hourly
#' @details The units of `climarray` must follow those in the dataset `climdata`.
#' It must be a list with each component of the list an array, named using the same
#' names as the column headers in weather (e.g. temp for temperature), excluding `obs_time`.
#' Dimensions 1 and 2 of the array must be the same as `r` and dimension 3 must have
#' the same length as `tme`. If `r` has a different resolution to `dtm` the climate
#' data are resampled to match the resolution of `dtm`. The array of Plant Area index values in `vegp` must
#' of the same x and y dims as `dtm` but can contain any number of repeated
#' measures up to the number of entries in `tme`. Data are interpolated to the
#' time increment of `tme`. Other vegetation paramaters, including vegetation
#' height are assumed time-invarient. The RasterLayer datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined. If `altcorrect`>0,
#' and the dimenions of `r` are not identical to those of `dtm`, the elevation
#' difference between each pixel of the dtm and the dtm coarsed to the resolution of
#' `r` is calaculated and an elevational lapse rate correction is applied to the
#' temperature data to accoutn for these elevation differences. If `altcorrect`=1,
#' a fixed lapse rate of 5 degrees per 100m is applied. If `altcorrect`=2, humidity-dependent
#' lapse rates are calaculate and applied.
#'
#' @seealso [modelin()], [modelina_dy()], [nctoarray()]
#'
#' @import raster
#' @export
modelina<-function(climarray,rainarray,tme,r,altcorrect = 0, vegp, soilc, dtm, merid = 0, dst = 0, runchecks = FALSE, daily = FALSE) {
  # Create weather and rainfall dataset
  weather<-.catoweather(climarray)
  rainfall<-apply(rainarray,3,mean,na.rm=T)
  # Run checks
  if (runchecks) {
    rc<-checkinputs(weather,rainfall,vegp,soilc,dtm,merid,dst,daily)
    weather<-rc$weather
    rainfall<-rc$rainfall
    vegp<-rc$vegp
    soilc<-rc$soilc
  }
  # Resample climdata
  tc<-suppressWarnings(.resa(climarray$temp,r,dtm))
  difr<-suppressWarnings(.resa(climarray$difrad,r,dtm))
  # ~~ Calculate dni
  di<-climarray$swrad-climarray$difrad
  jd<-.jday(tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  ll<-.latslonsfromr(r)
  n<-length(tme)
  sa<-.solalt(.vta(lt,r),.rta(raster(ll$lats),n),.rta(raster(ll$lons),n),.vta(jd,r),merid,dst)
  ze<-90-sa
  si<-cos(ze*(pi/180))
  si[si<0]<-0
  dirr<-suppressWarnings(.resa(di/si,r,dtm))
  dirr[is.na(dirr)]<-0
  dirr[dirr>1352]<-1352; dirr[dirr<0]<-0
  dp<-climarray$difrad/climarray$swrad
  dp[is.na(dp)]<-0.5; dp[dp<0]<-0; dp[dp>1]<-1
  dp<-suppressWarnings(.resa(dp,r,dtm))
  skyem<-suppressWarnings(.resa(climarray$skyem,r,dtm))
  pk<-suppressWarnings(.resa(climarray$pres,r,dtm))
  # ~~ Derived variables
  estl<-.satvap(climarray$temp)
  ea<-(climarray$relhum/100)*estl
  tdew<-.dewpoint(ea,climarray$temp)
  estl<-suppressWarnings(.resa(estl,r,dtm))
  ea<-suppressWarnings(.resa(ea,r,dtm))
  tdew<-suppressWarnings(.resa(tdew,r,dtm))
  vegx<-.rta(vegp$x,n)
  lref<-.rta(vegp$leafr,n)
  gref<-.rta(soilc$groundr,n)
  pai<-.unpackpai(vegp$pai,n)
  clump<-.rta(vegp$clump,n)
  soilp<-.soilinit(soilc)
  # Elevation correction
  # ~~ Fix lapse rate
  if (altcorrect>0) {
    dtmc<-resample(dtm,r)
    dtmc<-resample(dtmc,dtm)
    elevd<-dtmc-dtm
  }
  if (altcorrect==1) {
    tcdif<-elevd*(5/1000)
    tc<-.rta(tcdif,n)+tc
  }
  if (altcorrect==2) {
    lr<-.lapserate(climarray$temp,climarray$relhum,climarray$pres)
    lr<-suppressWarnings(.resa(lr,r,dtm))
    tcdif<-.rta(elevd,n)*lr
    tc<-tcdif+tc
  }
  out<-list(tme=tme,tc=tc,difr=difr,dirr=dirr,dp=dp,skyem=skyem,
            estl=estl,ea=ea,tdew=tdew,pk=pk,pai=pai,vegx=vegx,lref=lref,veghgt=vegp$hgt,
            gsmax=vegp$gsmax,clump=clump,gref=gref,rho=soilp$rho,Vm=soilp$Vm,leafd=vegp$leafd,
            Vq=soilp$Vq,Mc=soilp$Mc,soilb=soilp$soilb,psi_e=soilp$psi_e,Smax=soilp$Smax,
            dtm=dtm,lat=ll$lat,long=ll$long,merid=merid, dst=dst,
            climdata=weather,prec=rainfall,soilc=soilc)
  class(out) <-"microin"
  out
}
#' Create object of class microindaily
#'
#' @description The function `modelin` creates an object of class microindaily
#' which unpacks various component inputs and reformats as required
#' for running the model in daily timesteps
#'
#' @param weather a data.frame of hourly weather variables (see details)
#' @param rainfall a vector of daily rainfall
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a RasterLayer onject of elevations (see details)
#' @param merid optionally, longitude of local time zone meridian (decimal degrees)
#' @param dst optionally, numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#'
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. The array of Plant Area index values in `vegp` must
#' of the same x and y dims as `dtm` but can contain any number of repeated
#' measures up to the number of entries in `weather`. Data are interpolated to the
#' time increment of `weather`. Other vegetation paramaters, including vegetation
#' height are assumed time-invarient. The RasterLayer datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined.
#'
#' @seealso [inputchecks()], [modelin()]
#'
#' @import raster
#' @export
modelin_dy <- function(weather, rainfall, vegp, soilc, dtm, windhgt = 2, merid = 0, dst = 0, runchecks = TRUE) {
  climd<-.climtodaily(weather)
  micro_mn<-modelin(climd$climn,rainfall,vegp,soilc,dtm,windhgt,merid,dst,runchecks,daily=TRUE)
  micro_mx<-modelin(climd$climx,rainfall,vegp,soilc,dtm,windhgt,merid,dst,runchecks,daily=TRUE)
  out<-list(micro_mn=micro_mn,micro_mx=micro_mx,climdata=weather)
  class(out)<-"microindaily"
  return(out)
}
#' Create object of class microindaily with weather data as an array
#'
#' @description The function `modelin` creates an object of class microin
#' which unpacks various component inputs and reformats as required
#' for running the model in hourly timesteps. Here it is assumed that the input
#' weather data are as arrays - i.e. variable in space
#'
#' @param climarray a list of arrays of weather variables (see details). See also [nctoarray()]
#' @param rainarray an array of daily rainfall (see details)
#' @param tme an object of class POSIXlt giving the dates and times for each weather variable stroed in the array
#' @param r a raster object giving with the resolution, spatial extent, and projection of the weather data (see details)
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction, see details)
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a RasterLayer onject of elevations (see details)
#' @param merid optionally, longitude of local time zone meridian (decimal degrees)
#' @param dst optionally, numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' checks on format and units of input data.
#' @details The units of `climarray` must follow those in the dataset `climdata`.
#' It must be a list with each component of the list an array, named using the same
#' names as the column headers in weather (e.g. temp for temperature), excluding `obs_time`.
#' Dimensions 1 and 2 of the array must be the same as `r` and dimension 3 must have
#' the same length as `tme`. If `r` has a different resolution to `dtm` the climate
#' data are resampled to match the resolution of `dtm`. The array of Plant Area index values in `vegp` must
#' of the same x and y dims as `dtm` but can contain any number of repeated
#' measures up to the number of entries in `tme`. Data are interpolated to the
#' time increment of `tme`. Other vegetation paramaters, including vegetation
#' height are assumed time-invarient. The RasterLayer datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined. If `altcorrect`>0,
#' and the dimenions of `r` are not identical to those of `dtm`, the elevation
#' difference between each pixel of the dtm and the dtm coarsed to the resolution of
#' `r` is calaculated and an elevational lapse rate correction is applied to the
#' temperature data to accoutn for these elevation differences. If `altcorrect`=1,
#' a fixed lapse rate of 5 degrees per 100m is applied. If `altcorrect`=2, humidity-dependent
#' lapse rates are calaculate and applied.
#'
#' @seealso [modelin_dy()], [modelina()], [nctoarray()]
#'
#' @import raster
#' @export
modelina_dy <- function(climarray, rainarray, tme, r, altcorrect = 0, vegp, soilc, dtm, merid = 0, dst = 0, runchecks = FALSE) {
  climdata<-.catoweather(climarray)
  climd<-.climtodaily(climdata)
  tme2<-tme[climd$smn]
  micro_mn<-modelina(climd$climn,rainarray,tme2,r,altcorrect,vegp,soilc,dtm,merid,dst,runchecks,daily=TRUE)
  tme2<-tme[climd$smx]
  micro_mx<-modelina(climd$climn,rainarray,tme2,r,altcorrect,vegp,soilc,dtm,merid,dst,runchecks,daily=TRUE)
  out<-list(micro_mn=micro_mn,micro_mx=micro_mx,climdata=climdata)
  class(out)<-"microindaily"
  return(out)
}
#' Run microclimate model (hourly)
#'
#' @description The function `runmicro_hr` runs the microclimate model in hourly
#' time increments
#'
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above ground at which model outputs are needed (m).
#' @param pai_a an optional array of plant area index values above `reqhgt` (see details)
#' @param xyf optional input for called function [wind()]
#' @param zf optional input for called function [wind()]
#' @param soilinit initial soil moisture fractions in surface and subsurface layer (see [soilmpredict()])
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param surfwet an optional single numeric value of array of values specifying the proportion
#' of the canopy surface that should be treated as wet surface (see details)
#' @param slr an optional RasterLayer of slope values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param apr an optional RasterLayer of aspect values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from dtm if not supplied, but outer cells will be NA.
#' @param maxhgt an optional height (m) for which wind speed is needed. Determined
#' from height of tallest vegetation or as 2 m, whichever is greater, if not supplied.
#' @param twi optional raster of topographic wetness index values.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
#' @return an object of class microout with the following components:
#' @return `Tz` Array of air temperatures at height `reqhgt` (deg C). Identical to `T0`
#' if `reqhgt` = 0.
#' @return `tleaf` Array of leaf temperatures at height `reqhgt` (deg C).
#' NA if `reqhgt` greater than canopy height or `reqhgt` <= 0.
#' @return `T0` Array of ground surface temperatures (deg C)
#' @return `relhum` Array of relative humidities at height `reqhgt` (percentage).
#' NA if `reqhgt` <= 0.
#' @return `windspeed` Array of wind speeds at height `reqhgt` (m/s).
#' NA if `reqhgt` <= 0.
#' @return `raddir` Array of direct shortwave radiation received on horizontal
#' surface (W/m^2). NA if `reqhgt` <= 0.
#' @return `raddif` Array of direct shortwave radiation received on horizontal
#' surface (W/m^2). NA if `reqhgt` <= 0.
#' @return `radlw` Array of downward longwave radiation received on horizontal
#' surface (W/m^2). NA if `reqhgt` <= 0.
#'
#' @seealso [runmicro_dy()] for faster running microclimate model in daily time-steps,
#' including the option to expand daily outputs to hourly using the input diurnal
#' temperature cycle.
#'
#' @details `pai_a` is used to calaculate the radiation incercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calaculated from total plant area index by
#' assuming leaf density within the canopy is uniformly vertically distributed. If suuplied
#' it must have the same dimensions as micro$pai. I.e. with the same x and y dims as the
#' the supplied dtm and values for each hour as the z dimension. The paramater `surfwet`
#' determines how much of the canopy should be treated as wet surface when calaculating
#' latent heat fluxes. However, except when extremely droughted, the matric potential of leaves
#' is such that `surfwet` ~ 1.
#' @import raster
#' @export
#'
#' @examples
#' library(raster)
#' library(sp)
#' library(zoo)
#' library(abind)
#' # Create model input with inbuilt datasets
#' micro<-modelin(climdata,rainfall,vegp,soilc,dtmcaerth)
#' # Run model 5 cm above ground (takes ~ 3 mins to run on 50 x 50 x 8760 values)
#' mout<-runmicro_hr(micro, 0.05)
#' # Plot air temperatures on hottest hour
#' plot(raster(mout$Tz[,,4094]))
#' # Plot mean air temperatures
#' mairt<-apply(mout$Tz,c(1,2),mean)
#' plot(raster(mairt))
#' # Plot ground temperatures on hottest and coldest hour
#' plot(raster(mout$T0[,,991])) # coldest hour
#' plot(raster(mout$T0[,,4094])) # hottest hour
#' # Run model 20 cm below ground
#' mout<-runmicro_hr(micro, -0.2)
#' # Extract and plot mean soil temperatures
#' msoilt<-apply(mout$Tz,c(1,2),mean)
#' plot(raster(msoilt))
runmicro_hr <- function(micro, reqhgt, pai_a = NA, xyf = NA, zf = NA, soilinit = c(NA, NA),
                        tfact = 1.5, surfwet = 1, slr = NA, apr = NA, hor = NA, wsa = NA,
                        maxhgt = NA, twi = NA) {
  # Calculate soil surface temperature and soil moisture
  micro<-soiltemp_hr(micro,reqhgt,xyf,zf,soilinit,tfact,slr,apr,hor,wsa,maxhgt,twi)
  # Run above ground
  if (reqhgt > 0) {
    mout<-temphumE(micro,reqhgt,pai_a,xyf,zf,soilinit,tfact,surfwet,slr,apr,hor,wsa,maxhgt,twi)
  }
  # Run at ground level
  if (reqhgt == 0) {
    mout<-list(Tz=micro$T0,tleaf=NA,T0=micro$T0,soilm=micro$theta,
               relhum=NA,windspeed=NA,raddir=NA,raddif=NA,radlw=NA)
    class(mout)<-"microout"
  }
  # Run below ground
  if (reqhgt < 0) {
    Tz<-below_hr(micro,reqhgt,xyf,zf,soilinit,tfact,slr,apr,hor,wsa,maxhgt,twi)
    mout<-list(Tz=Tz,tleaf=NA,T0=micro$T0,soilm=micro$theta,
               relhum=NA,windspeed=NA,raddir=NA,raddif=NA,radlw=NA)
    class(mout)<-"microout"
  }
  return(mout)
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
#' @param xyf optional input for called function [wind()]
#' @param zf optional input for called function [wind()]
#' @param soilinit initial soil moisture fractions in surface and subsurface layer (see [soilmpredict()])
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param surfwet an optional single numeric value of array of values specifying the proportion
#' of the canopy surface that should be treated as wet surface (see details)
#' @param slr an optional RasterLayer of slope values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param apr an optional RasterLayer of aspect values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from dtm if not supplied, but outer cells will be NA.
#' @param maxhgt an optional height (m) for which wind speed is needed. Determined
#' from height of tallest vegetation or as 2 m, whichever is greater, if not supplied.
#' @param twi optional raster of topographic wetness index values.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
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
#' @return `raddir` Array of direct shortwave radiation received on horizontal
#' surface (W/m^2). NA if `reqhgt` <= 0.
#' @return `raddif` Array of direct shortwave radiation received on horizontal
#' surface (W/m^2). NA if `reqhgt` <= 0.
#' @return `radlw` Array of downward longwave radiation received on horizontal
#' surface (W/m^2). NA if `reqhgt` <= 0.
#' @return if expand = FALSE, an object of class microutdaily, list with the following components:
#' @return mout_mn an object of class microut for minimum daily temperatures
#' @return mout_mx an object of class microut for maximum daily temperatures
#' @seealso [runmicro_hr()] for running microclimate model in hourly time-steps
#'
#' @details
#' If expand = TRUE, daily minima and maxima are expanded to hourly using values in the hourly
#' weather dataset. The parmater `pai_a` is used to calaculate the radiation incercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calaculated from total plant area index by
#' assuming leaf density within the canopy is uniformly vertically distributed. If suuplied
#' it must have the same dimensions as micro$pai. I.e. with the same x and y dims as the
#' the supplied dtm and values for each hour as the z dimension. The paramater `surfwet`
#' determines how much of the canopy should be treated as wet surface when calaculating
#' latent heat fluxes. However, except when extremely droughted, the matric potential of leaves
#' is such that `surfwet` ~ 1.
#'
#' @import raster
#' @export
#'
#' @examples
#' library(raster)
#' library(sp)
#' library(zoo)
#' library(abind)
#' # Create model input with inbuilt datasets
#' microd<-modelin_dy(climdata,rainfall,vegp,soilc,dtmcaerth)
#' # Run model 5 cm above ground  performing hourly expansion (takes ~50 seconds to run)
#' mout<-runmicro_dy(microd, 0.05)
#' # Plot air temperatures on hottest hour
#' plot(raster(mout$Tz[,,4094]))
#' # Plot mean air temperatures
#' mairt<-apply(mout$Tz,c(1,2),mean)
#' plot(raster(mairt))
#' # Plot ground temperatures on hottest and coldest hour
#' plot(raster(mout$T0[,,991])) # coldest hour
#' plot(raster(mout$T0[,,4094])) # hottest hour
#' # Run model 20 cm below ground
#' mout<-runmicro_dy(microd, -0.2)
#' # Extract and plot mean soil temperatures
#' msoilt<-apply(mout$Tz,c(1,2),mean)
#' plot(raster(msoilt))
runmicro_dy <- function(microd, reqhgt, expand = TRUE, pai_a = NA, xyf = NA, zf = NA,
                        soilinit = c(NA, NA), tfact = 1.5, surfwet = 1, slr = NA,
                        apr = NA, hor = NA, wsa = NA, maxhgt = NA, twi = NA) {
  # Calculate soil surface temperature and soil moisture
  microd<-soiltemp_dy(microd,reqhgt,xyf,zf,soilinit,tfact,slr,apr,hor,wsa,maxhgt,twi)
  climdata<-microd$climdata
  climd<-.climtodaily(climdata)
  # Run above ground
  if (reqhgt > 0) {
    mout_mn<-temphumE(microd$micro_mn,reqhgt,pai_a,xyf,zf,soilinit,tfact,surfwet,slr,apr,hor,wsa,maxhgt,twi)
    mout_mx<-temphumE(microd$micro_mx,reqhgt,pai_a,xyf,zf,soilinit,tfact,surfwet,slr,apr,hor,wsa,maxhgt,twi)
    if (expand) {
      mout<-.expandclim(mout_mn,mout_mx,climdata)
      class(mout)<-"microout"
    } else {
      mout<-list(mout_mn=mout_mn,mout_mx=mout_mx)
      class(mout)<-"microoutdaily"
    }
  }
  mout_mn<-microd$micro_mn
  mout_mx<-microd$micro_mx
  # Run at ground level
  if (reqhgt == 0) {
    if (expand) {
      Tz<-.expandtohour(mout_mn$T0,mout_mx$T0,climd$smn,climd$smx,climdata$temp)
      mout<-list(Tz=Tz,tleaf=NA,T0=Tz,soilm=.ehr(mout_mn$theta),
                 relhum=NA,windspeed=NA,raddir=NA,raddif=NA,radlw=NA)
      class(mout)<-"microout"
    } else {
      mout_mn$trdf<-NULL
      mout_mx$trdf<-NULL
      mout<-list(mout_mn=mout_mn,mout_mx=mout_mx)
      class(mout)<-"microoutdaily"
    }
  }
  # Run below ground
  if (reqhgt < 0) {
    Tz<-below_dy(microd,reqhgt,expand,xyf,zf,soilinit,tfact,slr,apr,hor,wsa,maxhgt,twi)
    if (expand) {
      T0<-.expandtohour(mout_mn$T0,mout_mx$T0,climd$smn,climd$smx,climdata$temp)
      mout<-list(Tz=Tz,tleaf=NA,T0=T0,soilm=.ehr(mout_mn$theta),
                 relhum=NA,windspeed=NA,raddir=NA,raddif=NA,radlw=NA)
      class(mout)<-"microout"
    } else {
      micro_mn$Tz<-Tz
      micro_mx$Tz<-Tz
      mout_mn$trdf<-NULL
      mout_mx$trdf<-NULL
      mout<-list(mout_mn=mout_mn,mout_mx=mout_mx)
      class(mout)<-"microoutdaily"
    }
  }
  return(mout)
}
#' Writes model output as ncdf4 file
#'
#' @description The function `writetonc` writes hourly model outputs as an ncdf4 file to a specified file
#'
#' @param mout an object of class microut
#' @param fileout output filename
#' @param dtm a Rasterlayer covering the extent of the model outputs
#' @param weather a data.frame of hourly weather variables (used to get time stamp)
#' @param reqhgt at at which model was run (m)
#' @param merid longitude of local time zone meridian (decimal degrees)
#' @param dst numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#'
#' @import ncdf4 raster
#' @export
writetonc <- function(mout, fileout, dtm, weather, reqhgt, merid = 0, dst = 0) {
  atonc<-function(a,rd) {
    a<-apply(a,c(2,3),rev)
    a<-aperm(a,c(2,1,3))
    a <- round(a*rd,0)
    a <-array(as.integer(a),dim=dim(a))
    a
  }
  # Creat eastings and nothings sequence
  est<-seq(extent(dtm)@xmin+res(dtm)[1]/2,extent(dtm)@xmax-res(dtm)[1]/2,res(dtm)[1])
  nth<-seq(extent(dtm)@ymin+res(dtm)[2]/2,extent(dtm)@ymax-res(dtm)[2]/2,res(dtm)[2])
  east<-ncdim_def(name="east",units="metres",vals=est,longname="Eastings")
  north<-ncdim_def(name="north",units="metres",vals=nth,longname="Northings")
  # Create time variable
  tme<-as.POSIXct(weather$obs_time)
  times<-ncdim_def(name="Time",units="Decimal hours since 1970-01-01 00:00",vals=as.numeric(tme)/3600)
  # Variable names
  if (reqhgt > 0) {
    tname<-paste0("Air temperature at height ",reqhgt," m")
    lname<-paste0("Leaf temperature at height ",reqhgt," m")
    rname<-paste0("Relative humidity at height ",reqhgt," m")
    wname<-paste0("Wind speed at height ",reqhgt," m")
    # Define variables
    airtemp<-ncvar_def(name=tname,units="deg C x 100",dim=list(east,north,times),
                       missval=-9999,compression=9,prec="integer")
    leaftemp<-ncvar_def(name=lname,units="deg C x 100",dim=list(east,north,times),
                        missval=-9999,compression=9,prec="integer")
    relhum<-ncvar_def(name=rname,units="Percentage",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    windspeed<-ncvar_def(name=wname,units="m/s x 100",dim=list(east,north,times),
                         missval=-9999,compression=9,prec="integer")
    raddir<-ncvar_def(name="Downward direct shortwave radiation",units="W/m^2",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    raddif<-ncvar_def(name="Downward diffuse shortwave radiation",units="W/m^2",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    radlw<-ncvar_def(name="Downward longwave radiation",units="W/m^2",dim=list(east,north,times),
                     missval=-9999,compression=9,prec="integer")
    # Create nc file
    nc.name<-fileout
    ncnew<-nc_create(filename=nc.name,list(airtemp,leaftemp,relhum,windspeed,raddir,raddif,radlw))
    # Put variables in
    ncvar_put(ncnew,airtemp,vals=atonc(mout$Tz,100))
    ncvar_put(ncnew,leaftemp,vals=atonc(mout$tleaf,100))
    ncvar_put(ncnew,relhum,vals=atonc(mout$relhum,1))
    ncvar_put(ncnew,windspeed,vals=atonc(mout$windspeed,100))
    ncvar_put(ncnew,raddir,vals=atonc(mout$raddir,1))
    ncvar_put(ncnew,raddif,vals=atonc(mout$raddif,1))
    ncvar_put(ncnew,radlw,vals=atonc(mout$radlw,1))
    ncatt_put(ncnew,0,"Coordinate reference system",as.character(crs(dtm)))
    ncatt_put(ncnew,0,"Longitude of timezone meridian",merid)
    ncatt_put(ncnew,0,"Time difference from the timezone meridian",dst)
    nc_close(ncnew)
  }
  if (reqhgt == 0) {
    soiltemp<-ncvar_def(name="Soil surface temperature",units="deg C x 100",dim=list(east,north,times),
                        missval=-9999,compression=9,prec="integer")
    soilmoist<-ncvar_def(name="Soil surface moisture",units="Percentage volume",dim=list(east,north,times),
                         missval=-9999,compression=9,prec="integer")
    # Create nc file
    nc.name<-fileout
    ncnew<-nc_create(filename=nc.name,list(soiltemp,soilmoist))
    # Put variables in
    ncvar_put(ncnew,soiltemp,vals=atonc(mout$T0,100))
    ncvar_put(ncnew,leaftemp,vals=atonc(mout$soilm,100))
    ncatt_put(ncnew,0,"Coordinate reference system",as.character(crs(dtm)))
    ncatt_put(ncnew,0,"Longitude of timezone meridian",merid)
    ncatt_put(ncnew,0,"Time difference from the timezone meridian",dst)
    nc_close(ncnew)
  }
  if (reqhgt < 0) {
    tname<-paste0("Soil temperature at depth ",abs(reqhgt)," m")
    soiltemp<-ncvar_def(name=tname,units="deg C x 100",dim=list(east,north,times),
                        missval=-9999,compression=9,prec="integer")
    nc.name<-fileout
    ncnew<-nc_create(filename=nc.name,list(soiltemp))
    # Put variables in
    ncvar_put(ncnew,soiltemp,vals=atonc(mout$Tz,100))
    ncatt_put(ncnew,0,"Coordinate reference system",as.character(crs(dtm)))
    ncatt_put(ncnew,0,"Longitude of timezone meridian",merid)
    ncatt_put(ncnew,0,"Time difference from the timezone meridian",dst)
    nc_close(ncnew)
  }
}
#' runmicro on big areas
#'
#' The function `runmicro_big` tiles larger studies and saves outputs for each tile
#'
#' @param weather a data.frame of weather variables (see details)
#' @param rainfall a vector of daily rainfall
#' @param reqhgt height above ground at which model outputs are needed (m).
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a RasterLayer onject of elevations (see details)
#' @param pathout a file directory to which to save data
#' @param hourly optional logical indicating whether to expand model to hourly and write
#' outputs as ncdf4 files, or keep outputs at daily and save raw outputs
#' @param pai_a an optional array of plant area index values above `reqhgt` (see details)
#' @param xyf optional input for called function [wind()]
#' @param zf optional input for called function [wind()]
#' @param soilinit initial soil moisture fractions in surface and subsurface layer (see [soilmpredict()])
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param surfwet an optional single numeric value of array of values specifying the proportion
#' of the canopy surface that should be treated as wet surface (see details)
#' @param merid optionally, longitude of local time zone meridian (decimal degrees)
#' @param dst optionally, numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for BST if merid = 0).
#'
#' @return if `hourly` = TRUE, ncdf4 files of hourly values for each 100 x 100 grid cell tiles of
#' the study area, numbered by row and column, and saved in a folder `microut` in the
#' directory specified by path.
#' @return if `hourly` = FALSE, daily max and min values saved in files with a .R
#' extension for each 100 x 100 grid cell tiles of the study area, numbered by row
#' and column, and saved in a folder `microut` in the directory specified by path. This
#' saves disk space, is singnificantly faster and files can later be expanded to
#' hourly and written out as ncdf4 files using [expandtonc()]
#'
#' @import raster ncdf4
#' @export
#'
#' @details
#' The parmater `pai_a` is used to calaculate the radiation incercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calculated from total plant area index by
#' assuming leaf density within the canopy is uniformly vertically distributed. If suuplied
#' it must have the same dimensions as micro$pai. I.e. with the same x and y dims as the
#' the supplied dtm and values for each hour as the z dimension. The paramater `surfwet`
#' determines how much of the canopy should be treated as wet surface when calaculating
#' latent heat fluxes. However, except when extremely droughted, the matric potential of leaves
#' is such that `surfwet` ~ 1.
runmicro_big <- function(weather, rainfall, reqhgt, vegp, soilc, dtm, pathout, hourly = FALSE,
                         pai_a = NA, xyf = NA, zf = NA, soilinit = c(NA, NA), tfact = 1.5,
                         surfwet = 1, merid = 0, dst = 0, runchecks = TRUE) {
  # Calculate universal variables
  path2<-paste0(pathout,"microut/")
  dir.create(path2,showWarnings = FALSE)
  cat("Computing terrain variables over whole area\n")
  slr<-terrain(dtm,opt="slope")
  apr<-terrain(dtm,opt="aspect")
  twi<-.topidx(dtm)
  hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) hor[,,i]<-.horizon(dtm,(i-1)*15)
  dsm<-dtm+vegp$hgt
  maxhgt<-max(getValues(vegp$hgt),na.rm=T)
  maxhgt<-ifelse(maxhgt<2,2,maxhgt)
  if (is.na(xyf)) xyf<-trunc(pmax(10/res(dtm)[1],10))
  wsa<-.windsheltera(dsm,8,maxhgt,xyf)
  dms<-dim(dtm)[1:2]
  rws<-ceiling(dms[1]/100)
  cls<-ceiling(dms[2]/100)
  cat(paste0("Running model over ",rws*cls," tiles\n"))
  for (rw in 1:rws) {
    for (cl in 1:cls) {
      dtmi<-.cropraster(dtm,rw,cl)
      v<-getValues(dtmi)
      if (is.na(mean(v,na.rm=T))==F) {
        slri<-.cropraster(slr,rw,cl)
        apri<-.cropraster(apr,rw,cl)
        twii<-.cropraster(twi,rw,cl)
        hori<-.croparray(hor,rw,cl)
        wsai<-.croparray(wsa,rw,cl)
        vegpi<-.vegpcrop(vegp,rw,cl)
        soilci<-.soilccrop(soilc,rw,cl)
        microd<-suppressWarnings(modelin_dy(weather,rainfall,vegpi,soilci,dtmi,merid,dst,runchecks))
        cat(paste0("Running model for tile ",rw," ",cl,"\n"))
        if (hourly) {
          expand = TRUE
        } else expand = FALSE
        mout<-runmicro_dy(microd,reqhgt,expand,NA,xyf,zf,soilinit,tfact,surfwet,
                          slri,apri,hori,wsai,maxhgt,twii)
        rwt<-ifelse(rw<10,paste0("0",rw),paste0("",rw))
        clt<-ifelse(cl<10,paste0("0",cl),paste0("",cl))
        cat(paste0("Writing model outputs for tile ",rw," ",cl,"\n"))
        if (hourly) {
          fo<-paste0(path2,"area_",rwt,"_",clt,".nc")
          writetonc(mout,fo,dtmi,weather,reqhgt, merid, dst)
        } else {
          mout<-.outasint(mout,weather,dtmi,reqhgt,merid,dst)
          fo<-paste0(path2,"area_",rwt,"_",clt,".R")
          save(mout,file=fo)
        }
      }
    }
  }
}
#' expand runmicro_big daily output to hourly and write as ncff4 file
#'
#' @description the function `expandtonc` takes the daily output written by [runmicro_big()]
#' when hourly = F, expands this to hourly, and writes data out as an ncdf4 file.
#'
#' @param filein file to be loaded as written out by [runmicro_big()]
#' @param file to be written out. Must have .nc extension.
#' @import ncdf4 raster
#' @export
expandtonc <- function(filein, fileout) {
  mm<-load(filein)
  mout<-get(mm)
  dtm<-mout$dtm; mout$dtm<-NULL
  weather<-mout$weather; mout$weather<-NULL
  reqhgt<-mout$reqhgt; mout$reqhgt<-NULL
  merid<-mout$merid; mout$merid<-NULL
  dst<-mout$dst; mout$dst<-NULL
  mout_mn<-mout$mout_mn
  mout_mx<-mout$mout_mx
  mout_mn$Tz<-mout_mn$Tz/100
  mout_mx$Tz<-mout_mx$Tz/100
  mout_mn$tleaf<-mout_mn$tleaf/100
  mout_mx$tleaf<-mout_mx$tleaf/100
  mout_mn$windspeed<-mout_mn$windspeed/100
  mout_mx$windspeed<-mout_mx$windspeed/100
  cat("Expanding data to hourly \n")
  if (reqhgt > 0) {
    mout<-.expandclim2(mout_mn,mout_mx,weather)
    class(mout)<-"microout"
  }
  if (reqhgt == 0) {
    climd<-.climtodaily(weather)
    Tz<-.expandtohour(mout_mn$Tz,mout_mx$Tz,climd$smn,climd$smx,weather$temp)
    mout<-list(Tz=Tz,tleaf=NA,soilm=.ehr(mout_mn$theta),
               relhum=NA,windspeed=NA,raddir=NA,raddif=NA,radlw=NA)
    class(mout)<-"microout"
  }
  if (reqhgt < 0) {
    climd<-.climtodaily(weather)
    T0<-.expandtohour(mout_mn$T0,mout_mx$T0,climd$smn,climd$smx,climdata$temp)
    mout<-list(Tz=Tz,tleaf=NA,relhum=NA,windspeed=NA,raddir=NA,raddif=NA,radlw=NA)
    class(mout)<-"microout"
  }
  cat("Writing data as ncdf4 file \n")
  writetonc(mout, fileout, dtm, weather, reqhgt, merid = 0, dst = 0)
}
