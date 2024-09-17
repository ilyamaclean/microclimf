#' Subsets outputs from point microclimate model
#'
#' The function `subsetpointmodel` provides a means of selecting monthly or
#' yearly values from the ouputs of the point microclimate model
#'
#' @param pointmodel a list of model outputs from the point microclimate model as returned by [runpointmodel()].
#' @param tstep one of `year` or `month` (see details)
#' @param what one of `tmax`, `tmin` or `tmedian` (maximum, minimum or median temperature respectively - see details)
#' @param days optionally a vector of the days in the time sequence  to return data for (if provided `tstep` and other inputs are ignored)
#' @param Tc optionally a vector of canopy heat exchange surface temperatures (used to ensure same data subset when applying over arrays)
#' @seealso [runpointmodel()]
#' @return A list with the same format as `pointmodel` but with specified values
#' 'only selected.
#' @details setting 'tstep' to `year` identifies the day in each year with the e.g.
#' the hottest or coldest hourly temperature, and 'tstep' to `month` ideas the day
#' in each month in each year with e.g. the hottest or coldest hourly temperature.
#' Values for all hours of that day are returned, to ensure that the ground heat flux
#' in the grid microclimate model can be estimated. If `what` is set to `tmax` or `tmin`
#' the hottest or coldest hour within each month or year are identified. if If `what`
#' is set to `tmedian` hourly temperatures within the month or year are ranked and
#' the median hour identified.
#' @export
#' @examples
#' # Extract all hourly values for day in which hottest hour in each month occurs
#' micropoint <- runpointmodel(climdata, 0.05, dtmcaerth, vegp, soilc)
#' sub_micropoint_hr <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")
#' Tcanopy <- sub_micropoint_hr$dfo$Tc
#' tme <- as.POSIXct(sub_micropoint_hr$weather$obs_time, tz = "UTC")
#' # Plot
#' plot(Tcanopy ~ tme, type = "l")
subsetpointmodel <- function(pointmodel, tstep = "month", what = "tmax", days = NA, Tc = NA) {
  .extractday<-function(Tc,tme,sel,what) {
    if (what == "tmax") {
      s2<-which.max(Tc[sel])[1]
    } else if (what == "tmin") {
      s2<-which.min(Tc[sel])[1]
    } else if (what == "tmedian") {
      o<-order(Tc[sel])
      n<-trunc(length(o)/2)
      s2<-o[n]
    } else stop ("what must be one of tmax, tmin or tmedian")
    st<-tme[sel[s2]]
    i<-which(tme$year==st$year & tme$mon==st$mon & tme$mday==st$mday)
    return(i)
  }
  dfo<-pointmodel$dfo
  if (class(days) == "logical") {
    tme<-as.POSIXlt(pointmodel$weather$obs_time,tz="UTC")
    yrs<-unique(tme$year)
    if (class(Tc) == "logical") Tc<-dfo$Tc
    if (tstep == "year") {
      sel<-which(tme$year==yrs[1])
      ai<-.extractday(Tc,tme,sel,what)
      if (length(yrs) > 1) {
        for (y in 2:length(yrs)) {
          sel<-which(tme$year==yrs[y])
          i<-.extractday(Tc,tme,sel,what)
          ai<-c(ai,i)
        }
      }
    }
    if (tstep == "month") {
      sely<-which(tme$year==yrs[1])
      mths<-unique(tme$mon[sely])
      sel<-which(tme$mon[sely]==mths[1])
      ai<-.extractday(Tc,tme,sel,what)
      if (length(mths) > 1) {
        for (m in 2:length(mths)) {
          sel<-which(tme$mon[sely]==mths[m])
          i<-.extractday(Tc,tme,sel,what)
          ai<-c(ai,i)
        }
      }
      if (length(yrs) > 1) {
        for (y in 2:length(yrs)) {
          sely<-which(tme$year==yrs[y])
          mths<-unique(tme$mon[sely])
          sel<-which(tme$mon[sely]==mths[1])
          am<-.extractday(Tc,tme,sel,what)
          if (length(mths) > 1) {
            for (m in 2:length(mths)) {
              sel<-which(tme$mon[sely]==mths[m])
              i<-.extractday(Tc,tme,sel,what)
              am<-c(am,i)
            }
          }
          ai<-c(ai,am)
        }
      }
    }
  } else ai<-rep((days-1)*24,each=24)+rep(c(1:24),length(days))
  # Extract data
  microdf<-dfo[ai,]
  weather<-pointmodel$weather
  weatherout<-weather[ai,]
  pointmodel$weather<-weatherout
  pointmodel$dfo<-microdf
  pointmodel$subs<-pointmodel$subs[ai]
  if (class(pointmodel$Tbz) != "logical") pointmodel$Tbz<-pointmodel$Tbz[ai]
  class(pointmodel)<-"micropoint"
  return(pointmodel)
}
#' @title Subsets outputs from point microclimate model run over arrays
#' @description The function `subsetpointmodela` is the equivalent of [subsetpointmodel()] used
#' for outputs from [`runpointmodela()]
#' @param pointmodela a list of model outputs from the point microclimate model applied over arrays as returned by [runpointmodela()].
#' @param tstep one of `year` or `month` (see details)
#' @param what one of `tmax`, `tmin` or `tmedian` (maximum, minimum or median temperature repsectively - see details)
#' @param days optionally a vector of the days in the time sequence  to return data for (if provided tstep is ignored)
#' @seealso [subsetpointmodel()], [runpointmodela()]
#' @return A list with the same format as `pointmodela` but with specified values
#' 'only selected.
#' @export
subsetpointmodela <- function(pointmodela, tstep = "month", what = "tmax", days = NA) {
  # Calculate mean Tc
  Tc<-0
  n<-0
  for (i in 1:length(pointmodela)) {
    if (class(pointmodela[[i]]) != "logical") {
      pma<-pointmodela[[i]]
      Tc<-Tc+pma$dfo$Tc
      n<-n+1
    }
  }
  Tc<-Tc/n
  micropointa<-list()
  for (i in 1:length(pointmodela)) {
    if (class(pointmodela[[i]]) != "logical") {
      micropointa[[i]]<-subsetpointmodel(pointmodela[[i]],tstep,what,days,Tc=Tc)
    } else micropointa[[i]]<-NA
  }
  micropointa
}
#' @title Subsets outputs of snow model
#' @description The function `subsetsnowmodel` subsets snow model outputs as returned
#' by [runsnowmodel()]
#' @param snowmod a list of model outputs from the snow model as returned by [runsnowmodel()].
#' @param subs a vector of index values indicating which hours to select form the snow model
#' @export
#' @examples
#' # Run full snow model (takes ~90 seconds)
#' climdata$temp <- climdata$temp - 8 # Make it colder so there is snow
#' micropoint <- runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc) # Make it colder so there is snow
#' smod <- runsnowmodel(climdata, micropoint, vegp, soilc, dtmcaerth)
#' # Subset snow model using subset point model
#' smicropoint <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")
#' smods <- subsetsnowmodel(smod, smicropoint$subs)
subsetsnowmodel <- function(snowmod, subs) {
  snowmods<-list()
  snowmods$Tc<-snowmod$Tc[,,subs]
  snowmods$Tg<-snowmod$Tg[,,subs]
  snowmods$groundsnowdepth<-snowmod$groundsnowdepth[,,subs]
  snowmods$totalSWE<-snowmod$totalSWE[,,subs]
  snowmods$snowden<-snowmod$snowden[,,subs]
  if (is.array(snowmods$umu)) {
    snowmods$umu<-snowmod$umu[,,subs]
  } else snowmods$umu<-snowmod$umu[subs]
  return(snowmods)
}


#' Check format and values in model inputs
#'
#' The function `checkinputs` checks all the inputs into the model and returns errors
#' and warnings
#'
#' @param weather a data.frame of weather variables (see details)
#' @param precip a vector of daily precipitation
#' @param vegp an object of class vegparams as returned by [vegpfromhab()] (see details)
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a SpatRaster object of elevations (see details)
#' @param daily optional logical indicating whether input weather data are daily
#' @param windhgt height above ground (m) of wind speed measurement
#' @param tstep one of `hour` or `day` indicating whether data are hourly or daily
#' @return a list of checked inputs
#'
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. The array of Plant Area index values in `vegp` must
#' of the same x and y dims as `dtm` but can contain any number of repeated
#' measures up to the number of entries in `weather`. Data are interpolated to the
#' time increment of `weather`. Other vegetation paramaters, including vegetation
#' height are assumed time-invarient. The SpatRaster datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined.
#'
#' @export
#' @import terra
#'
#' @examples
#' # No warnings or errors given:
#' checks<-checkinputs(climdata, vegp, soilc, dtmcaerth)
#' # Warning given (not run)
#' # weather<-climdata
#' # weather$relhum[1]<-101
#' # checks<-checkinputs(weather, vegp, soilc, dtmcaerth)
#' # Error given (NB not run)
#' # weather<-climdata
#' # weather$pres<-weather$pres*1000
#' # checks<-checkinputs(weather, vegp, soilc, dtmcaerth)
#' # Error given for vegp (not run)
#' # vegp2<-vegp
#' # vegp2$clump<-1
#' # checks<-checkinputs(climdata, vegp2, soilc, dtmcaerth)
#' # Warning given for vegp (not run)
#' # vegp2<-vegp
#' # vegp2$pai<-vegp$pai*10
#' # checks<-checkinputs(climdata, vegp2, soilc, dtmcaerth)
checkinputs <- function(weather, vegp, soilc, dtm, windhgt = 2) {
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
  check.unpack<-function(r,nme="r") {
    if (class(r)[1] == "PackedSpatRaster" & class(r)[1] == "SpatRaster") {
      stop(paste0(nme," must be a SpatRaster produced using the terra package"))
    }
    if (class(r)[1] == "PackedSpatRaster") r<-rast(r)
    r
  }
  dtm<-check.unpack(dtm,"dtm")
  # check names
  nms<-names(weather)
  check.names(nms,"obs_time")
  check.names(nms,"temp")
  check.names(nms,"relhum")
  check.names(nms,"pres")
  check.names(nms,"swdown")
  check.names(nms,"difrad")
  check.names(nms,"lwdown")
  check.names(nms,"windspeed")
  check.names(nms,"winddir")
  check.names(nms,"precip")
  # check for NAS
  cd<-weather[,2:10]
  s<-which(is.na(cd))
  if (length(s) > 0) stop("weather contains NAs")
  # check timezone
  tz<-attr(weather$obs_time,"tzone")
  if (tz != "UTC") stop("timezone of obs_time in weather must be UTC")
  # get si
  if (is.na(crs(dtm))) stop("dtm must have a coordinate reference system specified")
  ll<-.latlongfromraster(dtm)
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  # check tme
  sel<-which(is.na(tme))
  if (length(sel)>0) stop("Cannot recognise all obs_time in weather")
  jd<-.jday(tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  sa<-.solalt(lt,ll$lat,ll$long,jd)
  ze<-90-sa
  si<-cos((90-sa)*(pi/180))
  si[si<0]<-0
  # Calculate pressure lims based on elevation
  mxelev<-max(.is(dtm),na.rm=T)
  mnelev<-min(.is(dtm),na.rm=T)
  mxp<-108.5*((293-0.0065*mnelev)/293)^5.26
  mnp<-87*((293-0.0065*mnelev)/293)^5.26
  # Check weather values
  check.vals(weather$temp,-50,65,"temperature","deg C")
  weather$relhum<-up.lim(weather$relhum,100,"relative humidity")
  check.vals(weather$relhum,0,100,"relative humidity","percentage (0-100)")
  check.mean(weather$relhum,5,100,"relative humidity","percentage (0-100)")
  check.vals(weather$pres,mnp,mxp,"pressure","kPa ~101.3")
  check.vals(weather$swdown,0,1350,"shortwave radiation","W / m^2")
  check.vals(weather$difrad,0,1350,"diffuse radiation","W / m^2")
  check.vals(weather$lwdown,0,600,"longwave radiation","W / m^2")
  check.vals(weather$windspeed,0,100,"wind speed","m/s")
  if (windhgt != 2) {
    ws<-(weather$windspeed*4.87)/log(67.8*windhgt-5.42)
  } else ws<-weather$windspeed
  if (max(ws)>30) {
    txt<-paste0("Maximum wind speed seems quite high. Check units are m/s and for ",windhgt," m above ground")
    warning(txt)
  }
  # Check direct radiation at low solar angles
  dirr<-weather$swdown-weather$difrad
  sel<-which(dirr<0)
  if (length(sel)>0) {
    weather$difrad[sel]<-weather$swdown[sel]
    warning("Diffuse radiation values higher than shortwave radiation, and so was set to shortwave radiation values")
  }
  # Calculate clear-sky radiation
  csr<-with(weather,.clearskyrad(tme,ll$lat,ll$long,temp,relhum,pres))
  csd<-csr+50
  sel<-which(dirr>csr)
  if (length(sel)>0) {
    warning("Direct radiation values higher than expected clear-sky radiation values. Assigning excess as diffuse radiation")
    dif<-ceiling((dirr[sel]-csr[sel])*100)/100
    weather$difrad[sel]<-round(weather$difrad[sel]+dif,3)
  }
  sel<-which(weather$swdown>csd)
  if (length(sel)>0) {
    warning("Short wave radiation values significantly higher than expected clear-sky radiation values")
  }
  # Wind direction
  mn<-min(weather$winddir)
  mx<-max(weather$winddir)
  if (mn<0 | mx>360) {
    weather$winddir<-weather$winddir%%360
    warning("wind direction adjusted to range 0-360 using modulo operation")
  }
  # Check other variables
  hrs<-dim(weather)[1]
  dys<-hrs
  if(dys != floor(dys)) stop ("weather needs to include data for entire days (24 hours)")
  # Check vegp data
  xy<-dim(dtm)[1:2]
  if (dim(vegp$pai)[1] != xy[1]) stop("y dimension of vegp$pai does not match dtm")
  if (dim(vegp$pai)[2] != xy[2]) stop("x dimension of vegp$pai does not match dtm")
  vegp$hgt<-check.unpack(vegp$hgt,"vegp$hgt")
  vegp$x<-check.unpack(vegp$x,"vegp$x")
  vegp$gsmax<-check.unpack(vegp$gsmax,"vegp$gsmax")
  vegp$leafr<-check.unpack(vegp$leafr,"vegp$leafr")
  vegp$leafd<-check.unpack(vegp$leafd,"vegp$leafd")
  vegp$leaft<-check.unpack(vegp$leaft,"vegp$leaft")
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
  if (dim(vegp$leaft)[1] != xy[1]) stop("y dimension of vegp$leaft does not match dtm")
  if (dim(vegp$leaft)[2] != xy[2]) stop("x dimension of vegp$leaft does not match dtm")
  if (dim(vegp$x)[3]>1) stop("time variant vegp$x not supported")
  if (dim(vegp$gsmax)[3]>1) stop("time variant vegp$gsmax not supported")
  if (dim(vegp$leafr)[3]>1) stop("time variant vegp$leafr not supported")
  if (dim(vegp$leafd)[3]>1) stop("time variant vegp$leafd not supported")
  if (length(vegp$clump) > 1) {
    if (length(vegp$clump) != length(vegp$pai)) stop("clump must be a single numeric value or have the same dimensions as vegp$pai")
  }
  # check leaf transmittance and reflectance not > 1
  totref<-suppressWarnings(.is(vegp$leafr+vegp$leaft))
  if (max(totref,na.rm=T) > 1) stop("leaf reflectance + transmittance cannot be greater than one")
  # Check soil data
  soilc$soiltype<-check.unpack(soilc$soiltype,"soilc$soiltype")
  soilc$groundr<-check.unpack(soilc$groundr,"soilc$groundr")
  xx<-unique(as.vector(soilc$soiltype))
  if (max(xx,na.rm=T) > 11) stop("Unrecognised soil type")
  if (min(xx,na.rm=T) < 1) stop("Unrecognised soil type")
  if (dim(soilc$soiltype)[1] != xy[1]) stop("y dimension of soilc$soiltype does not match dtm")
  if (dim(soilc$soiltype)[2] != xy[2]) stop("x dimension of soilc$soiltype does not match dtm")
  if (dim(soilc$groundr)[1] != xy[1]) stop("y dimension of soilc$groundr does not match dtm")
  if (dim(soilc$groundr)[2] != xy[2]) stop("x dimension of soilc$groundr does not match dtm")
  if (dim(soilc$soiltype)[3]>1) stop("time variant soilc$soiltype not supported")
  if (dim(soilc$groundr)[3]>1) stop("time variant soilc$groundr not supported")
  xx<-as.vector(soilc$groundr)
  xx<-xx[is.na(xx)==F]
  check.vals(xx,0,1,"soil reflectivity","range 0 to 1")
  xx<-as.vector(vegp$leafr)
  xx<-xx[is.na(xx)==F]
  check.vals(xx,0,1,"leaf reflectivity","range 0 to 1")
  xx<-as.vector(vegp$clump)
  xx<-xx[is.na(xx)==F]
  check.vals(xx,0,1,"vegetation clumping factor","range 0 to 1")
  xx<-as.vector(vegp$leafd)
  xx<-xx[is.na(xx)==F]
  check.vals(xx,0,5,"leaf diamater","metres")
  if (mean(xx)>1) warning(paste0("Mean leaf diameter of ",mean(xx)," seems large. Check units are in metres"))
  xx<-as.vector(vegp$gsmax)
  xx<-xx[is.na(xx)==F]
  check.vals(xx,0,2,"maximum stomatal conductance","mol / m^2 /s")
  # PAI
  if (class(vegp$pai)[1] == "PackedSpatRaster") vegp$pai<-rast(vegp$pai)
  xx<-.is(vegp$pai)
  xx<-xx[is.na(xx)==F]
  if (min(xx)<0) stop("Minimum vegp$pai must be greater than or equal to zero")
  if (max(xx)>15) warning(paste0("Maximum vegp$pai of ",max(xx)," seems high"))
  return(list(weather=weather,precip=precip,vegp=vegp,soilc=soilc))
}
.PAIforayear <- function(habitat, lat, long, yr) {
  laigaus <- function(minlai, maxlai, pkday, dhalf, yr) {
    diy <- 365
    sdev <- 0.0082 * dhalf^2 + 0.0717 * dhalf + 13.285
    difv <- maxlai - minlai
    x<-c(-diy:diy)
    y <- 1 / (sdev * sqrt(2 * pi)) * exp(-0.5 * (((x - 0) / sdev) ^ 2))
    y[(diy + ceiling(0.5 * diy)):(2 * diy + 1)] <- y[(diy - ceiling(0.5 * diy)):diy]
    st <- diy + 1 - pkday
    y <- y[st:(st + diy - 1)]
    x <- c(1:diy)
    x <- c(0, x, c(366:375))
    y <- c(y[diy], y, y[1:10])
    sel <-c(0:15) * 25 + 1
    x<-x[sel]
    y<-y[sel]
    tme <- as.POSIXct((x * 24 * 3600), origin = paste0(yr - 1,"-12-31 12:00"), tz = "GMT")
    xy <- stats::spline(tme, y, n = diy * 24 + 241)
    tme2 <- as.POSIXlt(xy$x, origin = "1970-01-01 00:00", tz = "GMT")
    sel <- which(tme2$year + 1900 == yr)
    y <- xy$y[sel]
    dify <- max(y) - min(y)
    y <- y * (difv / dify)
    y <- y + minlai - min(y)
    return(y)
  }
  diy<-365
  if (yr%%4 == 0) diy<-366
  long<- ifelse(long > 180.9375, long - 360, long)
  long<- ifelse(long < -179.0625, long + 360, long)
  lat<- ifelse(lat< -89.49406, -89.49406, lat)
  lat<- ifelse(lat> 89.49406, 89.49406, lat)
  ll<-vect(cbind(long,lat))
  mmonth <-c(16, 45.5, 75, 105.5, 136, 166.5, 197, 228, 258.5, 289, 319.5, 350)
  e <- ext(c(-179.0625, 180.9375, -89.49406, 89.49406))
  clim <- rep(NA,5)
  for (i in 1:5) {
    r <- rast(globclim[,,i])
    ext(r) <- e
    clim[i] <- extract(r, ll)$lyr.1
  }
  wgts <- function(x1, x2, ll, lmn, lmx) {
    ll <- ifelse(ll < lmn, lmn, lat)
    ll <- ifelse(ll > lmx, lmx, lat)
    w <- 1 - (abs(ll - lmn)  / (abs(ll - lmn)  + abs(ll - lmx)))
    y <- w * x1 + (1 - w) * x2
    y
  }
  # By habitat type
  if (habitat == 1) { # Evergreen needleleaf forest
    h2 <- 74.02 + 5.35 * clim[1]
    h1 <-  203.22 - 35.63 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 50, 50, hperiod)
    p2 <- 216.71 - 2.65 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <- 2.33 + 0.0132 * clim[1]
    } else maxlai <- 2.33 + 0.0132 * 20
    minlai <- 1.01
  }
  if (habitat == 2) { # Evergreen broadleaf forest
    hperiod <-  154.505 + 2.040 * clim[1]
    hperiod <- ifelse(hperiod < 50, 50, hperiod)
    peakdoy <- peakdoy <- mmonth[round(clim[5], 0)]
    maxlai <- 1.83 + 0.22 * log(clim[3])
    minlai <- (-1.09) + 0.4030 * log(clim[3])
    minlai <- ifelse(minlai < 1, 1, minlai)
  }
  if (habitat == 3) { # Deciduous needleleaf forest
    h2 <- 51.18 + 3.77  * clim[1]
    h1 <- 152
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    p2 <- 204.97 - 1.08 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <-  2.62 + 0.05 * clim[1]
    } else maxlai <- 2.62 + 0.05 * 20
    minlai <- 0.39
  }
  if (habitat == 4) { # Deciduous broadleaf forest
    h2 <- 47.6380 + 2.9232 * clim[1]
    h1 <- 220.06 - 79.19 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 32.5, 32.5, hperiod)
    p2 <- 209.760 - 1.208 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <- 3.98795 + 0.03330 * clim[1]
    } else maxlai <- 3.98795 * 0.03330 * 20
    minlai <- 0.4808
  }
  if (habitat == 5) { # Mixed forest
    h2 <- 74.02 + 5.35 * clim[1]
    h1 <-  203.22 - 35.63 * clim[4]
    hperiod1 <- wgts(h1, h2, abs(lat), 0, 20)
    h2 <- 51.18 +  3.77  * clim[1]
    h1 <-  152
    hperiod2 <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- (hperiod1 + hperiod2) / 2
    hperiod <- ifelse(hperiod < 30.5, 30.5, hperiod)
    p2 <- 216.71 - 2.65 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy1 <- wgts(p1, p2, abs(lat), 0, 30)
    p2 <-  204.97 - -1.08 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    peakdoy2 <- wgts(p1, p2, abs(lat), 0, 30)
    peakdoy <- (peakdoy1 + peakdoy2) / 2
    if (clim[1] <= 20) {
      maxlai1 <- 2.33 + 0.0132 * clim[1]
      maxlai2 <-  2.62 + 0.05 * clim[1]
      maxlai <- (maxlai1 + maxlai2) / 2
    } else maxlai <- 3.107
    minlai <- 0.7
  }
  if (habitat == 6) { # Closed shrublands
    h2 <- 33.867 + 6.324 * clim[1]
    h1 <-  284.20 - 102.51 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 30.5, 30.5, hperiod)
    p2 <- 223.55 - 3.125 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    maxlai <- 2.34
    minlai <- -0.4790 + 0.1450 * log(clim[3])
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
  }
  if (habitat == 7) { # Open shrublands
    h2 <- 8.908 + 4.907 * clim[1]
    h1 <-  210.09 - 28.62 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 38.3, 38.3, hperiod)
    p2 <- 211.7 - 4.085 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    maxlai <- -0.7206 + 0.272 * log(clim[3])
    minlai <- -0.146 +  0.059 * log(clim[3])
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
  }
  if (habitat == 8) { # Woody savannas
    hperiod1 <-  47.6380 + 2.9232 * clim[1]
    hperiod1 <- ifelse(hperiod1 < 32.5, 32.5, hperiod1)
    hperiod2 <- 71.72 + 3.012 * clim[1]
    h1 <- (hperiod1 + hperiod2) / 2
    h2 <- 282.04 - 92.28 * clim[4]
    h2 <- ifelse(hperiod1 < 31.9, 31.9, hperiod1)
    hperiod <- wgts(h1, h2, abs(lat), 25, 35)
    peakdoy1 <- 209.760 - 1.208 * clim[1]
    peakdoy1 <- ifelse(peakdoy1 > 244, 244, peakdoy1)
    if (lat < 0)  peakdoy1 <- ( peakdoy1 + diy / 2)%%diy
    peakdoy2 <- 211.98 - 3.4371 * clim[1]
    peakdoy2 <- ifelse(peakdoy2 > 244, 244, peakdoy2)
    if (lat < 0)  peakdoy2 <- (peakdoy2 + diy / 2)%%diy
    p2 <- (peakdoy1 + peakdoy2) / 2
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 40)
    if (clim[1] <= 20) {
      maxlai1 <- 3.98795 + 0.03330 * clim[1]
      maxlai2 <- 1.0532 * 0.016 * clim[1]
    } else {
      maxlai1 <- 3.98795 + 0.03330 * 20
      maxlai2 <- 1.0532 * 0.016 * 20
    }
    mx2 <- (maxlai1 + maxlai2) / 2
    minlai1 <- 0.4808
    minlai2 <- 0.0725 * 0.011 * clim[1]
    mn2 <- (minlai1 + minlai2) / 2
    mx1 <- 1.298 + 0.171 * log(clim[3])
    mn1 <- -2.9458 + 0.5889 * log(clim[3])
    maxlai <- wgts(mx1, mx2, abs(lat), 10, 40)
    minlai <- wgts(mn1, mn2, abs(lat), 10, 40)
    minlai <- ifelse(minlai < 0.0362, 0.0362, minlai)
  }
  if (habitat == 9 | habitat == 10 | habitat == 11) {  # Grasslands
    h2 <- 71.72 + 3.012 * clim[1]
    h1 <- 269.22 -  89.79 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 31.9, 31.9, hperiod)
    p2 <- 211.98 - 3.4371 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    mx2 <- 1.48
    mn2 <- 0.0725 * 0.011 * clim[1]
    mx1 <- 0.1215 + 0.2662 * log(clim[3])
    mn1 <- 0.331 + 0.0575 * log(clim[3])
    maxlai <- wgts(mx1, mx2, abs(lat), 10, 40)
    minlai <- wgts(mn1, mn2, abs(lat), 10, 40)
    minlai <- ifelse(minlai < 0.762, 0.762, minlai)
  }
  if (habitat == 12) {  # Permanent wetlands
    h2 <- 76 + 4.617 * clim[1]
    h1 <- 246.68 - 66.82 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 40, 40, hperiod)
    p2 <- 219.64 - 2.793 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    maxlai <- -0.1782 + 0.2608 * log(clim[3])
    maxlai <- ifelse(maxlai < 1.12, 1.12, maxlai)
    minlai <-  -0.1450 + 0.1440 * log(clim[3])
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
  }
  if (habitat == 13) { # Croplands
    h2 <- 54.893 +  1.785 * clim[1]
    h1 <- 243 - 112.18 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 10, 30)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 212.95 - 5.627 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <- 3.124 - 0.0886 * clim[1]
    } else maxlai <- 3.124 - 0.0886 * 20
    maxlai <- ifelse(maxlai > 3.14, 3.14, maxlai)
    minlai <- 0.13
  }
  if (habitat == 14) { # Urban and built-up
    h2 <- 66.669 +  5.618 * clim[1]
    h1 <- 283.44 - 86.11 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 215.998 - 4.2806 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 30)
    if (clim[1] <= 20) {
      maxlai <- 1.135 - 0.0244 * clim[1]
    } else maxlai <- 1.135 - 0.0244 * 20
    maxlai <- ifelse(maxlai > 1.15, 1.15, maxlai)
    minlai <- 0.28
  }
  if (habitat == 15) { # Cropland/Natural vegetation mosaic
    h2 <- 29.490 +  8.260 * clim[1]
    h1 <- 326.46 - 161.70 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 10, 30)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 210.867 - 3.5464 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 30)
    if (clim[1] <= 20) {
      maxlai <- 3.5485 - 0.09481 * clim[1]
    } else maxlai <- 3.5485 - 0.09481 * 20
    maxlai <- ifelse(maxlai > 3.14, 3.14, maxlai)
    if (clim[1] <= 20) {
      minlai <- -0.072815 - 0.044546 * clim[1]
    } else minlai <- -0.072815 - 0.044546 * 20
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
  }
  if (habitat == 16) { # Barren or sparsely vegetated
    h2 <- 80.557 +  6.440 * clim[1]
    h1 <- 344.65 -  -191.94 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 236.0143 - 3.4726 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 30)
    maxlai <- -0.05491 + 0.05991 * log(clim[4])
    maxlai <- ifelse(maxlai < 0.81, 0.81, maxlai)
    minlai <- 0.08
  }
  lai <- laigaus(minlai, maxlai, peakdoy, hperiod, 2000)
  if (habitat ==  "Short grasslands" | habitat == 10) lai <- lai / 2
  yhr<-round(mmonth*24,0)
  lai<-lai[yhr]
  return(lai)
}
.onehab<-function(habitat) {
  if (habitat == 1) {  # Evergreen needleleaf forest
    hgt<-15 # Vegetation height
    x<-0.4  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.25 # Leaf reflectivity (shortwave)
    leafd<-0.01 # Leaf width (m)
  }
  if (habitat == 2) { # Evergreen broadleaf forest
    hgt<-20 # Vegetation height
    x<-1.2  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.35 # Leaf width (m)
  }
  if (habitat == 3) { # Deciduous needleleaf forest
    hgt<-10 # Vegetation height
    x<-0.4  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.01  # Leaf width (m)
  }
  if (habitat == 4) { # Deciduous broadleaf forest
    hgt<-15 # Vegetation height
    x<-1.2  # Campbell x
    gsmax<-0.23 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.07 # Leaf width (m)
  }
  if (habitat == 5) { # Mixed forest
    hgt<-10 # Vegetation height
    x<-0.8  # Campbell x
    gsmax<-0.28 # Maximum stomatal conductance
    leafr<-0.28 # Leaf reflectivity (shortwave)
    leafd<-0.04 # Leaf width (m)
  }
  if (habitat == 6) { # Closed shrublands
    hgt<-2 # Vegetation height
    x<-1  # Campbell x
    gsmax<-0.35 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.04 # Leaf width (m)
  }
  if (habitat == 7) { # Open shrublands
    hgt<-1.5 # Vegetation height
    x<-0.7  # Campbell x
    gsmax<-0.35 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.04 # Leaf width (m)
  }
  if (habitat == 8) { # Woody savannas
    hgt<-3 # Vegetation height
    x<-0.7  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.35 # Leaf reflectivity (shortwave)
    leafd<-0.03 # Leaf width (m)
  }
  if (habitat == 9) { #Savannas
    hgt<-1.5 # Vegetation height
    x<-0.15  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.35 # Leaf reflectivity (shortwave)
    leafd<-0.01 # Leaf width (m)
  }
  if (habitat == 10) { # Short grasslands
    hgt<-0.25 # Vegetation height
    x<-0.15  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.35 # Leaf reflectivity (shortwave)
    leafd<-0.01 # Leaf width (m)
  }
  if (habitat == 11) { # Tall grasslands
    hgt<-1.5 # Vegetation height
    x<-0.15  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.35 # Leaf reflectivity (shortwave)
    leafd<-0.01 # Leaf width (m)
  }
  if (habitat == 12) { # Permanent wetlands
    hgt<-0.5 # Vegetation height
    x<-1.4  # Campbell x
    gsmax<-0.55 # Maximum stomatal conductance
    leafr<-0.5 # Leaf reflectivity (shortwave)
    leafd<-0.09 # Leaf width (m)
  }
  if (habitat == 13) { # Croplands
    hgt<-0.5 # Vegetation height
    x<-0.2  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.02 # Leaf width (m)
  }
  if (habitat == 14) { # Urban and built-up
    hgt<-1.5 # Vegetation height
    x<-1  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.04 # Leaf width (m)
  }
  if (habitat == 15) { # Cropland/Natural vegetation mosaic
    hgt<-1 # Vegetation height
    x<-0.5  # Campbell x
    gsmax<-0.3 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.03 # Leaf width (m)
  }
  if (habitat == 16) { # Barren or sparsely vegetated
    hgt<-0.15 # Vegetation height
    x<-0.6  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.015 # Leaf width (m)
  }
  return(list(hgt=hgt,x=x,gsmax=gsmax,leafr=leafr,leafd=leafd))
}
.paifromhabitat <- function(habitat, lat, long, tme) {
  yr<-unique(tme$year+1900)
  pai<-0
  for (i in yr) {
    lai<-.PAIforayear(habitat, lat, long, i)
    sel<-which(tme$year+1900==i)
    tme2<-tme[sel]
    mth<-unique(tme2$mon+1)
    pai1<-lai[mth]
    pai<-c(pai,pai1)
  }
  return(pai[-1])
}
#' Generate vegpetation paramaters form habitat type
#'
#' The function `vegpfromhab` generates an object of class vegparams from
#' a SpatRaster object of habitat types
#'
#' @param habitats a SpatRaster object of habitat types expressed as integers (see details)
#' @param hgts an optional SpatRaster object of vegetation heights. Estimated from habitat type if not provided.
#' @param pai an optional array of plant area index values. Estimated at monthly intervals
#' from habitat type, with seasonal variation dtermined form location and date if not provided.
#' @param lat latitude in decimal degrees. Only needed if `pai` not provided.
#' @param long longitude in decimal degrees. Only needed if `pai` not provided.
#' @param tme POSIXlt object of dates. Only needed if `pai` not provided (see details).
#' @param clump0 optional logical, which if TRUE sets the canopy clumping factor to 0, and if false, estimates it using [clumpestimate()]
#' @return an object of class vegparams - a list with the following elements:
#' @return `pai` an array of monthly plant area index values (see details).
#' @return `hgt` a SpatRaster object if vegetation heights (m)
#' @return `x` a SpatRaster object of ratios of vertical to horizontal projections of leaf foliage
#' @return `gsmax` a SpatRaster object of maximum stomatal conductances (mol / m^2 / s)
#' @return `leafr` a SpatRaster object of leaf reflectance values (to shortwave radiation)
#' @return `clump` an array of monthly values indicating the degree of canopy clumpiness, by default set to 0 (vegetation not clumped)
#' @return `leafd` a SpatRaster object of mean leaf widths (m)
#' @return `leaft` a SpatRaster object of mean leaf transmittance (m)
#'
#' @details
#' This function estimates the vegetation parameters needed to run microclimf from habitat types.
#' Plant area index values represent the combined one sided woody and green vegetation
#' plant area per unit ground area. If not provided, then approximated
#' from habitat type, location and date. The procedure is based on calibration
#' against MODIS-derived estimates, accounting for regional climate. An inbuilt dataset
#' of regional rainfall and temperature is included with the package. Monthly
#' values for each unique month in `tme` are returned. Note that values are
#' assumed spatially constant across a given habitat type, which is unlikely
#' to be the case in reality. Likewise, if vegetation height values are not provided,
#' these are estimated form habitat type and assumed constant within that habitat type.
#' Habitat types should be expressed as integers as follows:
#' (1) for Evergreen needleleaf forest,
#' (2) for Evergreen broadleaf forest,
#' (3) for Deciduous needleleaf forest,
#' (4) for Deciduous broadleaf forest,
#' (5) for Mixed forest,
#' (6) for Closed shrubland,
#' (7) for Open shrubland,
#' (8) for Woody savanna,
#' (9) for Savanna,
#' (10) for Short grassland,
#' (11) for Tall grassland,
#' (12) for Permanent wetland,
#' (13) for Cropland,
#' (14) for Urban and built-up,
#' (15) for Cropland / Natural vegetation mosaic and
#' (16) for Barren or sparsely vegetated

#' @export
#' @import terra sp
#'
#' @examples
#' library(terra)
#' tme<-as.POSIXlt(c(0:8783)*3600,origin="2000-01-01 00:00", tz = "GMT")
#' veg<-vegpfromhab(habitats,lat=50,long=-5,tme=tme)
#' plot(rast(veg$pai[,,1]), main = "Jan PAI")
#' plot(veg$hgt, main = "Vegetation height")
#' plot(veg$x, main = "Leaf angle coefficient")
#' plot(veg$gsmax, main = "Maximum stomatal conductance")
#' plot(veg$leafr, main = "Leaf reflectance")
#' plot(rast(veg$clump[,,1]), main = "Canopy clumping factor")
#' plot(veg$leafd, main = "Leaf diameter")
#' plot(veg$leaft, main = "Leaf transmittance")
vegpfromhab <- function(habitats, hgts = NA, pai = NA, lat, long, tme, clump0 = TRUE) {
  .poparray<-function(a,sel,v) {
    for (i in 1:length(v)) {
      m<-a[,,i]
      m[sel]<-v[i]
      a[,,i]<-m
    }
    a
  }
  if (class(habitats)[1] == "PackedSpatRaster") habitats<-rast(habitats)
  # unique habitats
  m<-.is(habitats)
  uh<-unique(as.vector(m))
  uh<-uh[is.na(uh)==F]
  # pai test
  pte<-base::mean(pai,na.rm=T)
  # Create blank array for pai
  if (is.na(pte)) {
    paii<-.paifromhabitat(1, lat, long, tme)
    pai<-array(NA,dim=c(dim(m),length(paii)))
  }
  # Create blank rasters
  x<-m; gsmax<-m; leafr<-m; leafd<-m; hgt<-m
  for (i in uh) {
    sel<-which(m==i)
    if (is.na(pte)) {
      paii<-.paifromhabitat(i, lat, long, tme)
      pai<-.poparray(pai,sel,paii)
    }
    vegi<-.onehab(i)
    x[sel]<-vegi$x
    gsmax[sel]<-vegi$gsmax
    leafr[sel]<-vegi$leafr
    leafd[sel]<-vegi$leafd
    hgt[sel]<-vegi$hgt
  }
  leaft<-0.5*leafr
  clump<-pai*0
  if (clump0 == F) {
    for (i in 1:dim(pai)[3]) clump[,,i]<-clumpestimate(hgt, leafd, pai[,,i])
  }
  # Convert to rasters
  if (class(hgts)=="logical") {
    hgt<-.rast(hgt,habitats)
  } else hgt<-hgts
  x<-.rast(x,habitats)
  gsmax<-.rast(gsmax,habitats)
  leafr<-.rast(leafr,habitats)
  leaft<-.rast(leaft,habitats)
  leafd<-.rast(leafd,habitats)
  vegp<-list(pai=pai,hgt=hgt,x=x,gsmax=gsmax,leafr=leafr,clump=clump,leafd=leafd,leaft=leaft)
  class(vegp)<-"vegparams"
  return(vegp)
}
#' Estimate canopy clumpiness
#'
#' The function `clumpestimate` estimates how clumpy a canopy is based on height and
#' plant area index
#'
#' @param hgt vegetation height (m).
#' @param pai plant area index value(s).
#' @param leafd mean leaf width (m).
#' @param maxclump optional parameter indicating maximum clumpiness
#' @return `clump` parameter indicating the fraction of radiation passing directly through larger gaps in the canopy.
#' @export
#' @details `hgt`, `pai` and `leafd` must all be single numeric values or arrays with identical dimensions
clumpestimate <- function(hgt, leafd, pai, maxclump = 0.95) {
  pai[pai>1]<-1
  sel<-which(leafd>hgt)
  leafd[sel]<-hgt[sel]
  clump<-(1-pai)^(hgt/leafd)
  sel<-which(clump>maxclump)
  clump[sel]<-maxclump
  clump
}
#' Derives leaf and ground reflectance from albedo
#' @param pai a SpatRaster of plant area index values
#' @param x a SpatRaster of the ratio of vertical to horizontal projections of leaf foliage
#' @param alb a SpatRaster of white-sky albedo
#' @param ltrr an optional numeric value giving an approximate estimate of the ratio of leaf transmittance to leaf reflectance (e.g. value of 1 makes leaf transmittance equal to reflectance). See details
#' @return a list of the following
#' \describe{
#'   \item{leafr}{Leaf reflectance (range 0 - 1)}
#'   \item{leaft}{Leaf transmittance (range 0 - 1)}
#'   \item{gref}{Ground reflectance (range 0 - 1)}
#' }
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib microclimf, .registration = TRUE
#' @export
#' @details the microclimate model is not unduly sensitive to `lttr` so if unknown, an approximate
#' value or the default can be used.
#' @examples
#' pai <- rast(vegp$pai)[[9]] # Plant Area Index in Sep (month in which albedo image was flown)
#' x <- rast(vegp$x)
#' alb <- rast(albedo)
#' lgr <- leafrfromalb(pai, x, alb)
#' plot(lgr$leafr)
#' plot(lgr$gref)
leafrfromalb<-function(pai, x, alb, ltrr = 0.5) {
  .Solveforlref <- function(pai,albin,x=1,gref=0.15,ltrr=0.5) {
    uniroot(leafrcpp, c(0.0001,0.9999),pai,gref,albin,x,ltrr,f.lower=-1,f.upper=1)$root
  }
  .Solveforlref2<-function(pai,albin,x,gref,ltrr) {
    out<-tryCatch(.Solveforlref(pai,albin,x,gref,ltrr),error=function(cond) -999)
  }
  cat("Computing ground reflectance \n")
  # calculate ground reflectance
  canc<-1-exp(-.is(pai))
  gref<-(.is(alb)-canc*0.6)/(1-canc)
  s<-which(canc>0.25)
  gref[s]<-NA
  gref<-.fillr(.rast(gref,pai),pai)
  gref[gref<0.01]<-0.01
  gref<-mask(gref,pai)
  # Convert to matricesz
  lai<-.is(pai)
  gref<-.is(gref)
  x<-.is(x)
  albi<-.is(alb)
  lref<-array(NA,dim=dim(lai))
  cat("Computing leaf reflectance \n")
  # Calculate leaf reflectance
  nn<-dim(gref)[1]
  pb <- txtProgressBar(min = 0, max = nn, style = 3)
  for (i in 1:dim(gref)[1]) {
    for (j in 1:dim(gref)[2]) {
      if (lai[i,j]>0 & is.na(lai[i,j]) == FALSE) {
        lref[i,j]<-.Solveforlref2(lai[i,j],albi[i,j],x[i,j],gref[i,j],ltrr)
      }
    }
    setTxtProgressBar(pb, i)
  }
  lref<-.rast(lref,alb)
  me<-mean(as.vector(leafr),na.rm=TRUE)
  lref[is.na(lref)]<-me
  mx<-0.95/(1+ltrr)
  lref[lref>mx]<-mx
  lref<-mask(lref,pai)
  return(list(leafr=lref,leaft=ltrr*lref,gref=.rast(gref,pai)))
}
#' Create climate arrays for inputting to microclimate model from a netCDF4 file
#'
#' @description The function `nctoclimarray` converts data in a netCDF4 file returned
#' by [mcera5::request_era5()] to the correct formal required for running [modelina()] or
#' [modelina_dy()].
#'
#' @param ncfile character vector containing the path and filename of the nc file
#' @param dtm a SpatRaster object of elevations covering the extent of the study area (see details)
#' @param dtr_cor_fac numeric value to be used in the diurnal temperature range
#' correction of coastal grid cells. Default = 1.285, based on calibration against UK Met Office
#' observations. If set to zero, no correction is applied.
#' @return a list of the following:
#' \describe{
#'   \item{tme}{POSIXlt object of times corresponding to climate observations}
#'   \item{climarray}{a list of arrays of hourly weather variables - same format as for [modelina()].}
#'   \item{precarray}{an array of daily precipitation values - same format as for [modelina()].}
#'   \item{dtmc}{a coarse resolution digital elevation dataset matching the resolution of input
#'   climate data, but with a coordinate reference system and extent matching `dtm`}
#' }
#' @export
#' @import ncdf4
#' @details the model requires that input climate data are projected using a coordinate reference
#' system in which x and y are in metres. Since values returned by [mcera5::request_era5()]
#' are in lat long, the output data are reprojected using the coordinate reference system and
#' extent of dtm (but retain the approximate original grid resolution of the input climate data).
#' Returned climate data match the resolution, coordinate reference system and extent of `dtmc`.
nctoclimarray <- function(ncfile, dtm, dtr_cor_fac = 1.285)  {
  t2m<-rast(ncfile,subds = "t2m") # Air temperature (K)
  d2m<-rast(ncfile,subds = "d2m") # Dewpoint temperature (K)
  sp<-rast(ncfile,subds = "sp") # Surface pressure (Pa)
  u10<-rast(ncfile,subds = "u10") # U-wind at 10m (m/s)
  v10<-rast(ncfile,subds = "v10") # V-wind at 10m (m/s)
  tp<-rast(ncfile,subds = "tp") # Total precipitation (m)
  msnlwrf<-rast(ncfile,subds = "msnlwrf")   # Mean surface net long-wave radiation flux (W/m^2)
  msdwlwrf<-rast(ncfile,subds = "msdwlwrf") # Mean surface downward long-wave radiation flux (W/m^2)
  fdir<-rast(ncfile,subds = "fdir") #  Total sky direct solar radiation at surface (W/m^2)
  ssrd<-rast(ncfile,subds = "ssrd") # Surface short-wave (solar) radiation downwards (W/m^2)
  lsm<-rast(ncfile,subds = "lsm") # Land sea mask
  # Create coarse-resolution dtm to use as template for resampling
  te<-terra::project(t2m[[1]],crs(dtm))
  agf<-res(te)[1]/res(dtm)[1]
  dtmc<-aggregate(dtm,fact=agf,fun=mean,na.rm=T)
  # Apply coastal correction to temperature data
  tmn<-.ehr(.hourtoday(as.array(t2m)-273.15,min))
  mu<-(1-as.array(lsm))*dtr_cor_fac+1
  tc<-.rast(((as.array(t2m)-273.15)-tmn)*mu+tmn,t2m)
  # Calculate vapour pressure
  ea<-.rast(.satvap(as.array(d2m)-273.15),t2m)
  # Resample all variables to match dtmc
  tc<-terra::project(tc,dtmc)
  ea<-terra::project(ea,dtmc)
  sp<-terra::project(sp,dtmc)
  u10<-terra::project(u10,dtmc)
  v10<-terra::project(v10,dtmc)
  tp<-terra::project(tp,dtmc)
  msnlwrf<-terra::project(msnlwrf,dtmc)
  msdwlwrf<-terra::project(msdwlwrf,dtmc)
  fdir<-terra::project(fdir,dtmc)
  ssrd<-terra::project(ssrd,dtmc)
  # Derive varies
  temp<-as.array(tc) # Temperature (deg c)
  relhum<-(as.array(ea)/.satvap(temp))*100 # Relative humidity (%)
  pres<-as.array(sp)/1000  # Surface pressure (kPa)
  swrad<-as.array(ssrd)/3600 # Downward shortwave radiation (W/m^2)
  difrad<-swrad-as.array(fdir)/3600 # Downward diffuse radiation (W/m^2)
  skyem<-as.array(msdwlwrf)/as.array(msdwlwrf-msnlwrf) # sky emissivity (0-1)
  windspeed<-sqrt(as.array(u10)^2+as.array(v10)^2)*log(67.8*2-5.42)/log(67.8*10-5.42) # Wind speed (m/s)
  winddir<-as.array((terra::atan2(u10,v10)*180/pi+180)%%360) # Wind direction (deg from N - from)
  precarray<-.hourtoday(as.array(tp)*1000,sum)
  # Save lists
  climarray<-list(temp=temp,relhum=relhum,pres=pres,swrad=swrad,difrad=difrad,
                  skyem=skyem,windspeed=windspeed,winddir=winddir)
  # Generate POSIXlt object of times
  tme<-as.POSIXlt(time(t2m),tz="UTC")
  # Output for returning
  out<-list(tme=tme,climarray=climarray,precarray=precarray,dtmc=dtmc)
  return(out)
}
#' Writes model output as ncdf4 file
#'
#' @description The function `writetonc` writes hourly model outputs as an ncdf4 file to a specified file
#'
#' @param mout an object of class microut as returned by runmicro_hr or runmicro_dy with expand = TRUE
#' @param fileout output filename
#' @param dtm a SpatRast covering the extent of the model outputs
#' @param reqhgt at at which model was run (m)
#' @import ncdf4 terra
#' @export
writetonc <- function(mout, fileout, dtm, reqhgt) {
  atonc<-function(a,rd) {
    a<-apply(a,c(2,3),rev)
    a<-aperm(a,c(2,1,3))
    a <- round(a*rd,0)
    a <-array(as.integer(a),dim=dim(a))
    a
  }
  # Create eastings and nothings sequence
  if (class(dtm)[1] == "PackedSpatRaster") dtm<-rast(dtm)
  est<-seq(ext(dtm)$xmin+res(dtm)[1]/2,ext(dtm)$xmax-res(dtm)[1]/2,res(dtm)[1])
  nth<-seq(ext(dtm)$ymin+res(dtm)[2]/2,ext(dtm)$ymax-res(dtm)[2]/2,res(dtm)[2])
  east<-ncdim_def(name="east",units="metres",vals=est,longname="Eastings")
  north<-ncdim_def(name="north",units="metres",vals=nth,longname="Northings")
  # Create time variable
  tme<-as.POSIXct(mout$tme)
  times<-ncdim_def(name="Time",units="Decimal hours since 1970-01-01 00:00",vals=as.numeric(tme)/3600)
  # Variable names
  if (reqhgt > 0) {
    tname<-paste0("Air temperature at height ",reqhgt," m")
    lname<-paste0("Leaf temperature at height ",reqhgt," m")
    rname<-paste0("Relative humidity at height ",reqhgt," m")
    wname<-paste0("Wind speed at height ",reqhgt," m")
    # Define variables
    airtemp<-ncvar_def(name="Tz",longname=tname,units="deg C x 100",dim=list(east,north,times),
                       missval=-9999,compression=9,prec="integer")
    leaftemp<-ncvar_def(name="tleaf",longname=lname,units="deg C x 100",dim=list(east,north,times),
                        missval=-9999,compression=9,prec="integer")
    relhum<-ncvar_def(name="relhum",longname=rname,units="Percentage",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    windspeed<-ncvar_def(name="windspeed",longname=wname,units="m/s x 100",dim=list(east,north,times),
                         missval=-9999,compression=9,prec="integer")
    raddir<-ncvar_def(name="Rdirdown",longname="Downward direct shortwave radiation",units="W/m^2",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    raddif<-ncvar_def(name="Rdifdown",longname="Downward diffuse shortwave radiation",units="W/m^2",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    radlw<-ncvar_def(name="Rlwdown",longname="Downward longwave radiation",units="W/m^2",dim=list(east,north,times),
                     missval=-9999,compression=9,prec="integer")
    radusw<-ncvar_def(name="Rswup",longname="Upward shortwave radiation",units="W/m^2",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    radulw<-ncvar_def(name="Rlwup",longname="Upward longwave radiation",units="W/m^2",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    # Create nc file
    nc.name<-fileout
    ncnew<-nc_create(filename=nc.name,list(airtemp,leaftemp,relhum,windspeed,raddir,raddif,radlw,radusw,radulw))
    # Put variables in
    ncvar_put(ncnew,airtemp,vals=atonc(mout$Tz,100))
    ncvar_put(ncnew,leaftemp,vals=atonc(mout$tleaf,100))
    ncvar_put(ncnew,relhum,vals=atonc(mout$relhum,1))
    ncvar_put(ncnew,windspeed,vals=atonc(mout$windspeed,100))
    ncvar_put(ncnew,raddir,vals=atonc(mout$Rdirdown,1))
    ncvar_put(ncnew,raddif,vals=atonc(mout$Rdifdown,1))
    ncvar_put(ncnew,radlw,vals=atonc(mout$Rlwdown,1))
    ncvar_put(ncnew,radusw,vals=atonc(mout$Rswup,1))
    ncvar_put(ncnew,radulw,vals=atonc(mout$Rlwup,1))
    ncatt_put(ncnew,0,"Coordinate reference system",as.character(crs(dtm)))
    nc_close(ncnew)
  }
  if (reqhgt == 0) {
    soiltemp<-ncvar_def(name="Tz",longname="Soil surface temperature",units="deg C x 100",dim=list(east,north,times),
                        missval=-9999,compression=9,prec="integer")
    soilmoist<-ncvar_def(name="soilm",longname="Soil surface moisture",units="Volume percentage soil moisture in top 10 cm of soil",
                         dim=list(east,north,times),
                         missval=-9999,compression=9,prec="integer")
    raddir<-ncvar_def(name="Rdirdown",longname="Downward direct shortwave radiation",units="W/m^2",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    raddif<-ncvar_def(name="Rdifdown",longname="Downward diffuse shortwave radiation",units="W/m^2",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    radlw<-ncvar_def(name="Rlwdown",longname="Downward longwave radiation",units="W/m^2",dim=list(east,north,times),
                     missval=-9999,compression=9,prec="integer")
    radusw<-ncvar_def(name="Rswup",longname="Upward shortwave radiation",units="W/m^2",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    radulw<-ncvar_def(name="Rlwup",longname="Upward longwave radiation",units="W/m^2",dim=list(east,north,times),
                      missval=-9999,compression=9,prec="integer")
    # Create nc file
    nc.name<-fileout
    ncnew<-nc_create(filename=nc.name,list(soiltemp,soilmoist,raddir,raddif,radlw,radusw,radulw))
    # Put variables in
    ncvar_put(ncnew,soiltemp,vals=atonc(mout$Tz,100))
    ncvar_put(ncnew,soilmoist,vals=atonc(mout$soilm,100))
    ncvar_put(ncnew,raddir,vals=atonc(mout$Rdirdown,1))
    ncvar_put(ncnew,raddif,vals=atonc(mout$Rdifdown,1))
    ncvar_put(ncnew,radlw,vals=atonc(mout$Rlwdown,1))
    ncvar_put(ncnew,radusw,vals=atonc(mout$Rswup,1))
    ncvar_put(ncnew,radulw,vals=atonc(mout$Rlwup,1))
    ncatt_put(ncnew,0,"Coordinate reference system",as.character(crs(dtm)))
    nc_close(ncnew)
  }
  if (reqhgt < 0) {
    tname<-paste0("Soil temperature at depth ",abs(reqhgt)," m")
    soiltemp<-ncvar_def(name="Tz",longname=tname,units="deg C x 100",dim=list(east,north,times),
                        missval=-9999,compression=9,prec="integer")
    soilmoist<-ncvar_def(name="soilm",longname="Soil surface moisture",units="Percentage volume",dim=list(east,north,times),
                         missval=-9999,compression=9,prec="integer")
    nc.name<-fileout
    ncnew<-nc_create(filename=nc.name,list(soiltemp,soilmoist))
    # Put variables in
    ncvar_put(ncnew,soiltemp,vals=atonc(mout$Tz,100))
    ncvar_put(ncnew,soilmoist,vals=atonc(mout$soilm,100))
    ncatt_put(ncnew,0,"Coordinate reference system",as.character(crs(dtm)))
    nc_close(ncnew)
  }
}
#' @title Mosaics a list of overlapping SpatRasters blending overlap areas
#' @description Mosaics a list of overlapping SpatRasters blending
#' the areas of overlap using a distance weighting to eliminate tiling effects
#' @param rlist a list of SpatRasters
#' @details
#' If rlist contains SpatRasters that are not
#' overlapping the conventional terra::moasic function is used.
#' If rlist contains SpatRasters that do overlap, they should comprise
#' a list of adjacent rasters in a single row or column.
#' @import terra
#' @export
mosaicblend <- function(rlist) {
  # order by row and then by column
  xmn<-0
  ymn<-0
  for (i in 1:length(rlist)) {
    e<-ext(rlist[[i]])
    xmn[i]<-e$xmin
    ymn[i]<-e$ymin
  }
  xmn1<-unique(xmn)
  ymn1<-unique(ymn)
  le<-min(length(xmn1),length(ymn1))
  if (le > 1) warning("rlist not a row or column. Blended mosaicing may not work")
  if (length(xmn1) > length(ymn1)) {
    o<-order(xmn)
  } else o<-order(ymn)
  rlist2<-list()
  for (i in 1:length(o)) rlist2[[i]]<-rlist[[o[i]]]
  rlist<-NULL
  rma<-rlist2[[1]]
  for (i in 2:length(rlist2)) {
    r<-rlist2[[i]]
    it<-intersect(ext(rma),ext(r))
    a<-as.numeric((it$xmax-it$xmin)*(it$ymax-it$ymin))
    if (a>0) {
      rma<-.blendmosaic(rma, r)
    } else rma<-mosaic(rma,r)
  }
  return(rma)
}


