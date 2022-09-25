#' Predicts soil moisture from rainfall and radiation
#'
#' @description The function `soilmpredict` implements a simple two layer
#' soil model to predict the daily soil moisture fraction in the surface and
#' sub-surface layers
#'
#' @param rainfall a vector of daily rainfall (mm/day)
#' @param rnet a vector of hourly net radiation (W/m^2)
#' @param soiltype a soil type (as listed in [soilparameters()])
#' @param soilinit initial soil moisture fractions in surface and sub-surface layers
#' @param soilmcoefs an optional list of coefficients for running the soil moisture model
#' as returned by [fitsoilm()]. Uses inbuilt parameters if not supplied.
#' @return a list of the following:
#' @return `soilm1` soil moisture fraction in surface layer
#' @return `soilm2` soil moisture fraction in sub-surface layer
#' @importFrom graphics par
#' @export
#' @examples
#' # Calculate net radiation from inbuilt climate dataset
#' swrad <- (1 - 0.15) * climdata$swrad
#' lwout <- 5.67e-8 * 0.95 * (climdata$temp + 273.15)^4
#' lwnet <- (1 - climdata$skyem) * lwout
#' rnet <- swrad - lwnet
#' # Run soil moisture model using inbuilt rainfall dataset
#' sm<-soilmpredict(rainfall, rnet, "Loam", c(0.4, 0.4))
#' # Plot results
#' plot(sm$soilm1, type="l", ylim=c(0.1, 0.5), xlab="Day",
#'      ylab="Soil moisture fraction")
#' par(new = TRUE)
#' plot(sm$soilm2, type="l", ylim=c(0.1,0.5), col="blue", lwd=2,
#' xlab="", ylab="")
soilmpredict <- function(rainfall, rnet, soiltype, soilinit=c(NA, NA), soilmcoefs = NA) {
  # get soil coefficients
  ii<-which(soilparameters$Soil.type==soiltype)
  if (is.na(soilmcoefs[1])) {
    soilmcoefs<-list(mult=soilparameters$mult[ii],
                     rmu=soilparameters$rmu[ii],
                     a=soilparameters$a[ii],
                     pwr=soilparameters$pwr[ii])
  }
  # Get net radiation
  rnet[rnet<0]<-0
  rnetd<-matrix(rnet,ncol=24,byrow=T)
  rnetd<-apply(rnetd,1,mean)
  if (is.na(soilinit[1])) {
    soilinit[1]<-(soilparameters$Smax[ii]+soilparameters$Smin[ii])/2
  }
  if (is.na(soilinit[2])) {
    soilinit[2]<-(soilparameters$Smax[ii]+soilparameters$Smin[ii])/2
  }
  s<-soilinit[1]
  s2<-soilinit[2]
  for (i in 2:length(rnetd)) {
    s[i]<-s[i-1]+soilmcoefs$rmu*rainfall[i]-soilmcoefs$mult*rnetd[i]
    sav<-(s[i-1]+s2[i-1])/2
    k<-soilparameters$Ksat[ii]*(sav/soilparameters$Smax[ii])^soilmcoefs$pwr
    dif<-s2[i-1]-s[i-1]
    s[i]<-s[i]+soilmcoefs$a*k*dif
    s2[i]<-s2[i-1]-((soilmcoefs$a*k*dif)/10)
    s[i]<-ifelse(s[i]>soilparameters$Smax[ii],soilparameters$Smax[ii],s[i])
    s[i]<-ifelse(s[i]<soilparameters$Smin[ii],soilparameters$Smin[ii],s[i])
    s2[i]<-ifelse(s2[i]>soilparameters$Smax[ii],soilparameters$Smax[ii],s2[i])
    s2[i]<-ifelse(s2[i]<soilparameters$Smin[ii],soilparameters$Smin[ii],s2[i])
  }
  return(list(soilm1=s,soilm2=s2))
}

#' Spatially distribute soil moisture by topographic wetness index
#'
#' The function `soilmdistribute` spatially distrubutes soil moisture by the Bevan and Kirkby
#' topographic wetness index
#'
#' @param soilm a vector of soil moisture fractions
#' @param a digital terrain model (used for calaculating topographic wetness)
#' @param Smin minimum fractional soil moisture content
#' @param Smax maximum fractional soil moisture content
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
#' @return a 3D array of soil moistures with the same x and y dims as `dtm` and z
#' equivelent to length(soilm)
#' @import terra
#' @export
#' @examples
#' # Calculate vector of soil moistures
#' library(terra)
#' swrad <- (1 - 0.15) * climdata$swrad
#' lwout <- 5.67e-8 * 0.95 * (climdata$temp + 273.15)^4
#' lwnet <- (1 - climdata$skyem) * lwout
#' rnet <- swrad - lwnet
#' soilm <- soilmpredict(rainfall, rnet, "Loam", c(0.4, 0.4))$soilm1
#' sma<- soilmdistribute(soilm, dtmcaerth)
#' plot(rast(sma[,,365])) # soil moisture on last day of year
#' plot(rast(sma[,,188])) # soil moisture on driest day of year
#'
soilmdistribute <- function(soilm, dtm, Smin = 0.074, Smax = 0.422, tfact = 1.5, twi = NA) {
  if (class(dtm)[1] == "PackedSpatRaster") dtm<-rast(dtm)
  if (class(twi)[1] == "logical") twi<-.topidx(dtm)
  rge<-Smax-Smin
  soilm<-ifelse(soilm<=Smin,Smin+0.001,soilm)
  soilm<-ifelse(soilm>=Smax,Smax-0.001,soilm)
  theta<-(soilm-Smin)/rge
  lt<-log(theta/(1-theta))
  lt<-.vta(lt,twi)
  ltwi<-log(.is(twi))/tfact
  me<-mean(ltwi,na.rm=T)
  add<-ltwi-me
  add<-.rta(rast(add),length(soilm))
  smout<-lt+add
  smout<-1/(1+exp(-smout))
  smout<-smout*rge+Smin
  smout
}
#' Applies two-stream radiation model
#'
#' @description The function `twostream` applies a varient of the Sellers two-stream
#' radiation model to calaculate long and showtwave radiation absorbed by the ground and
#' canopy and shortwave radiation absorbed by leafs.
#'
#' @param micro an object of class `micro` as returned by [modelin()]
#' @param reqhgt height at which temperatures are required (m). Negative if below ground
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied, but outer cells will be NA.
#' @return an object of class micro with additional terms added for subsequent modelling
#' @import terra
#' @export
#'
twostream<-function(micro, reqhgt = 0.05, pai_a = NA, slr = NA, apr = NA, hor = NA) {
  dtm<-micro$dtm
  dtm[is.na(dtm)]<-0
  # === (1a) Calculate solar altitude
  tme<-micro$tme
  jd<-.jday(tme=tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  salt<-.solalt(lt,micro$lat,micro$long,jd,micro$merid,micro$dst)
  sazi<-.solazi(lt,micro$lat,micro$long,jd,micro$merid,micro$dst)
  alt<-.vta(.ar(salt),dtm)
  # === (1b) Calculate solar index
  si<-.solarindex(dtm,alt,sazi,slr,apr)
  # === (1c) Calculate horizon angles
  if (class(hor)[1]=="logical") {
    hor<-array(NA,dim=c(dim(dtm)[1:2],24))
    for (i in 1:24) hor[,,i]<-.horizon(dtm,(i-1)*15)
  }
  i<-round(sazi/15,0)+1; i[i==25]<-1
  hora<-hor[,,i]
  # === (1d) Calculate terrain shading
  shadowmask<-hora*0+1
  shadowmask[hora>tan(alt)]<-0
  si<-si*shadowmask
  # === (1e) Calculate sky view
  msl<-tan(apply(atan(hor),c(1,2),mean))
  svf<-0.5*cos(2 *msl)+0.5
  # === (1f) Calculate canopy k
  kkd<-.cank(micro$vegx,alt,si)
  Kc<-kkd$kd/kkd$k0
  # === (1g) Calculate pai a
  micro$vha<-.rta(micro$veghgt,length(jd))
  if (class(pai_a)[1]=="logical") {
    if (reqhgt > 0) {
      fd<-foliageden(reqhgt,micro$vha,micro$pai)
      micro$pai_a<-fd$pai_a
      micro$leafden<-fd$leafden
    } else micro$pai_a<-0
  }
  # === (1g) Adjust paramaters for gap fraction and inclined surface
  n<-reqhgt/micro$vha
  n[n<1]<-1
  pai_a<-micro$pai_a/(1-micro$clump^n)
  pai_t<-micro$pai/(1-micro$clump)
  sk<-si/sin(alt)
  sk[sk>1]<-1
  sk[sk==0]<-1
  gref2<-micro$gref/sk
  # === (1f) Calculate two stream base parameters
  om<-micro$lref+micro$ltra
  a<-1-om
  del<-micro$lref-micro$ltra
  mla<-(9.65*(3+micro$vegx)^(-1.65))
  mla[mla>pi/2]<-pi/2
  J<-cos(mla)^2
  gma<-0.5*(om+J*del)
  s<-0.5*(om+J*del/kkd$k)*kkd$k
  sstr<-om*kkd$k-s
  # === (1g) Calculate two stream base parameters
  h<-sqrt(a^2+2*a*gma)
  sig<-kkd$kd^2+gma^2-(a+gma)^2
  S1<-exp(-h*pai_t)
  S2<-exp(-kkd$kd*pai_t)
  u1<-a+gma*(1-1/micro$gref)
  u2<-a+gma*(1-micro$gref)
  D1<-(a+gma+h)*(u1-h)*1/S1-(a+gma-h)*(u1+h)*S1
  D2<-(u2+h)*1/S1-(u2-h)*S1
  # === (1h) Calculate Diffuse radiation parameters
  p1<-(gma/(D1*S1))*(u1-h)
  p2<-(-gma*S1/D1)*(u1+h)
  p3<-(1/(D2*S1))*(u2+h)
  p4<-(-S1/D2)*(u2- h)
  # === (1i) Calculate Direct radiation parameters
  u1<-a+gma*(1-1/gref2)
  u2<-a+gma*(1-gref2)
  D1<-(a+gma+h)*(u1-h)*1/S1-(a+gma-h)*(u1+h)*S1
  D2<-(u2+h)*1/S1-(u2-h)*S1
  p5<- -s*(a+gma-kkd$kd)-gma*sstr
  v1<-s-(p5*(a+gma+kkd$kd))/sig
  v2<-s-gma-(p5/sig)*(u1+kkd$kd)
  p6<-(1/D1)*((v1/S1)*(u1-h)-(a+gma-h)*S2*v2)
  p7<-(-1/D1)*((v1*S1)*(u1+h)-(a+gma+h)*S2*v2)
  p8<-sstr*(a+gma+kkd$kd)-gma*s
  v3<-(sstr+gma*gref2-(p8/sig)*(u2-kkd$kd))*S2
  p9<-(-1/D2)*((p8/(sig*S1))*(u2+h)+v3)
  p10<-(1/D2)*(((p8*S1)/sig)*(u2-h)+v3)
  # === (1j) Calculate albedo
  albd<-p1+p2+micro$clump^2*micro$gref
  sel<-which(is.na(albd))
  albd[sel]<-micro$gref[sel]
  albb<-p5/sig+p6+p7+micro$clump^Kc*gref2
  sel<-which(is.na(albb))
  albb[sel]<-gref2[sel]
  albb[albb>0.95]<-0.95
  albedo<-(micro$dirr*sin(alt)*albb+micro$difr*albd)/(micro$dirr*si+micro$difr)
  albedo[is.na(albedo)]<-0.23
  albedo[albedo>0.99]<-0.99
  albedo[albedo<0.01]<-0.01
  # === (1k) Calculate ground absorbed radiation
  # Shortwave
  svfa<-.rta(rast(svf),length(jd))
  Rbgm<-(1-micro$clump^Kc)*exp(-kkd$kd*pai_t)+micro$clump^Kc
  Rdbm<-(1-micro$clump^2)*((p8/sig)*exp(-kkd$kd*pai_t)+p9*exp(-h*pai_t)+p10*exp(h*pai_t))
  sel<-which(is.na(Rdbm) | Rdbm<0)
  Rdbm[sel]<-0
  sel<-which(Rdbm>1)
  Rdbm[sel]<-1
  Rddm<-(1-micro$clump^2)*(p3*exp(-h*pai_t)+p4*exp(h*pai_t))+micro$clump^2
  sel<-which(is.na(Rddm))
  Rddm[sel]<-1
  micro$radGsw<-(1-micro$gref)*(micro$dirr*si*Rbgm+micro$dirr*sin(alt)*Rdbm+micro$difr*svfa*Rddm)
  # Longwave
  trd<-(1-micro$clump^2)*exp(-pai_t)+micro$clump^2
  micro$lwout<-0.97*5.67*10^-8*(micro$tc+273.15)^4 # Longwave emitted
  lwsky<-micro$skyem*micro$lwout # Longwave radiation down from sky
  micro$radGlw<-0.97*(trd*svfa*lwsky+(1-trd)*micro$lwout+(1-svfa)*micro$lwout)
  # === (1l) Calculate canopy and ground combined absorbed radiation
  trb<-(1-micro$clump^Kc)*exp(-kkd$kd*pai_t)+micro$clump^Kc
  micro$radCsw<-(1-albedo)*(micro$dirr*sin(alt)*(1-trb)+svfa*micro$difr*(1-trd))+micro$radGsw
  micro$radClw<-svfa*lwsky
  if (reqhgt > 0) {
    # === (1m) Calculate up and downstream at reqhgt
    # Downward
    Rbgm<-(1-micro$clump^(Kc*n))*exp(-kkd$kd*pai_a)+micro$clump^(Kc*n)
    Rdbm<-(1-micro$clump^(2*n))*((p8/sig)*exp(-kkd$kd*pai_a)+p9*exp(-h*pai_a)+p10*exp(h*pai_a))
    sel<-which(is.na(Rdbm) | Rdbm<0)
    Rdbm[sel]<-0
    sel<-which(Rdbm>1)
    Rdbm[sel]<-1
    Rddm<-(1-micro$clump^(2*n))*(p3*exp(-h*pai_a)+p4*exp(h*pai_t))+micro$clump^(2*n)
    sel<-which(is.na(Rddm))
    Rddm[sel]<-1
    micro$Rbdown<-micro$dirr*sin(alt)*Rbgm
    micro$Rddown<-micro$difr*svfa*Rddm+micro$dirr*sin(alt)*Rdbm
    # Upward
    Rdbm<-(1-micro$clump^(2*n))*((p5/sig)*exp(-kkd$kd*pai_a)+p6*exp(-h*pai_a)+p7*exp(h*pai_a))
    sel<-which(is.na(Rdbm) | Rdbm<0)
    Rdbm[sel]<-0
    sel<-which(Rdbm>1)
    Rdbm[sel]<-1
    Rddm<-(1-micro$clump^(2*n))*(p1*exp(-h*pai_a)+p2*exp(h*pai_t))+micro$clump^(2*n)
    sel<-which(is.na(Rddm))
    Rddm[sel]<-1
    micro$Rdup<-micro$difr*svfa*Rddm+micro$dirr*sin(alt)*Rdbm
    # === (1n) Calculate leaf swabs
    a[is.na(a)]<-mean(a,na.rm=T)
    micro$radLsw<-0.5*a*(micro$Rddown+micro$Rdup+kkd$k*sin(alt)*micro$Rbdown)
  }
  micro$k<-kkd$k
  micro$kc<-Kc
  micro$trb<-trb
  micro$trd<-trd
  micro$progress<-1
  # Clean micro
  micro$lat<-NULL
  micro$long<-NULL
  micro$merid<-NULL
  micro$dst<-NULL
  return(micro)
}
#' Downscales wind
#'
#' @description The function `wind` downscales wind speed accounting for
#' vegetation and terrain
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above ground for which wind speeds are wanted. If negative (below ground) wind friction velocity only is returned
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param xyf optional spatial smoothing factor applied in calculation of wind speeds (see details)
#' @param zf optional integer used to calculate how frequently to sample the plant area index array
#' when calculating smoothed wind speed (see details)
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from dtm if not supplied, but outer cells will be NA.
#' @param maxhgt an optional height (m) for which wind speed is needed. Determined
#' from height of tallest vegetation or as 2 m, whichever is greater, if not supplied.
#' @return an object of class microin with the following components added:
#' @return `uf` wind fricion velocities (m/s)
#' @return `uz` optionally wind speed at height `reqhgt` (m/s)
#' @return aditional terms required for subsequent modelling
#' @import terra zoo
#' @export
#' @details In downscaling wind, two processes are accounted for. Firstly the drag effects
#' of vegetation on wind , which ultimately dictate the wind height profile. In
#' so doing, it is necessary to accomodate the fact the wind speed is not just affected by
#' the surface roughness at each point location, but also by vegetation surrounding the
#' location. This is accomodated for by applying `xyf` which effectively smooths though surface roughness
#' using `terra::aggregate` where `xyf` is the aggregation factor. Because surface roughness
#' also depends on the plant area index (PAI), which is time-variant, and is computationally intensive
#' to perform aggregations on PAI at each hourly time increment, `zf` determines the time at
#' what frequency to aggregate PAI. E.g. where `zf` is 100, PAI at every 100th time increment is aggregated.
#' By default 'xyf' is set to 10 m or the pixel resolution, which ever is greater and `zf`
#' is set to aggregate PAI at monthly intervals. A second complication arises in that the maximum
#' vegetation height within the study area may exceed the height at which input wind speeds are
#' measured (by default assumed measured at a height of 2 m). For this reason, the wind
#' speed at the maximum height of vegetation within the study is calculated using a standard height
#' profile adjustment for an open area typical of those at which weather stations are sighted
#' prior to adjusting for the effects of local vegetation and terrain on wind speed. Terrain effects
#' are calculated by applying the topographic shelter coefficient described in Maclean et al
#' (2019) Methods in Ecology and Evolution, 10:280-290.
wind <- function(micro, reqhgt, pai_a = NA, xyf = 1, zf = NA, slr = NA,
                 apr = NA, hor = NA, wsa = NA, maxhgt = NA) {
  .sortuf<-function(ln,uf,uz)  {
    ln[ln<0.0001]<-0.0001
    uf<-(0.4*uz)/ln
    uf[uf<0.037]<-0.037
    sel<-which(uf>uz)
    uf[sel]<-uz[sel]
    uf
  }
  tme<-micro$tme
  ti<-trunc((as.numeric(tme[2])-as.numeric(tme[1]))/3600)
  if (micro$progress<1) {
    micro<-twostream(micro,reqhgt,pai_a,slr,apr,hor)
  }
  # calculate xyf and z factor
  mdm<-trunc(min(dim(micro$dtm)[1:2])/2)
  if (is.na(xyf)) {
    xyf<-trunc(pmax(10/res(micro$dtm)[1],10))
    xyf<-ifelse(xyf>mdm,mdm,xyf)
  }
  if (xyf > mdm) {
    warning(paste0("xyf too large. Setting to ",mdm))
    xyf<-mdm
  }
  if (is.na(zf)) {
    if(ti==1) {
      zf<-30*24
    } else zf<-30
  }
  # Calculate roughness lengths
  if (xyf > 1) {
    hgt<-.smr(micro$veghgt,xyf)
    pai2<-.sma(micro$pai,xyf,zf,crs(micro$dtm))
  } else {
    hgt<-micro$veghgt
    pai2<-micro$pai
  }
  hgt[hgt<0.004*10]<-0.004*10
  pai2[pai2<0.001]<-0.001
  le<-length(tme)
  d<-.zeroplanedis(.rta(hgt,le),pai2)
  d[d<(0.004*6.5)]<-0.004*6.5
  zm<-.roughlength(.rta(hgt,le),pai2,0.004)
  # Calculate wind frictian velocity
  if (class(maxhgt)[1] == "logical") {
    maxhgt<-max(as.vector(micro$veghgt),na.rm=T)
    maxhgt<-ifelse(maxhgt<2,2,maxhgt)
  }
  dsm<-micro$dtm+micro$veghgt
  climdata<-micro$climdata
  uz<-.windshelter(climdata$windspeed,climdata$winddir,dsm,maxhgt,8,xyf,wsa=wsa,cors=crs(micro$dtm))
  ln<-log((maxhgt-d)/zm)
  uf<-.sortuf(ln,uf,uz)
  # Recalculate with diabatic correction
  Hest<-0.5*(micro$radCsw+micro$radClw-micro$lwout)
  dbm<-.diabatic(micro$tc,uf,d,zm,maxhgt,Hest)
  micro$psi_m<-dbm$psi_m
  micro$psi_h<-dbm$psi_h
  ln<-log((maxhgt-d)/zm)+micro$psi_m
  uf<-.sortuf(ln,uf,uz)
  micro$uf<-uf
  micro$d<-d
  micro$zm<-zm
  micro$hgt<-hgt
  micro$maxhgt<-maxhgt
  micro$Hest<-Hest
  if (reqhgt > 0) {
    # Calculate wind speed above canopy
    le<-length(tme)
    mr<-suppressWarnings(log((reqhgt-micro$d)/micro$zm)/log((micro$maxhgt-micro$d)/micro$zm))
    ur<-suppressWarnings((micro$uf/0.4)*(log((reqhgt-micro$d)/micro$zm)+micro$psi_m*mr))
    ur<-.minwind(ur,micro$uf)
    # wind speed below canopy
    selb<-which(reqhgt<.rta(micro$veghgt,le))
    # == Wind speed at top of canopy
    hgt<-.rta(micro$hgt,le)
    mr<-suppressWarnings(log((reqhgt-micro$d)/micro$zm)/log((hgt-micro$d)/micro$zm))
    uh<-suppressWarnings((micro$uf/0.4)*(log(hgt-micro$d)/micro$zm))
    uh<-.minwind(uh,micro$uf)
    # == Calculate attenuation coefficient
    vha<-.rta(micro$veghgt,le)
    vha[vha<0.001]<-0.001
    pai<-micro$pai
    pai[pai<0.0001]<-0.0001
    l_m <- .mixinglength(vha,pai)
    a<-((0.2*pai*vha)/l_m)^0.5
    urb<-micro$uf*exp((a*(reqhgt/vha-1)))
    ur[selb]<-urb[selb]
    # == Check nothing higher than wind at reference height
    cda<-micro$climdata
    u2<-cda$windspeed
    u2<-.vta(u2,micro$dtm)
    sel<-which(ur>u2)
    ur[sel]<-u2[sel]
    micro$uz<-ur
    micro$vha<-vha
    micro$l_m<-l_m
    micro$a<-a
    micro$uh<-uh
  }
  micro$progress<-2
  return(micro)
}
#' Calculates ground surface temperature (hourly)
#'
#' @description The function `soiltemp_hr` estimates ground surface temperature
#'
#' @param micro an object of class `microin` as returned by [modelin()] (see details)
#' @param reqhgt height above ground at which model outputs are needed (m).
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param xyf optional input for called function [wind()]
#' @param zf optional input for called function [wind()]
#' @param soilinit initial soil moisture fractions in surface and subsurface layer (see [soilmpredict()])
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param slr an optional SpatRast object of slope values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRast object of aspect values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from dtm if not supplied, but outer cells will be NA.
#' @param maxhgt an optional height (m) for which wind speed is needed. Determined
#' from height of tallest vegetation or as 2 m, whichever is greater, if not supplied.
#' @param twi optional SpatRast object of topographic wetness index values.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
#' @param soilmcoefs optional list of soil moisture model coefficients as returned by [fitsoilm()]
#' @param soiltcoefs optional list of soil moisture model coefficients as returned by [fitsoilt()]
#' @return an object of class micro with the following components added:
#' @return `T0` ground surface temperature (deg C)
#' @return `theta` surface soil moisture fraction
#' @return additional terms needed for subsequent modelling.
#' @import terra
#' @export
#' @seealso [groundrad()], [soiltemp_dy()]
#' @details if [groundrad()] has not been run to add additional elements to `micro`
#' it is automatically called.
soiltemp_hr  <- function(micro, reqhgt = 0.05, pai_a = NA, xyf = 1, zf = NA, soilinit = c(NA, NA),
                         tfact = 1.5, slr = NA, apr = NA, hor = NA, wsa = NA, maxhgt = NA,
                         twi = NA, soilmcoefs = NA, soiltcoefs = NA) {
  # run two-stream and wind functions if not run
  if (micro$progress<2) {
    micro<-wind(micro,reqhgt,pai_a,xyf,zf,slr,apr,hor,wsa,maxhgt)
  }
  # Calculate net radiation at ground
  rnet<-micro$radGsw+micro$radGlw-micro$lwout
  # get summed radiation
  rsu<-suppressWarnings(apply(micro$radGsw,c(1,2),sum,na.rm=T))
  me<-mean(rsu,na.rm=T)
  rsu<-rsu/me
  mx<-max(rsu,na.rm=T)
  mn<-min(rsu,na.rm=T)
  # Calculate max and min soil moisture
  if (is.na(soilmcoefs[1])) {
    thetamx<-.multisoil(micro,soilinit,tfact,mx,twi)
    thetamn<-.multisoil(micro,soilinit,tfact,mn,twi)
  } else {
    clim<-micro$climdata
    soilc<-micro$soilc
    alb<-mean(.is(soilc$groundr),na.rm=T)
    swrad<-(1-alb)*clim$swrad
    lwout<- 5.67e-8*0.95*(clim$temp+273.15)^4
    lwnet<-(1-clim$skyem)*lwout
    rnet<-swrad-lwnet
    rnet1<-rnet*mx
    rnet2<-rnet*mn
    sm<-soilmpredict(micro$prec,rnet1,"Loam",soilinit,soilmcoefs)
    thetamx<-soilmdistribute(sm$soilm1,micro$dtm)
    sm<-soilmpredict(micro$prec,rnet2,"Loam",soilinit,soilmcoefs)
    thetamn<-soilmdistribute(sm$soilm1,micro$dtm)
  }
  # Adjust soil moisture
  rsua<-.rta(rast(rsu),dim(thetamn)[3])
  theta<-(rsua/mx)*(thetamx-thetamn)+thetamn
  thetam<-apply(theta,c(1,2),mean,na.rm=TRUE)
  theta<-.ehr(theta)
  hiy<-length(micro$tme)
  # Correct theta
  theta[theta<0.0002]<-0.0002
  sm<-log(theta/(1-theta))+8.516993
  # Get paramaters
  if (is.na(soiltcoefs[1])) {
    scfs<-.soilcoefs(micro$soilc)
  } else {
    d<-dim(micro$dtm)[1:2]
    scfs<-list(int=array(soiltcoefs$int,dim=d),
               t1=array(soiltcoefs$t1,dim=d),
               t2=array(soiltcoefs$t2,dim=d),
               t3=array(soiltcoefs$t3,dim=d),
               t4=array(soiltcoefs$t4,dim=d),
               t5=array(soiltcoefs$t5,dim=d),
               t6=array(soiltcoefs$t6,dim=d),
               t7=array(soiltcoefs$t7,dim=d))
  }
  # Calculate wind
  w2<-(micro$uf/0.4)*(log((micro$maxhgt-micro$d)/micro$zm)+micro$psi_m)
  w2<-log(w2+1)
  # Predict soil surface temperature
  T0<-.rta(rast(scfs$int),hiy)+
    .rta(rast(scfs$t1),hiy)*rnet+
    .rta(rast(scfs$t2),hiy)*sm+
    .rta(rast(scfs$t3),hiy)*w2+
    .rta(rast(scfs$t4),hiy)*sm*rnet+
    .rta(rast(scfs$t5),hiy)*sm*w2+
    .rta(rast(scfs$t6),hiy)*rnet*w2+
    .rta(rast(scfs$t7),hiy)*rnet*sm*w2
  T0<-T0+micro$tc
  T0<-.lim(T0,micro$tdew)
  # Calculate soil conductivity and specific heat capacity
  cs<-(2400*micro$rho/2.64+4180*thetam)
  phs<-(micro$rho*(1-thetam)+thetam)*1000
  pcs<-cs*phs # Total specific heat capacity
  frs<-micro$Vm+micro$Vq
  c1<-(0.57+1.73*micro$Vq+0.93*micro$Vm)/(1-0.74*micro$Vq-0.49*micro$Vm)-2.8*frs*(1-frs)
  c2<-1.06*micro$rho*thetam
  c3<-1+2.6*micro$Mc^-0.5
  c4<-0.03+0.7*frs^2
  k<-c1+c2*thetam-(c1-c4)*exp(-(c3*thetam)^4)
  # Calculate damping depth etc
  om<-(2*pi)/(24*3600)
  ka<-k/pcs
  DD<-sqrt((2*ka)/om)
  # Calculate G
  A0<-aperm(apply(T0,c(1,2),.A0f),c(2,3,1))
  r<-rast(A0[,,1])
  cda<-micro$climdata
  t0<-.vta(.t0f(cda$temp),r)
  tt<-(c(1:length(cda$temp))-1)%%24
  tt<-.vta(tt*3600,r)
  G<-(sqrt(2)*A0*.rta(rast(k),hiy)*sin(om*(tt-t0)+(pi/4)))/.rta(rast(DD),hiy)
  # Save elements to micro
  micro$T0<-T0
  micro$theta<-theta
  micro$ka<-ka
  micro$DD<-DD
  micro$G<-G
  micro$progress<-3
  # Clean micro
  micro$prec<-NULL
  micro$rho<-NULL
  micro$Vm<-NULL
  micro$Vq<-NULL
  micro$Mc<-NULL
  micro$soilc<-NULL
  micro$radGsw<-NULL
  micro$radGlw<-NULL
  return(micro)
}
#' Calculates ground surface temperature (daily)
#'
#' @description The function `soiltemp_dy` estimates ground surface temperature
#'
#' @param microd an object of class `microindaily` as returned by [modelin_dy()]
#' @param reqhgt height above ground at which model outputs are needed (m).
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param xyf optional input for called function [wind()]
#' @param zf optional input for called function [wind()]
#' @param soilinit initial soil moisture fractions in surface and subsurface layer (see [soilmpredict()])
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()]).
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' dtm if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from dtm if not supplied, but outer cells will be NA.
#' @param maxhgt an optional height (m) for which wind speed is needed. Determined
#' from height of tallest vegetation or as 2 m, whichever is greater, if not supplied.
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
#' @param soilmcoefs optional list of soil moisture model coefficients as returned by [fitsoilm()]
#' @param soiltcoefs optional list of soil moisture model coefficients as returned by [fitsoilt()]
#' @return an object of class micro with the following components added:
#' @return `T0` ground surface temperature (deg C)
#' @return `theta` surface soil moisture fraction
#' @return additional terms needed for subsequent modelling.
#' @import terra
#' @export
#' @seealso [groundrad()], [soiltemp_dy()]
#' @details if [groundrad()] has not been run to add additional elements to `micro`
#' it is automatically called.
soiltemp_dy  <- function(microd, reqhgt = 0.05, pai_a = NA, xyf = 1, zf = NA, soilinit = c(NA, NA), tfact = 1.5,
                         slr = NA, apr = NA, hor = NA, wsa = NA, maxhgt = NA, twi = NA,
                         soilmcoefs = NA, soiltcoefs = NA) {
  # run two-stream and wind functions if not run
  micro_mn<-microd$micro_mn
  micro_mx<-microd$micro_mx
  if (micro_mn$progress<2) {
    micro_mn<-wind(micro_mn,reqhgt,pai_a,xyf,zf,slr,apr,hor,wsa,maxhgt)
    micro_mx<-wind(micro_mx,reqhgt,pai_a,xyf,zf,slr,apr,hor,wsa,maxhgt)
  }
  # Calculate net radiation at ground
  rnet1<-micro_mn$radGsw+micro_mn$radGlw-micro_mn$lwout
  rnet2<-micro_mx$radGsw+micro_mx$radGlw-micro_mx$lwout
  # get summed radiation
  swabs<-micro_mn$radGsw+micro_mx$radGsw
  rsu<-suppressWarnings(apply(swabs,c(1,2),sum,na.rm=T))
  me<-mean(rsu,na.rm=T)
  rsu<-rsu/me
  mx<-max(rsu,na.rm=T)
  mn<-min(rsu,na.rm=T)
  # Calculate max and min soil moisture
  climd<-micro_mn$climdata
  micro_mn$climdata<-microd$climdata
  if (is.na(soilmcoefs[1])) {
    thetamx<-.multisoil(micro_mn,soilinit,tfact,mx,twi)
    thetamn<-.multisoil(micro_mn,soilinit,tfact,mn,twi)
  } else {
    clim<-micro_mn$climdata
    soilc<-micro_mn$soilc
    alb<-mean(.is(soilc$groundr),na.rm=T)
    swrad<-(1-alb)*clim$swrad
    lwout<- 5.67e-8*0.95*(clim$temp+273.15)^4
    lwnet<-(1-clim$skyem)*lwout
    rnet<-swrad-lwnet
    rnet1<-rnet*mx
    rnet2<-rnet*mn
    sm<-soilmpredict(micro_mn$prec,rnet1,"Loam",soilinit,soilmcoefs)
    thetamx<-soilmdistribute(sm$soilm1,micro_mn$dtm)
    sm<-soilmpredict(micro_mn$prec,rnet2,"Loam",soilinit,soilmcoefs)
    thetamn<-soilmdistribute(sm$soilm1,micro_mn$dtm)
  }
  # Adjust soil moisture
  rsua<-.rta(rast(rsu),dim(thetamn)[3])
  theta<-(rsua/mx)*(thetamx-thetamn)+thetamn
  thetam<-apply(theta,c(1,2),mean,na.rm=TRUE)
  micro_mn$climdata<-climd
  diy<-dim(micro_mn$tc)[3]
  # Correct theta
  theta[theta<0.0002]<-0.0002
  sm<-log(theta/(1-theta))+8.516993
  # Get paramaters
  if (is.na(soiltcoefs[1])) {
    scfs<-.soilcoefs(micro_mn$soilc)
  } else {
    d<-dim(micro_mn$dtm)[1:2]
    scfs<-list(int=array(soiltcoefs$int,dim=d),
               t1=array(soiltcoefs$t1,dim=d),
               t2=array(soiltcoefs$t2,dim=d),
               t3=array(soiltcoefs$t3,dim=d),
               t4=array(soiltcoefs$t3,dim=d),
               t5=array(soiltcoefs$t3,dim=d),
               t6=array(soiltcoefs$t3,dim=d),
               t7=array(soiltcoefs$t3,dim=d))
  }
  # Calculate wind speed
  w2_mn<-(micro_mn$uf/0.4)*(log((micro_mn$maxhgt-micro_mn$d)/micro_mn$zm)+micro_mn$psi_m)
  w2_mx<-(micro_mx$uf/0.4)*(log((micro_mn$maxhgt-micro_mx$d)/micro_mx$zm)+micro_mx$psi_m)
  w2_mn[w2_mn<0]<-0
  w2_mx[w2_mx<0]<-0
  w2_mn<-log(w2_mn+1)
  w2_mx<-log(w2_mx+1)
  # Predict soil surface temperature
  T0mn<-.rta(rast(scfs$int),diy)+
    .rta(rast(scfs$t1),diy)*rnet1+
    .rta(rast(scfs$t2),diy)*sm+
    .rta(rast(scfs$t3),diy)*w2_mn+
    .rta(rast(scfs$t4),diy)*sm*rnet1+
    .rta(rast(scfs$t5),diy)*sm*w2_mn+
    .rta(rast(scfs$t6),diy)*rnet1*w2_mn+
    .rta(rast(scfs$t7),diy)*rnet1*sm*w2_mn
  T0mx<-.rta(rast(scfs$int),diy)+
    .rta(rast(scfs$t1),diy)*rnet2+
    .rta(rast(scfs$t2),diy)*sm+
    .rta(rast(scfs$t3),diy)*w2_mx+
    .rta(rast(scfs$t4),diy)*sm*rnet2+
    .rta(rast(scfs$t5),diy)*sm*w2_mx+
    .rta(rast(scfs$t6),diy)*rnet2*w2_mx+
    .rta(rast(scfs$t7),diy)*rnet2*sm*w2_mx
  T0mn<-T0mn+micro_mn$tc
  T0mx<-T0mx+micro_mx$tc
  T0mn<-.lim(T0mn,micro_mn$tdew)
  T0mx<-.lim(T0mx,micro_mx$tdew)
  # Calculate soil conductivity and specific heat capacity
  cs<-(2400*micro_mn$rho/2.64+4180*thetam)
  phs<-(micro_mn$rho*(1-thetam)+thetam)*1000
  pcs<-cs*phs # Total specific heat capacity
  frs<-micro_mn$Vm+micro_mn$Vq
  c1<-(0.57+1.73*micro_mn$Vq+0.93*micro_mn$Vm)/
    (1-0.74*micro_mn$Vq-0.49*micro_mn$Vm)-2.8*frs*(1-frs)
  c2<-1.06*micro_mn$rho*thetam
  c3<-1+2.6*micro_mn$Mc^-0.5
  c4<-0.03+0.7*frs^2
  k<-c1+c2*thetam-(c1-c4)*exp(-(c3*thetam)^4)
  # Calculate damping depth etc
  om<-(2*pi)/(24*3600)
  ka<-k/pcs
  DD<-sqrt((2*ka)/om)
  # Calculate G
  A0<-(T0mx-T0mn)/2
  r<-rast(A0[,,1])
  cda<-microd$climdata
  t0<-.vta(.t0fd(cda$temp),r)
  tme_mn<-micro_mn$tme
  tme_mx<-micro_mx$tme
  tt_mn<-.vta(tme_mn$hour*3600,r)
  tt_mx<-.vta(tme_mx$hour*3600,r)
  G_mn<-(sqrt(2)*A0*.rta(rast(k),diy)*sin(om*(tt_mn-t0)+(pi/4)))/.rta(rast(DD),diy)
  G_mx<-(sqrt(2)*A0*.rta(rast(k),diy)*sin(om*(tt_mx-t0)+(pi/4)))/.rta(rast(DD),diy)
  # Save elements to micro
  micro_mn$T0<-T0mn
  micro_mn$theta<-theta
  micro_mn$ka<-ka
  micro_mn$DD<-DD
  micro_mn$G<-G_mn
  micro_mn$progress<-5
  micro_mx$T0<-T0mx
  micro_mx$theta<-theta
  micro_mx$ka<-ka
  micro_mx$DD<-DD
  micro_mx$G<-G_mx
  micro_mx$progress<-5
  # Clean micro
  micro_mn$prec<-NULL
  micro_mn$rho<-NULL
  micro_mn$Vm<-NULL
  micro_mn$Vq<-NULL
  micro_mn$Mc<-NULL
  micro_mn$soilc<-NULL
  micro_mn$radGsw<-NULL
  micro_mn$radGlw<-NULL
  micro_mx$prec<-NULL
  micro_mx$rho<-NULL
  micro_mx$Vm<-NULL
  micro_mx$Vq<-NULL
  micro_mx$Mc<-NULL
  micro_mx$soilc<-NULL
  micro_mx$radGsw<-NULL
  micro_mx$radGlw<-NULL
  out<-list(micro_mn=micro_mn,micro_mx=micro_mx,climdata=microd$climdata)
  return(out)
}
#' Calculates Latent heat exchange
#'
#' @description The function `PenMont` applies a varient of the Penman-Monteith
#' equation, with G calculated from ground temperatures, to calculate latent
#' heat exchange
#' @param tc air temperature (deg c)
#' @param pk atmospheric pressure (kPa)
#' @param ea vapour pressure of air (kPa)
#' @param radabs absorbed long and shortwave radiation (W/m2)
#' @param gHa boundary layer conductances for heat (mol/m^2/s)
#' @param gs stomatal conductance (mol/m^2/s)
#' @param G rate of ground heat storage (W/m2)
#' @param es saturated vapour pressure (kPa) at temperature tc (calculated if NA)
#' @param tdew dewpoint temperature (deg C) (calculated if NA, see details)
#' @param surfwet proportion of leaf or canopy acting as wet surface
#' @param allout optional logical indicating whether to return Latent, Sensible and
#' Ground heat flux in addition to canopy temperature
#' @return a list with the following components:
#' @return `tcan` the leaf or effective canopy temperature (deg C)
#' @return `H` the sensible heat flux (w/m^2)
#' @return `HR` hte fractional of net radiation that is sensible heat
#' @return Optionally `L` The latent heat flux (W/m^2)
#' @return Optionally `G` The ground heat flux (W/m^2)
#' @details `Latent` can be used to calaculate either the latent heat of
#' evapotranspiration from individual leaves or bulk surface latent heat of
#' evapotranspiration from the canopy as a whole depending on how `gHa` and `gs`
#' are defined. If for the canopy as a whole, heat storage `G` should also be
#' considered and `L` is the Latent heat per unit leaf area and should thus be
#' multiplied by the total leaf area per unit ground area to derive a total. If
#' for a leaf in isolation, set g0 to zero and T0 to an arbitrary value. The
#' paramater `surfwet` determines how much of the canopy should be treated as
#' wet surface. However, except when extremely droughted, the matric potential of
#' leaves is such that `surfwet` ~ 1.
#' @export
PenMont <- function(tc,pk,ea,radabs,gHa,gs,G,es=NA,tdew=NA,surfwet=1,allout=TRUE) {
  if (is.na(es[1])) es<-.satvap(tc)
  if (is.na(tdew[1])) tdew<-.dewpoint(ea,tc)
  delta<-.delta(tc,radabs)
  gv<-1/(1/gHa+1/gs)
  sb<-5.67*10^-8
  gHr<-gHa+(4*0.97*sb*(tc+273.15)^3)/29.3
  Rem<-0.97*sb*(tc+273.15)^4
  m<-44526*(gv/pk)
  L<-m*(es-ea)
  tcan<-(radabs-Rem-L-G)/(29.3*gHr+m*delta)+tc
  tcan<-.lim(tcan,tdew)
  H<-29.3*gHa*(tcan-tc)
  HR<-H/(radabs-Rem)
  HR[HR<0]<-0
  HR[HR>1]<-1
  if (allout) {
    estl<-.satvap(tcan)*surfwet
    L<-m*(estl-ea)
    L<-.lim(L,0)
    out<-list(tcan=tcan,H=H,L=L,G=G,HR=HR)
  } else out<-list(tcan=tcan,H=H,HR=HR)
  return(out)
}
#' Estimate foliage density and plant area index above z
#'
#' @description The function `foliageden` applies a mirrored gamma distribution
#' to estimate, for a given height `z` below canopy, the foliage density at
#' height z and the plant area index above `z`
#'
#' @param z height above ground (m)
#' @param hgt the height of the canopy (m)
#' @param pai total plant area index of canopy
#' @param shape optional shape parameter for gamma distribution (see details)
#' @param rate optional rate parameter for gamm distribution (see details)
#' @return a list with the following components:
#' @return `leafden` the foliage density (one sided leaf area per m^3) at height z
#' @return `paia` the plant area index above `z`
#' @details a broadly realistic estimate of foliage density is generated
#' using a mirrored gamma distribution (i.e. right (top) rather than left-skewed).
#' in applying the gamma distribution `z / hgt` is re-scaled to the range 0-10.
#' @importFrom stats dgamma pgamma
#' @export
#' @return a list of leaf densities and plant area index values
#'
#' @examples
#' z<-c(1:1000)/100
#' fdp<-foliageden(z,10,10)
#' plot(z~fdp$leafden, xlab = "Foliage density")
#' plot(z~fdp$pai_a, xlab = "PAI above z")
foliageden<-function(z,hgt,pai,shape=1.5,rate=shape/7) {
  # reverse and rescale z
  x<-((hgt-z)/hgt)*10
  # totAL DENS
  td<-pgamma(10,shape,rate) # total density
  rfd<-dgamma(x,shape,rate)/td # relative foliage density
  tdf<-(pai/hgt)*rfd*10 # total foliage density
  paia<-pgamma(x,shape,rate)*(pai/td) # pai above
  return(list(leafden=tdf,pai_a=paia))
}
#' Estimate temperature and humidity at specified height above ground
#'
#' @description The function `temphumE` runs the above ground component of
#' the microclimate model.
#'
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above ground at which model outputs are needed (m). Must be positive
#' @param pai_a an optional array of plant area index values above `reqhgt` (see details)
#' @param folden an optional array of foliage densities at height `reqhgt` (m^3/m^3)
#' @param surfwet an optional single numeric value of array of values specifying the proportion
#' of the canopy surface that should be treated as wet surface (see details)
#' @param xyf optional input for called function [wind()]
#' @param zf optional input for called function [wind()]
#' @param soilinit initial soil moisture fractions in surface and subsurface layer (see [soilmpredict()])
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
#' @param maxhgt an optional height (m) for which wind speed is needed. Determined
#' from height of tallest vegetation or as 2 m, whichever is greater, if not supplied.
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
#' @return an object of class microout with the following components:
#' @return `Tz` Array of air temperatures at height `reqhgt` (deg C)
#' @return `tleaf` Array of leaf temperatures at height `reqhgt` (deg C). if `reqhgt` greater than canopy height, the returned temperatures are those that would be attained by an unshaded leaf with average reflectance
#' @return `T0` Array of ground surface temperatures (deg C)
#' @return `relhum` Array of relative humidities at height `reqhgt` (percentage)
#' @return `windspeed` Array of wind speeds at height `reqhgt` (m/s)
#' @return `Rdirdown` Array of downward direct shortwave radiation received on horizontal surface at height `reqhgt` (W/m^2 - per unit one-sided area)
#' @return `Rdifdown` Array of downward diffuse shortwave radiation received on horizontal surface at height `reqhgt` (W/m^2 - per unit one-sided area)
#' @return `Rlwdown` Array of downward longwave radiation received on horizontal surface at height `reqhgt` (W/m^2 - per unit one-sided area)
#' @return `Rswup` Array of upward shortwave radiation received on the underside of a horizontal surface at height `reqhgt` (W/m^2 - per unit one-sided area). Upward shortwave radiation is assumed to be diffuse.
#' @return `Rlwup` Array of upward longwave radiation received on the underside of a horizontal surface at height `reqhgt` (W/m^2 - per unit one-sided area).
#' @seealso [below_hr()] and [below_day()] for running the below ground components of the microclimate model
#' @details `pai_a` is used to calculate the radiation intercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calculated from total plant area index by
#' assuming leaf density within the canopy is uniformly vertically distributed. If supplied
#' it must have the same dimensions as micro$pai. I.e. with the same x and y dims as the
#' the supplied dtm and values for each hour as the z dimension. The parameter `surfwet`
#' determines how much of the canopy should be treated as wet surface when calaculating
#' latent heat fluxes. However, except when extremely droughted, the matric potential of leaves
#' is such that `surfwet` ~ 1.
#'
#' @import terra zoo
#' @export
temphumE<-function(micro, reqhgt, pai_a = NA, folden = NA, xyf = 1, zf = NA, soilinit = c(NA, NA),
                   tfact = 1.5, surfwet = 1, slr = NA, apr = NA, hor = NA, wsa = NA,
                   maxhgt = NA, twi = NA) {
  # run soiltemp function if not run
  if (micro$progress<3) {
      micro<-soiltemp_hr(micro,reqhgt,pai_a,xyf,zf,soilinit,tfact,slr,apr,hor,wsa,maxhgt,twi)
  }
  # Calculate boundary layer conductivity
  pai<-micro$pai
  pai[pai<1]<-1
  lwi<-.rta(micro$leafd,dim(pai)[3])
  gmin<-.gfree(lwi,micro$Hest)*2*pai
  gHa<-.gturb(micro$uf,micro$d,micro$zm,micro$maxhgt,psi_h=micro$psi_h,gmin=gmin)
  # Calculate vapour conductivity
  climdata<-micro$climdata
  fsun<-.psunlit(micro$pai,micro$k,micro$kc,micro$clump)
  fsund<-.psunlitd(micro$pai,micro$clump)
  lais<-micro$pai*(micro$dp*fsund+(1-micro$dp)*fsun)
  gs<-.layercond(.vta(climdata$swrad,micro$dtm),micro$gsmax,100)
  gBs<-lais*gs
  # Calculate T0
  canabs<-micro$radCsw+micro$radClw
  # Calculate Latent heat flux
  TH<-PenMont(micro$tc,micro$pk,micro$ea,canabs,gHa,gBs,micro$G,micro$estl,
              micro$tdew,surfwet)
  # Set limits in TH
  TH$tcan<-.lim(TH$tcan,(micro$tc+40),up=TRUE)
  TH$tcan<-.lim(TH$tcan,80,up=TRUE)
  # Calculate temperature and vapour pressure above canopy, setting height to canopy top of reghgt < hgt
  z<-micro$vha
  sel<-which(micro$vha<reqhgt)
  z[sel]<-reqhgt
  Tzv<-.TVabove(TH,micro,z)
  micro$Tz<-Tzv$Tz
  ez<-Tzv$ez
  # Computes below canopy temperatures (Tz set to above canopy if reqhgt>hgt)
  Tzp<-.LangrangianSimT(reqhgt,micro,TH)
  micro$Tz<-Tzp$To
  # Compute leaf temperature
  micro<-.leaftemp(micro,gs,reqhgt,TH$tcan)
  # Compute relative humidity
  rh<-.LangrangianSimV(reqhgt,micro,ez,surfwet)
  # Return values needed
  micro<-list(Tz=micro$Tz,tleaf=micro$tleaf,T0=micro$T0,soilm=micro$theta,relhum=rh,windspeed=micro$uz,
              Rdirdown=micro$Rbdown,Rdifdown=micro$Rddown,Rlwdown=micro$Rldown,
              Rswup=micro$Rdup,Rlwup=micro$Rlup)
  class(micro)<-"microout"
  return(micro)
}
#' Estimate temperature below ground (hourly)
#'
#' @description The function `below_hr` runs the below ground component of
#' the microclimate model at hourly time increments.
#'
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above ground at which model outputs are needed (m). Must be negative.
#' @param xyf optional input for called function [wind()]
#' @param zf optional input for called function [wind()]
#' @param soilinit initial soil moisture fractions in surface and subsurface layer (see [soilmpredict()])
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
#' @param maxhgt an optional height (m) for which wind speed is needed. Determined
#' from height of tallest vegetation or as 2 m, whichever is greater, if not supplied.
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
#' @return `Tz` Array of air temperatures at height `reqhgt` below ground (deg C)
#' @seealso [temphumE()] for running the above ground component of the microclimate model
#'
#' @import terra zoo abind
#' @importFrom stats filter
#' @export
below_hr<-function(micro, reqhgt, pai_a = NA, xyf = 1, zf = NA, soilinit = c(NA, NA), tfact =1.5,
                   slr = NA, apr = NA, hor = NA, wsa = NA, maxhgt = NA, twi = NA) {
  xx <- function(x,y) as.numeric(stats::filter(c(x,x),rep(1/y,y), sides = 2))
  if (is.null(micro$T0[1])) {
    micro<-soiltemp_hr(micro,reqhgt,pai_a,xyf,zf,soilinit,tfact,slr,apr,hor,wsa,maxhgt,twi)
  }
  hiy<-length(micro$tme)
  Tav<-apply(micro$T0,c(1,2),mean)
  Tan<-micro$T0-.rta(rast(Tav),hiy)
  hr_sm<-(2*pi*abs(reqhgt)^2)/(7200*mean(micro$ka,na.rm=TRUE)) # calcaulate effective periodicity of soil fluctuations
  hr_ps<-round((24*abs(reqhgt))/(2*pi*mean(micro$DD,na.rm=TRUE)),0)
  if (hr_sm < hiy) {
    Tz<-aperm(apply(Tan,c(1,2),xx,hr_sm),c(2,3,1))
    Tz<-pmax(Tz[,,1:hiy],Tz[,,(hiy+1):dim(Tz)[3]],na.rm=T)
    if (hr_ps > 0) Tz<-abind(Tz[,,(dim(Tz)[3]-(hr_ps-1)):dim(Tz)[3]],Tz)
    Tz<-Tz[,,1:hiy]+.rta(rast(Tav),hiy)
  } else Tz<-.rta(rast(Tav),hiy)
  # Return values needed
  return(Tz)
}
#' Estimate temperature below ground (daily)
#'
#' @description The function `below_dy` runs the below ground component of
#' the microclimate model in daily time increments.
#'
#' @param microd object of class microin as returned by [modelin()]
#' @param reqhgt height above ground at which model outputs are needed (m). Must be negative.
#' @param expand optional logical indicating whethet to expand daily values to hourly
#' @param xyf optional input for called function [wind()]
#' @param zf optional input for called function [wind()]
#' @param soilinit initial soil moisture fractions in surface and subsurface layer (see [soilmpredict()])
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
#' @param maxhgt an optional height (m) for which wind speed is needed. Determined
#' from height of tallest vegetation or as 2 m, whichever is greater, if not supplied.
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from `dtm` of not supplied, but outer cells will be NA.
#' @return `Tz` Array of air temperatures at height `reqhgt` below ground (deg C)
#' @seealso [below_hr()] for running the below ground component of the microclimate model in hourly timesteps
#'
#' @import terra zoo abind
#' @importFrom stats filter
#' @export
below_dy<-function(microd, reqhgt, expand = TRUE, pai_a = NA, xyf = 1, zf = NA, soilinit = c(NA, NA), tfact = 1.5,
                   slr = NA, apr = NA, hor = NA, wsa = NA, maxhgt = NA, twi = NA) {
  xx <- function(x,y) as.numeric(stats::filter(c(x,x),rep(1/y,y), sides = 2))
  micro_mn<-microd$micro_mn
  if (micro_mn$progress<5) {
    microd<-soiltemp_dy(microd,reqhgt,pai_a,xyf,zf,soilinit,tfact,slr,apr,hor,wsa,maxhgt,twi)
    micro_mn<-microd$micro_mn
  }
  micro_mx<-microd$micro_mx
  # (1) Calculate average temperature grid
  T0d<-(micro_mn$T0+micro_mx$T0)/2
  Tav<-apply(T0d,c(1,2),mean)
  #(2) Calculate periodicity and phase shift
  l<-dim(T0d)[3]
  Tan<-T0d-.rta(rast(Tav),l)
  day_sm<-(2*pi*abs(reqhgt)^2)/(7200*24*mean(micro_mn$ka,na.rm=TRUE)) # calcaulate effective periodicity of soil fluctuations
  day_ps<-round(abs(reqhgt)/(2*pi*mean(micro_mn$DD,na.rm=TRUE)),0)
  # (3) Apply periodicity and phase shift
  if (day_sm < l) {
    if (day_sm > 1) {
      Tzd<-aperm(apply(Tan,c(1,2),xx,day_sm),c(2,3,1))
      Tzd<-pmax(Tzd[,,1:l],Tzd[,,(l+1):dim(Tzd)[3]],na.rm=T)
      x1<-(dim(Tzd)[3]-(day_ps-1))
      x2<-dim(Tzd)[3]
      if (day_ps > 0) Tzd<-abind(Tzd[,,(dim(Tzd)[3]-(day_ps-1)):dim(Tzd)[3]],Tzd)
      Tzd<-Tzd[,,1:l]+.rta(rast(Tav),l)
    } else {
      if (day_ps > 0) {
        Tzd<-abind(T0d[,,(dim(T0d)[3]-(day_ps-1)):dim(T0d)[3]],T0d)
        Tzd<-Tzd[,,1:l]
      } else Tzd<-T0d
    }
  } else Tzd<-.rta(rast(Tav),l)
  if (expand) {
    # (4) Calculate wgts
    d<-dim(Tzd)[3]
    r<-rast(Tzd[,,1])
    wg2<-rep(c(0:23)/23,d+1)
    wg1<-1-wg2
    # (5) Apply weights
    Tzd1<-.ehr(abind(Tzd[,,1],Tzd))
    Tzd2<-.ehr(abind(Tzd,Tzd[,,d]))
    Tzh<-.vta(wg1,r)*Tzd1+.vta(wg2,r)*Tzd2
    Tz<-Tzh[,,13:((d+1)*24-12)]
  } else {
    Tz<-Tzd
  }
  return(Tz)
}
#' Create an object of type `soilcharac`
#' @description The function `soilcfromtype` creates an object of class `soilcharac`
#' from a SpatRaster object of soil types and specified ground reflectance
#'
#' @param a SpatRaster object of soiltypes numbered numerically as integers as in `soilparameters`
#' @param groundr a single numeric value, matrix or SpatRaster object of soil reflectance (in range 0-1)
#' @return an object of class `soilcharac`, a list of two SpatRaster object of `soiltype` and
#' `groundr` converted to a SpatRaster object.
#' @import terra
#' @export
soilcfromtype <- function(soiltype, groundr = 0.15) {
  if (class(soiltype)[1] == "PackedSpatRaster") soiltype=rast(soiltype)
  if (class(groundr)[1] == "PackedSpatRaster") groundr=rast(groundr)
  r<-soiltype
  if (class(groundr)[1] != "SpatRaster") {
    groundr<-.is(soiltype)*0+groundr
    groundr<-.rast(groundr,r)
  }
  soilc<-list(soiltype=soiltype,groundr=groundr)
  class(soilc) <- "soilcharac"
  return(soilc)
}
