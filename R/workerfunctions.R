# =========================================================================== #
# ******************* True worker functions ********************************* #
# =========================================================================== #
#' Check if input is a SpatRaster or PackedSpatRaster and convert to matrix if it is
#' @import terra
.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r<-rast(r)
  if (class(r)[1] != "matrix") {
    if (dim(r)[3] > 1) {
      y<-as.array(r)
    } else y<-as.matrix(r,wide=TRUE)
  } else y<-r
  y
}
#' Convert matrix to array
.rta <- function(r,n) {
  m<-.is(r)
  a<-array(rep(m,n),dim=c(dim(r)[1:2],n))
  a
}
#' Convert vector to array
.vta <- function(v,r) {
  m<-.is(r)
  va<-rep(v,each=dim(m)[1]*dim(m)[2])
  a<-array(va,dim=c(dim(m),length(v)))
  a
}
#' Create SpatRaster object using a template
#' @import terra
.rast <- function(m,tem) {
  r<-rast(m)
  ext(r)<-ext(tem)
  crs(r)<-crs(tem)
  r
}
#' Calculates latitude and longitude from SpatRaster object
#' @import terra
.latlongfromraster<-function (r) {
  e <- ext(r)
  xy <- data.frame(x = (e$xmin + e$xmax)/2, y = (e$ymin + e$ymax)/2)
  xy <- sf::st_as_sf(xy, coords = c("x", "y"),
                     crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[2], long = sf::st_coordinates(ll)[1])
  return(ll)
}
#' Get mean, min, max etc from raster or packaged raster
#' @import terra
.mfr <- function(r,fun=mean) {
  if (class(r)[1] == "PackedSpatRaster") r<-rast(r)
  m<-fun(as.vector(r),na.rm=T)
  m
}
#' Calculates mode
.getmode <- function(v) {
  v<-v[is.na(v)==FALSE]
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#' Calculate moving average / rolling mean
.ma <- function(x, n = 5) {
  filter(x, rep(1 / n, n), sides = 1, circular = TRUE)
}
#' Cap values
.lim<-function(x,l,up=FALSE) {
  if (length(l)==1) l<-x*0+l
  if (up) {
    s<-which(x>l)
    x[s]<-l[s]
  } else {
    s<-which(x<l)
    x[s]<-l[s]
  }
  x
}
#' convert degrees to radians
.ar <- function(d) {
  d*(pi/180)
}
#' expand daily array to hourly array
.ehr<-function(a) {
  n<-dim(a)[1]*dim(a)[2]
  o1<-rep(c(1:n),24*dim(a)[3])
  o2<-rep(c(1:dim(a)[3]),each=24*n)-1
  o2<-o2*max(o1,na.rm=T)
  o<-o1+o2
  ah<-rep(a,24)
  ah<-ah[o]
  ah<-array(ah,dim=c(dim(a)[1:2],dim(a)[3]*24))
  ah
}
#' quicker implementation of apply used by .hourtoday
.applynotna<-function(a,fun) {
  m<-matrix(a,ncol=dim(a)[1]*dim(a)[2],byrow=T)
  sel<-which(is.na(m[1,])==F)
  r<-apply(m[,sel],2,fun)
  n<-dim(r)[1]
  ao<-array(NA,dim=c(dim(a)[1:2],n))
  sel<-which(is.na(a[,,1:n])==F)
  ao[sel]<-aperm(r,c(2,1))
  ao
}
#' Extract daily min, mean or max from array of hourly data
.hourtoday<-function(a,fun=mean) {
  .htd<-function(x) {
    y<-matrix(x,ncol=24,byrow=T)
    apply(y,1,fun,na.rm=T)
  }
  d<-.applynotna(a,.htd)
  d
}
# =========================================================================== #
# ********** Functions associated with modelin and steps before that ******** #
# =========================================================================== #
#' Calculates the astronomical Julian day
.jday <- function(tme) {
  yr<-tme$year+1900
  mth<-tme$mon+1
  dd<-tme$mday+(tme$hour+(tme$min+tme$sec/60)/60)/24
  madj<-mth+(mth<3)*12
  yadj<-yr+(mth<3)*-1
  jd<-trunc(365.25*(yadj+4716))+trunc(30.6001*(madj+1))+dd-1524.5
  b<-(2-trunc(yadj/100)+trunc(trunc(yadj/100)/4))
  jd<-jd+(jd>2299160)*b
  jd
}
#' Calculates solar time
.soltime <- function(localtime, long, jd, merid = 0, dst = 0) {
  m<-6.24004077+0.01720197*(jd-2451545)
  eot<- -7.659*sin(m)+9.863*sin(2*m+3.5932)
  st<-localtime+(4*(long-merid)+eot)/60-dst
  st
}
#' Calculates the solar altitude
.solalt <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st<-.soltime(localtime,long,jd,merid,dst)
  tt<-0.261799*(st-12)
  d<-(pi*23.5/180)*cos(2*pi*((jd-159.5)/365.25))
  sh<-sin(d)*sin(lat*pi/180)+cos(d)*cos(lat*pi/180)*cos(tt)
  sa<-(180*atan(sh/sqrt(1-sh^2)))/pi
  sa
}
#' Calculate clear sky radiation
.clearskyrad <- function(tme, lat, long, tc, rh, pk) {
  jd<-.jday(tme)
  lt <- tme$hour+tme$min/60+tme$sec/3600
  sa<-.solalt(lt,lat,long,jd)*pi/180
  m<-35*sin(sa)*((1224*sin(sa)^2+1)^(-0.5))
  TrTpg<-1.021-0.084*(m*0.00949*pk+0.051)^0.5
  xx<-log(rh/100)+((17.27*tc)/(237.3+tc))
  Td<-(237.3*xx)/(17.27-xx)
  u<-exp(0.1133-log(3.78)+0.0393*Td)
  Tw<-1-0.077*(u*m)^0.3
  Ta<-0.935*m
  od<-TrTpg*Tw*Ta
  Ic<-1352.778*sin(sa)*TrTpg*Tw*Ta
  Ic[is.na(Ic)]<-0
  Ic
}
#' unpack model inputs
.unpack<-function(dtm,vegp,soilc) {
  if (class(dtm)[1] == "PackedSpatRaster") dtm<-rast(dtm)
  if (class(vegp$hgt)[1] == "PackedSpatRaster") vegp$hgt<-rast(vegp$hgt)
  if (class(vegp$x)[1] == "PackedSpatRaster") vegp$x<-rast(vegp$x)
  if (class(vegp$gsmax)[1] == "PackedSpatRaster") vegp$gsmax<-rast(vegp$gsmax)
  if (class(vegp$leafr)[1] == "PackedSpatRaster") vegp$leafr<-rast(vegp$leafr)
  if (class(vegp$leaft)[1] == "PackedSpatRaster") vegp$leaft<-rast(vegp$leaft)
  if (class(vegp$leafd)[1] == "PackedSpatRaster") vegp$leafd<-rast(vegp$leafd)
  if (class(soilc$soiltype)[1] == "PackedSpatRaster") soilc$soiltype<-rast(soilc$soiltype)
  if (class(soilc$groundr)[1] == "PackedSpatRaster") soilc$groundr<-rast(soilc$groundr)
  crs(vegp$hgt)<-crs(dtm)
  crs(vegp$x)<-crs(dtm)
  crs(vegp$gsmax)<-crs(dtm)
  crs(vegp$leafr)<-crs(dtm)
  crs(vegp$leafd)<-crs(dtm)
  crs(soilc$soiltype)<-crs(dtm)
  crs(soilc$groundr)<-crs(dtm)
  return(list(dtm=dtm,vegp=vegp,soilc=soilc))
}
#' Calculate point surface temperature without diabatic coefficients and G set to 0
.point0<-function(micropoint) {
  d<-with(micropoint$vegp_p,.zeroplanedis(h,pai))
  zm<-with(micropoint$vegp_p,.roughlength(h,pai,d))
  uf<-(0.4*micropoint$weather$windspeed)/log((micropoint$zref-d)/zm)
  gHa<-(0.4*43*uf)/log((micropoint$zref-d)/(0.2*zm))
  es<-.satvap(micropoint$weather$temp)
  ea<-es*(micropoint$weather$relhum/100)
  T0p<-PenMont(micropoint$weather$temp,micropoint$weather$pres,ea,
               micropoint$microp$RabsG,gHa,gHa,es,NA)
  return(T0p$T0)
}
#' Calculate zero plane displacement
.zeroplanedis<-function(h,pai) {
  pai[pai<0.001]<-0.001
  d<-(1-(1-exp(-sqrt(7.5*pai)))/sqrt(7.5*pai))*h
  d
}
#' Calculate roughness length
.roughlength<-function(h,pai,d=NA,psi_h=0) {
  if (class(d)[1]=="logical") d<-.zeroplanedis(h,pai)
  Be<-sqrt(0.003+(0.2*pai)/2)
  zm<-(h-d)*exp(-0.4/Be)*exp(psi_h)
  zm[zm<0.0005]<-0.0005
  zm
}
#' Calculate saturated vapour pressure
.satvap <- function(tc) {
  es<-0.61078*exp(17.27*tc/(tc+237.3))
  ei<-0.61078*exp(21.875*tc/(tc+265.5))
  s<-which(tc<0)
  es[s]<-ei[s]
  es
}
#' Calculate dewpoint
.dewpoint <- function(ea, tc) {
  e0<-611.2/1000
  L<-(2.501*10^6)-(2340*tc)
  it<-1/273.15-(461.5/L)*log(ea/e0)
  Tdew<-1/it-273.15
  e0<-610.78/1000
  L<-2.834*10^6
  it<-1/273.15-(461.5/L)*log(ea/e0)
  Tfrost<-1/it-273.15
  sel<-which(Tdew<0)
  Tdew[sel]<-Tfrost[sel]
  Tdew
}
#' Calculates slope of the saturated vapour pressure curve
.delta <- function(tc) {
  es1<-.satvap(tc-0.5)
  es2<-.satvap(tc+0.5)
  delta<-es2-es1
  delta
}
#' Smooth SpatRaster
#' @import terra
.smr <- function(r,f,cors=NA) {
  r2<-r
  r2[is.na(r2)]<-0
  r3<-aggregate(r2,fact=f)
  if (class(cors)!="logical") {
    crs(r3)<-cors
    crs(r2)<-cors
    crs(r)<-cors
  }
  r3<-resample(r3,r2)
  r3<-mask(r3,r)
  r3
}
#' Smooth array
#' @import terra
.sma <- function(a,xyf,zf,cors=NA) {
  a3<-a
  a[is.na(a)]<-0
  ns<-floor(c(1:(dim(a)[3]/zf)))*zf
  ns<-unique(c(1,ns,dim(a)[3]))
  a2<-array(NA,dim=dim(a))
  for (i in 1:length(ns)) {
    a2[,,ns[i]]<-.is(.smr(rast(a[,,ns[i]]),xyf,cors))
  }
  a2<-apply(a2,c(1,2),na.approx)
  a2<-aperm(a2,c(2,3,1))
  a2[is.na(a3)]<-NA
  a2
}
#' Sort out smoothing of zero-plane displacement height and roughness length
.sortrough<-function(vegp,xyf) {
  # Calculate unsmoothed zero plane displacement height and roughness length
  hp<-.rta(vegp$hgt,dim(vegp$pai)[3])
  d<-.zeroplanedis(hp,vegp$pai)
  zm<-.roughlength(hp,vegp$pai,d)
  if (is.na(xyf) == FALSE) {
    # Smooth if xyf > 1
    if (xyf > 1) {
      d<-.sma(d,xyf,1)
      zm<-.sma(zm,xyf,1)
    }
  } else {
    # Average across study area if NA
    dp<-apply(d,3,mean,na.rm=T)
    zp<-apply(zm,3,mean,na.rm=T)
    d<-.vta(dp,vegp$hgt)
    zm<-.vta(dp,vegp$hgt)
  }
  zm[zm<0.001]<-0.001
  d[d<0.0065]<-0.0065
  return(list(d=d,zm=zm))
}
#' Unpacks pai values
.unpackpai <- function(pai,n) {
  x<-dim(pai)[3]
  i<-floor((x/n)*c(0:(n-1))+1)
  pai<-pai[,,i]
  pai
}
#' initialises soil parameters
.soilinit <- function(soilc) {
  soiltype<-soilc$soiltype
  u<-unique(as.vector(.is(soiltype)))
  u<-u[is.na(u)==F]
  rho<-array(NA,dim=dim(soiltype)[1:2])
  Vm<-rho; Vq<-rho; Mc<-rho;
  Smax<-rho; Smin<-rho
  psi_e<-rho; soilb<-rho
  for (i in 1:length(u)) {
    sel<-which(soilparameters$Number==u[i])
    sop<-soilparameters[sel,]
    sel<-which(.is(soiltype)==u[i])
    rho[sel]<-sop$rho
    Vm[sel]<-sop$Vm
    Vq[sel]<-sop$Vq
    Mc[sel]<-sop$Mc
    psi_e[sel]<-sop$psi_e
    soilb[sel]<-sop$b
    Smax[sel]<-sop$Smax
    Smin[sel]<-sop$Smin
  }
  return(list(rho=rho,Vm=Vm,Vq=Vq,Mc=Mc,psi_e=psi_e,soilb=soilb,Smax=Smax,Smin=Smin))
}
# =========================================================================== #
# ********** Functions associated with running the model ******************** #
# =========================================================================== #
# *** soilmdistribute()
#' Calculates flow direction
.flowdir <- function(md) {
  fd <- md * 0
  md2 <- array(NA, dim = c(dim(md)[1] + 2, dim(md)[2] + 2))
  md2[2:(dim(md)[1] + 1), 2:(dim(md)[2] + 1)] <- md
  v <- c(1:length(md))
  v <- v[is.na(md) == F]
  x <- arrayInd(v, dim(md))[, 1]
  y <- arrayInd(v, dim(md))[, 2]
  for (i in 1:length(x)) {
    md9 <- md2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
    fd[x[i], y[i]] <- round(mean(which(md9 == min(md9, na.rm = T)),na.rm=T), 0)
  }
  fd
}
#' Caclulates flow accumulation
.flowacc <- function (dtm) {
  dm<-.is(dtm)
  dm[is.na(dm)]<-0
  nx<-dim(dm)[1]
  ny<-dim(dm)[2]
  # stick buffer around
  dm3<-array(0,dim=c(nx+4,ny+4))
  dm3[3:(nx+2),3:(ny+2)]<-dm
  fd<-.flowdir(dm3)
  fa<-fd*0+1
  o<-order(dm3,decreasing=T,na.last=NA)
  for (i in 1:length(o)) {
    x <- arrayInd(o[i], dim(dm3))[1]
    y <- arrayInd(o[i], dim(dm3))[2]
    f <- fd[x, y]
    x2 <- x + (f - 1)%%3 - 1
    y2 <- y + (f - 1)%/%3 - 1
    if (x2 > 0 & x2 < dim(dm3)[1] & y2 > 0 & y2 < dim(dm3)[2])
      fa[x2, y2] <- fa[x, y] + 1
  }
  fa<-fa[3:(nx+2),3:(ny+2)]
  fa
}
#' Calculates topographic wetness index
.topidx <- function(dtm) {
  minslope<-atan(0.02/mean(res(dtm)))
  slope<-terrain(dtm,unit="radians")
  B<-.is(slope)
  B[B<minslope]<-minslope
  B[is.na(B)]<-median(B,na.rm=T)
  a<-.flowacc(dtm)
  a<-a*res(dtm)[1]* res(dtm)[2]
  tpx<- a/tan(B)
  r<-.rast(tpx,dtm)
  r
}
#' Calculates the solar azimuth
.solazi <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st<-.soltime(localtime,long,jd,merid,dst)
  tt<-0.261799*(st-12)
  d<-(pi*23.5/180)*cos(2*pi*((jd-159.5)/365.25))
  sh<-sin(d)*sin(lat*pi/180)+cos(d)*cos(lat*pi/180)*cos(tt)
  hh<-0.5*pi-acos(sh)
  sazi<-cos(d)*sin(tt)/cos(hh)
  cazi<-(sin(.ar(lat))*cos(d)*cos(tt)-cos(.ar(lat))*sin(d))/
    sqrt((cos(d)*sin(tt))^2+(sin(.ar(lat))*cos(d)*cos(tt)-cos(.ar(lat))*sin(d))^2)
  sqt <- 1-sazi^2
  sqt[sqt<0]<-0
  solz<-180+(180*atan(sazi/sqrt(sqt)))/pi
  solz[cazi<0 & sazi<0]<-180-solz[cazi<0 & sazi<0]
  solz[cazi<0 & sazi>=0]<-540-solz[cazi<0 & sazi>=0]
  solz
}
#' Calculates the solar coefficient
#' @import terra
.solarindex<- function(dtm,alt,azi,slr=NA,apr=NA) {
  if (class(slr)[1]=="logical") slr<-terrain(dtm,v="slope",unit="radians")
  if (class(apr)[1]=="logical") apr<-terrain(dtm,v="aspect",unit="radians")
  slr[is.na(slr)]<-0
  apr[is.na(apr)]<-pi/2
  sl<-.rta(slr,length(azi))
  ap<-.rta(apr,length(azi))
  zen<-pi/2-alt
  sazi<-.ar(.vta(azi,dtm))
  i<-cos(zen)*cos(sl)+sin(zen)*sin(sl)*cos(sazi-ap)
  i[i<0]<-0
  i
}
#' Calculates tangent of horizon angle
.horizon <- function(dtm, azimuth) {
  reso<-res(dtm)[1]
  dtm<-.is(dtm)
  dtm[is.na(dtm)]<-0
  dtm<-dtm/reso
  azi<-.ar(azimuth)
  horizon<-array(0,dim(dtm))
  dtm3<-array(0,dim(dtm)+200)
  x<-dim(dtm)[1]
  y<-dim(dtm)[2]
  dtm3[101:(x+100),101:(y+100)]<-dtm
  for (step in 1:10) {
    horizon[1:x,1:y]<-pmax(horizon[1:x,1:y],(dtm3[(101-cos(azi)*step^2):(x+100-cos(azi)*step^2),
                                                  (101+sin(azi)*step^2):(y+100+sin(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(step^2),na.rm=T)
  }
  horizon
}
#' Calculates radiation extinction coefficient for canopy
.cank <- function(x,sa,si) {
  # Raw k
  zen<-pi/2-sa
  k<-sqrt((x^2+(tan(zen)^2)))/(x+1.774*(x+1.182)^(-0.733))
  k0<-sqrt(x^2)/(x+1.774*(x+1.182)^(-0.733))
  # k dash
  kd<-k*cos(zen)/si
  sel<-which(si==0)
  kd[sel]<-1
  return(list(k=k,kd=kd,k0=k0))
}
# *** wind()
#' Calculates wind shelter coefficient in specified direction
.windcoef <- function(dsm, direction, hgt = 1, res = 1) {
  reso<-res(dsm)[1]
  dsm<-.is(dsm)
  dsm[is.na(dsm)]<-0
  dtm<-dsm/reso
  hgt<-hgt/reso
  azi<-direction*(pi/180)
  horizon <- array(0, dim(dtm))
  dtm3 <- array(0, dim(dtm) + 200)
  x <- dim(dtm)[1]
  y <- dim(dtm)[2]
  dtm3[101:(x + 100), 101:(y + 100)] <- dtm
  for (step in 1:10) {
    horizon[1:x,1:y]<-pmax(horizon[1:x,1:y],(dtm3[(101-cos(azi)*step^2):(x+100-cos(azi)*step^2),
                                                  (101+sin(azi)*step^2):(y+100+sin(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(step^2),na.rm=T)
    horizon[1:x,1:y]<-ifelse(horizon[1:x,1:y]<(hgt/step^2),0,horizon[1:x,1:y])
  }
  index<-1-atan(0.17*100*horizon)/1.65
  index
}
#' Calculates array of wind shelter coefficients in wdct directions
.windsheltera<-function(dtm, whgt, s) {
  if (is.na(s)) {
    s<-min(dim(dtm)[1:2])
    s<-ifelse(s>10,10,s)
  }
  dct2<-16
  drs<-c(0:(dct2-1))*360/dct2
  a<-array(NA,dim=c(dim(dtm)[1:2],dct2))
  for (i in 1:dct2) {
    r<-.rast(.windcoef(dtm,drs[i],whgt,res(dtm)[1]),dtm)
    r<-resample(aggregate(r,fact=s,fun="mean"),dtm)
    a[,,i]<-as.matrix(r,wide=T)
  }
  a2<-array(NA,dim=c(dim(dtm)[1:2],8))
  for (i in 1:8) {
    if (i == 1) {
      m<-0.5*a[,,1]+0.25*a[,,2]+0.25*a[,,dct2]
    } else m<-0.5*a[,,(i*2-1)]+0.25*a[,,i*2]+0.25*a[,,(i*2)-2]
    a2[,,i]<-m
  }
  a2
}
#' Calculate shelter coefficient for all wind directions
#' @import terra
.windshelter <- function(wdir,dtm,whgt,s,wsa = NA) {
  # Create array of wind shelter coefficients
  if (class(wsa)[1] == "logical") wsa<-.windsheltera(dtm,whgt,s)
  # Apply array of wind shelter coefficients
  i<-360/dim(wsa)[3]
  # if wind direction is an array
  if (length(dim(wdir)) == 3) {
    wa<-wdir*0
    sel<-(round(wdir/i,0))+1
    sel[sel==(dim(wsa)[3])+1]<-1
    sel[sel==0]<-dim(wsa)[3]
    for (ii in 1:dim(wa)[1]) {
      for (jj in 1:dim(wa)[2]) {
        wa[ii,jj,]<-wsa[ii,jj,sel[ii,jj,]]
      }
    }
  # if wind direction is a vector
  } else {
    sel<-(round(wdir/i,0))+1
    sel[sel==(dim(wsa)[3])+1]<-1
    sel[sel==0]<-dim(wsa)[3]
    wa<-wsa[,,sel]
  }
  wa
}
#' Calculates turbulent molar conductivity above canopy
.gturb<-function(uf,d,zm,z,psi_h=0,gmin) {
  zh<-0.2*zm
  ln<-log((z-d)/zh)
  g<-(0.4*43*uf)/(ln+psi_h)
  g<-.lim(g,gmin)
  g<-.lim(g,0.0001)
  g
}
# *** soiltemp_hr()
#' Calculate soil relative humidity
.soilrh <- function(theta, b, Psie, Smax, tc) {
  matric <- -Psie*(theta/Smax)^-b
  hr<-exp((0.018*matric)/(8.31*(tc+273.15)))
  hr[hr>1]<-1
  hr
}
#' Calculate soil conductivity
.soilcond<-function(rho,Vm,Vq,Mc,soilm) {
  # Find soil diffusivity
  cs<-(2400*rho/2.64+4180*soilm) # specific heat of soil in J/kg/K
  ph<-(rho*(1-soilm)+soilm)*1000   # bulk density in kg/m3
  frs<-Vm+Vq
  c1<-(0.57+1.73*Vq+0.93*Vm)/(1-0.74*Vq-0.49*Vm)-2.8*frs*(1-frs)
  c2<-1.06*rho*soilm
  c3<-1+2.6*Mc^-0.5
  c4<-0.03+0.7*frs^2
  k<-c1+c2*soilm-(c1-c4)*exp(-(c3*soilm)^4) # Thermal conductivity   W/m/K
  ka<-k/(cs*ph)
  omdy<-(2*pi)/(24*3600)
  DD<-sqrt(2*ka/omdy)
  return(list(k=k,DD=DD))
}
# *** temphumE()
#' Calculate effective ground  surface wetness from vapour pressure difference
.predwet<-function(eT) {
  eT[eT<0.001]<-0.001
  let<-log(eT)
  plf<-0.8753-1.7126*let
  soilwet<-1/(1+exp(-plf))
  soilwet
}
#' Calculates bulk surface stomatal conductivity
.layercond <- function(Rsw, gsmax, q50 = 300) {
  rpar<-Rsw*4.6
  gs<-(gsmax * rpar) / (rpar + q50)
  gs
}
#' Calculates temperature and vapour above canopy
.TVabove<-function(TH,micro,z,surfwet) {
  # Log ratio
  hgt<-micro$vha
  hgt[hgt<0.001]<-0.001
  micro$d<-.zeroplanedis(hgt,micro$pai)
  micro$zm<-.roughlength(hgt,micro$pai,micro$d)
  zh<-0.2*micro$zm
  lnr<-suppressWarnings(log((z-micro$d)/zh)/log((micro$zref-micro$d)/zh))
  # Temperature
  micro$Tz<-micro$tc+(TH$T0-micro$tc)*(1-lnr)
  # Vapour pressure
  micro$ez<-micro$ea+(.satvap(TH$T0)*surfwet-micro$ea)*(1-lnr)
  return(micro)
}
#' Calculate leaf energy balance
.leafbalance<-function(micro,TH,surfwet) {
  # longwave radiation absorbed by leaf
  lwsky<-micro$skyem*0.97*5.67*10^-8*(micro$tc+273.15)^4
  lwcan<-0.97*5.67*10^-8*(TH$T0+273.15)^4
  lwgro<-0.97*5.67*10^-8*(micro$T0+273.15)^4
  pai_g<-with(micro,pai-pai_a)
  lwup<-exp(-pai_g)*lwgro+(1-exp(-pai_g))*lwcan
  lwdn<-exp(-micro$pai_a)*lwgro+(1-exp(-micro$pai_a))*lwcan
  lwabs<-0.97*0.5*(lwup+lwdn)
  # Total radiation absorbed by leaf
  leafabs<-micro$radLsw+lwabs
  # Conductances
  gh<-with(micro,0.135*sqrt(uz/.rta(leafd,dim(uz)[3]))*1.4)
  rsw<-with(micro,Rddown+sin(alt)*Rbdown)
  gv<-.layercond(rsw,.rta(micro$gsmax,dim(micro$pai)[3]), q50 = 100)
  gV<-1/(1/gh+1/gv)
  Thl<-with(micro,PenMont(tc,pk,ea,leafabs,gh,gV,estl,tdew,surfwet))
  Thl$lwup<-lwup
  Thl$lwdn<-lwdn
  return(Thl)
}
#' Calculate inverse thermal diffusivity below canopy
.rhcanopy<-function(uf,h,d,z,phih=1) {
  dF<-d/h
  a2<-0.4*(1-dF)/(phih*1.25^2)
  int<-(2*h*((48*atan((sqrt(5)*sin((pi*z)/h))/(cos((pi*z)/h)+1)))/5^(3/2)+
               (32*sin((pi*z)/h))/((cos((pi*z)/h)+1)*((25*sin((pi*z)/h)^2)/
                                                        (cos((pi*z)/h)+1)^2+5))))/pi
  inth<-4.293251*h
  s<-which(z==h)
  int[s]<-inth[s]
  mu<-((uf/(a2*h))*(1/uf^2))
  rHa<-int*mu
  rHa[rHa<0.001]<-0.001
  rHa
}
#' Langrangian approximation function
.Langapprox<-function(micro,reqhgt,Flux,Fluxz,SH,SG,mu) {
  #  Calculate far-field independent of ground
  a2<-0.4*(1-micro$d/micro$vha)/1.5625
  int<-(11/16)*a2*10^2*micro$uf
  farm<-with(micro,(Flux/int)*vha)
  farp<-farm*(1-reqhgt/micro$vha)+SH
  # Weight by ground temperature
  rHa<-suppressWarnings(with(micro,.rhcanopy(uf,vha,d,reqhgt)))
  rh<- suppressWarnings(with(micro,.rhcanopy(uf,vha,d,vha)))
  rh[is.na(rh)]<-mean(rh,na.rm=T)
  wgt1<-(1-rHa/rh)
  wgt1[wgt1<0]<-0
  wgt2<-abs(farp-SG)/abs(farm+SH-SG)
  wgt2[wgt2>1]<-1
  wgt2[is.na(wgt2)]<-1
  wgt<-wgt1*wgt2
  farg<-wgt*SG+(1-wgt)*farp
  # Near-field
  # Calculate leaf temperature and flux
  SN<-Fluxz*micro$leafden
  near<-(3.047519+0.128642*log(micro$pai))*SN
  near[is.na(near)]<-0
  SC<-near+farg
  return(SC)
}
# =========================================================================== #
# ********** Functions associated with daily data *************************** #
# =========================================================================== #
#' Function for selecting daily min and max from point model
.subsetpointtoday<-function(micropoint) {
  microdf<-with(micropoint,data.frame(tme=microp$tme,Tc=microp$Tc,Tg=microp$Tg,
                      H=microp$H,G=microp$G,psih=microp$psih,
                      psim=microp$psim,OL=microp$OL,uf=microp$uf,
                      RabsG=microp$RabsG))
  Tc<-matrix( microdf$Tc,ncol=24,byrow=T)
  ta<-c(1:dim(Tc)[1])*24-24
  imn<-apply(Tc,1,which.min)+ta
  imx<-apply(Tc,1,which.max)+ta
  pointmodel_mn<-list()
  pointmodel_mx<-list()
  pointmodel_mn$microp<-as.list(microdf[imn,])
  pointmodel_mx$microp<-as.list(microdf[imx,])
  pointmodel_mn$weather<-micropoint$weather[imn,]
  pointmodel_mx$weather<-micropoint$weather[imx,]
  pointmodel_mn$soilm<-micropoint$soilm[imn]
  pointmodel_mx$soilm<-micropoint$soilm[imx]
  if (class(micropoint$Tbz) != "logical") {
    pointmodel_mn$Tbz<-micropoint$Tbz[imn]
    pointmodel_mx$Tbz<-micropoint$Tbz[imx]
    pointmodel_mn$DD<-micropoint$DD[imn]
    pointmodel_mx$DD<-micropoint$DD[imx]
  } else{
    pointmodel_mn$Tbz<-NA
    pointmodel_mx$Tbz<-NA
    pointmodel_mn$DD<-NA
    pointmodel_mx$DD<-NA
  }
  pointmodel_mn$precip<-micropoint$precip
  pointmodel_mx$precip<-micropoint$precip
  pointmodel_mn$vegp_p<-micropoint$vegp_p
  pointmodel_mx$vegp_p<-micropoint$vegp_p
  pointmodel_mn$groundp_p<-micropoint$groundp_p
  pointmodel_mx$groundp_p<-micropoint$groundp_p
  pointmodel_mn$lat<-micropoint$lat
  pointmodel_mx$lat<-micropoint$lat
  pointmodel_mn$long<-micropoint$long
  pointmodel_mx$long<-micropoint$long
  pointmodel_mn$zref<-micropoint$zref
  pointmodel_mx$zref<-micropoint$zref
  pointmodel_mn$tstep="day"
  pointmodel_mx$tstep="day"
  return(list(micropoint_mn_mn=pointmodel_mn,micropoint_mx=pointmodel_mx))
}
#' Identify daily min and max
.vartodaily<-function(temp) {
  td<-matrix(temp,ncol=24,byrow=TRUE)
  smn<-apply(td,1,which.min)
  smx<-apply(td,1,which.max)
  ta<-(c(1:dim(td)[1])-1)*24
  return(list(smx=smx+ta,smn=smn+ta))
}
#' Expand daily data to hourly
.expandtohour<-function(amn,amx,smn,smx,v) {
  r<-rast(amn[,,1])
  dtm<-rep(v[smn],each=24)
  dtx<-rep(v[smx],each=24)
  dtr<-dtx-dtm
  hfr<-.vta((v-dtm)/dtr,r)
  adtr<-.ehr(amx-amn)
  amnh<-.ehr(amn)
  ao<-hfr*adtr+amnh
  ao
}
#' Expand daily climate variables to hourly (all)
.expandclim<-function(mout_mn,mout_mx,micropoint,reqhgt) {
  # Temperature max & min hours
  scd<-.vartodaily(micropoint$microp$Tc)
  sgd<-.vartodaily(micropoint$microp$Tg)
  if (reqhgt<0) sbd<-.vartodaily(micropoint$microp$Tbz)
  # Radiation max and min hours
  rdir<-with(micropoint$weather,swrad-difrad)
  rlwd<-with(micropoint$weather,0.97*skyem*5.67*10^-8*(temp+273.15)^4)
  rlwu<-0.97*5.67*10^-8*(micropoint$microp$Tg+273.15)^4
  srdir<-.vartodaily(rdir)
  srdif<-.vartodaily(micropoint$weather$difrad)
  srsw<-.vartodaily(micropoint$weather$swrad)
  srlwd<-.vartodaily(rlwd)
  srlwu<-.vartodaily(rlwu)
  # Wind
  swi<-.vartodaily(micropoint$weather$windspeed)
  # Expand values regardless of height
  T0<-.expandtohour(mout_mn$T0,mout_mx$T0,sgd$smn,sgd$smx,micropoint$microp$Tg)
  soilm<-.ehr(mout_mn$soilm)
  tleaf<-NA
  relhum<-NA
  windspeed<-NA
  Rdirdown<-NA
  Rdifdown<-NA
  Rswup<-NA
  Rlwup<-NA
  Rlwdown<-NA
  if (reqhgt < 0) {
    Tz<-.expandtohour(mout_mn$Tz,mout_mx$Tz,sbd$smn,sbd$smx,micropoint$microp$Tbz)
  } else {
    # Expand shortwave radiation
    Rdirdown<-.expandtohour(mout_mn$Rdirdown,mout_mx$Rdirdown,srdir$smn,srdir$smx,rdir)
    Rdifdown<-.expandtohour(mout_mn$Rdifdown,mout_mx$Rdifdown,srdif$smn,srdif$smx,micropoint$weather$difrad)
    Rswup<-.expandtohour(mout_mn$Rswup,mout_mx$Rswup,srsw$smn,srsw$smx,micropoint$weather$swrad)
    # Expand longwave radiation
    Rlwup<-.expandtohour(mout_mn$Rlwup,mout_mx$Rlwup,srlwu$smn,srlwu$smx,rlwu)
    Rlwdown<-.expandtohour(mout_mn$Rlwdown,mout_mx$Rlwdown,srlwd$smn,srlwd$smx,rlwd)
  }
  if (reqhgt == 0) Tz<-T0
  if (reqhgt > 0) {
    Tz<-.expandtohour(mout_mn$Tz,mout_mx$Tz,scd$smn,scd$smx,micropoint$microp$Tc)
    tleaf<-.expandtohour(mout_mn$tleaf,mout_mx$tleaf,scd$smn,scd$smx,micropoint$microp$Tc)
    # Expand relhum
    ea_mn<-(mout_mn$relhum/100)*.satvap(mout_mn$Tz)
    ea_mx<-(mout_mx$relhum/100)*.satvap(mout_mx$Tz)
    ea<-(ea_mn+ea_mx)/2
    ea<-.ehr(ea)
    relhum<-(ea/.satvap(Tz))*100
    relhum[relhum>100]<-100
    # Expand wind speed
    windspeed<-.expandtohour(mout_mn$windspeed,mout_mx$windspeed,swi$smn,swi$smx,micropoint$weather$windspeed)
  }
  out<-list(Tz=Tz,tleaf=tleaf,T0=T0,soilm=soilm,relhum=relhum,windspeed=windspeed,
            Rdirdown=Rdirdown,Rdifdown=Rdifdown,Rlwdown=Rlwdown,Rswup=Rswup,Rlwup=Rlwup)
  return(out)
}
# =========================================================================== #
# ********** Functions associated with arrays ******************************* #
# =========================================================================== #
#' Latitudes from SpatRaster object
.latsfromr <- function(r) {
  e <- ext(r)
  lts <- rep(seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}
#' Longitudes from SpatRaster object
.lonsfromr <- function(r) {
  e <- ext(r)
  lns <- rep(seq(e$xmin + res(r)[1] / 2, e$xmax - res(r)[1] / 2, length.out = dim(r)[2]), dim(r)[1])
  lns <- lns[order(lns)]
  lns <- array(lns, dim = dim(r)[1:2])
  lns
}
#' lats and longs from SpatRaster object, including reprojection
.latslonsfromr <- function(r) {
  lats<-.latsfromr(r)
  lons<-.lonsfromr(r)
  xy<-data.frame(x=as.vector(lons),y=as.vector(lats))
  xy <- sf::st_as_sf(xy, coords = c('x', 'y'), crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[,2],
                   long = sf::st_coordinates(ll)[,1])
  lons<-array(ll$long,dim=dim(lons))
  lats<-array(ll$lat,dim=dim(lats))
  return(list(lats=lats,lons=lons))
}
#' Lape rates
.lapserate <- function(tc, ea, pk) {
  rv<-0.622*ea/(pk-ea)
  lr<-9.8076*(1+(2501000*rv)/(287*(tc+273.15)))/
    (1003.5+(0.622*2501000^2*rv)/(287*(tc+273.15)^2))
  lr
}
# =========================================================================== #
# ********** Functions associated with runmicrobig ************************** #
# =========================================================================== #
#' crop dem by row and column
.croprast<-function(r,rw,cl,tilesize) {
  e<-ext(r)
  xmn<-e$xmin+(cl-1)*tilesize*res(r)[1]
  xmx<-e$xmin+cl*tilesize*res(r)[1]
  ymn<-e$ymax-rw*tilesize*res(r)[2]
  ymx<-e$ymax-(rw-1)*tilesize*res(r)[2]
  xmn<-ifelse(xmn<e$xmin,e$xmin,xmn)
  xmx<-ifelse(xmx>e$xmax,e$xmax,xmx)
  ymn<-ifelse(ymn<e$ymin,e$ymin,ymn)
  ymx<-ifelse(ymx>e$ymax,e$ymax,ymx)
  e2<-ext(c(xmn,xmx,ymn,ymx))
  r2<-terra::crop(r,e2)
  r2
}
#' crop dem by row and column with overlap
.croprasto<-function(r,rw,cl,tilesize,ovlap=10) {
  e<-ext(r)
  xmn<-e$xmin+(cl-1)*tilesize*res(r)[1]-ovlap*res(r)[1]
  xmx<-e$xmin+cl*tilesize*res(r)[1]+ovlap*res(r)[1]
  ymn<-e$ymax-rw*tilesize*res(r)[2]-ovlap*res(r)[2]
  ymx<-e$ymax-(rw-1)*tilesize*res(r)[2]+ovlap*res(r)[2]
  xmn<-ifelse(xmn<e$xmin,e$xmin,xmn)
  xmx<-ifelse(xmx>e$xmax,e$xmax,xmx)
  ymn<-ifelse(ymn<e$ymin,e$ymin,ymn)
  ymx<-ifelse(ymx>e$ymax,e$ymax,ymx)
  e2<-ext(c(xmn,xmx,ymn,ymx))
  r2<-terra::crop(r,e2)
  r2
}
#' Crop vegetation parameters
.vegcrop<-function(vegp,dtmi) {
  pai<-crop(.rast(vegp$pai,vegp$x),ext(dtmi))
  hgt<-crop(vegp$hgt,ext(dtmi))
  x<-crop(vegp$x,ext(dtmi))
  gsmax<-crop(vegp$gsmax,ext(dtmi))
  clump<-crop(.rast(vegp$clump,vegp$x),ext(dtmi))
  leafr<-crop(vegp$leafr,ext(dtmi))
  leafd<-crop(vegp$leafd,ext(dtmi))
  leaft<-crop(vegp$leaft,ext(dtmi))
  vegpi<-list(pai=as.array(pai),hgt=hgt,x=x,gsmax=gsmax,
              clump=as.array(clump),leafr=leafr,
              leafd=leafd,leaft=leaft)
  return(vegpi)
}
#' Crop soil parameters
.soilcrop<-function(soilc,dtmi) {
  soiltype<-crop(soilc$soiltype,ext(dtmi))
  groundr<-crop(soilc$groundr,ext(dtmi))
  soilci<-list(soiltype=soiltype,groundr=groundr)
  return(soilci)
}
# =========================================================================== #
# ********** Functions associated with post-processing data ***************** #
# =========================================================================== #
#' Calculates reference time for phase of diurnal temperature fluctuation
.t0f<-function(tc) {
  tc<-matrix(tc,ncol=24,byrow=T)
  tmx<-apply(tc,1,which.max)-1
  t0<-(tmx-6)%%24
  rep(t0*3600,each=24)
}
#' Calculates reference time for phase of diurnal temperature fluctuation (daily)
.t0fd<-function(tc) {
  tc<-matrix(tc,ncol=24,byrow=T)
  tmx<-apply(tc,1,which.max)-1
  t0<-(tmx-6)%%24
  t0<-t0*3600
  t0
}
#' Selects relevent daily min (or max) data from hourly climate arrays
.climds<-function(climarray,sel) {
  climo<-list(temp=climarray$temp[,,sel],
              relhum=climarray$relhum[,,sel],
              pres=climarray$pres[,,sel],
              swrad=climarray$swrad[,,sel],
              difrad=climarray$difrad[,,sel],
              skyem=climarray$skyem[,,sel],
              windspeed=climarray$windspeed[,,sel],
              winddir=climarray$winddir[,,sel])
  climo
}
# =========================================================================== #
# ************* Functions associated with bioclim variables ***************** #
# =========================================================================== #
#' identifies which entries of weather time series should be selected
.biosel<-function(tme, tc) {
  tch<-matrix(tc,ncol=24,byrow=TRUE)
  tcd<-apply(tch,1,mean)
  tmd<-matrix(as.numeric(tme),ncol=24,byrow=TRUE)
  tmd<-as.POSIXlt(apply(tmd,1,mean),origin="1970-01-01 00:00",tz="UTC")
  # Median temperature of the month
  sel_med<-0
  for (mth in 1:12) {
    s<-which(tmd$mon+1==mth)
    o<-order(tcd[s])
    n<-trunc(length(o)/2)
    s2<-s[1]-1+o[n]
    sel_med<-c(sel_med,s2)
  }
  sel_med<-sel_med[-1]
  # Maximum and minimum  (median of multiple years)
  yrs<-unique(tmd$year)
  sel_max<-0
  sel_min<-0
  for (y in 1:length(yrs)) {
    s<-which(tmd$year==yrs[y])
    sel_max<-c(sel_max,which.max(tcd[s])+s[1]-1)
    sel_min<-c(sel_min,which.min(tcd[s])+s[1]-1)
  }
  sel_max<-sel_max[-1]
  sel_min<-sel_min[-1]
  o1<-order(tcd[sel_max])
  o2<-order(tcd[sel_min])
  sel_max<-sel_max[o1]
  sel_min<-sel_min[o2]
  n<-trunc(length(sel_max)/2)+1
  sel_max<-sel_max[n]
  sel_min<-sel_min[n]
  # Expand out to hour
  seld<-c(sel_med,sel_max,sel_min)
  selh<-rep((seld-1)*24,each=24)+rep(c(1:24),length(seld))
  return(list(selh=selh,seld=seld))
}
#' identifies which entries of bioclim time series should be selected for a given quarter
.selquarter<-function(tmeh,mon)  {
  mths<-c(mon-1,mon,mon+1)%%12
  mths[mths==0]<-12
  sel<-0
  for (i in 1:3) {
    seli<-which(tmeh$mon+1==mths[i])
    sel<-c(sel,seli)
  }
  return(sel[-1])
}
#' Run micropoint for bioclim calculations
.biomicropoint<-function(weather,precip,tme,vegp,soilc,reqhgt,windhgt,soilm,dTmx,maxiter,lat=NA,long=NA) {
  soilparams<-micropoint::soilparams # dirty fix
  if (class(weather) == "data.frame") {
    # Subset right values
    tme<-as.POSIXlt(weather$obs_time,tz="UTC")
    tc<-weather$temp
    sel<-.biosel(tme, tc)
    # Run full soil model
    if (class(soilc$soiltype) == "PackedSpatRaster") soilc$soiltype<-rast(soilc$soiltype)
    sn<-.getmode(as.vector(soilc$soiltype))
    soiltype<-soilparameters$Soil.type[sn]
    # Convert weather data format
    w2<-weather
    w2$lwdown<-weather$skyem*5.67*10^-8*(weather$temp+273.15)^4
    w2$precip<-rep(precip,each=24)/24
    w2$swdown<-weather$swrad
    w2$swrad<-NULL
    w2$skyem<-NULL
    if (class(soilm)=="logical") soilm<-micropoint::soilmmodel(w2, soiltype)
    # Subset weather and precipitation
    weather2<-weather[sel$selh,]
    precip2<-precip[sel$seld]
    soilm2<-soilm[sel$selh]
    # Run point model
    micropoint<-runpointmodel(weather2,precip2,reqhgt,vegp,soilc,windhgt,soilm2,dTmx,maxiter,yearG=FALSE,lat,long)
  } else {
    tc<-apply(weather$temp,3,mean, na.rm = TRUE)
    sel<-.biosel(tme, tc)
    # Run full soil model
    if (class(soilc$soiltype) == "PackedSpatRaster") soilc$soiltype<-rast(soilc$soiltype)
    sn<-.getmode(as.vector(soilc$soiltype))
    soiltype<-soilparameters$Soil.type[sn]
    # Create weather data.frame
    ea<-.satvap(weather$temp)*weather$relhum/100
    rh<-(apply(ea,3,mean, na.rm = TRUE)/.satvap(tc))*100
    rh[rh>100]<-100
    w2<-data.frame(obs_time=tme,
                   temp=tc,
                   relhum=rh,
                   pres=apply(weather$pres,3,mean,na.rm=TRUE),
                   swdown=apply(weather$swrad,3,mean,na.rm=TRUE),
                   difrad=apply(weather$difrad,3,mean,na.rm=TRUE),
                   lwdown=apply(weather$skyem,3,mean,na.rm=TRUE)*5.67*10^-8*(tc+273.15)^4,
                   windspeed=apply(weather$windspeed,3,mean,na.rm=TRUE),
                   winddir=0,
                   precip=rep(apply(precip,3,mean,na.rm=TRUE),each=24))
    if (class(soilm)=="logical") soilm<-micropoint::soilmmodel(w2, soiltype)
    # Subset weather and precipitation
    weather2<-list(temp=weather$temp[,,sel$selh],
                   relhum=weather$relhum[,,sel$selh],
                   pres=weather$pres[,,sel$selh],
                   swrad=weather$swrad[,,sel$selh],
                   difrad=weather$difrad[,,sel$selh],
                   skyem=weather$skyem[,,sel$selh],
                   windspeed=weather$windspeed[,,sel$selh],
                   winddir=weather$winddir[,,sel$selh])
    precip2<-precip[,,sel$seld]
    soilm<-soilm[sel$selh]
    soilmo<-array(rep(soilm,each=dim(precip)[1]*dim(precip)[2]),dim=c(dim(precip)[1],dim(precip)[2],length(soilm)))
    micropoint <- runpointmodela(weather2,precip2,tme[sel$selh],reqhgt,vegp,soilc,windhgt,soilmo,dTmx,maxiter)
  }
  return(micropoint)
}
#' Check SpatRasts of model input
.checkbiginputs<-function(dtm,vegp,soilc) {
  .cr<-function(r,dtm) {
    sel<-which(is.na(.is(dtm))==FALSE & is.na(.is(r)))
    length(sel)
  }
  pr<-.rast(vegp$pai[,,1],dtm)
  cr<-.rast(vegp$clump[,,1],dtm)
  le<-c(.cr(pr,dtm),.cr(vegp$hgt,dtm),.cr(vegp$x,dtm),.cr(vegp$gsmax,dtm),
        .cr(vegp$leafr,dtm),.cr(cr,dtm),.cr(vegp$leafd,dtm),.cr(vegp$leaft,dtm))
  nms<-c("pai","hgt","x","gsmax","leafr","clump","leafd","leaft")
  if (max(le)>0) {
    sel<-which(le>0)
    stop(paste0("The following layers of vegp contain NA that are not NA in dtm: ",nms[sel],"\n"))
  }
  le<-c(.cr(soilc$soiltype,dtm),.cr(soilc$groundr,dtm))
  if (max(le)>0) {
    sel<-which(le>0)
    stop(paste0("The following layers of soilc contain NA that are not NA in dtm: ",nms[sel],"\n"))
  }
}
# =========================================================================== #
# *************** Functions associated with snow modelling ****************** #
# =========================================================================== #
#' Get snow environment form lat and long
.extractsnowenv <- function(lat, long) {
  xy<-as.matrix(data.frame(x=long,y=lat))
  int<-as.numeric(extract(rast(snowenv),xy))
  env<-"Alpine"
  if (int == 1) env<-"Maritime"
  if (int == 2) env<-"Prairie"
  if (int == 3) env<-"Tundra"
  if (int == 4) env<-"Taiga"
  return(env)
}
#' Calculates snow density parameters from snow environment
.snowdens<-function(snowenv="Tundra") {
  densfun<-c(0.5975,0.2237,0.0012,0.0038)
  if (snowenv == "Maritime") densfun<-c(0.5979,0.2578,0.001,0.0038)
  if (snowenv == "Prairie") densfun<-c(0.594,0.2332,0.0016,0.0031)
  if (snowenv == "Tundra") densfun<-c(0.363,0.2425,0.0029,0.0049)
  if (snowenv == "Taiga") densfun<-c(0.217,0.217,0,0)
  densfun
}
#' Calculate zero-plane displacement for snow
.szeroplanedis<-function(h,pai,snowdepth) {
  mu<-(h-snowdepth)/h
  mu[mu<0]<-0
  pai<-pai*mu
  pai[pai<0.001]<-0.001
  hs<-h-snowdepth
  hs[hs<0]<-0
  d<-snowdepth+(1-(1-exp(-sqrt(7.5*pai)))/sqrt(7.5*pai))*hs
  return(d)
}
#' Smooth matrix by xyf
.applyxyf<-function(r,xyf) {
  if (is.na(xyf) == FALSE) {
    if (xyf > 1) {
      r<-as.matrix(.smr(rast(r),xyf),wide=TRUE)
    }
  } else {
    r<-r*0+mean(r,na.rm=T)
  }
  return(r)
}
#' Calculate roughness length for snow
.sroughlength<-function(h,pai,snowdepth,d=NA,psi_h=0) {
  if (class(d)[1]=="logical") d<-.szeroplanedis(h,pai,snowdepth)
  mu<-(h-snowdepth)/h
  mu[mu<0]<-0
  pai<-pai*mu
  pai[pai<0.001]<-0.001
  Be<-sqrt(0.003+(0.2*pai)/2)
  zm<-(h-d)*exp(-0.4/Be)*exp(psi_h)
  zm[zm<0.0005]<-0.0005
  return(zm)
}
#' Calculate snow albedo
.snowalb <- function(age) {
  snowalb<-(-9.8740*log(age/24)+78.3434)/100
  snowalb[snowalb>0.95]<-0.95
  return(snowalb)
}
#' Calculate integrated diabatic correction coefficient for momentum
.dpsim<-function(ze) {
  # unstable
  x<-(1-15*ze)^0.25
  psi_m<-log(((1+x)/2)^2*(1+x^2)/2)-2*atan(x)+pi/2
  # stable
  s<-which(ze>=0)
  p2<- -4.7*ze
  psi_m[s]<-p2[s]
  psi_m[psi_m < -4] <- -4
  psi_m[psi_m > 3] <- 3
  psi_m
}
#' Calculate integrated diabatic correction coefficient for heat
.dpsih<-function(ze) {
  # unstable
  y<-(1-9*ze)^0.5
  psi_h<-log(((1+y)/2)^2)
  # stable
  s<-which(ze>=0)
  p2<- -(4.7*ze)/0.74
  psi_h[s]<-p2[s]
  psi_h[psi_h < -4] <- -4
  psi_h[psi_h > 3] <- 3
  psi_h
}
#' Calculate pai above snow surface where snow depth is an array or matrix
.epaif<-function(pai,hgt,snowd) {
  epai<-((hgt-snowd)/hgt)*pai
  epai[epai<0]<-0
  epai[is.na(epai)]<-0
  epai
}

#' Calculates wind speed within canopy given vertical profile within the canopy
.meancanopywind<-function(uz,pai,hgt,d,zm,snowd=0,zu=2) {
  # Apply shelter coefficient
  uf<-suppressWarnings((0.4*uz)/log((zu-d)))
  ehgt<-pmax(hgt,snowd)
  uh<-suppressWarnings((uf/0.4)*log((ehgt-d)/zm))
  s<-which(uh<uf)
  uh[s]<-uf[s]
  Be<-0.205*pai^0.445+0.1
  a<-pai/hgt
  Lc<-(0.25*a)^-1
  Lm<-2*Be^3*Lc
  a<-(Be*hgt)/Lm
  ef<-exp(a*((snowd/hgt)-1))
  ha<-hgt/a
  int<-ha*(1-ef)
  ur<-int/(hgt-snowd)
  ur[is.na(ur)]<-1
  ur[ur>1]<-1
  u<-uh*ur
  sel<-which(is.na(u))
  u[sel]<-uf[sel]
  sel<-which(u<uf)
  u[sel]<-uf[sel]
  u
}
#' Calculates SWE (mm) of intercepted canopy snow (inputs as follows:)
#' prec: daily precipitation (mm),
#' uz: wind speed at height zu  (m/s),
#' gsnowd: ground snow depth (cm),
#' Lt: snow load in previous time step (mm SWE)
#' pai: plant area index
#' hgt: canopy height (m)
#' Sh: branch snow load coefficient (kg/m^2). 6.6 for pine and 5.9 kg for spruce.
#' zu: height above ground of wind speed measurement uz (m)
.canopysnowint<-function(precd,uz,pai,hgt,gsnowd,Lt,Sh=6.6,zu=2) {
  # Calculate mean canopy wind speed
  d<-.szeroplanedis(hgt,pai,gsnowd)
  zm<-.roughlength(hgt,pai,gsnowd,d)
  uc<-.meancanopywind(uz,dtm,pai,hgt,d,zm,gsnowd,zu)
  # Calculate parameters for canopy snow interception
  epai<-.epaif(pai,hgt,gsnowd) # pai above snow surface
  Cc<-1-exp(-0.8*epai) # Snow leaf contact
  S<-Sh*(0.26+46/375) # Branch snow load
  Lstr<-S*0.8*epai  # Max canopy load (mm SWE)
  ehgt<-hgt-gsnowd
  ehgt[ehgt<0]<-0
  btm<-1-Cc*((uc*hgt)/(800))
  Cp<-Cc/btm  # maximum plant area of the snow±leaf contact per unit area of ground
  k<-Cp/Lstr
  k[is.na(k)]<-1
  Iinit<-(Lstr-Lt)*(1-exp(-k*precd))
  Int<-Iinit*0.678
  sel<-which(Int>precd)
  Int[sel]<-precd[sel]
  return(Int)
}
