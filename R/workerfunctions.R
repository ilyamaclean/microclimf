#' Check if input is a SpatRaster or PackedSpatRaster and convert to matrix if it is
#' @import terra
.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r<-rast(r)
  as.matrix(r,wide=TRUE)
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
#' convert radians to degrees
.ad <-function(r) {
  r*(180/pi)
}
#' convert degrees to radians
.ar <- function(d) {
  d*(pi/180)
}
#' Create SpatRaster object using a template
#' @import terra
.rast <- function(m,tem) {
  r<-rast(m)
  ext(r)<-ext(tem)
  crs(r)<-crs(tem)
  r
}
#' unpack model inputs
.unpack<-function(dtm,vegp,soilc) {
  if (class(dtm)[1] == "PackedSpatRaster") dtm<-rast(dtm)
  if (class(vegp$hgt)[1] == "PackedSpatRaster") vegp$hgt<-rast(vegp$hgt)
  if (class(vegp$x)[1] == "PackedSpatRaster") vegp$x<-rast(vegp$x)
  if (class(vegp$gsmax)[1] == "PackedSpatRaster") vegp$gsmax<-rast(vegp$gsmax)
  if (class(vegp$leafr)[1] == "PackedSpatRaster") vegp$leafr<-rast(vegp$leafr)
  if (class(vegp$clump)[1] == "PackedSpatRaster") vegp$clump<-rast(vegp$clump)
  if (class(vegp$leafd)[1] == "PackedSpatRaster") vegp$leafd<-rast(vegp$leafd)
  if (class(soilc$soiltype)[1] == "PackedSpatRaster") soilc$soiltype<-rast(soilc$soiltype)
  if (class(soilc$groundr)[1] == "PackedSpatRaster") soilc$groundr<-rast(soilc$groundr)
  crs(vegp$hgt)<-crs(dtm)
  crs(vegp$x)<-crs(dtm)
  crs(vegp$gsmax)<-crs(dtm)
  crs(vegp$leafr)<-crs(dtm)
  crs(vegp$clump)<-crs(dtm)
  crs(vegp$leafd)<-crs(dtm)
  crs(soilc$soiltype)<-crs(dtm)
  crs(soilc$groundr)<-crs(dtm)
  return(list(dtm=dtm,vegp=vegp,soilc=soilc))
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
.unpackpai <- function(pai,n) {
  x<-dim(pai)[3]
  i<-floor((x/n)*c(0:(n-1))+1)
  pai<-pai[,,i]
  pai
}
.soilinit <- function(soilc) {
  soiltype<-soilc$soiltype
  u<-unique(as.vector(.is(soiltype)))
  u<-u[is.na(u)==F]
  rho<-array(NA,dim=dim(soiltype)[1:2])
  Vm<-rho; Vq<-rho; Mc<-rho;
  Smax<-rho; Smax<-rho
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
  }
  return(list(rho=rho,Vm=Vm,Vq=Vq,Mc=Mc,psi_e=psi_e,soilb=soilb,Smax=Smax))
}
# Get soil coefficients for temperature model
.soilcoefs <- function(soilc) {
  soiltype<-soilc$soiltype
  u<-unique(as.vector(.is(soiltype)))
  u<-u[is.na(u)==F]
  int<-array(NA,dim=dim(soiltype)[1:2])
  t1<-int; t2<-int;	t3<-int; t4<-int; t5<-int; t6<-int; t7<-int
  for (i in 1:length(u)) {
    sel<-which(soilparameters$Number==u[i])
    sop<-soilparameters[sel,]
    sel<-which(.is(soiltype)==u[i])
    int[sel]<-sop$int
    t1[sel]<-sop$t1
    t2[sel]<-sop$t2
    t3[sel]<-sop$t3
    t4[sel]<-sop$t4
    t5[sel]<-sop$t5
    t6[sel]<-sop$t6
    t7[sel]<-sop$t7
  }
  return(list(int=int,t1=t1,t2=t2,t3=t3,t4=t4,t5=t5,t6=t6,t7=t7))
}
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
                      (101+sin(azi)*step^2):(y+100+sin(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(step^2))
  }
  horizon
}
#' Calculates radiation extinction coefficient for canopy
.cank <- function(x,sa,si) {
  # Raw k
  zen<-pi/2-sa
  k<-sqrt((x^2+(tan(zen)^2)))/(x+1.774*(x+1.182)^(-0.733))
  # k dash
  kd<-k*cos(zen)/si
  sel<-which(si==0)
  kd[sel]<-1
  return(list(k=k,kd=kd))
}
#' Calculates direct radiation transmission through canopy
.cantransdir <- function(l, k, ref = 0.23, clump = 0) {
  f<-1/(1-clump)
  s<-sqrt(1-ref)
  ks<-k*s
  tr<-(1-clump)*exp(-ks*l*f)+clump
  tr
}
#' Calculates diffuse radiation transmission through canopy
.cantransdif <- function(l, ref = 0.23, clump = 0) {
  f<-1/(1-clump)
  s<-sqrt(1-ref)
  tr<-(1-clump)*exp(-s*l*f)+clump
  tr
}
#' Calculates proportion of sunlight leaves in canopy
.psunlit <- function(l, k, clump = 0) {
  f <- 1/(1-clump)
  Ls<-(1-exp(-k*l*f))/k
  Lp <- Ls/l
  Lp<-(1-clump)*Lp+clump
  Lp[is.na(Lp)]<-1
  Lp
}
#' Calculates equivelent of proportion of sunlight leaves for diffuse radiation
.psunlitd <- function(l, clump = 0) {
  f <- 1/(1-clump)
  Ls<-(1-exp(-l*f))
  Lp <- Ls/l
  Lp<-(1-clump)*Lp+clump
  Lp[is.na(Lp)]<-1
  Lp
}
#' Calculate saturated vapour pressure
.satvap <- function(tc) {
  e0<-(tc<0)*610.78/1000+(tc>=0)*611.2/1000
  L <- (tc<0)*2.834*10^6+(tc>=0)*((2.501*10^6)-(2340*tc))
  T0<-(tc<0)*273.15+(tc>=0)*273.15
  estl<-e0*exp((L/461.5)*(1/T0-1/(tc+273.15)))
  estl
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
#' Calculates bulk surface stomatal conductivity
.layercond <- function(Rsw, gsmax, q50 = 100) {
  rpar<-Rsw*4.6
  gs <- (.rta(gsmax,dim(Rsw)[3])*rpar)/(rpar+q50)
  gs
}
#' Calculates slope of the saturated vapour pressure curve
.delta <- function(tc, radabs) {
  ta<-2.685+0.0006664*radabs
  tc<-tc+0.5*ta
  es1<-.satvap(tc-0.5)
  es2<-.satvap(tc+0.5)
  delta<-es2-es1
  delta
}
#' Calculate zero plane displacement
.zeroplanedis <- function(hgt, pai) {
  d<-(1-(1-exp(-sqrt(7.5*pai)))/sqrt(7.5*pai))*hgt
  d
}
#' Calculate roughness length governing momentum transfer
.roughlength <- function(hgt, pai, zm0 = 0.003) {
  d<-.zeroplanedis(hgt,pai)
  ur<-sqrt(zm0+(0.3*pai)/2)
  ur[ur>0.3]<-0.3
  zm<-(hgt-d)*exp(-0.4*(1/ur)-0.193)
  zm[zm<zm0]<-zm0
  zm
}
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
                      (101+sin(azi)*step^2):(y+100+sin(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(step^2))
    horizon[1:x,1:y]<-ifelse(horizon[1:x,1:y]<(hgt/step^2),0,horizon[1:x,1:y])
  }
  index<-1-atan(0.17*100*horizon)/1.65
  index
}
#' Calculates array of wind shelter coefficients in wdct directions
.windsheltera <- function(dsm, wdct = 8, whgt = 2, xyf = 10, hs = TRUE, cors=NA) {
  if (hs) {
    dct2<-2*wdct
  } else dct2<-wdct
  drs <- c(0:(dct2-1))*360/dct2
  a<-array(NA,dim=c(dim(dsm)[1:2],dct2))
  for (i in 1:dct2) {
    r<-.windcoef(dsm,drs[i],whgt,res(dsm)[1])
    if (xyf > 1) {
      a[,,i]<-.is(.smr(rast(r),xyf,cors))
    } else a[,,i]<-r
  }
  if (hs) {
    a2<-array(NA,dim=c(dim(dsm)[1:2],wdct))
    for (i in 1:wdct) {
      if (i == 1) {
        m<-0.5*a[,,1]+0.25*a[,,2]+0.25*a[,,dct2]
      } else m<-0.5*a[,,(i*2-1)]+0.25*a[,,i*2]+0.25*a[,,(i*2)-2]
      a2[,,i]<-m
    }
  } else a2<-a
  a2
}
#' Calculate wind speeds accounting for topographic sheltering
.windshelter <- function(u2, wdir, dsm, whgt = 2, wdct = 8, xyf = 10, hs = TRUE, wsa = NA, cors = NA) {
  # Create array of wind shelter coefficients
  if (class(wsa)[1] == "logical") wsa<-.windsheltera(dsm,wdct,whgt,xyf,hs,cors)
  # Apply array of wind shelter coefficients
  i<-360/dim(wsa)[3]
  sel<-(round(wdir/i,0))+1
  sel[sel==(dim(wsa)[3])+1]<-1
  sel[sel==0]<-dim(wsa)[3]
  wa<-wsa[,,sel]
  # Calculate wind speed at two metres above tallest vegetation
  if (whgt!=2) u2 <- u2*log(67.8*whgt-5.42)/4.87
  # Apply shelter coefficient
  uz<-.vta(u2,dsm)
  uz<-uz*wa
  uz
}
#' Calculates diabatic correction factors
.diabatic<-function(tc,uf,d,zm,maxhgt,H) {
  hes<-d+zm
  r<-rast(H[,,1])
  Tk<-tc+273.15
  st<- -(0.4*9.81*(2+maxhgt-d)*H)/(1241*Tk*uf^3)
  # H>0
  psi_h <- -2*log(1+((1-16*st)^0.5)/2)
  psi_m <- 0.6*psi_h
  # H < 0
  sel<-which(H<0)
  psi_h[sel]<-6*log(1+st[sel])
  psi_m[sel]<-psi_h[sel]
  psi_m<-.lim(psi_m,4,up=TRUE)
  psi_m<-.lim(psi_m,-4)
  psi_h<-.lim(psi_h,4,up=TRUE)
  psi_h<-.lim(psi_h,-4)
  return(list(psi_h=psi_h,psi_m=psi_m))
}
# Calculates below canopy mixing length
.mixinglength <- function(hgt, pai, zm0 = 0.003, mnlm = 0.04) {
  d<-.zeroplanedis(hgt,pai)
  zm<-.roughlength(hgt,pai,zm0)
  l_m<-(0.32*(hgt-d))/log((hgt-d)/zm)
  l_m[l_m<mnlm]<-mnlm
  l_m
}
# Ensures uz above canopy cannot drop below friction velocity
.minwind <- function(uz, uf) {
  sel<-which(uz<uf)
  uz[sel]<-uf[sel]
  sel<-which(is.na(uz))
  uz[sel]<-uf[sel]
  uz
}
# Calculate forced convection (used for calculating minimum conductivity)
.gfree<-function(leafd,H) {
  d<-0.71*leafd
  dT<-0.7045388*(d*H^4)^(1/5)
  gha<-0.0375*(dT/d)^0.25
  gha[gha<0.1]<-0.1
  gha
}
#' Calculates turbulent molar conductivity above canopy
.gturb<-function(uf,d,zm,z1,z0=NA,psi_h=0,gmin) {
  zh<-0.2*zm
  if (is.na(z0)[1]) {
    z0<-d+zh
  }
  xx<-(z1-d)/(z0-d)
  xx[xx<0.001]<-0.001
  ln<-log(xx)
  lnr<-ln*log((z1-d)/zh)
  psx<-lnr*psi_h
  g<-(0.4*43*uf)/(ln+psx)
  sel<-which(g<gmin)
  g[sel]<-gmin[sel]
  sel<-which(g>1e10)
  g[sel]<-1e10
  g
}
#' Calculates temperature and vapour above canopy
.TVabove<-function(TH,micro,z) {
  # Temperature
  lnr<-suppressWarnings(log((micro$maxhgt-micro$d)/(z-micro$d)))
  sel<-which(is.na(lnr))
  Ta<-micro$tc+(TH$H/(0.4*43*29.3))*lnr
  Ta[sel]<-TH$tcan[sel]
  # Vapour pressure
  ez<-micro$ea+((TH$L*micro$pk)/(765847.2*micro$uf))*lnr
  ez[sel]<-.satvap(TH$tcan[sel])
  return(list(Tz=Ta,ez=ez))
}
#' Calculate temperature profile below canopy
.tleaf <- function(tair,tground,es,eref,pk,gtt,gt0,gHa,gv,gL,Rabsl,leafdens,surfwet,soilrh=1) {
  # Air temperature expressed as leaf temperature
  aL<-(gtt*(tair+273.15)+gt0*(tground+273.15))/(gtt+gt0)
  bL<-(leafdens*gL)/(gtt+gt0)
  # Vapour pressures
  esoil<-.satvap(tground)*soilrh
  sb<-5.67*10^-8
  Rnet<-Rabsl-0.97*sb*(tair+273.15)^4
  delta<-.delta(tair,Rabsl)
  ae<-(gtt*eref+gt0*esoil+gv*es)/(gtt+gt0+gv)
  be<-(gv*delta)/(gtt+gt0+gv)
  # Sensible heat
  bH<-gHa*29.3
  # Latent heat
  aX<-((44526*gv)/pk)*(surfwet*es-ae)
  bX<-((44526*gv)/pk)*(surfwet*delta-be)
  aX[aX<0]<-0
  bX[bX<0]<-0
  # Emmited radiation
  aR<-sb*0.97*aL^4
  bR<-4*0.97*sb*(aL^3*bL+(tair+273.15)^3)
  # Leaf temperature
  dTL<-(Rabsl-aR-aX)/(1+bR+bX+bH)
  # tz pass 1
  tn<-aL-273.15+bL*dTL
  tleaf<-tn+dTL
  # new vapour pressure
  eanew<-ae+be*dTL
  tmin<-.dewpoint(eanew,tn)
  esnew<-.satvap(tn)
  eanew<-.lim(eanew,esnew,up=TRUE)
  eanew<-.lim(eanew,0.01)
  rh<-(eanew/esnew)*100
  # Set both tair and tleaf so as not to drop below dewpoint
  tleaf<-.lim(tleaf,tmin)
  tleaf<-.lim(tleaf,tair+20,up=TRUE)
  tleaf<-.lim(tleaf,80,up=TRUE)
  tn<-.lim(tn,tmin)
  tn<-.lim(tn,tair+20,up=TRUE)
  tn<-.lim(tn,80,up=TRUE)
  return(list(tleaf=tleaf,tn=tn,rh=rh))
}
#' Caclulates flow direction
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
    fd[x[i], y[i]] <- round(mean(which(md9 == min(md9, na.rm = T))), 0)
  }
  fd
}
#' Caclulates flow accumulation
.flowacc <- function (dtm) {
  dm<-.is(dtm)
  fd<-.flowdir(dm)
  fa<-fd*0+1
  o<-order(dm,decreasing=T,na.last=NA)
  for (i in 1:length(o)) {
    x <- arrayInd(o[i], dim(dm))[1]
    y <- arrayInd(o[i], dim(dm))[2]
    f <- fd[x, y]
    x2 <- x + (f - 1)%%3 - 1
    y2 <- y + (f - 1)%/%3 - 1
    if (x2 > 0 & x2 < dim(dm)[1] & y2 > 0 & y2 < dim(dm)[2])
      fa[x2, y2] <- fa[x, y] + 1
  }
  fa
}
#' Calculates topographic wetness index
.topidx <- function(dtm) {
  minslope<-atan(0.005/mean(res(dtm)))
  slope<-terrain(dtm,unit="radians")
  B<-.is(slope)
  B[B<minslope]<-minslope
  a<-.flowacc(dtm)
  a<-a*res(dtm)[1]* res(dtm)[2]
  tpx<- a/tan(B)
  r<-.rast(tpx,dtm)
  r
}
.multisoil <- function(micro,soilinit,tfact,mult,twi) {
  soilc<-micro$soilc
  soiltype<-soilc$soiltype
  u<-unique(as.vector(.is(soiltype)))
  u<-u[is.na(u)==F]
  dtm<-micro$dtm
  clim<-micro$climdata
  smout<-array(NA,dim=c(dim(dtm)[1:2],length(micro$prec)))
  for (i in 1:length(u)) {
    stp<-soilparameters$Soil.type[u[i]]
    st<-soilparameters$Number[u[i]]
    sel<-which(.is(soiltype)==st)
    alb<-mean(.is(soilc$groundr)[sel],na.rm=T)
    swrad<-(1-alb)*clim$swrad
    lwout<- 5.67e-8*0.95*(clim$temp+273.15)^4
    lwnet<-(1-clim$skyem)*lwout
    rnet<-swrad-lwnet
    rnet<-rnet*mult
    soilm<-soilmpredict(micro$prec,rnet,stp,soilinit)$soilm1
    Smin<-soilparameters$Smin[u[i]]
    Smax<-soilparameters$Smax[u[i]]
    sm<-soilmdistribute(soilm,dtm,Smin,Smax,tfact,twi)
    sar<-.rta(soiltype,length(soilm))
    sela<-which(sar==st)
    smout[sela]<-sm[sela]
  }
  smout
}
#' Calculates canopy boundary layer conductivity and conductivity to ground
.conductivityE<-function(micro,reqhgt,xyf=1,zf=NA,slr=NA,apr=NA,hor=NA,wsa=NA,maxhgt=NA) {
  if (micro$progress<2) {
    micro<-canopyrad(micro,slr,apr,hor)
  }
  if (micro$progress<3) {
    micro<-wind(micro,xyf=xyf,zf=zf,slr=slr,apr=apr,hor=hor,wsa=wsa,maxhgt=maxhgt)
  }
  # Calculate approximate H
  Hest<-0.5*(micro$canabs-micro$lwout)
  # Calculate approximate diabatic correction
  dbm<-.diabatic(micro$tc,micro$uf,micro$d,micro$zm,micro$maxhgt,Hest)
  # Calculate boundary layer conductivity
  pai<-micro$pai
  pai[pai<1]<-1
  lwi<-.rta(micro$leafd,dim(pai)[3])
  gmin<-.gfree(lwi,Hest)*2*pai
  gHa<-.gturb(micro$uf,micro$d,micro$zm,micro$maxhgt,psi_h=dbm$psi_h,gmin=gmin)
  # Calculate g0
  micro<-wind(micro,xyf,zf,dbm$psi_m,reqhgt,slr=slr,apr=apr,hor=hor,wsa=wsa,maxhgt=maxhgt)
  micro$gHa<-gHa
  micro$dbm<-dbm
  micro$progress<-4
  # Clean micro
  micro$vegx<-NULL
  micro$veghgt<-NULL
  micro$si<-NULL
  return(micro)
}
#' expand daily array to hourly array
.ehr<-function(a) {
  n<-dim(a)[1]*dim(a)[2]
  o1<-rep(c(1:n),24*dim(a)[3])
  o2<-rep(c(1:dim(a)[3]),each=24*n)-1
  o2<-o2*max(o1)
  o<-o1+o2
  ah<-rep(a,24)
  ah<-ah[o]
  ah<-array(ah,dim=c(dim(a)[1:2],dim(a)[3]*24))
  ah
}
#' Calculate soil relative humidity
.soilrh <- function(theta, b, Psie, Smax, tc) {
  matric <- -Psie*(theta/Smax)^-b
  hr<-exp((0.018*matric)/(8.31*(tc+273.15)))
  hr[hr>1]<-1
  hr
}
# Calculate absorbed reflected radiation
.radref<-function(micro,pai_a,pai_b,reqhgt) {
  # canopy reflected
  a<- (-sqrt(1-micro$lref)*micro$pai*(micro$k*micro$dirr+micro$difr))/
    (micro$vha*(micro$dirr+micro$difr))
  b<-micro$lref*(micro$pai/micro$vha)*(1-micro$lref)*0.5*micro$lref*(micro$dirr+micro$difr)
  rcan<-b*(exp(a*(micro$vha-reqhgt))*((exp(2*a*reqhgt-1)/(2*a))+(micro$vha-reqhgt)))
  rcan[is.na(rcan)]<-0
  rcan[rcan<0]<-0
  # ground reflected
  trdi<-.cantransdir(pai_b,micro$k,micro$lref,micro$clump)
  trdf<-.cantransdif(pai_b,micro$lref,micro$clump)
  rgr<-(1-micro$lref)*trdf*micro$gref*(trdi*micro$radm*micro$dirr+trdf*micro$svfa*micro$difr)
  # both
  rref<-rcan+rgr
  rref
}
# Calculate leaf absorbed radiation
.radabs<-function(micro,pai_a,reqhgt,tcan) {
  pai_b<-micro$pai-pai_a
  #' Calculate radiation absorbed below PAI
  trdi<-.cantransdir(pai_a,micro$k,micro$lref,micro$clump)
  trdf<-.cantransdif(pai_a,micro$lref,micro$clump)
  trlw<-.cantransdif(pai_a,0.03,micro$clump)
  # SW radiation
  rad_dir<-trdi*micro$dirr*micro$radm
  rad_dif<-trdf*micro$difr*micro$svfa
  rad_ref<-.radref(micro,pai_a,pai_b,reqhgt)
  radzsw<-0.5*((1-micro$lref)*(rad_dir+rad_dif)+rad_ref)
  # LW radiation
  radzlw<-0.5*0.97*(trlw*micro$lwsky+(1-trlw)*0.97*5.67*10^-8*(tcan+273.15)^4)+ # above
          0.5*0.97*5.67*10^-8*(micro$T0+273.15)^4  # ground
  return(list(radzsw=radzsw,radzlw=radzlw,rad_dir=rad_dir,rad_dif=rad_dif,rad_ref=rad_ref,trdf=trdf))
}
# Convert hourly weather to daily
.climtodaily<-function(climdata) {
  td<-matrix(climdata$temp,ncol=24,byrow=TRUE)
  smn<-apply(td,1,which.min)
  smx<-apply(td,1,which.max)
  ta<-(c(1:dim(td)[1])-1)*24
  climn<-climdata[(smn+ta),]
  climx<-climdata[(smx+ta),]
  return(list(climn=climn,climx=climx,smx=smx+ta,smn=smn+ta))
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
.expandtohour2<-function(amn,amx,v) {
  r<-rast(amn[,,1])
  mev<-matrix(v,ncol=24,byrow=T)
  mev<-rep(apply(mev,1,mean),each=24)
  mu<-v/mev
  ame<-.ehr((amn+amx)/2)
  ao<-.vta(mu,r)*ame
  ao
}
#' Expand daily climate variables to hourly (all)
.expandclim<-function(mout_mn,mout_mx,climdata) {
  climd<-.climtodaily(climdata)
  # Expand Tz
  Tz<-.expandtohour(mout_mn$Tz,mout_mx$Tz,climd$smn,climd$smx,climdata$temp)
  # Expand tleaf
  tleaf<-.expandtohour(mout_mn$tleaf,mout_mx$tleaf,climd$smn,climd$smx,climdata$temp)
  # Expand T0
  T0<-.expandtohour(mout_mn$T0,mout_mx$T0,climd$smn,climd$smx,climdata$temp)
  # Expand soilm
  soilm<-.ehr(mout_mn$soilm)
  # Expand relhum
  ea_mn<-(mout_mn$relhum/100)*.satvap(mout_mn$Tz)
  ea_mx<-(mout_mx$relhum/100)*.satvap(mout_mx$Tz)
  ea<-(ea_mn+ea_mx)/2
  ea<-.ehr(ea)
  relhum<-(ea/.satvap(T0))*100
  relhum[relhum>100]<-100
  # Expand windspeed
  windspeed<-.expandtohour2(mout_mn$windspeed,mout_mx$windspeed,climdata$windspeed)
  # Expand shortwave radiation
  tr<-.ehr(mout_mn$trdf)
  r<-rast(tr[,,1])
  rdi<-climdata$swrad-climdata$difrad
  raddir<-.vta(rdi,r)*tr
  raddif<-.vta(climdata$difrad,r)*tr
  # Expand longwave radiation
  sb<-5.67*10^-8
  em_mn<-mout_mn$radlw/(sb*(mout_mn$Tz+273.15)^4)
  em_mx<-mout_mx$radlw/(sb*(mout_mx$Tz+273.15)^4)
  em<-.expandtohour2(em_mn,em_mx,climdata$skyem)
  em[em>1]<-1
  radlw<-em*sb*(Tz+273.15)^4
  out<-list(Tz=Tz,tleaf=tleaf,T0=T0,soilm=soilm,relhum=relhum,windspeed=windspeed,
            raddir=raddir,raddif=raddif,radlw=radlw)
  return(out)
}
# Crop SpatRaster object for running runmicro_big
#' @import terra
.cropraster<-function(r,rw,cl) {
  e<-ext(r)
  xmn<-e$xmin+(cl-1)*100*res(r)[1]
  xmx<-e$xmin+cl*100*res(r)[1]
  ymn<-e$ymax-rw*100*res(r)[2]
  ymx<-e$ymax-(rw-1)*100*res(r)[2]
  xmn<-ifelse(xmn<e$xmin,e$xmin,xmn)
  xmx<-ifelse(xmx>e$xmax,e$xmax,xmx)
  ymn<-ifelse(ymn<e$ymin,e$ymin,ymn)
  ymx<-ifelse(ymx>e$ymax,e$ymax,ymx)
  e2<-ext(c(xmn,xmx,ymn,ymx))
  r2<-terra::crop(r,e2)
  r2
}
# Crop array for running runmicro_big
.croparray<-function(a,rw,cl) {
  dms<-dim(a)
  xmn<-(cl-1)*100+1
  xmx<-xmn+99
  xmx<-ifelse(xmx>dms[2],dms[2],xmx)
  ymn<-(rw-1)*100+1
  ymx<-ymn+99
  ymx<-ifelse(ymx>dms[2],dms[2],ymx)
  ao<-a[ymn:ymx,xmn:xmx,]
  ao
}
# Crop vegp for running runmicro_big
.vegpcrop<-function(vegp,rw,cl) {
  pai<-.croparray(vegp$pai,rw,cl)
  hgt<-.cropraster(vegp$hgt,rw,cl)
  x<-.cropraster(vegp$x,rw,cl)
  gsmax<-.cropraster(vegp$gsmax,rw,cl)
  leafr<-.cropraster(vegp$leafr,rw,cl)
  clump<-.cropraster(vegp$clump,rw,cl)
  leafd<-.cropraster(vegp$leafd,rw,cl)
  vegp2<-list(pai=pai,hgt=hgt,x=x,gsmax=gsmax,leafr=leafr,clump=clump,leafd=leafd)
  class(vegp2)<-"vegparams"
  vegp2
}
# Crop soilc for running runmicro_big
.soilccrop<-function(soilc,rw,cl) {
  soiltype<-.cropraster(soilc$soiltype,rw,cl)
  groundr<-.cropraster(soilc$groundr,rw,cl)
  soilc2<-list(soiltype=soiltype,groundr=groundr)
  class(soilc2)<-"soilcharac"
  soilc2
}
#' Converts daily outputs to integers
.outasint <- function(mout, weather, dtm, reqhgt, merid, dst) {
  asa<-function(a,r) {
    a<-round(a*r,0)
    a<-array(as.integer(a),dim=dim(a))
    a
  }
  mout_mn<-mout$mout_mn
  mout_mx<-mout$mout_mx
  mout_mn$Tz<-asa(mout_mn$Tz,100)
  mout_mx$Tz<-asa(mout_mx$Tz,100)
  mout_mn$T0<-NULL
  mout_mx$T0<-NULL
  if (reqhgt > 0) {
    mout_mn$tleaf<-asa(mout_mn$tleaf,100)
    mout_mx$tleaf<-asa(mout_mx$tleaf,100)
    mout_mn$relhum<-asa(mout_mn$relhum,1)
    mout_mx$relhum<-asa(mout_mx$relhum,1)
    mout_mn$windspeed<-asa(mout_mn$windspeed,100)
    mout_mx$windspeed<-asa(mout_mx$windspeed,100)
    mout_mn$raddir<-asa(mout_mn$raddir,1)
    mout_mx$raddir<-asa(mout_mx$raddir,1)
    mout_mn$raddif<-asa(mout_mn$raddif,1)
    mout_mx$raddif<-asa(mout_mx$raddif,1)
    mout_mn$radlw<-asa(mout_mn$radlw,1)
    mout_mx$radlw<-asa(mout_mx$radlw,1)
    mout_mn$soilm<-NULL
    mout_mx$soilm<-NULL
  }
  if (reqhgt == 0) {
    mout_mn$soilm<-asa(mout_mn$soilm,100)
    mout_mx$soilm<-asa(mout_mx$soilm,100)
  }
  mout<-list(mout_mn=mout_mn,mout_mx=mout_mx,weather=weather,dtm=dtm,reqhgt=reqhgt,
             merid=merid,dst=dst)
  mout
}
# Expansion of climate data to hourly ignoring T0 and soilm
.expandclim2<-function(mout_mn,mout_mx,climdata) {
  climd<-.climtodaily(climdata)
  # Expand Tz
  Tz<-.expandtohour(mout_mn$Tz,mout_mx$Tz,climd$smn,climd$smx,climdata$temp)
  # Expand tleaf
  tleaf<-.expandtohour(mout_mn$tleaf,mout_mx$tleaf,climd$smn,climd$smx,climdata$temp)
  # Expand relhum
  ea_mn<-(mout_mn$relhum/100)*.satvap(mout_mn$Tz)
  ea_mx<-(mout_mx$relhum/100)*.satvap(mout_mx$Tz)
  ea<-(ea_mn+ea_mx)/2
  ea<-.ehr(ea)
  relhum<-(ea/.satvap(Tz))*100
  relhum[relhum>100]<-100
  # Expand windspeed
  windspeed<-.expandtohour2(mout_mn$windspeed,mout_mx$windspeed,climdata$windspeed)
  # Expand shortwave radiation
  tr<-.ehr(mout_mn$trdf)
  r<-rast(tr[,,1])
  rdi<-climdata$swrad-climdata$difrad
  raddir<-.vta(rdi,r)*tr
  raddif<-.vta(climdata$difrad,r)*tr
  # Expand longwave radiation
  sb<-5.67*10^-8
  em_mn<-mout_mn$radlw/(sb*(mout_mn$Tz+273.15)^4)
  em_mx<-mout_mx$radlw/(sb*(mout_mx$Tz+273.15)^4)
  em<-.expandtohour2(em_mn,em_mx,climdata$skyem)
  em[em>1]<-1
  radlw<-em*sb*(Tz+273.15)^4
  out<-list(Tz=Tz,tleaf=tleaf,relhum=relhum,windspeed=windspeed,
            raddir=raddir,raddif=raddif,radlw=radlw)
  return(out)
}
#' Resample array
#' @import terra
.resa<-function(varin,r,dtm) {
  if (dim(varin)[1]!=dim(dtm)[1] | dim(varin)[2]!=dim(dtm)[2]) {
    b<-rast(varin)
    ext(b)<-ext(r)
    crs(b) <- crs(r)
    if (as.character(crs(b)) != as.character(crs(dtm))) {
      b <- project(b,dtm, method = 'near')
      b <- projectRaster(from = b, crs = crs(dtm), method = 'ngb')
    }
    b<-resample(b,dtm)
    a<-as.array(b)
    a<-round(a,3)
  } else a<-varin
  a
}
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
#' Lats and longs from SpatRaster object, including reprojection
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
.lapserate <- function(tc, rh, pk) {
  e0 <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
  ws <- 0.622 * e0 / pk
  ea <- e0 * (rh / 100)
  rv <- 0.622 * ea / (pk - ea)
  lr <- 9.8076 * (1 + (2501000 * rv) / (287 * (tc + 273.15))) / (1003.5 +
                                                                   (0.622 * 2501000 ^ 2 * rv) / (287 * (tc + 273.15) ^ 2))
  lr
}
#' convert climate array to data.frame of weather
.catoweather<-function(climarray,tme) {
  # ~~~~ Wind
  uwind<-climarray$windspeed*cos(climarray$winddir*(pi/180))
  vwind<-climarray$windspeed*sin(climarray$winddir*(pi/180))
  uw<-apply(uwind,3,mean,na.rm=T)
  vw<-apply(vwind,3,mean,na.rm=T)
  ws=sqrt(uw^2+vw^2)
  wd=(atan2(vw,uw)*(180/pi))%%360
  weather<-data.frame(obs_time=tme,
                      temp=apply(climarray$temp,3,mean,na.rm=T),
                      relhum=apply(climarray$relhum,3,mean,na.rm=T),
                      pres=apply(climarray$pres,3,mean,na.rm=T),
                      swrad=apply(climarray$swrad,3,mean,na.rm=T),
                      difrad=apply(climarray$difrad,3,mean,na.rm=T),
                      skyem=apply(climarray$skyem,3,mean,na.rm=T),
                      windspeed=ws,
                      winddir=wd)
  return(weather)
}
#' Corrects wind profile
.windcorrect <- function(uz, windhgt, maxhgt) {
  d<-0.08
  zm<-0.01476
  zh<-0.1*zm
  uf<-(0.4*uz)/log((windhgt-d)/zm)
  uo<-(uf/0.4)*log((maxhgt+2-d)/zm)
  uo
}
#' Calculates diurnal temperature fluctuation
.A0f<-function(tc) {
  tc<-matrix(tc,ncol=24,byrow=T)
  mn<-apply(tc,1,min)
  mx<-apply(tc,1,max)
  A0<-(mx-mn)/2
  rep(A0,each=24)
}
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
