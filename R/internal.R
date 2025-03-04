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
#' Fill NAs in raster
#' @import terra
.fillr<-function(r,msk) {
  dd<-min(dim(r)[1:2])
  af<-c(2,4,8,16,32,64,128,256)
  af<-af[af<dd]
  ro<-list()
  for (i in 1:length(af)) {
    ro[[i]]<-resample(aggregate(r,af[i],fun="mean",na.rm=TRUE),r)
  }
  m<-.is(r)
  me<-mean(m,na.rm=TRUE)
  for (i in 1:length(af)) {
    s<-which(is.na(m))
    m2<-.is(ro[[i]])
    m[s]<-m2[s]
  }
  s<-which(is.na(m))
  m[s]<-me
  ro<-.rast(m,r)
  ro<-mask(ro,msk)
  ro
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
#' @import terra
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
#' convert degrees to radians
.ar <- function(d) {
  d*(pi/180)
}
# =========================================================================== #
# ******************* Used by point model *********************************** #
# =========================================================================== #
#' unpack model inputs
.unpack<-function(dtm,vegp,soilc) {
  if (class(dtm)[1] == "PackedSpatRaster") dtm<-rast(dtm)
  if (class(vegp$pai)[1] == "PackedSpatRaster") {
    vegp$pai<-rast(vegp$pai)
    crs(vegp$pai)<-crs(dtm)
  }
  if (class(vegp$hgt)[1] == "PackedSpatRaster") vegp$hgt<-rast(vegp$hgt)
  if (class(vegp$x)[1] == "PackedSpatRaster") vegp$x<-rast(vegp$x)
  if (class(vegp$gsmax)[1] == "PackedSpatRaster") vegp$gsmax<-rast(vegp$gsmax)
  if (class(vegp$leafr)[1] == "PackedSpatRaster") vegp$leafr<-rast(vegp$leafr)
  if (class(vegp$clump)[1] == "PackedSpatRaster") {
    vegp$clump<-rast(vegp$clump)
    crs(vegp$clump)<-crs(dtm)
  }
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
#' point model soil temperature below
.soilbelowT<-function(dfo,reqhgt) {
  # Calculate n
  n<- -118.35*reqhgt/dfo$DDp
  nmn<-floor(min(n))
  nmx<-ceiling(max(n))
  # Calculate temperatures
  Tnmn<-manCpp(dfo$Tg,nmn)
  Tnmx<-manCpp(dfo$Tg,nmx)
  # Calculate wgt (from Tnmx)
  wgt<-(n-nmn)/(nmx-nmn)
  Tb<-wgt*Tnmx+(1-wgt)*Tnmn
  # Calculate tmean weight
  wgt<- 0.041596*(reqhgt/mean(dfo$DDp))+0.87142
  Tbz<-wgt*Tb+(1-wgt)*mean(dfo$Tg)
  return(Tbz)
}
# replicate required number of layers to return array
.intr<-function(r,n,subs) {
  nr<-dim(r)[3]
  s<-round(seq(0.50001,nr+0.5,length.out=n),0)
  s[s>nr]<-nr
  s[s<1]<-1
  s<-s[subs]
  a<-as.array(r)
  ao<-a[,,s]
}
# max dims of vegp
.vegpdmx<-function(vegp) {
  dmx<-1
  if (dim(vegp$pai)[3] > dmx) dmx<-dim(vegp$pai)[3]
  if (dim(vegp$hgt)[3] > dmx) dmx<-dim(vegp$hgt)[3]
  if (dim(vegp$x)[3] > dmx) dmx<-dim(vegp$x)[3]
  if (dim(vegp$gsmax)[3] > dmx) dmx<-dim(vegp$gsmax)[3]
  if (dim(vegp$leafr)[3] > dmx) dmx<-dim(vegp$leafr)[3]
  if (dim(vegp$clump)[3] > dmx) dmx<-dim(vegp$clump)[3]
  if (dim(vegp$leafd)[3] > dmx) dmx<-dim(vegp$leafd)[3]
  if (dim(vegp$leaft)[3] > dmx) dmx<-dim(vegp$leaft)[3]
  return(dmx)
}
#' Sort out temporal nature of vegetation
.sortvegp<-function(vegp,method="R",n=8760,subs=c(1:8760)) {
  # check whether packed etc
  .unpack<-function(r,rte) {
    if (class(r)=="PackedSpatRaster") r<-rast(r)
    if (class(vegp$hgt) != "SpatRaster") r<-.rast(rte)
    return(r)
  }
  # test
  if (class(vegp$hgt)=="PackedSpatRaster") vegp$hgt<-rast(vegp$hgt)
  if (class(vegp$hgt) != "SpatRaster") stop("vegp$hgt must be a PackedSpatRaster or a SpatRaster")
  # unpack if Packed
  vegp$pai<-.unpack(vegp$pai,vegp$hgt)
  vegp$x<-.unpack(vegp$x,vegp$hgt)
  vegp$gsmax<-.unpack(vegp$gsmax,vegp$hgt)
  vegp$leafr<-.unpack(vegp$leafr,vegp$hgt)
  vegp$clump<-.unpack(vegp$clump,vegp$hgt)
  vegp$leafd<-.unpack(vegp$leafd,vegp$hgt)
  vegp$leaft<-.unpack(vegp$leaft,vegp$hgt)
  if (method == "P") {
    vegp<-c(mean(.is(vegp$hgt),na.rm=TRUE),
            mean(.is(vegp$pai),na.rm=TRUE),
            mean(.is(vegp$x),na.rm=TRUE),
            mean(.is(vegp$clump),na.rm=TRUE),
            mean(.is(vegp$leafr),na.rm=TRUE),
            mean(.is(vegp$leaft),na.rm=TRUE),
            mean(.is(vegp$leafd),na.rm=TRUE),
            0.97,
            mean(.is(vegp$gsmax),na.rm=TRUE),
            100)
  } else if (method == "R") {
    vegp$hgt<-.intr(vegp$hgt,n,subs)
    vegp$pai<-.intr(vegp$pai,n,subs)
    vegp$x<-.intr(vegp$x,n,subs)
    vegp$gsmax<-.intr(vegp$gsmax,n,subs)
    vegp$leafr<-.intr(vegp$leafr,n,subs)
    vegp$clump<-.intr(vegp$clump,n,subs)
    vegp$leafd<-.intr(vegp$leafd,n,subs)
    vegp$leaft<-.intr(vegp$leaft,n,subs)
  }  else {
    # Make all the layers as expanded as the biggest layer
    # Calculate max dims
    dmx<-.vegpdmx(vegp)
    # create subset
    s<-round(seq(0.50001,dmx+0.5,length.out=n),0)
    s[s>dmx]<-dmx
    s[s<1]<-1
    subs2<-unique(s[subs])
    vegp$hgt<-.intr(vegp$hgt,dmx,subs2)
    vegp$pai<-.intr(vegp$pai,dmx,subs2)
    vegp$x<-.intr(vegp$x,dmx,subs2)
    vegp$gsmax<-.intr(vegp$gsmax,dmx,subs2)
    vegp$leafr<-.intr(vegp$leafr,dmx,subs2)
    vegp$clump<-.intr(vegp$clump,dmx,subs2)
    vegp$leafd<-.intr(vegp$leafd,dmx,subs2)
    vegp$leaft<-.intr(vegp$leaft,dmx,subs2)
    # identify which layer is associated with which time sequence
    s<-round(seq(0.50001,dmx+0.5,length.out=n),0)
    s<-s[subs]
    vegp$lsubs<-s
  }
  return(vegp)
}
# sort out vegetation for bioclim
.sortvegp2<-function(vegp,selw,vegpisannual,n) {
  if (vegpisannual) {
    seld<-selw$seld%%365
    nd<-365
  } else {
    seld<-selw$seld
    nd<-round(n/24,0)
  }
  vegp2<-list()
  for (i in 1:length(vegp)) {
    dmx<-dim(vegp[[i]])[3]
    if (dmx == 1) {
      m<-.is(vegp[[i]])
      vegp2[[i]]<-array(m,dim=c(dim(m),14))
    } else {
      s<-round(seq(0.50001,dmx+0.5,length.out=nd),0)
      s[s>dmx]<-dmx
      s[s<1]<-1
      s<-c(s,s[length(s)])
      subs<-s[seld]
      a<-.is(vegp[[i]])
      vegp2[[i]]<-a[,,subs]
    }
  }
  names(vegp2)<-names(vegp)
  return(vegp2)
}
#' initialises soil parameters
.soilinit <- function(soilc) {
  .sortc<-function(varn,invar,soilc,u) {
    nms1<-names(soilc)
    s<-which(nms1 == varn)
    if (length(s) > 0) {
      invar<-.is(soilc[[s]])
    } else {
      for (i in 1:length(u)) {
        sel<-which(soilparameters$Number==u[i])
        sop<-soilparameters[sel,]
        nms2<-names(sop)
        s<-which(nms2 == varn)
        soiltype<-soilc$soiltype
        sel<-which(.is(soiltype)==u[i])
        invar[sel]<-as.numeric(sop[s])
      }
    }
    return(invar)
  }
  soiltype<-soilc$soiltype
  u<-unique(as.vector(.is(soiltype)))
  u<-u[is.na(u)==F]
  invar<-array(NA,dim=dim(soiltype)[1:2])
  rho<-.sortc("rho",invar,soilc,u)
  Vm<-.sortc("Vm",invar,soilc,u)
  Vq<-.sortc("Vq",invar,soilc,u)
  Mc<-.sortc("Mc",invar,soilc,u)
  psi_e<-.sortc("psi_e",invar,soilc,u)
  soilb<-.sortc("b",invar,soilc,u)
  Smax<-.sortc("Smax",invar,soilc,u)
  Smin<-.sortc("Smin",invar,soilc,u)
  return(list(rho=rho,Vm=Vm,Vq=Vq,Mc=Mc,psi_e=psi_e,soilb=soilb,Smax=Smax,Smin=Smin))
}
# Sort out soil parameters
.sortsoilc<-function(soilc,method="P") {
  soilcl<-.soilinit(soilc)
  if (method == "P") {
    sn<-.getmode(as.vector(soilc$soiltype))
    groundp_p<-c(.getmode(.is(soilc$groundr)),
                 0,180,0.97,
                 .getmode(soilcl$rho),
                 .getmode(soilcl$Vm),
                 .getmode(soilcl$Vq),
                 .getmode(soilcl$Mc),
                 .getmode(soilcl$soilb),
                 .getmode(soilcl$psi_e),
                 .getmode(soilcl$Smax),
                 .getmode(soilcl$Smin),
                 soilparamsp$alpha[sn],
                 soilparamsp$n[sn],
                 soilparamsp$Ksat[sn])
    soilcl<-groundp_p
  }
  return(soilcl)
}
.resamplev<-function(vegp,r) {
  af<-res(r)[1]/res(vegp$pai)[1]
  vegp$pai<-aggregate(vegp$pai,af,na.rm=TRUE)
  vegp$hgt<-aggregate(vegp$hgt,af,na.rm=TRUE)
  vegp$x<-aggregate(vegp$x,af,na.rm=TRUE)
  vegp$gsmax<-aggregate(vegp$gsmax,af,na.rm=TRUE)
  vegp$leafr<-aggregate(vegp$leafr,af,na.rm=TRUE)
  vegp$clump<-aggregate(vegp$clump,af,na.rm=TRUE)
  vegp$leafd<-aggregate(vegp$leafd,af,na.rm=TRUE)
  vegp$leaft<-aggregate(vegp$leaft,af,na.rm=TRUE)
  vegp$pai<-resample(vegp$pai,r)
  vegp$hgt<-resample(vegp$hgt,r)
  vegp$x<-resample(vegp$x,r)
  vegp$gsmax<-resample(vegp$gsmax,r)
  vegp$leafr<-resample(vegp$leafr,r)
  vegp$clump<-resample(vegp$clump,r)
  vegp$leafd<-resample(vegp$leafd,r)
  vegp$leaft<-resample(vegp$leaft,r)
  vegp
}
#' Converts array of climate data to a dataframe for a given cell
.todf<-function(climarray,i,j,tme,winddir) {
  climdf<-with(climarray,data.frame(obs_time=tme,
                                    temp=.is(temp)[i,j,],
                                    relhum=.is(relhum)[i,j,],
                                    pres=.is(pres)[i,j,],
                                    swdown=.is(swdown)[i,j,],
                                    difrad=.is(difrad)[i,j,],
                                    lwdown=.is(lwdown)[i,j,],
                                    windspeed=.is(windspeed)[i,j,],
                                    precip=.is(precip)[i,j,]))
  climdf$winddir<-winddir
  climdf
}
#' used by runpointmodela
.tovp<-function(vegp,i,j,vmx) {
  vegp_p<-c(mean(as.array(vegp$hgt)[i,j,]),
            mean(as.array(vegp$pai)[i,j,]),
            mean(as.array(vegp$x)[i,j,]),
            mean(as.array(vegp$clump)[i,j,]),
            mean(as.array(vegp$leafr)[i,j,]),
            mean(as.array(vegp$leaft)[i,j,]),
            mean(as.array(vegp$leafd)[i,j,]),
            0.97,
            mean(as.array(vegp$gsmax)[i,j,]),
            100)
}
#' used by runpointmodela
.togp<-function(soilc,i,j) {
  sn<-.is(soilc$soiltype)[i,j]
  soiltype<-soilparameters$Soil.type[sn]
  # Create list of soil parameters
  groundp_p<-c(.mfr(soilc$groundr),0,180,0.97,
               soilparameters$rho[sn],
               soilparameters$Vm[sn],
               soilparameters$Vq[sn],
               soilparameters$Mc[sn],
               soilparameters$b[sn],
               soilparameters$psi_e[sn],
               soilparameters$Smax[sn],
               soilparameters$Smin[sn],
               soilparameters$Smin[sn],
               soilparameters$Smin[sn],
               soilparameters$Smin[sn])
  return(list(groundp_p=groundp_p,soiltype=soiltype,sn=sn))
}
# =========================================================================== #
# ************ Used to run checks and prepare model inputs  ***************** #
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
#' Converts variable in a list of data.frame to an array as used by modelina
.cca<-function(weather,varn,h,r,rfi) {
  a<-array(NA,dim=c(dim(r)[1:2],h))
  k<-1
  for (i in 1:dim(a)[1]) {
    for (j in 1:dim(a)[2]) {
      onew<-weather[[k]]
      if (class(onew) != "logical") {
        s<-which(names(onew)==varn)
        a[i,j,]<-onew[,s]
      } else a[i,j,]<-NA
      k<-k+1
    }
  }
  ro<-.rast(a,r)
  res1<-res(r)[1]
  res2<-res(rfi)[1]
  if (res1 != res2) ro<-resample(ro,rfi)
  ro<-mask(ro,rfi)
  as.array(ro)
}

#' Lape rates
.lapserate <- function(tc, ea, pk) {
  rv<-0.622*ea/(pk-ea)
  lr<-9.8076*(1+(2501000*rv)/(287*(tc+273.15)))/
    (1003.5+(0.622*2501000^2*rv)/(287*(tc+273.15)^2))
  lr
}
# ================================================================ #
# ~~~~~~~ internal model in functions for data frames and arrays
# ================================================================ #
.modelin <- function(micropoint, vegp, soilc, dtm, runchecks = TRUE) {
  # ======== Unpack variables and run checks ======================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  weather<-micropoint$weather
  if (runchecks) {
    rc<-checkinputs(weather,vegp,soilc,dtm)
    weather<-rc$weather
    precip<-rc$precip
    vegp<-rc$vegp
    soilc<-rc$soilc
  }
  r<-dtm
  micro<-list()
  micro$tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  micro$zref<-micropoint$zref
  h<-length(micro$tme)
  # ============ Turn input climate variables into arrays ============= #
  micro$tc<-.vta(weather$temp,r)
  micro$difr<-.vta(weather$difrad,r)
  di<-weather$swdown-weather$difrad
  di[di<0]<-0
  jd<-.jday(micro$tme)
  lt<-with(micro$tme,hour+min/60+sec/3600)
  sa<-.solalt(lt,micropoint$lat,micropoint$long,jd)
  ze<-90-sa
  si<-cos(ze*(pi/180))
  si[si<0]<-0
  dni<-di/si
  dni[is.na(dni)]<-0
  dni[dni>1352]<-1352
  micro$dirr<-.vta(dni,r)
  micro$lwdown<-.vta(weather$lwdown,r)
  micro$u2<-.vta(weather$windspeed,r)
  micro$pk<-.vta(weather$pres,r)
  # ============ Turn point model variables into arrays ============= #
  dfo<-micropoint$dfo
  micro$umu<-.vta(dfo$umu,r)
  micro$kp<-.vta(dfo$kp,r)
  micro$muGp<-.vta(dfo$muGp,r)
  micro$DDp<-.vta(dfo$DDp,r)
  micro$T0p<-.vta(dfo$T0p,r)
  micro$dtrp<-.vta(dfo$dtrp,r)
  micro$Gp<-.vta(dfo$G,r)
  micro$soilm<-.vta(dfo$soilm,r)
  micro$Tgp<-.vta(dfo$Tg,r)
  # ============ Get derived climate variables ============= #
  # Obtain derived variables
  estl<-.satvap(weather$temp)
  ea<-(weather$relhum/100)*estl
  tdew<-.dewpoint(ea,weather$temp)
  micro$estl<-.vta(estl,r)
  micro$ea<-.vta(ea,r)
  micro$tdew<-.vta(tdew,r)
  # ============ Convert vegetation variables to arrays  ============= #
  n<-length(micropoint$tmeorig)
  subs<-micropoint$subs
  vegl<-.sortvegp(vegp,method="R",n,subs)
  micro$pai<-vegl$pai
  micro$vegx<-vegl$x
  micro$lref<-vegl$leafr
  micro$ltra<-vegl$leaft
  micro$veghgt<-vegl$hgt
  micro$gsmax<-vegl$gsmax
  micro$clump<-vegl$clump
  micro$leafd<-vegl$leafd
  # ============ Convert soil variables to arrays  ============= #
  micro$gref<-.rta(soilc$groundr,h)
  soilp<-.sortsoilc(soilc,method="R")
  micro$rho<-soilp$rho
  micro$Vm<-soilp$Vm
  micro$Vq<-soilp$Vq
  micro$Mc<-soilp$Mc
  micro$soilb<-soilp$soilb
  micro$psi_e<-soilp$psi_e
  micro$Smax<-soilp$Smax
  micro$Smin<-soilp$Smin
  # ============ Additional variables  ======================= #
  if (class(micropoint$Tbz) != "logical") micro$Tbp<-.vta(micropoint$Tbz,r)
  micro$winddir<-weather$winddir
  micro$dtm<-dtm
  ll<-.latslonsfromr(dtm)
  micro$lats<-ll$lats
  micro$lons<-ll$lons
  micro$matemp<-micropoint$matemp
  micro$progress<-0
  class(micro) <-"microin"
  return(micro)
}
.modelina <- function(micropointa, vegp, soilc, dtm, dtmc, altcorrect = 0, runchecks = TRUE) {
  # =========== Unpack variables and run checks ======================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  weather<-list()
  dfo<-list()
  Tbz<-list()
  mat<-rep(NA,length(micropointa))
  ii<-0
  for (i in 1:length(micropointa)) {
    onepoint<-micropointa[[i]]
    if (class(onepoint) != "logical") {
      weather[[i]]<-onepoint$weather
      dfo[[i]]<-onepoint$dfo
      tme<-as.POSIXlt(weather[[i]]$obs_time,tz="UTC")
      tmeorig<-onepoint$tmeorig
      subs<-onepoint$subs
      zref=onepoint$zref
      mat[i]<-onepoint$matemp
      if (class(onepoint$Tbz) != "logical") {
        Tbz[[i]]<-data.frame(Tbz=onepoint$Tbz)
        ii<-i
      }
      if (runchecks) {
        rc<-checkinputs(weather[[i]],vegp,soilc,dtm,onepoint$zref)
        weather[[i]]<-rc$weather
        vegp<-rc$vegp
        soilc<-rc$soilc
      }
    } else {
      weather[[i]]<-NA
      dfo[[i]]<-NA
      Tbz[[i]]<-NA
    }
  }
  if (class(dtmc)[1] == "PackedSpatRaster") {
    r<-rast(dtmc)
  } else r<-dtmc
  # ================== Turn climate variables into arrays ============= #
  matemp<-mean(mat,na.rm=TRUE)
  micro<-list()
  micro$tme<-tme
  rfi<-dtm
  h<-length(tme)
  micro$tc<-.cca(weather,"temp",h,r,rfi)
  pk<-.cca(weather,"pres",h,r,r)
  relhum<-.cca(weather,"relhum",h,r,rfi)
  micro$estl<-.satvap(micro$tc)
  micro$ea<-micro$estl*relhum/100
  micro$tdew<-.dewpoint(micro$ea,micro$tc)
  # Apply altitudinal correction
  if (altcorrect == 0) {  # No altitudinal correction
    micro$pk<-as.array(resample(.rast(pk,r),rfi))
  } else { # Altitudinal correction applied
    dtmc[is.na(dtmc)]<-0
    psl<-pk/(((293-0.0065*.rta(dtmc,h))/293)^5.26)
    psl<-as.array(resample(.rast(psl,r),rfi))
    micro$pk<-psl*(((293-0.0065*.rta(rfi,h))/293)^5.26)
    dc<-resample(dtmc,rfi)
    elevd<-.rta(dc-dtm,h)
    if (altcorrect==1) {  # Fixed lapse rate
      tcdif<-elevd*(5/1000)
    } else { # Humidity-dependent lapse rate
      lr<-.lapserate(micro$tc,micro$ea,micro$pk)
      tcdif<-lr*elevd
    }
    micro$tc<-tcdif+micro$tc
  }
  micro$difr<-.cca(weather,"difrad",h,r,rfi)
  swrad<-.cca(weather,"swdown",h,r,rfi)
  di<-swrad-micro$difr
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
  micro$dirr<-di/si
  micro$dirr[is.na(micro$dirr)]<-0
  micro$dirr[micro$dirr>1352]<-1352
  micro$lwdown<-.cca(weather,"lwdown",h,r,rfi)
  # Wind speed and direction
  u2<-.cca(weather,"windspeed",h,r,r)
  wd<-.cca(weather,"winddir",h,r,r)
  wu<-u2*cos(wd*pi/180)
  wv<-u2*sin(wd*pi/180)
  wuv<-apply(wu,3,mean,na.rm=TRUE)
  wvv<-apply(wv,3,mean,na.rm=TRUE)
  wu<-as.array(resample(.rast(wu,r),rfi))
  wv<-as.array(resample(.rast(wv,r),rfi))
  micro$u2<-sqrt(wu^2+wv^2)
  micro$winddir<-(atan2(wvv,wuv)*180/pi)%%360
  # =============== Turn point model variables into arrays ============= #
  micro$umu<-.cca(dfo,"umu",h,r,rfi)
  micro$kp<-.cca(dfo,"kp",h,r,rfi)
  micro$muGp<-.cca(dfo,"muGp",h,r,rfi)
  micro$DDp<-.cca(dfo,"DDp",h,r,rfi)
  micro$T0p<-.cca(dfo,"T0p",h,r,rfi)
  micro$dtrp<-.cca(dfo,"dtrp",h,r,rfi)
  micro$Gp<-.cca(dfo,"G",h,r,rfi)
  micro$soilm<-.cca(dfo,"soilm",h,r,rfi)
  micro$Tgp<-.cca(dfo,"Tg",h,r,rfi)
  # ============ Convert vegetation variables to arrays  =============== #
  n<-length(tmeorig)
  vegl<-.sortvegp(vegp,method="R",n,subs)
  micro$pai<-vegl$pai
  micro$vegx<-vegl$x
  micro$lref<-vegl$leafr
  micro$ltra<-vegl$leaft
  micro$veghgt<-vegl$hgt
  micro$gsmax<-vegl$gsmax
  micro$clump<-vegl$clump
  micro$leafd<-vegl$leafd
  # ================ Convert soil variables to arrays  ================= #
  micro$gref<-.rta(soilc$groundr,h)
  soilp<-.sortsoilc(soilc,method="R")
  micro$rho<-soilp$rho
  micro$Vm<-soilp$Vm
  micro$Vq<-soilp$Vq
  micro$Mc<-soilp$Mc
  micro$soilb<-soilp$soilb
  micro$psi_e<-soilp$psi_e
  micro$Smax<-soilp$Smax
  micro$Smin<-soilp$Smin
  # ===================== Additional variables  ======================== #
  if (ii>0) micro$Tbp<-.cca(Tbz,"Tbz",h,r,rfi)
  micro$dtm<-dtm
  micro$zref<-zref
  ll<-.latslonsfromr(dtm)
  micro$lats<-ll$lats
  micro$lons<-ll$lons
  micro$matemp<-matemp
  micro$progress<-0
  class(micro) <-"microin"
  return(micro)
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
#' Caclulates flow accumulation (NB - probably needs modifying wiht below)
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
#' Function above needs replaces with this:
flowacc<-function (dtm, basins = NA) {
  dm <- .is(dtm)
  fd <- .flowdir(dm)
  fa <- fd * 0 + 1
  if (class(basins) != "logical")
    ba <- .is(basins)
  o <- order(dm, decreasing = T, na.last = NA)
  for (i in 1:(length(o)-1)) {
    x<-arrayInd(o[i],dim(dm))[2]
    y<-arrayInd(o[i],dim(dm))[1]
    f<-fd[y,x]
    y2<-y+(f-1)%%3-1
    x2<-x+(f-1)%/%3-1
    if (class(basins) != "logical" & x2 > 0 & y2 > 0 & x2 <=
        dim(dm)[2] & y2 <= dim(dm)[1]) {
      b1 <- ba[y, x]
      b2 <- ba[y2, x2]
      if (!is.na(b1) && !is.na(b2)) {
        if (b1 == b2 & x2 > 0 & x2 < dim(dm)[2] & y2 >
            0 & y2 < dim(dm)[1])
          fa[y2,x2]<-fa[y2,x2]+fa[y,x]
      }
    }
    else if (x2 > 0 & x2 < dim(dm)[2] & y2 > 0 & y2 < dim(dm)[1])
      fa[y2,x2]<-fa[y2,x2]+fa[y,x]
  }
  fa <- .rast(fa, dtm)
  return(fa)
}
#' Calculates topographic wetness index
.topidx <- function(dtm) {
  minslope<-atan(0.02/mean(res(dtm)))
  slope<-terrain(dtm,unit="radians")
  B<-.is(slope)
  B[B<minslope]<-minslope
  B[is.na(B)]<-median(B,na.rm=T)
  a<-flowaccCpp(.is(dtm))+1
  a<-a*res(dtm)[1]* res(dtm)[2]
  tpx<- a/tan(B)
  r<-.rast(tpx,dtm)
  r<-mask(r,dtm)
  r
}
# *** twostreammodel()
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
.foliageden<-function(z,hgt,pai,shape=1.5,rate=shape/7) {
  # reverse and rescale z
  x<-((hgt-z)/hgt)*10
  # totAL DENS
  td<-pgamma(10,shape,rate) # total density
  rfd<-dgamma(x,shape,rate)/td # relative foliage density
  tdf<-(pai/hgt)*rfd*10 # total foliage density
  paia<-pgamma(x,shape,rate)*(pai/td) # pai above
  return(list(leafden=tdf,pai_a=paia))
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
# Sorts out NAs and zeros in inputs
.cleanvars<-function(vegp,soilc,dtm) {
  .cleanr<-function(r,s,v=NA) {
    d<-dim(r)[3]
    for (i in 1:d) {
      m<-.is(r[[i]])
      m[s]<-v
      ri<-.rast(m,r[[1]])
      r[[i]]<-ri
    }
    return(r)
  }
  # set NAs
  s1<-which(is.na(.is(vegp$pai[[1]])))
  s2<-which(is.na(.is(vegp$hgt[[1]])))
  s3<-which(is.na(.is(soilc$soiltype)))
  s4<-which(is.na(.is(soilc$groundr)))
  s5<-which(is.na(.is(dtm)))
  s<-unique(c(s1,s2,s3,s4,s5))
  vegp$pai<-.cleanr(vegp$pai,s)
  vegp$hgt<-.cleanr(vegp$hgt,s)
  vegp$x<-.cleanr(vegp$x,s)
  vegp$gsmax<-.cleanr(vegp$gsmax,s)
  vegp$leafr<-.cleanr(vegp$leafr,s)
  vegp$clump<-.cleanr(vegp$clump,s)
  vegp$leafd<-.cleanr(vegp$leafd,s)
  vegp$leaft<-.cleanr(vegp$leaft,s)
  soilc$soiltype<-.cleanr(soilc$soiltype,s)
  soilc$groundr<-.cleanr(soilc$groundr,s)
  dtm<-.cleanr(dtm,s)
  # sort out zeros
  s1<-which(.is(vegp$pai[[1]])==0)
  s2<-which(.is(vegp$hgt[[1]])==0)
  s3<-which(is.na(.is(vegp$gsmax[[1]])))
  s4<-which(is.na(.is(vegp$leafr[[1]])))
  s5<-which(is.na(.is(vegp$clump[[1]])))
  s6<-which(is.na(.is(vegp$leafd[[1]])))
  s7<-which(is.na(.is(vegp$leaft[[1]])))
  s<-unique(c(s1,s2,s3,s4,s5,s6,s7))
  vegp$pai<-.cleanr(vegp$pai,s,0)
  vegp$hgt<-.cleanr(vegp$hgt,0)
  vegp$pai<-mask(vegp$pai,dtm)
  vegp$hgt<-mask(vegp$hgt,dtm)
  return(list(vegp=vegp,dtm=dtm,soilc=soilc))
}
#' Run C++ version of model with time-invarient vegetation and data.frame weather input
.runmodel1Cpp<-function(micropoint,vegp,soilc,dtm, reqhgt=0.05, runchecks = TRUE, pai_a = NA, tfact = 1.5, out = rep(TRUE,10), slr = NA,
                        apr = NA, hor = NA, twi = NA, wsa = NA, svf = NA) {
  # Unpack and check variables
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  cv<-.cleanvars(vegp,soilc,dtm) # sorts out NAs and zeros
  dtm<-cv$dtm
  vegp<-cv$vegp
  soilc<-cv$soilc
  weather<-micropoint$weather
  if (runchecks) {
    rc<-checkinputs(weather,vegp,soilc,dtm)
    weather<-rc$weather
    vegp<-rc$vegp
    soilc<-rc$soilc
  }
  # derive additional climate variables
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                      hour=tme$hour+tme$min/60+tme$sec/3600)
  weather$es<-.satvap(weather$temp)
  weather$ea<-weather$es*weather$relhum/100
  weather$tdew<-.dewpoint(weather$ea,weather$temp)
  weather$relhum<-NULL
  weather$obs_time<-NULL
  # Sort out point model variables
  pointm<-micropoint$dfo
  if (reqhgt<0) {
    pointm$Tbp<-micropoint$Tbz
  } else pointm$Tbp<-0
  # Sort out vegp
  n<-length(micropoint$tmeorig)
  subs<-micropoint$subs
  vegp<-.sortvegp(vegp,method="C",n,subs)
  s<-which(vegp$pai==0)
  vegp$hgt[s]<-0
  s<-which(vegp$hgt==0)
  vegp$pai[s]<-0
  # add additional terms to vegp
  fd<-.foliageden(reqhgt,vegp$hgt,vegp$pai)
  if (class(pai_a) == "logical") vegp$paia<-fd$pai_a
  vegp$leafden<-fd$leafden
  # sort out soilc
  soilp<-.sortsoilc(soilc,method="R")
  soilc$gref<-.is(soilc$groundr)
  soilc$groundr<-NULL
  soilc$soiltype<-NULL
  soilc$Smin<-soilp$Smin
  soilc$Smax<-soilp$Smax
  soilc$soilb<-soilp$soilb
  soilc$Psie<-soilp$psi_e
  soilc$Vq<-soilp$Vq
  soilc$Vm<-soilp$Vm
  soilc$Mc<-soilp$Mc
  soilc$rho<-soilp$rho
  # Calculate slope, aspect and topographic wetness index
  if (class(slr)[1] == "logical") {
    soilc$slope<-terrain(dtm, v="slope")
  } else soilc$slope<-slr
  if (class(apr)[1] == "logical") {
    soilc$aspect<-terrain(dtm, v="aspect")
  }  else soilc$aspect<-apr
  if (class(twi)[1] == "logical") {
    soilc$twi<-.topidx(dtm)
  } else soilc$twi<-twi
  soilc$slope[is.na(soilc$slope)]<-0
  soilc$aspect[is.na(soilc$aspect)]<-0
  soilc$slope<-.is(mask(soilc$slope,dtm))
  soilc$aspect<-.is(mask(soilc$aspect,dtm))
  soilc$twi[is.na(soilc$twi)]<-1
  soilc$twi<-.is(mask(soilc$twi,dtm))
  # Calculate sky view factor and wind shelter coefficient  etc
  sazi<-0
  ll<-.latlongfromraster(dtm)
  if (class(hor)[1] == "logical") {
    soilc$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
    for (i in 1:24) soilc$hor[,,i]<-.horizon(dtm,(i-1)*15)
  } else soilc$hor<-hor
  if (class(svf)  == "logical") {
    msl<-tan(apply(atan(soilc$hor),c(1,2),mean))
    soilc$svfa<-0.5*cos(2*msl)+0.5
  } else soilc$svfa<-.is(svf)
  s<-1
  if (res(dtm)[1]<=100) s<-10
  if (class(wsa)[1] == "logical") {
    soilc$wsa<-.windsheltera(dtm,micropoint$zref,s)
  } else soilc$wsa<-wsa
  Sminp<-.getmode(soilc$Smin)
  Smaxp<-.getmode(soilc$Smax)
  complete<-FALSE
  if (length(subs)==length(micropoint$tmeorig)) complete<-TRUE
  if (reqhgt == 0) {
    out2<-c(1,0,0,1,0,1,1,1,1,1)
    out<-out2*out
  }
  if (reqhgt < 0) {
    out2<-c(1,0,0,1,0,0,0,0,0,0)
    out<-out2*out
  }
  matemp<-micropoint$matemp
  mout<-runmicro1Cpp(obstime,weather,pointm,vegp,soilc,reqhgt,micropoint$zref,ll$lat,ll$long,Sminp,Smaxp,tfact,complete,matemp,out)
  return(mout)
}
#' Run C++ version of model with time-invarient vegetation and array weather input
.runmodel2Cpp<-function(micropointa,vegp,soilc,dtm,dtmc,reqhgt=0.05, runchecks = TRUE, altcorrect=0,
                       pai_a = NA, tfact = 1.5, out = rep(TRUE,10), slr = NA,
                       apr = NA, hor = NA, twi = NA, wsa = NA, svf = NA) {
  # =========== Unpack variables and run checks ======================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  cv<-.cleanvars(vegp,soilc,dtm) # sorts out NAs and zeros
  dtm<-cv$dtm
  vegp<-cv$vegp
  soilc<-cv$soilc
  weather<-list()
  dfo<-list()
  Tbz<-list()
  mat<-rep(NA,length(micropointa))
  ii<-0
  for (i in 1:length(micropointa)) {
    onepoint<-micropointa[[i]]
    if (class(onepoint) != "logical") {
      weather[[i]]<-onepoint$weather
      dfo[[i]]<-onepoint$dfo
      tme<-as.POSIXlt(weather[[i]]$obs_time,tz="UTC")
      tmeorig<-onepoint$tmeorig
      subs<-onepoint$subs
      zref=onepoint$zref
      mat[i]<-onepoint$matemp
      if (class(onepoint$Tbz) != "logical") {
        Tbz[[i]]<-data.frame(Tbz=onepoint$Tbz)
        ii<-i
      }
      if (runchecks) {
        rc<-checkinputs(weather[[i]],vegp,soilc,dtm,onepoint$zref)
        weather[[i]]<-rc$weather
        vegp<-rc$vegp
        soilc<-rc$soilc
      }
    } else {
      weather[[i]]<-NA
      dfo[[i]]<-NA
      Tbz[[i]]<-NA
    }
  }
  if (class(dtmc)[1] == "PackedSpatRaster") {
    r<-rast(dtmc)
  } else r<-dtmc
  # derive additional climate variables
  matemp<-mean(mat,na.rm=TRUE)
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                      hour=tme$hour+tme$min/60+tme$sec/3600)
  # ================== Turn climate variables into arrays ============= #
  climdata<-list()
  h<-length(tme)
  climdata$tc<-.cca(weather,"temp",h,dtmc,dtm)
  pk<-.cca(weather,"pres",h,r,r)
  relhum<-.cca(weather,"relhum",h,dtmc,dtm)
  climdata$es<-.satvap(climdata$tc)
  climdata$ea<-climdata$es*relhum/100
  climdata$tdew<-.dewpoint(climdata$ea,climdata$tc)
  # Apply altitudinal correction
  if (altcorrect == 0) {  # No altitudinal correction
    climdata$pk<-as.array(resample(.rast(pk,dtmc),dtm))
  } else { # Altitudinal correction applied
    dtmc[is.na(dtmc)]<-0
    psl<-pk/(((293-0.0065*.rta(dtmc,h))/293)^5.26)
    psl<-as.array(resample(.rast(psl,dtmc),dtm))
    climdata$pk<-psl*(((293-0.0065*.rta(dtm,h))/293)^5.26)
    dc<-resample(dtmc,dtm)
    elevd<-.rta(dc-dtm,h)
    if (altcorrect==1) {  # Fixed lapse rate
      tcdif<-elevd*(5/1000)
    } else { # Humidity-dependent lapse rate
      lr<-.lapserate(climdata$tc,climdata$ea,climdata$pk)
      tcdif<-lr*elevd
    }
    climdata$tc<-tcdif+climdata$tc
  }
  climdata$difrad<-.cca(weather,"difrad",h,dtmc,dtm)
  climdata$swdown<-.cca(weather,"swdown",h,dtmc,dtm)
  climdata$lwdown<-.cca(weather,"lwdown",h,dtmc,dtm)
  # Wind speed and direction
  u2<-.cca(weather,"windspeed",h,dtmc,dtmc)
  wd<-.cca(weather,"winddir",h,dtmc,dtmc)
  wu<-u2*cos(wd*pi/180)
  wv<-u2*sin(wd*pi/180)
  wuv<-apply(wu,3,mean,na.rm=TRUE)
  wvv<-apply(wv,3,mean,na.rm=TRUE)
  wu<-as.array(resample(.rast(wu,dtmc),dtm))
  wv<-as.array(resample(.rast(wv,dtmc),dtm))
  climdata$windspeed<-sqrt(wu^2+wv^2)
  climdata$winddir<-(atan2(wvv,wuv)*180/pi)%%360
  # Point model variables
  # =============== Turn point model variables into arrays ============= #
  pointm<-list()
  pointm$umu<-.cca(dfo,"umu",h,dtmc,dtm)
  pointm$kp<-.cca(dfo,"kp",h,dtmc,dtm)
  pointm$muGp<-.cca(dfo,"muGp",h,dtmc,dtm)
  pointm$DDp<-.cca(dfo,"DDp",h,dtmc,dtm)
  pointm$T0p<-.cca(dfo,"T0p",h,dtmc,dtm)
  pointm$dtrp<-.cca(dfo,"dtrp",h,dtmc,dtm)
  pointm$Gp<-.cca(dfo,"G",h,dtmc,dtm)
  pointm$soilm<-.cca(dfo,"soilm",h,dtmc,dtm)
  pointm$Tg<-.cca(dfo,"Tg",h,dtmc,dtm)
  if (reqhgt<0)  {
    pointm$Tbp<-.cca(Tbz,"Tbz",h,dtmc,dtm)
  } else  pointm$Tbp<-pointm$soilm*0
  # ===================== Sort out vegp ================================== #
  n<-length(tmeorig)
  vegp<-.sortvegp(vegp,method="C",n,subs)
  # add additional terms to vegp
  fd<-.foliageden(reqhgt,vegp$hgt,vegp$pai)
  if (class(pai_a) == "logical") vegp$paia<-fd$pai_a
  vegp$leafden<-fd$leafden
  # ===================== Sort out soilc ================================== #
  soilp<-.sortsoilc(soilc,method="R")
  soilc$gref<-.is(soilc$groundr)
  soilc$groundr<-NULL
  soilc$soiltype<-NULL
  soilc$Smin<-soilp$Smin
  soilc$Smax<-soilp$Smax
  soilc$soilb<-soilp$soilb
  soilc$Psie<-soilp$psi_e
  soilc$Vq<-soilp$Vq
  soilc$Vm<-soilp$Vm
  soilc$Mc<-soilp$Mc
  soilc$rho<-soilp$rho
  # Calculate slope, aspect and topographic wetness index
  if (class(slr)[1] == "logical") {
    soilc$slope<-terrain(dtm, v="slope")
  } else soilc$slope<-slr
  if (class(apr)[1] == "logical") {
    soilc$aspect<-terrain(dtm, v="aspect")
  }  else soilc$aspect<-apr
  if (class(twi)[1] == "logical") {
    soilc$twi<-.topidx(dtm)
  } else soilc$twi<-twi
  soilc$slope[is.na(soilc$slope)]<-0
  soilc$aspect[is.na(soilc$aspect)]<-0
  soilc$slope<-.is(mask(soilc$slope,dtm))
  soilc$aspect<-.is(mask(soilc$aspect,dtm))
  soilc$twi[is.na(soilc$twi)]<-1
  soilc$twi<-.is(mask(soilc$twi,dtm))
  # Calculate sky view factor and wind shelter coefficient  etc
  ll<-.latslonsfromr(dtm)
  if (class(hor)[1] == "logical") {
    soilc$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
    for (i in 1:24) soilc$hor[,,i]<-.horizon(dtm,(i-1)*15)
  } else soilc$hor<-hor
  if (class(svf)  == "logical") {
    msl<-tan(apply(atan(soilc$hor),c(1,2),mean))
    soilc$svfa<-0.5*cos(2*msl)+0.5
  } else soilc$svfa<-.is(svf)
  s<-1
  if (res(dtm)[1]<=100) s<-10
  if (class(wsa)[1] == "logical") {
    soilc$wsa<-.windsheltera(dtm,zref,s)
  } else soilc$wsa<-wsa
  Sminp<-.getmode(soilc$Smin)
  Smaxp<-.getmode(soilc$Smax)
  complete<-FALSE
  if (length(subs)==length(tmeorig)) complete<-TRUE
  if (reqhgt == 0) {
    out2<-c(1,0,0,1,0,1,1,1,1,1)
    out<-out2*out
  }
  if (reqhgt < 0) {
    out2<-c(1,0,0,1,0,0,0,0,0,0)
    out<-out2*out
  }
  mout<-runmicro2Cpp(obstime,climdata,pointm,vegp,soilc,reqhgt,zref,ll$lats,ll$lons,Sminp,Smaxp,tfact,complete,matemp,out)
  return(mout)
}
.runmodel3Cpp<-function(micropoint,vegp,soilc,dtm, reqhgt=0.05, runchecks = TRUE, pai_a = NA, tfact = 1.5, out = rep(TRUE,10), slr = NA,
                        apr = NA, hor = NA, twi = NA, wsa = NA, svf = NA) {
  # Unpack and check variables
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  cv<-.cleanvars(vegp,soilc,dtm) # sorts out NAs and zeros
  dtm<-cv$dtm
  vegp<-cv$vegp
  soilc<-cv$soilc
  weather<-micropoint$weather
  if (runchecks) {
    rc<-checkinputs(weather,vegp,soilc,dtm)
    weather<-rc$weather
    precip<-rc$precip
    vegp<-rc$vegp
    soilc<-rc$soilc
  }
  # derive additional climate variables
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                      hour=tme$hour+tme$min/60+tme$sec/3600)
  weather$es<-.satvap(weather$temp)
  weather$ea<-weather$es*weather$relhum/100
  weather$tdew<-.dewpoint(weather$ea,weather$temp)
  weather$relhum<-NULL
  weather$obs_time<-NULL
  # Sort out point model variables
  pointm<-micropoint$dfo
  if (reqhgt<0) {
    pointm$Tbp<-micropoint$Tbz
  } else pointm$Tbp<-0
  # Sort out vegp
  n<-length(micropoint$tmeorig)
  subs<-micropoint$subs
  vegp<-.sortvegp(vegp,method="C",n,subs)
  s<-which(vegp$pai==0)
  vegp$hgt[s]<-0
  s<-which(vegp$hgt==0)
  vegp$pai[s]<-0
  # add additional terms to vegp
  fd<-.foliageden(reqhgt,vegp$hgt,vegp$pai)
  if (class(pai_a) == "logical") vegp$paia<-fd$pai_a
  vegp$leafden<-fd$leafden
  # Create dfsel
  ss<-vegp$lsubs
  dflyr<-data.frame(lyr=unique(ss),st=0,ed=0)
  for (i in 1:length(dflyr$lyr)) {
    s<-which(ss==dflyr$lyr[i])
    dflyr$st[i]<-floor(s[1]/24)*24
    dflyr$ed[i]<-floor(s[length(s)]/24)*24-1
  }
  dflyr$lyr<-c(1:length(dflyr$lyr))
  # sort out soilc
  soilp<-.sortsoilc(soilc,method="R")
  soilc$gref<-.is(soilc$groundr)
  soilc$groundr<-NULL
  soilc$soiltype<-NULL
  soilc$Smin<-soilp$Smin
  soilc$Smax<-soilp$Smax
  soilc$soilb<-soilp$soilb
  soilc$Psie<-soilp$psi_e
  soilc$Vq<-soilp$Vq
  soilc$Vm<-soilp$Vm
  soilc$Mc<-soilp$Mc
  soilc$rho<-soilp$rho
  # Calculate slope, aspect and topographic wetness index
  if (class(slr)[1] == "logical") {
    soilc$slope<-terrain(dtm, v="slope")
  } else soilc$slope<-slr
  if (class(apr)[1] == "logical") {
    soilc$aspect<-terrain(dtm, v="aspect")
  }  else soilc$aspect<-apr
  if (class(twi)[1] == "logical") {
    soilc$twi<-.topidx(dtm)
  } else soilc$twi<-twi
  soilc$slope[is.na(soilc$slope)]<-0
  soilc$aspect[is.na(soilc$aspect)]<-0
  soilc$slope<-.is(mask(soilc$slope,dtm))
  soilc$aspect<-.is(mask(soilc$aspect,dtm))
  soilc$twi[is.na(soilc$twi)]<-1
  soilc$twi<-.is(mask(soilc$twi,dtm))
  # Calculate sky view factor and wind shelter coefficient  etc
  sazi<-0
  ll<-.latlongfromraster(dtm)
  if (class(hor)[1] == "logical") {
    soilc$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
    for (i in 1:24) soilc$hor[,,i]<-.horizon(dtm,(i-1)*15)
  } else soilc$hor<-hor
  if (class(svf)  == "logical") {
    msl<-tan(apply(atan(soilc$hor),c(1,2),mean))
    soilc$svfa<-0.5*cos(2*msl)+0.5
  } else soilc$svfa<-.is(svf)
  s<-1
  if (res(dtm)[1]<=100) s<-10
  if (class(wsa)[1] == "logical") {
    soilc$wsa<-.windsheltera(dtm,micropoint$zref,s)
  } else soilc$wsa<-wsa
  Sminp<-.getmode(soilc$Smin)
  Smaxp<-.getmode(soilc$Smax)
  complete<-FALSE
  if (length(subs)==length(micropoint$tmeorig)) complete<-TRUE
  if (reqhgt == 0) {
    out2<-c(1,0,0,1,0,1,1,1,1,1)
    out<-out2*out
  }
  if (reqhgt < 0) {
    out2<-c(1,0,0,1,0,0,0,0,0,0)
    out<-out2*out
  }
  matemp<-micropoint$matemp
  mout<-runmicro3Cpp(dflyr,obstime,weather,pointm,vegp,soilc,reqhgt,micropoint$zref,ll$lat,ll$long,Sminp,Smaxp,tfact,complete,matemp,out)
  return(mout)
}
.runmodel4Cpp<-function(micropointa,vegp,soilc,dtm,dtmc,reqhgt=0.05, runchecks = TRUE, altcorrect=0,
                        pai_a = NA, tfact = 1.5, out = rep(TRUE,10), slr = NA,
                        apr = NA, hor = NA, twi = NA, wsa = NA, svf = NA) {
  # =========== Unpack variables and run checks ======================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  cv<-.cleanvars(vegp,soilc,dtm) # sorts out NAs and zeros
  dtm<-cv$dtm
  vegp<-cv$vegp
  soilc<-cv$soilc
  weather<-list()
  dfo<-list()
  Tbz<-list()
  mat<-rep(NA,length(micropointa))
  ii<-0
  for (i in 1:length(micropointa)) {
    onepoint<-micropointa[[i]]
    if (class(onepoint) != "logical") {
      weather[[i]]<-onepoint$weather
      dfo[[i]]<-onepoint$dfo
      tme<-as.POSIXlt(weather[[i]]$obs_time,tz="UTC")
      tmeorig<-onepoint$tmeorig
      subs<-onepoint$subs
      zref=onepoint$zref
      mat[i]<-onepoint$matemp
      if (class(onepoint$Tbz) != "logical") {
        Tbz[[i]]<-data.frame(Tbz=onepoint$Tbz)
        ii<-i
      }
      if (runchecks) {
        rc<-checkinputs(weather[[i]],vegp,soilc,dtm,onepoint$zref)
        weather[[i]]<-rc$weather
        vegp<-rc$vegp
        soilc<-rc$soilc
      }
    } else {
      weather[[i]]<-NA
      dfo[[i]]<-NA
      Tbz[[i]]<-NA
    }
  }
  if (class(dtmc)[1] == "PackedSpatRaster") {
    r<-rast(dtmc)
  } else r<-dtmc
  matemp<-mean(mat,na.rm=TRUE)
  # derive additional climate variables
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                      hour=tme$hour+tme$min/60+tme$sec/3600)
  # ================== Turn climate variables into arrays ============= #
  climdata<-list()
  h<-length(tme)
  climdata$tc<-.cca(weather,"temp",h,dtmc,dtm)
  pk<-.cca(weather,"pres",h,r,r)
  relhum<-.cca(weather,"relhum",h,dtmc,dtm)
  climdata$es<-.satvap(climdata$tc)
  climdata$ea<-climdata$es*relhum/100
  climdata$tdew<-.dewpoint(climdata$ea,climdata$tc)
  # Apply altitudinal correction
  if (altcorrect == 0) {  # No altitudinal correction
    climdata$pk<-as.array(resample(.rast(pk,dtmc),dtm))
  } else { # Altitudinal correction applied
    dtmc[is.na(dtmc)]<-0
    psl<-pk/(((293-0.0065*.rta(dtmc,h))/293)^5.26)
    psl<-as.array(resample(.rast(psl,dtmc),dtm))
    climdata$pk<-psl*(((293-0.0065*.rta(dtm,h))/293)^5.26)
    dc<-resample(dtmc,dtm)
    elevd<-.rta(dc-dtm,h)
    if (altcorrect==1) {  # Fixed lapse rate
      tcdif<-elevd*(5/1000)
    } else { # Humidity-dependent lapse rate
      lr<-.lapserate(climdata$tc,climdata$ea,climdata$pk)
      tcdif<-lr*elevd
    }
    climdata$tc<-tcdif+climdata$tc
  }
  climdata$difrad<-.cca(weather,"difrad",h,dtmc,dtm)
  climdata$swdown<-.cca(weather,"swdown",h,dtmc,dtm)
  climdata$lwdown<-.cca(weather,"lwdown",h,dtmc,dtm)
  # Wind speed and direction
  u2<-.cca(weather,"windspeed",h,dtmc,dtmc)
  wd<-.cca(weather,"winddir",h,dtmc,dtmc)
  wu<-u2*cos(wd*pi/180)
  wv<-u2*sin(wd*pi/180)
  wuv<-apply(wu,3,mean,na.rm=TRUE)
  wvv<-apply(wv,3,mean,na.rm=TRUE)
  wu<-as.array(resample(.rast(wu,dtmc),dtm))
  wv<-as.array(resample(.rast(wv,dtmc),dtm))
  climdata$windspeed<-sqrt(wu^2+wv^2)
  climdata$winddir<-(atan2(wvv,wuv)*180/pi)%%360
  # Point model variables
  # =============== Turn point model variables into arrays ============= #
  pointm<-list()
  pointm$umu<-.cca(dfo,"umu",h,dtmc,dtm)
  pointm$kp<-.cca(dfo,"kp",h,dtmc,dtm)
  pointm$muGp<-.cca(dfo,"muGp",h,dtmc,dtm)
  pointm$DDp<-.cca(dfo,"DDp",h,dtmc,dtm)
  pointm$T0p<-.cca(dfo,"T0p",h,dtmc,dtm)
  pointm$dtrp<-.cca(dfo,"dtrp",h,dtmc,dtm)
  pointm$Gp<-.cca(dfo,"G",h,dtmc,dtm)
  pointm$soilm<-.cca(dfo,"soilm",h,dtmc,dtm)
  pointm$Tg<-.cca(dfo,"Tg",h,dtmc,dtm)
  if (reqhgt<0)  {
    pointm$Tbp<-.cca(Tbz,"Tbz",h,dtmc,dtm)
  } else  pointm$Tbp<-pointm$soilm*0
  # ===================== Sort out vegp ================================== #
  n<-length(tmeorig)
  vegp<-.sortvegp(vegp,method="C",n,subs)
  # add additional terms to vegp
  fd<-.foliageden(reqhgt,vegp$hgt,vegp$pai)
  if (class(pai_a) == "logical") vegp$paia<-fd$pai_a
  vegp$leafden<-fd$leafden

  # ===================== Create dfsel ================================== #
  ss<-vegp$lsubs
  dflyr<-data.frame(lyr=unique(ss),st=0,ed=0)
  for (i in 1:length(dflyr$lyr)) {
    s<-which(ss==dflyr$lyr[i])
    dflyr$st[i]<-floor(s[1]/24)*24
    dflyr$ed[i]<-floor(s[length(s)]/24)*24-1
  }
  dflyr$lyr<-c(1:length(dflyr$lyr))
  # ===================== Sort out soilc ================================== #
  soilp<-.sortsoilc(soilc,method="R")
  soilc$gref<-.is(soilc$groundr)
  soilc$groundr<-NULL
  soilc$soiltype<-NULL
  soilc$Smin<-soilp$Smin
  soilc$Smax<-soilp$Smax
  soilc$soilb<-soilp$soilb
  soilc$Psie<-soilp$psi_e
  soilc$Vq<-soilp$Vq
  soilc$Vm<-soilp$Vm
  soilc$Mc<-soilp$Mc
  soilc$rho<-soilp$rho
  # Calculate slope, aspect and topographic wetness index
  if (class(slr)[1] == "logical") {
    soilc$slope<-terrain(dtm, v="slope")
  } else soilc$slope<-slr
  if (class(apr)[1] == "logical") {
    soilc$aspect<-terrain(dtm, v="aspect")
  }  else soilc$aspect<-apr
  if (class(twi)[1] == "logical") {
    soilc$twi<-.topidx(dtm)
  } else soilc$twi<-twi
  soilc$slope[is.na(soilc$slope)]<-0
  soilc$aspect[is.na(soilc$aspect)]<-0
  soilc$slope<-.is(mask(soilc$slope,dtm))
  soilc$aspect<-.is(mask(soilc$aspect,dtm))
  soilc$twi[is.na(soilc$twi)]<-1
  soilc$twi<-.is(mask(soilc$twi,dtm))
  # Calculate sky view factor and wind shelter coefficient  etc
  ll<-.latslonsfromr(dtm)
  if (class(hor)[1] == "logical") {
    soilc$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
    for (i in 1:24) soilc$hor[,,i]<-.horizon(dtm,(i-1)*15)
  } else soilc$hor<-hor
  if (class(svf)  == "logical") {
    msl<-tan(apply(atan(soilc$hor),c(1,2),mean))
    soilc$svfa<-0.5*cos(2*msl)+0.5
  } else soilc$svfa<-.is(svf)
  s<-1
  if (res(dtm)[1]<=100) s<-10
  if (class(wsa)[1] == "logical") {
    soilc$wsa<-.windsheltera(dtm,zref,s)
  } else soilc$wsa<-wsa
  Sminp<-.getmode(soilc$Smin)
  Smaxp<-.getmode(soilc$Smax)
  complete<-FALSE
  if (length(subs)==length(tmeorig)) complete<-TRUE
  if (reqhgt == 0) {
    out2<-c(1,0,0,1,0,1,1,1,1,1)
    out<-out2*out
  }
  if (reqhgt < 0) {
    out2<-c(1,0,0,1,0,0,0,0,0,0)
    out<-out2*out
  }
  mout<-runmicro4Cpp(dflyr,obstime,climdata,pointm,vegp,soilc,reqhgt,zref,ll$lats,ll$lons,Sminp,Smaxp,tfact,complete,matemp,out)
  return(mout)
}
#' Check SpatRasts of model input
.checkbiginputs<-function(dtm,vegp,soilc) {
  .cr<-function(r,dtm) {
    sel<-which(is.na(.is(dtm))==FALSE & is.na(.is(r[[1]])))
    length(sel)
  }
  le<-c(.cr(vegp$pai,dtm),.cr(vegp$hgt,dtm),.cr(vegp$x,dtm),.cr(vegp$gsmax,dtm),
        .cr(vegp$leafr,dtm),.cr(vegp$clump,dtm),.cr(vegp$leafd,dtm),.cr(vegp$leaft,dtm))
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
#' Crop raster for running runmicro_big
.croprast<-function(r,rw,cl,tilesize,toverlap) {
  e<-ext(r)
  xmn<-e$xmin+(cl-1)*tilesize*res(r)[1]-toverlap
  xmx<-e$xmin+cl*tilesize*res(r)[1]+toverlap
  ymn<-e$ymax-rw*tilesize*res(r)[2]-toverlap
  ymx<-e$ymax-(rw-1)*tilesize*res(r)[2]+toverlap
  xmn<-ifelse(xmn<e$xmin,e$xmin,xmn)
  xmx<-ifelse(xmx>e$xmax,e$xmax,xmx)
  ymn<-ifelse(ymn<e$ymin,e$ymin,ymn)
  ymx<-ifelse(ymx>e$ymax,e$ymax,ymx)
  e2<-ext(c(xmn,xmx,ymn,ymx))
  r2<-terra::crop(r,e2)
  r2
}
#' Crop list of rasters for running runmicro_big
.listcrop<-function(rlst,e) {
  rlsto<-list()
  for (ii in 1:length(rlst)) {
    rlsto[[ii]]<-crop(rlst[[ii]],e)
  }
  names(rlsto)<-names(rlst)
  return(rlsto)
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
.biomicropoint<-function(weather,tme,reqhgt,vegp,soilc,dtm,zref,windhgt,soilm) {
  if (class(weather) == "data.frame") {
    # Subset right values
    tc<-weather$temp
    sel<-.biosel(tme, tc)
    weather2<-weather[sel$selh,]
    micropoint<-runpointmodel(weather2,reqhgt,dtm,vegp,soilc,FALSE,zref,windhgt,soilm,NA,yearG=FALSE)
  } else {
    # Subset right values
    tc<-apply(.is(weather$temp),3,mean,na.rm=TRUE)
    sel<-.biosel(tme, tc)
    r<-weather$temp[[1]]
    climarray2<-list()
    for (i in 1:length(weather)) {
      a<-.is(weather[[i]])
      a<-a[,,sel$selh]
      climarray2[[i]]<-.rast(a,r)
    }
    names(climarray2)<-names(weather)
    micropoint<-runpointmodela(climarray2,tme[sel$selh],reqhgt,dtm,vegp,soilc,FALSE,zref,windhgt,soilm,NA,yearG=FALSE)
  }
  return(list(micropoint=micropoint,sel=sel))
}
#' Create Vector of quarter entries
.getselq<-function(iq,tme) {
  imn<-iq-1
  imx<-iq+1
  if (imn == 0) imn = 12
  if (imx == 13) imx = 1
  s1<-which(tme$mon+1==imn)
  s2<-which(tme$mon+1==iq)
  s3<-which(tme$mon+1==imx)
  sel<-c(s1,s2,s3)
  sel<-sel[order(sel)]
}
#' Run bioclim with data.frame weather and static veg
.runbioclim1<-function(weather, reqhgt = 0.05, vegp, soilc, dtm,
                       temp = "air", zref = 2, windhgt = zref, soilm = NA, runchecks = TRUE,
                       pai_a = NA, tfact = 1.5, out = rep(TRUE, 19)) {
  # ====================== Unpack and check variables ======================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  cv<-.cleanvars(vegp,soilc,dtm) # sorts out NAs and zeros
  dtm<-cv$dtm
  vegp<-cv$vegp
  soilc<-cv$soilc
  if (runchecks) {
    rc<-checkinputs(weather,vegp,soilc,dtm)
    weather<-rc$weather
    vegp<-rc$vegp
    soilc<-rc$soilc
  }
  # ========================== Calculate quarters ============================ #
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  agg<-stats::aggregate(weather$precip,by=list(tme$mon),mean,na.rm=TRUE)$x
  wq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
  dq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # driest quarter
  agg<-stats::aggregate(weather$temp,by=list(tme$mon),sum,na.rm=TRUE)$x
  hq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # hottest quarter
  cq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
  # =================== Subset time-series and run point model =============== #
  micropoint<-.biomicropoint(weather,tme,reqhgt,vegp,soilc,dtm,zref,windhgt,soilm)$micropoint
  matemp<-micropoint$matemp
  weather<-micropoint$weather
  # ==================== derive additional climate variables ================= #
  # derive additional climate variables
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                      hour=tme$hour+tme$min/60+tme$sec/3600)
  weather$es<-.satvap(weather$temp)
  weather$ea<-weather$es*weather$relhum/100
  weather$tdew<-.dewpoint(weather$ea,weather$temp)
  weather$relhum<-NULL
  weather$obs_time<-NULL
  # ======================== Sort out point model variables ================== #
  pointm<-micropoint$dfo
  if (reqhgt<0) {
    pointm$Tbp<-micropoint$Tbz
  } else pointm$Tbp<-0
  # ===================== Sort out vegp and add additional terms ============= #
  n<-length(micropoint$tmeorig)
  subs<-micropoint$subs
  vegp<-.sortvegp(vegp,method="C",n,subs)
  s<-which(vegp$pai==0)
  vegp$hgt[s]<-0
  s<-which(vegp$hgt==0)
  vegp$pai[s]<-0
  # add additional terms to vegp
  fd<-.foliageden(reqhgt,vegp$hgt,vegp$pai)
  if (class(pai_a) == "logical") vegp$paia<-fd$pai_a
  vegp$leafden<-fd$leafden
  # ========================== sort out soilc ================================ #
  soilp<-.sortsoilc(soilc,method="R")
  soilc$gref<-.is(soilc$groundr)
  soilc$groundr<-NULL
  soilc$soiltype<-NULL
  soilc$Smin<-soilp$Smin
  soilc$Smax<-soilp$Smax
  soilc$soilb<-soilp$soilb
  soilc$Psie<-soilp$psi_e
  soilc$Vq<-soilp$Vq
  soilc$Vm<-soilp$Vm
  soilc$Mc<-soilp$Mc
  soilc$rho<-soilp$rho
  # ============== Calculate slope, aspect and topographic wetness index ===== #
  soilc$slope<-terrain(dtm, v="slope")
  soilc$aspect<-terrain(dtm, v="aspect")
  soilc$twi<-.topidx(dtm)
  soilc$slope[is.na(soilc$slope)]<-0
  soilc$aspect[is.na(soilc$aspect)]<-0
  soilc$slope<-.is(mask(soilc$slope,dtm))
  soilc$aspect<-.is(mask(soilc$aspect,dtm))
  soilc$twi[is.na(soilc$twi)]<-1
  soilc$twi<-.is(mask(soilc$twi,dtm))
  # ======== Calculate sky view factor and wind shelter coefficient  etc ===== #
  sazi<-0
  ll<-.latlongfromraster(dtm)
  soilc$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) soilc$hor[,,i]<-.horizon(dtm,(i-1)*15)
  msl<-tan(apply(atan(soilc$hor),c(1,2),mean))
  soilc$svfa<-0.5*cos(2*msl)+0.5
  s<-1
  if (res(dtm)[1]<=100) s<-10
  soilc$wsa<-.windsheltera(dtm,micropoint$zref,s)
  Sminp<-.getmode(soilc$Smin)
  Smaxp<-.getmode(soilc$Smax)
  # ========================== Calculate selqs ============================ #
  wetq<-.getselq(wq,tme)-1
  dryq<-.getselq(dq,tme)-1
  hotq<-.getselq(hq,tme)-1
  colq<-.getselq(cq,tme)-1
  #' Run C++ version of model with time-invarient vegetation and data.frame weather input
  air<-FALSE
  if (temp == "air") air<-TRUE
  bioclim<-runbioclim1Cpp(obstime,weather,pointm,vegp,soilc,reqhgt,micropoint$zref,ll$lat,ll$long,
                          Sminp,Smaxp,tfact,matemp,out,wetq,dryq,hotq,colq,air)
  bior<-dtm
  index<-1
  for (i in 1:19) {
    if (out[i]) {
      bior<-c(bior,.rast(bioclim[[index]],dtm))
      index<-index+1
    }
  }
  bior<-bior[[-1]]
  bior<-mask(bior,dtm)
  nms<-paste0("bio",c(1:19))
  s<-which(out)
  nms<-nms[s]
  names(bior)<-nms
  return(bior)
}
# Array climate input, static veg
.runbioclim2<-function(climarray,tme,reqhgt=0.05,vegp,soilc,dtm,dtmc,temp="air",
                       zref = 2, windhgt = zref, soilm = NA, runchecks = TRUE, altcorrect=0,
                       pai_a = NA, tfact = 1.5, out = rep(TRUE,19)) {
  # ============================== Unpack variables  ========================= #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  cv<-.cleanvars(vegp,soilc,dtm) # sorts out NAs and zeros
  dtm<-cv$dtm
  vegp<-cv$vegp
  soilc<-cv$soilc
  # ========================== Calculate quarters ============================ #
  prech<-apply(.is(climarray$precip),3,mean,na.rm=TRUE)
  tc<-apply(.is(climarray$temp),3,mean,na.rm=TRUE)
  agg<-stats::aggregate(prech,by=list(tme$mon),mean,na.rm=TRUE)$x
  wq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
  dq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # driest quarter
  agg<-stats::aggregate(tc,by=list(tme$mon),sum,na.rm=TRUE)$x
  hq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # hottest quarter
  cq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
  # =============== subset climate array and run point models  =============== #
  micropointa<-.biomicropoint(climarray,tme,reqhgt,vegp,soilc,dtm,zref,windhgt,soilm)$micropoint
  weather<-list()
  # ==================== extract data from micropointa ======================= #
  dfo<-list()
  Tbz<-list()
  mat<-rep(NA,length(micropointa))
  ii<-0
  for (i in 1:length(micropointa)) {
    onepoint<-micropointa[[i]]
    if (class(onepoint) != "logical") {
      weather[[i]]<-onepoint$weather
      dfo[[i]]<-onepoint$dfo
      tme<-as.POSIXlt(weather[[i]]$obs_time,tz="UTC")
      tmeorig<-onepoint$tmeorig
      subs<-onepoint$subs
      zre<-onepoint$zref
      mat[i]<-onepoint$matemp
      if (class(onepoint$Tbz) != "logical") {
        Tbz[[i]]<-data.frame(Tbz=onepoint$Tbz)
        ii<-i
      }
      if (runchecks) {
        rc<-checkinputs(weather[[i]],vegp,soilc,dtm,onepoint$zref)
        weather[[i]]<-rc$weather
        vegp<-rc$vegp
        soilc<-rc$soilc
      }
    } else {
      weather[[i]]<-NA
      dfo[[i]]<-NA
      Tbz[[i]]<-NA
    }
  }
  if (class(dtmc)[1] == "PackedSpatRaster") {
    r<-rast(dtmc)
  } else r<-dtmc
  matemp<-mean(mat,na.rm=TRUE)
  # derive additional climate variables
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                      hour=tme$hour+tme$min/60+tme$sec/3600)
  # ================== Turn climate variables into arrays ============= #
  climdata<-list()
  h<-length(tme)
  climdata$tc<-.cca(weather,"temp",h,dtmc,dtm)
  pk<-.cca(weather,"pres",h,dtmc,dtmc)
  relhum<-.cca(weather,"relhum",h,dtmc,dtm)
  climdata$es<-.satvap(climdata$tc)
  climdata$ea<-climdata$es*relhum/100
  climdata$tdew<-.dewpoint(climdata$ea,climdata$tc)
  # Apply altitudinal correction
  if (altcorrect == 0) {  # No altitudinal correction
    climdata$pk<-as.array(resample(.rast(pk,dtmc),dtm))
  } else { # Altitudinal correction applied
    dtmc[is.na(dtmc)]<-0
    psl<-pk/(((293-0.0065*.rta(dtmc,h))/293)^5.26)
    psl<-as.array(resample(.rast(psl,dtmc),dtm))
    climdata$pk<-psl*(((293-0.0065*.rta(dtm,h))/293)^5.26)
    dc<-resample(dtmc,dtm)
    elevd<-.rta(dc-dtm,h)
    if (altcorrect==1) {  # Fixed lapse rate
      tcdif<-elevd*(5/1000)
    } else { # Humidity-dependent lapse rate
      lr<-.lapserate(climdata$tc,climdata$ea,climdata$pk)
      tcdif<-lr*elevd
    }
    climdata$tc<-tcdif+climdata$tc
  }
  climdata$difrad<-.cca(weather,"difrad",h,dtmc,dtm)
  climdata$swdown<-.cca(weather,"swdown",h,dtmc,dtm)
  climdata$lwdown<-.cca(weather,"lwdown",h,dtmc,dtm)
  # Wind speed and direction
  u2<-.cca(weather,"windspeed",h,dtmc,dtmc)
  wd<-.cca(weather,"winddir",h,dtmc,dtmc)
  wu<-u2*cos(wd*pi/180)
  wv<-u2*sin(wd*pi/180)
  wuv<-apply(wu,3,mean,na.rm=TRUE)
  wvv<-apply(wv,3,mean,na.rm=TRUE)
  wu<-as.array(resample(.rast(wu,dtmc),dtm))
  wv<-as.array(resample(.rast(wv,dtmc),dtm))
  climdata$windspeed<-sqrt(wu^2+wv^2)
  climdata$winddir<-(atan2(wvv,wuv)*180/pi)%%360
  # Point model variables
  # =============== Turn point model variables into arrays ============= #
  pointm<-list()
  pointm$umu<-.cca(dfo,"umu",h,dtmc,dtm)
  pointm$kp<-.cca(dfo,"kp",h,dtmc,dtm)
  pointm$muGp<-.cca(dfo,"muGp",h,dtmc,dtm)
  pointm$DDp<-.cca(dfo,"DDp",h,dtmc,dtm)
  pointm$T0p<-.cca(dfo,"T0p",h,dtmc,dtm)
  pointm$dtrp<-.cca(dfo,"dtrp",h,dtmc,dtm)
  pointm$Gp<-.cca(dfo,"G",h,dtmc,dtm)
  pointm$soilm<-.cca(dfo,"soilm",h,dtmc,dtm)
  pointm$Tg<-.cca(dfo,"Tg",h,dtmc,dtm)
  if (reqhgt<0)  {
    pointm$Tbp<-.cca(Tbz,"Tbz",h,dtmc,dtm)
  } else  pointm$Tbp<-pointm$soilm*0
  # ===================== Sort out vegp ================================== #
  n<-length(tmeorig)
  vegp<-.sortvegp(vegp,method="C",n,subs)
  # add additional terms to vegp
  fd<-.foliageden(reqhgt,vegp$hgt,vegp$pai)
  if (class(pai_a) == "logical") vegp$paia<-fd$pai_a
  vegp$leafden<-fd$leafden
  # ======================== Sort out soilc ================================== #
  soilp<-.sortsoilc(soilc,method="R")
  soilc$gref<-.is(soilc$groundr)
  soilc$groundr<-NULL
  soilc$soiltype<-NULL
  soilc$Smin<-soilp$Smin
  soilc$Smax<-soilp$Smax
  soilc$soilb<-soilp$soilb
  soilc$Psie<-soilp$psi_e
  soilc$Vq<-soilp$Vq
  soilc$Vm<-soilp$Vm
  soilc$Mc<-soilp$Mc
  soilc$rho<-soilp$rho
  # ============= Calculate slope, aspect and topographic wetness index ====== #
  soilc$slope<-terrain(dtm, v="slope")
  soilc$aspect<-terrain(dtm, v="aspect")
  soilc$twi<-.topidx(dtm)
  soilc$slope[is.na(soilc$slope)]<-0
  soilc$aspect[is.na(soilc$aspect)]<-0
  soilc$slope<-.is(mask(soilc$slope,dtm))
  soilc$aspect<-.is(mask(soilc$aspect,dtm))
  soilc$twi[is.na(soilc$twi)]<-1
  soilc$twi<-.is(mask(soilc$twi,dtm))
  # ====== Calculate sky view factor and wind shelter coefficient etc ======== #
  ll<-.latslonsfromr(dtm)
  soilc$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) soilc$hor[,,i]<-.horizon(dtm,(i-1)*15)
  msl<-tan(apply(atan(soilc$hor),c(1,2),mean))
  soilc$svfa<-0.5*cos(2*msl)+0.5
  s<-1
  if (res(dtm)[1]<=100) s<-10
  soilc$wsa<-.windsheltera(dtm,zref,s)
  Sminp<-.getmode(soilc$Smin)
  Smaxp<-.getmode(soilc$Smax)
  # ========================== Calculate selqs ============================ #
  wetq<-.getselq(wq,tme)-1
  dryq<-.getselq(dq,tme)-1
  hotq<-.getselq(hq,tme)-1
  colq<-.getselq(cq,tme)-1
  #' Run C++ version of model with time-invarient vegetation and data.frame weather input
  air<-FALSE
  if (temp == "air") air<-TRUE
  bioclim<-runbioclim2Cpp(obstime,climdata,pointm,vegp,soilc,reqhgt,zref,ll$lats,ll$lons,
                          Sminp,Smaxp,tfact,matemp,out,wetq,dryq,hotq,colq,air)
  bior<-dtm
  index<-1
  for (i in 1:19) {
    if (out[i]) {
      bior<-c(bior,.rast(bioclim[[index]],dtm))
      index<-index+1
    }
  }
  bior<-bior[[-1]]
  bior<-mask(bior,dtm)
  nms<-paste0("bio",c(1:19))
  s<-which(out)
  nms<-nms[s]
  names(bior)<-nms
  return(bior)
}
# Dataframe climate, variable vegp
.runbioclim3<-function(weather, reqhgt = 0.05, vegp, soilc, dtm,
                       temp = "air", zref = 2, windhgt = zref, soilm = NA, runchecks = TRUE,
                       pai_a = NA, tfact = 1.5, out = rep(TRUE, 19), vegpisannual = TRUE) {
  # ====================== Unpack and check variables ======================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  cv<-.cleanvars(vegp,soilc,dtm) # sorts out NAs and zeros
  dtm<-cv$dtm
  vegp<-cv$vegp
  soilc<-cv$soilc
  if (runchecks) {
    rc<-checkinputs(weather,vegp,soilc,dtm)
    weather<-rc$weather
    vegp<-rc$vegp
    soilc<-rc$soilc
  }
  n<-dim(weather)[1]
  # ========================== Calculate quarters ============================ #
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  agg<-stats::aggregate(weather$precip,by=list(tme$mon),mean,na.rm=TRUE)$x
  wq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
  dq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # driest quarter
  agg<-stats::aggregate(weather$temp,by=list(tme$mon),sum,na.rm=TRUE)$x
  hq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # hottest quarter
  cq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
  # =================== Subset time-series and run point model =============== #
  bmpt<-.biomicropoint(weather,tme,reqhgt,vegp,soilc,dtm,zref,windhgt,soilm)
  micropoint<-bmpt$micropoint
  selw<-bmpt$sel
  weather<-micropoint$weather
  matemp<-micropoint$matemp
  # ==================== derive additional climate variables ================= #
  # derive additional climate variables
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                      hour=tme$hour+tme$min/60+tme$sec/3600)
  weather$es<-.satvap(weather$temp)
  weather$ea<-weather$es*weather$relhum/100
  weather$tdew<-.dewpoint(weather$ea,weather$temp)
  weather$relhum<-NULL
  weather$obs_time<-NULL
  # ======================== Sort out point model variables ================== #
  pointm<-micropoint$dfo
  if (reqhgt<0) {
    pointm$Tbp<-micropoint$Tbz
  } else pointm$Tbp<-0
  # ===================== Sort out vegp and add attidional terms ============= #
  vegp<-.sortvegp2(vegp,selw,vegpisannual,n)
  # add additional terms to vegp
  fd<-.foliageden(reqhgt,vegp$hgt,vegp$pai)
  if (class(pai_a) == "logical") vegp$paia<-fd$pai_a
  vegp$leafden<-fd$leafden
  # ========================== sort out soilc ================================ #
  soilp<-.sortsoilc(soilc,method="R")
  soilc$gref<-.is(soilc$groundr)
  soilc$groundr<-NULL
  soilc$soiltype<-NULL
  soilc$Smin<-soilp$Smin
  soilc$Smax<-soilp$Smax
  soilc$soilb<-soilp$soilb
  soilc$Psie<-soilp$psi_e
  soilc$Vq<-soilp$Vq
  soilc$Vm<-soilp$Vm
  soilc$Mc<-soilp$Mc
  soilc$rho<-soilp$rho
  # ============== Calculate slope, aspect and topographic wetness index ===== #
  soilc$slope<-terrain(dtm, v="slope")
  soilc$aspect<-terrain(dtm, v="aspect")
  soilc$twi<-.topidx(dtm)
  soilc$slope[is.na(soilc$slope)]<-0
  soilc$aspect[is.na(soilc$aspect)]<-0
  soilc$slope<-.is(mask(soilc$slope,dtm))
  soilc$aspect<-.is(mask(soilc$aspect,dtm))
  soilc$twi[is.na(soilc$twi)]<-1
  soilc$twi<-.is(mask(soilc$twi,dtm))
  # ======== Calculate sky view factor and wind shelter coefficient  etc ===== #
  sazi<-0
  ll<-.latlongfromraster(dtm)
  soilc$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) soilc$hor[,,i]<-.horizon(dtm,(i-1)*15)
  msl<-tan(apply(atan(soilc$hor),c(1,2),mean))
  soilc$svfa<-0.5*cos(2*msl)+0.5
  s<-1
  if (res(dtm)[1]<=100) s<-10
  soilc$wsa<-.windsheltera(dtm,micropoint$zref,s)
  Sminp<-.getmode(soilc$Smin)
  Smaxp<-.getmode(soilc$Smax)
  # ========================== Calculate selqs ============================ #
  wetq<-.getselq(wq,tme)-1
  dryq<-.getselq(dq,tme)-1
  hotq<-.getselq(hq,tme)-1
  colq<-.getselq(cq,tme)-1
  #' Run C++ version of model with time-invarient vegetation and data.frame weather input
  air<-FALSE
  if (temp == "air") air<-TRUE
  bioclim<-runbioclim3Cpp(obstime,weather,pointm,vegp,soilc,reqhgt,micropoint$zref,ll$lat,ll$long,
                          Sminp,Smaxp,tfact,matemp,out,wetq,dryq,hotq,colq,air)
  bior<-dtm
  index<-1
  for (i in 1:19) {
    if (out[i]) {
      bior<-c(bior,.rast(bioclim[[index]],dtm))
      index<-index+1
    }
  }
  bior<-bior[[-1]]
  bior<-mask(bior,dtm)
  nms<-paste0("bio",c(1:19))
  s<-which(out)
  nms<-nms[s]
  names(bior)<-nms
  return(bior)
}
#  Array climate, variable vegp
.runbioclim4<-function(climarray,tme,reqhgt=0.05,vegp,soilc,dtm,dtmc,temp="air",
                       zref = 2, windhgt = zref, soilm = NA, runchecks = TRUE, altcorrect=0,
                       pai_a = NA, tfact = 1.5, out = rep(TRUE,19),
                       vegpisannual = TRUE) {
  # ============================== Unpack variables  ========================= #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  cv<-.cleanvars(vegp,soilc,dtm) # sorts out NAs and zeros
  dtm<-cv$dtm
  vegp<-cv$vegp
  soilc<-cv$soilc
  # ========================== Calculate quarters ============================ #
  prech<-apply(.is(climarray$precip),3,mean,na.rm=TRUE)
  tc<-apply(.is(climarray$temp),3,mean,na.rm=TRUE)
  agg<-stats::aggregate(prech,by=list(tme$mon),mean,na.rm=TRUE)$x
  wq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
  dq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # driest quarter
  agg<-stats::aggregate(tc,by=list(tme$mon),sum,na.rm=TRUE)$x
  hq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # hottest quarter
  cq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
  # =============== subset climate array and run point models  =============== #
  bmpt<-.biomicropoint(climarray,tme,reqhgt,vegp,soilc,dtm,zref,windhgt,soilm)
  selw<-bmpt$sel
  micropointa<-bmpt$micropoint
  weather<-list()
  # ==================== extract data from micropointa ======================= #
  dfo<-list()
  Tbz<-list()
  mat<-rep(NA,lengthmicropointa())
  ii<-0
  for (i in 1:length(micropointa)) {
    onepoint<-micropointa[[i]]
    if (class(onepoint) != "logical") {
      weather[[i]]<-onepoint$weather
      dfo[[i]]<-onepoint$dfo
      tme<-as.POSIXlt(weather[[i]]$obs_time,tz="UTC")
      tmeorig<-onepoint$tmeorig
      subs<-onepoint$subs
      zref<-onepoint$zref
      mat[i]<-onepoint$matemp
      if (class(onepoint$Tbz) != "logical") {
        Tbz[[i]]<-data.frame(Tbz=onepoint$Tbz)
        ii<-i
      }
      if (runchecks) {
        rc<-checkinputs(weather[[i]],vegp,soilc,dtm,onepoint$zref)
        weather[[i]]<-rc$weather
        vegp<-rc$vegp
        soilc<-rc$soilc
      }
    } else {
      weather[[i]]<-NA
      dfo[[i]]<-NA
      Tbz[[i]]<-NA
    }
  }
  if (class(dtmc)[1] == "PackedSpatRaster") {
    r<-rast(dtmc)
  } else r<-dtmc
  # derive additional climate variables
  matemp<-mean(mat,na.rm=TRUE)
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                      hour=tme$hour+tme$min/60+tme$sec/3600)
  # ================== Turn climate variables into arrays ============= #
  climdata<-list()
  h<-length(tme)
  climdata$tc<-.cca(weather,"temp",h,dtmc,dtm)
  pk<-.cca(weather,"pres",h,dtmc,dtmc)
  relhum<-.cca(weather,"relhum",h,dtmc,dtm)
  climdata$es<-.satvap(climdata$tc)
  climdata$ea<-climdata$es*relhum/100
  climdata$tdew<-.dewpoint(climdata$ea,climdata$tc)
  # Apply altitudinal correction
  if (altcorrect == 0) {  # No altitudinal correction
    climdata$pk<-as.array(resample(.rast(pk,dtmc),dtm))
  } else { # Altitudinal correction applied
    dtmc[is.na(dtmc)]<-0
    psl<-pk/(((293-0.0065*.rta(dtmc,h))/293)^5.26)
    psl<-as.array(resample(.rast(psl,dtmc),dtm))
    climdata$pk<-psl*(((293-0.0065*.rta(dtm,h))/293)^5.26)
    dc<-resample(dtmc,dtm)
    elevd<-.rta(dc-dtm,h)
    if (altcorrect==1) {  # Fixed lapse rate
      tcdif<-elevd*(5/1000)
    } else { # Humidity-dependent lapse rate
      lr<-.lapserate(climdata$tc,climdata$ea,climdata$pk)
      tcdif<-lr*elevd
    }
    climdata$tc<-tcdif+climdata$tc
  }
  climdata$difrad<-.cca(weather,"difrad",h,dtmc,dtm)
  climdata$swdown<-.cca(weather,"swdown",h,dtmc,dtm)
  climdata$lwdown<-.cca(weather,"lwdown",h,dtmc,dtm)
  # Wind speed and direction
  u2<-.cca(weather,"windspeed",h,dtmc,dtmc)
  wd<-.cca(weather,"winddir",h,dtmc,dtmc)
  wu<-u2*cos(wd*pi/180)
  wv<-u2*sin(wd*pi/180)
  wuv<-apply(wu,3,mean,na.rm=TRUE)
  wvv<-apply(wv,3,mean,na.rm=TRUE)
  wu<-as.array(resample(.rast(wu,dtmc),dtm))
  wv<-as.array(resample(.rast(wv,dtmc),dtm))
  climdata$windspeed<-sqrt(wu^2+wv^2)
  climdata$winddir<-(atan2(wvv,wuv)*180/pi)%%360
  # Point model variables
  # =============== Turn point model variables into arrays ============= #
  pointm<-list()
  pointm$umu<-.cca(dfo,"umu",h,dtmc,dtm)
  pointm$kp<-.cca(dfo,"kp",h,dtmc,dtm)
  pointm$muGp<-.cca(dfo,"muGp",h,dtmc,dtm)
  pointm$DDp<-.cca(dfo,"DDp",h,dtmc,dtm)
  pointm$T0p<-.cca(dfo,"T0p",h,dtmc,dtm)
  pointm$dtrp<-.cca(dfo,"dtrp",h,dtmc,dtm)
  pointm$Gp<-.cca(dfo,"G",h,dtmc,dtm)
  pointm$soilm<-.cca(dfo,"soilm",h,dtmc,dtm)
  pointm$Tg<-.cca(dfo,"Tg",h,dtmc,dtm)
  if (reqhgt<0)  {
    pointm$Tbp<-.cca(Tbz,"Tbz",h,dtmc,dtm)
  } else  pointm$Tbp<-pointm$soilm*0
  # ===================== Sort out vegp ================================== #
  vegp<-.sortvegp2(vegp,selw,vegpisannual,n)
  # add additional terms to vegp
  fd<-.foliageden(reqhgt,vegp$hgt,vegp$pai)
  if (class(pai_a) == "logical") vegp$paia<-fd$pai_a
  vegp$leafden<-fd$leafden
  # ======================== Sort out soilc ================================== #
  soilp<-.sortsoilc(soilc,method="R")
  soilc$gref<-.is(soilc$groundr)
  soilc$groundr<-NULL
  soilc$soiltype<-NULL
  soilc$Smin<-soilp$Smin
  soilc$Smax<-soilp$Smax
  soilc$soilb<-soilp$soilb
  soilc$Psie<-soilp$psi_e
  soilc$Vq<-soilp$Vq
  soilc$Vm<-soilp$Vm
  soilc$Mc<-soilp$Mc
  soilc$rho<-soilp$rho
  # ============= Calculate slope, aspect and topographic wetness index ====== #
  soilc$slope<-terrain(dtm, v="slope")
  soilc$aspect<-terrain(dtm, v="aspect")
  soilc$twi<-.topidx(dtm)
  soilc$slope[is.na(soilc$slope)]<-0
  soilc$aspect[is.na(soilc$aspect)]<-0
  soilc$slope<-.is(mask(soilc$slope,dtm))
  soilc$aspect<-.is(mask(soilc$aspect,dtm))
  soilc$twi[is.na(soilc$twi)]<-1
  soilc$twi<-.is(mask(soilc$twi,dtm))
  # ====== Calculate sky view factor and wind shelter coefficient etc ======== #
  ll<-.latslonsfromr(dtm)
  soilc$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) soilc$hor[,,i]<-.horizon(dtm,(i-1)*15)
  msl<-tan(apply(atan(soilc$hor),c(1,2),mean))
  soilc$svfa<-0.5*cos(2*msl)+0.5
  s<-1
  if (res(dtm)[1]<=100) s<-10
  soilc$wsa<-.windsheltera(dtm,zref,s)
  Sminp<-.getmode(soilc$Smin)
  Smaxp<-.getmode(soilc$Smax)
  # ========================== Calculate selqs ============================ #
  wetq<-.getselq(wq,tme)-1
  dryq<-.getselq(dq,tme)-1
  hotq<-.getselq(hq,tme)-1
  colq<-.getselq(cq,tme)-1
  #' Run C++ version of model with time-invarient vegetation and data.frame weather input
  air<-FALSE
  if (temp == "air") air<-TRUE
  bioclim<-runbioclim4Cpp(obstime,climdata,pointm,vegp,soilc,reqhgt,zref,ll$lats,ll$lons,
                          Sminp,Smaxp,tfact,matemp,out,wetq,dryq,hotq,colq,air)
  bior<-dtm
  index<-1
  for (i in 1:19) {
    if (out[i]) {
      bior<-c(bior,.rast(bioclim[[index]],dtm))
      index<-index+1
    }
  }
  bior<-bior[[-1]]
  bior<-mask(bior,dtm)
  nms<-paste0("bio",c(1:19))
  s<-which(out)
  nms<-nms[s]
  names(bior)<-nms
  return(bior)
}
#' calculate average vegetation value for hours with snow
.sortl<-function(vegp,sdep) {
  nms<-c("pai","hgt","leaft","clump")
  vegp2<-list()
  for (i in 1:4) {
    nn<-names(vegp)
    ss<-which(nn==nms[i])
    a<-.is(vegp[[ss]])
    d<-dim(a)
    if (length(d) > 2) {
      dmx<-dim(a)[3]
    } else dmx<-1
    if (dmx > 1) {
      n<-length(sdep)
      s<-round(seq(0.50001,dmx+0.5,length.out=n),0)
      s[s>dmx]<-dmx
      s[s<1]<-1
      sel<-which(sdep>0)
      if (length(sel) > 0) {
        s<-s[sel]
      } else s<-s[1]
      tb<-table(s)
      num<-as.numeric(names(tb))
      fre<-as.numeric(tb)
      m<-a[,,1]*0
      for (j in 1:length(num)) m<-m+a[,,j]*fre[j]
      m<-m/sum(fre)
      vegp2[[i]]<-m
    } else vegp2[[i]]<-.is(vegp[[ss]])
  }
  names(vegp2)<-nms
  return(vegp2)
}
#' calculate average vegetation value for hours with snow with paia and foliage density
.sortl2<-function(vegp,sdep,reqhgt,paia) {
  nms<-c("pai","hgt","leaft","clump","leafd")
  vegp2<-list()
  for (i in 1:5) {
    nn<-names(vegp)
    ss<-which(nn==nms[i])
    a<-.is(vegp[[ss]])
    d<-dim(a)
    if (length(d) > 2) {
      dmx<-dim(a)[3]
    } else dmx<-1
    if (dmx > 1) {
      n<-length(sdep)
      s<-round(seq(0.50001,dmx+0.5,length.out=n),0)
      s[s>dmx]<-dmx
      s[s<1]<-1
      sel<-which(sdep>0)
      if (length(sel) > 0) {
        s<-s[sel]
      } else s<-s[1]
      tb<-table(s)
      num<-as.numeric(names(tb))
      fre<-as.numeric(tb)
      m<-a[,,1]*0
      for (j in 1:length(num)) m<-m+a[,,j]*fre[j]
      m<-m/sum(fre)
      vegp2[[i]]<-m
    } else vegp2[[i]]<-.is(vegp[[ss]])
  }
  names(vegp2)<-nms
  fde<-.foliageden(reqhgt,vegp2$hgt,vegp2$pai)
  if (class(paia)[1] == "logical") {
    vegp2$paia<-fde$pai_a
  } else vegp2$paia<-paia
  vegp2$leafden<-fde$leafden
  return(vegp2)
}
#' Calculates topographic positioning index
.tpicalc<-function(af,me,dtm,tfact) {
  # Calculate topographic positioning index
  if (af < me/2) {
    dtmc<-aggregate(dtm,af,na.rm=T)
    dtmc<-resample(dtmc,dtm)
  } else dtmc<-dtm*0+mean(.is(dtm),na.rm=T)
  tpi<-dtmc-dtm
  tpic<-exp(.is(tpi)*tfact)
  tpic[tpic<0.05]<-0.1
  tpic[tpic>10]<-10
  me<-mean(tpic,na.rm=TRUE)
  tpic<-tpic/me
  tpic
}
#' full snow model with data.frame weather inputs
.snowmodel1<-function(weather,dtm,vegp,soilc,snowenv="Taiga",snowinitd = 0, snowinita = 0, zref = 2, windhgt = zref, tfact = 0.02) {
  # ================== unpack variables ==================================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  # =========== Perform wind height adjustment if necessary ===================== #
  if (zref != windhgt) {
    weather$windspeed<-weather$windspeed*log(67.8*zref-5.42)/log(67.8*windhgt-5.42)
  }
  # =========== Perform height adjustment if necessary ===================== #
  if (max(.is(vegp$hgt),na.rm=TRUE) > zref) {
    tme<-as.POSIXlt(weather$obs_time,tz="UTC")
    obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                        hour=tme$hour+tme$min/60+tme$sec/3600)
    zout<-max(.is(vegp$hgt))
    ll<-.latlongfromraster(dtm)
    tst<-1
    ctr<-0
    while (tst > 0) {
      weather2<-weatherhgtCpp(obstime,weather,zref,zref,zout,ll$lat,ll$long)
      tt<-mean(weather2$temp)
      if (is.na(tt) == FALSE) tst<-0
      ctr<-ctr+1
      if (ctr > 5) tst<-0
    }
    weather<-weather2
    weather$obs_time<-as.POSIXlt(tme)
    zref<-zout
  }
  # ============ Get variables for pointmodel =============================== #
  vegpp<-.sortvegp(vegp,method="P")
  vegpp<-c(vegpp[2],vegpp[1],vegpp[6],vegpp[4])
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,hour=tme$hour)
  up<-.unpack
  if (class(dtm) == "PackedSpatRaster") dtm<-rast(dtm)
  sdep<-.is(dtm)*0+snowinitd
  sage<-.is(dtm)*0+snowinita
  ll<-.latlongfromraster(dtm)
  other<-c(0,0,ll$lat,ll$long,zref,mean(sdep,na.rm=T),mean(sage,na.rm=T))
  # ========== Run point model and convert to input data.frame =============== #
  pmod<-pointmodelsnow(obstime,weather,vegpp,other,snowenv)
  pointm<-data.frame(Gp=pmod$G,Tc=pmod$Tc,RswabsG=pmod$RswabsG,RlwabsG=pmod$RlwabsG,
                     umu=pmod$umu,tr=pmod$tr)
  n<-dim(weather)[1]
  sdept<-pmod$sdepc[1:n]
  # ========== Prepare vegetation inputs for grid model =============== #
  vegp<-.sortl(vegp,sdept)
  # ========== Prepare other inputs for grid model =============== #
  other<-list()
  other$zref<-zref
  other$lat<-ll$lat
  other$lon<-ll$long
  other$isnowdc<-sdep
  other$isnowac<-sage
  other$isnowdg<-sdep*0.5
  other$isnowag<-sage
  # ========== Run grid model in daily timesteps =============== #
  h<-length(tme)
  n5days<-h/(24*5)
  pb <- txtProgressBar(min = 0, max = n5days, style = 3)
  Tc<-array(NA,dim=c(dim(dtm)[1:2],h))
  Tg<-Tc
  snowdepg<-Tc
  swe<-Tc
  sden<-Tc
  smx<-length(tme)
  dtms<-dtm+.rast(sdep*0.5,dtm)
  for (day in 1:n5days) {
    other$slope<-terra::terrain(dtms,v="slope")
    other$slope[is.na(other$slope)]<-0
    other$slope<-.is(mask(other$slope,dtm))
    other$aspect<-terra::terrain(dtms,v="aspect")
    other$aspect[is.na(other$aspect)]<-180
    other$aspect<-.is(mask(other$aspect,dtm))
    other$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
    for (i in 1:24) other$hor[,,i]<-.horizon(dtms,(i-1)*15)
    msl<-tan(apply(atan(other$hor),c(1,2),mean))
    other$skyview<-0.5*cos(2*msl)+0.5
    ss<-1
    if (res(dtm)[1]<=100) ss<-10
    other$wsa<-.windsheltera(dtms,zref,ss)
    # calculate s
    st<-(day-1)*(5*24)+1
    ed<-st+(5*24)-1
    if (ed > h) ed<-h
    s<-c(st:ed)
    smod<- gridmodelsnow1(obstime[s,],weather[s,],pointm[s,],vegp,other,snowenv)
    # Topographic snow re-distribution
    tpr<-10*mean(weather$windspeed[s])^0.5
    af<-round(tpr/res(dtm)[1],0)
    me<-min(dim(dtm)[1:2])
    tpi<-.tpicalc(af,me,dtms,tfact)
    # ground snow change
    asd<-.rta(other$isnowdg,length(s))
    dsnow<-smod$sdepg-asd
    dsnow2<-dsnow*.rta(tpi,length(s)) # topographically adjusted snow change
    sss<-which(dsnow<0)
    dsnow2[sss]<-dsnow[sss]
    # canopy only snow change
    asc<-.rta(other$isnowdc,length(s))
    cdsnow<-smod$sdepc-asc-dsnow
    # Assign to variables
    Tc[,,s]<-smod$Tc
    Tg[,,s]<-smod$Tg
    swe[,,s]<-(asc+cdsnow+dsnow2)*smod$sden
    snowdepg[,,s]<-asd+dsnow2
    sden[,,s]<-smod$sden
    other$isnowdc<-(asc+cdsnow+dsnow2)[,,length(s)]
    other$isnowac<-(asd+dsnow2)[,,length(s)]
    other$isnowac<-smod$agec
    other$isnowag<-smod$ageg
    nnn<-s[length(s)]
    dtms<-dtm+.rast(snowdepg[,,nnn],dtm)
    setTxtProgressBar(pb, day)
    ## NB need to return snow age from c++ code
  }
  setTxtProgressBar(pb, n5days)
  return(list(Tc=Tc,Tg=Tg,groundsnowdepth=snowdepg,totalSWE=swe,snowden=sden,umu=pointm$umu))
}
#' resample melt
.resamplemelt<-function(melt, sbtn, dtmc, dtm) {
  melts<-apply(melt[,,sbtn],c(1,2),sum)
  melts<-resample(.rast(melts,dtmc),dtm)
  return(.is(melts))
}
#' quick snow model with data.frame weather inputs
.snowmodelq1<-function(weather,dtm,vegp,soilc,subs,snowenv="Taiga",snowinitd = 0,
                       snowinita = 0, zref = 2, windhgt = zref, tfact = 0.02) {
  # ================== unpack variables ==================================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  # =========== Perform wind height adjustment if necessary ===================== #
  if (zref != windhgt) {
    weather$windspeed<-weather$windspeed*log(67.8*zref-5.42)/log(67.8*windhgt-5.42)
  }
  # =========== Perform height adjustment if necessary ===================== #
  if (max(.is(vegp$hgt),na.rm=TRUE) > zref) {
    tme<-as.POSIXlt(weather$obs_time,tz="UTC")
    obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                        hour=tme$hour+tme$min/60+tme$sec/3600)
    zout<-max(.is(vegp$hgt))
    ll<-.latlongfromraster(dtm)
    tst<-1
    ctr<-0
    while (tst > 0) {
      weather2<-weatherhgtCpp(obstime,weather,zref,zref,zout,ll$lat,ll$long)
      tt<-mean(weather2$temp)
      if (is.na(tt) == FALSE) tst<-0
      ctr<-ctr+1
      if (ctr > 5) tst<-0
    }
    weather<-weather2
    weather$obs_time<-as.POSIXlt(tme)
    zref<-zout
  }
  # ============ Get variables for pointmodel =============================== #
  vegpp<-.sortvegp(vegp,method="P")
  vegpp<-c(vegpp[2],vegpp[1],vegpp[6],vegpp[4])
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,hour=tme$hour)
  up<-.unpack
  if (class(dtm) == "PackedSpatRaster") dtm<-rast(dtm)
  sdep<-.is(dtm)*0+snowinitd
  sage<-.is(dtm)*0+snowinita
  ll<-.latlongfromraster(dtm)
  other<-c(0,0,ll$lat,ll$long,zref,mean(sdep,na.rm=T),mean(sage,na.rm=T))
  # ========== Run point model and convert to input data.frame =============== #
  pmod<-pointmodelsnow(obstime,weather,vegpp,other,snowenv,maxiter=20)
  s<-c(1:length(tme))+1
  pointm<-data.frame(Gp=pmod$G,Tc=pmod$Tc,RswabsG=pmod$RswabsG,RlwabsG=pmod$RlwabsG,
                     umu=pmod$umu,tr=pmod$tr,sdepc=pmod$sdepc[s],sdepg=pmod$sdepg[s],
                     sublmelt=pmod$sublmelt,tempmelt=pmod$tempmelt,rainmelt=pmod$rainmelt,
                     sstemp=pmod$sstemp)
  # ========== Calculate cumulative melt etc =============================== #
  snow<-ifelse(weather$temp>2,0,weather$prec)
  # ========== Prepare vegetation inputs for grid model =============== #
  n<-dim(weather)[1]
  sdept<-pmod$sdepc[1:n]
  vegp<-.sortl(vegp,sdept)
  vegp$leaft[is.na(vegp$leaft)]<-0.01
  # ===================== subset point models ================================ #
  weathero<-weather
  pointm<-pointm[subs,]
  weather<-weather[subs,]
  obstime<-obstime[subs,]
  # ========== Prepare other inputs for grid model =========================== #
  other<-list()
  other$zref<-zref
  other$lat<-ll$lat
  other$lon<-ll$long
  other$isnowdc<-snowinitd*.is(dtm)
  other$isnowac<-sage
  other$isnowag<-sage
  other$slope<-terra::terrain(dtm,v="slope")
  other$slope[is.na(other$slope)]<-0
  other$slope<-.is(mask(other$slope,dtm))
  other$aspect<-terra::terrain(dtm,v="aspect")
  other$aspect[is.na(other$aspect)]<-180
  other$aspect<-.is(mask(other$aspect,dtm))
  other$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) other$hor[,,i]<-.horizon(dtm,(i-1)*15)
  msl<-tan(apply(atan(other$hor),c(1,2),mean))
  other$skyview<-0.5*cos(2*msl)+0.5
  ss<-1
  if (res(dtm)[1]<=100) ss<-10
  other$wsa<-.windsheltera(dtm,zref,ss)
  weather$obs_time<-NULL
  # ========== Create arrays for storing data ============================== #
  n<-dim(weather)[1]
  Tc<-array(NA,dim=c(dim(dtm)[1:2],n))
  Tg<-Tc
  sdepc<-Tc
  sdepg<-array(0,dim=dim(Tc))
  sden<-Tc
  # Calculate typical canopy interception
  msnow<-mean(snow[snow>0])
  mtemp<-mean(weather$temp)
  intfrac<-canintfrac(vegp$hgt,vegp$pai,2,msnow,mtemp,0)
  other$isnowdg<-(1-intfrac)*other$isnowdc
  # ========== Loop through each day ======================================= #
  ndays<-length(weather$temp)/24
  mup<-.is(dtm)*0+1
  ped<-0
  pb <- txtProgressBar(min = 0, max = ndays, style = 3)
  for (day in 1:ndays) {
    st<-(day-1)*24+1
    ed<-st+23
    s<-c(st:ed)
    # work out snow balance in between current and previous model run
    if ((subs[st]-1) > 1) {
      sbtn<-c((ped+1):(subs[st]-1))
      mu<-meltmu(other$skyview,pmod$sstemp[sbtn],weathero$temp[sbtn])
      melt<-sum(pmod$sublmelt[sbtn])+sum(pmod$rainmelt[sbtn])+mu*sum(pmod$tempmelt[sbtn])
      balancec<-sum(snow[sbtn]/1000)-melt  # m SWE
      balanceg<-(1-intfrac)*sum(snow[sbtn]/1000)-exp(-vegp$pai)*melt
    } else {
      balancec<-0
      balanceg<-0
    }
    # adjust initial snow depths accordingly
    other$isnowdc<-other$isnowdc+balancec*(1000/mean(pmod$sdenc[sbtn]))
    other$isnowdg<-other$isnowdg+balanceg*(1000/mean(pmod$sdeng[sbtn]))
    other$isnowdc[other$isnowdc<0]<-0
    other$isnowdg[other$isnowdg<0]<-0
    # Run snow model
    smod<-gridmodelsnow1(obstime[s,],weather[s,],pointm[s,],vegp,other,snowenv)
    # Redistribute snow
    # Calculate change in snow depth
    dsnow<-smod$sdepc-.rta(other$isnowdc,24)  # ground and canopy
    dsnowg<-smod$sdepg-.rta(other$isnowdg,24) # ground only
    dsnowc<-dsnow-dsnowg # canopy only
    # Redistribute snow according to terrain position index
    dtms<-dtm+.rast(sdepg[,,ed],dtm)
    tpr<-10*mean(weather$windspeed[s])^0.5
    af<-round(tpr/res(dtm)[1],0)
    me<-min(dim(dtm)[1:2])
    tpi<-.tpicalc(af,me,dtms,tfact)
    dsnowg2<-dsnowg*.rta(tpi,length(s))
    dsnowc2<-dsnowc+dsnowg2
    # Assign variables
    Tc[,,s]<-smod$Tc
    Tg[,,s]<-smod$Tg
    sden[,,s]<-smod$sden
    sdc<-dsnowc2+.rta(other$isnowdc,24)
    sdg<-dsnowg2+.rta(other$isnowdg,24)
    sdc[sdc<0]<-0
    sdg[sdg<0]<-0
    sdepc[,,s]<-sdc
    sdepg[,,s]<-sdg
    # Recalculate ped
    pred<-subs[ed]
    utils::setTxtProgressBar(pb, day)
  }
  swe<-sdepc*sden
  return(list(Tc=Tc,Tg=Tg,groundsnowdepth=sdepg,totalSWE=swe,snowden=sden,umu=pointm$umu))
}
#' full snow model with array weather inputs
.snowmodel2<-function(weather,tme,dtm,dtmc,vegp,soilc,altcorrect = 0,snowenv="Taiga",
                      snowinitd = 0, snowinita = 0, zref = 2, windhgt = zref, tfact = 0.02) {
  # ================== unpack variables ==================================== #
  n5days<-length(tme)/(24*5)
  pb <- utils::txtProgressBar(min = 0, max = n5days*2, style = 3)
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  # ================= Sort out all the rasters ============================== #
  if (crs(weather$temp) != crs(dtm))  weather$temp<-project(weather$temp,crs(dtm))
  if (crs(weather$relhum) != crs(dtm))  weather$relhum<-project(weather$relhum,crs(dtm))
  if (crs(weather$pres) != crs(dtm))  weather$pres<-project(weather$pres,crs(dtm))
  if (crs(weather$swdown) != crs(dtm))  weather$swdown<-project(weather$swdown,crs(dtm))
  if (crs(weather$difrad) != crs(dtm))  weather$difrad<-project(weather$difrad,crs(dtm))
  if (crs(weather$lwdown) != crs(dtm))  weather$lwdown<-project(weather$lwdown,crs(dtm))
  if (crs(weather$windspeed) != crs(dtm))  weather$windspeed<-project(weather$windspeed,crs(dtm))
  if (crs(weather$precip) != crs(dtm))  weather$precip<-project(weather$precip,crs(dtm))
  wdir<-apply(.is(weather$winddir),3,.getmode)
  # ============= Run point model for each grid cell of clim array ========= #
  vegpc<-.resamplev(vegp,weather$temp)
  dms<-dim(weather$temp)
  n<-dms[1]*dms[2]
  ii<-0.5*n5days/n # for progress bar
  pointms<-list()
  clims<-list()
  mxhgt<-max(.is(vegp$hgt),na.rm=T)
  if (class(dtmc) == "PackedSpatRaster") dtmc<-rast(dtm)
  ll<-.latslonsfromr(dtmc)
  k<-1
  for (i in 1:dms[1]) {
    for (j in 1:dms[2]) {
      climdf<-.todf(weather,i,j,tme,wdir)
      if (is.na(climdf$temp[1]) == FALSE & is.na(.is(vegpc$hgt)[i,j]) == FALSE) {
        # ======== Perform wind height adjustment if necessary ============== #
        if (zref != windhgt) {
          climdf$windspeed<-climdf$windspeed*log(67.8*zref-5.42)/log(67.8*windhgt-5.42)
        }
        # =========== Perform height adjustment if necessary ================ #
        obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                            hour=tme$hour+tme$min/60+tme$sec/3600)
        if (mxhgt > zref) {
          tst<-1
          ctr<-0
          while (tst > 0) {
            weather2<-weatherhgtCpp(obstime,climdf,zref,zref,mxhgt,ll$lats[i,j],ll$lons[i,j])
            tt<-mean(weather2$temp)
            if (is.na(tt) == FALSE) tst<-0
            ctr<-ctr+1
            if (ctr > 5) tst<-0
          }
          climdf<-weather2
          climdfr$obs_time<-as.POSIXlt(tme)
          zref<-mxhgt
        } # end height adjust
        # ============ Get variables for pointmodel =========== ============== #
        vegpp<-.tovp(vegpc,i,j)
        vegpp<-c(vegpp[2],vegpp[1],vegpp[6],vegpp[4])
        other<-c(0,0,ll$lats[i,j],ll$lons[i,j],zref,snowinitd,snowinita)
        # ========== Run point model and convert to input data.frame ========= #
        pmod<-pointmodelsnow(obstime,climdf,vegpp,other,snowenv,maxiter=10)
        pointms[[k]]<-data.frame(Gp=pmod$G,Tc=pmod$Tc,RswabsG=pmod$RswabsG,
                                 RlwabsG=pmod$RlwabsG,umu=pmod$umu,tr=pmod$tr,sdepc=pmod$sdepc[1:length(tme)])
        clims[[k]]<-climdf
      } else {
        pointms[[k]]<-NA
        clims[[k]]<-NA
      }
      xx<-0.5*k*n5days/n
      utils::setTxtProgressBar(pb, k*ii)
      k<-k+1
    } # end j
  } # end i
  # ================ Turn climate data into arrays =========================== #
  climdata<-list()
  h<-length(tme)
  ii<-0.5*n5days/15
  climdata$temp<-.cca(clims,"temp",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 1*ii+n5days/2)
  pk<-.cca(clims,"pres",h,dtmc,dtmc)
  climdata$relhum<-.cca(clims,"relhum",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 2*ii+n5days/2)
  # Apply altitudinal correction
  if (altcorrect == 0) {  # No altitudinal correction
    climdata$pres<-as.array(resample(.rast(pk,dtmc),dtm))
  } else { # Altitudinal correction applied
    es<-.satvap(climdata$temp)
    ea<-es*climdata$relhum/100
    dtmc[is.na(dtmc)]<-0
    psl<-pk/(((293-0.0065*.rta(dtmc,h))/293)^5.26)
    psl<-as.array(resample(.rast(psl,dtmc),dtm))
    climdata$pres<-psl*(((293-0.0065*.rta(dtm,h))/293)^5.26)
    dc<-resample(dtmc,dtm)
    elevd<-.rta(dc-dtm,h)
    if (altcorrect==1) {  # Fixed lapse rate
      tcdif<-elevd*(5/1000)
    } else { # Humidity-dependent lapse rate
      lr<-.lapserate(climdata$temp,ea,climdata$pres)
      tcdif<-lr*elevd
    }
    climdata$temp<-tcdif+climdata$temp
    climdata$relhum<-(ea/.satvap(climdata$temp))*100
  }
  utils::setTxtProgressBar(pb, 3*ii+n5days/2)
  climdata$relhum[climdata$relhum>100]<-100
  utils::setTxtProgressBar(pb, 4*ii+n5days/2)
  climdata$difrad<-.cca(clims,"difrad",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 5*ii+n5days/2)
  climdata$swdown<-.cca(clims,"swdown",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 6*ii+n5days/2)
  climdata$lwdown<-.cca(clims,"lwdown",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 7*ii+n5days/2)
  climdata$precip<-.cca(clims,"precip",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 8*ii+n5days/2)
  # Wind speed and direction
  u2<-.cca(clims,"windspeed",h,dtmc,dtmc)
  wd<-.cca(clims,"winddir",h,dtmc,dtmc)
  wu<-u2*cos(wd*pi/180)
  wv<-u2*sin(wd*pi/180)
  wuv<-apply(wu,3,mean,na.rm=TRUE)
  wvv<-apply(wv,3,mean,na.rm=TRUE)
  wu<-as.array(resample(.rast(wu,dtmc),dtm))
  wv<-as.array(resample(.rast(wv,dtmc),dtm))
  climdata$windspeed<-sqrt(wu^2+wv^2)
  climdata$winddir<-(atan2(wvv,wuv)*180/pi)%%360
  utils::setTxtProgressBar(pb, 9*ii+n5days/2)
  # ================ Turn snow point variables into arrays ================== #
  pointm<-list()
  pointm$Gp<-.cca(pointms,"Gp",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 10*ii+n5days/2)
  pointm$Tc<-.cca(pointms,"Tc",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 11*ii+n5days/2)
  pointm$RswabsG<-.cca(pointms,"RswabsG",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 12*ii+n5days/2)
  pointm$RlwabsG<-.cca(pointms,"RlwabsG",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 13*ii+n5days/2)
  pointm$umu<-.cca(pointms,"umu",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 14*ii+n5days/2)
  pointm$tr<-.cca(pointms,"tr",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 15*ii+n5days/2)
  # ========== Prepare vegetation inputs for grid model =============== #
  sdepc<-.cca(pointms,"sdepc",h,dtmc,dtmc)
  sdept<-apply(sdepc,3,max,na.rm=TRUE)
  vegp<-.sortl(vegp,sdept)
  # ========== Prepare other inputs for grid model =============== #
  ll<-.latslonsfromr(dtm)
  other<-list()
  other$zref<-zref
  other$lats<-ll$lats
  other$lons<-ll$lons
  sdep<-.is(dtm)*0+snowinitd
  other$isnowdc<-sdep
  other$isnowac<-.is(dtm)*0+snowinita
  other$isnowdg<-other$isnowdc*0.5
  other$isnowag<-.is(dtm)*0+snowinita
  # ========== Run grid model in daily timesteps =============== #
  Tc<-array(NA,dim=c(dim(dtm)[1:2],h))
  Tg<-Tc
  snowdepg<-Tc
  swe<-Tc
  sden<-Tc
  smx<-length(tme)
  dtms<-dtm+.rast(sdep*0.5,dtm)
  for (day in 1:n5days) {
    other$slope<-terra::terrain(dtms,v="slope")
    other$slope[is.na(other$slope)]<-0
    other$slope<-.is(mask(other$slope,dtm))
    other$aspect<-terra::terrain(dtms,v="aspect")
    other$aspect[is.na(other$aspect)]<-180
    other$aspect<-.is(mask(other$aspect,dtm))
    other$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
    for (i in 1:24) other$hor[,,i]<-.horizon(dtms,(i-1)*15)
    msl<-tan(apply(atan(other$hor),c(1,2),mean))
    other$skyview<-0.5*cos(2*msl)+0.5
    ss<-1
    if (res(dtm)[1]<=100) ss<-10
    other$wsa<-.windsheltera(dtms,zref,ss)
    # calaculate s
    st<-(day-1)*(5*24)+1
    ed<-st+(5*24)-1
    s<-c(st:ed)
    s<-s[s<smx]
    # subset weather
    climdata1<-list()
    for (i in 1:8) climdata1[[i]]<-climdata[[i]][,,s]
    climdata1[[9]]<-climdata[[9]][s]
    names(climdata1)<-names(climdata)
    pointm1<-list()
    for (i in 1:6) pointm1[[i]]<-pointm[[i]][,,s]
    names(pointm1)<-names(pointm)
    vegp$leaft[is.na(vegp$leaft)]<-0.001
    smod<- gridmodelsnow2(obstime[s,],climdata1,pointm1,vegp,other,snowenv)
    # Topographic snow re-distribution
    wss<-sqrt(wuv[s]^2+wvv[s]^2)
    tpr<-10*mean(wss)^0.5
    af<-round(tpr/res(dtm)[1],0)
    me<-min(dim(dtm)[1:2])
    tpi<-.tpicalc(af,me,dtms,tfact)
    # ground snow change
    asd<-.rta(other$isnowdg,length(s))
    dsnow<-smod$sdepg-asd
    dsnow2<-dsnow*.rta(tpi,length(s)) # topographically adjusted snow change
    sss<-which(dsnow<0)
    dsnow2[sss]<-dsnow[sss]
    # canopy only snow change
    asc<-.rta(other$isnowdc,length(s))
    cdsnow<-smod$sdepc-asc-dsnow
    # Assign to variables
    Tc[,,s]<-smod$Tc
    Tg[,,s]<-smod$Tg
    swe[,,s]<-(asc+cdsnow+dsnow2)*smod$sden
    snowdepg[,,s]<-asd+dsnow2
    sden[,,s]<-smod$sden
    other$isnowdc<-(asc+cdsnow+dsnow2)[,,length(s)]
    other$isnowac<-(asd+dsnow2)[,,length(s)]
    other$isnowac<-smod$agec
    other$isnowag<-smod$ageg
    nnn<-s[length(s)]
    dtms<-dtm+.rast(snowdepg[,,nnn],dtm)
    utils::setTxtProgressBar(pb, day+n5days)
  }
  setTxtProgressBar(pb, n5days*2)
  return(list(Tc=Tc,Tg=Tg,groundsnowdepth=snowdepg,totalSWE=swe,snowden=sden,umu=pointm$umu))
}
#' quick snow model with array weather inputs
.snowmodelq2<-function(weather,tme,dtm,dtmc,vegp,soilc,subs,altcorrect=0,snowenv="Taiga",
                       snowinitd = 0, snowinita = 0, zref = 2, windhgt = zref, tfact = 0.02) {
  ndays<-length(subs)/24
  pb <- utils::txtProgressBar(min = 0, max = ndays*2, style = 3)
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  # ================= Sort out all the rasters ============================== #
  if (crs(weather$temp) != crs(dtm))  weather$temp<-project(weather$temp,crs(dtm))
  if (crs(weather$relhum) != crs(dtm))  weather$relhum<-project(weather$relhum,crs(dtm))
  if (crs(weather$pres) != crs(dtm))  weather$pres<-project(weather$pres,crs(dtm))
  if (crs(weather$swdown) != crs(dtm))  weather$swdown<-project(weather$swdown,crs(dtm))
  if (crs(weather$difrad) != crs(dtm))  weather$difrad<-project(weather$difrad,crs(dtm))
  if (crs(weather$lwdown) != crs(dtm))  weather$lwdown<-project(weather$lwdown,crs(dtm))
  if (crs(weather$windspeed) != crs(dtm))  weather$windspeed<-project(weather$windspeed,crs(dtm))
  if (crs(weather$precip) != crs(dtm))  weather$precip<-project(weather$precip,crs(dtm))
  wdir<-apply(.is(weather$winddir),3,.getmode)
  # ============= Run point model for each grid cell of clim array ========= #
  vegpc<-.resamplev(vegp,weather$temp)
  dms<-dim(weather$temp)
  n<-dms[1]*dms[2]
  ij<-ndays/n # for progress bar
  pointms<-list()
  pointms2<-list()
  clims<-list()
  mxhgt<-max(.is(vegp$hgt),na.rm=T)
  if (class(dtmc) == "PackedSpatRaster") dtmc<-rast(dtm)
  ll<-.latslonsfromr(dtmc)
  k<-1
  vegpp<-.sortvegp(vegp,method="P")
  vegpp<-c(vegpp[2],vegpp[1],vegpp[6],vegpp[4])
  for (i in 1:dms[1]) {
    for (j in 1:dms[2]) {
      climdf<-.todf(weather,i,j,tme,wdir)
      if (is.na(climdf$temp[1]) == FALSE & is.na(.is(vegpc$hgt)[i,j]) == FALSE) {
        # ======== Perform wind height adjustment if necessary ============== #
        if (zref != windhgt) {
          climdf$windspeed<-climdf$windspeed*log(67.8*zref-5.42)/log(67.8*windhgt-5.42)
        }
        # =========== Perform height adjustment if necessary ================ #
        obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                            hour=tme$hour+tme$min/60+tme$sec/3600)
        if (mxhgt > zref) {
          tst<-1
          ctr<-0
          while (tst > 0) {
            weather2<-weatherhgtCpp(obstime,climdf,zref,zref,mxhgt,ll$lats[i,j],ll$lons[i,j])
            tt<-mean(weather2$temp)
            if (is.na(tt) == FALSE) tst<-0
            ctr<-ctr+1
            if (ctr > 5) tst<-0
          }
          climdf<-weather2
          climdfr$obs_time<-as.POSIXlt(tme)
          zref<-mxhgt
        } # end height adjust
        # ============ Get variables for pointmodel =========== ============== #
        other<-c(0,0,ll$lats[i,j],ll$lons[i,j],zref,snowinitd,snowinita)
        # ========== Run point model and convert to input data.frame ========= #
        pmod<-pointmodelsnow(obstime,climdf,vegpp,other,snowenv,maxiter=10)
        s<-c(1:length(climdf$temp))+1
        pm<-data.frame(Gp=pmod$G,Tc=pmod$Tc,RswabsG=pmod$RswabsG,RlwabsG=pmod$RlwabsG,
                       umu=pmod$umu,tr=pmod$tr,sdepc=pmod$sdepc[s],sdenc=pmod$sdenc,
                       sdepg=pmod$sdepg[s],sdeng=pmod$sdeng,sublmelt=pmod$sublmelt,
                       tempmelt=pmod$tempmelt,rainmelt=pmod$rainmelt,sstemp=pmod$sstemp)
        pm2<-data.frame(sublmelt=pmod$sublmelt,tempmelt=pmod$tempmelt,rainmelt=pmod$rainmelt,sstemp=pmod$sstemp,
                        tc=climdf$temp,sdenc=pmod$sdenc,sdeng=pmod$sdeng)
        pm2$snow<-ifelse(climdf$temp>2,0,climdf$prec)
        pointms[[k]]<-pm[subs,]
        pointms2[[k]]<-pm2
        clims[[k]]<-climdf[subs,]
        h2<-length(pm2$snow)
      } else {
        pointms[[k]]<-NA
        pointms2[[k]]<-NA
        clims[[k]]<-NA
      }
      utils::setTxtProgressBar(pb, k*ij)
      k<-k+1
    } # end j
  } # end i
  # ================ Turn climate data into arrays =========================== #
  climdata<-list()
  h<-length(subs)
  ii<-0.25*ndays/5
  climdata$temp<-.cca(clims,"temp",h,dtmc,dtm)
  pk<-.cca(clims,"pres",h,dtmc,dtmc)
  climdata$relhum<-.cca(clims,"relhum",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 1*ii+ndays)
  # Apply altitudinal correction
  if (altcorrect == 0) {  # No altitudinal correction
    climdata$pres<-as.array(resample(.rast(pk,dtmc),dtm))
  } else { # Altitudinal correction applied
    es<-.satvap(climdata$temp)
    ea<-es*climdata$relhum/100
    dtmc[is.na(dtmc)]<-0
    psl<-pk/(((293-0.0065*.rta(dtmc,h))/293)^5.26)
    psl<-as.array(resample(.rast(psl,dtmc),dtm))
    climdata$pres<-psl*(((293-0.0065*.rta(dtm,h))/293)^5.26)
    dc<-resample(dtmc,dtm)
    elevd<-.rta(dc-dtm,h)
    if (altcorrect==1) {  # Fixed lapse rate
      tcdif<-elevd*(5/1000)
    } else { # Humidity-dependent lapse rate
      lr<-.lapserate(climdata$temp,ea,climdata$pres)
      tcdif<-lr*elevd
    }
    climdata$temp<-tcdif+climdata$temp
    climdata$relhum<-(ea/.satvap(climdata$temp))*100
  }
  utils::setTxtProgressBar(pb, 2*ii+ndays)
  climdata$relhum[climdata$relhum>100]<-100
  climdata$difrad<-.cca(clims,"difrad",h,dtmc,dtm)
  climdata$swdown<-.cca(clims,"swdown",h,dtmc,dtm)
  climdata$lwdown<-.cca(clims,"lwdown",h,dtmc,dtm)
  climdata$precip<-.cca(clims,"precip",h,dtmc,dtm)
  # Wind speed and direction
  u2<-.cca(clims,"windspeed",h,dtmc,dtmc)
  wd<-.cca(clims,"winddir",h,dtmc,dtmc)
  wu<-u2*cos(wd*pi/180)
  wv<-u2*sin(wd*pi/180)
  wuv<-apply(wu,3,mean,na.rm=TRUE)
  wvv<-apply(wv,3,mean,na.rm=TRUE)
  wu<-as.array(resample(.rast(wu,dtmc),dtm))
  wv<-as.array(resample(.rast(wv,dtmc),dtm))
  climdata$windspeed<-sqrt(wu^2+wv^2)
  climdata$winddir<-(atan2(wvv,wuv)*180/pi)%%360
  # ================ Turn snow point variables into arrays ================== #
  pointm<-list()
  pointm$Gp<-.cca(pointms,"Gp",h,dtmc,dtm)
  pointm$Tc<-.cca(pointms,"Tc",h,dtmc,dtm)
  pointm$RswabsG<-.cca(pointms,"RswabsG",h,dtmc,dtm)
  pointm$RlwabsG<-.cca(pointms,"RlwabsG",h,dtmc,dtm)
  pointm$umu<-.cca(pointms,"umu",h,dtmc,dtm)
  pointm$tr<-.cca(pointms,"tr",h,dtmc,dtm)
  pointm$sdepc<-.cca(pointms,"sdepc",h,dtmc,dtm)
  pointm$sdenc<-.cca(pointms,"sdenc",h,dtmc,dtm)
  pointm$sdeng<-.cca(pointms,"sdeng",h,dtmc,dtm)
  utils::setTxtProgressBar(pb, 3*ii+ndays)
  # ======== Turn snow point variables into arrays: entire time series ====== #
  pointm2<-list()
  pointm2$sublmelt<-.cca(pointms2,"sublmelt",h2,dtmc,dtmc)
  pointm2$tempmelt<-.cca(pointms2,"tempmelt",h2,dtmc,dtmc)
  pointm2$rainmelt<-.cca(pointms2,"rainmelt",h2,dtmc,dtmc)
  pointm2$snow<-.cca(pointms2,"snow",h2,dtmc,dtmc)
  pointm2$sstemp<-.cca(pointms2,"sstemp",h2,dtmc,dtm)
  utils::setTxtProgressBar(pb, 4*ii+ndays)
  pointm2$tc<-.cca(pointms2,"tc",h2,dtmc,dtm)
  pointm2$sdenc<-.cca(pointms2,"sdenc",h2,dtmc,dtmc)
  pointm2$sdeng<-.cca(pointms2,"sdeng",h2,dtmc,dtmc)
  utils::setTxtProgressBar(pb, 5*ii+ndays)
  # ========== Prepare vegetation inputs for grid model =============== #
  sdepc<-.cca(pointms,"sdepc",h,dtmc,dtmc)
  sdept<-apply(sdepc,3,max,na.rm=TRUE)
  vegp<-.sortl(vegp,sdept)
  vegp$leaft[is.na(vegp$leaft)]<-0.01
  ss<-subs[1]-1
  # =========================== subset models ================================ #
  obstime<-obstime[subs,]
  # ========== Prepare other inputs for grid model =========================== #
  ll<-.latslonsfromr(dtm)
  other<-list()
  other$zref<-zref
  other$lats<-ll$lats
  other$lons<-ll$lons
  other$isnowdc<-.is(dtm)*0+snowinitd
  other$isnowac<-.is(dtm)*0+snowinita
  other$isnowag<-.is(dtm)*0+snowinita
  other$slope<-terra::terrain(dtm,v="slope")
  other$slope[is.na(other$slope)]<-0
  other$slope<-.is(mask(other$slope,dtm))
  other$aspect<-terra::terrain(dtm,v="aspect")
  other$aspect[is.na(other$aspect)]<-180
  other$aspect<-.is(mask(other$aspect,dtm))
  other$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) other$hor[,,i]<-.horizon(dtm,(i-1)*15)
  msl<-tan(apply(atan(other$hor),c(1,2),mean))
  other$skyview<-0.5*cos(2*msl)+0.5
  ss<-1
  if (res(dtm)[1]<=100) ss<-10
  other$wsa<-.windsheltera(dtm,zref,ss)
  # ========== Create arrays for storing data ============================== #
  n<-length(subs)
  Tc<-array(NA,dim=c(dim(dtm)[1:2],n))
  Tg<-Tc
  sdepc<-Tc
  sdepg<-array(0,dim=dim(Tc))
  sden<-Tc
  msnow<-mean(pointm2$snow[pointm2$snow>0],na.rm=TRUE)
  mtemp<-mean(pointm2$tc,na.rm=TRUE)
  intfrac<-canintfrac(vegp$hgt,vegp$pai,2,msnow,mtemp,0)
  other$isnowdg<-(1-intfrac)*other$isnowdc
  # ========== Loop through each day ======================================= #
  mup<-.is(dtm)*0+1
  ped<-0
  for (day in 1:ndays) {
    st<-(day-1)*24+1
    ed<-st+23
    s<-c(st:ed)
    # work out snow balance in between current and previous model run
    if ((subs[st]-1) > 1) {
      sbtn<-c((ped+1):(subs[st]-1))
      mu<-meltmu2(other$skyview,pointm2$sstemp[,,sbtn],pointm2$tc[,,sbtn])
      sublmelt<-.resamplemelt(pointm2$sublmelt, sbtn, dtmc, dtm)
      rainmelt<-.resamplemelt(pointm2$rainmelt, sbtn, dtmc, dtm)
      tempmelt<-.resamplemelt(pointm2$tempmelt, sbtn, dtmc, dtm)
      snowsum<-.resamplemelt(pointm2$snow, sbtn, dtmc, dtm)
      melt<-sublmelt+rainmelt+mu*tempmelt
      balancec<-snowsum/1000-melt  # m SWE
      balanceg<-(1-intfrac)*snowsum/1000-exp(-vegp$pai)*melt
      # adjust initial snow depths accordingly
      sdec<-.resamplemelt(pointm2$sdenc, sbtn, dtmc, dtm)/length(sbtn)
      sdeg<-.resamplemelt(pointm2$sdeng, sbtn, dtmc, dtm)/length(sbtn)
      other$isnowdc<-other$isnowdc+balancec*(1000/sdec)
      other$isnowdg<-other$isnowdg+balanceg*(1000/sdeg)
    }
    other$isnowdc[other$isnowdc<0]<-0
    other$isnowdg[other$isnowdg<0]<-0
    # subset weather
    climdata1<-list()
    for (i in 1:8) climdata1[[i]]<-climdata[[i]][,,s]
    climdata1[[9]]<-climdata[[9]][s]
    names(climdata1)<-names(climdata)
    pointm1<-list()
    for (i in 1:9) pointm1[[i]]<-pointm[[i]][,,s]
    names(pointm1)<-names(pointm)
    smod<-gridmodelsnow2(obstime[s,],climdata1,pointm1,vegp,other,snowenv)
    # Calculate change in snow depth
    dsnow<-smod$sdepc-.rta(other$isnowdc,24)  # ground and canopy
    dsnowg<-smod$sdepg-.rta(other$isnowdg,24) # ground only
    dsnowc<-dsnow-dsnowg # canopy only
    # Redistribute snow according to terrain position index
    dtms<-dtm+.rast(sdepg[,,ed],dtm)
    wss<-sqrt(wuv[s]^2+wvv[s]^2)
    tpr<-10*mean(wss)^0.5
    af<-round(tpr/res(dtm)[1],0)
    me<-min(dim(dtm)[1:2])
    tpi<-.tpicalc(af,me,dtms,tfact)
    dsnowg2<-dsnowg*.rta(tpi,length(s))
    dsnowc2<-dsnowc+dsnowg2
    # Assign variables
    Tc[,,s]<-smod$Tc
    Tg[,,s]<-smod$Tg
    sden[,,s]<-smod$sden
    sdc<-dsnowc2+.rta(other$isnowdc,24)
    sdg<-dsnowg2+.rta(other$isnowdg,24)
    sdc[sdc<0]<-0
    sdg[sdg<0]<-0
    sdepc[,,s]<-sdc
    sdepg[,,s]<-sdg
    # Recalculate ped
    pred<-subs[ed]
    utils::setTxtProgressBar(pb, 0.5*day+ndays*1.5)
  }
  # Recalculate as mm snow water equivelent
  swe<-sdepc*sden
  return(list(Tc=Tc,Tg=Tg,groundsnowdepth=sdepg,totalSWE=swe,snowden=sden,umu=pointm$umu))
}
#' Runs microclimate model without snow
.runmicronosnow <- function(micropoint, reqhgt, vegp, soilc, dtm, dtmc = NA, altcorrect = 0,
                     runchecks = TRUE, pai_a = NA, tfact = 1.5,
                     out = rep(TRUE, 10), slr = NA, apr = NA, hor = NA,
                     twi = NA, wsa = NA, svf = NA, method = "Cpp") {
  if (method == "R") {
    micro<-modelin(micropoint,vegp,soilc,dtm,dtmc,altcorrect,runchecks)
    if (reqhgt == 0) {
      micro<-soiltemp(micro,reqhgt,pai_a,tfact)
      Tz<-micro$Tg
      mout<-aboveground(micro,reqhgt,pai_a,tfact)
      mout$Tz<-Tz
      if (out[1] == FALSE) mout$Tz<-NULL
      mout$tleaf<-NULL
      mout$relhum<-NULL
      if (out[4] == FALSE) mout$soilm<-NULL
      mout$windspeed<-NULL
      if (out[6] == FALSE) mout$Rdirdown<-NULL
      if (out[7] == FALSE) mout$Rdifdown<-NULL
      if (out[8] == FALSE) mout$Rlwdown<-NULL
      if (out[9] == FALSE) mout$Rswup<-NULL
      if (out[10] == FALSE) mout$Rlwup<-NULL
    }
    if (reqhgt > 0) {
      mout<-aboveground(micro,reqhgt,pai_a,tfact)
      if (out[1] == FALSE) mout$Tz<-NULL
      if (out[2] == FALSE) mout$tleaf<-NULL
      if (out[3] == FALSE) mout$relhum<-NULL
      if (out[4] == FALSE) mout$soilm<-NULL
      if (out[5] == FALSE) mout$windspeed<-NULL
      if (out[6] == FALSE) mout$Rdirdown<-NULL
      if (out[7] == FALSE) mout$Rdifdown<-NULL
      if (out[8] == FALSE) mout$Rlwdown<-NULL
      if (out[9] == FALSE) mout$Rswup<-NULL
      if (out[10] == FALSE) mout$Rlwup<-NULL
    }
    if (reqhgt < 0) {
      mout<-belowground(micro,reqhgt,pai_a,tfact)
      mout$T0<-NULL
      if (out[1] == FALSE) mout$Tz<-NULL
      if (out[4] == FALSE) mout$soilm<-NULL
    }
  } else {
    up<-.unpack(dtm,vegp,soilc)
    dmx<-.vegpdmx(up$vegp)
    if (dmx == 1) {  # temporally invarient vegetation
      if (class(micropoint) == "micropoint") { # data.frame climate input
        mout<-.runmodel1Cpp(micropoint,vegp,soilc,dtm,reqhgt,runchecks,pai_a,tfact,out,slr,apr,hor,twi,wsa)
      } else {  # array climate input
        mout<-.runmodel2Cpp(micropoint,vegp,soilc,dtm,dtmc,reqhgt,runchecks,altcorrect,pai_a,tfact,out,slr,apr,hor,twi,wsa)
      }
    } else { # time variant vegetation
      if (class(micropoint) == "micropoint") { # data.frame climate input
        mout<-.runmodel3Cpp(micropoint,vegp,soilc,dtm,reqhgt,runchecks,pai_a,tfact,out,slr,apr,hor,twi,wsa)
      } else { # array climate input
        mout<-.runmodel4Cpp(micropoint,vegp,soilc,dtm,dtmc,reqhgt,runchecks,altcorrect,pai_a,tfact,out,slr,apr,hor,twi,wsa)
      } # end if array
    }  # end if time variant
  } # end if R/C++
  return(mout)
}
#' prepare inputs for snow microclimate model (data.frame climate)
.prepsnowinputs1<-function(reqhgt,dtm,vegp,soilc,micropoints,snowdays,nosnowdays,smod,moutn,tmeorig,subs,runchecks,slr=NA,apr=NA,hor=NA,svf=NA,wsa=NA,paia) {
  # ================= unpack and run checks ================================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  cv<-.cleanvars(vegp,soilc,dtm) # sorts out NAs and zeros
  dtm<-cv$dtm
  vegp<-cv$vegp
  soilc<-cv$soilc
  weather<-micropoints$weather
  if (runchecks) {
    rc<-checkinputs(weather,vegp,soilc,dtm)
    weather<-rc$weather
    vegp<-rc$vegp
    soilc<-rc$soilc
  }
  # =========================  Create subset for snow model ================== #
  ai<-rep((snowdays-1)*24,each=24)+rep(c(1:24),length(snowdays))
  # =========================  Prep obstime ================================== #
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                      hour=tme$hour+tme$min/60+tme$sec/3600)
  # ======================== Prepare climate data ============================ #
  weather$obs_time<-NULL
  weather$umu<-smod$umu[ai]
  # ================= Prepare vegetation inputs for grid model =============== #
  sdept<-rep(0,length(tmeorig))
  sdept[micropoints$subs]<-1
  vegp<-.sortl2(vegp,sdept,reqhgt,paia)
  # ====================== Prepare other inputs for grid model =============== #
  other<-list()
  if (class(slr)[1] == "logical") {
    other$slope<-.is(terrain(dtm, v="slope"))
  } else other$slope<-.is(slr)
  if (class(apr)[1] == "logical") {
    other$aspect<-.is(terrain(dtm, v="aspect"))
  }  else other$aspect<-.is(apr)
  ll<-.latlongfromraster(dtm)
  if (class(hor)[1] == "logical") {
    other$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
    for (i in 1:24) other$hor[,,i]<-.horizon(dtm,(i-1)*15)
  } else other$hor<-hor
  if (class(svf)  == "logical") {
    msl<-tan(apply(atan(other$hor),c(1,2),mean))
    other$skyview<-0.5*cos(2*msl)+0.5
  } else other$skyview<-.is(svf)
  s<-1
  if (res(dtm)[1]<=100) s<-10
  if (class(wsa)[1] == "logical") {
    other$wsa<-.windsheltera(dtm,micropoint$zref,s)
  } else other$wsa<-wsa
  other$lat<-ll$lat
  other$lon<-ll$lon
  other$zref<-micropoints$zref
  soilp<-.sortsoilc(soilc,method="R")
  other$Smax<-soilp$Smax
  # ========================== Prepare microinput ============================ #
  t1<-length(snowdays)*24 # number of timesteps for days with snow
  # Days common to both datasets
  s<-c(1:t1)
  s1<-s[rep(snowdays %in%  nosnowdays,each=24)]  # entries in blank micro that need to be replaced by micron outputs
  s<-c(1:(length(nosnowdays)*24))
  s2<-s[rep(nosnowdays %in%  snowdays,each=24)]  # entries in micron that are being swapped in (length(s1) == length(s2))
  # Create a blank dataset
  a<-array(NA,dim=c(dim(dtm)[1:2],t1))
  nms<-names(moutn)
  micros<-list()
  for (i in 1:length(nms)) {
    v<-a
    v2<-moutn[[i]]
    v[,,s1]<-v2[,,s2]
    micros[[i]]<-v
  }
  names(micros)<-nms
  return(list(obstime=obstime,weather=weather,micro=micros,vegp=vegp,other=other))
}
#' prepare inputs for snow microclimate model (array climate)
.prepsnowinputs2<-function(reqhgt,dtm,dtmc,vegp,soilc,micropoints,altcorrect,runchecks,snowdays,
                          nosnowdays,moutn,slr=NA,apr=NA,hor=NA,svf=NA,wsa=NA,paia) {
  # ==================== unpack variables ==================================== #
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  weather<-list()
  # ==================== Extract weather and point data ======================= #
  weather<-list()
  dfo<-list()
  for (i in 1:length(micropoints)) {
    onepoint<-micropoints[[i]]
    if (class(onepoint) != "logical") {
      weather[[i]]<-onepoint$weather
      dfo[[i]]<-onepoint$dfo
      tme<-as.POSIXlt(weather[[i]]$obs_time,tz="UTC")
      tmeorig<-onepoint$tmeorig
      subs<-onepoint$subs
      zref=onepoint$zref
      if (runchecks) {
        rc<-checkinputs(weather[[i]],vegp,soilc,dtm,onepoint$zref)
        weather[[i]]<-rc$weather
        vegp<-rc$vegp
        soilc<-rc$soilc
      }
    } else {
      weather[[i]]<-NA
      dfo[[i]]<-NA
    }
  }
  if (class(dtmc)[1] == "PackedSpatRaster") dtmc<-rast(dtmc)
  # =========================  Prep obstime ================================== #
  obstime<-data.frame(year=tme$year+1900,month=tme$mon+1,day=tme$mday,
                      hour=tme$hour+tme$min/60+tme$sec/3600)
  # ================== Turn climate variables into arrays ============= #
  climdata<-list()
  h<-length(tme)
  tc<-.cca(weather,"temp",h,dtmc,dtmc)
  pk<-.cca(weather,"pres",h,dtmc,dtmc)
  rh<-.cca(weather,"relhum",h,dtmc,dtmc)
  ea<-.is(resample(.rast(.satvap(tc)*(rh/100),dtmc),dtm))
  climdata$temp<-as.array(resample(.rast(tc,dtmc),dtm))
  # Apply altitudinal correction
  if (altcorrect == 0) {  # No altitudinal correction
    climdata$pres<-as.array(resample(.rast(pk,dtmc),dtm))
  } else { # Altitudinal correction applied
    dtmc[is.na(dtmc)]<-0
    psl<-pk/(((293-0.0065*.rta(dtmc,h))/293)^5.26)
    psl<-as.array(resample(.rast(psl,dtmc),dtm))
    climdata$pres<-psl*(((293-0.0065*.rta(dtm,h))/293)^5.26)
    dc<-resample(dtmc,dtm)
    elevd<-.rta(dc-dtm,h)
    if (altcorrect==1) {  # Fixed lapse rate
      tcdif<-elevd*(5/1000)
    } else { # Humidity-dependent lapse rate
      lr<-.lapserate(climdata$temp,ea,climdata$pres)
      tcdif<-lr*elevd
    }
    climdata$temp<-tcdif+climdata$temp
  }
  climdata$relhum<-(ea/.satvap(climdata$temp))*100
  climdata$relhum[climdata$relhum>100]<-100
  climdata$relhum[climdata$relhum<20]<-20
  climdata$swdown<-.cca(weather,"swdown",h,dtmc,dtm)
  climdata$difrad<-.cca(weather,"difrad",h,dtmc,dtm)
  climdata$lwdown<-.cca(weather,"lwdown",h,dtmc,dtm)
  # Wind speed and direction
  u2<-.cca(weather,"windspeed",h,dtmc,dtmc)
  wd<-.cca(weather,"winddir",h,dtmc,dtmc)
  wu<-u2*cos(wd*pi/180)
  wv<-u2*sin(wd*pi/180)
  wuv<-apply(wu,3,mean,na.rm=TRUE)
  wvv<-apply(wv,3,mean,na.rm=TRUE)
  wu<-as.array(resample(.rast(wu,dtmc),dtm))
  wv<-as.array(resample(.rast(wv,dtmc),dtm))
  climdata$windspeed<-sqrt(wu^2+wv^2)
  climdata$winddir<-(atan2(wvv,wuv)*180/pi)%%360
  # precipitation and umu
  climdata$prec<-.cca(weather,"precip",h,dtmc,dtm)
  climdata$umu<-.cca(dfo,"umu",h,dtmc,dtm)
  # =========================  Create subset for snow model ================== #
  ai<-rep((snowdays-1)*24,each=24)+rep(c(1:24),length(snowdays))
  # ================= Prepare vegetation inputs for grid model =============== #
  sdept<-rep(0,length(tmeorig))
  sdept[micropoints$subs]<-1
  vegp<-.sortl2(vegp,sdept,reqhgt,paia)
  # ====================== Prepare other inputs for grid model =============== #
  other<-list()
  if (class(slr)[1] == "logical") {
    other$slope<-.is(terrain(dtm, v="slope"))
  } else other$slope<-.is(slr)
  if (class(apr)[1] == "logical") {
    other$aspect<-.is(terrain(dtm, v="aspect"))
  }  else other$aspect<-.is(apr)
  ll<-.latlongfromraster(dtm)
  if (class(hor)[1] == "logical") {
    other$hor<-array(NA,dim=c(dim(dtm)[1:2],24))
    for (i in 1:24) other$hor[,,i]<-.horizon(dtm,(i-1)*15)
  } else other$hor<-hor
  if (class(svf)  == "logical") {
    msl<-tan(apply(atan(other$hor),c(1,2),mean))
    other$skyview<-0.5*cos(2*msl)+0.5
  } else other$skyview<-.is(svf)
  s<-1
  if (res(dtm)[1]<=100) s<-10
  if (class(wsa)[1] == "logical") {
    other$wsa<-.windsheltera(dtm,zref,s)
  } else other$wsa<-wsa
  ll<-.latslonsfromr(dtm)
  other$lat<-ll$lats
  other$lon<-ll$lons
  other$zref<-zref
  soilp<-.sortsoilc(soilc,method="R")
  other$Smax<-soilp$Smax
  # ========================== Prepare microinput ============================ #
  t1<-length(snowdays)*24 # number of timesteps for days with snow
  # Days common to both datasets
  s<-c(1:t1)
  s1<-s[rep(snowdays %in%  nosnowdays,each=24)]  # entries in blank micro that need to be replaced by micron outputs
  s<-c(1:(length(nosnowdays)*24))
  s2<-s[rep(nosnowdays %in%  snowdays,each=24)]  # entries in micron that are being swapped in (length(s1) == length(s2))
  # Create a blank dataset
  a<-array(NA,dim=c(dim(dtm)[1:2],t1))
  nms<-names(moutn)
  micros<-list()
  for (i in 1:length(nms)) {
    v<-a
    v2<-moutn[[i]]
    v[,,s1]<-v2[,,s2]
    micros[[i]]<-v
  }
  names(micros)<-nms
  return(list(obstime=obstime,weather=climdata,micro=micros,vegp=vegp,other=other))
}
#' Run microclimate grid model with snow (data.frame weather)
.runmicrosnow1 <- function(micropoint,reqhgt,vegp,soilc,dtm,smod,runchecks=TRUE,pai_a=NA,tfact=1.5,
                           out=rep(TRUE,10),slr=NA,apr=NA,hor=NA,twi=NA,wsa=NA,svf=NA)  {
  pb <- utils::txtProgressBar(min = 0, max = 5, style = 3)
  if (class(dtm) == "PackedSpatRaster") dtm<-rast(dtm)
  # (1) Figure out snow and no snow days
  minsnow<-applycpp3(smod$totalSWE,"min")
  maxsnow<-applycpp3(smod$totalSWE,"max")
  snowd<-snowdaysfun(maxsnow, minsnow)
  v<-c(1:length(snowd$snowdays))
  snowdays<-v[snowd$snowdays==1]
  nosnowdays<-v[snowd$nosnowdays==1]
  # (2) Subset micropoint for snow and no snow days
  micropoints<-subsetpointmodel(micropoint,days=snowdays)
  micropointn<-subsetpointmodel(micropoint,days=nosnowdays)
  utils::setTxtProgressBar(pb,1)
  # (3) Run no snow microclimate model for no snow days
  moutn<-.runmicronosnow(micropointn,reqhgt,vegp,soilc,dtm,dtmc,altcorrect,runchecks,
                         pai_a,tfact,out,slr,apr,hor,twi,wsa,svf)
  utils::setTxtProgressBar(pb,3)
  # (4) Run snow microclimate model for snow days
  tmeorig<-micropoint$tmeorig
  subs<-micropoints$subs
  snowin<-.prepsnowinputs1(reqhgt,dtm,vegp,soilc,micropoints,snowdays,nosnowdays,
                           smod,moutn,tmeorig,subs,runchecks,slr,apr,hor,svf,wsa,pai_a)
  ss<-rep((snowdays-1)*24,each=24)+rep(c(1:24),length(snowdays))
  smods<-subsetsnowmodel(smod, ss)
  mouts<-gridmicrosnow1(reqhgt,snowin$obstime,snowin$weather,smods,snowin$micro,snowin$vegp,snowin$other,out)
  utils::setTxtProgressBar(pb,4)
  # (5) Merge the two model outputs back together to produce one output
  # (5a) - calculate days without any snow
  tdays<-unique(c(snowdays,nosnowdays))
  tdays<-tdays[order(tdays)]
  nosnow <- setdiff(tdays, snowdays)
  # (5b) get subset for no snow microclimate model to include only those days with no snow at all
  s<-c(1:dim(moutn[[1]])[3])
  s1<-s[rep(nosnowdays %in%  nosnow, each=24)]
  # (5c) produce list of entries for each dataset
  nosnowh<-rep((nosnow-1)*24,each=24)+rep(c(1:24),length(nosnow))
  snowh<-rep((snowdays-1)*24,each=24)+rep(c(1:24),length(snowdays))
  # (5c) merge the two datasets
  mout<-list()
  n<-length(nosnowh)+length(snowh)
  for (i in 1:length(moutn))  {
    # Create blank entry for storin g data
    mout[[i]]<-array(NA,dim=c(dim(dtm)[1:2],n))
    # perform subset of nosnow
    xx<-moutn[[i]][,,s1]
    mout[[i]][,,nosnowh]<-xx
    mout[[i]][,,snowh]<-mouts[[i]]
  }
  names(mout)<-names(moutn)
  utils::setTxtProgressBar(pb,5)
  return(mout)
}
#' Run microclimate grid model with snow (array weather)
.runmicrosnow2 <- function(micropoint,reqhgt,vegp,soilc,dtm,dtmc,smod,altcorrect = 0,runchecks=TRUE,pai_a=NA,tfact=1.5,
                           out=rep(TRUE,10),slr=NA,apr=NA,hor=NA,twi=NA,wsa=NA,svf=NA)  {
  pb <- utils::txtProgressBar(min = 0, max = 6, style = 3)
  if (class(dtm) == "PackedSpatRaster") dtm<-rast(dtm)
  # (1) Figure out snow and no snow days
  minsnow<-applycpp3(smod$totalSWE,"min")
  maxsnow<-applycpp3(smod$totalSWE,"max")
  snowd<-snowdaysfun(maxsnow, minsnow)
  v<-c(1:length(snowd$snowdays))
  snowdays<-v[snowd$snowdays==1]
  nosnowdays<-v[snowd$nosnowdays==1]
  # (2) Subset micropoint for snow and no snow days
  micropoints<-subsetpointmodela(micropoint,days=snowdays)
  micropointn<-subsetpointmodela(micropoint,days=nosnowdays)
  utils::setTxtProgressBar(pb,1)
  # (3) Run no snow microclimate model for no snow days
  moutn<-.runmicronosnow(micropointn,reqhgt,vegp,soilc,dtm,dtmc,altcorrect,runchecks,
                         pai_a,tfact,out,slr,apr,hor,twi,wsa,svf)
  utils::setTxtProgressBar(pb,3)
  # (4) Run snow microclimate model for snow days
  snowin<-.prepsnowinputs2(reqhgt,dtm,dtmc,vegp,soilc,micropoints,altcorrect,runchecks,snowdays,
                           nosnowdays,moutn,slr,apr,hor,svf,wsa,pai_a)
  utils::setTxtProgressBar(pb,4)
  ss<-rep((snowdays-1)*24,each=24)+rep(c(1:24),length(snowdays))
  smods<-subsetsnowmodel(smod, ss)
  mouts<-gridmicrosnow2(reqhgt,snowin$obstime,snowin$weather,smods,snowin$micro,snowin$vegp,snowin$other,out)
  utils::setTxtProgressBar(pb,5)
  # (5) Merge the two model outputs back together to produce one output
  # (5a) - calculate days without any snow
  tdays<-unique(c(snowdays,nosnowdays))
  tdays<-tdays[order(tdays)]
  nosnow <- setdiff(tdays, snowdays)
  # (5b) get subset for no snow microclimate model to include only those days with no snow at all
  s<-c(1:dim(moutn[[1]])[3])
  s1<-s[rep(nosnowdays %in%  nosnow, each=24)]
  # (5c) produce list of entries for each dataset
  nosnowh<-rep((nosnow-1)*24,each=24)+rep(c(1:24),length(nosnow))
  snowh<-rep((snowdays-1)*24,each=24)+rep(c(1:24),length(snowdays))
  # (5c) merge the two datasets
  mout<-list()
  n<-length(nosnowh)+length(snowh)
  for (i in 1:length(moutn))  {
    # Create blank entry for storin g data
    mout[[i]]<-array(NA,dim=c(dim(dtm)[1:2],n))
    # perform subset of nosnow
    xx<-moutn[[i]][,,s1]
    mout[[i]][,,nosnowh]<-xx
    mout[[i]][,,snowh]<-mouts[[i]]
  }
  names(mout)<-names(moutn)
  utils::setTxtProgressBar(pb,6)
  return(mout)
}
#' blend two adjacent rasters that have overlap
.blendmosaic<-function(r1, r2) {
  # run checks
  reso1<-res(r1)
  reso2<-res(r2)
  if (reso1[1] != reso2[1]) stop("resolutions must match")
  if (reso1[2] != reso2[2]) stop("resolutions must match")
  # Find whether r2 is TT, TR, RR, BR, BB, BL, LL or TL
  e1<-ext(r1)
  e2<-ext(r2)
  corner<-"ID"
  if (e2$ymax > e1$ymax & e2$xmax == e1$xmax) corner<-"TT"
  if (e2$ymax > e1$ymax & e2$xmax > e1$xmax) corner<-"TR"
  if (e2$ymax == e1$ymax & e2$xmax > e1$xmax) corner<-"RR"
  if (e2$ymax < e1$ymax & e2$xmax > e1$xmax) corner<-"BR"
  if (e2$ymax < e1$ymax & e2$xmax == e1$xmax) corner<-"BB"
  if (e2$ymax < e1$ymax & e2$xmax < e1$xmax) corner<-"BL"
  if (e2$ymax == e1$ymax & e2$xmax < e1$xmax) corner<-"LL"
  if (e2$ymax > e1$ymax & e2$xmax < e1$xmax) corner<-"TL"
  if (corner == "ID") {
    ro<-mosaic(r1,r2,fun="mean")
  } else {
    # Calculate overlap area
    if (corner == "TR" || corner == "TT") eo<-ext(e2$xmin,e1$xmax,e2$ymin,e1$ymax)
    if (corner == "BR" || corner == "RR") eo<-ext(e2$xmin,e1$xmax,e1$ymin,e2$ymax)
    if (corner == "BL" || corner == "BB") eo<-ext(e1$xmin,e2$xmax,e1$ymin,e2$ymax)
    if (corner == "TL" || corner == "LL") eo<-ext(e2$xmin,e1$xmax,e2$ymin,e1$ymax)
    # Calculate weights
    nx<-as.numeric((eo$xmax-eo$xmin)/res(r1)[1])
    ny<-as.numeric((eo$ymax-eo$ymin)/res(r1)[2])
    wx<-matrix(rep(seq(0,1,length.out=nx),each=ny),ncol=nx,nrow=ny)
    wy<-matrix(rep(seq(0,1,length.out=ny),nx),ncol=nx,nrow=ny)
    # Create a SpatRast of the blended area
    if (corner == "TT") w1<-wy
    if (corner == "TR") w1<-sqrt((1-wx)^2+wy^2)/sqrt(2)
    if (corner == "RR") w1<-1-wx
    if (corner == "BR") w1<-sqrt((1-wx)^2+(1-wy)^2)/sqrt(2)
    if (corner == "BB") w1<-1-wy
    if (corner == "BL") w1<-sqrt(wx^2+(1-wy)^2)/sqrt(2)
    if (corner == "LL") w1<-wx
    if (corner == "TL") w1<-sqrt(wx^2+wy^2)/sqrt(2)
    # Apply weights to raster
    nn<-dim(r1)[3]
    r1c<-crop(r1,eo)
    r2c<-crop(r2,eo)
    w1<-.rast(.rta(w1,nn),r1c[[1]])
    rb<-r1c*w1+r2c*(1-w1)
    # Clip out the overlap area
    ro<-mosaic(r1,r2)
    re<-w1*0-9999
    ro<-mosaic(ro,re,fun="min")
    ro[ro == -9999]<-NA
    # Mosaic with belnded data
    ro<-mosaic(ro,rb,fun="mean")
  }
  return(ro)
}
