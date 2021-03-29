#' generate sequence of soil data
.soilgen<-function(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr) {
  s<-soilm[1]
  s2<-0.2
  Smin<-min(soilm)
  Smax<-max(soilm)
  for (i in 2:length(rnetp)) {
    s[i]<-s[i-1]+rmu*prec[i]-mult*rnetp[i]
    sav<-(s[i-1]+s2[i-1])/2
    k<-ks*(sav/smx2)^pwr
    dif<-s2[i-1]-s[i-1]
    s[i]<-s[i]+a*k*dif
    s2[i]<-s2[i-1]-((a*k*dif)/10)
    s[i]<-ifelse(s[i]>Smax,Smax,s[i])
    s[i]<-ifelse(s[i]<Smin,Smin,s[i])
    s2[i]<-ifelse(s2[i]>smx2,smx2,s2[i])
    s2[i]<-ifelse(s2[i]<Smin,Smin,s2[i])
  }
  rms<-sqrt(sum((s-soilm)^2)/length(s))
  return(list(s=s,s2=s2,rms=rms))
}
#' fit radiation parameters
.fitrad <- function(prec,rnetp,rmu,a,soilm,ks,smx2,pwr) {
  rms<-0
  for (i in 1:20) {
    m<-i/1e4
    rms[i]<-.soilgen(prec,rnetp,m,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd1<-(sel[1]-1)/1e4
  for (i in 1:20) {
    m<-(i/1e5)+toadd1
    rms[i]<-.soilgen(prec,rnetp,m,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd2<-(sel[1]-1)/1e5+toadd1
  for (i in 1:20) {
    m<-(i/1e6)+toadd2
    rms[i]<-.soilgen(prec,rnetp,m,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd3<-(sel[1]-1)/1e6+toadd2
  for (i in 1:20) {
    m<-(i/1e7)+toadd3
    rms[i]<-.soilgen(prec,rnetp,m,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd4<-(sel[1]-1)/1e7+toadd3
  for (i in 1:20) {
    m<-(i/1e8)+toadd4
    rms[i]<-.soilgen(prec,rnetp,m,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd5<-(sel[1]-1)/1e8+toadd4
  for (i in 1:20) {
    m<-(i/1e9)+toadd5
    rms[i]<-.soilgen(prec,rnetp,m,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  mult<-sel[1]/1e9+toadd5
  mult
}
#' fit rain parameter
.fitrain <- function(prec,rnetp,mult,a,soilm,ks,smx2,pwr) {
  rms<-0
  for (i in 1:20) {
    rmu<-i
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd0<-(sel[1]-1)
  for (i in 1:20) {
    rmu<-(i/10)+toadd0
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd1<-(sel[1]-1)/10+toadd0
  for (i in 1:20) {
    rmu<-(i/100)+toadd1
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd2<-(sel[1]-1)/100+toadd1
  for (i in 1:20) {
    rmu<-(i/1000)+toadd2
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd3<-(sel[1]-1)/1000+toadd2
  for (i in 1:20) {
    rmu<-(i/10000)+toadd3
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd4<-(sel[1]-1)/10000+toadd3
  for (i in 1:20) {
    rmu<-(i/100000)+toadd4
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd5<-(sel[1]-1)/100000+toadd4
  for (i in 1:20) {
    rmu<-(i/1000000)+toadd5
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  rd<-sel[1]/1000000+toadd5
  rd
}
#' fit water exchange constant
.fita <- function(prec,rnetp,mult,rmu,soilm,ks,smx2,pwr) {
  rms<-0
  for (i in 1:20) {
    a<-i
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd0<-(sel[1]-1)
  for (i in 1:20) {
    a<-(i/10)+toadd0
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd1<-(sel[1]-1)/10+toadd0
  for (i in 1:20) {
    a<-(i/100)+toadd1
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd2<-(sel[1]-1)/100+toadd1
  for (i in 1:20) {
    a<-(i/1000)+toadd2
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd3<-(sel[1]-1)/1000+toadd2
  for (i in 1:20) {
    a<-(i/10000)+toadd3
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd4<-(sel[1]-1)/10000+toadd3
  for (i in 1:20) {
    a<-(i/100000)+toadd4
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd5<-(sel[1]-1)/100000+toadd4
  for (i in 1:20) {
    a<-(i/1000000)+toadd5
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  a<-sel[1]/1000000+toadd5
  a
}
#' fit water exchange power
.fitpwr <- function(prec,rnetp,mult,rmu,a,soilm,ks,smx2) {
  rms<-0
  for (i in 1:20) {
    pwr<-i
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd0<-(sel[1]-1)
  for (i in 1:20) {
    pwr<-(i/10)+toadd0
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd1<-(sel[1]-1)/10+toadd0
  for (i in 1:20) {
    pwr<-(i/100)+toadd1
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd2<-(sel[1]-1)/100+toadd1
  for (i in 1:20) {
    pwr<-(i/1000)+toadd2
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd3<-(sel[1]-1)/1000+toadd2
  for (i in 1:20) {
    pwr<-(i/10000)+toadd3
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd4<-(sel[1]-1)/10000+toadd3
  for (i in 1:20) {
    pwr<-(i/100000)+toadd4
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  toadd5<-(sel[1]-1)/100000+toadd4
  for (i in 1:20) {
    pwr<-(i/1000000)+toadd5
    rms[i]<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)$rms
  }
  sel<-which(rms==min(rms))
  pwr<-sel[1]/1000000+toadd5
  pwr
}
.fitall <- function(prec,rnetp,out,soilm,ii) {
  ks<-soilparameters$Ksat[ii]
  smx2<-soilparameters$Smax[ii]
  mult<-out$mult
  rmu<-out$rmu
  a<-out$a
  pwr<-out$pwr
  mult<-.fitrad(prec,rnetp,rmu,a,soilm,ks,smx2,pwr)
  rmu<-.fitrain(prec,rnetp,mult,a,soilm,ks,smx2,pwr)
  a<-.fita(prec,rnetp,mult,rmu,soilm,ks,smx2,pwr)
  pwr<-.fitpwr(prec,rnetp,mult,rmu,a,soilm,ks,smx2)
  sg<-.soilgen(prec,rnetp,mult,rmu,a,soilm,ks,smx2,pwr)
  return(list(mult=mult,rmu=rmu,a=a,pwr=pwr,rms=sg$rms))
}
#' Function for deriving parameters for soil moisture model
#'
#' The function `fitsoilm` derives parameters for fitting the soil moisture model
#' where daily moisture data are available
#'
#' @param weather a data.frame of weather variables (see details).
#' @param rainfall a vector of daily rainfall.
#' @param soilm a vector of daily soil moisture values.
#' @param soiltype a soil type as listed in [soilparameters()].
#' @param plotout optional logical indicating whether to plot fit on eahc iteration.
#' @param max.iter maximum number of iterations used for fitting the model.
#' @return a list wiht the following components:
#' @return `mult` multiplier for net radiation.
#' @return `rmu`  multiplier for rainfall.
#' @return `a` constant controlling sub-surface and surface water exchange.
#' @return `pwr` exponent controlling sub-surface and surface water exchange.
#' @return `rms` RMS Error of model fit.
#' @export
#'
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. In this version of the model, all precipitation is assumed to
#' fall as rain.
#'
#' @examples
#' require(NicheMapR) # - see https://github.com/mrke/NicheMapR
#' require(microclimc) # - see https://github.com/ilyamaclean/microclimc
#' # Run NicheMapR complex multilayer soil model
#' microout<-runNMR(climdata, rainfall, 50.2178, -5.32656, 0.05, 0.01, PAI = 0, soiltype="Clay loam")
#' # Extract surface soil moisture and covert to daily
#' sm<-microout$soilmoist
#' soilm<-sm[,3]
#' soilm<-matrix(soilm,ncol=24,byrow=T)
#' soilm<-apply(soilm,1,mean)
#' # Fit simple model
#' soilmcoefs<-fitsoilm(weather, rainfall, soilm, "Clay")
#' soilmcoefs
fitsoilm <- function(weather, rainfall, soilm, soiltype, plotout = TRUE, max.iter=20) {
  # Calculate net radiation
  swrad<-(1-0.15)*weather$swrad
  lwout<-5.67*10^-8*0.95*(climdata$temp+273.15)^4
  lwnet<-(1-climdata$skyem)*lwout
  rnet<-swrad-lwnet
  rnetp<-rnet
  rnetp[rnetp<0]<-0
  rnetp<-matrix(rnetp,ncol=24,byrow=T)
  rnetp<-apply(rnetp,1,mean)
  # soil types
  ii<-which(soilparameters$Soil.type==soiltype)
  out<-list(mult=1,rmu=0.04,a=0,pwr=1,rms=100)
  rmsa<-3
  tst<-1
  i<-2
  while (tst == 1) {
    out<-.fitall(rainfall,rnetp,out,soilm,ii)
    rmsa[i]<-out$rms
    dif<-rmsa[i-1]-rmsa[i]
    tst<-ifelse(dif<0.00005,0,1)
    tst<-ifelse(i>max.iter,0,tst)
    if (plotout) {
      ymx<-soilparameters$Smax[ii]
      sg<-.soilgen(rainfall,rnetp,out$mult,out$rmu,out$a,soilm,
                   soilparameters$Ksat[ii],ymx,out$pwr)
      plot(sg$s,type="l",ylim=c(0,(ymx+0.1)),xlab="Day of year",ylab="Soil moisture fraction")
      par(new=T)
      plot(soilm,type="l",ylim=c(0,(ymx+0.1)),col="red",xlab="",ylab="",main=paste("iter:",i-1))
    }
    i<-i+1
  }
  if (plotout) {
    par(mar=c(5,5,1,1))
    ymx<-max(sg$s,soilm)
    ymx<-ceiling(ymx*10)/10
    plot(sg$s,type="l",ylim=c(0,(ymx)),xlab="Day of year",ylab="Soil moisture fraction",cex.axis=2,cex.lab=2)
    par(new=T)
    plot(soilm,type="l",ylim=c(0,(ymx)),col="red",xlab="",ylab="",cex.axis=2,cex.lab=2)

  }
  return(out)
}
#' Function for deriving parameters for soil temperature model
#'
#' The function `fitsoilt` derives parameters for fitting the soil temperature model
#' where hourly temperature data are available
#'
#' @param weather a data.frame of weather variables (see details).
#' @param rainfall a vector of daily rainfall.
#' @param soiltemp a vector of hourly soil temperatures
#' @param soiltype a soil type as listed in [soilparameters()].
#' @param theta optional vector of hourly soil moistures. If NA estimated using [soilmpredict()]
#' @param plotout optional logical indicating whether to plot fit on eahc iteration.
#' @param soilmcoefs an optional list of soil moisture model coefficients as returned by [fitsoilm()] (see details)
#' @return a list with the following components:
#' @return `int` model intercept
#' @return `t1`  multiplier for net radiation.
#' @return `t2` multiplier for soil moisture.
#' @return `t3` multiplier for interaction term.
#' @return `rms` RMS Error of model fit.
#'
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. If `soilmcoefs` is not supplied, the inbuilt parameters in
#' [soilparameters()] are used to derive soil moisture.
#'
#' @examples
#' require(NicheMapR) # - see https://github.com/mrke/NicheMapR
#' require(microclimc) # - see https://github.com/ilyamaclean/microclimc
#' # Run NicheMapR complex multilayer soil model
#' microout<-runNMR(climdata, rainfall, 50.2178, -5.32656, 0.05, 0.01, PAI = 0, soiltype="Clay loam")
#' # Extract surface soil temperature
#' st<-microout$soiltemps
#' soiltemp<-st$D0cm
#' # Fit simple model
#' soiltcoefs<-fitsoilt(climdata, rainfall, soiltemp, "Clay loam")
#' soiltcoefs
fitsoilt <- function(weather, rainfall, soiltemp, soiltype, theta = NA, plotout = TRUE, soilmcoefs = NA) {
  reftemp<-weather$temp
  # Calculate net radiation
  swrad<-(1-0.15)*weather$swrad
  lwout<-5.67*10^-8*0.95*(climdata$temp+273.15)^4
  lwnet<-(1-climdata$skyem)*lwout
  rnet<-swrad-lwnet
  # Get soil moistures
  if (is.na(theta[1])) {
    theta<-soilmpredict(rainfall, rnet, soiltype)$soilm1
    theta[theta<0.0002]<-0.0002
    theta<-rep(theta,each=24)
  }
  sm<-log(theta/(1-theta))+8.516993
  # Get parameters
  Tanom<-soiltemp-reftemp
  m1<-summary(lm(Tanom~sm*rnet))
  pred<-m1$coef[1,1]+m1$coef[3,1]*rnet+m1$coef[2,1]*sm+m1$coef[4,1]*rnet*sm
  pred<-pred+reftemp
  ymn<-floor(min(pred,soiltemp))
  ymx<-ceiling(max(pred,soiltemp))
  tme<-as.POSIXct(weather$obs_time)
  if (plotout) {
    par(mar=c(5,5,1,1))
    plot(soiltemp~tme,type="l",col=rgb(1,0,0,0.5),ylim=c(ymn,ymx),xlab="Month",ylab="Temperature",cex.axis=2,cex.lab=2)
    par(new=T)
    plot(pred~tme,type="l",col=rgb(0,0,0,0.5),ylim=c(ymn,ymx),xlab="",ylab="",cex.axis=2,cex.lab=2)
  }
  rms<-sqrt(sum((pred-soiltemp)^2)/length(pred))
  return(list(int=m1$coef[1,1],t1=m1$coef[3,1],t2=m1$coef[2,1],t3=m1$coef[4,1],rms=rms))
}




