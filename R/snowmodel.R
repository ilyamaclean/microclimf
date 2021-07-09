#' Computes snow energy balance
.SnowEnergyBalance<-function(tc,u2,rh,pk,Rsw,skyem,snowalb,snowtemp,snowem=0.99,zm=0.002,umin=0.5) {
  # Variables set in function
  zpd=6.5*zm
  zh<-0.2*zm
  sb<-5.67*10^-8
  ph<-43.3
  cp<-29.3
  lambdaS<-51135.14
  lambdaW<-44526
  lambda<-ifelse(snowtemp>0,lambdaW,lambdaS)
  # Radiation emitted
  Rem<-snowem*sb*(snowtemp+273.15)^4
  # Radiation absorbed
  Rlwd<-skyem*sb*(tc+273.15)^4
  Rabs<-(1-snowalb)*Rsw+snowem*Rlwd
  # Sensible heat
  u2[u2<umin]<-umin
  gHa<-(0.4^2*ph*u2)/(log((2-zpd)/zm)*log((2-zpd)/zh))
  H<-cp*gHa*(snowtemp-tc)
  # Latent heat
  ests<-ifelse(snowtemp>0,0.61078*exp((17.27*snowtemp)/(snowtemp+237.3)),0.61078*exp((21.875*snowtemp)/(snowtemp+265.5)))
  esta<-ifelse(tc>0,0.61078*exp((17.27*tc)/(tc+237.3)),0.61078*exp((21.875*tc)/(tc+265.5)))
  ea<-(rh/100)*esta
  L<-((lambda*gHa)/pk)*(ests-ea)
  # Snow heat storage
  return(list(Rabs=Rabs,Rem=Rem,L=L,H=H))
}
#' Fits snow model (one step)
.FitSnowOne<-function(weather,snowdepth,precd,meltfact=6.54,snowem=0.99,zm=0.002,umin=0.5,astc=1.5) {
  predsnowd<-pSnow(weather,precd,meltfact,snowem,zm,umin,astc)$snowdepth
  RMS<-sqrt(sum((predsnowd-snowdepth)^2)/length(snowdepth))
  RMS
}
#' Estimates snow temperature
#'
#' The function `SnowTemp` estimates snow temperatures
#'
#' @param tc a vector of air temperature at reference height (deg C)
#' @param u2 a vector of wind at reference height (m/s)
#' @param rh a vector of relative humidities (percentage)
#' @param pk a vector of atmospheric pressure (kPa)
#' @param Rsw a vector of global shortwave radiation values (W/m^2)
#' @param snowalb a vector of snow albedos (as returned by [pSnow()])
#' @param snowem optionally, numeric value of snow emissivity
#' @param zm optionally, numeric value of roughness length for momentum transfer of snow surface (m)
#' @param umin optionally, numeric value indicating minimum wind speed used in conductivity calculations (m/s)
#' @return a vector of snow temperatures (deg C)
#' @export
SnowTemp<-function(tc,u2,rh,pk,Rsw,skyem,snowalb,snowem=0.99,zm=0.002,umin=0.5) {
  zpd=6.5*zm
  zh<-0.2*zm
  sb<-5.67*10^-8
  ph<-43.3
  cp<-29.3
  lambdaS<-51135.14
  lambda<-44526
  # Rabs
  Rabs<-(1-snowalb)*Rsw+skyem*sb*snowem*(tc+273.15)^4
  # gHr
  gHa<-(0.4^2*ph*u2)/(log((2-zpd)/zm)*log((2-zpd)/zh))
  gr<-4*snowem*sb*(tc+273.15)^3
  gHr<-gHa+gr
  # Deficit
  esta<-ifelse(tc>0,0.61078*exp((17.27*tc)/(tc+237.3)),0.61078*exp((21.875*tc)/(tc+265.5)))
  ea<-(rh/100)*esta
  DD<-esta-ea
  # s
  delta<-(4098*esta)/(tc-35.85)^2
  s<-delta/pk
  lam<-ifelse(tc<0,lambdaS,lambda)
  top<-Rabs-snowem*sb*(tc+273.15)^4-((lam*gHa*DD)/pk)
  btm<-gHr*cp+lam*s*gHa
  tsnow<-tc+top/btm
  tsnow
}
#' Estimates snow depth and other snow variables
#'
#' The function `pSnow` estimates snow depth, albedo, density and temperature
#'
#' @param weather a data.frame of weather variables (see details).
#' @param precd a vector of daily precipitation (mm).
#' @param meltfact a snow melt coefficient, as returned by [FitSnow()]
#' @param snowem optionally, numeric value of snow emissivity
#' @param zm optionally, numeric value of roughness length for momentum transfer of snow surface (m)
#' @param umin optionally, numeric value indicating minimum wind speed used in conductivity calculations (m/s)
#' @param astc optionally, numeric value indicating the temperature at which precipitation falls as snow (deg C)
#' @return a list of the following:
#' @return `snowdepth` snow depth (cm)
#' @return `snowdens` snow density (g/cm^3)
#' @return `snowtemp` snow temperature (deg C)
#' @return `snowalb` snow albedo (range 0 to 1)
#' @export
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`.
pSnow <- function(weather,precd,meltfact=6.54,snowem=0.99,zm=0.002,umin=0.5,astc=1.5) {
  # Calculate time since snowfall
  zpd=6.5*zm
  prech<-rep(0,length(precd)*24)
  sel<-(c(1:length(precd))-1)*24+1
  prech[sel]<-precd
  prech<-ifelse(weather$temp>astc,0,prech)
  snowyn<-ifelse(prech>0,1,0)
  age<-0
  for (i in 2:length(snowyn)) {
    age[i]<-ifelse(snowyn[i]==1,0,age[i-1]+1)
  }
  aged<-matrix(age,ncol=24,byrow=T)
  aged<-apply(aged,1,mean)
  aged<-aged/24
  # Calculate snow albedo
  snowalb<-(-9.8740*log(aged)+78.3434)/100
  snowalb[snowalb>0.95]<-0.95
  snowalb<-rep(snowalb,each=24)
  # Predict snow temperature
  psnowtemp<-SnowTemp(weather$temp,weather$windspeed,weather$relhum,weather$pres,
                      weather$swrad,weather$skyem,snowalb,snowem,zm,umin)
  # Energy balance
  SEB<-.SnowEnergyBalance(weather$temp,weather$windspeed,weather$relhum,weather$pres,
                          weather$swrad,weather$skyem,snowalb,psnowtemp,snowem,zm,umin)
  # Calculate net heat flux of snow
  G<-SEB$Rabs-SEB$Rem-SEB$H-SEB$L  # W / m2
  G[G<0]<-0
  # Calculate snow melt
  SHF<-334000
  meltkg<-(G*3600)/SHF
  meltcm<-(100*meltkg)/(0.2178*1000)
  # melt adjusted
  meltcm[psnowtemp<0]<-0
  meltcm<-meltcm*meltfact
  # Calculate rain melt in cm
  tcd<-matrix(weather$temp,ncol=24,byrow=T)
  tcd<-apply(tcd,1,mean)
  rmeltd<-ifelse(tcd>0,tcd*rainfall*0.0125,0)
  rmeltcm<-rep(rmeltd/24,each=24)/10
  # Calculate snow fall (cm)
  fallcm<-(prech*0.997)/(0.2178*10)
  # Calculate depth
  predsnowd<-0
  for (i in 2:(length(fallcm)+1)) {
    predsnowd[i]<-predsnowd[i-1]+fallcm[i-1]-meltcm[i-1]-rmeltcm[i-1]
    predsnowd[i]<-ifelse(predsnowd[i]<0,0,predsnowd[i])
  }
  predsnowd<-predsnowd[-1]
  # calculate snow pack age
  age<-0
  for (i in 1:length(predsnowd)) {
    age[i+1]<-ifelse(predsnowd[i]>0,age[i]+1/24,0)
  }
  age<-age[-1]
  # Calculate snow density
  pdensity<-(0.5979-0.2178)*(1-exp(-0.001*predsnowd-0.0038*age))+0.2178
  predsnowd<-predsnowd*(0.2178/pdensity)
  return(list(snowdepth=predsnowd,snowdens=pdensity,snowtemp=psnowtemp,
              snowalb=snowalb))

}
#' Derives snow melt coefficient from empirical data
#'
#' The function `fitsnow` estimates the snow melt coefficient from empirical data (or
#' the outputs of a more complex model) and returns this along with
#' estimates of snow depth and the RMS error
#'
#' @param weather a data.frame of weather variables (see details).
#' @param vector of empirically-derived hourly snow depths (cm)
#' @param precd a vector of daily precipitation (mm).
#' @param plotout optional logical indicating whether to plot observed and predicted data
#' @param snowem optionally, numeric value of snow emissivity
#' @param zm optionally, numeric value of roughness length for momentum transfer of snow surface (m)
#' @param umin optionally, numeric value indicating minimum wind speed used in conductivity calculations (m/s)
#' @param astc optionally, numeric value indicating the temperature at which precipitation falls as snow (deg C)
#' @return a list of the following:
#' @return `psnowdepth` predicted snow depth (cm)
#' @return `RMS` Root-mean-square error of predicted snow depths (cm)
#' @return `meltfact` a snow melt factor for use with [pSnow()]
#' @export
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`.
fitsnow<-function(weather,snowdepth,precd,plotout=T,snowem=0.99,zm=0.002,umin=0.5,astc=1.5) {
  rms<-0
  for (i in 1:20) {
    mf<-i
    rms[i]<-.FitSnowOne(weather,snowdepth,rainfall,mf,snowem,zm,umin,astc)
  }
  sel1<-which(rms==min(rms))[1]
  for (i in 1:20) {
    mf<-sel1-1+i/10
    rms[i]<-.FitSnowOne(weather,snowdepth,rainfall,mf,snowem,zm,umin,astc)
  }
  sel2<-which(rms==min(rms))[1]
  for (i in 1:20) {
    mf<-sel1-1+sel2/10-1/10+i/100
    rms[i]<-.FitSnowOne(weather,snowdepth,rainfall,mf,snowem,zm,umin,astc)
  }
  sel3<-which(rms==min(rms))[1]
  for (i in 1:20) {
    mf<-sel1-1+sel2/10-1/10+sel3/100-1/100+i/1000
    rms[i]<-.FitSnowOne(weather,snowdepth,rainfall,mf,snowem,zm,umin,astc)
  }
  sel4<-which(rms==min(rms))[1]
  mf<-sel1-1+sel2/10-1/10+sel3/100-1/100+sel4/1000
  # Generate snow predictions
  psnowd<-pSnow(weather,precd,mf,snowem,zm,umin,astc)$snowdepth
  if (plotout) {
    mx<-max(psnowd,snowdepth)
    plot(psnowd,type="l",col=rgb(0,0,0,0.5),ylim=c(0,mx),xlab="Hour of Year",ylab="Snow depth",main="")
    par(new=T)
    plot(snowdepth,col=rgb(1,0,0,0.5),type="l",ylim=c(0,mx),xlab="",ylab="",main=paste("RMS:",round(rms[sel4],3)))
  }
  return(list(psnowdepth=psnowd,RMS=rms[sel4],meltfact=mf))
}
