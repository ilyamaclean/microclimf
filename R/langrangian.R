.SourceD<-function(reqhgt,swabs,tcan,tground,tref,skyem,pai,pai_a,clump,Hratio,folden) {
  paie<-pai/(1-clump)
  pai_ae<-pai_a/(1-clump)
  # Lw down
  sb<-5.67*10^-8
  trd<-exp(-pai_ae)+clump
  lwdown<-trd*skyem*sb*(tref+273.15)^4+
    (1-trd)*0.97*sb*(tcan+273.15)^4
  # Lw up
  tru<-exp(-(paie-pai_ae))+clump
  lwup<-tru*sb*0.97*(tground+273.15)^4+
    (1-tru)*0.97*sb*(tcan+273.15)^4
  lwabs<-0.97*0.5*(lwdown+lwup)
  # H
  Rem<-0.97*sb*(tcan+273.15)^4
  H<-Hratio*(swabs+lwabs-Rem)
  S<-folden*H
  return(list(S=S,radabs=swabs+lwabs,Rldown=lwdown,Rlwup=lwup))
}
#' Computes near field concentration
.Cnear<-function(S,uf,hgt) {
  p1<-sqrt(uf)*(hgt^(2/3))
  slp<-0.42458*p1^0.85443
  Tn<-S/(slp*1000)
  Tn
}
#' Performs Langrangian simulation (temperature)
.LangrangianSimT<-function(reqhgt,micro,TH) {
  # Calculate d and zh of ground-layer
  d<-0.075*micro$vha
  zh<-0.001506*micro$vha
  lnr<-suppressWarnings(log((reqhgt-d)/zh)/log((micro$vha-d)/zh))
  sel<-which(reqhgt<(d+zh))
  # Calculate log-linear gradient
  Tzfo<-micro$T0-lnr*(micro$T0-micro$Tz)
  # Calculate linear gradient
  Tzfl<-micro$T0-(reqhgt/micro$vha)*(micro$T0-micro$Tz)
  Tzfo[sel]<-Tzfl[sel]
  # Partition according to canopy cover
  tr<-exp(-micro$pai)
  Tf<-tr*Tzfo+(1-tr)*Tzfl
  Tmx<-pmax(micro$T0,micro$Tz)
  Tmn<-pmin(micro$T0,micro$Tz)
  s1<-which(Tf>Tmx)
  s2<-which(Tf<Tmn)
  Tf[s1]<-Tmx[s1]
  Tf[s2]<-Tmn[s2]
  # Compute Source Density
  tref<-.vta(micro$climdata$temp,micro$dtm)
  skyem<-.vta(micro$climdata$skyem,micro$dtm)
  SR<-.SourceD(reqhgt,micro$radLs,TH$tcan,micro$T0,tref,skyem,micro$pai,
               micro$pai_a,micro$clump,TH$HR,micro$leafden)
  S<-SR$S
  leafabs<-SR$radabs
  # Compute near field temperature
  Tn<-.Cnear(S,micro$uf,micro$vha)
  To<-Tf+Tn
  # Replace with above canopy temperature if reqhgt > canopy height
  sel<-which(micro$vha<=reqhgt)
  To[sel]<-micro$Tz[sel]
  # Set limits
  dT<-To-micro$tc
  dTmx<- -0.6273*max(micro$tc,na.rm=T)+49.79
  dT[dT>dTmx]<-dTmx
  To<-micro$tc+dT
  To[To>72]<-72
  To<-.lim(To,micro$tdew)
  return(list(To=To,leafabs=leafabs,Rldown=SR$Rldown,Rlwup=SR$Rlwup))
}
# Calculates leaf temperature
.leaftemp<-function(micro,gs,reqhgt,tcan,leafabs) {
  Rem<-0.97*5.67*10^-8*(micro$tc+273.15^4)
  # Calculate conductivity
  n<-dim(Rem)[3]
  ld<-0.71*.rta(micro$leafd,n)
  gh<-0.023256*sqrt(micro$uz/ld)
  He<-0.7*(leafabs-Rem)
  gmn<-.gfree(ld/0.71,abs(He))
  sel<-which(gh<gmn)
  gh[sel]<-gmn[sel]
  TH<-PenMont(micro$Tz,micro$pk,micro$ea,leafabs,gh,
              gs,0.01,NA,NA,micro$tdew,1,T_est=tcan)
  micro$tleaf<-TH$tcan
  sel<-which(reqhgt>micro$vha)
  micro$tleaf[sel]<-tcan[sel]
  micro$gh<-gh
  return(micro)
}
.LangrangianSimV<-function(reqhgt,micro,ez,surfwet) {
  # Calculate soil surface effective vapour pressure
  n<-dim(ez)[3]
  Smin<-.rta(rast(micro$Smin),n)
  Smax<-.rta(rast(micro$Smax),n)
  rhs<-(micro$theta-Smin)/(Smax-Smin)
  e0<-.satvap(micro$T0)*rhs
  # Calculate d and zh of ground-layer
  d<-0.075*micro$vha
  zh<-0.001506*micro$vha
  # Calculate far field vapour pressure
  lnr<-suppressWarnings(log((reqhgt-d)/zh)/log((micro$vha-d)/zh))
  sel<-which(reqhgt<(d+zh))
  # Calculate log-linear gradient
  efo<-e0-lnr*(e0-ez)
  # Calculate linear gradient
  efl<-e0-(reqhgt/micro$vha)*(e0-ez)
  efo[sel]<-efl[sel]
  tr<-exp(-micro$pai)
  ef<-tr*efo+(1-tr)*efl
  # Compute near field vapour pressure
  gs<-.layercond(micro$Rbdown+micro$Rddown,micro$gsmax)
  gv<-1/(1/gs+1/micro$gh)
  es<-.satvap(micro$tleaf)*surfwet
  L<-44526*gv/micro$pk*(es-micro$ea)*surfwet
  S<-L*micro$leafden
  Cn<-.Cnear(S,micro$uf,micro$vha)*29.3*43
  Cn[Cn>2500]<-2500
  # Vapour pressure
  eo<-ef+micro$pk/(44526*43)*Cn
  sel<-which(micro$vha<=reqhgt)
  eo[sel]<-ez[sel]
  eas<-.satvap(micro$Tz)
  # RElative humidity
  rh<-(eo/eas)*100
  rh[rh>100]<-100
  rh[rh<30]<-30
  return(rh)
}
