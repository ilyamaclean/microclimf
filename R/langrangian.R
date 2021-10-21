#' Computes source concentration
.SourceD<-function(z,R0,tref,skyem,HRatio,pai_a,folden,hgt,lref,difprop,k=1) {
  ks<-k*sqrt(1-lref)
  # Direct radiation:
  Rb<-R0*(1-difprop)*(1-lref)*exp(-ks*pai_a)
  # Diffuse radiation:
  ks<-sqrt(1-lref)
  Rd<-R0*difprop*(1-lref)*exp(-ks*pai_a)
  # Longwave radiation
  Rem<-0.97*5.67*10^-8*(tref+273.15)^4
  tr<-exp(-pai_a)
  Rlwd<-Rem*((1-tr)+skyem*tr)
  H<-HRatio*(Rb+Rd+Rlwd-Rem)
  S<-folden*H
  S
}
#' Computes near field concentration
.Cnear<-function(S,uf,hgt) {
  p1<-sqrt(uf)*(hgt^(2/3))
  slp<-0.42458*p1^0.85443
  Tn<-S/(slp*1000)
  Tn
}
#' Performs Langrangian simulation (temperature)
.LangrangianSimT<-function(reqhgt,micro,TH,pai_a=NA,folden=NA) {
  # Calculate d and zh of ground-layer
  n<-dim(micro$tc)[3]
  hgt<-.rta(micro$hgt,n)
  d<-0.075*hgt
  zh<-0.001506*hgt
  lnr<-suppressWarnings(log((reqhgt-d)/zh)/log((hgt-d)/zh))
  sel<-which(reqhgt<(d+zh))
  # Calculate log-linear gradient
  Tzfo<-micro$T0-lnr*(micro$T0-micro$Tz)
  # Calculate linear gradient
  Tzfl<-micro$T0-(reqhgt/hgt)*(micro$T0-micro$Tz)
  Tzfo[sel]<-Tzfl[sel]
  # Partition according to canopy cover
  Tf<-micro$trlw*Tzfo+(1-micro$trdi)*Tzfl
  # Compute Source Density
  ldp<-foliageden(reqhgt,micro$vha,micro$pai)
  if (class(pai_a)[1]=="logical") pai_a<-ldp$pai_a
  if (class(folden)[1]=="logical") folden<-ldp$leafden
  cd<-micro$climdata
  R0<-.vta(cd$swrad,micro$gsmax)
  S<-.SourceD(reqhgt,R0,micro$tc,micro$skyem,TH$HR,pai_a,folden,
              micro$vha,micro$lref,micro$dp,micro$k)
  # Compute near field temperature
  Tn<-.Cnear(S,micro$uf,micro$vha)
  To<-Tf+Tn
  # Set limits
  sel<-which(To<micro$tdew)
  To[sel]<-micro$tdew[sel]
  tmx<-pmax(micro$T0,TH$tcan,micro$tc)
  sel<-which(To>tmx)
  To[sel]<-tmx[sel]
  sel<-which(micro$vha<=reqhgt)
  To[sel]<-micro$Tz[sel]
  return(list(To=To,pai_a=pai_a,folden=folden,R0=R0))
}
# Calculates leaf temperature
.leaftemp<-function(micro,pai_a,gs) {
  # Calculate radiation absorbption
  rad<-.radabs(micro,pai_a)
  radabs<-rad$radzsw+rad$radzlw
  Rem<-0.97*5.67*10^-8*(micro$tc+273.15^4)
  # Calculate conductivity
  n<-dim(Rem)[3]
  ld<-0.71*.rta(micro$leafd,n)
  gh<-0.023256*sqrt(micro$uz/ld)
  He<-0.7*(radabs-Rem)
  gmn<-.gfree(ld/0.71,abs(He))
  sel<-which(gh<gmn)
  gh[sel]<-gmn[sel]
  TH<-PenMont(micro$Tz,micro$pk,micro$ea,radabs,gh,
              gs,0,NA,micro$tdew,1)
  tleaf<-TH$tcan
  return(list(tleaf=tleaf,gh=gh,rad=rad))
}
# Calculates relative humidity
.LangrangianSimV<-function(reqhgt,micro,folden,ez,tlg) {
  # Calculate soil surface effective vapour pressure
  n<-dim(ez)[3]
  rhs<-.soilrh(micro$theta,
               .rta(raster(micro$soilb),n),
               .rta(raster(micro$psi_e),n),
               .rta(raster(micro$Smax),n),
               micro$T0)
  e0<-.satvap(micro$T0)*rhs
  # Calculate d and zh of ground-layer
  hgt<-.rta(micro$hgt,n)
  d<-0.075*hgt
  zh<-0.001506*hgt
  # Calculate far field vapour pressure
  lnr<-suppressWarnings(log((reqhgt-d)/zh)/log((hgt-d)/zh))
  sel<-which(reqhgt<(d+zh))
  # Calculate log-linear gradient
  efo<-e0-lnr*(e0-ez)
  # Calculate linear gradient
  efl<-e0-(reqhgt/hgt)*(e0-ez)
  efo[sel]<-efl[sel]
  # Partition according to canopy cover
  ef<-micro$trlw*efo+(1-micro$trlw)*efl
  # Compute near field vapour pressure
  gh<-tlg$gh
  radp<-tlg$rad
  gs<-.layercond(radp$rad_dir+radp$rad_dif,micro$gsmax)
  gv<-1/(1/gs+1/gh)
  es<-.satvap(tlg$tleaf)
  L<-44526*gv/micro$pk*(es-micro$ea)
  S<-L*folden
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
