#' Spatially distribute soil moisture by topographic wetness index
#'
#' The function `soilmdistribute` spatially distrubutes soil moisture by the Bevan and Kirkby
#' topographic wetness index
#'
#' @param micro an object of class microin as returned by e.g. modelin
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from he dtm in micro of not supplied, but outer cells will be NA.
#' @return a 3D array of soil moistures with the same x and y dims as `dtm` and z
#' equivelent to length(soilm)
#' @import terra
#' @export
soilmdistribute <- function(micro, tfact = 1.5, twi = NA) {
  soilm<-micro$soilm
  dtm<-micro$dtm
  Smin<-micro$groundp_p$Smin
  Smax<-micro$groundp_p$Smax
  if (class(dtm)[1] == "PackedSpatRaster") dtm<-rast(dtm)
  if (class(twi)[1] == "logical") twi<-.topidx(dtm)
  rge<-Smax-Smin
  soilm[soilm<=Smin]<-Smin+0.001
  soilm[soilm>=Smax]<-Smax-0.001
  theta<-(soilm-Smin)/rge
  lt<-log(theta/(1-theta))
  ltwi<-log(.is(twi))/tfact
  me<-mean(ltwi,na.rm=T)
  add<-ltwi-me
  add<-.rta(rast(add),dim(soilm)[3])
  smout<-lt+add
  smout<-1/(1+exp(-smout))
  smout<-smout*rge+Smin
  smout<-.rast(smout, micro$dtm)
  smout<-mask(smout,micro$dtm)
  micro$soild<-as.array(smout)
  micro$progress<-1
  return(micro)
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
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' the dtm in micro if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' the dtm in micro if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from he dtm in micro if not supplied, but outer cells will be NA.
#' @return an object of class micro with additional terms added for subsequent modelling
#' @import terra
#' @export
twostream<-function(micro, reqhgt = 0.05, pai_a = NA, tfact=1.5, slr = NA, apr = NA, hor = NA, twi = NA) {
  if (micro$progress<1) micro<-soilmdistribute(micro,tfact,twi)
  dtm<-micro$dtm
  dtm[is.na(dtm)]<-0
  # === (1a) Calculate solar altitude
  tme<-micro$tme
  jd<-.jday(tme=tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  salt<-.solalt(lt,micro$lat,micro$long,jd)
  sazi<-.solazi(lt,micro$lat,micro$long,jd)
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
  svf<-0.5*cos(2*msl)+0.5
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
  # === (1g) Adjust parameters for gap fraction and inclined surface
  n<-(micro$vha-reqhgt)/micro$vha
  n[n<0]<-0
  pai_a<-micro$pai_a/(1-micro$clump*n)
  pai_t<-micro$pai/(1-micro$clump)
  gref2<-1-(((1-micro$gref)*sin(alt))/si)
  s<-which(is.infinite(gref2))
  gref2[s]<-micro$gref[s]
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
  albedo<-(micro$dirr*sin(alt)*albb+micro$difr*albd)/(micro$dirr*sin(alt)+micro$difr)
  albedo[is.na(albedo)]<-0.23
  albedo[albedo>0.99]<-0.99
  albedo[albedo<0.01]<-0.01
  # === (1k) Calculate ground absorbed radiation
  # Shortwave
  svfa<-.rta(rast(svf),length(jd))
  # ~~ Contribution of direct to downward diffuse
  Rdbm<-(1-micro$clump^Kc)*((p8/sig)*exp(-kkd$kd*pai_t)+p9*exp(-h*pai_t)+p10*exp(h*pai_t))+micro$clump^Kc
  s<-which(is.na(Rdbm) | Rdbm<0)
  Rdbm[s]<-0
  Rdbm[Rdbm>1]<-1
  # ~~ Downward diffuse
  Rddm<-(p3*exp(-h*pai_t)+p4*exp(h*pai_t))+micro$clump^2
  s<-which(is.na(Rddm))
  Rddm[s]<-1
  # ~~ Downward direct
  Rbgm<-exp(-kkd$kd*pai_t)+micro$clump^Kc
  # ~~ Ground absorbed
  micro$radGdir<-micro$dirr*si*Rbgm
  micro$radGdif<-micro$difr*svfa*Rddm+micro$dirr*sin(alt)*Rdbm
  micro$radGsw<-(1-micro$gref)*(micro$radGdir+micro$radGdif)
  # Longwave
  trd<-exp(-pai_t)+micro$clump^2
  micro$lwout<-0.97*5.67*10^-8*(micro$tc+273.15)^4 # Longwave emitted
  lwsky<-micro$skyem*micro$lwout # Longwave radiation down from sky
  micro$radGlw<-0.97*(trd*svfa*lwsky+(1-trd)*(1-svfa)*micro$lwout)
  # === (1l) Calculate canopy and ground combined absorbed radiation
  trb<-exp(-kkd$kd*pai_t)+micro$clump^Kc
  micro$radCsw<-(1-albedo)*(micro$dirr*sin(alt)*(1-trb)+svfa*micro$difr*(1-trd))+micro$radGsw
  micro$radClw<-svfa*lwsky
  if (reqhgt > 0) {
    # === (1m) Calculate up and downstream at reqhgt
    # Downward
    trcb<-micro$clump^(Kc*n) # transmission through canopy gaps (direct, reqhgt)
    trcd<-micro$clump^(2*n)  # transmission through canopy gaps (diffuse)
    # ~~ Contribution of direct to downward diffuse
    Rdbm<-(1-trcb)*((p8/sig)*exp(-kkd$kd*pai_a)+p9*exp(-h*pai_a)+p10*exp(h*pai_a))+trcb
    s<-which(is.na(Rdbm) | Rdbm<0)
    Rdbm[s]<-0
    s<-which(Rdbm>1)
    Rdbm[s]<-1
    # ~~ Contribution of direct to upward diffuse
    Rubm<-(1-trcb)*((p5/sig)*exp(-kkd$kd*pai_a)+p6*exp(-h*pai_a)+p7*exp(h*pai_a))+trcb
    s<-which(is.na(Rubm) | Rubm<0)
    Rubm[s]<-0
    s<-which(Rubm>1)
    Rubm[s]<-1
    # ~~ Downward diffuse
    Rddm<-(p3*exp(-h*pai_a)+p4*exp(h*pai_a))+trcd
    Rddm[is.na(Rddm)]<-1
    # ~~ Upward diffuse
    Rudm<-(p1*exp(-h*pai_a)+p2*exp(h*pai_a))+trcd
    s<-which(is.na(Rudm))
    Rudm[s]<-1-micro$gref[s]
    # Downward direct
    Rbgm<-exp(-kkd$kd*pai_a)+trcb
    Rbgm[is.na(Rbgm)]<-1
    # Calculate actual fluxes
    micro$Rbdown<-micro$dirr*sin(alt)*Rbgm # Direct down
    micro$Rddown<-micro$difr*svfa*Rddm+micro$dirr*sin(alt)*Rdbm
    micro$Rdup<-micro$difr*svfa*Rudm+micro$dirr*sin(alt)*Rubm
    # === (1n) Calculate leaf swabs
    a[is.na(a)]<-mean(a,na.rm=T)
    micro$radLsw<-0.5*a*(micro$Rddown+micro$Rdup+kkd$k*sin(alt)*micro$Rbdown)
    micro$si<-si
    micro$alt<-alt
  }
  micro$progress<-2
  # Clean micro
  micro$lat<-NULL
  micro$long<-NULL
  micro$vegx<-NULL
  micro$lref<-NULL
  micro$ltra<-NULL
  micro$clump<-NULL
  return(micro)
}
#' Downscales wind
#'
#' @description The function `wind` downscales wind speed accounting for
#' vegetation and terrain
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above ground for which wind speeds are wanted. If negative (below ground) wind friction velocity only is returned
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from the
#' dtm in micro if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' the dtm ion micro if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @param twi optional SpatRast object of topographic wetness index values.
#' Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @return an object of class microin with the following components added:
#' @return `uf` wind friction velocity (m/s)
#' @return `uz` if `reqhgt > 0` wind speed at height `reqhgt` (m/s)
#' @return `gHa` convective conductance (mol/m^2/s)
#' @import terra zoo
#' @export
#' @details In downscaling wind, two processes are accounted for. Firstly the drag effects
#' of vegetation on wind , which ultimately dictate the wind height profile. Secondly,
#' the effects of local vegetation and terrain on wind speed. Terrain effects are calculated
#' by applying the topographic shelter coefficient described in Maclean et al (2019) Methods
#' in Ecology and Evolution, 10:280-290.
wind <- function(micro, reqhgt = 0.05, pai_a = NA, tfact = 1.5, slr = NA, apr = NA,
                 hor = NA, twi = NA, wsa = NA) {
  # Run two-stream model if not run
  if (micro$progress < 2) micro<-twostream(micro,reqhgt,pai_a,tfact,slr,apr,hor,twi)
  # Calculate ratio of u* ref to u* pixel (no diabnatic correction)
  dp<-with(micro$vegp_p,.zeroplanedis(h,pai))
  zmp<-with(micro$vegp_p,.roughlength(h,pai,dp))
  ufp<-with(micro,(0.4*u2)/log((zref-dp)/zmp))
  ufa<-with(micro,(0.4*u2)/log((zref-d)/zm))
  umu<-ufp/ufa
  # Calculate wind terrain shelter coefficient
  dsm<-micro$dtm+micro$veghgt
  s<-1
  if (res(micro$dtm)[1]<=100) s<-10
  if (class(wsa) == "logical") wsa<-.windsheltera(micro$dtm,micro$zref,s)
  ws<-.windshelter(micro$winddir,dsm,2,s,wsa)
  # Calculate friction velocity
  micro$uf<-micro$ufp*umu*ws
  # Calculate wind speed at height z
  if (reqhgt > 0) {
    s1<-which(reqhgt >= micro$vha)
    uza<-suppressWarnings(with(micro,(uf[s1]/0.4)*log((reqhgt-d[s1])/zm[s1])))
    sel<-which(is.na(uza) | uza<micro$uf[s1])
    uza[sel]<-micro$uf[s1[sel]]
    # Below canopy
    s2<-which(reqhgt < micro$vha)
    uh<-suppressWarnings(with(micro,(uf[s2]/0.4)*log((vha[s2]-d[s2])/zm[s2])))
    sel<-which(is.na(uh) | uh<micro$uf[s2])
    uh[sel]<-micro$uf[s2[sel]]
    sel<-which(uh>micro$u2[s2])
    uh[sel]<-micro$uf[s2[sel]]
    Be<-(micro$uf[s2])/uh
    Be[Be>1]<-1
    Be[Be<0.001]<-0.001
    a<-micro$pai[s2]/micro$vha[s2]
    Lc<-(0.25*a)^-1
    Lm<-2*Be^3*Lc
    uzb<-uh*exp(Be*(reqhgt-micro$vha[s2])/Lm)
    # Combine below cand above canopy
    uz<-array(NA,dim=dim(micro$vha))
    uz[s1]<-uza
    uz[s2]<-uzb
    # Cap at max wind speed
    sel<-which(uz>micro$u2)
    uz[sel]<-micro$u2[sel]
    micro$uz<-uz
  }
  # Calculate convective conductance
  micro$gHa<-with(micro,.gturb(uf,d,zm,zref,psi_h=0,0.01))
  micro$progress<-3
  # Clean micro
  micro$d<-NULL
  micro$zm<-NULL
  micro$u2<-NULL
  return(micro)
}
#' Applies Penman-Monteith equation
#'
#' @description The function `PenMont` applies a varient of the Penman-Monteith
#' equation, to calculate surface temperature and the energy fluxes at the surface
#' @param tc air temperature (deg c)
#' @param pk atmospheric pressure (kPa)
#' @param ea vapour pressure of air (kPa)
#' @param radabs absorbed long and shortwave radiation (W/m2)
#' @param gHa boundary layer conductances for heat (mol/m^2/s)
#' @param gV conductance for vapour (mol/m^2/s)
#' @param es optional saturated vapour pressure (kPa) at temperature tc (calculated if NA)
#' @param tdew optional dewpoint temperature (deg C) (calculated if NA, see details)
#' @param surfwet proportion of surface acting as wet surface
#' @param G optional ground heat flux (W/m^2)
#' @param allout optional logical indicating whether to return Latent, Sensible and
#' Emitted radiation flux in addition to surface temperature
#' @return a list with the following components:
#' @return `T0` surface temperature (deg C)
#' @return `H` the sensible heat flux (w/m^2) (NA if allout set to FALSE)
#' @return `L` the latent heat flux (w/m^2) (NA if allout set to FALSE)
#' @return `Rem` the emitted radiation flux (w/m^2) (NA if allout set to FALSE)
#' @details The paramater `surfwet` determines how much of the canopy should be treated as
#' wet surface. However, except when extremely droughted, the matric potential of
#' of soil leaves is such that `surfwet` ~ 1.
#' @export
PenMont <- function(tc,pk,ea,radabs,gHa,gV,es=NA,tdew=NA,surfwet=1,G=0,allout=TRUE) {
  if (is.na(es[1])) es<-.satvap(tc)
  if (is.na(tdew[1])) tdew<-.dewpoint(ea,tc)
  delta<-.delta(tc)
  sb<-5.67*10^-8
  gHr<-gHa+(4*0.97*sb*(tc+273.15)^3)/29.3
  Rem<-0.97*sb*(tc+273.15)^4
  la<-45068.7-42.8428*tc
  sel<-which(tc<0)
  la[sel]<-51078.69-4.338*tc[sel]-0.06367*tc[sel]^2
  m<-la*(gV/pk)
  L<-m*(es-ea)*surfwet
  dT<-(radabs-Rem-L-G)/(29.3*gHr+m*delta)
  dTmx<- -0.6273*max(tc,na.rm=T)+49.79
  dT[dT>dTmx]<-dTmx
  T0<-dT+tc
  T0<-.lim(T0,tdew)
  T0[T0>80]<-80
  out<-list(T0=T0,H=NA,L=NA,Rem=NA)
  if (allout) {
    out$H<-29.3*gHa*(T0-tc)
    out$L<-m*(.satvap(T0)-ea)*surfwet
    out$Rem<-0.97*sb*(T0+273.15)^4
  }
  return(out)
}
#' Calculates ground surface temperature (hourly)
#'
#' @description The function `soiltemp_hr` estimates ground surface temperature
#'
#' @param micro an object of class `microin` as returned by [modelin()] (see details)
#' @param reqhgt height above ground at which model outputs are needed (m).
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param slr an optional SpatRast object of slope values (Radians). Calculated from
#' the dtm in micro if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRast object of aspect values (Radians). Calculated from
#' the dtm in micro if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @param twi optional SpatRast object of topographic wetness index values.
#' Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @return an object of class micro with the following components added:
#' @return `T0` ground surface temperature (deg C)
#' @return `G` ground heat flux (W/m^2)
#' @return additional terms needed for subsequent modelling.
#' @import terra
#' @export
#' @seealso [wind()], [soiltemp_dy()]
#' @details if [wind()] has not been run to add additional elements to `micro`
#' it is automatically called.
soiltemp_hr  <- function(micro, reqhgt = 0.05, pai_a = NA, tfact = 1.5, slr = NA,
                         apr = NA, hor = NA, twi = NA, wsa = NA) {
  # run two-stream and wind functions if not run
  if (micro$progress<3) {
    micro<-wind(micro,reqhgt,pai_a,tfact,slr,apr,hor,twi,wsa)
  }
  # Estimate diurnal range in grid ground surface temperature with G set to zero
  radabs<-micro$radGsw+micro$radGlw  # Absorbed radiation
  # * Calculate soil relative humidity
  n<-dim(radabs)[3]
  surfwet<-with(micro,.soilrh(soild,.rta(soilb,n),.rta(psi_e,n),.rta(Smax,n),tc))
  # * Calculate soil surface temperature
  T0<-with(micro,PenMont(tc,pk,ea,radabs,gHa,gHa,estl,tdew,surfwet,G=0))
  # * Calculate diurnal range
  dtr<-.hourtoday(T0$T0,max)-.hourtoday(T0$T0,min)
  # * Estimate point surface temperature with G set to zero
  dtp<-.hourtoday(micro$T0p,max)-.hourtoday(micro$T0p,min)
  # * Calculate ratio of diurnal temperature fluxes
  dtR<-.ehr(dtr/dtp)
  # Thermal conductance and damping depth (grid)
  kDDg<-with(micro,.soilcond(.rta(rho,n),.rta(Vm,n),.rta(Vq,n),.rta(Mc,n),soild))
  kDDp<-with(micro$groundp_p,.soilcond(rho,Vm,Vq,Mc,micro$soilm))
  # Calculate Gmu
  Gmu<-dtR*(kDDg$k*kDDp$DD)/(kDDp$k*kDDg$DD)
  # Compute ground heat flux
  G<-micro$Gp*Gmu
  # Set limits to G:
  # Calculate daily Rnet cycle
  Rnet<-radabs-T0$Rem
  Rd<-.ehr(pmax(.hourtoday(Rnet,max),-.hourtoday(Rnet,min)))
  G<-pmin(G,0.6*Rd)
  micro$G<-pmax(G,-0.6*Rd)
  # Calculate soil surface temperature
  micro$T0<-with(micro,PenMont(tc,pk,ea,radabs,gHa,gHa,estl,tdew,surfwet,G,allout=FALSE))$T0
  micro$DDg<-kDDg$DD
  micro$DDp<-kDDp$DD
  micro$progress<-4
  # Clean micro
  micro$rho<-NULL
  micro$Vm<-NULL
  micro$Vq<-NULL
  micro$Mc<-NULL
  micro$soilb<-NULL
  micro$psi_e<-NULL
  return(micro)
}
#' Estimate temperature and humidity at specified height above ground
#'
#' @description The function `temphumE` runs the above ground component of
#' the microclimate model.
#'
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above ground for which wind speeds are wanted. If negative (below ground) wind friction velocity only is returned
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param surfwet an optional coefficient describing the fraction of the
#' vegetation surface that acts like a free water surface. Calculated from rainfall if not supplied.
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' the dtm in micro if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' the dtm in micro if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @param twi optional SpatRast object of topographic wetness index values.
#' Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from the dtm in micro if not supplied, but outer cells will be NA.
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
temphumE<-function(micro, reqhgt = 0.05, pai_a = NA, tfact = 1.5, surfwet = NA,
                   slr = NA, apr = NA, hor = NA, twi = NA, wsa = NA) {
  # run soiltemp function if not run
  if (micro$progress<4) {
    micro<-soiltemp_hr(micro,reqhgt,pai_a,tfact,slr,apr,hor,twi,wsa)
  }
  # Calculate vapour conductivity
  Rsw<-with(micro,difr+dirr*si)
  n<-dim(Rsw)[3]
  # Calculate ground wetness
  eT<-.satvap(micro$T0)-micro$ea
  gwet<-.predwet(eT)
  if (class(surfwet) == "logical") {
    surfwet<-with(micro,(soild-.rta(Smin,n))/.rta(Smax-Smin,n))
  }
  gwet<-pmax(gwet,surfwet)
  gv<-.layercond(Rsw, .rta(micro$gsmax,n))*3
  gV<-1/(1/micro$gHa+1/gv)
  # Calculate canopy temperature
  canabs<-micro$radCsw+micro$radClw
  TH<-with(micro,PenMont(tc,pk,ea,canabs,gHa,gV,estl,tdew,surfwet,G))
  # Temperature and vapour above canopy
  # Set z to canopy height if reqhgt < canopy height
  z<-micro$vha*0+reqhgt
  sel<-which(micro$vha>z)
  z[sel]<-micro$vha[sel]
  micro<-.TVabove(TH,micro,z,surfwet)
  # Leaf temperature
  Thl<-.leafbalance(micro,TH,surfwet)
  micro$tleaf<-Thl$T0
  # Temperature below canopy
  Flux<-TH$H*(1-exp(-micro$pai))
  Fluxz<-Thl$H
  mu<-29.3*43
  SH<-micro$Tz*mu
  SG<-micro$T0*mu
  Tzb<-.Langapprox(micro,reqhgt,Flux,Fluxz,SH,SG)/mu
  selb<-which(reqhgt<micro$vha)
  micro$Tz[selb]<-Tzb[selb]
  # Set limits
  tmx<-with(micro,pmax(tleaf,tc,T0,TH$T0))
  tmn<-with(micro,pmin(tleaf,tc,T0,TH$T0))
  Tz<-pmin(micro$Tz,tmx)
  Tz<-pmax(Tz,tmn)
  # Humidity below canopy
  Flux<-TH$L*(1-exp(-micro$pai))
  Fluxz<-Thl$L
  mu<-44525.68*43/micro$pk
  SH<-micro$ez*mu
  SG<-.satvap(micro$T0)*gwet*mu
  estb<-.Langapprox(micro,reqhgt,Flux,Fluxz,SH,SG)/mu
  micro$ez[selb]<-estb[selb]
  rh<-(micro$ez/.satvap(micro$Tz))*100
  rh[rh>100]<-100
  # Return values needed
  micro<-list(Tz=Tz,tleaf=micro$tleaf,T0=micro$T0,soilm=micro$soild,relhum=rh,windspeed=micro$uz,
              Rdirdown=micro$Rbdown,Rdifdown=micro$Rddown,Rlwdown=Thl$lwdn,
              Rswup=micro$Rdup,Rlwup=Thl$lwup)
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
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' the dtm in micro if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' the dtm in micro if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from the dtm in micro if not supplied, but outer cells will be NA.
#' @return `Tz` Array of air temperatures at height `reqhgt` below ground (deg C)
#' @seealso [temphumE()] for running the above ground component of the microclimate
#' model and [below_dy()] for running the model in daily time increments.
#'
#' @import terra zoo abind
#' @importFrom stats filter
#' @export
below_hr<-function(micro, reqhgt, pai_a = NA, tfact = 1.5,
                   slr = NA, apr = NA, hor = NA, twi = NA, wsa = NA) {
  if (is.null(micro$T0[1])) {
    micro<-soiltemp_hr(micro,reqhgt,pai_a,tfact,slr,apr,hor,twi,wsa)
  }
  # Calculate n
  n<- -118.35*reqhgt/mean(micro$DDp)
  hiy<-ifelse((micro$tme$year[1]+1900)%%4==0,366*24,365*24)
  # Test whether comple time sequence
  tst<-as.numeric(micro$tme[25])-as.numeric(micro$tme[24])
  if (tst > 3600) { # Not complete time sequence
    if (n >= 24) {
      micro$Tz<-.rta(rast(apply(micro$T0,c(1,2),mean)),dim(micro$T0)[3])
      dTa<-micro$Tz-mean(micro$Tg,na.rm=T)
    }
    if (n < hiy) {
      # Calculate daily annual wgts
      w1<-24/n
      w2<-n/hiy
      wgt<-w1/(w1+w2)
      # Temperature difference daily
      Tdp<-.ehr(.hourtoday(micro$Tg,mean))
      Tdg<-.ehr(.hourtoday(micro$T0,mean))
      dTd<-Tdg-Tdp
      micro$Tz<-wgt*dTd+(1-wgt)*dTa+micro$Tbp
    }
    if (n < 24) {
      # Calculate hourly & daily wgts
      w1<-1/n
      w2<-n/24
      wgt<-w1/(w1+w2)
      # Temperature difference hourly
      dTh<-micro$T0-micro$Tg
      micro$Tz<-wgt*dTh+(1-wgt)*dTd+micro$Tbp
    }
  } else { # Complete time sequence
    if (n < length(micro$tme)) {
      micro$Tz<-aperm(apply(micro$T0,c(1,2),.ma,n),c(2,3,1))
    } else {
      micro$Tz<-.rta(rast(apply(micro$T0,c(1,2),mean)),dim(micro$T0)[3])
    }
  }
  class(micro)<-"microout"
  return(micro)
}
#' Calculates ground surface temperature (daily)
#'
#' @description The function `soiltemp_dy` estimates ground surface temperature
#'
#' @param micro_dy an object of class `microindaily` as returned by [modelin_dy()] (see details)
#' @param reqhgt height above ground at which model outputs are needed (m).
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param slr an optional SpatRast object of slope values (Radians). Calculated from
#' the dtm in micro_dy if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRast object of aspect values (Radians). Calculated from
#' the dtm in micro_dy if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from the dtm in micro_dy if not supplied, but outer cells will be NA.
#' @param twi optional SpatRast object of topographic wetness index values.
#' Calculated from the dtm in micro_dy if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from the dtm in micro_dy if not supplied, but outer cells will be NA.
#' @return an object of class microdaily with the following components added:
#' @return `T0` ground surface temperature (deg C)
#' @return `G` ground heat flux (W/m^2)
#' @return additional terms needed for subsequent modelling.
#' @import terra
#' @export
#' @seealso [wind()], [soiltemp_hr()]
#' @details if [wind()] has not been run to add additional elements to `micro`
#' it is automatically called.
soiltemp_dy  <- function(micro_dy, reqhgt = 0.05, pai_a = NA, tfact = 1.5,
                         slr = NA, apr = NA, hor = NA, twi = NA, wsa = NA) {
  # run two-stream and wind functions if not run
  micro_mn<-micro_dy$micro_mn
  micro_mx<-micro_dy$micro_mx
  if (micro_mn$progress<3) {
    micro_mn<-wind(micro_mn,reqhgt,pai_a,tfact,slr,apr,hor,twi,wsa)
    micro_mx<-wind(micro_mx,reqhgt,pai_a,tfact,slr,apr,hor,twi,wsa)
  }
  # Estimate diurnal range in grid ground surface temperature with G set to zero
  radabs_mn<-micro_mn$radGsw+micro_mn$radGlw  # Absorbed radiation (min)
  radabs_mx<-micro_mx$radGsw+micro_mx$radGlw  # Absorbed radiation (max)
  # * Calculate soil relative humidity
  n<-dim(radabs_mn)[3]
  surfwet_mn<-with(micro_mn,.soilrh(soild,.rta(soilb,n),.rta(psi_e,n),.rta(Smax,n),tc))
  surfwet_mx<-with(micro_mx,.soilrh(soild,.rta(soilb,n),.rta(psi_e,n),.rta(Smax,n),tc))
  # * Calculate soil surface temperature
  T0_mn<-with(micro_mn,PenMont(tc,pk,ea,radabs_mn,gHa,gHa,estl,tdew,surfwet_mn,G=0))
  T0_mx<-with(micro_mx,PenMont(tc,pk,ea,radabs_mx,gHa,gHa,estl,tdew,surfwet_mx,G=0))
  # * Calculate diurnal range
  dtr<-T0_mx$T0-T0_mn$T0
  # * Estimate point surface temperature with G set to zero
  dtp<-micro_mx$T0p-micro_mn$T0p
  # * Calculate ratio of diurnal temperature fluxes
  dtR<-dtr/dtp
  # Thermal conductance and damping depth (grid)
  kDDg_mn<-with(micro_mn,.soilcond(.rta(rho,n),.rta(Vm,n),.rta(Vq,n),.rta(Mc,n),soild))
  kDDg_mx<-with(micro_mx,.soilcond(.rta(rho,n),.rta(Vm,n),.rta(Vq,n),.rta(Mc,n),soild))
  kDDp_mn<-with(micro_mn$groundp_p,.soilcond(rho,Vm,Vq,Mc,micro_mn$soilm))
  kDDp_mx<-with(micro_mx$groundp_p,.soilcond(rho,Vm,Vq,Mc,micro_mx$soilm))
  # Calculate Gmu
  Gmu_mn<-dtR*(kDDg_mn$k*kDDp_mn$DD)/(kDDp_mn$k*kDDg_mn$DD)
  Gmu_mx<-dtR*(kDDg_mx$k*kDDp_mx$DD)/(kDDp_mx$k*kDDg_mx$DD)
  # Compute ground heat flux
  G_mn<-micro_mn$Gp*Gmu_mn
  G_mx<-micro_mx$Gp*Gmu_mx
  # Set limits to G:
  # Calculate daily Rnet cycle
  Rnet_mn<-radabs_mn-T0_mn$Rem
  Rnet_mx<-radabs_mx-T0_mx$Rem
  Rd<-pmax(Rnet_mx,-Rnet_mn)
  G_mn<-pmin(G_mn,0.6*Rd)
  G_mx<-pmin(G_mx,0.6*Rd)
  micro_mn$G<-pmax(G_mn,-0.6*Rd)
  micro_mx$G<-pmax(G_mx,-0.6*Rd)
  # Calculate soil surface temperature
  micro_mn$T0<-with(micro_mn,PenMont(tc,pk,ea,radabs_mn,gHa,gHa,estl,tdew,surfwet_mn,G,allout=FALSE))$T0
  micro_mx$T0<-with(micro_mx,PenMont(tc,pk,ea,radabs_mx,gHa,gHa,estl,tdew,surfwet_mx,G,allout=FALSE))$T0
  micro_mn$DDg<-kDDg_mn$DD
  micro_mx$DDg<-kDDg_mx$DD
  micro_mn$DDp<-kDDp_mn$DD
  micro_mx$DDp<-kDDp_mx$DD
  micro_mn$progress<-4
  micro_mx$progress<-4
  # Clean micromin
  micro_mn$rho<-NULL
  micro_mn$Vm<-NULL
  micro_mn$Vq<-NULL
  micro_mn$Mc<-NULL
  micro_mn$soilb<-NULL
  micro_mn$psi_e<-NULL
  # Clean micromax
  micro_mx$rho<-NULL
  micro_mx$Vm<-NULL
  micro_mx$Vq<-NULL
  micro_mx$Mc<-NULL
  micro_mx$soilb<-NULL
  micro_mx$psi_e<-NULL
  micro_dy$micro_mn<-micro_mn
  micro_dy$micro_mx<-micro_mx
  return(micro_dy)
}
#' Estimate temperature below ground (daily)
#'
#' @description The function `below_dy` runs the below ground component of
#' the microclimate model at hourly time increments.
#'
#' @param micro object of class microin as returned by [modelin()]
#' @param reqhgt height above ground at which model outputs are needed (m). Must be negative.
#' @param pai_a plant area index above `reqhgt`. Determined from total `pai` if not supplied.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#' in topographic wetness (see [soilmdistribute()])
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' from the dtm in micro_dy if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' the dtm in micro_dy if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from the dtm micro_dy if not supplied, but outer cells will be NA.
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from the dtm in micro_dy if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from the dtm in micro_dy if not supplied, but outer cells will be NA.
#' @return `Tz` Array of air temperatures at height `reqhgt` below ground (deg C)
#' @seealso [temphumE()] for running the above ground component of the microclimate
#' model and [below_hr()] for running the model in hourly time increments.
#'
#' @import terra zoo abind
#' @importFrom stats filter
#' @export
below_dy<-function(micro_dy, reqhgt, pai_a = NA, tfact = 1.5,
                   slr = NA, apr = NA, hor = NA, twi = NA, wsa = NA) {

  if (is.null(micro_dy$micro_mn$T0[1])) {
    micro_dy<-soiltemp_dy(micro_dy,reqhgt,pai_a,tfact,slr,apr,hor,twi,wsa)
  }
  micro_mn<-micro_dy$micro_mn
  micro_mx<-micro_dy$micro_mx
  # Calculate n
  DDp<-(micro_mn$DDp+micro_mx$DDp)/2
  n<- -118.35*reqhgt/mean(DDp)
  hiy<-ifelse((micro_mn$tme$year[1]+1900)%%4==0,366*24,365*24)
  if (n >= 24) {
    T0<-(micro_mn$T0+micro_mx$T0)/2
    T0<-rast(apply((micro_mn$T0+micro_mx$T0)/2,c(1,2),mean))
    micro_mn$Tz<-.rta(T0,dim(micro_mn$T0)[3])
    micro_mx$Tz<-micro_mn$Tz
    Tg<-(micro_mn$Tg+micro_mx$Tg)/2
    dTa<-micro_mn$Tz-mean(Tg,na.rm=T)
  }
  if (n < hiy) {
    # Calculate daily annual wgts
    w1<-24/n
    w2<-n/hiy
    wgt<-w1/(w1+w2)
    # Temperature difference daily
    Tdp<-(micro_mn$Tg+micro_mx$Tg)/2
    Tdg<-(micro_mn$T0+micro_mx$T0)/2
    dTd<-Tdg-Tdp
    micro_mn$Tz<-wgt*dTd+(1-wgt)*dTa+micro_mn$Tbp
    micro_mx$Tz<-wgt*dTd+(1-wgt)*dTa+micro_mx$Tbp
  }
  if (n < 24) {
    # Calculate hourly & daily wgts
    w1<-1/n
    w2<-n/24
    wgt<-w1/(w1+w2)
    # Temperature difference hourly
    dTh_mn<-micro_mn$T0-micro_mn$Tg
    dTh_mx<-micro_mx$T0-micro_mx$Tg
    micro_mn$Tz<-wgt*dTh_mn+(1-wgt)*dTd+micro_mn$Tbp
    micro_mx$Tz<-wgt*dTh_mx+(1-wgt)*dTd+micro_mx$Tbp
  }
  micro_dy$micro_mn<-micro_mn
  micro_dy$micro_mx<-micro_mx
  class(micro_dy)<-"microoutdaily"
  return(micro_dy)
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
    windspeed<-ncvar_def(name="uz",longname=wname,units="m/s x 100",dim=list(east,north,times),
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
    soilmoist<-ncvar_def(name="soilm",longname="Soil surface moisture",units="Percentage volume",dim=list(east,north,times),
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
    ncnew<-nc_create(filename=nc.name,list(soiltemp,soilmoist))
    # Put variables in
    ncvar_put(ncnew,soiltemp,vals=atonc(mout$T0,100))
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
    ncnew<-nc_create(filename=nc.name,list(soiltemp))
    # Put variables in
    ncvar_put(ncnew,soiltemp,vals=atonc(mout$Tz,100))
    ncvar_put(ncnew,soilmoist,vals=atonc(mout$soilm,100))
    ncatt_put(ncnew,0,"Coordinate reference system",as.character(crs(dtm)))
    nc_close(ncnew)
  }
}


#' expand runmicro_big daily output to hourly and write as ncdf4 file
#'
#' @description the function `expandtonc` takes the daily output written by [runmicro_big()]
#' when hourly = F, expands this to hourly, and writes data out as an ncdf4 file.
#'
#' @param filein file to be loaded as written out by [runmicro_big()]
#' @param file to be written out. Must have .nc extension.
#' @import ncdf4 terra
#' @export
expandtonc <- function(filein, fileout) {
  mout<-readRDS(filein)
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
