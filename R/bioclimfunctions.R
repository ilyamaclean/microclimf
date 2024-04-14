#' Runs microclimate model to produce bioclim variables
#'
#' The function `runbioclim` runs the microclimate model and produces microclimate
#' equivalents of the 19 Worldclim bioclimate variables
#'
#' @param weather a data.frame or list of arrays of weather variables (as for [runpointmodel() or [runpointmodela()])
#' @param precip a vector or array of daily precipitation variables ((as for [runpointmodel() or [runpointmodela()]))
#' @param tme If weather is not a data.frame, an object of class POSIXlt giving the dates and times for each weather variable stored in the array. Set to NA if weather is a data,frame
#' @param reqhgt height for which temperatures are needed
#' @param micropoint optional subset output of point microclimate model for days with monthly monthly median, maximum and minimum temperatures. Calculated if not supplied.
#' @param vegp an object of class vegparams as returned by [vegpfromhab()]
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a SpatRaster object of elevations in metres
#' @param dtmc If weather is s not a data.frame, a SpatRaster object giving the resolution,
#' spatial extent, and projection of the weather data. Must give elevations in meters if
#' `altcorrect` > 0.
#' @param temp one of "air" or "leaf" indicating whether bioclim variables are constructed using air or leaf
#' temperatures. If `temp = leaf`, for grid cells where `reqhgt` is above vegetation, vertically averaged canopy
#' temperature is used.
#' @param hourly optional logical indicating whether to run the model in hourly or daily mode.
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction)
#' @param windhgt height above ground of wind speed data in weather
#' @param soilm optional vector of soil moisture values in upper 10 cm of the soil (calculated if not supplied)
#' @param dTmx optional maximum amount by which canopy or ground surface temperatures can exceed air temperatures.
#' Included to ensure model convergence.
#' @param maxiter optional integer indicating the maximum number of iterations (see details).
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' @param xyf optional spatial smoothing factor applied in calculation of surface
#' roughness and zero-plane displacement heights (see [modelin()])
#' @param pai_a an optional array of plant area index values above `reqhgt`.
#' @param tfact coefficient determining sensitivity of soil moisture to variation in topographic wetness (see [soilmdistribute()]).
#' @param surfwet an optional single numeric value of array of values specifying the proportion
#' of the canopy surface that should be treated as wet surface (modelled if not supplied)
#' @param slr slr an optional SpatRaster object of slope values (Radians). Calculated from
#' `dtm` if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' `dtm` if not supplied, but outer cells will be NA.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from `dtm` if not supplied, but outer cells will be NA.
#' @param twi optional SpatRaster object of topographic wetness index values.
#' Calculated from `dtm` if not supplied, but outer cells will be NA.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from `dtm` if not supplied, but outer cells will be NA.
#' @param lat optional central latitude of study area (removes tile effects when running in tiles)
#' @param long optional central longitude of study area (removes tile effects when running in tiles)
#' @return a multilayer SpatRast of the following:
#' \describe{
#'  \item{BIO1}{Mean of monthy median temperatures (degrees C)}
#'  \item{BIO2}{Mean diurnal temperature range (degrees C)}
#'  \item{BIO3}{Isothermality (BIO2/BIO7) (×100)}
#'  \item{BIO4}{Temperature Seasonality (standard deviation of monthly median temperatures × 100)}
#'  \item{BIO5}{Maximum temperature (degrees C)}
#'  \item{BIO6}{Minimum temperature (degrees C)}
#'  \item{BIO7}{Temperature Annual Range (BIO5-BIO6) (degrees C)}
#'  \item{BIO8}{Mean of monthly median temperatures in wettest three months (degrees C)}
#'  \item{BIO9}{Mean of monthly median temperatures in driest three months (degrees C)}
#'  \item{BIO10}{Mean of monthly median temperatures in warmest three months (degrees C)}
#'  \item{BIO11}{Mean of monthly median temperatures in coldest three months (degrees C)}
#'  \item{BIO12}{Mean of monthly soil moistures on the hottest, coldest and median temperature day of each month (m^3 / m^3)}
#'  \item{BIO13}{Maximum of monthly soil moistures on the hottest, coldest and median temperature day of each month (m^3 / m^3)}
#'  \item{BIO14}{Minimum of monthly soil moistures on the hottest, coldest and median temperature day of each month (m^3 / m^3)}
#'  \item{BIO15}{Soil moisture seasonality (Coefficient of Variation) on day in each month with median temperature}
#'  \item{BIO16}{Mean of monthly soil moistures on the hottest, coldest and median temperature day of the wettest three months (m^3 / m^3)}
#'  \item{BIO17}{Mean of monthly soil moistures on the hottest, coldest and median temperature day of the driest three months (m^3 / m^3)}
#'  \item{BIO18}{Mean of monthly soil moistures on the hottest, coldest and median temperature day of the warmest three months (m^3 / m^3)}
#'  \item{BIO19}{Mean of monthly soil moistures on the hottest, coldest and median temperature day of the coldest three months (m^3 / m^3)}
#' }
#' @details
#' To enhance computational efficiency the microclimate model is run for selected days only. Thus, to compute
#' mean annual temperature, the mean ambient temperature of each day in the input weather data is calculated,
#' the day with median temperatures in each month selected and the mean across months calculated. This is
#' not, strictly speaking, the same as the mean temperature, but differences are likely to be minor, and for each year of
#' data supplied, there is an approximately 30-fold gain in computational efficiency by calculating
#' BIO1 in this way. Similarly, to calculate maximum temperature (BIO5), the day of the year with the hottest ambient
#' temperature is selected, and microclimate temperatures calculated on this day only. This ignores
#' the possibility that on a slightly cooler, but sunnier day, microclimate temperatures may be hotter
#' hotter at certain locations. If `hourly = TRUE` all hours within a given day are selected and calculations
#' performed on hourly data. If If `hourly = FALSE` only the hours corresponding to times when hourly temperatures
#' are at their daily maximum and minimum and selected. This results in a c. 10-fold increase in computational
#' efficiency, but cannot pick out areas where terrain results in near-ground temperatures reaching
#' a maximum later in the afternoon than the peak in ambient temperature. If weather data for more than one year are supplied,
#' only one set of median, maximum and minimum monthly temperature data are selected representing an average
#' across years. Resultant, there is little computational penalty if providing data for multiple years
#' in comparison to one year of data.
#'
#' @import terra
#' @export
#' @examples
#' # Takes ~20 seconds to run
#' bioclim<-runbioclim(climdata, rainfall, tme = NA, reqhgt = 0.05, micropoint = NA,
#'                     vegp, soilc, dtmcaerth)
#' plot(bioclim[[1]]) # BIO1
runbioclim<-function(weather, precip, tme = NA, reqhgt = 0.05, micropoint = NA, vegp, soilc, dtm, dtmc = NA,
                     temp = "air", hourly = TRUE, altcorrect = 1, windhgt = 2, soilm = NA, dTmx = 25,
                     maxiter = 100, runchecks = TRUE, xyf = 1, pai_a = NA, tfact = 1.5, surfwet = NA,
                     slr = NA, apr = NA, hor = NA, twi = NA, wsa = NA, lat = NA, long = NA) {
  if (class(dtm) == "PackedSpatRaster") dtm<-rast(dtm) # raster template
  ll<-.latlongfromraster(dtm)
  if (class(weather) == "data.frame") {
    if (class(micropoint) == "logical") micropoint<-.biomicropoint(weather,precip,tme,vegp,soilc,reqhgt,windhgt,soilm=NA,dTmx,maxiter)
    # Calculate relevant quarters
    tme<-as.POSIXlt(weather$obs_time,tz="UTC")
    tc<-weather$temp
    prech<-rep(precip,each=24)/24
    agg<-stats::aggregate(prech,by=list(tme$mon),sum)$x
    wq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
    dq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # driest quarter
    agg<-stats::aggregate(tc,by=list(tme$mon),sum)$x
    hq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # hottest quarter
    cq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
    # Prepare model input
    if (hourly) {
      micro<-modelin(micropoint,vegp,soilc,dtm,runchecks,xyf)
      mout<-runmicro_hr(micro,reqhgt,pai_a,tfact,surfwet,slr,apr,hor,twi,wsa)
    } else {
      microd<-modelin_dy(micropoint,vegp,soilc,dtm,runchecks,xyf)
      mout<-runmicro_dy(microd,reqhgt,TRUE,pai_a,tfact,surfwet,slr,apr,hor,twi,wsa)
    }
  } else {
    tc<-apply(weather$temp,3,mean)
    sel<-.biosel(tme, tc)
    if (class(micropoint) == "logical") micropoint<-.biomicropoint(weather,precip,tme,vegp,soilc,reqhgt,windhgt,soilm=NA,dTmx,maxiter)
    # Calculate relevent quarters
    prech<-rep(apply(precip,3,mean),each=24)/24
    agg<-stats::aggregate(prech,by=list(tme$mon),sum)$x
    wq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
    dq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # driest quarter
    agg<-stats::aggregate(tc,by=list(tme$mon),sum)$x
    hq<-which.max(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # hottest quarter
    cq<-which.min(filter(agg, rep(1/3, 3), sides = 2, circular = TRUE)) # wettest quarter
    # Prepare model input
    if (hourly) {
      micro<-modelina(micropoint,vegp,soilc,dtm,dtmc,altcorrect,runchecks,xyf)
      mout<-runmicro_hr(micro,reqhgt,pai_a,tfact,surfwet,slr,apr,hor,twi,wsa)

    } else {
      microd<-modelina_dy(micropoint,vegp,soilc,dtm,dtmc,altcorrect,runchecks,xyf)
      mout<-runmicro_dy(microd,reqhgt,TRUE,pai_a,tfact,surfwet,slr,apr,hor,twi,wsa)
    }
  }
  if (temp == "air") {
    Tz<-mout$Tz
  } else Tz<-mout$tleaf
  # BIO1:  Annual Mean Temperature
  bio1<-.rast(apply(Tz[,,1:288],c(1,2),mean),dtm)
  # BIO2:  Mean Diurnal Range
  th<-array(Tz[,,1:288],dim=c(dim(Tz)[1:2],24,12))
  dtr<-apply(th,c(1,2,4),max)-apply(th,c(1,2,4),min)
  bio2<-.rast(apply(dtr,c(1,2),mean),dtm)
  # BIO4: Temperature Seasonality (standard deviation of monthly means ×100)
  tmth<-apply(th,c(1,2,4),mean)
  bio4<-.rast(apply(tmth,c(1,2),sd)*100,dtm)
  # BIO5: Max Temperature
  bio5<-.rast(apply(Tz,c(1,2),max),dtm)
  # BIO6: Min Temperature
  bio6<-.rast(apply(Tz,c(1,2),min),dtm)
  # BIO7: Temperature Annual Range (BIO5-BIO6)
  bio7<-bio5-bio6
  # BIO3: Isothermality (BIO2/BIO7) (×100)
  bio3<-(bio2/bio7)*100
  # BIO8: Mean Temperature of Wettest Quarter
  if (class(weather) == "data.frame") {
    tmeh<-as.POSIXlt(micropoint$weather$obs_time,tz="UTC")[1:288]
  } else {
    s<-which(is.na(micropoint) == FALSE)[1]
    xx<-micropoint[[s]]
    tmeh<-as.POSIXlt(xx$weather$obs_time,tz="UTC")[1:288]
  }
  sq<-.selquarter(tmeh,wq)
  bio8<-.rast(apply(Tz[,,sq],c(1,2),mean),dtm)
  # BIO9: Mean Temperature of Driest Quarter
  sq<-.selquarter(tmeh,dq)
  bio9<-.rast(apply(Tz[,,sq],c(1,2),mean),dtm)
  # BIO10: Mean Temperature of Warmest Quarter
  sq<-.selquarter(tmeh,hq)
  bio10<-.rast(apply(Tz[,,sq],c(1,2),mean),dtm)
  # BIO11: Mean Temperature of Coldest Quarter
  sq<-.selquarter(tmeh,cq)
  bio11<-.rast(apply(Tz[,,sq],c(1,2),mean),dtm)
  # BIO12: Mean soil moisture
  bio12<-.rast(apply(mout$soilm[,,1:288],c(1,2),mean),dtm)
  # BIO13: Max soil moisture
  bio13<-.rast(apply(mout$soilm,c(1,2),max),dtm)
  # BIO14: Min soil moisture
  bio14<-.rast(apply(mout$soilm,c(1,2),min),dtm)
  # BIO15: Soil moisture seasonality (Coefficient of Variation)
  sh<-array(mout$soilm[,,1:288],dim=c(dim(Tz)[1:2],24,12))
  smth<-apply(sh,c(1,2,4),mean)
  bio15<-.rast(apply(smth,c(1,2),sd)/apply(smth,c(1,2),mean),dtm)
  # BIO16: Mean soil moisture of Wettest Quarter
  sq<-.selquarter(tmeh,wq)
  bio16<-.rast(apply(mout$soilm[,,sq],c(1,2),mean),dtm)
  # BIO17: Mean soil moisture of Driest Quarter
  sq<-.selquarter(tmeh,dq)
  bio17<-.rast(apply(mout$soilm[,,sq],c(1,2),mean),dtm)
  # BIO18: Mean soil moisture of Warmest Quarter
  sq<-.selquarter(tmeh,hq)
  bio18<-.rast(apply(mout$soilm[,,sq],c(1,2),mean),dtm)
  # BIO19: Mean soil moisture of Coldest Quarter
  sq<-.selquarter(tmeh,cq)
  bio19<-.rast(apply(mout$soilm[,,sq],c(1,2),mean),dtm)
  # Stack bioclim variables
  bioclim<-c(bio1,bio2,bio3,bio4,bio5,bio6,bio7,bio8,bio9,bio10,
             bio11,bio12,bio13,bio14,bio15,bio16,bio17,bio18,bio19)
  nms<-paste0("BIO",c(1:19))
  names(bioclim)<-nms
  return(bioclim)
}
#' The function `runbioclim_big` runs the microclimate model and produces microclimate
#' equivalents of the 19 Worldclim bioclimate variables in tiles over large areas
#'
#' @param weather a data.frame or list of arrays of weather variables (as for [runpointmodel() or [runpointmodela()])
#' @param precip a vector or array of daily precipitation variables ((as for [runpointmodel() or [runpointmodela()]))
#' @param tme if weather is not a data.frame, an object of class POSIXlt giving the dates and times for each weather variable stored in the array. Set to NA if weather is a data,frame
#' @param reqhgt height for which temperatures are needed
#' @param vegp an object of class vegparams as returned by [vegpfromhab()]
#' @param soilc an object of class soilcharac as returned by [soilcfromtype()]
#' @param dtm a SpatRaster object of elevations in metres
#' @param dtmc If weather is s not a data.frame, a SpatRaster object giving the resolution,
#' spatial extent, and projection of the weather data. Must give elevations in meters if
#' `altcorrect` > 0.
#' @param temp one of "air" or "leaf" indicating whether bioclim variables are constructed using air or leaf
#' temperatures. If `temp = leaf`, for grid cells where `reqhgt` is above vegetation, vertically averaged canopy
#' temperature is used.
#' @param hourly optional logical indicating whether to run the model in hourly or daily mode.
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction)
#' @param windhgt height above ground of wind speed data in weather
#' @param soilm optional vector of soil moisture values in upper 10 cm of the soil (calculated if not supplied)
#' @param dTmx optional maximum amount by which canopy or ground surface temperatures can exceed air temperatures.
#' Included to ensure model convergence.
#' @param maxiter optional integer indicating the maximum number of iterations (see details).
#' @param runchecks optional logical indicating whether to call [checkinputs()] to run
#' @param xyf optional spatial smoothing factor applied in calculation of surface
#' roughness and zero-plane displacement heights (see [modelin()])
#' @param pai_a an optional array of plant area index values above `reqhgt`.
#' @param tfact coefficient determining sensitivity of soil moisture to variation in topographic wetness (see [soilmdistribute()]).
#' @param surfwet an optional single numeric value of array of values specifying the proportion
#' of the canopy surface that should be treated as wet surface (modelled if not supplied)
#' @param silent optional logical indicating whether to report on progress (default FALSE - progress reported).
#' @return a multilayer SpatRast of the following:
#' \describe{
#'  \item{BIO1}{Mean of monthy median temperatures (degrees C)}
#'  \item{BIO2}{Mean diurnal temperature range (degrees C)}
#'  \item{BIO3}{Isothermality (BIO2/BIO7) (×100)}
#'  \item{BIO4}{Temperature Seasonality (standard deviation of monthly median temperatures × 100)}
#'  \item{BIO5}{Maximum temperature (degrees C)}
#'  \item{BIO6}{Minimum temperature (degrees C)}
#'  \item{BIO7}{Temperature Annual Range (BIO5-BIO6) (degrees C)}
#'  \item{BIO8}{Mean of monthly median temperatures in wettest three months (degrees C)}
#'  \item{BIO9}{Mean of monthly median temperatures in driest three months (degrees C)}
#'  \item{BIO10}{Mean of monthly median temperatures in warmest three months (degrees C)}
#'  \item{BIO11}{Mean of monthly median temperatures in coldest three months (degrees C)}
#'  \item{BIO12}{Mean of monthly soil moistures on the hottest, coldest and median temperature day of each month (m^3 / m^3)}
#'  \item{BIO13}{Maximum of monthly soil moistures on the hottest, coldest and median temperature day of each month (m^3 / m^3)}
#'  \item{BIO14}{Minimum of monthly soil moistures on the hottest, coldest and median temperature day of each month (m^3 / m^3)}
#'  \item{BIO15}{Soil moisture seasonality (Coefficient of Variation) on day in each month with median temperature}
#'  \item{BIO16}{Mean of monthly soil moistures on the hottest, coldest and median temperature day of the wettest three months (m^3 / m^3)}
#'  \item{BIO17}{Mean of monthly soil moistures on the hottest, coldest and median temperature day of the driest three months (m^3 / m^3)}
#'  \item{BIO18}{Mean of monthly soil moistures on the hottest, coldest and median temperature day of the warmest three months (m^3 / m^3)}
#'  \item{BIO19}{Mean of monthly soil moistures on the hottest, coldest and median temperature day of the coldest three months (m^3 / m^3)}
#' }
#' @seealso [runbioclim()]
#' @details
#' The model is run in 100 pixel by 100 pixel tiles and the tiles are automatically merged back-togther
#' to create a single merged output. When climate data are provided as arrays tiling effects can be
#' present, and to avoid this, the tiles are run with a 10 pixel overlap on each side, and distance-weighted
#' blending undertaken to derive the final output.
#' @import terra
#' @export
runbioclim_big<-function(weather, precip, tme = NA, reqhgt = 0.05, vegp, soilc, dtm, dtmc = NA,
                         temp = "air", hourly = TRUE, altcorrect = 1, windhgt = 2, soilm = NA,
                         dTmx = 25, maxiter = 100, runchecks = TRUE, xyf = 1,
                         pai_a = NA, tfact = 1.5, surfwet = NA, silent = FALSE) {
  soilparams<-micropoint::soilparams # dirty fix
  # Unpack data
  up<-.unpack(dtm,vegp,soilc)
  dtm<-up$dtm
  vegp<-up$vegp
  soilc<-up$soilc
  .checkbiginputs(dtm,vegp,soilc) # check no odd NAs
  # Calculate universal variables
  if (silent == FALSE) cat("Computing terrain variables over whole area\n")
  slr<-terrain(dtm,v="slope",unit="radians")
  apr<-terrain(dtm,v="aspect",unit="radians")
  twi<-.topidx(dtm)
  hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) hor[,,i]<-.horizon(dtm,(i-1)*15)
  dsm<-dtm+vegp$hgt
  wsa<-.windsheltera(dsm,8,NA)
  # Determine rows and columns
  dms<-dim(dtm)[1:2]
  rws<-ceiling(dms[1]/100)
  cls<-ceiling(dms[2]/100)
  if (silent == FALSE) cat(paste0("Running model over ",rws*cls," tiles\n"))
  # Run point model
  ll<-.latlongfromraster(dtm)
  micropoint<-.biomicropoint(weather,precip,tme,vegp,soilc,reqhgt,windhgt,soilm=NA,dTmx,maxiter,ll$lat,ll$long)
  # Create list for storing data
  blst<-list()
  i<-1
  for (rw in 1:rws) {
    for (cl in 1:cls) {
      dtmi<-.croprast(dtm,rw,cl,100) # crop with overlap
      dc<-crop(dtmc,ext(dtmi))
      v<-as.vector(dtmi)
      v<-v[is.na(v)==F]
      v2<-as.vector(dc)
      v2<-v2[is.na(v2)==F]
      if (length(v)>1 & length(v2) > 1) {
        slri<-crop(slr,ext(dtmi))
        apri<-crop(apr,ext(dtmi))
        twii<-crop(twi,ext(dtmi))
        hori<-as.array(crop(.rast(hor,dtm),ext(dtmi)))
        wsai<-as.array(crop(.rast(wsa,dtm),ext(dtmi)))
        vegpi<-.vegcrop(vegp,dtmi)
        soilci<-.soilcrop(soilc,dtmi)
        if (silent==FALSE) cat(paste0("Running model for tile ",rw," ",cl,"\n"))
        blst[[i]]<-runbioclim(weather,precip,tme,reqhgt,micropoint,vegpi,soilci,dtmi,dtmc,temp,hourly,
                              altcorrect,windhgt,soilm,dTmx,maxiter,runchecks,xyf,
                              pai_a,tfact,surfwet,slri,apri,hori,twii,wsai,ll$lat,ll$long)
        i<-i+1
      }
    } # end column
  } # end row
  bsrc<-sprc(blst)
  mos<-mosaic(bsrc)
  return(mos)
}

