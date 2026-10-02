## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = FALSE
)


## ----eval=FALSE---------------------------------------------------------------
# require(devtools)
# install_github("ilyamaclean/microclimf", ref = "dev")


## ----eval=FALSE---------------------------------------------------------------
# library(microclimf)
# library(terra)
# # Runs point microclimate model with inbuilt datasets
# micropoint <- runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc)
# # Subset point model outputs
# micropoint_mx <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")
# micropoint_mn <- subsetpointmodel(micropoint, tstep = "month", what = "tmin")
# # Run grid model 5 cm above ground with subset values and inbuilt datasets (takes ~20 seconds)
# mout_mx <- rungridmodel(micropoint_mx, reqhgt = 0.05, dtmcaerth, vegp, soilc)
# mout_mn <- rungridmodel(micropoint_mn, reqhgt = 0.05, dtmcaerth, vegp, soilc)
# # Plot air temperatures on hottest hour in micropoint (2017-06-20 13:00:00 UTC)
# plot(mout_mx$Tz[[134]])
# # Plot air temperatures on coldest hour in micropoint (2017-01-03 08:00:00 UTC)
# plot(mout_mn$Tz[[9]])
# # Plot mean of monthly max and min
# mairt <- mean((mout_mn$Tz + mout_mx$Tz) / 2)
# plot(mairt)


## -----------------------------------------------------------------------------
library(microclimf)
head(climdata)


## ----eval = FALSE-------------------------------------------------------------
# ?PFTparams


## ----fig.show='hold'----------------------------------------------------------
library(terra)

names(vegp)

plot(rast(vegp$habitat), main = "Habitat class")
plot(rast(vegp$pai)[[1]], main = "Jan PAI")

paiarray <- as.array(rast(vegp$pai))
plot(apply(paiarray, 3, mean, na.rm = TRUE),
     type = "l", main = "Seasonal variation in PAI")

plot(rast(vegp$hgt), main = "Vegetation height")
plot(rast(vegp$x), main = "Leaf angle distribution")
plot(rast(vegp$gsmax), main = "Max. stomatal conductance")
plot(rast(vegp$clump)[[1]], main = "Jan canopy clumping factor")
plot(rast(vegp$leafr), col = gray(0:255/255), main = "Leaf reflectance")
plot(rast(vegp$leaft), col = gray(0:255/255), main = "Leaf transmittance")


## ----fig.show='hold'----------------------------------------------------------
attributes(soilc)
plot(rast(soilc$soiltype), main = "Soiltype") # Clay loam throughout
plot(rast(soilc$groundr), col=gray(0:255/255), main = "Soil reflectance")


## -----------------------------------------------------------------------------
soilparamstable


## ----eval=FALSE---------------------------------------------------------------
# ?soilparamstable


## ----eval=FALSE---------------------------------------------------------------
# # Run the point model
# micropoint <- runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc)


## ----eval=FALSE---------------------------------------------------------------
# attributes(micropoint)


## ----eval=FALSE---------------------------------------------------------------
# micropoint <- runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc)
# mod <- micropoint$model
# tme <- as.POSIXct(micropoint$weather$obs_time)
# par(mar=c(5,5,3,3))
# plot(mod$Tground ~ tme, type="l", ylim = c(-5, 45), col = rgb(1,0,0,0.5), xlab = "Month", ylab = "Temperature") # temperature of ground surface
# par(new = TRUE)
# plot(mod$Tcanopy ~ tme, type="l", ylim = c(-5, 45), col = rgb(0,0.5,0.5,0.5), xlab = "", ylab = "")


## ----eval=FALSE---------------------------------------------------------------
# micropoint <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")


## ----eval=FALSE---------------------------------------------------------------
# 


## ----eval=FALSE---------------------------------------------------------------
# # Takes ~10 seconds to run
# micro <- rungridmodel(micropoint, reqhgt = 0.0, dtmcaerth, vegp, soilc)
# # Surface soil moisture (2017-01-30 19:00:00 UTC)
# plot(micro$soilm[[20]], col = rev(map.pal("viridis")))


## ----eval=FALSE---------------------------------------------------------------
# # Downward direct and diffuse shortwave radiation (2017-06-20 13:00:00 UTC)
# plot(micro$Rdirdown[[134]])
# plot(micro$Rdifdown[[134]])
# # Upward shortwave radiation
# plot(micro$Rswup[[134]])


## ----eval=FALSE---------------------------------------------------------------
# # Downward and upward longwave radiation (2017-06-20 13:00:00 UTC)
# plot(micro$Rlwdown[[134]])
# plot(micro$Rlwup[[134]])


## ----eval=FALSE---------------------------------------------------------------
# # Takes ~10 seconds to run
# micro <- rungridmodel(micropoint, reqhgt = 0.1, dtmcaerth, vegp, soilc)
# # Wind speed 10 cm above ground (2017-01-30 19:00:00 UTC, direction from south)
# plot(micro$windspeed[[20]])


## ----eval=FALSE---------------------------------------------------------------
# # Run point model and select hottest and coldest days of each month
# micropoint <- runpointmodel(climdata,0.0,dtmcaerth,vegp,soilc)
# micropoint_mx <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")
# micropoint_mn <- subsetpointmodel(micropoint, tstep = "month", what = "tmin")
# # Run grid model and get absolute max and min
# mout_mx <- rungridmodel(micropoint_mx, reqhgt = 0.0, dtmcaerth, vegp, soilc)
# mout_mn <- rungridmodel(micropoint_mn, reqhgt = 0.0, dtmcaerth, vegp, soilc)
# # Plot ground surface temperatures
# plot(max(mout_mx$Tz, na.rm = TRUE)) # hottest
# plot(min(mout_mn$Tz, na.rm = TRUE)) # coldest


## ----eval=FALSE---------------------------------------------------------------
# # 10 cm above ground - micropoint_mx and micropoint_mn recycled
# mout_mx <- rungridmodel(micropoint_mx, reqhgt = 0.1, dtmcaerth, vegp, soilc)
# mout_mn <- rungridmodel(micropoint_mn, reqhgt = 0.1, dtmcaerth, vegp, soilc)
# # Plot air temperatures 10cm above ground
# plot(max(mout_mx$Tz, na.rm = TRUE)) # hottest
# plot(min(mout_mn$Tz, na.rm = TRUE)) # coldest


## ----eval=FALSE---------------------------------------------------------------
# # Plot leaf temperatures 10cm above ground
# plot(max(mout_mx$tleaf, na.rm = TRUE)) # hottest
# plot(min(mout_mn$tleaf, na.rm = TRUE)) # coldest


## ----eval=FALSE---------------------------------------------------------------
# # Run point model and select hottest days of each month
# micropoint10 <- runpointmodel(climdata, reqhgt = -0.1, dtmcaerth, vegp, soilc)
# micropoint50 <- runpointmodel(climdata, reqhgt = -0.5, dtmcaerth, vegp, soilc)
# micropoint10 <- subsetpointmodel(micropoint10, tstep = "month", what = "tmax")
# micropoint50 <- subsetpointmodel(micropoint50, tstep = "month", what = "tmax")
# # Run grid model and get temperatures 10 cm and 50 cm below ground
# mout10 <- rungridmodel(micropoint10, reqhgt = -0.1, dtmcaerth, vegp, soilc)
# mout50 <- rungridmodel(micropoint50, reqhgt = -0.5, dtmcaerth, vegp, soilc)
# # Plot below ground temperatures
# plot(mout10$Tz[[134]], range = c(13, 38)) # 10 cm below
# plot(mout50$Tz[[134]], range = c(13, 38)) # 50 cm below


## ----eval=FALSE---------------------------------------------------------------
# library(terra)
# # Create a dummy 5 x 5 array of the climate variables
# #  -- Convert array to raster --
# .rast <- function(x, template) {
#   r <- rast(x)
#   ext(r) <- ext(template)
#   crs(r) <- crs(template)
#   r
# }
# # -- replicate each variable 5 x 5 times --
# .ta <- function(x, template, xdim = 5, ydim = 5) {
#   a <- array(rep(x, each = xdim * ydim),
#              dim = c(ydim, xdim, length(x)))
#   .rast(a, template)
# }
# # -- Create dummy list of rasters
# dtm <- rast(dtmcaerth)
# vars <- names(climdata)[2:10]
# climarrayr <- lapply(climdata[vars], .ta, template = dtm)
# # Get times corresponding to each layer
# tme <- as.POSIXlt(climdata$obs_time, tz="UTC")
# # Create coarse-resolution dtm matching resolution of climate data
# dtmc <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE)
# # Run point model array (takes ~ one minute)
# micropointa <- runpointmodel(climarrayr, reqhgt = 0.05, dtm, vegp, soilc,
#                              tme, dtmc, cores = "auto")
# # Subset point model using defaults (monthly tmax)
# micropointa <- subsetpointmodel(micropointa)
# # Run microclimate model with no altitude correction
# mout <- rungridmodel(micropointa, reqhgt = 0.05, dtm, vegp, soilc, altcorrect = 0)
# # plot maximum temperature
# plot(mout$Tz[[134]])


## ----eval=FALSE---------------------------------------------------------------
# # Create coarse resolution vegp and soilc layers
# vegpc <- lapply(vegp, function(x) {
#   x <- unwrap(x)
#   aggregate(x, fact = 10, fun = mean, na.rm = TRUE)
# })
# vegpc$habitat <- aggregate(unwrap(vegp$habitat), fact = 10,
#                            fun = "modal", na.rm = TRUE)
# soilcc <- lapply(soilc, function(x) {
#   x <- unwrap(x)
#   aggregate(x, fact = 10, fun = mean, na.rm = TRUE)
# })
# soilcc$soiltype <- aggregate(unwrap(soilc$soiltype), fact = 10,
#                            fun = "modal", na.rm = TRUE)
# # Re-use tme, climarrayr and dtmc created in code above
# mout <- runpointmodelasgrid(climarrayr, reqhgt = 0.05, dtmc,
#                             vegpc, soilcc, tme, cores = "auto")
# # Find hottest day
# hotday <- which.max(climdata$temp)
# plot(mout$Tz[[hotday]])


## ----eval=FALSE---------------------------------------------------------------
# # Download example data from Zenodo
# url <- "https://zenodo.org/records/23091040/files/modeldatav2.zip"
# pathout<-"C:/Temp/tiles/"
# dir.create(pathout)
# setwd(pathout)
# download.file(url, "modeldatav2.zip")
# unzip("modeldatav2.zip")
# # Read in spatial data
# big_vegp <- readRDS("vegp_big.RDS")
# big_soilc <-readRDS("soilc_big.RDS")
# dtm <- rast("dem.tif")
# # Run and subset point model using inbuilt climate dataset (~15 seconds)
# micropoint <- runpointmodel(climdata, reqhgt = 0.05, dtm, big_vegp, big_soilc)
# micropoint <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")
# # Run the model in tiles
# rungridmodelbig(micropoint, reqhgt = 0.05, dtm, big_vegp, big_soilc,
#                 cores = "auto")


## ----eval=FALSE---------------------------------------------------------------
# # Use inbuilt datasets
# vegp <- microclimf::vegp
# soilc <- microclimf::soilc
# bioclim <- runbioclim(climdata, 0.05, dtmcaerth, vegp, soilc, temp = "air")
# # Mean temperature of the coldest quarter:
# plot(bioclim[[11]])
# # Volumetric soil water fraction of driest quarter
# plot(bioclim[[17]], col = rev(map.pal("viridis")))

