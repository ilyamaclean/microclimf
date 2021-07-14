globalVariables("globclim")
.PAIforayear <- function(habitat, lat, long) {
  laigaus <- function(minlai, maxlai, pkday, dhalf, yr) {
    diy <- 365
    sdev <- 0.0082 * dhalf^2 + 0.0717 * dhalf + 13.285
    difv <- maxlai - minlai
    x<-c(-diy:diy)
    y <- 1 / (sdev * sqrt(2 * pi)) * exp(-0.5 * (((x - 0) / sdev) ^ 2))
    y[(diy + ceiling(0.5 * diy)):(2 * diy + 1)] <- y[(diy - ceiling(0.5 * diy)):diy]
    st <- diy + 1 - pkday
    y <- y[st:(st + diy - 1)]
    x <- c(1:diy)
    x <- c(0, x, c(366:375))
    y <- c(y[diy], y, y[1:10])
    sel <-c(0:15) * 25 + 1
    x<-x[sel]
    y<-y[sel]
    tme <- as.POSIXct((x * 24 * 3600), origin = paste0(yr - 1,"-12-31 12:00"), tz = "GMT")
    xy <- stats::spline(tme, y, n = diy * 24 + 241)
    tme2 <- as.POSIXlt(xy$x, origin = "1970-01-01 00:00", tz = "GMT")
    sel <- which(tme2$year + 1900 == yr)
    y <- xy$y[sel]
    dify <- max(y) - min(y)
    y <- y * (difv / dify)
    y <- y + minlai - min(y)
    return(y)
  }
  long<- ifelse(long > 180.9375, long - 360, long)
  long<- ifelse(long < -179.0625, long + 360, long)
  lat<- ifelse(lat< -89.49406, -89.49406, lat)
  lat<- ifelse(lat> 89.49406, 89.49406, lat)
  ll<- SpatialPoints(data.frame(x = long, y = lat))
  mmonth <-c(16, 45.5, 75, 105.5, 136, 166.5, 197, 228, 258.5, 289, 319.5, 350)
  e <- extent(c(-179.0625, 180.9375, -89.49406, 89.49406))
  clim <- rep(NA,5)
  for (i in 1:5) {
    r <- raster(globclim[,,i])
    extent(r) <- e
    clim[i] <- extract(r, ll)
  }
  wgts <- function(x1, x2, ll, lmn, lmx) {
    ll <- ifelse(ll < lmn, lmn, lat)
    ll <- ifelse(ll > lmx, lmx, lat)
    w <- 1 - (abs(ll - lmn)  / (abs(ll - lmn)  + abs(ll - lmx)))
    y <- w * x1 + (1 - w) * x2
    y
  }
  # By habitat type
  if (habitat == 1) { # Evergreen needleleaf forest
    h2 <- 74.02 + 5.35 * clim[1]
    h1 <-  203.22 - 35.63 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 50, 50, hperiod)
    p2 <- 216.71 - 2.65 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <- 2.33 + 0.0132 * clim[1]
    } else maxlai <- 2.33 + 0.0132 * 20
    minlai <- 1.01
  }
  if (habitat == 2) { # Evergreen broadleaf forest
    hperiod <-  154.505 + 2.040 * clim[1]
    hperiod <- ifelse(hperiod < 50, 50, hperiod)
    peakdoy <- peakdoy <- mmonth[round(clim[5], 0)]
    maxlai <- 1.83 + 0.22 * log(clim[3])
    minlai <- (-1.09) + 0.4030 * log(clim[3])
    minlai <- ifelse(minlai < 1, 1, minlai)
  }
  if (habitat == 3) { # Deciduous needleleaf forest
    h2 <- 51.18 + 3.77  * clim[1]
    h1 <- 152
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    p2 <- 204.97 - 1.08 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <-  2.62 + 0.05 * clim[1]
    } else maxlai <- 2.62 + 0.05 * 20
    minlai <- 0.39
  }
  if (habitat == 4) { # Deciduous broadleaf forest
    h2 <- 47.6380 + 2.9232 * clim[1]
    h1 <- 220.06 - 79.19 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 32.5, 32.5, hperiod)
    p2 <- 209.760 - 1.208 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <- 3.98795 + 0.03330 * clim[1]
    } else maxlai <- 3.98795 * 0.03330 * 20
    minlai <- 0.4808
  }
  if (habitat == 5) { # Mixed forest
    h2 <- 74.02 + 5.35 * clim[1]
    h1 <-  203.22 - 35.63 * clim[4]
    hperiod1 <- wgts(h1, h2, abs(lat), 0, 20)
    h2 <- 51.18 +  3.77  * clim[1]
    h1 <-  152
    hperiod2 <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- (hperiod1 + hperiod2) / 2
    hperiod <- ifelse(hperiod < 30.5, 30.5, hperiod)
    p2 <- 216.71 - 2.65 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy1 <- wgts(p1, p2, abs(lat), 0, 30)
    p2 <-  204.97 - -1.08 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    peakdoy2 <- wgts(p1, p2, abs(lat), 0, 30)
    peakdoy <- (peakdoy1 + peakdoy2) / 2
    if (clim[1] <= 20) {
      maxlai1 <- 2.33 + 0.0132 * clim[1]
      maxlai2 <-  2.62 + 0.05 * clim[1]
      maxlai <- (maxlai1 + maxlai2) / 2
    } else maxlai <- 3.107
    minlai <- 0.7
  }
  if (habitat == 6) { # Closed shrublands
    h2 <- 33.867 + 6.324 * clim[1]
    h1 <-  284.20 - 102.51 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 30.5, 30.5, hperiod)
    p2 <- 223.55 - 3.125 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    maxlai <- 2.34
    minlai <- -0.4790 + 0.1450 * log(clim[3])
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
  }
  if (habitat == 7) { # Open shrublands
    h2 <- 8.908 + 4.907 * clim[1]
    h1 <-  210.09 - 28.62 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 38.3, 38.3, hperiod)
    p2 <- 211.7 - 4.085 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    maxlai <- -0.7206 + 0.272 * log(clim[3])
    minlai <- -0.146 +  0.059 * log(clim[3])
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
  }
  if (habitat == 8) { # Woody savannas
    hperiod1 <-  47.6380 + 2.9232 * clim[1]
    hperiod1 <- ifelse(hperiod1 < 32.5, 32.5, hperiod1)
    hperiod2 <- 71.72 + 3.012 * clim[1]
    h1 <- (hperiod1 + hperiod2) / 2
    h2 <- 282.04 - 92.28 * clim[4]
    h2 <- ifelse(hperiod1 < 31.9, 31.9, hperiod1)
    hperiod <- wgts(h1, h2, abs(lat), 25, 35)
    peakdoy1 <- 209.760 - 1.208 * clim[1]
    peakdoy1 <- ifelse(peakdoy1 > 244, 244, peakdoy1)
    if (lat < 0)  peakdoy1 <- ( peakdoy1 + diy / 2)%%diy
    peakdoy2 <- 211.98 - 3.4371 * clim[1]
    peakdoy2 <- ifelse(peakdoy2 > 244, 244, peakdoy2)
    if (lat < 0)  peakdoy2 <- (peakdoy2 + diy / 2)%%diy
    p2 <- (peakdoy1 + peakdoy2) / 2
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 40)
    if (clim[1] <= 20) {
      maxlai1 <- 3.98795 + 0.03330 * clim[1]
      maxlai2 <- 1.0532 * 0.016 * clim[1]
    } else {
      maxlai1 <- 3.98795 + 0.03330 * 20
      maxlai2 <- 1.0532 * 0.016 * 20
    }
    mx2 <- (maxlai1 + maxlai2) / 2
    minlai1 <- 0.4808
    minlai2 <- 0.0725 * 0.011 * clim[1]
    mn2 <- (minlai1 + minlai2) / 2
    mx1 <- 1.298 + 0.171 * log(clim[3])
    mn1 <- -2.9458 + 0.5889 * log(clim[3])
    maxlai <- wgts(mx1, mx2, abs(lat), 10, 40)
    minlai <- wgts(mn1, mn2, abs(lat), 10, 40)
    minlai <- ifelse(minlai < 0.0362, 0.0362, minlai)
  }
  if (habitat == 9 | habitat == 10 | habitat == 11) {  # Grasslands
    h2 <- 71.72 + 3.012 * clim[1]
    h1 <- 269.22 -  89.79 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 31.9, 31.9, hperiod)
    p2 <- 211.98 - 3.4371 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    mx2 <- 1.48
    mn2 <- 0.0725 * 0.011 * clim[1]
    mx1 <- 0.1215 + 0.2662 * log(clim[3])
    mn1 <- 0.331 + 0.0575 * log(clim[3])
    maxlai <- wgts(mx1, mx2, abs(lat), 10, 40)
    minlai <- wgts(mn1, mn2, abs(lat), 10, 40)
    minlai <- ifelse(minlai < 0.762, 0.762, minlai)
  }
  if (habitat == 12) {  # Permanent wetlands
    h2 <- 76 + 4.617 * clim[1]
    h1 <- 246.68 - 66.82 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 40, 40, hperiod)
    p2 <- 219.64 - 2.793 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    maxlai <- -0.1782 + 0.2608 * log(clim[3])
    maxlai <- ifelse(maxlai < 1.12, 1.12, maxlai)
    minlai <-  -0.1450 + 0.1440 * log(clim[3])
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
  }
  if (habitat == 13) { # Croplands
    h2 <- 54.893 +  1.785 * clim[1]
    h1 <- 243 - 112.18 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 10, 30)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 212.95 - 5.627 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 0, 30)
    if (clim[1] <= 20) {
      maxlai <- 3.124 - 0.0886 * clim[1]
    } else maxlai <- 3.124 - 0.0886 * 20
    maxlai <- ifelse(maxlai > 3.14, 3.14, maxlai)
    minlai <- 0.13
  }
  if (habitat == 14) { # Urban and built-up
    h2 <- 66.669 +  5.618 * clim[1]
    h1 <- 283.44 - 86.11 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 215.998 - 4.2806 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 30)
    if (clim[1] <= 20) {
      maxlai <- 1.135 - 0.0244 * clim[1]
    } else maxlai <- 1.135 - 0.0244 * 20
    maxlai <- ifelse(maxlai > 1.15, 1.15, maxlai)
    minlai <- 0.28
  }
  if (habitat == 15) { # Cropland/Natural vegetation mosaic
    h2 <- 29.490 +  8.260 * clim[1]
    h1 <- 326.46 - 161.70 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 10, 30)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 210.867 - 3.5464 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 30)
    if (clim[1] <= 20) {
      maxlai <- 3.5485 - 0.09481 * clim[1]
    } else maxlai <- 3.5485 - 0.09481 * 20
    maxlai <- ifelse(maxlai > 3.14, 3.14, maxlai)
    if (clim[1] <= 20) {
      minlai <- -0.072815 - 0.044546 * clim[1]
    } else minlai <- -0.072815 - 0.044546 * 20
    minlai <- ifelse(minlai < 0.001, 0.001, minlai)
  }
  if (habitat == 16) { # Barren or sparsely vegetated
    h2 <- 80.557 +  6.440 * clim[1]
    h1 <- 344.65 -  -191.94 * clim[4]
    hperiod <- wgts(h1, h2, abs(lat), 0, 20)
    hperiod <- ifelse(hperiod < 43.5, 43.5, hperiod)
    p2 <- 236.0143 - 3.4726 * clim[1]
    p2 <- ifelse(p2 > 244, 244, p2)
    if (lat < 0) p2 <- (p2 + diy / 2)%%diy
    p1 <- mmonth[round(clim[5], 0)]
    peakdoy <- wgts(p1, p2, abs(lat), 10, 30)
    maxlai <- -0.05491 + 0.05991 * log(clim[4])
    maxlai <- ifelse(maxlai < 0.81, 0.81, maxlai)
    minlai <- 0.08
  }
  lai <- laigaus(minlai, maxlai, peakdoy, hperiod, 2000)
  if (habitat ==  "Short grasslands" | habitat == 10) lai <- lai / 2
  yhr<-round(mmonth*24,0)
  lai<-lai[yhr]
  return(lai)
}

.paifromhabitat <- function(habitat, lat, long, tme) {
  yr<-unique(tme$year+1900)
  pai<-0
  lai<-.PAIforayear(habitat, lat, long)
  for (i in yr) {
    sel<-which(tme$year+1900==i)
    tme2<-tme[sel]
    mth<-unique(tme2$mon+1)
    pai1<-lai[mth]
    pai<-c(pai,pai1)
  }
  return(pai[-1])
}
.onehab<-function(habitat) {
  if (habitat == 1) {  # Evergreen needleleaf forest
    hgt<-15 # Vegetation height
    x<-0.4  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.25 # Leaf reflectivity (shortwave)
    leafd<-0.01 # Leaf width (m)
  }
  if (habitat == 2) { # Evergreen broadleaf forest
    hgt<-20 # Vegetation height
    x<-1.2  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.35 # Leaf width (m)
  }
  if (habitat == 3) { # Deciduous needleleaf forest
    hgt<-10 # Vegetation height
    x<-0.4  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.01  # Leaf width (m)
  }
  if (habitat == 4) { # Deciduous broadleaf forest
    hgt<-15 # Vegetation height
    x<-1.2  # Campbell x
    gsmax<-0.23 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.07 # Leaf width (m)
  }
  if (habitat == 5) { # Mixed forest
    hgt<-10 # Vegetation height
    x<-0.8  # Campbell x
    gsmax<-0.28 # Maximum stomatal conductance
    leafr<-0.28 # Leaf reflectivity (shortwave)
    leafd<-0.04 # Leaf width (m)
  }
  if (habitat == 6) { # Closed shrublands
    hgt<-2 # Vegetation height
    x<-1  # Campbell x
    gsmax<-0.35 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.04 # Leaf width (m)
  }
  if (habitat == 7) { # Open shrublands
    hgt<-1.5 # Vegetation height
    x<-0.7  # Campbell x
    gsmax<-0.35 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.04 # Leaf width (m)
  }
  if (habitat == 8) { # Woody savannas
    hgt<-3 # Vegetation height
    x<-0.7  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.35 # Leaf reflectivity (shortwave)
    leafd<-0.03 # Leaf width (m)
  }
  if (habitat == 9) { #Savannas
    hgt<-1.5 # Vegetation height
    x<-0.15  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.35 # Leaf reflectivity (shortwave)
    leafd<-0.01 # Leaf width (m)
  }
  if (habitat == 10) { # Short grasslands
    hgt<-0.25 # Vegetation height
    x<-0.15  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.35 # Leaf reflectivity (shortwave)
    leafd<-0.01 # Leaf width (m)
  }
  if (habitat == 11) { # Tall grasslands
    hgt<-1.5 # Vegetation height
    x<-0.15  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.35 # Leaf reflectivity (shortwave)
    leafd<-0.01 # Leaf width (m)
  }
  if (habitat == 12) { # Permanent wetlands
    hgt<-0.5 # Vegetation height
    x<-1.4  # Campbell x
    gsmax<-0.55 # Maximum stomatal conductance
    leafr<-0.5 # Leaf reflectivity (shortwave)
    leafd<-0.09 # Leaf width (m)
  }
  if (habitat == 13) { # Croplands
    hgt<-0.5 # Vegetation height
    x<-0.2  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.02 # Leaf width (m)
  }
  if (habitat == 14) { # Urban and built-up
    hgt<-1.5 # Vegetation height
    x<-1  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.04 # Leaf width (m)
  }
  if (habitat == 15) { # Cropland/Natural vegetation mosaic
    hgt<-1 # Vegetation height
    x<-0.5  # Campbell x
    gsmax<-0.3 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.03 # Leaf width (m)
  }
  if (habitat == 16) { # Barren or sparsely vegetated
    hgt<-0.15 # Vegetation height
    x<-0.6  # Campbell x
    gsmax<-0.33 # Maximum stomatal conductance
    leafr<-0.3 # Leaf reflectivity (shortwave)
    leafd<-0.015 # Leaf width (m)
  }
  return(list(hgt=hgt,x=x,gsmax=gsmax,leafr=leafr,leafd=leafd))
}
#' Generate vegpetation paramaters form habitat type
#'
#' The function `vegpfromhab` generates an object of class vegparams from
#' a RasterLayer object of habitat types
#'
#' @param habitats a data.frame of weather variables (see details)
#' @param hgts an optional raster of vegetation heights. Estimated from habitat type if not provided.
#' @param pai an optional array of plant area index values. Estimated at monthly intervals
#' from habitat type, with seasonal variation dtermined form location and date if not provided.
#' @param lat latitude in decimal degrees. Only needed if `pai` not provided.
#' @param long longitude in decimal degrees. Only needed if `pai` not provided.
#' @param tme POSIXlt object of dates. Only needed if `pai` not provided (see details).
#' @return an object of class vegparams - a list with the following elements:
#' @return `pai` an array of monthly plant area index values (see details).
#' @return `hgt` a raster if vegetation heights (m)
#' @return `x` a raster of ratios of vertical to horizontal projections of leaf foliage
#' @return `gsmax` a raster of maximum stomatal conductances (mol / m^2 / s)
#' @return `leafr` a raster of leaf reflectance values (to shortwave radiation)
#' @return `clump` a raster indicating the degree of canopy clumpiness, here set to 0 (vegetation not clumped)
#' @return `leafd` a raster of mean leaf widths (m)
#'
#' @details
#' This function estimates the vegetation parameters needed to run microclimf from habitat types.
#' Plant area index values represent the combined one sided woody and green vegetation
#' plant area per unit ground area. If not provided, then approximated
#' from habitat type, location and date. The procedure is based on calibration
#' against MODIS-derived estimates, accounting for regional climate. An inbuilt dataset
#' of regional rainfall and temperature is included with the package. Monthly
#' values for each unique month in `tme` are returned. Note that values are
#' assumed spatially constant across a given habitat type, which is unlikely
#' to be the case in reality. Likewise, if vegetation height values are not provided,
#' these are estimated form habitat type and assumed constant within that habitat type.
#' Habitat types should be expressed as integers as follows:
#' (1) for Evergreen needleleaf forest,
#' (2) for Evergreen broadleaf forest,
#' (3) for Deciduous needleleaf forest,
#' (4) for Deciduous broadleaf forest,
#' (5) for Mixed forest,
#' (6) for Closed shrubland,
#' (7) for Open shrubland,
#' (8) for Woody savanna,
#' (9) for Savanna,
#' (10) for Short grassland,
#' (11) for Tall grassland,
#' (12) for Permanent wetland,
#' (13) for Cropland,
#' (14) for Urban and built-up,
#' (15) for Cropland / Natural vegetation mosaic and
#' (16) for Barren or sparsely vegetated

#' @export
#' @import raster sp
#'
#' @examples
#' library(raster)
#' tme<-as.POSIXlt(c(0:8783)*3600,origin="2000-01-01 00:00", tz = "GMT")
#' veg<-vegpfromhab(habitats,lat=50,long=-5,tme=tme)
#' plot(raster(veg$pai[,,1]), main = "Jan PAI")
#' plot(veg$hgt, main = "Vegetation height")
#' plot(veg$x, main = "Leaf angle coefficient")
#' plot(vegp$gsmax, main = "Maximum stomatal conductance")
#' plot(veg$leafr, main = "Leaf reflectance")
#' plot(veg$clump, main = "Canopy clumping factor")
#' plot(veg$leafd, main = "Leaf diamater")
vegpfromhab <- function(habitats, hgts = NA, pai = NA, lat, long, tme) {
  .poparray<-function(a,sel,v) {
    for (i in 1:length(v)) {
      m<-a[,,i]
      m[sel]<-v[i]
      a[,,i]<-m
    }
    a
  }
  # unique habitats
  m<-.is(habitats)
  uh<-unique(as.vector(m))
  uh<-uh[is.na(uh)==F]
  # pai test
  pte<-mean(pai,na.rm=T)
  hte<-mean(hgts,na.rm=T)
  # Create blank array for pai
  if (is.na(pte)) {
    paii<-.paifromhabitat(1, lat, long, tme)
    pai<-array(NA,dim=c(dim(m),length(paii)))
  }
  # Create blank rasters
  x<-m; gsmax<-m; leafr<-m; leafd<-m; hgt<-m
  clump<-habitats*0
  for (i in uh) {
    sel<-which(m==i)
    if (is.na(pte)) {
      paii<-.paifromhabitat(i, lat, long, tme)
      pai<-.poparray(pai,sel,paii)
    }
    vegi<-.onehab(i)
    x[sel]<-vegi$x
    gsmax[sel]<-vegi$gsmax
    leafr[sel]<-vegi$leafr
    leafd[sel]<-vegi$leafd
    hgt[sel]<-vegi$hgt
  }
  # Convert to rasters
  if (is.na(hte)) {
    hgt<-raster(hgt,template=habitats)
  } else hgt<-hgts
  x<-raster(x,template=habitats)
  gsmax<-raster(gsmax,template=habitats)
  leafr<-raster(leafr,template=habitats)
  clump<-habitats*0
  leafd<-raster(leafd,template=habitats)
  vegp<-list(pai=pai,hgt=hgt,x=x,gsmax=gsmax,leafr=leafr,clump=clump,leafd=leafd)
  class(vegp)<-"vegparams"
  return(vegp)
}
