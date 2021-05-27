paifromhab <- function(habitat_r, lat, long, year, meantemp = NA, cvtemp = NA,
                        rainfall = NA, cvrain = NA, wetmonth = NA) {
  warning("Users of `microclimf` are urged to express caution if relying upon vegetation parameters derived from this function. Best guess values such as those provided here are unlikely to capture real dynamics, and empirical measurements of vegetation are encouraged instead.")
  laigaus <- function(minlai, maxlai, pkday, dhalf, yr) {
    diy <- 365
    diy_array <- array(rep(365, dim(minlai)[1] *dim(minlai)[2]),
                       dim = c(dim(minlai)[1], dim(minlai)[2]))

    sdev <- 0.0082 * dhalf^2 + 0.0717 * dhalf + 13.285
    difv <- maxlai - minlai
    x <- aperm(array(rep(c(-diy:diy), dim(sdev)[1] * dim(sdev)[2]),
             dim = c(length(c(-diy:diy)), dim(sdev)[1], dim(sdev)[2])), c(2,3,1))
    # Generating y was one line in laifromhabitat(), but here I've broken it into
    # two commands because needs to be applied to arrays
    count <- 0
    x_sdev <- array(apply(x, c(3), function(i) {
      count <<- count + 1
      exp(-0.5 * (((i - 0) / sdev) ^ 2)) # 2nd half of y command
      # exp(-0.5 * (((x - 0) / sdev) ^ 2))
    }),
    dim = c(dim(x)[1], dim(x)[2], dim(x)[3]))

    y <- array(apply(x_sdev, c(3), function(i) {
      1 / (sdev * sqrt(2 * pi)) * i
    }),
    dim = c(dim(x)[1], dim(x)[2], dim(x)[3])) # 1st half of y command

    y <- aperm(apply(y, c(1,2), function(i) {
      i[(diy + ceiling(0.5 * diy)):(2 * diy + 1)] <- i[(diy - ceiling(0.5 * diy)):diy]
      i
    }), c(2,3,1))

    # y[(diy + ceiling(0.5 * diy)):(2 * diy + 1)] <- y[(diy - ceiling(0.5 * diy)):diy]

    st <- diy + 1 - pkday

    st_apply <- aperm(apply(st, c(1,2), function(i) {
      i:(i + diy - 1)
    }), c(2, 3, 1))

    # temporarily convert y and st_apply to 2D
    y2d <- array(y, dim = c(dim(y)[1] * dim(y)[2], dim(y)[3]))
    st_apply2d <- array(st_apply,
                        dim = c(dim(st_apply)[1] * dim(st_apply)[2], dim(st_apply)[3]))
    count <- 0
    y2d <- aperm(apply(y2d, c(1), function(i) {
      count <<- count + 1
      i[st_apply2d[count,]]
    }), c(2, 1))
    y <- array(y2d, dim = c(dim(y)[1], dim(y)[2], dim(y2d)[2]))

    # y <- y[st:(st + diy - 1)]
    x <- c(1:diy)
    x <- c(0, x, c(366:375))

    y <- aperm(apply(y, c(1,2), function(i) {
      c(i[diy], i, i[1:10])
    }), c(2,3,1))
    # y <- c(y[diy], y, y[1:10])

    sel <-c(0:15) * 25 + 1
    x<-x[sel]

    y <- aperm(apply(y, c(1,2), function(i) {
      i[sel]
    }), c(2,3,1))
    # y<-y[sel]
    tme <- as.POSIXct((x * 24 * 3600), origin = paste0(yr - 1,"-12-31 12:00"), tz = "GMT")

    xy <- apply(y, c(1,2), function(i) {
      spline(tme, i, n = diy * 24 + 241)
    })
    # xy <- spline(tme, y, n = diy * 24 + 241)
    tme2 <- as.POSIXlt(xy[[1]]$x, origin = "1970-01-01 00:00", tz = "GMT")
    sel <- which(tme2$year + 1900 == yr)
    y <- aperm(apply(xy, c(1,2), function(i) {
      i[[1]]$y[sel]
    }), c(2,3,1))
    # y <- xy$y[sel]
    dify <- apply(y, c(1,2), function(i) {
      max(i) - min(i)
    })
    ynew <- apply(y, c(3), function(i) {
      i * (difv / dify)
    })
    y <- array(ynew, dim = c(dim(difv)[1], dim(difv)[2], dim(y)[3]))
    # y <- y * (difv / dify)

    count <- 0
    y <- aperm(apply(y, c(1,2), function(i) {
      count <<- count + 1
      i + minlai[count] - min(i)
    }), c(2,3,1))
    # y <- y + minlai - min(y)
    return(y)
  }
  long <- ifelse(long > 180.9375, long - 360, long)
  long <- ifelse(long < -179.0625, long + 360, long)
  lat <- ifelse(lat< -89.49406, -89.49406, lat)
  lat <- ifelse(lat> 89.49406, 89.49406, lat)
  ll <- SpatialPoints(data.frame(x = long, y = lat))
  diy <- 366
  if (year%%4 == 0) diy <- 366
  if (year%%100 == 0 & year%%400 != 0) diy <- 365
  mmonth <-c(16, 45.5, 75, 105.5, 136, 166.5, 197, 228, 258.5, 289, 319.5, 350)
  if (diy == 365) mmonth[2:12] <- mmonth[2:12] + 0.5
  e <- extent(c(-179.0625, 180.9375, -89.49406, 89.49406))
  clim <- c(meantemp, cvtemp, rainfall, cvrain, wetmonth)
  for (i in 1:5) {
    if (is.na(clim[i])) {
      r <- raster(microctools::globalclimate[,,i])
      extent(r) <- e
      clim[i] <- raster::extract(r, ll)
    }
  }
  wgts <- function(x1, x2, ll, lmn, lmx) {
    ll <- ifelse(ll < lmn, lmn, lat)
    ll <- ifelse(ll > lmx, lmx, lat)
    w <- 1 - (abs(ll - lmn)  / (abs(ll - lmn)  + abs(ll - lmx)))
    y <- w * x1 + (1 - w) * x2
    y
  }

  habitat_array <- raster::as.array(habitat_r)

  # By habitat type
  habfun <- function(habitat, returnvar) {
    if (habitat == "Evergreen needleleaf forest" | habitat == 1) {
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
      x <- 0.4
      hgt <- 15
    }
    if (habitat == "Evergreen broadleaf forest" | habitat == 2) {
      hperiod <-  154.505 + 2.040 * clim[1]
      hperiod <- ifelse(hperiod < 50, 50, hperiod)
      peakdoy <- peakdoy <- mmonth[round(clim[5], 0)]
      maxlai <- 1.83 + 0.22 * log(clim[3])
      minlai <- (-1.09) + 0.4030 * log(clim[3])
      minlai <- ifelse(minlai < 1, 1, minlai)
      x <- 1.2
      hgt <- 20
    }
    if (habitat == "Deciduous needleleaf forest" | habitat == 3) {
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
      x <- 0.4
      hgt <- 10
    }
    if (habitat == "Deciduous broadleaf forest" | habitat == 4) {
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
      x <- 1.2
      hgt <- 15
    }
    if (habitat == "Mixed forest" | habitat == 5) {
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
      x <- 0.8
      hgt <- 10
    }
    if (habitat == "Closed shrublands" | habitat == 6) {
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
      x <- 1
      hgt <- 2
    }
    if (habitat == "Open shrublands" | habitat == 7) {
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
      x <- 0.7
      hgt <- 1.5
    }
    if (habitat == "Woody savannas" | habitat == 8) {
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
      x <- 0.7
      hgt <- 3
    }
    if (habitat == "Savannas" | habitat == 9 |
        habitat == "Short grasslands" | habitat == 10 |
        habitat == "Tall grasslands" | habitat == 11) {
      h2 <- 71.72 + 3.012 * clim[1]
      h1 <- 269.22 -  89.79 * clim[4]
      hperiod <- wgts(h1, h2, abs(lat), 0, 20)
      hperioid <- ifelse(hperiod < 31.9, 31.9, hperiod)
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
      x <- 0.15
      if (habitat == "Savannas" | habitat == 9) hgt <- 1.5
      if (habitat == "Short grasslands" | habitat == 10) hgt <- 0.25
      if (habitat == "Tall grasslands" | habitat == 11) hgt <- 1.5
    }
    if (habitat == "Permanent wetlands" | habitat == 12) {
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
      x <- 1.4
      hgt <- 0.5
    }
    if (habitat == "Croplands" | habitat == 13) {
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
      x <- 0.2
      hgt <- 0.5
    }
    if (habitat == "Urban and built-up" | habitat == 14) {
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
      x <- 1
      hgt <- 1.5
    }
    if (habitat == "Cropland/Natural vegetation mosaic" | habitat == 15) {
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
      x <- 0.5
      hgt <- 1
    }
    if (habitat == "Barren or sparsely vegetated" | habitat == 16) {
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
      x <- 0.6
      hgt <- 0.15
    }
    if (habitat == "Open water" | habitat == 17) {
      hperiod <- 100
      peakdoy <- 50
      maxlai <- 0
      minlai <- 0
      hgt <- 0
      x <- 0
    }
    return(eval(parse(text = returnvar)))
  }

  minlai_array <- apply(habitat_array, c(1,2), habfun, "minlai")
  maxlai_array <- apply(habitat_array, c(1,2), habfun, "maxlai")
  x_array <- apply(habitat_array, c(1,2), habfun, "x")
  hgt_array <- apply(habitat_array, c(1,2), habfun, "hgt")
  hperiod_array <- apply(habitat_array, c(1,2), habfun, "hperiod")
  peakdoy_array <- apply(habitat_array, c(1,2), habfun, "peakdoy")

  lai <- laigaus(minlai_array, maxlai_array, peakdoy_array, hperiod_array, year)

  count <- 0
  lai <- aperm(apply(lai, c(1,2), function(i) {
    count <<- count + 1
    if (habitat_array[count] ==  "Short grasslands" | habitat_array[count] == 10) {
      i <- i / 2
    }
    i
  }), c(2,3,1))

  # Time series will be identical for every pixel, so extract length from just
  # first pixel: lai[1,1,]
  tme <- c(1:length(lai[1,1,])) - 1
  tme <- as.POSIXlt(tme * 3600, origin = paste0(year,"-01-01 00:00"), tz = "UTC")
  return(list(lai = lai, x = x_array, height = hgt_array, obs_time = tme))
}
#' Returns PAI for multiple years
.PAI.sort_array <- function(habitat_array, lat, long, tme) {
  tme2<-.tme.sort(tme)
  lai<-array(0, dim = c(50,50,1)); ota<-as.POSIXlt(0,origin="1970-01-01 00:00",tz="UTC")
  yrs<-unique(tme2$year) +1900
  for (i in 1:length(yrs)) {
    pai <- paifromhab(habitat_array, lat, long, yrs[i])
    sel<-which(tme2$year+1900 == yrs[i])
    tme1<-tme2[sel]
    tmemn<-as.numeric(min(tme1))
    tmemx<-as.numeric(max(tme1))
    ot<-as.numeric(pai$obs_time)
    sel <- which(ot >= tmemn & ot <= tmemx)
    l <- aperm(apply(pai$lai, c(1,2), function(i) {
      i[sel]}), c(2,3,1))
    otx<- pai$obs_time[sel]
    lai <- abind::abind(lai, l, along = 3)
    ota<-c(ota,otx)
  }
  lai<-lai[,,-1]; ota<-ota[-1]
  ota<-as.POSIXlt(as.numeric(ota),origin="1970-01-01 00:00",tz="UTC")
  # spline to correct time interval
  if (length(tme) > 1) {
    int <- as.numeric(tme[2])-as.numeric(tme[1])
  } else int <- 3600
  n <- (length(ota)-1)*(3600/(int))+1
  n<-round(n,0)
  ota2<-stats::spline(ota~ota,n=n)$y
  lai2 <- aperm(apply(lai, c(1,2), function(i) {
    spline(i~ota, n = n)$y
  }), c(2,3,1))

  # lai2<-stats::spline(lai~ota,n=n)$y
  tz <- attr(tme,"tzone")[1]
  ota2<-as.POSIXlt(ota2,origin="1970-01-01 00:00",tz=tz)
  sel <- which(ota2 >= tme[1])[1:length(tme)]
  pai$lai <- aperm(apply(lai2, c(1,2), function(i) {
    i[sel]}), c(2,3,1))

  # pai$lai<-lai2[sel]
  pai$obs_time<-tme
  pai
}

vegpfromhab <- function(habitat_r, lat, long, tme, m = 1) {
  # Check class of tme to prevent future failure
  if (any(class(tme) != "POSIXlt")) {
    try_convert <- try(tme<-as.POSIXlt(tme))
    if (any(class(try_convert) == "try-error")) {
      stop("tme must be provided as a single value or vector of POSIXlt objects.")
    }
  }
  year <- unique(as.double(format(as.Date(tme, format="%d/%m/%Y"),"%Y")))

  habitat_array <- raster::as.array(habitat_r)

  if (length(year) == 1) {
    pai <- paifromhab(habitat_r, lat, long, year)
  }
  if (length(year) > 1) {
    pai <- .PAI.sort_array(habitat_array, lat, long, tme)
  }

  # By habitat type
  habfun2 <- function(habitat, returnvar) {
    if (habitat == "Evergreen needleleaf forest" | habitat == 1) {
      # m2 <- round((1/15)*m,0)
      # m2<-ifelse(m2<1,1,m2)
      # pl<-.habgen(1,lat,long,tme,PAIt,m,m2,7.5,70,0.75)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.4,0.7,0.1)
      # iw <- iwgeometry(m, iwmin = 0.8, iwmax = 1.2)
      refls = 0.25 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.3 #  Reflectivity of leaves to PAR
      lw <- 0.01
      gsmax <- 0.33
      leafd <- 0.04 # leaf diameter
      phw <- 500
      uhgt <- 1
      # wgt <- (m2 / m) * 0.25
      # if (any(class(PAI) == "matrix")) {
      #   sPAI<-apply(PAI,2,sum)*wgt
      # } else sPAI <- sum(PAI)*wgt
      # zm0 <- roughlength(uhgt, PAI = sPAI)
    }
    if (habitat == "Evergreen broadleaf forest" | habitat == 2) {
      # m2 <- round((4/20)*m,0)
      # m2<-ifelse(m2<1,1,m2)
      # pl<-.habgen(2,lat,long,tme,PAIt,m,m2,6.5,70,0.8)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.55, 1,0.2)
      # iw <- iwgeometry(m, iwmin = 0.8, iwmax = 1.2)
      refls = 0.3 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.35
      gsmax <- 0.33
      leafd <- 0.055 # leaf diameter
      phw <- 1100
      uhgt <- 4
      # wgt <- (m2 / m) * 0.25
      # if (any(class(PAI) == "matrix")) {
      #   sPAI<-apply(PAI,2,sum)*wgt
      # } else sPAI <- sum(PAI)*wgt
      # zm0 <- roughlength(uhgt, PAI = sPAI)
    }
    if (habitat == "Deciduous needleleaf forest" | habitat == 3) {
      # m2 <- round((1/10)*m,0)
      # m2<-ifelse(m2<1,1,m2)
      # pl<-.habgen(3,lat,long,tme,PAIt,m,m2,7.5,70,0.75)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.4,0.7,0.1)
      # iw <- iwgeometry(m, iwmin = 0.8, iwmax = 1.2)
      refls = 0.3 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.25 #  Reflectivity of leaves to PAR
      lw <- 0.01
      gsmax <- 0.33
      leafd <- 0.04 # leaf diameter
      phw <- 500
      uhgt <- 1
      # wgt <- (m2 / m) * 0.25
      # if (any(class(PAI) == "matrix")) {
      #   sPAI<-apply(PAI,2,sum)*wgt
      # } else sPAI <- sum(PAI)*wgt
      # zm0 <- roughlength(uhgt, PAI = sPAI)
    }
    if (habitat == "Deciduous broadleaf forest" | habitat == 4) {
      # m2 <- round((2/15)*m,0)
      # m2<-ifelse(m2<1,1,m2)
      # pl<-.habgen(4,lat,long,tme,PAIt,m,m2,6.5,70,0.775)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.5, 0.6,0.15)
      # iw <- iwgeometry(m, iwmin = 0.8, iwmax = 1.2)
      refls = 0.3 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.07
      gsmax <- 0.23
      leafd <- 0.05 # leaf diameter
      phw <- 700
      uhgt <- 2
      # wgt <- (m2 / m) * 0.25
      # if (any(class(PAI) == "matrix")) {
      #   sPAI<-apply(PAI,2,sum)*wgt
      # } else sPAI <- sum(PAI)*wgt
      # zm0 <- roughlength(uhgt, PAI = sPAI)
    }
    if (habitat == "Mixed forest" | habitat == 5) {
      # m2 <- round((1.5/10)*m,0)
      # m2<-ifelse(m2<1,1,m2)
      # pl<-.habgen(5,lat,long,tme,PAIt,m,m2,7,70,0.775)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.45,0.65,0.12)
      # iw <- iwgeometry(m, iwmin = 0.8, iwmax = 1.2)
      refls = 0.28 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.25 #  Reflectivity of leaves to PAR
      lw <- 0.04
      gsmax <- 0.28
      leafd <- 0.045 # leaf diameter
      phw <- 600
      uhgt <- 1.5
      # wgt <- (m2 / m) * 0.25
      # if (any(class(PAI) == "matrix")) {
      #   sPAI<-apply(PAI,2,sum)*wgt
      # } else sPAI <- sum(PAI)*wgt
      # zm0 <- roughlength(uhgt, PAI = sPAI)
    }
    if (habitat == "Closed shrublands" | habitat == 6) {
      # pl<-.habgen(6,lat,long,tme,PAIt,m,0,6,80,0.85,under=F)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.2,0.5,0.05)
      # iw <- iwgeometry(m, iwmin = 0.4, iwmax = 0.9)
      refls = 0.3 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.03
      gsmax <- 0.35
      leafd <- 0.02 # leaf diameter
      phw <- 500
      uhgt <- 0.05
      # zm0 <- 0.004
    }
    if (habitat == "Open shrublands" | habitat == 7) {
      # pl<-.habgen(7,lat,long,tme,PAIt,m,0,6,80,0.75,under=F)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.2,0.5,0.05)
      # iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
      refls = 0.3 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.03
      gsmax <- 0.35
      leafd <- 0.02 # leaf diameter
      phw <- 500
      uhgt <- 0.05
      # zm0 <- 0.004
    }
    if (habitat == "Woody savannas" | habitat == 8) {
      # m2<-round((0.75/3)*m,0)
      # m2<-ifelse(m2<1,1,m2)
      # pl<-.habgen(8,lat,long,tme,PAIt,m,m2,6.5,70,0.75)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.2,0.3,0.1)
      # iw <- iwgeometry(m, iwmin = 0.5, iwmax = 0.98)
      refls = 0.35 # reflectivity (shortwave radiation) of leaves
      refw = 0.2 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.01
      gsmax <- 0.33
      leafd <- 0.02 # leaf diameter
      phw <- 300
      uhgt <- 0.75
      # wgt <- (m2 / m) * 0.25
      # if (any(class(PAI) == "matrix")) {
      #   sPAI<-apply(PAI,2,sum)*wgt
      # } else sPAI <- sum(PAI)*wgt
      # zm0 <- roughlength(uhgt, PAI = sPAI)
    }
    if (habitat == "Savannas" | habitat == 9) {
      # pl<-.habgen(9,lat,long,tme,PAIt,m,0,1,50,0.7,under=F)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.2,0.3,0.02)
      # iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
      refls = 0.35 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.01
      gsmax <- 0.33
      leafd <- 0.02 # leaf diameter
      phw <- 300
      uhgt <- 0.05
      # zm0 <- 0.004
    }
    if (habitat == "Short grasslands" | habitat == 10) {
      # pl<-.habgen(10,lat,long,tme,PAIt,m,0,1,50,0.85,under=F)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.2,0.3,0.02)
      # iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
      refls = 0.35 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.01
      gsmax <- 0.33
      leafd <- 0.02 # leaf diameter
      phw <- 300
      uhgt <- 0.05
      # zm0 <- 0.004
    }
    if (habitat == "Tall grasslands" | habitat == 11) {
      # pl<-.habgen(11,lat,long,tme,PAIt,m,0,1,50,0.85,under=F)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.2,0.3,0.02)
      # iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
      refls = 0.35 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.01
      gsmax <- 0.33
      leafd <- 0.02 # leaf diameter
      phw <- 300
      uhgt <- 0.05
      # zm0 <- 0.004
    }
    if (habitat == "Permanent wetlands" | habitat == 12) {
      # pl<-.habgen(12,lat,long,tme,PAIt,m,0,1,50,0.95,under=F)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.2,0.3,0.02)
      # iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
      refls = 0.5 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.09
      gsmax <- 0.55
      leafd <- 0.02 # leaf diameter
      phw <- 600
      uhgt <- 0.02
      # zm0 <- 0.002
    }
    if (habitat == "Croplands" | habitat == 13) {
      # pl<-.habgen(13,lat,long,tme,PAIt,m,0,1,50,0.8,under=F)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.25,0.3,0.02)
      # iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
      refls = 0.3 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.02
      gsmax <- 0.33
      leafd <- 0.025 # leaf diameter
      phw <- 300
      uhgt <- 0.05
      # zm0 <- 0.004
    }
    if (habitat == "Urban and built-up" | habitat == 14) {
      # m2 <- round((0.5/1.5) * m, 0)
      # m2<-ifelse(m2<1,1,m2)
      # pl<-.habgen(14,lat,long,tme,PAIt,m,m2,6.5,70,0.75)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.5, 0.6,0.15)
      # iw <- iwgeometry(m, iwmin = 0.6, iwmax = 1.4)
      refls = 0.3 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.04
      gsmax <- 0.33
      leafd <- 0.05 # leaf diameter
      phw <- 600
      uhgt <- 0.1
      # wgt <- (m2 / m) * 0.25
      # if (any(class(PAI) == "matrix")) {
      #   sPAI<-apply(PAI,2,sum)*wgt
      # } else sPAI <- sum(PAI)*wgt
      # zm0 <- roughlength(uhgt, PAI = sPAI)
    }
    if (habitat == "Cropland/Natural vegetation mosaic" | habitat == 15) {
      # pl<-.habgen(15,lat,long,tme,PAIt,m,0,1,50,0.75,under=F)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.35,0.5,0.02)
      # iw <- iwgeometry(m, iwmin = 0.5, iwmax = 0.98)
      refls = 0.3 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.03
      gsmax <- 0.3
      leafd <- 0.035 # leaf diameter
      phw <- 400
      uhgt <- 0.05
      # zm0 <- 0.004
    }
    if (habitat == "Barren or sparsely vegetated" | habitat == 16) {
      # pl<-.habgen(16,lat,long,tme,PAIt,m,0,1,50,0.55,under=F)
      # PAI<-pl$PAI
      # pLAI<-pl$pLAI
      # thickw <- thickgeometry(m, 0.25,0.3,0.02)
      # iw <- iwgeometry(m, iwmin = 0.36, iwmax = 0.9)
      refls = 0.3 # reflectivity (shortwave radiation) of leaves
      refw = 0.1 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0.2 #  Reflectivity of leaves to PAR
      lw <- 0.015
      gsmax <- 0.33
      leafd <- 0.025 # leaf diameter
      phw <- 600
      uhgt <- 0.01
      # zm0 <- 0.001
    }
    if (habitat == "Open water" | habitat == 17) {
      # PAI<-0
      # pLAI<-0
      # thickw <- 0
      # iw <- 0
      refls = 0 # reflectivity (shortwave radiation) of leaves
      refw = 0 # reflectivity (shortwave radiation) of woody vegetation
      reflp = 0 #  Reflectivity of leaves to PAR
      lw <- 0
      gsmax <- 0
      leafd <- 0 # leaf diameter
      phw <- 0
      uhgt <- 0
      # zm0 <- 0
    }

    return(eval(parse(text = returnvar)))

  }

  # PAI_array <- apply(habitat_array, c(1,2), habfun2, "PAI")
  lw_array <- apply(habitat_array, c(1,2), habfun2, "lw")
  # iw_array <- apply(habitat_array, c(1,2), habfun2, "iw")
  uhgt_array <- apply(habitat_array, c(1,2), habfun2, "uhgt")
  # zm0_array <- apply(habitat_array, c(1,2), habfun2, "zm0")
  # pLAI_array <- apply(habitat_array, c(1,2), habfun2, "pLAI")
  refls_array <- apply(habitat_array, c(1,2), habfun2, "refls")
  refw_array <- apply(habitat_array, c(1,2), habfun2, "refw")
  reflp_array <- apply(habitat_array, c(1,2), habfun2, "reflp")
  gsmax_array <- apply(habitat_array, c(1,2), habfun2, "gsmax")
  leafd_array <- apply(habitat_array, c(1,2), habfun2, "leafd")
  # thickw_array <- apply(habitat_array, c(1,2), habfun2, "thickw")
  phw_array <- apply(habitat_array, c(1,2), habfun2, "phw")

  if (any(is.na(pai))) {
    warning("NA values returned in PAI. Possible reasons:\n  - tme time series covers multiple years (only one year allowed)\n  - tme time series is sub-daily but not at hourly resoultion")
  }

  # Coerce to raster
  hgt <- habitat_r
  values(hgt) <- pai$height
  x <- habitat_r
  values(x) <- pai$x
  gsmax <- habitat_r
  values(gsmax) <- gsmax_array
  leafr <- habitat_r
  values(leafr) <- refls_array
  clump <- habitat_r
  values(clump) <- 0
  leafd <- habitat_r
  values(leafd) <- leafd_array

  vegp <- list(pai = pai$lai, hgt = hgt, x = x, gsmax = gsmax,
               leafr = leafr, clump=clump, leafd = leafd)
  class(vegp) <- "vegparams"
  return(vegp)
}
