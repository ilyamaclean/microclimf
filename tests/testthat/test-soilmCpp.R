test_that("Test soil water model", {
  hrs   <- c(0:23, 0:23)
  n     <- length(hrs)
  obstime <- data.frame(
    year  = rep(2024L, n),
    month = rep(3L,    n),
    day   = rep(c(21, 22),   each = n/2),
    hour  = as.numeric(hrs)
  )
  # Temperatures in deg C (sinusoid), humidity %, pressure kPa
  Tair   <- 10 + 5 * sin((hrs - 8) / 24 * 2 * pi)          # °C
  ea <- 0.7 * satvapCpp(mean(Tair))
  RH <- 0; for (hr in 1:n) RH[hr] <- (ea / satvapCpp(Tair[hr])) * 100
  Pk     <- rep(101.3, n)                                  # kPa  (your code uses kPa)
  # Simulate radiation
  csr <- clearskyradCpp(obstime$year, obstime$month, obstime$day, obstime$hour,
                        lat = 50, lon = -5, Tair, RH, Pk)
  solvars <- solpositionvCpp(obstime$year, obstime$month, obstime$day, obstime$hour, 
                             lat = 50, lon = -5, slope = 0, aspect = 180)
  si <- solvars$si
  Rdir <- 0.3 * csr * si
  SWd <- 0.5 * csr                                         # W m-2
  Rdif <- SWd - Rdir                                       # W m-2
  LWd    <- rep(350, n)                                    # W m-2
  U2     <- rep(2, n)                                      # m s-1
  Wdir   <- rep(180, n)                                    # deg
  Prec   <- rep(0, n)                                      # mm h-1
  climdata <- data.frame(
    temp      = Tair,
    relhum    = RH,
    pres      = Pk,
    swdown    = SWd,
    difrad    = Rdif,
    lwdown    = LWd,
    windspeed = U2,
    winddir   = Wdir,
    precip    = Prec
  )
  soilm <- soilmCpp(climdata, rmu = 0.021303, mult = 0.000191202, pwr = 1.134773, 
                    Smax = 0.419, Smin = 0.091, Ksat = 5.89, a = 0.059765)
  expect_true(length(soilm) == 2)
  expect_true(all(is.finite(soilm)))
  expect_gte(min(soilm), 0.35); expect_lte(max(soilm), 0.419) 
})
