test_that("weatherhgtCpp runs on a 24 h test case and returns expected outputs", {
  # ---- Build minimal 24h inputs ----
  hrs   <- 0:23
  n     <- length(hrs)
  obstime <- data.frame(
    year  = rep(2024L, n),
    month = rep(3L,    n),
    day   = rep(21L,   n),
    hour  = as.numeric(hrs)
  )
  # Temperatures in deg C (sinusoid), humidity %, pressure kPa
  Tair   <- 10 + 5 * sin((hrs - 8) / 24 * 2 * pi)          # °C
  ea <- 0.7 * satvapCpp(mean(Tair))
  RH <- 0; for (hr in 1:24) RH[hr] <- (ea / satvapCpp(Tair[hr])) * 100
  Pk     <- rep(101.3, n)                                  # kPa  (your code uses kPa)
  SWd    <- pmax(0, 600 * sin((hrs - 6) / 12 * pi))        # W m-2
  Rdif   <- pmin(SWd, 0.3 * SWd)                           # W m-2 (simple split)
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

  # ---- Call the weather height adjustment function ----
  out <- weatherhgtCpp(
    obstime, climdata, zin = 2, uzin = 2, zout = 10, lat = 50, lon = -5
  )

  # ---- Basic structure ----
  expect_true(is.data.frame(out))

  # Pull known fields the function returns
  getv <- function(nm) if (!is.null(out[[nm]])) out[[nm]] else NULL

  temp    <- getv("temp")
  relhum    <- getv("relhum")
  pres     <- getv("pres")
  swdown     <- getv("swdown")
  difrad <- getv("difrad")
  lwdown  <- getv("lwdown")
  windspeed  <- getv("windspeed")
  winddir    <- getv("winddir")
  precip    <- getv("precip")

  # Calculate difference form inputs
  d_temp <- abs(Tair - temp)
  d_pres <- abs(pres - Pk)
  d_swdown <- abs(swdown - SWd)
  d_difrad <- abs(difrad - Rdif)
  d_lwdown <- abs(lwdown - LWd)
  d_winddir <- abs(winddir - Wdir)
  d_precip <- abs(precip - Prec)

  mu_windspeed <- windspeed / U2
  ea02 <- 0; for (hr in 1:24) ea02[hr] <- satvapCpp(Tair[hr]) * RH[hr] / 100 # vapour pressure (2m)
  ea10 <- 0; for (hr in 1:24) ea10[hr] <- satvapCpp(temp[hr]) * relhum[hr] / 100 # vapour pressure (2m)
  d_ea <- abs(ea - ea10)

  # perform difference tests
  if (!is.null(d_temp)) { expect_lte(max(d_temp), 4) }
  if (!is.null(d_pres)) { expect_lte(max(d_pres), 1) }
  if (!is.null(d_swdown)) { expect_lte(max(d_swdown), 1) }
  if (!is.null(d_difrad)) { expect_lte(max(d_difrad), 1) }
  if (!is.null(d_lwdown)) { expect_lte(max(d_lwdown), 1) }
  if (!is.null(d_winddir)) { expect_lte(max(d_winddir), 1) }
  if (!is.null(d_precip)) { expect_lte(max(d_precip), 1) }
  if (!is.null(mu_windspeed)) { expect_gte(min(mu_windspeed), 1.2); expect_lte(max(mu_windspeed), 1.4) }
  if (!is.null(d_ea)) { expect_lte(max(d_ea), 0.5) }
})
