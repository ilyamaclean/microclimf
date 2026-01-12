test_that("Test snow model", {
  hrs   <- 0:23
  n     <- length(hrs)
  obstime <- data.frame(
    year  = rep(2024L, n),
    month = rep(3L,    n),
    day   = rep(21L,   n),
    hour  = as.numeric(hrs)
  )
  # Temperatures in deg C (sinusoid), humidity %, pressure kPa
  Tair   <- -5 + 5 * sin((hrs - 8) / 24 * 2 * pi)          # °C
  ea <- 0.7 * satvapCpp(mean(Tair))
  RH <- 0; for (hr in 1:24) RH[hr] <- (ea / satvapCpp(Tair[hr])) * 100
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
  Prec   <- rep(1, n)                                      # mm h-1
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
  # vegp
  vegp <- c(pai = 2, hgt = 0.5, ltra = 0.05, clump = 0)
  # other
  other <- c(slope = 0, aspect = 180, lat = 50, lon = -5, 
             zref = 2, isnowd = 0, isnowa = 0)
  # Run snow model
  pmod<-pointmodelsnow(obstime, climdata, vegp, other, "Taiga")
  # Pull known fields the function returns
  getv <- function(nm) if (!is.null(pmod[[nm]])) pmod[[nm]] else NULL
  Tc    <- getv("Tc")
  Tg    <- getv("Tg")
  sdepc <- getv("sdepc")
  sdepg <- getv("sdepg")
  sdenc <- getv("sdenc")
  sdeng <- getv("sdeng")
  # ---- Finiteness & lengths ----
  nums <- unlist(Filter(is.numeric, list(Tc,Tg,sdepc,sdepg,sdenc,sdeng)))
  expect_true(length(nums) > 0)
  expect_true(all(is.finite(nums)))
  # Lengths should match the 24 h input where applicable
  check_len <- function(x) if (!is.null(x)) expect_equal(length(x), n)
  check_len2 <- function(x) if (!is.null(x)) expect_equal(length(x), n + 1)
  lapply(list(Tc,Tg,sdenc,sdeng), check_len)
  lapply(list(sdepc,sdepg), check_len2)
  # Calculate variables for checking
  Tcdif <- abs(Tc - Tair)
  Tgdif <- abs(Tg - Tair)
  depdif <- sdepc - sdepg
  snowacc <- (sdepc[n] - sdepc[1]) * sdenc[n] # snow water equivelent in mm
  precsum <- sum(Prec)
    # ---- Bounds / sanity checks ----
  if (!is.null(Tcdif)) { expect_lte(max(Tcdif), 1.5) }
  if (!is.null(Tgdif)) { expect_lte(max(Tgdif), 2.5) }
  if (!is.null(depdif)) { expect_gte(min(depdif), 0); expect_lte(max(depdif), 0.05) }
  if (!is.null(snowacc)) { expect_gte(snowacc, 0.8 * precsum); expect_lte(snowacc, precsum) }
  if (!is.null(sdenc)) { expect_gte(min(sdenc), 216); expect_lte(max(sdenc), 218) }
  if (!is.null(sdeng)) { expect_gte(min(sdeng), 216); expect_lte(max(sdeng), 218) }
})
