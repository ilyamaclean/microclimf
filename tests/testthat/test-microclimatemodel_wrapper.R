test_that("The full microclimate model runs on a 24h test case and returns expected outputs", {
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
  # vegp
  vegp <- c(
    h=0.5, pai=2, vegx=1.0, clump=0.1, lref=0.4,
    ltra=0.2, leafd=0.05, em = 0.97, gsmax=0.13
  )
  # groundp
  groundp <- c(
    gref=0.15, slope=0, aspect=180, groundem = 0.97, rho=1.53, Vm=0.509, Vq=0.06, Mc=0.5422,
    soilb = 5.2, psie = 2.6, Smin =0.091, Smax=0.419
  )
  # Run below canopy microclimate model
  BLout = BigLeafCpp(obstime, climdata, vegp, groundp, rep(0.3, 24),
    lat = 50, lon = -5, 25.0, 2, 100, 0.5, 0.5, 0.1, FALSE)
  out = microclimatemodel_wrapper(obstime, climdata, BLout, vegp, groundp,
             reqhgt = 0.05, zref = 2, lat = 50, lon = -5)
  # Pull known fields the function returns
  getv <- function(nm) if (!is.null(out[[nm]])) out[[nm]] else NULL
  Tz    <- getv("Tz")
  tleaf    <- getv("tleaf")
  rh     <- getv("rh")
  uz     <- getv("uz")
  Rdirdown <- getv("Rdirdown")
  Rdifdown  <- getv("Rdifdown")
  Rswup  <- getv("Rswup")
  Rlwdn  <- getv("Rlwdown")
  Rlwup    <- getv("Rlwup")
  # ---- Finiteness & lengths ----
  nums <- unlist(Filter(is.numeric, list(Tz,tleaf,rh,uz,Rdirdown,Rdifdown,Rswup,Rlwdn,Rlwup)))
  expect_true(length(nums) > 0)
  expect_true(all(is.finite(nums)))
  # Lengths should match the 24 h input where applicable
  check_len <- function(x) if (!is.null(x)) expect_equal(length(x), n)
  lapply(list(Tz,tleaf,rh,uz,Rdirdown,Rdifdown,Rswup,Rlwdn,Rlwup), check_len)
  # Calculate variables for checking
  Tair_dif <- abs(Tz - Tair)
  Tleaf_dif <- abs(tleaf - Tair)
  uz_rat <- uz / U2
  Rdir_rat <- (Rdirdown * si) / Rdir; Rdir_rat[Rdir == 0] <- 0.2
  Rdif_rat <- Rdifdown / Rdif; Rdif_rat[Rdif == 0] <- 0.3
  Rsw_rat <- Rswup / SWd; Rsw_rat[SWd == 0] <- 0.112
  Rlw_rat1 <- Rlwdn / LWd
  Rlw_rat2 <- Rlwup / LWd
  # ---- Bounds / sanity checks ----
  if (!is.null(Tair_dif)) { expect_lte(max(Tair_dif), 5) }
  if (!is.null(Tleaf_dif)) { expect_lte(max(Tleaf_dif), 2) }
  if (!is.null(rh)) { expect_gte(min(rh), min(RH) - 5); expect_lte(max(rh), 100) }
  if (!is.null(uz_rat)) { expect_gte(min(uz_rat), 0.08); expect_lte(max(uz_rat), 0.1) }
  if (!is.null(Rdir_rat)) { expect_gte(min(Rdir_rat), 0); expect_lte(max(Rdir_rat), 0.25) }
  if (!is.null(Rdif_rat)) { expect_gte(min(Rdif_rat), 0.27); expect_lte(max(Rdif_rat), 0.32) }
  if (!is.null(Rsw_rat)) { expect_gte(min(Rsw_rat), 0.03); expect_lte(max(Rsw_rat), 0.15) }
  if (!is.null(Rlw_rat1)) { expect_gte(min(Rlw_rat1), 0.94); expect_lte(max(Rlw_rat1), 1.2) }
  if (!is.null(Rlw_rat2)) { expect_gte(min(Rlw_rat2), 0.94); expect_lte(max(Rlw_rat2), 1.2) }
  # Run microclimate model at ground level
  out <- microclimatemodel_wrapper(obstime, climdata, BLout, vegp, groundp,
             reqhgt = 0.0, zref = 2, lat = 50, lon = -5)
  getv <- function(nm) if (!is.null(out[[nm]])) out[[nm]] else NULL
  Tz    <- getv("Tz")
  soilm    <- getv("soilm")
  Rdirdown <- getv("Rdirdown")
  Rdifdown  <- getv("Rdifdown")
  Rswup  <- getv("Rswup")
  Rlwdn  <- getv("Rlwdown")
  Rlwup    <- getv("Rlwup")
  # ---- Finiteness & lengths ----
  nums <- unlist(Filter(is.numeric, list(Tz,soilm,Rdirdown,Rdifdown,Rswup,Rlwdn,Rlwup)))
  expect_true(length(nums) > 0)
  expect_true(all(is.finite(nums)))
  # Lengths should match the 24 h input where applicable
  check_len <- function(x) if (!is.null(x)) expect_equal(length(x), n)
  lapply(list(Tz,soilm,Rdirdown,Rdifdown,Rswup,Rlwdn,Rlwup), check_len)
  # Calculate variables for checking
  Tz_dif <- abs(Tz - Tair)
  Rdir_rat <- (Rdirdown * si) / Rdir; Rdir_rat[Rdir == 0] <- 0.2
  Rdif_rat <- Rdifdown / Rdif; Rdif_rat[Rdif == 0] <- 0.238
  Rsw_rat <- Rswup / SWd; Rsw_rat[SWd == 0] <- 0.148
  Rlw_rat1 <- Rlwdn / LWd
  Rlw_rat2 <- Rlwup / LWd
  # ---- Bounds / sanity checks ----
  if (!is.null(Tz_dif)) { expect_lte(max(Tz_dif), 5) }
  if (!is.null(soilm)) { expect_gte(min(soilm), 0.299); expect_lte(max(soilm), 0.301) }
  if (!is.null(Rdir_rat)) { expect_gte(min(Rdir_rat), 0); expect_lte(max(Rdir_rat), 0.25) }
  if (!is.null(Rdif_rat)) { expect_gte(min(Rdif_rat), 0.2); expect_lte(max(Rdif_rat), 0.3) }
  if (!is.null(Rsw_rat)) { expect_gte(min(Rsw_rat), 0.03); expect_lte(max(Rsw_rat), 0.15) }
  if (!is.null(Rlw_rat1)) { expect_gte(min(Rlw_rat1), 0.9); expect_lte(max(Rlw_rat1), 1.1) }
  if (!is.null(Rlw_rat2)) { expect_gte(min(Rlw_rat2), 0.9); expect_lte(max(Rlw_rat2), 1.15) }
  # Run microclimate model below ground
  out <- microclimatemodel_wrapper(obstime, climdata, BLout, vegp, groundp,
              reqhgt = -0.05, zref = 2, lat = 50, lon = -5)
  getv <- function(nm) if (!is.null(out[[nm]])) out[[nm]] else NULL
  Tz    <- getv("Tz")
  soilm    <- getv("soilm")
  # ---- Finiteness & lengths ----
  nums <- unlist(Filter(is.numeric, list(Tz,soilm)))
  expect_true(length(nums) > 0)
  expect_true(all(is.finite(nums)))
  # Lengths should match the 24 h input where applicable
  check_len <- function(x) if (!is.null(x)) expect_equal(length(x), n)
  lapply(list(Tz,soilm), check_len)
  if (!is.null(Tz)) { expect_gte(min(Tz), 8); expect_lte(max(Tz), 9) }
  if (!is.null(soilm)) { expect_gte(min(soilm), 0.299); expect_lte(max(soilm), 0.301) }
})
