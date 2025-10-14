test_that("BigLeafCpp runs on a 24 h test case and returns expected outputs", {
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

  # vegp
  vegp <- c(
    h=0.5, pai=2.0, vegx=1.0, clump=0.1,
    lref=0.4, ltra=0.2, leafd=0.05, em=0.97,
    gsmax=0.33, q50=100
  )

  # groundp
  groundp <- c(
    gref=0.15, slope=0, aspect=180, groundem=0.97,
    rho=1.53, Vm=0.509, Vq=0.06, Mc=0.5422,
    b =5.2, Psie=-5.6, Smax=0.42, Smin=0.074
  )

  # soilm: length = n
  soilm <- rep(0.3, n)

  lat <- 50.0; lon <- -5.0

  # ---- Call the model (disable the yearly G tweak for speed/repro) ----
  out <- BigLeafCpp(
    obstime, climdata, vegp, groundp, soilm,
    lat, lon,
    dTmx = 25, zref = 2, maxiter = 50, bwgt = 0.5, tol = 0.5, gmn = 0.1, yearG = FALSE
  )

  # ---- Test basic structure ----
  expect_true(is.list(out))

  # Pull known fields the function returns
  # Tc, Tg (K); H, G, RabsG (W m-2); psih, psim, phih; OL; uf; err; albedo
  getv <- function(nm) if (!is.null(out[[nm]])) out[[nm]] else NULL

  Tc    <- getv("Tc")
  Tg    <- getv("Tg")
  H     <- getv("H")
  G     <- getv("G")
  RabsG <- getv("RabsG")
  psih  <- getv("psih")
  psim  <- getv("psim")
  phih  <- getv("phih")
  OL    <- getv("OL")
  uf    <- getv("uf")
  err   <- getv("err")
  alb   <- getv("albedo")

  # ---- Finiteness & lengths ----
  nums <- unlist(Filter(is.numeric, list(Tc,Tg,H,G,RabsG,psih,psim,phih,OL,uf,err,alb)))
  expect_true(length(nums) > 0)
  expect_true(all(is.finite(nums)))

  # Lengths should match the 24 h input where applicable
  check_len <- function(x) if (!is.null(x)) expect_equal(length(x), n)
  lapply(list(Tc,Tg,H,G,RabsG,psih,psim,phih,OL,uf,alb), check_len)

  # ---- Bounds / sanity checks ----
  if (!is.null(Tc)) { expect_gte(min(Tc), min(Tair) - 5); expect_lte(max(Tc), max(Tair)+10) }
  if (!is.null(Tg)) { expect_gte(min(Tg), min(Tair) - 5); expect_lte(max(Tg), max(Tair)+10) }
  if (!is.null(alb)) { expect_gte(min(alb), 0.01); expect_lte(max(alb), 0.99) }
  if (!is.null(uf))  { expect_gte(min(uf), 2e-4); expect_lte(max(uf), max(U2))}
  if (!is.null(Tc)) { expect_gte(min(Tc), min(Tair) - 5); expect_lte(max(Tc), max(Tair)+10) }
  if (!is.null(Tg)) { expect_gte(min(Tg), min(Tair) - 5); expect_lte(max(Tg), max(Tair)+10) }
  if (!is.null(alb)) { expect_gte(min(alb), 0.01); expect_lte(max(alb), 0.99) }

  # Energy balance limits
  Lwup <- 5.67e-8 * (Tc + 273.15)^4
  Rnet <- (1 - alb) * SWd + 0.97 * (LWd - Lwup)

  if (!is.null(H))  { expect_gte(min(H), min(Rnet + G - 20)); expect_lte(max(H), max(Rnet + G - 20))}
  if (!is.null(G))  { expect_gte(min(G), min(Rnet - 20)); expect_lte(max(G), max(Rnet + 20))}
  if (!is.null(RabsG))  { expect_gte(min(RabsG), 0.5 * min(LWd)); expect_lte(max(RabsG), 1000)}

  # Stability correction limits
  if (!is.null(psih)) { expect_gte(min(psih), -4); expect_lte(max(psih), 3) }
  if (!is.null(psim)) { expect_gte(min(psim), -4); expect_lte(max(psim), 3) }

  # ---- Convergence signal (err small on average) ----
  if (!is.null(err)) {
    expect_lt(err, 0.5)  # uses your default tol scale; adjust if needed
  }
})



