# The full tier: year-long and grid runs. Skipped unless
# MICROCLIMF_FULL_TESTS=true; see helper-microclimf.R.

test_that("the bundled year converges in every hour and stays physical", {
  skip_unless_full()
  r <- suppressWarnings(runpointmodel(climdata, reqhgt = -0.3, dtm = dtmcaerth, vegp = vegp,
                                      soilc = soilc))
  m <- r$model
  expect_true(all(m$iters < 100))
  for (f in c("Tcanopy", "Tground", "G", "Hout", "theta0", "Tzbelow", "thetazbelow")) {
    expect_true(all(is.finite(m[[f]])), info = f)
  }
  expect_true(all(m$theta0 > 0 & m$theta0 <= max(r$soilc$Smax)))
})

test_that("a part-year run with matemp reproduces the same part of the full year", {
  skip_unless_full()
  mat <- mean(climdata$temp)
  h <- seq_len(180 * 24)
  full <- suppressWarnings(runpointmodel(climdata, reqhgt = -0.3, dtm = dtmcaerth, vegp = vegp,
                                         soilc = soilc, matemp = mat))
  part <- suppressWarnings(runpointmodel(climdata[h, ], reqhgt = -0.3, dtm = dtmcaerth, vegp = vegp,
                                         soilc = soilc, matemp = mat))
  expect_identical(part$weather, full$weather[h, ])
  expect_identical(part$model$Tground, full$model$Tground[h])
})

test_that("a run restarted from soilend continues the unbroken run", {
  skip_unless_full()
  mat <- mean(climdata$temp)
  h1 <- seq_len(180 * 24); h2 <- (180 * 24 + 1):nrow(climdata)
  full <- suppressWarnings(runpointmodel(climdata, reqhgt = -0.3, dtm = dtmcaerth, vegp = vegp,
                                         soilc = soilc, matemp = mat))
  a <- suppressWarnings(runpointmodel(climdata[h1, ], reqhgt = -0.3, dtm = dtmcaerth, vegp = vegp,
                                      soilc = soilc, matemp = mat, soilend = TRUE))
  b <- suppressWarnings(runpointmodel(climdata[h2, ], reqhgt = -0.3, dtm = dtmcaerth, vegp = vegp,
                                      soilc = soilc, matemp = mat, soilinit = a$soilend))
  # Canopy water and any surface pool are not carried across, so allow small
  # differences at the surface; the soil at 30 cm should follow closely.
  first2d <- 1:48
  expect_lt(max(abs(b$model$Tzbelow[first2d] - full$model$Tzbelow[h2][first2d])), 0.02)
  expect_lt(max(abs(b$model$Tground[first2d] - full$model$Tground[h2][first2d])), 0.1)
})

test_that("tiled grid runs equal untiled ones", {
  skip_unless_full()
  pm <- suppressWarnings(runpointmodel(climdata, reqhgt = 0.05, dtm = dtmcaerth, vegp = vegp,
                                       soilc = soilc))
  pms <- subsetpointmodel(pm, days = c(15, 172))
  # A 50 x 50 corner in four tiles: enough for each tile to differ from the
  # whole, which is what the check needs.
  d <- terra::rast(dtmcaerth); e <- terra::ext(d)
  box <- terra::ext(terra::xmin(e), terra::xmin(e) + 50 * terra::res(d)[1],
                    terra::ymax(e) - 50 * terra::res(d)[2], terra::ymax(e))
  crp <- function(x) terra::crop(terra::rast(x), box)
  vg <- lapply(vegp, crp); sl <- lapply(soilc, crp); dm <- crp(dtmcaerth)
  u <- suppressWarnings(rungridmodel(pms, reqhgt = 0.05, dtm = dm, vegp = vg, soilc = sl))
  t <- suppressWarnings(rungridmodel(pms, reqhgt = 0.05, dtm = dm, vegp = vg, soilc = sl, cores = 2, tilesize = 25))
  expect_named(t, names(u))
  for (k in names(u)) {
    # A field with no value at this height is a single NA in both.
    if (!inherits(u[[k]], "SpatRaster")) {
      expect_identical(t[[k]], u[[k]], info = k)
      next
    }
    a <- as.vector(terra::values(u[[k]])); b <- as.vector(terra::values(t[[k]]))
    # Cells off the land are missing in both; the tiled mosaic stores them as
    # NaN rather than NA, which is the same missing value.
    expect_identical(is.na(a), is.na(b), info = k)
    expect_identical(a[!is.na(a)], b[!is.na(b)], info = k)
  }
})

test_that("runbioclim works below ground, on soil temperature whatever temp says", {
  skip_unless_full()
  d <- terra::rast(dtmcaerth); e <- terra::ext(d)
  box <- terra::ext(terra::xmin(e), terra::xmin(e) + 10 * terra::res(d)[1],
                    terra::ymax(e) - 10 * terra::res(d)[2], terra::ymax(e))
  crp <- function(x) terra::crop(terra::rast(x), box)
  vg <- lapply(vegp, crp); sl <- lapply(soilc, crp); dm <- crp(dtmcaerth)
  air <- suppressWarnings(runbioclim(climdata, reqhgt = -0.1, dtm = dm, vegp = vg, soilc = sl,
                                     temp = "air", matemp = 11))
  leaf <- suppressWarnings(runbioclim(climdata, reqhgt = -0.1, dtm = dm, vegp = vg, soilc = sl,
                                      temp = "leaf", matemp = 11))
  a <- terra::values(air)
  land <- is.finite(terra::values(dm))
  expect_true(all(is.finite(a[land, ])))
  expect_identical(terra::values(leaf), a)
})

test_that("a canopy of vanishing plant area stays close to bare ground", {
  skip_unless_full()
  # The two meet only approximately (a standing limitation: the canopy step's
  # approximation against the bare solve); measured at about 0.1 K RMS.
  wx <- climdata[1:(30 * 24), ]
  bare <- core_run(wx, 0, 0)
  thin <- core_run(wx, 0.1, 1e-8)
  expect_lt(sqrt(mean((thin$model$Tcanopy - bare$model$Tcanopy)^2)), 0.3)
})
