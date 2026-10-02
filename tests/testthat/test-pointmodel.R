# Short point-model runs: sanity, and identities established when the relevant
# parts of the model were settled.

wx <- climdata[1:72, ]

test_that("a short run is finite, bounded and converges within the pass cap", {
  r <- core_run(wx, 0.3, 2)
  m <- r$model
  for (f in c("Tcanopy", "Tground", "G", "Hout", "uf", "theta0", "psi_r")) {
    expect_true(all(is.finite(m[[f]])), info = f)
  }
  expect_true(all(m$theta0 > 0 & m$theta0 <= max(r$soilc$Smax)))
  expect_true(all(m$iters < 100))
  expect_true(all(m$Tground > -40 & m$Tground < 70))
})

test_that("zero height and zero plant area are the same bare ground", {
  a <- core_run(wx, 0, 0)
  b <- core_run(wx, 0.5, 0)
  c <- core_run(wx, 0, 2)
  expect_identical(a$model, b$model)
  expect_identical(a$model, c$model)
})

test_that("on bare ground the bulk surface is the soil surface", {
  m <- core_run(wx, 0, 0)$model
  expect_identical(m$Tcanopy, m$Tground)
})

test_that("the wind returned at the weather height is the measured wind", {
  # A short canopy, so the forcing is used at the height it was measured.
  r <- core_run(wx, 0.1, 1, reqhgt = 2)
  expect_equal(r$zref, 2)
  expect_equal(r$model$windzabove, wx$windspeed, tolerance = 1e-9)
})

test_that("the wind returned at a translated reference height is the translated forcing", {
  r <- core_run(wx, 15, 5, reqhgt = NA)
  expect_gt(r$zref, 2)
  r2 <- core_run(wx, 15, 5, reqhgt = r$zref)
  expect_equal(r2$model$windzabove, r2$weather$windspeed, tolerance = 1e-9)
})

test_that("the fields at reqhgt are returned only where that height has them", {
  soil <- c("Tzbelow", "thetazbelow"); air <- c("Tzabove", "windzabove", "RHzabove")
  cols <- function(reqhgt) names(core_run(wx[1:24, ], 0.3, 2, reqhgt = reqhgt)$model)
  below <- cols(-0.1); surface <- cols(0); within <- cols(0.1); top <- cols(0.3); above <- cols(2)
  expect_true(all(soil %in% below) && !any(air %in% below))
  expect_identical(surface, below)
  expect_false(any(c(soil, air) %in% within))
  expect_true(all(air %in% above) && !any(soil %in% above))
  expect_identical(top, above)
  # On bare ground the surface is still soil, and any height above it is above canopy.
  bare <- function(reqhgt) names(core_run(wx[1:24, ], 0, 0, reqhgt = reqhgt)$model)
  expect_identical(bare(0), below)
  expect_identical(bare(0.05), above)
})

test_that("the deepest soil node stays at matemp whatever starting profile is given", {
  prof <- data.frame(depth = c(0, 1, 4), temp = c(15, 15, 22))
  r <- suppressWarnings(runpointmodel(wx, reqhgt = -5, dtm = dtmcaerth, vegp = vegp, soilc = soilc,
                                      matemp = 12, soilinit = prof, runchecks = FALSE))
  expect_equal(r$model$Tzbelow, rep(12, nrow(wx)))
})

test_that("soilend is returned only when asked for, in the form soilinit takes", {
  r0 <- runpointmodel(wx, dtm = dtmcaerth, vegp = vegp, soilc = soilc, matemp = 11, runchecks = FALSE)
  expect_null(r0$soilend)
  r1 <- runpointmodel(wx, dtm = dtmcaerth, vegp = vegp, soilc = soilc, matemp = 11, runchecks = FALSE,
                      soilend = TRUE)
  expect_named(r1$soilend, c("depth", "temp", "theta"))
  expect_equal(r1$soilend$temp[nrow(r1$soilend)], 11)
  # Passing it back raises no warning and starts from it.
  expect_no_warning(runpointmodel(climdata[73:96, ], dtm = dtmcaerth, vegp = vegp, soilc = soilc,
                                  matemp = 11, runchecks = FALSE, soilinit = r1$soilend))
})

test_that("omitting soilinit leaves the run as it was", {
  a <- runpointmodel(wx, dtm = dtmcaerth, vegp = vegp, soilc = soilc, matemp = 11, runchecks = FALSE)
  sc <- a$soilc; n <- sc$nLayers
  zn <- microclimf:::.soilNodeDepths(sc)
  def <- data.frame(depth = zn[1:n], temp = 11, theta = (0.5 * (sc$Smin + sc$Smax))[1:n])
  b <- runpointmodel(wx, dtm = dtmcaerth, vegp = vegp, soilc = soilc, matemp = 11, runchecks = FALSE,
                     soilinit = def)
  expect_identical(a$model, b$model)
})
