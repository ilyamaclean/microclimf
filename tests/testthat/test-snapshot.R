# Regression snapshots: short runs compared with stored results in _snaps/.
# These fail after any change to the model's results, deliberate or not. After a
# deliberate change, review the differences, accept them with
# testthat::snapshot_accept("snapshot"), and record what changed.

fields <- c("Tcanopy", "Tground", "G", "Hout", "uf", "theta0", "psi_r", "Tzbelow")

test_that("point model: three bundled days", {
  r <- runpointmodel(climdata[4321:4392, ], reqhgt = -0.1, dtm = dtmcaerth, vegp = vegp,
                     soilc = soilc, matemp = 11, runchecks = FALSE)
  expect_snapshot_value(as.list(r$model[fields]), style = "serialize", tolerance = 1e-7)
})

test_that("point model: bare ground and a forest, two bundled days each", {
  wx <- climdata[169:216, ]
  bare <- core_run(wx, 0, 0, reqhgt = -0.1, matemp = 11)
  forest <- core_run(wx, 15, 5, reqhgt = -0.1, matemp = 11)
  expect_snapshot_value(list(bare = as.list(bare$model[fields]), forest = as.list(forest$model[fields]),
                             forestForcing = as.list(forest$weather[c("temp", "relhum", "windspeed")])),
                        style = "serialize", tolerance = 1e-7)
})

test_that("grid model: one bundled day on a corner of the landscape", {
  pm <- runpointmodel(climdata[(164 * 24 + 1):(174 * 24), ], reqhgt = 0.05, dtm = dtmcaerth,
                      vegp = vegp, soilc = soilc, matemp = 11, runchecks = FALSE)
  pms <- subsetpointmodel(pm, days = 8)
  d <- terra::rast(dtmcaerth)
  e <- terra::ext(d)
  box <- terra::ext(terra::xmin(e), terra::xmin(e) + 20 * terra::res(d)[1],
                    terra::ymax(e) - 20 * terra::res(d)[2], terra::ymax(e))
  crp <- function(x) terra::crop(terra::rast(x), box)
  g <- suppressWarnings(rungridmodel(pms, reqhgt = 0.05, dtm = crp(dtmcaerth),
                                     vegp = lapply(vegp, crp), soilc = lapply(soilc, crp), out = c(Tz = TRUE, relhum = TRUE)))
  expect_snapshot_value(lapply(g, function(x) as.vector(terra::values(x))),
                        style = "serialize", tolerance = 1e-7)
})
