# How topographic wetness distributes the reference run's soil moisture.

d <- terra::rast(dtmcaerth)
e <- terra::ext(d)
box <- terra::ext(terra::xmin(e), terra::xmin(e) + 8 * terra::res(d)[1],
                  terra::ymax(e) - 8 * terra::res(d)[2], terra::ymax(e))
crp <- function(x) terra::crop(terra::rast(x), box)
dm <- crp(dtmcaerth); vg <- lapply(vegp, crp); sl <- lapply(soilc, crp)

test_that("the anchor is the index of a flat cell with no contributing area", {
  anchor <- microclimf:::.referenceWetnessIndex(dm)
  expect_equal(anchor, log(terra::res(dm)[1] / tan(0.1 * pi / 180)))
  # An index equal to the anchor everywhere leaves the reference's moisture unshifted.
  w <- microclimf:::.topoWetness(dm, dm * 0 + anchor, tfact = 1.7, twiwet = Inf)
  expect_true(all(terra::values(w$shift) == 0, na.rm = TRUE))
  expect_true(all(terra::values(w$wet) == 0, na.rm = TRUE))
})

test_that("the wet weight is 1 at and above the threshold and falls off below it", {
  idx <- dm
  terra::values(idx) <- rep(c(6, 9, 9.5, 10, 12), length.out = terra::ncell(dm))
  w <- terra::values(microclimf:::.topoWetness(dm, idx, tfact = 1.7, twiwet = 10)$wet)[, 1]
  i <- terra::values(idx)[, 1]
  ok <- is.finite(i)
  expect_true(all(w[ok & i >= 10] == 1))
  expect_equal(unique(w[ok & i == 9.5]), exp(-1))
  expect_equal(unique(w[ok & i == 9]), exp(-2))
  expect_lt(max(w[ok & i == 6]), 1e-3)
})

test_that("ground at or above the wet threshold is saturated in every hour, and other ground is untouched", {
  wx <- climdata[(170 * 24 + 1):(172 * 24), ]
  pm <- subsetpointmodel(suppressWarnings(runpointmodel(wx, reqhgt = 0, dtm = dtmcaerth, vegp = vegp,
                                                        soilc = soilc, matemp = 11, runchecks = FALSE)), days = 2)
  # An index of 12 in the western half and 3 in the eastern half.
  idx <- dm * 0 + 3
  idx[, 1:4] <- 12
  idx <- terra::mask(idx, dm)
  run <- function(twiwet) suppressWarnings(rungridmodel(pm, reqhgt = 0, dtm = dm, vegp = vg, soilc = sl, twi = idx,
                                                        twiwet = twiwet, out = c(soilm = TRUE)))$soilm
  lim <- terra::values(run(10)); nolim <- terra::values(run(Inf))
  west <- which(terra::values(idx)[, 1] == 12); east <- which(terra::values(idx)[, 1] == 3)
  smax <- terra::values(terra::classify(sl$soiltype, cbind(soilparamstable$Number, soilparamstable$Smax)))[, 1]
  expect_equal(lim[west, ], matrix(smax[west], length(west), ncol(lim)), ignore_attr = TRUE)
  expect_true(all(nolim[west, ] < smax[west]))
  # Seven index units below the threshold the weight is exp(-14).
  expect_lt(max(abs(lim[east, ] - nolim[east, ])), 1e-6)
})
