# How topographic wetness distributes the reference run's soil moisture.

d <- terra::rast(dtmcaerth)
e <- terra::ext(d)
box <- terra::ext(terra::xmin(e), terra::xmin(e) + 8 * terra::res(d)[1],
                  terra::ymax(e) - 8 * terra::res(d)[2], terra::ymax(e))
crp <- function(x) terra::crop(terra::rast(x), box)
dm <- crp(dtmcaerth); vg <- lapply(vegp, crp); sl <- lapply(soilc, crp)

test_that("the anchor is the typical-landscape index for the cell size, and twiref replaces it", {
  anchor <- microclimf:::.referenceWetnessIndex(dm)
  expect_equal(anchor, 8 + 0.5 * log(terra::res(dm)[1] / 10))
  # An index equal to the anchor everywhere leaves the reference's moisture unshifted.
  w <- microclimf:::.topoWetness(dm, dm * 0 + anchor, tfact = 1.7, twiwet = Inf)
  expect_true(all(terra::values(w$shift) == 0, na.rm = TRUE))
  expect_true(all(terra::values(w$wet) == 0, na.rm = TRUE))
  w <- microclimf:::.topoWetness(dm, dm * 0 + 5, tfact = 2, twiwet = Inf, twiref = 7)
  expect_true(all(terra::values(w$shift) == -1, na.rm = TRUE))
})

test_that("unshifted cells keep the reference water content, and no cell is drier than air-dry soil or the reference", {
  m <- function(v) matrix(v, 1, 3)
  smax <- m(c(0.37, 0.45, 0.5)); psie <- m(c(0.73, 2.59, 3.73)); b <- m(c(1.69, 5.15, 7.63))
  theta0 <- c(0.2, 0.05)
  # The floor is set by the day's lowest relative humidity, at the air temperature of that hour.
  Ta <- c(20, 15); rh <- c(48, 90)
  psiAirDry <- 8.314 * (20 + 273.15) / 0.018015 * log(0.48)
  dry <- as.vector(smax) * (as.vector(psie) / -psiAirDry)^(1 / as.vector(b))
  # With no shift each cell takes the reference's water content.
  s0 <- microclimf:::soilmDistributeCpp(smax, psie, b, m(0), m(0), theta0, Ta, rh, 2)
  expect_equal(as.vector(s0[1, , 1]), rep(0.2, 3), tolerance = 1e-12)
  expect_equal(as.vector(s0[1, , 2]), rep(0.05, 3), tolerance = 1e-12)
  # A large negative shift stops at each soil's air-dry water content, or at the reference where drier.
  s1 <- microclimf:::soilmDistributeCpp(smax, psie, b, m(-30), m(0), theta0, Ta, rh, 2)
  expect_equal(as.vector(s1[1, , 1]), pmin(dry, 0.2), tolerance = 1e-12)
  expect_equal(as.vector(s1[1, , 2]), pmin(dry, 0.05), tolerance = 1e-12)
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
