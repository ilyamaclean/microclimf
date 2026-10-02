# Checks on the model's building blocks, each an identity or limit the model's
# formulation requires.

ns <- asNamespace("microclimf")

test_that("the soil layer grid starts at the surface, deepens and ends at the requested depth", {
  z <- ns$geometricCpp(7L, 3)
  expect_equal(z[1], 0)
  expect_true(all(diff(z) > 0))
  expect_equal(z[8], 3)
})

test_that("roughness sublayer limits hold", {
  # No canopy, no sublayer: the clear height is ground level.
  expect_equal(ns$rslClearHeightCpp(0, 1), 0)
  # A canopy's clear height is never below its own top.
  h <- c(0.05, 0.5, 2, 20); p <- c(0.1, 1, 3, 6)
  expect_true(all(ns$rslClearHeightCpp(h, p) >= h))
  # A sublayer of depth coefficient one has no influence.
  expect_equal(ns$sublayerInfluenceCpp(1), 0)
  expect_equal(ns$momentumInfluenceCpp(0), 0)
})

test_that("displacement height lies within the canopy", {
  h <- c(0.1, 1, 10); p <- c(0.5, 2, 5)
  d <- ns$zeroplanedisCpp(h, p)
  expect_true(all(d > 0 & d < h))
})

test_that("CO2 rises with calendar year", {
  expect_true(all(diff(sapply(1980:2020, Cafromyear)) > 0))
})

test_that("a cell with zero height or zero plant area becomes bare ground", {
  z <- ns$zeroBareVegCpp(matrix(c(0, 1, 2, NA)), matrix(c(2, 0, 3, NA)))
  expect_equal(as.vector(z$hgt), c(0, 0, 2, NA))
  expect_equal(as.vector(z$pai), c(0, 0, 3, NA))
  expect_equal(z$nBare, 2L)
  expect_equal(ns$zeroBareVegCpp(matrix(c(0, 1, 2, NA)), matrix(c(2, 0, 3, NA)), countOnly = TRUE)$nBare, 2L)
})

test_that("reference vegetation counts bare cells as grass and is bare only when all cells are", {
  bare <- ns$.referenceVegetation(c(0, 0), c(0, 0), c(1, 1), c(0, 0), c(NA, NA), c(NA, NA), c(NA, NA),
                                  habitat = c(16, 16), lat = 50)
  expect_true(is.na(bare$pft))
  mixed <- ns$.referenceVegetation(c(0, 0, 10), c(0, 0, 5), c(1, 1, 1), c(0, 0, 0),
                                   c(NA, NA, NA), c(NA, NA, NA), c(NA, NA, NA),
                                   habitat = c(16, 16, 4), lat = 50)
  expect_equal(mixed$pft, "C3")
})

test_that("bioclim's representative day for each month lies in that month over several years", {
  tme <- seq(as.POSIXct("2017-01-01", tz = "UTC"), as.POSIXct("2019-12-31 23:00", tz = "UTC"), by = "hour")
  tc <- 10 + 8 * sin(2 * pi * (seq_along(tme) %% 8766) / 8766) + stats::rnorm(length(tme))
  seld <- ns$.biosel(tme, tc)$seld
  dayStart <- tme[(seld[1:12] - 1) * 24 + 1]
  expect_equal(as.POSIXlt(dayStart)$mon + 1, 1:12)
})

test_that("a starting soil profile never sets the fixed bottom node", {
  sc <- ns$createsoilc("Loam", nlayers = 7, totalDepth = 2)
  n <- sc$nLayers
  zn <- ns$.soilNodeDepths(sc)
  prof <- data.frame(depth = c(0, 1, 5), temp = c(5, 8, 30))
  out <- suppressWarnings(ns$.soilInitNodes(prof, sc, matemp = 11, warn = FALSE))
  expect_equal(out$Te[n + 1], 11)
  expect_equal(out$Te[1], 5)
  expect_equal(out$Te[3], 5 + 3 * zn[3])
  expect_null(out$theta)
  # A profile given at the node depths is taken exactly.
  exact <- data.frame(depth = zn[1:n], temp = seq(4, 10, length.out = n))
  expect_equal(ns$.soilInitNodes(exact, sc, 11, warn = FALSE)$Te[1:n], exact$temp)
  # Water is kept between oven dryness and saturation.
  wet <- data.frame(depth = c(0, 1), theta = c(0.99, 0.2))
  th <- suppressWarnings(ns$.soilInitNodes(wet, sc, 11, warn = FALSE))$theta
  expect_true(all(th <= sc$Smax))
})

test_that("starting-profile warnings fire where they should and only there", {
  sc <- ns$createsoilc("Loam", nlayers = 7, totalDepth = 2)
  w <- function(p) testthat::capture_warnings(ns$.soilInitNodes(p, sc, 11))
  deep <- w(data.frame(depth = c(0, 5), temp = c(8, 30)))
  expect_true(any(grepl("differ from matemp", deep)))
  expect_true(any(grepl("are not used", deep)))
  expect_true(any(grepl("runs linearly to matemp", w(data.frame(depth = c(0, 0.5), temp = c(8, 9))))))
  expect_true(any(grepl("clipped", w(data.frame(depth = c(0, 1), theta = c(0.99, 0.2))))))
  expect_length(w(data.frame(depth = c(0, 2.5), temp = c(8, 11))), 0)
  expect_length(w(data.frame(depth = c(0, 1), theta = c(0.3, 0.2))), 0)
  expect_error(ns$.checkSoilinit(data.frame(d = 1)), "soilinit must be")
})
