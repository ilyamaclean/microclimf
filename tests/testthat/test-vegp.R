# The rules by which a supplied vegp is completed before any model reads it.

ns <- asNamespace("microclimf")
d <- terra::rast(dtmcaerth)
vals1 <- function(r) { v <- terra::values(ns$.unpackRaster(r))[, 1]; v[is.nan(v)] <- NA; v }

test_that("a complete vegp resolves to itself, and resolving twice changes nothing", {
  r <- ns$.resolveVegp(vegp, d)
  for (f in names(vegp)) expect_identical(r[[f]], vegp[[f]])
  expect_identical(ns$.resolveVegp(r, d), r)
})

test_that("vegp needs hgt and pai together, or habitat, and only accepted fields", {
  expect_error(ns$.resolveVegp(list(x = vegp$x), d), "both hgt and pai, or habitat")
  expect_error(ns$.resolveVegp(list(hgt = vegp$hgt), d), "both hgt and pai, or habitat")
  expect_error(ns$.resolveVegp(c(vegp, list(vegem = vegp$x)), d), "not an accepted")
  h <- terra::rast(vegp$habitat)
  expect_error(ns$.resolveVegp(c(vegp[names(vegp) != "habitat"], list(habitat = c(h, h))), d),
               "may vary through time")
})

test_that("habitat alone gives height and plant area from the table, with a warning", {
  expect_warning(r <- ns$.resolveVegp(list(habitat = vegp$habitat), d), "very approximate")
  land <- !is.na(vals1(vegp$habitat))
  expect_true(all(is.finite(vals1(r$hgt)[land])) && all(is.finite(vals1(r$pai)[land])))
  expect_true(all(vals1(r$hgt)[land & vals1(vegp$habitat) == 16] == 0))
})

test_that("habitat is inferred from structure where it is missing", {
  hh <- terra::rast(vegp$habitat); v <- vals1(hh); i <- which(!is.na(v))[1:50]; v[i] <- NA
  terra::values(hh) <- v
  r <- ns$.resolveVegp(c(vegp[names(vegp) != "habitat"], list(habitat = hh)), d)
  expect_false(anyNA(vals1(r$habitat)[i]))
  r0 <- ns$.resolveVegp(vegp[names(vegp) != "habitat"], d)
  expect_identical(is.na(vals1(r0$habitat)), is.na(vals1(vegp$habitat)))
})

test_that("inferred habitat follows the evergreen, grass and bare rules", {
  inf <- ns$.inferHabitat
  expect_equal(inf(matrix(20, 1, 1), matrix(rep(4, 12), 1)), 2L)
  expect_equal(inf(matrix(20, 1, 1), matrix(4, 1, 1)), 4L)
  expect_equal(inf(matrix(c(0.5, 0.9), 2, 1), matrix(1, 2, 1)), c(10L, 11L))
  expect_equal(inf(matrix(0, 1, 1), matrix(0, 1, 1)), 16L)
})

test_that("a supplied leaf-physiology value replaces the table's in that cell only", {
  pftM <- matrix(c("C3", "C3", "BDT", NA), 2, 2); veg <- !is.na(pftM)
  base <- ns$.leafPhysiologyMatrices(pftM, veg, 2, 2)
  r <- terra::rast(nrows = 2, ncols = 2, vals = c(NA, 80, NA, NA))
  over <- ns$.leafPhysiologyMatrices(pftM, veg, 2, 2, list(Vcmx25 = r))
  expect_equal(over$Vcmax25[1, 2], 80e-6)
  expect_identical(over$Vcmax25[-3], base$Vcmax25[-3])
  expect_identical(over[names(over) != "Vcmax25"], base[names(base) != "Vcmax25"])
})

test_that("barren land carrying vegetation is grass, split by latitude", {
  expect_equal(habitattoPFT(c(16, 16), c(50, 10)), c("C3", "C4"))
})
