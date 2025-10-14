test_that("Test soil water model", {
  obstime <- c(2024, 12, 15, 12)
  clim <- c(
    tc = -5,
    rh = 75,
    pk = 101.3,
    swdown = 100,
    difrad = 90,
    lwdown = 234,
    windspeed = 2,
    precip = 2
  )
  vegp = c(
    pai = 2,
    h = 0.5,
    ltra = 0.1,
    clump = 0.1
  )
  out <- snowoneBtest(obstime, clim, vegp, albedo = 0.95, initdepth = 0.1,
                      lat = 50, lon = -5, zref = 2)
  # format checks
  expect_true(length(out) == 15)
  expect_true(all(is.finite(out)))

  # Value checks: Combined canopy and ground
  expect_gte(out[1], -8); expect_lte(out[1], -4)    # snow surface temperature (deg C)
  expect_gte(out[2], 0); expect_lte(out[2], 0.001)  # sublimation (m SWE)
  expect_gte(out[3], 0); expect_lte(out[3], 0.0001)  # temperature melt (m SWE)
  expect_gte(out[4], 0); expect_lte(out[4], 0.0001)  # rain melt (m SWE)
  # Value checks: ground only
  expect_gte(out[5], -8); expect_lte(out[5], -4)    # Combined canopy and ground snow temperature
  expect_gte(out[6], 0); expect_lte(out[6], 0.001)  # sublimation (m SWE)
  expect_gte(out[7], 0); expect_lte(out[7], 0.0001)  # temperature melt (m SWE)
  expect_gte(out[8], 0); expect_lte(out[8], 0.0001)  # rain melt (m SWE)
  # Other variables
  expect_gte(out[9], 0); expect_lte(out[9], 1) # Canopy interception
  expect_gte(out[10], 0.2); expect_lte(out[10], 0.25) # wind friction velocity
  expect_gte(out[11], 0.1); expect_lte(out[11], 5) # Ground absorbed shortwave radiation
  expect_gte(out[12], 200); expect_lte(out[12], 300) # Ground absorbed longwave radiation
  expect_gte(out[13], 0.1); expect_lte(out[13], 0.15) # Canopy transmission
  expect_gte(out[14], 0.8); expect_lte(out[14], 1) # Convective conductance

})
