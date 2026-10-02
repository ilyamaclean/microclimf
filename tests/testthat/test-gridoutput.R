# The grid model's output fields: which hold values at each requested height,
# and what is returned for those that do not.

d <- terra::rast(dtmcaerth)
e <- terra::ext(d)
box <- terra::ext(terra::xmin(e), terra::xmin(e) + 8 * terra::res(d)[1],
                  terra::ymax(e) - 8 * terra::res(d)[2], terra::ymax(e))
crp <- function(x) terra::crop(terra::rast(x), box)
dm <- crp(dtmcaerth); vg <- lapply(vegp, crp); sl <- lapply(soilc, crp)
wx <- climdata[(170 * 24 + 1):(172 * 24), ]
radiation <- c("Rdirdown", "Rdifdown", "Rswup", "Rlwdown", "Rlwup")
allFields <- c("Tz", "tleaf", "soilm", "relhum", "windspeed", radiation)

pointrun <- function(reqhgt) {
  subsetpointmodel(suppressWarnings(runpointmodel(wx, reqhgt = reqhgt, dtm = dtmcaerth, vegp = vegp,
                                                  soilc = soilc, matemp = 11, runchecks = FALSE)), days = 2)
}
isRaster <- function(g) vapply(g, inherits, logical(1), what = "SpatRaster")
isSingleNA <- function(g) vapply(g, function(x) is.atomic(x) && length(x) == 1 && is.na(x), logical(1))

test_that("at the ground surface Tz is the ground temperature and the air fields are single NAs", {
  g <- suppressWarnings(rungridmodel(pointrun(0), reqhgt = 0, dtm = dm, vegp = vg, soilc = sl))
  expect_named(g, allFields)
  expect_identical(names(g)[isRaster(g)], c("Tz", "soilm", radiation))
  expect_identical(names(g)[isSingleNA(g)], c("tleaf", "relhum", "windspeed"))
  expect_true(all(is.finite(terra::values(g$Tz)[is.finite(terra::values(dm)), ])))
})

test_that("above ground only soil moisture is a single NA, and below ground all but Tz and soilm", {
  pm <- pointrun(0.05)
  above <- suppressWarnings(rungridmodel(pm, reqhgt = 0.05, dtm = dm, vegp = vg, soilc = sl))
  expect_identical(names(above)[isSingleNA(above)], "soilm")
  expect_identical(names(above)[isRaster(above)], setdiff(allFields, "soilm"))
  below <- suppressWarnings(rungridmodel(pointrun(-0.1), reqhgt = -0.1, dtm = dm, vegp = vg, soilc = sl))
  expect_identical(names(below)[isRaster(below)], c("Tz", "soilm"))
  expect_identical(names(below)[isSingleNA(below)], setdiff(allFields, c("Tz", "soilm")))
})

test_that("out selects fields, including ones with no value at the height, and T0 is not a field", {
  pm <- pointrun(0)
  g <- suppressWarnings(rungridmodel(pm, reqhgt = 0, dtm = dm, vegp = vg, soilc = sl,
                                     out = c(Tz = TRUE, windspeed = TRUE)))
  expect_named(g, c("Tz", "windspeed"))
  expect_s4_class(g$Tz, "SpatRaster")
  expect_identical(g$windspeed, NA)
  expect_error(rungridmodel(pm, reqhgt = 0, dtm = dm, vegp = vg, soilc = sl, out = c(T0 = TRUE)),
               "not outputs")
})
