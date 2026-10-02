# Two tiers of tests.
#
# The quick tier runs by default: checks on the model's functions and on runs of
# a few bundled days, in well under a minute. Run it after any change:
#
#   devtools::test()
#
# The full tier adds year-long and grid runs, a minute or two more. It is
# skipped unless switched on, which also keeps it off CRAN:
#
#   Sys.setenv(MICROCLIMF_FULL_TESTS = "true"); devtools::test()
#
# devtools::test() compiles a debug build (-O0), which is slow and leaves debug
# object files in src/. To test the installed build instead (NOT_CRAN = "true"
# is what devtools::test() sets; without it testthat skips the snapshots):
#
#   Sys.setenv(NOT_CRAN = "true")
#   testthat::test_dir("tests/testthat", package = "microclimf", load_package = "installed")
#
# The snapshot tests (test-snapshot.R) compare short runs with stored results
# in _snaps/. After a deliberate change to the model's results they fail; review
# the differences, then accept the new values with
#
#   testthat::snapshot_accept("snapshot")
#
# and record what changed, so that results saved before it are not compared with
# new ones.

full_tests <- function() identical(tolower(Sys.getenv("MICROCLIMF_FULL_TESTS")), "true")

skip_unless_full <- function() {
  testthat::skip_if_not(full_tests(), "full tier: set MICROCLIMF_FULL_TESTS=true to run")
}

# One point-model run of the core, at the bundled site, for a given canopy.
core_run <- function(weather, hgt, pai, soil = "Loam", reqhgt = 0, matemp = mean(climdata$temp), ...) {
  ll <- microclimf:::.latlongFromRaster(terra::rast(dtmcaerth))
  suppressWarnings(microclimf:::.runpointmodelCore(weather, ll$lat, ll$lon,
    meanhgt = hgt, meanpai = pai, meanx = 1, meanclump = 0, meangsmax = NA_real_,
    meanleafr = 0.3, meanleaft = 0.15, meangroundr = 0.15, soiltypename = soil,
    pft = if (hgt > 0 && pai > 0) "C3" else NA_character_, reqhgt = reqhgt,
    matemp = matemp, tolerance = 0.1, ...))
}
