# Solved soil states live at the soil nodes: every conversion between a node
# and a depth uses the node depths. Two bundled days, default column.

wx <- climdata[1:48, ]
run <- function(reqhgt, ...) {
  runpointmodel(wx, reqhgt = reqhgt, dtm = dtmcaerth, vegp = vegp, soilc = soilc, matemp = 11,
                runchecks = FALSE, ...)
}
ref <- run(-0.1, soilend = TRUE)
zn <- microclimf:::.soilNodeDepths(ref$soilc)
last <- nrow(wx)

test_that("soilend reports the node depths", {
  expect_identical(ref$soilend$depth, zn)
  expect_equal(zn[1], 0)
  expect_equal(max(zn), 1.5 * 2)
})

test_that("a node's own depth returns that node's state exactly", {
  for (i in 2:(length(zn) - 1)) {
    r <- run(-zn[i])
    expect_identical(r$model$Tzbelow[last], ref$soilend$temp[i], info = paste("node", i))
    expect_identical(r$model$thetazbelow[last], ref$soilend$theta[i], info = paste("node", i))
  }
})

test_that("between nodes the state is interpolated linearly, surface and bottom are the end nodes", {
  mid <- run(-(zn[2] + zn[3]) / 2)
  expect_equal(mid$model$Tzbelow[last], mean(ref$soilend$temp[2:3]), tolerance = 1e-12)
  surf <- run(0)
  expect_identical(surf$model$Tzbelow, surf$model$Tground)
  deep <- suppressWarnings(run(-5))
  expect_equal(deep$model$Tzbelow, rep(11, last))
})

test_that("a profile at the node depths lands on those nodes, and soilend goes back in unchanged", {
  sc <- ref$soilc
  back <- microclimf:::.soilInitNodes(ref$soilend, sc, matemp = 11, warn = FALSE)
  n <- length(zn)
  expect_equal(back$Te[1:(n - 1)], ref$soilend$temp[1:(n - 1)], tolerance = 1e-12)
  expect_equal(back$Te[n], 11)
  expect_equal(back$theta, ref$soilend$theta, tolerance = 1e-12)
})
