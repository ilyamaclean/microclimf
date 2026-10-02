# The convergence tolerance is a user-facing argument wherever a point model is
# solved. These checks confirm a supplied value reaches the solver, and that
# omitting it leaves the shipped default in force.

wx <- climdata[1:48, ]
matemp <- mean(climdata$temp)

solve_at <- function(...) {
  runpointmodel(wx, reqhgt = 0, dtm = dtmcaerth, vegp = vegp, soilc = soilc,
                zin = 2, matemp = matemp, runchecks = FALSE, ...)
}

test_that("the shipped defaults are the documented ones", {
  expect_null(formals(runpointmodel)$tolerance)
  expect_equal(eval(formals(runpointmodelasgrid)$tolerance), 0.1)
  expect_null(formals(runbioclim)$tolerance)
})

test_that("omitting tolerance reproduces the default of 0.1 exactly", {
  expect_equal(solve_at()$model$Tcanopy, solve_at(tolerance = 0.1)$model$Tcanopy)
})

test_that("a supplied tolerance reaches the solver", {
  loose <- solve_at(tolerance = 0.1)
  tight <- solve_at(tolerance = 1e-4, maxIter = 400)
  expect_gt(mean(tight$model$iters), mean(loose$model$iters))
  expect_gt(max(abs(tight$model$Tcanopy - loose$model$Tcanopy)), 1e-3)
})

test_that("runbioclim forwards tolerance to its reference run", {
  ref <- microclimf:::.runbioclim1Ref
  expect_true("tolerance" %in% names(formals(ref)))
  expect_true(any(grepl("tolerance = tolerance", deparse(body(ref)), fixed = TRUE)))
})
