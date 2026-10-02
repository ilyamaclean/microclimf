# bioclimate.R
# Bioclimatic summaries of modelled microclimate.
#
# BIO1-BIO11 describe the thermal environment using modelled air or leaf
# temperature; BIO12-BIO19 describe moisture availability using modelled soil
# moisture rather than precipitation. To avoid running the spatial model for
# every hour of a full year, temperature and soil moisture are evaluated for a
# small set of representative days, while the wet/dry/warm/cold quarters are
# identified from the full forcing record before that reduction.
#
# The point model itself is always run on the continuous weather sequence so
# that soil heat, soil water and canopy-water state evolve chronologically.
# Representative days are selected only after that stateful calculation.

#' Compute bioclimatic variables from modelled microclimate
#'
#' Calculates 19 WorldClim-style bioclimatic variables from modelled
#' microclimate. BIO1-BIO11 summarise air or leaf temperature and BIO12-BIO19
#' summarise soil moisture. The spatial model is evaluated for 14
#' representative days rather than every day of the year.
#'
#' @details
#' \code{weather} must span at least one complete annual cycle (all 12
#' calendar months present) and be a chronologically continuous record, not a
#' pre-subsetted or pre-run point-model result.
#'
#' Fourteen representative days are selected: one day per calendar month at the
#' centre of that month's ranked daily-mean forcing temperature, plus the
#' warmest and coldest days by daily mean. When several years are supplied,
#' annual warm/cold candidates are combined across years and the median
#' candidate is selected. For spatial
#' climate forcing, one domain-mean temperature series is used for day
#' selection so every grid cell is evaluated for the same dates.
#'
#' Wettest, driest, warmest and coldest quarters (BIO8-BIO11, BIO16-BIO19) are
#' identified from a 3-month moving average of the forcing's precipitation and
#' temperature, averaged by calendar month across the full forcing record, before
#' representative-day subsetting. The point model is also run over the full
#' continuous record before its output is subset, preserving the temporal
#' state of the soil and canopy water calculations.
#'
#' @param weather a data.frame of hourly weather or a multilayer raster of
#'   weather data with each raster layer named as one of the weather
#'   variables. See inbuilt dataset \code{climdata}.
#' @param reqhgt Height above ground (m) or depth (m, negative) for which model
#'   outputs are required.
#' @param dtm a raster of digital elevation data with heights in m. The
#'   coordinate reference system used must preserve x and y in metres.
#' @param vegp a list of raster layers each corresponding to vegetation inputs.
#'   See inbuilt dataset \code{vegp} and package vignette for details. The
#'   resolution and coordinate reference system of each layer must match
#'   \code{dtm}.
#' @param soilc a list of raster layers corresponding to soil inputs. See
#'   inbuilt dataset \code{soilc} and package vignette for details. The
#'   resolution and coordinate reference system of each layer must match
#'   \code{dtm}.
#' @param tme this is only used when climate data are provided as a
#'   multi-layer raster and is a POSIXlt/POSIXct time vector (UTC)
#'   corresponding to the time of each layer.
#' @param temp Temperature used for BIO1-BIO11: \code{"air"} or \code{"leaf"}.
#'   Below ground (negative \code{reqhgt}) soil temperature is used whatever
#'   \code{temp} is.
#' @param zin Height (m) of temperature and humidity data in \code{weather}.
#' @param uzin Height (m) of the wind speed data in \code{weather}.
#' @param tfact coefficient determining sensitivity of soil moisture to variation
#'   in topographic wetness.
#' @param twiwet topographic wetness index at and above which the ground is
#'   permanently wet: soil there is saturated at every time step. Set to
#'   \code{Inf} for no permanently wet ground.
#' @param hor an optional array of horizon angles (degrees) in 24 directions.
#'   Calculated automatically from \code{dtm} if not supplied.
#' @param twi optional raster object of topographic wetness index values.
#'   Calculated automatically from \code{dtm} if not supplied. A supplied index
#'   must be the natural logarithm of upslope area per metre of contour over
#'   the tangent of slope, as returned by \code{terravars::twi()}.
#' @param svf optional raster object of sky view factor. Calculated automatically from \code{dtm} if not supplied.
#' @param slr slope in degrees. Calculated automatically from \code{dtm} if not supplied.
#' @param apr aspect in degrees. Calculated automatically from \code{dtm} if not supplied.
#' @param out Length-19 logical vector selecting BIO1-BIO19 outputs (default:
#'   all \code{TRUE}).
#' @param runchecks if set to \code{TRUE} (the default) the model automatically
#'   performs checks on the weather data supplied to ensure values are within
#'   range (i.e. correct units are used) and that other inputs are provided in
#'   the right format.
#' @param matemp temperature to which the lowest soil layer is set as a lower
#'   boundary condition. By default, mean annual temperature, automatically
#'   computed from \code{weather}.
#' @param dtmc Optional coarse gridded digital elevation dataset. Used only for
#'   elevation adjustments to weather data when climate data are provided as a
#'   multi-layer raster.
#' @param cores Parallel execution, splitting the domain into tiles:
#'   \code{"off"} (default, sequential), \code{"auto"}, \code{"max"}, or a
#'   specific number of worker processes; see \code{\link{rungridmodel}}.
#'   Tiled and untiled runs give the same result.
#' @param maxmemGB Memory budget (GB) used to choose tile size.
#' @param tolerance Convergence tolerance (degrees C): the maximum allowed
#'   change in canopy or ground surface temperature between successive
#'   iterations for the solution to be considered converged.
#' @param tilesize Optional tile width/height in pixels.
#'
#' @return A \code{SpatRaster} with one layer for each requested bioclimatic
#'   variable, named \code{bio1} to \code{bio19}.
#'
#' @export
runbioclim <- function(weather, reqhgt = 0.05, dtm, vegp, soilc, tme = NULL,
                        temp = "air", zin = 2, uzin = zin, tfact = 1.7, twiwet = 10,
                        hor = NA, twi = NA, svf = NA, slr = NA, apr = NA,
                        out = rep(TRUE, 19), runchecks = TRUE, matemp = NA_real_,
                        dtmc = NULL, cores = "off", maxmemGB = 4, tilesize = NA,
                        tolerance = NULL) {
  .forceArgs()
  if (!identical(temp, "air") && !identical(temp, "leaf")) {
    stop("runbioclim(): `temp` must be \"air\" or \"leaf\"")
  }
  if (!is.logical(out) || length(out) != 19) {
    stop("runbioclim(): `out` must be a length-19 logical vector, one entry per BIO1-BIO19")
  }
  # Below ground there are no leaves: the temperature there is the soil's.
  # The reference point model is run at that depth, so that the grid has its
  # soil temperature and moisture to scale from; above ground it is run at
  # the surface.
  belowGround <- !is.na(reqhgt) && reqhgt < 0
  if (belowGround) temp <- "air"
  refhgt <- if (belowGround) reqhgt else 0
  if (is.data.frame(weather)) {
    # Single-location forcing has one domain-wide reference microclimate.
    # Build that reference from the whole domain before tiling so spatial
    # subdivision cannot change the climate state against which cells are
    # downscaled.
    vegp <- .resolveVegp(vegp, dtm)
    ref <- .runbioclim1Ref(weather, dtm, vegp, soilc, zin, uzin, matemp, runchecks,
                           tolerance, refhgt)
    nT <- nrow(ref$pm$weather)
    worker <- function(dtm, vegp, soilc, hor, twi, svf, slr, apr) {
      list(bio = .runbioclim1PerTile(ref, reqhgt, dtm, vegp, soilc, temp, tfact, twiwet,
                                      hor, twi, svf, slr, apr, out))
    }
    result <- .tiledDispatch(dtm, vegp, soilc, hor, twi, svf, slr, apr,
                              nT = nT, cores = cores, maxmemGB = maxmemGB,
                              tilesize = tilesize, workerFn = worker)
    result$bio
  } else {
    if (is.null(tme)) {
      stop("runbioclim(): `tme` is required when `weather` is array-mode climate data ",
           "(a named list of SpatRaster climate layers)")
    }
    # Spatial forcing has one reference point-model run per coarse climate
    # cell. Resolve those references for the full domain before tiling, then
    # use the same set for every fine-resolution tile.
    vegp <- .resolveVegp(vegp, dtm)
    ref <- .runbioclim2Ref(weather, tme, dtm, vegp, soilc, zin, uzin, matemp, runchecks, dtmc,
                           tolerance, refhgt)
    validIdx <- which(vapply(ref$pma$points, is.list, logical(1)))[1]
    nT <- nrow(ref$pma$points[[validIdx]]$weather)
    worker <- function(dtm, vegp, soilc, hor, twi, svf, slr, apr) {
      list(bio = .runbioclim2PerTile(ref, reqhgt, dtm, vegp, soilc, temp, tfact, twiwet,
                                      hor, twi, svf, slr, apr, out))
    }
    result <- .tiledDispatch(dtm, vegp, soilc, hor, twi, svf, slr, apr,
                              nT = nT, cores = cores, maxmemGB = maxmemGB,
                              tilesize = tilesize, workerFn = worker)
    result$bio
  }
}

# Select 14 representative days from the forcing temperature: one
# central-ranked day from each calendar month plus representative annual
# warmest and coldest days. The returned hour indices retain all 24 hours of
# each selected day, because diurnal range is itself part of the BIO metrics.
.biosel <- function(tme, tc) {
  # Reduce the hourly forcing to one mean temperature per complete day.
  tch <- matrix(tc, ncol = 24, byrow = TRUE)
  tcd <- apply(tch, 1, mean)
  tmd <- matrix(as.numeric(tme), ncol = 24, byrow = TRUE)
  tmd <- as.POSIXlt(apply(tmd, 1, mean), origin = "1970-01-01 00:00", tz = "UTC")

  # Represent each month by a day near the centre of that month's ranked
  # daily-temperature distribution.
  sel_med <- integer(0)
  for (mth in 1:12) {
    s <- which(tmd$mon + 1 == mth)
    if (length(s) == 0) {
      stop("runbioclim(): `weather` does not cover calendar month ", mth,
           " -- bioclim requires at least one full annual cycle (see runbioclim()'s own Details)")
    }
    o <- order(tcd[s])
    n <- trunc(length(o) / 2)
    if (n < 1) n <- 1
    sel_med <- c(sel_med, s[o[n]])
  }

  # For multi-year forcing, first identify the warmest and coldest day in
  # each year, then choose a central candidate across years so one anomalous
  # year does not define the whole bioclimatic surface.
  yrs <- unique(tmd$year)
  sel_max <- integer(0)
  sel_min <- integer(0)
  for (y in seq_along(yrs)) {
    s <- which(tmd$year == yrs[y])
    sel_max <- c(sel_max, which.max(tcd[s]) + s[1] - 1)
    sel_min <- c(sel_min, which.min(tcd[s]) + s[1] - 1)
  }
  sel_max <- sel_max[order(tcd[sel_max])]
  sel_min <- sel_min[order(tcd[sel_min])]
  n <- trunc(length(sel_max) / 2) + 1
  sel_max <- sel_max[n]
  sel_min <- sel_min[n]

  # Expand the 14 selected day numbers back to complete 24-hour blocks.
  seld <- c(sel_med, sel_max, sel_min)
  selh <- rep((seld - 1) * 24, each = 24) + rep(1:24, length(seld))
  list(selh = selh, seld = seld)
}

# Return the hours belonging to the three-month window centred on a chosen
# calendar month, wrapping correctly across the end of the year.
.getselq <- function(iq, tme) {
  imn <- iq - 1; if (imn == 0) imn <- 12
  imx <- iq + 1; if (imx == 13) imx <- 1
  sel <- which((tme$mon + 1) %in% c(imn, iq, imx))
  sel[order(sel)]
}

# Identify the centre months of the wettest, driest, warmest and coldest
# three-month periods from the forcing's precipitation and air temperature,
# each averaged within calendar month so months of different length, or months
# repeated in the record, weigh alike. Quarter identity is fixed before
# representative-day selection so it reflects the full seasonal cycle.
.bioclimQuarters <- function(precip, temp, mon0) {
  aggp <- stats::aggregate(precip, by = list(mon0), mean, na.rm = TRUE)$x
  aggt <- stats::aggregate(temp, by = list(mon0), mean, na.rm = TRUE)$x
  if (length(aggp) < 12 || length(aggt) < 12) {
    stop("runbioclim(): `weather` does not cover all 12 calendar months -- ",
         "bioclim requires at least one full annual cycle")
  }
  filt <- function(x) stats::filter(x, rep(1 / 3, 3), sides = 2, circular = TRUE)
  list(wq = which.max(filt(aggp)), dq = which.min(filt(aggp)),
       hq = which.max(filt(aggt)), cq = which.min(filt(aggt)))
}

# Attach the terrain grid's spatial geometry to the per-cell BIO summaries
# returned by the C++ reduction, producing one raster layer per requested
# variable.
.bioclimToRaster <- function(bout, out, dtm) {
  dtm <- .unpackRaster(dtm)
  nms <- paste0("bio", 1:19)[out]
  layers <- lapply(bout[nms], function(m) {
    r <- terra::rast(m)
    terra::ext(r) <- terra::ext(dtm)
    terra::crs(r) <- terra::crs(dtm)
    r
  })
  r <- terra::rast(layers)
  names(r) <- nms
  r
}

# Build the single-location reference state used by every fine-resolution
# cell. Quarter identity and representative days are derived from the full
# weather record, and the point model is run for the whole domain before any
# spatial tiling is applied.
.runbioclim1Ref <- function(weather, dtm, vegp, soilc, zin, uzin, matemp, runchecks,
                            tolerance = NULL, refhgt = 0) {
  if (is.null(weather$obs_time) || is.null(weather$temp) || is.null(weather$precip)) {
    stop("runbioclim(): `weather` must have obs_time/temp/precip columns ",
         "(a full runpointmodel()-style data.frame)")
  }
  if (nrow(weather) %% 24 != 0) {
    stop("runbioclim(): `weather` must span a whole number of days (nrow a multiple of 24)")
  }
  tme <- as.POSIXlt(weather$obs_time, tz = "UTC")

  q <- .bioclimQuarters(weather$precip, weather$temp, tme$mon)
  sel <- .biosel(tme, weather$temp)

  # Initialise the deep soil thermal state from the forcing climate when no
  # independent mean annual temperature has been supplied.
  if (is.na(matemp)) matemp <- mean(weather$temp, na.rm = TRUE)

  # Soil temperature, soil water and canopy interception are stateful: each
  # hour starts from the state left by the previous one. Run the point model
  # over the continuous record first, then extract the representative days;
  # concatenating isolated days before the solve would create artificial
  # jumps in those stores.
  pmFull <- runpointmodel(weather, reqhgt = refhgt, dtm = dtm, vegp = vegp, soilc = soilc,
                           runchecks = runchecks, zin = zin, uzin = uzin, matemp = matemp,
                           tolerance = tolerance)
  pm <- subsetpointmodel(pmFull, days = sel$seld)

  tmeSub <- as.POSIXlt(pm$weather$obs_time, tz = "UTC")
  list(pm = pm,
       wetq = .getselq(q$wq, tmeSub) - 1,
       dryq = .getselq(q$dq, tmeSub) - 1,
       hotq = .getselq(q$hq, tmeSub) - 1,
       colq = .getselq(q$cq, tmeSub) - 1)
}

# Downscale the shared reference state across one fine-resolution tile and
# reduce the resulting temperature/moisture series to BIO variables.
.runbioclim1PerTile <- function(ref, reqhgt, dtm, vegp, soilc, temp, tfact, twiwet,
                                 hor, twi, svf, slr, apr, out) {
  outT <- c(Tz = identical(temp, "air"), tleaf = identical(temp, "leaf"),
            soilm = (!is.na(reqhgt) && reqhgt <= 0),
            relhum = FALSE, windspeed = FALSE, Rdirdown = FALSE, Rdifdown = FALSE,
            Rswup = FALSE, Rlwdown = FALSE, Rlwup = FALSE)
  rawT <- .gridmodelCore1(ref$pm, vegp, soilc, dtm, reqhgt = reqhgt, tfact = tfact, twiwet = twiwet,
                           hor = hor, twi = twi, svf = svf, slr = slr, apr = apr, out = outT)
  Tz <- if (identical(temp, "air")) rawT$Tz else rawT$tleaf

  # BIO12-BIO19 describe surface moisture availability, irrespective of the
  # height used for the thermal variables. Reuse soil moisture when the
  # temperature call already evaluated the ground; otherwise obtain a
  # separate ground-level moisture field.
  soilm <- if (!is.na(reqhgt) && reqhgt <= 0) {
    rawT$soilm
  } else {
    rawS <- .gridmodelCore1(ref$pm, vegp, soilc, dtm, reqhgt = 0, tfact = tfact, twiwet = twiwet,
                             hor = hor, twi = twi, svf = svf, slr = slr, apr = apr,
                             out = c(soilm = TRUE))
    rawS$soilm
  }

  bout <- runbioclimCpp(Tz, soilm, out, ref$wetq, ref$dryq, ref$hotq, ref$colq)
  .bioclimToRaster(bout, out, dtm)
}

# Build the set of coarse-cell reference states for spatial climate forcing.
# Representative dates are selected from the domain-mean forcing so every
# coarse climate cell, and therefore every fine cell, refers to the same
# parts of the annual cycle.
.runbioclim2Ref <- function(weatherA, tme, dtm, vegp, soilc, zin, uzin, matemp, runchecks, dtmc,
                            tolerance = NULL, refhgt = 0) {
  climvars <- c("temp", "relhum", "pres", "swdown", "difrad", "lwdown",
                "windspeed", "winddir", "precip")
  missingvars <- setdiff(climvars, names(weatherA))
  if (length(missingvars) > 0) {
    stop("runbioclim(): `weather` is missing: ", paste(missingvars, collapse = ", "))
  }
  weatherA <- lapply(weatherA, .unpackRaster)
  nT <- terra::nlyr(weatherA$temp)
  if (length(tme) != nT) {
    stop("runbioclim(): length(tme) must match the number of layers in `weather` (", nT, ")")
  }
  if (nT %% 24 != 0) {
    stop("runbioclim(): `weather` must span a whole number of days (nlyr a multiple of 24)")
  }
  tmeL <- as.POSIXlt(tme, tz = "UTC")

  # Use one domain-mean climate series to choose representative dates and
  # seasonal quarters, preserving a common temporal basis across space.
  tcMean <- terra::global(weatherA$temp, "mean", na.rm = TRUE)[[1]]
  precipMean <- terra::global(weatherA$precip, "mean", na.rm = TRUE)[[1]]
  q <- .bioclimQuarters(precipMean, tcMean, tmeL$mon)
  sel <- .biosel(tmeL, tcMean)

  # Use the domain-mean forcing temperature to initialise deep soil
  # temperature when no independent value is supplied.
  if (is.na(matemp)) matemp <- mean(tcMean, na.rm = TRUE)

  # Preserve the chronological state of each coarse-cell point model by
  # running the complete climate sequence first. Apply the common
  # representative-day selection only to the solved output.
  pmaFull <- runpointmodel(weatherA, tme = tme, reqhgt = refhgt, dtm = dtm, vegp = vegp,
                            soilc = soilc, runchecks = runchecks, zin = zin, uzin = uzin,
                            matemp = matemp, dtmc = dtmc, tolerance = tolerance)
  pma <- subsetpointmodel(pmaFull, days = sel$seld)

  # All valid coarse cells were subset to the same dates, so one valid time
  # vector is sufficient for mapping the selected seasonal quarters onto
  # the reduced 336-hour sequence.
  validIdx <- which(vapply(pma$points, is.list, logical(1)))[1]
  tmeSubL <- as.POSIXlt(pma$points[[validIdx]]$weather$obs_time, tz = "UTC")
  list(pma = pma,
       wetq = .getselq(q$wq, tmeSubL) - 1,
       dryq = .getselq(q$dq, tmeSubL) - 1,
       hotq = .getselq(q$hq, tmeSubL) - 1,
       colq = .getselq(q$cq, tmeSubL) - 1)
}

# Downscale the coarse-cell reference states across one fine-resolution tile
# and reduce the resulting temperature/moisture series to BIO variables.
.runbioclim2PerTile <- function(ref, reqhgt, dtm, vegp, soilc, temp, tfact, twiwet,
                                 hor, twi, svf, slr, apr, out) {
  outT <- c(Tz = identical(temp, "air"), tleaf = identical(temp, "leaf"),
            soilm = (!is.na(reqhgt) && reqhgt <= 0),
            relhum = FALSE, windspeed = FALSE, Rdirdown = FALSE, Rdifdown = FALSE,
            Rswup = FALSE, Rlwdown = FALSE, Rlwup = FALSE)
  rawT <- .gridmodelCore2(ref$pma, vegp, soilc, dtm, reqhgt = reqhgt, tfact = tfact, twiwet = twiwet,
                           hor = hor, twi = twi, svf = svf, slr = slr, apr = apr, out = outT)
  Tz <- if (identical(temp, "air")) rawT$Tz else rawT$tleaf

  soilm <- if (!is.na(reqhgt) && reqhgt <= 0) {
    rawT$soilm
  } else {
    rawS <- .gridmodelCore2(ref$pma, vegp, soilc, dtm, reqhgt = 0, tfact = tfact, twiwet = twiwet,
                             hor = hor, twi = twi, svf = svf, slr = slr, apr = apr,
                             out = c(soilm = TRUE))
    rawS$soilm
  }

  bout <- runbioclimCpp(Tz, soilm, out, ref$wetq, ref$dryq, ref$hotq, ref$colq)
  .bioclimToRaster(bout, out, dtm)
}
