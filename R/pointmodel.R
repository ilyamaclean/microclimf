# pointmodel.R
# User-facing point-model entry points.
#
# The point model represents the vegetation/ground system at one representative
# location as a coupled surface energy and water balance. It can be driven by a
# single weather series or repeated across a coarser climate grid, with each
# climate cell parameterised from the finer vegetation, soil and terrain lying
# beneath it. Height-resolved spatial microclimate is handled by rungridmodel();
# this file provides the reference surface states that drive that calculation.

#' Check point-model inputs
#'
#' Checks that weather, vegetation, soil and terrain inputs contain the fields,
#' units and value ranges required by the point model. Implausible inputs that
#' would invalidate the calculation cause an error; unusual but possible values
#' may generate a warning.
#'
#' @param weather a data.frame of hourly weather or a multilayer raster of
#'   weather data with each raster layer named as one of the weather
#'   variables. See inbuilt dataset \code{climdata}.
#' @param vegp a list of raster layers each corresponding to vegetation inputs.
#'   See inbuilt dataset \code{vegp} and package vignette for details. The
#'   resolution and coordinate reference system of each layer must match
#'   \code{dtm}.
#' @param soilc a list of raster layers corresponding to soil inputs. See
#'   inbuilt dataset \code{soilc} and package vignette for details. The
#'   resolution and coordinate reference system of each layer must match
#'   \code{dtm}.
#' @param dtm a raster of digital elevation data with heights in m. The
#'   coordinate reference system used must preserve x and y in metres.
#' @param uzin Height (m) of the wind speed data in \code{weather}.
#'
#' @return Invisibly, \code{TRUE} if all checks pass.
#' @export
checkinputs <- function(weather, vegp, soilc, dtm, uzin = 2) {
  check.names <- function(nms, char) {
    if (!(char %in% nms)) stop("Cannot find ", char, " in weather")
  }
  check.vals <- function(x, mn, mx, char, unit) {
    if (anyNA(x)) stop("Missing values in weather$", char)
    if (any(x < mn)) stop(signif(min(x), 4), " outside range of typical ", char, " values. Units should be ", unit)
    if (any(x > mx)) stop(signif(max(x), 4), " outside range of typical ", char, " values. Units should be ", unit)
  }
  check.mean <- function(x, mn, mx, char, unit) {
    me <- mean(x, na.rm = TRUE)
    if (me < mn || me > mx) stop("Mean ", char, " of ", signif(me, 4), " implausible. Units should be ", unit)
  }
  checkRasterDims <- function(r, xy, nme) {
    d <- dim(.unpackRaster(r))
    if (d[1] != xy[1]) stop("y dimension of ", nme, " does not match dtm")
    if (d[2] != xy[2]) stop("x dimension of ", nme, " does not match dtm")
  }
  # Smallest and largest value over every cell and layer of a raster, computed
  # without copying its values into R; NULL if it holds no values.
  rasterRange <- function(r) {
    r <- .unpackRaster(r)
    lo <- terra::global(r, "min", na.rm = TRUE)[[1]]
    hi <- terra::global(r, "max", na.rm = TRUE)[[1]]
    if (all(is.na(lo))) return(NULL)
    c(min(lo, na.rm = TRUE), max(hi, na.rm = TRUE))
  }
  checkRange01 <- function(r, nme) {
    rng <- rasterRange(r)
    if (!is.null(rng) && (rng[1] < 0 || rng[2] > 1)) {
      stop(nme, " must lie in the range 0 to 1")
    }
  }

  dtm <- .unpackRaster(dtm)
  if (is.na(terra::crs(dtm)) || terra::crs(dtm) == "") {
    stop("dtm must have a coordinate reference system specified")
  }
  xy <- dim(dtm)[1:2]

  # ---- weather ----
  nms <- names(weather)
  for (v in c("obs_time", "temp", "relhum", "pres", "swdown", "difrad",
              "lwdown", "windspeed", "winddir", "precip")) {
    check.names(nms, v)
  }
  wvars <- c("temp", "relhum", "pres", "swdown", "difrad", "lwdown",
             "windspeed", "winddir", "precip")
  if (anyNA(weather[, wvars])) stop("weather contains NAs")
  tz <- attr(weather$obs_time, "tzone")
  if (is.null(tz) || tz != "UTC") stop("timezone of obs_time in weather must be UTC")
  tme <- as.POSIXlt(weather$obs_time, tz = "UTC")
  if (anyNA(tme)) stop("Cannot recognise all obs_time in weather")
  if (length(weather$temp) %% 24 != 0) {
    stop("weather needs to include data for entire days (24 hours)")
  }

  check.vals(weather$temp, -50, 65, "temperature", "deg C")
  if (anyNA(weather$relhum)) stop("Missing values in weather$relhum")
  if (any(weather$relhum < 0) || any(weather$relhum > 150)) {
    stop(signif(range(weather$relhum), 4), " outside range of typical relative humidity values. Units should be percentage (0-100)")
  }
  if (any(weather$relhum > 100)) {
    warning("weather$relhum exceeds 100 in places -- runpointmodel() does not correct this, fix the input data if this is a sensor artefact")
  }
  check.mean(weather$relhum, 5, 100, "relative humidity", "percentage (0-100)")
  check.vals(weather$pres, 30, 108.5, "pressure", "kPa ~101.3")
  check.vals(weather$swdown, 0, 1350, "shortwave radiation", "W/m^2")
  check.vals(weather$difrad, 0, 1350, "diffuse radiation", "W/m^2")
  check.vals(weather$lwdown, 0, 600, "longwave radiation", "W/m^2")
  check.vals(weather$windspeed, 0, 100, "wind speed", "m/s")
  if (uzin != 2) {
    ws10 <- weather$windspeed * log((10 - 0.08) / 0.01) / log((uzin - 0.08) / 0.01)
    if (max(ws10) > 30) {
      warning("Maximum wind speed seems quite high once adjusted to a nominal 10 m reference height. Check units are m/s and that uzin (",
              uzin, " m) is correct")
    }
  } else if (max(weather$windspeed) > 30) {
    warning("Maximum wind speed seems quite high. Check units are m/s")
  }
  mn <- min(weather$winddir); mx <- max(weather$winddir)
  if (mn < 0 || mx > 360) {
    warning("wind direction adjusted to range 0-360 using modulo operation")
  }

  # Direct radiation vs. clear-sky consistency (data-quality check, not a
  # physical constraint the model itself relies on)
  ll <- .latlongFromRaster(dtm)
  dirr <- weather$swdown - weather$difrad
  if (any(dirr < 0)) {
    warning("Diffuse radiation values higher than shortwave radiation in places; set to shortwave radiation values for the check below")
    weather$difrad[dirr < 0] <- weather$swdown[dirr < 0]
    dirr <- weather$swdown - weather$difrad
  }
  csr <- .clearskyrad(tme, ll$lat, ll$lon, weather$temp, weather$relhum, weather$pres)
  if (any(dirr > csr, na.rm = TRUE)) {
    warning("Direct radiation values higher than expected clear-sky radiation values in places")
  }
  if (any(weather$swdown > csr + 50, na.rm = TRUE)) {
    warning("Shortwave radiation values significantly higher than expected clear-sky radiation values in places")
  }

  # ---- vegp ----
  vegp <- lapply(vegp, .unpackRaster)
  .validateVegp(vegp, dtm)
  for (v in c("leafr", "leaft", "clump", "Lfrac")) {
    if (!is.null(vegp[[v]])) checkRange01(vegp[[v]], paste0("vegp$", v))
  }
  if (!is.null(vegp$leafr) && !is.null(vegp$leaft)) {
    totref <- .rasterMean(vegp$leafr) + .rasterMean(vegp$leaft)
    if (totref > 1) stop("Mean leaf reflectance + transmittance cannot be greater than one")
  }
  if (!is.null(vegp$pai)) {
    rng <- rasterRange(vegp$pai)
    if (!is.null(rng) && rng[1] < 0) stop("Minimum vegp$pai must be greater than or equal to zero")
    if (!is.null(rng) && rng[2] > 15) warning("Maximum vegp$pai of ", signif(rng[2], 4), " seems high")
  }
  if (!is.null(vegp$hgt) && !is.null(vegp$pai)) {
    nBare <- zeroBareVegCpp(terra::values(vegp$hgt, mat = TRUE),
                            terra::values(vegp$pai, mat = TRUE), countOnly = TRUE)$nBare
    if (nBare > 0) {
      warning(nBare, " cells have zero vegp$hgt or zero vegp$pai but not both; ",
              "they are treated as bare ground, with both set to zero")
    }
  }
  if (!is.null(vegp$gsmax)) {
    rng <- rasterRange(vegp$gsmax)
    if (!is.null(rng) && (rng[1] < 0 || rng[2] > 2)) {
      stop("vegp$gsmax should be in the range 0 to 2 mol/m^2/s")
    }
  }

  # ---- soilc ----
  checkRasterDims(soilc$soiltype, xy, "soilc$soiltype")
  checkRasterDims(soilc$groundr, xy, "soilc$groundr")
  checkRange01(soilc$groundr, "soilc$groundr")
  st <- terra::values(.unpackRaster(soilc$soiltype)); st <- st[!is.na(st)]
  if (length(st) > 0 && !all(st %in% microclimf::soilparamstable$Number)) {
    stop("soilc$soiltype contains a code not present in soilparamstable$Number")
  }

  invisible(TRUE)
}
#' Run the point microclimate model
#'
#' Solves the coupled radiation, turbulent exchange, vegetation energy balance,
#' soil heat and water balance, stomatal conductance and canopy interception for
#' a representative location. The vegetation is represented as a bulk exchange
#' surface rather than a vertically resolved canopy.
#'
#' The same calculation can be used in two scientifically distinct ways. With a
#' single weather data frame, one representative vegetation/soil state is
#' derived for the study area and one point simulation is run. With gridded
#' weather forcing, a separate point simulation is run for each coarse climate
#' cell using that cell's weather and the fine vegetation and soil beneath it
#' (the point simulation itself is flat). The resulting point simulations provide
#' the reference states used by
#' \code{\link{rungridmodel}} for fine-resolution spatial downscaling.
#'
#' @details
#' Each cell's plant functional type is determined from its habitat class,
#' supplied or inferred from vegetation structure, combined with latitude;
#' bare-ground cells count as grass (C3 at absolute latitude >= 30 degrees, C4
#' otherwise). The point simulation uses the most frequent type, and is bare
#' ground only if every cell is. Height, plant area index and clumping are
#' averaged over all cells. Leaf angle, leaf reflectance, transmittance and
#' maximum stomatal conductance are averaged over the modal type's cells only
#' and replace that type's defaults; its physiological parameters come from the
#' type itself.
#'
#' If the supplied weather is referenced below the vegetation height, the
#' forcing is first adjusted to a common height above the canopy by running a
#' short reference-grass proxy calculation and substituting its output into the
#' weather. This changes the reference height of the forcing, not the
#' requested soil depth.
#'
#' A \code{reqhgt} at or below ground (\code{0} = the surface) requests soil
#' temperature and water content at that depth from the solved soil profile.
#' A \code{reqhgt} at or above the resolved canopy height requests a real
#' above-canopy profile (\code{Tzabove}/\code{windzabove}/\code{RHzabove} in
#' \code{model}), evaluated at that single representative location. A
#' \code{reqhgt} strictly within the canopy (between \code{0} and canopy
#' height) is not resolved by this function -- that is a per-cell quantity, so
#' use \code{\link{runpointmodelasgrid}} for a within-canopy profile, or
#' \code{\link{rungridmodel}} for spatial microclimate at a specified height
#' above ground.
#'
#' The soil is solved as a column 1.5 times \code{totalDepth} deep (3 m by
#' default), in layers that thicken with depth. Its deepest node is not solved:
#' it is the lower boundary, held at \code{matemp} for the whole run. A
#' starting profile supplied through \code{soilinit} therefore sets the solved
#' layers only. Supplied temperatures at or below the deepest node are not
#' used, and below the deepest supplied temperature the starting profile runs
#' linearly to \code{matemp} at the deepest node. Supply \code{matemp} with a
#' temperature profile where the mean annual temperature is known. At the base
#' of the column, water drains freely (\code{FreeDrain = TRUE}) or is held
#' saturated. The early part of any run reflects its starting profiles, most
#' at depth: running a year with \code{soilend = TRUE} and passing the result
#' back as \code{soilinit} starts a run from a soil state consistent with its
#' climate.
#'
#' In gridded-weather mode, small missing coastal or masked climate cells are
#' filled by inverse-distance weighting from neighbouring valid climate cells
#' where possible. Coarse cells without a valid land footprint or with
#' unfillable forcing are returned as \code{NA}.
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
#' @param tme Used only when climate data are provided as a multi-layer raster.
#'   A POSIXlt/POSIXct time vector (UTC) corresponding to the time of each layer.
#' @param dtmc Optional coarse-gridded digital elevation dataset. Used only for
#'   elevation adjustments to weather data when climate data are provided as a
#'   multi-layer raster.
#' @param paiFlat Optional logical indicating whether plant area index is
#'   expressed per unit horizontal ground area (\code{TRUE}) or per unit local
#'   inclined ground-surface area (\code{FALSE}).
#' @param runchecks If set to \code{TRUE} (the default), the model
#'   automatically performs checks on the supplied weather data to ensure values
#'   are within range (i.e. the correct units are used) and that other inputs are
#'   provided in the correct format.
#' @param zin Height (m) of temperature and humidity data in \code{weather}.
#' @param uzin Height (m) of wind speed data in \code{weather}.
#' @param maxIter Maximum number of iterations allowed for model convergence.
#' @param tolerance Convergence tolerance (°C): the maximum allowed change in
#'   canopy or ground surface temperature between successive iterations for the
#'   solution to be considered converged.
#' @param nlayers Number of soil layers used for the soil temperature and water
#'   model.
#' @param totalDepth Total depth of the soil profile (m).
#' @param surface_organicmu Parameter controlling the concentration of soil
#'   organic matter towards the surface.
#' @param FreeDrain Optional logical indicating whether the lowest soil layer
#'   is free draining (\code{TRUE}) or not (\code{FALSE}).
#' @param matemp Temperature to which the lowest soil layer is set as a lower
#'   boundary condition. By default, mean annual temperature is automatically
#'   computed from \code{weather}.
#' @param pooling Logical indicating whether to retain water that temporarily
#'   exceeds infiltration capacity for later infiltration rather than treating it
#'   as runoff.
#' @param cores Used only when climate data are supplied as multi-layer
#'   rasters. Controls parallel processing across coarse cells: \code{off}
#'   (default) runs sequentially; \code{auto} uses one fewer than the number of
#'   available cores; \code{max} uses all available cores; alternatively, a
#'   positive integer can be supplied to specify the number of cores to use.
#' @param soilinit Optional data frame specifying the initial soil water and
#'   temperature conditions in the soil profile, with a column \code{depth}
#'   indicating depth (m), columns \code{temp} indicating soil temperature (deg
#'   C) and \code{theta} volumetric water content (m3/m3) at that depth. Depth
#'   values should be positive, and values are interpolated to the depths at
#'   which the model solves soil temperature and water. When omitted, every node
#'   is assumed to have temperature set to \code{matemp} and moisture value
#'   midway between its residual and saturated water content.
#' @param soilend Optional logical. When \code{TRUE}, returns a data.frame of
#'   the final soil water state at the end of the model run, in a form that can
#'   be passed as \code{soilinit} to initialise a subsequent run.
#'
#' @return For a single weather series, a \code{micropointv2} list containing
#'   the forcing actually used (\code{weather}) and the hourly solved surface
#'   state (\code{model}: \code{Tcanopy}/\code{Tground} (deg C), ground heat
#'   flux \code{G} (W/m2), friction velocity \code{uf} (m/s), Monin-Obukhov
#'   length \code{LL} (m), surface soil water content \code{theta0} (m3/m3),
#'   root-zone water potential \code{psi_r} (MPa), canopy water storage
#'   \code{swaterdepth} (mm)), plus \code{vegp}/\code{soilc} (the resolved parameters used),
#'   \code{pft}, \code{lat}/\code{lon}, \code{zref} (the height weather was
#'   actually run at) and \code{matemp} (the mean annual temperature used).
#'   When \code{reqhgt} requests a soil depth, \code{model} also includes
#'   \code{Tzbelow} (deg C) and \code{thetazbelow} (m3/m3) at that depth; when
#'   \code{reqhgt} is at or above the resolved canopy height, \code{model}
#'   instead includes \code{Tzabove} (deg C), \code{windzabove} (m/s) and
#'   \code{RHzabove} (percent) at that height. Within the canopy it includes
#'   neither set. For gridded weather, a
#'   list contains one such point result per valid coarse climate cell
#'   together with the grid defining those cells. With \code{soilend = TRUE},
#'   each point result also carries \code{soilend}.
#' @export
runpointmodel <- function(weather, reqhgt = 0, dtm, vegp, soilc, tme = NULL, dtmc = NULL,
                           paiFlat = FALSE, runchecks = TRUE, zin = 2, uzin = zin, maxIter = 100,
                           tolerance = NULL, nlayers = 7, totalDepth = 2,
                           surface_organicmu = 3, FreeDrain = TRUE,
                           matemp = NA_real_,
                           pooling = FALSE, cores = "off",
                           soilinit = NULL, soilend = FALSE) {
  .checkSoilinit(soilinit)
  if (is.data.frame(weather)) {
    args <- list(weather = weather, reqhgt = reqhgt, dtm = dtm, vegp = vegp, soilc = soilc,
                 runchecks = runchecks, zin = zin, uzin = uzin, maxIter = maxIter,
                 nlayers = nlayers, totalDepth = totalDepth,
                 surface_organicmu = surface_organicmu, FreeDrain = FreeDrain,
                 matemp = matemp, pooling = pooling, paiFlat = paiFlat,
                 soilinit = soilinit, soilend = soilend)
    if (!is.null(tolerance)) args$tolerance <- tolerance
    do.call(.runpointmodelp, args)
  } else if (is.list(weather)) {
    if (is.null(tme)) {
      stop("runpointmodel(): `tme` is required when `weather` is a named list of ",
           "climate rasters (array mode), not a single weather data frame")
    }
    args <- list(climarrayr = weather, tme = tme, reqhgt = reqhgt, dtm = dtm, vegp = vegp,
                 soilc = soilc, zin = zin, uzin = uzin, maxIter = maxIter, dtmc = dtmc,
                 nlayers = nlayers, totalDepth = totalDepth,
                 surface_organicmu = surface_organicmu, FreeDrain = FreeDrain,
                 matemp = matemp, runchecks = runchecks,
                 cores = cores, pooling = pooling, paiFlat = paiFlat,
                 soilinit = soilinit, soilend = soilend)
    if (!is.null(tolerance)) args$tolerance <- tolerance
    do.call(.runpointmodela, args)
  } else {
    stop("runpointmodel(): `weather` must be either a data frame of hourly weather ",
         "(single-location mode) or a named list of climate rasters (array mode, ",
         "requires `tme` too) -- got class ", paste(class(weather), collapse = "/"))
  }
}

# Collapse a spatial vegetation/soil domain to the representative state used
# by a single point simulation, then run the coupled point model at the domain
# centroid.
.runpointmodelp <- function(weather, reqhgt = 0, dtm, vegp, soilc, runchecks = TRUE,
                           zin = 2, uzin = zin, maxIter = 100,
                           nlayers = 7, totalDepth = 2,
                           surface_organicmu = 3, FreeDrain = TRUE,
                           matemp = NA_real_, tolerance = 0.1,
                           pooling = FALSE, paiFlat = FALSE, soilinit = NULL, soilend = FALSE) {
  dtm <- .unpackRaster(dtm)
  vegp <- lapply(vegp, .unpackRaster)
  soilc <- lapply(soilc, .unpackRaster)
  if (runchecks) checkinputs(weather, vegp, soilc, dtm, uzin = uzin)
  vegp <- .resolveVegp(vegp, dtm)
  ll <- .latlongFromRaster(dtm)
  lat <- ll$lat
  lon <- ll$lon

  # A point simulation needs one vegetation and soil state: the modal plant
  # functional type, with vegetation properties summarised as
  # .referenceVegetation() describes, the modal soil class and the mean ground
  # reflectance. For seasonal multi-layer inputs, use the middle layer as a
  # representative state rather than averaging phenologically distinct states
  # together.
  cellValues <- function(v) {
    if (is.null(v)) return(NULL)
    as.vector(terra::values(.unpackRaster(.middleLayerOf(v))))
  }
  refv <- .referenceVegetation(cellValues(vegp$hgt), cellValues(vegp$pai), cellValues(vegp$x),
                               cellValues(vegp$clump), cellValues(vegp$gsmax),
                               cellValues(vegp$leafr), cellValues(vegp$leaft),
                               cellValues(vegp$habitat), ll$lat,
                               phys = .physOverrides(vegp), Lfrac = cellValues(vegp$Lfrac))
  meangroundr <- .rasterMean(.middleLayerOf(soilc$groundr))
  modsoiltype <- .rasterMode(.middleLayerOf(soilc$soiltype))

  soiltab <- microclimf::soilparamstable
  s <- which(soiltab$Number == modsoiltype)
  if (length(s) == 0) stop("Unrecognised soilc$soiltype code: ", modsoiltype)
  soiltypename <- as.character(soiltab$Soil.type[s[1]])

  # Reference weather must sit above every vegetation state that may later
  # be encountered by the spatial model. Use the highest roughness-sublayer
  # clear height across the full vegetation input, rather than one from the
  # representative state used in this point simulation, as the forcing height.
  zoutOverride <- max(.maxClearHeightOf(vegp$hgt, vegp$pai), zin)

  result <- .runpointmodelCore(weather, lat, lon,
                     refv$hgt, refv$pai, refv$x, refv$clump, refv$gsmax,
                     refv$leafr, refv$leaft,
                     meangroundr, soiltypename, pft = refv$pft,
                     reqhgt = reqhgt, zin = zin, uzin = uzin,
                     nlayers = nlayers, totalDepth = totalDepth,
                     surface_organicmu = surface_organicmu, FreeDrain = FreeDrain,
                     Lfrac = refv$Lfrac, matemp = matemp, maxIter = maxIter, tolerance = tolerance,
                     pooling = pooling, zoutOverride = zoutOverride,
                     soilinit = soilinit, soilend = soilend, phys = refv$phys)
  # The point model represents a flat reference surface, so the distinction
  # between horizontal- and local-surface-area PAI affects only the later
  # slope-aware grid calculation. Preserve the user's convention on the
  # returned reference object for that calculation.
  result$paiFlat <- paiFlat
  result
}
#' Reduce a point-model run to representative days
#'
#' Selects complete 24-hour periods representing warm, cold or central thermal
#' conditions within each month or year. Keeping complete days preserves the
#' diurnal cycle needed when the reduced point-model output is used to drive a
#' spatial microclimate calculation.
#'
#' For gridded point-model output, all coarse cells are subset to the same
#' calendar days, chosen from the domain-mean canopy-temperature series. This
#' preserves a common time axis across the landscape.
#'
#' @param pointmodel Result from \code{\link{runpointmodel}}.
#' @param tstep Selection interval: \code{"month"} or \code{"year"}.
#' @param what Thermal criterion: \code{"tmax"}, \code{"tmin"} or
#'   \code{"tmedian"}.
#' @param days Optional day-of-run indices to select directly.
#' @param Tc Optional temperature series used to rank days in single-location
#'   mode. In gridded mode this is ignored: days are chosen from mean
#'   \code{Tcanopy} over valid cells.
#'
#' @return The same point-model structure with weather and model time series
#'   restricted to the selected complete days; \code{subs} records the original
#'   hourly rows retained.
#' @export
subsetpointmodel <- function(pointmodel, tstep = "month", what = "tmax", days = NA, Tc = NA) {
  if (!is.null(pointmodel$points) && !is.null(pointmodel$dtmc)) {
    .subsetpointmodela(pointmodel, tstep = tstep, what = what, days = days)
  } else if (!is.null(pointmodel$weather)) {
    .subsetpointmodelp(pointmodel, tstep = tstep, what = what, days = days, Tc = Tc)
  } else {
    stop("subsetpointmodel(): `pointmodel` must be either a single-location ",
         "\"micropointv2\" list (with $weather) or an array-mode runpointmodel() ",
         "result (with $points and $dtmc)")
  }
}

# Select complete representative days from one point-model time series and
# apply the same hourly subset to its weather and model state.
.subsetpointmodelp <- function(pointmodel, tstep = "month", what = "tmax", days = NA, Tc = NA) {
  .extractday <- function(Tc, tme, sel, what) {
    if (what == "tmax") {
      s2 <- which.max(Tc[sel])[1]
    } else if (what == "tmin") {
      s2 <- which.min(Tc[sel])[1]
    } else if (what == "tmedian") {
      o <- order(Tc[sel])
      n <- trunc(length(o) / 2)
      s2 <- o[n]
    } else stop("what must be one of tmax, tmin or tmedian")
    st <- tme[sel[s2]]
    which(tme$year == st$year & tme$mon == st$mon & tme$mday == st$mday)
  }
  model <- pointmodel$model
  if (is.logical(days)) {
    # Automatic day selection is defined only at monthly or annual scale.
    if (!(tstep %in% c("year", "month"))) {
      stop("tstep must be one of \"year\" or \"month\" when days is not supplied")
    }
    tme <- as.POSIXlt(pointmodel$weather$obs_time, tz = "UTC")
    yrs <- unique(tme$year)
    if (is.logical(Tc)) Tc <- model$Tcanopy
    if (tstep == "year") {
      sel <- which(tme$year == yrs[1])
      ai <- .extractday(Tc, tme, sel, what)
      if (length(yrs) > 1) {
        for (y in 2:length(yrs)) {
          sel <- which(tme$year == yrs[y])
          ai <- c(ai, .extractday(Tc, tme, sel, what))
        }
      }
    }
    if (tstep == "month") {
      sely <- which(tme$year == yrs[1])
      mths <- unique(tme$mon[sely])
      sel <- which(tme$mon[sely] == mths[1])
      ai <- .extractday(Tc, tme, sel, what)
      if (length(mths) > 1) {
        for (m in 2:length(mths)) {
          sel <- which(tme$mon[sely] == mths[m])
          ai <- c(ai, .extractday(Tc, tme, sel, what))
        }
      }
      if (length(yrs) > 1) {
        for (y in 2:length(yrs)) {
          sely <- which(tme$year == yrs[y])
          mths <- unique(tme$mon[sely])
          sel <- which(tme$mon[sely] == mths[1])
          am <- .extractday(Tc, tme, sel, what)
          if (length(mths) > 1) {
            for (m in 2:length(mths)) {
              sel <- which(tme$mon[sely] == mths[m])
              am <- c(am, .extractday(Tc, tme, sel, what))
            }
          }
          ai <- c(ai, am)
        }
      }
    }
  } else {
    ai <- rep((days - 1) * 24, each = 24) + rep(c(1:24), length(days))
  }
  # The length of the run the days were chosen from, over which seasonal
  # vegetation layers are spread; kept from the first subsetting.
  if (is.null(pointmodel$nhoursFull)) pointmodel$nhoursFull <- nrow(pointmodel$weather)
  pointmodel$weather <- pointmodel$weather[ai, ]
  pointmodel$model <- model[ai, ]
  pointmodel$subs <- ai
  class(pointmodel) <- "micropointv2"
  pointmodel
}
# Run the point model for spatially varying climate forcing. Each coarse climate
# cell receives its own weather series and a representative vegetation/soil state
# from the fine landscape beneath it, producing one physically independent reference
# simulation per valid climate cell.
.runpointmodela <- function(climarrayr, tme, reqhgt = 0, dtm, vegp, soilc,
                            zin = 2, uzin = zin, maxIter = 100,
                            dtmc = NULL, nlayers = 7, totalDepth = 2,
                            surface_organicmu = 3, FreeDrain = TRUE,
                            matemp = NA_real_, tolerance = 0.1,
                            runchecks = TRUE, cores = "off",
                            pooling = FALSE, paiFlat = FALSE, soilinit = NULL, soilend = FALSE) {
  dtm <- .unpackRaster(dtm)
  vegp <- lapply(vegp, .unpackRaster)
  soilc <- lapply(soilc, .unpackRaster)
  climarrayr <- lapply(climarrayr, .unpackRaster)

  climvars <- c("temp", "relhum", "pres", "swdown", "difrad",
                "lwdown", "windspeed", "winddir", "precip")
  missingvars <- setdiff(climvars, names(climarrayr))
  if (length(missingvars) > 0) {
    stop("climarrayr is missing: ", paste(missingvars, collapse = ", "))
  }
  nT <- terra::nlyr(climarrayr[[1]])
  if (length(tme) != nT) {
    stop("length(tme) must match the number of layers in climarrayr (", nT, ")")
  }

  # Each coarse climate cell receives one representative point simulation;
  # define that coarse grid explicitly or derive it from the forcing grid.
  rc <- if (!is.null(dtmc)) {
    .unpackRaster(dtmc)
  } else {
    .aggregateToCoarse(dtm, climarrayr[[1]][[1]])
  }
  if (terra::nlyr(rc) > 1) rc <- rc[[1]]
  if (isTRUE(all.equal(dim(dtm)[1:2], dim(rc)[1:2]))) {
    stop("dtm and the coarse grid (dtmc, or climarrayr's own resolution) ",
         "have the same dimensions -- runpointmodela collapses a fine dtm/",
         "vegp/soilc footprint into each coarse cell, so the two grids ",
         "must genuinely differ in resolution. For weather already at dtm's ",
         "resolution, use runpointmodelasgrid().")
  }

  if (runchecks) {
    checkinputs(.climarrayrDomainMeanWeather(climarrayr, tme), vegp, soilc, dtm, uzin = uzin)
  }
  vegp <- .resolveVegp(vegp, dtm)

  # Fill small forcing gaps over cells that have a real land footprint so a
  # coastal mask in the climate product does not automatically remove valid
  # terrestrial cells from the microclimate calculation.
  for (v in names(climarrayr)) {
    if (!identical(terra::crs(climarrayr[[v]]), terra::crs(rc))) {
      climarrayr[[v]] <- terra::project(climarrayr[[v]], terra::crs(rc))
    }
  }
  landMaskMat <- terra::as.matrix(rc, wide = TRUE)
  climarrayr <- lapply(climarrayr, function(r) {
    arr <- terra::as.array(r)
    filled <- .fillClimArrayNA(arr, landMaskMat)
    outr <- terra::rast(filled)
    terra::crs(outr) <- terra::crs(r)
    terra::ext(outr) <- terra::ext(r)
    outr
  })

  # Represent the fine vegetation/soil footprint beneath each climate cell by
  # its modal plant functional type, with vegetation properties summarised as
  # .referenceVegetation() describes, its modal soil class and its mean ground
  # reflectance. If a field contains seasonal layers, use its middle layer so
  # the reference state remains an actual vegetation state rather than an
  # average of phenologically distinct states.
  ll <- .latlonsFromCells(rc)
  physFields <- intersect(.VEGP_PHYSIOLOGY, names(vegp))
  vegFields <- c("hgt", "pai", "x", "clump", "gsmax", "leafr", "leaft", "habitat", "Lfrac", physFields)
  vegFields <- vegFields[!vapply(vegp[vegFields], is.null, logical(1))]
  vegStack <- terra::rast(lapply(vegFields, function(f) .unpackRaster(.middleLayerOf(vegp[[f]]))))
  names(vegStack) <- vegFields
  if (!identical(terra::crs(vegStack), terra::crs(rc))) {
    vegStack <- terra::project(vegStack, terra::crs(rc), method = "near")
  }
  vegVals <- terra::values(vegStack)
  byCoarse <- split(seq_len(nrow(vegVals)), .coarseCellOf(vegStack, rc))
  refV <- vector("list", terra::ncell(rc))
  for (k in names(byCoarse)) {
    i <- byCoarse[[k]]
    col <- function(f) if (f %in% vegFields) vegVals[i, f] else NULL
    refV[[as.integer(k)]] <- .referenceVegetation(col("hgt"), col("pai"), col("x"), col("clump"),
                                                  col("gsmax"), col("leafr"), col("leaft"),
                                                  col("habitat"), ll$lat[as.integer(k)],
                                                  phys = stats::setNames(lapply(physFields, function(f)
                                                    col(f) * .physMultiplier(f)), physFields),
                                                  Lfrac = col("Lfrac"))
  }
  refField <- function(f, na) vapply(refV, function(r) if (is.null(r)) na else r[[f]], na)
  meanhgtV     <- refField("hgt", NA_real_)
  pftV         <- refField("pft", NA_character_)
  meangroundrV <- terra::values(.aggregateToCoarse(.middleLayerOf(soilc$groundr), rc))[, 1]
  modsoiltypeV <- terra::values(.aggregateToCoarse(.middleLayerOf(soilc$soiltype), rc, .modeIgnoreNA))[, 1]

  # Use one forcing reference height above the tallest vegetation state in
  # the whole domain. A common height keeps the coarse reference simulations
  # physically comparable when they are later interpolated across the fine
  # grid.
  zoutOverride <- max(.maxClearHeightOf(vegp$hgt, vegp$pai), zin)

  # Resolve the dominant soil class beneath each coarse climate cell.
  soiltab <- microclimf::soilparamstable
  soiltypeNames <- vapply(modsoiltypeV, function(code) {
    if (is.na(code)) return(NA_character_)
    s <- which(soiltab$Number == code)
    if (length(s) == 0) return(NA_character_)
    as.character(soiltab$Soil.type[s[1]])
  }, character(1))

  # Only climate cells with both a vegetation footprint and a recognised soil
  # class can support a point-model reference simulation.
  validCells <- which(!is.na(meanhgtV) & !is.na(soiltypeNames))
  if (length(validCells) == 0) {
    warning("runpointmodela: no coarse cell has any valid (non-NA) dtm ",
            "coverage or resolvable soil type -- returning a list of NA")
    return(list(points = as.list(rep(NA, terra::ncell(rc))), dtmc = terra::wrap(rc), paiFlat = paiFlat))
  }

  # Extract each coarse cell's complete weather time series before running the
  # independent point simulations. Conversion to plain matrices here is an
  # execution detail needed for process-based parallelism, not part of the
  # physical calculation.
  climMat <- lapply(climarrayr, terra::values)
  tmeUTC <- as.POSIXct(tme, tz = "UTC")

  runOneCell <- function(k) {
    cellnum <- validCells[k]
    wx <- lapply(climvars, function(v) climMat[[v]][cellnum, ])
    names(wx) <- climvars
    wx <- as.data.frame(wx)
    wx$obs_time <- tmeUTC
    if (anyNA(wx[climvars])) return(NA) # unfillable gap at this cell
    .runpointmodelCore(
      wx, ll$lat[cellnum], ll$lon[cellnum],
      meanhgt = meanhgtV[cellnum], meanpai = refV[[cellnum]]$pai,
      meanx = refV[[cellnum]]$x, meanclump = refV[[cellnum]]$clump,
      meangsmax = refV[[cellnum]]$gsmax,
      meanleafr = refV[[cellnum]]$leafr, meanleaft = refV[[cellnum]]$leaft,
      meangroundr = meangroundrV[cellnum],
      soiltypename = soiltypeNames[cellnum], pft = pftV[cellnum],
      reqhgt = reqhgt, zin = zin, uzin = uzin,
      nlayers = nlayers, totalDepth = totalDepth,
      surface_organicmu = surface_organicmu, FreeDrain = FreeDrain,
      Lfrac = refV[[cellnum]]$Lfrac, matemp = matemp, maxIter = maxIter, tolerance = tolerance,
      pooling = pooling, zoutOverride = zoutOverride,
      soilinit = soilinit, soilend = soilend, soilinitWarn = (k == 1L),
      phys = refV[[cellnum]]$phys)
  }

  coresN <- .resolveCores(cores)
  oldplan <- future::plan()
  on.exit(future::plan(oldplan), add = TRUE)
  if (coresN > 1) {
    future::plan(future::multisession, workers = coresN)
  } else {
    future::plan(future::sequential)
  }
  results <- future.apply::future_lapply(seq_along(validCells), runOneCell,
                                          future.seed = TRUE)

  out <- as.list(rep(NA, terra::ncell(rc)))
  out[validCells] <- results
  # Return the coarse grid used to define the reference simulations together
  # with the PAI-area convention needed by the subsequent slope-aware spatial
  # calculation.
  list(points = out, dtmc = terra::wrap(rc), paiFlat = paiFlat)
}

# Apply one common representative-day subset across all valid coarse-cell
# point simulations, retaining the coarse spatial grid unchanged.
.subsetpointmodela <- function(pointmodela, tstep = "month", what = "tmax", days = NA) {
  if (is.null(pointmodela$points) || is.null(pointmodela$dtmc)) {
    stop("subsetpointmodel(): `pointmodel` must be an array-mode runpointmodel() ",
         "result -- a list with $points (one micropointv2 entry per coarse cell) ",
         "and $dtmc (the coarse grid those runs were assigned to)")
  }
  # Temporal subsetting does not alter the spatial grid or the PAI-area
  # convention attached to the reference simulations.
  list(points = .subsetpointmodelaCore(pointmodela$points, tstep, what, days),
       dtmc = pointmodela$dtmc, paiFlat = pointmodela$paiFlat)
}
