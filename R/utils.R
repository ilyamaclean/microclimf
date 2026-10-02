# utils.R
# Shared R-side model helpers.
#
# These functions translate spatial inputs into the physical quantities needed
# by the point and grid models: representative vegetation/soil states, solar
# geometry, plant functional types, vertical foliage profiles, topographic
# wetness and wind shelter, and coarse-to-fine climate fields. Functions near
# the end of the file handle tiling and parallel execution; those comments are
# intentionally kept separate from the scientific interpretation of the model.

# =========================================================================== #
# Generic physical and spatial helpers
# =========================================================================== #

#' Estimate atmospheric CO2 concentration from year
#'
#' Estimates atmospheric CO2 concentration from calendar year using piecewise
#' polynomial fits spanning the historical and projection periods used by the
#' model.
#'
#' @param year Calendar year (1750-2100).
#'
#' @return Estimated atmospheric CO2 concentration (ppm).
#' @export
Cafromyear <- function(year) {
  if (year > 2100) {
    stop("Cannot estimate CO2 for year > 2100\n")
  } else if (year > 2026) {
    Ca <- -0.011905 * year^2 + 51.571429 * year - 55191.666667
  } else if (year > 1962) {
    Ca <- 0.013401 * year^2 - 51.718345 * year + 50201.979441
  } else if (year >= 1750) {
    Ca <- 0.001040 * year^2 - 3.677483 * year + 3528.996720
  } else {
    stop("Cannot estimate CO2 for year < 1750\n")
  }
  return(Ca)
}

# =========================================================================== #
# Spatial summaries used to construct representative point-model states
# =========================================================================== #

# Converts a PackedSpatRaster to a SpatRaster if necessary; a plain
# SpatRaster is returned unchanged. terra's own raster functions only
# operate on SpatRaster, but model inputs are commonly stored packed
# (see terra::wrap) so they can be saved to or loaded from an R data file.
.unpackRaster <- function(r) {
  if (inherits(r, "PackedSpatRaster")) r <- terra::rast(r)
  r
}

# Mean of every value in a raster, across all cells and all layers (for a
# multi-layer, e.g. monthly, raster), ignoring NAs.
.rasterMean <- function(r) {
  r <- .unpackRaster(r)
  mean(terra::values(r), na.rm = TRUE)
}

# The highest roughness-sublayer clear height anywhere in the domain, over
# every cell and every seasonal layer. Height and plant area are paired cell by
# cell, since the clear height depends on both; a single-layer field is reused
# for each layer of a seasonal one.
.maxClearHeightOf <- function(hgt, pai) {
  maxClearHeightCpp(terra::values(.unpackRaster(hgt), mat = TRUE),
                    terra::values(.unpackRaster(pai), mat = TRUE))
}

# Select a representative layer from a seasonal vegetation/soil field.
# The middle supplied layer is used so the point reference represents an
# actual seasonal state rather than an arithmetic mixture of, for example,
# dormant and peak-canopy conditions. Single-layer inputs are unchanged.
.middleLayerOf <- function(v) {
  if (is.null(v)) return(v)
  nl <- .nLayersOf(v)
  if (nl <= 1L) return(v)
  if (inherits(v, "PackedSpatRaster")) v <- terra::rast(v)
  v[[ceiling(nl / 2)]]
}

# Modal (most frequent) value in a raster, ignoring NAs. Used for
# categorical rasters (soil type codes, habitat classes) where a mean
# would be meaningless.
.rasterMode <- function(r) {
  r <- .unpackRaster(r)
  v <- terra::values(r)
  v <- v[!is.na(v)]
  uv <- unique(v)
  uv[which.max(tabulate(match(v, uv)))]
}

# Latitude/longitude (decimal degrees, WGS84) of the centre of a raster's
# extent, reprojected from its own coordinate reference system.
.latlongFromRaster <- function(r) {
  r <- .unpackRaster(r)
  e <- terra::ext(r)
  xy <- data.frame(x = (e$xmin + e$xmax) / 2, y = (e$ymin + e$ymax) / 2)
  xy <- sf::st_as_sf(xy, coords = c("x", "y"), crs = terra::crs(r))
  ll <- sf::st_transform(xy, 4326)
  co <- sf::st_coordinates(ll)
  data.frame(lat = co[1, 2], lon = co[1, 1])
}

# Where a single weather series' radiation was measured, for solar geometry:
# the location its reference point-model run records, the centre of the whole
# study area, so that every tile of a tiled run shares it. The raster's own
# centre serves only for a reference that records none.
.referenceLatLon <- function(pointmodel, dtm) {
  if (is.numeric(pointmodel$lat) && is.numeric(pointmodel$lon) &&
      is.finite(pointmodel$lat) && is.finite(pointmodel$lon)) {
    return(data.frame(lat = pointmodel$lat, lon = pointmodel$lon))
  }
  .latlongFromRaster(dtm)
}

# Latitude/longitude (decimal degrees, WGS84) of raster-cell centres: every
# cell, or those numbered in cells. Spatial climate runs use these coordinates
# so solar geometry follows location across the climate grid rather than being
# represented by one domain centroid.
.latlonsFromCells <- function(r, cells = NULL) {
  r <- .unpackRaster(r)
  if (is.null(cells)) cells <- seq_len(terra::ncell(r))
  xy <- terra::xyFromCell(r, cells)
  ll <- sf::sf_project(sf::st_crs(terra::crs(r)), sf::st_crs(4326), xy)
  data.frame(lat = ll[, 2], lon = ll[, 1])
}

# Modal value of a numeric vector, ignoring missing values. Used when
# aggregating categorical environmental fields for which an arithmetic mean has
# no physical meaning.
.modeIgnoreNA <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# The coarse cell containing each fine cell's centre (NA outside the coarse
# grid), for two rasters in the same coordinate reference system.
.coarseCellOf <- function(fineR, coarseR) {
  xy <- terra::xyFromCell(fineR, seq_len(terra::ncell(fineR)))
  terra::cellFromXY(coarseR, xy)
}

# Summarise a fine-resolution environmental field within each coarse climate
# cell. Continuous properties normally use the mean; categorical fields can
# supply a modal summary instead. Cells with no valid underlying land remain
# NA, allowing them to be excluded from the point-reference calculations.
.aggregateToCoarse <- function(fineR, coarseR, fun = function(x) mean(x, na.rm = TRUE)) {
  fineR <- .unpackRaster(fineR)
  coarseR <- .unpackRaster(coarseR)
  if (!identical(terra::crs(fineR), terra::crs(coarseR))) {
    fineR <- terra::project(fineR, terra::crs(coarseR))
  }
  cellIdx <- .coarseCellOf(fineR, coarseR)
  vmat <- terra::values(fineR)
  nlyr <- ncol(vmat)
  # Every layer of a fine cell belongs to the same coarse cell, so the
  # cell-index vector is replicated once per layer to line up with
  # as.vector(vmat)'s column-major (cell varies fastest, then layer) order.
  cellIdxRep <- rep(cellIdx, times = nlyr)
  v <- as.vector(vmat)
  keep <- !is.na(cellIdxRep)
  agg <- tapply(v[keep], cellIdxRep[keep], fun)
  out <- rep(NA_real_, terra::ncell(coarseR))
  out[as.integer(names(agg))] <- as.numeric(agg)
  out[is.nan(out)] <- NA_real_
  outR <- coarseR
  if (terra::nlyr(outR) > 1) outR <- outR[[1]]
  terra::values(outR) <- out
  outR
}

# Fill missing climate values by inverse-distance weighting from neighbouring
# valid climate cells, but never interpolate into locations excluded by the
# land mask.
.fillClimArrayNA <- function(arr, landMask) {
  d <- dim(arr)
  filled <- fillNAIDWCpp(as.vector(arr), d[1], d[2], d[3], as.vector(landMask))
  array(filled, dim = d)
}

# =========================================================================== #
# Aerodynamic surface geometry
# =========================================================================== #

# Zero-plane displacement, roughness length and the height at which the flow
# clears the canopy's roughness sublayers are taken from the model's own C++
# (zeroplanedisCpp, roughlengthCpp, rslClearHeightCpp; src/utils.cpp), each
# vectorised over cells, so R and the model share one definition.

# =========================================================================== #
# Solar geometry and clear-sky radiation
# =========================================================================== #
# These calculations provide an independent clear-sky benchmark used to flag
# meteorological radiation inputs that are physically implausible.

# Julian day (including fractional day) from a POSIXlt time.
.jday <- function(tme) {
  yr <- tme$year + 1900
  mth <- tme$mon + 1
  dd <- tme$mday + (tme$hour + (tme$min + tme$sec / 60) / 60) / 24
  madj <- mth + (mth < 3) * 12
  yadj <- yr + (mth < 3) * -1
  jd <- trunc(365.25 * (yadj + 4716)) + trunc(30.6001 * (madj + 1)) + dd - 1524.5
  b <- 2 - trunc(yadj / 100) + trunc(trunc(yadj / 100) / 4)
  jd + (jd > 2299160) * b
}

# Solar time (decimal hours), accounting for longitude and the equation of time.
.soltime <- function(localtime, long, jd, merid = 0, dst = 0) {
  m <- 6.24004077 + 0.01720197 * (jd - 2451545)
  eot <- -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932)
  localtime + (4 * (long - merid) + eot) / 60 - dst
}

# Solar altitude (degrees above the horizon).
.solalt <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st <- .soltime(localtime, long, jd, merid, dst)
  tt <- 0.261799 * (st - 12)
  d <- (pi * 23.5 / 180) * cos(2 * pi * ((jd - 159.5) / 365.25))
  sh <- sin(d) * sin(lat * pi / 180) + cos(d) * cos(lat * pi / 180) * cos(tt)
  (180 * atan(sh / sqrt(1 - sh^2))) / pi
}

# Estimated clear-sky shortwave radiation (W/m^2) on a horizontal surface,
# from Crawford & Duchon (1999)-style atmospheric transmittance terms
# (Rayleigh/permanent-gas, water-vapour and aerosol transmittances).
.clearskyrad <- function(tme, lat, long, tc, rh, pk) {
  jd <- .jday(tme)
  lt <- tme$hour + tme$min / 60 + tme$sec / 3600
  sa <- .solalt(lt, lat, long, jd) * pi / 180
  m <- 35 * sin(sa) * ((1224 * sin(sa)^2 + 1)^(-0.5))
  TrTpg <- 1.021 - 0.084 * (m * 0.00949 * pk + 0.051)^0.5
  xx <- log(rh / 100) + ((17.27 * tc) / (237.3 + tc))
  Td <- (237.3 * xx) / (17.27 - xx)
  u <- exp(0.1133 - log(3.78) + 0.0393 * Td)
  Tw <- 1 - 0.077 * (u * m)^0.3
  Ta <- 0.935 * m
  Ic <- 1352.778 * sin(sa) * TrTpg * Tw * Ta
  Ic[is.na(Ic)] <- 0
  Ic
}

# =========================================================================== #
# Point-model biological and soil parameterisation
# =========================================================================== #

# Map the 16 habitat classes used by the package to the plant functional
# types that determine photosynthetic, stomatal and hydraulic parameters.
# Habitat class determines the plant functional type used for physiology.
# Evergreen broadleaf forest is split into tropical/temperate types at 25 degrees
# absolute latitude; savanna, grassland, cropland and mosaic classes are split
# between C4 and C3 grass at 30 degrees. The remaining habitat-to-PFT mapping is
# fixed by the lookup below; vegetation structure such as PAI is supplied
# independently and therefore still controls whether a cell is effectively bare.
#
# Several of the 16 classes have no clean one-to-one PFT match (this package
# has no "mixed forest" or "savanna" type) -- this is a first-pass best guess,
# not a definitive crosswalk. Two mappings are the weakest links: mixed forest
# defaults to needleleaf evergreen (no "mixed" PFT exists), and the closed/open
# shrubland split (evergreen/deciduous shrub) approximates what the underlying
# IGBP classes actually distinguish by cover fraction, not leaf habit.
habitattoPFT <- function(habitat, lat) {
  if (any(habitat < 1 | habitat > 16 | habitat != round(habitat))) {
    stop("habitat must be an integer between 1 and 16")
  }
  if (length(lat) == 1 && length(habitat) > 1) lat <- rep(lat, length(habitat))
  base <- c("NET", "BET.Tr", "NDT", "BDT", "NET", "ESh", "DSh", "ESh",
            "C4", "C4", "C4", "C3", "C4", "C3", "C4", "C4")
  pft <- base[habitat]
  tropicalBETthreshold <- 25
  tropicalGrassThreshold <- 30
  isBET <- habitat == 2
  pft[isBET & abs(lat) >= tropicalBETthreshold] <- "BET.Te"
  # Barren or sparsely vegetated land (16) that carries vegetation is grass.
  grassLike <- habitat %in% c(9, 10, 11, 13, 15, 16)
  pft[grassLike & abs(lat) >= tropicalGrassThreshold] <- "C3"
  return(pft)
}

# The plant functional types of PFTparams, and one variable's value for each of
# them, in that order and in the model's working units.
.PFT_NAMES <- c("BET.Tr", "BET.Te", "BDT", "NET", "NDT", "C3", "C4", "ESh", "DSh")
.pftTableRow <- function(varname) {
  tab <- microclimf::PFTparams
  s <- which(tab$varname == varname)
  as.numeric(tab[s, .PFT_NAMES]) * tab$multiplier[s]
}

# Fields a user may supply in vegp. The structure and optics fields and the
# leaf-physiology variables of PFTparams (in that table's own units) override
# the values the habitat's plant functional type would give. Only the seasonal
# fields may carry more than one layer.
.VEGP_STRUCTURE <- c("hgt", "pai", "x", "clump", "leafr", "leaft", "gsmax", "Lfrac", "habitat")
.VEGP_PHYSIOLOGY <- c("Vcmx25", "Tup", "Tlow", "Dcrit", "alpha", "rpmin", "f0", "fd",
                      "psi50", "apsi", "len", "wid")
.VEGP_SEASONAL <- c("hgt", "pai", "x", "clump", "leafr", "leaft", "Lfrac")

# Habitat class (1-16) of each cell from its canopy structure. Rows of hgtV and
# paiV are cells, columns seasonal layers; xV (leaf angle) may be NULL, and
# cells marked in skip are passed over. The plant functional type is the one
# whose reference height, plant area and leaf angle lie closest to the cell's
# (matchStructureCpp, C++ source), on the middle layer, or on the layer of
# largest plant area where the middle one is leafless. Grass is
# short grassland (10) below 0.75 m and tall grassland (11) above. Woody types
# are evergreen (needleleaf 1, broadleaf 2, closed shrubland 6) only where a
# seasonal plant area never falls below half its maximum, otherwise deciduous
# (3, 4, 7): structure alone cannot tell the two apart. A cell bare in every
# layer is barren (16). Cells without height or plant area, or passed over,
# are NA.
.inferHabitat <- function(hgtV, paiV, xV = NULL, skip = logical(0)) {
  hgtV <- as.matrix(hgtV); paiV <- as.matrix(paiV)
  n <- nrow(hgtV)
  xV <- if (is.null(xV)) matrix(NA_real_, n, 0) else as.matrix(xV)
  m <- matchStructureCpp(hgtV, paiV, xV, skip, .pftTableRow("h"), .pftTableRow("pai"),
                         .pftTableRow("x"))
  hab <- rep(NA_integer_, n)
  hab[which(m$pft == 0L)] <- 16L
  veg <- which(m$pft > 0L)
  if (length(veg) > 0) {
    pft <- .PFT_NAMES[m$pft[veg]]; h <- m$hgt[veg]; evergreen <- m$evergreen[veg]
    out <- ifelse(pft %in% c("C3", "C4"), ifelse(h < 0.75, 10L, 11L),
           ifelse(pft %in% c("BET.Tr", "BET.Te", "BDT"), ifelse(evergreen, 2L, 4L),
           ifelse(pft %in% c("NET", "NDT"), ifelse(evergreen, 1L, 3L),
                  ifelse(evergreen, 6L, 7L))))
    hab[veg] <- as.integer(out)
  }
  hab
}

# Complete a user-supplied vegp for the cells of dtm, so that every model reads
# the same resolved vegetation. vegp must hold either both hgt and pai, or
# habitat; any other field is optional. The result always holds:
#   habitat  supplied, or inferred from structure (.inferHabitat) for land
#            cells without one; single layer.
#   hgt, pai supplied, or taken from the PFT table through the cell's habitat,
#            with a warning that such values are very approximate. Barren land
#            (16) is bare unless a supplied height or plant area is positive.
#            A cell with zero height or zero plant area is bare ground.
#   x, leafr, leaft, clump, Lfrac
#            supplied, with any missing vegetated cell taken from the PFT table
#            (x, and leafr/leaft from lref/ltra), clump 0 and Lfrac 1.
# gsmax and the leaf-physiology overrides are kept as supplied. Each field keeps
# its layers and its class (packed or not). Resolving a resolved vegp returns it
# unchanged.
.resolveVegp <- function(vegp, dtm) {
  dtm <- .unpackRaster(dtm)
  vegp <- vegp[!vapply(vegp, is.null, logical(1))]
  .validateVegp(vegp, dtm)
  hasH <- !is.null(vegp$hgt); hasP <- !is.null(vegp$pai); hasHab <- !is.null(vegp$habitat)
  packed <- inherits(vegp[[1]], "PackedSpatRaster")
  R <- lapply(vegp, .unpackRaster)
  n <- terra::ncell(dtm)
  # Fields are worked on as cells x layers matrices and written back to rasters
  # only at the end. M holds height and plant area, read once, and any field
  # whose values have changed; a field set in full (rebuilt) becomes a new
  # raster on dtm, one changed in place (dirty) keeps its own.
  M <- list()
  fields <- names(R)
  touched <- character(0); rebuilt <- character(0); dirty <- character(0)
  has <- function(f) f %in% fields
  vals <- function(f) {
    if (!is.null(M[[f]])) return(M[[f]])
    v <- terra::values(R[[f]], mat = TRUE)
    if (f %in% c("hgt", "pai")) M[[f]] <<- v
    v
  }
  setField <- function(f, V) {
    M[[f]] <<- V
    if (!has(f)) fields <<- c(fields, f)
    rebuilt <<- union(rebuilt, f); touched <<- union(touched, f)
  }
  zeroBare <- function() {
    zb <- zeroBareVegCpp(vals("hgt"), vals("pai"))
    if (zb$nBare == 0) return(invisible())
    for (f in c("hgt", "pai")) if (!is.null(zb[[f]])) {
      M[[f]] <<- zb[[f]]; dirty <<- union(dirty, f)
    }
    touched <<- union(touched, c("hgt", "pai"))
  }

  # Habitat: supplied, or inferred where a land cell has structure but no class.
  hab <- rep(NA_real_, n)
  if (hasHab) { hab <- vals("habitat")[, 1]; hab[is.nan(hab)] <- NA }
  inferred <- integer(0)
  if (hasH && hasP) {
    zeroBare()
    if (anyNA(hab)) {
      new <- .inferHabitat(vals("hgt"), vals("pai"), if (has("x")) vals("x"), skip = !is.na(hab))
      inferred <- which(!is.na(new))
      hab[inferred] <- new[inferred]
    }
  }
  land <- !is.na(hab)
  if (!any(land)) return(vegp)
  if (!hasHab || length(inferred) > 0) setField("habitat", matrix(hab, ncol = 1))
  # Each land cell's plant functional type, as an index into .PFT_NAMES.
  type <- rep(NA_integer_, n)
  type[land] <- match(habitattoPFT(hab[land], .latlonsFromCells(dtm, which(land))$lat), .PFT_NAMES)

  # Height and plant area from the table where not supplied. Barren land stays
  # bare unless the supplied height or plant area says it carries vegetation.
  fillStructure <- function(f, varname, other) {
    V <- if (has(f)) vals(f) else matrix(NA_real_, n, 1)
    vegetated16 <- if (has(other)) rowAnyPositiveCpp(vals(other)) else rep(FALSE, n)
    barren <- length(.PFT_NAMES) + 1L
    typeOrBarren <- type
    typeOrBarren[land & hab == 16 & !vegetated16] <- barren
    filled <- fillByTypeCpp(V, typeOrBarren, c(.pftTableRow(varname), 0))
    if (filled$n > 0) setField(f, filled$values)
    filled$n
  }
  nH <- fillStructure("hgt", "h", "pai")
  nP <- fillStructure("pai", "pai", "hgt")
  if (nH + nP > 0) {
    warning("vegp: height or plant area index was not supplied for ", max(nH, nP), " cells ",
            "and has been taken from the plant functional type table through their habitat ",
            "class. Such values are very approximate; supply vegp$hgt and vegp$pai where ",
            "they are known.", call. = FALSE)
  }
  zeroBare()

  # Optional fields: supplied values, with missing vegetated cells from the table.
  vegetated <- land & rowAnyPositiveCpp(vals("hgt")) & rowAnyPositiveCpp(vals("pai"))
  typeVeg <- type
  typeVeg[!vegetated] <- NA
  nType <- length(.PFT_NAMES)
  defaults <- list(x = .pftTableRow("x"), leafr = .pftTableRow("lref"), leaft = .pftTableRow("ltra"),
                   clump = rep(0, nType), Lfrac = rep(1, nType))
  for (f in names(defaults)) {
    V <- if (has(f)) vals(f) else matrix(NA_real_, n, 1)
    filled <- fillByTypeCpp(V, typeVeg, defaults[[f]])
    if (filled$n == 0 && has(f)) next
    setField(f, if (filled$n > 0) filled$values else V)
  }

  out <- lapply(fields, function(f) {
    if (!(f %in% touched)) return(vegp[[f]])
    r <- R[[f]]
    if (f %in% rebuilt) {
      r <- terra::rast(dtm, nlyrs = ncol(M[[f]]))
      terra::values(r) <- M[[f]]
      names(r) <- rep(f, ncol(M[[f]]))
    } else if (f %in% dirty) {
      terra::values(r) <- M[[f]]
    }
    if (packed) terra::wrap(r) else r
  })
  names(out) <- fields
  out
}

# The structural rules a vegp must meet, whatever else it holds: accepted
# field names only, hgt and pai together or habitat, rasters with dtm's rows
# and columns, only the seasonal fields multi-layer, and habitat an integer
# class 1-16.
.validateVegp <- function(vegp, dtm) {
  dtm <- .unpackRaster(dtm)
  vegp <- vegp[!vapply(vegp, is.null, logical(1))]
  unknown <- setdiff(names(vegp), c(.VEGP_STRUCTURE, .VEGP_PHYSIOLOGY))
  if (length(unknown) > 0) {
    stop("vegp contains ", paste(unknown, collapse = ", "), ", which is not an accepted ",
         "vegetation input. Accepted: ", paste(c(.VEGP_STRUCTURE, .VEGP_PHYSIOLOGY), collapse = ", "),
         ".", call. = FALSE)
  }
  if (is.null(vegp$habitat) && (is.null(vegp$hgt) || is.null(vegp$pai))) {
    stop("vegp must contain both hgt and pai, or habitat (integer classes 1-16).", call. = FALSE)
  }
  for (f in names(vegp)) {
    r <- .unpackRaster(vegp[[f]])
    if (!inherits(r, "SpatRaster")) stop("vegp$", f, " must be a raster.", call. = FALSE)
    if (terra::nrow(r) != terra::nrow(dtm) || terra::ncol(r) != terra::ncol(dtm)) {
      stop("vegp$", f, " does not have the same rows and columns as dtm.", call. = FALSE)
    }
    if (terra::nlyr(r) > 1 && !(f %in% .VEGP_SEASONAL)) {
      stop("vegp$", f, " has ", terra::nlyr(r), " layers; only ",
           paste(.VEGP_SEASONAL, collapse = ", "), " may vary through time.", call. = FALSE)
    }
  }
  if (!is.null(vegp$habitat)) {
    hv <- terra::values(.unpackRaster(vegp$habitat))
    hv <- hv[!is.na(hv)]
    if (length(hv) > 0 && (min(hv) < 1 || max(hv) > 16 || any(hv != round(hv)))) {
      stop("vegp$habitat must contain integers between 1 and 16.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# =========================================================================== #
# Construct vegetation and soil states used by the point model
# =========================================================================== #

# Vegetation parameter list for a plant functional type, looked up from
# the inbuilt PFTparams table (see data.R) and scaled by
# PFTparams$multiplier to the model's working units. vegtype: one of
# "BET.Tr"/"BET.Te" (broadleaf evergreen tree, tropical/temperate), "BDT"
# (broadleaf deciduous tree), "NET"/"NDT" (needleleaf evergreen/deciduous
# tree), "C3"/"C4" grass, "ESh"/"DSh" (evergreen/deciduous shrub).
createvegp <- function(vegtype = "BET.Te") {
  tab <- microclimf::PFTparams
  s <- which(names(tab) == vegtype)
  v <- tab[, s] * tab$multiplier
  vegp <- as.list(v)
  names(vegp) <- tab$varname
  return(vegp)
}

# Factor converting a PFTparams variable from the table's units to the model's
# working units.
.physMultiplier <- function(varname) {
  tab <- microclimf::PFTparams
  tab$multiplier[tab$varname == varname]
}

# The leaf-physiology values supplied in vegp, one vector per variable over the
# cells in order, in the model's working units.
.physOverrides <- function(vegp) {
  f <- intersect(.VEGP_PHYSIOLOGY, names(vegp))
  out <- lapply(f, function(v) {
    x <- as.vector(terra::values(.unpackRaster(vegp[[v]])))
    x[is.nan(x)] <- NA
    x * .physMultiplier(v)
  })
  stats::setNames(out, f)
}

# Retrieve physiological and leaf-geometry parameters for a vector of PFTs,
# allowing the grid model to give each cell the photosynthetic, hydraulic and
# leaf-size properties associated with its own vegetation type.
.pftPhysiologyLookup <- function(pft, varnames) {
  type <- match(pft, .PFT_NAMES)
  if (anyNA(type)) {
    stop("Unrecognised plant functional type: ", paste(unique(pft[is.na(type)]), collapse = ", "))
  }
  out <- lapply(varnames, function(vn) .pftTableRow(vn)[type])
  names(out) <- varnames
  out
}

# Expand each cell's PFT into the physiological fields needed by the leaf
# energy/stomatal calculation. Leaf characteristic dimension is taken as the
# geometric mean of length and width -- a deliberately simpler stand-in for a
# full forced/free-convection-derived leaf geometry; all PFTs except C4 grass
# use the C3 photosynthetic pathway. PFTparams' own Vcmx25/Tlow columns are
# renamed Vcmax25/Tlw in the returned list to match toVegpstructCpp's own
# field names (C++ source).
.leafPhysiologyMatrices <- function(pftM, vegcell, rows, cols, vegp = NULL) {
  physVars <- c("Vcmx25", "Tup", "Tlow", "Dcrit", "alpha", "f0", "fd",
                "psi50", "apsi", "rpmin", "len", "wid")
  out <- stats::setNames(vector("list", length(physVars)), physVars)
  for (v in physVars) out[[v]] <- matrix(NA_real_, nrow = rows, ncol = cols)
  if (any(vegcell)) {
    looked <- .pftPhysiologyLookup(pftM[vegcell], physVars)
    for (v in physVars) out[[v]][vegcell] <- looked[[v]]
    # Values supplied in vegp replace the table's, cell by cell.
    for (v in intersect(physVars, names(vegp))) {
      m <- terra::as.matrix(.unpackRaster(vegp[[v]]), wide = TRUE) * .physMultiplier(v)
      sel <- vegcell & !is.na(m)
      out[[v]][sel] <- m[sel]
    }
  }
  leafd <- sqrt(out$len * out$wid)
  isC3 <- matrix(NA, nrow = rows, ncol = cols)
  isC3[vegcell] <- pftM[vegcell] != "C4"
  list(Vcmax25 = out$Vcmx25, Tup = out$Tup, Tlw = out$Tlow, Dcrit = out$Dcrit,
       alpha = out$alpha, f0 = out$f0, fd = out$fd, psi50 = out$psi50,
       apsi = out$apsi, rpmin = out$rpmin, leafd = leafd, isC3 = isC3)
}

# Construct the vertical soil profile for a named soil class. Campbell water-
# retention/conductivity parameters and mineral composition come from the soil
# table; organic matter is concentrated near the surface and bulk density increases
# with depth. The lower boundary either drains freely under gravity or is held
# saturated, depending on FreeDrain.
createsoilc <- function(soiltype = "Clay loam", nlayers = 15, totalDepth = 2, slope = 0, aspect = 180,
                        surface_organicmu = 3, FreeDrain = TRUE) {
  # The profile actually solved extends 1.5x beyond the caller's requested
  # totalDepth, to allow for a lower boundary layer below the depth range a
  # caller cares about.
  totalDepth <- 1.5 * totalDepth
  tab <- microclimf::soilparamstable
  s <- which(tab$Soil.type == soiltype)
  # Extract the Campbell retention/conductivity and thermal-composition
  # parameters for the selected soil type.
  campbellcols <- c("Smax", "Smin", "Ksat", "Vq", "Vm", "Vo", "Mc", "rho", "b", "psi_e")
  v <- as.numeric(tab[s, campbellcols])
  names(v) <- campbellcols
  soilc <- as.list(v)
  for (nm in campbellcols) {
    soilc[[nm]] <- rep(v[[nm]], nlayers + 1)
  }
  # Concentrate organic matter towards the surface, rescaling the other
  # solid constituents to preserve total solid volume
  mu <- rev(geometricCpp(nlayers, surface_organicmu))[2:(nlayers + 2)]
  soilc$Vo <- soilc$Vo * mu
  dif <- soilc$Vo - v[["Vo"]]
  sm <- v[["Vq"]] + v[["Vm"]]
  mu <- (sm - dif) / sm
  soilc$Vm <- soilc$Vm * mu
  soilc$Vq <- soilc$Vq * mu
  soilc$Mc <- soilc$Mc * mu
  # Bulk density increases slightly with depth
  mu <- seq(0.9, 1.1, length.out = nlayers + 1)
  soilc$rho <- soilc$rho * mu
  soilc$gref <- 0.15
  soilc$groundem <- 0.97
  soilc$grefPAR <- 0.75 * soilc$gref
  soilc$slope <- slope
  soilc$aspect <- aspect
  soilc$nLayers <- nlayers
  soilc$totalDepth <- totalDepth
  soilc$FreeDrain <- FreeDrain
  return(soilc)
}

# Run the coupled point-model physics through the complete weather sequence.
# At each timestep the model iterates radiation, atmospheric stability,
# stomatal/transpirational exchange, soil heat and water, and the surface
# energy balance until the surface temperatures converge. Soil thermal and
# hydraulic states persist between timesteps, so the sequence is genuinely
# stateful rather than a set of independent hourly equilibria.
#
# Atmospheric CO2 is represented by the concentration associated with the
# median calendar year of the forcing. matemp (mean annual temperature,
# deg C) sets the temperature of the soil model's lower boundary, the deepest
# soil node, held fixed for the whole run; any NA-supplied-by-caller fallback
# is resolved one level up, by this function's own caller, not here. A
# reqhgt at or below ground (0 = the surface) linearly interpolates this
# run's own solved soil temperature/moisture between the soil nodes that
# bracket the depth: the surface node at the surface, and the deepest node
# (matemp itself, for temperature) below the column. zAbove is a
# separate, unrelated parameter requesting a temporary above-canopy
# extrapolation, used when weather forcing has to be raised above a tall
# canopy.
#
# Returns a data frame with one row per row of climdata, including: uf
# (friction velocity, m/s), LL (Monin-Obukhov length, m), Tground/G (soil
# surface temperature and ground heat flux, from the full multi-layer soil
# heat solve), Tcanopy (the canopy/ground combined heat-exchange surface
# temperature), psi_r (root-zone mean water potential, MPa), theta0
# (surface soil water content), D_D (diurnal damping depth of the top soil
# layer) and rHa (ground-to-reference-height aerodynamic resistance) --
# D_D/rHa are returned every timestep regardless of reqhgt, since the grid
# model's ground-heat-flux shortcut needs them to evaluate a target cell's
# surface energy balance without running its own full soil solve. Plus
# Tzbelow/thetazbelow when reqhgt is at or below ground,
# Tzabove/windzabove/RHzabove when it is at or above canopy top, and
# Tabove/windAbove/RHabove when zAbove is supplied (NA otherwise).
.runbigleaf <- function(climdata, vegp, soilc, Lfrac = 1, lat, lon, zref = 2, clump = 0,
                        gsmaxCap = NA_real_,
                        reqhgt = NA_real_, matemp = NA_real_, zAbove = NA_real_,
                        maxIter = 100, tolerance = 1e-2, pooling = FALSE,
                        slope = 0, aspect = 0, svfa = 1, shelterc = 1,
                        Tinit = NULL, thetainit = NULL) {
  vegp$clump <- clump
  vegp$gsmaxCap <- ifelse(is.na(gsmaxCap), -1, gsmaxCap)
  # Forcing times are UTC. Reading them without a timezone would parse them in
  # the session's local time, where a daylight-saving gap makes R abandon the
  # time of day for every entry and read the whole record as midnight.
  obs_time <- as.POSIXlt(climdata$obs_time, tz = "UTC")
  year <- obs_time$year + 1900
  month <- obs_time$mon + 1
  day <- obs_time$mday
  hour <- obs_time$hour + obs_time$min / 60 + obs_time$sec / 3600
  Ca <- Cafromyear(stats::median(year))
  BigLeafCpp(year, month, day, hour,
             climdata$temp, climdata$relhum, climdata$pres,
             climdata$swdown, climdata$difrad, climdata$lwdown,
             climdata$windspeed, climdata$winddir, climdata$precip,
             vegp, soilc, Lfrac, lat, lon, zref, reqhgt, Ca, matemp, zAbove, maxIter, tolerance, pooling,
             slope, aspect, svfa, shelterc,
             if (is.null(Tinit)) numeric(0) else Tinit,
             if (is.null(thetainit)) numeric(0) else thetainit)
}

# Output fields of the grid models, in the order they are returned.
.GRID_OUT_FIELDS <- c("Tz", "tleaf", "soilm", "relhum", "windspeed",
                      "Rdirdown", "Rdifdown", "Rswup", "Rlwdown", "Rlwup")

# The output fields that have a value at a requested height. Below ground only
# soil temperature and moisture exist; at the ground surface there is no leaf,
# humidity or wind; above ground soil moisture is not at the requested height.
.fieldsAtHeight <- function(reqhgt) {
  radiation <- c("Rdirdown", "Rdifdown", "Rswup", "Rlwdown", "Rlwup")
  if (!is.na(reqhgt) && reqhgt < 0) return(c("Tz", "soilm"))
  if (!is.na(reqhgt) && reqhgt == 0) return(c("Tz", "soilm", radiation))
  c("Tz", "tleaf", "relhum", "windspeed", radiation)
}

# Georeference the output arrays of a grid run on the terrain grid, one raster
# layer per timestep. A field not among `live`, those with a value at the
# requested height, is returned as a single NA.
.gridOutputRasters <- function(result, dtm, live) {
  dtm <- .unpackRaster(dtm)
  out <- lapply(names(result), function(f) {
    if (!(f %in% live)) return(NA)
    r <- terra::rast(result[[f]])
    terra::ext(r) <- terra::ext(dtm)
    terra::crs(r) <- terra::crs(dtm)
    r
  })
  names(out) <- names(result)
  out
}

# Error messages for the grid model's input checks, worded for the user of
# rungridmodel() and rungridmodelbig(), whichever internal pathway raised them.
.incompletePointmodelMessage <- function() {
  paste0("`pointmodel` is incomplete: it is not an unmodified result of runpointmodel() ",
         "or subsetpointmodel(). If it was saved from an earlier version of ",
         "microclimf, run runpointmodel() again.")
}
.outNotNamedMessage <- function(valid) {
  paste0("`out` must be a named vector of TRUE or FALSE, for example ",
         "out = c(Tz = TRUE, soilm = TRUE). Valid names: ", paste(valid, collapse = ", "), ".")
}
.outUnknownMessage <- function(unknown, valid) {
  paste0("`out` contains names that are not outputs: ", paste(unknown, collapse = ", "),
         ". Valid names: ", paste(valid, collapse = ", "), ".")
}
.depthMismatchMessage <- function(reqhgt, runReqhgt) {
  ran <- if (is.null(runReqhgt) || is.na(runReqhgt)) "NA" else format(runReqhgt)
  paste0("reqhgt = ", format(reqhgt), " asks for soil conditions ", format(-reqhgt),
         " m below ground, but `pointmodel` was run with reqhgt = ", ran, ". For a depth ",
         "below ground, run runpointmodel() with the same reqhgt and pass that result.")
}
.noUsableRunMessage <- function() {
  paste0("`pointmodel` holds no usable point model run: no climate cell had land, ",
         "vegetation and soil data beneath it with complete weather. Check that the ",
         "climate data overlaps dtm, vegp and soilc.")
}
.sameResolutionMessage <- function() {
  paste0("`dtm` has the same resolution as the climate data `pointmodel` was run on. ",
         "rungridmodel() downscales climate data to a finer dtm. For climate data ",
         "already at the dtm's resolution, use runpointmodelasgrid().")
}
.notWholeDaysMessage <- function(nHours) {
  paste0("`pointmodel` covers ", nHours, " hours, which is not a whole number of days. ",
         "With vegetation that changes through the year, the point model result must ",
         "cover whole days: run runpointmodel() on whole days of weather, or use ",
         "subsetpointmodel(), which keeps whole days.")
}

# Depths (m) of the point model's soil nodes for the column `soilc` describes
# (createsoilc()): node 1 is the surface, node nLayers + 1 the fixed lower
# boundary at the column's full depth. The solved soil temperature and water
# content are the states at these depths.
.soilNodeDepths <- function(soilc) {
  n <- soilc$nLayers
  geometricCpp(n, soilc$totalDepth)[1:(n + 1)]
}

# Map a user-supplied starting soil profile onto the point model's soil nodes.
# `soilinit` is a data frame of depth (m, positive downward) with temp (deg C)
# and/or theta (m3/m3); `soilc` is the column the run solves (createsoilc()).
# Values are interpolated linearly in depth to the node depths, and held
# constant above the shallowest supplied depth.
# The deepest node is the fixed-temperature lower boundary for the whole run and
# always starts at matemp: supplied temperatures at or below it are not used, and
# below the deepest usable temperature the profile runs linearly to matemp there.
# Water below the deepest supplied depth takes that depth's value, and is kept
# between the solver's own dry bound and saturation. Returns list(Te, theta),
# either NULL where not supplied; `warn` controls the explanatory warnings.
.soilInitNodes <- function(soilinit, soilc, matemp, warn = TRUE) {
  out <- list(Te = NULL, theta = NULL)
  if (is.null(soilinit)) return(out)
  n <- soilc$nLayers
  zc <- .soilNodeDepths(soilc)
  zb <- zc[n + 1]
  dep <- soilinit$depth
  if (any(dep > zb) && warn) {
    warning("soilinit: depths below ", round(zb, 2), " m, the deepest soil node of this run, ",
            "are not used (the solved column is 1.5 x totalDepth deep).", call. = FALSE)
  }
  if (!is.null(soilinit$temp)) {
    ok <- !is.na(soilinit$temp)
    deep <- ok & dep >= zb
    if (warn && any(deep & abs(soilinit$temp - matemp) > 1)) {
      warning("soilinit: supplied temperatures at or below ", round(zb, 2), " m differ from ",
              "matemp (", round(matemp, 2), " deg C) by more than 1 K. They are not used: the ",
              "deepest soil node is the fixed lower boundary for the whole run and is held at ",
              "matemp.", call. = FALSE)
    }
    use <- ok & dep < zb
    if (any(use)) {
      d <- dep[use]; v <- soilinit$temp[use]
      Te <- if (length(unique(d)) > 1) stats::approx(d, v, xout = zc, rule = 2, ties = mean)$y
            else rep(mean(v), n + 1)
      dmax <- max(d); vmax <- mean(v[d == dmax])
      below <- zc > dmax
      Te[below] <- vmax + (matemp - vmax) * (zc[below] - dmax) / (zb - dmax)
      Te[n + 1] <- matemp
      if (warn && dmax < zc[n]) {
        warning("soilinit: the supplied temperature profile ends at ", round(dmax, 2), " m; ",
                "below that the starting temperature runs linearly to matemp (",
                round(matemp, 2), " deg C) at ", round(zb, 2), " m, the fixed lower boundary.",
                call. = FALSE)
      }
      out$Te <- Te
    } else if (warn) {
      warning("soilinit: no usable temperature above the deepest soil node; the default ",
              "starting temperature (matemp) is used.", call. = FALSE)
    }
  }
  if (!is.null(soilinit$theta)) {
    ok <- !is.na(soilinit$theta)
    if (any(ok)) {
      d <- dep[ok]; v <- soilinit$theta[ok]
      th <- if (length(unique(d)) > 1) stats::approx(d, v, xout = zc, rule = 2, ties = mean)$y
            else rep(mean(v), n + 1)
      thS <- soilc$Smax
      # Water content at oven dryness (-1e6 J/kg), the soil water solve's own lower bound.
      thDry <- thS * (1e6 / abs(soilc$psi_e))^(-1 / soilc$b)
      clipped <- th < thDry | th > thS
      th <- pmin(pmax(th, thDry), thS)
      if (warn && any(clipped)) {
        warning("soilinit: supplied water content lies outside the soil's own range (oven dry ",
                "to saturation) at some depths and was clipped into it.", call. = FALSE)
      }
      out$theta <- th
    }
  }
  out
}

# Check the form of a user-supplied starting soil profile.
.checkSoilinit <- function(soilinit) {
  if (is.null(soilinit)) return(invisible(NULL))
  if (!is.data.frame(soilinit) || is.null(soilinit$depth) ||
      (is.null(soilinit$temp) && is.null(soilinit$theta))) {
    stop("soilinit must be a data frame with a column depth (m, positive downward) and ",
         "at least one of temp (deg C) and theta (m3/m3)", call. = FALSE)
  }
  if (!is.numeric(soilinit$depth) || anyNA(soilinit$depth) || any(soilinit$depth < 0)) {
    stop("soilinit$depth must be non-negative numbers (m below the surface), without NA",
         call. = FALSE)
  }
  invisible(NULL)
}

# Representative vegetation for one point simulation, from per-cell values over
# the cells it represents. Each vegetated cell is given the plant functional type
# of its habitat class (.resolveVegp supplies a class for every land cell);
# a bare cell counts as grass (C3 poleward of 30 degrees, C4 equatorward), and
# the representative type is the most frequent one. Only if every cell is bare
# is the point bare ground. Leaf angle, leaf optical properties and maximum
# stomatal conductance are averaged over the cells of that type alone, so they
# describe the vegetation the type represents; where none of those cells has a
# value, the type's own defaults apply. Height, plant area and clumping vary
# within a type and are averaged over all cells. phys holds any supplied
# leaf-physiology values (.physOverrides), averaged like leaf angle, as is the
# living-leaf fraction Lfrac (1 where none is given).
.referenceVegetation <- function(hgt, pai, x, clump, gsmax, leafr, leaft, habitat, lat,
                                 phys = list(), Lfrac = NULL) {
  land <- !is.na(hgt)
  mn <- function(v, sel) {
    if (is.null(v)) return(NA_real_)
    m <- mean(v[sel], na.rm = TRUE)
    if (is.nan(m)) NA_real_ else m
  }
  out <- list(hgt = mn(hgt, land), pai = mn(pai, land), clump = mn(clump, land),
              pft = NA_character_, x = NA_real_, gsmax = NA_real_,
              leafr = NA_real_, leaft = NA_real_,
              phys = lapply(phys, function(v) NA_real_), Lfrac = 1)
  bare <- land & (hgt <= 0 | (!is.na(pai) & pai <= 0))
  veg <- land & !bare
  if (!any(veg)) return(out)
  pft <- rep(NA_character_, length(hgt))
  pft[bare] <- if (abs(lat) >= 30) "C3" else "C4"
  pft[veg] <- habitattoPFT(habitat[veg], lat)
  counts <- table(factor(pft[land], levels = unique(pft[land])))
  out$pft <- names(counts)[which.max(counts)]
  ofType <- veg & pft == out$pft
  out$x <- mn(x, ofType); out$gsmax <- mn(gsmax, ofType)
  out$leafr <- mn(leafr, ofType); out$leaft <- mn(leaft, ofType)
  out$phys <- lapply(phys, mn, ofType)
  if (!is.null(Lfrac)) out$Lfrac <- mn(Lfrac, ofType)
  if (is.na(out$Lfrac)) out$Lfrac <- 1
  out
}

# Construct and run one representative point-model state. The supplied
# vegetation/soil summaries describe one location or one coarse climate cell.
# This step resolves its plant functional type and physiological parameters,
# constructs the vertical soil profile, ensures the weather forcing is
# referenced above the vegetation, and then runs the coupled surface model.
.runpointmodelCore <- function(weather, lat, lon,
                                meanhgt, meanpai, meanx, meanclump, meangsmax,
                                meanleafr, meanleaft,
                                meangroundr, soiltypename, modhabitat = NA, pft = NA_character_,
                                reqhgt = 0, zin = 2, uzin = zin,
                                nlayers = 7, totalDepth = 2,
                                surface_organicmu = 3, FreeDrain = TRUE,
                                Lfrac = 1, matemp = NA_real_,
                                maxIter = 100, tolerance = 1e-2, pooling = FALSE,
                                zoutOverride = NA_real_,
                                slope = 0, aspect = 0, svfa = 1, shelterc = 1, shadow = 0,
                                soilinit = NULL, soilend = FALSE, soilinitWarn = TRUE,
                                phys = list()) {
  # Mean annual temperature anchors the deep soil thermal boundary as well as
  # the default initial temperature profile. A short weather record can give a
  # seasonally biased estimate, so warn when it is being used as a substitute
  # for an independently supplied annual mean.
  if (is.na(matemp) && soilinitWarn && !is.null(soilinit$temp)) {
    warning("soilinit supplies a starting temperature profile but matemp is not given: the ",
            "deepest soil node, the fixed lower boundary for the whole run, is held at the ",
            "mean of the supplied weather. Supply matemp if the mean annual temperature is ",
            "known.", call. = FALSE)
  }
  if (is.na(matemp)) {
    obs_time_all <- as.POSIXct(weather$obs_time, tz = "UTC")
    span_days <- as.numeric(difftime(max(obs_time_all), min(obs_time_all), units = "days"))
    if (span_days < 364) {
      warning("matemp not supplied and weather spans only ", round(span_days, 1),
              " days (less than a full year); using mean(weather$temp) as an ",
              "estimate of mean annual temperature, but this will be seasonally ",
              "biased since the supplied period doesn't cover a full annual ",
              "cycle. Supply matemp explicitly (or a full year of weather) for ",
              "a reliable below-ground deep-boundary temperature.", call. = FALSE)
    }
    matemp <- mean(weather$temp, na.rm = TRUE)
  }

  # Resolve vegetation physiology from a plant functional type supplied by the
  # caller, otherwise from the habitat class. Bare ground takes the C3
  # parameter set only so the common model structures can be constructed,
  # while its height and PAI remain zero and therefore remove vegetation
  # fluxes. A plant of zero height or zero plant area has no mass, so either
  # one zero is bare ground and both are set to zero, here as in the solver.
  bareGround <- !is.na(meanhgt) && (meanhgt <= 0 || (!is.na(meanpai) && meanpai <= 0))
  if (bareGround) {
    meanhgt <- 0
    meanpai <- 0
    pft <- "C3"
  } else if (is.na(pft)) {
    pft <- habitattoPFT(modhabitat, lat)
  }

  realvegp <- createvegp(pft)
  realvegp$h <- meanhgt
  realvegp$pai <- meanpai
  # Use the observed/supplied canopy structure and shortwave optical
  # properties for this location where they exist; PFT defaults provide them
  # otherwise, and provide the physiological parameters that cannot be
  # inferred directly from those structural inputs.
  if (!is.na(meanx)) realvegp$x <- meanx
  if (!is.na(meanleafr)) realvegp$lref <- meanleafr
  if (!is.na(meanleaft)) realvegp$ltra <- meanleaft
  # Supplied leaf-physiology values (working units) replace the type's.
  if (!bareGround) for (v in names(phys)) if (!is.na(phys[[v]])) realvegp[[v]] <- phys[[v]]

  realsoilc <- createsoilc(soiltypename, nlayers = nlayers, totalDepth = totalDepth,
                          surface_organicmu = surface_organicmu, FreeDrain = FreeDrain)
  realsoilc$gref <- meangroundr
  realsoilc$grefPAR <- 0.75 * meangroundr
  init <- .soilInitNodes(soilinit, realsoilc, matemp, warn = soilinitWarn)

  # Below totalDepth lies the column's lower boundary layer: temperature there
  # is drawn toward matemp, which the deepest node holds for the whole run, so
  # warn rather than let it be read as the soil's own behaviour at that depth.
  if (!is.na(reqhgt) && reqhgt < 0 && -reqhgt > totalDepth) {
    zb <- max(.soilNodeDepths(realsoilc))
    warning("reqhgt requests a depth of ", -reqhgt, " m, below totalDepth (", totalDepth,
            " m). Below totalDepth the soil temperature is drawn toward matemp, which ",
            "the deepest point of the column (", round(zb, 2), " m) holds for the whole ",
            "run", if (-reqhgt >= zb) ", and at this depth it is matemp" else "",
            ". Increase totalDepth if the soil's own behaviour at this depth is needed.",
            call. = FALSE)
  }

  # Meteorological forcing must be referenced above the canopy's roughness
  # sublayer. When the height at which the flow clears that sublayer (or a
  # caller-supplied common height) lies above the weather height, translate the
  # forcing up to it so all point references contributing to a spatial run share
  # a consistent atmospheric height.
  weather_adj <- weather
  zref <- zin
  zout <- if (!is.na(zoutOverride)) zoutOverride else
    max(rslClearHeightCpp(meanhgt, meanpai), zin)
  if (zout > zin) {
    grasshgt <- 0.12
    grasspai <- 1
    if (uzin != zin) {
      d0 <- zeroplanedisCpp(grasshgt, grasspai)
      zm0 <- roughlengthCpp(grasshgt, grasspai, d0)
      lnr <- log((zin - d0) / zm0) / log((uzin - d0) / zm0)
      weather_adj$windspeed <- weather_adj$windspeed * lnr
    }
    grassvegp <- createvegp("C3")
    grassvegp$h <- grasshgt
    grassvegp$pai <- grasspai
    grasssoilc <- createsoilc(soiltypename, nlayers = nlayers, totalDepth = totalDepth)
    # Use a short reference-grass surface calculation to obtain the
    # stability-dependent temperature, humidity and wind profiles needed to
    # translate the forcing from its measurement height to the common
    # above-canopy reference height.
    grassrun <- .runbigleaf(weather_adj, grassvegp, grasssoilc, Lfrac = 1,
                           lat = lat, lon = lon, zref = zin, zAbove = zout, matemp = matemp,
                           maxIter = maxIter, tolerance = tolerance)
    weather_adj$temp <- grassrun$Tabove
    weather_adj$relhum <- grassrun$RHabove
    weather_adj$windspeed <- grassrun$windAbove
    zref <- zout
  }

  # Run the actual vegetation/soil system once the atmospheric forcing has
  # been placed at its appropriate reference height. Terrain modifiers are
  # available for specialised callers; ordinary point-reference runs are flat
  # and unsheltered, with fine-scale terrain effects applied by the grid model.
  # In hours the terrain shades the surface its shortwave is diffuse alone; this
  # is applied after the forcing height adjustment, which describes open air.
  weather_run <- weather_adj
  shaded <- rep_len(as.numeric(shadow), nrow(weather_run)) > 0
  weather_run$swdown[shaded] <- weather_run$difrad[shaded]
  model <- .runbigleaf(weather_run, realvegp, realsoilc, Lfrac = Lfrac,
                      lat = lat, lon = lon, zref = zref, clump = meanclump,
                      gsmaxCap = meangsmax, reqhgt = reqhgt, matemp = matemp,
                      maxIter = maxIter, tolerance = tolerance, pooling = pooling,
                      slope = slope, aspect = aspect, svfa = svfa, shelterc = shelterc,
                      Tinit = init$Te, thetainit = init$theta)
  TeEnd <- attr(model, "TeEnd"); thetaEnd <- attr(model, "thetaEnd")
  attr(model, "TeEnd") <- NULL; attr(model, "thetaEnd") <- NULL

  # Above-reference profile fields are diagnostic products of the temporary
  # forcing-height calculation, not part of the physical output of the main
  # vegetation/soil simulation.
  model$Tabove <- NULL
  model$windAbove <- NULL
  model$RHabove <- NULL

  # Soil state at reqhgt exists only at or below the ground surface, and the
  # above-canopy profile only at or above canopy top. The fields of whichever
  # does not apply hold no values and are not returned; within the canopy
  # neither applies.
  belowGround <- !is.na(reqhgt) && reqhgt <= 0
  aboveCanopy <- !is.na(reqhgt) && !belowGround && reqhgt >= meanhgt
  if (!belowGround) model[c("Tzbelow", "thetazbelow")] <- NULL
  if (!aboveCanopy) model[c("Tzabove", "windzabove", "RHzabove")] <- NULL

  out <- list(weather = weather_adj, model = model, vegp = realvegp,
              soilc = realsoilc, pft = pft, lat = lat, lon = lon, zref = zref,
              matemp = matemp, reqhgt = reqhgt, Lfrac = Lfrac)
  if (isTRUE(soilend)) {
    # The soil state at the end of the run, in the form soilinit takes, so a
    # run can be started again from it, at the node depths.
    out$soilend <- data.frame(depth = .soilNodeDepths(realsoilc),
                              temp = TeEnd, theta = thetaEnd)
  }
  class(out) <- "micropointv2"
  out
}

# Domain-mean forcing used for broad plausibility checks on gridded weather.
# Spatial structure is checked separately; this summary is not used as the
# climate forcing for any modelled cell.
.climarrayrDomainMeanWeather <- function(climarrayr, tme) {
  vars <- names(climarrayr)
  out <- lapply(vars, function(v) {
    vmat <- terra::values(climarrayr[[v]])
    apply(vmat, 2, mean, na.rm = TRUE)
  })
  names(out) <- vars
  out <- as.data.frame(out)
  out$obs_time <- as.POSIXct(tme, tz = "UTC")
  out
}

# Choose representative days from the mean canopy-temperature signal across
# valid coarse climate cells, then apply those same dates to every cell. This
# preserves a common temporal axis across the landscape rather than allowing
# each location to represent a different calendar day.
.subsetpointmodelaCore <- function(points, tstep, what, days) {
  isValid <- vapply(points, is.list, logical(1))
  validCells <- which(isValid)
  if (length(validCells) == 0) {
    stop("subsetpointmodela(): `pointmodela$points` has no valid (non-NA) coarse ",
         "cells -- nothing to subset")
  }

  Tc <- 0
  for (i in validCells) Tc <- Tc + points[[i]]$model$Tcanopy
  Tc <- Tc / length(validCells)

  outPoints <- as.list(rep(NA, length(points)))
  for (i in validCells) {
    outPoints[[i]] <- .subsetpointmodelp(points[[i]], tstep = tstep, what = what, days = days, Tc = Tc)
  }
  outPoints
}

# =========================================================================== #
# Grid-model spatial parameterisation
# =========================================================================== #

# Plant area index above height z for the generic gamma-shaped vertical foliage
# profile used here. At ground level all plant area lies above the requested height;
# at or above canopy top none does. This quantity controls how much vegetation can
# intercept radiation between the requested height and the sky.
.paiaboveheight <- function(z, hgt, pai, shape = 1.5, rate = shape / 7) {
  zc <- pmin(pmax(z, 0), hgt) # clamp to the canopy's own physical extent
  xx <- ((hgt - zc) / hgt) * 10
  stats::pgamma(xx, shape, rate) * (pai / stats::pgamma(10, shape, rate))
}

# Vertical foliage distributions used to determine how much plant area lies
# above a requested within-canopy height. Broadleaf, needleleaf and shrub
# profiles are precomputed stand-level templates expressed against relative
# canopy height; each cell rescales the appropriate template by its own canopy
# height and total PAI. Grass uses a strongly basal power-law distribution
# with an analytic cumulative profile.
#
# The woody templates come from simulating a stand of many individual plants
# (the Lalic & Mihailovic leaf-density shape for trees, a softer reasoned
# version for shrubs, each plus its own wood/branching term peaked a little
# above the ground) with realistic height/allometry variance and a genuine
# understory population, then averaging -- the averaging itself, not any
# hand-tuned smoothing, is what makes a stand's aggregate profile smoother
# than any one plant's own sharp profile. That simulation is far too
# expensive to re-run per grid cell, so its result is precomputed once per
# group as a normalised template and reused via cheap interpolation +
# rescale by each cell's own height/PAI, not by re-simulating at runtime.
.PFT_GROUP_OF <- c(BET.Tr = "broadleaf", BET.Te = "broadleaf", BDT = "broadleaf",
                    NET = "needleleaf", NDT = "needleleaf",
                    C3 = "grass", C4 = "grass",
                    ESh = "shrub", DSh = "shrub")

# Relative canopy-height coordinates shared by the woody vegetation templates.
.PFT_TEMPLATE_ZFRAC <- c(0, 0.0101, 0.0202, 0.0303, 0.0404, 0.05051, 0.06061, 0.07071, 0.08081, 0.09091, 0.10101, 0.11111, 0.12121, 0.13131, 0.14141, 0.15152, 0.16162, 0.17172, 0.18182, 0.19192, 0.20202, 0.21212, 0.22222, 0.23232, 0.24242, 0.25253, 0.26263, 0.27273, 0.28283, 0.29293, 0.30303, 0.31313, 0.32323, 0.33333, 0.34343, 0.35354, 0.36364, 0.37374, 0.38384, 0.39394, 0.40404, 0.41414, 0.42424, 0.43434, 0.44444, 0.45455, 0.46465, 0.47475, 0.48485, 0.49495, 0.50505, 0.51515, 0.52525, 0.53535, 0.54545, 0.55556, 0.56566, 0.57576, 0.58586, 0.59596, 0.60606, 0.61616, 0.62626, 0.63636, 0.64646, 0.65657, 0.66667, 0.67677, 0.68687, 0.69697, 0.70707, 0.71717, 0.72727, 0.73737, 0.74747, 0.75758, 0.76768, 0.77778, 0.78788, 0.79798, 0.80808, 0.81818, 0.82828, 0.83838, 0.84848, 0.85859, 0.86869, 0.87879, 0.88889, 0.89899, 0.90909, 0.91919, 0.92929, 0.93939, 0.94949, 0.9596, 0.9697, 0.9798, 0.9899, 1)

# Fraction of total PAI above each relative height. Each template declines
# from the whole canopy at ground level to zero at canopy top and is
# interpolated onto the requested relative height for each grid cell.
# Computed once from a 600-tree-plus-understory Monte Carlo stand simulation
# at a representative reference height/PAI per group (broadleaf: h=20,
# pai=4; needleleaf: h=22.5, pai=3; shrub: h=1.25, pai=2), reused for any
# cell's own height/PAI via the normalised coordinates above.
.PFT_PAIA_TEMPLATE <- list(
  broadleaf = c(1, 0.993149, 0.986122, 0.978997, 0.971798, 0.964628, 0.957565, 0.950528, 0.94346, 0.93632, 0.929113, 0.921919, 0.914701, 0.907464, 0.900267, 0.893052, 0.885793, 0.878451, 0.871097, 0.863769, 0.856593, 0.849598, 0.84269, 0.835888, 0.829117, 0.822306, 0.815446, 0.808528, 0.801544, 0.7945, 0.787388, 0.780239, 0.773036, 0.765765, 0.758395, 0.75093, 0.743373, 0.735717, 0.727951, 0.720059, 0.712029, 0.703874, 0.69558, 0.687132, 0.678511, 0.669714, 0.660741, 0.651582, 0.642217, 0.632631, 0.62283, 0.612814, 0.602573, 0.592097, 0.581365, 0.570365, 0.559162, 0.547817, 0.536358, 0.524744, 0.512896, 0.500809, 0.488467, 0.475858, 0.462948, 0.449706, 0.436141, 0.422238, 0.408028, 0.393471, 0.37853, 0.363275, 0.347882, 0.332551, 0.317161, 0.301663, 0.286166, 0.270637, 0.255201, 0.239973, 0.22489, 0.209981, 0.195147, 0.180347, 0.165716, 0.15139, 0.137557, 0.124355, 0.111549, 0.0991144, 0.0872262, 0.0757946, 0.0648039, 0.0543387, 0.0444261, 0.0349423, 0.0258465, 0.0170637, 0.00846546, 0),
  needleleaf = c(1, 0.993569, 0.986983, 0.980317, 0.973602, 0.966934, 0.960294, 0.953573, 0.946769, 0.939864, 0.932842, 0.925766, 0.918602, 0.911369, 0.904112, 0.896766, 0.889305, 0.881698, 0.874, 0.866246, 0.858518, 0.850841, 0.843149, 0.835442, 0.827663, 0.819758, 0.811721, 0.803546, 0.795225, 0.786762, 0.778147, 0.769403, 0.760518, 0.751476, 0.742259, 0.73287, 0.723306, 0.713565, 0.70364, 0.693517, 0.683188, 0.672666, 0.661939, 0.650999, 0.639838, 0.628454, 0.616847, 0.605017, 0.592952, 0.580644, 0.568109, 0.555358, 0.542394, 0.529229, 0.515874, 0.502341, 0.488656, 0.474839, 0.460894, 0.446803, 0.43258, 0.418242, 0.4038, 0.38927, 0.374658, 0.359977, 0.345244, 0.330476, 0.315671, 0.300812, 0.28595, 0.271205, 0.256699, 0.242569, 0.228855, 0.215456, 0.202414, 0.189723, 0.177369, 0.165348, 0.153666, 0.142325, 0.131299, 0.120599, 0.110261, 0.100295, 0.0907121, 0.0815112, 0.0726858, 0.0642582, 0.0562664, 0.0486659, 0.0414412, 0.0346058, 0.0281597, 0.0220008, 0.0160912, 0.0104534, 0.00509001, 0),
  shrub = c(1, 0.989557, 0.979002, 0.96834, 0.957582, 0.946741, 0.935817, 0.924812, 0.913718, 0.902538, 0.891272, 0.879918, 0.868477, 0.85695, 0.845328, 0.833617, 0.821824, 0.809937, 0.79796, 0.78589, 0.773722, 0.761463, 0.74911, 0.736663, 0.724118, 0.711477, 0.698739, 0.685904, 0.672976, 0.659959, 0.646855, 0.633672, 0.620406, 0.607059, 0.593636, 0.580142, 0.566583, 0.552963, 0.53929, 0.525569, 0.511808, 0.498015, 0.484199, 0.470368, 0.456533, 0.442704, 0.428891, 0.415105, 0.401359, 0.387668, 0.374043, 0.360498, 0.347052, 0.333718, 0.320507, 0.307438, 0.29452, 0.281764, 0.269187, 0.256803, 0.244622, 0.232659, 0.220932, 0.209448, 0.198225, 0.187285, 0.176634, 0.166282, 0.156256, 0.146555, 0.137184, 0.12816, 0.11948, 0.111143, 0.103159, 0.0955278, 0.0882415, 0.081304, 0.0747168, 0.0684666, 0.0625505, 0.0569713, 0.0517084, 0.0467518, 0.0421021, 0.0377345, 0.0336365, 0.0298038, 0.0262174, 0.0228631, 0.0197347, 0.0168191, 0.014102, 0.0115776, 0.00923958, 0.00707517, 0.00507787, 0.00324274, 0.00155378, 0)
)

# Closed-form cumulative plant area above height z for the grass profile used
# by the model. Density is assumed proportional to (1/(z/h)^n) - 1, n = 0.2:
# a genuine, integrable (for n<1) singularity at z=0 tapering to exactly 0 at
# z=h, matching real grass structure (densest at the base, nothing at the
# tallest blade tips). The cumulative integral below was verified against
# numeric integration to ~1e-9.
.paiaAboveHeightGrass <- function(z, hgt, pai, n = 0.2) {
  t <- pmin(pmax(z, 0), hgt) / hgt
  frac <- (1 - t^(1 - n) - (1 - n) * (1 - t)) / n
  ifelse(t >= 1, 0, pai * frac)
}

# Plant area above a requested height for each cell, using the vertical
# foliage profile associated with that cell's structural vegetation group.
.paiaboveheightPFT <- function(z, hgt, pai, group) {
  out <- hgt
  out[] <- NA_real_
  for (g in names(.PFT_PAIA_TEMPLATE)) {
    sel <- group == g
    if (!any(sel)) next
    zfrac <- pmin(pmax(z, 0), hgt[sel]) / hgt[sel]
    frac <- stats::approx(.PFT_TEMPLATE_ZFRAC, .PFT_PAIA_TEMPLATE[[g]], xout = zfrac,
                           rule = 2)$y
    out[sel] <- pai[sel] * frac
  }
  selg <- group == "grass"
  if (any(selg)) out[selg] <- .paiaAboveHeightGrass(z, hgt[selg], pai[selg])
  out
}

# Map soil classes to their residual and saturated volumetric water contents
# so topographic redistribution of soil moisture respects each cell's own
# physically attainable range.
.soilcSminSmax <- function(soiltypeR) {
  if (inherits(soiltypeR, "PackedSpatRaster")) soiltypeR <- terra::rast(soiltypeR)
  lookup <- microclimf::soilparamstable[, c("Number", "Smin", "Smax")]
  list(Smin = terra::classify(soiltypeR, rcl = cbind(lookup$Number, lookup$Smin)),
       Smax = terra::classify(soiltypeR, rcl = cbind(lookup$Number, lookup$Smax)))
}

# Surface-layer Campbell soil parameters required by the grid model's local
# ground-energy approximation: solid composition controls thermal properties,
# while saturation, air-entry potential and the Campbell exponent determine
# moisture-dependent conductivity and surface humidity. psi_e is passed
# through here in the soil table's own positive-magnitude convention; the
# sign is applied downstream, at the point of use in gridmodel::miniSoilpCpp,
# not here.
.soilcCampbellParams <- function(soiltypeR) {
  if (inherits(soiltypeR, "PackedSpatRaster")) soiltypeR <- terra::rast(soiltypeR)
  lookup <- microclimf::soilparamstable[, c("Number", "Vq", "Vm", "Vo", "Mc", "Smax", "psi_e", "b")]
  list(Vq    = terra::classify(soiltypeR, rcl = cbind(lookup$Number, lookup$Vq)),
       Vm    = terra::classify(soiltypeR, rcl = cbind(lookup$Number, lookup$Vm)),
       Vo    = terra::classify(soiltypeR, rcl = cbind(lookup$Number, lookup$Vo)),
       Mc    = terra::classify(soiltypeR, rcl = cbind(lookup$Number, lookup$Mc)),
       thetaS = terra::classify(soiltypeR, rcl = cbind(lookup$Number, lookup$Smax)),
       psie  = terra::classify(soiltypeR, rcl = cbind(lookup$Number, lookup$psi_e)),
       b     = terra::classify(soiltypeR, rcl = cbind(lookup$Number, lookup$b)))
}

# Number of seasonal states represented by a vegetation field. Scalars,
# matrices and single-layer rasters represent one fixed state; multi-layer rasters
# represent a sequence of states through time.
.nLayersOf <- function(v) {
  if (inherits(v, "PackedSpatRaster")) v <- terra::rast(v)
  if (inherits(v, "SpatRaster")) return(as.integer(terra::nlyr(v)))
  1L
}

# Return one vegetation state as a rows x cols matrix. Fixed fields return
# the same values for every requested state; seasonal rasters return the selected
# layer.
.tomatLayer <- function(v, layer, rows, cols) {
  if (inherits(v, "PackedSpatRaster")) v <- terra::rast(v)
  if (inherits(v, "SpatRaster")) return(terra::as.matrix(v[[layer]], wide = TRUE))
  if (is.matrix(v)) return(v)
  matrix(v, nrow = rows, ncol = cols)
}

# Smallest slope (degrees) allowed in the topographic wetness index, which
# divides by the tangent of slope.
.TWI_MIN_SLOPE <- 0.1

# Topographic wetness index computed from the terrain. Sinks are filled before
# upslope area is accumulated, so that flow is not stopped at every pit and
# stream channels keep their upslope area; flow is shared among all downhill
# neighbours (multiple flow directions). Slope is taken from the unfilled
# surface, so that a filled hollow is not treated as flat ground.
.wetnessIndex <- function(dtm) {
  dtm <- .unpackRaster(dtm)
  fill_ext <- terravars:::.fillna_dtm(dtm)
  slope <- terra::crop(terra::terrain(fill_ext$extended, "slope", unit = "degrees"), dtm)
  terravars::twi(terravars::fillsinks(dtm), method = "modified", route = "mfd", slope = slope,
                 min_slope = .TWI_MIN_SLOPE, fillna = TRUE)
}

# Topographic wetness index of the reference run. The reference is a flat cell
# with no contributing area: its upslope area per metre of contour is one cell
# width, and its slope is the smallest the index allows. A cell with this
# index takes the reference run's soil moisture.
.referenceWetnessIndex <- function(dtm) {
  log(terra::res(.unpackRaster(dtm))[1] / tan(.TWI_MIN_SLOPE * pi / 180))
}

# Width, in units of the topographic wetness index, over which ground just
# below the wet threshold is joined to saturated ground at and above it.
.TWI_WET_JOIN <- 0.5

# The two terms by which topographic wetness redistributes reference soil
# moisture across the landscape.
#   shift: the log-scaled topographic wetness index less the reference run's
#     (.referenceWetnessIndex), divided by tfact. Positive values identify
#     cells topographically wetter than the reference run and negative values
#     drier ones.
#   wet: the share of the gap to saturation closed in every hour. It is 1 at
#     and above the wet threshold twiwet, where the ground is permanently
#     saturated, and falls off exponentially below it over .TWI_WET_JOIN, so
#     that cells either side of the threshold do not differ by a step.
# Supplied TWI is treated as raw and log-transformed when its magnitude
# indicates an unlogged index -- log-TWI rarely exceeds ~25-30 while raw TWI
# reaches into the hundreds or thousands, so a value above 50 is taken as raw.
# An index not supplied is computed from dtm (.wetnessIndex).
.topoWetness <- function(dtm, twi, tfact, twiwet) {
  if (is.logical(twi)) {
    ltwi <- .wetnessIndex(dtm)
  } else {
    if (inherits(twi, "PackedSpatRaster")) twi <- terra::rast(twi)
    maxtwi <- max(terra::values(twi), na.rm = TRUE)
    ltwi <- if (maxtwi > 50) log(twi) else twi
  }
  list(shift = (ltwi - .referenceWetnessIndex(dtm)) / tfact,
       wet = terra::clamp(exp((ltwi - twiwet) / .TWI_WET_JOIN), upper = 1))
}

# =========================================================================== #
# Wind shelter (grid model) -- derived from the horizon-angle array
# =========================================================================== #
#
# Wind shelter is derived from the upwind terrain horizon using the Ryan
# (1977) transformation. Horizon angles are available every 15 degrees; the
# shelter coefficient is interpolated between the two azimuths bracketing the
# actual wind direction so shelter varies smoothly as wind direction changes.
# Built directly from the horizon-angle array already computed for sky-view
# factor, rather than a separate terravars::wind_shelter() call, avoiding a
# redundant precompute; no separate neighbourhood-smoothing step is needed
# either, since the horizon array is already a full-neighbourhood search.

# Interpolate terrain shelter in the actual upwind direction. The function
# accepts either one direction per timestep for the whole domain or a fully
# spatial wind-direction field.
.windShelterInterp <- function(hora, wdir) {
  rows <- dim(hora)[1]; cols <- dim(hora)[2]
  raynCoef <- function(h) 1 - atan(17 * h) / 1.65

  if (is.null(dim(wdir))) {
    # ---- domain-uniform direction, one value per timestep ----
    wdirv <- as.vector(wdir) %% 360
    tsteps <- length(wdirv)
    idxLo <- (floor(wdirv / 15) %% 24) + 1L
    idxHi <- (idxLo %% 24) + 1L
    frac <- (wdirv %% 15) / 15
    shelterLo <- raynCoef(hora[, , idxLo, drop = FALSE])
    shelterHi <- raynCoef(hora[, , idxHi, drop = FALSE])
    frac3 <- rep(frac, each = rows * cols)
    result <- shelterLo * (1 - frac3) + shelterHi * frac3
    if (tsteps == 1) dim(result) <- c(rows, cols)
    result
  } else {
    # ---- genuinely per-fine-cell direction (array-mode weather) ----
    tsteps <- dim(wdir)[3]
    ncell <- rows * cols
    horaMat <- matrix(hora, nrow = ncell, ncol = 24)
    wdirMat <- matrix(wdir %% 360, nrow = ncell)
    idxLoMat <- (floor(wdirMat / 15) %% 24) + 1L
    idxHiMat <- (idxLoMat %% 24) + 1L
    fracMat <- (wdirMat %% 15) / 15
    cellIdx <- rep(seq_len(ncell), times = tsteps)
    loLin <- cellIdx + (as.vector(idxLoMat) - 1L) * ncell
    hiLin <- cellIdx + (as.vector(idxHiMat) - 1L) * ncell
    shelterLo <- raynCoef(horaMat[loLin])
    shelterHi <- raynCoef(horaMat[hiLin])
    result <- shelterLo * (1 - as.vector(fracMat)) + shelterHi * as.vector(fracMat)
    dim(result) <- c(rows, cols, tsteps)
    result
  }
}

# Resample a coarse environmental field to the fine terrain grid and mask it
# to the modelled land footprint. Plain arrays first inherit the spatial
# geometry of their coarse template; existing rasters already carry that
# geometry and retain their values directly.
.resampleToFine <- function(a, template, dtm) {
  r <- if (inherits(a, "SpatRaster")) {
    a
  } else {
    rr <- terra::rast(a)
    terra::ext(rr) <- terra::ext(template)
    terra::crs(rr) <- terra::crs(template)
    rr
  }
  if (!isTRUE(all.equal(terra::res(r), terra::res(dtm)))) {
    r <- terra::resample(r, dtm)
  }
  terra::mask(r, dtm)
}

# Assemble one variable from the list of coarse-cell point simulations into
# a coarse rows x cols x time array. Cells with no valid point simulation
# remain NA.
.coarseListToArray <- function(lst, varn, nT, dtmc) {
  rows <- dim(dtmc)[1]; cols <- dim(dtmc)[2]
  a <- array(NA_real_, dim = c(rows, cols, nT))
  for (k in seq_along(lst)) {
    if (is.null(lst[[k]])) next
    i <- (k - 1) %/% cols + 1
    j <- (k - 1) %% cols + 1
    a[i, j, ] <- lst[[k]][[varn]]
  }
  a
}

# Assemble a coarse-cell time series and interpolate it onto the fine terrain
# grid, retaining the fine land mask.
.coarseListToFineArray <- function(lst, varn, nT, dtmc, dtm) {
  a <- .coarseListToArray(lst, varn, nT, dtmc)
  terra::as.array(.resampleToFine(a, dtmc, dtm))
}

# Interpolate a time-invariant quantity defined once per coarse climate cell
# onto the fine terrain grid. This is used for reference properties such as
# measurement height or canopy geometry that do not vary through the run.
.coarseScalarToFineMatrix <- function(vals, dtmc, dtm) {
  r <- terra::rast(dtmc)
  terra::values(r) <- vals
  terra::as.matrix(.resampleToFine(r, dtmc, dtm), wide = TRUE)
}

# Saturation vapour pressure (kPa), using the Magnus relationship over water
# above freezing and over ice below freezing. Vectorised for spatial climate
# fields.
.satvap <- function(tc) {
  es <- 0.61078 * exp(17.27 * tc / (tc + 237.3))
  ei <- 0.61078 * exp(21.875 * tc / (tc + 265.5))
  ifelse(tc < 0, ei, es)
}

# Moist-air environmental lapse rate (degrees C/m) from temperature, vapour
# pressure and atmospheric pressure.
.lapserate <- function(tc, ea, pk) {
  rv <- 0.622 * ea / (pk - ea)
  9.8076 * (1 + (2501000 * rv) / (287 * (tc + 273.15))) /
    (1003.5 + (0.622 * 2501000^2 * rv) / (287 * (tc + 273.15)^2))
}

# Each cell's vapour network for the surface balances, from its own radiation,
# stomata, soil and resistances and the reference run's wet-canopy state
# (surfaceResistGridCpp): `rSurf`, `hSurf`, the foliage's share `hFol` and the
# soil humidity `hr`. `raw` is the grid radiation and wind state; climate fields
# may be hourly vectors or fine-grid arrays.
.cellSurfaceResistance <- function(raw, Rsw, Rdif, Ta, rh, pk, precip, year, soilm,
                                   campbellMats, hgt, pai, x, lref, ltra, svfa, gsmax,
                                   leafphys, Lfrac, wetShare, filmAvailable, paiRef) {
  nT <- dim(soilm)[3]
  veg <- c(list(hgt = hgt, pai = pai, x = x, leafr = lref, leaft = ltra, svfa = svfa,
                gsmax = gsmax, Lfrac = Lfrac), leafphys)
  surfaceResistGridCpp(
    rad = list(radCpar = raw$radCpar, radGpar = raw$radGpar, RabsCanopy = raw$radCsw + raw$radClw,
               zend = raw$zend, rHa = raw$rHa, rGz = raw$rGz, rGm = raw$rGm, emCanopy = raw$emCanopy),
    clim = list(Rsw = as.vector(Rsw), Rdif = as.vector(Rdif), Ta = as.vector(Ta), rh = as.vector(rh),
                pk = as.vector(pk), precip = as.vector(precip),
                Ca = rep(Cafromyear(stats::median(year)), nT)),
    soilm = soilm,
    soil = list(thetaS = campbellMats$thetaS, psie = campbellMats$psie, b = campbellMats$b),
    veg = veg,
    ref = list(wetShare = as.vector(wetShare), filmAvailable = as.vector(filmAvailable),
               pai = as.vector(paiRef)))
}

# Bare cells, as an array matching the ground pass's cell-hour layout.
.bareCells <- function(gh, hgt) {
  array(rep(as.vector(hgt <= 0), dim(gh$Ts_est)[3]), dim = dim(gh$Ts_est))
}

# Absorbed radiation of the bulk surface, kept only in the cells whose bulk
# temperature the ground pass has not formed and which are vegetated, so that
# canopyTempCpp, which skips missing values, estimates only those.
.oneSurfaceRabs <- function(RabsCanopy, gh, hgt) {
  RabsCanopy[is.finite(gh$Tcanopy_est) | is.finite(gh$Tcanopy_one) | .bareCells(gh, hgt)] <- NA_real_
  RabsCanopy
}

# The surface state every grid worker reports, from the ground pass and the
# one-surface bulk estimate. The ground pass gives the bulk temperature where it
# solves the bulk surface together with the ground, and where it forms the
# one-surface estimate itself; the separate estimate serves any other vegetated
# cell. On bare ground the bulk surface is the ground, and both are the ground
# pass's temperature.
# Below 10 mm the reported ground temperature is blended toward the bulk
# surface's with the point model's matching weight, 0 at the soil's heat
# roughness height (0.8 mm) and 1 at 10 mm, end points taken exactly; this is
# the one-shot form of the point model's matching constraint. The humidity
# factor handed to the humidity profiles is the one that reproduces the vapour
# network's flux at the reported temperatures, with the soil's source at the
# ground's own temperature: hFol + (1 - hFol) hr es(Ts)/es(Tc). Bare cells keep
# the soil's own humidity factor.
.coupledSurfaceState <- function(gh, TcOne, sr, hgt) {
  Tc <- ifelse(is.finite(gh$Tcanopy_est), gh$Tcanopy_est,
               ifelse(is.finite(gh$Tcanopy_one), gh$Tcanopy_one,
                      ifelse(.bareCells(gh, hgt), gh$Ts_est, TcOne)))
  Ts <- gh$Ts_est
  # Bare ground (zero height) is the ground's own balance, weight one.
  wM <- ifelse(hgt > 0, pmin(1, pmax(0, (hgt - 0.0008) / (0.010 - 0.0008))), 1)
  wArr <- array(rep(as.vector(wM), dim(Ts)[3]), dim = dim(Ts))
  Ts <- ifelse(wArr >= 1, Ts, ifelse(wArr <= 0, Tc, wArr * Ts + (1 - wArr) * Tc))
  veg <- is.finite(gh$Tcanopy_est) & is.finite(sr$hFol)
  hSurf <- ifelse(veg, sr$hFol + (1 - sr$hFol) * sr$hr * .satvap(Ts) / .satvap(Tc), sr$hSurf)
  list(Ts_est = Ts, Tcanopy_est = Tc, hSurf = hSurf)
}

# Stomatal conductance ceiling per cell from vegp$gsmax (first layer), or no
# ceiling where it is not supplied.
.gsmaxMatrix <- function(vegp, rows, cols) {
  g <- vegp$gsmax
  if (is.null(g)) return(matrix(NA_real_, rows, cols))
  if (inherits(g, "PackedSpatRaster")) g <- terra::rast(g)
  if (inherits(g, "SpatRaster")) return(terra::as.matrix(g[[1]], wide = TRUE))
  matrix(g, rows, cols)
}

# Mean over time of the elevation correction applied to each fine cell's air
# temperature, so that an annual mean temperature follows the same correction.
.meanAltitudeOffset <- function(tcCorrected, tcUncorrected) {
  apply(tcCorrected - tcUncorrected, c(1, 2), mean)
}

# Correct coarse climate forcing for elevation differences between the
# climate grid and the fine terrain. Pressure is corrected via sea-level
# pressure using the barometric relationship. With altcorrect = 0, pressure is
# only resampled and temperature is unchanged; altcorrect = 1 applies that
# pressure correction plus a fixed 5 degrees C/km temperature lapse rate; any
# other nonzero value uses the pressure correction plus a moist-air lapse rate
# derived from local temperature, humidity and pressure.
.altitudeCorrectPkTc <- function(pkCoarse, tcFine, relhumFine, dtmc, dtm, altcorrect) {
  if (altcorrect == 0) {
    return(list(pk = terra::as.array(.resampleToFine(pkCoarse, dtmc, dtm)), tc = tcFine))
  }

  dtmcZeroed <- dtmc
  dtmcZeroed[is.na(dtmcZeroed)] <- 0
  # Elevation is constant through time. Flatten the spatial elevation fields
  # so R recycles the same cell-specific correction across every timestep of
  # the three-dimensional climate arrays.
  dtmcMat <- as.vector(terra::as.matrix(dtmcZeroed, wide = TRUE))
  dtmMat <- as.vector(terra::as.matrix(dtm, wide = TRUE))

  # Reduce coarse pressure to sea level using the coarse grid's own
  # elevation (a smooth field, safe to resample directly), then re-inflate
  # using the fine grid's own real elevation -- recovers the orographic
  # pressure signal a direct coarse -> fine resample of raw pressure would
  # otherwise discard.
  pslCoarse <- pkCoarse / ((293 - 0.0065 * dtmcMat) / 293)^5.26
  pslFine <- terra::as.array(.resampleToFine(pslCoarse, dtmcZeroed, dtm))
  pkFine <- pslFine * ((293 - 0.0065 * dtmMat) / 293)^5.26

  # Fine-scale temperature correction depends on the elevation difference
  # between the actual terrain and the elevation represented by the resampled
  # coarse climate grid.
  dtmcOnDtm <- as.vector(terra::as.matrix(terra::resample(dtmcZeroed, dtm), wide = TRUE))
  elevd <- dtmcOnDtm - dtmMat

  if (altcorrect == 1) {
    tcdif <- elevd * (5 / 1000) # fixed environmental lapse rate, 5 degC/km
  } else {
    esFine <- .satvap(tcFine)
    eaFine <- esFine * relhumFine / 100
    tcdif <- .lapserate(tcFine, eaFine, pkFine) * elevd
  }
  list(pk = pkFine, tc = tcFine + tcdif)
}

# Interpolate wind from coarse to fine resolution in vector form. Speed and
# direction are converted to orthogonal components before interpolation so
# directions close to 0/360 degrees combine correctly, then reconstructed as
# fine-resolution wind speed and direction.
.coarseWindToFine <- function(lst, nT, dtmc, dtm) {
  u2 <- .coarseListToArray(lst, "windspeed", nT, dtmc)
  wd <- .coarseListToArray(lst, "winddir", nT, dtmc)
  wu <- u2 * cos(wd * pi / 180)
  wv <- u2 * sin(wd * pi / 180)

  wuFine <- terra::as.array(.resampleToFine(wu, dtmc, dtm))
  wvFine <- terra::as.array(.resampleToFine(wv, dtmc, dtm))
  winddirFine <- (atan2(wvFine, wuFine) * 180 / pi) %% 360
  list(windspeed = sqrt(wuFine^2 + wvFine^2), winddir = winddirFine)
}

# =========================================================================== #
# Tiled execution and raster mosaics
# =========================================================================== #

# Choose a conservative square tile size from the requested memory budget and
# number of timesteps. The estimate allows for the several simultaneous
# space-by-time arrays required by the grid calculation and keeps a safety margin
# rather than attempting to consume the full stated memory budget: only half
# the caller's stated budget is actually used, and the candidate tile size is
# the largest one that still fits under that half-budget, so the choice is
# always rounded down to a smaller, safer tile rather than up.
.autoTileSize <- function(rows, cols, tsteps, maxmemGB = 4) {
  arraysAtPeak <- 30
  bytesPerElement <- 8
  budgetBytes <- maxmemGB * 1e9 * 0.5
  maxCellsTimesteps <- budgetBytes / (arraysAtPeak * bytesPerElement)
  maxCells <- maxCellsTimesteps / tsteps
  osize <- sqrt(maxCells)
  sizeCandidates <- c(10, 20, 50, 100, 200, 500, 1000, 2000)
  fits <- sizeCandidates[sizeCandidates <= osize]
  tilesize <- if (length(fits) == 0) sizeCandidates[1] else max(fits)
  min(tilesize, max(rows, cols)) # no point exceeding the domain's own extent
}

# Spatial extent of one square tile, including an optional overlap measured
# in pixels and clipped to the domain boundary.
.tileExtent <- function(dtm, rw, cl, tilesize, toverlap) {
  e <- terra::ext(dtm)
  rs <- terra::res(dtm)
  xmn <- e$xmin + (cl - 1) * tilesize * rs[1] - toverlap * rs[1]
  xmx <- e$xmin + cl * tilesize * rs[1] + toverlap * rs[1]
  ymn <- e$ymax - rw * tilesize * rs[2] - toverlap * rs[2]
  ymx <- e$ymax - (rw - 1) * tilesize * rs[2] + toverlap * rs[2]
  xmn <- max(xmn, e$xmin); xmx <- min(xmx, e$xmax)
  ymn <- max(ymn, e$ymin); ymx <- min(ymx, e$ymax)
  terra::ext(c(xmn, xmx, ymn, ymx))
}

# Crop spatial vegetation or soil fields to one tile while leaving scalar
# parameters unchanged. Seasonal raster stacks retain all of their layers.
.cropVegSoilList <- function(lst, e) {
  out <- lapply(lst, function(r) {
    if (is.null(r)) return(r)
    if (inherits(r, "PackedSpatRaster")) r <- terra::rast(r)
    if (!inherits(r, "SpatRaster")) return(r)
    terra::crop(r, e)
  })
  names(out) <- names(lst)
  out
}

# Convert the user-facing parallel setting to a physical worker count.
# Sequential execution is the default; "auto" leaves one physical core free,
# "max" uses all physical cores, and an explicit count is capped at the number
# available to avoid spawning redundant memory-heavy R processes.
.resolveCores <- function(cores = "off") {
  physicalCores <- function() {
    p <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
    if (is.na(p) || p < 1) p <- 1L
    p
  }
  if (is.character(cores)) {
    cores <- match.arg(cores, c("off", "auto", "max"))
    if (cores == "off") return(1L)
    physical <- physicalCores()
    if (cores == "auto") return(max(1L, physical - 1L))
    # "max"
    if (physical > 1) {
      warning("cores = \"max\" requested: using all ", physical, " physical ",
               "cores. Your machine may be slow or unresponsive for other ",
               "work while this runs.", call. = FALSE)
    }
    return(physical)
  }
  coresN <- suppressWarnings(as.integer(cores))
  if (is.na(coresN) || coresN < 1) {
    stop(".resolveCores(): `cores` must be \"off\", \"auto\", \"max\", or a ",
         "positive whole number, not: ", format(cores))
  }
  physical <- physicalCores()
  if (coresN > physical) {
    message("Requested cores = ", coresN, " exceeds this machine's ", physical,
            " physical core(s) -- using ", physical, " instead.")
    coresN <- physical
  }
  coresN
}

# Convert spatial rasters to a serialisable form before they cross R process
# boundaries during parallel execution, preserving the surrounding list structure.
.wrapForWorker <- function(x) {
  if (inherits(x, "SpatRaster")) return(terra::wrap(x))
  if (is.list(x)) return(lapply(x, .wrapForWorker))
  x
}

# Restore serialised rasters to live SpatRaster objects inside the receiving
# process after parallel transfer.
.unwrapFromWorker <- function(x) {
  if (inherits(x, "PackedSpatRaster")) return(terra::rast(x))
  if (is.list(x)) return(lapply(x, .unwrapFromWorker))
  x
}

# Directions (degrees from north) in which the terrain horizon is searched.
.HORIZON_AZIMUTHS <- seq(0, 345, by = 15)

# Terrain horizon angle (degrees) in 24 directions, one layer each: every cell
# out to 20 cells, then more widely spaced. tick is called after each direction.
.horizonRaster <- function(dtm, tick = NULL) {
  layers <- lapply(.HORIZON_AZIMUTHS, function(az) {
    r <- terravars::horizon(dtm, azi = az, near = 20)
    if (!is.null(tick)) tick()
    r
  })
  terra::rast(layers)
}

# Tangent of the horizon angle (rows x cols x 24) and sky view factor for a grid
# run. A sky view factor not supplied is formed from the horizon.
.horizonAndSkyview <- function(dtm, hor, svf) {
  if (is.logical(hor)) {
    horR <- .horizonRaster(dtm)
    hor <- tan(terra::as.array(horR) * pi / 180)
  } else if (is.logical(svf)) {
    horR <- terra::rast(atan(hor) * 180 / pi, extent = terra::ext(dtm), crs = terra::crs(dtm))
  }
  if (is.logical(svf)) svf <- terravars::skyview(dtm, hor = horR)
  list(hora = hor, svf = svf)
}

# Text progress bar in terra's style, over n steps. finish() wipes it, leaving
# the console at the start of an empty line.
.progressBar <- function(n, show = TRUE) {
  if (!isTRUE(show) || n < 1) {
    return(list(tick = function() invisible(NULL), finish = function() invisible(NULL)))
  }
  scale <- "|---------|---------|---------|---------|"
  width <- nchar(scale)
  done <- 0L
  drawn <- 0L
  cat("\r", scale, "\r", sep = "")
  utils::flush.console()
  list(
    tick = function() {
      done <<- min(done + 1L, n)
      upto <- round(done * width / n)
      if (upto > drawn) {
        cat(strrep("=", upto - drawn))
        utils::flush.console()
        drawn <<- upto
      }
      invisible(NULL)
    },
    finish = function() {
      cat("\r", strrep(" ", width), "\r", sep = "")
      utils::flush.console()
      invisible(NULL)
    }
  )
}

# Compute neighbourhood-dependent terrain properties on the full landscape
# before tiling. Horizon, sky view, slope/aspect and topographic wetness must retain
# their domain context at tile boundaries. show = TRUE draws a progress bar.
.domainWideTerrain <- function(dtm, hor, twi, svf, slr, apr, show = FALSE) {
  bar <- .progressBar(if (is.logical(hor)) length(.HORIZON_AZIMUTHS) + 3L else 3L, show)
  on.exit(bar$finish(), add = TRUE)
  horR <- if (is.logical(hor)) {
    .horizonRaster(dtm, tick = bar$tick) # degrees, NOT yet tan-transformed
  } else .unpackRaster(hor)
  svfR <- if (is.logical(svf)) {
    terravars::skyview(dtm, hor = horR)
  } else .unpackRaster(svf)
  bar$tick()
  if (is.logical(slr) || is.logical(apr)) {
    fill_ext <- terravars:::.fillna_dtm(dtm)
  }
  slrR <- if (is.logical(slr)) {
    terra::crop(terra::terrain(fill_ext$extended, "slope", unit = "degrees"), dtm)
  } else .unpackRaster(slr)
  aprR <- if (is.logical(apr)) {
    terra::crop(terra::terrain(fill_ext$extended, "aspect", unit = "degrees"), dtm)
  } else .unpackRaster(apr)
  bar$tick()
  twiR <- if (is.logical(twi)) .wetnessIndex(dtm) else .unpackRaster(twi)
  bar$tick()

  list(svfR = svfR, horR = horR, slrR = slrR, aprR = aprR, twiR = twiR)
}

# Evaluate every argument of the calling function now. R defers evaluating an
# argument until it is first used, keeping only a reference to the caller's
# environment; a parallel worker receives that reference but not the caller's
# environment, so an argument passed as a variable would not be found there.
.forceArgs <- function(env = parent.frame()) {
  invisible(mget(names(formals(sys.function(sys.parent()))), envir = env))
}

# Execute a landscape calculation either as one domain or as independent,
# non-overlapping tiles. Neighbourhood-dependent terrain inputs are computed
# before splitting, then each tile is evaluated and merged back onto the
# original grid.
.tiledDispatch <- function(dtm, vegp, soilc, hor, twi, svf, slr, apr,
                            nT, cores, maxmemGB, tilesize, workerFn, extraArgs = list(),
                            silent = FALSE) {
  dtm <- .unpackRaster(dtm)
  coresN <- .resolveCores(cores)
  if (coresN <= 1) {
    return(do.call(workerFn, c(list(dtm = dtm, vegp = vegp, soilc = soilc,
                                     hor = hor, twi = twi, svf = svf, slr = slr, apr = apr),
                               extraArgs)))
  }

  dwt <- .domainWideTerrain(dtm, hor, twi, svf, slr, apr)
  rows <- dim(dtm)[1]; cols <- dim(dtm)[2]
  ts <- if (!is.na(tilesize)) tilesize else .autoTileSize(rows, cols, nT, maxmemGB = maxmemGB)
  rws <- ceiling(rows / ts); cls <- ceiling(cols / ts)

  jobs <- list()
  for (rw in seq_len(rws)) {
    for (cl in seq_len(cls)) {
      e <- .tileExtent(dtm, rw, cl, ts, 0)
      dtmi <- terra::crop(dtm, e)
      if (all(is.na(terra::values(dtmi)))) next
      jobs[[length(jobs) + 1]] <- list(
        dtmi = terra::wrap(dtmi),
        vegpi = .wrapForWorker(.cropVegSoilList(vegp, e)),
        soilci = .wrapForWorker(.cropVegSoilList(soilc, e)),
        svfi = terra::wrap(terra::crop(dwt$svfR, e)),
        hori = tan(terra::as.array(terra::crop(dwt$horR, e)) * pi / 180),
        slri = terra::wrap(terra::crop(dwt$slrR, e)),
        apri = terra::wrap(terra::crop(dwt$aprR, e)),
        twii = terra::wrap(terra::crop(dwt$twiR, e))
      )
    }
  }
  if (length(jobs) == 0) {
    stop(".tiledDispatch(): every tile of `dtm` is entirely NA -- nothing to compute")
  }

  runOneJob <- function(job, .extraArgs, .workerFn, p = NULL) {
    res <- do.call(.workerFn, c(list(
      dtm = terra::rast(job$dtmi),
      vegp = .unwrapFromWorker(job$vegpi),
      soilc = .unwrapFromWorker(job$soilci),
      hor = job$hori,
      twi = terra::rast(job$twii),
      svf = terra::rast(job$svfi),
      slr = terra::rast(job$slri),
      apr = terra::rast(job$apri)), .extraArgs))
    if (!is.null(p)) p()
    .wrapForWorker(res)
  }

  oldplan <- future::plan()
  on.exit(future::plan(oldplan), add = TRUE)
  future::plan(future::multisession, workers = coresN)

  haveProgressr <- requireNamespace("progressr", quietly = TRUE)
  if (haveProgressr) {
    results <- progressr::with_progress({
      p <- progressr::progressor(steps = length(jobs))
      future.apply::future_lapply(jobs, runOneJob, .extraArgs = extraArgs,
                                   .workerFn = workerFn, p = p, future.seed = TRUE)
    })
  } else {
    if (!isTRUE(silent)) {
      message(".tiledDispatch(): running ", length(jobs), " tile(s) on ", coresN,
              " cores. Install the \"progressr\" package for a live progress bar ",
              "during parallel runs (none shown otherwise).")
    }
    results <- future.apply::future_lapply(jobs, runOneJob, .extraArgs = extraArgs,
                                            .workerFn = workerFn, p = NULL, future.seed = TRUE)
  }

  # A field with no value at the requested height is a single NA in every tile.
  fieldNames <- names(results[[1]])
  merged <- lapply(fieldNames, function(f) {
    if (!inherits(results[[1]][[f]], "PackedSpatRaster")) return(NA)
    layers <- lapply(results, function(r) terra::unwrap(r[[f]]))
    if (length(layers) == 1) return(layers[[1]])
    do.call(terra::merge, layers)
  })
  names(merged) <- fieldNames
  merged
}

# Blend two neighbouring raster tiles across their shared overlap using a
# linear distance weighting. Each tile dominates the edge nearest its own
# non-overlap region, giving a smooth transition across the seam while leaving
# all non-overlapping cells unchanged.
.blendmosaic <- function(r1, r2, direction = c("right", "bottom")) {
  direction <- match.arg(direction)
  ov <- terra::intersect(terra::ext(r1), terra::ext(r2))
  if (is.null(ov) || terra::is.empty(ov)) return(terra::merge(r1, r2))

  r1o <- terra::crop(r1, ov)
  r2o <- terra::crop(r2, ov)
  nx <- dim(r1o)[2]; ny <- dim(r1o)[1]
  # weight for r1 -- ramps 1 -> 0 across the overlap in the direction r2
  # sits away from r1, so r1 dominates the r1-facing edge of the overlap
  # and r2 dominates the r2-facing edge, smoothly in between.
  w1 <- if (direction == "right") {
    matrix(rep(seq(1, 0, length.out = nx), each = ny), nrow = ny, ncol = nx)
  } else {
    matrix(rep(seq(1, 0, length.out = ny), times = nx), nrow = ny, ncol = nx)
  }
  w1r <- terra::rast(r1o[[1]]); terra::values(w1r) <- w1
  blended <- r1o * w1r + r2o * (1 - w1r)

  # Remove the original values only within the overlap before merging, so the
  # blended band is the sole source of values there and the rest of each tile
  # remains unchanged.
  ovPoly <- terra::as.polygons(ov)
  terra::crs(ovPoly) <- terra::crs(r1)
  r1rem <- terra::mask(r1, ovPoly, inverse = TRUE)
  r2rem <- terra::mask(r2, ovPoly, inverse = TRUE)
  terra::merge(terra::merge(r1rem, blended), r2rem)
}
