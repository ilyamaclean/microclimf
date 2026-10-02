# gridmodel.R
# Spatially resolves the point-model surface energy balance across a landscape.
# Terrain and vegetation modify radiation, aerodynamic exposure, soil wetness and
# heat exchange at each cell, while the point-model solution supplies the
# time-varying atmospheric and surface state that anchors those spatial
# adjustments. The requested height may lie below ground, at the surface, within
# vegetation or above the canopy; the physically meaningful outputs differ
# accordingly.
#
# Terrain sky exposure, horizon shading and wind shelter describe topographic
# obstruction only. Canopy interception and attenuation are handled separately
# by the vegetation radiation and turbulence calculations, so vegetation is not
# folded into those terrain factors a second time.

#' Model microclimate across a landscape
#'
#' Uses a solved point-model run as the time-varying physical reference for
#' spatial estimates of radiation, wind, temperature, humidity and soil moisture
#' across a terrain grid. Local differences arise from elevation, slope and
#' aspect, horizon shading, sky exposure, topographic wetness and vegetation
#' structure rather than from independently solving the full coupled surface
#' energy balance at every fine-resolution cell. In contrast,
#' \code{\link{runpointmodelasgrid}} solves a full point model independently
#' in every cell.
#'
#' The driving climate can be represented by one point-model run for the whole
#' domain or by a coarser grid of point-model runs. Vegetation can likewise be
#' fixed through time or supplied as seasonal layers. These choices change the
#' spatial and temporal information available to the calculation, but not the
#' meaning of the requested outputs.
#'
#' @details
#' Soil moisture is distributed from the point model's by the topographic
#' wetness index. The point model stands for a flat cell with no contributing
#' area: a cell with that cell's index takes the point model's soil moisture, a
#' cell with a higher index is wetter and one with a lower index drier, by an
#' amount set by \code{tfact}. Ground whose index is at or above \code{twiwet},
#' such as a stream channel, is saturated at every time step.
#'
#' For vegetation supplied as several layers, each field is treated as a
#' sequence of seasonal states. With \code{splineveg = FALSE}, days are assigned
#' to discrete states; with \code{splineveg = TRUE}, values are interpolated
#' smoothly between them. Plant functional type is held fixed for a cell during
#' the run, while canopy structure may vary seasonally.
#'
#' The interpretation of \code{reqhgt} depends on its position relative to the
#' ground and the local canopy. Below ground the model returns soil state at the
#' requested depth; at zero height it returns surface quantities; above ground,
#' temperature, humidity, wind and radiation are evaluated either within or
#' above the local canopy. Leaf temperature is defined only where the requested
#' height lies within vegetation.
#'
#' @param pointmodel Point-model output from \code{\link{runpointmodel}} or
#'   \code{\link{subsetpointmodel}}, in single-location or gridded-climate mode.
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
#' @param tfact coefficient determining sensitivity of soil moisture to
#'   variation in topographic wetness.
#' @param twiwet topographic wetness index at and above which the ground is
#'   permanently wet: soil there is saturated at every time step. Set to
#'   \code{Inf} for no permanently wet ground.
#' @param altcorrect a single numeric value indicating whether to apply an
#'   elevational lapse rate correction to temperatures (0 = no correction, 1 =
#'   fixed lapse rate correction, 2 = humidity-dependent variable lapse rate
#'   correction). Used only when climate data are provided as a multi-layer
#'   raster.
#' @param hor an optional array of horizon angles (degrees) in 24 directions.
#'   Calculated automatically from \code{dtm} if not supplied, but can also be
#'   provided as separate inputs to avoid edge effects.
#' @param twi optional raster object of topographic wetness index values.
#'   Calculated automatically from \code{dtm} if not supplied, but can also be
#'   provided as separate inputs to avoid edge effects. A supplied index must
#'   be the natural logarithm of upslope area per metre of contour over the
#'   tangent of slope, as returned by \code{terravars::twi()}.
#' @param svf optional raster object of sky view factor. Calculated
#'   automatically from \code{dtm} if not supplied, but can also be provided as
#'   separate inputs to avoid edge effects.
#' @param slr,apr slope and aspect in degrees. Calculated automatically from
#'   \code{dtm} if not supplied, but can also be provided as separate inputs to
#'   avoid edge effects.
#' @param out Optional named vector of logicals indicating which variables to
#'   return, e.g. \code{c(Tz = TRUE, soilm = TRUE)}.
#' @param splineveg Logical indicating whether to interpolate seasonal
#'   vegetation layers smoothly rather than assigning discrete layers to days.
#' @param saveout Logical; when set to \code{TRUE} the results are saved to
#'   \code{filename}.
#' @param filename File path used when \code{saveout = TRUE}.
#' @param cores Controls parallel processing across tiles: \code{off} (default)
#'   runs sequentially; \code{auto} uses one fewer than the number of
#'   available cores; \code{max} uses all available cores; alternatively, a
#'   positive integer can be supplied to specify the number of cores to use.
#' @param maxmemGB Memory budget (GB) used to choose tile size when
#'   parallelised. Used only when \code{cores} is set such that the model runs in
#'   parallel tiles; ignored with the default \code{cores = "off"}.
#' @param tilesize Optional tile width/height in pixels. Used only when
#'   \code{cores} is set to run the model in parallel mode.
#'
#' @return A named list of \code{SpatRaster}s, with one layer per timestep, for
#'   the requested subset of: \code{Tz} (air temperature above ground,
#'   ground-surface temperature at \code{reqhgt = 0}, soil temperature below
#'   ground, deg C), \code{tleaf} (leaf temperature, deg C), \code{soilm}
#'   (volumetric soil moisture), \code{relhum} (relative humidity, percent),
#'   \code{windspeed} (m/s), \code{Rdirdown}, \code{Rdifdown} and \code{Rswup}
#'   (shortwave radiation, W/m2), and \code{Rlwdown} and \code{Rlwup} (longwave
#'   radiation, W/m2). A quantity with no physical meaning at the requested
#'   height/depth is returned as a single \code{NA} rather than omitted, so
#'   the field set is fixed regardless of \code{reqhgt}. Which fields are
#'   meaningful depends on the regime: below ground, only
#'   \code{Tz}/\code{soilm} (soil temperature and moisture at that depth); at
#'   the ground surface, \code{Tz}/\code{soilm} (surface temperature and soil
#'   moisture) and the radiation fields; above ground,
#'   \code{Tz}/\code{relhum}/\code{windspeed}, the radiation fields and
#'   \code{tleaf}, but not \code{soilm}. \code{tleaf} holds values only in
#'   cells where \code{reqhgt} lies within the canopy.
#' @export
rungridmodel <- function(pointmodel, reqhgt = 0, dtm, vegp, soilc, tfact = 1.7, twiwet = 10,
                          altcorrect = 0, hor = NA, twi = NA, svf = NA,
                          slr = NA, apr = NA, out = NULL, splineveg = FALSE,
                          saveout = FALSE, filename = NULL,
                          cores = "off", maxmemGB = 4, tilesize = NA) {
  .forceArgs()
  isArray <- !is.null(pointmodel$points) && !is.null(pointmodel$dtmc)
  isSingle <- !is.null(pointmodel$weather)
  if (!isArray && !isSingle) {
    stop("`pointmodel` must be the result of runpointmodel() or subsetpointmodel().",
         call. = FALSE)
  }
  if (isTRUE(saveout) && is.null(filename)) {
    stop("rungridmodel(): `filename` is required when saveout = TRUE -- the ",
         "full path (directory + file name) to save the wrapped result to, ",
         "e.g. filename = \"C:/path/to/output.rds\"", call. = FALSE)
  }

  vegp <- .resolveVegp(vegp, dtm)
  vegFieldNames <- c("hgt", "pai", "x", "leafr", "leaft", "clump", "Lfrac")
  isTimeVariant <- any(vapply(vegFieldNames, function(f) .nLayersOf(vegp[[f]]), integer(1)) > 1)

  # Four combinations are possible: shared or spatially varying climate forcing,
  # crossed with fixed or seasonally varying vegetation. The appropriate
  # pathway preserves that distinction while applying the same landscape-scale
  # radiation, turbulence, soil and temperature calculations.
  dispatchOne <- function(dtm, vegp, soilc, hor, twi, svf, slr, apr) {
    if (!isArray && !isTimeVariant) {
      .rungridmodel1(pointmodel, vegp, soilc, dtm, reqhgt = reqhgt, tfact = tfact, twiwet = twiwet,
                     hor = hor, twi = twi, svf = svf, slr = slr, apr = apr, out = out)
    } else if (isArray && !isTimeVariant) {
      .rungridmodel2(pointmodel, vegp, soilc, dtm, reqhgt = reqhgt, tfact = tfact, twiwet = twiwet,
                     altcorrect = altcorrect, hor = hor, twi = twi, svf = svf,
                     slr = slr, apr = apr, out = out)
    } else if (!isArray && isTimeVariant) {
      .rungridmodel3(pointmodel, vegp, soilc, dtm, reqhgt = reqhgt, tfact = tfact, twiwet = twiwet,
                     hor = hor, twi = twi, svf = svf, slr = slr, apr = apr, out = out,
                     splineveg = splineveg)
    } else {
      .rungridmodel4(pointmodel, vegp, soilc, dtm, reqhgt = reqhgt, tfact = tfact, twiwet = twiwet,
                     altcorrect = altcorrect, hor = hor, twi = twi, svf = svf,
                     slr = slr, apr = apr, out = out, splineveg = splineveg)
    }
  }

  nT <- if (isSingle) {
    nrow(pointmodel$weather)
  } else {
    validCells <- which(vapply(pointmodel$points, is.list, logical(1)))
    if (length(validCells) == 0) stop(.noUsableRunMessage(), call. = FALSE)
    nrow(pointmodel$points[[validCells[1]]]$weather)
  }

  result <- .tiledDispatch(dtm, vegp, soilc, hor, twi, svf, slr, apr,
                            nT = nT, cores = cores, maxmemGB = maxmemGB,
                            tilesize = tilesize, workerFn = dispatchOne)

  # Keep the returned rasters directly usable in the current R session; when
  # saving, wrap a separate copy so terra's external pointers survive reload.
  if (isTRUE(saveout)) {
    saveRDS(.wrapForWorker(result), filename)
  }
  result
}

#' Model very large landscapes in tiles
#'
#' Runs \code{\link{rungridmodel}} over domains that are too large to process
#' comfortably as one raster. Terrain quantities whose values depend on the
#' surrounding landscape are calculated from the whole domain before it is
#' divided, so the physical terrain context is not reset independently for each
#' tile. Tile outputs are written to disk and can be recombined with
#' \code{\link{mosaicblend}}.
#'
#' Overlap can be retained around tile edges when a smoothly blended mosaic is
#' required. This affects how saved tiles are joined; it does not change the
#' underlying interpretation of the microclimate variables.
#'
#' @param pointmodel Point-model output used to drive the grid calculation.
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
#' @param tfact coefficient determining sensitivity of soil moisture to
#'   variation in topographic wetness.
#' @param twiwet topographic wetness index at and above which the ground is
#'   permanently wet: soil there is saturated at every time step. Set to
#'   \code{Inf} for no permanently wet ground.
#' @param altcorrect a single numeric value indicating whether to apply an
#'   elevational lapse rate correction to temperatures (0 = no correction, 1 =
#'   fixed lapse rate correction, 2 = humidity-dependent variable lapse rate
#'   correction). Used only when climate data are provided as a multi-layer
#'   raster.
#' @param out Optional named vector of logicals indicating which variables to
#'   return, e.g. \code{c(Tz = TRUE, soilm = TRUE)}.
#' @param splineveg Logical indicating whether to interpolate seasonal
#'   vegetation layers smoothly rather than assigning discrete layers to days.
#' @param hor an optional array of horizon angles (degrees) in 24 directions.
#'   Calculated automatically from \code{dtm} if not supplied, but can also be
#'   provided as separate inputs to avoid edge effects.
#' @param twi optional raster object of topographic wetness index values.
#'   Calculated automatically from \code{dtm} if not supplied, but can also be
#'   provided as separate inputs to avoid edge effects. A supplied index must
#'   be the natural logarithm of upslope area per metre of contour over the
#'   tangent of slope, as returned by \code{terravars::twi()}.
#' @param svf optional raster object of sky view factor. Calculated
#'   automatically from \code{dtm} if not supplied, but can also be provided as
#'   separate inputs to avoid edge effects.
#' @param slr,apr slope and aspect in degrees. Calculated automatically from
#'   \code{dtm} if not supplied, but can also be provided as separate inputs to
#'   avoid edge effects.
#' @param pathout Directory in which tile outputs are written.
#' @param tilesize Optional tile width/height in pixels.
#' @param toverlap Number of overlapping pixels retained around adjacent tiles.
#' @param maxmemGB Memory budget (GB) used to choose tile size.
#' @param writeformat Output format: \code{"rds"}, \code{"nc"}, or both.
#' @param scalefactor Optional multiplier used before integer storage of outputs.
#' @param cores Controls parallel processing across tiles: \code{off} (default)
#'   runs sequentially; \code{auto} uses one fewer than the number of
#'   available cores; \code{max} uses all available cores; alternatively, a
#'   positive integer can be supplied to specify the number of cores to use.
#' @param silent Logical; suppress progress reporting when \code{TRUE}.
#'
#' @return Invisibly, the paths of tile files written (empty if every tile
#'   was entirely non-land and so skipped -- not an error).
#' @export
rungridmodelbig <- function(pointmodel, reqhgt = 0, dtm, vegp, soilc, tfact = 1.7, twiwet = 10,
                             altcorrect = 0, out = NULL, splineveg = FALSE,
                             hor = NA, twi = NA, svf = NA, slr = NA, apr = NA,
                             pathout = getwd(), tilesize = NA, toverlap = 0,
                             maxmemGB = 4, writeformat = "rds", scalefactor = NULL,
                             cores = "off", silent = FALSE) {
  .forceArgs()
  writeformat <- match.arg(writeformat, c("rds", "nc"), several.ok = TRUE)
  if ("nc" %in% writeformat && !requireNamespace("ncdf4", quietly = TRUE)) {
    stop("rungridmodelbig(): writeformat includes \"nc\" but the ncdf4 package ",
         "is not installed -- install.packages(\"ncdf4\") or drop \"nc\" from ",
         "writeformat.", call. = FALSE)
  }

  isArray <- !is.null(pointmodel$points) && !is.null(pointmodel$dtmc)
  isSingle <- !is.null(pointmodel$weather)
  if (!isArray && !isSingle) {
    stop("`pointmodel` must be the result of runpointmodel() or subsetpointmodel().",
         call. = FALSE)
  }

  dtm <- .unpackRaster(dtm)
  vegp <- .resolveVegp(lapply(vegp, .unpackRaster), dtm)
  soilc <- lapply(soilc, .unpackRaster)

  # The number of timesteps contributes to the memory cost of each spatial tile.
  if (isSingle) {
    nT <- nrow(pointmodel$weather)
  } else {
    validCells <- which(vapply(pointmodel$points, is.list, logical(1)))
    if (length(validCells) == 0) stop(.noUsableRunMessage(), call. = FALSE)
    nT <- nrow(pointmodel$points[[validCells[1]]]$weather)
  }

  rows <- dim(dtm)[1]; cols <- dim(dtm)[2]

  if (is.na(tilesize)) {
    tilesize <- .autoTileSize(rows, cols, nT, maxmemGB = maxmemGB)
    if (!isTRUE(silent)) {
      message("rungridmodelbig(): auto-selected tilesize = ", tilesize,
              " pixels (maxmemGB = ", maxmemGB, ") -- pass `tilesize` directly to override.")
    }
  }

  outdir <- file.path(pathout, "microut")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Compute neighbourhood-dependent terrain properties on the complete landscape
  # before cropping tiles, so horizon, sky exposure and wetness retain the same
  # physical neighbourhood they would have in an untiled run.
  show <- !isTRUE(silent)
  if (show) message("Terrain: horizon, sky view, slope and wetness of the whole area")
  dwt <- .domainWideTerrain(dtm, hor, twi, svf, slr, apr, show = show)

  rws <- ceiling(rows / tilesize)
  cls <- ceiling(cols / tilesize)

  # The tiles that hold land. Overlap is retained where requested so
  # neighbouring predictions can subsequently be blended.
  tiles <- list()
  for (rw in seq_len(rws)) {
    for (cl in seq_len(cls)) {
      e <- .tileExtent(dtm, rw, cl, tilesize, toverlap)
      if (all(is.na(terra::values(terra::crop(dtm, e))))) next
      tiles[[length(tiles) + 1]] <- list(rw = rw, cl = cl, e = e)
    }
  }

  # One tile's inputs, cropped only when the tile is about to run.
  tileInputs <- function(tile) {
    e <- tile$e
    list(rw = tile$rw, cl = tile$cl,
         dtmi = terra::crop(dtm, e),
         vegpi = .cropVegSoilList(vegp, e),
         soilci = .cropVegSoilList(soilc, e),
         svfi = terra::crop(dwt$svfR, e),
         hori = tan(terra::as.array(terra::crop(dwt$horR, e)) * pi / 180),
         slri = terra::crop(dwt$slrR, e),
         apri = terra::crop(dwt$aprR, e),
         twii = terra::crop(dwt$twiR, e))
  }
  # What every tile shares.
  common <- list(pointmodel = pointmodel, reqhgt = reqhgt, tfact = tfact, twiwet = twiwet, altcorrect = altcorrect,
                 out = out, splineveg = splineveg,
                 scalefactor = scalefactor, writeformat = writeformat, outdir = outdir,
                 rowWidth = nchar(as.character(rws)), colWidth = nchar(as.character(cls)))

  coresN <- .resolveCores(cores)
  written <- vector("list", length(tiles))
  if (show) {
    message("Tiles: ", length(tiles), if (coresN > 1) paste0(", on ", coresN, " cores"))
  }
  bar <- .progressBar(length(tiles), show)
  on.exit(bar$finish(), add = TRUE)

  if (coresN <= 1) {
    for (i in seq_along(tiles)) {
      written[[i]] <- .runGridTile(tileInputs(tiles[[i]]), common)
      bar$tick()
    }
  } else {
    oldplan <- future::plan(future::multisession, workers = coresN)
    on.exit(future::plan(oldplan), add = TRUE)
    oldopt <- options(future.globals.maxSize = +Inf)
    on.exit(options(oldopt), add = TRUE)

    # The tile runner is sent to each worker under a name that is not one of the
    # package's own. A global named as a function the session has attached is
    # not sent, the worker being left to find it by attaching the package, and
    # a worker cannot see an internal function that way.
    runTile <- .runGridTile

    # A tile is sent whenever a worker is free.
    running <- list()
    nextTile <- 1L
    while (nextTile <= length(tiles) || length(running) > 0) {
      while (length(running) < coresN && nextTile <= length(tiles)) {
        job <- .wrapForWorker(tileInputs(tiles[[nextTile]]))
        fut <- future::future(runTile(job, common), seed = TRUE,
                              globals = list(job = job, common = common, runTile = runTile))
        running[[length(running) + 1]] <- list(i = nextTile, fut = fut)
        nextTile <- nextTile + 1L
      }
      done <- vapply(running, function(r) future::resolved(r$fut), logical(1))
      for (r in running[done]) {
        written[[r$i]] <- future::value(r$fut)
        bar$tick()
      }
      running <- running[!done]
      if (!any(done)) Sys.sleep(0.05)
    }
  }
  manifest <- unlist(written, use.names = FALSE)
  if (is.null(manifest)) manifest <- character(0)

  invisible(manifest)
}

# Run one tile of rungridmodelbig() and write it to disk; returns the file paths.
.runGridTile <- function(job, common) {
  dtmi <- .unwrapFromWorker(job$dtmi)

  # Each tile is itself solved sequentially; any requested parallelism is
  # across independent tiles rather than nested within them.
  resulti <- rungridmodel(common$pointmodel, reqhgt = common$reqhgt, dtm = dtmi,
                          vegp = .unwrapFromWorker(job$vegpi), soilc = .unwrapFromWorker(job$soilci),
                          tfact = common$tfact, twiwet = common$twiwet,
                          altcorrect = common$altcorrect, hor = job$hori,
                          twi = .unwrapFromWorker(job$twii), svf = .unwrapFromWorker(job$svfi),
                          slr = .unwrapFromWorker(job$slri), apr = .unwrapFromWorker(job$apri),
                          out = common$out, splineveg = common$splineveg, cores = "off")

  scalefactor <- common$scalefactor
  if (!is.null(scalefactor)) {
    resulti <- lapply(resulti, function(r) round(r * scalefactor))
  }

  rwt <- formatC(job$rw, width = common$rowWidth, flag = "0")
  clt <- formatC(job$cl, width = common$colWidth, flag = "0")
  basefn <- file.path(common$outdir, paste0("tile_", rwt, "_", clt))

  fns <- character(0)
  if ("rds" %in% common$writeformat) {
    fn <- paste0(basefn, ".rds")
    saveRDS(list(result = .wrapForWorker(resulti), scalefactor = scalefactor,
                 dtm = terra::wrap(dtmi)),
            fn)
    fns <- c(fns, fn)
  }
  # A NetCDF file holds only the fields with a value at the requested height.
  hasValues <- vapply(resulti, inherits, logical(1), what = "SpatRaster")
  if ("nc" %in% common$writeformat && any(hasValues)) {
    fn <- paste0(basefn, ".nc")
    sdsObj <- terra::sds(resulti[hasValues])
    names(sdsObj) <- names(resulti)[hasValues]
    if (!is.null(scalefactor)) {
      terra::writeCDF(sdsObj, fn, overwrite = TRUE, datatype = "INT4S")
    } else {
      terra::writeCDF(sdsObj, fn, overwrite = TRUE)
    }
    fns <- c(fns, fn)
  }
  fns
}

#' Blend overlapping grid-model tiles
#'
#' Recombines tiles written by \code{\link{rungridmodelbig}} and averages
#' smoothly across their overlapping margins. Use this when tiles were saved
#' with \code{toverlap > 0}; without overlap, adjacent tiles can simply abut and
#' there is no transition zone to blend.
#'
#' @param path Directory containing tile files.
#' @param field Output field to blend for \code{format = "rds"}.
#' @param format Tile format, \code{"rds"} or \code{"nc"}.
#'
#' @return A blended \code{SpatRaster} for RDS input, or a named list of
#'   blended \code{SpatRaster}s for NetCDF input. A field with no physical
#'   meaning at the height the tiles were run for is returned as a single
#'   \code{NA} for RDS input, and is absent from NetCDF tiles.
#' @export
mosaicblend <- function(path, field = NULL, format = c("rds", "nc")) {
  format <- match.arg(format)
  if (format == "rds" && is.null(field)) {
    stop("mosaicblend(): `field` is required for format = \"rds\" -- the ",
         "name of the single output field to blend, e.g. field = \"Tz\"", call. = FALSE)
  }

  pattern <- if (format == "rds") "^tile_.*\\.rds$" else "^tile_.*\\.nc$"
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("mosaicblend(): no tile files matching \"", pattern, "\" found in ", path,
         call. = FALSE)
  }

  # Recover the spatial arrangement of tiles from their row/column labels so
  # blending follows landscape position rather than filesystem ordering.
  base <- basename(files)
  m <- regmatches(base, regexec("^tile_([0-9]+)_([0-9]+)\\.", base))
  rw <- vapply(m, function(x) as.integer(x[2]), integer(1))
  cl <- vapply(m, function(x) as.integer(x[3]), integer(1))
  if (anyNA(rw) || anyNA(cl)) {
    stop("mosaicblend(): one or more file names in ", path, " don't match the ",
         "expected tile_<row>_<col> naming convention. The folder should hold only ",
         "the tiles written by rungridmodelbig().", call. = FALSE)
  }

  readOneField <- function(fn, fld) {
    if (format == "rds") {
      d <- readRDS(fn)
      terra::unwrap(d$result[[fld]])
    } else {
      terra::rast(fn, subds = fld)
    }
  }
  fieldNames <- if (format == "nc") {
    d0 <- terra::sds(files[1])
    names(d0)
  } else {
    field
  }

  blendOneField <- function(fld) {
    rowRasters <- vector("list", max(rw))
    for (r in sort(unique(rw))) {
      idx <- which(rw == r)
      idx <- idx[order(cl[idx])]
      rasters <- lapply(files[idx], readOneField, fld = fld)
      merged <- rasters[[1]]
      if (length(rasters) > 1) {
        for (k in 2:length(rasters)) merged <- .blendmosaic(merged, rasters[[k]], direction = "right")
      }
      rowRasters[[r]] <- merged
    }
    rowRasters <- rowRasters[!vapply(rowRasters, is.null, logical(1))]
    full <- rowRasters[[1]]
    if (length(rowRasters) > 1) {
      for (k in 2:length(rowRasters)) full <- .blendmosaic(full, rowRasters[[k]], direction = "bottom")
    }
    full
  }

  if (format == "rds") {
    # A field with no value at the height the tiles were run for is a single NA
    # in every tile, and is returned as one.
    first <- readRDS(files[1])$result[[field]]
    if (is.atomic(first) && length(first) == 1 && is.na(first)) return(NA)
    blendOneField(field)
  } else {
    out <- lapply(fieldNames, blendOneField)
    names(out) <- fieldNames
    out
  }
}

# Spatial microclimate from one shared climate reference and fixed vegetation.
# The point-model time series anchors the atmospheric and surface state for the
# whole landscape. Each grid cell then modifies that reference according to its
# own vegetation, radiation geometry, exposure, soil and topographic wetness.
# The calculation proceeds from radiation/wind, through soil moisture and ground
# temperature, to canopy temperature and finally the requested vertical profile.
.gridmodelCore1 <- function(pointmodel, vegp, soilc, dtm, reqhgt = 0, tfact = 1.7, twiwet = 10,
                           hor = NA, twi = NA, svf = NA, slr = NA, apr = NA, out = NULL) {
  if (is.null(pointmodel$weather) || is.null(pointmodel$weather$obs_time)) {
    stop(.incompletePointmodelMessage(), call. = FALSE)
  }
  weather <- pointmodel$weather
  tsteps <- nrow(weather)

  # Requested outputs determine which expensive physical stages are needed;
  # quantities required only as intermediate state are still calculated when a
  # downstream requested field depends on them.
  .outFieldNames <- .GRID_OUT_FIELDS
  if (is.null(out)) {
    out <- stats::setNames(rep(TRUE, length(.outFieldNames)), .outFieldNames)
  } else {
    if (is.null(names(out)) || any(!nzchar(names(out)))) {
      stop(.outNotNamedMessage(.outFieldNames), call. = FALSE)
    }
    unknown <- setdiff(names(out), .outFieldNames)
    if (length(unknown) > 0) {
      stop(.outUnknownMessage(unknown, .outFieldNames), call. = FALSE)
    }
    # Opt-in scoping: fields the caller didn't mention default to FALSE
    # here (the opposite of a partially-supplied vegp/soilc) -- the whole
    # point of `out` is that out = c(soilm = TRUE) means "soil moisture
    # only", not "soil moisture plus everything else".
    full <- stats::setNames(rep(FALSE, length(.outFieldNames)), .outFieldNames)
    full[names(out)] <- as.logical(out)
    out <- full
  }

  if (inherits(dtm, "PackedSpatRaster")) dtm <- terra::rast(dtm)
  rows <- dim(dtm)[1]; cols <- dim(dtm)[2]
  # A single non-spatial weather series has one solar-geometry reference: where its
  # radiation was measured, which the reference run records as the centre of the
  # whole study area. It is shared by every fine cell and every tile.
  ll <- .referenceLatLon(pointmodel, dtm)
  latM <- matrix(ll$lat, 1, 1)
  lonM <- matrix(ll$lon, 1, 1)
  naArr3 <- array(NA_real_, dim = c(rows, cols, tsteps)) # shared placeholder for any field `out` skips

  # Put vegetation properties onto the fine grid so every cell carries its own
  # canopy structure and optical properties; a seasonal field gives its middle
  # layer, as the reference run takes.
  tomat <- function(v) {
    if (inherits(v, "PackedSpatRaster")) v <- terra::rast(v)
    if (inherits(v, "SpatRaster")) return(terra::as.matrix(v[[ceiling(terra::nlyr(v) / 2)]], wide = TRUE))
    if (is.matrix(v)) return(v)
    matrix(v, nrow = rows, ncol = cols) # scalar -> constant matrix
  }
  hgt <- tomat(vegp$hgt)
  pai <- tomat(vegp$pai)
  x <- tomat(vegp$x)
  lref <- tomat(vegp$leafr)
  ltra <- tomat(vegp$leaft)
  clump <- tomat(vegp$clump)

  # Slope and aspect control direct-beam interception. Derivatives are evaluated
  # on a one-cell-extended terrain surface so valid edge cells are not implicitly
  # treated as flat simply because neighbours outside the supplied extent are absent.
  if (is.logical(slr) || is.logical(apr)) {
    fill_ext <- terravars:::.fillna_dtm(dtm)
  }
  slope <- if (is.logical(slr)) {
    terra::as.matrix(terra::crop(terra::terrain(fill_ext$extended, "slope", unit = "degrees"), dtm), wide = TRUE)
  } else {
    tomat(slr)
  }
  aspect <- if (is.logical(apr)) {
    terra::as.matrix(terra::crop(terra::terrain(fill_ext$extended, "aspect", unit = "degrees"), dtm), wide = TRUE)
  } else {
    tomat(apr)
  }

  # Resolve vegetation type for each vegetated cell and use its characteristic
  # vertical foliage profile to determine how much plant area lies above the
  # requested height. Bare ground has no foliage and therefore no PFT or overlying PAI.
  # A plant of zero height or zero plant area has no mass, so either one zero is
  # bare ground, and both are set to zero.
  bareground <- !is.na(hgt) & (hgt <= 0 | pai <= 0)
  hgt[bareground] <- 0
  pai[bareground] <- 0
  vegcell <- !is.na(hgt) & !bareground
  group <- matrix(NA_character_, nrow = rows, ncol = cols)
  pftM <- matrix(NA_character_, nrow = rows, ncol = cols)
  if (any(vegcell)) {
    pft <- habitattoPFT(as.vector(tomat(vegp$habitat))[vegcell], ll$lat)
    group[vegcell] <- .PFT_GROUP_OF[pft]
    pftM[vegcell] <- pft
  }
  # Leaf temperature requires PFT-specific photosynthetic, hydraulic and leaf
  # geometry parameters in addition to the canopy structural fields above.
  leafphys <- .leafPhysiologyMatrices(pftM, vegcell, rows, cols, vegp)

  # On sloping ground, convert PAI expressed per horizontal ground area to the
  # inclined-area convention used by the slope-aware canopy radiation geometry.
  # PFT classification above deliberately uses the unconverted structural PAI.
  if (isTRUE(pointmodel$paiFlat)) {
    pai <- pai * cos(slope * pi / 180)
  }

  paia <- matrix(NA_real_, nrow = rows, ncol = cols)
  paia[bareground] <- 0
  if (any(vegcell)) paia[vegcell] <- .paiaboveheightPFT(reqhgt, hgt[vegcell], pai[vegcell], group[vegcell])
  # Ground reflectance sets the lower boundary for shortwave exchange beneath the canopy.
  grefm <- tomat(soilc$groundr)

  # Soil type sets the local residual and saturated water contents. Topographic
  # wetness then redistributes the reference soil-water state so convergent terrain
  # is relatively wetter and divergent terrain relatively drier.
  sminsmax <- .soilcSminSmax(soilc$soiltype)
  Smin <- tomat(sminsmax$Smin)
  Smax <- tomat(sminsmax$Smax)
  wetness <- .topoWetness(dtm, twi, tfact, twiwet)
  tadd <- tomat(wetness$shift)
  wetw <- tomat(wetness$wet)

  # Terrain sky exposure controls diffuse and longwave exchange, while the local
  # horizon determines whether the direct solar beam is blocked at each timestep.
  sky <- .horizonAndSkyview(dtm, hor, svf)
  svfa <- tomat(sky$svf)
  hora <- sky$hora # tangents: runmicro1Cpp compares against tan(solar altitude)

  # Calendar date and time are retained explicitly because solar position and
  # terrain shadowing change through the day and year.
  tme <- as.POSIXlt(weather$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year + 1900, month = tme$mon + 1, day = tme$mday,
                         hour = tme$hour + tme$min / 60 + tme$sec / 3600)
  if (is.null(weather$windspeed) || is.null(weather$winddir)) {
    stop(.incompletePointmodelMessage(), call. = FALSE)
  }
  climdata <- data.frame(temp = weather$temp, swdown = weather$swdown,
                          difrad = weather$difrad, lwdown = weather$lwdown,
                          windspeed = weather$windspeed,
                          relhum = weather$relhum, pres = weather$pres)

  vegpl <- list(hgt = hgt, pai = pai, paia = paia, x = x, leafr = lref, leaft = ltra, clump = clump)
  soilcl <- list(gref = grefm, slope = slope, aspect = aspect, svfa = svfa, hor = hora)

  # Work backwards from the requested outputs to the physical state they require.
  # For example, leaf temperature needs the local air profile and longwave field;
  # canopy temperature needs ground heat flux; ground heat flux needs soil moisture.
  reqhgt_below0 <- !is.na(reqhgt) && reqhgt < 0
  reqhgt_eq0    <- !is.na(reqhgt) && reqhgt == 0
  reqhgt_above0 <- !is.na(reqhgt) && reqhgt > 0

  needLeafTemp     <- reqhgt_above0 && isTRUE(out[["tleaf"]])
  needLongwave     <- isTRUE(out[["Rlwdown"]]) || isTRUE(out[["Rlwup"]]) || needLeafTemp
  needProfile      <- reqhgt_above0 && (isTRUE(out[["Tz"]]) || isTRUE(out[["relhum"]]) || needLeafTemp)
  needCanopyTemp   <- isTRUE(out[["tleaf"]]) || needLongwave || needProfile
  needGroundTemp   <- (reqhgt_eq0 && isTRUE(out[["Tz"]])) || needCanopyTemp ||
    (reqhgt_below0 && isTRUE(out[["Tz"]])) || needLongwave
  needSurfaceSoilm <- needGroundTemp || (reqhgt_eq0 && isTRUE(out[["soilm"]]))
  needBelowSoilm   <- reqhgt_below0 && isTRUE(out[["soilm"]])
  # Radiation and wind form the common physical foundation for all above-ground
  # quantities and for the surface energy-balance estimates. They can be skipped
  # only when the requested result is soil moisture alone.
  needRaw <- any(unlist(out[setdiff(names(out), "soilm")]))
  zref <- pointmodel$zref
  if (needRaw) {
    if (is.null(pointmodel$model$uf)) stop(.incompletePointmodelMessage(), call. = FALSE)
    ufRef <- pointmodel$model$uf
    hRef <- pointmodel$vegp$h[1]
    paiRef <- pointmodel$vegp$pai[1]
    dRef <- zeroplanedisCpp(hRef, paiRef)
    zmRef <- roughlengthCpp(hRef, paiRef, dRef)
    sheltarr <- .windShelterInterp(hora, weather$winddir)

    raw <- runmicro1Cpp(obstime, climdata, vegpl, soilcl, latM, lonM,
                         zref, reqhgt, ufRef, pointmodel$model$Hout, dRef, zmRef, as.vector(sheltarr))
  } else {
    # `out` requests soil moisture alone -- skip the radiation/wind core
    # entirely (see needRaw's own comment above) and substitute NA
    # placeholders shaped like its real output for the field names
    # unconditionally read below.
    naMat <- matrix(NA_real_, rows, cols)
    raw <- list(radGsw = naArr3, radGlw = naArr3, radCsw = naArr3, radClw = naArr3,
                emGround = naArr3, emCanopy = naArr3, rGh = naArr3,
                rGreq = naArr3, zs = matrix(NA_real_, rows, cols),
                shapeR = matrix(NA_real_, rows, cols), shapeC = matrix(NA_real_, rows, cols),
                Rbdown = naArr3, Rddown = naArr3, Rdup = naArr3,
                uz = naArr3, uf = naArr3, rGz = naArr3, rGm = naArr3, rHa = naArr3, a2 = naArr3, L = naArr3,
                d = naMat, zm = naMat)
  }

  # Distribute the reference surface soil-water trajectory across cells according
  # to each soil's water-content limits and its topographic wetness anomaly.
  if (needSurfaceSoilm) {
    if (is.null(pointmodel$model$theta0)) stop(.incompletePointmodelMessage(), call. = FALSE)
    raw$soilm <- soilmDistributeCpp(Smin, Smax, tadd, wetw,pointmodel$model$theta0, tsteps)
  } else {
    raw$soilm <- naArr3
  }

  # Estimate each cell's ground heat flux and surface temperature by scaling the
  # reference thermal response for local soil properties and local absorbed energy,
  # while shifting the diurnal response to reflect slope/aspect effects on illumination.
  if (needGroundTemp) {
    campbell <- .soilcCampbellParams(soilc$soiltype)
    RabsGround <- raw$radGsw + raw$radGlw
    if (is.null(pointmodel$model$wetShare) || is.null(pointmodel$model$filmAvailable)) {
      stop(.incompletePointmodelMessage(), call. = FALSE)
    }
    campbellMats <- lapply(campbell, tomat)
    sr <- .cellSurfaceResistance(raw, weather$swdown, weather$difrad, weather$temp, weather$relhum,
                                 weather$pres, weather$precip, obstime$year, raw$soilm, campbellMats,
                                 hgt, pai, x, lref, ltra, svfa, .gsmaxMatrix(vegp, rows, cols),
                                 leafphys, tomat(vegp$Lfrac), pointmodel$model$wetShare,
                                 pointmodel$model$filmAvailable, paiRef)
    gh <- groundHeatFluxCpp(RabsGround, raw$rGz, raw$rGm,
                             raw$radCsw + raw$radClw, raw$rHa, sr$rSurf, sr$hSurf, sr$hFol,
                             raw$soilm,
                             weather$temp, weather$relhum, weather$pres,
                             tomat(campbell$Vq), tomat(campbell$Vm), tomat(campbell$Vo),
                             tomat(campbell$Mc), tomat(campbell$thetaS), tomat(campbell$psie),
                             tomat(campbell$b), svfa, raw$emGround, raw$emCanopy,
                             slope, aspect, weather$swdown, weather$difrad,
                             latM, lonM,
                             obstime$year, obstime$month, obstime$day, obstime$hour,
                             pointmodel$model$G, pointmodel$model$RabsGround,
                             referenceGroundResistCpp(weather$windspeed, pointmodel$model$uf,
                                                      pointmodel$model$Hout, weather$temp,
                                                      weather$pres, zref, hRef, paiRef),
                             pointmodel$model$theta0, pointmodel$model$emGround,
                             pointmodel$soilc$Vq[1], pointmodel$soilc$Vm[1],
                             pointmodel$soilc$Vo[1], pointmodel$soilc$Mc[1],
                             pointmodel$soilc$Smax[1], pointmodel$soilc$psi_e[1],
                             pointmodel$soilc$b[1])
    raw$G_est <- gh$G_est
    raw$Ts_est <- gh$Ts_est
  } else {
    raw$G_est <- naArr3
    raw$Ts_est <- naArr3
  }

  # Estimate the canopy/ground exchange-surface temperature from local absorbed
  # radiation, aerodynamic coupling and ground heat flux. Where foliage is present
  # the ground pass has formed it, jointly with the ground or as a one-surface
  # estimate, and that is used; the one-surface estimate here serves only
  # vegetated cells the pass leaves, and bare cells take the ground's own
  # temperature. The reported ground temperature is then blended toward it for
  # vegetation shorter than 10 mm (.coupledSurfaceState).
  if (needGroundTemp) {
    RabsCanopy <- .oneSurfaceRabs(raw$radCsw + raw$radClw, gh, hgt)
    TcOne <- canopyTempCpp(RabsCanopy, raw$rHa, raw$G_est,
                           weather$temp, weather$relhum, weather$pres,
                           svfa, raw$emCanopy, sr$rSurf, sr$hSurf)
    st <- .coupledSurfaceState(gh, TcOne, sr, hgt)
    raw$Tcanopy_est <- st$Tcanopy_est
    raw$Ts_est <- st$Ts_est
    sr$hSurf <- st$hSurf
  } else {
    raw$Tcanopy_est <- naArr3
  }

  # Assemble outputs according to the physical position of reqhgt. Below ground
  # only soil state is meaningful; at the surface the ground energy balance and
  # radiation apply; above ground each cell is classified as within or above its
  # own canopy before temperature, humidity, wind and radiation are selected.
  if (!is.na(reqhgt) && reqhgt < 0) {
    # Below-ground temperature retains the reference soil solution at the same
    # depth, adjusted according to how the local surface-temperature forcing differs
    # from the reference. The influence of surface differences is damped with depth.
    if (isTRUE(out[["Tz"]])) {
      if (is.null(pointmodel$model$Tzbelow)) {
        stop(.depthMismatchMessage(reqhgt, pointmodel$reqhgt), call. = FALSE)
      }
      if (!is.null(pointmodel$reqhgt) && !isTRUE(all.equal(pointmodel$reqhgt, reqhgt))) {
        stop(.depthMismatchMessage(reqhgt, pointmodel$reqhgt), call. = FALSE)
      }
      DD <- dampingDepthGridCpp(tomat(campbell$Vq), tomat(campbell$Vm), tomat(campbell$Vo),
                                tomat(campbell$Mc), tomat(campbell$thetaS), tomat(campbell$psie),
                                tomat(campbell$b), raw$soilm, raw$Ts_est, weather$pres)
      Tzbelow_est <- belowGroundShortcutCpp(raw$Ts_est, pointmodel$model$Tground,
                                             pointmodel$model$Tzbelow, DD,
                                             pointmodel$matemp, reqhgt, 8760)
    } else {
      Tzbelow_est <- naArr3
    }
    # Soil moisture at depth uses the reference moisture trajectory at that same
    # depth, redistributed spatially with the same soil/topographic constraints.
    soilm_below <- if (needBelowSoilm) {
      soilmDistributeCpp(Smin, Smax, tadd, wetw,pointmodel$model$thetazbelow, tsteps)
    } else {
      naArr3
    }

    result <- list(Tz = Tzbelow_est, tleaf = naArr3, soilm = soilm_below,
                    relhum = naArr3, windspeed = naArr3,
                    Rdirdown = naArr3, Rdifdown = naArr3, Rswup = naArr3,
                    Rlwdown = naArr3, Rlwup = naArr3)
  } else if (!is.na(reqhgt) && reqhgt == 0) {
    # At the ground surface, longwave radiation includes sky transmission through
    # the canopy and thermal emission from the canopy and ground.
    lw <- if (needLongwave) {
      longwaveGridCpp(raw$Ts_est, raw$Tcanopy_est, weather$temp, hgt, pai, paia, svfa, weather$lwdown, reqhgt)
    } else {
      list(Rlwdown = naArr3, Rlwup = naArr3)
    }
    result <- list(Tz = raw$Ts_est, tleaf = naArr3, soilm = raw$soilm,
                    relhum = naArr3, windspeed = naArr3,
                    Rdirdown = raw$Rbdown, Rdifdown = raw$Rddown, Rswup = raw$Rdup,
                    Rlwdown = lw$Rlwdown, Rlwup = lw$Rlwup)
  } else {
    # Above the canopy, Monin-Obukhov similarity links the local exchange surface to
    # air temperature and humidity at reqhgt. Within the canopy, a resistance-based
    # profile links the ground surface to the air above the canopy.
    ac <- if (needProfile) {
      aboveCanopyProfileGridCpp(raw$Tcanopy_est, raw$rHa, raw$L, raw$uf, raw$d, raw$zh,
                                 weather$temp, weather$relhum, sr$rSurf, sr$hSurf,
                                 raw$rGz, raw$rGreq, raw$zs,
                                 reqhgt, zref)
    } else {
      list(Tz = naArr3, RHz = naArr3)
    }

    bc <- if (needProfile) {
      belowCanopyProfileGridCpp(
        raw$Tcanopy_est, raw$Ts_est, raw$rHa, raw$rGh, raw$rGm, raw$rGz,
        raw$shapeR, raw$shapeC, hgt, pai,
        weather$temp, weather$relhum, weather$pres, sr$rSurf,
        sr$hSurf,
        reqhgt, zref)
    } else {
      list(Tz = naArr3, RHz = naArr3)
    }

    tsteps <- dim(raw$Tcanopy_est)[3]
    # The same requested height may be above short vegetation but within taller
    # vegetation, so select the appropriate vertical profile independently by cell.
    above_mask <- rep(as.vector(hgt <= reqhgt), tsteps)

    Tz_pc <- ifelse(above_mask, as.vector(ac$Tz), as.vector(bc$Tz))
    relhum_pc <- ifelse(above_mask, as.vector(ac$RHz), as.vector(bc$RHz))
    dim(Tz_pc) <- dim(raw$Tcanopy_est)
    dim(relhum_pc) <- dim(raw$Tcanopy_est)

    # Evaluate longwave exchange at the requested height using the amount of canopy
    # above and below that level.
    lw <- if (needLongwave) {
      longwaveGridCpp(raw$Ts_est, raw$Tcanopy_est, weather$temp, hgt, pai, paia, svfa, weather$lwdown, reqhgt)
    } else {
      list(Rlwdown = naArr3, Rlwup = naArr3)
    }

    # For leaves present at reqhgt, solve leaf temperature from local shortwave/PAR,
    # longwave exchange, wind, local air temperature/humidity and plant hydraulic/
    # photosynthetic properties. This is distinct from the bulk canopy exchange temperature.
    tleaf_pc <- if (needLeafTemp) {
      psiCell <- soilWaterPotentialGridCpp(tomat(campbell$thetaS), tomat(campbell$psie),
                                           tomat(campbell$b), raw$soilm)
      radLlw <- 0.5 * 0.97 * (lw$Rlwdown + lw$Rlwup) # 0.97 = mc::surfaceEmissivity (src/constants.h),
                                                       # the same shared canopy/ground emissivity value
                                                       # used throughout this package -- no R-side
                                                       # binding to the C++ constant exists, hardcoded
                                                       # here matching soilc$groundem's own convention
                                                       # (R/utils.R)
      Ca <- rep(Cafromyear(stats::median(obstime$year)), tsteps)
      # Use the local air state at reqhgt, so the leaf energy balance experiences
      # the same above- or below-canopy microclimate being returned for that cell.
      leafTempCpp(raw$radLpar, raw$radLsw, radLlw, raw$uz,
                  Tz_pc, relhum_pc, weather$pres, Ca, psiCell,
                  hgt, pai, leafphys$Vcmax25, leafphys$Tup, leafphys$Tlw, leafphys$Dcrit,
                  leafphys$alpha, leafphys$f0, leafphys$fd, leafphys$psi50, leafphys$apsi,
                  leafphys$rpmin, leafphys$leafd, leafphys$isC3, reqhgt)
    } else {
      naArr3
    }

    result <- list(Tz = Tz_pc, tleaf = tleaf_pc, soilm = naArr3,
                    relhum = relhum_pc, windspeed = raw$uzActual,
                    Rdirdown = raw$Rbdown, Rdifdown = raw$Rddown, Rswup = raw$Rdup,
                    Rlwdown = lw$Rlwdown, Rlwup = lw$Rlwup)
  }

  # Return only fields requested by the user; intermediate state used to obtain
  # those fields is not exposed.
  result <- result[out[names(result)]]
  result
}

# Convert the raw spatial/time arrays from the single-reference calculation into
# georeferenced raster time series matching the input terrain.
.rungridmodel1 <- function(pointmodel, vegp, soilc, dtm, reqhgt = 0, tfact = 1.7, twiwet = 10,
                           hor = NA, twi = NA, svf = NA, slr = NA, apr = NA, out = NULL) {
  result <- .gridmodelCore1(pointmodel, vegp, soilc, dtm, reqhgt, tfact, twiwet,
                             hor, twi, svf, slr, apr, out)
  # Restore the terrain extent and coordinate reference system on every output field.
  .gridOutputRasters(result, dtm, .fieldsAtHeight(reqhgt))
}

# Spatial microclimate from coarse-gridded climate references and fixed vegetation.
# Each fine cell inherits the time-varying atmospheric and surface state of its
# surrounding coarse climate field, with smooth resampling where appropriate, then
# receives the same local terrain, vegetation, soil and vertical-profile adjustments
# as the single-reference pathway.
.gridmodelCore2 <- function(pointmodela, vegp, soilc, dtm, reqhgt = 0, tfact = 1.7, twiwet = 10,
                            altcorrect = 0, hor = NA, twi = NA, svf = NA,
                            slr = NA, apr = NA, out = NULL) {
  # The coarse grid defines where independent reference energy-balance solutions
  # exist; fine cells are subsequently related to those spatial climate references.
  if (is.null(pointmodela$points) || is.null(pointmodela$dtmc)) {
    stop("`pointmodel` must be the result of runpointmodel() or subsetpointmodel().",
         call. = FALSE)
  }
  points <- pointmodela$points
  dtmc <- .unpackRaster(pointmodela$dtmc)

  # Requested outputs determine the minimum set of physical stages that must be run.
  .outFieldNames <- .GRID_OUT_FIELDS
  if (is.null(out)) {
    out <- stats::setNames(rep(TRUE, length(.outFieldNames)), .outFieldNames)
  } else {
    if (is.null(names(out)) || any(!nzchar(names(out)))) {
      stop(.outNotNamedMessage(.outFieldNames), call. = FALSE)
    }
    unknown <- setdiff(names(out), .outFieldNames)
    if (length(unknown) > 0) {
      stop(.outUnknownMessage(unknown, .outFieldNames), call. = FALSE)
    }
    full <- stats::setNames(rep(FALSE, length(.outFieldNames)), .outFieldNames)
    full[names(out)] <- as.logical(out)
    out <- full
  }

  if (inherits(dtm, "PackedSpatRaster")) dtm <- terra::rast(dtm)
  if (isTRUE(all.equal(dim(dtm)[1:2], dim(dtmc)[1:2]))) {
    stop(.sameResolutionMessage(), call. = FALSE)
  }

  # Coarse cells without land or usable climate forcing have no physical reference
  # solution and are excluded from the spatial interpolation.
  isValid <- vapply(points, is.list, logical(1))
  validCells <- which(isValid)
  if (length(validCells) == 0) stop(.noUsableRunMessage(), call. = FALSE)

  # Gather the weather, solved surface state, reference canopy geometry and soil
  # properties for every valid coarse cell. These are the quantities fine cells
  # inherit from the large-scale climate field before local corrections are applied.
  weatherList <- vector("list", length(points))
  modelList <- vector("list", length(points))
  matemp <- rep(NA_real_, length(points))
  zrefV <- rep(NA_real_, length(points))
  hRefV <- rep(NA_real_, length(points))
  paiRefV <- rep(NA_real_, length(points))
  # Retain each reference soil texture because local ground thermal response is
  # scaled relative to the soil represented by its coarse-cell point solution.
  VqRefV <- rep(NA_real_, length(points))
  VmRefV <- rep(NA_real_, length(points))
  VoRefV <- rep(NA_real_, length(points))
  McRefV <- rep(NA_real_, length(points))
  thetaSRefV <- rep(NA_real_, length(points))
  psieRefV <- rep(NA_real_, length(points))
  bRefV <- rep(NA_real_, length(points))
  for (i in validCells) {
    weatherList[[i]] <- points[[i]]$weather
    modelList[[i]] <- points[[i]]$model
    matemp[i] <- points[[i]]$matemp
    zrefV[i] <- points[[i]]$zref
    hRefV[i] <- points[[i]]$vegp$h
    paiRefV[i] <- points[[i]]$vegp$pai
    VqRefV[i] <- points[[i]]$soilc$Vq[1]
    VmRefV[i] <- points[[i]]$soilc$Vm[1]
    VoRefV[i] <- points[[i]]$soilc$Vo[1]
    McRefV[i] <- points[[i]]$soilc$Mc[1]
    thetaSRefV[i] <- points[[i]]$soilc$Smax[1]
    psieRefV[i] <- points[[i]]$soilc$psi_e[1]
    bRefV[i] <- points[[i]]$soilc$b[1]
    # The reference's ground-to-reference resistance as the grid evaluates it,
    # so that its side of the ground heat flux scaling matches a cell's.
    modelList[[i]]$rGzGrid <- referenceGroundResistCpp(points[[i]]$weather$windspeed,
                                                       points[[i]]$model$uf, points[[i]]$model$Hout,
                                                       points[[i]]$weather$temp, points[[i]]$weather$pres,
                                                       points[[i]]$zref,
                                                       points[[i]]$vegp$h, points[[i]]$vegp$pai)
  }

  # Interpolate the coarse atmospheric state to the terrain grid. Temperature and
  # pressure can additionally be adjusted for the fine cell's elevation.
  nT <- nrow(weatherList[[validCells[1]]])
  tcFine <- .coarseListToFineArray(weatherList, "temp", nT, dtmc, dtm)
  relhumFine <- .coarseListToFineArray(weatherList, "relhum", nT, dtmc, dtm)
  pkCoarse <- .coarseListToArray(weatherList, "pres", nT, dtmc)
  alt <- .altitudeCorrectPkTc(pkCoarse, tcFine, relhumFine, dtmc, dtm, altcorrect)

  # Radiation fields interpolate directly; wind is interpolated as vector components
  # so directional wrap-around cannot create spurious winds.
  swdownFine <- .coarseListToFineArray(weatherList, "swdown", nT, dtmc, dtm)
  difradFine <- .coarseListToFineArray(weatherList, "difrad", nT, dtmc, dtm)
  lwdownFine <- .coarseListToFineArray(weatherList, "lwdown", nT, dtmc, dtm)
  wind <- .coarseWindToFine(weatherList, nT, dtmc, dtm)

  climdata <- list(temp = alt$tc, swdown = swdownFine, difrad = difradFine,
                    lwdown = lwdownFine, windspeed = wind$windspeed,
                    relhum = relhumFine, pres = alt$pk)

  # A common calendar underlies every coarse reference, preserving one consistent
  # solar geometry through the spatially varying climate field.
  tme <- as.POSIXlt(weatherList[[validCells[1]]]$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year + 1900, month = tme$mon + 1, day = tme$mday,
                         hour = tme$hour + tme$min / 60 + tme$sec / 3600)

  # As in the single-reference pathway, derive which physical stages are required
  # from the requested outputs and height/depth regime.
  reqhgt_below0 <- !is.na(reqhgt) && reqhgt < 0
  reqhgt_eq0    <- !is.na(reqhgt) && reqhgt == 0
  reqhgt_above0 <- !is.na(reqhgt) && reqhgt > 0
  needLeafTemp     <- reqhgt_above0 && isTRUE(out[["tleaf"]])
  needLongwave     <- isTRUE(out[["Rlwdown"]]) || isTRUE(out[["Rlwup"]]) || needLeafTemp
  needProfile      <- reqhgt_above0 && (isTRUE(out[["Tz"]]) || isTRUE(out[["relhum"]]) || needLeafTemp)
  needCanopyTemp   <- isTRUE(out[["tleaf"]]) || needLongwave || needProfile
  needGroundTemp   <- (reqhgt_eq0 && isTRUE(out[["Tz"]])) || needCanopyTemp ||
    (reqhgt_below0 && isTRUE(out[["Tz"]])) || needLongwave
  needSurfaceSoilm <- needGroundTemp || (reqhgt_eq0 && isTRUE(out[["soilm"]]))
  needBelowSoilm   <- reqhgt_below0 && isTRUE(out[["soilm"]])
  needRaw <- any(unlist(out[setdiff(names(out), "soilm")]))

  # Fine-grid vegetation and terrain remain local properties; only the climate and
  # reference surface state originate on the coarser grid.
  rows <- dim(dtm)[1]; cols <- dim(dtm)[2]
  naArr3 <- array(NA_real_, dim = c(rows, cols, nT)) # shared placeholder for a field `out`/reqhgt regime skips
  # A seasonal field gives its middle layer, as the reference runs take.
  tomat <- function(v) {
    if (inherits(v, "PackedSpatRaster")) v <- terra::rast(v)
    if (inherits(v, "SpatRaster")) return(terra::as.matrix(v[[ceiling(terra::nlyr(v) / 2)]], wide = TRUE))
    if (is.matrix(v)) return(v)
    matrix(v, nrow = rows, ncol = cols)
  }
  # Solar position is evaluated at each fine cell because the climate forcing is
  # genuinely spatial here, avoiding artificial changes at coarse-cell boundaries.
  llFine <- .latlonsFromCells(dtm)
  latR <- dtm; terra::values(latR) <- llFine$lat
  lonR <- dtm; terra::values(lonR) <- llFine$lon
  latM <- tomat(latR)
  lonM <- tomat(lonR)
  hgt <- tomat(vegp$hgt)
  pai <- tomat(vegp$pai)
  x <- tomat(vegp$x)
  lref <- tomat(vegp$leafr)
  ltra <- tomat(vegp$leaft)
  clump <- tomat(vegp$clump)

  # Resolve slope and aspect on an edge-extended terrain surface so boundary cells
  # retain physically plausible direct-beam geometry.
  if (is.logical(slr) || is.logical(apr)) {
    fill_ext <- terravars:::.fillna_dtm(dtm)
  }
  slope <- if (is.logical(slr)) {
    terra::as.matrix(terra::crop(terra::terrain(fill_ext$extended, "slope", unit = "degrees"), dtm), wide = TRUE)
  } else {
    tomat(slr)
  }
  aspect <- if (is.logical(apr)) {
    terra::as.matrix(terra::crop(terra::terrain(fill_ext$extended, "aspect", unit = "degrees"), dtm), wide = TRUE)
  } else {
    tomat(apr)
  }

  # Classify vegetation from each cell's habitat before any slope-area conversion
  # of PAI, then derive PFT-dependent leaf and canopy-profile properties.
  # Either zero is bare ground (a plant of zero height or plant area has no mass).
  bareground <- !is.na(hgt) & (hgt <= 0 | pai <= 0)
  hgt[bareground] <- 0
  pai[bareground] <- 0
  vegcell <- !is.na(hgt) & !bareground
  group <- matrix(NA_character_, nrow = rows, ncol = cols)
  pftM <- matrix(NA_character_, nrow = rows, ncol = cols)
  if (any(vegcell)) {
    pft <- habitattoPFT(as.vector(tomat(vegp$habitat))[vegcell], as.vector(latM)[vegcell])
    group[vegcell] <- .PFT_GROUP_OF[pft]
    pftM[vegcell] <- pft
  }
  # Attach the physiological and leaf-geometry parameters needed for leaf energy balance.
  leafphys <- .leafPhysiologyMatrices(pftM, vegcell, rows, cols, vegp)

  # Convert horizontal-area PAI to the inclined-area convention used by the
  # slope-aware radiation calculation.
  if (isTRUE(pointmodela$paiFlat)) {
    pai <- pai * cos(slope * pi / 180)
  }

  paia <- matrix(NA_real_, nrow = rows, ncol = cols)
  paia[bareground] <- 0
  if (any(vegcell)) paia[vegcell] <- .paiaboveheightPFT(reqhgt, hgt[vegcell], pai[vegcell], group[vegcell])
  grefm <- tomat(soilc$groundr)

  # Local soil class sets the physical moisture bounds; topographic wetness
  # shifts the reference soil-water state within those bounds across the fine grid.
  sminsmax <- .soilcSminSmax(soilc$soiltype)
  Smin <- tomat(sminsmax$Smin)
  Smax <- tomat(sminsmax$Smax)
  wetness <- .topoWetness(dtm, twi, tfact, twiwet)
  tadd <- tomat(wetness$shift)
  wetw <- tomat(wetness$wet)

  # Terrain controls exposure to diffuse/longwave sky radiation and whether the
  # direct solar beam is blocked by the local horizon. Vegetation shading is handled
  # separately by the canopy radiation calculation.
  sky <- .horizonAndSkyview(dtm, hor, svf)
  svfa <- tomat(sky$svf)
  hora <- sky$hora # tangents: runmicro1Cpp compares against tan(solar altitude)

  vegpl <- list(hgt = hgt, pai = pai, paia = paia, x = x, leafr = lref, leaft = ltra, clump = clump)
  soilcl <- list(gref = grefm, slope = slope, aspect = aspect, svfa = svfa, hor = hora)

  # Fine cells inherit the reference canopy geometry and friction velocity of the
  # surrounding coarse climate field. Displacement and roughness are then derived
  # after interpolation so the nonlinear aerodynamic geometry is evaluated locally.
  zrefFine <- .coarseScalarToFineMatrix(zrefV, dtmc, dtm)
  if (needRaw) {
    hRefFine <- .coarseScalarToFineMatrix(hRefV, dtmc, dtm)
    paiRefFine <- .coarseScalarToFineMatrix(paiRefV, dtmc, dtm)
    dRef <- zeroplanedisCpp(hRefFine, paiRefFine)
    zmRef <- roughlengthCpp(hRefFine, paiRefFine, dRef)
    ufRef <- .coarseListToFineArray(modelList, "uf", nT, dtmc, dtm)

    # Apply direction-dependent terrain shelter to the locally interpolated wind.
    zrefRepresentative <- mean(zrefV, na.rm = TRUE)
    sheltarr <- .windShelterInterp(hora, wind$winddir)
  } else {
    dRef <- matrix(NA_real_, rows, cols); zmRef <- matrix(NA_real_, rows, cols)
    ufRef <- naArr3
    sheltarr <- rep(NA_real_, rows * cols * nT)
    zrefRepresentative <- NA_real_
  }

  # Solve shortwave/longwave radiation and wind on the fine grid using local
  # vegetation/terrain and the interpolated coarse atmospheric/aerodynamic state.
  if (needRaw) {
    climdataFlat <- data.frame(temp = as.vector(climdata$temp), swdown = as.vector(climdata$swdown),
                                difrad = as.vector(climdata$difrad), lwdown = as.vector(climdata$lwdown),
                                windspeed = as.vector(climdata$windspeed), relhum = as.vector(climdata$relhum),
                                pres = as.vector(climdata$pres))
    raw <- runmicro1Cpp(obstime, climdataFlat, vegpl, soilcl, latM, lonM,
                         as.vector(zrefFine), reqhgt, as.vector(ufRef),
                         as.vector(.coarseListToFineArray(modelList, "Hout", nT, dtmc, dtm)),
                         as.vector(dRef), as.vector(zmRef), as.vector(sheltarr))
  } else {
    naMat <- matrix(NA_real_, rows, cols)
    raw <- list(radGsw = naArr3, radGlw = naArr3, radCsw = naArr3, radClw = naArr3,
                emGround = naArr3, emCanopy = naArr3, rGh = naArr3,
                rGreq = naArr3, zs = matrix(NA_real_, rows, cols),
                shapeR = matrix(NA_real_, rows, cols), shapeC = matrix(NA_real_, rows, cols),
                Rbdown = naArr3, Rddown = naArr3, Rdup = naArr3,
                uz = naArr3, uf = naArr3, rGz = naArr3, rGm = naArr3, rHa = naArr3, a2 = naArr3, L = naArr3,
                d = naMat, zm = naMat)
  }

  # Interpolate the reference soil-water trajectory and redistribute it within each
  # fine cell according to local soil limits and topographic wetness.
  if (needSurfaceSoilm) {
    theta0Fine <- .coarseListToFineArray(modelList, "theta0", nT, dtmc, dtm)
    raw$soilm <- soilmDistributeCpp(Smin, Smax, tadd, wetw,theta0Fine, nT)
  } else {
    theta0Fine <- naArr3
    raw$soilm <- naArr3
  }

  # Reconstruct local ground heat flux and surface temperature using local absorbed
  # energy and soil properties, referenced to the corresponding coarse-cell thermal state.
  if (needGroundTemp) {
    campbell <- .soilcCampbellParams(soilc$soiltype)
    GFine <- .coarseListToFineArray(modelList, "G", nT, dtmc, dtm)
    RabsGFine <- .coarseListToFineArray(modelList, "RabsGround", nT, dtmc, dtm)
    rGzFine <- .coarseListToFineArray(modelList, "rGzGrid", nT, dtmc, dtm)
    emGFine <- .coarseListToFineArray(modelList, "emGround", nT, dtmc, dtm)
    thetaSRefFine <- .coarseScalarToFineMatrix(thetaSRefV, dtmc, dtm)
    psieRefFine <- .coarseScalarToFineMatrix(psieRefV, dtmc, dtm)
    bRefFine <- .coarseScalarToFineMatrix(bRefV, dtmc, dtm)
    VqRefFine <- .coarseScalarToFineMatrix(VqRefV, dtmc, dtm)
    VmRefFine <- .coarseScalarToFineMatrix(VmRefV, dtmc, dtm)
    VoRefFine <- .coarseScalarToFineMatrix(VoRefV, dtmc, dtm)
    McRefFine <- .coarseScalarToFineMatrix(McRefV, dtmc, dtm)
    RabsGround <- raw$radGsw + raw$radGlw
    sr <- .cellSurfaceResistance(raw, climdata$swdown, climdata$difrad, climdata$temp, climdata$relhum,
                                 climdata$pres, .coarseListToFineArray(weatherList, "precip", nT, dtmc, dtm),
                                 obstime$year, raw$soilm, lapply(campbell, tomat),
                                 hgt, pai, x, lref, ltra, svfa, .gsmaxMatrix(vegp, rows, cols), leafphys,
                                 tomat(vegp$Lfrac),
                                 .coarseListToFineArray(modelList, "wetShare", nT, dtmc, dtm),
                                 .coarseListToFineArray(modelList, "filmAvailable", nT, dtmc, dtm),
                                 paiRefFine)
    rSurfFine <- sr$rSurf
    hSurfFine <- sr$hSurf
    gh <- groundHeatFluxCpp(RabsGround, raw$rGz, raw$rGm,
                             raw$radCsw + raw$radClw, raw$rHa, rSurfFine, hSurfFine, sr$hFol,
                             raw$soilm,
                             climdata$temp, climdata$relhum, climdata$pres,
                             tomat(campbell$Vq), tomat(campbell$Vm), tomat(campbell$Vo),
                             tomat(campbell$Mc), tomat(campbell$thetaS), tomat(campbell$psie),
                             tomat(campbell$b), svfa, raw$emGround, raw$emCanopy,
                             slope, aspect, climdata$swdown, climdata$difrad,
                             latM, lonM,
                             obstime$year, obstime$month, obstime$day, obstime$hour,
                             GFine, RabsGFine, rGzFine, theta0Fine, emGFine,
                             VqRefFine, VmRefFine, VoRefFine, McRefFine,
                             thetaSRefFine, psieRefFine, bRefFine)
    raw$G_est <- gh$G_est
    raw$Ts_est <- gh$Ts_est
  } else {
    raw$G_est <- naArr3
    raw$Ts_est <- naArr3
  }

  # Estimate local canopy exchange temperature using the fine-cell radiation and
  # aerodynamics: where foliage is present, as the ground pass formed it or, in
  # cells it leaves, the one-surface estimate; on bare ground, the ground's own
  # temperature (see .coupledSurfaceState).
  if (needGroundTemp) {
    RabsCanopy <- .oneSurfaceRabs(raw$radCsw + raw$radClw, gh, hgt)
    TcOne <- canopyTempCpp(RabsCanopy, raw$rHa, raw$G_est,
                           climdata$temp, climdata$relhum, climdata$pres,
                           svfa, raw$emCanopy, rSurfFine, hSurfFine)
    st <- .coupledSurfaceState(gh, TcOne, sr, hgt)
    raw$Tcanopy_est <- st$Tcanopy_est
    raw$Ts_est <- st$Ts_est
    hSurfFine <- st$hSurf
  } else {
    rSurfFine <- naArr3
    hSurfFine <- naArr3
    raw$Tcanopy_est <- naArr3
  }

  # Select soil, surface, within-canopy or above-canopy outputs according to reqhgt,
  # using the same physical regimes as the single-reference pathway.
  if (!is.na(reqhgt) && reqhgt < 0) {
    # The soil profile must have been solved at the same depth as the requested
    # output; a profile from another depth cannot be rescaled interchangeably.
    if (!is.null(points[[validCells[1]]]$reqhgt) &&
        !isTRUE(all.equal(points[[validCells[1]]]$reqhgt, reqhgt))) {
      stop(.depthMismatchMessage(reqhgt, points[[validCells[1]]]$reqhgt), call. = FALSE)
    }
    # Below ground, carry the reference soil solution to the local cell using the
    # same surface-forcing and topographic adjustments as the single-reference pathway.
    Tzbelow_est <- if (isTRUE(out[["Tz"]])) {
      TgroundFine <- .coarseListToFineArray(modelList, "Tground", nT, dtmc, dtm)
      TzbelowFine <- .coarseListToFineArray(modelList, "Tzbelow", nT, dtmc, dtm)
      DD <- dampingDepthGridCpp(tomat(campbell$Vq), tomat(campbell$Vm), tomat(campbell$Vo),
                                tomat(campbell$Mc), tomat(campbell$thetaS), tomat(campbell$psie),
                                tomat(campbell$b), raw$soilm, raw$Ts_est, climdata$pres)
      matFine <- .coarseScalarToFineMatrix(matemp, dtmc, dtm) + .meanAltitudeOffset(alt$tc, tcFine)
      belowGroundShortcutCpp(raw$Ts_est, TgroundFine, TzbelowFine, DD, matFine, reqhgt, 8760)
    } else {
      naArr3
    }

    soilm_below <- if (needBelowSoilm) {
      thetazbelowFine <- .coarseListToFineArray(modelList, "thetazbelow", nT, dtmc, dtm)
      soilmDistributeCpp(Smin, Smax, tadd, wetw,thetazbelowFine, nT)
    } else {
      naArr3
    }

    result <- list(Tz = Tzbelow_est, tleaf = naArr3, soilm = soilm_below,
                    relhum = naArr3, windspeed = naArr3,
                    Rdirdown = naArr3, Rdifdown = naArr3, Rswup = naArr3,
                    Rlwdown = naArr3, Rlwup = naArr3)
  } else if (!is.na(reqhgt) && reqhgt == 0) {
    # At the surface, return the local ground state and radiation field.
    lw <- if (needLongwave) {
      longwaveGridCpp(raw$Ts_est, raw$Tcanopy_est, climdata$temp, hgt, pai, paia, svfa, climdata$lwdown, reqhgt)
    } else {
      list(Rlwdown = naArr3, Rlwup = naArr3)
    }
    result <- list(Tz = raw$Ts_est, tleaf = naArr3, soilm = raw$soilm,
                    relhum = naArr3, windspeed = naArr3,
                    Rdirdown = raw$Rbdown, Rdifdown = raw$Rddown, Rswup = raw$Rdup,
                    Rlwdown = lw$Rlwdown, Rlwup = lw$Rlwup)
  } else {
    # Above ground, determine separately for each cell whether reqhgt is inside or
    # above the canopy and evaluate the corresponding air profile.
    ac <- if (needProfile) {
      aboveCanopyProfileGridCpp(raw$Tcanopy_est, raw$rHa, raw$L, raw$uf, raw$d, raw$zh,
                                 climdata$temp, climdata$relhum, rSurfFine, hSurfFine,
                                 raw$rGz, raw$rGreq, raw$zs,
                                 reqhgt, zrefRepresentative)
    } else {
      list(Tz = naArr3, RHz = naArr3)
    }


    bc <- if (needProfile) {
      belowCanopyProfileGridCpp(
        raw$Tcanopy_est, raw$Ts_est, raw$rHa, raw$rGh, raw$rGm, raw$rGz,
        raw$shapeR, raw$shapeC, hgt, pai,
        climdata$temp, climdata$relhum, climdata$pres, rSurfFine, hSurfFine,
        reqhgt, zrefRepresentative)
    } else {
      list(Tz = naArr3, RHz = naArr3)
    }

    above_mask <- rep(as.vector(hgt <= reqhgt), nT)
    Tz_pc <- ifelse(above_mask, as.vector(ac$Tz), as.vector(bc$Tz))
    relhum_pc <- ifelse(above_mask, as.vector(ac$RHz), as.vector(bc$RHz))
    dim(Tz_pc) <- dim(raw$Tcanopy_est)
    dim(relhum_pc) <- dim(raw$Tcanopy_est)

    lw <- if (needLongwave) {
      longwaveGridCpp(raw$Ts_est, raw$Tcanopy_est, climdata$temp, hgt, pai, paia, svfa, climdata$lwdown, reqhgt)
    } else {
      list(Rlwdown = naArr3, Rlwup = naArr3)
    }

    # Leaf temperature at reqhgt, from the radiation, wind, air temperature and
    # humidity at that height. The leaf's stomata respond to soil drying through
    # the soil water potential, calculated from each cell's own soil moisture
    # and soil type.
    tleaf_pc <- if (needLeafTemp) {
      psiCell <- soilWaterPotentialGridCpp(tomat(campbell$thetaS), tomat(campbell$psie),
                                           tomat(campbell$b), raw$soilm)
      radLlw <- 0.5 * 0.97 * (lw$Rlwdown + lw$Rlwup) # 0.97 = mc::surfaceEmissivity (src/constants.h)
      Ca <- rep(Cafromyear(stats::median(obstime$year)), nT)
      # Use local microclimate at reqhgt rather than the unmodified coarse forcing.
      leafTempCpp(raw$radLpar, raw$radLsw, radLlw, raw$uz,
                  Tz_pc, relhum_pc, climdata$pres, Ca, psiCell,
                  hgt, pai, leafphys$Vcmax25, leafphys$Tup, leafphys$Tlw, leafphys$Dcrit,
                  leafphys$alpha, leafphys$f0, leafphys$fd, leafphys$psi50, leafphys$apsi,
                  leafphys$rpmin, leafphys$leafd, leafphys$isC3, reqhgt)
    } else {
      naArr3
    }

    result <- list(Tz = Tz_pc, tleaf = tleaf_pc, soilm = naArr3,
                    relhum = relhum_pc, windspeed = raw$uzActual,
                    Rdirdown = raw$Rbdown, Rdifdown = raw$Rddown, Rswup = raw$Rdup,
                    Rlwdown = lw$Rlwdown, Rlwup = lw$Rlwup)
  }

  # Return only requested output fields.
  result <- result[out[names(result)]]
  result
}

# Convert coarse-reference raw arrays into georeferenced fine-grid raster time series.
.rungridmodel2 <- function(pointmodela, vegp, soilc, dtm, reqhgt = 0, tfact = 1.7, twiwet = 10,
                            altcorrect = 0, hor = NA, twi = NA, svf = NA,
                            slr = NA, apr = NA, out = NULL) {
  result <- .gridmodelCore2(pointmodela, vegp, soilc, dtm, reqhgt, tfact, twiwet,
                             altcorrect, hor, twi, svf, slr, apr, out)
  # Restore the terrain extent and coordinate reference system on each field.
  .gridOutputRasters(result, dtm, .fieldsAtHeight(reqhgt))
}

# Single shared climate reference with seasonally varying vegetation.
# Fixed terrain and soil properties are prepared once, while canopy structure and
# optical properties are resolved for each modelled day before the radiation, wind
# and energy-balance shortcuts are applied.
.rungridmodel3 <- function(pointmodel, vegp, soilc, dtm, reqhgt = 0, tfact = 1.7, twiwet = 10,
                           hor = NA, twi = NA, svf = NA, slr = NA, apr = NA, out = NULL,
                           splineveg = FALSE) {
  if (is.null(pointmodel$weather) || is.null(pointmodel$weather$obs_time)) {
    stop(.incompletePointmodelMessage(), call. = FALSE)
  }
  weather <- pointmodel$weather
  tsteps <- nrow(weather)
  if (tsteps %% 24 != 0) stop(.notWholeDaysMessage(tsteps), call. = FALSE)
  nDaysNeeded <- tsteps %/% 24L

  # Requested fields determine which physical stages need to be evaluated.
  .outFieldNames <- .GRID_OUT_FIELDS
  if (is.null(out)) {
    out <- stats::setNames(rep(TRUE, length(.outFieldNames)), .outFieldNames)
  } else {
    if (is.null(names(out)) || any(!nzchar(names(out)))) {
      stop(.outNotNamedMessage(.outFieldNames), call. = FALSE)
    }
    unknown <- setdiff(names(out), .outFieldNames)
    if (length(unknown) > 0) {
      stop(.outUnknownMessage(unknown, .outFieldNames), call. = FALSE)
    }
    full <- stats::setNames(rep(FALSE, length(.outFieldNames)), .outFieldNames)
    full[names(out)] <- as.logical(out)
    out <- full
  }

  if (inherits(dtm, "PackedSpatRaster")) dtm <- terra::rast(dtm)
  rows <- dim(dtm)[1]; cols <- dim(dtm)[2]
  # The solar-geometry reference of the single weather series (see .referenceLatLon).
  ll <- .referenceLatLon(pointmodel, dtm)
  latM <- matrix(ll$lat, 1, 1)
  lonM <- matrix(ll$lon, 1, 1)
  naArr3 <- array(NA_real_, dim = c(rows, cols, tsteps))

  tomat <- function(v) {
    if (inherits(v, "PackedSpatRaster")) v <- terra::rast(v)
    if (inherits(v, "SpatRaster")) return(terra::as.matrix(v[[1]], wide = TRUE))
    if (is.matrix(v)) return(v)
    matrix(v, nrow = rows, ncol = cols)
  }

  # Map each modelled day back onto the full seasonal sequence represented by the
  # vegetation layers, including runs that have been temporally subsetted.
  if (!is.null(pointmodel$nhoursFull)) {
    nDaysFull <- ceiling(pointmodel$nhoursFull / 24)
  } else if (!is.null(pointmodel$subs)) {
    nDaysFull <- ceiling(max(pointmodel$subs) / 24)
  } else {
    nDaysFull <- nDaysNeeded
  }
  dayPos <- integer(nDaysNeeded)
  for (d in seq_len(nDaysNeeded)) {
    if (!is.null(pointmodel$subs)) {
      dayPos[d] <- ceiling(pointmodel$subs[(d - 1) * 24 + 1] / 24)
    } else {
      dayPos[d] <- d
    }
  }

  # Each vegetation property can be fixed, assigned as discrete seasonal states, or
  # smoothly interpolated through the supplied states independently of other fields.
  vegFieldNames <- c("hgt", "pai", "x", "leafr", "leaft", "clump", "Lfrac")
  nLayersV <- stats::setNames(vapply(vegFieldNames, function(f) .nLayersOf(vegp[[f]]), integer(1)), vegFieldNames)

  fieldSetup <- list()
  for (f in vegFieldNames) {
    nl <- nLayersV[[f]]
    v <- vegp[[f]]
    if (nl == 1L) {
      fieldSetup[[f]] <- list(invariant = TRUE, matrix = tomat(v))
    } else if (!splineveg) {
      dayLyr <- ceiling(seq_len(nDaysFull) * nl / nDaysFull)
      fieldSetup[[f]] <- list(invariant = FALSE, raster = v, dayLyr = dayLyr)
    } else {
      knotX <- (seq_len(nl) - 0.5) * nDaysFull / nl
      knotY <- vapply(seq_len(nl), function(i) as.vector(.tomatLayer(v, i, rows, cols)), numeric(rows * cols))
      fit <- splineFitCpp(knotX, knotY)
      fieldSetup[[f]] <- list(invariant = FALSE, fit = fit)
    }
  }

  # Resolve the vegetation state for a given day; slope-dependent PAI conversion is
  # applied only after the seasonal value itself has been obtained.
  resolveFieldForDay <- function(f, dpos) {
    fs <- fieldSetup[[f]]
    m <- if (fs$invariant) {
      fs$matrix
    } else if (!splineveg) {
      .tomatLayer(vegp[[f]], fs$dayLyr[dpos], rows, cols)
    } else {
      matrix(splineEvalCpp(fs$fit, dpos), nrow = rows, ncol = cols)
    }
    if (f == "pai" && isTRUE(pointmodel$paiFlat)) m <- m * cos(slope * pi / 180)
    m
  }

  # PFT describes the enduring vegetation type rather than seasonal canopy state:
  # it comes from each cell's habitat, and whether a cell is vegetated from its
  # middle-season structure.
  midHgt <- .tomatLayer(vegp$hgt, ceiling(nLayersV[["hgt"]] / 2), rows, cols)
  midPai <- .tomatLayer(vegp$pai, ceiling(nLayersV[["pai"]] / 2), rows, cols)
  midX   <- .tomatLayer(vegp$x,   ceiling(nLayersV[["x"]] / 2),   rows, cols)
  midBareground <- !is.na(midHgt) & (midHgt <= 0 | midPai <= 0)
  midVegcell <- !is.na(midHgt) & !midBareground
  group <- matrix(NA_character_, nrow = rows, ncol = cols)
  pftM <- matrix(NA_character_, nrow = rows, ncol = cols)
  if (any(midVegcell)) {
    pft <- habitattoPFT(as.vector(tomat(vegp$habitat))[midVegcell], ll$lat)
    group[midVegcell] <- .PFT_GROUP_OF[pft]
    pftM[midVegcell] <- pft
  }
  # PFT-specific physiology is likewise fixed while canopy structure changes seasonally.
  leafphys <- .leafPhysiologyMatrices(pftM, midVegcell, rows, cols, vegp)

  # Terrain geometry and soil physical properties do not vary with seasonal vegetation
  # and are therefore resolved once for the whole run.
  grefm <- tomat(soilc$groundr)
  sminsmax <- .soilcSminSmax(soilc$soiltype)
  Smin <- tomat(sminsmax$Smin)
  Smax <- tomat(sminsmax$Smax)
  wetness <- .topoWetness(dtm, twi, tfact, twiwet)
  tadd <- tomat(wetness$shift)
  wetw <- tomat(wetness$wet)
  sky <- .horizonAndSkyview(dtm, hor, svf)
  svfa <- tomat(sky$svf)
  hora <- sky$hora
  if (is.logical(slr) || is.logical(apr)) {
    fill_ext <- terravars:::.fillna_dtm(dtm)
  }
  slope <- if (is.logical(slr)) {
    terra::as.matrix(terra::crop(terra::terrain(fill_ext$extended, "slope", unit = "degrees"), dtm), wide = TRUE)
  } else {
    tomat(slr)
  }
  aspect <- if (is.logical(apr)) {
    terra::as.matrix(terra::crop(terra::terrain(fill_ext$extended, "aspect", unit = "degrees"), dtm), wide = TRUE)
  } else {
    tomat(apr)
  }
  campbell <- .soilcCampbellParams(soilc$soiltype)
  campbellMats <- list(Vq = tomat(campbell$Vq), Vm = tomat(campbell$Vm), Vo = tomat(campbell$Vo),
                       Mc = tomat(campbell$Mc), thetaS = tomat(campbell$thetaS),
                       psie = tomat(campbell$psie), b = tomat(campbell$b))

  # Preserve the calendar explicitly for solar position while the shared weather
  # time series supplies the atmospheric forcing for every fine cell.
  tme <- as.POSIXlt(weather$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year + 1900, month = tme$mon + 1, day = tme$mday,
                         hour = tme$hour + tme$min / 60 + tme$sec / 3600)
  if (is.null(weather$windspeed) || is.null(weather$winddir)) {
    stop(.incompletePointmodelMessage(), call. = FALSE)
  }
  climdata <- data.frame(temp = weather$temp, swdown = weather$swdown,
                          difrad = weather$difrad, lwdown = weather$lwdown,
                          windspeed = weather$windspeed,
                          relhum = weather$relhum, pres = weather$pres)

  # The requested outputs and height/depth regime determine which physical stages
  # are needed; this dependency structure is constant through the seasonal run.
  reqhgt_below0 <- !is.na(reqhgt) && reqhgt < 0
  reqhgt_eq0    <- !is.na(reqhgt) && reqhgt == 0
  reqhgt_above0 <- !is.na(reqhgt) && reqhgt > 0
  # Leaf temperature is meaningful only above ground and requires both the local
  # air profile and longwave radiation at the requested height for its energy balance.
  needLeafTemp     <- reqhgt_above0 && isTRUE(out[["tleaf"]])
  needLongwave     <- isTRUE(out[["Rlwdown"]]) || isTRUE(out[["Rlwup"]]) || needLeafTemp
  # Air temperature and humidity profiles are therefore also needed whenever leaf
  # temperature is requested, even if those profile variables are not returned.
  needProfile      <- reqhgt_above0 && (isTRUE(out[["Tz"]]) || isTRUE(out[["relhum"]]) || needLeafTemp)
  needCanopyTemp   <- isTRUE(out[["tleaf"]]) || needLongwave || needProfile
  needGroundTemp   <- (reqhgt_eq0 && isTRUE(out[["Tz"]])) || needCanopyTemp ||
    (reqhgt_below0 && isTRUE(out[["Tz"]])) || needLongwave
  needSurfaceSoilm <- needGroundTemp || (reqhgt_eq0 && isTRUE(out[["soilm"]]))
  needBelowSoilm   <- reqhgt_below0 && isTRUE(out[["soilm"]])
  needBelowTemp    <- reqhgt_below0 && isTRUE(out[["Tz"]])
  needRaw <- any(unlist(out[setdiff(names(out), "soilm")]))

  # The reference aerodynamic geometry remains fixed through time; seasonal changes
  # in the fine-grid canopy are handled locally in each day's radiation/wind solve.
  zref <- pointmodel$zref
  if (needRaw) {
    if (is.null(pointmodel$model$uf)) {
      stop(.incompletePointmodelMessage(), call. = FALSE)
    }
    hRef <- pointmodel$vegp$h[1]
    paiRef <- pointmodel$vegp$pai[1]
    dRef <- zeroplanedisCpp(hRef, paiRef)
    zmRef <- roughlengthCpp(hRef, paiRef, dRef)
    sheltarr <- as.vector(.windShelterInterp(hora, weather$winddir))
    # The reference's ground-to-reference resistance as the grid evaluates it,
    # so that its side of the ground heat flux scaling matches a cell's.
    rGzRefGrid <- referenceGroundResistCpp(weather$windspeed, pointmodel$model$uf, pointmodel$model$Hout,
                                           weather$temp, weather$pres, zref, hRef, paiRef)
  } else {
    dRef <- NA_real_; zmRef <- NA_real_
    sheltarr <- rep(NA_real_, rows * cols * tsteps)
    rGzRefGrid <- rep(NA_real_, tsteps)
  }

  # A below-ground request must refer to the same depth solved by the reference point
  # model; whole-run thermal quantities used to scale that soil profile are fixed.
  if (needBelowTemp) {
    if (is.null(pointmodel$model$Tzbelow)) {
      stop(.depthMismatchMessage(reqhgt, pointmodel$reqhgt), call. = FALSE)
    }
    if (!is.null(pointmodel$reqhgt) && !isTRUE(all.equal(pointmodel$reqhgt, reqhgt))) {
      stop(.depthMismatchMessage(reqhgt, pointmodel$reqhgt), call. = FALSE)
    }
  }

  # Only the requested fields with a value at this height are accumulated over
  # the run; any other requested field is returned as a single NA.
  requested <- names(out)[out]
  live <- intersect(requested, .fieldsAtHeight(reqhgt))
  resultFull <- stats::setNames(lapply(live, function(f) naArr3), live)

  fixed <- list(rows = rows, cols = cols, latM = latM, lonM = lonM, reqhgt = reqhgt,
                out = out, needRaw = needRaw, needSurfaceSoilm = needSurfaceSoilm,
                needBelowSoilm = needBelowSoilm, needBelowTemp = needBelowTemp,
                needGroundTemp = needGroundTemp, needCanopyTemp = needCanopyTemp,
                needProfile = needProfile, needLongwave = needLongwave,
                needLeafTemp = needLeafTemp,
                group = group, leafphys = leafphys,
                grefm = grefm, slope = slope, aspect = aspect, svfa = svfa, hora = hora,
                Smin = Smin, Smax = Smax, tadd = tadd, wetw = wetw,
                zref = zref, zrefRepresentative = zref, dRef = dRef, zmRef = zmRef,
                campbellMats = campbellMats,
                refVq = pointmodel$soilc$Vq[1], refVm = pointmodel$soilc$Vm[1],
                refVo = pointmodel$soilc$Vo[1], refMc = pointmodel$soilc$Mc[1],
                refThetaS = pointmodel$soilc$Smax[1], refPsie = pointmodel$soilc$psi_e[1],
                refB = pointmodel$soilc$b[1],
                gsmax = .gsmaxMatrix(vegp, rows, cols),
                paiRef = if (needRaw) paiRef else NA_real_,
                matemp = pointmodel$matemp)

  # Advance one calendar day at a time so the vegetation state is held constant
  # within a day but can change between days. Each day is inserted back into the
  # full hourly output sequence after the common physical calculation.
  for (d in seq_len(nDaysNeeded)) {
    idx <- ((d - 1) * 24 + 1):(d * 24)
    dpos <- dayPos[d]
    shidx <- ((d - 1) * rows * cols * 24 + 1):(d * rows * cols * 24)

    vegpl <- list(
      hgt = resolveFieldForDay("hgt", dpos), pai = resolveFieldForDay("pai", dpos),
      x = resolveFieldForDay("x", dpos), leafr = resolveFieldForDay("leafr", dpos),
      leaft = resolveFieldForDay("leaft", dpos), clump = resolveFieldForDay("clump", dpos),
      Lfrac = resolveFieldForDay("Lfrac", dpos)
    )

    day <- list(
      obstime = obstime[idx, ], climdata = climdata[idx, ], vegpl = vegpl,
      ufRef = if (needRaw) pointmodel$model$uf[idx] else rep(NA_real_, 24),
      HRef = if (needRaw) pointmodel$model$Hout[idx] else rep(NA_real_, 24),
      shelterc = sheltarr[shidx],
      refG = if (needGroundTemp) pointmodel$model$G[idx] else rep(NA_real_, 24),
      refRabsG = if (needGroundTemp) pointmodel$model$RabsGround[idx] else rep(NA_real_, 24),
      refEmG = if (needGroundTemp) pointmodel$model$emGround[idx] else rep(NA_real_, 24),
      refrGz = if (needGroundTemp) rGzRefGrid[idx] else rep(NA_real_, 24),
      reftheta0 = if (needSurfaceSoilm) pointmodel$model$theta0[idx] else rep(NA_real_, 24),
      precip = weather$precip[idx],
      refWet = if (needGroundTemp) pointmodel$model$wetShare[idx] else rep(NA_real_, 24),
      refAvail = if (needGroundTemp) pointmodel$model$filmAvailable[idx] else rep(NA_real_, 24),
      refTground = if (needBelowTemp) pointmodel$model$Tground[idx] else rep(NA_real_, 24),
      refTzbelow = if (needBelowTemp) pointmodel$model$Tzbelow[idx] else rep(NA_real_, 24),
      refthetazbelow = if (needBelowSoilm) pointmodel$model$thetazbelow[idx] else rep(NA_real_, 24)
    )

    dayResult <- .gridmodelDayCore(day, fixed)
    for (fld in live) {
      resultFull[[fld]][, , idx] <- dayResult[[fld]]
    }
  }

  # Restore georeferencing so each model field is returned as a directly usable
  # raster time series on the input terrain grid.
  result <- stats::setNames(lapply(requested, function(f) resultFull[[f]]), requested)
  .gridOutputRasters(result, dtm, live)
}

# Evaluate one day of a seasonal-vegetation grid run. The daily canopy state changes
# radiation interception, aerodynamic structure and foliage above reqhgt, while the
# terrain, soil and PFT identity remain fixed. The same sequence is then followed as
# in the time-invariant grid model: radiation/wind -> soil water -> ground temperature
# -> canopy temperature -> vertical air/leaf state.
.gridmodelDayCore <- function(day, fixed) {
  rows <- fixed$rows; cols <- fixed$cols
  hgt <- day$vegpl$hgt; pai <- day$vegpl$pai; x <- day$vegpl$x
  lref <- day$vegpl$leafr; ltra <- day$vegpl$leaft; clump <- day$vegpl$clump
  tsteps <- 24L
  naArr3 <- array(NA_real_, dim = c(rows, cols, tsteps))

  # Re-evaluate whether foliage is present and how much lies above reqhgt from this
  # day's canopy height and PAI; seasonal leaf-off can therefore change the local regime.
  # Either zero is bare ground (a plant of zero height or plant area has no mass).
  bareground <- !is.na(hgt) & (hgt <= 0 | pai <= 0)
  hgt[bareground] <- 0
  pai[bareground] <- 0
  vegcell <- !is.na(hgt) & !bareground
  paia <- matrix(NA_real_, nrow = rows, ncol = cols)
  paia[bareground] <- 0
  if (any(vegcell)) {
    paia[vegcell] <- .paiaboveheightPFT(fixed$reqhgt, hgt[vegcell], pai[vegcell], fixed$group[vegcell])
  }

  vegpl <- list(hgt = hgt, pai = pai, paia = paia, x = x, leafr = lref, leaft = ltra, clump = clump)
  soilcl <- list(gref = fixed$grefm, slope = fixed$slope, aspect = fixed$aspect,
                 svfa = fixed$svfa, hor = fixed$hora)

  if (fixed$needRaw) {
    raw <- runmicro1Cpp(day$obstime, day$climdata, vegpl, soilcl, fixed$latM, fixed$lonM,
                         fixed$zref, fixed$reqhgt, day$ufRef, day$HRef, fixed$dRef, fixed$zmRef, day$shelterc)
  } else {
    naMat <- matrix(NA_real_, rows, cols)
    raw <- list(radGsw = naArr3, radGlw = naArr3, radCsw = naArr3, radClw = naArr3,
                emGround = naArr3, emCanopy = naArr3, rGh = naArr3,
                rGreq = naArr3, zs = matrix(NA_real_, rows, cols),
                shapeR = matrix(NA_real_, rows, cols), shapeC = matrix(NA_real_, rows, cols),
                Rbdown = naArr3, Rddown = naArr3, Rdup = naArr3,
                uz = naArr3, uf = naArr3, rGz = naArr3, rGm = naArr3, rHa = naArr3, a2 = naArr3, L = naArr3,
                d = naMat, zm = naMat)
  }

  if (fixed$needSurfaceSoilm) {
    raw$soilm <- soilmDistributeCpp(fixed$Smin, fixed$Smax, fixed$tadd, fixed$wetw,day$reftheta0, tsteps)
  } else {
    raw$soilm <- naArr3
  }

  if (fixed$needGroundTemp) {
    sr <- .cellSurfaceResistance(raw, day$climdata$swdown, day$climdata$difrad, day$climdata$temp,
                                 day$climdata$relhum, day$climdata$pres, day$precip, day$obstime$year,
                                 raw$soilm, fixed$campbellMats, hgt, pai, x, lref, ltra, fixed$svfa,
                                 fixed$gsmax, fixed$leafphys, day$vegpl$Lfrac, day$refWet, day$refAvail,
                                 fixed$paiRef)
    RabsGround <- raw$radGsw + raw$radGlw
    gh <- groundHeatFluxCpp(RabsGround, raw$rGz, raw$rGm,
                             raw$radCsw + raw$radClw, raw$rHa, sr$rSurf, sr$hSurf, sr$hFol,
                             raw$soilm,
                             day$climdata$temp, day$climdata$relhum, day$climdata$pres,
                             fixed$campbellMats$Vq, fixed$campbellMats$Vm, fixed$campbellMats$Vo,
                             fixed$campbellMats$Mc, fixed$campbellMats$thetaS, fixed$campbellMats$psie,
                             fixed$campbellMats$b, fixed$svfa, raw$emGround, raw$emCanopy,
                             fixed$slope, fixed$aspect, day$climdata$swdown, day$climdata$difrad,
                             fixed$latM, fixed$lonM,
                             day$obstime$year, day$obstime$month, day$obstime$day, day$obstime$hour,
                             day$refG, day$refRabsG, day$refrGz, day$reftheta0, day$refEmG,
                             fixed$refVq, fixed$refVm, fixed$refVo, fixed$refMc,
                             fixed$refThetaS, fixed$refPsie, fixed$refB)
    raw$G_est <- gh$G_est
    raw$Ts_est <- gh$Ts_est
  } else {
    raw$G_est <- naArr3
    raw$Ts_est <- naArr3
  }

  # The bulk surface: as the ground pass formed it where foliage is present, the
  # one-surface estimate in vegetated cells it leaves, the ground's own
  # temperature on bare ground; see .coupledSurfaceState.
  if (fixed$needGroundTemp) {
    RabsCanopy <- .oneSurfaceRabs(raw$radCsw + raw$radClw, gh, hgt)
    TcOne <- canopyTempCpp(RabsCanopy, raw$rHa, raw$G_est,
                           day$climdata$temp, day$climdata$relhum, day$climdata$pres,
                           fixed$svfa, raw$emCanopy, sr$rSurf, sr$hSurf)
    st <- .coupledSurfaceState(gh, TcOne, sr, hgt)
    raw$Tcanopy_est <- st$Tcanopy_est
    raw$Ts_est <- st$Ts_est
    sr$hSurf <- st$hSurf
  } else {
    raw$Tcanopy_est <- naArr3
  }

  if (!is.na(fixed$reqhgt) && fixed$reqhgt < 0) {
    # ---- BELOW GROUND ----
    Tzbelow_est <- if (fixed$needBelowTemp) {
      DD <- dampingDepthGridCpp(fixed$campbellMats$Vq, fixed$campbellMats$Vm, fixed$campbellMats$Vo,
                                fixed$campbellMats$Mc, fixed$campbellMats$thetaS,
                                fixed$campbellMats$psie, fixed$campbellMats$b,
                                raw$soilm, raw$Ts_est, day$climdata$pres)
      belowGroundShortcutCpp(raw$Ts_est, day$refTground, day$refTzbelow,
                              DD, fixed$matemp, fixed$reqhgt, 8760)
    } else {
      naArr3
    }
    soilm_below <- if (fixed$needBelowSoilm) {
      soilmDistributeCpp(fixed$Smin, fixed$Smax, fixed$tadd, fixed$wetw,day$refthetazbelow, tsteps)
    } else {
      naArr3
    }
    result <- list(Tz = Tzbelow_est, tleaf = naArr3, soilm = soilm_below,
                    relhum = naArr3, windspeed = naArr3,
                    Rdirdown = naArr3, Rdifdown = naArr3, Rswup = naArr3,
                    Rlwdown = naArr3, Rlwup = naArr3)
  } else if (!is.na(fixed$reqhgt) && fixed$reqhgt == 0) {
    # ---- GROUND LEVEL ----
    lw <- if (fixed$needLongwave) {
      longwaveGridCpp(raw$Ts_est, raw$Tcanopy_est, day$climdata$temp, hgt, pai, paia, fixed$svfa, day$climdata$lwdown, fixed$reqhgt)
    } else {
      list(Rlwdown = naArr3, Rlwup = naArr3)
    }
    result <- list(Tz = raw$Ts_est, tleaf = naArr3, soilm = raw$soilm,
                    relhum = naArr3, windspeed = naArr3,
                    Rdirdown = raw$Rbdown, Rdifdown = raw$Rddown, Rswup = raw$Rdup,
                    Rlwdown = lw$Rlwdown, Rlwup = lw$Rlwup)
  } else {
    # ---- reqhgt > 0: per-cell above/below-canopy split ----
    # Air profiles are referenced to one scalar atmospheric measurement height for
    # the day. In gridded-climate mode this is represented by the mean reference
    # height of the contributing coarse cells.
    ac <- if (fixed$needProfile) {
      aboveCanopyProfileGridCpp(raw$Tcanopy_est, raw$rHa, raw$L, raw$uf, raw$d, raw$zh,
                                 day$climdata$temp, day$climdata$relhum, sr$rSurf, sr$hSurf,
                                 raw$rGz, raw$rGreq, raw$zs,
                                 fixed$reqhgt, fixed$zrefRepresentative)
    } else {
      list(Tz = naArr3, RHz = naArr3)
    }

    bc <- if (fixed$needProfile) {
      belowCanopyProfileGridCpp(
        raw$Tcanopy_est, raw$Ts_est, raw$rHa, raw$rGh, raw$rGm, raw$rGz,
        raw$shapeR, raw$shapeC, hgt, pai,
        day$climdata$temp, day$climdata$relhum, day$climdata$pres, sr$rSurf,
        sr$hSurf,
        fixed$reqhgt, fixed$zrefRepresentative)
    } else {
      list(Tz = naArr3, RHz = naArr3)
    }

    above_mask <- rep(as.vector(hgt <= fixed$reqhgt), tsteps)
    Tz_pc <- ifelse(above_mask, as.vector(ac$Tz), as.vector(bc$Tz))
    relhum_pc <- ifelse(above_mask, as.vector(ac$RHz), as.vector(bc$RHz))
    dim(Tz_pc) <- dim(raw$Tcanopy_est)
    dim(relhum_pc) <- dim(raw$Tcanopy_est)

    lw <- if (fixed$needLongwave) {
      longwaveGridCpp(raw$Ts_est, raw$Tcanopy_est, day$climdata$temp, hgt, pai, paia, fixed$svfa, day$climdata$lwdown, fixed$reqhgt)
    } else {
      list(Rlwdown = naArr3, Rlwup = naArr3)
    }

    # Where foliage occupies reqhgt, solve leaf energy balance using that day's local
    # air state, radiation, wind and the cell's fixed PFT physiology.
    tleaf_pc <- if (fixed$needLeafTemp) {
      radLlw <- 0.5 * 0.97 * (lw$Rlwdown + lw$Rlwup) # 0.97 = mc::surfaceEmissivity
      Ca <- rep(Cafromyear(stats::median(day$obstime$year)), tsteps)
      # Leaf forcing uses the local microclimate at reqhgt.
      leafTempCpp(raw$radLpar, raw$radLsw, radLlw, raw$uz,
                  Tz_pc, relhum_pc, day$climdata$pres, Ca,
                  soilWaterPotentialGridCpp(fixed$campbellMats$thetaS, fixed$campbellMats$psie,
                                            fixed$campbellMats$b, raw$soilm),
                  hgt, pai, fixed$leafphys$Vcmax25, fixed$leafphys$Tup, fixed$leafphys$Tlw,
                  fixed$leafphys$Dcrit, fixed$leafphys$alpha, fixed$leafphys$f0, fixed$leafphys$fd,
                  fixed$leafphys$psi50, fixed$leafphys$apsi, fixed$leafphys$rpmin,
                  fixed$leafphys$leafd, fixed$leafphys$isC3, fixed$reqhgt)
    } else {
      naArr3
    }

    result <- list(Tz = Tz_pc, tleaf = tleaf_pc, soilm = naArr3,
                    relhum = relhum_pc, windspeed = raw$uzActual,
                    Rdirdown = raw$Rbdown, Rdifdown = raw$Rddown, Rswup = raw$Rdup,
                    Rlwdown = lw$Rlwdown, Rlwup = lw$Rlwup)
  }
  result
}

# Coarse-gridded climate references with seasonally varying vegetation.
# This is the most spatially and temporally resolved shortcut pathway: atmospheric
# and reference surface state vary across the coarse climate grid, while canopy
# structure can also change by day on the fine terrain grid.
.rungridmodel4 <- function(pointmodela, vegp, soilc, dtm, reqhgt = 0, tfact = 1.7, twiwet = 10,
                           altcorrect = 0, hor = NA, twi = NA, svf = NA,
                           slr = NA, apr = NA, out = NULL, splineveg = FALSE) {
  if (is.null(pointmodela$points) || is.null(pointmodela$dtmc)) {
    stop("`pointmodel` must be the result of runpointmodel() or subsetpointmodel().",
         call. = FALSE)
  }
  points <- pointmodela$points
  dtmc <- .unpackRaster(pointmodela$dtmc)

  # Requested outputs determine which physical stages need to be evaluated.
  .outFieldNames <- .GRID_OUT_FIELDS
  if (is.null(out)) {
    out <- stats::setNames(rep(TRUE, length(.outFieldNames)), .outFieldNames)
  } else {
    if (is.null(names(out)) || any(!nzchar(names(out)))) {
      stop(.outNotNamedMessage(.outFieldNames), call. = FALSE)
    }
    unknown <- setdiff(names(out), .outFieldNames)
    if (length(unknown) > 0) {
      stop(.outUnknownMessage(unknown, .outFieldNames), call. = FALSE)
    }
    full <- stats::setNames(rep(FALSE, length(.outFieldNames)), .outFieldNames)
    full[names(out)] <- as.logical(out)
    out <- full
  }

  if (inherits(dtm, "PackedSpatRaster")) dtm <- terra::rast(dtm)
  if (isTRUE(all.equal(dim(dtm)[1:2], dim(dtmc)[1:2]))) {
    stop(.sameResolutionMessage(), call. = FALSE)
  }

  isValid <- vapply(points, is.list, logical(1))
  validCells <- which(isValid)
  if (length(validCells) == 0) {
    stop(.noUsableRunMessage(), call. = FALSE)
  }

  # Gather the atmospheric forcing, solved surface state, aerodynamic reference
  # geometry and reference soil properties from every valid coarse climate cell.
  weatherList <- vector("list", length(points))
  modelList <- vector("list", length(points))
  matemp <- rep(NA_real_, length(points))
  zrefV <- rep(NA_real_, length(points))
  hRefV <- rep(NA_real_, length(points))
  paiRefV <- rep(NA_real_, length(points))
  VqRefV <- rep(NA_real_, length(points))
  VmRefV <- rep(NA_real_, length(points))
  VoRefV <- rep(NA_real_, length(points))
  McRefV <- rep(NA_real_, length(points))
  thetaSRefV <- rep(NA_real_, length(points))
  psieRefV <- rep(NA_real_, length(points))
  bRefV <- rep(NA_real_, length(points))
  for (i in validCells) {
    weatherList[[i]] <- points[[i]]$weather
    modelList[[i]] <- points[[i]]$model
    matemp[i] <- points[[i]]$matemp
    zrefV[i] <- points[[i]]$zref
    hRefV[i] <- points[[i]]$vegp$h
    paiRefV[i] <- points[[i]]$vegp$pai
    VqRefV[i] <- points[[i]]$soilc$Vq[1]
    VmRefV[i] <- points[[i]]$soilc$Vm[1]
    VoRefV[i] <- points[[i]]$soilc$Vo[1]
    McRefV[i] <- points[[i]]$soilc$Mc[1]
    thetaSRefV[i] <- points[[i]]$soilc$Smax[1]
    psieRefV[i] <- points[[i]]$soilc$psi_e[1]
    bRefV[i] <- points[[i]]$soilc$b[1]
    # The reference's ground-to-reference resistance as the grid evaluates it,
    # so that its side of the ground heat flux scaling matches a cell's.
    modelList[[i]]$rGzGrid <- referenceGroundResistCpp(points[[i]]$weather$windspeed,
                                                       points[[i]]$model$uf, points[[i]]$model$Hout,
                                                       points[[i]]$weather$temp, points[[i]]$weather$pres,
                                                       points[[i]]$zref,
                                                       points[[i]]$vegp$h, points[[i]]$vegp$pai)
  }

  # Map modelled days onto the full seasonal sequence represented by the vegetation
  # layers; all coarse references share the same calendar.
  refPoint <- points[[validCells[1]]]
  nDaysNeeded <- nrow(refPoint$weather) %/% 24L
  if (nrow(refPoint$weather) %% 24 != 0) {
    stop(.notWholeDaysMessage(nrow(refPoint$weather)), call. = FALSE)
  }
  if (!is.null(refPoint$nhoursFull)) {
    nDaysFull <- ceiling(refPoint$nhoursFull / 24)
  } else if (!is.null(refPoint$subs)) {
    nDaysFull <- ceiling(max(refPoint$subs) / 24)
  } else {
    nDaysFull <- nDaysNeeded
  }
  dayPos <- integer(nDaysNeeded)
  for (d in seq_len(nDaysNeeded)) {
    if (!is.null(refPoint$subs)) {
      dayPos[d] <- ceiling(refPoint$subs[(d - 1) * 24 + 1] / 24)
    } else {
      dayPos[d] <- d
    }
  }
  tsteps <- nDaysNeeded * 24L

  # Interpolate coarse atmospheric forcing to the fine terrain grid and optionally
  # adjust temperature and pressure for local elevation. Radiation interpolates
  # directly, while wind is interpolated as vector components.
  nT <- nrow(weatherList[[validCells[1]]])
  tcFine <- .coarseListToFineArray(weatherList, "temp", nT, dtmc, dtm)
  relhumFine <- .coarseListToFineArray(weatherList, "relhum", nT, dtmc, dtm)
  pkCoarse <- .coarseListToArray(weatherList, "pres", nT, dtmc)
  alt <- .altitudeCorrectPkTc(pkCoarse, tcFine, relhumFine, dtmc, dtm, altcorrect)
  swdownFine <- .coarseListToFineArray(weatherList, "swdown", nT, dtmc, dtm)
  difradFine <- .coarseListToFineArray(weatherList, "difrad", nT, dtmc, dtm)
  lwdownFine <- .coarseListToFineArray(weatherList, "lwdown", nT, dtmc, dtm)
  wind <- .coarseWindToFine(weatherList, nT, dtmc, dtm)
  climdataFineList <- list(temp = alt$tc, swdown = swdownFine, difrad = difradFine,
                            lwdown = lwdownFine, windspeed = wind$windspeed,
                            relhum = relhumFine, pres = alt$pk)

  tme <- as.POSIXlt(weatherList[[validCells[1]]]$obs_time, tz = "UTC")
  obstime <- data.frame(year = tme$year + 1900, month = tme$mon + 1, day = tme$mday,
                         hour = tme$hour + tme$min / 60 + tme$sec / 3600)

  rows <- dim(dtm)[1]; cols <- dim(dtm)[2]
  tomat <- function(v) {
    if (inherits(v, "PackedSpatRaster")) v <- terra::rast(v)
    if (inherits(v, "SpatRaster")) return(terra::as.matrix(v[[1]], wide = TRUE))
    if (is.matrix(v)) return(v)
    matrix(v, nrow = rows, ncol = cols)
  }

  # Preserve fine-cell geographic position for solar geometry across the domain.
  llFine <- .latlonsFromCells(dtm)
  latR <- dtm; terra::values(latR) <- llFine$lat
  lonR <- dtm; terra::values(lonR) <- llFine$lon
  latM <- tomat(latR)
  lonM <- tomat(lonR)

  # Resolve each vegetation property from fixed or seasonal layers independently,
  # using discrete day assignment or smooth interpolation as requested.
  vegFieldNames <- c("hgt", "pai", "x", "leafr", "leaft", "clump", "Lfrac")
  nLayersV <- stats::setNames(vapply(vegFieldNames, function(f) .nLayersOf(vegp[[f]]), integer(1)), vegFieldNames)

  fieldSetup <- list()
  for (f in vegFieldNames) {
    nl <- nLayersV[[f]]
    v <- vegp[[f]]
    if (nl == 1L) {
      fieldSetup[[f]] <- list(invariant = TRUE, matrix = tomat(v))
    } else if (!splineveg) {
      dayLyr <- ceiling(seq_len(nDaysFull) * nl / nDaysFull)
      fieldSetup[[f]] <- list(invariant = FALSE, raster = v, dayLyr = dayLyr)
    } else {
      knotX <- (seq_len(nl) - 0.5) * nDaysFull / nl
      knotY <- vapply(seq_len(nl), function(i) as.vector(.tomatLayer(v, i, rows, cols)), numeric(rows * cols))
      fit <- splineFitCpp(knotX, knotY)
      fieldSetup[[f]] <- list(invariant = FALSE, fit = fit)
    }
  }

  # Obtain each day's vegetation state before converting horizontal-area PAI to
  # the slope-aware inclined-area convention where required.
  resolveFieldForDay <- function(f, dpos) {
    fs <- fieldSetup[[f]]
    m <- if (fs$invariant) {
      fs$matrix
    } else if (!splineveg) {
      .tomatLayer(vegp[[f]], fs$dayLyr[dpos], rows, cols)
    } else {
      matrix(splineEvalCpp(fs$fit, dpos), nrow = rows, ncol = cols)
    }
    if (f == "pai" && isTRUE(pointmodela$paiFlat)) m <- m * cos(slope * pi / 180)
    m
  }

  # PFT and physiology represent the persistent vegetation type, taken from each
  # cell's habitat, while canopy state varies by day; whether a cell is vegetated
  # is judged from its middle-season structure.
  midHgt <- .tomatLayer(vegp$hgt, ceiling(nLayersV[["hgt"]] / 2), rows, cols)
  midPai <- .tomatLayer(vegp$pai, ceiling(nLayersV[["pai"]] / 2), rows, cols)
  midX   <- .tomatLayer(vegp$x,   ceiling(nLayersV[["x"]] / 2),   rows, cols)
  midBareground <- !is.na(midHgt) & (midHgt <= 0 | midPai <= 0)
  midVegcell <- !is.na(midHgt) & !midBareground
  group <- matrix(NA_character_, nrow = rows, ncol = cols)
  pftM <- matrix(NA_character_, nrow = rows, ncol = cols)
  if (any(midVegcell)) {
    pft <- habitattoPFT(as.vector(tomat(vegp$habitat))[midVegcell], as.vector(latM)[midVegcell])
    group[midVegcell] <- .PFT_GROUP_OF[pft]
    pftM[midVegcell] <- pft
  }
  # PFT-specific physiology and leaf geometry remain fixed while seasonal canopy
  # height, area and optical properties vary through the year.
  leafphys <- .leafPhysiologyMatrices(pftM, midVegcell, rows, cols, vegp)

  # Soil properties and terrain geometry are fixed landscape attributes and are
  # resolved once before the daily vegetation cycle is evaluated.
  grefm <- tomat(soilc$groundr)
  sminsmax <- .soilcSminSmax(soilc$soiltype)
  Smin <- tomat(sminsmax$Smin)
  Smax <- tomat(sminsmax$Smax)
  wetness <- .topoWetness(dtm, twi, tfact, twiwet)
  tadd <- tomat(wetness$shift)
  wetw <- tomat(wetness$wet)
  sky <- .horizonAndSkyview(dtm, hor, svf)
  svfa <- tomat(sky$svf)
  hora <- sky$hora
  if (is.logical(slr) || is.logical(apr)) {
    fill_ext <- terravars:::.fillna_dtm(dtm)
  }
  slope <- if (is.logical(slr)) {
    terra::as.matrix(terra::crop(terra::terrain(fill_ext$extended, "slope", unit = "degrees"), dtm), wide = TRUE)
  } else {
    tomat(slr)
  }
  aspect <- if (is.logical(apr)) {
    terra::as.matrix(terra::crop(terra::terrain(fill_ext$extended, "aspect", unit = "degrees"), dtm), wide = TRUE)
  } else {
    tomat(apr)
  }
  campbell <- .soilcCampbellParams(soilc$soiltype)
  campbellMats <- list(Vq = tomat(campbell$Vq), Vm = tomat(campbell$Vm), Vo = tomat(campbell$Vo),
                       Mc = tomat(campbell$Mc), thetaS = tomat(campbell$thetaS),
                       psie = tomat(campbell$psie), b = tomat(campbell$b))

  # Determine once which physical quantities are required by the selected outputs
  # and requested height/depth.
  reqhgt_below0 <- !is.na(reqhgt) && reqhgt < 0
  reqhgt_eq0    <- !is.na(reqhgt) && reqhgt == 0
  reqhgt_above0 <- !is.na(reqhgt) && reqhgt > 0
  # Leaf temperature is meaningful only above ground and requires both the local
  # air profile and longwave radiation at the requested height for its energy balance.
  needLeafTemp     <- reqhgt_above0 && isTRUE(out[["tleaf"]])
  needLongwave     <- isTRUE(out[["Rlwdown"]]) || isTRUE(out[["Rlwup"]]) || needLeafTemp
  # Air temperature and humidity profiles are therefore also needed whenever leaf
  # temperature is requested, even if those profile variables are not returned.
  needProfile      <- reqhgt_above0 && (isTRUE(out[["Tz"]]) || isTRUE(out[["relhum"]]) || needLeafTemp)
  needCanopyTemp   <- isTRUE(out[["tleaf"]]) || needLongwave || needProfile
  needGroundTemp   <- (reqhgt_eq0 && isTRUE(out[["Tz"]])) || needCanopyTemp ||
    (reqhgt_below0 && isTRUE(out[["Tz"]])) || needLongwave
  needSurfaceSoilm <- needGroundTemp || (reqhgt_eq0 && isTRUE(out[["soilm"]]))
  needBelowSoilm   <- reqhgt_below0 && isTRUE(out[["soilm"]])
  needBelowTemp    <- reqhgt_below0 && isTRUE(out[["Tz"]])
  needRaw <- any(unlist(out[setdiff(names(out), "soilm")]))

  # Interpolate the reference aerodynamic state from coarse climate cells to the
  # fine grid; local seasonal canopy effects are applied subsequently day by day.
  zrefFine <- .coarseScalarToFineMatrix(zrefV, dtmc, dtm)
  if (needRaw) {
    hRefFine <- .coarseScalarToFineMatrix(hRefV, dtmc, dtm)
    paiRefFine <- .coarseScalarToFineMatrix(paiRefV, dtmc, dtm)
    dRefFine <- zeroplanedisCpp(hRefFine, paiRefFine)
    zmRefFine <- roughlengthCpp(hRefFine, paiRefFine, dRefFine)
    ufRefFine <- .coarseListToFineArray(modelList, "uf", nT, dtmc, dtm)
    HRefFine <- .coarseListToFineArray(modelList, "Hout", nT, dtmc, dtm)
    zrefRepresentative <- mean(zrefV, na.rm = TRUE)
    sheltarrFine <- as.vector(.windShelterInterp(hora, wind$winddir))
  } else {
    dRefFine <- matrix(NA_real_, rows, cols); zmRefFine <- matrix(NA_real_, rows, cols)
    ufRefFine <- array(NA_real_, dim = c(rows, cols, nT))
    HRefFine <- ufRefFine
    sheltarrFine <- rep(NA_real_, rows * cols * nT)
    zrefRepresentative <- NA_real_
  }

  # Interpolate the coarse reference surface and soil state needed to anchor local
  # ground and canopy energy-balance estimates throughout the run.
  if (needGroundTemp) {
    GFine <- .coarseListToFineArray(modelList, "G", nT, dtmc, dtm)
    RabsGFine <- .coarseListToFineArray(modelList, "RabsGround", nT, dtmc, dtm)
    rGzFine <- .coarseListToFineArray(modelList, "rGzGrid", nT, dtmc, dtm)
    emGFine <- .coarseListToFineArray(modelList, "emGround", nT, dtmc, dtm)
    thetaSRefFine <- .coarseScalarToFineMatrix(thetaSRefV, dtmc, dtm)
    psieRefFine <- .coarseScalarToFineMatrix(psieRefV, dtmc, dtm)
    bRefFine <- .coarseScalarToFineMatrix(bRefV, dtmc, dtm)
    VqRefFine <- .coarseScalarToFineMatrix(VqRefV, dtmc, dtm)
    VmRefFine <- .coarseScalarToFineMatrix(VmRefV, dtmc, dtm)
    VoRefFine <- .coarseScalarToFineMatrix(VoRefV, dtmc, dtm)
    McRefFine <- .coarseScalarToFineMatrix(McRefV, dtmc, dtm)
  } else {
    GFine <- array(NA_real_, dim = c(rows, cols, nT)); RabsGFine <- GFine; rGzFine <- GFine
    emGFine <- GFine
    VqRefFine <- matrix(NA_real_, rows, cols); VmRefFine <- VqRefFine
    thetaSRefFine <- VqRefFine; psieRefFine <- VqRefFine; bRefFine <- VqRefFine
    VoRefFine <- VqRefFine; McRefFine <- VqRefFine
  }
  if (needSurfaceSoilm) {
    theta0Fine <- .coarseListToFineArray(modelList, "theta0", nT, dtmc, dtm)
  } else {
    theta0Fine <- array(NA_real_, dim = c(rows, cols, nT))
  }
  if (needGroundTemp) {
    wetFine <- .coarseListToFineArray(modelList, "wetShare", nT, dtmc, dtm)
    availFine <- .coarseListToFineArray(modelList, "filmAvailable", nT, dtmc, dtm)
    precipFine <- .coarseListToFineArray(weatherList, "precip", nT, dtmc, dtm)
  } else {
    wetFine <- array(NA_real_, dim = c(rows, cols, nT))
    availFine <- wetFine; precipFine <- wetFine
  }

  # For below-ground output, interpolate the reference soil profiles and their
  # damping/mean-temperature context to the fine grid at the requested depth.
  if (needBelowTemp) {
    if (is.null(modelList[[validCells[1]]]$Tzbelow)) {
      stop(.depthMismatchMessage(reqhgt, refPoint$reqhgt), call. = FALSE)
    }
    if (!is.null(refPoint$reqhgt) && !isTRUE(all.equal(refPoint$reqhgt, reqhgt))) {
      stop(.depthMismatchMessage(reqhgt, refPoint$reqhgt), call. = FALSE)
    }
    TgroundFine <- .coarseListToFineArray(modelList, "Tground", nT, dtmc, dtm)
    TzbelowFine <- .coarseListToFineArray(modelList, "Tzbelow", nT, dtmc, dtm)
    matFine <- .coarseScalarToFineMatrix(matemp, dtmc, dtm) + .meanAltitudeOffset(alt$tc, tcFine)
  } else {
    TgroundFine <- array(NA_real_, dim = c(rows, cols, nT)); TzbelowFine <- TgroundFine
    matFine <- matrix(NA_real_, rows, cols)
  }
  if (needBelowSoilm) {
    thetazbelowFine <- .coarseListToFineArray(modelList, "thetazbelow", nT, dtmc, dtm)
  } else {
    thetazbelowFine <- array(NA_real_, dim = c(rows, cols, nT))
  }

  # Only the requested fields with a value at this height are accumulated over
  # the run; any other requested field is returned as a single NA.
  requested <- names(out)[out]
  live <- intersect(requested, .fieldsAtHeight(reqhgt))
  resultFull <- stats::setNames(
    lapply(live, function(f) array(NA_real_, dim = c(rows, cols, tsteps))), live)

  fixed <- list(rows = rows, cols = cols, latM = latM, lonM = lonM, reqhgt = reqhgt,
                out = out, needRaw = needRaw, needSurfaceSoilm = needSurfaceSoilm,
                needBelowSoilm = needBelowSoilm, needBelowTemp = needBelowTemp,
                needGroundTemp = needGroundTemp, needCanopyTemp = needCanopyTemp,
                needProfile = needProfile, needLongwave = needLongwave,
                needLeafTemp = needLeafTemp,
                group = group, leafphys = leafphys,
                grefm = grefm, slope = slope, aspect = aspect, svfa = svfa, hora = hora,
                Smin = Smin, Smax = Smax, tadd = tadd, wetw = wetw,
                zref = as.vector(zrefFine), zrefRepresentative = zrefRepresentative,
                dRef = as.vector(dRefFine), zmRef = as.vector(zmRefFine),
                campbellMats = campbellMats,
                refVq = as.vector(VqRefFine), refVm = as.vector(VmRefFine),
                refVo = as.vector(VoRefFine), refMc = as.vector(McRefFine),
                refThetaS = as.vector(thetaSRefFine), refPsie = as.vector(psieRefFine),
                refB = as.vector(bRefFine),
                gsmax = .gsmaxMatrix(vegp, rows, cols),
                paiRef = if (needRaw) as.vector(paiRefFine) else NA_real_,
                matemp = as.vector(matFine))

  # For each day, combine that day's fine-grid canopy state with the corresponding
  # fine-grid atmospheric/reference state, then apply the common daily physics.
  for (d in seq_len(nDaysNeeded)) {
    idx <- ((d - 1) * 24 + 1):(d * 24)
    dpos <- dayPos[d]
    shidx <- ((d - 1) * rows * cols * 24 + 1):(d * rows * cols * 24)

    vegpl <- list(
      hgt = resolveFieldForDay("hgt", dpos), pai = resolveFieldForDay("pai", dpos),
      x = resolveFieldForDay("x", dpos), leafr = resolveFieldForDay("leafr", dpos),
      leaft = resolveFieldForDay("leaft", dpos), clump = resolveFieldForDay("clump", dpos),
      Lfrac = resolveFieldForDay("Lfrac", dpos)
    )

    day <- list(
      obstime = obstime[idx, ],
      climdata = data.frame(temp = as.vector(climdataFineList$temp)[shidx],
                             swdown = as.vector(climdataFineList$swdown)[shidx],
                             difrad = as.vector(climdataFineList$difrad)[shidx],
                             lwdown = as.vector(climdataFineList$lwdown)[shidx],
                             windspeed = as.vector(climdataFineList$windspeed)[shidx],
                             relhum = as.vector(climdataFineList$relhum)[shidx],
                             pres = as.vector(climdataFineList$pres)[shidx]),
      vegpl = vegpl,
      ufRef = as.vector(ufRefFine)[shidx],
      HRef = as.vector(HRefFine)[shidx],
      shelterc = sheltarrFine[shidx],
      refG = as.vector(GFine)[shidx],
      refRabsG = as.vector(RabsGFine)[shidx],
      refrGz = as.vector(rGzFine)[shidx],
      refEmG = as.vector(emGFine)[shidx],
      reftheta0 = as.vector(theta0Fine)[shidx],
      precip = as.vector(precipFine)[shidx],
      refWet = as.vector(wetFine)[shidx],
      refAvail = as.vector(availFine)[shidx],
      refTground = as.vector(TgroundFine)[shidx],
      refTzbelow = as.vector(TzbelowFine)[shidx],
      refthetazbelow = as.vector(thetazbelowFine)[shidx]
    )

    dayResult <- .gridmodelDayCore(day, fixed)
    for (fld in live) {
      resultFull[[fld]][, , idx] <- dayResult[[fld]]
    }
  }

  # Restore georeferencing so each model field is returned as a directly usable
  # raster time series on the input terrain grid.
  result <- stats::setNames(lapply(requested, function(f) resultFull[[f]]), requested)
  .gridOutputRasters(result, dtm, live)
}


# Full per-cell microclimate solve. This pathway is used when the weather grid
# already matches the terrain resolution: each land cell receives its own weather,
# vegetation, soil and terrain and solves an independent coupled point-model energy
# and water balance. Terrain slope/aspect, sky exposure, horizon shadowing and wind
# shelter are represented locally. Soil columns do not exchange water laterally, so
# no additional topographic-wetness redistribution is imposed on their solved state.


#' Solve the full point model independently for every grid cell
#'
#' Runs a complete coupled point-model energy and water balance at every land cell
#' when the weather data already have the same spatial resolution as the terrain.
#' Unlike \code{\link{rungridmodel}}, no coarse or shared point solution is used as
#' a reference: each cell experiences its own weather, vegetation, soil and terrain.
#'
#' @details
#' Each cell independently resolves canopy and ground radiation, aerodynamic
#' exchange, stomatal conductance, soil heat and water, and canopy interception.
#' Above- and below-canopy air conditions and below-ground temperatures therefore
#' derive from that cell's own solved state. Terrain slope, aspect, sky exposure and
#' wind shelter are local to the cell. Lateral redistribution of soil water through
#' topographic wetness is not included because each soil column is solved
#' independently.
#'
#' Terrain-horizon shadowing is applied to direct solar radiation. Vegetation is
#' time-invariant within a run; seasonal vegetation layers are not supported here.
#' If a vegetation field has multiple layers, only its middle layer is used,
#' silently.
#'
#' A \code{reqhgt} at or below ground (\code{0} = the surface) draws on each
#' cell's own solved soil profile: \code{Tz} and \code{soilm} are the soil
#' temperature and moisture at that depth, the ground surface's at
#' \code{reqhgt = 0}.
#'
#' @param weather Named raster time series of temperature, humidity, pressure,
#'   shortwave/diffuse/longwave radiation, wind speed/direction and precipitation.
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
#' @param tme A POSIXlt/POSIXct time vector (UTC) corresponding to the time of
#'   each layer.
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
#' @param cores Controls parallel processing across cells: \code{off} (default)
#'   runs sequentially; \code{auto} uses one fewer than the number of
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
#' @param hor an optional array of horizon angles (degrees) in 24 directions.
#'   Calculated automatically from \code{dtm} if not supplied, but can also be
#'   provided as separate inputs to avoid edge effects.
#' @param svf optional raster object of sky view factor. Calculated
#'   automatically from \code{dtm} if not supplied, but can also be provided as
#'   separate inputs to avoid edge effects.
#' @param slr,apr slope and aspect in degrees. Calculated automatically from
#'   \code{dtm} if not supplied, but can also be provided as separate inputs to
#'   avoid edge effects.
#' @param out Optional named vector of logicals indicating which variables to
#'   return, e.g. \code{c(Tz = TRUE, soilm = TRUE)}.
#' @param saveout Logical; when set to \code{TRUE} the results are saved to
#'   \code{filename}.
#' @param filename File path used when \code{saveout = TRUE}.
#'
#' @return A named list of \code{SpatRaster} time series, on the \code{dtm}
#'   grid, for the requested subset of: \code{Tz} (air temperature above
#'   ground, ground-surface temperature at \code{reqhgt = 0}, soil temperature
#'   below ground, deg C), \code{tleaf} (leaf temperature, deg C), \code{soilm}
#'   (volumetric soil moisture: at the requested depth, or at the surface when
#'   \code{reqhgt} is above ground), \code{relhum} (relative humidity,
#'   percent), \code{windspeed} (m/s), \code{Rdirdown}, \code{Rdifdown},
#'   \code{Rswup}, \code{Rlwdown} and \code{Rlwup} (radiation flux densities,
#'   W/m2). A quantity with no physical meaning at the requested height/depth
#'   is returned as a single \code{NA}: everything but \code{Tz} and
#'   \code{soilm} below ground, and \code{tleaf}, \code{relhum} and
#'   \code{windspeed} at the ground surface.
#' @export
runpointmodelasgrid <- function(weather, reqhgt = 0, dtm, vegp, soilc, tme,
                                 paiFlat = FALSE, runchecks = TRUE, zin = 2, uzin = zin,
                                 maxIter = 100, tolerance = 0.1, nlayers = 7, totalDepth = 2,
                                 surface_organicmu = 3, FreeDrain = TRUE, matemp = NA_real_,
                                 pooling = FALSE, cores = "off", soilinit = NULL,
                                 hor = NA, svf = NA, slr = NA, apr = NA, out = NULL,
                                 saveout = FALSE, filename = NULL) {
  .checkSoilinit(soilinit)

  if (isTRUE(saveout) && is.null(filename)) {
    stop("runpointmodelasgrid(): saveout = TRUE requires `filename`", call. = FALSE)
  }

  dtm <- .unpackRaster(dtm)
  vegp <- lapply(vegp, .unpackRaster)
  soilc <- lapply(soilc, .unpackRaster)
  weather <- lapply(weather, .unpackRaster)

  climvars <- c("temp", "relhum", "pres", "swdown", "difrad", "lwdown",
                "windspeed", "winddir", "precip")
  missingvars <- setdiff(climvars, names(weather))
  if (length(missingvars) > 0) {
    stop("runpointmodelasgrid(): weather is missing: ", paste(missingvars, collapse = ", "),
         call. = FALSE)
  }
  nT <- terra::nlyr(weather[[1]])
  if (length(tme) != nT) {
    stop("runpointmodelasgrid(): length(tme) must match the number of layers in weather (",
         nT, ")", call. = FALSE)
  }

  if (inherits(dtm, "PackedSpatRaster")) dtm <- terra::rast(dtm)
  if (!isTRUE(all.equal(dim(dtm)[1:2], dim(weather[[1]])[1:2]))) {
    stop("runpointmodelasgrid(): weather and dtm must share the same spatial ",
         "resolution -- this function is for the case where climate array ",
         "data already matches dtm's own resolution (no coarse/fine split). ",
         "If they differ, use runpointmodel() (array mode) with rungridmodel() instead.",
         call. = FALSE)
  }
  rows <- dim(dtm)[1]; cols <- dim(dtm)[2]

  if (runchecks) {
    checkinputs(.climarrayrDomainMeanWeather(weather, tme), vegp, soilc, dtm, uzin = uzin)
  }
  vegp <- .resolveVegp(vegp, dtm)

  # Put every climate field onto the terrain coordinate system and fill isolated
  # missing climate cells over land so each land cell has a complete forcing series.
  for (v in names(weather)) {
    if (!identical(terra::crs(weather[[v]]), terra::crs(dtm))) {
      weather[[v]] <- terra::project(weather[[v]], terra::crs(dtm))
    }
  }
  landMaskMat <- terra::as.matrix(dtm, wide = TRUE)
  weatherArr <- lapply(weather, function(r) {
    arr <- terra::as.array(r)
    .fillClimArrayNA(arr, landMaskMat)
  })
  # Extract the fixed vegetation and soil state for each independently solved
  # cell; a seasonal field gives its middle layer, as the reference run takes.
  tomat <- function(v) {
    if (inherits(v, "PackedSpatRaster")) v <- terra::rast(v)
    if (inherits(v, "SpatRaster")) return(terra::as.matrix(v[[ceiling(terra::nlyr(v) / 2)]], wide = TRUE))
    if (is.matrix(v)) return(v)
    matrix(v, nrow = rows, ncol = cols)
  }
  hgt <- tomat(vegp$hgt)
  pai <- tomat(vegp$pai)
  xv <- tomat(vegp$x)
  lref <- tomat(vegp$leafr)
  ltra <- tomat(vegp$leaft)
  clumpm <- tomat(vegp$clump)
  grefm <- tomat(soilc$groundr)
  soiltypeM <- tomat(soilc$soiltype)
  # Where habitat is supplied, use it consistently for the PFT identity that controls
  # both physiology and vertical canopy structure.
  modhabitatM <- tomat(vegp$habitat)
  physM <- lapply(.physOverrides(vegp), function(v) matrix(v, nrow = rows, ncol = cols, byrow = TRUE))
  LfracM <- tomat(vegp$Lfrac)
  gsmaxM <- if (is.null(vegp$gsmax)) matrix(NA_real_, rows, cols) else tomat(vegp$gsmax)

  llCells <- .latlonsFromCells(dtm) # per-cell lat/lon, rows*cols length, terra cell-number order
  # Preserve each cell's own geographic position for solar geometry.
  latM <- matrix(llCells$lat, nrow = rows, ncol = cols, byrow = TRUE)
  lonM <- matrix(llCells$lon, nrow = rows, ncol = cols, byrow = TRUE)

  # Resolve local slope, aspect and sky exposure used by each cell's energy balance.
  if (is.logical(slr) || is.logical(apr)) {
    fill_ext <- terravars:::.fillna_dtm(dtm)
  }
  slope <- if (is.logical(slr)) {
    terra::as.matrix(terra::crop(terra::terrain(fill_ext$extended, "slope", unit = "degrees"), dtm), wide = TRUE)
  } else {
    tomat(slr)
  }
  aspect <- if (is.logical(apr)) {
    terra::as.matrix(terra::crop(terra::terrain(fill_ext$extended, "aspect", unit = "degrees"), dtm), wide = TRUE)
  } else {
    tomat(apr)
  }
  horR <- if (is.logical(hor)) .horizonRaster(dtm) else hor
  svfa <- tomat(if (is.logical(svf)) terravars::skyview(dtm, hor = horR) else svf)
  hora <- terra::as.array(horR)

  # Reduce each cell's wind according to the surrounding terrain in the current
  # wind direction, retaining the full temporal variation in directional shelter.
  shelterArr <- .windShelterInterp(hora, weatherArr$winddir)

  # Resolve PFT for vegetated cells from their habitat. This identity controls
  # physiology and within-canopy foliage profiles; bare ground has no PFT.
  # Either zero is bare ground (a plant of zero height or plant area has no mass).
  bareground <- !is.na(hgt) & (hgt <= 0 | pai <= 0)
  hgt[bareground] <- 0
  pai[bareground] <- 0
  vegcell <- !is.na(hgt) & !bareground
  group <- matrix(NA_character_, nrow = rows, ncol = cols)
  pftM <- matrix(NA_character_, nrow = rows, ncol = cols)
  if (any(vegcell)) {
    pft <- habitattoPFT(as.vector(tomat(vegp$habitat))[vegcell], latM[vegcell])
    group[vegcell] <- .PFT_GROUP_OF[pft]
    pftM[vegcell] <- pft
  }
  # PFT-specific physiology and leaf geometry provide the parameters required
  # for the leaf energy balance in vegetated cells.
  leafphys <- .leafPhysiologyMatrices(pftM, vegcell, rows, cols, vegp)

  # Convert PAI from horizontal ground area to the inclined surface convention used
  # by the slope-aware radiation calculation, after PFT classification.
  if (isTRUE(paiFlat)) {
    pai <- pai * cos(slope * pi / 180)
  }

  # Derive each cell's aerodynamic displacement and roughness from its own canopy
  # structure for the above- and within-canopy wind/temperature profiles.
  landcell <- !is.na(hgt)
  dMat <- matrix(NA_real_, rows, cols)
  zmMat <- matrix(NA_real_, rows, cols)
  if (any(landcell)) {
    dMat[landcell] <- zeroplanedisCpp(hgt[landcell], pai[landcell])
    zmMat[landcell] <- roughlengthCpp(hgt[landcell], pai[landcell], dMat[landcell])
  }

  # Map each soil code to the physical soil class used by the independent column solve.
  soiltab <- microclimf::soilparamstable
  soiltypeNameV <- vapply(as.vector(soiltypeM), function(code) {
    if (is.na(code)) return(NA_character_)
    s <- which(soiltab$Number == code)
    if (length(s) == 0) return(NA_character_)
    as.character(soiltab$Soil.type[s[1]])
  }, character(1))
  soiltypeNameM <- matrix(soiltypeNameV, nrow = rows, ncol = cols)

  # Only cells with both land and a recognised soil type can support a complete
  # point-model energy and water balance.
  validMask <- landcell & !is.na(soiltypeNameM)
  validIJ <- which(validMask, arr.ind = TRUE)
  if (nrow(validIJ) == 0) {
    stop("runpointmodelasgrid(): no cell has both vegetation and a recognised soil type: ",
         "check that vegp and soilc cover dtm.", call. = FALSE)
  }

  # The full cell solve is performed once; `out` controls which microclimate
  # fields are assembled and returned.
  .outFieldNames <- .GRID_OUT_FIELDS
  if (is.null(out)) {
    out <- stats::setNames(rep(TRUE, length(.outFieldNames)), .outFieldNames)
  } else {
    if (is.null(names(out)) || any(!nzchar(names(out)))) {
      stop("runpointmodelasgrid(): ", .outNotNamedMessage(.outFieldNames), call. = FALSE)
    }
    unknown <- setdiff(names(out), .outFieldNames)
    if (length(unknown) > 0) {
      stop("runpointmodelasgrid(): ", .outUnknownMessage(unknown, .outFieldNames), call. = FALSE)
    }
    full <- stats::setNames(rep(FALSE, length(.outFieldNames)), .outFieldNames)
    full[names(out)] <- as.logical(out)
    out <- full
  }
  # Only the requested fields with a value at this height are kept from each
  # cell; any other requested field is returned as a single NA. Soil moisture
  # is each cell's own solved value at every height: at the surface when reqhgt
  # is above ground.
  requested <- names(out)[out]
  live <- intersect(requested, union(.fieldsAtHeight(reqhgt), "soilm"))

  tmeUTC <- as.POSIXct(tme, tz = "UTC")
  # Calendar time together with each cell's latitude/longitude determines solar
  # position and therefore both radiation geometry and terrain shadowing.
  tmeLt <- as.POSIXlt(tmeUTC, tz = "UTC")
  solYear <- tmeLt$year + 1900
  solMonth <- tmeLt$mon + 1
  solDay <- tmeLt$mday
  solHour <- tmeLt$hour + tmeLt$min / 60 + tmeLt$sec / 3600
  # The shadow test compares horizon angle with solar altitude; retain the original
  # horizon angles separately because wind shelter uses their angular values directly.
  horaTan <- tan(hora * pi / 180)

  # Solve each valid land cell independently using its own weather, terrain,
  # vegetation and soil. Cells can therefore be evaluated in parallel without
  # exchanging energy or water laterally.
  runOneCell <- function(k) {
    i <- validIJ[k, 1]; j <- validIJ[k, 2]

    wx <- as.data.frame(lapply(weatherArr, function(a) a[i, j, ]))
    if (anyNA(wx)) return(NA) # unfillable gap at this cell
    wx$obs_time <- tmeUTC

    # Direct sunlight is removed when the sun is below the astronomical horizon or
    # below the terrain horizon in its current azimuthal sector.
    solp <- computeSolarSeriesRCpp(latM[i, j], lonM[i, j], solYear, solMonth, solDay, solHour)
    solaralt <- (pi / 2) - solp$zenr
    ha <- horaTan[cbind(i, j, solp$sindex + 1L)]
    shadowVec <- as.numeric((solaralt <= 0) | (ha > tan(solaralt)))

    res <- .runpointmodelCore(
      wx, lat = latM[i, j], lon = lonM[i, j],
      meanhgt = hgt[i, j], meanpai = pai[i, j], meanx = xv[i, j],
      meanclump = clumpm[i, j], meangsmax = gsmaxM[i, j],
      meanleafr = lref[i, j], meanleaft = ltra[i, j],
      meangroundr = grefm[i, j], soiltypename = soiltypeNameM[i, j],
      modhabitat = modhabitatM[i, j], reqhgt = reqhgt, zin = zin, uzin = uzin,
      nlayers = nlayers, totalDepth = totalDepth,
      surface_organicmu = surface_organicmu, FreeDrain = FreeDrain,
      Lfrac = if (is.na(LfracM[i, j])) 1 else LfracM[i, j], matemp = matemp, maxIter = maxIter,
      tolerance = tolerance, pooling = pooling, slope = slope[i, j], aspect = aspect[i, j],
      svfa = svfa[i, j], shelterc = shelterArr[i, j, ], shadow = shadowVec,
      soilinit = soilinit, soilinitWarn = (k == 1L),
      phys = lapply(physM, function(m) m[i, j]))

    m <- res$model
    Tz <- m$Tzabove; RHz <- m$RHzabove; windz <- m$windzabove
    windzExchange <- windz
    # For non-positive reqhgt this pathway reports the solved soil-profile state at
    # that depth; above ground it reports the surface soil-water state. At exactly
    # zero height this convention differs from rungridmodel(), which treats zero as
    # the ground-surface regime.
    soilm <- if (!is.na(reqhgt) && reqhgt <= 0.0) m$thetazbelow else m$theta0

    if (!is.na(reqhgt) && reqhgt <= 0.0) {
      Tz <- m$Tzbelow; RHz <- rep(NA_real_, nT); windz <- rep(NA_real_, nT)
    } else if (!is.na(reqhgt) && reqhgt > 0.0 && hgt[i, j] > 0.0 && reqhgt < hgt[i, j]) {
      {
        bc <- belowCanopyProfilePointCpp(
          Tcanopy = m$Tcanopy, Tground = m$Tground, groundhr = m$groundhr,
          rBL = m$rBL, L = m$LL, uf = m$uf,
          d = dMat[i, j], zm = zmMat[i, j], hgt = hgt[i, j], pai = pai[i, j],
          Ta = res$weather$temp, rh = res$weather$relhum, pk = wx$pres,
          rSurf = m$rSurf, hSurf = m$hSurf,
          reqhgt = reqhgt, zref = res$zref)
        Tz <- bc$Tz; RHz <- bc$RHz
        # Within the canopy, attenuate the solved canopy-top wind exponentially with
        # depth using the cell's own friction velocity, canopy height and PAI.
        uhFloored <- pmax(m$uh, 1e-6)
        Be <- m$uf / uhFloored
        Lc <- 1 / (0.25 * pai[i, j] / hgt[i, j])
        Lm <- 2 * Be^3 * Lc
        windz <- uhFloored * exp(Be * (reqhgt - hgt[i, j]) / Lm)
        # The leaf boundary layer uses this exchange-scale wind; the reported
        # wind is scaled to the actual wind.
        windzExchange <- windz
        windz <- windz * m$windScale
      }
    }

    # Determine the amount of foliage above reqhgt once for use by both shortwave
    # and longwave exchange at that height.
    hgtCell <- hgt[i, j]; paiCell <- pai[i, j]
    paiaCell <- if (!is.na(reqhgt) && reqhgt >= 0.0 && paiCell > 0.0 && !is.na(group[i, j])) {
      .paiaboveheightPFT(reqhgt, hgtCell, paiCell, group[i, j])
    } else {
      0.0
    }

    # Reconstruct upward and downward longwave fluxes at reqhgt from this cell's
    # solved ground/canopy temperatures and its canopy geometry. Longwave is not
    # defined below the soil surface.
    Rlwdown <- rep(NA_real_, nT); Rlwup <- rep(NA_real_, nT)
    if (!is.na(reqhgt) && reqhgt >= 0.0) {
      lw <- longwaveGridCpp(
        Ts_est = array(m$Tground, dim = c(1L, 1L, nT)),
        Tcanopy_est = array(m$Tcanopy, dim = c(1L, 1L, nT)), Ta = wx$temp,
        hgt = matrix(hgtCell, 1, 1), pai = matrix(paiCell, 1, 1),
        paia = matrix(paiaCell, 1, 1), svfa = matrix(svfa[i, j], 1, 1),
        lwdown = wx$lwdown, reqhgt = reqhgt)
      Rlwdown <- as.vector(lw$Rlwdown); Rlwup <- as.vector(lw$Rlwup)
    }

    # Resolve direct, diffuse and upward shortwave radiation at reqhgt from this
    # cell's solar geometry, terrain exposure and canopy optical properties. The
    # same calculation supplies absorbed shortwave/PAR for the leaf energy balance.
    Rdirdown <- rep(NA_real_, nT); Rdifdown <- rep(NA_real_, nT); Rswup <- rep(NA_real_, nT)
    radLsw <- rep(NA_real_, nT); radLpar <- rep(NA_real_, nT)
    if (!is.na(reqhgt) && reqhgt >= 0.0) {
      rad <- radiationPointCpp(
        pai = paiCell, paia = paiaCell, x = xv[i, j], lref = lref[i, j], ltra = ltra[i, j],
        clump = clumpm[i, j], gref = grefm[i, j], svfa = svfa[i, j],
        slope = slope[i, j], aspect = aspect[i, j],
        zend = solp$zend, azid = solp$azid, shadow = shadowVec,
        Rsw = wx$swdown, Rdif = wx$difrad)
      Rdirdown <- rad$Rbdown; Rdifdown <- rad$Rddown; Rswup <- rad$Rdup
      radLsw <- rad$radLsw; radLpar <- rad$radLpar
    }

    # Leaf temperature is the resolved leaf energy-balance temperature at reqhgt,
    # forced by the local air temperature and humidity from this cell's own vertical profile.
    tleafOut <- rep(NA_real_, nT)
    if (!is.na(reqhgt) && reqhgt > 0.0 && hgtCell > 0.0 && reqhgt < hgtCell) {
      radLlw <- 0.5 * 0.97 * (Rlwdown + Rlwup) # 0.97 = mc::surfaceEmissivity
      Ca <- rep(Cafromyear(stats::median(solYear)), nT)
      radLparV <- radLpar; dim(radLparV) <- c(1L, 1L, nT)
      radLswV <- radLsw; dim(radLswV) <- c(1L, 1L, nT)
      radLlwV <- radLlw; dim(radLlwV) <- c(1L, 1L, nT)
      lp <- leafphys
      tleafArr <- leafTempCpp(
        radLpar = radLparV, radLsw = radLswV, radLlw = radLlwV, uz = windzExchange,
        Ta = Tz, rh = RHz, pk = wx$pres, Ca = Ca, psi_r = m$psi_r,
        hgt = matrix(hgtCell, 1, 1), pai = matrix(paiCell, 1, 1),
        Vcmax25 = matrix(lp$Vcmax25[i, j], 1, 1), Tup = matrix(lp$Tup[i, j], 1, 1),
        Tlw = matrix(lp$Tlw[i, j], 1, 1), Dcrit = matrix(lp$Dcrit[i, j], 1, 1),
        alpha = matrix(lp$alpha[i, j], 1, 1), f0 = matrix(lp$f0[i, j], 1, 1),
        fd = matrix(lp$fd[i, j], 1, 1), psi50 = matrix(lp$psi50[i, j], 1, 1),
        apsi = matrix(lp$apsi[i, j], 1, 1), rpmin = matrix(lp$rpmin[i, j], 1, 1),
        leafd = matrix(lp$leafd[i, j], 1, 1), isC3 = matrix(lp$isC3[i, j], 1, 1),
        reqhgt = reqhgt)
      tleafOut <- as.vector(tleafArr)
    }

    list(Tz = Tz, tleaf = tleafOut, soilm = soilm,
         relhum = RHz, windspeed = windz,
         Rdirdown = Rdirdown, Rdifdown = Rdifdown, Rswup = Rswup,
         Rlwdown = Rlwdown, Rlwup = Rlwup)[live]
  }

  coresN <- .resolveCores(cores)
  oldplan <- future::plan()
  on.exit(future::plan(oldplan), add = TRUE)
  if (coresN > 1) {
    future::plan(future::multisession, workers = coresN)
  } else {
    future::plan(future::sequential)
  }
  results <- future.apply::future_lapply(seq_len(nrow(validIJ)), runOneCell,
                                          future.seed = TRUE)

  # Reassemble the independently solved cells onto the original terrain grid and
  # preserve the input coordinate system in the returned raster time series.
  resultFull <- stats::setNames(
    lapply(live, function(f) array(NA_real_, dim = c(rows, cols, nT))), live)
  for (k in seq_len(nrow(validIJ))) {
    r <- results[[k]]
    if (!is.list(r)) next # a cell with an unfillable gap in its weather
    i <- validIJ[k, 1]; j <- validIJ[k, 2]
    for (fn in live) resultFull[[fn]][i, j, ] <- r[[fn]]
  }

  result <- stats::setNames(lapply(requested, function(f) resultFull[[f]]), requested)
  outList <- .gridOutputRasters(result, dtm, live)

  if (isTRUE(saveout)) {
    saveRDS(.wrapForWorker(outList), filename)
  }
  outList
}
