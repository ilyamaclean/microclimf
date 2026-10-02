# microclimf (development version)

`microclimf` is an R package for fast, mechanistic modelling of microclimate above, within and
below vegetation canopies and below ground. From standard hourly weather data, a digital
elevation model, and vegetation and soil properties, it predicts gridded temperature, humidity,
wind speed, radiation and soil moisture across real landscapes.

## About this branch

This is the `dev` branch: version 3 of the package, a substantial rebuild of the model.

- The `main` branch holds version 2, the version described in Maclean (2026) *Methods in Ecology
  and Evolution* 17: 1112-1123. It remains the default until snow is integrated here.
- **Snow is not yet represented in this version.** If you need the snow model, use `main`.
- Function names, inputs and outputs differ in places from version 2. The vignette
  "Running microclimf" lists what has changed.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ilyamaclean/microclimf", ref = "dev")
```

The package contains C++ code, so a compiler is needed: Rtools on Windows, Xcode command line
tools on macOS. Its dependency `terravars` is installed from GitHub automatically; the others
are on CRAN. The vignettes come already built, so nothing extra is needed to read them.

To install version 2 instead:

```r
remotes::install_github("ilyamaclean/microclimf")
```

## Quick start

```r
library(microclimf)
library(terra)
# Run the point model for a reference location, using the inbuilt datasets
micropoint <- runpointmodel(climdata, reqhgt = 0.05, dtmcaerth, vegp, soilc)
# Keep the hottest day of each month
micropoint <- subsetpointmodel(micropoint, tstep = "month", what = "tmax")
# Run the grid model for 5 cm above ground
mout <- rungridmodel(micropoint, reqhgt = 0.05, dtmcaerth, vegp, soilc)
# Air temperature at 13:00 on 20 June 2017
plot(mout$Tz[[134]])
```

## Documentation

Two vignettes are included:

```r
vignette("running-microclimf")  # how to prepare inputs and run the model
vignette("modelequations")      # the equations the model solves
```

Each function also has its own help page, for example `?rungridmodel`.

## Main functions

| Function | Purpose |
|---|---|
| `runpointmodel()` | Solves the full model for a reference location |
| `subsetpointmodel()` | Selects representative days from the point model |
| `rungridmodel()` | Derives gridded microclimate from the point model |
| `rungridmodelbig()`, `mosaicblend()` | Run large areas in tiles and recombine them |
| `runpointmodelasgrid()` | Runs the full model independently for every grid cell |
| `runbioclim()` | Microclimate equivalents of the 19 bioclimatic variables |

## Problems and questions

Please report problems at <https://github.com/ilyamaclean/microclimf/issues>, saying that you
are using the `dev` branch.
