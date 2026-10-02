// bioclimate.h
// BIO1-BIO19: the standard "WorldClim-style" set of annual, seasonal and
// extreme temperature/moisture summary variables widely used in species
// distribution and climate-niche modelling. Here they are computed from
// this model's own microclimate output (air/leaf temperature and soil
// moisture) rather than from coarse regional climate data, giving
// fine-scale bioclimatic variables at the scale a plant or animal actually
// experiences. Following this model's own convention, BIO12-BIO19 (the
// variables that describe moisture rather than temperature) are computed
// from soil moisture rather than precipitation -- the more physiologically
// relevant quantity for what most organisms actually respond to at
// microclimate scale, and one this model already estimates directly.
//
// To keep the cost of this reasonable, these variables are computed from a
// small set of 14 representative days (336 hours) -- the temperature-
// median day of each calendar month, plus representative hottest and coldest
// days -- rather than from a full year run at every grid cell. With multiple
// years of weather, the latter are the central-ranked days among the annual
// hottest/coldest days, not the record-wide extremes; see R/bioclimate.R.
//
// The C++ layer receives temperature and soil-moisture arrays already
// reduced to these representative hours and converts each grid cell to the
// requested BIO summaries; spatial raster handling remains on the R side.
#include <Rcpp.h>

namespace bioclimate {

    // Sample standard deviation (n-1 denominator), NA_REAL if n <= 1.
    double calc_std_dev(Rcpp::NumericVector vec);

    // BIO1-BIO19 formulas, each operating on one grid cell's own already-
    // subsetted 336-hour (14-day) temperature or soil-moisture time
    // series:
    //   BIO1  Annual Mean Temperature
    //   BIO2  Mean Diurnal Temperature Range (mean of each day's own max
    //         minus min temperature)
    //   BIO4  Temperature Seasonality (standard deviation of monthly mean
    //         temperature, x100 by WorldClim convention)
    //   BIO5  Maximum Temperature of the Warmest Month
    //   BIO6  Minimum Temperature of the Coldest Month
    //   BIO8/9/10/11  Mean Temperature of the Wettest/Driest/Warmest/
    //         Coldest Quarter. Wet/dry quarters are defined by forcing
    //         precipitation and warm/cold quarters by forcing air temperature,
    //         never by the modelled temperature or soil moisture being averaged;
    //         see R/bioclimate.R's quarter-selection.
    //   BIO12/13/14/15  Annual mean soil moisture, wettest/driest soil
    //         moisture, and the package's soil-moisture seasonality index
    //         (mean / SD; larger values therefore mean less variability),
    //         using soil moisture rather than precipitation
    //   BIO16/17/18/19  Mean soil moisture of the Wettest/Driest/Warmest/
    //         Coldest Quarter, using those same forcing-defined quarter windows.
    // BIO3 (Isothermality, how large the diurnal temperature swing is
    // relative to the total annual swing) and BIO7 (Temperature Annual
    // Range) are derived quantities, not independent measurements --
    // BIO7 = BIO5 - BIO6 and BIO3 = BIO2 / BIO7 -- so they are computed
    // inline in runbioclimCpp() once BIO2/BIO5/BIO6 are available for
    // that cell, rather than as their own separate functions.
    double bioclim1(Rcpp::NumericVector Tz);
    double bioclim2(Rcpp::NumericVector Tz);
    double bioclim4(Rcpp::NumericVector Tz);
    double bioclim5(Rcpp::NumericVector Tz);
    double bioclim6(Rcpp::NumericVector Tz);
    double bioclim8(Rcpp::NumericVector Tz, Rcpp::IntegerVector wetq);
    double bioclim9(Rcpp::NumericVector Tz, Rcpp::IntegerVector dryq);
    double bioclim10(Rcpp::NumericVector Tz, Rcpp::IntegerVector hotq);
    double bioclim11(Rcpp::NumericVector Tz, Rcpp::IntegerVector colq);
    double bioclim12(Rcpp::NumericVector soilm);
    double bioclim13(Rcpp::NumericVector soilm);
    double bioclim14(Rcpp::NumericVector soilm);
    double bioclim15(Rcpp::NumericVector soilm);
    double bioclim16(Rcpp::NumericVector soilm, Rcpp::IntegerVector wetq);
    double bioclim17(Rcpp::NumericVector soilm, Rcpp::IntegerVector dryq);
    double bioclim18(Rcpp::NumericVector soilm, Rcpp::IntegerVector hotq);
    double bioclim19(Rcpp::NumericVector soilm, Rcpp::IntegerVector colq);

    // rows x cols matrix pre-filled with NA -- used to allocate an output
    // field's memory only for the BIO variables actually requested.
    Rcpp::NumericMatrix bioclimfill(int rows, int cols);

    // Convert each grid cell's representative temperature and soil-moisture
    // record into the requested annual, seasonal and extreme BIO summaries.
    // Quarter index vectors identify the representative hours belonging to
    // the wettest/driest three-month periods defined by forcing precipitation
    // and the warmest/coldest periods defined by forcing air temperature.
    Rcpp::List runbioclimCpp(Rcpp::NumericVector Tz, Rcpp::NumericVector soilm,
        std::vector<bool> out, Rcpp::IntegerVector wetq, Rcpp::IntegerVector dryq,
        Rcpp::IntegerVector hotq, Rcpp::IntegerVector colq);

} // namespace bioclimate
