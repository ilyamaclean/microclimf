// bioclimate.cpp
// Microclimate-derived BIO1-BIO19 summaries. Temperature metrics are
// calculated from representative hourly microclimate temperatures; the
// moisture metrics use modelled soil moisture rather than precipitation.
// See bioclimate.h for the scientific definitions and sampling convention.
#include "bioclimate.h"
#include <cmath>
#include <numeric>

using namespace Rcpp;

namespace bioclimate {

    // Sample standard deviation used by the seasonality metrics.
    double calc_std_dev(NumericVector vec) {
        int n = vec.size();
        if (n <= 1) {
            return NA_REAL;
        }
        double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / n;
        double sum_squared_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            sum_squared_diff += std::pow(vec[i] - mean, 2.0);
        }
        double variance = sum_squared_diff / (n - 1);
        return std::sqrt(variance);
    }

    // Every BIO1-19 function below operates on a single cell's own
    // 336-hour representative-day time series, laid out as: the 12
    // monthly-median days first (hours 0-287, one representative day per
    // calendar month), followed by representative hottest and coldest days
    // (hours 288-311 and 312-335). With multiple years of weather, each is
    // the central-ranked day among the annual hottest/coldest days rather
    // than the record-wide extreme; see R/bioclimate.R for the selection.

    double bioclim1(NumericVector Tz) {
        // Annual Mean Temperature: the mean over the 12 monthly-median
        // representative days, a broad summary of a location's overall
        // thermal regime.
        double out = 0.0;
        for (int i = 0; i < 288; ++i) out += Tz[i];
        return out / 288.0;
    }

    double bioclim2(NumericVector Tz) {
        // Mean Diurnal Temperature Range: the average, across the 12
        // monthly-median days, of each day's own difference between its
        // warmest and coolest hour -- how much temperature typically
        // swings within a single day, as distinct from how much it swings
        // across the year (BIO7).
        NumericVector dtr(12);
        int index = 0;
        for (int day = 0; day < 12; day++) {
            double tmx = -273.15;
            double tmn = 273.15;
            for (int hr = 0; hr < 24; hr++) {
                if (Tz[index] > tmx) tmx = Tz[index];
                if (Tz[index] < tmn) tmn = Tz[index];
                index++;
            }
            dtr[day] = tmx - tmn;
        }
        double out = 0.0;
        for (int day = 0; day < 12; day++) out += dtr[day];
        return out / 12;
    }

    double bioclim4(NumericVector Tz) {
        // Temperature Seasonality: how much the mean temperature varies
        // from month to month across the year, expressed as the standard
        // deviation of the 12 monthly-median days' own mean temperature
        // (x100, the standard WorldClim scaling convention).
        NumericVector monmean(12);
        int index = 0;
        for (int mth = 0; mth < 12; mth++) {
            monmean[mth] = 0.0;
            for (int hr = 0; hr < 24; hr++) {
                monmean[mth] += Tz[index];
                index++;
            }
            monmean[mth] = monmean[mth] / 24;
        }
        return calc_std_dev(monmean) * 100.0;
    }

    double bioclim5(NumericVector Tz) {
        // Maximum Temperature of the Warmest Month: the single hottest
        // hour of the year's own representative day -- an extreme, not an
        // average, capturing peak thermal exposure.
        double tmx = -273.15;
        for (int i = 288; i < 312; i++) if (Tz[i] > tmx) tmx = Tz[i];
        return tmx;
    }

    double bioclim6(NumericVector Tz) {
        // Minimum Temperature of the Coldest Month: the single coldest
        // hour of the year's own representative day -- the opposite
        // extreme to BIO5, capturing peak cold exposure.
        double tmn = 273.15;
        for (int i = 312; i < 336; i++) if (Tz[i] < tmn) tmn = Tz[i];
        return tmn;
    }

    // BIO8-11/16-19 (quarter means, below) each average their variable
    // over whichever representative hours fall within the wettest/driest/
    // warmest/coldest three-calendar-month window. Only the centre month of
    // that window is identified from the full weather record: calendar-month
    // means of forcing precipitation define wet/dry quarters and calendar-
    // month means of forcing air temperature define warm/cold quarters
    // (R/bioclimate.R). The indices are then the representative hours whose
    // months fall in that window: normally 72, but more if the separately
    // selected hottest or coldest representative day falls inside it. The
    // true number of matching hours is always used as the divisor.

    double bioclim8(NumericVector Tz, IntegerVector wetq) {
        // Mean temperature during the three-month period identified from
        // the full weather record as the wettest quarter (highest rainfall).
        // Soil moisture replaces precipitation only for BIO12-BIO19; it is
        // not what defines the wet/dry quarter used here.
        double out = 0.0;
        for (R_xlen_t i = 0; i < wetq.size(); i++) out += Tz[wetq[i]];
        return out / (double)wetq.size();
    }

    double bioclim9(NumericVector Tz, IntegerVector dryq) {
        // Mean Temperature of the Driest Quarter.
        double out = 0.0;
        for (R_xlen_t i = 0; i < dryq.size(); i++) out += Tz[dryq[i]];
        return out / (double)dryq.size();
    }

    double bioclim10(NumericVector Tz, IntegerVector hotq) {
        // Mean Temperature of the Warmest Quarter.
        double out = 0.0;
        for (R_xlen_t i = 0; i < hotq.size(); i++) out += Tz[hotq[i]];
        return out / (double)hotq.size();
    }

    double bioclim11(NumericVector Tz, IntegerVector colq) {
        // Mean Temperature of the Coldest Quarter.
        double out = 0.0;
        for (R_xlen_t i = 0; i < colq.size(); i++) out += Tz[colq[i]];
        return out / (double)colq.size();
    }

    double bioclim12(NumericVector soilm) {
        // Annual mean soil moisture -- the moisture-availability
        // counterpart to BIO1, standing in for WorldClim's Annual
        // Precipitation using the actual moisture a plant's roots or an
        // organism at the surface would experience, rather than rainfall
        // input alone.
        double out = 0.0;
        for (int i = 0; i < 288; i++) out += soilm[i];
        return out / 288.0;
    }

    double bioclim13(NumericVector soilm) {
        // Soil moisture of the Wettest Month: the single highest soil
        // moisture value across the whole representative-day record.
        double out = 0.0;
        for (R_xlen_t i = 0; i < soilm.size(); i++) if (soilm[i] > out) out = soilm[i];
        return out;
    }

    double bioclim14(NumericVector soilm) {
        // Soil moisture of the Driest Month: the single lowest soil
        // moisture value across the whole representative-day record.
        double out = 1.0;
        for (R_xlen_t i = 0; i < soilm.size(); i++) if (soilm[i] < out) out = soilm[i];
        return out;
    }

    double bioclim15(NumericVector soilm) {
        // BIO15 analogue: mean soil moisture over the 12 monthly
        // representative days divided by the standard deviation over all
        // 336 representative hours. This is the inverse of WorldClim's
        // coefficient-of-variation form, so larger values mean less relative
        // variability.
        double me = 0.0;
        for (int i = 0; i < 288; i++) me += soilm[i];
        me = me / 288.0;
        double sd = calc_std_dev(soilm);
        return me / sd;
    }

    double bioclim16(NumericVector soilm, IntegerVector wetq) {
        // Mean soil moisture of the Wettest Quarter.
        double me = 0.0;
        for (R_xlen_t i = 0; i < wetq.size(); i++) me += soilm[wetq[i]];
        return me / (double)wetq.size();
    }

    double bioclim17(NumericVector soilm, IntegerVector dryq) {
        // Mean soil moisture of the Driest Quarter.
        double me = 0.0;
        for (R_xlen_t i = 0; i < dryq.size(); i++) me += soilm[dryq[i]];
        return me / (double)dryq.size();
    }

    double bioclim18(NumericVector soilm, IntegerVector hotq) {
        // Mean soil moisture of the Warmest Quarter.
        double me = 0.0;
        for (R_xlen_t i = 0; i < hotq.size(); i++) me += soilm[hotq[i]];
        return me / (double)hotq.size();
    }

    double bioclim19(NumericVector soilm, IntegerVector colq) {
        // Mean soil moisture of the Coldest Quarter.
        double me = 0.0;
        for (R_xlen_t i = 0; i < colq.size(); i++) me += soilm[colq[i]];
        return me / (double)colq.size();
    }

    NumericMatrix bioclimfill(int rows, int cols) {
        NumericMatrix bio(rows, cols);
        std::fill(bio.begin(), bio.end(), NA_REAL);
        return bio;
    }

    List runbioclimCpp(NumericVector Tz, NumericVector soilm, std::vector<bool> out,
        IntegerVector wetq, IntegerVector dryq, IntegerVector hotq, IntegerVector colq)
    {
        IntegerVector dims = Tz.attr("dim");
        int rows = dims[0];
        int cols = dims[1];
        int tsteps = dims[2];

        // BIO3 and BIO7 are derived from the primary temperature metrics:
        // annual range is BIO5 - BIO6, and isothermality is BIO2 / BIO7.
        // The latter is returned as a ratio (not multiplied by 100). These
        // dependencies are calculated even when they are needed only as
        // intermediates for another requested BIO variable.
        bool needBio7 = out[6] || out[2];
        bool needBio6 = out[5] || needBio7;
        bool needBio5 = out[4] || needBio7;
        bool needBio2 = out[1] || out[2];

        NumericMatrix bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10,
            bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19;

        if (out[0])   bio1  = bioclimfill(rows, cols);
        if (needBio2) bio2  = bioclimfill(rows, cols);
        if (out[2])   bio3  = bioclimfill(rows, cols);
        if (out[3])   bio4  = bioclimfill(rows, cols);
        if (needBio5) bio5  = bioclimfill(rows, cols);
        if (needBio6) bio6  = bioclimfill(rows, cols);
        if (needBio7) bio7  = bioclimfill(rows, cols);
        if (out[7])  bio8  = bioclimfill(rows, cols);
        if (out[8])  bio9  = bioclimfill(rows, cols);
        if (out[9])  bio10 = bioclimfill(rows, cols);
        if (out[10]) bio11 = bioclimfill(rows, cols);
        if (out[11]) bio12 = bioclimfill(rows, cols);
        if (out[12]) bio13 = bioclimfill(rows, cols);
        if (out[13]) bio14 = bioclimfill(rows, cols);
        if (out[14]) bio15 = bioclimfill(rows, cols);
        if (out[15]) bio16 = bioclimfill(rows, cols);
        if (out[16]) bio17 = bioclimfill(rows, cols);
        if (out[17]) bio18 = bioclimfill(rows, cols);
        if (out[18]) bio19 = bioclimfill(rows, cols);

        // Processed one grid cell at a time: each cell's own 336-hour
        // temperature and soil-moisture series is extracted, reduced to
        // BIO1-19, and discarded before moving to the next cell, so only
        // one cell's time series is ever held in memory at once.
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double val = Tz[i + rows * j];
                if (!NumericVector::is_na(val)) {
                    NumericVector Tzv(tsteps);
                    NumericVector soilmv(tsteps);
                    for (int k = 0; k < tsteps; ++k) {
                        int idx = i + rows * j + cols * rows * k;
                        Tzv[k] = Tz[idx];
                        soilmv[k] = soilm[idx];
                    }
                    if (out[0])   bio1(i, j)  = bioclim1(Tzv);
                    if (needBio2) bio2(i, j)  = bioclim2(Tzv);
                    if (out[3])   bio4(i, j)  = bioclim4(Tzv);
                    if (needBio5) bio5(i, j)  = bioclim5(Tzv);
                    if (needBio6) bio6(i, j)  = bioclim6(Tzv);
                    if (out[7])  bio8(i, j)  = bioclim8(Tzv, wetq);
                    if (out[8])  bio9(i, j)  = bioclim9(Tzv, dryq);
                    if (out[9])  bio10(i, j) = bioclim10(Tzv, hotq);
                    if (out[10]) bio11(i, j) = bioclim11(Tzv, colq);
                    if (out[11]) bio12(i, j) = bioclim12(soilmv);
                    if (out[12]) bio13(i, j) = bioclim13(soilmv);
                    if (out[13]) bio14(i, j) = bioclim14(soilmv);
                    if (out[14]) bio15(i, j) = bioclim15(soilmv);
                    if (out[15]) bio16(i, j) = bioclim16(soilmv, wetq);
                    if (out[16]) bio17(i, j) = bioclim17(soilmv, dryq);
                    if (out[17]) bio18(i, j) = bioclim18(soilmv, hotq);
                    if (out[18]) bio19(i, j) = bioclim19(soilmv, colq);
                    if (needBio7) bio7(i, j)  = bio5(i, j) - bio6(i, j);
                    if (out[2])   bio3(i, j)  = bio2(i, j) / bio7(i, j);
                }
            }
        }

        List outp;
        if (out[0])  outp["bio1"]  = bio1;
        if (out[1])  outp["bio2"]  = bio2;
        if (out[2])  outp["bio3"]  = bio3;
        if (out[3])  outp["bio4"]  = bio4;
        if (out[4])  outp["bio5"]  = bio5;
        if (out[5])  outp["bio6"]  = bio6;
        if (out[6])  outp["bio7"]  = bio7;
        if (out[7])  outp["bio8"]  = bio8;
        if (out[8])  outp["bio9"]  = bio9;
        if (out[9])  outp["bio10"] = bio10;
        if (out[10]) outp["bio11"] = bio11;
        if (out[11]) outp["bio12"] = bio12;
        if (out[12]) outp["bio13"] = bio13;
        if (out[13]) outp["bio14"] = bio14;
        if (out[14]) outp["bio15"] = bio15;
        if (out[15]) outp["bio16"] = bio16;
        if (out[16]) outp["bio17"] = bio17;
        if (out[17]) outp["bio18"] = bio18;
        if (out[18]) outp["bio19"] = bio19;
        return outp;
    }

} // namespace bioclimate

// [[Rcpp::export]]
Rcpp::List runbioclimCpp(Rcpp::NumericVector Tz, Rcpp::NumericVector soilm,
    std::vector<bool> out, Rcpp::IntegerVector wetq, Rcpp::IntegerVector dryq,
    Rcpp::IntegerVector hotq, Rcpp::IntegerVector colq)
{
    return bioclimate::runbioclimCpp(Tz, soilm, out, wetq, dryq, hotq, colq);
}
