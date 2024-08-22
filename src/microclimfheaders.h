// microclimfheaders.h
// Radiation model
#include <Rcpp.h>
struct radmodel {
    std::vector<double> ground; // ground absorbed radiation (W/m^2)
    std::vector<double> canopy; // canopy absorbed radiation (W/m^2)
    std::vector<double> albedo; // albedo
};
struct radmodel2 {
    std::vector<double> radGsw;
    std::vector<double> radGlw;
    std::vector<double> radCsw;
    std::vector<double> radClw;
    std::vector<double> Rbdown;
    std::vector<double> Rddown;
    std::vector<double> Rdup;
    std::vector<double> radLsw;
    std::vector<double> lwout;
};
struct Gmodel {
    std::vector<double> G;
    std::vector<double> Gmin;
    std::vector<double> Gmax;
};
struct windmodel {
    std::vector<double> uf;
    std::vector<double> uz;
    std::vector<double> gHa;
};
struct soilmodel {
    std::vector<double> Tg;
    std::vector<double> G;
    std::vector<double> kDDg;
    std::vector<double> kDDp;
};
struct abovemodel {
    std::vector<double> Tz;
    std::vector<double> tleaf;
    std::vector<double> rh;
    std::vector<double> lwdn;
    std::vector<double> lwup;
};
struct microm {
    Rcpp::NumericVector Tz;
    Rcpp::NumericVector tleaf;
    Rcpp::NumericVector relhum;
    Rcpp::NumericVector soilm;
    Rcpp::NumericVector uz;
    Rcpp::NumericVector Rdirdown;
    Rcpp::NumericVector Rdifdown;
    Rcpp::NumericVector Rlwdown;
    Rcpp::NumericVector Rswup;
    Rcpp::NumericVector Rlwup;
};
