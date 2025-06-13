#include <Rcpp.h>
#include "microclimfheaders.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <tuple>
#include <cfloat> 
using namespace Rcpp;
// ** Calculates Astronomical Julian day ** //
int juldayCpp(int year, int month, int day)
{
    double dd = day + 0.5;
    int madj = month + (month < 3) * 12;
    int yadj = year + (month < 3) * -1;
    double j = std::trunc(365.25 * (yadj + 4716)) + std::trunc(30.6001 * (madj + 1)) + dd - 1524.5;
    int b = 2 - std::trunc(yadj / 100) + std::trunc(std::trunc(yadj / 100) / 4);
    int jd = static_cast<int>(j + (j > 2299160) * b);
    return jd;
}
// ** Calculates solar time ** //
double soltimeCpp(int jd, double lt, double lon)
{

    double m = 6.24004077 + 0.01720197 * (jd - 2451545);
    double eot = -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932);
    double st = lt + (4 * lon + eot) / 60;
    return st;
}
// ** Calculates solar position ** //
// [[Rcpp::export]]
std::vector<double> solpositionCpp(double lat, double lon, int year, int month, int day, double lt)
{
    int jd = juldayCpp(year, month, day);
    double st = soltimeCpp(jd, lt, lon);
    // Calculate solar zenith (degrees)
    double latr = lat * M_PI / 180;
    double tt = 0.261799 * (st - 12);
    double dec = (M_PI * 23.5 / 180) * cos(2 * M_PI * ((jd - 159.5) / 365.25));
    double coh = sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt);
    double z = acos(coh) * (180 / M_PI);
    // Calculate solar azimuth (degrees)
    double sh = sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt);
    double hh = atan(sh / sqrt(1 - sh * sh));
    double sazi = cos(dec) * sin(tt) / cos(hh);
    double cazi = (sin(latr) * cos(dec) * cos(tt) - cos(latr) * sin(dec)) /
        sqrt(pow(cos(dec) * sin(tt), 2) + pow(sin(latr) * cos(dec) * cos(tt) - cos(latr) * sin(dec), 2));
    double sqt = 1 - sazi * sazi;
    if (sqt < 0) sqt = 0;
    double azi = 180 + (180 * atan(sazi / sqrt(sqt))) / M_PI;
    if (cazi < 0) {
        if (sazi < 0) {
            azi = 180 - azi;
        }
        else {
            azi = 540 - azi;
        }
    }
    // Define and return output variable
    std::vector<double> solpos(2, 0.0);
    solpos[0] = z;
    solpos[1] = azi;
    return solpos;
}
// ** Calculate solar index ** //
// [[Rcpp::export]]
double solarindexCpp(double slope, double aspect, double zen, double azi, bool shadowmask = false)
{
    double si;
    if (zen > 90.0 && !shadowmask) {
        si = 0;
    }
    else {
        if (slope == 0.0) {
            si = cos(zen * M_PI / 180);
        }
        else {
            si = cos(zen * M_PI / 180) * cos(slope * M_PI / 180) + sin(zen * M_PI / 180) *
                sin(slope * M_PI / 180) * cos((azi - aspect) * M_PI / 180);
        }
    }
    if (si < 0.0) si = 0.0;
    return si;
}
// ** Calculate clear sky radiation ** //
// [[Rcpp::export]]
std::vector<double> clearskyradCpp(std::vector<int> year, std::vector<int> month, std::vector<int> day,
    std::vector<double> lt, double lat, double lon, std::vector<double> tc, std::vector<double> rh, 
    std::vector<double> pk)
{
    std::vector<double> Ic(year.size());
    for (size_t i = 0; i < Ic.size(); ++i) {
        std::vector<double> sp = solpositionCpp(lat, lon, year[i], month[i], day[i], lt[i]);
        double zen = sp[0];
        double z = zen * M_PI / 180.0;
        if (zen <= 90.0) {
            double m = 35 * cos(z) * pow(1224 * cos(z) * cos(z) + 1, -0.5);
            double TrTpg = 1.021 - 0.084 * sqrt(m * 0.00949 * pk[i] + 0.051);
            double xx = log(rh[i] / 100) + ((17.27 * tc[i]) / (237.3 + tc[i]));
            double Td = (237.3 * xx) / (17.27 - xx);
            double u = exp(0.1133 - log(3.78) + 0.0393 * Td);
            double Tw = 1 - 0.077 * pow(u * m, 0.3);
            double Ta = 0.935 * m;
            double od = TrTpg * Tw * Ta;
            Ic[i] = 1352.778 * cos(z) * od;
        }
    }
    return Ic;
}
// ** Calculate canopy extinction coefficient for sloped ground surfaces ** //
std::vector<double> cankCpp(double zen, double x, double si) {
    double k;
    if (zen > 90.0) zen = 90.0;
    if (si < 0.0) si = 0.0;
    double Z = zen * M_PI / 180.0;
    // Calculate normal canopy extinction coefficient
    if (x == 1.0) {
        k = 1 / (2 * cos(Z));
    }
    else if (std::isinf(x)) {
        k = 1.0;
    }
    else if (x == 0.0) {
        k = tan(Z);
    }
    else {
        k = sqrt(x * x + (tan(Z) * tan(Z))) / (x + 1.774 * pow((x + 1.182), -0.733));
    }
    if (k > 6000.0) k = 6000.0;
    // Calculate adjusted k
    double kd = k * cos(Z) / si;
    if (si == 0) kd = 1.0;
    double Kc = 1.0 / si;
    if (si == 0.0) Kc = 600.0;
    std::vector<double> kparams(3, 0.0);
    kparams[0] = k;
    kparams[1] = kd;
    kparams[2] = Kc;
    return kparams;
}
// ** Calculates parameters for diffuse radiation using two-stream model ** //
std::vector<double> twostreamdifCpp(double pait, double om, double a, double gma, double h, double gref)
{
    // Adjust leaf area for clumping factor
    // Calculate base parameters
    double S1 = exp(-h * pait);
    double u1 = a + gma * (1 - 1 / gref);
    double u2 = a + gma * (1 - gref);
    double D1 = (a + gma + h) * (u1 - h) * 1 / S1 - (a + gma - h) * (u1 + h) * S1;
    double D2 = (u2 + h) * 1 / S1 - (u2 - h) * S1;
    // Calculate parameters
    double p1 = (gma / (D1 * S1)) * (u1 - h);
    double p2 = (-gma * S1 / D1) * (u1 + h);
    double p3 = (1 / (D2 * S1)) * (u2 + h);
    double p4 = (-S1 / D2) * (u2 - h);
    // Define and return output variable
    std::vector<double> params(8, 0.0);
    params[0] = p1;
    params[1] = p2;
    params[2] = p3;
    params[3] = p4;
    params[4] = u1;
    params[5] = S1;
    params[6] = D1;
    params[7] = D2;
    return params;
}
// ** Calculates parameters for direct radiation using two-stream model ** //
std::vector<double> twostreamdirCpp(double pait, double om, double a, double gma, double J, double del, double h, double gref,
    double k, double sig, double u1, double S1, double D1, double D2)
{
    // Calculate base parameters
    double ss = 0.5 * (om + J * del / k) * k;
    double sstr = om * k - ss;
    double S2 = exp(-k * pait);
    double u2 = a + gma * (1 - gref);
    double p5 = -ss * (a + gma - k) - gma * sstr;
    double v1 = ss - (p5 * (a + gma + k)) / sig;
    double v2 = ss - gma - (p5 / sig) * (u1 + k);
    double p6 = (1 / D1) * ((v1 / S1) * (u1 - h) - (a + gma - h) * S2 * v2);
    double p7 = (-1 / D1) * ((v1 * S1) * (u1 + h) - (a + gma + h) * S2 * v2);
    sig = -sig;
    double p8 = sstr * (a + gma + k) - gma * ss;
    double v3 = (sstr + gma * gref - (p8 / sig) * (u2 - k)) * S2;
    double p9 = (-1 / D2) * ((p8 / (sig * S1)) * (u2 + h) + v3);
    double p10 = (1 / D2) * (((p8 * S1) / sig) * (u2 - h) + v3);
    // Define and return output variable
    std::vector<double> params(6, 0.0);
    params[0] = p5;
    params[1] = p6;
    params[2] = p7;
    params[3] = p8;
    params[4] = p9;
    params[5] = p10;
    return params;
}
// ** Calculate absorbed shortwave radiation ** //
radmodel RadswabsCpp(double pai, double x, double lref, double ltra, double clump, double gref,
    double slope, double aspect, double lat, double lon,
    std::vector<int> year, std::vector<int> month, std::vector<int> day, std::vector<double> lt,
    std::vector<double> Rsw, std::vector<double> Rdif)
{
    // ** Define variables
    std::vector<double> radGsw(Rsw.size());
    std::vector<double> radCsw(Rsw.size());
    std::vector<double> albedo(Rsw.size());
    if (pai > 0.0) {
        // Calculate time-invariant variables
        double pait = pai;
        if (clump > 0.0) pait = pai / (1 - clump);
        double om = lref + ltra;
        double a = 1 - om;
        double del = lref - ltra;
        double J = 1.0 / 3.0;
        if (x != 1.0) {
            double mla = 9.65 * pow((3 + x), -1.65);
            if (mla > M_PI / 2) mla = M_PI / 2;
            J = cos(mla) * cos(mla);
        }
        double gma = 0.5 * (om + J * del);
        double h = sqrt(a * a + 2 * a * gma);
        // Calculate two-stream parameters (diffuse)
        std::vector<double> tspdif = twostreamdifCpp(pait, om, a, gma, h, gref);
        double p1 = tspdif[0];
        double p2 = tspdif[1];
        double p3 = tspdif[2];
        double p4 = tspdif[3];
        double u1 = tspdif[4];
        double S1 = tspdif[5];
        double D1 = tspdif[6];
        double D2 = tspdif[7];
        double trd = clump * clump;
        // Calculate albedo and ground flux
        double amx = gref;
        if (amx < lref) amx = lref;
        double albd = gref * trd + (1.0 - trd) * (p1 + p2);
        if (albd > amx) albd = amx;
        if (albd < 0.01) albd = 0.01;
        double groundRdd = trd + (1.0 - trd) * (p3 * exp(-h * pait) + p4 * exp(h * pait));
        // Calculate time-variant two-stream parameters (direct)
        for (size_t i = 0; i < Rsw.size(); ++i) {
            if (Rsw[i] > 0.0) {
                // Calculate solar variables
                std::vector<double> solp = solpositionCpp(lat, lon, year[i], month[i], day[i], lt[i]);
                double zen = solp[0];
                double azi = solp[1];
                double si = solarindexCpp(slope, aspect, zen, azi);
                if (zen > 90.0) zen = 90.0;
                // Calculate canopy extinction coefficient
                double cosz = cos(zen * M_PI / 180);
                std::vector<double> kp = cankCpp(zen, x, si);
                double kd = kp[1];
                double Kc = kp[2];
                // Calculate two-stream parameters (direct)      
                double sig = (kd * kd + gma * gma - pow((a + gma), 2.0));
                std::vector<double> tspdir = twostreamdirCpp(pait, om, a, gma, J, del, h, gref, kd, sig, u1, S1, D1, D2);
                double p5 = tspdir[0];
                double p6 = tspdir[1];
                double p7 = tspdir[2];
                double p8 = tspdir[3];
                double p9 = tspdir[4];
                double p10 = tspdir[5];
                // Calculate beam normalisations
                double Rbeam = (Rsw[i] - Rdif[i]) / cosz;
                if (Rbeam > 1352.0) Rbeam = 1352.0;
                double trb = pow(clump, Kc);
                if (trb > 0.999) trb = 0.999;
                if (trb < 0.0) trb = 0.0;
                double Rb = Rbeam * cosz;
                double trg = trb + exp(-kd * pait); // tranmission to ground though gaps and leaves
                double Rbc = (trg * si + (1 - trg) * cosz) * Rbeam;
                // Calculate albedo and ground flux
                double albb = trd * gref + (1.0 - trd) * (p5 / sig + p6 + p7);
                if (albb > amx) albb = amx;
                if (albb < 0.01) albb = 0.01;
                double groundRbdd = trb + (1.0 - trb) * ((p8 / -sig) * exp(-kd * pait) +
                    p9 * exp(-h * pait) + p10 * exp(h * pait));
                if (groundRbdd > amx) groundRbdd = amx;
                if (groundRbdd < 0.0) groundRbdd = 0.0;
                // Calculate canopy absorbed
                radCsw[i] = (1.0 - albd) * Rdif[i] + (1.0 - albb) * Rbc;
                // Calculate ground absorbed
                double Rgdif = groundRdd * Rdif[i] + groundRbdd * Rb;
                radGsw[i] = (1.0 - gref) * (Rgdif + exp(-kd * pait) * Rbeam * si);
                // Calculate albedo
                albedo[i] = 1.0 - (radCsw[i] / (Rdif[i] + Rb));
                if (albedo[i] > amx) albedo[i] = amx;
                if (albedo[i] < 0.01) albedo[i] = 0.01;
            }
            else {
                radGsw[i] = 0;
                radCsw[i] = 0;
                albedo[i] = lref;
            }
        }
    }
    else {
        for (size_t i = 0; i < Rsw.size(); ++i) {
            albedo[i] = gref;
            if (Rsw[i] > 0) {
                std::vector<double> solp = solpositionCpp(lat, lon, year[i], month[i], day[i], lt[i]);
                double zen = solp[0];
                double azi = solp[1];
                double si = solarindexCpp(slope, aspect, zen, azi);
                if (zen > 90.0) zen = 90.0;
                double dirr = (Rsw[i] - Rdif[i]) / cos(zen * M_PI / 180.0);
                radGsw[i] = (1 - gref) * (Rdif[i] + si * dirr);
                radCsw[i] = radGsw[i];
            }
            else {
                radGsw[i] = 0;
                radCsw[i] = 0;
            }
        }
    }
    radmodel out;
    out.ground = radGsw;
    out.canopy = radCsw;
    out.albedo = albedo;
    return out;
}
// ** Calculate molar density of air ** //
double phairCpp(double tc, double pk)
{
    double tk = tc + 273.15;
    double ph = 44.6 * (pk / 101.3) * (273.15 / tk);
    return ph;
}
// ** Calculate specific heat of air at constant pressure ** //
double cpairCpp(double tc)
{
    double cp = 2e-05 * pow(tc, 2) + 0.0002 * tc + 29.119;
    return cp;
}
// ** Calculate zero-plane displacement ** //
// [[Rcpp::export]]
double zeroplanedisCpp(double h, double pai)
{
    if (pai < 0.001) pai = 0.001;
    double d = (1 - (1 - exp(-sqrt(7.5 * pai))) / sqrt(7.5 * pai)) * h;
    return d;
}
// ** Calculate roughness length ** //
// [[Rcpp::export]]
double roughlengthCpp(double h, double pai, double d, double psi_h)
{
    double Be = sqrt(0.003 + (0.2 * pai) / 2);
    double zm = (h - d) * exp(-0.4 / Be) * exp(psi_h);
    if (zm < 0.0005) zm = 0.0005;
    return zm;
}
// **  Calculate integrated diabatic correction coefficient for momentum ** //
// [[Rcpp::export]]
double dpsimCpp(double ze)
{
    double psim;
    // unstable
    if (ze < 0) {
        double x = pow((1 - 15 * ze), 0.25);
        psim = log(pow((1 + x) / 2, 2) * (1 + pow(x, 2)) / 2) - 2 * atan(x) + M_PI / 2;
    }
    // stable
    else {
        psim = -4.7 * ze;
    }
    if (psim < -4) psim = -4;
    if (psim > 3) psim = 3;
    return psim;
}
// **  Calculate integrated diabatic correction coefficient for heat ** //
// [[Rcpp::export]]
double dpsihCpp(double ze)
{
    double psih;
    // unstable
    if (ze < 0) {
        double y = sqrt(1 - 9 * ze);
        psih = log(pow((1 + y) / 2, 2));
    }
    // stable
    else {
        psih = -(4.7 * ze) / 0.74;
    }
    if (psih < -4) psih = -4;
    if (psih > 3) psih = 3;
    return psih;
}
// **  Calculate diabatic influencing factor for heat ** //  
// [[Rcpp::export]]
double dphihCpp(double ze)
{
    double phih;
    // unstable
    if (ze < 0) {
        double phim = 1 / pow((1.0 - 16.0 * ze), 0.25);
        phih = pow(phim, 2.0);
    }
    // stable
    else {
        phih = 1 + ((6.0 * ze) / (1.0 + ze));
    }
    if (phih > 1.5) phih = 1.5;
    if (phih < 0.5) phih = 0.5;
    return phih;
}
// **  Calculate free convection ** //  
double gfreeCpp(double leafd, double H)
{
    double d = 0.71 * leafd;
    double dT = 0.7045388 * pow((d * pow(H, 4)), 0.2);
    double gha = 0.0375 * pow(dT / d, 0.25);
    if (gha < 0.1) gha = 0.1;
    return gha;
}
// **  Calculate molar conductance above canopy ** //  
double gturbCpp(double uf, double d, double zm, double zref, double ph, double psi_h, double gmin)
{
    double z0 = 0.2 * zm + d; // heat exchange surface height
    double ln = log((zref - d) / (z0 - d));
    double g = (0.4 * ph * uf) / (ln + psi_h);
    if (g < gmin) g = gmin;
    return g;
}
// **  Stomatal conductance ** //  
double stomcondCpp(double Rsw, double gsmax, double q50)
{
    double rpar = Rsw * 4.6;
    double gs = (gsmax * rpar) / (rpar + q50);
    return gs;
}
// **  Saturated vapour pressure ** //
// [[Rcpp::export]]  
double satvapCpp(double tc)
{
    double es;
    if (tc > 0) {
        es = 0.61078 * exp(17.27 * tc / (tc + 237.3));
    }
    else {
        es = 0.61078 * exp(21.875 * tc / (tc + 265.5));
    }
    return es;
}
// **  Dewpoint temperature ** // 
// [[Rcpp::export]]   
double dewpointCpp(double tc, double ea)
{
    double e0;
    double L;
    double it;
    // Dew point
    if (tc >= 0) {
        e0 = 611.2 / 1000;
        L = (2.501 * pow(10, 6)) - (2340 * tc);
        it = 1 / 273.15 - (461.5 / L) * log(ea / e0);
    }
    // Frost point
    else {
        e0 = 610.78 / 1000;
        L = 2.834 * pow(10, 6);
        it = 1 / 273.16 - (461.5 / L) * log(ea / e0);
    }
    double Tdew = 1 / it - 273.15;
    return Tdew;
}
// **  Penman-Monteith equation ** //
double PenmanMonteithCpp(double Rabs, double gHa, double gV, double tc, double te, double pk, double ea, double em, double G, double erh)
{
    double sb = 5.67 * pow(10, -8);
    double Rema = em * sb * pow(tc + 273.15, 4);
    double la = 0;
    if (te >= 0) {
        la = 45068.7 - 42.8428 * te;
    }
    else {
        la = 51078.69 - 4.338 * te - 0.06367 * te * te;
    }
    double cp = cpairCpp(te);
    double Da = satvapCpp(tc) - ea;
    double gR = (4 * em * sb * pow(te + 273.15, 3)) / cp;
    double De = satvapCpp(te + 0.5) - satvapCpp(te - 0.5);
    double Ts = tc + ((Rabs - Rema - la * (gV / pk) * Da * erh - G) / (cp * (gHa + gR) + la * (gV / pk) * De * erh));
    return Ts;
}
// **  Function to compute daily from hourly** // 
// [[Rcpp::export]]   
std::vector<double> hourtodayCpp(std::vector<double> hourly, std::string stat) {
    int numDays = hourly.size() / 24;
    std::vector<double> daily(numDays * 24, 0.0);
    for (int i = 0; i < numDays; ++i) {
        double dailyStat = hourly[i * 24];
        if (stat == "max") {
            for (int j = 1; j < 24; ++j) {
                dailyStat = std::max(dailyStat, hourly[i * 24 + j]);
            }
        }
        else if (stat == "min") {
            for (int j = 1; j < 24; ++j) {
                dailyStat = std::min(dailyStat, hourly[i * 24 + j]);
            }
        }
        else if (stat == "mean") {
            dailyStat = 0.0;
            for (int j = 0; j < 24; ++j) {
                dailyStat += hourly[i * 24 + j];
            }
            dailyStat /= 24;
        }
        else {
            throw std::invalid_argument("Invalid statistic. Please use 'max', 'min', or 'mean'.");
        }
        // Fill daily with replicated values for each day
        for (int j = 0; j < 24; ++j) {
            daily[i * 24 + j] = dailyStat;
        }
    }
    return daily;
}
// **  Function to compute daily from hourly without replicating each value 24 times** //  
std::vector<double> hourtodayCpp2(std::vector<double> hourly, std::string stat) {
    int numDays = hourly.size() / 24;
    std::vector<double> daily(numDays, 0.0);
    for (int i = 0; i < numDays; ++i) {
        double dailyStat = hourly[i * 24];
        if (stat == "max") {
            for (int j = 1; j < 24; ++j) {
                dailyStat = std::max(dailyStat, hourly[i * 24 + j]);
            }
        }
        else if (stat == "min") {
            for (int j = 1; j < 24; ++j) {
                dailyStat = std::min(dailyStat, hourly[i * 24 + j]);
            }
        }
        else if (stat == "mean") {
            dailyStat = 0.0;
            for (int j = 0; j < 24; ++j) {
                dailyStat += hourly[i * 24 + j];
            }
            dailyStat /= 24;
        }
        else if (stat == "sum") {
            dailyStat = 0.0;
            for (int j = 0; j < 24; ++j) {
                dailyStat += hourly[i * 24 + j];
            }
        }
        else {
            throw std::invalid_argument("Invalid statistic. Please use 'max', 'min', 'mean' or 'sum'.");
        }
        // Fill daily with replicated values for each day
        daily[i] = dailyStat;
    }
    return daily;
}
// **  Function to compute rolling mean temp ** // 
// [[Rcpp::export]] 
std::vector<double> maCpp(std::vector<double> x, int n) {
    std::vector<double> y(x.size());
    int m = x.size();
    for (int i = 0; i < m; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += x[(i - j + m) % m];
        }
        y[i] = sum / n;
    }
    return y;
}
// **  Function to compute rolling mean yearly ** // 
// [[Rcpp::export]]  
std::vector<double> mayCpp(std::vector<double> x) {
    // Calculate daily mean
    int numDays = x.size() / 24;
    std::vector<double> d(numDays, 0.0);
    for (int i = 0; i < numDays; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 24; ++j) {
            sum += x[i * 24 + j];
        }
        d[i] = sum / 24.0;
    }
    // Calculate rolling mean
    std::vector<double> y = maCpp(d, 91);
    std::vector<double> z;
    for (double val : y) {
        for (int i = 0; i < 24; ++i) {
            z.push_back(val);
        }
    }
    return z;
}
// **  Function to compute ground heat flux ** //  
Gmodel GFluxCpp(std::vector<double> Tg, std::vector<double> soilm, double rho, double Vm, double Vq, double Mc,
    std::vector<double> Gmax, std::vector<double> Gmin, int iter, bool yearG = true) {
    // Time invariant variables
    double frs = Vm + Vq;
    double c1 = (0.57 + 1.73 * Vq + 0.93 * Vm) / (1 - 0.74 * Vq - 0.49 * Vm) - 2.8 * frs * (1 - frs);
    double c3 = 1 + 2.6 * pow(Mc, -0.5);
    double c4 = 0.03 + 0.7 * frs * frs;
    double mu1 = 2400 * rho / 2.64;
    double mu2 = 1.06 * rho;
    // Calculate daily mean soil surface temperature
    std::vector<double> Td = hourtodayCpp(Tg, "mean");
    // Initalise variables that need retaining
    std::vector<double> Gmu(Tg.size());
    std::vector<double> dT(Tg.size());
    std::vector<double> k(Tg.size());
    std::vector<double> ka(Tg.size());
    for (size_t i = 0; i < Tg.size(); ++i) {
        // Find soil diffusivity
        double cs = mu1 + 4180 * soilm[i];
        double ph = (rho * (1 - soilm[i]) + soilm[i]) * 1000;
        double c2 = mu2 * soilm[i];
        k[i] = c1 + c2 * soilm[i] - (c1 - c4) * exp(-pow(c3 * soilm[i], 4));
        ka[i] = k[i] / (cs * ph);
        double omdy = (2 * M_PI) / (24 * 3600);
        double DD = sqrt(2 * ka[i] / omdy);
        Gmu[i] = sqrt(2) * (k[i] / DD) * 0.5;
        // Calculate T fluctuation from daily mean
        dT[i] = Tg[i] - Td[i];
    }
    // Calculate 6 hour back rolling mean of Gmu and dT to get 3 hour lag
    std::vector<double> Gmud = maCpp(Gmu, 6);
    std::vector<double> G = maCpp(dT, 6);
    // remove this line **
    for (size_t i = 0; i < G.size(); ++i) G[i] = G[i] * Gmud[i] * 1.1171;
    if (iter == 0) {
        Gmin = hourtodayCpp(G, "min");
        Gmax = hourtodayCpp(G, "max");
    }
    for (size_t i = 0; i < G.size(); ++i) {
        if (G[i] < Gmin[i]) G[i] = Gmin[i];
        if (G[i] > Gmax[i]) G[i] = Gmax[i];
    }
    // Calculate yearly ground flux
    if (yearG) {
        // Calculate moving average of k
        std::vector<double> kma = mayCpp(k);
        std::vector<double> kama = mayCpp(ka);
        std::vector<double> Gmuy(Tg.size());
        for (size_t i = 0; i < G.size(); ++i) {
            double omyr = (2 * M_PI) / (G.size() * 3600);
            Gmuy[i] = sqrt(2) * kma[i] / sqrt(2 * kama[i] / omyr);
        }
        // Calculate mean Td
        double sumTd = 0.0;
        for (double num : Td) sumTd += num;
        std::vector<double> dTy(Tg.size());
        for (size_t i = 0; i < dTy.size(); ++i) dTy[i] = Td[i] - sumTd / Td.size();
        std::vector<double> madTy = mayCpp(dTy);
        for (size_t i = 0; i < G.size(); ++i) {
            G[i] = G[i] + madTy[i] * Gmuy[i] * 1.1171;
        }
    }
    Gmodel out;
    out.G = G;
    out.Gmin = Gmin;
    out.Gmax = Gmax;
    return out;
}
// Run point model in progress
// [[Rcpp::export]]
Rcpp::List BigLeafCpp(DataFrame obstime, DataFrame climdata, std::vector<double> vegp, std::vector<double> groundp, std::vector<double> soilm,
    double lat, double lon, double dTmx = 25.0, double zref = 2.0, int maxiter = 100, double bwgt = 0.5, double tol = 0.5, double gmn = 0.1, bool yearG = true)
{
    // Access items of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access columns of climdata
    std::vector<double> tc = climdata["temp"];
    std::vector<double> rh = climdata["relhum"];
    std::vector<double> pk = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> wspeed = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> prec = climdata["precip"];
    // Access items of vegp
    double h = vegp[0];
    double pai = vegp[1];
    double vegx = vegp[2];
    double clump = vegp[3];
    double lref = vegp[4];
    double ltra = vegp[5];
    double leafd = vegp[6];
    double em = vegp[7];
    double gsmax = vegp[8];
    double q50 = vegp[9];
    // Access items of groundp
    double gref = groundp[0];
    double slope = groundp[1];
    double aspect = groundp[2];
    double groundem = groundp[3];
    double rho = groundp[4];
    double Vm = groundp[5];
    double Vq = groundp[6];
    double Mc = groundp[7];
    double Smin = groundp[11];
    double Smax = groundp[10];
    // Calculate SW radiation
    radmodel swabs = RadswabsCpp(pai, vegx, lref, ltra, clump, gref, slope, aspect, lat, lon, year, month, day, hour, Rsw, Rdif);
    // Calculate time-invarient variables
    double pait = pai / (1 - clump);
    double trd = (1 - clump * clump) * exp(-pait) + clump * clump;
    double sb = 5.67 * pow(10, -8);
    double d = zeroplanedisCpp(h, pai);
    // Used to avoid (h-d)/zm being less than one, meaning log((h-d)/zm) becomes negative
    double Belim = 0.4 / sqrt(0.003 + (0.2 * pai) / 2);
    // Initialise variables set to zero on first run
    std::vector<double> Tg = tc;
    std::vector<double> Tc = tc;
    std::vector<double> tcc = tc;
    std::vector<double> tcg = tc;
    std::vector<double> psim(Rsw.size(), 0.0);
    std::vector<double> psih(Rsw.size(), 0.0);
    std::vector<double> phih(Rsw.size(), 0.0);
    std::vector<double> LL(Rsw.size(), 0.0);
    std::vector<double> G(Rsw.size(), 0.0);
    std::vector<double> Gmin(Rsw.size(), -999.0);
    std::vector<double> Gmax(Rsw.size(), 999.0);
    std::vector<double> uf(Rsw.size(), 999.0);
    std::vector<double> RabsG(Rsw.size(), 999.0);
    // New variables for storing
    // Initalise H
    std::vector<double> H(Rsw.size());
    for (size_t i = 0; i < H.size(); ++i) H[i] = 0.5 * Rsw[i] - em * sb * pow(tc[i] + 273.15, 4);
    double tstf = tol * 2;
    int iter = 0;
    double tst = 0;
    while (tstf > tol) {
        tst = 0;
        for (size_t i = 0; i < Rsw.size(); ++i) {
            // Calculate longwave radiation
            double RemC = em * sb * pow(Tc[i] + 273.15, 4);
            double radClw = em * Rlw[i]; // Longwave radiation down from sky
            double radGlw = groundem * (trd * radClw + (1 - trd) * RemC);
            // Calculate absorbed radiation
            RabsG[i] = swabs.ground[i] + radGlw;
            double RabsC = swabs.canopy[i] + radClw;
            // Calculate canopy temperature
            double zm = roughlengthCpp(h, pai, d, psih[i]);
            uf[i] = (0.4 * wspeed[i]) / (log((zref - d) / zm) + psim[i]);
            if (uf[i] < 0.0002) uf[i] = 0.0002;
            double gmin = gfreeCpp(leafd, abs(H[i])) * 2 * pai;
            double ph = phairCpp(tcc[i], pk[i]);
            double gHa = gturbCpp(uf[i], d, zm, zref, ph, psih[i], gmin);
            double gC = stomcondCpp(Rsw[i], gsmax * 3, q50 * 3);
            double gV = 1 / (1 / gHa + 1 / gC);
            if (gC == 0) gV = 0;
            double ea = satvapCpp(tc[i]) * rh[i] / 100;
            double Tcn = PenmanMonteithCpp(RabsC, gHa, gV, tc[i], tcc[i], pk[i], ea, em, G[i], 1);
            double tdew = dewpointCpp(tc[i], ea);
            if (Tcn < tdew) Tcn = tdew;
            // Calculate ground surface temperature
            double srh = (soilm[i] - Smin) / (Smax - Smin);
            double Tgn = PenmanMonteithCpp(RabsG[i], gHa, gHa, tcg[i], tcc[i], pk[i], ea, em, G[i], srh);
            if (Tgn < tdew) Tgn = tdew;
            // Cap values
            double dTc = Tcn - tc[i];
            double dTg = Tgn - tc[i];
            if (dTc > dTmx) dTc = dTmx;
            if (dTg > dTmx) dTg = dTmx;
            Tcn = tc[i] + dTc;
            Tgn = tc[i] + dTg;
            // Run tests for convergence
            double tst2 = abs(Tcn - Tc[i]);
            double tst3 = abs(Tgn - Tg[i]);
            if (tst2 > tst) tst = tst2;
            if (tst3 > tst) tst = tst3;
            // Reassign Tc and Tg using bwgt
            Tc[i] = bwgt * Tc[i] + (1 - bwgt) * Tcn;
            Tg[i] = bwgt * Tg[i] + (1 - bwgt) * Tgn;
            // Recalculate variables
            tcc[i] = (Tc[i] + tc[i]) / 2;
            tcg[i] = (Tg[i] + tc[i]) / 2;
            double Tk = 273.15 + tcc[i];
            ph = phairCpp(tcc[i], pk[i]);
            double cp = cpairCpp(tcc[i]);
            // Calculate H
            H[i] = bwgt * H[i] + (1 - bwgt) * (gHa * cp * (Tcn - tc[i]));
            // Set limits to H
            double Rnet = RabsC - sb * em * pow(Tc[i] + 273.15, 4);
            if (Rnet > 0 && H[i] > Rnet) H[i] = Rnet;
            // Recalculate stablity variables
            // Stability
            LL[i] = (ph * cp * pow(uf[i], 3) * Tk) / (-0.4 * 9.81 * H[i]);
            //if (LL[i] > 10000.0) LL[i] = 10000.0;
            //if (LL[i] < -10000.0) LL[i] = -10000.0;
            psim[i] = dpsimCpp(zm / LL[i]) - dpsimCpp((zref - d) / LL[i]);
            psih[i] = dpsihCpp((0.2 * zm) / LL[i]) - dpsihCpp((zref - d) / LL[i]);
            phih[i] = dphihCpp((zref - d) / LL[i]);
            // Set limits to diabatic coefficients
            double ln1 = log((zref - d) / zm);
            double ln2 = log((zref - d) / (0.2 * zm));
            if (psim[i] < -0.9 * ln1) psim[i] = -0.9 * ln1;
            if (psih[i] < -0.9 * ln2) psih[i] = -0.9 * ln2;
            if (psim[i] > 0.9 * ln1) psim[i] = 0.9 * ln1;
            if (psih[i] > 0.9 * ln2) psih[i] = 0.9 * ln2;
            if (psih[i] > 0.9 * Belim) psih[i] = 0.9 * Belim;
        }
        // Recalculate Ground heat flux
        Gmodel GG = GFluxCpp(Tg, soilm, rho, Vm, Vq, Mc, Gmax, Gmin, iter, yearG);
        G = GG.G;
        Gmin = GG.Gmin;
        Gmax = GG.Gmax;
        tstf = tst;
        ++iter;
        if (iter >= maxiter) tstf = 0;
    }
    // Return outputs
    Rcpp::List out;
    out["Tc"] = Rcpp::wrap(Tc);
    out["Tg"] = Rcpp::wrap(Tg);
    out["H"] = Rcpp::wrap(H);
    out["G"] = Rcpp::wrap(G);
    out["psih"] = Rcpp::wrap(psih);
    out["psim"] = Rcpp::wrap(psim);
    out["phih"] = Rcpp::wrap(phih);
    out["OL"] = Rcpp::wrap(LL);
    out["uf"] = Rcpp::wrap(uf);
    out["RabsG"] = Rcpp::wrap(RabsG);
    out["err"] = Rcpp::wrap(tst);
    out["albedo"] = Rcpp::wrap(swabs.albedo);
    return out;
}
// Perform weather height adjustment
// [[Rcpp::export]]
DataFrame weatherhgtCpp(DataFrame obstime, DataFrame climdata, double zin, double uzin, double zout, double lat, double lon)
{
    std::vector<double> tc = climdata["temp"];
    std::vector<double> rh = climdata["relhum"];
    std::vector<double> ws = climdata["windspeed"];
    std::vector<double> vegp({ 0.12,1,1,0.1,0.4,0.2,0.05,0.97,0.33,100.0 });
    std::vector<double> groundp({ 0.15,0.0, 180.0,0.97,1.529643,0.509,0.06,0.5422,5.2,-5.6,0.419,0.074 });
    std::vector<double> soilm(tc.size(), 0.2); // Initialize soilm with size and value 0.2
    // Run point model
    Rcpp::List bigleafp = BigLeafCpp(obstime, climdata, vegp, groundp, soilm, lat, lon, 25, 2, 20, 0.5, 0.5, 0.1);
    // Extract things needed from list
    std::vector<double> psih = Rcpp::as<std::vector<double>>(bigleafp["psih"]);
    std::vector<double> Tc = Rcpp::as<std::vector<double>>(bigleafp["Tc"]);
    // Define variables
    std::vector<double> Tz(tc.size());
    std::vector<double> Rh(tc.size());
    std::vector<double> Uz(tc.size());
    // temporary variables to diagnose
    double d = zeroplanedisCpp(0.12, 1);
    for (size_t i = 0; i < tc.size(); ++i) {
        double zm = roughlengthCpp(0.12, 1, d, psih[i]);
        double zh = 0.2 * zm;
        double lnr = log((zout - d) / zh) / log((zin - d) / zh);
        // Temperature
        Tz[i] = (Tc[i] - tc[i]) * (1 - lnr) + tc[i];
        // Humidity
        double ea = satvapCpp(tc[i]) * rh[i] / 100;
        double es = satvapCpp(Tc[i]) * sqrt(rh[i] / 100);
        double ez = ea + (es - ea) * (1 - lnr);
        es = satvapCpp(Tz[i]);
        Rh[i] = (ez / es) * 100;
        // Cap Rh
        if (Rh[i] < 0.25 * rh[i]) Rh[i] = 0.25 * rh[i];
        if (Rh[i] > 100.0) Rh[i] = 100.0;
        // Wind speed
        double lnru = log((zout - d) / zm) / log((uzin - d) / zm);
        Uz[i] = ws[i] * lnru;
    }
    // Clone climdata to stop original input being over-written in R
    DataFrame climdata_copy = clone(climdata); // Make a copy of climdata
    climdata_copy["temp"] = Tz;
    climdata_copy["relhum"] = Rh;
    climdata_copy["windspeed"] = Uz;
    return climdata_copy;
}
// Soil moisture model
// [[Rcpp::export]]
std::vector<double> soilmCpp(DataFrame climdata, double rmu, double mult, double pwr, double Smax, double Smin, double Ksat, double a)
{
    // Extract data from columns
    std::vector<double> temp = climdata["temp"];
    std::vector<double> swdown = climdata["swdown"];
    std::vector<double> lwdown = climdata["lwdown"];
    std::vector<double> rainh = climdata["precip"];
    // Calculate net radiation
    std::vector<double> rnet(temp.size());
    for (size_t i = 0; i < temp.size(); ++i) {
        double swrad = (1 - 0.15) * swdown[i];
        double lwout = 5.67 * pow(10, -8) * 0.95 * pow(temp[i] + 273.15, 4);
        double lwnet = lwout - lwdown[i];
        rnet[i] = swrad - lwnet;
        if (rnet[i] < 0) rnet[i] = 0;
    }
    // Convert to daily
    std::vector<double> rnetd = hourtodayCpp2(rnet, "mean");
    std::vector<double> rain = hourtodayCpp2(rainh, "sum");
    // Run soil model
    std::vector<double> soilm(rain.size());
    double s1 = Smax;
    double s2 = Smax;
    soilm[0] = Smax;
    for (size_t i = 1; i < rain.size(); ++i) {
        double sav = (s1 + s2) / 2;
        double dif = s2 - s1;
        s1 = s1 + rmu * rain[i] - mult * rnetd[i];
        double k = Ksat * pow(sav / Smax, pwr);
        s1 = s1 + a * k * dif;
        s2 = s2 - ((a * k * dif) / 10);
        if (s1 > Smax) s1 = Smax;
        if (s2 > Smax) s2 = Smax;
        if (s1 < Smin) s1 = Smin;
        if (s2 < Smin) s2 = Smin;
        soilm[i] = (s1 + s2) / 2;
    }
    return soilm;
}
// Bits needed for microclimf model
// [[Rcpp::export]]
DataFrame pointmprocess(DataFrame pointvars, double zref, double h, double pai, 
    double rho, double Vm, double Vq, double Mc)
{
    // Extract required variables form point model
    std::vector<double> u2 = pointvars["windspeed"];
    std::vector<double> tc = pointvars["tc"];
    std::vector<double> rh = pointvars["rh"];
    std::vector<double> pk = pointvars["pk"];
    std::vector<double> uf = pointvars["uf"];
    std::vector<double> soilm = pointvars["soilm"];
    std::vector<double> RabsG = pointvars["RabsG"];
    // Calculate time invariant variables
    double dp = zeroplanedisCpp(h, pai);
    double zmp = roughlengthCpp(h, pai, dp, 0);
    double frs = Vm + Vq;
    double c1 = (0.57 + 1.73 * Vq + 0.93 * Vm) / (1 - 0.74 * Vq - 0.49 * Vm) - 2.8 * frs * (1 - frs);
    double c3 = 1 + 2.6 * pow(Mc, -0.5);
    double c4 = 0.03 + 0.7 * pow(frs, 2);
    // Create variables
    std::vector<double> umu(u2.size());
    std::vector<double> kp(u2.size());
    std::vector<double> DDp(u2.size());
    std::vector<double> muGp(u2.size());
    std::vector<double> dtrp(u2.size());
    std::vector<double> T0p(u2.size());
    for (size_t i = 0; i < u2.size(); ++i) {
        // Calculate umu
        double ufps = (0.4 * u2[i]) / log((zref - dp) / zmp);
        umu[i] = uf[i] / ufps;
        // Calculate conductivity
        double cs = (2400 * rho / 2.64 + 4180 * soilm[i]); // specific heat of soil in J / kg / K
        double ph = (rho * (1 - soilm[i]) + soilm[i]) * 1000; // bulk density in kg / m3
        double c2 = 1.06 * rho * soilm[i];
        kp[i] = c1 + c2 * soilm[i] - (c1 - c4) * exp(-pow(c3 * soilm[i], 4));
        // Calculate muGp
        double ka = kp[i] / (cs * ph);
        double omdy = (2 * M_PI) / (24 * 3600);
        muGp[i] = pow(2 * ka / omdy, 0.5);
        // Calculate DDp
        DDp[i] = pow(2.0 * ka / omdy, 0.5);
        // Calculate T0p
        double gHa = (0.4 * 43.0 * ufps) / log((zref - dp) / zmp);
        double es = satvapCpp(tc[i]);
        double ea = es * rh[i] / 100.0;
        T0p[i] = PenmanMonteithCpp(RabsG[i], gHa, gHa, tc[i], tc[i], pk[i], ea, 0.97, 0.0, 1.0);
    }
    // Calculate dtrp
    std::vector<double> T0pmx = hourtodayCpp(T0p, "max");
    std::vector<double> T0pmn = hourtodayCpp(T0p, "min");
    for (size_t i = 0; i < T0pmx.size(); ++i) {
        dtrp[i] = T0pmx[i] - T0pmn[i];
    }
    DataFrame out;
    out["umu"] = umu;
    out["kp"] = kp;
    out["muGp"] = muGp;
    out["DDp"] = DDp;
    out["T0p"] = T0p;
    out["dtrp"] = dtrp;
    return out;
}
// ********************************************************************* //
// ~~~~~~~~~~~~~~~~~~ microclimf model from here ~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************************************************** //
// Worker functions
// Function to permute dimensions of a 3D array
NumericVector aperm3D(NumericVector Tz, int rows, int cols, int tsteps) {
    // Initialize a new vector to store the permuted array
    NumericVector permuted(Tz.size());

    // Permute the dimensions
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < tsteps; ++k) {
                // Calculate the index in the original array
                int index_original = k + tsteps * (j + cols * i);
                // Calculate the index in the permuted array
                int index_permuted = i + rows * (j + cols * k);
                // Copy the element from the original array to the permuted array
                permuted[index_permuted] = Tz[index_original];
            }
        }
    }
    // Set the dimensions of the permuted array
    permuted.attr("dim") = NumericVector::create(rows, cols, tsteps);
    return permuted;
}
// Function to permute dimensions of a 3D array and convert to a NumericVector (the other way)
// [[Rcpp::export]]
NumericVector aperm3D2(NumericVector Tz, int rows, int cols, int tsteps)
{
    // Initialize a new vector to store the permuted array
    NumericVector permuted(Tz.size());
    // Permute the dimensions
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int k = 0; k < tsteps; ++k) {
                // Calculate the index in the original array
                int index_original = i + rows * (j + cols * k);
                // Calculate the index in the permuted array
                int index_permuted = k + tsteps * (j + cols * i);
                // Copy the element from the original array to the permuted array
                permuted[index_permuted] = Tz[index_original];
            }
        }
    }
    // Set the dimensions of the permuted array
    return permuted;
}
// Function to slice a 3D array along the 3rd dimension for a fixed k (2D slice)
// [[Rcpp::export]]
NumericVector slice_2d(const NumericVector& a, int k, IntegerVector dim) {
    int nrow = dim[0];
    int ncol = dim[1];
    int total_elements = nrow * ncol;
    NumericVector result(total_elements);
    int index = 0;
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            double value = a[i + j * nrow + k * nrow * ncol];
            if (!NumericVector::is_na(value)) {
                result[index++] = value;
            }
        }
    }
    result = result[seq(0, index - 1)];  // Trim the result to remove unused elements
    return result;
}
// Apply the function (mean, sum, min, max) over the third dimension
// [[Rcpp::export]]
NumericVector applycpp3(NumericVector a, std::string fun_name) {
    IntegerVector dim = a.attr("dim");
    int nlayer = dim[2];  // Third dimension (number of slices)
    NumericVector result(nlayer);  // Result vector with size equal to number of slices
    for (int k = 0; k < nlayer; k++) {
        NumericVector slice = slice_2d(a, k, dim);
        if (slice.size() == 0) {
            // If the cleaned slice is empty (i.e., all values were NA), assign NA to the result
            result[k] = NA_REAL;
        }
        else {
            // Apply the desired function on the cleaned vector
            if (fun_name == "mean") {
                result[k] = mean(slice);
            }
            else if (fun_name == "sum") {
                result[k] = sum(slice);
            }
            else if (fun_name == "min") {
                result[k] = min(slice);
            }
            else if (fun_name == "max") {
                result[k] = max(slice);
            }
            else {
                stop("Unknown function name");
            }
        }
    }
    return result;
}
// Converts NumericVector to std::vector<double>
std::vector<double> tosvd(NumericVector x)
{
    std::vector<double> y(x.size());
    for (int i = 0; i < x.size(); ++i) y[i] = x[i];
    return y;
}
// **  Function to compute rolling mean yearly ** // 
// [[Rcpp::export]]  
std::vector<double> manCpp(std::vector<double> x, int n) {
    std::vector<double> z(x.size());
    // Calculate rolling mean if n < 48
    if (n <= 48) {
        z = maCpp(x, n);
    }
    // Average to daily first if n > 48
    else {
        // Calculate daily mean
        int numDays = x.size() / 24;
        std::vector<double> d(numDays, 0.0);
        for (int i = 0; i < numDays; ++i) {
            double sum = 0.0;
            for (int j = 0; j < 24; ++j) {
                sum += x[i * 24 + j];
            }
            d[i] = sum / 24.0;
        }
        int n2 = n / 24;
        // Calculate rolling mean of daily
        std::vector<double> y = maCpp(d, n2);
        // Expand back out to hourly
        for (int i = 0; i < numDays; ++i) {
            for (int j = 0; j < 24; ++j) {
                z[i * 24 + j] = y[i];
            }
        }
        z = maCpp(z, 24);
    }
    return z;
}
NumericMatrix subsetArray(NumericVector a, int layer) {
    // Extract dimensions
    IntegerVector dims = a.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int layers = dims[2];
    // Check if the layer is within bounds
    if (layer >= layers) {
        stop("Layer index out of bounds");
    }
    // Create the matrix to store the desired layer
    NumericMatrix layer_mat(rows, cols);
    // Extract the desired layer
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Corrected index calculation
            layer_mat(i, j) = a[layer * rows * cols + i + j * rows];
        }
    }
    return layer_mat;
}
int getn(NumericVector a)
{
    // Extract dimensions
    SEXP dims = a.attr("dim");
    int out = INTEGER(dims)[2];
    return out;
}
DataFrame ssdf(DataFrame df, IntegerVector rows) {
    int n = rows.size();
    DataFrame output = DataFrame::create(); // Create an empty DataFrame to populate
    CharacterVector names = df.names();     // Get names to maintain structure
    // Loop over each column
    for (int i = 0; i < df.size(); i++) {
        SEXP column = df[i];
        switch (TYPEOF(column)) {
        case INTSXP: {
            IntegerVector col = as<IntegerVector>(column);
            IntegerVector sub(n);
            for (int j = 0; j < n; j++) {
                sub[j] = col[rows[j]]; // Use 0-based index directly
            }
            output.push_back(sub, as<std::string>(names[i]));
            break;
        }
        case REALSXP: {
            NumericVector col = as<NumericVector>(column);
            NumericVector sub(n);
            for (int j = 0; j < n; j++) {
                sub[j] = col[rows[j]]; // Use 0-based index directly
            }
            output.push_back(sub, as<std::string>(names[i]));
            break;
        }
        case STRSXP: {
            StringVector col = as<StringVector>(column);
            StringVector sub(n);
            for (int j = 0; j < n; j++) {
                sub[j] = col[rows[j]]; // Use 0-based index directly
            }
            output.push_back(sub, as<std::string>(names[i]));
            break;
        }
        default:
            stop("Unhandled column type in DataFrame");
        }
    }
    return output;
}
NumericVector abind3D(NumericVector a1, NumericVector a2) {
    // Extract dimensions of a1
    SEXP dims1 = a1.attr("dim");
    int rows1 = INTEGER(dims1)[0];
    int cols1 = INTEGER(dims1)[1];
    int nlyr1 = INTEGER(dims1)[2];
    // Extract dimensions of a2
    SEXP dims2 = a2.attr("dim");
    int nlyr2 = INTEGER(dims2)[2];
    NumericVector abind(rows1 * cols1 * (nlyr1 + nlyr2));
    abind.attr("dim") = Dimension(rows1, cols1, nlyr1 + nlyr2);
    // Copy elements from a1 and a2 to abind
    std::copy(a1.begin(), a1.end(), abind.begin());
    std::copy(a2.begin(), a2.end(), abind.begin() + nlyr1 * rows1 * cols1);
    return abind;
}
NumericVector subset3DArrayBySlices(NumericVector array, IntegerVector slices) {
    // Ensure that array has dimensions attribute
    if (!array.hasAttribute("dim")) {
        stop("Input array does not have dimensions.");
    }
    // Get dimensions of the original array
    IntegerVector dims = array.attr("dim");
    if (dims.size() != 3) {
        stop("Input array is not a 3D array.");
    }
    int n1 = dims[0];
    int n2 = dims[1];
    int n3 = dims[2];
    // Exclude any invalid slice indices
    int numSlices = slices.size();
    for (int i = 0; i < numSlices; ++i) {
        if (slices[i] < 0 || slices[i] >= n3) {
            stop("Slice index out of bounds.");
        }
    }
    // Dimensions of the resulting array
    NumericVector result(n1 * n2 * numSlices);
    result.attr("dim") = IntegerVector::create(n1, n2, numSlices);
    // Populate the resulting array
    for (int k = 0; k < numSlices; k++) {
        int sliceIndex = slices[k];
        for (int j = 0; j < n2; j++) {
            for (int i = 0; i < n1; i++) {
                result[i + n1 * (j + n2 * k)] = array[i + n1 * (j + n2 * sliceIndex)];
            }
        }
    }
    return result;
}
NumericVector slicea(NumericVector array, IntegerVector slices) {
    // Ensure that array has dimensions attribute
    if (!array.hasAttribute("dim")) {
        stop("Input array does not have dimensions.");
    }
    // Get dimensions of the original array
    IntegerVector dims = array.attr("dim");
    if (dims.size() != 3) {
        stop("Input array is not a 3D array.");
    }
    int n1 = dims[0];
    int n2 = dims[1];
    int n3 = dims[2];
    // Exclude any invalid slice indices
    int numSlices = slices.size();
    for (int i = 0; i < numSlices; ++i) {
        if (slices[i] < 0 || slices[i] >= n3) {
            stop("Error slice index out of bounds.");
        }
    }
    // Dimensions of the resulting array
    NumericVector result(n1 * n2 * numSlices);
    result.attr("dim") = IntegerVector::create(n1, n2, numSlices);
    // Populate the resulting array
    for (int k = 0; k < numSlices; k++) {
        int sliceIndex = slices[k];
        for (int j = 0; j < n2; j++) {
            for (int i = 0; i < n1; i++) {
                result[i + n1 * (j + n2 * k)] = array[i + n1 * (j + n2 * sliceIndex)];
            }
        }
    }
    return result;
}
// ************************************************************* //
// ~ Functions needed for calaculated topographic wetness index
// ************************************************************* //
// [[Rcpp::export]]
IntegerMatrix flowdirCpp(NumericMatrix md) {
    int nrow = md.nrow();
    int ncol = md.ncol();
    // Create a padded matrix with NA around the edges
    NumericMatrix md2(nrow + 2, ncol + 2);
    std::fill(md2.begin(), md2.end(), NA_REAL);
    // Copy original matrix into the center of md2
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            md2(i + 1, j + 1) = md(i, j);
        }
    }
    // Initialize outputs
    IntegerMatrix fd(nrow, ncol);
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            // Check whether focal cell is NA
            double val = md(i, j);
            if (!R_IsNA(val)) {
                double minval = 9999.99;
                int indx = 1;
                for (int jj = 0; jj < 3; ++jj) {
                    for (int ii = 0; ii < 3; ++ii) {
                        double val2 = md2(i + ii, j + jj);
                        if (!R_IsNA(val2)) {
                            if (val2 < minval) {
                                minval = val2;
                                fd(i, j) = indx;
                            }
                        }
                        ++indx;
                    } // end ii
                } // end jj
            } // end NA check
            else {
                fd(i, j) = NA_REAL;
            }
        }
    }
    return fd;
}
// [[Rcpp::export]]
NumericMatrix flowaccCpp(NumericMatrix dm) {
    int nrow = dm.nrow();
    int ncol = dm.ncol();
    // Get flow direction
    IntegerMatrix fd = flowdirCpp(dm);
    // Initialize flow accumulation matrix with 1s
    NumericMatrix fa(nrow, ncol);
    std::fill(fa.begin(), fa.end(), 1); // Set all elements to 1 initially
    // Set NA values in `fa` where `dm` is NA
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            if (R_IsNA(dm(i, j))) {
                fa(i, j) = NA_INTEGER;  // Propagate NA to output
            }
        }
    }
    // Create an index vector for ordering `dm` in decreasing order
    std::vector<std::pair<double, int>> orderVec;
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            if (!R_IsNA(dm(i, j))) {
                orderVec.push_back({ dm(i, j), i * ncol + j }); // Store value and 1D index
            }
        }
    }
    // Sort in decreasing order
    std::sort(orderVec.begin(), orderVec.end(), std::greater<std::pair<double, int>>());
    // Process ordered indices
    for (size_t i = 0; i < orderVec.size() - 1; i++) {
        int index = orderVec[i].second;
        int y = index / ncol; // Row index
        int x = index % ncol; // Column index
        if (R_IsNA(dm(y, x))) continue; // Skip NA cells
        int f = fd(y, x); // Flow direction from fd
        if (f < 1 || f > 9) continue; // Skip invalid flow direction values
        // Compute new coordinates based on flow direction (1-based indexing)
        int y2 = y + (f - 1) % 3 - 1;
        int x2 = x + (f - 1) / 3 - 1;
        // Ensure the new coordinates are within bounds
        if (x2 >= 0 && x2 < ncol && y2 >= 0 && y2 < nrow) {
            if (fa(y2, x2) != NA_INTEGER) { // Only update if the target is not NA
                fa(y2, x2) += fa(y, x);
            }
        }
    }
    return fa;
}
// ************************************************************* //
// ~~~~~~~~~~~~~~~  Grid model from here ~~~~~~~~~~~~~~~~~~~~~~~
// ************************************************************* //
// soilmdistribute code (vector)
// [[Rcpp::export]]
std::vector<double> soildCppv(std::vector<double> soilm, double Smin, double Smax, double tadd)
{
    double rge = Smax - Smin;
    std::vector<double> smout(soilm.size());
    for (size_t i = 0; i < soilm.size(); ++i) {
        double theta = (soilm[i] - Smin) / rge;
        if (theta > 0.9999) theta = 0.9999;
        if (theta < 0.0001) theta = 0.0001;
        double lt = log(theta / (1 - theta));
        double sm = lt + tadd;
        sm = 1 / (1 + exp(-sm));
        smout[i] = sm * rge + Smin;
    }
    return smout;
}
// soilmdistribute code(matrix - i.e.returns tadd)
// [[Rcpp::export]]
NumericMatrix soildCppm(NumericMatrix twi, double Smin, double Smax, double tfact)
{
    int rows = twi.nrow();
    int cols = twi.ncol();
    // Calculate ltwi
    NumericMatrix ltwi(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = twi(i, j);
            if (!NumericMatrix::is_na(val)) {
                ltwi(i, j) = log(val) / tfact;
            }
            else {
                ltwi(i, j) = NumericVector::get_na(); // Set NA if twi is NA
            }
        }
    }
    // Calculate me
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = ltwi(i, j);
            if (!NumericMatrix::is_na(val)) {
                sum += val;
                count++;
            }
        }
    }
    double me = sum / count;
    // Calculate tadd
    NumericMatrix tadd(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = ltwi(i, j);
            if (!NumericMatrix::is_na(val)) {
                tadd(i, j) = val - me;
            }
            else {
                tadd(i, j) = NumericVector::get_na();
            }
        }
    }
    return tadd;
}
// ** Calculate two-stream coefficients vector ** //
radmodel2 twostreamvCpp(double reqhgt, double hgt, double pai, double paia, double x, double lref, double ltra, double clump, double gref,
    double svfa, std::vector<double> si, std::vector<double> zen, std::vector<double> Rsw, std::vector<double> Rdif,
    std::vector<double> tc, std::vector<double> lwdown)
{
    // ** Define variables
    std::vector<double> radGsw(Rsw.size());
    std::vector<double> radCsw(Rsw.size());
    std::vector<double> radGlw(Rsw.size());
    std::vector<double> radClw(Rsw.size());
    std::vector<double> Rbdown(Rsw.size());
    std::vector<double> Rddown(Rsw.size());
    std::vector<double> Rdup(Rsw.size());
    std::vector<double> radLsw(Rsw.size());
    std::vector<double> lwout(Rsw.size());
    if (pai == 0.0) {
        for (size_t i = 0; i < Rsw.size(); ++i) {
            Rbdown[i] = (Rsw[i] - Rdif[i]) / cos(zen[i]);
            radGsw[i] = (1 - gref) * (svfa * Rdif[i] + si[i] * Rbdown[i]);
            radCsw[i] = radGsw[i];
            radGlw[i] = 0.97 * svfa * lwdown[i];
            radClw[i] = radGlw[i];
            Rddown[i] = Rdif[i];
            Rdup[i] = gref * Rsw[i];
            radLsw[i] = 0;
            lwout[i] = 0.97 * 5.67 * pow(10.0, -8.0) * pow(tc[i] + 273.15, 4.0);
        }
    }
    else {
        // Calculate time-invariant variables
        // ** Perform adjustment to total leaf area for clump
        double pait = pai / (1 - clump);
        // ** Calculate base parameters
        double om = lref + ltra;
        double a = 1 - om;
        double del = lref - ltra;
        double J = 1.0 / 3.0;
        if (x != 1.0) {
            double mla = 9.65 * pow((3 + x), -1.65);
            if (mla > M_PI / 2) mla = M_PI / 2;
            J = cos(mla) * cos(mla);
        }
        double gma = 0.5 * (om + J * del);
        double h = sqrt(a * a + 2 * a * gma);
        // Calculate two-stream parameters (diffuse)
        std::vector<double> tspdif = twostreamdifCpp(pait, om, a, gma, h, gref);
        double p1 = tspdif[0];
        double p2 = tspdif[1];
        double p3 = tspdif[2];
        double p4 = tspdif[3];
        double u1 = tspdif[4];
        double S1 = tspdif[5];
        double D1 = tspdif[6];
        double D2 = tspdif[7];
        // ** calculate tranmissions through gaps
        double trdn = pow(clump, 2.0); // downward at ground
        double gi = pow(clump, paia / pai);
        double giu = pow(clump, (pai - paia) / pai);
        if (gi > 0.99) gi = 0.99;
        if (giu > 0.99) giu = 0.99;
        double trd = gi * gi; // diffuse transmission downward
        double trdu = giu * giu; // diffuse transmission upward
        double paiaa = paia / (1.0 - gi); // paia adjusted for gap fraction
        // ** white-sky albedo
        double amx = gref;
        if (amx < lref) amx = lref;
        double albd = (1.0 - trdn) * (p1 + p2) + trdn * gref;
        if (albd > amx) albd = amx;
        if (albd < 0.01) albd = 0.01;
        // ** normalised downward diffuse flux at ground
        double Rddn_g = (1.0 - trdn) * (p3 * exp(-h * pait) + p4 * exp(h * pait)) + trdn;
        if (Rddn_g > 1.0) Rddn_g = 1.0;
        if (Rddn_g < 0.0) Rddn_g = 0.0;
        // ** normalised upward diffuse flux at z
        double Rdup_z = (1.0 - trdu) * (p1 * exp(-h * paiaa) + p2 * exp(h * paiaa)) 
            + trdu * gref;
        if (Rdup_z > 1.0) Rdup_z = 1.0;
        if (Rdup_z < 0.0) Rdup_z = 0.0;
        // ** normalised downward diffuse flux at z
        double Rddn_z = (1.0 - trd) * (p3 * exp(-h * paiaa) + p4 * exp(h * paiaa))
            + trd;
        if (Rddn_z > 1.0) Rddn_z = 1.0;
        if (Rddn_z < 0.0) Rddn_z = 0.0;
        // Calculate time-variant two-stream parameters (direct)
        for (size_t i = 0; i < Rsw.size(); ++i) {
            if (Rsw[i] > 0) {
                // Calculate canopy extinction coefficient
                std::vector<double> kp = cankCpp(zen[i] * 180 / M_PI, x, si[i]);
                double k = kp[0];
                double kd = kp[1];
                double Kc = kp[2];
                // Calculate two-stream parameters (direct)      
                double sig = kd * kd + gma * gma - pow((a + gma), 2);
                std::vector<double> tspdir = twostreamdirCpp(pait, om, a, gma, J, del, h, gref, kd, sig, u1, S1, D1, D2);
                double p5 = tspdir[0];
                double p6 = tspdir[1];
                double p7 = tspdir[2];
                double p8 = tspdir[3];
                double p9 = tspdir[4];
                double p10 = tspdir[5];
                // Calculate gap tranmissions
                double trbn = pow(clump, Kc);
                double trb = pow(gi, Kc); // diffuse transmission donward
                if (trb > 0.999) trb = 0.999;
                if (trb < 0.0) trb = 0.0;
                if (trbn > 0.999) trbn = 0.999;
                if (trbn < 0.0) trbn = 0.0;
                // ** black-sky albedo
                double albb = (1.0 - trdn) * ((p5 / sig) + p6 + p7) + trdn * gref;
                if (albb > amx) albb = amx;
                if (albb < 0.01) albb = 0.01;
                // ** normalised contribution of direct to downward diffuse flux at ground
                double Rdbdn_g = (1.0 - trbn) * ((p8 / -sig) * exp(-kd * pait) +
                    p9 * exp(-h * pait) + p10 * exp(h * pait));
                if (Rdbdn_g > amx) Rdbdn_g = amx;
                if (Rdbdn_g < 0.0) Rdbdn_g = 0.0;
                // ** normalised contribution of direct to upward diffuse flux at z
                double Rdbup_z = (1.0 - trdu) * ((p5 / sig) * exp(-kd * pait) +
                    p6 * exp(-h * paiaa) + p7 * exp(h * paiaa))
                    + trdu * gref;
                if (Rdbup_z > amx) Rdbup_z = amx;
                if (Rdbup_z < 0.0) Rdbup_z = 0.0;
                // ** normalised contribution of direct to downward diffuse flux at z
                double Rdbdn_z = (1.0 - trdu) * ((p8 / -sig) * exp(-kd * paiaa) +
                    p9 * exp(-h * paiaa) + p10 * exp(h * paiaa));
                if (Rdbdn_z > amx) Rdbdn_z = amx;
                if (Rdbdn_z < 0.0) Rdbdn_z = 0.0;
                // Calculate incident flux
                double Rbeam = (Rsw[i] - Rdif[i]) / cos(zen[i]);
                if (Rbeam > 1352.0) Rbeam = 1352.0;
                double Rb = Rbeam * cos(zen[i]);
                double trg = trb + exp(-kd * pait); // tranmission to ground though gaps and leaves
                double Rbc = (trg * si[i] + (1 - trg) * cos(zen[i])) * Rbeam;
                // Calculate ground absorbed radiation
                double Rbdn_g = trbn + (1.0 - trbn) * exp(-kd * pait);
                if (Rbdn_g > 1.0) Rbdn_g = 1.0;
                if (Rbdn_g < 0.0) Rbdn_g = 0.0;
                radGsw[i] = (1.0 - gref) * (Rddn_g * Rdif[i] * svfa + Rdbdn_g * Rb +
                    Rbdn_g * Rbeam * si[i]);
                double maxg = radGsw[i] = (1 - gref) * (Rdif[i] * svfa + Rbeam * si[i]);
                if (radGsw[i] > maxg) radGsw[i] = maxg;
                // Calculate canopy absorbed radiation
                radCsw[i] = (1.0 - albd) * Rdif[i] * svfa + (1.0 - albb) * Rbc;
                // Calculate upward and downward fluxes
                Rbdown[i] = (trb + (1.0 - trb) * exp(-kd * paiaa)) * Rbeam;
                Rddown[i] = Rddn_z * Rdif[i] * svfa + Rdbdn_z * Rb;
                Rdup[i] = Rdup_z * Rdif[i] * svfa + Rdbup_z * Rb;
                // Leaf absorbed
                radLsw[i] = 0.5 * a * (Rddown[i] + Rdup[i] + k * cos(zen[i]) * Rbdown[i]);
            }
            else {
                radGsw[i] = 0;
                radCsw[i] = 0;
                Rbdown[i] = 0;
                Rddown[i] = 0;
                Rdup[i] = 0;
                radLsw[i] = 0;
            }
            // Calculate longwave radiation
            lwout[i] = 0.97 * 5.67 * pow(10, -8) * pow(tc[i] + 273.15, 4); // Longwave emitted
            radGlw[i] = 0.97 * (trd * svfa * lwdown[i] + (1 - trd) * (1 - svfa) * lwout[i]);
            radClw[i] = 0.97 * svfa * lwdown[i];
        } // end for i
    } // end pai > 0
    radmodel2 out;
    out.radGsw = radGsw;
    out.radGlw = radGlw;
    out.radCsw = radCsw;
    out.radClw = radClw;
    out.Rbdown = Rbdown;
    out.Rddown = Rddown;
    out.Rdup = Rdup;
    out.radLsw = radLsw;
    out.lwout = lwout;
    return out;
}
// Calculates a 3D array zenith angles and solar index values
// [[Rcpp::export]]
List solargrid(NumericMatrix slope, NumericMatrix aspect, DataFrame obstime, List micro)
{
    // Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access latlong from micro
    NumericMatrix lats = micro["lats"];
    NumericMatrix lons = micro["lons"];
    // Get dimensions
    int rows = slope.nrow();
    int cols = slope.ncol();
    int tsteps = year.size();
    // Initialise output variables
    NumericVector zen(rows * cols * tsteps);
    NumericVector si(rows * cols * tsteps);
    int index = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double val = slope(i, j);
            if (!NumericMatrix::is_na(val)) {
                for (int k = 0; k < tsteps; k++) {
                    std::vector<double> sp = solpositionCpp(lats(i, j), lons(i, j), year[k], 
                        month[k], day[k], hour[k]);
                    double zend = sp[0];
                    double azi = sp[1];
                    si[index] = solarindexCpp(slope(i, j), aspect(i, j), zend, azi);
                    zen[index] = zend * M_PI / 180;
                    index++;
                } // end k
            } // na check
            else {
                for (int k = 0; k < tsteps; k++) {
                    si[index] = NA_REAL;
                    zen[index] = NA_REAL;
                    index++;
                }
            }
        } // end j
    } // end i
    // Shape results
    zen.attr("dim") = NumericVector::create(tsteps, rows, cols);
    si.attr("dim") = NumericVector::create(tsteps, rows, cols);
    // Reshape results
    zen = aperm3D(zen, rows, cols, tsteps);
    si = aperm3D(si, rows, cols, tsteps);
    // Calculate vector of solar azimuths
    double lat = 0.0;
    double lon = 0.0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            lat = lat + lats(i, j);
            lon = lon + lons(i, j);
        }
    }
    lat = lat / (rows * cols);
    lon = lat / (rows * cols);
    NumericVector sazi(tsteps);
    for (int k = 0; k < tsteps; k++) {
        std::vector<double> sp = solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
        sazi[k] = sp[1] * M_PI / 180.0;
    }
    // Assign to micro
    micro["zen"] = Rcpp::wrap(zen);
    micro["si"] = Rcpp::wrap(si);
    micro["sazi"]= Rcpp::wrap(sazi);
    return micro;
}
// Grid version of two-stream radiation model
// [[Rcpp::export]]
List twostreamgrid(double reqhgt, List micro)
{
    // Access items of micro
    NumericVector hgt = micro["veghgt"];
    NumericVector pai = micro["pai"];
    NumericVector paia = micro["paia"];
    NumericVector x = micro["vegx"];
    NumericVector lref = micro["lref"];
    NumericVector ltra = micro["ltra"];
    NumericVector clump = micro["clump"];
    NumericVector gref = micro["gref"];
    NumericVector dirr = micro["dirr"];
    NumericVector difr = micro["difr"];
    NumericVector tc = micro["tc"];
    NumericVector lwdown = micro["lwdown"];
    NumericVector zen = micro["zen"];
    NumericVector si = micro["si"];
    NumericVector svfa = micro["svfa"];
    NumericVector lats = micro["lats"];
    NumericVector lons = micro["lons"];
    // Create output variables for storing data
    NumericVector radGsw(hgt.size());
    NumericVector radGlw(hgt.size());
    NumericVector radCsw(hgt.size());
    NumericVector radClw(hgt.size());
    NumericVector Rbdown(hgt.size());
    NumericVector Rddown(hgt.size());
    NumericVector Rdup(hgt.size());
    NumericVector radLsw(hgt.size());
    NumericVector lwout(hgt.size());
    for (R_xlen_t i = 0; i < hgt.size(); i++) {
        double val = hgt[i];
        if (!NumericVector::is_na(val)) {
            // ** Perform adjustments for gap fraction
            double pait = pai[i] / (1 - clump[i]);
            double rad = dirr[i] + difr[i];
            if (rad > 0.0) {
                double chck = lref[i];
                if (pai[i] < chck) chck = pai[i];
                if (chck > 0.0) { 
                    // ** Calculate base parameters
                    double om = lref[i] + ltra[i];
                    double a = 1 - om;
                    double del = lref[i] - ltra[i];
                    double J = 1.0 / 3.0;
                    if (x[i] != 1.0) {
                        double mla = 9.65 * pow(3 + x[i], -1.65);
                        if (mla > M_PI / 2) mla = M_PI / 2;
                        J = cos(mla) * cos(mla);
                    }
                    double gma = 0.5 * (om + J * del);
                    double h = sqrt(a * a + 2 * a * gma);
                    // Calculate two-stream parameters (diffuse)
                    std::vector<double> tspdif = twostreamdifCpp(pait, om, a, gma, h, gref[i]);
                    double p1 = tspdif[0];
                    double p2 = tspdif[1];
                    double p3 = tspdif[2];
                    double p4 = tspdif[3];
                    double u1 = tspdif[4];
                    double S1 = tspdif[5];
                    double D1 = tspdif[6];
                    double D2 = tspdif[7];
                    // Calculate canopy extinction coefficient
                    std::vector<double> kp = cankCpp(zen[i] * 180.0 / M_PI, x[i], si[i]);
                    double k = kp[0];
                    double kd = kp[1];
                    double Kc = kp[2];
                    // Calculate two-stream parameters (direct)      
                    double sig = kd * kd + gma * gma - pow((a + gma), 2);
                    std::vector<double> tspdir = twostreamdirCpp(pait, om, a, gma, J, del, h, gref[i], kd, sig, u1, S1, D1, D2);
                    double p5 = tspdir[0];
                    double p6 = tspdir[1];
                    double p7 = tspdir[2];
                    double p8 = tspdir[3];
                    double p9 = tspdir[4];
                    double p10 = tspdir[5];
                    // Calculate transmissions
                    double trdn = pow(clump[i], 2.0);
                    double trbn = pow(clump[i], Kc);
                    if (trbn > 0.999) trbn = 0.999;
                    if (trbn < 0.0) trbn = 0.0;
                    double gi = pow(clump[i], paia[i] / pai[i]);
                    double giu = pow(clump[i], (pai[i] - paia[i]) / pai[i]);
                    if (gi > 0.99) gi = 0.99;
                    if (giu > 0.99) giu = 0.99;
                    double trd = gi * gi; // diffuse transmission downard
                    double trb = pow(gi, Kc); // direct tranmission downward
                    if (trb > 0.999) trb = 0.999;
                    if (trb < 0.0) trb = 0.0;
                    double trdu = giu * giu; // diffuse transmission upward
                    double paiaa = paia[i] / (1.0 - gi); // paia adjusted for gap fraction
                    // white- and black-sky albedo
                    double amx = gref[i];
                    if (amx < lref[i]) amx = lref[i];
                    double albd = gref[i] * trdn + (1.0 - trdn) * (p1 + p2);
                    if (albd > amx) albd = amx;
                    if (albd < 0.01) albd = 0.01;
                    double albb = gref[i] * trdn + (1.0 - trdn) * ((p5 / sig) + p6 + p7);
                    if (albb > amx) albb = amx;
                    if (albb < 0.01) albb = 0.01;
                    // Normalised downward diffuse only at ground
                    double Rddn_g = (1.0 - trdn) * (p3 * exp(-h * pait) + p4 * exp(h * pait)) + trdn;
                    if (Rddn_g > 1.0) Rddn_g = 1.0;
                    if (Rddn_g < 0.0) Rddn_g = 0.0;
                    // Normalised upward diffuse only at z
                    double Rddn_z = (1.0 - trdu) * (p1 * exp(-h * paiaa) + p2 * exp(h * paiaa)) + trdu * gref[i];
                    if (Rddn_z > 1.0) Rddn_z = 1.0;
                    if (Rddn_z < 0.0) Rddn_z = 0.0;
                    // Normalised donward diffuse only at z
                    double Rdup_z = (1.0 - trd) * (p3 * exp(-h * paiaa) + p4 * exp(h * paiaa)) + trd;
                    if (Rdup_z > 1.0) Rdup_z = 1.0;
                    if (Rdup_z < 0.0) Rdup_z = 0.0;
                    // Normalised contribution of direct to downward diffuse at ground
                    double Rdbdn_g = (1.0 - trbn) * ((p8 / -sig) * exp(-kd * pait) +
                        p9 * exp(-h * pait) + p10 * exp(h * pait));
                    if (Rdbdn_g > amx) Rdbdn_g = amx;
                    if (Rdbdn_g < 0.0) Rdbdn_g = 0.0;
                    // Normalised contribution of direct to upward diffuse at z
                    double Rdbdn_z = (1.0 - trdu) * ((p5 / sig) * exp(-kd * paiaa) +
                        p6 * exp(-h * paiaa) + p7 * exp(h * paiaa)) + trdu * gref[i];
                    if (Rdbdn_z > amx) Rdbdn_z = amx;
                    if (Rdbdn_z < 0.0) Rdbdn_z = 0.0;
                    // Normalised contribution of direct to downward diffuse at z
                    double Rdbup_z = (1.0 - trd) * ((p8 / -sig) * exp(-kd * paiaa) + 
                        p9 * exp(-h * paiaa) + p10 * exp(h * paiaa)) + trd;
                    if (Rdbup_z > amx) Rdbup_z = amx;
                    if (Rdbup_z < 0.0) Rdbup_z = 0.0;
                    // Calculate incident flux
                    double Rb = dirr[i] * cos(zen[i]);
                    double trg = trb + exp(-kd * pait); // tranmission to ground though gaps and leaves
                    double Rbc = (trg * si[i] + (1 - trg) * cos(zen[i])) * dirr[i];
                    // Calculate ground absorbed radiation
                    double Rdirg = (trbn + (1.0 - trbn) * exp(-kd * pait)) * dirr[i] * si[i];
                    radGsw[i] = (Rdirg + Rddn_g * difr[i] * svfa[i] +
                        Rdbdn_g * Rb) * (1.0 - gref[i]);
                    double maxg = radGsw[i] = (1 - gref[i]) * (difr[i] * svfa[i] + dirr[i] * si[i]); 
                    if (radGsw[i] > maxg) radGsw[i] = maxg;
                    // Calculate canopy and ground absorbed radiation
                    radCsw[i] = (1.0 - albd) * difr[i] * svfa[i] +
                        (1.0 - albb) * Rbc;
                    // Calculate fluxes
                    Rbdown[i] = ((1 - trb) * exp(-kd * paiaa) + trb) * dirr[i];
                    Rddown[i] = Rddn_z * difr[i] * svfa[i] + Rdbdn_z * Rb;
                    Rdup[i] = Rdup_z * difr[i] * svfa[i] + Rdbup_z * Rb;
                    // Calculate leaf absorbed
                    radLsw[i] = 0.5 * a * (Rddown[i] + Rdup[i] + k * cos(zen[i]) * Rbdown[i]);
                } // end pait > 0
                else {
                    radGsw[i] = (1 - gref[i]) * (difr[i] * svfa[i] + dirr[i] * si[i]); // Ground absorbed
                    radCsw[i] = radGsw[i];
                    Rbdown[i] = dirr[i] * cos(zen[i]);
                    Rddown[i] = difr[i] * svfa[i];
                    Rdup[i] = gref[i] * (difr[i] * svfa[i] + dirr[i] * cos(zen[i]));
                }
            } // end rad > 0
            // Calculate longwave radiation
            double trd = (1 - pow(clump[i], 2)) * exp(-pait) + pow(clump[i], 2);
            lwout[i] = 0.97 * 5.67 * pow(10, -8) * pow(tc[i] + 273.15, 4); // Longwave emitted
            radGlw[i] = 0.97 * (trd * svfa[i] * lwdown[i] + (1 - trd) * (1 - svfa[i]) * lwout[i]);
            radClw[i] = svfa[i] * lwdown[i];
        } // end NA check
        else {
            radGsw[i] = NA_REAL;
            radGlw[i] = NA_REAL;
            radCsw[i] = NA_REAL;
            radClw[i] = NA_REAL;
            Rbdown[i] = NA_REAL;
            Rddown[i] = NA_REAL;
            Rdup[i] = NA_REAL;
            radLsw[i] = NA_REAL;
            lwout[i] = NA_REAL;
        } // end NA check
    } // end i
    // Get dimensions of arrays
    IntegerVector dims = hgt.attr("dim");
    int n1 = dims[0];
    int n2 = dims[1];
    int n3 = dims[2];
    // reshape results
    radGsw.attr("dim") = NumericVector::create(n1, n2, n3);
    radGlw.attr("dim") = NumericVector::create(n1, n2, n3);
    radCsw.attr("dim") = NumericVector::create(n1, n2, n3);
    radClw.attr("dim") = NumericVector::create(n1, n2, n3);
    Rbdown.attr("dim") = NumericVector::create(n1, n2, n3);
    Rddown.attr("dim") = NumericVector::create(n1, n2, n3);
    Rdup.attr("dim") = NumericVector::create(n1, n2, n3);
    radLsw.attr("dim") = NumericVector::create(n1, n2, n3);
    lwout.attr("dim") = NumericVector::create(n1, n2, n3);
    // Create output
    micro["radGsw"] = Rcpp::wrap(radGsw);
    micro["radGlw"] = Rcpp::wrap(radGlw);
    micro["radCsw"] = Rcpp::wrap(radCsw);
    micro["radClw"] = Rcpp::wrap(radClw);
    micro["Rbdown"] = Rcpp::wrap(Rbdown);
    micro["Rddown"] = Rcpp::wrap(Rddown);
    micro["Rdup"] = Rcpp::wrap(Rdup);
    micro["radLsw"] = Rcpp::wrap(radLsw);
    micro["lwout"] = Rcpp::wrap(lwout);
    return micro;
}
// Point version of two-stream radiation model with snow
NumericVector twostreampoint(double reqhgt, double hgt, double pai, double paia, double lref, double ltra, double clump,
    double albg, double albc, double Rsw, double Rdif, double shadowmask, double tc, double lwdown, double zenr, 
    double si, double svfa)
{
    // ** Perform adjustments for gap fraction
    double pait = pai / (1 - clump);
    double Rbdown = 0.0;
    double Rddown = 0.0;
    double Rdup = 0.0;
    double radLsw = 0.0;
    if (Rsw > 0.0) {
        double Rbeam = (Rsw - Rdif) / cos(zenr);
        if (Rbeam > 1352.0) Rbeam = 1352.0;
        if (pai > 0.0) {
            // ** Calculate base parameters
            double om = lref + ltra;
            double a = 1 - om;
            double del = lref - ltra;
            double J = 1.0 / 3.0;
            double gma = 0.5 * (om + J * del);
            double h = sqrt(a * a + 2 * a * gma);
            // Calculate two-stream parameters (diffuse)
            std::vector<double> tspdif = twostreamdifCpp(pait, om, a, gma, h, albg);
            double p1 = tspdif[0];
            double p2 = tspdif[1];
            double p3 = tspdif[2];
            double p4 = tspdif[3];
            double u1 = tspdif[4];
            double S1 = tspdif[5];
            double D1 = tspdif[6];
            double D2 = tspdif[7];
            // Calculate canopy extinction coefficient
            std::vector<double> kp = cankCpp(zenr * 180.0 / M_PI, 1.0, si);
            double kd = kp[1];
            double Kc = kp[2];
            // Calculate two-stream parameters (direct)      
            double sig = kd * kd + gma * gma - pow((a + gma), 2);
            std::vector<double> tspdir = twostreamdirCpp(pait, om, a, gma, J, del, h, albg, kd, sig, u1, S1, D1, D2);
            double p5 = tspdir[0];
            double p6 = tspdir[1];
            double p7 = tspdir[2];
            double p8 = tspdir[3];
            double p9 = tspdir[4];
            double p10 = tspdir[5];
            // Calculate tranmissions
            // ~~ Gap fraction
            double gi = pow(clump, paia / pai);
            double giu = pow(clump, (pai - paia) / pai);
            if (gi > 0.99) gi = 0.99;
            if (giu > 0.99) giu = 0.99;
            double paiaa = paia / (1.0 - gi); // paia adjusted for gap fraction
            // ~~ transmissions
            double trd = gi * gi;
            double trdu = giu * giu;
            double trb = pow(gi, Kc);
            if (trb > 0.999) trb = 0.999;
            if (trb < 0.0) trb = 0.0;
            // Normalised upward diffuse only at z
            double Rddn_z = (1.0 - trdu) * (p1 * exp(-h * paiaa) + p2 * exp(h * paiaa)) + trdu * albg;
            if (Rddn_z > 1.0) Rddn_z = 1.0;
            if (Rddn_z < 0.0) Rddn_z = 0.0;
            // Normalised downward diffuse only at z
            double Rdup_z = (1.0 - trd) * (p3 * exp(-h * paiaa) + p4 * exp(h * paiaa)) + trd;
            if (Rdup_z > 1.0) Rdup_z = 1.0;
            if (Rdup_z < 0.0) Rdup_z = 0.0;
            // Normalised contribution of direct to upward diffuse at z
            double Rdbdn_z = (1.0 - trdu) * ((p5 / sig) * exp(-kd * paiaa) +
                p6 * exp(-h * paiaa) + p7 * exp(h * paiaa)) + trdu * albg;
            if (Rdbdn_z > 1.0) Rdbdn_z = 1.0;
            if (Rdbdn_z < 0.0) Rdbdn_z = 0.0;
            // Normalised contribution of direct to downward diffuse at z
            double Rdbup_z = (1.0 - trd) * ((p8 / -sig) * exp(-kd * paiaa) +
                p9 * exp(-h * paiaa) + p10 * exp(h * paiaa)) + trd;
            if (Rdbup_z > 1.0) Rdbup_z = 1.0;
            if (Rdbup_z < 0.0) Rdbup_z = 0.0;
            // Calculate incident flux
            double Rbeam = (Rsw - Rdif) / cos(zenr);
            double Rb = Rsw - Rdif;
            Rbdown = (trb + (1.0 - trb) * exp(-kd * paiaa)) * Rbeam * shadowmask;
            Rddown = Rddn_z * Rdif * svfa + Rdbdn_z * Rb * shadowmask;
            Rdup = Rdup_z * Rdif * svfa + Rdbup_z * Rb * shadowmask;
            radLsw = 0.5 * a * (Rddown + Rdup + cos(zenr) * 0.5 * Rbdown);
        } // end pait > 0
        else {
            Rbdown = Rbeam * shadowmask;
            Rddown = Rdif * svfa;
            Rdup = albg * (Rdif * svfa + Rbeam * shadowmask * cos(zenr));
        }
    } // end rad > 0
    // Create output
    NumericVector out(4);
    out[0] = Rbdown;
    out[1] = Rddown;
    out[2] = Rdup;
    out[3] = radLsw;
    return out;
}
// ** Calculate wind speed vector ** //
windmodel windvCpp(double reqhgt, double zref, std::vector<double> u2, std::vector<double> umu,
    double h, double pai, std::vector<double> ws)
{
    // Calculate roughness length and zero-plane displacement
    double d = zeroplanedisCpp(h, pai);
    double zm = roughlengthCpp(h, pai, d, 0);
    // Create variables
    std::vector<double> uf(u2.size());
    std::vector<double> uz(u2.size());
    std::vector<double> gHa(u2.size());
    for (size_t i = 0; i < u2.size(); ++i) {
        double ufs = (0.4 * u2[i]) / log((zref - d) / zm);
        uf[i] = ufs * umu[i] * ws[i];
        uz[i] = uf[i];
        if (reqhgt > 0) {
            if (reqhgt >= h) {
                uz[i] = (uf[i] / 0.4) * log((reqhgt - d) / zm);
            }
            else {
                // Calculate wind speed at height z below canopy
                double uh = (uf[i] / 0.4) * log((h - d) / zm);
                if (uh < uf[i]) uh = uf[i];
                double Be = uf[i] / uh;
                if (Be < 0.001) Be = 0.001;
                double a = pai / h;
                double Lc = pow(0.25 * a, -1);
                double Lm = 2 * pow(Be, 3) * Lc;
                uz[i] = uh * exp(Be * (reqhgt - h) / Lm);
            }
            if (uz[i] > u2[i]) uz[i] = u2[i];
        }
        // Calculate gHa
        gHa[i] = gturbCpp(uf[i], d, zm, zref, 43, 0, 0.0001);
    }
    windmodel out;
    out.uf = uf;
    out.uz = uz;
    out.gHa = gHa;
    return(out);
}
// Grid version of wind downscale model
// [[Rcpp::export]]
List windgrid(double reqhgt, List micro)
{
    // Access items of micro
    NumericVector u2 = micro["u2"]; // wind speed at zref
    NumericVector h = micro["veghgt"]; // vegetation height
    NumericVector pai = micro["pai"]; // vegetation plant area index
    NumericVector ws = micro["ws"]; // wind shelter coefficient
    NumericVector umu = micro["umu"]; // Stability multiplier
    double zref = micro["zref"]; // reference height
    // Create output variables
    NumericVector uf(u2.size());
    NumericVector uz(u2.size());
    NumericVector gHa(u2.size());
    for (R_xlen_t i = 0; i < h.size(); i++) {
        // Calculate point model variables
        double d = zeroplanedisCpp(h[i], pai[i]);
        double zm = roughlengthCpp(h[i], pai[i], d, 0);
        double ufs = (0.4 * u2[i]) / log((zref - d) / zm);
        // Calculate wind friction velocity
        uf[i] = ufs * umu[i] * ws[i];
        uz[i] = uf[i];
        if (reqhgt > 0) {
            // Calculate wind speed at height z above canopy
            if (reqhgt >= h[i]) {
                uz[i] = (uf[i] / 0.4) * log((reqhgt - d) / zm);
            }
            // Calculate wind speed at height z below canopy
            else {
                // Calculate wind speed at height z below canopy
                double uh = (uf[i] / 0.4) * log((h[i] - d) / zm);
                if (uh < uf[i]) uh = uf[i];
                double Be = uf[i] / uh;
                if (Be < 0.001) Be = 0.001;
                double a = pai[i] / h[i];
                double Lc = pow(0.25 * a, -1.0);
                double Lm = 2 * pow(Be, 3.0) * Lc;
                uz[i] = uh * exp(Be * (reqhgt - h[i]) / Lm);
            }
            if (uz[i] > u2[i]) uz[i] = u2[i];
        }
        // Calculate gHa
        gHa[i] = gturbCpp(uf[i], d, zm, zref, 43, 0, 0.0001);
    }
    // Get dimensions of arrays
    IntegerVector dims = h.attr("dim");
    int n1 = dims[0];
    int n2 = dims[1];
    int n3 = dims[2];
    // reshape results
    uf.attr("dim") = NumericVector::create(n1, n2, n3);
    uz.attr("dim") = NumericVector::create(n1, n2, n3);
    gHa.attr("dim") = NumericVector::create(n1, n2, n3);
    // Create output
    micro["uf"] = Rcpp::wrap(uf);
    micro["uz"] = Rcpp::wrap(uz);
    micro["gHa"] = Rcpp::wrap(gHa);
    return micro;
}
NumericVector windpoint(double reqhgt, double u2, double h, double pai, double ws, double umu, double zref)
{
    // Calculate point model variables
    double d = zeroplanedisCpp(h, pai);
    double zm = roughlengthCpp(h, pai, d, 0);
    double ufs = (0.4 * u2) / log((zref - d) / zm);
    // Calculate wind friction velocity
    double uf = ufs * umu * ws;
    double uz = uf;
    // Calculate wind speed at height z above canopy
    if (reqhgt >= h) {
        uz = (uf / 0.4) * log((reqhgt - d) / zm);
    }
    // Calculate wind speed at height z below canopy
    else {
        // Calculate wind speed at height z below canopy
        double uh = (uf / 0.4) * log((h - d) / zm);
        if (uh < uf) uh = uf;
        double Be = uf / uh;
        if (Be < 0.001) Be = 0.001;
        double a = pai / h;
        double Lc = pow(0.25 * a, -1.0);
        double Lm = 2 * pow(Be, 3.0) * Lc;
        uz = uh * exp(Be * (reqhgt - h) / Lm);
    }
    if (uz > u2) uz = u2;
    double gHa = gturbCpp(uf, d, zm, zref, 43, 0, 0.0001);
    NumericVector out(3);
    out[0] = uz;
    out[1] = uf;
    out[2] = gHa;
    return out;
}
// Simple PenmanMonteith function
// [[Rcpp::export]]
std::vector<double> PenmanMonteith2Cpp(double Rabs, double gHa, double gV, double tc, double mxtc, double pk, double ea, double es,
    double G, double surfwet, double tdew)
{
    double De = satvapCpp(tc + 0.5) - satvapCpp(tc - 0.5);
    double sb = 5.67 * pow(10, -8);
    double gHr = gHa + (4 * 0.97 * sb * pow(tc + 273.15, 3)) / 29.3;
    double Rem = 0.97 * sb * pow(tc + 273.15, 4);
    double la = 0;
    if (tc >= 0) {
        la = 45068.7 - 42.8428 * tc;
    }
    else {
        la = 51078.69 - 4.338 * tc - 0.06367 * tc * tc;
    }
    double m = la * (gV / pk);
    double L = m * (es - ea) * surfwet;
    double dT = (Rabs - Rem - L - G) / (29.3 * gHr + m * De);
    double dTmx = -0.6273 * mxtc + 49.79;
    if (dT > dTmx) dT = dTmx;
    if (dT > 80.0) dT = 80.0;
    double T0 = dT + tc;
    if (T0 < tdew) T0 = tdew;
    std::vector<double> out(5);
    out[0] = T0;
    out[1] = 29.3 * gHa * (T0 - tc);
    out[2] = m * (satvapCpp(T0) - ea) * surfwet;
    out[3] = 0.97 * sb * pow(T0 + 273.15, 4);
    out[4] = la * (43 / pk);
    return out;
}
// Calculate soil conductance
NumericMatrix soilcondCpp(double rho, double Vm, double Vq, double Mc, std::vector<double> soilm)
{
    // Create output variable
    int n = soilm.size();
    NumericMatrix out(n, 3); // 0 = k, 1 = DD
    // Time invariant variables
    double frs = Vm + Vq;
    double c1 = (0.57 + 1.73 * Vq + 0.93 * Vm) / (1 - 0.74 * Vq - 0.49 * Vm) - 2.8 * frs * (1 - frs);
    double c3 = 1 + 2.6 * pow(Mc, -0.5);
    double c4 = 0.03 + 0.7 * pow(frs, 2);
    for (size_t i = 0; i < soilm.size(); ++i) {
        // Find soil diffusivity
        double cs = (2400 * rho / 2.64 + 4180 * soilm[i]); // specific heat of soil in J / kg / K
        double ph = (rho * (1 - soilm[i]) + soilm[i]) * 1000; // bulk density in kg / m3
        double c2 = 1.06 * rho * soilm[i];
        out(i, 0) = c1 + c2 * soilm[i] - (c1 - c4) * exp(-pow(c3 * soilm[i], 4)); // Thermal conductivity W / m / K
        double ka = out(i, 0) / (cs * ph);
        double omdy = (2 * M_PI) / (24 * 3600);
        out(i, 1) = pow(2 * ka / omdy, 0.5); // Damping depth
    }
    return out;
}
// ** Calculate ground surface temperature using hourly data  ** //
soilmodel soiltemp_hrCpp(std::vector<double> tc, std::vector<double> ea, std::vector<double> es, 
    std::vector<double> tdew, std::vector<double> pk, std::vector<double> radGsw, 
    std::vector<double> radGlw, std::vector<double> gHa,
    std::vector<double> soilparamsm, std::vector<double> Gp,
    std::vector<double> kp, std::vector<double> dtrp, std::vector<double> muGp,
    std::vector<double> soilm, double mxtc)
{
    // Create variables
    std::vector<double> T0(tc.size());
    std::vector<double> Rnet(tc.size());
    std::vector<double> surfwet(tc.size());
    std::vector<double> radabs(tc.size());
    for (size_t i = 0; i < tc.size(); ++i) {
        // Radiation absorbed by ground
        radabs[i] = radGsw[i] + radGlw[i];
        // Calculate soil surface effective relative humidity
        double matric = -soilparamsm[3] * pow(soilm[i] / soilparamsm[0], -soilparamsm[2]);
        surfwet[i] = exp((0.018 * matric) / (8.31 * (tc[i] + 273.15)));
        if (surfwet[i] > 1) surfwet[i] = 1;
        // Calculate soil surface temperature assuming G to be zero
        std::vector<double> T0v = PenmanMonteith2Cpp(radabs[i], gHa[i], gHa[i], tc[i], mxtc, pk[i], ea[i], es[i], 0.0, surfwet[i], tdew[i]);
        T0[i] = T0v[0];
        Rnet[i] = radabs[i] - T0v[3];
    }
    // Calculate diurnal range
    std::vector<double> T0mx = hourtodayCpp(T0, "max");
    std::vector<double> T0mn = hourtodayCpp(T0, "min");
    std::vector<double> dtR(T0mx.size());
    for (size_t i = 0; i < dtR.size(); ++i) {
        double dtr = T0mx[i] - T0mn[i];
        dtR[i] = dtr / dtrp[i];
    }
    // Calculate thermal conductance and damping depth
    NumericMatrix kDDg = soilcondCpp(soilparamsm[7], soilparamsm[5], soilparamsm[4], soilparamsm[6], soilm);
    // Calculate diurnal radiation cycle for setting limits to G
    std::vector<double> Rdmx = hourtodayCpp(Rnet, "max");
    std::vector<double> Rdmn = hourtodayCpp(Rnet, "min");
    std::vector<double> G(tc.size());
    for (size_t i = 0; i < G.size(); ++i) {
        // Calculate Gmu
        double Gmu = dtR[i] * (kDDg(i, 0) * muGp[i]) / (kp[i] * kDDg(i, 1));
        // Compute ground heat flux
        G[i] = Gp[i] * Gmu;
        double Rd = Rdmx[i];
        if (-Rdmn[i] > Rd) Rd = -Rdmn[i];
        if (G[i] > 0.6 * Rd) G[i] = 0.6 * Rd;
        if (G[i] < -0.6 * Rd) G[i] = -0.6 * Rd;
        // Re-calculate soil surface temperature
        std::vector<double> T0v = PenmanMonteith2Cpp(radabs[i], gHa[i], gHa[i], tc[i], mxtc, pk[i], ea[i], es[i], G[i], surfwet[i], tdew[i]);
        T0[i] = T0v[0];
    }
    soilmodel out;
    out.Tg = T0;
    out.G = G;
    for (size_t i = 0; i < G.size(); ++i) {
        out.kDDg.push_back(kDDg(i, 1)); // kDDg is the damping depth
    }
    return out;
}
// Calculate gorund surface temperature: grid
// [[Rcpp::export]]
List soiltempgrid(List micro)
{
    // Access items of micro
    // weather
    NumericVector tc = micro["tc"];
    NumericVector ea = micro["ea"];
    NumericVector es = micro["estl"];
    NumericVector tdew = micro["tdew"];
    NumericVector pk = micro["pk"];
    NumericVector radGsw = micro["radGsw"];
    NumericVector radGlw = micro["radGlw"];
    NumericVector gHa = micro["gHa"];
    // soil variables
    NumericMatrix rho = micro["rho"];
    NumericMatrix Vm = micro["Vm"];
    NumericMatrix Vq = micro["Vq"];
    NumericMatrix Mc = micro["Mc"];
    NumericMatrix soilb = micro["soilb"];
    NumericMatrix psi_e = micro["psi_e"];
    NumericMatrix Smax = micro["Smax"];
    NumericMatrix Smin = micro["Smin"];
    // point model soil
    NumericVector Gp = micro["Gp"]; 
    NumericVector kp = micro["kp"];
    NumericVector muGp = micro["muGp"]; 
    NumericVector dtrp = micro["dtrp"];
    // Soil misture
    NumericVector soilm = micro["soilm"];
    // Get dimensions of arrays
    IntegerVector dims = tc.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    // repermutate variables
    // weather
    tc = aperm3D2(tc, rows, cols, tsteps);
    ea = aperm3D2(ea, rows, cols, tsteps);
    es = aperm3D2(es, rows, cols, tsteps);
    tdew = aperm3D2(tdew, rows, cols, tsteps);
    pk = aperm3D2(pk, rows, cols, tsteps);
    radGsw = aperm3D2(radGsw, rows, cols, tsteps);
    radGlw = aperm3D2(radGlw, rows, cols, tsteps);
    gHa = aperm3D2(gHa, rows, cols, tsteps);
    // point model and soilm
    Gp = aperm3D2(Gp, rows, cols, tsteps);
    kp = aperm3D2(kp, rows, cols, tsteps);
    muGp = aperm3D2(muGp, rows, cols, tsteps);
    dtrp = aperm3D2(dtrp, rows, cols, tsteps);
    soilm = aperm3D2(soilm, rows, cols, tsteps);
    // Create output variables
    NumericVector Tg(tc.size());
    NumericVector G(tc.size());
    NumericVector kDDg(tc.size());
    // set indices for storing results
    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = rho(i, j);
            if (!NumericMatrix::is_na(val)) {
                // subset climate variables
                int st = j * tsteps + i * tsteps * cols;
                int ed = st + tsteps - 1;
                // subset to vector
                std::vector<double> tcv = tosvd(tc[Range(st, ed)]);
                std::vector<double> eav = tosvd(ea[Range(st, ed)]);
                std::vector<double> esv = tosvd(es[Range(st, ed)]);
                std::vector<double> tdewv = tosvd(tdew[Range(st, ed)]);
                std::vector<double> pkv = tosvd(pk[Range(st, ed)]);
                std::vector<double> radGswv = tosvd(radGsw[Range(st, ed)]);
                std::vector<double> radGlwv = tosvd(radGlw[Range(st, ed)]);
                std::vector<double> gHav = tosvd(gHa[Range(st, ed)]);
                // soil point model
                std::vector<double> Gpv = tosvd(Gp[Range(st, ed)]);
                std::vector<double> kpv = tosvd(kp[Range(st, ed)]);
                std::vector<double> muGpv = tosvd(muGp[Range(st, ed)]);
                std::vector<double> dtrpv = tosvd(dtrp[Range(st, ed)]);
                std::vector<double> soilmv = tosvd(soilm[Range(st, ed)]);
                // create soilparamsm
                std::vector<double> soilparamsm = { Smax(i, j), Smin(i, j), soilb(i, j),
                    psi_e(i, j), Vq(i, j), Vm(i, j), Mc(i, j), rho(i, j) };
                // maxtc calc
                double mxtc = -273.15;
                for (int k = 0; k < tsteps; ++k) {
                    if (tcv[k] > mxtc) mxtc = tcv[k];
                }
                soilmodel sout = soiltemp_hrCpp(tcv, eav, esv, tdewv, pkv, radGswv, radGlwv,
                    gHav, soilparamsm, Gpv, kpv, dtrpv, muGpv, soilmv, mxtc);
                for (int k = 0; k < tsteps; ++k) {
                    Tg[index] = sout.Tg[k];
                    G[index] = sout.G[k];
                    kDDg[index] = sout.kDDg[k];
                    index++;
                }
            }
            else {
                for (int k = 0; k < tsteps; ++k) {
                    Tg[index] = NA_REAL;
                    G[index] = NA_REAL;
                    kDDg[index] = NA_REAL;
                    index++;
                }
            } // end NA test
        } // end j
    } // end i
    // Shape results
    Tg.attr("dim") = NumericVector::create(tsteps, rows, cols);
    G.attr("dim") = NumericVector::create(tsteps, rows, cols);
    kDDg.attr("dim") = NumericVector::create(tsteps, rows, cols);
    // Reshape results
    Tg = aperm3D(Tg, rows, cols, tsteps);
    G = aperm3D(G, rows, cols, tsteps);
    kDDg = aperm3D(kDDg, rows, cols, tsteps);
    // Assign to micro
    // Create output
    micro["Tg"] = Rcpp::wrap(Tg);
    micro["G"] = Rcpp::wrap(G);
    micro["kDDg"] = Rcpp::wrap(kDDg);
    return micro;
}
// Calculate temperature and vapour pressure above canopy
std::vector<double> TVabove(double reqhgt, double zref, double h, double pai, double T0, double tc, double ea, double surfwet = 1)
{
    if (h < 0.001) h = 0.001;
    double d = zeroplanedisCpp(h, pai);
    double zm = roughlengthCpp(h, pai, d, 0);
    double zh = 0.2 * zm;
    double lnr = log((reqhgt - d) / zh) / log((zref - d) / zh);
    std::vector<double> out(2);
    // Temperature
    out[0] = tc + (T0 - tc) * (1 - lnr);
    double estl = satvapCpp(T0);
    out[1] = ea + (estl - ea) * surfwet * (1 - lnr);
    //out[1] = lnr;
    return out;
}
std::vector<double> leaftemp(double Tcan, double T0, double tc, double mxtc, double pk, double ea, double es, double uz, double tdew,
    double surfwet, double zen, double radLsw, double Rddown, double Rbdown, double Rlw,
    double pai, double paia, double leafd, double gsmax)
{
    // radiation absorbed by leaf
    double lwcan = 0.97 * 5.67 * pow(10, -8) * pow(Tcan + 273.15, 4);
    double lwgro = 0.97 * 5.67 * pow(10, -8) * pow(T0 + 273.15, 4);
    double paig = pai - paia;
    double lwup = exp(-paig) * lwgro + (1 - exp(-paig)) * lwcan;
    double lwdn = exp(-paia) * Rlw + (1 - exp(-paia)) * lwcan;
    double lwabs = 0.97 * 0.5 * (lwup + lwdn);
    double leafabs = radLsw + lwabs;
    // ** Conductances
    double gh = 0.135 * sqrt(uz / leafd) * 1.4;
    double rsw = Rddown + cos(zen) * Rbdown;
    double gV = gh;
    if (gsmax < 999.99) {
        double gv = stomcondCpp(rsw, gsmax, 100);
        gV = 1 / (1 / gh + 1 / gv);
    }
    // Temperature
    std::vector<double> tl = PenmanMonteith2Cpp(leafabs, gh, gV, tc, mxtc, pk, ea, es, 0.0, surfwet, tdew);
    std::vector<double> out(5);
    out[0] = tl[0];
    out[1] = tl[1];
    out[2] = tl[2];
    out[3] = lwdn;
    out[4] = lwup;
    return out;
}
// [[Rcpp::export]]
double rhcanopy(double uf, double h, double d, double z)
{
    double a2 = 0.4 * (1.0 - (d / h)) / std::pow(1.25, 2);
    double inth = 4.293251 * h;
    if (z != h) {
        inth = (2.0 * h * ((48 * atan((sqrt(5.0) * sin((M_PI * z) / h)) / (cos((M_PI * z) / h) + 1))) /
            pow(5.0, 1.5) + (32 * sin((M_PI * z) / h)) / ((cos((M_PI * z) / h) + 1) * ((25 * pow(sin((M_PI * z) / h), 2)) /
                pow((cos((M_PI * z) / h) + 1), 2) + 5)))) / M_PI;
    }
    double mu = uf / (a2 * h) * 1.0 / (uf * uf);
    double rHa = inth * mu;
    if (rHa < 0.001) rHa = 0.001;
    return rHa;
}
double TVbelow(double zref, double z, double d, double h, double pai, double uf,
    double leafden, double Flux, double Fluxz, double SH, double SG)
{
    // Calculate thermal diffusivities
    double Rc = rhcanopy(uf, h, d, h);
    double Kc = h / Rc;
    double Kg = 1.0 / rhcanopy(uf, h, d, z);
    double Kh = 1.0 / (Rc - rhcanopy(uf, h, d, z));
    Kg = Kg / z;
    Kh = Kh / (h - z);
    // Calculate canopy source concentration
    double SC = SH + Flux / Kc;
    // Calculate far-field source concentration
    double farg = (Kg * SG + Kh * SH + Kc * SC) / (Kg + Kh + Kc);
    // Calculate leaf temperature and flux
    double SN = Fluxz * leafden;
    // Calculate near-field
    double near = (3.047519 + 0.128642 * log(pai)) * SN;
    if (std::isnan(near)) near = 0;
    double SCFN = near + farg;
    return SCFN;
}
// Calculate temperature or vapour pressure above ground (zen in radiations)
abovemodel TVabovegroundv(double reqhgt, double zref, std::vector<double> tc, std::vector<double> pk, std::vector<double> ea,
    std::vector<double> es, std::vector<double> tdew, std::vector<double> Rsw, std::vector<double> lwdn,
    std::vector<double> T0, std::vector<double> uz, std::vector<double> uf, std::vector<double> gHa, std::vector<double> zen,
    std::vector<double> radCsw, std::vector<double> radClw, std::vector<double> radLsw, std::vector<double> Rddown, std::vector<double> Rbdown,
    std::vector<double> soild, std::vector<double> G,
    double hgt, double pai, double paia, double leafd, double leafden,
    double Smin, double Smax, double gsmax, double mxtc)
{
    // Create variable
    std::vector<double> Tz(tc.size());
    std::vector<double> tleaf(tc.size());
    std::vector<double> rh(tc.size());
    std::vector<double> lwup(tc.size());
    for (size_t i = 0; i < tc.size(); ++i) {
        // Calculate ground and surface wetness
        double eT = satvapCpp(T0[i]) - ea[i];
        if (eT < 0.001) eT = 0.001;
        double plf = 0.8753 - 1.7126 * log(eT);
        double gwet = 1 / (1 + exp(-plf));
        double surfwet = (soild[i] - Smin) / (Smax - Smin);
        if (surfwet > gwet) gwet = surfwet;
        // Calculate vapour conductivity
        double gv = stomcondCpp(Rsw[i], gsmax, 300) * 3;
        double gV = 1 / (1 / gHa[i] + 1 / gv);
        // Calculate canopy temperature
        double Rabs = radCsw[i] + radClw[i];
        std::vector<double> TH = PenmanMonteith2Cpp(Rabs, gHa[i], gV, tc[i], mxtc, pk[i], ea[i], es[i], G[i], surfwet, tdew[i]);
        double Tcan = TH[0];
        double ez = 0;
        if (reqhgt >= hgt) {
            std::vector<double> tv = TVabove(reqhgt, zref, hgt, pai, Tcan, tc[i], ea[i], surfwet);
            Tz[i] = tv[0];
            tleaf[i] = Tcan;
            ez = tv[1];
            lwup[i] = 0.97 * 5.67 * pow(10, -8) * pow(Tcan + 273.15, 4);
        }
        else {
            // Calculate Leaf temperature
            std::vector<double> tvl = leaftemp(Tcan, T0[i], tc[i], mxtc, pk[i], ea[i], es[i], uz[i], tdew[i], surfwet, zen[i],
                radLsw[i], Rddown[i], Rbdown[i], lwdn[i], pai, paia, leafd, gsmax);
            tleaf[i] = tvl[0];
            // Calculate temperature below canopy
            double Flux = TH[1] * (1 - exp(-pai));
            double Fluxz = tvl[1];
            std::vector<double> tv = TVabove(hgt, zref, hgt, pai, Tcan, tc[i], ea[i], surfwet);
            double SH = tv[0] * 29.3 * 43;
            double SG = T0[i] * 29.3 * 43;
            double d = zeroplanedisCpp(hgt, pai);
            Tz[i] = TVbelow(zref, reqhgt, d, hgt, pai, uf[i], leafden, Flux, Fluxz, SH, SG) / (29.3 * 43);
            // Calculate humidity below canopy
            Flux = TH[2] * (1 - exp(-pai));
            Fluxz = tvl[2];
            SH = tv[1] * TH[4];
            SG = satvapCpp(T0[i]) * gwet * TH[4];
            ez = TVbelow(zref, reqhgt, d, hgt, pai, uf[i], leafden, Flux, Fluxz, SH, SG) / TH[4];
            lwdn[i] = tvl[3];
            lwup[i] = tvl[4];
        }
        // Set value limits
        rh[i] = (ez / satvapCpp(Tz[i])) * 100;
        if (rh[i] > 100) rh[i] = 100;
        double tmx = std::max({ tleaf[i], tc[i], T0[i], Tcan }) + 2.0;
        double tmn = std::min({ tleaf[i], tc[i], T0[i], Tcan }) - 2.0;
        if (Tz[i] > tmx) Tz[i] = tmx;
        if (Tz[i] < tmn) Tz[i] = tmn;
    }
    abovemodel out;
    out.Tz = Tz;
    out.tleaf = tleaf;
    out.rh = rh;
    out.lwdn = lwdn;
    out.lwup = lwup;
    return out;
}
NumericVector mxtccalc(NumericVector tc) {
    IntegerVector dims = tc.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    NumericVector mxtc(tc.size());
    tc = aperm3D2(tc, rows, cols, tsteps);
    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // subset climate variables
            int st = j * tsteps + i * tsteps * cols;
            int ed = st + tsteps - 1;
            // subset to vector
            std::vector<double> tcv = tosvd(tc[Range(st, ed)]);
            double mx = -273.15;
            for (int k = 0; k < tsteps; ++k) if (tcv[k] > mx) mx = tcv[k];
            for (int k = 0; k < tsteps; ++k) {
                mxtc[index] = mx;
                index++;
            }
        }
    }
    mxtc.attr("dim") = NumericVector::create(tsteps, rows, cols);
    // Reshape results
    mxtc = aperm3D(mxtc, rows, cols, tsteps);
    return mxtc;
}
// microclimate model as grid above ground
// [[Rcpp::export]]
List abovegrid(double reqhgt, List micro)
{
    // Access items of micro
    double zref = micro["zref"];
    NumericVector zen = micro["zen"];
    // weather
    NumericVector tc = micro["tc"];
    NumericVector pk = micro["pk"];
    NumericVector ea = micro["ea"];
    NumericVector es = micro["estl"];
    NumericVector tdew = micro["tdew"];
    NumericVector difr = micro["difr"];
    NumericVector dirr = micro["dirr"];
    NumericVector lwdn = micro["lwdown"];
    // microclimate
    NumericVector T0 = micro["Tg"];
    NumericVector uz = micro["uz"];
    NumericVector uf = micro["uf"];
    NumericVector gHa = micro["gHa"];
    NumericVector radCsw = micro["radCsw"];
    NumericVector radClw = micro["radClw"];
    NumericVector radLsw = micro["radLsw"];
    NumericVector Rddown = micro["Rddown"];
    NumericVector Rbdown = micro["Rbdown"];
    // soil variables
    NumericVector soilm = micro["soilm"];
    NumericVector G = micro["G"];
    // Habitat variables
    NumericVector hgt = micro["veghgt"];
    NumericVector pai = micro["pai"];
    NumericVector paia = micro["paia"];
    NumericVector leafd = micro["leafd"];
    NumericVector leafden = micro["leafden"];
    NumericVector gsmax = micro["gsmax"];
    // Soil variables
    NumericVector Smin = micro["Smin"];
    NumericVector Smax = micro["Smax"];
    // Create variable
    NumericVector Tz(tc.size());
    NumericVector tleaf(tc.size());
    NumericVector rh(tc.size());
    NumericVector lwup(tc.size());
    // Calculate mxtc
    NumericVector mxtc = mxtccalc(tc);
    // Run model
    for (R_xlen_t i = 0; i < tc.size(); ++i) {
        // Calculate ground and surface wetness
        double eT = satvapCpp(T0[i]) - ea[i];
        if (eT < 0.001) eT = 0.001;
        double plf = 0.8753 - 1.7126 * log(eT);
        double gwet = 1 / (1 + exp(-plf));
        double surfwet = (soilm[i] - Smin[i]) / (Smax[i] - Smin[i]);
        if (surfwet > gwet) gwet = surfwet;
        // Calculate vapour conductivity
        double Rsw = difr[i] + dirr[i] * cos(zen[i]);
        double gv = stomcondCpp(Rsw, gsmax[i], 300) * 3;
        double gV = 1 / (1 / gHa[i] + 1 / gv);
        // Calculate canopy temperature
        double Rabs = radCsw[i] + radClw[i];
        std::vector<double> TH = PenmanMonteith2Cpp(Rabs, gHa[i], gV, tc[i], mxtc[i], pk[i], ea[i], es[i], G[i], surfwet, tdew[i]);
        double Tcan = TH[0];
        double ez = 0;
        if (reqhgt >= hgt[i]) {
            std::vector<double> tv = TVabove(reqhgt, zref, hgt[i], pai[i], Tcan, tc[i], ea[i], surfwet);
            Tz[i] = tv[0];
            tleaf[i] = Tcan;
            ez = tv[1];
            lwup[i] = 0.97 * 5.67 * pow(10, -8) * pow(Tcan + 273.15, 4);
        }
        else {
            // Calculate Leaf temperature
            std::vector<double> tvl = leaftemp(Tcan, T0[i], tc[i], mxtc[i], pk[i], ea[i], es[i], uz[i], tdew[i], surfwet, zen[i],
                radLsw[i], Rddown[i], Rbdown[i], lwdn[i], pai[i], paia[i], leafd[i], gsmax[i]);
            tleaf[i] = tvl[0];
            // Calculate temperature below canopy
            double Flux = TH[1] * (1 - exp(-pai[i]));
            double Fluxz = tvl[1];
            std::vector<double> tv = TVabove(hgt[i], zref, hgt[i], pai[i], Tcan, tc[i], ea[i], surfwet);
            double SH = tv[0] * 29.3 * 43;
            double SG = T0[i] * 29.3 * 43;
            double d = zeroplanedisCpp(hgt[i], pai[i]);
            Tz[i] = TVbelow(zref, reqhgt, d, hgt[i], pai[i], uf[i], leafden[i], Flux, Fluxz, SH, SG) / (29.3 * 43);
            // Calculate humidity below canopy
            Flux = TH[2] * (1 - exp(-pai[i]));
            Fluxz = tvl[2];
            SH = tv[1] * TH[4];
            SG = satvapCpp(T0[i]) * gwet * TH[4];
            ez = TVbelow(zref, reqhgt, d, hgt[i], pai[i], uf[i], leafden[i], Flux, Fluxz, SH, SG) / TH[4];
            lwdn[i] = tvl[3];
            lwup[i] = tvl[4];
        }
        // Set value limits
        rh[i] = (ez / satvapCpp(Tz[i])) * 100;
        if (rh[i] > 100) rh[i] = 100;
        double tmx = std::max({ tleaf[i], tc[i], T0[i], Tcan }) + 2.0;
        double tmn = std::min({ tleaf[i], tc[i], T0[i], Tcan }) - 2.0;
        if (Tz[i] > tmx) Tz[i] = tmx;
        if (Tz[i] < tmn) Tz[i] = tmn;
    }
    // Get dimensions of arrays
    IntegerVector dims = tc.attr("dim");
    int n1 = dims[0];
    int n2 = dims[1];
    int n3 = dims[2];
    // reshape results
    Tz.attr("dim") = NumericVector::create(n1, n2, n3);
    tleaf.attr("dim") = NumericVector::create(n1, n2, n3);
    rh.attr("dim") = NumericVector::create(n1, n2, n3);
    lwdn.attr("dim") = NumericVector::create(n1, n2, n3);
    lwup.attr("dim") = NumericVector::create(n1, n2, n3);
    // Create output
    List mout;
    mout["Tz"] = Rcpp::wrap(Tz);
    mout["tleaf"] = Rcpp::wrap(tleaf);
    mout["T0"] = Rcpp::wrap(T0);
    mout["soilm"] = Rcpp::wrap(soilm);
    mout["relhum"] = Rcpp::wrap(rh);
    mout["windspeed"] = Rcpp::wrap(uz);
    mout["Rdirdown"] = Rcpp::wrap(Rbdown);
    mout["Rdifdown"] = Rcpp::wrap(Rddown);
    mout["Rlwdown"] = Rcpp::wrap(lwdn);
    mout["Rswup"] = Rcpp::wrap(micro["Rdup"]);
    mout["Rlwup"] = Rcpp::wrap(lwup);
    return mout;
}
// microclimate model as point above ground
NumericVector abovepoint(double reqhgt, double zref, double tc, double pk, double relhum, double u2, 
    double Rsw, double Rdif, double Rlw, double hgt, double pai, double paia, double lref, double ltra, 
    double clump, double leafd, double leafden, double albg, double albc, double snowtempg, double snowtempc, double zenr,
    double shadowmask, double si, double svfa, double ws, double umu, double mxtc)
{
    // Calculate base variables
    double es = satvapCpp(tc);
    double ea = es * relhum / 100.0;
    double tdew = dewpointCpp(tc, ea);
    // Calculate two-stream variables
    NumericVector tsv = twostreampoint(reqhgt, hgt, pai, paia, lref, ltra, clump, albg, albc, Rsw, Rdif,
        shadowmask, tc, Rlw, zenr, si, svfa);
    double Rbdown = tsv[0];
    double Rddown = tsv[1];
    double Rdup = tsv[2];
    double radLsw = tsv[3];
    // Calculate wind speed
    NumericVector uzv = windpoint(reqhgt, u2, hgt, pai, ws, umu, zref);
    double uz = uzv[0];
    double uf = uzv[1];
    double gHa = uzv[2];
    // Initalise variables
    double Tz = 0.0;
    double tleaf = snowtempc;
    double ez = ea;
    double lwup = 0.97 * 5.67 * pow(10, -8) * pow(snowtempc + 273.15, 4);
    double lwdn = Rlw;
    // Above canopy
    if (reqhgt >= hgt) {
        std::vector<double> tv = TVabove(reqhgt, zref, hgt, pai, snowtempc, tc, ea, 1.0);
        Tz = tv[0];
        ez = tv[1];
    }
    else {
        // Calculate Leaf temperature
        std::vector<double> tvl = leaftemp(snowtempc, snowtempc, tc, mxtc, pk, ea, es, uz, tdew, 1.0, zenr,
            radLsw, Rddown, Rbdown, Rlw, pai, paia, leafd, 1000.0);
        tleaf = tvl[0];
        // Calculate temperature below canopy
        double H = 29.3 * gHa * (snowtempc - tc);
        double la = 45068.7 - 42.8428 * tc;
        if (tc < 0) {
            la = 51078.69 - 4.338 * tc - 0.06367 * tc * tc;
        }
        double m = la * (gHa / pk);
        double L = m * (es - ea);
        double Flux = H * (1 - exp(-pai));
        double Fluxz = tvl[1];
        std::vector<double> tv = TVabove(hgt, zref, hgt, pai, snowtempc, tc, ea, 1.0);
        double SH = tv[0] * 29.3 * 43;
        double SG = snowtempg * 29.3 * 43;
        double d = zeroplanedisCpp(hgt, pai);
        Tz = TVbelow(zref, reqhgt, d, hgt, pai, uf, leafden, Flux, Fluxz, SH, SG) / (29.3 * 43);
        // Calculate humidity below canopy
        Flux = L * (1 - exp(-pai));
        Fluxz = tvl[2];
        double TH4 = la * (43 / pk);
        SH = tv[1] * TH4;
        SG = satvapCpp(snowtempg) * TH4;
        ez = TVbelow(zref, reqhgt, d, hgt, pai, uf, leafden, Flux, Fluxz, SH, SG) / TH4;
        lwdn = tvl[3];
        lwup = tvl[4];
    }
    // Set value limits
    double rh = (ez / satvapCpp(Tz)) * 100.0;
    if (rh > 100) rh = 100.0;
    double tmx = std::max({ tleaf, tc, snowtempg, snowtempc }) + 2.0;
    double tmn = std::min({ tleaf, tc, snowtempg, snowtempc }) - 2.0;
    if (Tz > tmx) Tz = tmx;
    if (Tz < tmn) Tz = tmn;
    // Create output
    NumericVector mout(9);
    mout[0] = Tz;
    mout[1] = tleaf;
    mout[2] = rh;
    mout[3] = uz;
    mout[4] = Rbdown;
    mout[5] = Rddown;
    mout[6] = lwdn;
    mout[7] = Rdup;
    mout[8] = lwup;
    return mout;
}
// Calculate temperature below ground
std::vector<double> Tbelowgroundv(double reqhgt, std::vector<double> Tg, std::vector<double> Tgp, 
    std::vector<double> Tbp, double meanD, double mat, int hiy, bool complete)
{
    std::vector<double> Tz = Tg;
    if (reqhgt < 0) {
        // Calculate n (smoothing parameter)
        double nb = -118.35 * reqhgt / meanD;
        int n = static_cast<int>(std::round(nb));
        if (complete) {
            if (n < static_cast<int>(Tg.size())) {
                Tz = manCpp(Tg, n);
            }
            else {
                double sumT = 0;
                for (size_t i = 0; i < Tg.size(); ++i) sumT = sumT + Tg[i];
                double meanT = sumT / Tg.size();
                for (size_t i = 0; i < Tg.size(); ++i) Tz[i] = meanT;
            }
        }
        // Time sequence incomplete
        else {
            // Calculate Tz at daily damping depth
            // ** Create variable
            std::vector<double> Tzd(Tg.size());
            // ** Calculate mean daily Tbp and difference from this
            std::vector<double> Tbpd = hourtodayCpp(Tbp, "mean");
            std::vector<double> Tbpa(Tg.size());
            for (size_t i = 0; i < Tbp.size(); ++i) Tbpa[i] = Tbp[i] - Tbpd[i];
            // ** Calculate grid to point ratio and difference
            std::vector<double> dTggmx = hourtodayCpp(Tg, "max");
            std::vector<double> dTggmn = hourtodayCpp(Tg, "min");
            std::vector<double> dTggme = hourtodayCpp(Tg, "mean");
            std::vector<double> dTgpmx = hourtodayCpp(Tgp, "max");
            std::vector<double> dTgpmn = hourtodayCpp(Tgp, "min");
            std::vector<double> dTgpme = hourtodayCpp(Tgp, "mean");
            double sat = 0.0;
            for (size_t i = 0; i < Tg.size(); ++i) {
                double rat = (dTggmx[i] - dTggmn[i]) / (dTgpmx[i] - dTgpmn[i]);
                double dif = dTggme[i] - dTgpme[i];
                Tzd[i] = rat * Tbpa[i] + Tbpd[i] + dif;
                sat = sat + Tzd[i];
            }
            if (nb > 1.0 && nb <= 24.0) {
                double w1 = 1.0 / nb;
                double w2 = nb / 24.0;
                double wgt = w1 / (w1 + w2);
                for (size_t i = 0; i < Tg.size(); ++i) Tz[i] = wgt * Tg[i] + (1 - wgt) * Tzd[i];
            }
            if (nb > 24.0) {
                // Calculate Tz at annual damping depth
                if (nb < hiy) {
                    // Calculate daily to annual weights
                    double w1 = 24.0 / nb;
                    double w2 = nb / hiy;
                    double wgt = w1 / (w1 + w2);
                    for (size_t i = 0; i < Tg.size(); ++i) Tz[i] = wgt * Tzd[i] + (1 - wgt) * mat;
                }
                else {
                    for (size_t i = 0; i < Tg.size(); ++i) Tz[i] = mat;
                } // end nb >= hiy
            } // end nb > 24.0
        } // end none complete
    }
    return Tz;
}
// Calculate below ground temperature: grid
// [[Rcpp::export]]
List belowgrid(double reqhgt, List micro, int hiy, bool complete)
{
    // Access required items from micro
    // weather
    NumericVector Tg = micro["Tg"];
    NumericVector Tgp = micro["Tgp"];
    NumericVector Tbp = micro["Tbp"];
    NumericVector DD = micro["kDDg"];
    double mat = micro["matemp"];
    // Get dimensions of arrays
    IntegerVector dims = Tg.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    // repermutate variables
    Tg = aperm3D2(Tg, rows, cols, tsteps);
    Tgp = aperm3D2(Tgp, rows, cols, tsteps);
    Tbp = aperm3D2(Tbp, rows, cols, tsteps);
    DD = aperm3D2(DD, rows, cols, tsteps);
    // Create output variables
    NumericVector Tz(Tg.size());
    // set indices for storing results
    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = Tg[index];
            if (!NumericVector::is_na(val)) {
                // subset climate variables
                int st = j * tsteps + i * tsteps * cols;
                int ed = st + tsteps - 1;
                // subset to vector
                std::vector<double> Tgv = tosvd(Tg[Range(st, ed)]);
                std::vector<double> Tgpv = tosvd(Tgp[Range(st, ed)]);
                std::vector<double> Tbpv = tosvd(Tbp[Range(st, ed)]);
                std::vector<double> DDv = tosvd(DD[Range(st, ed)]);
                // Calculate mean damping depth
                double sumD = 0.0;
                for (int k = 0; k < tsteps; k++) sumD = sumD + DDv[i];
                double meanD = sumD / static_cast<double>(DDv.size());
                std::vector<double> Tzv = Tbelowgroundv(reqhgt, Tgv, Tgpv, Tbpv,
                    meanD, mat, hiy, complete);
                for (int k = 0; k < tsteps; k++) Tz[index++] = Tzv[k];
            }
            else {
                for (int k = 0; k < tsteps; k++) Tz[index++] = NA_REAL;
            }
        }
    }
    // Shape results
    Tz.attr("dim") = NumericVector::create(tsteps, rows, cols);
    // Reshape results
    Tz = aperm3D(Tz, rows, cols, tsteps);
    // Create output
    List mout;
    mout["Tz"] = Rcpp::wrap(Tz);
    mout["T0"] = Rcpp::wrap(micro["Tg"]);
    mout["soilm"] = Rcpp::wrap(micro["soilm"]);
    return mout;
}
double belowpoint(double reqhgt, double meanD, double snowtempg, double Tzd, double Tza, double hiy) {
    double nb = -118.35 * reqhgt / meanD;
    double Tz = snowtempg;
    if (nb > 1.0) {
        if (nb <= 24.0) {
            double w1 = 1.0 / nb;
            double w2 = nb / 24.0;
            double wgt = w1 / (w1 + w2);
            Tz = wgt * snowtempg + (1 - wgt) * Tzd;
        }
        else {
            if (nb <= hiy) {
                double w1 = 24.0 / nb;
                double w2 = nb / hiy;
                double wgt = w1 / (w1 + w2);
                Tz = wgt * Tzd + (1 - wgt) * Tza;
            }
            else {
                Tz = Tza;
            } // end nb > hiy
        } // end nb > 24
    } // end nb > 1.0
    return Tz;
}
// Run microclimate model (vector, hourly)
microm runmicrovhCpp(double reqhgt, double zref, std::vector<double> si, std::vector<double> zen,
    std::vector<double> tc, std::vector<double> pk, std::vector<double> u2, 
    std::vector<double> Rsw, std::vector<double> Rdif, std::vector<double> lwdown, 
    std::vector<double> ea, std::vector<double> es, std::vector<double> tdew, 
    std::vector<double> soilm, std::vector<double> umu, 
    std::vector<double> T0p, std::vector<double> Tgp, std::vector<double> Tbp, 
    std::vector<double> Gp, std::vector<double> kp, std::vector<double> dtrp, std::vector<double> muGp, 
    double hgt, double pai, double x, double lref, double ltra, double clump, double leafd, 
    double gsmax, double Smin, double Smax, double gref, double soilb, double Psie, double Vq,
    double Vm, double Mc, double rho, std::vector<double> ws, double tadd, double paia, double leafden, 
    double mxtc, double svfa, int hiy, bool complete, double mat)
{
    double reqhgt2 = reqhgt;
    if (reqhgt2 == 0) reqhgt2 = 0.00001;
    // Calculate distributed soil moisture
    std::vector<double> soild = soildCppv(soilm, Smin, Smax, tadd);
    // Calculate radiation fluxes
    // ** Calculate two-stream coefficients vector ** //
    radmodel2 radm = twostreamvCpp(reqhgt2, hgt, pai, paia, x, lref, ltra, clump, gref, svfa, si, zen, Rsw, Rdif, tc, lwdown);
    // Calculate wind speed
    windmodel windm = windvCpp(reqhgt2, zref, u2, umu, hgt, pai, ws);
    // Calculate ground surface temperature
    std::vector<double> soilparamsm = { Smax, Smin, soilb, Psie, Vq, Vm, Mc, rho };
    soilmodel smod = soiltemp_hrCpp(tc, ea, es, tdew, pk, radm.radGsw, radm.radGlw, windm.gHa,
        soilparamsm, Gp, kp, dtrp, muGp, soild, mxtc);
    // Initalise output variables
    microm out;
    out.soilm = soild;
    // Calculate ground surface temperature
    if (reqhgt2 > 0) {
        abovemodel tva = TVabovegroundv(reqhgt, zref, tc, pk, ea, es, tdew, Rsw, lwdown, smod.Tg,
            windm.uz, windm.uf, windm.gHa, zen, radm.radCsw, radm.radClw, radm.radLsw, radm.Rddown,
            radm.Rbdown, soild, smod.G, hgt, pai, paia, leafd, leafden, Smin, Smax, gsmax, mxtc);
        if (reqhgt == 0) {
            out.Tz = smod.Tg;
            out.tleaf = { 0,0 };
            out.relhum = { 0,0 };
            out.uz = { 0,0 };
        }
        else {
            out.Tz = tva.Tz;
            out.tleaf = tva.tleaf;
            out.relhum = tva.rh;
            out.uz = windm.uz;
        }
        out.Rlwdown = tva.lwdn;
        out.Rlwup = tva.lwup;
        out.Rdirdown = radm.Rbdown;
        out.Rdifdown = radm.Rddown;
        out.Rswup = radm.Rdup;
    }
    else {
        // Calculate meanD
        double SumD = 0.0;
        for (size_t i = 0; i < smod.kDDg.size(); i++) SumD = SumD + smod.kDDg[i];
        double meanD = SumD / static_cast<double>(smod.kDDg.size());
        out.Tz = Tbelowgroundv(reqhgt, smod.Tg, Tgp, Tbp, meanD, mat, hiy, complete);
        out.tleaf = { 0,0 };
        out.relhum = { 0,0 };
        out.Rlwdown = { 0,0 };
        out.Rlwup = { 0,0 };
        out.Rdirdown = { 0,0 };
        out.Rdifdown = { 0,0 };
        out.Rswup = { 0,0 };
    }
    return out;
}
// Run microclimate model(hourly, static vegetation, data.frame climate input)
// [[Rcpp::export]]
List runmicro1Cpp(DataFrame obstime, DataFrame climdata, DataFrame pointm, List vegp, List soilc,
    double reqhgt, double zref, double lat, double lon, double Sminp, double Smaxp, double tfact, 
    bool complete, double mat, std::vector<bool> out)
{
    // Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access columns of climdata
    std::vector<double> tc = climdata["temp"];
    std::vector<double> es = climdata["es"];
    std::vector<double> ea = climdata["ea"];
    std::vector<double> tdew = climdata["tdew"];
    std::vector<double> pk = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> lwdown = climdata["lwdown"];
    std::vector<double> u2 = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    // Access columns of pointm
    std::vector<double> soilmp = pointm["soilm"];
    std::vector<double> Tgp = pointm["Tg"];
    std::vector<double> T0p = pointm["T0p"];
    std::vector<double> Tbp = pointm["Tbp"];
    std::vector<double> Gp = pointm["G"];
    std::vector<double> DDp = pointm["DDp"];
    std::vector<double> umu = pointm["umu"];
    std::vector<double> kp = pointm["kp"];
    std::vector<double> muGp = pointm["muGp"];
    std::vector<double> dtrp = pointm["dtrp"];
    // Access items from vegp
    NumericMatrix hgt = vegp["hgt"];
    NumericMatrix pai = vegp["pai"];
    NumericMatrix x = vegp["x"];
    NumericMatrix gsmax = vegp["gsmax"];
    NumericMatrix lref = vegp["leafr"];
    NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"];
    NumericMatrix leafd = vegp["leafd"];
    // ** additional
    NumericMatrix paia = vegp["paia"];
    NumericMatrix leafden = vegp["leafden"];
    // Access items from soilc
    NumericMatrix Smin = soilc["Smin"]; 
    NumericMatrix Smax = soilc["Smax"];
    NumericMatrix gref = soilc["gref"];
    NumericMatrix soilb = soilc["soilb"];
    NumericMatrix Psie = soilc["Psie"]; 
    NumericMatrix Vq = soilc["Vq"];
    NumericMatrix Vm = soilc["Vm"]; 
    NumericMatrix Mc = soilc["Mc"];
    NumericMatrix rho = soilc["rho"]; 
    // *** additional
    NumericMatrix slope = soilc["slope"];
    NumericMatrix aspect = soilc["aspect"]; 
    NumericMatrix twi = soilc["twi"];
    NumericMatrix svfa = soilc["svfa"];
    NumericVector wsa = soilc["wsa"]; 
    NumericVector hor = soilc["hor"];
    // Get dimensions
    int rows = hgt.nrow();
    int cols = hgt.ncol();
    int tsteps = tc.size();
    // Calculate zenith and humidity variables
    std::vector<double> zend(tc.size());
    std::vector<double> azid(tc.size());
    std::vector<int> sindex(tsteps);
    for (int k = 0; k < tsteps; ++k) {
        std::vector<double> sp = solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
        zend[k] = sp[0];
        azid[k] = sp[1];
        sindex[k] = static_cast<int>(std::round(azid[k] / 15)) % 24;
    }
    // Calculate other odds and sods
    int hiy = 365 * 24;
    if (year[0] % 4 == 0) hiy = 366 * 24;
    double mxtc = -273.15;
    for (int k = 0; k < tsteps; ++k) if (tc[k] > mxtc) mxtc = tc[k];
    // Distribute soil moisture
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    int index = 0;
    // Compute wind index
    std::vector<int> windex(tsteps);
    for (int k = 0; k < tsteps; ++k) windex[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    // Create output variables
    NumericVector Tz;
    NumericVector tleaf;
    NumericVector relhum;
    NumericVector soilm;
    NumericVector uz;
    NumericVector Rdirdown;
    NumericVector Rdifdown;
    NumericVector Rlwdown;
    NumericVector Rswup;
    NumericVector Rlwup;
    // Conditionally assign memory to output variables
    if (out[0]) Tz = NumericVector(rows * cols * tsteps);
    if (out[1]) tleaf = NumericVector(rows * cols * tsteps);
    if (out[2]) relhum = NumericVector(rows * cols * tsteps);
    if (out[3]) soilm = NumericVector(rows * cols * tsteps);
    if (out[4]) uz = NumericVector(rows * cols * tsteps);
    if (out[5]) Rdirdown = NumericVector(rows * cols * tsteps);
    if (out[6]) Rdifdown = NumericVector(rows * cols * tsteps);
    if (out[7]) Rlwdown = NumericVector(rows * cols * tsteps);
    if (out[8]) Rswup = NumericVector(rows * cols * tsteps);
    if (out[9]) Rlwup = NumericVector(rows * cols * tsteps);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                std::vector<double> si(tsteps);
                std::vector<double> zen(tsteps);
                std::vector<double> ws(tsteps);
                // Calculate solar index etc 
                for (int k = 0; k < tsteps; k++) {
                    si[k] = solarindexCpp(slope(i, j), aspect(i, j), zend[k], azid[k], true);
                    if (si[k] < 0.0) si[k] = 0.0;
                    zen[k] = zend[k] * M_PI / 180;
                    ws[k] = wsa[windex[k] * rows * cols + j * rows + i];
                    double ha = hor[sindex[k] * rows * cols + j * rows + i];
                    double sa = 90 - zend[k];
                    if (ha > tan(sa * M_PI / 180)) si[k] = 0.0;
                }
                microm micro = runmicrovhCpp(reqhgt, zref, si, zen, tc, pk, u2, Rsw, Rdif, lwdown,
                    ea, es, tdew, soilmp, umu, T0p, Tgp, Tbp, Gp, kp, dtrp, muGp,
                    hgt(i, j), pai(i, j), x(i, j), lref(i, j), ltra(i, j), clump(i, j), leafd(i, j),
                    gsmax(i, j), Smin(i, j), Smax(i, j), gref(i, j), soilb(i, j), Psie(i, j), Vq(i, j),
                    Vm(i, j), Mc(i, j), rho(i, j), ws, tadd(i, j), paia(i, j), leafden(i, j),
                    mxtc, svfa(i, j), hiy, complete, mat);
                for (int k = 0; k < tsteps; k++) {
                    if (out[0]) Tz[index] = micro.Tz[k];
                    if (out[1]) tleaf[index] = micro.tleaf[k];
                    if (out[2]) relhum[index] = micro.relhum[k];
                    if (out[3]) soilm[index] = micro.soilm[k];
                    if (out[4]) uz[index] = micro.uz[k];
                    if (out[5]) Rdirdown[index] = micro.Rdirdown[k];
                    if (out[6]) Rdifdown[index] = micro.Rdifdown[k];
                    if (out[7]) Rlwdown[index] = micro.Rlwdown[k];
                    if (out[8]) Rswup[index] = micro.Rswup[k];
                    if (out[9]) Rlwup[index] = micro.Rlwup[k];
                    index++;
                }
            }
            else {
                for (int k = 0; k < tsteps; k++) {
                    if (out[0]) Tz[index] = NA_REAL;
                    if (out[1]) tleaf[index] = NA_REAL;
                    if (out[2]) relhum[index] = NA_REAL;
                    if (out[3]) soilm[index] = NA_REAL;
                    if (out[4]) uz[index] = NA_REAL;
                    if (out[5]) Rdirdown[index] = NA_REAL;
                    if (out[6]) Rdifdown[index] = NA_REAL;
                    if (out[7]) Rlwdown[index] = NA_REAL;
                    if (out[8]) Rswup[index] = NA_REAL;
                    if (out[9]) Rlwup[index] = NA_REAL;
                    index++;
                } // end for
            } // end else
        } // end col
    } // end row
    // reshape results
    if (out[0]) Tz.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[1]) tleaf.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[2]) relhum.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[3]) soilm.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[4]) uz.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[5]) Rdirdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[6]) Rdifdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[7]) Rlwdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[8]) Rswup.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[9]) Rlwup.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[0]) Tz = aperm3D(Tz, rows, cols, tsteps);
    if (out[1]) tleaf = aperm3D(tleaf, rows, cols, tsteps);
    if (out[2]) relhum = aperm3D(relhum, rows, cols, tsteps);
    if (out[3]) soilm = aperm3D(soilm, rows, cols, tsteps);
    if (out[4]) uz = aperm3D(uz, rows, cols, tsteps);
    if (out[5]) Rdirdown = aperm3D(Rdirdown, rows, cols, tsteps);
    if (out[6]) Rdifdown = aperm3D(Rdifdown, rows, cols, tsteps);
    if (out[7]) Rlwdown = aperm3D(Rlwdown, rows, cols, tsteps);
    if (out[8]) Rswup = aperm3D(Rswup, rows, cols, tsteps);
    if (out[9]) Rlwup = aperm3D(Rlwup, rows, cols, tsteps);
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Rcpp::wrap(Tz);
    if (out[1]) outp["tleaf"] = Rcpp::wrap(tleaf);
    if (out[2]) outp["relhum"] = Rcpp::wrap(relhum);
    if (out[3]) outp["soilm"] = Rcpp::wrap(soilm);
    if (out[4]) outp["windspeed"] = Rcpp::wrap(uz);
    if (out[5]) outp["Rdirdown"] = Rcpp::wrap(Rdirdown);
    if (out[6]) outp["Rdifdown"] = Rcpp::wrap(Rdifdown);
    if (out[7]) outp["Rlwdown"] = Rcpp::wrap(Rlwdown);
    if (out[8]) outp["Rswup"] = Rcpp::wrap(Rswup);
    if (out[9]) outp["Rlwup"] = Rcpp::wrap(Rlwup);
    return outp;
}
// Run microclimate model(hourly, static vegetation, array climate inputs)
// [[Rcpp::export]]
List runmicro2Cpp(DataFrame obstime, List climdata, List pointm, List vegp, List soilc,
    double reqhgt, double zref, NumericMatrix lats, NumericMatrix lons, double Sminp, double Smaxp, double tfact,
    bool complete, double mat, std::vector<bool> out)
{
    // Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access vectors of climdata
    NumericVector tc = climdata["tc"];
    NumericVector es = climdata["es"];
    NumericVector ea = climdata["ea"];
    NumericVector tdew = climdata["tdew"];
    NumericVector pk = climdata["pk"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector lwdown = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"];
    NumericVector wdir = climdata["winddir"]; // not an array
    // Access vectors of pointm
    NumericVector soilmp = pointm["soilm"];
    NumericVector Tgp = pointm["Tg"];
    NumericVector T0p = pointm["T0p"];
    NumericVector Tbp = pointm["Tbp"];
    NumericVector Gp = pointm["Gp"];
    NumericVector DDp = pointm["DDp"];
    NumericVector umu = pointm["umu"];
    NumericVector kp = pointm["kp"];
    NumericVector muGp = pointm["muGp"];
    NumericVector dtrp = pointm["dtrp"];
    // Access items from vegp
    NumericMatrix hgt = vegp["hgt"];
    NumericMatrix pai = vegp["pai"];
    NumericMatrix x = vegp["x"];
    NumericMatrix gsmax = vegp["gsmax"];
    NumericMatrix lref = vegp["leafr"];
    NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"];
    NumericMatrix leafd = vegp["leafd"];
    // ** additional
    NumericMatrix paia = vegp["paia"];
    NumericMatrix leafden = vegp["leafden"];
    // Access items from soilc
    NumericMatrix Smin = soilc["Smin"];
    NumericMatrix Smax = soilc["Smax"];
    NumericMatrix gref = soilc["gref"];
    NumericMatrix soilb = soilc["soilb"];
    NumericMatrix Psie = soilc["Psie"];
    NumericMatrix Vq = soilc["Vq"];
    NumericMatrix Vm = soilc["Vm"];
    NumericMatrix Mc = soilc["Mc"];
    NumericMatrix rho = soilc["rho"];
    // *** additional
    NumericMatrix slope = soilc["slope"];
    NumericMatrix aspect = soilc["aspect"];
    NumericMatrix twi = soilc["twi"];
    NumericMatrix svfa = soilc["svfa"];
    NumericVector wsa = soilc["wsa"];
    NumericVector hor = soilc["hor"];
    // Get dimensions
    int rows = hgt.nrow();
    int cols = hgt.ncol();
    int tsteps = year.size();
    // Calculate hours in year
    int hiy = 365 * 24;
    if (year[0] % 4 == 0) hiy = 366 * 24;
    // Distribute soil moisture
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    // Create output variables
    NumericVector Tz;
    NumericVector tleaf;
    NumericVector relhum;
    NumericVector soilm;
    NumericVector uz;
    NumericVector Rdirdown;
    NumericVector Rdifdown;
    NumericVector Rlwdown;
    NumericVector Rswup;
    NumericVector Rlwup;
    // Conditionally assign memory to output variables
    if (out[0]) Tz = NumericVector(rows * cols * tsteps);
    if (out[1]) tleaf = NumericVector(rows * cols * tsteps);
    if (out[2]) relhum = NumericVector(rows * cols * tsteps);
    if (out[3]) soilm = NumericVector(rows * cols * tsteps);
    if (out[4]) uz = NumericVector(rows * cols * tsteps);
    if (out[5]) Rdirdown = NumericVector(rows * cols * tsteps);
    if (out[6]) Rdifdown = NumericVector(rows * cols * tsteps);
    if (out[7]) Rlwdown = NumericVector(rows * cols * tsteps);
    if (out[8]) Rswup = NumericVector(rows * cols * tsteps);
    if (out[9]) Rlwup = NumericVector(rows * cols * tsteps);
    // repermutate climate variables
    tc = aperm3D2(tc, rows, cols, tsteps);
    es = aperm3D2(es, rows, cols, tsteps);
    ea = aperm3D2(ea, rows, cols, tsteps);
    tdew = aperm3D2(tdew, rows, cols, tsteps);
    pk = aperm3D2(pk, rows, cols, tsteps);
    Rsw = aperm3D2(Rsw, rows, cols, tsteps);
    Rdif = aperm3D2(Rdif, rows, cols, tsteps);
    lwdown = aperm3D2(lwdown, rows, cols, tsteps);
    u2 = aperm3D2(u2, rows, cols, tsteps);
    // repermutate point model variables
    soilmp = aperm3D2(soilmp, rows, cols, tsteps);
    Tgp = aperm3D2(Tgp, rows, cols, tsteps);
    T0p = aperm3D2(T0p, rows, cols, tsteps);
    Tbp = aperm3D2(Tbp, rows, cols, tsteps);
    Gp = aperm3D2(Gp, rows, cols, tsteps);
    DDp = aperm3D2(DDp, rows, cols, tsteps);
    umu = aperm3D2(umu, rows, cols, tsteps);
    kp = aperm3D2(kp, rows, cols, tsteps);
    muGp = aperm3D2(muGp, rows, cols, tsteps);
    dtrp = aperm3D2(dtrp, rows, cols, tsteps);
    // set indices for storing results
    int index = 0;
    // Compute wind index
    std::vector<int> windex(tsteps);
    for (int k = 0; k < tsteps; ++k) windex[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                // Calculate solar and wind shelter variables
                std::vector<double> zen(tsteps);
                std::vector<double> si(tsteps);
                std::vector<double> ws(tsteps);
                for (int k = 0; k < tsteps; ++k) {
                    ws[k] = wsa[windex[k] * rows * cols + j * rows + i];
                    std::vector<double> sp = solpositionCpp(lats(i, j), lons(i, j), year[k], month[k], day[k], hour[k]);
                    double zend = sp[0];
                    double azid = sp[1];
                    si[k] = solarindexCpp(slope(i, j), aspect(i, j), zend, azid);
                    if (si[k] < 0.0) si[k] = 0.0;
                    zen[k] = zend * M_PI / 180;
                    int sindex = static_cast<int>(std::round(azid / 15)) % 24;
                    double ha = hor[sindex * rows * cols + j * rows + i];
                    double sa = 90 - zend;
                    if (ha > tan(sa * M_PI / 180)) si[k] = 0.0;
                }
                int st = j * tsteps + i * tsteps * cols;
                int ed = st + tsteps - 1;
                // subset climate variables
                std::vector<double> tcv = tosvd(tc[Range(st, ed)]);
                std::vector<double> esv = tosvd(es[Range(st, ed)]);
                std::vector<double> eav = tosvd(ea[Range(st, ed)]);
                std::vector<double> tdewv = tosvd(tdew[Range(st, ed)]);
                std::vector<double> pkv = tosvd(pk[Range(st, ed)]);
                std::vector<double> Rswv = tosvd(Rsw[Range(st, ed)]);
                std::vector<double> Rdifv = tosvd(Rdif[Range(st, ed)]);
                std::vector<double> lwdownv = tosvd(lwdown[Range(st, ed)]);
                std::vector<double> u2v = tosvd(u2[Range(st, ed)]);
                // subset point model variables
                std::vector<double> soilmpv = tosvd(soilmp[Range(st, ed)]);
                std::vector<double> Tgpv = tosvd(Tgp[Range(st, ed)]);
                std::vector<double> T0pv = tosvd(T0p[Range(st, ed)]);
                std::vector<double> Tbpv = tosvd(Tbp[Range(st, ed)]);
                std::vector<double> Gpv = tosvd(Gp[Range(st, ed)]);
                std::vector<double> DDpv = tosvd(DDp[Range(st, ed)]);
                std::vector<double> umuv = tosvd(umu[Range(st, ed)]);
                std::vector<double> kpv = tosvd(kp[Range(st, ed)]);
                std::vector<double> muGpv = tosvd(muGp[Range(st, ed)]);
                std::vector<double> dtrpv = tosvd(dtrp[Range(st, ed)]);
                // calculate other odds and sods
                double mxtc = -273.15;
                for (int k = 0; k < tsteps; ++k) if (tcv[k] > mxtc) mxtc = tcv[k];
                microm micro = runmicrovhCpp(reqhgt, zref, si, zen, tcv, pkv, u2v, Rswv, Rdifv, lwdownv,
                    eav, esv, tdewv, soilmpv, umuv, T0pv, Tgpv, Tbpv, Gpv, kpv, dtrpv, muGpv,
                    hgt(i, j), pai(i, j), x(i, j), lref(i, j), ltra(i, j), clump(i, j), leafd(i, j),
                    gsmax(i, j), Smin(i, j), Smax(i, j), gref(i, j), soilb(i, j), Psie(i, j), Vq(i, j),
                    Vm(i, j), Mc(i, j), rho(i, j), ws, tadd(i, j), paia(i, j), leafden(i, j),
                    mxtc, svfa(i, j), hiy, complete, mat);
                for (int k = 0; k < tsteps; k++) {
                    if (out[0]) Tz[index] = micro.Tz[k];
                    if (out[1]) tleaf[index] = micro.tleaf[k];
                    if (out[2]) relhum[index] = micro.relhum[k];
                    if (out[3]) soilm[index] = micro.soilm[k];
                    if (out[4]) uz[index] = micro.uz[k];
                    if (out[5]) Rdirdown[index] = micro.Rdirdown[k];
                    if (out[6]) Rdifdown[index] = micro.Rdifdown[k];
                    if (out[7]) Rlwdown[index] = micro.Rlwdown[k];
                    if (out[8]) Rswup[index] = micro.Rswup[k];
                    if (out[9]) Rlwup[index] = micro.Rlwup[k];
                    index++;
                }
            }
            else {
                for (int k = 0; k < tsteps; k++) {
                    if (out[0]) Tz[index] = NA_REAL;
                    if (out[1]) tleaf[index] = NA_REAL;
                    if (out[2]) relhum[index] = NA_REAL;
                    if (out[3]) soilm[index] = NA_REAL;
                    if (out[4]) uz[index] = NA_REAL;
                    if (out[5]) Rdirdown[index] = NA_REAL;
                    if (out[6]) Rdifdown[index] = NA_REAL;
                    if (out[7]) Rlwdown[index] = NA_REAL;
                    if (out[8]) Rswup[index] = NA_REAL;
                    if (out[9]) Rlwup[index] = NA_REAL;
                    index++;
                } // end for
            } // end else
        } // end col
    } // end row
    // reshape results
    if (out[0]) Tz.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[1]) tleaf.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[2]) relhum.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[3]) soilm.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[4]) uz.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[5]) Rdirdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[6]) Rdifdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[7]) Rlwdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[8]) Rswup.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[9]) Rlwup.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[0]) Tz = aperm3D(Tz, rows, cols, tsteps);
    if (out[1]) tleaf = aperm3D(tleaf, rows, cols, tsteps);
    if (out[2]) relhum = aperm3D(relhum, rows, cols, tsteps);
    if (out[3]) soilm = aperm3D(soilm, rows, cols, tsteps);
    if (out[4]) uz = aperm3D(uz, rows, cols, tsteps);
    if (out[5]) Rdirdown = aperm3D(Rdirdown, rows, cols, tsteps);
    if (out[6]) Rdifdown = aperm3D(Rdifdown, rows, cols, tsteps);
    if (out[7]) Rlwdown = aperm3D(Rlwdown, rows, cols, tsteps);
    if (out[8]) Rswup = aperm3D(Rswup, rows, cols, tsteps);
    if (out[9]) Rlwup = aperm3D(Rlwup, rows, cols, tsteps);
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Rcpp::wrap(Tz);
    if (out[1]) outp["tleaf"] = Rcpp::wrap(tleaf);
    if (out[2]) outp["relhum"] = Rcpp::wrap(relhum);
    if (out[3]) outp["soilm"] = Rcpp::wrap(soilm);
    if (out[4]) outp["windspeed"] = Rcpp::wrap(uz);
    if (out[5]) outp["Rdirdown"] = Rcpp::wrap(Rdirdown);
    if (out[6]) outp["Rdifdown"] = Rcpp::wrap(Rdifdown);
    if (out[7]) outp["Rlwdown"] = Rcpp::wrap(Rlwdown);
    if (out[8]) outp["Rswup"] = Rcpp::wrap(Rswup);
    if (out[9]) outp["Rlwup"] = Rcpp::wrap(Rlwup);
    return outp;
}
// Run microclimate model(hourly, changing vegetation, dataframe climate inputs)
// [[Rcpp::export]]
List runmicro3Cpp(DataFrame dfsel, DataFrame obstime, DataFrame climdata, DataFrame pointm,
    List vegp, List soilc, double reqhgt, double zref, double lat, double lon, double Sminp,
    double Smaxp, double tfact, bool complete, double mat, std::vector<bool> out)
{
    // Extract data from dfsel
    std::vector<int> lyr = dfsel["lyr"];
    std::vector<int> st = dfsel["st"];
    std::vector<int> ed = dfsel["ed"];
    // First increment of vegetation layers
    // ** create a subset vector
    int n = ed[0] - st[0] + 1;
    IntegerVector s(n);
    for (int i = 0; i < n; ++i) s[i] = st[0] + i;
    // ** subset the data.frames
    DataFrame obstimeo = ssdf(obstime, s);
    DataFrame climdatao = ssdf(climdata, s);
    DataFrame pointmo = ssdf(pointm, s);
    // ** subset the vegetation layers
    NumericMatrix hgto = subsetArray(vegp["hgt"], 0);
    NumericMatrix paio = subsetArray(vegp["pai"], 0);
    NumericMatrix xo = subsetArray(vegp["x"], 0);
    NumericMatrix gsmaxo = subsetArray(vegp["gsmax"], 0);
    NumericMatrix lrefo = subsetArray(vegp["leafr"], 0);
    NumericMatrix ltrao = subsetArray(vegp["leaft"], 0);
    NumericMatrix clumpo = subsetArray(vegp["clump"], 0);
    NumericMatrix leafdo = subsetArray(vegp["leafd"], 0);
    // ** additional
    NumericMatrix paiao = subsetArray(vegp["paia"], 0);
    NumericMatrix leafdeno = subsetArray(vegp["leafden"], 0);
    // wrap into list
    List vegpl;
    vegpl["hgt"] = Rcpp::wrap(hgto);
    vegpl["pai"] = Rcpp::wrap(paio);
    vegpl["x"] = Rcpp::wrap(xo);
    vegpl["gsmax"] = Rcpp::wrap(gsmaxo);
    vegpl["leafr"] = Rcpp::wrap(lrefo);
    vegpl["leaft"] = Rcpp::wrap(ltrao);
    vegpl["clump"] = Rcpp::wrap(clumpo);
    vegpl["leafd"] = Rcpp::wrap(leafdo);
    // ** additional
    vegpl["paia"] = Rcpp::wrap(paiao);
    vegpl["leafden"] = Rcpp::wrap(leafdeno);
    // Create output variables
    NumericVector Tz;
    NumericVector tleaf;
    NumericVector relhum;
    NumericVector soilm;
    NumericVector windspeed;
    NumericVector Rdirdown;
    NumericVector Rdifdown;
    NumericVector Rlwdown;
    NumericVector Rswup;
    NumericVector Rlwup;
    // ** run the model for first increment of vegetation layer
    // ** run the model for first increment of vegetation layer
    List micro = runmicro1Cpp(obstimeo, climdatao, pointmo, vegpl, soilc, reqhgt, zref, lat, lon,
        Sminp, Smaxp, tfact, complete, mat, out);
    // Extract data
    if (out[0]) Tz = micro["Tz"];
    if (out[1]) tleaf = micro["tleaf"];
    if (out[2]) relhum = micro["relhum"];
    if (out[3]) soilm = micro["soilm"];
    if (out[4]) windspeed = micro["windspeed"];
    if (out[5]) Rdirdown = micro["Rdirdown"];
    if (out[6]) Rdifdown = micro["Rdifdown"];
    if (out[7]) Rlwdown = micro["Rlwdown"];
    if (out[8]) Rswup = micro["Rswup"];
    if (out[9]) Rlwup = micro["Rlwup"];
    for (size_t ilyr = 1; ilyr < lyr.size(); ++ilyr) {
        // ** create a subset vector
        int n = ed[ilyr] - st[ilyr] + 1;
        IntegerVector s(n);
        for (int i = 0; i < n; ++i) s[i] = st[ilyr] + i;
        // ** subset the data.frames
        obstimeo = ssdf(obstime, s);
        climdatao = ssdf(climdata, s);
        pointmo = ssdf(pointm, s);
        // ** subset the vegetation layers
        hgto = subsetArray(vegp["hgt"], ilyr);
        paio = subsetArray(vegp["pai"], ilyr);
        xo = subsetArray(vegp["x"], ilyr);
        gsmaxo = subsetArray(vegp["gsmax"], ilyr);
        lrefo = subsetArray(vegp["leafr"], ilyr);
        ltrao = subsetArray(vegp["leaft"], ilyr);
        clumpo = subsetArray(vegp["clump"], ilyr);
        leafdo = subsetArray(vegp["leafd"], ilyr);
        // ** additional
        paiao = subsetArray(vegp["paia"], ilyr);
        leafdeno = subsetArray(vegp["leafden"], ilyr);
        // Turn into list
        vegpl["hgt"] = Rcpp::wrap(hgto);
        vegpl["pai"] = Rcpp::wrap(paio);
        vegpl["x"] = Rcpp::wrap(xo);
        vegpl["gsmax"] = Rcpp::wrap(gsmaxo);
        vegpl["leafr"] = Rcpp::wrap(lrefo);
        vegpl["leaft"] = Rcpp::wrap(ltrao);
        vegpl["clump"] = Rcpp::wrap(clumpo);
        vegpl["leafd"] = Rcpp::wrap(leafdo);
        // ** additional
        vegpl["paia"] = Rcpp::wrap(paiao);
        vegpl["leafden"] = Rcpp::wrap(leafdeno);
        micro = runmicro1Cpp(obstimeo, climdatao, pointmo, vegpl, soilc, reqhgt, zref, lat, lon,
            Sminp, Smaxp, tfact, complete, mat, out);
        // abind layers
        if (out[0]) Tz = abind3D(Tz, micro["Tz"]);
        if (out[1]) tleaf = abind3D(tleaf, micro["tleaf"]);
        if (out[2]) relhum = abind3D(relhum, micro["relhum"]);
        if (out[3]) soilm = abind3D(soilm, micro["soilm"]);
        if (out[4]) windspeed = abind3D(windspeed, micro["windspeed"]);
        if (out[5]) Rdirdown = abind3D(Rdirdown, micro["Rdirdown"]);
        if (out[6]) Rdifdown = abind3D(Rdifdown, micro["Rdifdown"]);
        if (out[7]) Rlwdown = abind3D(Rlwdown, micro["Rlwdown"]);
        if (out[8]) Rswup = abind3D(Rswup, micro["Rswup"]);
        if (out[9]) Rlwup = abind3D(Rlwup, micro["Rlwup"]);
    }
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Rcpp::wrap(Tz);
    if (out[1]) outp["tleaf"] = Rcpp::wrap(tleaf);
    if (out[2]) outp["relhum"] = Rcpp::wrap(relhum);
    if (out[3]) outp["soilm"] = Rcpp::wrap(soilm);
    if (out[4]) outp["windspeed"] = Rcpp::wrap(windspeed);
    if (out[5]) outp["Rdirdown"] = Rcpp::wrap(Rdirdown);
    if (out[6]) outp["Rdifdown"] = Rcpp::wrap(Rdifdown);
    if (out[7]) outp["Rlwdown"] = Rcpp::wrap(Rlwdown);
    if (out[8]) outp["Rswup"] = Rcpp::wrap(Rswup);
    if (out[9]) outp["Rlwup"] = Rcpp::wrap(Rlwup);
    return outp;
}
// Run microclimate model(hourly, variable vegetation, array climate inputs)
// [[Rcpp::export]]
List runmicro4Cpp(DataFrame dfsel, DataFrame obstime, List climdata, List pointm, List vegp,
    List soilc, double reqhgt, double zref, NumericMatrix lats, NumericMatrix lons, double Sminp,
    double Smaxp, double tfact, bool complete, double mat, std::vector<bool> out)
{
    // Extract data from dfsel
    std::vector<int> lyr = dfsel["lyr"];
    std::vector<int> st = dfsel["st"];
    std::vector<int> ed = dfsel["ed"];
    // First increment of vegetation layers
    // ** create a subset vector
    int n = ed[0] - st[0] + 1;
    IntegerVector s(n);
    for (int i = 0; i < n; ++i) s[i] = st[0] + i;
    // ** subset the data.frame
    DataFrame obstimeo = ssdf(obstime, s);
    // Extract and subset climate arrays
    NumericVector tco = slicea(climdata["tc"], s);
    NumericVector eso = slicea(climdata["es"], s);
    NumericVector eao = slicea(climdata["ea"], s);
    NumericVector tdewo = slicea(climdata["tdew"], s);
    NumericVector pko = slicea(climdata["pk"], s);
    NumericVector Rswo = slicea(climdata["swdown"], s);
    NumericVector Rdifo = slicea(climdata["difrad"], s);
    NumericVector lwdowno = slicea(climdata["lwdown"], s);
    NumericVector u2o = slicea(climdata["windspeed"], s);
    NumericVector wdir = climdata["winddir"]; // not an array
    NumericVector wdiro(n);
    for (int i = 0; i < n; ++i) wdiro[i] = wdir[i + st[0]];
    // wrap into list
    List climdatao;
    climdatao["tc"] = Rcpp::wrap(tco);
    climdatao["es"] = Rcpp::wrap(eso);
    climdatao["ea"] = Rcpp::wrap(eao);
    climdatao["tdew"] = Rcpp::wrap(tdewo);
    climdatao["pk"] = Rcpp::wrap(pko);
    climdatao["swdown"] = Rcpp::wrap(Rswo);
    climdatao["difrad"] = Rcpp::wrap(Rdifo);
    climdatao["lwdown"] = Rcpp::wrap(lwdowno);
    climdatao["windspeed"] = Rcpp::wrap(u2o);
    climdatao["winddir"] = Rcpp::wrap(wdiro);
    // Extract and subset the point model arrays
    NumericVector soilmo = slicea(pointm["soilm"], s);
    NumericVector Tgo = slicea(pointm["Tg"], s);
    NumericVector T0po = slicea(pointm["T0p"], s);
    NumericVector Tbpo = slicea(pointm["Tbp"], s);
    NumericVector Gpo = slicea(pointm["Gp"], s);
    NumericVector DDpo = slicea(pointm["DDp"], s);
    NumericVector umuo = slicea(pointm["umu"], s);
    NumericVector kpo = slicea(pointm["kp"], s);
    NumericVector muGpo = slicea(pointm["muGp"], s);
    NumericVector dtrpo = slicea(pointm["dtrp"], s);
    List pointmo;
    pointmo["soilm"] = Rcpp::wrap(soilmo);
    pointmo["Tg"] = Rcpp::wrap(Tgo);
    pointmo["T0p"] = Rcpp::wrap(T0po);
    pointmo["Tbp"] = Rcpp::wrap(Tbpo);
    pointmo["Gp"] = Rcpp::wrap(Gpo);
    pointmo["DDp"] = Rcpp::wrap(DDpo);
    pointmo["umu"] = Rcpp::wrap(umuo);
    pointmo["kp"] = Rcpp::wrap(kpo);
    pointmo["muGp"] = Rcpp::wrap(muGpo);
    pointmo["dtrp"] = Rcpp::wrap(dtrpo);
    // ** subset the vegetation layers
    NumericMatrix hgto = subsetArray(vegp["hgt"], 0);
    NumericMatrix paio = subsetArray(vegp["pai"], 0);
    NumericMatrix xo = subsetArray(vegp["x"], 0);
    NumericMatrix gsmaxo = subsetArray(vegp["gsmax"], 0);
    NumericMatrix lrefo = subsetArray(vegp["leafr"], 0);
    NumericMatrix ltrao = subsetArray(vegp["leaft"], 0);
    NumericMatrix clumpo = subsetArray(vegp["clump"], 0);
    NumericMatrix leafdo = subsetArray(vegp["leafd"], 0);
    // ** additional
    NumericMatrix paiao = subsetArray(vegp["paia"], 0);
    NumericMatrix leafdeno = subsetArray(vegp["leafden"], 0);
    // wrap into list
    List vegpl;
    vegpl["hgt"] = Rcpp::wrap(hgto);
    vegpl["pai"] = Rcpp::wrap(paio);
    vegpl["x"] = Rcpp::wrap(xo);
    vegpl["gsmax"] = Rcpp::wrap(gsmaxo);
    vegpl["leafr"] = Rcpp::wrap(lrefo);
    vegpl["leaft"] = Rcpp::wrap(ltrao);
    vegpl["clump"] = Rcpp::wrap(clumpo);
    vegpl["leafd"] = Rcpp::wrap(leafdo);
    // ** additional
    vegpl["paia"] = Rcpp::wrap(paiao);
    vegpl["leafden"] = Rcpp::wrap(leafdeno);
    List micro = runmicro2Cpp(obstimeo, climdatao, pointmo, vegpl, soilc, reqhgt, zref, lats, lons, Sminp,
        Smaxp, tfact, complete, mat, out);
    // Create vectors for storing data
    // Create output variables
    NumericVector Tz;
    NumericVector tleaf;
    NumericVector relhum;
    NumericVector soilm;
    NumericVector windspeed;
    NumericVector Rdirdown;
    NumericVector Rdifdown;
    NumericVector Rlwdown;
    NumericVector Rswup;
    NumericVector Rlwup;
    if (out[0]) Tz = micro["Tz"];
    if (out[1]) tleaf = micro["tleaf"];
    if (out[2]) relhum = micro["relhum"];
    if (out[3]) soilm = micro["soilm"];
    if (out[4]) windspeed = micro["windspeed"];
    if (out[5]) Rdirdown = micro["Rdirdown"];
    if (out[6]) Rdifdown = micro["Rdifdown"];
    if (out[7]) Rlwdown = micro["Rlwdown"];
    if (out[8]) Rswup = micro["Rswup"];
    if (out[9]) Rlwup = micro["Rlwup"];
    for (size_t ilyr = 1; ilyr < lyr.size(); ++ilyr) {
        // ** create a subset vector
        int n = ed[ilyr] - st[ilyr] + 1;
        IntegerVector s(n);
        for (int i = 0; i < n; ++i) s[i] = st[ilyr] + i;
        // ** subset the data.frame
        obstimeo = ssdf(obstime, s);
        // Extract and subset climate arrays
        tco = slicea(climdata["tc"], s);
        eso = slicea(climdata["es"], s);
        eao = slicea(climdata["ea"], s);
        tdewo = slicea(climdata["tdew"], s);
        pko = slicea(climdata["pk"], s);
        Rswo = slicea(climdata["swdown"], s);
        Rdifo = slicea(climdata["difrad"], s);
        lwdowno = slicea(climdata["lwdown"], s);
        u2o = slicea(climdata["windspeed"], s);
        wdir = climdata["winddir"]; // not an array
        NumericVector wdiro(n);
        for (int i = 0; i < n; ++i) wdiro[i] = wdir[i + st[ilyr]];
        // wrap into list
        climdatao["tc"] = Rcpp::wrap(tco);
        climdatao["es"] = Rcpp::wrap(eso);
        climdatao["ea"] = Rcpp::wrap(eao);
        climdatao["tdew"] = Rcpp::wrap(tdewo);
        climdatao["pk"] = Rcpp::wrap(pko);
        climdatao["swdown"] = Rcpp::wrap(Rswo);
        climdatao["difrad"] = Rcpp::wrap(Rdifo);
        climdatao["lwdown"] = Rcpp::wrap(lwdowno);
        climdatao["windspeed"] = Rcpp::wrap(u2o);
        climdatao["winddir"] = Rcpp::wrap(wdiro);
        // Extract and subset the point model arrays
        soilmo = slicea(pointm["soilm"], s);
        Tgo = slicea(pointm["Tg"], s);
        T0po = slicea(pointm["T0p"], s);
        Tbpo = slicea(pointm["Tbp"], s);
        Gpo = slicea(pointm["Gp"], s);
        DDpo = slicea(pointm["DDp"], s);
        umuo = slicea(pointm["umu"], s);
        kpo = slicea(pointm["kp"], s);
        muGpo = slicea(pointm["muGp"], s);
        dtrpo = slicea(pointm["dtrp"], s);
        // Convert back into list
        pointmo["soilm"] = Rcpp::wrap(soilmo);
        pointmo["Tg"] = Rcpp::wrap(Tgo);
        pointmo["T0p"] = Rcpp::wrap(T0po);
        pointmo["Tbp"] = Rcpp::wrap(Tbpo);
        pointmo["Gp"] = Rcpp::wrap(Gpo);
        pointmo["DDp"] = Rcpp::wrap(DDpo);
        pointmo["umu"] = Rcpp::wrap(umuo);
        pointmo["kp"] = Rcpp::wrap(kpo);
        pointmo["muGp"] = Rcpp::wrap(muGpo);
        pointmo["dtrp"] = Rcpp::wrap(dtrpo);
        // ** subset the vegetation layers
        hgto = subsetArray(vegp["hgt"], ilyr);
        paio = subsetArray(vegp["pai"], ilyr);
        xo = subsetArray(vegp["x"], ilyr);
        gsmaxo = subsetArray(vegp["gsmax"], ilyr);
        lrefo = subsetArray(vegp["leafr"], ilyr);
        ltrao = subsetArray(vegp["leaft"], ilyr);
        clumpo = subsetArray(vegp["clump"], ilyr);
        leafdo = subsetArray(vegp["leafd"], ilyr);
        // ** additional
        paiao = subsetArray(vegp["paia"], ilyr);
        leafdeno = subsetArray(vegp["leafden"], ilyr);
        // wrap into list
        vegpl["hgt"] = Rcpp::wrap(hgto);
        vegpl["pai"] = Rcpp::wrap(paio);
        vegpl["x"] = Rcpp::wrap(xo);
        vegpl["gsmax"] = Rcpp::wrap(gsmaxo);
        vegpl["leafr"] = Rcpp::wrap(lrefo);
        vegpl["leaft"] = Rcpp::wrap(ltrao);
        vegpl["clump"] = Rcpp::wrap(clumpo);
        vegpl["leafd"] = Rcpp::wrap(leafdo);
        // ** additional
        vegpl["paia"] = Rcpp::wrap(paiao);
        vegpl["leafden"] = Rcpp::wrap(leafdeno);
        micro = runmicro2Cpp(obstimeo, climdatao, pointmo, vegpl, soilc, reqhgt, zref, lats, lons, Sminp,
            Smaxp, tfact, complete, mat, out);
        // abind layers
        if (out[0]) Tz = abind3D(Tz, micro["Tz"]);
        if (out[1]) tleaf = abind3D(tleaf, micro["tleaf"]);
        if (out[2]) relhum = abind3D(relhum, micro["relhum"]);
        if (out[3]) soilm = abind3D(soilm, micro["soilm"]);
        if (out[4]) windspeed = abind3D(windspeed, micro["windspeed"]);
        if (out[5]) Rdirdown = abind3D(Rdirdown, micro["Rdirdown"]);
        if (out[6]) Rdifdown = abind3D(Rdifdown, micro["Rdifdown"]);
        if (out[7]) Rlwdown = abind3D(Rlwdown, micro["Rlwdown"]);
        if (out[8]) Rswup = abind3D(Rswup, micro["Rswup"]);
        if (out[9]) Rlwup = abind3D(Rlwup, micro["Rlwup"]);
    }
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Rcpp::wrap(Tz);
    if (out[1]) outp["tleaf"] = Rcpp::wrap(tleaf);
    if (out[2]) outp["relhum"] = Rcpp::wrap(relhum);
    if (out[3]) outp["soilm"] = Rcpp::wrap(soilm);
    if (out[4]) outp["windspeed"] = Rcpp::wrap(windspeed);
    if (out[5]) outp["Rdirdown"] = Rcpp::wrap(Rdirdown);
    if (out[6]) outp["Rdifdown"] = Rcpp::wrap(Rdifdown);
    if (out[7]) outp["Rlwdown"] = Rcpp::wrap(Rlwdown);
    if (out[8]) outp["Rswup"] = Rcpp::wrap(Rswup);
    if (out[9]) outp["Rlwup"] = Rcpp::wrap(Rlwup);
    return outp;
}
// ********************************************************************** //
// ~~~~~~~~~~~~~~~~~~~~~ bioclim model from here ~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************************************************** //
double calc_std_dev(NumericVector vec) {
    int n = vec.size();
    if (n <= 1) {
        return NA_REAL; // Not enough data to compute standard deviation
    }
    // Calculate the mean
    double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / n;
    // Calculate the sum of squared differences from the mean
    double sum_squared_diff = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_squared_diff += std::pow(vec[i] - mean, 2);
    }
    // Calculate the variance
    double variance = sum_squared_diff / (n - 1);
    // Calculate the standard deviation
    double std_dev = std::sqrt(variance);
    return std_dev;
}
double bioclim1(NumericVector Tz)
{
    // BIO1 = Annual Mean Temperature
    // intepreted as mean of first 12 x 24 entries supplied
    double out = 0.0;
    for (int i = 0; i < 288; ++i) {
        out = out + Tz[i];
    }
    out = out / 288.0;
    return out;
}
double bioclim2(NumericVector Tz)
{
    // BIO2 = Mean Diurnal Range
    // intepreted as diurnal temperature range of first 12 x 24 entries supplied
    NumericVector dtr(12);
    int index = 0;
    // Calculate diurnal temperature range
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
    // Calculate mean diurnal temperature range
    double out = 0.0;
    for (int day = 0; day < 12; day++) out = out + dtr[day];
    out = out / 12;
    return out;
}
double bioclim4(NumericVector Tz)
{
    // BIO4 = Temperature Seasonality (standard deviation  times 100)
    // intepreted as standard deviation of the 12 mean daily temperatures for each month
    NumericVector monmean(12);
    int index = 0;
    // Calculate monthly means
    for (int mth = 0; mth < 12; mth++) {
        monmean[mth] = 0.0;
        for (int hr = 0; hr < 24; hr++) {
            monmean[mth] = monmean[mth] + Tz[index];
            index++;
        }
        monmean[mth] = monmean[mth] / 24;
    }
    double out = calc_std_dev(monmean) * 100.0;
    return out;
}
double bioclim5(NumericVector Tz)
{
    // BIO5 = Max Temperature of Warmest Month
    // intepreted as maximum temperature in the time sequence of 1 x 24 hours
    double tmx = -273.15;
    for (int i = 288; i < 312; i++) {
        if (Tz[i] > tmx) tmx = Tz[i];
    }
    return tmx;
}
double bioclim6(NumericVector Tz)
{
    // BIO6 = Min Temperature of Coldest Month
    // intepreted as minimum temperature in the time sequence of 1 x 24 hours
    double tmn = 273.15;
    for (int i = 312; i < 336; i++) {
        if (Tz[i] < tmn) tmn = Tz[i];
    }
    return tmn;
}
double bioclim8(NumericVector Tz, IntegerVector wetq)
{
    // BIO8 = Mean Temperature of Wettest Quarter
    // intepreted as mean temperature of all entries of wetq (wetq indexed to start on zero)
    double out = 0.0;
    for (R_xlen_t i = 0; i < wetq.size(); i++) {
        out = out + Tz[wetq[i]];
    }
    out = out / 72.0;
    return out;
}
double bioclim9(NumericVector Tz, IntegerVector dryq)
{
    // BIO9 =  Mean Temperature of Driest Quarter
    // intepreted as mean temperature of all entries of dryq (dryq indexed to start on zero)
    double out = 0.0;
    for (R_xlen_t i = 0; i < dryq.size(); i++) {
        out = out + Tz[dryq[i]];
    }
    out = out / 72.0;
    return out;
}
double bioclim10(NumericVector Tz, IntegerVector hotq)
{
    // BIO10 =  Mean Temperature of Warmest Quarter
    // intepreted as mean temperature of all entries of hotq (hotq indexed to start on zero)
    double out = 0.0;
    for (R_xlen_t i = 0; i < hotq.size(); i++) {
        out = out + Tz[hotq[i]];
    }
    out = out / 72.0;
    return out;
}
double bioclim11(NumericVector Tz, IntegerVector colq)
{
    // BIO11 =  Mean Temperature of Coldest Quarter
    // intepreted as mean temperature of all entries of colq (colq indexed to start on zero)
    double out = 0.0;
    for (R_xlen_t i = 0; i < colq.size(); i++) {
        out = out + Tz[colq[i]];
    }
    out = out / 72.0;
    return out;
}
double bioclim12(NumericVector soilm)
{
    // BIO12 =  Annual precipitation
    // intepreted as mean soil moisture over 12 months
    double out = 0.0;
    for (int i = 0; i < 288; i++) {
        out = out + soilm[i];
    }
    out = out / 288.0;
    return out;
}
double bioclim13(NumericVector soilm)
{
    // BIO13 =  Precipitation of Wettest Month
    // intepreted as wettest soil moisture in entire time sequence
    double out = 0.0;
    for (R_xlen_t i = 0; i < soilm.size(); i++) {
        if (soilm[i] > out) out = soilm[i];
    }
    return out;
}
double bioclim14(NumericVector soilm)
{
    // BIO14 =  Precipitation of Driest Month
    // intepreted as driest soil moisture in entire time sequence
    double out = 1.0;
    for (R_xlen_t i = 0; i < soilm.size(); i++) {
        if (soilm[i] < out) out = soilm[i];
    }
    return out;
}
double bioclim15(NumericVector soilm)
{
    // BIO15 =  Precipitation Seasonality (Coefficient of Variation)
    // intepreted as standard deviation / mean of 12 x 24 values
    double me = 0.0;
    for (int i = 0; i < 288; i++) {
        me = me + soilm[i];
    }
    me = me / 288.0;
    double sd = calc_std_dev(soilm);
    double cv = me / sd;
    return cv;
}
double bioclim16(NumericVector soilm, IntegerVector wetq)
{
    // BIO16 =  Precipitation of wettest quarter
    // intepreted as mean soilm in wetq
    double me = 0.0;
    for (R_xlen_t i = 0; i < wetq.size(); i++) {
        me = me + soilm[wetq[i]];
    }
    me = me / 72.0;
    return me;
}
double bioclim17(NumericVector soilm, IntegerVector dryq)
{
    // BIO17 =  Precipitation of driest quarter
    // intepreted as mean soilm in dryq
    double me = 0.0;
    for (R_xlen_t i = 0; i < dryq.size(); i++) {
        me = me + soilm[dryq[i]];
    }
    me = me / 72.0;
    return me;
}
double bioclim18(NumericVector soilm, IntegerVector hotq)
{
    // BIO17 =  Precipitation of warmest quarter
    // intepreted as mean soilm in hotq
    double me = 0.0;
    for (R_xlen_t i = 0; i < hotq.size(); i++) {
        me = me + soilm[hotq[i]];
    }
    me = me / 72.0;
    return me;
}
double bioclim19(NumericVector soilm, IntegerVector colq)
{
    // BIO17 =  Precipitation of coldest quarter
    // intepreted as mean soilm in colq
    double me = 0.0;
    for (R_xlen_t i = 0; i < colq.size(); i++) {
        me = me + soilm[colq[i]];
    }
    me = me / 72.0;
    return me;
}
// Converts NumericVector of hourly data to 14 x 24 matrix
std::vector<std::vector<double>> asdailym(DataFrame df, std::string coln)
{
    NumericVector hr = df[coln];
    // Check that the input has exactly 336 values
    if (hr.size() != 336) {
        Rcpp::stop("Input vector must have exactly 336 elements.");
    }
    // Initialize the 14x24 matrix
    std::vector<std::vector<double>> m(14, std::vector<double>(24));
    // Fill the matrix
    for (size_t i = 0; i < 14; ++i) {
        for (size_t j = 0; j < 24; ++j) {
            m[i][j] = hr[i * 24 + j];
        }
    }
    return m;
}
// Converts IntegerVector of hourly data to 14 x 24 matrix
std::vector<std::vector<int>> asdailymi(DataFrame df, std::string coln)
{
    IntegerVector hr = df[coln];
    // Check that the input has exactly 336 values
    if (hr.size() != 336) {
        Rcpp::stop("Input vector must have exactly 336 elements.");
    }
    // Initialize the 14x24 matrix
    std::vector<std::vector<int>> m(14, std::vector<int>(24));
    // Fill the matrix
    for (size_t i = 0; i < 14; ++i) {
        for (size_t j = 0; j < 24; ++j) {
            m[i][j] = hr[i * 24 + j];
        }
    }
    return m;
}
// Run microclimate model(hourly, static vegetation, data.frame climate input)
// [[Rcpp::export]]
List runbioclim1Cpp(DataFrame obstime, DataFrame climdata, DataFrame pointm, List vegp, List soilc,
    double reqhgt, double zref, double lat, double lon, double Sminp, double Smaxp, double tfact,
    double mat, std::vector<bool> out, IntegerVector wetq, IntegerVector dryq, IntegerVector hotq, 
    IntegerVector colq, bool air)
{
    // Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access columns of climdata
    std::vector<double> tc = climdata["temp"];
    std::vector<double> es = climdata["es"];
    std::vector<double> ea = climdata["ea"];
    std::vector<double> tdew = climdata["tdew"];
    std::vector<double> pk = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> lwdown = climdata["lwdown"];
    std::vector<double> u2 = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    // Access columns of pointm
    std::vector<double> soilmp = pointm["soilm"];
    std::vector<double> Tgp = pointm["Tg"];
    std::vector<double> T0p = pointm["T0p"];
    std::vector<double> Tbp = pointm["Tbp"];
    std::vector<double> Gp = pointm["G"];
    std::vector<double> DDp = pointm["DDp"];
    std::vector<double> umu = pointm["umu"];
    std::vector<double> kp = pointm["kp"];
    std::vector<double> muGp = pointm["muGp"];
    std::vector<double> dtrp = pointm["dtrp"];
    // Access items from vegp
    NumericMatrix hgt = vegp["hgt"];
    NumericMatrix pai = vegp["pai"];
    NumericMatrix x = vegp["x"];
    NumericMatrix gsmax = vegp["gsmax"];
    NumericMatrix lref = vegp["leafr"];
    NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"];
    NumericMatrix leafd = vegp["leafd"];
    // ** additional
    NumericMatrix paia = vegp["paia"];
    NumericMatrix leafden = vegp["leafden"];
    // Access items from soilc
    NumericMatrix Smin = soilc["Smin"];
    NumericMatrix Smax = soilc["Smax"];
    NumericMatrix gref = soilc["gref"];
    NumericMatrix soilb = soilc["soilb"];
    NumericMatrix Psie = soilc["Psie"];
    NumericMatrix Vq = soilc["Vq"];
    NumericMatrix Vm = soilc["Vm"];
    NumericMatrix Mc = soilc["Mc"];
    NumericMatrix rho = soilc["rho"];
    // *** additional
    NumericMatrix slope = soilc["slope"];
    NumericMatrix aspect = soilc["aspect"];
    NumericMatrix twi = soilc["twi"];
    NumericMatrix svfa = soilc["svfa"];
    NumericVector wsa = soilc["wsa"];
    NumericVector hor = soilc["hor"];
    // Get dimensions
    int rows = hgt.nrow();
    int cols = hgt.ncol();
    int tsteps = tc.size();
    // Calculate zenith and humidity variables
    std::vector<double> zend(tc.size());
    std::vector<double> azid(tc.size());
    std::vector<int> sindex(tsteps);
    for (int k = 0; k < tsteps; ++k) {
        std::vector<double> sp = solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
        zend[k] = sp[0];
        azid[k] = sp[1];
        sindex[k] = static_cast<int>(std::round(azid[k] / 15)) % 24;
    }
    // Calculate other odds and sods
    int hiy = 365 * 24;
    if (year[0] % 4 == 0) hiy = 366 * 24;
    double mxtc = -273.15;
    for (int k = 0; k < tsteps; ++k) if (tc[k] > mxtc) mxtc = tc[k];
    // Distribute soil moisture
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    // Compute wind index
    std::vector<int> windex(tsteps);
    for (int k = 0; k < tsteps; ++k) windex[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    // Initialise
    NumericMatrix bio1;
    NumericMatrix bio2;
    NumericMatrix bio3;
    NumericMatrix bio4;
    NumericMatrix bio5;
    NumericMatrix bio6;
    NumericMatrix bio7;
    NumericMatrix bio8;
    NumericMatrix bio9;
    NumericMatrix bio10;
    NumericMatrix bio11;
    NumericMatrix bio12;
    NumericMatrix bio13;
    NumericMatrix bio14;
    NumericMatrix bio15;
    NumericMatrix bio16;
    NumericMatrix bio17;
    NumericMatrix bio18;
    NumericMatrix bio19;
    // Conditionally assign memory to output variables
    if (out[0]) bio1 = NumericMatrix(rows, cols);
    if (out[1]) bio2 = NumericMatrix(rows, cols);
    if (out[2]) bio3 = NumericMatrix(rows, cols);
    if (out[3]) bio4 = NumericMatrix(rows, cols);
    if (out[4]) bio5 = NumericMatrix(rows, cols);
    if (out[5]) bio6 = NumericMatrix(rows, cols);
    if (out[6]) bio7 = NumericMatrix(rows, cols);
    if (out[7]) bio8 = NumericMatrix(rows, cols);
    if (out[8]) bio9 = NumericMatrix(rows, cols);
    if (out[9]) bio10 = NumericMatrix(rows, cols);
    if (out[10]) bio11 = NumericMatrix(rows, cols);
    if (out[11]) bio12 = NumericMatrix(rows, cols);
    if (out[12]) bio13 = NumericMatrix(rows, cols);
    if (out[13]) bio14 = NumericMatrix(rows, cols);
    if (out[14]) bio15 = NumericMatrix(rows, cols);
    if (out[15]) bio16 = NumericMatrix(rows, cols);
    if (out[16]) bio17 = NumericMatrix(rows, cols);
    if (out[17]) bio18 = NumericMatrix(rows, cols);
    if (out[18]) bio19 = NumericMatrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                std::vector<double> si(tsteps);
                std::vector<double> zen(tsteps);
                std::vector<double> ws(tsteps);
                // Calculate solar index etc 
                for (int k = 0; k < tsteps; k++) {
                    si[k] = solarindexCpp(slope(i, j), aspect(i, j), zend[k], azid[k], true);
                    if (si[k] < 0.0) si[k] = 0.0;
                    zen[k] = zend[k] * M_PI / 180;
                    ws[k] = wsa[windex[k] * rows * cols + j * rows + i];
                    double ha = hor[sindex[k] * rows * cols + j * rows + i];
                    double sa = 90 - zend[k];
                    if (ha > tan(sa * M_PI / 180)) si[k] = 0.0;
                }
                microm micro = runmicrovhCpp(reqhgt, zref, si, zen, tc, pk, u2, Rsw, Rdif, lwdown,
                    ea, es, tdew, soilmp, umu, T0p, Tgp, Tbp, Gp, kp, dtrp, muGp,
                    hgt(i, j), pai(i, j), x(i, j), lref(i, j), ltra(i, j), clump(i, j), leafd(i, j),
                    gsmax(i, j), Smin(i, j), Smax(i, j), gref(i, j), soilb(i, j), Psie(i, j), Vq(i, j),
                    Vm(i, j), Mc(i, j), rho(i, j), ws, tadd(i, j), paia(i, j), leafden(i, j),
                    mxtc, svfa(i, j), hiy, false, mat);
                if (out[0]) {
                    if (air) {
                        bio1(i, j) = bioclim1(micro.Tz);
                    }
                    else {
                        bio1(i, j) = bioclim1(micro.tleaf);
                    }
                }
                if (out[1]) {
                    if (air) {
                        bio2(i, j) = bioclim2(micro.Tz);
                    }
                    else {
                        bio2(i, j) = bioclim2(micro.tleaf);
                    }
                }
                if (out[3]) {
                    if (air) {
                        bio4(i, j) = bioclim4(micro.Tz);
                    }
                    else {
                        bio4(i, j) = bioclim4(micro.tleaf);
                    }
                }
                if (out[4]) {
                    if (air) {
                        bio5(i, j) = bioclim5(micro.Tz);
                    }
                    else {
                        bio5(i, j) = bioclim5(micro.tleaf);
                    }
                }
                if (out[5]) {
                    if (air) {
                        bio6(i, j) = bioclim6(micro.Tz);
                    }
                    else {
                        bio6(i, j) = bioclim6(micro.tleaf);
                    }
                }
                if (out[7]) {
                    if (air) {
                        bio8(i, j) = bioclim8(micro.Tz, wetq);
                    }
                    else {
                        bio8(i, j) = bioclim8(micro.tleaf, wetq);
                    }
                }
                if (out[8]) {
                    if (air) {
                        bio9(i, j) = bioclim9(micro.Tz, dryq);
                    }
                    else {
                        bio9(i, j) = bioclim9(micro.tleaf, dryq);
                    }
                }
                if (out[9]) {
                    if (air) {
                        bio10(i, j) = bioclim10(micro.Tz, hotq);
                    }
                    else {
                        bio10(i, j) = bioclim10(micro.tleaf, hotq);
                    }
                }
                if (out[10]) {
                    if (air) {
                        bio11(i, j) = bioclim11(micro.Tz, colq);
                    }
                    else {
                        bio11(i, j) = bioclim11(micro.tleaf, colq);
                    }
                }
                if (out[6]) bio7(i, j) = bio5(i, j) - bio6(i, j);
                if (out[2]) bio3(i, j) = bio2(i, j) / bio7(i, j);
                if (out[11]) bio12(i, j) = bioclim12(micro.soilm);
                if (out[12]) bio13(i, j) = bioclim13(micro.soilm);
                if (out[13]) bio14(i, j) = bioclim14(micro.soilm);
                if (out[14]) bio15(i, j) = bioclim15(micro.soilm);
                if (out[15]) bio16(i, j) = bioclim16(micro.soilm, wetq);
                if (out[16]) bio17(i, j) = bioclim17(micro.soilm, dryq);
                if (out[17]) bio18(i, j) = bioclim18(micro.soilm, hotq);
                if (out[18]) bio19(i, j) = bioclim19(micro.soilm, colq);
            }
            else {
                if (out[0]) bio1(i, j) = NA_REAL;
                if (out[1]) bio2(i, j) = NA_REAL;
                if (out[3]) bio4(i, j) = NA_REAL;
                if (out[4]) bio5(i, j) = NA_REAL;
                if (out[5]) bio6(i, j) = NA_REAL;
                if (out[7]) bio8(i, j) = NA_REAL;
                if (out[8]) bio9(i, j) = NA_REAL;
                if (out[9]) bio10(i, j) = NA_REAL;
                if (out[10]) bio11(i, j) = NA_REAL;
                if (out[6]) bio7(i, j) = NA_REAL;
                if (out[2]) bio3(i, j) = NA_REAL;
                if (out[11]) bio12(i, j) = NA_REAL;
                if (out[12]) bio13(i, j) = NA_REAL;
                if (out[13]) bio14(i, j) = NA_REAL;
                if (out[14]) bio15(i, j) = NA_REAL;
                if (out[15]) bio16(i, j) = NA_REAL;
                if (out[16]) bio17(i, j) = NA_REAL;
                if (out[17]) bio18(i, j) = NA_REAL;
                if (out[18]) bio19(i, j) = NA_REAL;
            } // end else
        } // end col
    } // end row
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["bio1"] = Rcpp::wrap(bio1);
    if (out[1]) outp["bio2"] = Rcpp::wrap(bio2);
    if (out[2]) outp["bio3"] = Rcpp::wrap(bio3);
    if (out[3]) outp["bio4"] = Rcpp::wrap(bio4);
    if (out[4]) outp["bio5"] = Rcpp::wrap(bio5);
    if (out[5]) outp["bio6"] = Rcpp::wrap(bio6);
    if (out[6]) outp["bio7"] = Rcpp::wrap(bio7);
    if (out[7]) outp["bio8"] = Rcpp::wrap(bio8);
    if (out[8]) outp["bio9"] = Rcpp::wrap(bio9);
    if (out[9]) outp["bio10"] = Rcpp::wrap(bio10);
    if (out[10]) outp["bio11"] = Rcpp::wrap(bio11);
    if (out[1]) outp["bio12"] = Rcpp::wrap(bio12);
    if (out[12]) outp["bio13"] = Rcpp::wrap(bio13);
    if (out[13]) outp["bio14"] = Rcpp::wrap(bio14);
    if (out[14]) outp["bio15"] = Rcpp::wrap(bio15);
    if (out[15]) outp["bio16"] = Rcpp::wrap(bio16);
    if (out[16]) outp["bio17"] = Rcpp::wrap(bio17);
    if (out[17]) outp["bio18"] = Rcpp::wrap(bio18);
    if (out[18]) outp["bio19"] = Rcpp::wrap(bio19);
    return outp;
}
// Run microclimate model(hourly, static vegetation, array climate inputs)
// [[Rcpp::export]]
List runbioclim2Cpp(DataFrame obstime, List climdata, List pointm, List vegp, List soilc,
    double reqhgt, double zref, NumericMatrix lats, NumericMatrix lons, double Sminp, double Smaxp, 
    double tfact, double mat, std::vector<bool> out, IntegerVector wetq, IntegerVector dryq, 
    IntegerVector hotq, IntegerVector colq, bool air)
{
    // Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access vectors of climdata
    NumericVector tc = climdata["tc"];
    NumericVector es = climdata["es"];
    NumericVector ea = climdata["ea"];
    NumericVector tdew = climdata["tdew"];
    NumericVector pk = climdata["pk"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector lwdown = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"];
    NumericVector wdir = climdata["winddir"]; // not an array
    // Access vectors of pointm
    NumericVector soilmp = pointm["soilm"];
    NumericVector Tgp = pointm["Tg"];
    NumericVector T0p = pointm["T0p"];
    NumericVector Tbp = pointm["Tbp"];
    NumericVector Gp = pointm["Gp"];
    NumericVector DDp = pointm["DDp"];
    NumericVector umu = pointm["umu"];
    NumericVector kp = pointm["kp"];
    NumericVector muGp = pointm["muGp"];
    NumericVector dtrp = pointm["dtrp"];
    // Access items from vegp
    NumericMatrix hgt = vegp["hgt"];
    NumericMatrix pai = vegp["pai"];
    NumericMatrix x = vegp["x"];
    NumericMatrix gsmax = vegp["gsmax"];
    NumericMatrix lref = vegp["leafr"];
    NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"];
    NumericMatrix leafd = vegp["leafd"];
    // ** additional
    NumericMatrix paia = vegp["paia"];
    NumericMatrix leafden = vegp["leafden"];
    // Access items from soilc
    NumericMatrix Smin = soilc["Smin"];
    NumericMatrix Smax = soilc["Smax"];
    NumericMatrix gref = soilc["gref"];
    NumericMatrix soilb = soilc["soilb"];
    NumericMatrix Psie = soilc["Psie"];
    NumericMatrix Vq = soilc["Vq"];
    NumericMatrix Vm = soilc["Vm"];
    NumericMatrix Mc = soilc["Mc"];
    NumericMatrix rho = soilc["rho"];
    // *** additional
    NumericMatrix slope = soilc["slope"];
    NumericMatrix aspect = soilc["aspect"];
    NumericMatrix twi = soilc["twi"];
    NumericMatrix svfa = soilc["svfa"];
    NumericVector wsa = soilc["wsa"];
    NumericVector hor = soilc["hor"];
    // Get dimensions
    int rows = hgt.nrow();
    int cols = hgt.ncol();
    int tsteps = year.size();
    // Calculate hours in year
    int hiy = 365 * 24;
    if (year[0] % 4 == 0) hiy = 366 * 24;
    // Distribute soil moisture
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    // Initialise
    NumericMatrix bio1;
    NumericMatrix bio2;
    NumericMatrix bio3;
    NumericMatrix bio4;
    NumericMatrix bio5;
    NumericMatrix bio6;
    NumericMatrix bio7;
    NumericMatrix bio8;
    NumericMatrix bio9;
    NumericMatrix bio10;
    NumericMatrix bio11;
    NumericMatrix bio12;
    NumericMatrix bio13;
    NumericMatrix bio14;
    NumericMatrix bio15;
    NumericMatrix bio16;
    NumericMatrix bio17;
    NumericMatrix bio18;
    NumericMatrix bio19;
    // Conditionally assign memory to output variables
    if (out[0]) bio1 = NumericMatrix(rows, cols);
    if (out[1]) bio2 = NumericMatrix(rows, cols);
    if (out[2]) bio3 = NumericMatrix(rows, cols);
    if (out[3]) bio4 = NumericMatrix(rows, cols);
    if (out[4]) bio5 = NumericMatrix(rows, cols);
    if (out[5]) bio6 = NumericMatrix(rows, cols);
    if (out[6]) bio7 = NumericMatrix(rows, cols);
    if (out[7]) bio8 = NumericMatrix(rows, cols);
    if (out[8]) bio9 = NumericMatrix(rows, cols);
    if (out[9]) bio10 = NumericMatrix(rows, cols);
    if (out[10]) bio11 = NumericMatrix(rows, cols);
    if (out[11]) bio12 = NumericMatrix(rows, cols);
    if (out[12]) bio13 = NumericMatrix(rows, cols);
    if (out[13]) bio14 = NumericMatrix(rows, cols);
    if (out[14]) bio15 = NumericMatrix(rows, cols);
    if (out[15]) bio16 = NumericMatrix(rows, cols);
    if (out[16]) bio17 = NumericMatrix(rows, cols);
    if (out[17]) bio18 = NumericMatrix(rows, cols);
    if (out[18]) bio19 = NumericMatrix(rows, cols);
    // repermutate climate variables
    tc = aperm3D2(tc, rows, cols, tsteps);
    es = aperm3D2(es, rows, cols, tsteps);
    ea = aperm3D2(ea, rows, cols, tsteps);
    tdew = aperm3D2(tdew, rows, cols, tsteps);
    pk = aperm3D2(pk, rows, cols, tsteps);
    Rsw = aperm3D2(Rsw, rows, cols, tsteps);
    Rdif = aperm3D2(Rdif, rows, cols, tsteps);
    lwdown = aperm3D2(lwdown, rows, cols, tsteps);
    u2 = aperm3D2(u2, rows, cols, tsteps);
    // repermutate point model variables
    soilmp = aperm3D2(soilmp, rows, cols, tsteps);
    Tgp = aperm3D2(Tgp, rows, cols, tsteps);
    T0p = aperm3D2(T0p, rows, cols, tsteps);
    Tbp = aperm3D2(Tbp, rows, cols, tsteps);
    Gp = aperm3D2(Gp, rows, cols, tsteps);
    DDp = aperm3D2(DDp, rows, cols, tsteps);
    umu = aperm3D2(umu, rows, cols, tsteps);
    kp = aperm3D2(kp, rows, cols, tsteps);
    muGp = aperm3D2(muGp, rows, cols, tsteps);
    dtrp = aperm3D2(dtrp, rows, cols, tsteps);
    // set indices for storing results
    // Compute wind index
    std::vector<int> windex(tsteps);
    for (int k = 0; k < tsteps; ++k) windex[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                // Calculate solar and wind shelter variables
                std::vector<double> zen(tsteps);
                std::vector<double> si(tsteps);
                std::vector<double> ws(tsteps);
                for (int k = 0; k < tsteps; ++k) {
                    ws[k] = wsa[windex[k] * rows * cols + j * rows + i];
                    std::vector<double> sp = solpositionCpp(lats(i, j), lons(i, j), year[k], month[k], day[k], hour[k]);
                    double zend = sp[0];
                    double azid = sp[1];
                    si[k] = solarindexCpp(slope(i, j), aspect(i, j), zend, azid);
                    if (si[k] < 0.0) si[k] = 0.0;
                    zen[k] = zend * M_PI / 180;
                    int sindex = static_cast<int>(std::round(azid / 15)) % 24;
                    double ha = hor[sindex * rows * cols + j * rows + i];
                    double sa = 90 - zend;
                    if (ha > tan(sa * M_PI / 180)) si[k] = 0.0;
                }
                int st = j * tsteps + i * tsteps * cols;
                int ed = st + tsteps - 1;
                // subset climate variables
                std::vector<double> tcv = tosvd(tc[Range(st, ed)]);
                std::vector<double> esv = tosvd(es[Range(st, ed)]);
                std::vector<double> eav = tosvd(ea[Range(st, ed)]);
                std::vector<double> tdewv = tosvd(tdew[Range(st, ed)]);
                std::vector<double> pkv = tosvd(pk[Range(st, ed)]);
                std::vector<double> Rswv = tosvd(Rsw[Range(st, ed)]);
                std::vector<double> Rdifv = tosvd(Rdif[Range(st, ed)]);
                std::vector<double> lwdownv = tosvd(lwdown[Range(st, ed)]);
                std::vector<double> u2v = tosvd(u2[Range(st, ed)]);
                // subset point model variables
                std::vector<double> soilmpv = tosvd(soilmp[Range(st, ed)]);
                std::vector<double> Tgpv = tosvd(Tgp[Range(st, ed)]);
                std::vector<double> T0pv = tosvd(T0p[Range(st, ed)]);
                std::vector<double> Tbpv = tosvd(Tbp[Range(st, ed)]);
                std::vector<double> Gpv = tosvd(Gp[Range(st, ed)]);
                std::vector<double> DDpv = tosvd(DDp[Range(st, ed)]);
                std::vector<double> umuv = tosvd(umu[Range(st, ed)]);
                std::vector<double> kpv = tosvd(kp[Range(st, ed)]);
                std::vector<double> muGpv = tosvd(muGp[Range(st, ed)]);
                std::vector<double> dtrpv = tosvd(dtrp[Range(st, ed)]);
                // calculate other odds and sods
                double mxtc = -273.15;
                for (int k = 0; k < tsteps; ++k) if (tcv[k] > mxtc) mxtc = tcv[k];
                microm micro = runmicrovhCpp(reqhgt, zref, si, zen, tcv, pkv, u2v, Rswv, Rdifv, lwdownv,
                    eav, esv, tdewv, soilmpv, umuv, T0pv, Tgpv, Tbpv, Gpv, kpv, dtrpv, muGpv,
                    hgt(i, j), pai(i, j), x(i, j), lref(i, j), ltra(i, j), clump(i, j), leafd(i, j),
                    gsmax(i, j), Smin(i, j), Smax(i, j), gref(i, j), soilb(i, j), Psie(i, j), Vq(i, j),
                    Vm(i, j), Mc(i, j), rho(i, j), ws, tadd(i, j), paia(i, j), leafden(i, j),
                    mxtc, svfa(i, j), hiy, false, mat);
                if (out[0]) {
                    if (air) {
                        bio1(i, j) = bioclim1(micro.Tz);
                    }
                    else {
                        bio1(i, j) = bioclim1(micro.tleaf);
                    }
                }
                if (out[1]) {
                    if (air) {
                        bio2(i, j) = bioclim2(micro.Tz);
                    }
                    else {
                        bio2(i, j) = bioclim2(micro.tleaf);
                    }
                }
                if (out[3]) {
                    if (air) {
                        bio4(i, j) = bioclim4(micro.Tz);
                    }
                    else {
                        bio4(i, j) = bioclim4(micro.tleaf);
                    }
                }
                if (out[4]) {
                    if (air) {
                        bio5(i, j) = bioclim5(micro.Tz);
                    }
                    else {
                        bio5(i, j) = bioclim5(micro.tleaf);
                    }
                }
                if (out[5]) {
                    if (air) {
                        bio6(i, j) = bioclim6(micro.Tz);
                    }
                    else {
                        bio6(i, j) = bioclim6(micro.tleaf);
                    }
                }
                if (out[7]) {
                    if (air) {
                        bio8(i, j) = bioclim8(micro.Tz, wetq);
                    }
                    else {
                        bio8(i, j) = bioclim8(micro.tleaf, wetq);
                    }
                }
                if (out[8]) {
                    if (air) {
                        bio9(i, j) = bioclim9(micro.Tz, dryq);
                    }
                    else {
                        bio9(i, j) = bioclim9(micro.tleaf, dryq);
                    }
                }
                if (out[9]) {
                    if (air) {
                        bio10(i, j) = bioclim10(micro.Tz, hotq);
                    }
                    else {
                        bio10(i, j) = bioclim10(micro.tleaf, hotq);
                    }
                }
                if (out[10]) {
                    if (air) {
                        bio11(i, j) = bioclim11(micro.Tz, colq);
                    }
                    else {
                        bio11(i, j) = bioclim11(micro.tleaf, colq);
                    }
                }
                if (out[6]) bio7(i, j) = bio5(i, j) - bio6(i, j);
                if (out[2]) bio3(i, j) = bio2(i, j) / bio7(i, j);
                if (out[11]) bio12(i, j) = bioclim12(micro.soilm);
                if (out[12]) bio13(i, j) = bioclim13(micro.soilm);
                if (out[13]) bio14(i, j) = bioclim14(micro.soilm);
                if (out[14]) bio15(i, j) = bioclim15(micro.soilm);
                if (out[15]) bio16(i, j) = bioclim16(micro.soilm, wetq);
                if (out[16]) bio17(i, j) = bioclim17(micro.soilm, dryq);
                if (out[17]) bio18(i, j) = bioclim18(micro.soilm, hotq);
                if (out[18]) bio19(i, j) = bioclim19(micro.soilm, colq);
            }
            else {
                if (out[0]) bio1(i, j) = NA_REAL;
                if (out[1]) bio2(i, j) = NA_REAL;
                if (out[3]) bio4(i, j) = NA_REAL;
                if (out[4]) bio5(i, j) = NA_REAL;
                if (out[5]) bio6(i, j) = NA_REAL;
                if (out[7]) bio8(i, j) = NA_REAL;
                if (out[8]) bio9(i, j) = NA_REAL;
                if (out[9]) bio10(i, j) = NA_REAL;
                if (out[10]) bio11(i, j) = NA_REAL;
                if (out[6]) bio7(i, j) = NA_REAL;
                if (out[2]) bio3(i, j) = NA_REAL;
                if (out[11]) bio12(i, j) = NA_REAL;
                if (out[12]) bio13(i, j) = NA_REAL;
                if (out[13]) bio14(i, j) = NA_REAL;
                if (out[14]) bio15(i, j) = NA_REAL;
                if (out[15]) bio16(i, j) = NA_REAL;
                if (out[16]) bio17(i, j) = NA_REAL;
                if (out[17]) bio18(i, j) = NA_REAL;
                if (out[18]) bio19(i, j) = NA_REAL;
            } // end else
        } // end col
    } // end row
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["bio1"] = Rcpp::wrap(bio1);
    if (out[1]) outp["bio2"] = Rcpp::wrap(bio2);
    if (out[2]) outp["bio3"] = Rcpp::wrap(bio3);
    if (out[3]) outp["bio4"] = Rcpp::wrap(bio4);
    if (out[4]) outp["bio5"] = Rcpp::wrap(bio5);
    if (out[5]) outp["bio6"] = Rcpp::wrap(bio6);
    if (out[6]) outp["bio7"] = Rcpp::wrap(bio7);
    if (out[7]) outp["bio8"] = Rcpp::wrap(bio8);
    if (out[8]) outp["bio9"] = Rcpp::wrap(bio9);
    if (out[9]) outp["bio10"] = Rcpp::wrap(bio10);
    if (out[10]) outp["bio11"] = Rcpp::wrap(bio11);
    if (out[1]) outp["bio12"] = Rcpp::wrap(bio12);
    if (out[12]) outp["bio13"] = Rcpp::wrap(bio13);
    if (out[13]) outp["bio14"] = Rcpp::wrap(bio14);
    if (out[14]) outp["bio15"] = Rcpp::wrap(bio15);
    if (out[15]) outp["bio16"] = Rcpp::wrap(bio16);
    if (out[16]) outp["bio17"] = Rcpp::wrap(bio17);
    if (out[17]) outp["bio18"] = Rcpp::wrap(bio18);
    if (out[18]) outp["bio19"] = Rcpp::wrap(bio19);
    return outp;
}

// Run bioclim model(hourly, changing vegetation, data.frame climate input)
// [[Rcpp::export]]
List runbioclim3Cpp(DataFrame obstime, DataFrame climdata, DataFrame pointm, List vegp, List soilc,
    double reqhgt, double zref, double lat, double lon, double Sminp, double Smaxp, double tfact,
    double mat, std::vector<bool> out, IntegerVector wetq, IntegerVector dryq, IntegerVector hotq,
    IntegerVector colq, bool air)
{
    // Access columns of obstime
    std::vector<std::vector<int>> year = asdailymi(obstime, "year");
    std::vector<std::vector<int>> month = asdailymi(obstime, "month");
    std::vector<std::vector<int>> day = asdailymi(obstime, "day");
    // Access columns of climdata
    std::vector<std::vector<double>> tc = asdailym(climdata, "temp");
    std::vector<std::vector<double>> es = asdailym(climdata, "es");
    std::vector<std::vector<double>> ea = asdailym(climdata, "ea");
    std::vector<std::vector<double>> tdew = asdailym(climdata, "tdew");
    std::vector<std::vector<double>> pk = asdailym(climdata, "pres");
    std::vector<std::vector<double>> Rsw = asdailym(climdata, "swdown");
    std::vector<std::vector<double>> Rdif = asdailym(climdata, "difrad");
    std::vector<std::vector<double>> lwdown = asdailym(climdata, "lwdown");
    std::vector<std::vector<double>> u2 = asdailym(climdata, "windspeed");
    std::vector<std::vector<double>> wdir = asdailym(climdata, "winddir");
    // Access columns of pointm
    std::vector<std::vector<double>> soilmp = asdailym(pointm, "soilm");
    std::vector<std::vector<double>> Tgp = asdailym(pointm, "Tg");
    std::vector<std::vector<double>> T0p = asdailym(pointm, "T0p");
    std::vector<std::vector<double>> Tbp = asdailym(pointm, "Tbp");
    std::vector<std::vector<double>> Gp = asdailym(pointm, "G");
    std::vector<std::vector<double>> DDp = asdailym(pointm, "DDp");
    std::vector<std::vector<double>> umu = asdailym(pointm, "umu");
    std::vector<std::vector<double>> kp = asdailym(pointm, "kp");
    std::vector<std::vector<double>> muGp = asdailym(pointm, "muGp");
    std::vector<std::vector<double>> dtrp = asdailym(pointm, "dtrp");
    // Access items from vegp
    NumericVector hgt = vegp["hgt"];
    NumericVector pai = vegp["pai"];
    NumericVector x = vegp["x"];
    NumericVector gsmax = vegp["gsmax"];
    NumericVector lref = vegp["leafr"];
    NumericVector ltra = vegp["leaft"];
    NumericVector clump = vegp["clump"];
    NumericVector leafd = vegp["leafd"];
    // ** additional
    NumericVector paia = vegp["paia"];
    NumericVector leafden = vegp["leafden"];
    // Access items from soilc
    NumericMatrix Smin = soilc["Smin"];
    NumericMatrix Smax = soilc["Smax"];
    NumericMatrix gref = soilc["gref"];
    NumericMatrix soilb = soilc["soilb"];
    NumericMatrix Psie = soilc["Psie"];
    NumericMatrix Vq = soilc["Vq"];
    NumericMatrix Vm = soilc["Vm"];
    NumericMatrix Mc = soilc["Mc"];
    NumericMatrix rho = soilc["rho"];
    // *** additional
    NumericMatrix slope = soilc["slope"];
    NumericMatrix aspect = soilc["aspect"];
    NumericMatrix twi = soilc["twi"];
    NumericMatrix svfa = soilc["svfa"];
    NumericVector wsa = soilc["wsa"];
    NumericVector hor = soilc["hor"];
    // Get dimensions
    int rows = Smin.nrow();
    int cols = Smin.ncol();
    // Calculate zenith and humidity variables
    std::vector<double> zend(336);
    std::vector<double> azid(336);
    std::vector<int> sindex(336);
    int idx = 0;
    for (int d = 0; d < 14; ++d) {
        for (int h = 0; h < 24; ++h) {
            double hr = static_cast<double>(h);
            std::vector<double> sp = solpositionCpp(lat, lon, year[d][h], month[d][h], day[d][h], hr);
            zend[idx] = sp[0];
            azid[idx] = sp[1];
            sindex[idx] = static_cast<int>(std::round(azid[idx] / 15)) % 24;
            ++idx;
        }
    }
    // Calculate other odds and sods
    int hiy = 365 * 24;
    if (year[0][0] % 4 == 0) hiy = 366 * 24;
    double mxtc = -273.15;
    for (int d = 0; d < 14; ++d) for (int h = 0; h < 24; ++h) if (tc[d][h] > mxtc) mxtc = tc[d][h];
    // Distribute soil moisture
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    // Compute wind index
    idx = 0;
    std::vector<int> windex(336);
    for (int d = 0; d < 14; ++d) {
        for (int h = 0; h < 24; ++h) {
            windex[idx] = static_cast<int>(std::round(wdir[d][h] / 45)) % 8;
            ++idx;
        }
    }
    // Initialise
    NumericMatrix bio1;
    NumericMatrix bio2;
    NumericMatrix bio3;
    NumericMatrix bio4;
    NumericMatrix bio5;
    NumericMatrix bio6;
    NumericMatrix bio7;
    NumericMatrix bio8;
    NumericMatrix bio9;
    NumericMatrix bio10;
    NumericMatrix bio11;
    NumericMatrix bio12;
    NumericMatrix bio13;
    NumericMatrix bio14;
    NumericMatrix bio15;
    NumericMatrix bio16;
    NumericMatrix bio17;
    NumericMatrix bio18;
    NumericMatrix bio19;
    // Conditionally assign memory to output variables
    if (out[0]) bio1 = NumericMatrix(rows, cols);
    if (out[1]) bio2 = NumericMatrix(rows, cols);
    if (out[2]) bio3 = NumericMatrix(rows, cols);
    if (out[3]) bio4 = NumericMatrix(rows, cols);
    if (out[4]) bio5 = NumericMatrix(rows, cols);
    if (out[5]) bio6 = NumericMatrix(rows, cols);
    if (out[6]) bio7 = NumericMatrix(rows, cols);
    if (out[7]) bio8 = NumericMatrix(rows, cols);
    if (out[8]) bio9 = NumericMatrix(rows, cols);
    if (out[9]) bio10 = NumericMatrix(rows, cols);
    if (out[10]) bio11 = NumericMatrix(rows, cols);
    if (out[11]) bio12 = NumericMatrix(rows, cols);
    if (out[12]) bio13 = NumericMatrix(rows, cols);
    if (out[13]) bio14 = NumericMatrix(rows, cols);
    if (out[14]) bio15 = NumericMatrix(rows, cols);
    if (out[15]) bio16 = NumericMatrix(rows, cols);
    if (out[16]) bio17 = NumericMatrix(rows, cols);
    if (out[17]) bio18 = NumericMatrix(rows, cols);
    if (out[18]) bio19 = NumericMatrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = Smin(i, j);
            if (!NumericMatrix::is_na(val)) {
                NumericVector Tz(336);
                NumericVector soilm(336);
                int idx = 0;
                int idx2 = 0;
                for (int d = 0; d < 14; ++d) {
                    std::vector<double> si(24);
                    std::vector<double> zen(24);
                    std::vector<double> ws(24);
                    // Calculate solar index etc 
                    for (int h = 0; h < 24; ++h) {
                        si[h] = solarindexCpp(slope(i, j), aspect(i, j), zend[idx], azid[idx], true);
                        if (si[h] < 0.0) si[h] = 0.0;
                        zen[h] = zend[idx] * M_PI / 180;
                        ws[h] = wsa[windex[idx] * rows * cols + j * rows + i];
                        double ha = hor[sindex[idx] * rows * cols + j * rows + i];
                        double sa = 90 - zend[idx];
                        if (ha > tan(sa * M_PI / 180)) si[h] = 0.0;
                        ++idx;
                    }
                    int idxv = i + rows * (j + cols * d);
                    microm micro = runmicrovhCpp(reqhgt, zref, si, zen, tc[d], pk[d], u2[d], Rsw[d], Rdif[d], 
                        lwdown[d], ea[d], es[d], tdew[d], soilmp[d], umu[d], T0p[d], Tgp[d], Tbp[d], Gp[d], kp[d], dtrp[d], muGp[d],
                        hgt[idxv], pai[idxv], x[idxv], lref[idxv], ltra[idxv], clump[idxv], leafd[idxv],
                        gsmax[idxv], Smin(i, j), Smax(i, j), gref(i, j), soilb(i, j), Psie(i, j), Vq(i, j),
                        Vm(i, j), Mc(i, j), rho(i, j), ws, tadd(i, j), paia[idxv], leafden[idxv],
                        mxtc, svfa(i, j), hiy, false, mat);
                    for (int h = 0; h < 24; ++h) {
                        if (air) {
                            Tz[idx2] = micro.Tz[h];
                        }
                        else {
                            Tz[idx2] = micro.tleaf[h];
                        }
                        soilm[idx2] = micro.soilm[h];
                        ++idx2;
                    }
                }
                if (out[0]) bio1(i, j) = bioclim1(Tz);
                if (out[1]) bio2(i, j) = bioclim2(Tz);
                if (out[3]) bio4(i, j) = bioclim4(Tz);
                if (out[4]) bio5(i, j) = bioclim5(Tz);
                if (out[5]) bio6(i, j) = bioclim6(Tz);
                if (out[7]) bio8(i, j) = bioclim8(Tz, wetq);
                if (out[8]) bio9(i, j) = bioclim9(Tz, dryq);
                if (out[9]) bio10(i, j) = bioclim10(Tz, hotq); 
                if (out[10]) bio11(i, j) = bioclim11(Tz, colq);
                if (out[6]) bio7(i, j) = bio5(i, j) - bio6(i, j);
                if (out[2]) bio3(i, j) = bio2(i, j) / bio7(i, j);
                if (out[11]) bio12(i, j) = bioclim12(soilm);
                if (out[12]) bio13(i, j) = bioclim13(soilm);
                if (out[13]) bio14(i, j) = bioclim14(soilm);
                if (out[14]) bio15(i, j) = bioclim15(soilm);
                if (out[15]) bio16(i, j) = bioclim16(soilm, wetq);
                if (out[16]) bio17(i, j) = bioclim17(soilm, dryq);
                if (out[17]) bio18(i, j) = bioclim18(soilm, hotq);
                if (out[18]) bio19(i, j) = bioclim19(soilm, colq);
            }
            else {
                if (out[0]) bio1(i, j) = NA_REAL;
                if (out[1]) bio2(i, j) = NA_REAL;
                if (out[3]) bio4(i, j) = NA_REAL;
                if (out[4]) bio5(i, j) = NA_REAL;
                if (out[5]) bio6(i, j) = NA_REAL;
                if (out[7]) bio8(i, j) = NA_REAL;
                if (out[8]) bio9(i, j) = NA_REAL;
                if (out[9]) bio10(i, j) = NA_REAL;
                if (out[10]) bio11(i, j) = NA_REAL;
                if (out[6]) bio7(i, j) = NA_REAL;
                if (out[2]) bio3(i, j) = NA_REAL;
                if (out[11]) bio12(i, j) = NA_REAL;
                if (out[12]) bio13(i, j) = NA_REAL;
                if (out[13]) bio14(i, j) = NA_REAL;
                if (out[14]) bio15(i, j) = NA_REAL;
                if (out[15]) bio16(i, j) = NA_REAL;
                if (out[16]) bio17(i, j) = NA_REAL;
                if (out[17]) bio18(i, j) = NA_REAL;
                if (out[18]) bio19(i, j) = NA_REAL;
            } // end else
        } // end col
    } // end row
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["bio1"] = Rcpp::wrap(bio1);
    if (out[1]) outp["bio2"] = Rcpp::wrap(bio2);
    if (out[2]) outp["bio3"] = Rcpp::wrap(bio3);
    if (out[3]) outp["bio4"] = Rcpp::wrap(bio4);
    if (out[4]) outp["bio5"] = Rcpp::wrap(bio5);
    if (out[5]) outp["bio6"] = Rcpp::wrap(bio6);
    if (out[6]) outp["bio7"] = Rcpp::wrap(bio7);
    if (out[7]) outp["bio8"] = Rcpp::wrap(bio8);
    if (out[8]) outp["bio9"] = Rcpp::wrap(bio9);
    if (out[9]) outp["bio10"] = Rcpp::wrap(bio10);
    if (out[10]) outp["bio11"] = Rcpp::wrap(bio11);
    if (out[1]) outp["bio12"] = Rcpp::wrap(bio12);
    if (out[12]) outp["bio13"] = Rcpp::wrap(bio13);
    if (out[13]) outp["bio14"] = Rcpp::wrap(bio14);
    if (out[14]) outp["bio15"] = Rcpp::wrap(bio15);
    if (out[15]) outp["bio16"] = Rcpp::wrap(bio16);
    if (out[16]) outp["bio17"] = Rcpp::wrap(bio17);
    if (out[17]) outp["bio18"] = Rcpp::wrap(bio18);
    if (out[18]) outp["bio19"] = Rcpp::wrap(bio19);
    return outp;
}
// Run bioclim model (hourly, changing vegetation, array climate input)
// [[Rcpp::export]]
List runbioclim4Cpp(DataFrame obstime, List climdata, List pointm, List vegp,
    List soilc, double reqhgt, double zref, NumericMatrix lats, NumericMatrix lons, double Sminp,
    double Smaxp, double tfact, double mat, std::vector<bool> out, IntegerVector wetq, IntegerVector dryq,
    IntegerVector hotq, IntegerVector colq, bool air)
{
    // Access columns of obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access vectors of climdata
    NumericVector tc = climdata["tc"];
    NumericVector es = climdata["es"];
    NumericVector ea = climdata["ea"];
    NumericVector tdew = climdata["tdew"];
    NumericVector pk = climdata["pk"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector lwdown = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"];
    NumericVector wdir = climdata["winddir"]; // not an array
    // Access vectors of pointm
    NumericVector soilmp = pointm["soilm"];
    NumericVector Tgp = pointm["Tg"];
    NumericVector T0p = pointm["T0p"];
    NumericVector Tbp = pointm["Tbp"];
    NumericVector Gp = pointm["Gp"];
    NumericVector DDp = pointm["DDp"];
    NumericVector umu = pointm["umu"];
    NumericVector kp = pointm["kp"];
    NumericVector muGp = pointm["muGp"];
    NumericVector dtrp = pointm["dtrp"];
    // Access items from vegp
    NumericVector hgt = vegp["hgt"];
    NumericVector pai = vegp["pai"];
    NumericVector x = vegp["x"];
    NumericVector gsmax = vegp["gsmax"];
    NumericVector lref = vegp["leafr"];
    NumericVector ltra = vegp["leaft"];
    NumericVector clump = vegp["clump"];
    NumericVector leafd = vegp["leafd"];
    // ** additional
    NumericVector paia = vegp["paia"];
    NumericVector leafden = vegp["leafden"];
    // Access items from soilc
    NumericMatrix Smin = soilc["Smin"];
    NumericMatrix Smax = soilc["Smax"];
    NumericMatrix gref = soilc["gref"];
    NumericMatrix soilb = soilc["soilb"];
    NumericMatrix Psie = soilc["Psie"];
    NumericMatrix Vq = soilc["Vq"];
    NumericMatrix Vm = soilc["Vm"];
    NumericMatrix Mc = soilc["Mc"];
    NumericMatrix rho = soilc["rho"];
    // *** additional
    NumericMatrix slope = soilc["slope"];
    NumericMatrix aspect = soilc["aspect"];
    NumericMatrix twi = soilc["twi"];
    NumericMatrix svfa = soilc["svfa"];
    NumericVector wsa = soilc["wsa"];
    NumericVector hor = soilc["hor"];
    // Distribute soil moisture
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    // Get dimensions
    int rows = Smin.nrow();
    int cols = Smin.ncol();
    // Calculate hours in year
    int hiy = 365 * 24;
    if (year[0] % 4 == 0) hiy = 366 * 24;
    // Initialise
    NumericMatrix bio1;
    NumericMatrix bio2;
    NumericMatrix bio3;
    NumericMatrix bio4;
    NumericMatrix bio5;
    NumericMatrix bio6;
    NumericMatrix bio7;
    NumericMatrix bio8;
    NumericMatrix bio9;
    NumericMatrix bio10;
    NumericMatrix bio11;
    NumericMatrix bio12;
    NumericMatrix bio13;
    NumericMatrix bio14;
    NumericMatrix bio15;
    NumericMatrix bio16;
    NumericMatrix bio17;
    NumericMatrix bio18;
    NumericMatrix bio19;
    // Conditionally assign memory to output variables
    if (out[0]) bio1 = NumericMatrix(rows, cols);
    if (out[1]) bio2 = NumericMatrix(rows, cols);
    if (out[2]) bio3 = NumericMatrix(rows, cols);
    if (out[3]) bio4 = NumericMatrix(rows, cols);
    if (out[4]) bio5 = NumericMatrix(rows, cols);
    if (out[5]) bio6 = NumericMatrix(rows, cols);
    if (out[6]) bio7 = NumericMatrix(rows, cols);
    if (out[7]) bio8 = NumericMatrix(rows, cols);
    if (out[8]) bio9 = NumericMatrix(rows, cols);
    if (out[9]) bio10 = NumericMatrix(rows, cols);
    if (out[10]) bio11 = NumericMatrix(rows, cols);
    if (out[11]) bio12 = NumericMatrix(rows, cols);
    if (out[12]) bio13 = NumericMatrix(rows, cols);
    if (out[13]) bio14 = NumericMatrix(rows, cols);
    if (out[14]) bio15 = NumericMatrix(rows, cols);
    if (out[15]) bio16 = NumericMatrix(rows, cols);
    if (out[16]) bio17 = NumericMatrix(rows, cols);
    if (out[17]) bio18 = NumericMatrix(rows, cols);
    if (out[18]) bio19 = NumericMatrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = Smin(i, j);
            if (!NumericMatrix::is_na(val)) {
                // Calculate mxtc
                double mxtc = -273.15;
                for (int d = 0; d < 14; ++d) {
                    for (int h = 0; h < 24; ++h) {
                        int t = d * 24 + h;
                        int s = i + rows * (j + cols * t);
                        if (tc[s] > mxtc) mxtc = tc[s];
                    }
                }
                NumericVector Tz(336);
                NumericVector soilm(336);
                int idx = 0;
                int idx2 = 0;
                // Loop through days
                for (int d = 0; d < 14; ++d) {
                    // Create variables
                    // ** solar variables
                    std::vector<double> si(24);
                    std::vector<double> zen(24);
                    std::vector<double> ws(24);
                    // ** standard climate variables
                    std::vector<double> tcv(24);
                    std::vector<double> pkv(24);
                    std::vector<double> u2v(24);
                    std::vector<double> Rswv(24);
                    std::vector<double> Rdifv(24);
                    std::vector<double> lwdownv(24);
                    // ** extra climate variables
                    std::vector<double> eav(24);
                    std::vector<double> esv(24);
                    std::vector<double> tdewv(24);
                    // ** point model variables
                    std::vector<double> soilmpv(24);
                    std::vector<double> umuv(24);
                    std::vector<double> Tgpv(24);
                    std::vector<double> T0pv(24);
                    std::vector<double> Tbpv(24);
                    std::vector<double> Gpv(24);
                    std::vector<double> kpv(24);
                    std::vector<double> dtrpv(24);
                    std::vector<double> muGpv(24);
                    for (int h = 0; h < 24; ++h) {
                        // Calculate solar position
                        std::vector<double> sp = solpositionCpp(lats(i, j), lons(i, j), year[idx], month[idx],
                            day[idx], hour[idx]);
                        // Calculate wind index
                        int windex = static_cast<int>(std::round(wdir[idx] / 45)) % 8;
                        ws[h] = wsa[windex * rows * cols + j * rows + i];
                        // Calculate solar index
                        si[h] = solarindexCpp(slope(i, j), aspect(i, j), sp[0], sp[1], true);
                        if (si[h] < 0.0) si[h] = 0.0;
                        zen[h] = sp[0] * M_PI / 180;
                        int sindex = static_cast<int>(std::round(sp[1] / 15)) % 24;
                        double ha = hor[sindex * rows * cols + j * rows + i];
                        if (ha > tan(M_PI / 2 - sp[0])) si[h] = 0.0;
                        // Extract climate and point model data
                        int t = d * 24 + h;
                        int s = i + rows * (j + cols * t);
                        tcv[h] = tc[s];
                        pkv[h] = pk[s];;
                        u2v[h] = u2[s];;
                        Rswv[h] = Rsw[s];;
                        Rdifv[h] = Rdif[s];;
                        lwdownv[h] = lwdown[s];;
                        // ** extra climate variables
                        eav[h] = ea[s];
                        esv[h] = es[s];
                        tdewv[h] = tdew[s];
                        // ** point model variables
                        soilmpv[h] = soilmp[s];
                        umuv[h] = umu[s];
                        Tgpv[h] = Tgp[s];
                        T0pv[h] = T0p[s];
                        Tbpv[h] = Tbp[s];
                        Gpv[h] = Gp[s];
                        kpv[h] = kp[s];
                        dtrpv[h] = dtrp[s];
                        muGpv[h] = muGp[s];
                        ++idx;
                    }
                    int idxv = i + rows * (j + cols * d);
                    microm micro = runmicrovhCpp(reqhgt, zref, si, zen, tcv, pkv, u2v, Rswv, Rdifv,
                        lwdownv, eav, esv, tdewv, soilmpv, umuv, T0pv, Tgpv, Tbpv, Gpv, kpv, dtrpv, muGpv,
                        hgt[idxv], pai[idxv], x[idxv], lref[idxv], ltra[idxv], clump[idxv], leafd[idxv],
                        gsmax[idxv], Smin(i, j), Smax(i, j), gref(i, j), soilb(i, j), Psie(i, j), Vq(i, j),
                        Vm(i, j), Mc(i, j), rho(i, j), ws, tadd(i, j), paia[idxv], leafden[idxv],
                        mxtc, svfa(i, j), hiy, false, mat);
                    for (int h = 0; h < 24; ++h) {
                        if (air) {
                            Tz[idx2] = micro.Tz[h];
                        }
                        else {
                            Tz[idx2] = micro.tleaf[h];
                        }
                        soilm[idx2] = micro.soilm[h];
                        ++idx2;
                    } // end hour
                } // end day
                if (out[0]) bio1(i, j) = bioclim1(Tz);
                if (out[1]) bio2(i, j) = bioclim2(Tz);
                if (out[3]) bio4(i, j) = bioclim4(Tz);
                if (out[4]) bio5(i, j) = bioclim5(Tz);
                if (out[5]) bio6(i, j) = bioclim6(Tz);
                if (out[7]) bio8(i, j) = bioclim8(Tz, wetq);
                if (out[8]) bio9(i, j) = bioclim9(Tz, dryq);
                if (out[9]) bio10(i, j) = bioclim10(Tz, hotq);
                if (out[10]) bio11(i, j) = bioclim11(Tz, colq);
                if (out[6]) bio7(i, j) = bio5(i, j) - bio6(i, j);
                if (out[2]) bio3(i, j) = bio2(i, j) / bio7(i, j);
                if (out[11]) bio12(i, j) = bioclim12(soilm);
                if (out[12]) bio13(i, j) = bioclim13(soilm);
                if (out[13]) bio14(i, j) = bioclim14(soilm);
                if (out[14]) bio15(i, j) = bioclim15(soilm);
                if (out[15]) bio16(i, j) = bioclim16(soilm, wetq);
                if (out[16]) bio17(i, j) = bioclim17(soilm, dryq);
                if (out[17]) bio18(i, j) = bioclim18(soilm, hotq);
                if (out[18]) bio19(i, j) = bioclim19(soilm, colq);
            }
            else {
                if (out[0]) bio1(i, j) = NA_REAL;
                if (out[1]) bio2(i, j) = NA_REAL;
                if (out[3]) bio4(i, j) = NA_REAL;
                if (out[4]) bio5(i, j) = NA_REAL;
                if (out[5]) bio6(i, j) = NA_REAL;
                if (out[7]) bio8(i, j) = NA_REAL;
                if (out[8]) bio9(i, j) = NA_REAL;
                if (out[9]) bio10(i, j) = NA_REAL;
                if (out[10]) bio11(i, j) = NA_REAL;
                if (out[6]) bio7(i, j) = NA_REAL;
                if (out[2]) bio3(i, j) = NA_REAL;
                if (out[11]) bio12(i, j) = NA_REAL;
                if (out[12]) bio13(i, j) = NA_REAL;
                if (out[13]) bio14(i, j) = NA_REAL;
                if (out[14]) bio15(i, j) = NA_REAL;
                if (out[15]) bio16(i, j) = NA_REAL;
                if (out[16]) bio17(i, j) = NA_REAL;
                if (out[17]) bio18(i, j) = NA_REAL;
                if (out[18]) bio19(i, j) = NA_REAL;
            } // end else
        }
    }
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["bio1"] = Rcpp::wrap(bio1);
    if (out[1]) outp["bio2"] = Rcpp::wrap(bio2);
    if (out[2]) outp["bio3"] = Rcpp::wrap(bio3);
    if (out[3]) outp["bio4"] = Rcpp::wrap(bio4);
    if (out[4]) outp["bio5"] = Rcpp::wrap(bio5);
    if (out[5]) outp["bio6"] = Rcpp::wrap(bio6);
    if (out[6]) outp["bio7"] = Rcpp::wrap(bio7);
    if (out[7]) outp["bio8"] = Rcpp::wrap(bio8);
    if (out[8]) outp["bio9"] = Rcpp::wrap(bio9);
    if (out[9]) outp["bio10"] = Rcpp::wrap(bio10);
    if (out[10]) outp["bio11"] = Rcpp::wrap(bio11);
    if (out[1]) outp["bio12"] = Rcpp::wrap(bio12);
    if (out[12]) outp["bio13"] = Rcpp::wrap(bio13);
    if (out[13]) outp["bio14"] = Rcpp::wrap(bio14);
    if (out[14]) outp["bio15"] = Rcpp::wrap(bio15);
    if (out[15]) outp["bio16"] = Rcpp::wrap(bio16);
    if (out[16]) outp["bio17"] = Rcpp::wrap(bio17);
    if (out[17]) outp["bio18"] = Rcpp::wrap(bio18);
    if (out[18]) outp["bio19"] = Rcpp::wrap(bio19);
    return outp;
}
// ============================================================================== #
// ~~~~~~~~~~~~~~~~~~~~~~~~~ Snow model from here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
// ============================================================================== #
// Function to calculate canopy interception of snow
// ** h = canopy height (m)
// ** pai = plant area index (dimensionless)
// ** snowdepth = snow depth (m) 
// ** uf = wind friction velocity (m/s)
// ** prec = hourly snowfall (mm)
// ** tc = air temperature (deg C)
// ** Li = snow load in previous time-step (mm snow water equivelent)
// ** Sh = snow load per unit branch area coefficient (kg/m^3)
double canopysnowintCpp(double hgt, double pai, double uf, double prec,
    double tc, double Li, double Sh = 6.2)
{
    if (hgt < 0.001) hgt = 0.001;
    if (pai < 0.001) pai = 0.001;
    // Calculate mean canopy wind
    double Be = sqrt(0.003 + (0.2 * pai) / 2);
    double uh = uf / Be;
    double a = pai / hgt;
    double Lc = pow(0.25 * a, -1.0);
    double Lm = 2.0 * pow(Be, 3) * Lc;
    double k1 = Be / Lm;
    double uzm = (uh / (hgt * k1)) * (1 - exp(-k1 * hgt));
    if (uzm < uf) uzm = uf;
    // Calculate snow interception
    double rhos = 67.92 + 51.25 * exp(tc / 2.59); // fresh snow density (kg/m^3)
    double S = Sh * (0.26 + 46 / rhos); // maximum snow load per unit branch area (kg/m^3)
    double Lstr = S * pai;  // maximum canopy snow load (mm SWE)
    double Z = atan(uzm / 0.8); // zenith angle of snow fall direction (0.8m/s = terminal velocity of snow flake)
    double kc = 1.0 / (2.0 * cos(Z)); // extinction coefficient(from Campbell assuming spherical distribution)
    double Cp = 1.0 - exp(-kc * pai); // effective canopy cover perpendicular to direction of snow flake
    double k2 = Cp / Lstr; // dimensionless proportionality factor
    double I1 = (Lstr - Li) * (1.0 - exp(-k2 * prec)); // the intercepted snow load at the start of unloading(mm SWE)
    double cis = I1 * 0.678; // Canopy snow interception(mm SWE)
    if (cis > prec) cis = prec;
    return cis;
}
std::vector<double> radoneB(std::vector<double> obstime, std::vector<double> clim, std::vector<double> vegp, std::vector<double> snow,
    std::vector<double> other)
{
    // Extract variables obstime
    int year = static_cast<int>(obstime[0]);
    int month = static_cast<int>(obstime[1]);
    int day = static_cast<int>(obstime[2]);
    double hour = obstime[3];
    // Extract climate
    double Rsw = clim[3];
    double Rdif = clim[4];
    double Rlw = clim[5];
    double Tci = clim[8];
    // Extract variables: vegp
    double pai = vegp[0];
    double hgt = vegp[1];
    double ltra = vegp[2];
    double clump = vegp[3];
    // Extract snow
    double alb = snow[0];
    // Extract other
    double slope = other[0];
    double aspect = other[1];
    double lat = other[2];
    double lon = other[3];
    // ************ Calculate longwave radiation ************* //
    double RlwabsC = 0.97 * Rlw;
    double RlwabsG = RlwabsC;
    double cld = clump * clump;
    double pait = pai / (1 - clump);
    double tr = (1 - cld) * exp(-pait) + cld;
    if (hgt > 0) {
        double Rsky = tr * Rlw;
        double Rcan = (1 - tr) * 0.97 * 5.67 * pow(10.0, -8.0) * pow(Tci + 273.15, 4.0);
        RlwabsG = 0.97 * (Rsky + Rcan);
    }
    // ************ Calculate shortwave radiation ************* //
    double RabsC = RlwabsC;
    double RswabsG = 0.0;
    if (Rsw > 0) {
        // *** Calculate Solar variables ********************* //
        std::vector<double> solp = solpositionCpp(lat, lon, year, month, day, hour);
        double zen = solp[0];
        double azi = solp[1];
        double si = solarindexCpp(slope, aspect, zen, azi);
        if (zen > 90.0) zen = 90.0;
        if (si < 0.0) si = 0.0;
        // *** Calculate radiation absorbed by canopy *** //
        double cosz = cos(zen * M_PI / 180.0);
        double Rbeam = (Rsw - Rdif) / cosz;
        if (Rbeam > 1352.2) Rbeam = 1352.2;
        double RswabsC = (1 - alb) * (Rdif + Rbeam * cosz);
        RabsC = RswabsC + RlwabsC;
        // *** Calculate shortwave radiation absorbed by ground *** //
        RswabsG = RswabsC;
        if (hgt > 0) {
            // base parameters
            double om = alb + ltra;
            if (om > 0.99) om = 0.99;
            double a = 1 - om;
            double del = alb - ltra;
            if (del < 0.01) del = 0.01;
            double J = 1.0 / 3.0;
            double gma = 0.5 * (om + J * del);
            double h = sqrt(a * a + 2 * a * gma);
            // two-stream parameters (diffuse)
            std::vector<double> tspdif = twostreamdifCpp(pait, om, a, gma, h, alb);
            double p3 = tspdif[2];
            double p4 = tspdif[3];
            double u1 = tspdif[4];
            double S1 = tspdif[5];
            double D1 = tspdif[6];
            double D2 = tspdif[7];
            // Calculate canopy extinction coefficient
            std::vector<double> kp = cankCpp(zen, 1, si);
            double kd = kp[1];
            double Kc = kp[2];
            // Calculate two-stream parameters (direct)      
            double sig = kd * kd + gma * gma - pow((a + gma), 2);
            std::vector<double> tspdir = twostreamdirCpp(pait, om, a, gma, J, del, h, alb, kd, sig, u1, S1, D1, D2);
            double p8 = tspdir[3];
            double p9 = tspdir[4];
            double p10 = tspdir[5];
            // Downward diffuse stream
            double clb = pow(clump, Kc);
            double  Rddm = (1 - cld) * (p3 * exp(-h * pait) + p4 * exp(h * pait)) + cld;
            // Contribution of direct to downward diffuse stream
            double Rdbm = (1 - clb) * ((p8 / -sig) * exp(-kd * pait) + p9 * exp(-h * pait) + p10 * exp(h * pait)); 
            // Downward direct stream
            double Rbgm = (1 - clb) * exp(-kd * pait) + clb;
            // Radiation absorbed by ground
            double RdifG = (1 - alb) * (Rdbm * Rbeam * cos(zen) + Rddm * Rdif);
            double RdirG = (1 - alb) * (Rbgm * Rbeam * 0.5);
            RswabsG = RdifG + RdirG;
        } // end if not snow covered
    } // end if daytime
    std::vector<double> out(4);
    out[0] = RabsC;
    out[1] = RswabsG;
    out[2] = RlwabsG;
    out[3] = tr;
    return out;
}
// One point, once cell, bigleaf.
std::vector<double> snowoneB(std::vector<double> obstime, std::vector<double> clim, std::vector<double> vegp, std::vector<double> snow,
    std::vector<double> other, double umu = 1.0)
{
    // Run radiation model
    std::vector<double> rad = radoneB(obstime, clim, vegp, snow, other);
    double RabsC = rad[0];
    double RswabsG = rad[1];
    double RlwabsG = rad[2];
    double RabsG = RswabsG + RlwabsG;
    // Extract climate
    double tc = clim[0];
    double ea = clim[1];
    double pk = clim[2];
    double u2 = clim[6];
    double prec = clim[7];
    double te = clim[8];
    // Extract variables: vegp
    double pai = vegp[0];
    double hgt = vegp[1];
    // Extract snow
    double sdenc = snow[1];
    double sdeng = snow[2];
    double sdepc = snow[3];
    double sdepg = snow[4];
    // Extract other
    double psim = other[4];
    double psih = other[5];
    double G = other[6];
    double zref = other[7];
    // *** Calculate temperature of ground and canopy snow surfaces *** //
    // Calculate convective conductivity
    double d = 0.0;
    double zm = 0.005; 
    if (hgt > 0.0) {
        d = zeroplanedisCpp(hgt, pai);
        zm = roughlengthCpp(hgt, pai, d, psih);
    }
    if (zm < 0.005) zm = 0.005;
    double uf = (0.4 * u2) / (log((zref - d) / zm) + psim);
    uf = uf * umu;
    double ph = phairCpp(tc, pk);
    double gHa = gturbCpp(uf, d, zm, zref, ph, psih, 0.03);
    // Calculate temperatures
    double Tc = PenmanMonteithCpp(RabsC, gHa, gHa, tc, te, pk, ea, 0.97, G, 1.0);
    double Tg = PenmanMonteithCpp(RabsG, gHa, gHa, tc, te, pk, ea, 0.97, G, 1.0);
    double tdew = dewpointCpp(tc, ea);
    if (Tc < tdew) Tc = tdew;
    if (Tg < tdew) Tg = tdew;
    // ******* Calculate mass balance of snowpack (canopy + ground) ******
    // Sublimation
    double la = 45068.7 - 42.8428 * Tc;
    if (Tc < 0) {
        la = 51078.69 - 4.338 * Tc - 0.06367 * Tc * Tc;
    }
    double L = la * (gHa / pk) * (satvapCpp(Tc) - ea);
    la = la / 0.018015; // Conversion to J/kg
    double mSc = (L / la) * 3.6; // m SWE sublimation
    // Melt if snowpack temperature above zero (m SWE)
    double mMc = 0.0;
    double Tcp = Tc;
    if (Tc > 0.0) {
        double S = sdepc * (sdenc / 1000);
        double Fm = 583.3 * Tc * S;
        mMc = (Fm / 334000.0) * 3.6;
        if (sdepc > 0.0) Tc = 0.0;
    }
    // Rain melt (m SWE)
    double mRc = 0.0;
    if (tc > 0.0) {
        mRc = 0.0125 * tc * prec / 1000;
    }
    // ******* Calculate mass balance of snowpack (ground only) ******
    // Sublimation
    la = 45068.7 - 42.8428 * Tg;
    if (Tg < 0) {
        la = 51078.69 - 4.338 * Tg - 0.06367 * Tg * Tg;
    }
    double mu = exp(-pai);
    if (mu > 1.0) mu = 1.0;
    L = la * (gHa / pk) * (satvapCpp(Tg) - ea) * mu;
    la = la / 0.018015; // Conversion to J/kg
    double mSg = (L / la) * 3.6; // m SWE sublimation
    // Melt if snowpack temperature above zero (m SWE)
    double mMg = 0.0;
    if (Tg > 0.0) {
        double S = sdepg * (sdeng / 1000);
        double Fm = 583.3 * Tg * S;
        mMg = (Fm / 334000.0) * 3.6;
        if (sdepg > 0.0) Tg = 0.0;
    }
    // Canopy interception
    double Li = 0.0;
    if (sdepc > 0.0) {
        double wgtg = sdepg / sdepc;
        double sdencc = wgtg * sdeng + (1 - wgtg) * sdenc;
        Li = (sdepc - sdepg) * sdencc;
    }
    if (Li < 0.0) Li = 0.0;
    double cis = canopysnowintCpp(hgt, pai, uf, prec, tc, Li);
    if (cis > prec) cis = prec;
    // Rain melt (m SWE)
    double mRg = 0.0;
    if (tc > 0.0) {
        mRg = 0.0125 * tc * (prec - cis) / 1000;
    }
    // ********************** return outputs ******************
    std::vector<double> out(15);
    // Canopy + ground
    out[0] = Tc;
    out[1] = mSc; // sublimation (m SWE)
    out[2] = mMc; // temperature melt (m SWE)
    out[3] = mRc; // rain melt (m SWE)
    // Ground
    out[4] = Tg;
    out[5] = mSg; // sublimation (m SWE)
    out[6] = mMg; // temperature melt (m SWE)
    out[7] = mRg; // rain melt (m SWE)
    // canopy interception (mm) and uf
    out[8] = cis;
    out[9] = uf;
    out[10] = RswabsG;
    out[11] = RlwabsG;
    out[12] = rad[3];
    out[13] = gHa;
    out[14] = Tcp;
    return out;
}
// **  Function to compute rate of heat storage by snow ** //  
std::vector<double> GFluxCppsnow(std::vector<double> snowt, std::vector<double> snowden)
{
    // Initalise variables that need retaining
    std::vector<double> Gmu(snowt.size());
    std::vector<double> dT(snowt.size());
    std::vector<double> k(snowt.size());
    std::vector<double> ka(snowt.size());
    // Calculate daily mean soil surface temperature
    std::vector<double> Td = hourtodayCpp(snowt, "mean");
    for (size_t i = 0; i < snowt.size(); ++i) {
        // Dive by 1000 to convert snow density in kg/m^3 to g/cm^3
        k[i] = 0.0442 * exp(5.181 * snowden[i] / 1000);
        ka[i] = k[i] / (snowden[i] * 2090);
        double omdy = (2 * M_PI) / (24 * 3600);
        double DD = sqrt(2 * ka[i] / omdy);
        Gmu[i] = sqrt(2) * (k[i] / DD) * 0.5;
        // Calculate T fluctuation from daily mean
        dT[i] = snowt[i] - Td[i];
    }
    // Calculate 6 hour back rolling mean of Gmu and dT to get 3 hour lag
    std::vector<double> Gmud = maCpp(Gmu, 6);
    std::vector<double> G = maCpp(dT, 6);
    for (size_t i = 0; i < G.size(); ++i) G[i] = G[i] * Gmud[i] * 1.1171;
    return G;
}
// ** Function to calculate snow density parameters
std::vector<double> snowdenp(std::string snowenv)
{
    std::vector<double> densfun = { 0.5975, 0.2237, 0.0012, 0.0038 };
    if (snowenv == "Maritime") densfun = { 0.5979, 0.2578, 0.001, 0.0038 };
    else if (snowenv == "Prairie") densfun = { 0.594, 0.2332, 0.0016, 0.0031 };
    else if (snowenv == "Tundra") densfun = { 0.363, 0.2425, 0.0029, 0.0049 };
    else if (snowenv == "Taiga") densfun = { 0.217, 0.217, 0.0, 0.0 };
    return densfun;
}
// Function to calculate snow albedo from age
// ** hs = hours since fresh snow
std::vector<double> snowalbCpp(std::vector<double> prec) {
    std::vector<int> hs(prec.size());
    hs[0] = 0;
    for (size_t i = 1; i < prec.size(); ++i) {
        if (prec[i] > 0) {
            hs[i] = 0;
        }
        else {
            hs[i] = hs[i - 1] + 1;
        }
    }
    std::vector<double> alb(hs.size());  // snow albedo (0-1)
    for (size_t i = 0; i < hs.size(); ++i) {
        alb[i] = (-9.8740 * log(hs[i] / 24) + 78.3434) / 100.0;
        if (alb[i] > 0.95) alb[i] = 0.95;
    }
    return alb;
}
// One point, once cell, bigleaf.
// [[Rcpp::export]]
List pointmodelsnow(DataFrame obstime, DataFrame climdata, std::vector<double> vegp,
    std::vector<double> other, std::string snowenv, double tol = 0.5, double maxiter = 100)
{
    // Access columns of obstime
    std::vector<double> year = obstime["year"];
    std::vector<double> month = obstime["month"];
    std::vector<double> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Access columns of climdata
    std::vector<double> tc = climdata["temp"];
    std::vector<double> rh = climdata["relhum"];
    std::vector<double> pk = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> u2 = climdata["windspeed"];
    std::vector<double> prec = climdata["precip"];
    std::vector<double> ea(tc.size());
    for (size_t i = 0; i < tc.size(); ++i) ea[i] = satvapCpp(tc[i]) * rh[i] / 100;
    std::vector<double> te = tc;
    // Extract other
    double slope = other[0];
    double aspect = other[1];
    double lat = other[2];
    double lon = other[3];
    double zref = other[4];
    double isnowd = other[5]; // initial snow depth
    double isnowa = other[6]; // initial snow age
    // Calculate snow albedo
    std::vector<double> salb = snowalbCpp(prec);
    // Have a first stab at guessing H
    std::vector<double> H(tc.size());
    // Sensible heat flux
    for (size_t i = 0; i < tc.size(); ++i) {
        double Rabs = (1 - salb[i]) * Rsw[i] + 0.97 * Rlw[i];
        H[i] = 0.5 * Rabs;
    }
    // Have a first stab at guessing snow density
    std::vector<double> sdp = snowdenp(snowenv);
    std::vector<double> sdenc(tc.size());
    for (size_t i = 0; i < tc.size(); ++i) {
        sdenc[i] = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * isnowd / 100 -
            sdp[3] * 0)) + sdp[1]) * 1000;
    }
    std::vector<double> sdeng = sdenc;
    // Have a first stab at guessing snow G
    std::vector<double> G = GFluxCppsnow(tc, sdenc);
    // Have a first stab at guessing diabatic coefficients
    std::vector<double> psih(tc.size());
    std::vector<double> psim(tc.size());
    std::vector<double> phih(tc.size(), 1.0);
    // ******************** Initial step ************************* //
    // Initalize variables
    std::vector<double> Tc(tc.size(), tc[1]);
    std::vector<double> Tg(tc.size());
    std::vector<double> sdepc(tc.size() + 1, isnowd);
    std::vector<double> sdepg(tc.size() + 1, 0.5 * isnowd);
    std::vector<double> RswabsG(tc.size());
    std::vector<double> RlwabsG(tc.size());
    std::vector<double> tr(tc.size());
    std::vector<double> umu(tc.size());
    std::vector<double> sublmelt(tc.size());
    std::vector<double> rainmelt(tc.size());
    std::vector<double> tempmelt(tc.size());
    std::vector<double> sstemp(tc.size());
    double tst = 100.0;
    // **** Iterate until convergene
    int iter = 0;
    double mxdif = 0.0;
    while (tst > tol) {
        // **** Iterate through all time steps
        int snowagec = isnowa;
        int snowageg = isnowa;
        // **** Extract antecident Tc and Tg
        std::vector<double> Tco = Tc;
        std::vector<double> Tgo = Tg;
        mxdif = 0.0;
        for (size_t i = 0; i < tc.size(); ++i) {
            // Adjust veg parameters for presence of snow
            double pai = 0.0;
            if (vegp[1] > sdepg[i]) {
                pai = vegp[0] * (vegp[1] - sdepg[i]) / vegp[1];
            }
            double hgt = vegp[1] - sdepg[i];
            double zi = 0.0;
            if (sdepc[i] > 0.0) zi = ((sdepc[i] - sdepg[i]) * sdenc[i]) / (hgt * 1000);
            double ltra = vegp[2] * exp(-10.1 * zi);
            // Get model inputs
            std::vector<double> obstimeo = { year[i],month[i],day[i],hour[i] };
            std::vector<double> climo = { tc[i],ea[i],pk[i],Rsw[i],Rdif[i],Rlw[i],u2[i],prec[i],te[i],Tc[i] };
            std::vector<double> vegpo = { pai,hgt,ltra,vegp[3] };
            std::vector<double> snowo = { salb[i],sdenc[i],sdeng[i],sdepc[i],sdepg[i] };
            std::vector<double> othero = { slope,aspect,lat,lon,psim[i],psih[i],G[i],zref };
            // Run model
            std::vector<double> smod = snowoneB(obstimeo, climo, vegpo, snowo, othero, 1.0);
            // Recalculate change in SWE
            Tc[i] = smod[0];
            Tg[i] = smod[4];
            double snow = prec[i];
            double snowg = prec[i] - smod[8];
            if (tc[i] > 2.0) {
                snow = 0.0;
                snowg = 0.0;
            }
            double swec = snow / 1000.0 - smod[1] - smod[2] - smod[3];
            double sweg = snowg / 1000.0 - smod[5] - smod[6] - smod[7];
            // Recalculate snow density
            snowagec = snowagec + 1;
            snowageg = snowageg + 1;
            sdenc[i] = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * sdepc[i - 1] / 100 -
                sdp[3] * snowagec / 24)) + sdp[1]) * 1000;
            sdeng[i] = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * sdepg[i - 1] / 100 -
                sdp[3] * snowageg / 24)) + sdp[1]) * 1000;
            sdepc[i + 1] = sdepc[i] + (swec * 1000.0) / sdenc[i];
            sdepg[i + 1] = sdepg[i] + (sweg * 1000.0) / sdeng[i];
            if (sdepc[i + 1] < 0.0) {
                sdepc[i + 1] = 0.0;
                snowagec = 0;
            }
            if (sdepg[i + 1] < 0.0) {
                sdepg[i + 1] = 0.0;
                snowageg = 0;
            }
            // Combine the data
            Tc[i] = 0.5 * Tco[i] + 0.5 * Tc[i];
            Tg[i] = 0.5 * Tgo[i] + 0.5 * Tg[i];
            double abs1 = abs(Tc[i] - Tco[i]);
            double abs2 = abs(Tg[i] - Tgo[i]);
            if (mxdif < abs1) mxdif = abs1;
            if (mxdif < abs2) mxdif = abs2;
            // Recalculate H
            double cp = cpairCpp(tc[i]);
            double ph = phairCpp(tc[i], pk[i]);
            double uf = smod[9];
            double gHa = smod[13];
            H[i] = cp * gHa * (Tc[i] - tc[i]);
            // Recalculate diabatic coefficients
            double d = zeroplanedisCpp(hgt, pai);
            double zm = roughlengthCpp(hgt, pai, d, psih[i]);
            double Tk = tc[i] + 273.15;
            double LL = (ph * cp * pow(uf, 3) * Tk) / (-0.4 * 9.81 * H[i]);
            psim[i] = dpsimCpp(zm / LL) - dpsimCpp((zref - d) / LL);
            psih[i] = dpsihCpp((0.2 * zm) / LL) - dpsihCpp((zref - d) / LL);
            phih[i] = dphihCpp((zref - d) / LL);
            // Set limits to diabatic coefficients
            double Belim = 0.4 / sqrt(0.003 + (0.2 * pai) / 2);
            double ln1 = log((zref - d) / zm);
            double ln2 = log((zref - d) / (0.2 * zm));
            if (psim[i] < -0.9 * ln1) psim[i] = -0.9 * ln1;
            if (psih[i] < -0.9 * ln2) psih[i] = -0.9 * ln2;
            if (psim[i] > 0.9 * ln1) psim[i] = 0.9 * ln1;
            if (psih[i] > 0.9 * ln2) psih[i] = 0.9 * ln2;
            if (psih[i] > 0.9 * Belim) psih[i] = 0.9 * Belim;
            RswabsG[i] = smod[10];
            RlwabsG[i] = smod[11];
            tr[i] = smod[12];
            double ufps = (0.4 * u2[i]) / log((zref - d) / zm);
            umu[i] = uf / ufps;
            te[i] = (Tc[i] + tc[i]) / 2.0;
            // Add melt
            sublmelt[i] = smod[1];
            tempmelt[i] = smod[2];
            rainmelt[i] = smod[3];
            sstemp[i] = smod[14];
        }
        // Recalculate G
        G = GFluxCppsnow(Tg, sdenc);
        tst = mxdif;
        ++iter;
        if (iter > maxiter) tst = 0;
    }
    List out;
    out["Tc"] = Rcpp::wrap(Tc);
    out["Tg"] = Rcpp::wrap(Tg);
    out["sdepc"] = Rcpp::wrap(sdepc);
    out["sdepg"] = Rcpp::wrap(sdepg);
    out["sdenc"] = Rcpp::wrap(sdenc);
    out["sdeng"] = Rcpp::wrap(sdeng);
    out["G"] = Rcpp::wrap(G);
    out["RswabsG"] = Rcpp::wrap(RswabsG);
    out["RlwabsG"] = Rcpp::wrap(RlwabsG);
    out["tr"] = Rcpp::wrap(tr);
    out["umu"] = Rcpp::wrap(umu);
    out["mxdif"] = Rcpp::wrap(mxdif);
    out["sublmelt"] = Rcpp::wrap(sublmelt);
    out["tempmelt"] = Rcpp::wrap(tempmelt);
    out["rainmelt"] = Rcpp::wrap(rainmelt);
    out["sstemp"] = Rcpp::wrap(sstemp);
    return out;
}
// snow model - data.frame climate
// [[Rcpp::export]]
List gridmodelsnow1(DataFrame obstime, DataFrame climdata, DataFrame pointm, List vegp,
    List other, std::string snowenv)
{
    // EXtract obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Convert to double
    std::vector<double> yeard = obstime["year"];
    std::vector<double> monthd = obstime["month"];
    std::vector<double> dayd = obstime["day"];
    // Extract pointm
    std::vector<double> Gp = pointm["Gp"];
    std::vector<double> Tcp = pointm["Tc"];
    std::vector<double> RswabsG = pointm["RswabsG"];
    std::vector<double> RlwabsG = pointm["RlwabsG"];
    std::vector<double> umu = pointm["umu"];
    std::vector<double> tr = pointm["tr"];
    // Extract variables: vegp
    NumericMatrix pai = vegp["pai"];
    NumericMatrix hgt = vegp["hgt"];
    NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"];
    // Extract other
    NumericMatrix slope = other["slope"];
    NumericMatrix aspect = other["aspect"];
    NumericMatrix skyview = other["skyview"];
    NumericVector wsa = other["wsa"];
    NumericVector hor = other["hor"];
    double lat = other["lat"];
    double lon = other["lon"];
    double zref = other["zref"];
    NumericMatrix isnowdc = other["isnowdc"]; // initial snow depth (canopy + ground)
    NumericMatrix isnowdg = other["isnowdg"]; // initial snow depth (ground only)
    IntegerMatrix isnowac = other["isnowac"]; // initial snow age (in hours - canopy and ground)
    IntegerMatrix isnowag = other["isnowag"]; // initial snow age (in hours - gorund only)
    // Calculate stuff needed to calculate Gmu
    std::vector<double> tc = climdata["temp"];
    std::vector<double> rh = climdata["relhum"];
    std::vector<double> pk = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> u2 = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> prec = climdata["precip"];
    std::vector<double> salb = snowalbCpp(prec);
    std::vector<double> ea(tc.size());
    std::vector<double> te(tc.size());
    std::vector<double> Rnet(tc.size());
    for (size_t i = 0; i < tc.size(); ++i) {
        ea[i] = satvapCpp(tc[i]) * rh[i] / 100.0;
        te[i] = (Tcp[i] + tc[i]) / 2;
        double Rem = 0.97 * 5.67 * pow(10.0, -8.0) * pow(tc[i] + 273.15, 4);
        Rnet[i] = RswabsG[i] + RlwabsG[i] - Rem;
    }
    int ndays = tr.size() / 24;
    std::vector<double> Rmxd(ndays, -1352.0);
    std::vector<double> Rmnd(ndays, 1352.0);
    std::vector<double> Rswmnd(ndays);
    std::vector<double> Rlwmnd(ndays);
    std::vector<double> Rswmxd(ndays);
    std::vector<double> Rlwmxd(ndays);
    int idx = 0;
    for (int i = 0; i < ndays; ++i) {
        for (int h = 0; h < 24; ++h) {
            if (Rmxd[i] < Rnet[idx]) {
                Rmxd[i] = Rnet[idx];
                Rswmxd[i] = Rsw[idx];
                Rlwmxd[i] = Rlw[idx];
            }
            if (Rmnd[i] > Rnet[idx]) {
                Rmnd[i] = Rnet[idx];
                Rswmnd[i] = Rsw[idx];
                Rlwmnd[i] = Rlw[idx];
            }
            ++idx;
        }
    }
    idx = 0;
    std::vector<double> Rmx(tc.size());
    std::vector<double> Rmn(tc.size());
    std::vector<double> Rswmn(tc.size());
    std::vector<double> Rlwmn(tc.size());
    std::vector<double> Rswmx(tc.size());
    std::vector<double> Rlwmx(tc.size());
    for (int d = 0; d < ndays; ++d) {
        for (int h = 0; h < 24; ++h) {
            Rmx[idx] = Rmxd[d];
            Rmn[idx] = Rmnd[d];
            Rswmn[idx] = Rswmnd[d];
            Rlwmn[idx] = Rlwmnd[d];
            Rswmx[idx] = Rswmxd[d];
            Rlwmx[idx] = Rlwmxd[d];
            ++idx;
        }
    }
    // Calculate Gmx
    idx = 0;
    std::vector<double> Gmxd(ndays,0.0);
    for (int d = 0; d < ndays; ++d) {
        for (int h = 0; h < 24; ++h) {
            if (abs(Rnet[idx]) > Gmxd[d]) Gmxd[d] = abs(Rnet[idx]);
            ++idx;
        }
    }
    idx = 0;
    std::vector<double> Gmx(tc.size());
    for (int d = 0; d < ndays; ++d) {
        for (int h = 0; h < 24; ++h) {
            Gmx[idx] = Gmxd[d];
            ++idx;
        }
    }
    // set dims
    int rows = pai.nrow();
    int cols = pai.ncol();
    int tsteps = tc.size();
    // Calculate solar variables
    std::vector<double> zend(tc.size());
    std::vector<double> azid(tc.size());
    std::vector<int> sindex(tsteps);
    for (int i = 0; i < tsteps; ++i) {
        std::vector<double> sp = solpositionCpp(lat, lon, year[i], month[i], day[i], hour[i]);
        zend[i] = sp[0];
        azid[i] = sp[1];
        sindex[i] = static_cast<int>(std::round(azid[i] / 15)) % 24;
    }
    // Compute wind index
    std::vector<int> windex(tsteps);
    for (int i = 0; i < tsteps; ++i) windex[i] = static_cast<int>(std::round(wdir[i] / 45)) % 8;
    NumericVector Tc(tsteps * rows * cols);
    NumericVector Tg(tsteps * rows * cols);
    NumericVector sdepg(tsteps * rows * cols);
    NumericVector sdepc(tsteps * rows * cols);
    NumericVector sdenca(tsteps * rows * cols);
    NumericMatrix ageg(rows, cols);
    NumericMatrix agec(rows, cols);
    NumericMatrix meltg(rows, cols);
    NumericMatrix meltc(rows, cols);
    /// Snow density in here
    std::vector<double> sdp = snowdenp(snowenv);
    idx = 0;
    // **** Loop through each grid cell and time step
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!Rcpp::NumericMatrix::is_na(val)) {
                int snowagec = isnowac(i, j);
                int snowageg = isnowag(i, j);
                for (int k = 0; k < tsteps; ++k) {
                    // Do snow model test
                    double sdepca = isnowdc(i, j);
                    double sdepga = isnowdg(i, j) * 0.5;
                    if (k > 0) {
                        sdepca = sdepc[idx - 1];
                        sdepga = sdepg[idx - 1];
                    }
                    int snowtest = 0; // whether to run snow model
                    if (sdepca > 0.0) snowtest = 1;
                    if (tc[k] < 2.0 && prec[k] > 0.0) snowtest = 1;
                    if (snowtest > 0) {
                        // Calculate snow density
                        double sdenc = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * sdepca / 100 -
                            sdp[3] * snowagec / 24)) + sdp[1]) * 1000;
                        double sdeng = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * sdepga / 100 -
                            sdp[3] * snowageg / 24)) + sdp[1]) * 1000;
                        // Adjust veg parameters for presence of snow
                        double paip = pai(i, j);
                        if (hgt(i, j) > sdepga) {
                            paip = paip * (hgt(i, j) - sdepga) / hgt(i, j);
                        }
                        double hgtp = hgt(i, j) - sdepga;
                        double zi = 0.0;
                        if (sdepca > 0.0) zi = ((sdepca - sdepga) * sdenc) / (hgtp * 1000);
                        double ltrap = ltra(i, j) * exp(-10.1 * zi);
                        // Calculate Ground heat flux
                        double dtR = Rmx[k] - Rmn[k];
                        double trS = skyview(i, j) * exp(-paip);
                        double Rem = 0.97 * 5.67 * pow(10.0, -8.0) * pow(tc[k] + 273.15, 4);
                        double dmxS = trS * Rswmx[k] + trS * Rlwmx[k] + (1 - trS) * Rem - Rem;
                        double dmnS = trS * Rswmn[k] + trS * Rlwmn[k] + (1 - trS) * Rem - Rem;
                        double Gmu = (dmxS - dmnS) / dtR;
                        double G = Gp[k] * Gmu;
                        if (G > Gmx[k]) G = Gmx[k];
                        if (G < -Gmx[k]) G = -Gmx[k];
                        // Calculate solar multiplier
                        double ha = hor[sindex[k] * rows * cols + j * rows + i];
                        double sa = 90 - zend[k];
                        double smu = 1.0;
                        if (ha > tan(sa * M_PI / 180)) smu = 0.0;
                        // Calculate wind shelter
                        double ws = wsa[windex[k] * rows * cols + j * rows + i];
                        // Get model inputs
                        double u2p = umu[k] * ws;
                        double Rdifp = Rdif[k] * skyview(i, j);
                        double Rdirp = (Rsw[k] - Rdif[k]) * smu;
                        double Rswp = Rdirp + Rdifp;
                        double Rlwp = Rlw[k] * skyview(i, j);
                        // Create model inputs
                        std::vector<double> obstimeo = { yeard[k],monthd[k],dayd[k],hour[k] };
                        std::vector<double> climo = { tc[k],ea[k],pk[k],Rswp,Rdifp,Rlwp,u2p,prec[k],te[k],Tcp[k] };
                        std::vector<double> vegpo = { paip,hgtp,ltrap,clump(i, j) };
                        std::vector<double> snowo = { salb[k],sdenc,sdeng,sdepca,sdepga };
                        std::vector<double> othero = { slope(i,j),aspect(i,j),lat,lon,0.0,0.0,G,zref };
                        // Run snow model for one grid cell
                        std::vector<double> smod = snowoneB(obstimeo, climo, vegpo, snowo, othero, 1.0);
                        // Recalculate change in SWE
                        Tc[idx] = smod[0];
                        Tg[idx] = smod[4];
                        double melc = smod[1] + smod[2] + smod[3];
                        double melg = smod[5] + smod[6] + smod[7];
                        double snoc = prec[k];
                        double snog = prec[k] - smod[8];
                        if (tc[k] > 2.0) {
                            snoc = 0.0;
                            snog = 0.0;
                        }
                        double swec = snoc / 1000.0 - melc;
                        double sweg = snog / 1000.0 - melg;
                        // Recalculate snow density
                        snowagec = snowagec + 1;
                        snowageg = snowageg + 1;
                        sdenc = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * sdepca / 100 -
                            sdp[3] * snowagec / 24)) + sdp[1]) * 1000;
                        sdeng = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * sdepga / 100 -
                            sdp[3] * snowageg / 24)) + sdp[1]) * 1000;
                        sdepc[idx] = sdepca + (swec * 1000.0) / sdenc;
                        sdepg[idx] = sdepga + (sweg * 1000.0) / sdeng;
                        // Calculate cumulative snow melt in m
                        meltc(i, j) = meltc(i, j) + (melc * 1000) / sdenc;
                        meltg(i, j) = meltg(i, j) + (melg * 1000) / sdeng;
                        if (sdepc[idx] < 0.0) {
                            sdepc[idx] = 0.0;
                            snowagec = 0;
                        }
                        if (sdepg[idx] < 0.0) {
                            sdepg[idx] = 0.0;
                            snowageg = 0;
                        }
                        if (k == (tsteps - 1)) {
                            agec(i, j) = snowagec;
                            ageg(i, j) = snowageg;
                        }
                        sdenca[idx] = sdenc;
                    } // end snowtest
                    ++idx;
                } // end k
            } // end NA
            else {
                for (int k = 0; k < tsteps; ++k) {
                    Tc[idx] = NA_REAL;
                    Tg[idx] = NA_REAL;
                    sdepc[idx] = NA_REAL;
                    sdepg[idx] = NA_REAL;
                    sdenca[idx] = NA_REAL;
                    if (k == (tsteps - 1)) {
                        agec(i, j) = NA_REAL;
                        ageg(i, j) = NA_REAL;
                    }
                    ++idx;
                } // and k
            } // end NA test
        } // end col
    } // end row
    // Reshape results
    Tc.attr("dim") = NumericVector::create(tsteps, rows, cols);
    Tg.attr("dim") = NumericVector::create(tsteps, rows, cols);
    sdepc.attr("dim") = NumericVector::create(tsteps, rows, cols);
    sdepg.attr("dim") = NumericVector::create(tsteps, rows, cols);
    sdenca.attr("dim") = NumericVector::create(tsteps, rows, cols);
    Tc = aperm3D(Tc, rows, cols, tsteps);
    Tg = aperm3D(Tg, rows, cols, tsteps);
    sdepc = aperm3D(sdepc, rows, cols, tsteps);
    sdepg = aperm3D(sdepg, rows, cols, tsteps);
    sdenca = aperm3D(sdenca, rows, cols, tsteps);
    // return output
    List out;
    out["Tc"] = Rcpp::wrap(Tc);
    out["Tg"] = Rcpp::wrap(Tg);
    out["sdepc"] = Rcpp::wrap(sdepc);
    out["sdepg"] = Rcpp::wrap(sdepg);
    out["sden"] = Rcpp::wrap(sdenca);
    out["agec"] = Rcpp::wrap(agec);
    out["ageg"] = Rcpp::wrap(ageg);
    out["meltc"] = Rcpp::wrap(meltc);
    out["meltg"] = Rcpp::wrap(meltg);
    return out;
}
// Calculate daily and annual mean snow temperature where n is the number of entries for one year
// [[Rcpp::export]]
List snowdayan(NumericVector stempg)
{
    // Get dimensions of arrays
    IntegerVector dims = stempg.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    // repermutate variables
    stempg = aperm3D2(stempg, rows, cols, tsteps);
    // Create output variables
    NumericVector snowtempd(stempg.size());
    NumericMatrix snowtempa(rows, cols);
    // set indices for storing results
    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = stempg[index];
            if (!NumericVector::is_na(val)) {
                // subset climate variables
                int st = j * tsteps + i * tsteps * cols;
                int ed = st + tsteps - 1;
                // subset to vector
                std::vector<double> stempgv = tosvd(stempg[Range(st, ed)]);
                std::vector<double> stempdv = hourtodayCpp(stempgv, "mean");
                double sumst = 0.0;
                for (int k = 0; k < tsteps; ++k) {
                    snowtempd[index] = stempdv[k];
                    sumst = sumst + stempgv[k];
                    ++index;
                }
                snowtempa(i, j) = sumst / tsteps;
            }
            else {
                for (int k = 0; k < tsteps; ++k) {
                    snowtempd[index] = NA_REAL;
                    ++index;
                }
                snowtempa(i, j) = NA_REAL;
            }
        }
    }
    // Shape results
    snowtempd.attr("dim") = NumericVector::create(tsteps, rows, cols);
    // Reshape results
    snowtempd = aperm3D(snowtempd, rows, cols, tsteps);
    // Create output
    List mout;
    mout["snowtempd"] = Rcpp::wrap(snowtempd);
    mout["snowtempa"] = Rcpp::wrap(snowtempa);
    return mout;
}
NumericMatrix meanDsnow(NumericVector snowden)
{
    IntegerVector dims = snowden.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    // repermutate variables
    snowden = aperm3D2(snowden, rows, cols, tsteps);
    double omdy = (2 * M_PI) / (24 * 3600);
    NumericMatrix meanD(rows, cols);
    int index = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = snowden[index];
            if (!NumericVector::is_na(val)) {
                double sumD = 0.0;
                for (int k = 0; k < tsteps; ++k) {
                    double co = 0.0442 * exp(5.181 * snowden[index] / 1000);
                    double ka = co / (snowden[index] * 2090);
                    sumD = sumD + sqrt(2 * ka / omdy);
                    ++index;
                }
                meanD(i, j) = sumD / tsteps;
            }
            else {
                meanD(i, j) = NA_REAL;
            }
        }
    }
    return meanD;
}
// Calculate Matrix canopy interception fraction
// [[Rcpp::export]]
NumericMatrix canintfrac(NumericMatrix hgt, NumericMatrix pai, double uf,
    double prec, double tc, double Li) {
    int rows = hgt.nrow();
    int cols = hgt.ncol();
    NumericMatrix frac(rows, cols);
    if (prec > 0.0) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double val = hgt(i, j);
                if (!Rcpp::NumericMatrix::is_na(val)) {
                    double ci = canopysnowintCpp(hgt(i, j), pai(i, j), uf, prec, tc, Li);
                    frac(i, j) = ci / prec;
                } 
                else {
                    frac(i, j) = NA_REAL;
                }
            }
        }
    }
    else {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double val = hgt(i, j);
                if (!Rcpp::NumericMatrix::is_na(val)) {
                    frac(i, j) = 0.5;
                }
                else {
                    frac(i, j) = NA_REAL;
                }
            }
        }
    }
    return frac;
}
// Calculate melt mu
// [[Rcpp::export]]
NumericMatrix meltmu(NumericMatrix mu, NumericVector stemp, NumericVector tc)
{
    double dhp = 0.0;
    for (R_xlen_t i = 0; i < stemp.size(); ++i) {
        if (stemp[i] > 0.0) dhp += stemp[i];
    }
    // Calculate multiplier for temperature melt
    int rows = mu.nrow();
    int cols = mu.ncol();
    NumericMatrix mumelt(rows, cols);
    if (dhp > 0.0) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double val = mu(i, j);
                if (!Rcpp::NumericMatrix::is_na(val)) {
                    double dhm = 0.0;
                    for (R_xlen_t k = 0; k < stemp.size(); ++k) {
                        double dif = stemp[k] - tc[k];
                        double stemp2 = dif * mu(i, j) + tc[k];
                        if (stemp2 > 0.0) dhm += stemp2;
                    }
                    mumelt(i, j) = dhm / dhp;
                }
                else {
                    mumelt(i, j) = NA_REAL;
                }
            }
        }
    }
    else {
        // Set all values to 1
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                mumelt(i, j) = 1.0;
            }
        }
    }
    return mumelt;
}
// Calculate melt mu with data provided as arrays
// [[Rcpp::export]]
NumericMatrix meltmu2(NumericMatrix mu, NumericVector stemp, NumericVector tc)
{
    IntegerVector dims = stemp.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    // repermutate variables
    stemp = aperm3D2(stemp, rows, cols, tsteps);
    tc = aperm3D2(tc, rows, cols, tsteps);
    int idx = 0;
    NumericMatrix mumelt(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = mu(i, j);
            if (!Rcpp::NumericMatrix::is_na(val)) {
                double dhp = 0.0;
                double dhm = 0.0;
                for (int k = 0; k < tsteps; ++k) {
                    if (stemp[idx] > 0.0) dhp += stemp[idx];
                    double dif = stemp[idx] - tc[idx];
                    double stemp2 = dif * mu(i, j) + tc[idx];
                    if (stemp2 > 0.0) dhm += stemp2;
                    ++idx;
                }
                if (dhp > 0.0) {
                    mumelt(i, j) = dhm / dhp;
                }
                else {
                    mumelt(i, j) = 0.5;
                }
            }
            else {
                for (int k = 0; k > tsteps; ++k) ++idx;
                mumelt(i, j) = NA_REAL;
            }
        }
    }
    return mumelt;
}
// Calculate snow and no snow days
// [[Rcpp::export]]
List snowdaysfun(NumericVector maxsnowdepth, NumericVector minsnowdepth)
{
    int days = maxsnowdepth.size() / 24;
    IntegerVector snowdays(days);
    IntegerVector nosnowdays(days);
    int index = 0;
    for (int d = 0; d < days; ++d) {
        snowdays[d] = 0;
        nosnowdays[d] = 0;
        for (int h = 0; h < 24; ++h) {
            if (maxsnowdepth[index] > 0.0) snowdays[d] = 1; // 1 indicates snow
            if (minsnowdepth[index] == 0.0) nosnowdays[d] = 1; // 1 indicates no snow
            ++index;
        }
    }
    List out;
    out["snowdays"] = Rcpp::wrap(snowdays);
    out["nosnowdays"] = Rcpp::wrap(nosnowdays);
    return out;
}
// snow microclimate model - data.frame climate
// [[Rcpp::export]]
List gridmicrosnow1(double reqhgt, DataFrame obstime, DataFrame climdata, List snowm, List micro, List vegp, List other,
    std::vector<bool> out) {
    // Extract obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Extract climdata
    std::vector<double> tc = climdata["temp"];
    std::vector<double> rh = climdata["relhum"];
    std::vector<double> pk = climdata["pres"];
    std::vector<double> Rsw = climdata["swdown"];
    std::vector<double> Rdif = climdata["difrad"];
    std::vector<double> Rlw = climdata["lwdown"];
    std::vector<double> u2 = climdata["windspeed"];
    std::vector<double> wdir = climdata["winddir"];
    std::vector<double> prec = climdata["precip"];
    std::vector<double> umu = climdata["umu"];
    // Calculate mxtc
    double mxtc = -273.15;
    for (size_t i = 0; i < tc.size(); ++i) {
        if (tc[i] > mxtc) mxtc = tc[i];
    }
    // Calculate snow albedo
    std::vector<double> salb = snowalbCpp(prec);
    // Extract variables: vegp
    NumericMatrix pai = vegp["pai"];
    NumericMatrix paia = vegp["paia"];
    NumericMatrix hgt = vegp["hgt"];
    NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"];
    NumericMatrix leafd = vegp["leafd"];
    NumericMatrix leafden = vegp["leafden"];
    // Extract other
    NumericMatrix slope = other["slope"];
    NumericMatrix aspect = other["aspect"];
    NumericMatrix skyview = other["skyview"];
    NumericVector wsa = other["wsa"];
    NumericVector hor = other["hor"];
    double lat = other["lat"];
    double lon = other["lon"];
    double zref = other["zref"];
    NumericMatrix Smax = other["Smax"];
    // set dims
    int rows = pai.nrow();
    int cols = pai.ncol();
    int tsteps = tc.size();
    // Extract snow variables
    NumericVector snowtempc = aperm3D2(snowm["Tc"], rows, cols, tsteps);
    NumericVector snowtempg = aperm3D2(snowm["Tg"], rows, cols, tsteps);
    NumericVector swe = aperm3D2(snowm["totalSWE"], rows, cols, tsteps);
    NumericVector sdepg = aperm3D2(snowm["groundsnowdepth"], rows, cols, tsteps);
    NumericVector sden = aperm3D2(snowm["snowden"], rows, cols, tsteps);
    // Calculate snow variables
    List stda = snowdayan(snowm["Tg"]);
    NumericVector Tzd = stda["snowtempd"];
    NumericMatrix Tza = stda["snowtempa"];
    NumericMatrix meanD = meanDsnow(snowm["snowden"]);
    // Extract existing microclimate model variables
    // Initialise
    NumericVector Tz;
    NumericVector tleaf;
    NumericVector relhum;
    NumericVector soilm;
    NumericVector uz;
    NumericVector Rdirdown;
    NumericVector Rdifdown;
    NumericVector Rlwdown;
    NumericVector Rswup;
    NumericVector Rlwup;
    if (out[0]) Tz = aperm3D2(micro["Tz"], rows, cols, tsteps);
    if (out[1]) tleaf = aperm3D2(micro["tleaf"], rows, cols, tsteps);
    if (out[2]) relhum = aperm3D2(micro["relhum"], rows, cols, tsteps);
    if (out[3]) soilm = aperm3D2(micro["soilm"], rows, cols, tsteps);
    if (out[4]) uz = aperm3D2(micro["windspeed"], rows, cols, tsteps);
    if (out[5]) Rdirdown = aperm3D2(micro["Rdirdown"], rows, cols, tsteps);
    if (out[6]) Rdifdown = aperm3D2(micro["Rdifdown"], rows, cols, tsteps);
    if (out[7]) Rlwdown = aperm3D2(micro["Rlwdown"], rows, cols, tsteps);
    if (out[8]) Rswup = aperm3D2(micro["Rswup"], rows, cols, tsteps);
    if (out[9]) Rlwup = aperm3D2(micro["Rlwup"], rows, cols, tsteps);
    // Calculate solar variables
    std::vector<double> zend(tc.size());
    std::vector<double> azid(tc.size());
    std::vector<int> sindex(tsteps);
    for (int i = 0; i < tsteps; ++i) {
        std::vector<double> sp = solpositionCpp(lat, lon, year[i], month[i], day[i], hour[i]);
        zend[i] = sp[0];
        azid[i] = sp[1];
        sindex[i] = static_cast<int>(std::round(azid[i] / 15)) % 24;
    }
    // Compute wind index
    std::vector<int> windex(tsteps);
    for (int i = 0; i < tsteps; ++i) windex[i] = static_cast<int>(std::round(wdir[i] / 45)) % 8;
    // Calculate hours in year
    int hiy = (year[0] % 4 == 0 && (year[0] % 100 != 0 || year[0] % 400 == 0)) ? 366 * 24 : 365 * 24;
    // Go through each grid cell
    int idx = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                for (int k = 0; k < tsteps; ++k) {
                    // check whether to run model
                    if (swe[idx] > 0.0) {
                        // Calculate adjusted reqhgt
                        double reqhgts = reqhgt - sdepg[idx];
                        if (reqhgts >= 0.0) {
                            // Calculate shadowmask
                            double shadowmask = 1.0;
                            double zenr = zend[k] * M_PI / 180.0;
                            double ha = hor[sindex[k] * rows * cols + j * rows + i];
                            double sa = 90 - zend[k];
                            double si = solarindexCpp(slope(i, j), aspect(i, j), zend[k], azid[k], true);
                            // Need to add in sky view
                            if (ha > tan(sa * M_PI / 180)) shadowmask = 0.0;
                            // Calculate wind shelter coefficient
                            double ws = wsa[windex[k] * rows * cols + j * rows + i];
                            // Adjust veg parameters for presence of snow
                            double hgts = hgt(i, j) - sdepg[idx];
                            double pais = 0.0;
                            double paias = 0.0;
                            if (hgts > 0.0) {
                                pais = pai(i, j) * hgts / hgt(i, j);
                                paias = paia(i, j) * hgts / hgt(i, j);
                            }
                            double zi = 0.0;
                            double sdepc = swe[idx] / sden[idx];
                            if (swe[idx] > 0.0) zi = ((sdepc - sdepg[idx]) * sden[idx]) / (hgts * 1000);
                            double ltras = ltra(i, j) * exp(-10.1 * zi);
                            // run microclimate model
                            NumericVector apv = abovepoint(reqhgts, zref, tc[k], pk[k], rh[k], u2[k],
                                Rsw[k], Rdif[k], Rlw[k], hgts, pais, paias, salb[k], ltras,
                                clump(i, j), leafd(i, j), leafden(i, j), salb[k], salb[k], snowtempg[idx],
                                snowtempc[idx], zenr, shadowmask, si, skyview(i, j), ws, umu[k], mxtc);
                            if (out[0]) Tz[idx] = apv[0];
                            if (out[1]) tleaf[idx] = apv[1];
                            if (out[2]) relhum[idx] = apv[2];
                            if (out[4]) uz[idx] = apv[3];
                            if (out[5]) Rdirdown[idx] = apv[4];
                            if (out[6]) Rdifdown[idx] = apv[5];
                            if (out[7]) Rlwdown[idx] = apv[6];
                            if (out[8]) Rswup[idx] = apv[7];
                            if (out[9]) Rlwup[idx] = apv[8];
                        } // end above snow
                        else {
                            double bpv = belowpoint(reqhgts, meanD(i, j), snowtempg[idx], Tzd[idx], Tza(i, j), hiy);
                            if (out[0]) Tz[idx] = bpv;
                            if (out[1]) tleaf[idx] = bpv;
                            if (out[2]) relhum[idx] = 100.0;
                            if (out[4]) uz[idx] = 0.0;
                            if (out[5]) Rdirdown[idx] = 0.0;
                            if (out[6]) Rdifdown[idx] = 0.0;
                            if (out[7]) Rlwdown[idx] = 0.0;
                            if (out[8]) Rswup[idx] = 0.0;
                            if (out[9]) Rlwup[idx] = 0.0;
                        }
                        if (out[3]) soilm[idx] = Smax(i, j);
                    } // end check whether to run model
                    ++idx;
                } // end k
            } // end NA check
            else {
                for (int k = 0; k < tsteps; ++k) ++idx;
            } // end NA
        } // end j
    } // end i
    // Reshape results
    if (out[0]) Tz.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[1]) tleaf.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[2]) relhum.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[3]) soilm.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[4]) uz.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[5]) Rdirdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[6]) Rdifdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[7]) Rlwdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[8]) Rswup.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[9]) Rlwup.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[0]) Tz = aperm3D(Tz, rows, cols, tsteps);
    if (out[1]) tleaf = aperm3D(tleaf, rows, cols, tsteps);
    if (out[2]) relhum = aperm3D(relhum, rows, cols, tsteps);
    if (out[3]) soilm = aperm3D(soilm, rows, cols, tsteps);
    if (out[4]) uz = aperm3D(uz, rows, cols, tsteps);
    if (out[5]) Rdirdown = aperm3D(Rdirdown, rows, cols, tsteps);
    if (out[6]) Rdifdown = aperm3D(Rdifdown, rows, cols, tsteps);
    if (out[7]) Rlwdown = aperm3D(Rlwdown, rows, cols, tsteps);
    if (out[8]) Rswup = aperm3D(Rswup, rows, cols, tsteps);
    if (out[9]) Rlwup = aperm3D(Rlwup, rows, cols, tsteps);
    // return output
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Rcpp::wrap(Tz);
    if (out[1]) outp["tleaf"] = Rcpp::wrap(tleaf);
    if (out[2]) outp["relhum"] = Rcpp::wrap(relhum);
    if (out[3]) outp["soilm"] = Rcpp::wrap(soilm);
    if (out[4]) outp["windspeed"] = Rcpp::wrap(uz);
    if (out[5]) outp["Rdirdown"] = Rcpp::wrap(Rdirdown);
    if (out[6]) outp["Rdifdown"] = Rcpp::wrap(Rdifdown);
    if (out[7]) outp["Rlwdown"] = Rcpp::wrap(Rlwdown);
    if (out[8]) outp["Rswup"] = Rcpp::wrap(Rswup);
    if (out[9]) outp["Rlwup"] = Rcpp::wrap(Rlwup);
    return outp;
}
// snow model - array climate
// [[Rcpp::export]]
List gridmodelsnow2(DataFrame obstime, List climdata, List pointm, List vegp,
    List other, std::string snowenv)
{
    // EXtract obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Convert to double
    std::vector<double> yeard = obstime["year"];
    std::vector<double> monthd = obstime["month"];
    std::vector<double> dayd = obstime["day"];
    // Extract variables: vegp
    NumericMatrix pai = vegp["pai"];
    NumericMatrix hgt = vegp["hgt"];
    NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"];
    // set dims
    int rows = pai.nrow();
    int cols = pai.ncol();
    int tsteps = year.size();
    // Extract pointm
    NumericVector Gp = aperm3D2(pointm["Gp"], rows, cols, tsteps);
    NumericVector Tcp = aperm3D2(pointm["Tc"], rows, cols, tsteps);
    NumericVector RswabsG = aperm3D2(pointm["RswabsG"], rows, cols, tsteps);
    NumericVector RlwabsG = aperm3D2(pointm["RlwabsG"], rows, cols, tsteps);
    NumericVector umu = aperm3D2(pointm["umu"], rows, cols, tsteps);
    NumericVector tr = aperm3D2(pointm["tr"], rows, cols, tsteps);
    // Extract other
    NumericMatrix slope = other["slope"];
    NumericMatrix aspect = other["aspect"];
    NumericMatrix skyview = other["skyview"];
    NumericVector wsa = other["wsa"];
    NumericVector hor = other["hor"];
    NumericMatrix lats = other["lats"];
    NumericMatrix lons = other["lons"];
    double zref = other["zref"];
    NumericMatrix isnowdc = other["isnowdc"]; // initial snow depth (canopy + ground)
    NumericMatrix isnowdg = other["isnowdg"]; // initial snow depth (ground only)
    IntegerMatrix isnowac = other["isnowac"]; // initial snow age (in hours - canopy and ground)
    IntegerMatrix isnowag = other["isnowag"]; // initial snow age (in hours - gorund only)
    // Extract climate data
    NumericVector tc = aperm3D2(climdata["temp"], rows, cols, tsteps);
    NumericVector rh = aperm3D2(climdata["relhum"], rows, cols, tsteps);
    NumericVector pk = aperm3D2(climdata["pres"], rows, cols, tsteps);
    NumericVector Rsw = aperm3D2(climdata["swdown"], rows, cols, tsteps);
    NumericVector Rdif = aperm3D2(climdata["difrad"], rows, cols, tsteps);
    NumericVector Rlw = aperm3D2(climdata["lwdown"], rows, cols, tsteps);
    NumericVector u2 = aperm3D2(climdata["windspeed"], rows, cols, tsteps);
    std::vector<double> wdir = climdata["winddir"];
    NumericVector prec = aperm3D2(climdata["precip"], rows, cols, tsteps);
    // Compute wind index
    std::vector<int> windex(tsteps);
    for (int i = 0; i < tsteps; ++i) windex[i] = static_cast<int>(std::round(wdir[i] / 45)) % 8;
    // Snow density parameters
    std::vector<double> sdp = snowdenp(snowenv);
    // Initialize variables
    NumericVector Tc(tsteps * rows * cols);
    NumericVector Tg(tsteps * rows * cols);
    NumericVector sdepg(tsteps * rows * cols);
    NumericVector sdepc(tsteps * rows * cols);
    NumericVector sdenca(tsteps * rows * cols);
    NumericMatrix ageg(rows, cols);
    NumericMatrix agec(rows, cols);
    NumericMatrix meltg(rows, cols);
    NumericMatrix meltc(rows, cols);
    // **** Loop through each grid cell and time step
    int idx1 = 0;
    int idx2 = 0;
    int idx3 = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!Rcpp::NumericMatrix::is_na(val)) {
                // Calculate base variables
                std::vector<double> precv(tsteps);
                std::vector<double> ea(tsteps);
                std::vector<double> te(tsteps);
                std::vector<double> Rnet(tsteps);
                for (int k = 0; k < tsteps; ++k) {
                    precv[k] = prec[idx1];
                    ea[k] = satvapCpp(tc[idx1]) * rh[idx1] / 100.0;
                    te[k] = (Tcp[idx1] + tc[idx1]) / 2;
                    double Rem = 0.97 * 5.67 * pow(10.0, -8.0) * pow(tc[idx1] + 273.15, 4);
                    Rnet[k] = RswabsG[idx1] + RlwabsG[idx1] - Rem;
                    ++idx1;
                }
                std::vector<double> salb = snowalbCpp(precv);
                // Calculate Gmu variables
                int ndays = tsteps / 24;
                std::vector<double> Rmxd(ndays, -1352.0);
                std::vector<double> Rmnd(ndays, 1352.0);
                std::vector<double> Rswmnd(ndays);
                std::vector<double> Rlwmnd(ndays);
                std::vector<double> Rswmxd(ndays);
                std::vector<double> Rlwmxd(ndays);
                int idx = 0;
                for (int d = 0; d < ndays; ++d) {
                    for (int h = 0; h < 24; ++h) {
                        if (Rmxd[d] < Rnet[idx]) {
                            Rmxd[d] = Rnet[idx];
                            Rswmxd[d] = Rsw[idx2];
                            Rlwmxd[d] = Rlw[idx2];
                        }
                        if (Rmnd[d] > Rnet[idx]) {
                            Rmnd[d] = Rnet[idx];
                            Rswmnd[d] = Rsw[idx2];
                            Rlwmnd[d] = Rlw[idx2];
                        }
                        ++idx;
                        ++idx2;
                    }
                }
                std::vector<double> Rmx(tsteps);
                std::vector<double> Rmn(tsteps);
                std::vector<double> Rswmn(tsteps);
                std::vector<double> Rlwmn(tsteps);
                std::vector<double> Rswmx(tsteps);
                std::vector<double> Rlwmx(tsteps);
                idx = 0;
                for (int d = 0; d < ndays; ++d) {
                    for (int h = 0; h < 24; ++h) {
                        Rmx[idx] = Rmxd[d];
                        Rmn[idx] = Rmnd[d];
                        Rswmn[idx] = Rswmnd[d];
                        Rlwmn[idx] = Rlwmnd[d];
                        Rswmx[idx] = Rswmxd[d];
                        Rlwmx[idx] = Rlwmxd[d];
                        ++idx;
                    }
                }
                // Calculate Gmx
                idx = 0;
                std::vector<double> Gmxd(ndays, 0.0);
                for (int d = 0; d < ndays; ++d) {
                    for (int h = 0; h < 24; ++h) {
                        if (abs(Rnet[idx]) > Gmxd[d]) Gmxd[d] = abs(Rnet[idx]);
                        ++idx;
                    }
                }
                idx = 0;
                std::vector<double> Gmx(tsteps);
                for (int d = 0; d < ndays; ++d) {
                    for (int h = 0; h < 24; ++h) {
                        Gmx[idx] = Gmxd[d];
                        ++idx;
                    }
                }
                // Calculate solar variables
                std::vector<double> zend(tsteps);
                std::vector<double> azid(tsteps);
                std::vector<int> sindex(tsteps);
                for (int k = 0; k < tsteps; ++k) {
                    std::vector<double> sp = solpositionCpp(lats(i, j), lons(i, j), year[k], month[k], day[k], hour[k]);
                    zend[k] = sp[0];
                    azid[k] = sp[1];
                    sindex[k] = static_cast<int>(std::round(azid[k] / 15)) % 24;
                }
                // snow density
                int snowagec = isnowac(i, j);
                int snowageg = isnowag(i, j);
                for (int k = 0; k < tsteps; ++k) {
                    // Do snow model test
                    double sdepca = isnowdc(i, j);
                    double sdepga = isnowdg(i, j);
                    if (k > 0) {
                        sdepca = sdepc[idx3 - 1];
                        sdepga = sdepg[idx3 - 1];
                    }
                    int snowtest = 0; // whether to run snow model
                    if (sdepca > 0.0) snowtest = 1;
                    if (tc[idx3] < 2.0 && prec[idx3] > 0.0) snowtest = 1;
                    if (snowtest > 0) {
                        // Calculate snow density
                        double sdenc = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * sdepca / 100 -
                            sdp[3] * snowagec / 24)) + sdp[1]) * 1000;
                        double sdeng = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * sdepga / 100 -
                            sdp[3] * snowageg / 24)) + sdp[1]) * 1000;
                        // Adjust veg parameters for presence of snow
                        double paip = pai(i, j);
                        if (hgt(i, j) > sdepga) {
                            paip = paip * (hgt(i, j) - sdepga) / hgt(i, j);
                        }
                        double hgtp = hgt(i, j) - sdepga;
                        double zi = 0.0;
                        if (sdepca > 0.0) zi = ((sdepca - sdepga) * sdenc) / (hgtp * 1000);
                        double ltrap = ltra(i, j) * exp(-10.1 * zi);
                        // Calculate Ground heat flux
                        double dtR = Rmx[k] - Rmn[k];
                        double trS = skyview(i, j) * exp(-paip);
                        double Rem = 0.97 * 5.67 * pow(10.0, -8.0) * pow(tc[idx3] + 273.15, 4);
                        double dmxS = trS * Rswmx[k] + trS * Rlwmx[k] + (1 - trS) * Rem - Rem;
                        double dmnS = trS * Rswmn[k] + trS * Rlwmn[k] + (1 - trS) * Rem - Rem;
                        double Gmu = (dmxS - dmnS) / dtR;
                        double G = Gp[idx3] * Gmu;
                        if (G > Gmx[k]) G = Gmx[k];
                        if (G < -Gmx[k]) G = -Gmx[k];
                        // Calculate solar multiplier
                        double ha = hor[sindex[k] * rows * cols + j * rows + i];
                        double sa = 90 - zend[k];
                        double smu = 1.0;
                        if (ha > tan(sa * M_PI / 180)) smu = 0.0;
                        // Calculate wind shelter
                        double ws = wsa[windex[k] * rows * cols + j * rows + i];
                        // Get model inputs
                        double u2p = umu[idx3] * ws;
                        double Rdifp = Rdif[idx3] * skyview(i, j);
                        double Rdirp = (Rsw[idx3] - Rdif[idx3]) * smu;
                        double Rswp = Rdirp + Rdifp;
                        double Rlwp = Rlw[idx3] * skyview(i, j);
                        // Create model inputs
                        std::vector<double> obstimeo = { yeard[k],monthd[k],dayd[k],hour[k] };
                        std::vector<double> climo = { tc[idx3],ea[k],pk[idx3],Rswp,Rdifp,Rlwp,u2p,prec[idx3],te[k],Tcp[idx3] };
                        std::vector<double> vegpo = { paip,hgtp,ltrap,clump(i, j) };
                        std::vector<double> snowo = { salb[k],sdenc,sdeng,sdepca,sdepga };
                        std::vector<double> othero = { slope(i,j),aspect(i,j),lats(i,j),lons(i,j),0.0,0.0,G,zref };
                        // Run snow model for one grid cell
                        std::vector<double> smod = snowoneB(obstimeo, climo, vegpo, snowo, othero, 1.0);
                        // Recalculate change in SWE
                        Tc[idx3] = smod[0];
                        Tg[idx3] = smod[4];
                        double melc = smod[1] + smod[2] + smod[3];
                        double melg = smod[5] + smod[6] + smod[7];
                        double snoc = prec[idx3];
                        double snog = prec[idx3] - smod[8];
                        if (tc[idx3] > 2.0) {
                            snoc = 0.0;
                            snog = 0.0;
                        }
                        double swec = snoc / 1000.0 - melc;
                        double sweg = snog / 1000.0 - melg;
                        // Recalculate snow density
                        snowagec = snowagec + 1;
                        snowageg = snowageg + 1;
                        sdenc = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * sdepca / 100 -
                            sdp[3] * snowagec / 24)) + sdp[1]) * 1000;
                        sdeng = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * sdepga / 100 -
                            sdp[3] * snowageg / 24)) + sdp[1]) * 1000;
                        sdepc[idx3] = sdepca + (swec * 1000.0) / sdenc;
                        sdepg[idx3] = sdepga + (sweg * 1000.0) / sdeng;
                        // Calculate cumulative snow melt in m
                        meltc(i, j) = meltc(i, j) + (melc * 1000) / sdenc;
                        meltg(i, j) = meltg(i, j) + (melg * 1000) / sdeng;
                        if (sdepc[idx3] < 0.0) {
                            sdepc[idx3] = 0.0;
                            snowagec = 0;
                        }
                        if (sdepg[idx3] < 0.0) {
                            sdepg[idx3] = 0.0;
                            snowageg = 0;
                        }
                        if (k == (tsteps - 1)) {
                            agec(i, j) = snowagec;
                            ageg(i, j) = snowageg;
                        }
                        sdenca[idx3] = sdenc;
                    } // end snowtest
                    ++idx3;
                } // end k


            } // end NA
            else {
                for (int k = 0; k < tsteps; ++k) {
                    Tc[idx1] = NA_REAL;
                    Tg[idx1] = NA_REAL;
                    sdepc[idx1] = NA_REAL;
                    sdepg[idx1] = NA_REAL;
                    sdenca[idx1] = NA_REAL;
                    if (k == (tsteps - 1)) {
                        agec(i, j) = NA_REAL;
                        ageg(i, j) = NA_REAL;
                    }
                    ++idx1;
                    ++idx2;
                    ++idx3;
                } // end k
            } // end NA
        } // end j
    } // end k
    // Reshape results
    Tc.attr("dim") = NumericVector::create(tsteps, rows, cols);
    Tg.attr("dim") = NumericVector::create(tsteps, rows, cols);
    sdepc.attr("dim") = NumericVector::create(tsteps, rows, cols);
    sdepg.attr("dim") = NumericVector::create(tsteps, rows, cols);
    sdenca.attr("dim") = NumericVector::create(tsteps, rows, cols);
    Tc = aperm3D(Tc, rows, cols, tsteps);
    Tg = aperm3D(Tg, rows, cols, tsteps);
    sdepc = aperm3D(sdepc, rows, cols, tsteps);
    sdepg = aperm3D(sdepg, rows, cols, tsteps);
    sdenca = aperm3D(sdenca, rows, cols, tsteps);
    // return output
    List out;
    out["Tc"] = Rcpp::wrap(Tc);
    out["Tg"] = Rcpp::wrap(Tg);
    out["sdepc"] = Rcpp::wrap(sdepc);
    out["sdepg"] = Rcpp::wrap(sdepg);
    out["sden"] = Rcpp::wrap(sdenca);
    out["agec"] = Rcpp::wrap(agec);
    out["ageg"] = Rcpp::wrap(ageg);
    out["meltc"] = Rcpp::wrap(meltc);
    out["meltg"] = Rcpp::wrap(meltg);
    return out;
}
// snow microclimate model -array climate
// [[Rcpp::export]]
List gridmicrosnow2(double reqhgt, DataFrame obstime, List climdata, List snowm, List micro, List vegp, List other,
    std::vector<bool> out) {
    // Extract obstime
    std::vector<int> year = obstime["year"];
    std::vector<int> month = obstime["month"];
    std::vector<int> day = obstime["day"];
    std::vector<double> hour = obstime["hour"];
    // Extract variables: vegp
    NumericMatrix pai = vegp["pai"];
    NumericMatrix paia = vegp["paia"];
    NumericMatrix hgt = vegp["hgt"];
    NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"];
    NumericMatrix leafd = vegp["leafd"];
    NumericMatrix leafden = vegp["leafden"];
    // set dims
    int rows = pai.nrow();
    int cols = pai.ncol();
    int tsteps = hour.size();
    // Extract climdata
    NumericVector tc = aperm3D2(climdata["temp"], rows, cols, tsteps);
    NumericVector rh = aperm3D2(climdata["relhum"], rows, cols, tsteps);
    NumericVector pk = aperm3D2(climdata["pres"], rows, cols, tsteps);
    NumericVector Rsw = aperm3D2(climdata["swdown"], rows, cols, tsteps);
    NumericVector Rdif = aperm3D2(climdata["difrad"], rows, cols, tsteps);
    NumericVector Rlw = aperm3D2(climdata["lwdown"], rows, cols, tsteps);
    NumericVector u2 = aperm3D2(climdata["windspeed"], rows, cols, tsteps);
    std::vector<double> wdir = climdata["winddir"];
    NumericVector prec = aperm3D2(climdata["prec"], rows, cols, tsteps);
    NumericVector umu = aperm3D2(climdata["umu"], rows, cols, tsteps);
    // Calculate snow albedo
    NumericVector salb(prec.size());
    int id = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                int st = j * tsteps + i * tsteps * cols;
                int ed = st + tsteps - 1;
                // subset precipitation
                std::vector<double> precv = tosvd(prec[Range(st, ed)]);
                std::vector<double> alb = snowalbCpp(precv);
                for (int k = 0; k < tsteps; ++k) {
                    salb[id] = alb[k];
                    ++id;
                } // end k
            } // end NA check
            else {
                for (int k = 0; k < tsteps; ++k) {
                    salb[id] = NA_REAL;
                    ++id;
                } // end k
            } // end NA check
        } // end j
    } // end i
    // Extract other
    NumericMatrix slope = other["slope"];
    NumericMatrix aspect = other["aspect"];
    NumericMatrix skyview = other["skyview"];
    NumericVector wsa = other["wsa"];
    NumericVector hor = other["hor"];
    NumericMatrix lat = other["lat"];
    NumericMatrix lon = other["lon"];
    NumericMatrix Smax = other["Smax"];
    double zref = other["zref"];
    // Extract snow variables
    NumericVector snowtempc = aperm3D2(snowm["Tc"], rows, cols, tsteps);
    NumericVector snowtempg = aperm3D2(snowm["Tg"], rows, cols, tsteps);
    NumericVector swe = aperm3D2(snowm["totalSWE"], rows, cols, tsteps);
    NumericVector sdepg = aperm3D2(snowm["groundsnowdepth"], rows, cols, tsteps);
    NumericVector sden = aperm3D2(snowm["snowden"], rows, cols, tsteps);
    // Calculate snow variables
    List stda = snowdayan(snowm["Tg"]);
    NumericVector Tzd = stda["snowtempd"];
    NumericMatrix Tza = stda["snowtempa"];
    NumericMatrix meanD = meanDsnow(snowm["snowden"]);
    // Extract existing microclimate model variables
    // Initialise
    NumericVector Tz;
    NumericVector tleaf;
    NumericVector relhum;
    NumericVector soilm;
    NumericVector uz;
    NumericVector Rdirdown;
    NumericVector Rdifdown;
    NumericVector Rlwdown;
    NumericVector Rswup;
    NumericVector Rlwup;
    if (out[0]) Tz = aperm3D2(micro["Tz"], rows, cols, tsteps);
    if (out[1]) tleaf = aperm3D2(micro["tleaf"], rows, cols, tsteps);
    if (out[2]) relhum = aperm3D2(micro["relhum"], rows, cols, tsteps);
    if (out[3]) soilm = aperm3D2(micro["soilm"], rows, cols, tsteps);
    if (out[4]) uz = aperm3D2(micro["windspeed"], rows, cols, tsteps);
    if (out[5]) Rdirdown = aperm3D2(micro["Rdirdown"], rows, cols, tsteps);
    if (out[6]) Rdifdown = aperm3D2(micro["Rdifdown"], rows, cols, tsteps);
    if (out[7]) Rlwdown = aperm3D2(micro["Rlwdown"], rows, cols, tsteps);
    if (out[8]) Rswup = aperm3D2(micro["Rswup"], rows, cols, tsteps);
    if (out[9]) Rlwup = aperm3D2(micro["Rlwup"], rows, cols, tsteps);
    // Compute wind index
    std::vector<int> windex(tsteps);
    for (int i = 0; i < tsteps; ++i) windex[i] = static_cast<int>(std::round(wdir[i] / 45)) % 8;
    // Calculate hours in year
    int hiy = (year[0] % 4 == 0 && (year[0] % 100 != 0 || year[0] % 400 == 0)) ? 366 * 24 : 365 * 24;
    // Go through each grid cell
    int idx0 = 0;
    int idx = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                double mxtc = -273.15;
                for (int k = 0; k < tsteps; ++k) {
                    if (tc[idx0] > mxtc) mxtc = tc[idx0];
                    ++idx0;
                }
                for (int k = 0; k < tsteps; ++k) {
                    // check whether to run model
                    if (swe[idx] > 0.0) {
                        // Calculate adjusted reqhgt
                        double reqhgts = reqhgt - sdepg[idx];
                        if (reqhgts >= 0.0) {
                            // Calculate solar variables
                            std::vector<double> sp = solpositionCpp(lat(i, j), lon(i, j), year[k], month[k], day[k], hour[k]);
                            double zend = sp[0];
                            double azid = sp[1];
                            double sindex = static_cast<int>(std::round(azid / 15)) % 24;
                            // Calculate shadowmask
                            double shadowmask = 1.0;
                            double zenr = zend * M_PI / 180.0;
                            double ha = hor[sindex * rows * cols + j * rows + i];
                            double sa = 90 - zend;
                            double si = solarindexCpp(slope(i, j), aspect(i, j), zend, azid, true);
                            // Need to add in sky view
                            if (ha > tan(sa * M_PI / 180)) shadowmask = 0.0;
                            // Calculate wind shelter coefficient
                            double ws = wsa[windex[k] * rows * cols + j * rows + i];
                            // Adjust veg parameters for presence of snow
                            double hgts = hgt(i, j) - sdepg[idx];
                            double pais = 0.0;
                            double paias = 0.0;
                            if (hgts > 0.0) {
                                pais = pai(i, j) * hgts / hgt(i, j);
                                paias = paia(i, j) * hgts / hgt(i, j);
                            }
                            double sdepc = swe[idx] / sden[idx];
                            double zi = 0.0;
                            if (swe[idx] > 0.0) zi = ((sdepc - sdepg[idx]) * sden[idx]) / (hgts * 1000);
                            double ltras = ltra(i, j) * exp(-10.1 * zi);
                            // run microclimate model
                            NumericVector apv = abovepoint(reqhgts, zref, tc[idx], pk[idx], rh[idx], u2[idx],
                                Rsw[idx], Rdif[idx], Rlw[idx], hgts, pais, paias, salb[idx], ltras,
                                clump(i, j), leafd(i, j), leafden(i, j), salb[idx], salb[idx], snowtempg[idx],
                                snowtempc[idx], zenr, shadowmask, si, skyview(i, j), ws, umu[idx], mxtc);
                            if (out[0]) Tz[idx] = apv[0];
                            if (out[1]) tleaf[idx] = apv[1];
                            if (out[2]) relhum[idx] = apv[2];
                            if (out[4]) uz[idx] = apv[3];
                            if (out[5]) Rdirdown[idx] = apv[4];
                            if (out[6]) Rdifdown[idx] = apv[5];
                            if (out[7]) Rlwdown[idx] = apv[6];
                            if (out[8]) Rswup[idx] = apv[7];
                            if (out[9]) Rlwup[idx] = apv[8];
                        } // end above snow
                        else {
                            double bpv = belowpoint(reqhgts, meanD(i, j), snowtempg[idx], Tzd[idx], Tza(i, j), hiy);
                            if (out[0]) Tz[idx] = bpv;
                            if (out[1]) tleaf[idx] = bpv;
                            if (out[2]) relhum[idx] = 100.0;
                            if (out[4]) uz[idx] = 0.0;
                            if (out[5]) Rdirdown[idx] = 0.0;
                            if (out[6]) Rdifdown[idx] = 0.0;
                            if (out[7]) Rlwdown[idx] = 0.0;
                            if (out[8]) Rswup[idx] = 0.0;
                            if (out[9]) Rlwup[idx] = 0.0;
                        }
                        if (out[3]) soilm[idx] = Smax(i, j);
                    } // end check whether to run model
                    ++idx;
                } // end k
            } // end NA check
            else {
                for (int k = 0; k < tsteps; ++k) ++idx;
            } // end NA
        } // end j
    } // end i
    // Reshape results
    if (out[0]) Tz.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[1]) tleaf.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[2]) relhum.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[3]) soilm.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[4]) uz.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[5]) Rdirdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[6]) Rdifdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[7]) Rlwdown.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[8]) Rswup.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[9]) Rlwup.attr("dim") = NumericVector::create(tsteps, rows, cols);
    if (out[0]) Tz = aperm3D(Tz, rows, cols, tsteps);
    if (out[1]) tleaf = aperm3D(tleaf, rows, cols, tsteps);
    if (out[2]) relhum = aperm3D(relhum, rows, cols, tsteps);
    if (out[3]) soilm = aperm3D(soilm, rows, cols, tsteps);
    if (out[4]) uz = aperm3D(uz, rows, cols, tsteps);
    if (out[5]) Rdirdown = aperm3D(Rdirdown, rows, cols, tsteps);
    if (out[6]) Rdifdown = aperm3D(Rdifdown, rows, cols, tsteps);
    if (out[7]) Rlwdown = aperm3D(Rlwdown, rows, cols, tsteps);
    if (out[8]) Rswup = aperm3D(Rswup, rows, cols, tsteps);
    if (out[9]) Rlwup = aperm3D(Rlwup, rows, cols, tsteps);
    // return output
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Rcpp::wrap(Tz);
    if (out[1]) outp["tleaf"] = Rcpp::wrap(tleaf);
    if (out[2]) outp["relhum"] = Rcpp::wrap(relhum);
    if (out[3]) outp["soilm"] = Rcpp::wrap(soilm);
    if (out[4]) outp["windspeed"] = Rcpp::wrap(uz);
    if (out[5]) outp["Rdirdown"] = Rcpp::wrap(Rdirdown);
    if (out[6]) outp["Rdifdown"] = Rcpp::wrap(Rdifdown);
    if (out[7]) outp["Rlwdown"] = Rcpp::wrap(Rlwdown);
    if (out[8]) outp["Rswup"] = Rcpp::wrap(Rswup);
    if (out[9]) outp["Rlwup"] = Rcpp::wrap(Rlwup);
    return outp;
}
// Function used to calculate leaf reflectance
// [[Rcpp::export]]
double leafrcpp(double om, double pai, double gref, double albin, double x, double ltrr)
{
    // Base parameters
    double lr = om / (ltrr + 1);
    double ltr = ltrr * lr;
    double lref = om - ltr;
    double del = lref - ltr;
    double J = 1.0 / 3.0;
    if (x != 1.0) {
        double mla = 9.65 * pow((3 + x), -1.65);
        if (mla > M_PI / 2) mla = M_PI / 2;
        J = cos(mla) * cos(mla);
    }
    // Two two - stream parameters
    double a = 1 - om;
    double gma = 0.5 * (om + J * del);
    // Intermediate parameters
    double h = sqrt(a * a + 2 * a * gma);
    double S1 = exp(-h * pai);
    double u1 = a + gma * (1 - 1 / gref);
    double D1 = (a + gma + h) * (u1 - h) * 1 / S1 - (a + gma - h) * (u1 + h) * S1;
    // p parameters
    double p1 = gma / (D1 * S1) * (u1 - h);
    double p2 = (-gma * S1 / D1) * (u1 + h);
    //  albedo
    double albw = p1 + p2;
    double out = albw - albin;
    return out;
}