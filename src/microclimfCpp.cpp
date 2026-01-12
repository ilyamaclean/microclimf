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
#include <numeric>
using namespace Rcpp;
constexpr double pi = 3.14159265358979323846;
constexpr double torad = 3.14159265358979323846 / 180.0;
constexpr double sb = 5.67e-8;
constexpr double thetam = 0.365;
constexpr double ka = 0.4;
constexpr double omdy = (2.0 * 3.14159265358979323846) / (24.0 * 3600.0);
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************** Worker functions ****************************************************************** //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ** Derive emitted radiation
double radem(double tc) {
    return std::pow(tc + 273.15, 4.0);
}
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

    double m = 6.24004077 + 0.01720197 * (jd - 2451545.0);
    double eot = -7.659 * std::sin(m) + 9.863 * std::sin(2 * m + 3.5932);
    double st = lt + (4.0 * lon + eot) / 60.0;
    return st;
}
// ** Calculates solar position ** //
solmodel solpositionCpp(double lat, double lon, int year, int month, int day, double lt)
{
    int jd = juldayCpp(year, month, day);
    double st = soltimeCpp(jd, lt, lon);
    // Calculate solar zenith (degrees)
    double latr = lat * pi / 180.0;
    double tt = 0.261799 * (st - 12);
    double dec = (pi * 23.5 / 180) * std::cos(2 * pi * ((jd - 159.5) / 365.25));
    double coh = std::sin(dec) * std::sin(latr) + std::cos(dec) * std::cos(latr) * std::cos(tt);
    double z = std::acos(coh) * (180 / pi);
    // Calculate solar azimuth (degrees)
    double sh = std::sin(dec) * std::sin(latr) + std::cos(dec) * std::cos(latr) * std::cos(tt);
    double hh = std::atan(sh / std::sqrt(1 - sh * sh));
    double sazi = std::cos(dec) * std::sin(tt) / std::cos(hh);
    double cazi = (std::sin(latr) * std::cos(dec) * std::cos(tt) - std::cos(latr) * std::sin(dec)) /
        std::sqrt(std::pow(std::cos(dec) * std::sin(tt), 2) + std::pow(std::sin(latr) *
            std::cos(dec) * std::cos(tt) - std::cos(latr) * std::sin(dec), 2));
    double sqt = 1 - sazi * sazi;
    if (sqt < 0) sqt = 0;
    double azi = 180 + (180 * std::atan(sazi / std::sqrt(sqt))) / pi;
    if (cazi < 0) {
        if (sazi < 0) {
            azi = 180 - azi;
        }
        else {
            azi = 540 - azi;
        }
    }
    // Define and return output variable
    solmodel solpos;
    solpos.zend = z;
    solpos.zenr = z * torad;
    solpos.azid = azi;
    solpos.azir = azi * torad;
    return solpos;
}
// ** Calculates solar index ** //
double solarindexCpp(double slope, double aspect, double zend, double azid, bool shadowmask = false)
{
    double si;
    if (zend > 90.0 && !shadowmask) {
        si = 0;
    }
    else {
        if (slope == 0.0) {
            si = std::cos(zend * torad);
        }
        else {
            si = std::cos(zend * torad) * std::cos(slope * torad) + std::sin(zend * torad) *
                std::sin(slope * torad) * std::cos((azid - aspect) * torad);
        }
    }
    if (si < 0.0) si = 0.0;
    return si;
}
// ** Calculate canopy extinction coefficient for sloped ground surfaces ** //
kstruct cankCpp(double zenr, double x, double si) {
    double k;
    if (zenr > (pi / 2.0)) zenr = pi / 2.0;
    if (si < 0.0) si = 0.0;
    // Calculate normal canopy extinction coefficient
    if (x == 1.0) {
        k = 1.0 / (2.0 * std::cos(zenr));
    }
    else if (std::isinf(x)) {
        k = 1.0;
    }
    else if (x == 0.0) {
        k = std::tan(zenr);
    }
    else {
        k = std::sqrt(x * x + (std::tan(zenr) * std::tan(zenr))) / (x + 1.774 * std::pow((x + 1.182), -0.733));
    }
    if (k > 6000.0) k = 6000.0;
    // Calculate adjusted k
    double kd = k * std::cos(zenr) / si;
    if (si == 0) kd = 1.0;
    double Kc = 1.0 / si;
    if (si == 0.0) Kc = 600.0;
    kstruct kparams;
    kparams.k = k;
    kparams.kd = kd;
    kparams.Kc = Kc;
    return kparams;
}
// ** Calculates parameters for diffuse radiation using two-stream model ** //
tsdifstruct twostreamdifCpp(double pait, double x, double lref, double ltra, double gref)
{
    tsdifstruct params;
    params.om = lref + ltra;
    params.a = 1.0 - params.om;
    params.del = lref - ltra;
    params.J = 1.0 / 3.0;
    if (x != 1.0) {
        double mla = 9.65 * std::pow((3.0 + x), -1.65);
        if (mla > pi / 2.0) mla = pi / 2.0;
        params.J = std::cos(mla) * std::cos(mla);
    }
    params.gma = 0.5 * (params.om + params.J * params.del);
    params.h = std::sqrt(params.a * params.a + 2.0 * params.a * params.gma);
    // Calculate base parameters
    params.S1 = std::exp(-params.h * pait);
    params.u1 = params.a + params.gma * (1.0 - 1.0 / gref);
    double u2 = params.a + params.gma * (1.0 - gref);
    params.D1 = (params.a + params.gma + params.h) * (params.u1 - params.h) * 
        1.0 / params.S1 - (params.a + params.gma - params.h) * (params.u1 + params.h) * params.S1;
    params.D2 = (u2 + params.h) * 1.0 / params.S1 - (u2 - params.h) * params.S1;
    // Calculate parameters
    params.p1 = (params.gma / (params.D1 * params.S1)) * (params.u1 - params.h);
    params.p2 = (-params.gma * params.S1 / params.D1) * (params.u1 + params.h);
    params.p3 = (1.0 / (params.D2 * params.S1)) * (u2 + params.h);
    params.p4 = (-params.S1 / params.D2) * (u2 - params.h);
    // Define and return output 
    return params;
}
// ** Calculates parameters for direct radiation using two-stream model ** //
tsdirstruct twostreamdirCpp(double pait, double om, double a, double gma, double J, double del, double h, double gref,
    double kd, double u1, double S1, double D1, double D2)
{
    tsdirstruct params;
    // Calculate base parameters
    double sig = kd * kd + gma * gma - std::pow((a + gma), 2.0);
    double ss = 0.5 * (om + J * del / kd) * kd;
    double sstr = om * kd - ss;
    double S2 = std::exp(-kd * pait);
    double u2 = a + gma * (1.0 - gref);
    params.p5 = -ss * (a + gma - kd) - gma * sstr;
    double v1 = ss - (params.p5 * (a + gma + kd)) / sig;
    double v2 = ss - gma - (params.p5 / sig) * (u1 + kd);
    params.p6 = (1.0 / D1) * ((v1 / S1) * (u1 - h) - (a + gma - h) * S2 * v2);
    params.p7 = (-1.0 / D1) * ((v1 * S1) * (u1 + h) - (a + gma + h) * S2 * v2);
    params.sig = -sig;
    params.p8 = sstr * (a + gma + kd) - gma * ss;
    double v3 = (sstr + gma * gref - (params.p8 / params.sig) * (u2 - kd)) * S2;
    params.p9 = (-1 / D2) * ((params.p8 / (params.sig * S1)) * (u2 + h) + v3);
    params.p10 = (1 / D2) * (((params.p8 * S1) / params.sig) * (u2 - h) + v3);
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
        tsdifstruct tspdif = twostreamdifCpp(pait, x, lref, ltra, gref);
        double trd = clump * clump;
        // Calculate albedo and ground flux
        double amx = gref;
        if (amx < lref) amx = lref;
        double albd = gref * (trd * trd) + (1.0 - trd * trd) * (tspdif.p1 + tspdif.p2);
        if (albd > amx) albd = amx;
        if (albd < 0.01) albd = 0.01;
        double groundRdd = trd + (1.0 - trd) * (tspdif.p3 * std::exp(-tspdif.h * pait) + tspdif.p4 * std::exp(tspdif.h * pait));
        // Calculate time-variant two-stream parameters (direct)
        for (size_t i = 0; i < Rsw.size(); ++i) {
            if (Rsw[i] > 0.0) {
                // Calculate solar variables
                solmodel solp = solpositionCpp(lat, lon, year[i], month[i], day[i], lt[i]);
                double si = solarindexCpp(slope, aspect, solp.zend, solp.azid);
                if (solp.zenr > pi / 2.0) solp.zenr = pi / 2.0;
                // Calculate canopy extinction coefficient
                double cosz = std::cos(solp.zenr);
                kstruct kp = cankCpp(solp.zenr, x, si);
                // Calculate two-stream parameters (direct)  
                tsdirstruct tspdir = twostreamdirCpp(pait, tspdif.om, tspdif.a, tspdif.gma, tspdif.J, tspdif.del, tspdif.h, gref,
                    kp.kd, tspdif.u1, tspdif.S1, tspdif.D1, tspdif.D2);
                // Calculate beam normalisations
                double Rbeam = (Rsw[i] - Rdif[i]) / cosz;
                if (Rbeam > 1352.0) Rbeam = 1352.0;
                double trb = std::pow(clump, kp.Kc);
                if (trb > 0.999) trb = 0.999;
                if (trb < 0.0) trb = 0.0;
                double Rb = Rbeam * cosz;
                double trg = trb + (1 - trb) * std::exp(-kp.kd * pait); // tranmission to ground though gaps and leaves
                double Rbc = (trg * si + (1 - trg) * cosz) * Rbeam;
                // Calculate albedo and ground flux
                double albb = trd * trb * gref + (1.0 - trd * trb) * (tspdir.p5 / -tspdir.sig + tspdir.p6 + tspdir.p7);
                if (albb > amx) albb = amx;
                if (albb < 0.01) albb = 0.01;
                double groundRbdd = trb + (1.0 - trb) * ((tspdir.p8 / tspdir.sig) * std::exp(-kp.kd * pait) +
                    tspdir.p9 * std::exp(-tspdif.h * pait) + tspdir.p10 * std::exp(tspdif.h * pait));
                if (groundRbdd > amx) groundRbdd = amx;
                if (groundRbdd < 0.0) groundRbdd = 0.0;
                // Calculate canopy absorbed
                radCsw[i] = (1.0 - albd) * Rdif[i] + (1.0 - albb) * Rbc;
                // Calculate ground absorbed
                double Rgdif = groundRdd * Rdif[i] + groundRbdd * Rb;
                radGsw[i] = (1.0 - gref) * (Rgdif + std::exp(-kp.kd * pait) * Rbeam * si);
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
                solmodel solp = solpositionCpp(lat, lon, year[i], month[i], day[i], lt[i]);
                double si = solarindexCpp(slope, aspect, solp.zend, solp.azid);
                if (solp.zenr > pi / 2.0) solp.zenr = pi / 2.0;
                double dirr = (Rsw[i] - Rdif[i]) / std::cos(solp.zenr);
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
    double cp = 2e-05 * std::pow(tc, 2.0) + 0.0002 * tc + 29.119;
    return cp;
}
// ** Calculate zero-plane displacement ** //
// [[Rcpp::export]]
double zeroplanedisCpp(double h, double pai)
{
    if (pai < 0.001) pai = 0.001;
    double d = (1.0 - (1.0 - std::exp(-std::sqrt(7.5 * pai))) / std::sqrt(7.5 * pai)) * h;
    return d;
}
// ** Calculate roughness length ** //
// [[Rcpp::export]]
double roughlengthCpp(double h, double pai, double d, double psi_h)
{
    double Be = std::sqrt(0.003 + (0.2 * pai) / 2);
    double zm = (h - d) * std::exp(-ka / Be) * std::exp(ka * psi_h);
    // safety check to stop diabatic coefficient reversing profile
    if (zm > (0.9 * (h - d))) zm = 0.9 * (h - d);
    if (zm < 0.0005) zm = 0.0005; // sets a minimum
    return zm;
}
// **  Calculate integrated diabatic correction coefficient for momentum ** //
double dpsimCpp(double ze)
{
    double psim;
    // unstable
    if (ze < 0) {
        double x = std::pow((1.0 - 15.0 * ze), 0.25);
        psim = std::log(std::pow((1.0 + x) / 2.0, 2.0) * (1 + std::pow(x, 2.0)) / 2.0) -
            2.0 * std::atan(x) + pi / 2.0;
    }
    // stable
    else {
        psim = -4.7 * ze;
    }
    if (psim < -4.0) psim = -4.0;
    if (psim > 3.0) psim = 3.0;
    return psim;
}
// **  Calculate integrated diabatic correction coefficient for heat ** //
double dpsihCpp(double ze)
{
    double psih;
    // unstable
    if (ze < 0) {
        double y = std::sqrt(1.0 - 9.0 * ze);
        psih = std::log(std::pow((1.0 + y) / 2.0, 2.0));
    }
    // stable
    else {
        psih = -(4.7 * ze) / 0.74;
    }
    if (psih < -4.0) psih = -4.0;
    if (psih > 3.0) psih = 3.0;
    return psih;
}
// **  Calculate diabatic influencing factor for heat ** //  
double dphihCpp(double ze)
{
    double phih;
    // unstable
    if (ze < 0) {
        double phim = 1 / std::pow((1.0 - 16.0 * ze), 0.25);
        phih = std::pow(phim, 2.0);
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
    double dT = 0.7045388 * std::pow((d * std::pow(H, 4.0)), 0.2);
    double gha = 0.0375 * std::pow(dT / d, 0.25);
    if (gha < 0.1) gha = 0.1;
    return gha;
}
// **  Calculate molar conductance above canopy ** //  
double gturbCpp(double uf, double d, double zm, double zref, double ph, double psi_h, double gmin)
{
    double z0 = 0.2 * zm + d; // heat exchange surface height
    double ln = std::log((zref - d) / (z0 - d));
    double g = (ka * ph * uf) / (ln + psi_h);
    if (g < gmin) g = gmin;
    return g;
}
// **  Convert soil water fraction to water potential ** //
double psiwfromthetaCpp(double theta, double Smax, double psi_e, double b)
{
    psi_e = std::abs(psi_e);
    double Se = theta / Smax;
    if (Se > 1.0) Se = 1.0;
    double psiw = -psi_e * std::pow(Se, -b) * 0.01;
    return psiw;
}
// **  Calculate stomatal paramaters ** //
stompstruct stomparamsCpp(double hgt, double lat, double x)
{
    stompstruct out;
    // Initially assume C3 grass
    out.Rsmx = 420.0;
    out.psiw0 = -3.1;
    out.kk = 0.34;
    out.rat = 0.9;
    if (hgt < 1.0 && std::abs(lat) < 22.5) {
        // C4 grass
        out.Rsmx = 450.0;
        out.psiw0 = -2.7;
        out.kk = 0.39;
        out.rat = 0.9;
    }
    if (hgt >= 1.0 && hgt < 7.0) {
        // Shrub
        out.Rsmx = 430.0;
        out.psiw0 = -4.0;
        out.kk = 0.28;
        out.rat = 0.75;
    }
    if (hgt >= 7.0) {
        // tree
        if (std::abs(lat) < 22.5) {
            // Broadleaf evergreen tropical trees
            out.Rsmx = 500.0;
            out.psiw0 = -1.75;
            out.kk = 0.67;
            out.rat = 0.4;
        }
        else {
            if (x < 0.8 || std::abs(lat) > 58.0) {
                // Needleleaf forest
                out.Rsmx = 420.0;
                out.psiw0 = -4.09;
                out.kk = 0.29;
                out.rat = 0.6;
            }
            else {
                // Decidious  forest
                out.Rsmx = 500.0;
                out.psiw0 = -2.51;
                out.kk = 0.46;
                out.rat = 0.45;
            }
        }  // end temperate
    } // end forest
    return out;
}
// ** Calculate stomatal conductance using new function
double stomcondCpp(double Rswabs, double theta, double gsmax, double Smax,
    double psi_e, double b, stompstruct stomp)
{
    // Compute absorbed radiation
    if (Rswabs <= 0.0) return 0.0;
    // Model PAR response
    if (Rswabs > stomp.Rsmx) Rswabs = stomp.Rsmx;
    double gs = gsmax * std::pow(2.0, -(stomp.Rsmx - Rswabs) / (0.2 * stomp.Rsmx));
    // Compute psiw
    double thetan = stomp.rat * theta + (1 - stomp.rat) * thetam;
    double psiw = psiwfromthetaCpp(thetan, Smax, psi_e, b);
    if (psiw < stomp.psiw0) psiw = stomp.psiw0;
    double mu = 1.0 - (std::exp(-stomp.kk * psiw) - 1.0) / (std::exp(-stomp.kk * stomp.psiw0) - 1.0);
    double gs2 = mu * gsmax;
    if (gs > gs2) gs = gs2;
    return gs;
}
// ** Calculate canopy conductance using new function
double canopycondCpp(double Rsw, double Rdif, double k, double om, double theta, double gsmax, double PAI,
    double Smax, double psi_e, double b, stompstruct stomp)
{
    double Gs = 9999.99;
    if (!std::isnan(om)) {
        // Sunlit leaf area
        double P_sun = (1.0 - std::exp(-k * PAI)) / k;
        double P_shade = PAI - P_sun;
        // Absorbed - sun and shaded parts
        double Rshade_abs = Rdif * ((1.0 - std::exp(-PAI)) / PAI) * (1.0 - om);
        double Rsun_abs = (Rsw - Rdif) * k * (1 - om) + Rshade_abs;
        // Calculate stomatal conductance
        double gs_sun = stomcondCpp(Rsun_abs, theta, gsmax, Smax, psi_e, b, stomp);
        double gs_shade = stomcondCpp(Rshade_abs, theta, gsmax, Smax, psi_e, b, stomp);
        Gs = gs_sun * P_sun + gs_shade * P_shade;
    }
    return Gs;
}
// **  Saturated vapour pressure ** //
// [[Rcpp::export]]
double satvapCpp(double tc)
{
    double es;
    if (tc > 0) {
        es = 0.61078 * std::exp(17.27 * tc / (tc + 237.3));
    }
    else {
        es = 0.61078 * std::exp(21.875 * tc / (tc + 265.5));
    }
    return es;
}
// **  Dewpoint temperature ** // 
// [[Rcpp::export]]
double dewpointCpp(double ea) {
    return 243.5 * std::log(ea / 0.6112) /
        (17.67 - std::log(ea / 0.6112));
}
// **  Penman-Monteith equation ** //
double PenmanMonteithCpp(double Rabs, double gHa, double gV, double tc, double te, double pk, double ea, double em, double G, double erh)
{
    double Rema = em * sb * radem(tc);
    double la = 0;
    if (te >= 0) {
        la = 45068.7 - 42.8428 * te;
    }
    else {
        la = 51078.69 - 4.338 * te - 0.06367 * te * te;
    }
    double cp = cpairCpp(te);
    double Da = satvapCpp(tc) - ea;
    double gR = (4.0 * em * sb * std::pow(te + 273.15, 3.0)) / cp;
    double De = satvapCpp(te + 0.5) - satvapCpp(te - 0.5);
    double Ts = tc + ((Rabs - Rema - la * (gV / pk) * Da * erh - G) / (cp * (gHa + gR) + la * (gV / pk) * De * erh));
    return Ts;
}
// **  Function to compute daily from hourly - each value replicated 24 times ** // 
// [[Rcpp::export]]
std::vector<double> hourtodayCpp(std::vector<double> hourly, std::string stat, bool rephour = true) {
    int numDays = hourly.size() / 24;
    std::vector<double> daily(rephour ? numDays * 24 : numDays);
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
        if (rephour) {
            // Fill daily with replicated values for each day
            for (int j = 0; j < 24; ++j) {
                daily[i * 24 + j] = dailyStat;
            }
        } 
        else {
            daily[i] = dailyStat;
        }
    }
    return daily;
}
// **  Function to compute rolling mean temp ** // 
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
// **  Microclimf function to compute rolling mean yearly ** //
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
soilstruct soilpfun(double Vm, double Vq, double Mc, double rho)
{
    soilstruct out;
    double frs = Vm + Vq;
    out.c1 = (0.57 + 1.73 * Vq + 0.93 * Vm) / (1.0 - 0.74 * Vq - 0.49 * Vm) - 2.8 * frs * (1.0 - frs);
    out.c3 = 1.0 + 2.6 * std::pow(Mc, -0.5);
    out.c4 = 0.03 + 0.7 * frs * frs;
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************** Big leaf point model here ********************************************************* //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// **  Function to compute ground heat flux ** //  
Gmodel GFluxCpp(std::vector<double> Tg, std::vector<double> soilm, double rho, double Vm, double Vq, double Mc,
    std::vector<double> Gmax, std::vector<double> Gmin, int iter, bool yearG = true) {
    // Time invariant variables
    double frs = Vm + Vq;
    double c1 = (0.57 + 1.73 * Vq + 0.93 * Vm) / (1.0 - 0.74 * Vq - 0.49 * Vm) - 2.8 * frs * (1.0 - frs);
    double c3 = 1.0 + 2.6 * std::pow(Mc, -0.5);
    double c4 = 0.03 + 0.7 * frs * frs;
    double mu1 = 2400.0 * rho / 2.64;
    double mu2 = 1.06 * rho;
    // Calculate daily mean soil surface temperature
    std::vector<double> Td = hourtodayCpp(Tg, "mean");
    // Initalise variables that need retaining
    std::vector<double> Gmu(Tg.size());
    std::vector<double> dT(Tg.size());
    std::vector<double> k(Tg.size());
    std::vector<double> kap(Tg.size());
    for (size_t i = 0; i < Tg.size(); ++i) {
        // Find soil diffusivity
        double cs = mu1 + 4180 * soilm[i];
        double ph = (rho * (1.0 - soilm[i]) + soilm[i]) * 1000;
        double c2 = mu2 * soilm[i];
        k[i] = c1 + c2 * soilm[i] - (c1 - c4) * std::exp(-pow(c3 * soilm[i], 4.0));
        kap[i] = k[i] / (cs * ph);
        double DD = std::sqrt(2 * kap[i] / omdy);
        Gmu[i] = std::sqrt(2) * (k[i] / DD) * 0.5;
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
        std::vector<double> kama = mayCpp(kap);
        std::vector<double> Gmuy(Tg.size());
        for (size_t i = 0; i < G.size(); ++i) {
            double omyr = (2 * pi) / (G.size() * 3600);
            Gmuy[i] = std::sqrt(2) * kma[i] / std::sqrt(2 * kama[i] / omyr);
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
    // Access items of groundp
    double gref = groundp[0];
    double slope = groundp[1];
    double aspect = groundp[2];
    double groundem = groundp[3];
    double rho = groundp[4];
    double Vm = groundp[5];
    double Vq = groundp[6];
    double Mc = groundp[7];
    double soilb = groundp[8];
    double psie = groundp[9];
    double Smax = groundp[10];
    double Smin = groundp[11];
    // Calculate SW radiation
    radmodel swabs = RadswabsCpp(pai, vegx, lref, ltra, clump, gref, slope, aspect, lat, lon, year, month, day, hour, Rsw, Rdif);
    // Calculate time-invarient variables
    double pait = pai / (1 - clump);
    double trd = (1 - clump * clump) * std::exp(-pait) + clump * clump;
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
    // Create stomatal parameters
    stompstruct stomp = stomparamsCpp(h, lat, vegx);
    double om = 0.5 * (lref + ltra);
    // New variables for storing
    // Initalise H
    std::vector<double> H(Rsw.size());
    for (size_t i = 0; i < H.size(); ++i) H[i] = 0.5 * Rsw[i] - em * sb * radem(tc[i]);
    double tstf = tol * 2;
    int iter = 0;
    double tst = 0;
    while (tstf > tol) {
        tst = 0;
        for (size_t i = 0; i < Rsw.size(); ++i) {
            // Calculate longwave radiation
            double RemC = em * sb * radem(Tc[i]);
            double radClw = em * Rlw[i]; // Longwave radiation down from sky
            double radGlw = groundem * (trd * radClw + (1 - trd) * RemC);
            // Calculate absorbed radiation
            RabsG[i] = swabs.ground[i] + radGlw;
            double RabsC = swabs.canopy[i] + radClw;
            // Calculate canopy temperature
            double zm = roughlengthCpp(h, pai, d, psih[i]);
            uf[i] = (ka * wspeed[i]) / (std::log((zref - d) / zm) + psim[i]);
            if (uf[i] < 0.0002) uf[i] = 0.0002;
            double gmin = gfreeCpp(leafd, std::abs(H[i])) * 2 * pai;
            double ph = phairCpp(tcc[i], pk[i]);
            double gHa = gturbCpp(uf[i], d, zm, zref, ph, psih[i], gmin);
            solmodel solp = solpositionCpp(lat, lon, year[i], month[i], day[i], hour[i]);
            kstruct kp = cankCpp(solp.zenr, vegx, std::cos(solp.zenr));
            double gC = canopycondCpp(Rsw[i], Rdif[i], kp.k, om, soilm[i], gsmax, pai, Smax, psie, soilb, stomp);
            double gV = 1 / (1 / gHa + 1 / gC);
            if (gC == 0) gV = 0;
            double ea = satvapCpp(tc[i]) * rh[i] / 100;
            double Tcn = PenmanMonteithCpp(RabsC, gHa, gV, tc[i], tcc[i], pk[i], ea, em, G[i], 1);
            double tdew = dewpointCpp(ea);
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
            double tst2 = std::abs(Tcn - Tc[i]);
            double tst3 = std::abs(Tgn - Tg[i]);
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
            double Rnet = RabsC - sb * em * radem(Tc[i]);
            if (Rnet > 0 && H[i] > Rnet) H[i] = Rnet;
            // Recalculate stablity variables
            // Stability
            LL[i] = (ph * cp * std::pow(uf[i], 3.0) * Tk) / (-0.4 * 9.81 * H[i]);
            //if (LL[i] > 10000.0) LL[i] = 10000.0;
            //if (LL[i] < -10000.0) LL[i] = -10000.0;
            psim[i] = dpsimCpp(zm / LL[i]) - dpsimCpp((zref - d) / LL[i]);
            psih[i] = dpsihCpp((0.2 * zm) / LL[i]) - dpsihCpp((zref - d) / LL[i]);
            phih[i] = dphihCpp((zref - d) / LL[i]);
            // Set limits to diabatic coefficients
            double ln1 = std::log((zref - d) / zm);
            double ln2 = std::log((zref - d) / (0.2 * zm));
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
    std::vector<double> groundp({ 0.15,0.0,180.0,0.97,1.529643,0.509,0.06,0.5422,5.2,2.6,0.419,0.074 });
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
        double lnr = std::log((zout - d) / zh) / std::log((zin - d) / zh);
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
        double lnru = std::log((zout - d) / zm) / std::log((uzin - d) / zm);
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
        double lwout = sb * 0.95 * radem(temp[i]);
        double lwnet = lwout - lwdown[i];
        rnet[i] = swrad - lwnet;
        if (rnet[i] < 0) rnet[i] = 0;
    }
    // Convert to daily
    std::vector<double> rnetd = hourtodayCpp(rnet, "mean", false);
    std::vector<double> rain = hourtodayCpp(rainh, "sum", false);
    // Run soil model
    std::vector<double> soilm(rain.size());
    double s1 = Smax;
    double s2 = Smax;
    soilm[0] = Smax;
    for (size_t i = 1; i < rain.size(); ++i) {
        double sav = (s1 + s2) / 2;
        double dif = s2 - s1;
        s1 = s1 + rmu * rain[i] - mult * rnetd[i];
        double k = Ksat * std::pow(sav / Smax, pwr);
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
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************* Functions used to drive microclimf ************************************************* //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
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
                ltwi(i, j) = std::log(val) / tfact;
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
// soilmdistribute code (single value)
double soildCpp(double soilm, double Smin, double Smax, double tadd)
{
    double rge = Smax - Smin;
    double theta = (soilm - Smin) / rge;
    if (theta > 0.9999) theta = 0.9999;
    if (theta < 0.0001) theta = 0.0001;
    double lt = std::log(theta / (1 - theta));
    double sm = lt + tadd;
    sm = 1 / (1 + std::exp(-sm));
    double smout = sm * rge + Smin;
    return smout;
}
// Calculate time-invarient vegetation parameters for radiation model
tirstruct twostreamdif(double pai, double paia, double x, double lref, double ltra, double clump, double gref)
{
    tirstruct out;
    // ** Derive diffuse parameters
    out.pait = pai / (1.0 - clump);
    tsdifstruct tspdif = twostreamdifCpp(out.pait, x, lref, ltra, gref);
    out.om = tspdif.om;
    out.omp = 0.5 * out.om;
    out.a = tspdif.a;
    out.gma = tspdif.gma;
    out.J = tspdif.J;
    out.del = tspdif.del;
    out.h = tspdif.h;
    out.u1 = tspdif.u1;
    out.S1 = tspdif.S1;
    out.D1 = tspdif.D1;
    out.D2 = tspdif.D2;
    // ** Calculate things needed prior to computing fluxes
    out.gi = 0.0;
    if (clump > 0.0) out.gi = std::pow(clump, paia / pai); // gaps from canopy top
    if (out.gi > 0.99) out.gi = 0.99;
    double giu = 0.0;
    if (clump > 0.0) giu = std::pow(clump, (pai - paia) / pai);  // gaps upward from ground
    if (giu > 0.99) giu = 0.99;
    double trd = out.gi * out.gi;
    out.trdn = std::pow(clump, 2.0); //  // transmission downward to ground
    out.trdu = giu * giu;
    // paia adjusted for gap fraction
    out.paiaa = paia / (1.0 - out.gi);
    // transmission downward and z
    // ** white-sky albedo
    out.amx = gref;
    if (out.amx < lref) out.amx = lref;
    out.albd = (1.0 - out.trdn * out.trdn) * (tspdif.p1 + tspdif.p2) + out.trdn * out.trdn * gref;
    if (out.albd > out.amx) out.albd = out.amx;
    if (out.albd < 0.01) out.albd = 0.01;
    // ** normalised downward diffuse flux at ground
    out.Rddn_g = (1.0 - out.trdn) * (tspdif.p3 * std::exp(-tspdif.h * out.pait) + tspdif.p4 * std::exp(tspdif.h * out.pait)) + out.trdn;
    if (out.Rddn_g > 1.0)  out.Rddn_g = 1.0;
    if (out.Rddn_g < 0.0)  out.Rddn_g = 0.0;
    // ** normalised upward diffuse flux at z
    out.Rdup_z = (1.0 - out.trdu * out.trdn) * (tspdif.p1 * std::exp(-tspdif.h * out.paiaa) + tspdif.p2 * std::exp(tspdif.h * out.paiaa))
        + out.trdu * out.trdn * gref;
    if (out.Rdup_z > 1.0) out.Rdup_z = 1.0;
    if (out.Rdup_z < 0.0) out.Rdup_z = 0.0;
    // ** normalised downward diffuse flux at z
    out.Rddn_z = (1.0 - trd) * (tspdif.p3 * std::exp(-tspdif.h * out.paiaa) + tspdif.p4 * std::exp(tspdif.h * out.paiaa)) + trd;
    if (out.Rddn_z > 1.0) out.Rddn_z = 1.0;
    if (out.Rddn_z < 0.0) out.Rddn_z = 0.0;
    return out;
}
// ** Run two-stream radiation model ** //
radmodel2 twostreamCpp(double pai, double clump, double gref, double svfa, double si,
    double tc, double Rsw, double Rdif, double Rlw, solmodel solp, kstruct kp, 
    tsdirstruct tspdir, tirstruct tir)
{
    radmodel2 out;
    if (Rsw > 0.0) {
        double cosz = std::cos(solp.zenr);
        if (pai > 0.0) {
            // Calculate gap tranmissions
            double trbn = std::pow(clump, kp.Kc); // direct transmission downward at ground
            if (trbn > 0.999) trbn = 0.999;
            if (trbn < 0.0) trbn = 0.0;
            double trb = std::pow(tir.gi, kp.Kc); // direct transmission downward at z
            if (trb > 0.999) trb = 0.999;
            if (trb < 0.0) trb = 0.0;
            // ** black-sky albedo
            double albb = (1.0 - tir.trdn * trbn) * ((tspdir.p5 / -tspdir.sig) + tspdir.p6 + tspdir.p7) + tir.trdn * trbn * gref;
            if (albb > tir.amx) albb = tir.amx;
            if (albb < 0.01) albb = 0.01;
            // ** normalised contribution of direct to downward diffuse flux at ground
            double Rdbdn_g = (1.0 - trbn) * ((tspdir.p8 / tspdir.sig) * std::exp(-kp.kd * tir.pait) +
                tspdir.p9 * std::exp(-tir.h * tir.pait) + tspdir.p10 * std::exp(tir.h * tir.pait));
            if (Rdbdn_g > tir.amx) Rdbdn_g = tir.amx;
            if (Rdbdn_g < 0.0) Rdbdn_g = 0.0;
            // ** normalised contribution of direct to upward diffuse flux at z
            double Rdbup_z = (1.0 - tir.trdu * trbn) * ((tspdir.p5 / -tspdir.sig) * exp(-kp.kd * tir.paiaa) +
                tspdir.p6 * std::exp(-tir.h * tir.paiaa) + tspdir.p7 * std::exp(tir.h * tir.paiaa))
                + tir.trdu * trbn * gref;
            if (Rdbup_z > tir.amx) Rdbup_z = tir.amx;
            if (Rdbup_z < 0.0) Rdbup_z = 0.0;
            // ** normalised contribution of direct to downward diffuse flux at z
            double Rdbdn_z = (1.0 - trb) * ((tspdir.p8 / tspdir.sig) * std::exp(-kp.kd * tir.paiaa) +
                tspdir.p9 * std::exp(-tir.h * tir.paiaa) + tspdir.p10 * std::exp(tir.h * tir.paiaa));
            if (Rdbdn_z > tir.amx) Rdbdn_z = tir.amx;
            if (Rdbdn_z < 0.0) Rdbdn_z = 0.0;
            // Calculate incident flux
            double Rbeam = (Rsw - Rdif) / cosz;
            if (Rbeam > 1352.0) Rbeam = 1352.0;
            double Rb = Rbeam * cosz;
            double trg = trb + (1 - trb) * std::exp(-kp.kd * tir.pait); // tranmission to ground though gaps and leaves
            double Rbc = (trg * si + (1 - trg) * cosz) * Rbeam;
            // Calculate ground absorbed radiation
            double Rbdn_g = trbn + (1.0 - trbn) * std::exp(-kp.kd * tir.pait);
            if (Rbdn_g > 1.0) Rbdn_g = 1.0;
            if (Rbdn_g < 0.0) Rbdn_g = 0.0;
            out.radGsw = (1.0 - gref) * (tir.Rddn_g * Rdif * svfa + Rdbdn_g * Rb +
                Rbdn_g * Rbeam * si);
            double maxg = (1.0 - gref) * (Rdif * svfa + Rbeam * si);
            if (out.radGsw > maxg) out.radGsw = maxg;
            // Calculate canopy absorbed radiation
            out.radCsw = (1.0 - tir.albd) * Rdif * svfa + (1.0 - albb) * Rbc;
            // Calculate upward and downward fluxes
            out.Rbdown = (trb + (1.0 - trb) * std::exp(-kp.kd * tir.paiaa)) * Rbeam;
            out.Rddown = tir.Rddn_z * Rdif * svfa + Rdbdn_z * Rb;
            out.Rdup = tir.Rdup_z * Rdif * svfa + Rdbup_z * Rb;
            // Leaf absorbed
            out.radLsw = 0.5 * (1.0 - tir.om) * (out.Rddown + out.Rdup + kp.k * cosz * out.Rbdown);
            out.radLpar = 0.5 * (1.0 - tir.omp) * (out.Rddown + out.Rdup + kp.k * cosz * out.Rbdown);
        } // end if pai > 0
        else {
            out.Rbdown = (Rsw - Rdif) / cosz;
            out.Rddown = Rdif * svfa;
            out.Rdup = gref * (Rdif * svfa + (Rsw - Rdif));
            out.radGsw = (1.0 - gref) * (svfa * Rdif + si * out.Rbdown);
            out.radCsw = out.radGsw;
            out.radLsw = 0.0;
            out.radLpar = 0.0;
        } // end if pai = 0
    } // end if Rsw > 0
    else {
        out.Rbdown = 0.0;
        out.Rddown = 0.0;
        out.Rdup = 0.0;
        out.radGsw = 0.0;
        out.radCsw = 0.0;
        out.radLsw = 0.0;
        out.radLpar = 0.0;
    }
    // Longwave
    if (pai > 0.0) {
        double trdif = (1.0 - tir.trdn) * std::exp(-tir.pait) + tir.trdn;
        out.lwout = 0.97 * sb * radem(tc); // Longwave emitted
        out.radGlw = 0.97 * (trdif * svfa * Rlw + (1.0 - trdif) * out.lwout);
        out.radClw = 0.97 * svfa * Rlw;
    }
    else {
        out.lwout = 0.97 * sb * radem(tc);
        out.radGlw = 0.97 * svfa * Rlw;
        out.radClw = out.radGlw;
    }
    out.zend = solp.zend;
    return out;
}
tiwstruct windtiCpp(double h, double pai)
{
    tiwstruct out;
    out.d = zeroplanedisCpp(h, pai);
    out.zm = roughlengthCpp(h, pai, out.d, 0.0);
    if (out.zm < 1e-6) out.zm = 1e-6;
    out.a = pai / h;
    return out;
}
// ** Calculate wind speed (single value) ** //
windmodel windCpp(double reqhgt, double zref, double h, double pai, double uref, double umu, double ws, tiwstruct tiw)
{
    
    windmodel out;
    double ufs = (ka * uref) / std::log((zref - tiw.d) / tiw.zm);
    out.uf = ufs * umu * ws;
    out.uz = out.uf;
    if (reqhgt > 0) {
        if (reqhgt >= h) {
            out.uz = (out.uf / ka) * std::log((reqhgt - tiw.d) / tiw.zm);
        }
        else {
            // Calculate wind speed at height z below canopy
            double uh = (out.uf / ka) * std::log((h - tiw.d) / tiw.zm);
            if (uh < out.uf) uh = out.uf;
            double Be = out.uf / uh;
            if (Be < 0.001) Be = 0.001;
            double Lc = std::pow(0.25 * tiw.a, -1.0);
            double Lm = 2 * std::pow(Be, 3.0) * Lc;
            out.uz = uh * std::exp(Be * (reqhgt - h) / Lm);
        }
        if (out.uz > uref) out.uz = uref;
    }
    // Calculate gHa
    out.gHa = gturbCpp(out.uf, tiw.d, tiw.zm, zref, 43, 0, 0.0001);
    return out;
}
// Simple PenmanMonteith function
penmonstruct PenmanMonteith2Cpp(double Rabs, double gHa, double gV, double tc, double mxtc, double pk, double ea, double es,
    double G, double surfwet, double tdew)
{
    double De = satvapCpp(tc + 0.5) - satvapCpp(tc - 0.5);
    double gHr = gHa + (4 * 0.97 * sb * std::pow(tc + 273.15, 3.0)) / 29.3;
    double Rem = 0.97 * sb * radem(tc);
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
    penmonstruct out;
    out.Ts = dT + tc;
    if (out.Ts < tdew) out.Ts = tdew;
    out.H = 29.3 * gHa * (out.Ts - tc);
    out.L = m * (satvapCpp(out.Ts) - ea) * surfwet;
    out.Rem = 0.97 * sb * radem(out.Ts);
    out.mu = la * (43.0 / pk);
    return out;
}
// Calculate soil conductance
soilkstruct soilcondCpp(double rho, double soilm, soilstruct sp)
{
    // Find soil diffusivity
    double cs = (2400 * rho / 2.64 + 4180.0 * soilm); // specific heat of soil in J / kg / K
    double ph = (rho * (1.0 - soilm) + soilm) * 1000.0; // bulk density in kg / m3
    double c2 = 1.06 * rho * soilm;
    soilkstruct out;
    out.k = sp.c1 + c2 * soilm - (sp.c1 - sp.c4) * exp(-std::pow(sp.c3 * soilm, 4.0)); // Thermal conductivity W / m / K
    double kap = out.k / (cs * ph);
    out.DD = std::pow(2.0 * kap / omdy, 0.5); // Damping depth
    return out;
}
// Calculate soil surface temperature assuming G to be zero
soilmodelG0 soiltempG0(double tc, double es, double ea, double pk, double radGsw, double radGlw, 
    double tdew, double gHa, double soilm, double mxtc, soilpstruct soilparams)
{
    soilmodelG0 out;
    out.radabs = radGsw + radGlw;
    double matric = -std::abs(soilparams.psi_e) * std::pow(soilm / soilparams.Smax, -soilparams.soilb);
    out.surfwet = std::exp((0.018 * matric) / (8.31 * (tc + 273.15)));
    if (out.surfwet > 1.0) out.surfwet = 1.0;
    // Calculate soil surface temperature assuming G to be zero
    penmonstruct pm = PenmanMonteith2Cpp(out.radabs, gHa, gHa, tc, mxtc, pk, ea, es, 0.0, out.surfwet, tdew);
    out.Tg = pm.Ts;
    out.Rnet = out.radabs - pm.Rem;
    return out;
}
// ** Calculate ground surface temperature using hourly data  ** //
soilmodel soiltemp_hrCpp(double tc, double es, double ea, double pk, double radabs, double surfwet, double tdew,
    double gHa, double soilm, double mxtc, double Gp, double dtr, double dtrp, double muGp, double kp, 
    double Rdmx, soilstruct sp, soilpstruct soilparams)
{
    // Diurnal temperature range in soil temperature
    double dtR = dtr / dtrp;
    // Calculate thermal conductance and damping depth
    soilkstruct kDDg = soilcondCpp(soilparams.rho, soilm, sp);
    // Calculate Gmu
    double Gmu = dtR * (kDDg.k * muGp) / (kp * kDDg.DD);
    soilmodel out;
    // Compute ground heat flux
    out.G = Gp * Gmu;
    if (out.G > 0.6 * Rdmx) out.G = 0.6 * Rdmx;
    if (out.G < -0.6 * Rdmx) out.G = -0.6 * Rdmx;
    penmonstruct pm = PenmanMonteith2Cpp(radabs, gHa, gHa, tc, mxtc, pk, ea, es, out.G, surfwet, tdew);
    out.Tg = pm.Ts;
    out.DD = kDDg.DD;
    return out;
}
// Calculate temperature and vapour pressure above canopy
abovecanstruct TVabove(double reqhgt, double zref, double h, double d, double zm, double T0, double tc, double ea, double surfwet = 1.0)
{
    double zh = 0.2 * zm;
    double estl = satvapCpp(T0);
    abovecanstruct out;
    if (reqhgt > (d + zh)) {
        double lnr = log((reqhgt - d) / zh) / log((zref - d) / zh);
        out.Tz = tc + (T0 - tc) * (1 - lnr);
        out.ez = ea + (estl - ea) * surfwet * (1 - lnr);
    }
    else {
        out.Tz = T0;
        out.ez = ea + (estl - ea) * surfwet;
    }
    return out;
}
// Calculate leaf temperature and fluxes
leaftempstruct leaftemp(double Tcan, double Tg, double tc, double mxtc, double pk, double ea, double es, double uz, double tdew,
    double surfwet, double radLsw, double Rddown, double Rbdown, double Rlw, double pai, double paia, double leafd, 
    double gsmax, double PARabs, double theta, double Smax, double psi_e, double soilb, stompstruct stomp)
{
    leaftempstruct out;
    // radiation absorbed by leaf
    double lwcan = 0.97 * sb * radem(Tcan);
    double lwgro = 0.97 * sb * radem(Tg);
    double paig = pai - paia;
    out.lwup = std::exp(-paig) * lwgro + (1 - std::exp(-paig)) * lwcan;
    out.lwdn = std::exp(-paia) * Rlw + (1 - std::exp(-paia)) * lwcan;
    double lwabs = 0.97 * 0.5 * (out.lwup + out.lwdn);
    double leafabs = radLsw + lwabs;
    // ** Conductances
    double gh = 0.135 * std::sqrt(uz / leafd) * 1.4;
    double gV = gh;
    if (gsmax < 999.99) {
        gV = 0.0;
        double gs = stomcondCpp(PARabs, theta, gsmax, Smax, psi_e, soilb, stomp);
        if (gs > 0.0) gV = 1 / (1 / gh + 1 / gs);
    }
    // Temperature
    penmonstruct pm = PenmanMonteith2Cpp(leafabs, gh, gV, tc, mxtc, pk, ea, es, 0.0, surfwet, tdew);
    out.tleaf = pm.Ts;
    out.H = pm.H;
    out.L = pm.L;
    return out;
}
double rhcanopy(double uf, double h, double d, double z)
{
    double a2 = 0.4 * (1.0 - (d / h)) / std::pow(1.25, 2);
    double inth = 4.293251 * h;
    if (z != h) {
        inth = (2.0 * h * ((48 * std::atan((std::sqrt(5.0) * std::sin((pi * z) / h)) /
            (std::cos((pi * z) / h) + 1))) / std::pow(5.0, 1.5) +
            (32.0 * std::sin((pi * z) / h)) / ((std::cos((pi * z) / h) + 1) *
                ((25.0 * std::pow(std::sin((pi * z) / h), 2.0)) /
                    std::pow((std::cos((pi * z) / h) + 1.0), 2.0) + 5.0)))) / pi;
    }
    double mu = uf / (a2 * h) * 1.0 / (uf * uf);
    double rHa = inth * mu;
    if (rHa < 0.001) rHa = 0.001;
    return rHa;
}
double TVbelow(double zref, double z, double d, double h, double pai, double uf,
    double leafden, double Flux, double Fluxz, double SH, double SG, double mxnear)
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
    double near = (3.047519 + 0.128642 * std::log(pai)) * SN;
    if (std::abs(near) > mxnear) {
        if (near > 0.0) {
            near = mxnear;
        }
        else {
            near = -mxnear;
        }
    }
    if (std::isnan(near)) near = 0;
    return near + farg;
}
// Calculate temperature or vapour pressure above ground
abovemodel TVaboveground(double reqhgt, double zref, double tc, double pk, double ea, double es, double tdew, 
    double Rsw, double Rdif, double Rlw, double soilm, double hgt, double pai, double paia, double vegx, double leafd, double leafden,
    double Smin, double Smax, double psi_e, double soilb, double gsmax, double mxtc, stompstruct stomp, tirstruct tir, radmodel2 rvars, 
    tiwstruct tiw, windmodel wvars, soilmodel Gvars)
{
    abovemodel out;
    // Calculate ground and surface wetness
    double eT = satvapCpp(Gvars.Tg) - ea;
    if (eT < 0.001) eT = 0.001;
    double plf = 0.8753 - 1.7126 * std::log(eT);
    double gwet = 1.0 / (1.0 + std::exp(-plf));
    double surfwet = (soilm - Smin) / (Smax - Smin);
    if (surfwet > gwet) gwet = surfwet;
    // Calculate vapour conductivity
    kstruct kp = cankCpp(rvars.zend, vegx, std::cos(rvars.zend * torad));
    double gS = canopycondCpp(Rsw, Rdif, kp.k, tir.omp, soilm, gsmax, pai, Smax, psi_e, soilb, stomp);
    double gV = 0.0;
    if (gS > 0.0) gV = 1.0 / (1.0 / wvars.gHa + 1 / gS);
    // Calculate canopy temperature
    double Rabs = rvars.radCsw + rvars.radClw;
    penmonstruct pm = PenmanMonteith2Cpp(Rabs, wvars.gHa, gV, tc, mxtc, pk, ea, es, Gvars.G, surfwet, tdew);
    double Tcan = pm.Ts;
    double ez = 0;
    if (reqhgt >= hgt) {
        abovecanstruct tv = TVabove(reqhgt, zref, hgt, tiw.d, tiw.zm, Tcan, tc, ea, surfwet);
        out.Tz = tv.Tz;
        out.tleaf = Tcan;
        out.lwup = 0.97 * sb * radem(Tcan);
        out.lwdn = Rlw;
        ez = tv.ez;
    }
    else {
        leaftempstruct tvl = leaftemp(Tcan, Gvars.Tg, tc, mxtc, pk, ea, es, wvars.uz, tdew, surfwet, rvars.radLsw,
            rvars.Rddown, rvars.Rbdown, Rlw, pai, paia, leafd, gsmax, rvars.radLpar, soilm, Smax, psi_e, soilb, stomp);
        out.tleaf = tvl.tleaf;
        // Calculate temperature below canopy
        double Flux = pm.H * (1.0 - std::exp(-pai));
        double Fluxz = tvl.H;
        abovecanstruct tv = TVabove(hgt, zref, hgt, tiw.d, tiw.zm, Tcan, tc, ea, surfwet);
        double SH = tv.Tz * 29.3 * 43.0;
        double SG = Gvars.Tg * 29.3 * 43.0;
        double mxnear = std::abs(out.tleaf - tv.Tz) * 29.3 * 43.0;
        out.Tz = TVbelow(zref, reqhgt, tiw.d, hgt, pai, wvars.uf, leafden, Flux, Fluxz, SH, SG, mxnear) / (29.3 * 43.0);
        // Calculate humidity below canopy
        Flux = pm.L * (1.0 - std::exp(-pai));
        Fluxz = tvl.L;
        SH = tv.ez * pm.mu;
        SG = satvapCpp(Gvars.Tg) * gwet * pm.mu;
        mxnear = std::abs(satvapCpp(out.tleaf) - tv.ez) * pm.mu;
        ez = TVbelow(zref, reqhgt, tiw.d, hgt, pai, wvars.uf, leafden, Flux, Fluxz, SH, SG, mxnear) / pm.mu;
        out.lwdn = tvl.lwdn;
        out.lwup = tvl.lwup;
    }
    out.rh = (ez / satvapCpp(out.Tz)) * 100.0;
    // Set limits on variables
    if (out.rh > 100.0) out.rh = 100.0;
    double tmx = std::max({ out.tleaf, tc, Gvars.Tg, Tcan }) + 2.0;
    double tmn = std::min({ out.tleaf, tc, Gvars.Tg, Tcan }) - 2.0;
    if (out.Tz > tmx) out.Tz = tmx;
    if (out.Tz < tmn) out.Tz = tmn;
    return out;
}
// Calculate temperature below ground (has to be a vector)
std::vector<double> Tbelowgroundv(double reqhgt, std::vector<double> Tg, std::vector<double> Tgp,
    std::vector<double> Tbp, double meanD, double mat, int hiy, bool complete)
{
    std::vector<double> Tz = Tg;
    int tsteps = Tg.size();
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
                for (int i = 0; i < tsteps; ++i) sumT = sumT + Tg[i];
                double meanT = sumT / Tg.size();
                for (int i = 0; i < tsteps; ++i) Tz[i] = meanT;
            }
        }
        // Time sequence incomplete
        else {
            // Calculate Tz at daily damping depth
            // ** Create variable
            std::vector<double> Tzd(Tg.size());
            // ** Calculate mean daily Tbp and difference from this
            std::vector<double> Tbpd = hourtodayCpp(Tbp, "mean", true);
            std::vector<double> Tbpa(tsteps);
            for (int i = 0; i < tsteps; ++i) Tbpa[i] = Tbp[i] - Tbpd[i];
            // ** Calculate grid to point ratio and difference
            std::vector<double> dTggmx = hourtodayCpp(Tg, "max", true);
            std::vector<double> dTggmn = hourtodayCpp(Tg, "min", true);
            std::vector<double> dTggme = hourtodayCpp(Tg, "mean", true);
            std::vector<double> dTgpmx = hourtodayCpp(Tgp, "max", true);
            std::vector<double> dTgpmn = hourtodayCpp(Tgp, "min", true);
            std::vector<double> dTgpme = hourtodayCpp(Tgp, "mean", true);
            double sat = 0.0;
            for (int i = 0; i < tsteps; ++i) {
                double rat = (dTggmx[i] - dTggmn[i]) / (dTgpmx[i] - dTgpmn[i]);
                double dif = dTggme[i] - dTgpme[i];
                Tzd[i] = rat * Tbpa[i] + Tbpd[i] + dif;
                sat = sat + Tzd[i];
            }
            if (nb > 1.0 && nb <= 24.0) {
                double w1 = 1.0 / nb;
                double w2 = nb / 24.0;
                double wgt = w1 / (w1 + w2);
                for (int i = 0; i < tsteps; ++i) Tz[i] = wgt * Tg[i] + (1 - wgt) * Tzd[i];
            }
            if (nb > 24.0) {
                // Calculate Tz at annual damping depth
                if (nb < hiy) {
                    // Calculate daily to annual weights
                    double w1 = 24.0 / nb;
                    double w2 = nb / hiy;
                    double wgt = w1 / (w1 + w2);
                    for (int i = 0; i < tsteps; ++i) Tz[i] = wgt * Tzd[i] + (1 - wgt) * mat;
                }
                else {
                    for (int i = 0; i < tsteps; ++i) Tz[i] = mat;
                } // end nb >= hiy
            } // end nb > 24.0
        } // end none complete
    }
    return Tz;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************* Functions used to drive R version of microclimf **************************************** //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Calculates a 3D array zenith angles (radians) and solar index values
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
    IntegerVector dim = IntegerVector::create(rows, cols, tsteps);
    // Initialise output variables
    NumericVector zen(rows * cols * tsteps);
    NumericVector si(rows * cols * tsteps);
    zen.attr("dim") = dim;
    si.attr("dim") = dim;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double val = slope(i, j);
            if (!NumericMatrix::is_na(val)) {
                for (int k = 0; k < tsteps; k++) {
                    int idx = i + rows * j + cols * rows * k;
                    solmodel sp = solpositionCpp(lats(i, j), lons(i, j), year[k],
                        month[k], day[k], hour[k]);
                    zen[idx] = sp.zenr;
                    if (sp.zend < 90.0) {
                        si[idx] = solarindexCpp(slope(i, j), aspect(i, j), sp.zend, sp.azid);
                    }
                    else {
                        si[idx] = 0.0;
                    }
                } // end k
            } // na check
            else {
                for (int k = 0; k < tsteps; k++) {
                    int idx = i + rows * j + cols * rows * k;
                    si[idx] = NA_REAL;
                    zen[idx] = NA_REAL;
                }
            }
        } // end j
    } // end i
    // Calculate vector of solar azimuths
    int mrows = rows / 2;
    int mcols = cols / 2;
    double lat = lats(mrows, mcols);
    double lon = lons(mrows, mcols);
    NumericVector sazi(tsteps);
    for (int k = 0; k < tsteps; k++) {
        solmodel sp = solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
        sazi[k] = sp.azir;
    }
    // Assign to micro
    micro["zen"] = zen;
    micro["si"] = si;
    micro["sazi"] = sazi;
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
    // Get dimensions
    IntegerVector dim = hgt.attr("dim");
    int rows = dim[0]; // row
    int cols = dim[1]; // column
    int tsteps = dim[2]; // time
    // Create output variables for storing data
    NumericVector radGsw(hgt.size(), NA_REAL);
    NumericVector radGlw(hgt.size(), NA_REAL);
    NumericVector radCsw(hgt.size(), NA_REAL);
    NumericVector radClw(hgt.size(), NA_REAL);
    NumericVector Rbdown(hgt.size(), NA_REAL);
    NumericVector Rddown(hgt.size(), NA_REAL);
    NumericVector Rdup(hgt.size(), NA_REAL);
    NumericVector radLsw(hgt.size(), NA_REAL);
    NumericVector radLpar(hgt.size(), NA_REAL);
    NumericVector lwout(hgt.size(), NA_REAL);
    // Shape outputs
    radGsw.attr("dim") = dim;
    radGlw.attr("dim") = dim;
    radCsw.attr("dim") = dim;
    radClw.attr("dim") = dim;
    Rbdown.attr("dim") = dim;
    Rddown.attr("dim") = dim;
    Rdup.attr("dim") = dim;
    radLsw.attr("dim") = dim;
    radLpar.attr("dim") = dim;
    lwout.attr("dim") = dim;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double val = hgt[i + rows * j + cols * rows];
            if (!NumericVector::is_na(val)) {
                for (int k = 0; k < tsteps; k++) {
                    int idx = i + rows * j + cols * rows * k;
                    double Rsw = difr[idx] + dirr[idx] * std::cos(zen[idx]);
                    solmodel solp;
                    solp.zenr = zen[idx];
                    kstruct kp = cankCpp(zen[idx], x[idx], si[idx]);
                    tirstruct tir = twostreamdif(pai[idx], paia[idx], x[idx], lref[idx], ltra[idx], clump[idx], gref[idx]);
                    tsdirstruct tspdir = twostreamdirCpp(tir.pait, tir.om, tir.a, tir.gma, tir.J, tir.del, tir.h, gref[idx],
                        kp.kd, tir.u1, tir.S1, tir.D1, tir.D2);
                    radmodel2 radm = twostreamCpp(pai[idx], clump[idx], gref[idx], svfa[idx], si[idx], tc[idx],
                        Rsw, difr[idx], lwdown[idx], solp, kp, tspdir, tir);
                    radGsw[idx] = radm.radGsw;
                    radGlw[idx] = radm.radGlw;
                    radCsw[idx] = radm.radCsw;
                    radClw[idx] = radm.radClw;
                    Rbdown[idx] = radm.Rbdown;
                    Rddown[idx] = radm.Rddown;
                    Rdup[idx] = radm.Rdup;
                    radLsw[idx] = radm.radLsw;
                    radLpar[idx] = radm.radLpar;
                    lwout[idx] = radm.lwout;
                } // end time steps
            } // end NA check
        } // end cols
    } // end rows
    // Create output
    micro["radGsw"] = radGsw;
    micro["radGlw"] = radGlw;
    micro["radCsw"] = radCsw;
    micro["radClw"] = radClw;
    micro["Rbdown"] = Rbdown;
    micro["Rddown"] = Rddown;
    micro["Rdup"] = Rdup;
    micro["radLsw"] = radLsw;
    micro["radLpar"] = radLpar;
    micro["lwout"] = lwout;
    return micro;
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
    // Get dimensions
    IntegerVector dim = h.attr("dim");
    int rows = dim[0]; // row
    int cols = dim[1]; // column
    int tsteps = dim[2]; // time
    // Create output variables for storing data
    NumericVector uf(h.size(), NA_REAL);
    NumericVector uz(h.size(), NA_REAL);
    NumericVector gHa(h.size(), NA_REAL);
    // Shape outputs
    uf.attr("dim") = dim;
    uz.attr("dim") = dim;
    gHa.attr("dim") = dim;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double val = h[i + rows * j + cols * rows];
            if (!NumericVector::is_na(val)) {
                for (int k = 0; k < tsteps; k++) {
                    int idx = i + rows * j + cols * rows * k;
                    tiwstruct tiw = windtiCpp(h[idx], pai[idx]);
                    windmodel wmod = windCpp(reqhgt, zref, h[idx], pai[idx], u2[idx], umu[idx], ws[idx], tiw);
                    uf[idx] = wmod.uf;
                    uz[idx] = wmod.uz;
                    gHa[idx] = wmod.gHa;
                } // end time steps
            } // end NA check
        } // end cols
    } // end rows
    // Create output
    micro["uf"] = uf;
    micro["uz"] = uz;
    micro["gHa"] = gHa;
    return micro;
}
// Calculate ground surface temperature: grid
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
    IntegerVector dim = tc.attr("dim");
    int rows = dim[0];
    int cols = dim[1];
    int tsteps = dim[2];
    int ndays = tsteps / 24;
    // Create output variables
    NumericVector Tg(tc.size(), NA_REAL);
    NumericVector G(tc.size(), NA_REAL);
    NumericVector kDDg(tc.size(), NA_REAL);
    // Shape outputs
    Tg.attr("dim") = dim;
    G.attr("dim") = dim;
    kDDg.attr("dim") = dim;
    soilpstruct spa;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = rho(i, j);
            if (!NumericMatrix::is_na(val)) {
                // Calculate mxtc
                double mxtc = -273.15;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    if (tc[idx] > mxtc) mxtc = tc[idx];
                }
                spa.Smax = Smax(i, j); spa.Smin = Smin(i, j); spa.soilb = soilb(i, j); spa.psi_e = psi_e(i, j);
                spa.Vq = Vq(i, j); spa.Vm = Vm(i, j); spa.Mc = Mc(i, j); spa.rho = rho(i, j);
                soilstruct sp = soilpfun(Vm(i, j), Vq(i, j), Mc(i, j), rho(i, j));
                for (int day = 0; day < ndays; ++day) {
                    // Calculate limits
                    double Rmx = -999.9;
                    double tmx = -999.0;
                    double tmn = 999.0;
                    NumericVector surfwet(24);
                    NumericVector radabs(24);
                    for (int hr = 0; hr < 24; ++hr) {
                        int k = day * 24 + hr;
                        int idx = i + rows * j + cols * rows * k;
                        soilmodelG0 sG0 = soiltempG0(tc[idx], es[idx], ea[idx], pk[idx], radGsw[idx], radGlw[idx], tdew[idx],
                            gHa[idx], soilm[idx], mxtc, spa);
                        double Rval = std::abs(sG0.Rnet);
                        if (Rmx < Rval) Rmx = Rval;
                        if (tmx < sG0.Tg) tmx = sG0.Tg;
                        if (tmn > sG0.Tg) tmn = sG0.Tg;
                        surfwet[hr] = sG0.surfwet;
                        radabs[hr] = sG0.radabs;
                    }
                    double dtr = tmx - tmn;
                    for (int hr = 0; hr < 24; ++hr) {
                        int k = day * 24 + hr;
                        int idx = i + rows * j + cols * rows * k;
                        soilmodel stemp = soiltemp_hrCpp(tc[idx], es[idx], ea[idx], pk[idx], radabs[hr], surfwet[hr], tdew[idx],
                            gHa[idx], soilm[idx], mxtc, Gp[idx], dtr, dtrp[idx], muGp[idx], kp[idx], Rmx, sp, spa);
                        Tg[idx] = stemp.Tg;
                        G[idx] = stemp.G;
                        kDDg[idx] = stemp.DD;
                    } // end hour
                } // end day
            } // end NA check
        } // end cols
    } // end rows
    micro["Tg"] = Tg;
    micro["G"] = G;
    micro["kDDg"] = kDDg;
    return micro;
}
// Calculates mxtc - needed by function abovegrid
NumericVector mxtccalc(NumericVector tc) {
    IntegerVector dim = tc.attr("dim");
    int rows = dim[0];
    int cols = dim[1];
    int tsteps = dim[2];
    // Cceate and shape output
    NumericVector mxtc(tc.size(), NA_REAL);
    mxtc.attr("dim") = dim;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = tc[i + rows * j + cols * rows];
            if (!NumericVector::is_na(val)) {
                double mx = -273.15;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    if (tc[idx] > mx) mx = tc[idx];
                }
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    mxtc[idx] = mx;
                } // end k
            } // end NA check
        } // end cols
    } // end row
    return mxtc;
}
// microclimate model as grid above ground
// [[Rcpp::export]]
List abovegrid(double reqhgt, List micro)
{
    // Access items of micro
    double zref = micro["zref"];
    NumericVector zen = micro["zen"];
    NumericVector lat = micro["lats"];
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
    NumericVector radLpar = micro["radLpar"];
    NumericVector Rddown = micro["Rddown"];
    NumericVector Rbdown = micro["Rbdown"];
    // soil variables
    NumericVector soilm = micro["soilm"];
    NumericVector G = micro["G"];
    // Habitat variables
    NumericVector hgt = micro["veghgt"];
    NumericVector vegx = micro["vegx"];
    NumericVector pai = micro["pai"];
    NumericVector paia = micro["paia"];
    NumericVector leafd = micro["leafd"];
    NumericVector leafden = micro["leafden"];
    NumericVector gsmax = micro["gsmax"];
    NumericVector lref = micro["lref"];
    NumericVector ltra = micro["ltra"];
    // Soil variables
    NumericVector Smin = micro["Smin"];
    NumericVector Smax = micro["Smax"];
    NumericVector soilb = micro["soilb"];
    NumericVector psi_e = micro["psi_e"];
    IntegerVector dim = tc.attr("dim");
    int rows = dim[0];
    int cols = dim[1];
    int tsteps = dim[2];
    // Create variable
    NumericVector Tz(tc.size(), NA_REAL);
    NumericVector tleaf(tc.size(), NA_REAL);
    NumericVector rh(tc.size(), NA_REAL);
    NumericVector lwup(tc.size(), NA_REAL);
    // Shape variables
    Tz.attr("dim") = dim;
    tleaf.attr("dim") = dim;
    rh.attr("dim") = dim;
    lwup.attr("dim") = dim;
    tirstruct tir;
    radmodel2 rvars;
    windmodel wvars;
    soilmodel Gvars;
    // Calculate mxtc
    NumericVector mxtc = mxtccalc(tc);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt[i + rows * j + cols * rows];
            if (!NumericVector::is_na(val)) {
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    stompstruct stomp = stomparamsCpp(hgt[idx], lat[idx], vegx[idx]);
                    tir.omp = 0.5 * (lref[idx] + ltra[idx]);
                    //rvars.radGsw = radGsw[idx];
                    //rvars.radGlw = radGlw[idx];
                    rvars.radCsw = radCsw[idx];
                    rvars.radClw = radClw[idx];
                    rvars.Rbdown = Rbdown[idx];
                    rvars.Rddown = Rddown[idx];
                    //rvars.Rdup = Rdup[idx];
                    rvars.radLsw = radLsw[idx];
                    rvars.radLpar = radLpar[idx];
                    //rvars.lwout = lwout[idx];
                    rvars.zend = zen[idx] * 180.0 / pi;
                    tiwstruct tiw = windtiCpp(hgt[idx], pai[idx]);
                    wvars.uf = uf[idx];
                    wvars.uz = uz[idx];
                    wvars.gHa = gHa[idx];
                    Gvars.Tg = T0[idx];
                    Gvars.G = G[idx];
                    double Rsw = difr[idx] + dirr[idx] * cos(zen[idx]);
                    abovemodel tv = TVaboveground(reqhgt, zref, tc[idx], pk[idx], ea[idx], es[idx], tdew[idx], Rsw, difr[idx],
                        lwdn[idx], soilm[idx], hgt[idx], pai[idx], paia[idx], vegx[idx], leafd[idx], leafden[idx],
                        Smin[idx], Smax[idx], psi_e[idx], soilb[idx], gsmax[idx], mxtc[idx], stomp, tir, rvars,
                        tiw, wvars, Gvars);
                    Tz[idx] = tv.Tz;
                    tleaf[idx] = tv.tleaf;
                    rh[idx] = tv.rh;
                    lwup[idx] = tv.lwup;
                    lwdn[idx] = tv.lwdn;
                } // end time step
            } // end na check
            else {
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    lwdn[idx] = NA_REAL;
                }
            }
        } // end cols
    } // end rows;
    // Create output
    List mout;
    mout["Tz"] = Tz;
    mout["tleaf"] = tleaf;
    mout["T0"] = T0;
    mout["soilm"] = soilm;
    mout["relhum"] = rh;
    mout["windspeed"] = uz;
    mout["Rdirdown"] = Rbdown;
    mout["Rdifdown"] = Rddown;
    mout["Rlwdown"] = lwdn;
    mout["Rswup"] = micro["Rdup"];
    mout["Rlwup"] = lwup;
    return mout;
}
// Calculate below ground temperature : grid
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
    IntegerVector dim = Tg.attr("dim");
    int rows = dim[0];
    int cols = dim[1];
    int tsteps = dim[2];
    // Create and reshape output variables
    NumericVector Tz(Tg.size(), NA_REAL);
    Tz.attr("dim") = dim;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = Tg[i + rows * j + cols * rows];
            if (!NumericVector::is_na(val)) {
                // Subset input variables to vector
                std::vector<double> Tgv(tsteps);
                std::vector<double> Tgpv(tsteps);
                std::vector<double> Tbpv(tsteps);
                double sumD = 0.0;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    Tgv[k] = Tg[idx];
                    Tgpv[k] = Tgp[idx];
                    Tbpv[k] = Tbp[idx];
                    sumD += DD[idx];
                }
                double meanD = sumD / static_cast<double>(tsteps);
                // Run below ground vector model
                std::vector<double> Tzv = Tbelowgroundv(reqhgt, Tgv, Tgpv, Tbpv,
                    meanD, mat, hiy, complete);
                // Extract values to Tz
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    Tz[idx] = Tzv[k];
                } // end timestep
            } // end NA check
        } // end cols
    } // end rows
    // Create output
    List mout;
    mout["Tz"] = Tz;
    mout["T0"] = micro["Tg"];
    mout["soilm"] = micro["soilm"];
    return mout;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ********************************* Functions used to drive C++ version of microclimf **************************************** //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Run microclimate model (hourly, static vegetation, data.frame climate input)
// [[Rcpp::export]]
List runmicro1Cpp(DataFrame obstime, DataFrame climdata, DataFrame pointm, List vegp, List soilc,
    double reqhgt, double zref, double lat, double lon, double Sminp, double Smaxp, double tfact,
    bool complete, double mat, std::vector<bool> out)
{
    // Access columns of obstime
    IntegerVector year = obstime["year"];
    IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];
    NumericVector hour = obstime["hour"];
    // Access columns of climdata
    NumericVector tc = climdata["temp"];
    NumericVector es = climdata["es"];
    NumericVector ea = climdata["ea"];
    NumericVector tdew = climdata["tdew"];
    NumericVector pk = climdata["pres"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector lwdown = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"];
    NumericVector wdir = climdata["winddir"];
    // Access columns of pointm
    NumericVector soilmp = pointm["soilm"];
    NumericVector Tgp = pointm["Tg"];
    NumericVector T0p = pointm["T0p"];
    NumericVector Tbp = pointm["Tbp"];
    NumericVector Gp = pointm["G"];
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
    int tsteps = tc.size();
    int ndays = tsteps / 24;
    IntegerVector dim = { rows, cols, tsteps };
    int n = rows * cols * tsteps;
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
    if (out[0]) Tz = NumericVector(n, NA_REAL);
    if (out[1]) tleaf = NumericVector(n, NA_REAL);
    if (out[2]) relhum = NumericVector(n, NA_REAL);
    if (out[3]) soilm = NumericVector(n, NA_REAL);
    if (out[4]) uz = NumericVector(n, NA_REAL);
    if (out[5]) Rdirdown = NumericVector(n, NA_REAL);
    if (out[6]) Rdifdown = NumericVector(n, NA_REAL);
    if (out[7]) Rlwdown = NumericVector(n, NA_REAL);
    if (out[8]) Rswup = NumericVector(n, NA_REAL);
    if (out[9]) Rlwup = NumericVector(n, NA_REAL);
    // Shape output variables
    if (out[0]) Tz.attr("dim") = dim;
    if (out[1]) tleaf.attr("dim") = dim;
    if (out[2]) relhum.attr("dim") = dim;
    if (out[3]) soilm.attr("dim") = dim;
    if (out[4]) uz.attr("dim") = dim;
    if (out[5]) Rdirdown.attr("dim") = dim;
    if (out[6]) Rdifdown.attr("dim") = dim;
    if (out[7]) Rlwdown.attr("dim") = dim;
    if (out[8]) Rswup.attr("dim") = dim;
    if (out[9]) Rlwup.attr("dim") = dim;
    // Calculate things that vary in time, but not space
    IntegerVector sindex(tsteps);
    IntegerVector windex(tsteps);
    NumericVector zend(tsteps);
    NumericVector zenr(tsteps);
    NumericVector azid(tsteps);
    NumericVector azir(tsteps);
    double mxtc = -273.15;
    for (int k = 0; k < tsteps; ++k) {
        solmodel sp = solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
        zend[k] = sp.zend;
        zenr[k] = sp.zenr;
        azid[k] = sp.azid;
        azir[k] = sp.azir;
        sindex[k] = static_cast<int>(std::round(sp.azid / 15.0)) % 24;
        windex[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
        if (tc[k] > mxtc) mxtc = tc[k];
    }
    // Calculate other odds and sods
    int hiy = 365 * 24;
    if (year[0] % 4 == 0) hiy = 366 * 24;
    // Distribute soil moisture
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    // Loop through and run microclimate model
    soilpstruct spa;
    solmodel solp;
    std::vector<double> Tgp2 = as<std::vector<double>>(Tgp);
    std::vector<double> Tbp2 = as<std::vector<double>>(Tbp);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                // Calculate variables that don't vary temporally
                tirstruct tir = twostreamdif(pai(i, j), paia(i, j), x(i, j), lref(i, j), ltra(i, j), clump(i, j), gref(i, j));
                stompstruct stomp = stomparamsCpp(hgt(i, j), lat, x(i, j));
                tiwstruct tiw = windtiCpp(hgt(i, j), pai(i, j));
                spa.Smax = Smax(i, j); spa.Smin = Smin(i, j); spa.soilb = soilb(i, j); spa.psi_e = Psie(i, j);
                spa.Vq = Vq(i, j); spa.Vm = Vm(i, j); spa.Mc = Mc(i, j); spa.rho = rho(i, j);
                soilstruct sp = soilpfun(Vm(i, j), Vq(i, j), Mc(i, j), rho(i, j));
                // Soil variables needed as vectors
                std::vector<double> Tg(tsteps);
                std::vector<double> DD(tsteps);
                for (int dy = 0; dy < ndays; ++dy) {
                    // Limits for soil model
                    double Rmx = -999.9;
                    double tmx = -999.0;
                    double tmn = 999.0;
                    // Daily vectors needed (general)
                    NumericVector surfwet(24);
                    NumericVector radabs(24);
                    NumericVector soilmday(24);
                    // Daily vectors needed (rad model)
                    NumericVector radCsw(24);
                    NumericVector radClw(24);
                    NumericVector Rddown(24);
                    NumericVector Rbdown(24);
                    NumericVector radLsw(24);
                    NumericVector radLpar(24);
                    // Daily vectors needed (wind model)
                    NumericVector uf(24);
                    NumericVector uzday(24);
                    NumericVector gHa(24);
                    for (int hr = 0; hr < 24; ++hr) {
                        int k = dy * 24 + hr;
                        int idx = i + rows * j + cols * rows * k;
                        // Calculate terrain adjusted solar index and wind shelter coefficient
                        double si = solarindexCpp(slope(i, j), aspect(i, j), zend[k], azid[k], true);
                        if (si < 0.0) si = 0.0;
                        double ws = wsa[windex[k] * rows * cols + j * rows + i];
                        double ha = hor[sindex[k] * rows * cols + j * rows + i];
                        double sa = (pi / 2.0) - zenr[k];
                        if (ha > std::tan(sa)) si = 0.0;
                        // Calculate distributed soil moisture
                        double soild = soildCpp(soilmp[k], Smin(i, j), Smax(i, j), tadd(i, j));
                        soilmday[hr] = soild;
                        if (out[3]) soilm[idx] = soild;
                        // Calculate radiation 
                        solp.zend = zend[k];
                        solp.zenr = zenr[k];
                        kstruct kpp = cankCpp(zenr[k], x(i, j), si);
                        tsdirstruct tspdir = twostreamdirCpp(tir.pait, tir.om, tir.a, tir.gma, tir.J, tir.del, tir.h,
                            gref(i, j), kpp.kd, tir.u1, tir.S1, tir.D1, tir.D2);
                        radmodel2 radm = twostreamCpp(pai(i, j), clump(i, j), gref(i, j), svfa(i, j), si, tc[k], Rsw[k], Rdif[k],
                            lwdown[k], solp, kpp, tspdir, tir);
                        radCsw[hr] = radm.radCsw;
                        radClw[hr] = radm.radClw;
                        Rddown[hr] = radm.Rddown;
                        Rbdown[hr] = radm.Rbdown;
                        radLsw[hr] = radm.radLsw;
                        radLpar[hr] = radm.radLpar;
                        if (out[5]) Rdirdown[idx] = radm.Rbdown;
                        if (out[6]) Rdifdown[idx] = radm.Rddown;
                        if (out[8]) Rswup[idx] = radm.Rdup;
                        // Calculate wind speed
                        double reqhgt2 = reqhgt;
                        if (reqhgt2 < 0.00001) reqhgt2 = 0.00001;
                        windmodel windm = windCpp(reqhgt2, zref, hgt(i, j), pai(i, j), u2[k], umu[k], ws, tiw);
                        // Calculate ground surface temperature assuming G to be zero
                        uf[hr] = windm.uf;
                        uzday[hr] = windm.uz;
                        gHa[hr] = windm.gHa;
                        if (out[4]) uz[idx] = windm.uz;
                        soilmodelG0 sG0 = soiltempG0(tc[k], es[k], ea[k], pk[k], radm.radGsw, radm.radGlw, tdew[k],
                            windm.gHa, soild, mxtc, spa);
                        double Rval = std::abs(sG0.Rnet);
                        if (Rmx < Rval) Rmx = Rval;
                        if (tmx < sG0.Tg) tmx = sG0.Tg;
                        if (tmn > sG0.Tg) tmn = sG0.Tg;
                        surfwet[hr] = sG0.surfwet;
                        radabs[hr] = sG0.radabs;
                    } // end hour loop for soil
                    double dtr = tmx - tmn;
                    for (int hr = 0; hr < 24; ++hr) {
                        int k = dy * 24 + hr;
                        int idx = i + rows * j + cols * rows * k;
                        // Derive ground surface temperature
                        soilmodel Gvars = soiltemp_hrCpp(tc[k], es[k], ea[k], pk[k], radabs[hr], surfwet[hr], tdew[k],
                            gHa[hr], soilmday[hr], mxtc, Gp[k], dtr, dtrp[k], muGp[k], kp[k], Rmx, sp, spa);
                        Tg[k] = Gvars.Tg;
                        DD[k] = Gvars.DD;
                        if (reqhgt >= 0.0) {
                            radmodel2 rvars;
                            rvars.zend = zend[k];
                            rvars.radCsw = radCsw[hr];
                            rvars.radClw = radClw[hr];
                            rvars.Rddown = Rddown[hr];
                            rvars.Rbdown = Rbdown[hr];
                            rvars.radLsw = radLsw[hr];
                            rvars.radLpar = radLpar[hr];
                            windmodel wvars;
                            wvars.uf = uf[hr];
                            wvars.uz = uzday[hr];
                            wvars.gHa = gHa[hr];
                            double reqhgt2 = reqhgt;
                            if (reqhgt2 < 0.00001) reqhgt2 = 0.00001;
                            abovemodel tv = TVaboveground(reqhgt2, zref, tc[k], pk[k], ea[k], es[k], tdew[k], Rsw[k], Rdif[k],
                                lwdown[k], soilmday[hr], hgt(i, j), pai(i, j), paia(i, j), x(i, j), leafd(i, j), leafden(i, j),
                                Smin(i, j), Smax(i, j), Psie(i, j), soilb(i, j), gsmax(i, j), mxtc, stomp, tir, rvars,
                                tiw, wvars, Gvars);
                            // Return outputs
                            if (reqhgt > 0.0) {
                                if (out[0]) Tz[idx] = tv.Tz;
                            }
                            else {
                                if (out[0]) Tz[idx] = Tg[k];
                            }
                            if (out[7]) Rlwdown[idx] = tv.lwdn;
                            if (out[9]) Rlwup[idx] = tv.lwup;
                            if (reqhgt > 0.0) {
                                if (out[1]) tleaf[idx] = tv.tleaf;
                                if (out[2]) relhum[idx] = tv.rh;
                            } // end reqhgt > 0
                        } // end reqhgt >= 0
                    } // end hr
                } // end day
                if (reqhgt < 0.0 && out[0]) {
                    // Run below ground vector model
                    double sumD = 0.0;
                    for (int k = 0; k < tsteps; ++k) {
                        sumD += DD[k];
                    }
                    double meanD = sumD / static_cast<double>(tsteps);
                    std::vector<double> Tzv = Tbelowgroundv(reqhgt, Tg, Tgp2, Tbp2,
                        meanD, mat, hiy, complete);
                    for (int k = 0; k < tsteps; ++k) {
                        int idx = i + rows * j + cols * rows * k;
                        Tz[idx] = Tzv[k];
                    } // end k
                } // end if reqhgt < 0
            } // end NA check
        } // end col
    } // end row
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Tz;
    if (out[1]) outp["tleaf"] = tleaf;
    if (out[2]) outp["relhum"] = relhum;
    if (out[3]) outp["soilm"] = soilm;
    if (out[4]) outp["windspeed"] = uz;
    if (out[5]) outp["Rdirdown"] = Rdirdown;
    if (out[6]) outp["Rdifdown"] = Rdifdown;
    if (out[7]) outp["Rlwdown"] = Rlwdown;
    if (out[8]) outp["Rswup"] = Rswup;
    if (out[9]) outp["Rlwup"] = Rlwup;
    return outp;
}
// Run microclimate model(hourly, static vegetation, array climate inputs)
// [[Rcpp::export]]
List runmicro2Cpp(DataFrame obstime, List climdata, List pointm, List vegp, List soilc,
    double reqhgt, double zref, NumericMatrix lats, NumericMatrix lons, double Sminp, double Smaxp, double tfact,
    bool complete, double mat, std::vector<bool> out)
{
    // Access columns of obstime
    IntegerVector year = obstime["year"];
    IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];
    NumericVector hour = obstime["hour"];
    // Accesss climdata
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
    // Access columns of pointm
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
    int ndays = tsteps / 24;
    IntegerVector dim = { rows, cols, tsteps };
    int n = rows * cols * tsteps;
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
    if (out[0]) Tz = NumericVector(n, NA_REAL);
    if (out[1]) tleaf = NumericVector(n, NA_REAL);
    if (out[2]) relhum = NumericVector(n, NA_REAL);
    if (out[3]) soilm = NumericVector(n, NA_REAL);
    if (out[4]) uz = NumericVector(n, NA_REAL);
    if (out[5]) Rdirdown = NumericVector(n, NA_REAL);
    if (out[6]) Rdifdown = NumericVector(n, NA_REAL);
    if (out[7]) Rlwdown = NumericVector(n, NA_REAL);
    if (out[8]) Rswup = NumericVector(n, NA_REAL);
    if (out[9]) Rlwup = NumericVector(n, NA_REAL);
    // Shape output variables
    if (out[0]) Tz.attr("dim") = dim;
    if (out[1]) tleaf.attr("dim") = dim;
    if (out[2]) relhum.attr("dim") = dim;
    if (out[3]) soilm.attr("dim") = dim;
    if (out[4]) uz.attr("dim") = dim;
    if (out[5]) Rdirdown.attr("dim") = dim;
    if (out[6]) Rdifdown.attr("dim") = dim;
    if (out[7]) Rlwdown.attr("dim") = dim;
    if (out[8]) Rswup.attr("dim") = dim;
    if (out[9]) Rlwup.attr("dim") = dim;
    // Calculate things that vary in time, but not space
    IntegerVector windex(tsteps);
    for (int k = 0; k < tsteps; ++k) {
        windex[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    }
    // Calculate other odds and sods
    int hiy = 365 * 24;
    if (year[0] % 4 == 0) hiy = 366 * 24;
    // Distribute soil moisture
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    // Loop through and run microclimate model
    soilpstruct spa;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                // Calculate variables that don't vary temporally
                tirstruct tir = twostreamdif(pai(i, j), paia(i, j), x(i, j), lref(i, j), ltra(i, j), clump(i, j), gref(i, j));
                stompstruct stomp = stomparamsCpp(hgt(i, j), lats(i, j), x(i, j));
                tiwstruct tiw = windtiCpp(hgt(i, j), pai(i, j));
                spa.Smax = Smax(i, j); spa.Smin = Smin(i, j); spa.soilb = soilb(i, j); spa.psi_e = Psie(i, j);
                spa.Vq = Vq(i, j); spa.Vm = Vm(i, j); spa.Mc = Mc(i, j); spa.rho = rho(i, j);
                soilstruct sp = soilpfun(Vm(i, j), Vq(i, j), Mc(i, j), rho(i, j));
                // Soil variables needed as vectors
                std::vector<double> Tg(tsteps);
                std::vector<double> DD(tsteps);
                // Calculate mxtc
                double mxtc = -273.15;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    if (tc[idx] > mxtc) mxtc = tc[idx];
                }
                for (int dy = 0; dy < ndays; ++dy) {
                    // Limits for soil model
                    double Rmx = -999.9;
                    double tmx = -999.0;
                    double tmn = 999.0;
                    // Daily vectors needed (general)
                    NumericVector surfwet(24);
                    NumericVector radabs(24);
                    NumericVector soilmday(24);
                    // Daily vectors needed (rad model)
                    NumericVector radCsw(24);
                    NumericVector radClw(24);
                    NumericVector Rddown(24);
                    NumericVector Rbdown(24);
                    NumericVector radLsw(24);
                    NumericVector radLpar(24);
                    // Daily vectors needed (wind model)
                    NumericVector uf(24);
                    NumericVector uzday(24);
                    NumericVector gHa(24);
                    NumericVector zendday(24);
                    for (int hr = 0; hr < 24; ++hr) {
                        int k = dy * 24 + hr;
                        int idx = i + rows * j + cols * rows * k;
                        // Calculate solar index with terrain shadowing
                        solmodel solp = solpositionCpp(lats(i, j), lons(i, j), year[k], month[k], day[k], hour[k]);
                        zendday[hr] = solp.zend;
                        double si = solarindexCpp(slope(i, j), aspect(i, j), solp.zend, solp.azid);
                        if (si < 0.0) si = 0.0;
                        int sindex = static_cast<int>(std::round(solp.azid / 15)) % 24;
                        double ha = hor[sindex * rows * cols + j * rows + i];
                        double sa = (pi / 2.0) - solp.zenr;
                        if (ha > std::tan(sa)) si = 0.0;
                        double ws = wsa[windex[k] * rows * cols + j * rows + i];
                        // Calculate distributed soil moisture
                        double soild = soildCpp(soilmp[idx], Smin(i, j), Smax(i, j), tadd(i, j));
                        soilmday[hr] = soild;
                        if (out[3]) soilm[idx] = soild;
                        // Calculate radiation 
                        kstruct kpp = cankCpp(solp.zenr, x(i, j), si);
                        tsdirstruct tspdir = twostreamdirCpp(tir.pait, tir.om, tir.a, tir.gma, tir.J, tir.del, tir.h,
                            gref(i, j), kpp.kd, tir.u1, tir.S1, tir.D1, tir.D2);
                        radmodel2 radm = twostreamCpp(pai(i, j), clump(i, j), gref(i, j), svfa(i, j), si, tc[idx], Rsw[idx], Rdif[idx],
                            lwdown[idx], solp, kpp, tspdir, tir);
                        radCsw[hr] = radm.radCsw;
                        radClw[hr] = radm.radClw;
                        Rddown[hr] = radm.Rddown;
                        Rbdown[hr] = radm.Rbdown;
                        radLsw[hr] = radm.radLsw;
                        radLpar[hr] = radm.radLpar;
                        if (out[5]) Rdirdown[idx] = radm.Rbdown;
                        if (out[6]) Rdifdown[idx] = radm.Rddown;
                        if (out[8]) Rswup[idx] = radm.Rdup;
                        // Calculate wind speed
                        double reqhgt2 = reqhgt;
                        if (reqhgt2 < 0.00001) reqhgt2 = 0.00001;
                        windmodel windm = windCpp(reqhgt2, zref, hgt(i, j), pai(i, j), u2[idx], umu[idx], ws, tiw);
                        // Calculate ground surface temperature assuming G to be zero
                        uf[hr] = windm.uf;
                        uzday[hr] = windm.uz;
                        gHa[hr] = windm.gHa;
                        if (out[4]) uz[idx] = windm.uz;
                        soilmodelG0 sG0 = soiltempG0(tc[idx], es[idx], ea[idx], pk[idx], radm.radGsw, radm.radGlw, tdew[idx],
                            windm.gHa, soild, mxtc, spa);
                        double Rval = std::abs(sG0.Rnet);
                        if (Rmx < Rval) Rmx = Rval;
                        if (tmx < sG0.Tg) tmx = sG0.Tg;
                        if (tmn > sG0.Tg) tmn = sG0.Tg;
                        surfwet[hr] = sG0.surfwet;
                        radabs[hr] = sG0.radabs;
                    } // end hour loop for soil
                    double dtr = tmx - tmn;
                    for (int hr = 0; hr < 24; ++hr) {
                        int k = dy * 24 + hr;
                        int idx = i + rows * j + cols * rows * k;
                        // Derive ground surface temperature
                        soilmodel Gvars = soiltemp_hrCpp(tc[idx], es[idx], ea[idx], pk[idx], radabs[hr], surfwet[hr], tdew[idx],
                            gHa[hr], soilmday[hr], mxtc, Gp[idx], dtr, dtrp[idx], muGp[idx], kp[idx], Rmx, sp, spa);
                        Tg[k] = Gvars.Tg;
                        DD[k] = Gvars.DD;
                        if (reqhgt >= 0.0) {
                            radmodel2 rvars;
                            rvars.zend = zendday[hr];
                            rvars.radCsw = radCsw[hr];
                            rvars.radClw = radClw[hr];
                            rvars.Rddown = Rddown[hr];
                            rvars.Rbdown = Rbdown[hr];
                            rvars.radLsw = radLsw[hr];
                            rvars.radLpar = radLpar[hr];
                            windmodel wvars;
                            wvars.uf = uf[hr];
                            wvars.uz = uzday[hr];
                            wvars.gHa = gHa[hr];
                            double reqhgt2 = reqhgt;
                            if (reqhgt2 < 0.00001) reqhgt2 = 0.00001;
                            abovemodel tv = TVaboveground(reqhgt2, zref, tc[idx], pk[idx], ea[idx], es[idx], tdew[idx], Rsw[idx], Rdif[idx],
                                lwdown[idx], soilmday[hr], hgt(i, j), pai(i, j), paia(i, j), x(i, j), leafd(i, j), leafden(i, j),
                                Smin(i, j), Smax(i, j), Psie(i, j), soilb(i, j), gsmax(i, j), mxtc, stomp, tir, rvars,
                                tiw, wvars, Gvars);
                            // Return outputs
                            if (reqhgt > 0.0) {
                                if (out[0]) Tz[idx] = tv.Tz;
                            }
                            else {
                                if (out[0]) Tz[idx] = Tg[k];
                            }
                            if (out[7]) Rlwdown[idx] = tv.lwdn;
                            if (out[9]) Rlwup[idx] = tv.lwup;
                            if (reqhgt > 0.0) {
                                if (out[1]) tleaf[idx] = tv.tleaf;
                                if (out[2]) relhum[idx] = tv.rh;
                            } // end reqhgt > 0
                        } // end reqhgt >= 0
                    } // end hr
                } // end day
                if (reqhgt < 0.0 && out[0]) {
                    // Run below ground vector model
                    std::vector<double> Tgpv(tsteps);
                    std::vector<double> Tbpv(tsteps);
                    double sumD = 0.0;
                    for (int k = 0; k < tsteps; ++k) {
                        int idx = i + rows * j + cols * rows * k;
                        Tgpv[k] = Tgp[idx];
                        Tbpv[k] = Tbp[idx];
                        sumD += DD[k];
                    }
                    double meanD = sumD / static_cast<double>(tsteps);
                    std::vector<double> Tzv = Tbelowgroundv(reqhgt, Tg, Tgpv, Tbpv, meanD, mat, hiy, complete);
                    for (int k = 0; k < tsteps; ++k) {
                        int idx = i + rows * j + cols * rows * k;
                        Tz[idx] = Tzv[k];
                    } // end k
                } // end if reqhgt < 0
            } // end NA check
        } // end col
    } // end row
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Tz;
    if (out[1]) outp["tleaf"] = tleaf;
    if (out[2]) outp["relhum"] = relhum;
    if (out[3]) outp["soilm"] = soilm;
    if (out[4]) outp["windspeed"] = uz;
    if (out[5]) outp["Rdirdown"] = Rdirdown;
    if (out[6]) outp["Rdifdown"] = Rdifdown;
    if (out[7]) outp["Rlwdown"] = Rlwdown;
    if (out[8]) outp["Rswup"] = Rswup;
    if (out[9]) outp["Rlwup"] = Rlwup;
    return outp;
}
// Run microclimate model (hourly, changing vegetation, dataframe climate inputs)
// [[Rcpp::export]]
List runmicro3Cpp(DataFrame dfsel, DataFrame obstime, DataFrame climdata, DataFrame pointm,
    List vegp, List soilc, double reqhgt, double zref, double lat, double lon, double Sminp,
    double Smaxp, double tfact, bool complete, double mat, std::vector<bool> out)
{
    // Extract data from dfsel
    IntegerVector lyr = dfsel["lyr"];
    IntegerVector st = dfsel["st"];
    IntegerVector ed = dfsel["ed"];
    int nlyrs = lyr.size();
    IntegerVector ndays(nlyrs);
    for (int i = 0; i < nlyrs; ++i) {
        int span = ed[i] - st[i] + 1;
        if (span < 24)
            Rcpp::stop("Too many layers in vegp. Max layers must be <= max days");
        ndays[i] = span / 24;
    }
    // Access columns of obstime
    IntegerVector year = obstime["year"];
    IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];
    NumericVector hour = obstime["hour"];
    // Access columns of climdata
    NumericVector tc = climdata["temp"];
    NumericVector es = climdata["es"];
    NumericVector ea = climdata["ea"];
    NumericVector tdew = climdata["tdew"];
    NumericVector pk = climdata["pres"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector lwdown = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"];
    NumericVector wdir = climdata["winddir"]; // not an array
    // Access columns of pointm
    NumericVector soilmp = pointm["soilm"];
    NumericVector Tgp = pointm["Tg"];
    NumericVector T0p = pointm["T0p"];
    NumericVector Tbp = pointm["Tbp"];
    NumericVector Gp = pointm["G"];
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
    // Get dimensions
    int rows = Vq.nrow();
    int cols = Vq.ncol();
    int tsteps = year.size();
    IntegerVector dim = { rows, cols, tsteps };
    int n = rows * cols * tsteps;
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
    if (out[0]) Tz = NumericVector(n, NA_REAL);
    if (out[1]) tleaf = NumericVector(n, NA_REAL);
    if (out[2]) relhum = NumericVector(n, NA_REAL);
    if (out[3]) soilm = NumericVector(n, NA_REAL);
    if (out[4]) uz = NumericVector(n, NA_REAL);
    if (out[5]) Rdirdown = NumericVector(n, NA_REAL);
    if (out[6]) Rdifdown = NumericVector(n, NA_REAL);
    if (out[7]) Rlwdown = NumericVector(n, NA_REAL);
    if (out[8]) Rswup = NumericVector(n, NA_REAL);
    if (out[9]) Rlwup = NumericVector(n, NA_REAL);
    // Shape output variables
    if (out[0]) Tz.attr("dim") = dim;
    if (out[1]) tleaf.attr("dim") = dim;
    if (out[2]) relhum.attr("dim") = dim;
    if (out[3]) soilm.attr("dim") = dim;
    if (out[4]) uz.attr("dim") = dim;
    if (out[5]) Rdirdown.attr("dim") = dim;
    if (out[6]) Rdifdown.attr("dim") = dim;
    if (out[7]) Rlwdown.attr("dim") = dim;
    if (out[8]) Rswup.attr("dim") = dim;
    if (out[9]) Rlwup.attr("dim") = dim;
    // Calculate things that vary in time, but not space
    IntegerVector sindex(tsteps);
    IntegerVector windex(tsteps);
    NumericVector zend(tsteps);
    NumericVector zenr(tsteps);
    NumericVector azid(tsteps);
    NumericVector azir(tsteps);
    double mxtc = -273.15;
    for (int k = 0; k < tsteps; ++k) {
        solmodel sp = solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
        zend[k] = sp.zend;
        zenr[k] = sp.zenr;
        azid[k] = sp.azid;
        azir[k] = sp.azir;
        sindex[k] = static_cast<int>(std::round(sp.azid / 15.0)) % 24;
        windex[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
        if (tc[k] > mxtc) mxtc = tc[k];
    }
    // Calculate other odds and sods
    int hiy = 365 * 24;
    if (year[0] % 4 == 0) hiy = 366 * 24;
    // Distribute soil moisture
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    // Loop through and run microclimate model
    soilpstruct spa;
    solmodel solp;
    std::vector<double> Tgp2 = as<std::vector<double>>(Tgp);
    std::vector<double> Tbp2 = as<std::vector<double>>(Tbp);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt[i + rows * j];
            if (!NumericVector::is_na(val)) {
                // Soil variables needed as vectors for complete time series
                std::vector<double> Tg(tsteps);
                std::vector<double> DD(tsteps);
                for (int lyr = 0; lyr < nlyrs; ++lyr) {
                    int idxl = i + rows * j + cols * rows * lyr;
                    // Calculate variables that don't vary temporally except by lyr
                    tirstruct tir = twostreamdif(pai[idxl], paia[idxl], x[idxl], lref[idxl], ltra[idxl], clump[idxl], gref(i, j));
                    stompstruct stomp = stomparamsCpp(hgt[idxl], lat, x[idxl]);
                    tiwstruct tiw = windtiCpp(hgt[idxl], pai[idxl]);
                    spa.Smax = Smax(i, j); spa.Smin = Smin(i, j); spa.soilb = soilb(i, j); spa.psi_e = Psie(i, j);
                    spa.Vq = Vq(i, j); spa.Vm = Vm(i, j); spa.Mc = Mc(i, j); spa.rho = rho(i, j);
                    soilstruct sp = soilpfun(Vm(i, j), Vq(i, j), Mc(i, j), rho(i, j));
                    for (int dy = 0; dy < ndays[lyr]; ++dy) {
                        // Limits for soil model
                        double Rmx = -999.9;
                        double tmx = -999.0;
                        double tmn = 999.0;
                        // Daily vectors needed (general)
                        NumericVector surfwet(24);
                        NumericVector radabs(24);
                        NumericVector soilmday(24);
                        // Daily vectors needed (rad model)
                        NumericVector radCsw(24);
                        NumericVector radClw(24);
                        NumericVector Rddown(24);
                        NumericVector Rbdown(24);
                        NumericVector radLsw(24);
                        NumericVector radLpar(24);
                        // Daily vectors needed (wind model)
                        NumericVector uf(24);
                        NumericVector uzday(24);
                        NumericVector gHa(24);
                        for (int hr = 0; hr < 24; ++hr) {
                            int k = dy * 24 + hr + st[lyr];
                            int idx = i + rows * j + cols * rows * k;
                            // Calculate terrain adjusted solar index and wind shelter coefficient
                            double si = solarindexCpp(slope(i, j), aspect(i, j), zend[k], azid[k], true);
                            if (si < 0.0) si = 0.0;
                            double ws = wsa[windex[k] * rows * cols + j * rows + i];
                            double ha = hor[sindex[k] * rows * cols + j * rows + i];
                            double sa = (pi / 2.0) - zenr[k];
                            if (ha > std::tan(sa)) si = 0.0;
                            // Calculate distributed soil moisture
                            double soild = soildCpp(soilmp[k], Smin(i, j), Smax(i, j), tadd(i, j));
                            soilmday[hr] = soild;
                            if (out[3]) soilm[idx] = soild;
                            // Calculate radiation 
                            solp.zend = zend[k];
                            solp.zenr = zenr[k];
                            kstruct kpp = cankCpp(zenr[k], x[idxl], si);
                            tsdirstruct tspdir = twostreamdirCpp(tir.pait, tir.om, tir.a, tir.gma, tir.J, tir.del, tir.h,
                                gref(i, j), kpp.kd, tir.u1, tir.S1, tir.D1, tir.D2);
                            radmodel2 radm = twostreamCpp(pai[idxl], clump[idxl], gref(i, j), svfa(i, j), si, tc[k], Rsw[k], Rdif[k],
                                lwdown[k], solp, kpp, tspdir, tir);
                            radCsw[hr] = radm.radCsw;
                            radClw[hr] = radm.radClw;
                            Rddown[hr] = radm.Rddown;
                            Rbdown[hr] = radm.Rbdown;
                            radLsw[hr] = radm.radLsw;
                            radLpar[hr] = radm.radLpar;
                            if (out[5]) Rdirdown[idx] = radm.Rbdown;
                            if (out[6]) Rdifdown[idx] = radm.Rddown;
                            if (out[8]) Rswup[idx] = radm.Rdup;
                            // Calculate wind speed
                            double reqhgt2 = reqhgt;
                            if (reqhgt2 < 0.00001) reqhgt2 = 0.00001;
                            windmodel windm = windCpp(reqhgt2, zref, hgt[idxl], pai[idxl], u2[k], umu[k], ws, tiw);
                            // Calculate ground surface temperature assuming G to be zero
                            uf[hr] = windm.uf;
                            uzday[hr] = windm.uz;
                            gHa[hr] = windm.gHa;
                            if (out[4]) uz[idx] = windm.uz;
                            soilmodelG0 sG0 = soiltempG0(tc[k], es[k], ea[k], pk[k], radm.radGsw, radm.radGlw, tdew[k],
                                windm.gHa, soild, mxtc, spa);
                            double Rval = std::abs(sG0.Rnet);
                            if (Rmx < Rval) Rmx = Rval;
                            if (tmx < sG0.Tg) tmx = sG0.Tg;
                            if (tmn > sG0.Tg) tmn = sG0.Tg;
                            surfwet[hr] = sG0.surfwet;
                            radabs[hr] = sG0.radabs;
                        } // end hour loop for soil
                        double dtr = tmx - tmn;
                        for (int hr = 0; hr < 24; ++hr) {
                            int k = dy * 24 + hr + st[lyr];
                            int idx = i + rows * j + cols * rows * k;
                            // Derive ground surface temperature
                            soilmodel Gvars = soiltemp_hrCpp(tc[k], es[k], ea[k], pk[k], radabs[hr], surfwet[hr], tdew[k],
                                gHa[hr], soilmday[hr], mxtc, Gp[k], dtr, dtrp[k], muGp[k], kp[k], Rmx, sp, spa);
                            Tg[k] = Gvars.Tg;
                            DD[k] = Gvars.DD;
                            if (reqhgt >= 0.0) {
                                radmodel2 rvars;
                                rvars.zend = zend[k];
                                rvars.radCsw = radCsw[hr];
                                rvars.radClw = radClw[hr];
                                rvars.Rddown = Rddown[hr];
                                rvars.Rbdown = Rbdown[hr];
                                rvars.radLsw = radLsw[hr];
                                rvars.radLpar = radLpar[hr];
                                windmodel wvars;
                                wvars.uf = uf[hr];
                                wvars.uz = uzday[hr];
                                wvars.gHa = gHa[hr];
                                double reqhgt2 = reqhgt;
                                if (reqhgt2 < 0.00001) reqhgt2 = 0.00001;
                                abovemodel tv = TVaboveground(reqhgt2, zref, tc[k], pk[k], ea[k], es[k], tdew[k], Rsw[k], Rdif[k],
                                    lwdown[k], soilmday[hr], hgt[idxl], pai[idxl], paia[idxl], x[idxl], leafd[idxl], leafden[idxl],
                                    Smin(i, j), Smax(i, j), Psie(i, j), soilb(i, j), gsmax[idxl], mxtc, stomp, tir, rvars,
                                    tiw, wvars, Gvars);
                                // Return outputs
                                if (reqhgt > 0.0) {
                                    if (out[0]) Tz[idx] = tv.Tz;
                                }
                                else {
                                    if (out[0]) Tz[idx] = Tg[k];
                                }
                                if (out[7]) Rlwdown[idx] = tv.lwdn;
                                if (out[9]) Rlwup[idx] = tv.lwup;
                                if (reqhgt > 0.0) {
                                    if (out[1]) tleaf[idx] = tv.tleaf;
                                    if (out[2]) relhum[idx] = tv.rh;
                                } // end reqhgt > 0
                            } // end reqhgt >= 0
                        } // end hr
                    } // end day
                } // end lyr
                if (reqhgt < 0.0 && out[0]) {
                    // Run below ground vector model
                    double sumD = 0.0;
                    for (int k = 0; k < tsteps; ++k) {
                        sumD += DD[k];
                    }
                    double meanD = sumD / static_cast<double>(tsteps);
                    std::vector<double> Tzv = Tbelowgroundv(reqhgt, Tg, Tgp2, Tbp2,
                        meanD, mat, hiy, complete);
                    for (int k = 0; k < tsteps; ++k) {
                        int idx = i + rows * j + cols * rows * k;
                        Tz[idx] = Tzv[k];
                    } // end k
                } // end if reqhgt < 0
            } // end NA check
        } // end col
    } // end row
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Tz;
    if (out[1]) outp["tleaf"] = tleaf;
    if (out[2]) outp["relhum"] = relhum;
    if (out[3]) outp["soilm"] = soilm;
    if (out[4]) outp["windspeed"] = uz;
    if (out[5]) outp["Rdirdown"] = Rdirdown;
    if (out[6]) outp["Rdifdown"] = Rdifdown;
    if (out[7]) outp["Rlwdown"] = Rlwdown;
    if (out[8]) outp["Rswup"] = Rswup;
    if (out[9]) outp["Rlwup"] = Rlwup;
    return outp;
}
// Run microclimate model(hourly, changing vegetation, dataframe climate inputs)
// [[Rcpp::export]]
List runmicro4Cpp(DataFrame dfsel, DataFrame obstime, List climdata, List pointm, List vegp,
    List soilc, double reqhgt, double zref, NumericMatrix lats, NumericMatrix lons, double Sminp,
    double Smaxp, double tfact, bool complete, double mat, std::vector<bool> out)
{
    // Extract data from dfsel
    IntegerVector lyr = dfsel["lyr"];
    IntegerVector st = dfsel["st"];
    IntegerVector ed = dfsel["ed"];
    int nlyrs = lyr.size();
    IntegerVector ndays(nlyrs);
    for (int i = 0; i < nlyrs; ++i) {
        int span = ed[i] - st[i] + 1;
        if (span < 24)
            Rcpp::stop("Too many layers in vegp. Max layers must be <= max days");
        ndays[i] = span / 24;
    }
    // Access columns of obstime
    IntegerVector year = obstime["year"];
    IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];
    NumericVector hour = obstime["hour"];
    // Accesss climdata
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
    // Access columns of pointm
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
    // Get dimensions
    int rows = Smin.nrow();
    int cols = Smin.ncol();
    int tsteps = year.size();
    IntegerVector dim = { rows, cols, tsteps };
    int n = rows * cols * tsteps;
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
    if (out[0]) Tz = NumericVector(n, NA_REAL);
    if (out[1]) tleaf = NumericVector(n, NA_REAL);
    if (out[2]) relhum = NumericVector(n, NA_REAL);
    if (out[3]) soilm = NumericVector(n, NA_REAL);
    if (out[4]) uz = NumericVector(n, NA_REAL);
    if (out[5]) Rdirdown = NumericVector(n, NA_REAL);
    if (out[6]) Rdifdown = NumericVector(n, NA_REAL);
    if (out[7]) Rlwdown = NumericVector(n, NA_REAL);
    if (out[8]) Rswup = NumericVector(n, NA_REAL);
    if (out[9]) Rlwup = NumericVector(n, NA_REAL);
    // Shape output variables
    if (out[0]) Tz.attr("dim") = dim;
    if (out[1]) tleaf.attr("dim") = dim;
    if (out[2]) relhum.attr("dim") = dim;
    if (out[3]) soilm.attr("dim") = dim;
    if (out[4]) uz.attr("dim") = dim;
    if (out[5]) Rdirdown.attr("dim") = dim;
    if (out[6]) Rdifdown.attr("dim") = dim;
    if (out[7]) Rlwdown.attr("dim") = dim;
    if (out[8]) Rswup.attr("dim") = dim;
    if (out[9]) Rlwup.attr("dim") = dim;
    // Calculate things that vary in time, but not space
    IntegerVector windex(tsteps);
    for (int k = 0; k < tsteps; ++k) {
        windex[k] = static_cast<int>(std::round(wdir[k] / 45)) % 8;
    }
    // Calculate other odds and sods
    int hiy = 365 * 24;
    if (year[0] % 4 == 0) hiy = 366 * 24;
    // Distribute soil moisture
    NumericMatrix tadd = soildCppm(twi, Sminp, Smaxp, tfact);
    // Loop through and run microclimate model
    soilpstruct spa;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt[i + rows * j];
            if (!NumericVector::is_na(val)) {
                // Calculate mxtc
                double mxtc = -273.15;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    if (tc[idx] > mxtc) mxtc = tc[idx];
                }
                // Soil variables needed as vectors for complete time series
                std::vector<double> Tg(tsteps);
                std::vector<double> DD(tsteps);
                for (int lyr = 0; lyr < nlyrs; ++lyr) {
                    int idxl = i + rows * j + cols * rows * lyr;
                    // Calculate variables that don't vary temporally except by lyr
                    tirstruct tir = twostreamdif(pai[idxl], paia[idxl], x[idxl], lref[idxl], ltra[idxl], clump[idxl], gref(i, j));
                    stompstruct stomp = stomparamsCpp(hgt[idxl], lats(i, j), x[idxl]);
                    tiwstruct tiw = windtiCpp(hgt[idxl], pai[idxl]);
                    spa.Smax = Smax(i, j); spa.Smin = Smin(i, j); spa.soilb = soilb(i, j); spa.psi_e = Psie(i, j);
                    spa.Vq = Vq(i, j); spa.Vm = Vm(i, j); spa.Mc = Mc(i, j); spa.rho = rho(i, j);
                    soilstruct sp = soilpfun(Vm(i, j), Vq(i, j), Mc(i, j), rho(i, j));
                    for (int dy = 0; dy < ndays[lyr]; ++dy) {
                        // Limits for soil model
                        double Rmx = -999.9;
                        double tmx = -999.0;
                        double tmn = 999.0;
                        // Daily vectors needed (general)
                        NumericVector surfwet(24);
                        NumericVector radabs(24);
                        NumericVector soilmday(24);
                        // Daily vectors needed (rad model)
                        NumericVector radCsw(24);
                        NumericVector radClw(24);
                        NumericVector Rddown(24);
                        NumericVector Rbdown(24);
                        NumericVector radLsw(24);
                        NumericVector radLpar(24);
                        // Daily vectors needed (wind model)
                        NumericVector uf(24);
                        NumericVector uzday(24);
                        NumericVector gHa(24);
                        NumericVector zendday(24);
                        for (int hr = 0; hr < 24; ++hr) {
                            int k = dy * 24 + hr + st[lyr];
                            int idx = i + rows * j + cols * rows * k;
                            // Calculate solar model
                            solmodel solp = solpositionCpp(lats(i, j), lons(i, j), year[k], month[k], day[k], hour[k]);
                            // Calculate terrain adjusted solar index and wind shelter coefficient
                            zendday[hr] = solp.zend;
                            double si = solarindexCpp(slope(i, j), aspect(i, j), solp.zend, solp.azid);
                            if (si < 0.0) si = 0.0;
                            int sindex = static_cast<int>(std::round(solp.azid / 15)) % 24;
                            double ws = wsa[windex[k] * rows * cols + j * rows + i];
                            double ha = hor[sindex * rows * cols + j * rows + i];
                            double sa = (pi / 2.0) - solp.zenr;
                            if (ha > std::tan(sa)) si = 0.0;
                            // Calculate distributed soil moisture
                            double soild = soildCpp(soilmp[idx], Smin(i, j), Smax(i, j), tadd(i, j));
                            soilmday[hr] = soild;
                            if (out[3]) soilm[idx] = soild;
                            // Calculate radiation 
                            kstruct kpp = cankCpp(solp.zenr, x[idxl], si);
                            tsdirstruct tspdir = twostreamdirCpp(tir.pait, tir.om, tir.a, tir.gma, tir.J, tir.del, tir.h,
                                gref(i, j), kpp.kd, tir.u1, tir.S1, tir.D1, tir.D2);
                            radmodel2 radm = twostreamCpp(pai[idxl], clump[idxl], gref(i, j), svfa(i, j), si, tc[idx], Rsw[idx], Rdif[idx],
                                lwdown[idx], solp, kpp, tspdir, tir);
                            radCsw[hr] = radm.radCsw;
                            radClw[hr] = radm.radClw;
                            Rddown[hr] = radm.Rddown;
                            Rbdown[hr] = radm.Rbdown;
                            radLsw[hr] = radm.radLsw;
                            radLpar[hr] = radm.radLpar;
                            if (out[5]) Rdirdown[idx] = radm.Rbdown;
                            if (out[6]) Rdifdown[idx] = radm.Rddown;
                            if (out[8]) Rswup[idx] = radm.Rdup;
                            // Calculate wind speed
                            double reqhgt2 = reqhgt;
                            if (reqhgt2 < 0.00001) reqhgt2 = 0.00001;
                            windmodel windm = windCpp(reqhgt2, zref, hgt[idxl], pai[idxl], u2[idx], umu[idx], ws, tiw);
                            // Calculate ground surface temperature assuming G to be zero
                            uf[hr] = windm.uf;
                            uzday[hr] = windm.uz;
                            gHa[hr] = windm.gHa;
                            if (out[4]) uz[idx] = windm.uz;
                            soilmodelG0 sG0 = soiltempG0(tc[idx], es[idx], ea[idx], pk[idx], radm.radGsw, radm.radGlw, tdew[idx],
                                windm.gHa, soild, mxtc, spa);
                            double Rval = std::abs(sG0.Rnet);
                            if (Rmx < Rval) Rmx = Rval;
                            if (tmx < sG0.Tg) tmx = sG0.Tg;
                            if (tmn > sG0.Tg) tmn = sG0.Tg;
                            surfwet[hr] = sG0.surfwet;
                            radabs[hr] = sG0.radabs;
                        } // end hour loop for soil
                        double dtr = tmx - tmn;
                        for (int hr = 0; hr < 24; ++hr) {
                            int k = dy * 24 + hr + st[lyr];
                            int idx = i + rows * j + cols * rows * k;
                            // Derive ground surface temperature
                            soilmodel Gvars = soiltemp_hrCpp(tc[idx], es[idx], ea[idx], pk[idx], radabs[hr], surfwet[hr], tdew[idx],
                                gHa[hr], soilmday[hr], mxtc, Gp[idx], dtr, dtrp[idx], muGp[idx], kp[idx], Rmx, sp, spa);
                            Tg[k] = Gvars.Tg;
                            DD[k] = Gvars.DD;
                            if (reqhgt >= 0.0) {
                                radmodel2 rvars;
                                rvars.zend = zendday[hr];
                                rvars.radCsw = radCsw[hr];
                                rvars.radClw = radClw[hr];
                                rvars.Rddown = Rddown[hr];
                                rvars.Rbdown = Rbdown[hr];
                                rvars.radLsw = radLsw[hr];
                                rvars.radLpar = radLpar[hr];
                                windmodel wvars;
                                wvars.uf = uf[hr];
                                wvars.uz = uzday[hr];
                                wvars.gHa = gHa[hr];
                                double reqhgt2 = reqhgt;
                                if (reqhgt2 < 0.00001) reqhgt2 = 0.00001;
                                abovemodel tv = TVaboveground(reqhgt2, zref, tc[idx], pk[idx], ea[idx], es[idx], tdew[idx], Rsw[idx], 
                                    Rdif[idx], lwdown[idx], soilmday[hr], hgt[idxl], pai[idxl], paia[idxl], x[idxl], leafd[idxl], 
                                    leafden[idxl], Smin(i, j), Smax(i, j), Psie(i, j), soilb(i, j), gsmax[idxl], mxtc, stomp, tir, rvars,
                                    tiw, wvars, Gvars);
                                // Return outputs
                                if (reqhgt > 0.0) {
                                    if (out[0]) Tz[idx] = tv.Tz;
                                }
                                else {
                                    if (out[0]) Tz[idx] = Tg[k];
                                }
                                if (out[7]) Rlwdown[idx] = tv.lwdn;
                                if (out[9]) Rlwup[idx] = tv.lwup;
                                if (reqhgt > 0.0) {
                                    if (out[1]) tleaf[idx] = tv.tleaf;
                                    if (out[2]) relhum[idx] = tv.rh;
                                } // end reqhgt > 0
                            } // end reqhgt >= 0
                        } // end hr
                    } // end day
                } // end lyr
                if (reqhgt < 0.0 && out[0]) {
                    std::vector<double> Tgp2(tsteps);
                    std::vector<double> Tbp2(tsteps);
                    // Run below ground vector model
                    double sumD = 0.0;
                    for (int k = 0; k < tsteps; ++k) {
                        int idx = i + rows * j + cols * rows * k;
                        sumD += DD[k];
                        Tgp2[k] = Tgp[idx];
                        Tbp2[k] = Tbp[idx];
                    }
                    double meanD = sumD / static_cast<double>(tsteps);
                    std::vector<double> Tzv = Tbelowgroundv(reqhgt, Tg, Tgp2, Tbp2,
                        meanD, mat, hiy, complete);
                    for (int k = 0; k < tsteps; ++k) {
                        int idx = i + rows * j + cols * rows * k;
                        Tz[idx] = Tzv[k];
                    } // end k
                } // end if reqhgt < 0
            } // end NA check
        } // end col
    } // end row
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Tz;
    if (out[1]) outp["tleaf"] = tleaf;
    if (out[2]) outp["relhum"] = relhum;
    if (out[3]) outp["soilm"] = soilm;
    if (out[4]) outp["windspeed"] = uz;
    if (out[5]) outp["Rdirdown"] = Rdirdown;
    if (out[6]) outp["Rdifdown"] = Rdifdown;
    if (out[7]) outp["Rlwdown"] = Rlwdown;
    if (out[8]) outp["Rswup"] = Rswup;
    if (out[9]) outp["Rlwup"] = Rlwup;
    return outp;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************** BIOCLIM model from here *********************************************************** //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
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
        sum_squared_diff += std::pow(vec[i] - mean, 2.0);
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
// Used to pad a NumericMatrix with NAs
NumericMatrix bioclimfill(int rows, int cols)
{
    NumericMatrix bio(rows, cols);
    std::fill(bio.begin(), bio.end(), NA_REAL);
    return bio;
}
// Extract bioclim
List runbioclimCpp(NumericVector Tz, NumericVector soilm, std::vector<bool> out, 
    IntegerVector wetq, IntegerVector dryq, IntegerVector hotq, IntegerVector colq)
{
    // Get dims
    IntegerVector dims = Tz.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    // Initialize variables
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
    if (out[0]) bio1 = bioclimfill(rows, cols);
    if (out[1]) bio2 = bioclimfill(rows, cols);
    if (out[2]) bio3 = bioclimfill(rows, cols);
    if (out[3]) bio4 = bioclimfill(rows, cols);
    if (out[4]) bio5 = bioclimfill(rows, cols);
    if (out[5]) bio6 = bioclimfill(rows, cols);
    if (out[6]) bio7 = bioclimfill(rows, cols);
    if (out[7]) bio8 = bioclimfill(rows, cols);
    if (out[8]) bio9 = bioclimfill(rows, cols);
    if (out[9]) bio10 = bioclimfill(rows, cols);
    if (out[10]) bio11 = bioclimfill(rows, cols);
    if (out[11]) bio12 = bioclimfill(rows, cols);
    if (out[12]) bio13 = bioclimfill(rows, cols);
    if (out[13]) bio14 = bioclimfill(rows, cols);
    if (out[14]) bio15 = bioclimfill(rows, cols);
    if (out[15]) bio16 = bioclimfill(rows, cols);
    if (out[16]) bio17 = bioclimfill(rows, cols);
    if (out[17]) bio18 = bioclimfill(rows, cols);
    if (out[18]) bio19 = bioclimfill(rows, cols);
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
                if (out[0]) bio1(i, j) = bioclim1(Tzv);
                if (out[1]) bio2(i, j) = bioclim2(Tzv);
                if (out[3]) bio4(i, j) = bioclim4(Tzv);
                if (out[4]) bio5(i, j) = bioclim5(Tzv);
                if (out[5]) bio6(i, j) = bioclim6(Tzv);
                if (out[7]) bio8(i, j) = bioclim8(Tzv, wetq);
                if (out[8]) bio9(i, j) = bioclim9(Tzv, dryq);
                if (out[9]) bio10(i, j) = bioclim10(Tzv, hotq);
                if (out[10]) bio11(i, j) = bioclim11(Tzv, colq);
                if (out[11]) bio12(i, j) = bioclim12(soilmv);
                if (out[12]) bio13(i, j) = bioclim13(soilmv);
                if (out[13]) bio14(i, j) = bioclim14(soilmv);
                if (out[14]) bio15(i, j) = bioclim15(soilmv);
                if (out[15]) bio16(i, j) = bioclim16(soilmv, wetq);
                if (out[16]) bio17(i, j) = bioclim17(soilmv, dryq);
                if (out[17]) bio18(i, j) = bioclim18(soilmv, hotq);
                if (out[18]) bio19(i, j) = bioclim19(soilmv, colq);
                if (out[6]) bio7(i, j) = bio5(i, j) - bio6(i, j);
                if (out[2]) bio3(i, j) = bio2(i, j) / bio7(i, j);
            } // end NA check
        } // end col
    } // end row
    // Assign to list
    Rcpp::List outp;
    if (out[0]) outp["bio1"] = bio1;
    if (out[1]) outp["bio2"] = bio2;
    if (out[2]) outp["bio3"] = bio3;
    if (out[3]) outp["bio4"] = bio4;
    if (out[4]) outp["bio5"] = bio5;
    if (out[5]) outp["bio6"] = bio6;
    if (out[6]) outp["bio7"] = bio7;
    if (out[7]) outp["bio8"] = bio8;
    if (out[8]) outp["bio9"] = bio9;
    if (out[9]) outp["bio10"] = bio10;
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

// Run microclimate model(hourly, static vegetation, data.frame climate input)
// [[Rcpp::export]]
List runbioclim1Cpp(DataFrame obstime, DataFrame climdata, DataFrame pointm, List vegp, List soilc,
    double reqhgt, double zref, double lat, double lon, double Sminp, double Smaxp, double tfact,
    double mat, std::vector<bool> out, IntegerVector wetq, IntegerVector dryq, IntegerVector hotq,
    IntegerVector colq, bool air) 
{
    std::vector<bool> outm;
    if (air) {
        outm = { true, false, false, true, false, false, false, false, false, false };
    }
    else {
        outm = { false, true, false, true, false, false, false, false, false, false };
    }
    // Run microclimate model
    List mout = runmicro1Cpp(obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lat, lon, Sminp, Smaxp, tfact,
        true, mat, outm);
    NumericVector Tz;
    if (air) {
        Tz = mout["Tz"];
    }
    else {
        Tz = mout["tleaf"];
    }
    NumericVector soilm = mout["soilm"];
    List bout = runbioclimCpp(Tz, soilm, out, wetq, dryq, hotq, colq);
    return bout;
}
// Run microclimate model(hourly, static vegetation, array climate inputs)
// [[Rcpp::export]]
List runbioclim2Cpp(DataFrame obstime, List climdata, List pointm, List vegp, List soilc,
    double reqhgt, double zref, NumericMatrix lats, NumericMatrix lons, double Sminp, double Smaxp,
    double tfact, double mat, std::vector<bool> out, IntegerVector wetq, IntegerVector dryq,
    IntegerVector hotq, IntegerVector colq, bool air)
{
    // Create output
    std::vector<bool> outm;
    if (air) {
        outm = { true, false, false, true, false, false, false, false, false, false };
    }
    else {
        outm = { false, true, false, true, false, false, false, false, false, false };
    }
    // Run microclimate model
    List mout = runmicro2Cpp(obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lats, lons, Sminp, Smaxp, tfact,
        true, mat, outm);
    NumericVector Tz;
    if (air) {
        Tz = mout["Tz"];
    }
    else {
        Tz = mout["tleaf"];
    }
    NumericVector soilm = mout["soilm"];
    List bout = runbioclimCpp(Tz, soilm, out, wetq, dryq, hotq, colq);
    return bout;
}
// Run bioclim model(hourly, changing vegetation, data.frame climate input)
// [[Rcpp::export]]
List runbioclim3Cpp(DataFrame obstime, DataFrame climdata, DataFrame pointm, List vegp, List soilc,
    double reqhgt, double zref, double lat, double lon, double Sminp, double Smaxp, double tfact,
    double mat, std::vector<bool> out, IntegerVector wetq, IntegerVector dryq, IntegerVector hotq,
    IntegerVector colq, bool air)
{
    // Create output
    std::vector<bool> outm;
    if (air) {
        outm = { true, false, false, true, false, false, false, false, false, false };
    }
    else {
        outm = { false, true, false, true, false, false, false, false, false, false };
    }
    // Create dfsel
    IntegerVector lyr(14);
    IntegerVector st(14);
    IntegerVector ed(14);
    for (int i = 0; i < 14; ++i) {
        lyr[i] = i + 1;
        st[i] = i * 24;
        ed[i] = i * 24 + 23;
    }
    DataFrame dfsel;
    dfsel["lyr"] = lyr;
    dfsel["st"] = st;
    dfsel["ed"] = ed;
    NumericVector Tz;
    List mout = runmicro3Cpp(dfsel, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lat, lon, Sminp,
        Smaxp, tfact, true, mat, outm);
    if (air) {
        Tz = mout["Tz"];
    }
    else {
        Tz = mout["tleaf"];
    }
    NumericVector soilm = mout["soilm"];
    List bout = runbioclimCpp(Tz, soilm, out, wetq, dryq, hotq, colq);
    return bout;
}
// Run bioclim model (hourly, changing vegetation, array climate input)
// [[Rcpp::export]]
List runbioclim4Cpp(DataFrame obstime, List climdata, List pointm, List vegp,
    List soilc, double reqhgt, double zref, NumericMatrix lats, NumericMatrix lons, double Sminp,
    double Smaxp, double tfact, double mat, std::vector<bool> out, IntegerVector wetq, IntegerVector dryq,
    IntegerVector hotq, IntegerVector colq, bool air)
{
    // Create output
    std::vector<bool> outm;
    if (air) {
        outm = { true, false, false, true, false, false, false, false, false, false };
    }
    else {
        outm = { false, true, false, true, false, false, false, false, false, false };
    }
    // Create dfsel
    IntegerVector lyr(14);
    IntegerVector st(14);
    IntegerVector ed(14);
    for (int i = 0; i < 14; ++i) {
        lyr[i] = i + 1;
        st[i] = i * 24;
        ed[i] = i * 24 + 23;
    }
    DataFrame dfsel;
    dfsel["lyr"] = lyr;
    dfsel["st"] = st;
    dfsel["ed"] = ed;
    List mout = runmicro4Cpp(dfsel, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lats, lons, Sminp,
        Smaxp, tfact, true, mat, outm);
    NumericVector Tz;
    if (air) {
        Tz = mout["Tz"];
    }
    else {
        Tz = mout["tleaf"];
    }
    NumericVector soilm = mout["soilm"];
    List bout = runbioclimCpp(Tz, soilm, out, wetq, dryq, hotq, colq);
    return bout;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************* Snow energy and mass balance model ************************************************* //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
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
    double Be = std::sqrt(0.003 + (0.2 * pai) / 2.0);
    double uh = uf / Be;
    double a = pai / hgt;
    double Lc = std::pow(0.25 * a, -1.0);
    double Lm = 2.0 * std::pow(Be, 3.0) * Lc;
    double k1 = Be / Lm;
    double uzm = (uh / (hgt * k1)) * (1 - exp(-k1 * hgt));
    if (uzm < uf) uzm = uf;
    // Calculate snow interception
    double rhos = 67.92 + 51.25 * exp(tc / 2.59); // fresh snow density (kg/m^3)
    double S = Sh * (0.26 + 46 / rhos); // maximum snow load per unit branch area (kg/m^3)
    double Lstr = S * pai;  // maximum canopy snow load (mm SWE)
    double Z = std::atan(uzm / 0.8); // zenith angle of snow fall direction (0.8m/s = terminal velocity of snow flake)
    double kc = 1.0 / (2.0 * std::cos(Z)); // extinction coefficient(from Campbell assuming spherical distribution)
    double Cp = 1.0 - std::exp(-kc * pai); // effective canopy cover perpendicular to direction of snow flake
    double k2 = Cp / Lstr; // dimensionless proportionality factor
    double I1 = (Lstr - Li) * (1.0 - std::exp(-k2 * prec)); // the intercepted snow load at the start of unloading(mm SWE)
    double cis = I1 * 0.678; // Canopy snow interception(mm SWE)
    if (cis > prec) cis = prec;
    return cis;
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
NumericVector snowalbCpp(NumericVector prec) {
    int tsteps = prec.size();
    IntegerVector hs(tsteps);
    hs[0] = 0;
    for (int i = 1; i < tsteps; ++i) {
        if (prec[i] > 0) {
            hs[i] = 0;
        }
        else {
            hs[i] = hs[i - 1] + 1;
        }
    }
    NumericVector alb(tsteps);  // snow albedo (0-1)
    for (int i = 0; i < tsteps; ++i) {
        alb[i] = (-9.8740 * std::log(hs[i] / 24) + 78.3434) / 100.0;
        if (alb[i] > 0.95) alb[i] = 0.95;
        if (alb[i] < 0.1) alb[i] = 0.1;
    }
    return alb;
}
// Calculates absorbed radiation for snow
snowrad radoneB(obspoint obstime, climpoint clim, vegpoint vegp, snowpoint snow,
    otherpoint other)
{
    snowrad out;
    // ************ Calculate longwave radiation ************* //
    double RlwabsC = 0.97 * clim.Rlw;
    out.RlwabsG = RlwabsC;
    double cld = vegp.clump * vegp.clump;
    double pait = vegp.pai / (1 - vegp.clump);
    out.tr = (1 - cld) * std::exp(-pait) + cld;
    if (vegp.hgt > 0) {
        double Rsky = out.tr * clim.Rlw;
        double Rcan = (1.0 - out.tr) * 0.97 * sb * radem(clim.Tci);
        out.RlwabsG = 0.97 * (Rsky + Rcan);
    }
    // ************ Calculate shortwave radiation ************* //
    out.RabsC = RlwabsC;
    out.RswabsG = 0.0;
    if (clim.Rsw > 0.0) {
        // *** Calculate Solar variables ********************* //
        solmodel solp = solpositionCpp(other.lat, other.lon, obstime.year, obstime.month, obstime.day, obstime.hour);
        double si = solarindexCpp(other.slope, other.aspect, solp.zend, solp.azid);
        if (solp.zend > 90.0) solp.zend = 90.0;
        if (si < 0.0) si = 0.0;
        // *** Calculate radiation absorbed by canopy *** //
        double cosz = std::cos(solp.zenr);
        double Rbeam = (clim.Rsw - clim.Rdif) / cosz;
        if (Rbeam > 1352.2) Rbeam = 1352.2;
        double RswabsC = (1.0 - snow.alb) * (clim.Rdif + Rbeam * cosz);
        out.RabsC = RswabsC + RlwabsC;
        // *** Calculate shortwave radiation absorbed by ground *** //
        out.RswabsG = RswabsC;
        if (vegp.hgt > 0.0) {
            if ((snow.alb + vegp.ltra) > 0.999) vegp.ltra = 0.999 - snow.alb;
            tsdifstruct tspdif = twostreamdifCpp(pait, 1.0, snow.alb, vegp.ltra, snow.alb);
            // Calculate canopy extinction coefficient
            kstruct kp = cankCpp(solp.zenr, 1.0, si);
            tsdirstruct tspdir = twostreamdirCpp(pait, tspdif.om, tspdif.a, tspdif.gma, tspdif.J, tspdif.del,
                tspdif.h, snow.alb, kp.kd, tspdif.u1, tspdif.S1, tspdif.D1, tspdif.D2);
            // Downward diffuse stream
            double clb = std::pow(vegp.clump, kp.Kc);
            double  Rddm = (1.0 - cld) * (tspdif.p3 * std::exp(-tspdif.h * pait) + tspdif.p4 * std::exp(tspdif.h * pait)) + cld;
            if (Rddm > 1.0) Rddm = 1.0;
            if (Rddm < 0.0) Rddm = 0.0;
            // Contribution of direct to downward diffuse stream
            double Rdbm = (1.0 - clb) * ((tspdir.p8 / tspdir.sig) * std::exp(-kp.kd * pait) +
                tspdir.p9 * std::exp(-tspdif.h * pait) + tspdir.p10 * std::exp(tspdif.h * pait));
            if (Rdbm > 1.0) Rdbm = 1.0;
            if (Rdbm < 0.0) Rdbm = 0.0;
            // Downward direct stream
            double Rbgm = (1.0 - clb) * std::exp(-kp.kd * pait) + clb;
            if (Rbgm > 1.0) Rbgm = 1.0;
            if (Rbgm < 0.0) Rbgm = 0.0;
            // Radiation absorbed by ground
            double RdifG = (1.0 - snow.alb) * (Rdbm * Rbeam * cosz) + Rddm * clim.Rdif;
            double RdirG = (1.0 - snow.alb) * (Rbgm * Rbeam * 0.5);
            out.RswabsG = RdifG + RdirG;
        } // end if not snow covered
    } // end if daytime
    return out;
}
// One point in time, once cell, bigleaf.
snowmodpoint snowoneB(obspoint obstime, climpoint clim, vegpoint vegp, snowpoint snow,
    otherpoint other, std::vector<double> sdp, double umu = 1.0)
{
    snowmodpoint out;
    // *** Adjust veg parameters for presence of snow
    double pai = 0.0;
    if (vegp.hgt > snow.sdepg) {
        pai = vegp.pai * (vegp.hgt - snow.sdepg) / vegp.hgt;
    }
    double hgt = vegp.hgt - snow.sdepg;
    if (hgt < 0.0) hgt = 0.0;
    double zi = 0.0;
    if (snow.sdepg > 0.0) zi = ((snow.sdepc - snow.sdepg) * snow.sdenc) / (hgt * 1000.0);
    double ltra = vegp.ltra * std::exp(-10.1 * zi);
    // Run radiation model
    vegp.hgt = hgt;
    vegp.ltra = ltra;
    vegp.pai = pai;
    snowrad rad = radoneB(obstime, clim, vegp, snow, other);
    double RabsG = rad.RswabsG + rad.RlwabsG;
    // *** Calculate temperature of ground and canopy snow surfaces *** //
    // Calculate convective conductivity
    double d = 0.0;
    double zm = 0.005;
    if (vegp.hgt > 0.0) {
        d = zeroplanedisCpp(hgt, pai);
        zm = roughlengthCpp(hgt, pai, d, other.psih);
    }
    if (zm < 0.0009) zm = 0.0009;
    out.hgt = hgt;
    out.pai = pai;
    out.uf = (ka * clim.u2) / (std::log((other.zref - d) / zm) + other.psim);
    out.uf = out.uf * umu;
    double ph = phairCpp(clim.tc, clim.pk);
    out.gHa = gturbCpp(out.uf, d, zm, other.zref, ph, other.psih, 0.03);
    // Calculate temperatures
    out.Tc = PenmanMonteithCpp(rad.RabsC, out.gHa, out.gHa, clim.tc, clim.te, clim.pk, clim.ea, 0.97, other.G, 1.0);
    out.Tg = PenmanMonteithCpp(RabsG, out.gHa, out.gHa, clim.tc, clim.te, clim.pk, clim.ea, 0.97, other.G, 1.0);
    double tdew = dewpointCpp(clim.ea);
    if (out.Tc < tdew) out.Tc = tdew;
    if (out.Tg < tdew) out.Tg = tdew;
    // ******* Calculate mass balance of snowpack (canopy + ground) ******
    // Sublimation
    double la;
    if (out.Tc < 0.0) {
        la = 51078.69 - 4.338 * out.Tc - 0.06367 * out.Tc * out.Tc;
    }
    else {
        la = 45068.7 - 42.8428 * out.Tc;
    }
    double L = la * (out.gHa / clim.pk) * (satvapCpp(out.Tc) - clim.ea);
    la = la / 0.018015; // Conversion to J/kg
    out.mSc = (L / la) * 3.6; // m SWE sublimation
    // Melt if snowpack temperature above zero (m SWE)
    double Tcp = out.Tc;
    out.mMc = 0.0;
    if (out.Tc > 0.0) {
        double S = snow.sdepc * (snow.sdenc / 1000);
        double Fm = 583.3 * out.Tc * S;
        out.mMc = (Fm / 334000.0) * 3.6;
        if (snow.sdepc > 0.0) out.Tc = 0.0;
    }
    // Rain melt (m SWE)
    out.mRc = 0.0;
    if (clim.tc > 0.0) {
        out.mRc = 0.0125 * clim.tc * clim.prec / 1000;
    }
    // ******* Calculate mass balance of snowpack (ground only) ******
    // Sublimation 
    if (out.Tg < 0.0) {
        la = 51078.69 - 4.338 * out.Tg - 0.06367 * out.Tg * out.Tg;
    }
    else {
        la = 45068.7 - 42.8428 * out.Tg;
    }
    double mu = std::exp(-vegp.pai);
    if (mu > 1.0) mu = 1.0;
    L = la * (out.gHa / clim.pk) * (satvapCpp(out.Tg) - clim.ea) * mu;
    la = la / 0.018015; // Conversion to J/kg
    out.mSg = (L / la) * 3.6; // m SWE sublimation
    // Melt if snowpack temperature above zero (m SWE)
    out.mMg = 0.0;
    if (out.Tg > 0.0) {
        double S = snow.sdepg * (snow.sdeng / 1000.0);
        double Fm = 583.3 * out.Tg * S;
        out.mMg = (Fm / 334000.0) * 3.6;
        if (snow.sdepg > 0.0) out.Tg = 0.0;
    }
    // Canopy interception
    double Li = 0.0;
    if (snow.sdepc > 0.0) {
        double wgtg = snow.sdepg / snow.sdepc;
        double sdencc = wgtg * snow.sdeng + (1.0 - wgtg) * snow.sdenc;
        Li = (snow.sdepc - snow.sdepg) * sdencc;
    }
    if (Li < 0.0) Li = 0.0;
    out.cis = canopysnowintCpp(vegp.hgt, vegp.pai, out.uf, clim.prec, clim.tc, Li);
    if (out.cis > clim.prec) out.cis = clim.prec;
    // Rain melt (m SWE)
    out.mRg = 0.0;
    if (clim.tc > 0.0) {
        out.mRg = 0.0125 * clim.tc * (clim.prec - out.cis) / 1000.0;
    }
    // Re-calculate change in SWE
    double snowc = clim.prec;
    double snowg = clim.prec - out.cis;
    if (clim.tc > 2.0) {
        snowc = 0.0;
        snowg = 0.0;
    }
    double swec = snowc / 1000.0 - out.mSc - out.mMc - out.mRc; // units are m
    double sweg = snowg / 1000.0 - out.mSg - out.mMg - out.mRg; // units are m
    // Recalculate snow density and age
    out.snowagec = snow.snowagec + 1.0;
    out.snowageg = snow.snowageg + 1.0;
    out.sdenc = ((sdp[0] - sdp[1]) * (1.0 - std::exp(-sdp[2] * snow.sdepc / 100.0 -
        sdp[3] * out.snowagec / 24.0)) + sdp[1]) * 1000.0;
    out.sdeng = ((sdp[0] - sdp[1]) * (1.0 - std::exp(-sdp[2] * snow.sdepg / 100.0 -
        sdp[3] * out.snowageg / 24.0)) + sdp[1]) * 1000.0;
    out.sdepc = snow.sdepc + (swec * 1000.0) / out.sdenc; // units are metres
    out.sdepg = snow.sdepg + (sweg * 1000.0) / out.sdeng; // units are metres
    if (out.sdepc < 0.0) {
        out.sdepc = 0.0;
        out.snowagec = 0.0;
    }
    if (out.sdepg < 0.0) {
        out.sdepg = 0.0;
        out.snowageg = 0.0;
    }
    // ********************** return outputs ******************
    out.RswabsG = rad.RswabsG;
    out.RlwabsG = rad.RlwabsG;
    out.tr = rad.tr;
    out.Tcp = Tcp;
    return out;
}
// **  Function to compute rate of heat storage by snow ** //  
NumericVector GFluxCppsnow(NumericVector snowt, NumericVector snowden)
{
    // Initalise variables that need retaining
    int tsteps = snowt.size();
    std::vector<double> Gmu(tsteps);
    std::vector<double> dT(tsteps);
    // Calculate daily mean soil surface temperature
    std::vector<double> Td = hourtodayCpp(Rcpp::as<std::vector<double>>(snowt), "mean");
    for (int i = 0; i < tsteps; ++i) {
        // Dive by 1000 to convert snow density in kg/m^3 to g/cm^3
        double k = 0.0442 * std::exp(5.181 * snowden[i] / 1000);
        double kap = k / (snowden[i] * 2090);
        double DD = std::sqrt(2.0 * kap / omdy);
        Gmu[i] = std::sqrt(2.0) * (k / DD) * 0.5;
        // Calculate T fluctuation from daily mean
        dT[i] = snowt[i] - Td[i];
    }
    // Calculate 6 hour back rolling mean of Gmu and dT to get 3 hour lag
    std::vector<double> Gmud = maCpp(Gmu, 6);
    std::vector<double> G = maCpp(dT, 6);
    for (int i = 0; i < tsteps; ++i) G[i] = G[i] * Gmud[i] * 1.1171;
    NumericVector Gv = Rcpp::wrap(G);
    return Gv;
}
// One point, once cell, bigleaf.
// [[Rcpp::export]]
List pointmodelsnow(DataFrame obstime, DataFrame climdata, NumericVector vegp,
    NumericVector other, std::string snowenv, double tol = 0.5, double maxiter = 100)
{
    // Access columns of obstime
    IntegerVector year = obstime["year"];
    IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];
    NumericVector hour = obstime["hour"];
    // Access columns of climdata
    NumericVector tc = climdata["temp"];
    NumericVector rh = climdata["relhum"];
    NumericVector pk = climdata["pres"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector Rlw = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"];
    NumericVector prec = climdata["precip"];
    NumericVector ea(tc.size());
    int tsteps = tc.size();
    for (int i = 0; i < tsteps; ++i) ea[i] = satvapCpp(tc[i]) * rh[i] / 100.0;
    NumericVector te = tc;
    // Extract other
    double slope = other[0];
    double aspect = other[1];
    double lat = other[2];
    double lon = other[3];
    double zref = other[4];
    double isnowd = other[5]; // initial snow depth
    double isnowa = other[6]; // initial snow age
    // Calculate snow albedo
    NumericVector salb = snowalbCpp(prec);
    // Have a first stab at guessing H
    NumericVector H(tsteps);
    // Sensible heat flux
    for (int i = 0; i < tsteps; ++i) {
        double Rabs = (1 - salb[i]) * Rsw[i] + 0.97 * Rlw[i];
        H[i] = 0.5 * Rabs;
    }
    // Have a first stab at guessing snow density
    std::vector<double> sdp = snowdenp(snowenv);
    NumericVector sdenc(tsteps);
    for (int i = 0; i < tsteps; ++i) {
        sdenc[i] = ((sdp[0] - sdp[1]) * (1 - exp(-sdp[2] * isnowd / 100.0 -
            sdp[3] * 0)) + sdp[1]) * 1000.0;
    }
    NumericVector sdeng = sdenc;
    // Have a first stab at guessing snow G
    NumericVector G = GFluxCppsnow(tc, sdenc);
    // Have a first stab at guessing diabatic coefficients
    NumericVector psih(tsteps);
    NumericVector psim(tsteps);
    NumericVector phih(tsteps, 1.0);
    // ******************** Initial step ************************* //
    // Initalize variables
    NumericVector Tc(tsteps, tc[1]);
    NumericVector Tg(tsteps);
    NumericVector sdepc(tsteps + 1, isnowd);
    NumericVector sdepg(tsteps + 1, 0.5 * isnowd);
    NumericVector RswabsG(tsteps);
    NumericVector RlwabsG(tsteps);
    NumericVector tr(tsteps);
    NumericVector umu(tsteps);
    NumericVector sublmelt(tsteps);
    NumericVector rainmelt(tsteps);
    NumericVector tempmelt(tsteps);
    NumericVector sstemp(tsteps);
    double tst = 100.0;
    // **** Iterate until convergene
    int iter = 0;
    double mxdif = 0.0;
    // Initalize variables that are passed to snowoneB
    obspoint obstimeo;
    climpoint climo;
    vegpoint vegpo; vegpo.pai = vegp[0]; vegpo.hgt = vegp[1];  vegpo.clump = vegp[3]; vegpo.ltra = vegp[2];
    otherpoint othero; othero.slope = slope; othero.aspect = aspect; othero.lat = lat; othero.lon = lon; othero.zref = zref;
    snowpoint snowo;
    while (tst > tol) {
        // **** Iterate through all time steps
        int snowagec = isnowa;
        int snowageg = isnowa;
        // **** Extract antecident Tc and Tg
        NumericVector Tco = Tc;
        NumericVector Tgo = Tg;
        mxdif = 0.0;
        for (int i = 0; i < tsteps; ++i) {
            // Get model inputs
            obstimeo.year = year[i]; obstimeo.month = month[i]; obstimeo.day = day[i]; obstimeo.hour = hour[i];
            climo.tc = tc[i]; climo.ea = ea[i]; climo.pk = pk[i]; climo.u2 = u2[i]; climo.prec = prec[i];
            climo.Rsw = Rsw[i]; climo.Rdif = Rdif[i]; climo.Rlw = Rlw[i]; climo.Tci = Tc[i]; climo.te = te[i];
            othero.psim = psim[i]; othero.psih = psih[i]; othero.G = G[i];
            snowo.alb = salb[i]; snowo.sdenc = sdenc[i]; snowo.sdeng = sdeng[i]; snowo.sdepc = sdepc[i]; snowo.sdepg = sdepg[i];
            snowo.snowagec = snowagec; snowo.snowageg = snowageg;
            // Run model
            snowmodpoint smod = snowoneB(obstimeo, climo, vegpo, snowo, othero, sdp, 1.0);
            Tc[i] = smod.Tc;
            Tg[i] = smod.Tg;
            // Re-assign snow density and age
            snowagec = smod.snowagec;
            snowageg = smod.snowageg;
            sdepc[i + 1] = smod.sdepc;
            sdepg[i + 1] = smod.sdepg;
            // Combine the data
            Tc[i] = 0.5 * Tco[i] + 0.5 * Tc[i];
            Tg[i] = 0.5 * Tgo[i] + 0.5 * Tg[i];
            double abs1 = std::abs(Tc[i] - Tco[i]);
            double abs2 = std::abs(Tg[i] - Tgo[i]);
            if (mxdif < abs1) mxdif = abs1;
            if (mxdif < abs2) mxdif = abs2;
            // Recalculate H
            double cp = cpairCpp(tc[i]);
            double ph = phairCpp(tc[i], pk[i]);
            H[i] = cp * smod.gHa * (Tc[i] - tc[i]);
            // Recalculate diabatic coefficients
            double d = zeroplanedisCpp(smod.hgt, smod.pai);
            double zm = roughlengthCpp(smod.hgt, smod.pai, d, psih[i]);
            if (zm < 0.001) zm = 0.001;
            double Tk = tc[i] + 273.15;
            double LL = (ph * cp * std::pow(smod.uf, 3.0) * Tk) / (-ka * 9.81 * H[i]);
            psim[i] = dpsimCpp(zm / LL) - dpsimCpp((zref - d) / LL);
            psih[i] = dpsihCpp((0.2 * zm) / LL) - dpsihCpp((zref - d) / LL);
            phih[i] = dphihCpp((zref - d) / LL);
            // Set limits to diabatic coefficients
            double Belim = 0.4 / std::sqrt(0.003 + (0.2 * smod.pai) / 2.0);
            double ln1 = std::log((zref - d) / zm);
            double ln2 = std::log((zref - d) / (0.2 * zm));
            if (psim[i] < -0.9 * ln1) psim[i] = -0.9 * ln1;
            if (psih[i] < -0.9 * ln2) psih[i] = -0.9 * ln2;
            if (psim[i] > 0.9 * ln1) psim[i] = 0.9 * ln1;
            if (psih[i] > 0.9 * ln2) psih[i] = 0.9 * ln2;
            if (psih[i] > 0.9 * Belim) psih[i] = 0.9 * Belim;
            RswabsG[i] = smod.RswabsG;
            RlwabsG[i] = smod.RlwabsG;
            tr[i] = smod.tr;
            double ufps = (0.4 * u2[i]) / std::log((zref - d) / zm);
            umu[i] = smod.uf / ufps;
            te[i] = (Tc[i] + tc[i]) / 2.0;
            // Add melt
            sublmelt[i] = smod.mSc;
            tempmelt[i] = smod.mMc;
            rainmelt[i] = smod.mRc;
            sstemp[i] = smod.Tcp;
        }
        // Recalculate G
        G = GFluxCppsnow(Tg, sdenc);
        tst = mxdif;
        ++iter;
        if (iter > maxiter) tst = 0;
    }
    List out;
    out["Tc"] = Tc;
    out["Tg"] = Tg;
    out["sdepc"] = sdepc;
    out["sdepg"] = sdepg;
    out["sdenc"] = sdenc;
    out["sdeng"] = sdeng;
    out["G"] = G;
    out["RswabsG"] = RswabsG;
    out["RlwabsG"] = RlwabsG;
    out["tr"] = tr;
    out["umu"] = umu;
    out["mxdif"] = mxdif;
    out["sublmelt"] = sublmelt;
    out["tempmelt"] = tempmelt;
    out["rainmelt"] = rainmelt;
    out["sstemp"] = sstemp;
    return out;
}
// snow energy and mass balance model - data.frame climate
// [[Rcpp::export]]
List gridmodelsnow1(DataFrame obstime, DataFrame climdata, DataFrame pointm, List vegp,
    List other, std::string snowenv)
{
    // EXtract obstime
    IntegerVector year = obstime["year"];
    IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];
    NumericVector hour = obstime["hour"];
    // Extract pointm
    NumericVector Gp = pointm["Gp"];
    NumericVector Tcp = pointm["Tc"];
    NumericVector RswabsG = pointm["RswabsG"];
    NumericVector RlwabsG = pointm["RlwabsG"];
    NumericVector umu = pointm["umu"];
    NumericVector tr = pointm["tr"];
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
    // Extract climate data
    NumericVector tc = climdata["temp"];
    NumericVector rh = climdata["relhum"];
    NumericVector pk = climdata["pres"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector Rlw = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"];
    NumericVector wdir = climdata["winddir"];
    NumericVector prec = climdata["precip"];
    NumericVector salb = snowalbCpp(prec);
    // set dims
    int rows = pai.nrow();
    int cols = pai.ncol();
    int tsteps = tc.size();
    IntegerVector dim = { rows, cols, tsteps };
    // Calculate extras
    NumericVector ea(tsteps);
    NumericVector te(tsteps);
    NumericVector Rnet(tsteps);
    for (int i = 0; i < tsteps; ++i) {
        ea[i] = satvapCpp(tc[i]) * rh[i] / 100.0;
        te[i] = (Tcp[i] + tc[i]) / 2.0;
        double Rem = 0.97 * sb * radem(tc[i]);
        Rnet[i] = RswabsG[i] + RlwabsG[i] - Rem;
    }
    int ndays = tr.size() / 24;
    NumericVector Rmxd(ndays, -1352.0);
    NumericVector Rmnd(ndays, 1352.0);
    NumericVector Rswmnd(ndays);
    NumericVector Rlwmnd(ndays);
    NumericVector Rswmxd(ndays);
    NumericVector Rlwmxd(ndays);
    for (int d = 0; d < ndays; ++d) {
        for (int h = 0; h < 24; ++h) {
            int idx = d * 24 + h;
            if (Rmxd[d] < Rnet[idx]) {
                Rmxd[d] = Rnet[idx];
                Rswmxd[d] = Rsw[idx];
                Rlwmxd[d] = Rlw[idx];
            }
            if (Rmnd[d] > Rnet[idx]) {
                Rmnd[d] = Rnet[idx];
                Rswmnd[d] = Rsw[idx];
                Rlwmnd[d] = Rlw[idx];
            }
        }
    }
    NumericVector Rmx(tsteps);
    NumericVector Rmn(tsteps);
    NumericVector Rswmn(tsteps);
    NumericVector Rlwmn(tsteps);
    NumericVector Rswmx(tsteps);
    NumericVector Rlwmx(tsteps);
    for (int d = 0; d < ndays; ++d) {
        for (int h = 0; h < 24; ++h) {
            int idx = d * 24 + h;
            Rmx[idx] = Rmxd[d];
            Rmn[idx] = Rmnd[d];
            Rswmn[idx] = Rswmnd[d];
            Rlwmn[idx] = Rlwmnd[d];
            Rswmx[idx] = Rswmxd[d];
            Rlwmx[idx] = Rlwmxd[d];
        }
    }
    // Calculate Gmx
    NumericVector Gmxd(ndays, 0.0);
    NumericVector Gmx(tc.size());
    for (int d = 0; d < ndays; ++d) {
        for (int h = 0; h < 24; ++h) {
            int idx = d * 24 + h;
            if (std::abs(Rnet[idx]) > Gmxd[d]) Gmxd[d] = std::abs(Rnet[idx]);
        }
        for (int h = 0; h < 24; ++h) {
            int idx = d * 24 + h;
            Gmx[idx] = Gmxd[d];
        }
    }
    // Calculate location invarient variables
    NumericVector zend(tsteps);
    NumericVector azid(tsteps);
    IntegerVector sindex(tsteps);
    IntegerVector windex(tsteps);
    for (int i = 0; i < tsteps; ++i) {
        solmodel sp = solpositionCpp(lat, lon, year[i], month[i], day[i], hour[i]);
        zend[i] = sp.zend;
        azid[i] = sp.azid;
        sindex[i] = static_cast<int>(std::round(azid[i] / 15)) % 24;
        windex[i] = static_cast<int>(std::round(wdir[i] / 45)) % 8;
    }
    // Create and shape variables
    NumericVector Tc(tsteps * rows * cols, NA_REAL);
    NumericVector Tg(tsteps * rows * cols, NA_REAL);
    NumericVector sdepg(tsteps * rows * cols, NA_REAL);
    NumericVector sdepc(tsteps * rows * cols, NA_REAL);
    NumericVector sdenc(tsteps * rows * cols, NA_REAL);
    Tc.attr("dim") = dim;
    Tg.attr("dim") = dim;
    sdepg.attr("dim") = dim;
    sdepc.attr("dim") = dim;
    sdenc.attr("dim") = dim;
    NumericMatrix ageg = bioclimfill(rows, cols);
    NumericMatrix agec = bioclimfill(rows, cols);
    NumericMatrix meltg = bioclimfill(rows, cols);
    NumericMatrix meltc = bioclimfill(rows, cols);
    // Initalize variables that are passed to snowoneB
    obspoint obstimeo;
    climpoint climo;
    vegpoint vegpo;
    otherpoint othero; othero.lat = lat; othero.lon = lon; othero.zref = zref;  othero.psim = 0.0; othero.psih = 0.0;
    snowpoint snowo;
    std::vector<double> sdp = snowdenp(snowenv);  /// Snow density parameters
    // **** Loop through each grid cell and time step
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!Rcpp::NumericMatrix::is_na(val)) {
                vegpo.pai = pai(i, j); vegpo.hgt = hgt(i, j);  vegpo.clump = clump(i, j); vegpo.ltra = ltra(i, j);
                othero.slope = slope(i, j); othero.aspect = aspect(i, j);
                // initialize snow density, age and depth
                int snowagec = isnowac(i, j);
                int snowageg = isnowag(i, j);
                double sdencp = ((sdp[0] - sdp[1]) * (1 - std::exp(-sdp[2] * isnowdc(i, j) / 100.0 -
                    sdp[3] * snowagec / 24.0)) + sdp[1]) * 1000.0;
                double sdengp = ((sdp[0] - sdp[1]) * (1 - std::exp(-sdp[2] * isnowdg(i, j) * 0.5 / 100.0 -
                    sdp[3] * snowageg / 24.0)) + sdp[1]) * 1000.0;
                double sdepcp = isnowdc(i, j);
                double sdepgp = isnowdg(i, j);
                meltc(i, j) = 0.0;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k; 
                    int snowtest = 0; // whether to run snow model
                    if (sdepcp > 0.0) snowtest = 1;
                    if (tc[k] < 2.0 && prec[k] > 0.0) snowtest = 1;
                    if (snowtest > 0) {
                        // Adjust pai for presence of snow
                        double paip = pai(i, j);
                        if (hgt(i, j) > sdepgp) {
                            paip = paip * (hgt(i, j) - sdepgp) / hgt(i, j);
                        }
                        // Calculate Ground heat flux
                        double dtR = Rmx[k] - Rmn[k];
                        double trS = skyview(i, j) * std::exp(-paip);
                        double Rem = 0.97 * sb * radem(tc[k]);
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
                        if (ha > std::tan(sa * torad)) smu = 0.0;
                        // Calculate wind shelter
                        double ws = wsa[windex[k] * rows * cols + j * rows + i];
                        // Get model inputs
                        double u2p = umu[k] * ws;
                        double Rdifp = Rdif[k] * skyview(i, j);
                        double Rdirp = (Rsw[k] - Rdif[k]) * smu;
                        double Rswp = Rdirp + Rdifp;
                        double Rlwp = Rlw[k] * skyview(i, j);
                        // Create model inputs
                        obstimeo.year = year[k]; obstimeo.month = month[k]; obstimeo.day = day[k]; obstimeo.hour = hour[k];
                        climo.tc = tc[k]; climo.ea = ea[k]; climo.pk = pk[k]; climo.u2 = u2p; climo.prec = prec[k];
                        climo.Rsw = Rswp; climo.Rdif = Rdifp; climo.Rlw = Rlwp; climo.Tci = Tcp[k]; climo.te = te[k];
                        othero.G = G;
                        snowo.alb = salb[k]; snowo.sdenc = sdencp; snowo.sdeng = sdengp;
                        snowo.sdepc = sdepcp; snowo.sdepg = sdepgp;
                        snowo.snowagec = static_cast<double>(snowagec); snowo.snowageg = static_cast<double>(snowageg);
                        // Run snow model for one grid cell
                        snowmodpoint smod = snowoneB(obstimeo, climo, vegpo, snowo, othero, sdp, 1.0);
                        // Assign variables
                        Tc[idx] = smod.Tc;
                        Tg[idx] = smod.Tg;
                        sdepc[idx] = smod.sdepc;
                        sdepg[idx] = smod.sdepg;
                        sdenc[idx] = smod.sdenc;
                        // reassign variables here
                        sdencp = smod.sdenc;
                        sdengp = smod.sdeng;
                        sdepcp = smod.sdepc;
                        sdepgp = smod.sdepg;
                        snowagec = smod.snowagec;
                        snowageg = smod.snowageg;
                        // Calculate cumulative snow melt in m
                        double melc = smod.mSc + smod.mMc + smod.mRc;
                        double melg = smod.mSg + smod.mMg + smod.mRg;
                        meltc(i, j) = meltc(i, j) + melc;
                        meltc(i, j) = meltc(i, j) + (melc * 1000.0) / smod.sdenc;
                        meltg(i, j) = meltg(i, j) + (melg * 1000.0) / smod.sdeng;
                    } // end snowtest
                    else {
                        Tc[idx] = 0.0;
                        Tg[idx] = 0.0;
                        sdepc[idx] = 0.0;
                        sdepg[idx] = 0.0;
                        sdenc[idx] = sdp[1] * 1000.0;
                    }
                } // end k
                agec(i, j) = snowagec;
                ageg(i, j) = snowageg;
            } // end NA test
        } // end col
    } // end row
    // return output
    List out;
    out["Tc"] = Tc;
    out["Tg"] = Tg;
    out["sdepc"] = sdepc;
    out["sdepg"] = sdepg;
    out["sden"] = sdenc;
    out["agec"] = agec;
    out["ageg"] = ageg;
    out["meltc"] = meltc;
    out["meltg"] = meltg;
    return out;
}
// snow model - array climate
// [[Rcpp::export]]
List gridmodelsnow2(DataFrame obstime, List climdata, List pointm, List vegp,
    List other, std::string snowenv)
{
    // Extract obstime
    IntegerVector year = obstime["year"];
    IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];
    NumericVector hour = obstime["hour"];
    // Extract variables: vegp
    NumericMatrix pai = vegp["pai"];
    NumericMatrix hgt = vegp["hgt"];
    NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"];
    // set dims
    int rows = pai.nrow();
    int cols = pai.ncol();
    int tsteps = year.size();
    IntegerVector dim = { rows, cols, tsteps };
    // Extract pointm
    NumericVector Gp = pointm["Gp"];
    NumericVector Tcp = pointm["Tc"];
    NumericVector RswabsG = pointm["RswabsG"];
    NumericVector RlwabsG = pointm["RlwabsG"];
    NumericVector umu = pointm["umu"];
    NumericVector tr = pointm["tr"];
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
    NumericVector tc = climdata["temp"];
    NumericVector rh = climdata["relhum"];
    NumericVector pk = climdata["pres"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector Rlw = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"];
    NumericVector wdir = climdata["winddir"];
    NumericVector prec = climdata["precip"];
    // Compute wind index
    IntegerVector windex(tsteps);
    for (int i = 0; i < tsteps; ++i) windex[i] = static_cast<int>(std::round(wdir[i] / 45)) % 8;
    // Initialize variables
    NumericVector Tc(tsteps * rows * cols, NA_REAL);
    NumericVector Tg(tsteps * rows * cols, NA_REAL);
    NumericVector sdepg(tsteps * rows * cols, NA_REAL);
    NumericVector sdepc(tsteps * rows * cols, NA_REAL);
    NumericVector sdenc(tsteps * rows * cols, NA_REAL);
    Tc.attr("dim") = dim;
    Tg.attr("dim") = dim;
    sdepg.attr("dim") = dim;
    sdepc.attr("dim") = dim;
    sdenc.attr("dim") = dim;
    NumericMatrix ageg = bioclimfill(rows, cols);
    NumericMatrix agec = bioclimfill(rows, cols);
    NumericMatrix meltg = bioclimfill(rows, cols);
    NumericMatrix meltc = bioclimfill(rows, cols);
    // Initalize variables that are passed to snowoneB
    obspoint obstimeo;
    climpoint climo;
    vegpoint vegpo;
    otherpoint othero; othero.zref = zref;  othero.psim = 0.0; othero.psih = 0.0;
    snowpoint snowo;
    std::vector<double> sdp = snowdenp(snowenv);  /// Snow density parameters
    // **** Loop through each grid cell and time step
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!Rcpp::NumericMatrix::is_na(val)) {
                // Calculate ground heat flux 
                NumericVector Rnet(tsteps);
                NumericVector precv(tsteps);
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    double Rem = 0.97 * sb * radem(tc[idx]);
                    Rnet[k] = RswabsG[idx] + RlwabsG[idx] - Rem;
                    precv[k] = prec[idx];
                }
                NumericVector salb = snowalbCpp(precv);
                int ndays = tsteps / 24;
                NumericVector Rmxd(ndays, -1352.0);
                NumericVector Rmnd(ndays, 1352.0);
                NumericVector Rswmnd(ndays);
                NumericVector Rlwmnd(ndays);
                NumericVector Rswmxd(ndays);
                NumericVector Rlwmxd(ndays);
                for (int d = 0; d < ndays; ++d) {
                    for (int h = 0; h < 24; ++h) {
                        int k = d * 24 + h;
                        int idx = i + rows * j + cols * rows * k;
                        if (Rmxd[d] < Rnet[k]) {
                            Rmxd[d] = Rnet[k];
                            Rswmxd[d] = Rsw[idx];
                            Rlwmxd[d] = Rlw[idx];
                        }
                        if (Rmnd[d] > Rnet[k]) {
                            Rmnd[d] = Rnet[k];
                            Rswmnd[d] = Rsw[idx];
                            Rlwmnd[d] = Rlw[idx];
                        }
                    }
                }
                NumericVector Rmx(tsteps);
                NumericVector Rmn(tsteps);
                NumericVector Rswmn(tsteps);
                NumericVector Rlwmn(tsteps);
                NumericVector Rswmx(tsteps);
                NumericVector Rlwmx(tsteps);
                for (int d = 0; d < ndays; ++d) {
                    for (int h = 0; h < 24; ++h) {
                        int k = d * 24 + h;
                        Rmx[k] = Rmxd[d];
                        Rmn[k] = Rmnd[d];
                        Rswmn[k] = Rswmnd[d];
                        Rlwmn[k] = Rlwmnd[d];
                        Rswmx[k] = Rswmxd[d];
                        Rlwmx[k] = Rlwmxd[d];
                    }
                }
                // Calculate Gmx
                NumericVector Gmxd(ndays, 0.0);
                NumericVector Gmx(tsteps);
                for (int d = 0; d < ndays; ++d) {
                    for (int h = 0; h < 24; ++h) {
                        int k = d * 24 + h;
                        if (std::abs(Rnet[k]) > Gmxd[d]) Gmxd[d] = std::abs(Rnet[k]);
                    }
                    for (int h = 0; h < 24; ++h) {
                        int k = d * 24 + h;
                        Gmx[k] = Gmxd[d];
                    }
                }
                // create vegetation variables
                vegpo.pai = pai(i, j); vegpo.hgt = hgt(i, j);  vegpo.clump = clump(i, j); vegpo.ltra = ltra(i, j);
                othero.slope = slope(i, j); othero.aspect = aspect(i, j); othero.lat = lats(i, j); othero.lon = lons(i, j);
                // initialize snow density, age and depth
                int snowagec = isnowac(i, j);
                int snowageg = isnowag(i, j);
                double sdencp = ((sdp[0] - sdp[1]) * (1 - std::exp(-sdp[2] * isnowdc(i, j) / 100.0 -
                    sdp[3] * snowagec / 24.0)) + sdp[1]) * 1000.0;
                double sdengp = ((sdp[0] - sdp[1]) * (1 - std::exp(-sdp[2] * isnowdg(i, j) * 0.5 / 100.0 -
                    sdp[3] * snowageg / 24.0)) + sdp[1]) * 1000.0;
                double sdepcp = isnowdc(i, j);
                double sdepgp = isnowdg(i, j);
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    int snowtest = 0; // whether to run snow model
                    if (sdepcp > 0.0) snowtest = 1;
                    if (tc[idx] < 2.0 && prec[idx] > 0.0) snowtest = 1;
                    if (snowtest > 0) {
                        // Adjust pai for presence of snow
                        double paip = pai(i, j);
                        if (hgt(i, j) > sdepgp) {
                            paip = paip * (hgt(i, j) - sdepgp) / hgt(i, j);
                        }
                        // Calculate Ground heat flux
                        double dtR = Rmx[k] - Rmn[k];
                        double trS = skyview(i, j) * std::exp(-paip);
                        double Rem = 0.97 * sb * radem(tc[idx]);
                        double dmxS = trS * Rswmx[k] + trS * Rlwmx[k] + (1 - trS) * Rem - Rem;
                        double dmnS = trS * Rswmn[k] + trS * Rlwmn[k] + (1 - trS) * Rem - Rem;
                        double Gmu = (dmxS - dmnS) / dtR;
                        double G = Gp[idx] * Gmu;
                        if (G > Gmx[k]) G = Gmx[k];
                        if (G < -Gmx[k]) G = -Gmx[k];
                        // Calculate solar multiplier
                        solmodel solp = solpositionCpp(lats(i, j), lons(i, j), year[k], month[k], day[k], hour[k]);
                        int sindex = static_cast<int>(std::round(solp.azid / 15.0)) % 24;
                        double ha = hor[sindex * rows * cols + j * rows + i];
                        double sa = pi/2.0 - solp.zenr;
                        double smu = 1.0;
                        if (ha > std::tan(sa)) smu = 0.0;
                        // Calculate wind shelter
                        double ws = wsa[windex[k] * rows * cols + j * rows + i];
                        // Get model inputs
                        double u2p = umu[idx] * ws;
                        double Rdifp = Rdif[idx] * skyview(i, j);
                        double Rdirp = (Rsw[idx] - Rdif[idx]) * smu;
                        double Rswp = Rdirp + Rdifp;
                        double Rlwp = Rlw[idx] * skyview(i, j);
                        // Calculate additional variables
                        double ea = satvapCpp(tc[idx]) * rh[idx] / 100.0;
                        double te = (Tcp[idx] + tc[idx]) / 2.0;
                        // Create model inputs
                        obstimeo.year = year[k]; obstimeo.month = month[k]; obstimeo.day = day[k]; obstimeo.hour = hour[k];
                        climo.tc = tc[idx]; climo.ea = ea; climo.pk = pk[idx]; climo.u2 = u2p; climo.prec = prec[idx];
                        climo.Rsw = Rswp; climo.Rdif = Rdifp; climo.Rlw = Rlwp; climo.Tci = Tcp[idx]; climo.te = te;
                        othero.G = G;
                        snowo.alb = salb[k]; snowo.sdenc = sdencp; snowo.sdeng = sdengp;
                        snowo.sdepc = sdepcp; snowo.sdepg = sdepgp;
                        snowo.snowagec = snowagec; snowo.snowageg = snowageg;
                        // Run snow model for one grid cell
                        snowmodpoint smod = snowoneB(obstimeo, climo, vegpo, snowo, othero, sdp, 1.0);
                        // Assign variables
                        Tc[idx] = smod.Tc;
                        Tg[idx] = smod.Tg;
                        sdepc[idx] = smod.sdepc;
                        sdepg[idx] = smod.sdepg;
                        sdenc[idx] = smod.sdenc;
                        // reassign variables here
                        sdencp = smod.sdenc;
                        sdengp = smod.sdeng;
                        sdepcp = smod.sdepc;
                        sdepgp = smod.sdepg;
                        snowagec = smod.snowagec;
                        snowageg = smod.snowageg;
                        // Calculate cumulative snow melt in m
                        double melc = smod.mSc + smod.mMc + smod.mRc;
                        double melg = smod.mSg + smod.mMg + smod.mRg;
                        meltc(i, j) = meltc(i, j) + melc;
                        meltc(i, j) = meltc(i, j) + (melc * 1000.0) / smod.sdenc;
                        meltg(i, j) = meltg(i, j) + (melg * 1000.0) / smod.sdeng;
                    } // end snowtest
                    else {
                        Tc[idx] = 0.0;
                        Tg[idx] = 0.0;
                        sdepc[idx] = 0.0;
                        sdepg[idx] = 0.0;
                        sdenc[idx] = sdp[1] * 1000.0;
                    }
                } // end k
                agec(i, j) = snowagec;
                ageg(i, j) = snowageg;
            } // end NA test
        } // end col
    } // end row
    // return output
    List out;
    out["Tc"] = Tc;
    out["Tg"] = Tg;
    out["sdepc"] = sdepc;
    out["sdepg"] = sdepg;
    out["sden"] = sdenc;
    out["agec"] = agec;
    out["ageg"] = ageg;
    out["meltc"] = meltc;
    out["meltg"] = meltg;
    return out;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************************ Snow microclimate model ************************************************* //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Calculate daily mean snow temperature
// [[Rcpp::export]]
NumericVector snowdayan(NumericVector stempg)
{
    // Get dimensions of arrays
    IntegerVector dims = stempg.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    int ndays = tsteps / 24;
    // Create output variables
    NumericVector snowtempd(stempg.size(), NA_REAL);
    snowtempd.attr("dim") = dims;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = stempg[i + rows * j];
            if (!NumericVector::is_na(val)) {
                for (int d = 0; d < ndays; ++d) {
                    double sumd = 0.0;
                    for (int h = 0; h < 24; ++h) {
                        int k = d * 24 + h;
                        int idx = i + rows * j + cols * rows * k;
                        sumd += stempg[idx];
                    }
                    double meand = sumd / 24.0;
                    for (int h = 0; h < 24; ++h) {
                        int k = d * 24 + h;
                        int idx = i + rows * j + cols * rows * k;
                        snowtempd[idx] = meand;
                    } // end hour
                } // end day
            } // end NA check
        } // end j
    } // end i
    return snowtempd;
}
NumericMatrix meanDsnow(NumericVector snowden)
{
    IntegerVector dims = snowden.attr("dim");
    int rows = dims[0];
    int cols = dims[1];
    int tsteps = dims[2];
    // repermutate variables
    NumericMatrix meanD = bioclimfill(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = snowden[i + rows * j];
            if (!NumericVector::is_na(val)) {
                double sumD = 0.0;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    double co = 0.0442 * std::exp(5.181 * snowden[idx] / 1000.0);
                    double kap = co / (snowden[idx] * 2090.0);
                    sumD += std::sqrt(2.0 * kap / omdy);
                } // end k
                meanD(i, j) = sumD / static_cast<double>(tsteps);
            } // end NA check
        } // end j
    } // end i
    return meanD;
}
// Run full model as point (snow)
snowmicro snowabovepoint(double reqhgt, double zref, double tc, double relhum, double pk, double u2, 
    double Rsw, double Rdif, double Rlw, 
    double hgt, double pai, double paia, double leafd, double clump, double ltra, double leafden,
    solmodel solp, double si, double svfa, int shadowmask, double ws, double umu, double mxtc, snowpoint2 snowp)
{
    if (reqhgt == 0.0) reqhgt = 0.001;
    snowmicro out;
    // Calculate base variables
    double es = satvapCpp(tc);
    double ea = es * relhum / 100.0;
    double tdew = dewpointCpp(ea);
    // Adjust veg hgt and pai for presence of snow
    double hgts = hgt - snowp.sdepg;
    if (hgts < 0.0) hgts = 0.0;
    double pais = 0.0;
    tiwstruct tiw;
    if (hgts > 0.0) {
        pais = pai * hgts / hgt;
        tiw = windtiCpp(hgts, pais);
    }
    else {
        tiw.d = 0.0;
        tiw.zm = 1e-5;
    }
    // Calculate wind speed
    windmodel wind = windCpp(reqhgt, zref, hgts, pais, u2, umu, ws, tiw);
    out.uz = wind.uz;
    // Above canopy
    double ez;
    if (reqhgt >= hgts) {
        if (Rsw > 0.0) {
            out.Rddown = Rdif * svfa;
            if (si > 0.0) {
                if (shadowmask > 0) {
                    out.Rbdown = (Rsw - Rdif) / si;
                    if (out.Rbdown > 1352.0) out.Rbdown = 1352.0;
                    out.Rdup = snowp.albc * Rsw * svfa;
                }
                else {
                    out.Rbdown = 0.0;
                    out.Rdup = snowp.albc * Rdif * svfa;
                }
            }
            else {
                out.Rbdown = 0.0;
                out.Rdup = snowp.albc * Rdif * svfa;
            }
        }
        else {
            out.Rbdown = 0.0;
            out.Rddown = 0.0;
            out.Rdup = 0.0;
        }
        out.Rlwdn = svfa * Rlw;
        out.Rlwup = svfa * 0.97 * sb * radem(snowp.snowtempc);
        abovecanstruct tv = TVabove(reqhgt, zref, hgts, tiw.d, tiw.zm, snowp.snowtempc, tc, ea, 1.0);
        out.Tz = tv.Tz;
        out.tleaf = snowp.snowtempc;
        ez = tv.ez;
    }
    else {
        double paias = 0.0;
        if (hgts > 0.0) {
            paias = paia * hgts / hgt;
        }
        double zi = 0.0;
        if (snowp.sdepg > 0.0) zi = ((snowp.sdepc - snowp.sdepg) * snowp.sdenc) / (hgts * 1000.0);
        double ltras = ltra * std::exp(-10.1 * zi);
        double sm = ltras + snowp.albc;
        if (sm > 0.999) ltras = 0.999 - snowp.albc;
        double clumps = clump;
        if (clump > 0.0) clumps = std::pow(clump, pais / pai);
        double pait = pais;
        if (clump > 0.0) pait = pais / (1.0 - clumps);
        // Calculate two-stream parameters
        tsdifstruct tspdif = twostreamdifCpp(pait, 1.0, snowp.albc, ltras, snowp.albg);
        tirstruct tir = twostreamdif(pais, paias, 1.0, snowp.albc, ltras, clumps, snowp.albg);
        kstruct kp = cankCpp(solp.zenr, 1.0, si);
        tsdirstruct tspdir = twostreamdirCpp(pait, tspdif.om, tspdif.a, tspdif.gma, tspdif.J, tspdif.del, tspdif.h, snowp.albg,
            kp.kd, tspdif.u1, tspdif.S1, tspdif.D1, tspdif.D2);
        // Calculate radiation
        radmodel2 rad = twostreamCpp(pais, clumps, snowp.albg, svfa, si, tc, Rsw, Rdif, Rlw, solp, kp, tspdir, tir);
        if (shadowmask == 0) rad.Rbdown = 0.0;
        // Calculate temperature below canopy
        stompstruct stomp;
        leaftempstruct tvl = leaftemp(snowp.snowtempc, snowp.snowtempg, tc, mxtc, pk, ea, es, wind.uz, tdew, 1.0, rad.radLsw,
            rad.Rddown, rad.Rbdown, Rlw, pais, paias, leafd, 999.999, rad.radLpar, 0.4, 0.4, 2.6, 5.2, stomp);
        out.tleaf = tvl.tleaf;
        double H = 29.3 * wind.gHa * (snowp.snowtempc - tc);
        double Flux = H * (1.0 - std::exp(-pais));
        double Fluxz = tvl.H;
        abovecanstruct tv = TVabove(hgts, zref, hgts, tiw.d, tiw.zm, snowp.snowtempc, tc, ea, 1.0);
        double SH = tv.Tz * 29.3 * 43.0;
        double SG = snowp.snowtempg * 29.3 * 43.0;
        double mxnear = std::abs(out.tleaf - tv.Tz) * 29.3 * 43.0;
        out.Tz = TVbelow(zref, reqhgt, tiw.d, hgts, pais, wind.uf, leafden, Flux, Fluxz, SH, SG, mxnear) / (29.3 * 43);
        // Calculative vapour pressure below canopy
        double la;
        if (tc < 0) {
            la = 51078.69 - 4.338 * tc - 0.06367 * tc * tc;
        }
        else {
            la = 45068.7 - 42.8428 * tc;
        }
        double m = la * (wind.gHa / pk);
        double L = m * (es - ea);
        Flux = L * (1.0 - std::exp(-pais));
        Fluxz = tvl.L;
        double mu = la * (43 / pk);
        SH = tv.ez * mu;
        SG = satvapCpp(snowp.snowtempg) * mu;
        mxnear = std::abs(satvapCpp(out.tleaf) - tv.ez) * mu;
        ez = TVbelow(zref, reqhgt, tiw.d, hgts, pais, wind.uf, leafden, Flux, Fluxz, SH, SG, mxnear) / mu;
        out.Rbdown = rad.Rbdown;
        out.Rddown = rad.Rddown;
        out.Rdup = rad.Rdup;
        out.Rlwdn = tvl.lwdn;
        out.Rlwup = tvl.lwup;
    }
    // Set value limits
    out.rh = (ez / satvapCpp(out.Tz)) * 100.0;
    if (out.rh > 100.0) out.rh = 100.0;
    double tmx = std::max({ out.tleaf, tc, snowp.snowtempg, snowp.snowtempc }) + 2.0;
    double tmn = std::min({ out.tleaf, tc, snowp.snowtempg, snowp.snowtempc }) - 2.0;
    if (out.Tz > tmx) out.Tz = tmx;
    if (out.Tz < tmn) out.Tz = tmn;
    return out;
}
// [[Rcpp::export]]
double belowpointsnow(double reqhgt, double meanD, double snowtempg, double Tzd, double Tza, double hiy) {
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
// snow microclimate model - data.frame climate
// [[Rcpp::export]]
List gridmicrosnow1(double reqhgt, DataFrame obstime, DataFrame climdata, List snowm, List micro, List vegp, List other,
    double mat, std::vector<bool> out) {
    // Extract obstime
    IntegerVector year = obstime["year"];
    IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];
    NumericVector hour = obstime["hour"];
    // Extract climdata
    NumericVector tc = climdata["temp"];
    NumericVector rh = climdata["relhum"];
    NumericVector pk = climdata["pres"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector Rlw = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"];
    NumericVector wdir = climdata["winddir"];
    NumericVector prec = climdata["precip"];
    NumericVector umu = climdata["umu"];
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
    // Extract snow variables (3D)
    NumericVector snowtempc = snowm["Tc"];
    NumericVector snowtempg = snowm["Tg"];
    NumericVector swe = snowm["totalSWE"];
    NumericVector sdepg = snowm["groundsnowdepth"];
    NumericVector sden = snowm["snowden"];
    // Calculate additional snow variables
    NumericVector Tzd = snowdayan(snowm["Tg"]);
    NumericMatrix meanD = meanDsnow(snowm["snowden"]);
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
    if (out[0]) Tz = micro["Tz"];
    if (out[1]) tleaf = micro["tleaf"];
    if (out[2]) relhum = micro["relhum"];
    if (out[3]) soilm = micro["soilm"];
    if (out[4]) uz = micro["windspeed"];
    if (out[5]) Rdirdown = micro["Rdirdown"];
    if (out[6]) Rdifdown = micro["Rdifdown"];
    if (out[7]) Rlwdown = micro["Rlwdown"];
    if (out[8]) Rswup = micro["Rswup"];
    if (out[9]) Rlwup = micro["Rlwup"];
    // Calculate things that don't vary spatially
    NumericVector salb = snowalbCpp(prec);
    double mxtc = -273.15;
    NumericVector zend(tsteps);
    NumericVector azid(tsteps);
    NumericVector zenr(tsteps);
    NumericVector azir(tsteps);
    IntegerVector sindex(tsteps);
    IntegerVector windex(tsteps);
    for (int i = 0; i < tsteps; ++i) {
        if (tc[i] > mxtc) mxtc = tc[i];
        solmodel sp = solpositionCpp(lat, lon, year[i], month[i], day[i], hour[i]);
        zend[i] = sp.zend;
        azid[i] = sp.azid;
        zenr[i] = sp.zenr;
        azir[i] = sp.azir;
        sindex[i] = static_cast<int>(std::round(azid[i] / 15)) % 24;
        windex[i] = static_cast<int>(std::round(wdir[i] / 45)) % 8;
    }
    // Calculate hours in year
    int hiy = (year[0] % 4 == 0 && (year[0] % 100 != 0 || year[0] % 400 == 0)) ? 366 * 24 : 365 * 24;
    // Go through each grid cell
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    // check whether to run model
                    if (swe[idx] > 0.0) {
                        // Calculate adjusted reqhgt
                        double reqhgts = reqhgt - sdepg[idx];
                        if (reqhgts >= 0.0) {
                            // Calculate shadowmask
                            int shadowmask = 1;
                            double ha = hor[sindex[k] * rows * cols + j * rows + i];
                            double sa = (pi / 2.0) - zenr[k];
                            double si = solarindexCpp(slope(i, j), aspect(i, j), zend[k], azid[k], true);
                            if (Rcpp::NumericVector::is_na(si)) si = std::cos(zend[k] * torad);
                            if (ha > std::tan(sa)) shadowmask = 0;
                            // Calculate wind shelter coefficient
                            double ws = wsa[windex[k] * rows * cols + j * rows + i];
                            // call snow model
                            solmodel solp; solp.zend = zend[k]; solp.zenr = zenr[k]; solp.azid = azid[k]; solp.azir = azir[k];
                            double sdepc = swe[idx] / sden[idx];
                            snowpoint2 snowp; snowp.snowtempg = snowtempg[idx]; snowp.snowtempc = snowtempc[idx];
                            snowp.sdepc = sdepc; snowp.sdepg = sdepg[idx]; snowp.sdenc = sden[idx]; 
                            snowp.albc = salb[k]; snowp.albg = salb[k];
                            snowmicro apv = snowabovepoint(reqhgts, zref, tc[k], rh[k], pk[k], u2[k], Rsw[k], Rdif[k], Rlw[k],
                                hgt(i, j), pai(i, j), paia(i, j), leafd(i, j), clump(i, j), ltra(i, j), leafden(i, j),
                                solp, si, skyview(i, j), shadowmask, ws, umu[k], mxtc, snowp);
                            if (out[0]) Tz[idx] = apv.Tz;
                            if (out[1]) tleaf[idx] = apv.tleaf;
                            if (out[2]) relhum[idx] = apv.rh;
                            if (out[4]) uz[idx] = apv.uz;
                            if (out[5]) Rdirdown[idx] = apv.Rbdown;
                            if (out[6]) Rdifdown[idx] = apv.Rddown;
                            if (out[7]) Rlwdown[idx] = apv.Rlwdn;
                            if (out[8]) Rswup[idx] = apv.Rdup;
                            if (out[9]) Rlwup[idx] = apv.Rlwup;
                        } // end above snow
                        else {
                            double bpv = belowpointsnow(reqhgts, meanD(i, j), snowtempg[idx], Tzd[idx], mat, hiy);
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
                } // end k
            } // end NA check
        } // end j
    } // end i
    // return output
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Tz;
    if (out[1]) outp["tleaf"] = tleaf;
    if (out[2]) outp["relhum"] = relhum;
    if (out[3]) outp["soilm"] = soilm;
    if (out[4]) outp["windspeed"] = uz;
    if (out[5]) outp["Rdirdown"] = Rdirdown;
    if (out[6]) outp["Rdifdown"] = Rdifdown;
    if (out[7]) outp["Rlwdown"] = Rlwdown;
    if (out[8]) outp["Rswup"] = Rswup;
    if (out[9]) outp["Rlwup"] = Rlwup;
    return outp;
}
// snow microclimate model -array climate
// [[Rcpp::export]]
List gridmicrosnow2(double reqhgt, DataFrame obstime, List climdata, List snowm, List micro, List vegp, List other,
    double mat, std::vector<bool> out) {
    // Extract obstime
    IntegerVector year = obstime["year"];
    IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];
    NumericVector hour = obstime["hour"];
    // Extract variables: vegp
    NumericMatrix pai = vegp["pai"];
    NumericMatrix paia = vegp["paia"];
    NumericMatrix hgt = vegp["hgt"];
    NumericMatrix ltra = vegp["leaft"];
    NumericMatrix clump = vegp["clump"];
    NumericMatrix leafd = vegp["leafd"];
    NumericMatrix leafden = vegp["leafden"];
    // Extract climdata
    NumericVector tc = climdata["temp"];
    NumericVector rh = climdata["relhum"];
    NumericVector pk = climdata["pres"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector Rlw = climdata["lwdown"];
    NumericVector u2 = climdata["windspeed"];
    NumericVector wdir = climdata["winddir"]; // not 3D
    NumericVector prec = climdata["prec"];
    NumericVector umu = climdata["umu"];
    // Extract other
    NumericMatrix slope = other["slope"];
    NumericMatrix aspect = other["aspect"];
    NumericMatrix skyview = other["skyview"];
    NumericVector wsa = other["wsa"];
    NumericVector hor = other["hor"];
    NumericMatrix lats = other["lat"];
    NumericMatrix lons = other["lon"];
    NumericMatrix Smax = other["Smax"];
    double zref = other["zref"];
    // Extract snow variables
    NumericVector snowtempc = snowm["Tc"];
    NumericVector snowtempg = snowm["Tg"];
    NumericVector swe = snowm["totalSWE"];
    NumericVector sdepg = snowm["groundsnowdepth"];
    NumericVector sden = snowm["snowden"];
    // Calculate snow variables
    NumericVector Tzd = snowdayan(snowm["Tg"]);
    NumericMatrix meanD = meanDsnow(snowm["snowden"]);
    // Extract existing microclimate model variables
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
    if (out[0]) Tz = micro["Tz"];
    if (out[1]) tleaf = micro["tleaf"];
    if (out[2]) relhum = micro["relhum"];
    if (out[3]) soilm = micro["soilm"];
    if (out[4]) uz = micro["windspeed"];
    if (out[5]) Rdirdown = micro["Rdirdown"];
    if (out[6]) Rdifdown = micro["Rdifdown"];
    if (out[7]) Rlwdown = micro["Rlwdown"];
    if (out[8]) Rswup = micro["Rswup"];
    if (out[9]) Rlwup = micro["Rlwup"];
    // set dims
    int rows = pai.nrow();
    int cols = pai.ncol();
    int tsteps = hour.size();
    // Compute wind index
    IntegerVector windex(tsteps);
    for (int i = 0; i < tsteps; ++i) windex[i] = static_cast<int>(std::round(wdir[i] / 45)) % 8;
    // Calculate hours in year
    int hiy = (year[0] % 4 == 0 && (year[0] % 100 != 0 || year[0] % 400 == 0)) ? 366 * 24 : 365 * 24;
    // Loop through to calculate microclimate
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = hgt(i, j);
            if (!NumericMatrix::is_na(val)) {
                double mxtc = -273.15;
                NumericVector precv(tsteps);
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    if (tc[idx] > mxtc) mxtc = tc[idx];
                    precv[k] = prec[idx];
                }
                NumericVector salb = snowalbCpp(precv);
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    // check whether to run model
                    if (swe[idx] > 0.0) {
                        // Calculate adjusted reqhgt
                        double reqhgts = reqhgt - sdepg[idx];
                        if (reqhgts >= 0.0) {
                            solmodel solp = solpositionCpp(lats(i, j), lons(i, j), year[k], month[k], day[k], hour[k]);
                            int sindex = static_cast<int>(std::round(solp.azid / 15)) % 24;
                            // Calculate shadowmask
                            int shadowmask = 1;
                            double ha = hor[sindex * rows * cols + j * rows + i];
                            double sa = (pi / 2.0) - solp.zenr;
                            double si = solarindexCpp(slope(i, j), aspect(i, j), solp.zend, solp.azid, true);
                            if (Rcpp::NumericVector::is_na(si)) si = std::cos(solp.zenr);
                            if (ha > std::tan(sa)) shadowmask = 0;
                            // Calculate wind shelter coefficient
                            double ws = wsa[windex[k] * rows * cols + j * rows + i];
                            // call snow model
                            double sdepc = swe[idx] / sden[idx];
                            snowpoint2 snowp; snowp.snowtempg = snowtempg[idx]; snowp.snowtempc = snowtempc[idx];
                            snowp.sdepc = sdepc; snowp.sdepg = sdepg[idx]; snowp.sdenc = sden[idx];
                            snowp.albc = salb[k]; snowp.albg = salb[k];
                            snowmicro apv = snowabovepoint(reqhgts, zref, tc[idx], rh[idx], pk[idx], u2[idx], Rsw[idx], Rdif[idx], Rlw[idx],
                                hgt(i, j), pai(i, j), paia(i, j), leafd(i, j), clump(i, j), ltra(i, j), leafden(i, j),
                                solp, si, skyview(i, j), shadowmask, ws, umu[idx], mxtc, snowp);
                            if (out[0]) Tz[idx] = apv.Tz;
                            if (out[1]) tleaf[idx] = apv.tleaf;
                            if (out[2]) relhum[idx] = apv.rh;
                            if (out[4]) uz[idx] = apv.uz;
                            if (out[5]) Rdirdown[idx] = apv.Rbdown;
                            if (out[6]) Rdifdown[idx] = apv.Rddown;
                            if (out[7]) Rlwdown[idx] = apv.Rlwdn;
                            if (out[8]) Rswup[idx] = apv.Rdup;
                            if (out[9]) Rlwup[idx] = apv.Rlwup;
                        } // end above snow
                        else {
                            double bpv = belowpointsnow(reqhgts, meanD(i, j), snowtempg[idx], Tzd[idx], mat, hiy);
                            if (out[0]) Tz[idx] = bpv;
                            if (out[1]) tleaf[idx] = bpv;
                            if (out[2]) relhum[idx] = 100.0;
                            if (out[4]) uz[idx] = 0.0;
                            if (out[5]) Rdirdown[idx] = 0.0;
                            if (out[6]) Rdifdown[idx] = 0.0;
                            if (out[7]) Rlwdown[idx] = 0.0;
                            if (out[8]) Rswup[idx] = 0.0;
                            if (out[9]) Rlwup[idx] = 0.0;
                        } // end below snow
                        if (out[3]) soilm[idx] = Smax(i, j);
                    } // end check whether to run model
                } // end k
            } // end NA check
        } // end j
    } // end i
    // return output
    Rcpp::List outp;
    if (out[0]) outp["Tz"] = Tz;
    if (out[1]) outp["tleaf"] = tleaf;
    if (out[2]) outp["relhum"] = relhum;
    if (out[3]) outp["soilm"] = soilm;
    if (out[4]) outp["windspeed"] = uz;
    if (out[5]) outp["Rdirdown"] = Rdirdown;
    if (out[6]) outp["Rdifdown"] = Rdifdown;
    if (out[7]) outp["Rlwdown"] = Rlwdown;
    if (out[8]) outp["Rswup"] = Rswup;
    if (out[9]) outp["Rlwup"] = Rlwup;
    return outp;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************** Functions used by R *************************************************************** //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ** Calculate clear sky radiation - used by test that etc ** //
// [[Rcpp::export]]
NumericVector clearskyradCpp(IntegerVector year, IntegerVector month, IntegerVector day,
    NumericVector lt, double lat, double lon, NumericVector tc, NumericVector rh,
    NumericVector pk)
{
    int n = year.size();
    NumericVector Ic(n);
    for (int i = 0; i < n; ++i) {
        solmodel sp = solpositionCpp(lat, lon, year[i], month[i], day[i], lt[i]);
        if (sp.zend <= 90.0) {
            double m = 35 * std::cos(sp.zenr) * std::pow(1224.0 * std::cos(sp.zenr) * std::cos(sp.zenr) + 1.0, -0.5);
            double TrTpg = 1.021 - 0.084 * std::sqrt(m * 0.00949 * pk[i] + 0.051);
            double xx = std::log(rh[i] / 100.0) + ((17.27 * tc[i]) / (237.3 + tc[i]));
            double Td = (237.3 * xx) / (17.27 - xx);
            double u = std::exp(0.1133 - std::log(3.78) + 0.0393 * Td);
            double Tw = 1 - 0.077 * std::pow(u * m, 0.3);
            double Ta = 0.935 * m;
            double od = TrTpg * Tw * Ta;
            Ic[i] = 1352.778 * std::cos(sp.zenr) * od;
        }
    }
    return Ic;
}
// ** Calculate solar variables as data.frame of vectors
// [[Rcpp::export]]
DataFrame solpositionvCpp(IntegerVector year, IntegerVector month, IntegerVector day,
    NumericVector lt, double lat, double lon, double slope, double aspect)
{
    int n = year.size();
    NumericVector zend(n);
    NumericVector azid(n);
    NumericVector si(n);
    for (int i = 0; i < n; ++i) {
        solmodel sp = solpositionCpp(lat, lon, year[i], month[i], day[i], lt[i]);
        zend[i] = sp.zend;
        azid[i] = sp.azid;
        si[i] = solarindexCpp(slope, aspect, sp.zend, sp.azid, true);
    }
    DataFrame df = DataFrame::create(
        _["zen"] = zend,
        _["azi"] = azid,
        _["si"] = si);
    return df;
}
// Needed to process point model
// [[Rcpp::export]]
DataFrame pointmprocess(DataFrame pointvars, double zref, double h, double pai,
    double rho, double Vm, double Vq, double Mc)
{
    // Extract required variables form point model
    NumericVector u2 = pointvars["windspeed"];
    NumericVector tc = pointvars["tc"];
    NumericVector rh = pointvars["rh"];
    NumericVector pk = pointvars["pk"];
    NumericVector uf = pointvars["uf"];
    NumericVector soilm = pointvars["soilm"];
    NumericVector RabsG = pointvars["RabsG"];
    // Calculate time invariant variables
    double dp = zeroplanedisCpp(h, pai);
    double zmp = roughlengthCpp(h, pai, dp, 0);
    soilstruct sp = soilpfun(Vm, Vq, Mc, rho);
    // Create variables
    int tsteps = u2.size();
    NumericVector umu(tsteps);
    NumericVector kp(tsteps);
    NumericVector DDp(tsteps);
    NumericVector muGp(tsteps);
    NumericVector dtrp(tsteps);
    std::vector<double> T0p(tsteps);
    for (int i = 0; i < tsteps; ++i) {
        // Calculate umu
        double ufps = (ka * u2[i]) / std::log((zref - dp) / zmp);
        umu[i] = uf[i] / ufps;
        // Calculate conductivity
        double cs = (2400 * rho / 2.64 + 4180 * soilm[i]); // specific heat of soil in J / kg / K
        double ph = (rho * (1.0 - soilm[i]) + soilm[i]) * 1000; // bulk density in kg / m3
        double c2 = 1.06 * rho * soilm[i];
        kp[i] = sp.c1 + c2 * soilm[i] - (sp.c1 - sp.c4) * std::exp(-std::pow(sp.c3 * soilm[i], 4.0));
        // Calculate muGp
        double ka = kp[i] / (cs * ph);
        muGp[i] = std::pow(2.0 * ka / omdy, 0.5);
        // Calculate DDp
        DDp[i] = std::pow(2.0 * ka / omdy, 0.5);
        // Calculate T0p
        double gHa = (0.4 * 43.0 * ufps) / std::log((zref - dp) / zmp);
        double es = satvapCpp(tc[i]);
        double ea = es * rh[i] / 100.0;
        T0p[i] = PenmanMonteithCpp(RabsG[i], gHa, gHa, tc[i], tc[i], pk[i], ea, 0.97, 0.0, 1.0);
    }
    // Calculate dtrp
    std::vector<double> T0pmx = hourtodayCpp(T0p, "max");
    std::vector<double> T0pmn = hourtodayCpp(T0p, "min");
    int tsteps2 = T0pmx.size();
    for (int i = 0; i < tsteps2; ++i) {
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
//  ** Functions needed for calculated topographic wetness index
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
// Calculate melt mu - this function derives a temperature melt factor based on sky view 
// for all temperatures that exceed zero
// [[Rcpp::export]]
NumericMatrix meltmu(NumericMatrix skyview, NumericVector stemp, NumericVector tc)
{
    double dhp = 0.0;
    for (R_xlen_t i = 0; i < stemp.size(); ++i) {
        if (stemp[i] > 0.0) dhp += stemp[i];
    }
    // Calculate multiplier for temperature melt
    int rows = skyview.nrow();
    int cols = skyview.ncol();
    NumericMatrix mumelt(rows, cols);
    if (dhp > 0.0) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double val = skyview(i, j);
                if (!Rcpp::NumericMatrix::is_na(val)) {
                    double dhm = 0.0;
                    for (R_xlen_t k = 0; k < stemp.size(); ++k) {
                        double dif = stemp[k] - tc[k];
                        double stemp2 = dif * skyview(i, j) + tc[k];
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
    NumericMatrix mumelt(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = mu(i, j);
            if (!Rcpp::NumericMatrix::is_na(val)) {
                double dhp = 0.0;
                double dhm = 0.0;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + cols * rows * k;
                    if (stemp[idx] > 0.0) dhp += stemp[idx];
                    double dif = stemp[idx] - tc[idx];
                    double stemp2 = dif * mu(i, j) + tc[idx];
                    if (stemp2 > 0.0) dhm += stemp2;
                }
                if (dhp > 0.0) {
                    mumelt(i, j) = dhm / dhp;
                }
                else {
                    mumelt(i, j) = 0.5;
                }
            }
            else {
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
// Apply the function (mean, sum, min, max) over the third dimension
// [[Rcpp::export]]
NumericVector applycpp3(NumericVector a, std::string fun_name) {
    IntegerVector dim = a.attr("dim");
    int rows = dim[0], cols = dim[1], tsteps = dim[2];
    int n = rows * cols;
    int fun;
    if (fun_name == "mean") fun = 0;
    else if (fun_name == "sum") fun = 1;
    else if (fun_name == "max") fun = 2;
    else if (fun_name == "min") fun = 3;
    else stop("Unknown function name");
    NumericVector result(tsteps);
    for (int k = 0; k < tsteps; ++k) {
        int base = k * n;
        double out, cnt = 0;
        if (fun < 2) {
            out = 0.0;
            for (int i = 0; i < n; ++i) {
                double v = a[base + i];
                if (!NumericVector::is_na(v)) { out += v; ++cnt; }
            }
            result[k] = (fun == 0) ? (cnt > 0 ? out / cnt : R_NaN) : out;
        }
        else {
            out = (fun == 2) ? R_NegInf : R_PosInf;
            for (int i = 0; i < n; ++i) {
                double v = a[base + i];
                if (!NumericVector::is_na(v)) {
                    if (fun == 2) { if (v > out) out = v; }
                    else { if (v < out) out = v; }
                }
            }
            result[k] = out;
        }
    }
    return result;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ************************************** Functions used by Testthat ******************************************************** //
// Test wrappers: microclimate model
// [[Rcpp::export]]
DataFrame microclimatemodel_wrapper(DataFrame obstime, DataFrame climdata, List BLout,
    std::vector<double> vegp, std::vector<double> groundp,
    double reqhgt, double zref, double lat, double lon)
{
    // Access items of obstime
    IntegerVector year = obstime["year"];
    IntegerVector month = obstime["month"];
    IntegerVector day = obstime["day"];
    NumericVector hour = obstime["hour"];
    // Access columns of climdata
    NumericVector tc = climdata["temp"];
    NumericVector rh = climdata["relhum"];
    NumericVector pk = climdata["pres"];
    NumericVector Rsw = climdata["swdown"];
    NumericVector Rdif = climdata["difrad"];
    NumericVector Rlw = climdata["lwdown"];
    NumericVector wspeed = climdata["windspeed"];
    NumericVector wdir = climdata["winddir"];
    NumericVector prec = climdata["precip"];
    // Access items of vegp
    double hgt = vegp[0];
    double pai = vegp[1];
    double vegx = vegp[2];
    double clump = vegp[3];
    double lref = vegp[4];
    double ltra = vegp[5];
    double leafd = vegp[6];
    // double em = vegp[7];
    double gsmax = vegp[8];
    // Access items of groundp
    double gref = groundp[0];
    double slope = groundp[1];
    double aspect = groundp[2];
    // double groundem = groundp[3];
    double rho = groundp[4];
    double Vm = groundp[5];
    double Vq = groundp[6];
    double Mc = groundp[7];
    double soilb = groundp[8];
    double psie = groundp[9];
    double Smax = groundp[10];
    double Smin = groundp[11];
    // Compute additional climate variables and solar zenith and si
    int n = tc.size();
    double mxtc = -999.99;
    std::vector<double> soilm(n);
    for (int i = 0; i < n; ++i) {
        soilm[i] = 0.3;
        if (tc[i] > mxtc) mxtc = tc[i];
    }
    // Extract from Bigleaf model
    NumericVector varTg = BLout["Tg"];
    NumericVector G = BLout["G"];
    NumericVector uf = BLout["uf"];
    if (reqhgt >= 0) {
        // Create additional output variables
        NumericVector Rdirdown(n);
        NumericVector Rdifdown(n);
        NumericVector Rlwdown(n);
        NumericVector Rswup(n);
        NumericVector Rlwup(n);
        NumericVector Tz(n);
        NumericVector tleaf(n);
        NumericVector rhz(n);
        NumericVector uz(n);
        // Compute radiation variables
        double paia = 0;
        if (reqhgt < hgt) paia = (1.0 - reqhgt / hgt) * pai;
        // Compute time-invarent valiables
        double pait = pai;
        if (clump > 0.0) pait = pai / (1.0 - clump);
        tsdifstruct tspdif = twostreamdifCpp(pait, vegx, lref, ltra, gref);
        tirstruct tir = twostreamdif(pai, paia, vegx, lref, ltra, clump, gref);
        stompstruct stomp;
        tiwstruct tiw;
        double dp;
        double zmp;
        if (reqhgt > 0.0) {
            stomp = stomparamsCpp(hgt, lat, vegx);
            tiw = windtiCpp(hgt, pai);
            dp = zeroplanedisCpp(hgt, pai);
            zmp = roughlengthCpp(hgt, pai, dp, 0);
        }
        for (int i = 0; i < n; ++i) {
            solmodel solp = solpositionCpp(lat, lon, year[i], month[i], day[i], hour[i]);
            double si = solarindexCpp(slope, aspect, solp.zend, solp.azid);
            kstruct kp = cankCpp(solp.zenr, vegx, si);
            tsdirstruct tspdir = twostreamdirCpp(pait, tspdif.om, tspdif.a, tspdif.gma, tspdif.J, tspdif.del, tspdif.h, gref,
                kp.kd, tspdif.u1, tspdif.S1, tspdif.D1, tspdif.D2);
            radmodel2 radvars = twostreamCpp(pai, clump, gref, 1.0, si, tc[i], Rsw[i], Rdif[i], Rlw[i],
                solp, kp, tspdir, tir);
            // Save radiation variables
            Rdirdown[i] = radvars.Rbdown;
            Rdifdown[i] = radvars.Rddown;
            Rswup[i] = radvars.Rdup;
            if (reqhgt == 0) {
                Rlwdown[i] = radvars.radGlw;
                Rlwup[i] = 0.97 * sb * radem(varTg[i]);
                Tz[i] = varTg[i];
            }
            else {
                // Compute wind
                double ufps = (ka * wspeed[i]) / std::log((zref - dp) / zmp);
                double umu = uf[i] / ufps;
                windmodel wvars = windCpp(reqhgt, zref, hgt, pai, zref, umu, 1.0, tiw);
                uz[i] = wvars.uz;
                // Compute additional variabales needed
                double es = satvapCpp(tc[i]);
                double ea = es * (rh[i] / 100.0);
                double tdew = dewpointCpp(ea);
                double leafden = pai / hgt;
                // Compute ground
                soilmodel Gvars;
                Gvars.Tg = varTg[i];
                Gvars.G = G[i];
                abovemodel out = TVaboveground(reqhgt, zref, tc[i], pk[i], ea, es, tdew, Rsw[i], Rdif[i], Rlw[i], soilm[i], hgt, pai, paia, vegx,
                    leafd, leafden, Smin, Smax, psie, soilb, gsmax, mxtc, stomp, tir, radvars, tiw, wvars, Gvars);
                Tz[i] = out.Tz;
                tleaf[i] = out.tleaf;
                rhz[i] = out.rh;
                Rlwdown[i] = out.lwdn;
                Rlwup[i] = out.lwup;
            } // end above gorund check
        } // end time step
        if (reqhgt == 0.0) {
            DataFrame df = DataFrame::create(
                _["Tz"] = Tz,
                _["Rdirdown"] = Rdirdown,
                _["Rdifdown"] = Rdifdown,
                _["Rswup"] = Rswup,
                _["Rlwdown"] = Rlwdown,
                _["Rlwup"] = Rlwup,
                _["soilm"] = Rcpp::wrap(soilm));
            return df;
        }
        else {
            DataFrame df = DataFrame::create(
                _["Tz"] = Tz,
                _["tleaf"] = tleaf,
                _["rh"] = rhz,
                _["uz"] = uz,
                _["Rdirdown"] = Rdirdown,
                _["Rdifdown"] = Rdifdown,
                _["Rswup"] = Rswup,
                _["Rlwdown"] = Rlwdown,
                _["Rlwup"] = Rlwup,
                _["soilm"] = Rcpp::wrap(soilm));
            return df;
        }
    }
    else { // Below ground
        // Compute meanD
        soilstruct sp = soilpfun(Vm, Vq, Mc, rho);
        double sumD = 0.0;
        for (int i = 0; i < n; ++i) {
            soilkstruct skd = soilcondCpp(rho, soilm[i], sp);
            sumD += skd.DD;
        }
        double meanD = sumD / static_cast<double>(n);
        // Compute temperature below ground
        double nb = -118.35 * reqhgt / meanD;
        int nn = static_cast<int>(std::round(nb));
        std::vector<double> Tg = BLout["Tg"];;
        std::vector<double> Tz = manCpp(Tg, nn);
        DataFrame df = DataFrame::create(
            _["Tz"] = Tz,
            _["soilm"] = Rcpp::wrap(soilm));
        return df;
    }
}
