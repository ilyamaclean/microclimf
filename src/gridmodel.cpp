// gridmodel.cpp
// Landscape-scale microclimate calculations. The grid model uses the
// physically resolved point-model solution as a reference state, then
// represents how terrain, vegetation and soil alter radiation, aerodynamic
// exchange, surface temperature, soil state and within-/above-canopy profiles
// from cell to cell without repeating the full iterative point solve everywhere.
#include "gridmodel.h"
#include "constants.h"
#include "utils.h"
#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using mc::pi;
using mc::sb;
using mc::surfaceEmissivity;

using namespace Rcpp;

namespace gridmodel {

    // ========================================================================
    // Radiation: two-stream shortwave/longwave absorption, ground/canopy
    // totals and the at-height direct/diffuse/reflected streams
    // ========================================================================

    // Precompute the canopy optical state for one grid cell. Total PAI controls
    // whole-canopy transmission and albedo, while PAI above the requested
    // height defines the radiation field seen at that height. Canopy clumping
    // is represented as a mixture of open gaps and the homogeneous two-stream
    // canopy solution.
    tirstruct GridRadswabsSetupCpp(double pai, double paia, double x, double lref, double ltra,
        double clump, double gref)
    {
        tirstruct tir;
        tir.vegetated = (pai > 0.0);
        tir.amx = gref;
        if (tir.amx < lref) tir.amx = lref;
        if (tir.vegetated) {
            // Whole-canopy diffuse transfer and albedo. The effective PAI is
            // increased to describe the vegetation-bearing part of a clumped
            // canopy; explicit gap fractions then restore the open pathways.
            tir.pait = (clump > 0.0) ? pai / (1.0 - clump) : pai;
            tir.tsd = utils::twostreamdifCpp(tir.pait, x, lref, ltra, gref);
            tir.trdn = clump * clump;
            tir.albd = utils::canopyGapMixCpp(tir.trdn * tir.trdn, gref, tir.tsd.p1 + tir.tsd.p2);
            if (tir.albd > tir.amx) tir.albd = tir.amx;
            if (tir.albd < 0.01) tir.albd = 0.01;
            tir.Rddn_g = utils::canopyGapMixCpp(tir.trdn, 1.0,
                tir.tsd.p3 * std::exp(-tir.tsd.h * tir.pait) + tir.tsd.p4 * std::exp(tir.tsd.h * tir.pait));
            if (tir.Rddn_g > 1.0) tir.Rddn_g = 1.0;
            if (tir.Rddn_g < 0.0) tir.Rddn_g = 0.0;

            // Partition the same clumped canopy at the requested height.
            // `gi` and `giu` are gap fractions above and below that height,
            // allowing radiation there to depend on foliage on both sides.
            tir.gi = (clump > 0.0) ? std::pow(clump, paia / pai) : 0.0;
            if (tir.gi > 0.99) tir.gi = 0.99;
            double giu = (clump > 0.0) ? std::pow(clump, (pai - paia) / pai) : 0.0;
            if (giu > 0.99) giu = 0.99;
            double trdz = tir.gi * tir.gi;
            tir.trdu = giu * giu;
            tir.paiaa = paia / (1.0 - tir.gi);

            tir.Rddn_z = utils::canopyGapMixCpp(trdz, 1.0,
                tir.tsd.p3 * std::exp(-tir.tsd.h * tir.paiaa) + tir.tsd.p4 * std::exp(tir.tsd.h * tir.paiaa));
            if (tir.Rddn_z > 1.0) tir.Rddn_z = 1.0;
            if (tir.Rddn_z < 0.0) tir.Rddn_z = 0.0;

            tir.Rdup_z = utils::canopyGapMixCpp(tir.trdu * tir.trdn, gref,
                tir.tsd.p1 * std::exp(-tir.tsd.h * tir.paiaa) + tir.tsd.p2 * std::exp(tir.tsd.h * tir.paiaa));
            if (tir.Rdup_z > 1.0) tir.Rdup_z = 1.0;
            if (tir.Rdup_z < 0.0) tir.Rdup_z = 0.0;
        } else {
            // With no vegetation, radiation at any positive height sees the
            // ground directly: no canopy attenuation and ground reflectance
            // supplies the upward shortwave stream.
            tir.pait = 0.0;
            tir.trdn = 0.0;
            tir.albd = gref;
            tir.Rddn_g = 0.0;
            tir.gi = 0.0;
            tir.trdu = 0.0;
            tir.paiaa = 0.0;
            tir.Rddn_z = 1.0;
            tir.Rdup_z = gref;
        }
        return tir;
    }

    // Apply the timestep's solar geometry and incoming shortwave to the cell's
    // precomputed optical state. Returns absorption by the ground and the
    // combined canopy-ground system, plus direct, diffuse and upward-reflected
    // radiation at the requested height and the shortwave/PAR absorbed by a
    // representative leaf there.
    radmodel2 GridRadswabsStepCpp(const tirstruct& tir, double pai, double clump, double gref,
        double svfa, double si, bool shadow, double zenr, double x, double Rsw, double Rdif)
    {
        radmodel2 out;
        if (Rsw <= 0.0) return out;
        if (zenr > pi / 2.0) zenr = pi / 2.0;
        double cosz = std::cos(zenr);
        // With the sun at or below the horizon there is no direct beam, but the
        // diffuse sky radiation the weather records still arrives.
        const bool beam = !shadow && cosz > 1e-6;

        if (!tir.vegetated) {
            // Bare ground: terrain can remove the direct beam and restrict the
            // diffuse sky, but there is no vegetation to scatter or absorb it.
            out.Rbdown = beam ? std::min((Rsw - Rdif) / cosz, 1352.0) : 0.0;
            out.Rddown = Rdif * svfa;
            out.Rdup = gref * (Rdif * svfa + si * out.Rbdown);
            out.radGsw = (1.0 - gref) * (svfa * Rdif + si * out.Rbdown);
            out.radCsw = out.radGsw;
            return out;
        }

        // The direct beam's contribution to each stream, zero without a beam.
        // Direct-beam extinction changes with solar angle and leaf-angle
        // distribution; the two-stream direct solution describes how that beam
        // is scattered upward/downward as it crosses the canopy.
        double Rbeam = 0.0, Rb = 0.0, kcos = 0.0, albb = 0.0, Rdbdn_g = 0.0, Rbdn_g = 0.0;
        double trgz = 0.0, Rdbdn_z = 0.0, Rdbup_z = 0.0;
        if (beam) {
            utils::kstruct kp = utils::cankCpp(zenr, x, si);
            utils::tsdirstruct tspdir = utils::twostreamdirCpp(tir.pait, tir.tsd.om, tir.tsd.a, tir.tsd.gma,
                tir.tsd.J, tir.tsd.del, tir.tsd.h, gref, kp.kd, tir.tsd.u1, tir.tsd.S1, tir.tsd.D1, tir.tsd.D2);
            Rbeam = (Rsw - Rdif) / cosz;
            if (Rbeam > 1352.0) Rbeam = 1352.0;
            // Beam flux crossing the canopy plane, the quantity every beam
            // fraction from the two-stream solution is expressed against (see
            // the same step in the point model). On level ground si equals
            // cos(zenith).
            Rb = Rbeam * si;
            kcos = kp.k * cosz;

            // Whole-canopy beam transmission, reflection and absorption.
            // Clumping mixes radiation travelling through large gaps with
            // radiation crossing the foliage-bearing portion of the canopy.
            double trbn = std::pow(clump, kp.Kc);
            if (trbn > 0.999) trbn = 0.999;
            if (trbn < 0.0) trbn = 0.0;
            double trg = utils::canopyGapMixCpp(trbn, 1.0, std::exp(-kp.kd * tir.pait));
            albb = utils::canopyGapMixCpp(tir.trdn * trbn, gref,
                tspdir.p5 / -tspdir.sig + tspdir.p6 + tspdir.p7);
            if (albb > tir.amx) albb = tir.amx;
            if (albb < 0.01) albb = 0.01;
            Rdbdn_g = utils::canopyGapMixCpp(trbn, 0.0,
                (tspdir.p8 / tspdir.sig) * std::exp(-kp.kd * tir.pait) +
                tspdir.p9 * std::exp(-tir.tsd.h * tir.pait) + tspdir.p10 * std::exp(tir.tsd.h * tir.pait));
            if (Rdbdn_g > tir.amx) Rdbdn_g = tir.amx;
            if (Rdbdn_g < 0.0) Rdbdn_g = 0.0;
            Rbdn_g = trg;

            // Repeat the beam/scattering calculation only for the canopy
            // above the requested height, then combine it with reflection from
            // below to obtain the local radiation field experienced by a leaf
            // at that height.
            double trbz = std::pow(tir.gi, kp.Kc);
            if (trbz > 0.999) trbz = 0.999;
            if (trbz < 0.0) trbz = 0.0;
            trgz = utils::canopyGapMixCpp(trbz, 1.0, std::exp(-kp.kd * tir.paiaa));
            Rdbdn_z = utils::canopyGapMixCpp(trbz, 0.0,
                (tspdir.p8 / tspdir.sig) * std::exp(-kp.kd * tir.paiaa) +
                tspdir.p9 * std::exp(-tir.tsd.h * tir.paiaa) + tspdir.p10 * std::exp(tir.tsd.h * tir.paiaa));
            if (Rdbdn_z > tir.amx) Rdbdn_z = tir.amx;
            if (Rdbdn_z < 0.0) Rdbdn_z = 0.0;
            Rdbup_z = utils::canopyGapMixCpp(tir.trdu * trbn, gref,
                (tspdir.p5 / -tspdir.sig) * std::exp(-kp.kd * tir.paiaa) +
                tspdir.p6 * std::exp(-tir.tsd.h * tir.paiaa) + tspdir.p7 * std::exp(tir.tsd.h * tir.paiaa));
            if (Rdbup_z > tir.amx) Rdbup_z = tir.amx;
            if (Rdbup_z < 0.0) Rdbup_z = 0.0;
        }

        out.radCsw = (1.0 - tir.albd) * Rdif * svfa + (1.0 - albb) * Rb;
        out.radGsw = (1.0 - gref) * (tir.Rddn_g * Rdif * svfa + Rdbdn_g * Rb + Rbdn_g * Rb);
        double maxg = (1.0 - gref) * (Rdif * svfa + Rb);
        if (out.radGsw > maxg) out.radGsw = maxg;

        // Radiation streams at the requested height. Sunlit leaf absorption
        // additionally includes interception of the direct beam; diffuse and
        // reflected radiation contribute to all leaves.
        out.Rbdown = trgz * Rbeam;
        out.Rddown = tir.Rddn_z * Rdif * svfa + Rdbdn_z * Rb;
        out.Rdup = tir.Rdup_z * Rdif * svfa + Rdbup_z * Rb;
        out.radLsw = 0.5 * (1.0 - tir.tsd.om) * (out.Rddown + out.Rdup + kcos * out.Rbdown);
        out.radLpar = 0.5 * (1.0 - 0.5 * tir.tsd.om) * (out.Rddown + out.Rdup + kcos * out.Rbdown);
        return out;
    }

    // Longwave exchange, as pointmodel::RadlwabsStepCpp sets it out. The
    // ground's hemisphere is canopy over `1 - trdif`, open sky over
    // `trdif*svfa` and surrounding terrain over the rest; the canopy+ground
    // surface sees sky over `svfa` and terrain over the rest. Terrain radiates
    // at air temperature. Canopy and terrain return `1 - em` of what reaches
    // them, so exchange with either carries `em^2`, which `emGround` and
    // `emCanopy` fold into each surface's own emission. At this radiation-only
    // stage no local canopy temperature has been solved, so canopy emission is
    // also approximated at the driving air temperature supplied as `tc`.
    void GridRadlwabsStepCpp(const tirstruct& tir, double svfa, double tc, double Rlw, radmodel2& out)
    {
        const double em2 = surfaceEmissivity * surfaceEmissivity;
        double Rem = sb * utils::rademCpp(tc);  // canopy and terrain emission, before emissivity
        out.lwout = surfaceEmissivity * Rem;
        double trdif = tir.vegetated ? utils::canopyGapMixCpp(tir.trdn, 1.0, std::exp(-tir.pait)) : 1.0;
        double skyG = trdif * svfa;
        out.radGlw = surfaceEmissivity * skyG * Rlw + em2 * (1.0 - skyG) * Rem;
        out.radClw = surfaceEmissivity * svfa * Rlw + em2 * (1.0 - svfa) * Rem;
        out.emGround = surfaceEmissivity * skyG + em2 * (1.0 - skyG);
        // With no canopy `trdif` is 1, so the ground's view and the whole
        // surface's coincide and so do both pairs of quantities above.
        out.emCanopy = surfaceEmissivity * svfa + em2 * (1.0 - svfa);
    }

    // Solar position through time for one location. `sindex` maps azimuth to
    // one of the 24 horizon sectors used for terrain shadowing.
    struct SolarSeries {
        std::vector<double> zend, zenr, azid;
        std::vector<int> sindex;
    };

    SolarSeries computeSolarSeriesCpp(double lat, double lon,
        const IntegerVector& year, const IntegerVector& month, const IntegerVector& day,
        const NumericVector& hour, int tsteps)
    {
        SolarSeries s;
        s.zend.resize(tsteps);
        s.zenr.resize(tsteps);
        s.azid.resize(tsteps);
        s.sindex.resize(tsteps);
        for (int k = 0; k < tsteps; ++k) {
            pointmodel::solmodel sp = pointmodel::solpositionCpp(lat, lon, year[k], month[k], day[k], hour[k]);
            s.zend[k] = sp.zend;
            s.zenr[k] = sp.zenr;
            s.azid[k] = sp.azid;
            s.sindex[k] = static_cast<int>(std::round(sp.azid / 15.0)) % 24;
        }
        return s;
    }

    // Expose solar geometry for callers that need to combine the same solar
    // position with independently calculated terrain horizon angles.
    List computeSolarSeriesRCpp(double lat, double lon, IntegerVector year, IntegerVector month,
        IntegerVector day, NumericVector hour)
    {
        int tsteps = year.size();
        SolarSeries s = computeSolarSeriesCpp(lat, lon, year, month, day, hour, tsteps);
        IntegerVector sindexR(s.sindex.begin(), s.sindex.end());
        return List::create(Named("zend") = s.zend, Named("zenr") = s.zenr,
            Named("azid") = s.azid, Named("sindex") = sindexR);
    }

    // Evaluate the grid radiation formulation for one location through time.
    // This is useful when terrain-aware direct/diffuse/reflected radiation at a
    // specific height is needed without constructing a full spatial grid.
    List radiationPointCpp(double pai, double paia, double x, double lref, double ltra,
        double clump, double gref, double svfa, double slope, double aspect,
        NumericVector zend, NumericVector azid, NumericVector shadow,
        NumericVector Rsw, NumericVector Rdif)
    {
        int tsteps = Rsw.size();
        NumericVector Rbdown(tsteps, NA_REAL), Rddown(tsteps, NA_REAL), Rdup(tsteps, NA_REAL),
            radLsw(tsteps, NA_REAL), radLpar(tsteps, NA_REAL);

        tirstruct tir = GridRadswabsSetupCpp(pai, paia, x, lref, ltra, clump, gref);
        for (int k = 0; k < tsteps; ++k) {
            bool shadow_k = (shadow[k] != 0.0);
            double si = utils::solarindexCpp(slope, aspect, zend[k], azid[k], true, shadow_k);
            double zenr_k = zend[k] * pi / 180.0;
            radmodel2 radm = GridRadswabsStepCpp(tir, pai, clump, gref, svfa, si, shadow_k,
                zenr_k, x, Rsw[k], Rdif[k]);
            Rbdown[k] = radm.Rbdown; Rddown[k] = radm.Rddown; Rdup[k] = radm.Rdup;
            radLsw[k] = radm.radLsw; radLpar[k] = radm.radLpar;
        }
        return List::create(Named("Rbdown") = Rbdown, Named("Rddown") = Rddown,
            Named("Rdup") = Rdup, Named("radLsw") = radLsw, Named("radLpar") = radLpar);
    }

    // ========================================================================
    // Grid radiation and wind core
    // ========================================================================
    // Evaluate each land cell's shortwave, PAR and longwave radiation
    // absorption, together with its aerodynamic state, for every timestep.
    // Climate forcing may be common to the whole domain or vary by
    // coarse reference cell; terrain and vegetation always vary at the fine-grid
    // scale.

    List runmicro1Cpp(DataFrame obstime, DataFrame climdata, List vegp, List soilc,
        NumericMatrix lats, NumericMatrix lons,
        NumericVector zref, double z, NumericVector ufRef, NumericVector HRef, NumericVector dRef,
        NumericVector zmRef, NumericVector shelterc)
    {
        IntegerVector year = obstime["year"];
        IntegerVector month = obstime["month"];
        IntegerVector day = obstime["day"];
        NumericVector hour = obstime["hour"];
        NumericVector tc = climdata["temp"];
        NumericVector Rsw = climdata["swdown"];
        NumericVector Rdif = climdata["difrad"];
        NumericVector Rlw = climdata["lwdown"];
        NumericVector uref = climdata["windspeed"];
        NumericVector pk = climdata["pres"];

        NumericMatrix hgt = vegp["hgt"];
        NumericMatrix pai = vegp["pai"];
        NumericMatrix paia = vegp["paia"];
        NumericMatrix x = vegp["x"];
        NumericMatrix lref = vegp["leafr"];
        NumericMatrix ltra = vegp["leaft"];
        NumericMatrix clump = vegp["clump"];
        NumericMatrix gref = soilc["gref"];
        NumericMatrix slope = soilc["slope"];
        NumericMatrix aspect = soilc["aspect"];
        NumericMatrix svfa = soilc["svfa"];
        NumericVector hor = soilc["hor"];

        int rows = hgt.nrow();
        int cols = hgt.ncol();
        // The time axis is shared even when the climate values themselves vary
        // spatially.
        int tsteps = year.size();
        int n = rows * cols * tsteps;
        IntegerVector dim = { rows, cols, tsteps };

        // Reference climate and friction velocity can represent one common
        // forcing series or spatially varying coarse-cell reference states.
        bool climSpatial = (tc.size() != tsteps);
        bool ufRefSpatial = (ufRef.size() != tsteps);
        bool HRefSpatial = (HRef.size() != tsteps);

        NumericVector radGsw(n, NA_REAL), radGlw(n, NA_REAL), radCsw(n, NA_REAL), radClw(n, NA_REAL);
        NumericVector emGround(n, NA_REAL), emCanopy(n, NA_REAL);
        NumericVector rGh(n, NA_REAL);
        // Ground-to-requested-height resistance, and the top of the canopy's
        // roughness sublayer, for the above-canopy profile.
        NumericVector rGreq(n, NA_REAL);
        NumericMatrix zsOut(rows, cols);
        std::fill(zsOut.begin(), zsOut.end(), NA_REAL);
        // The below-canopy profile's vertical shape at the requested height:
        // fixed by canopy geometry, so one pair of numbers per cell.
        NumericMatrix shapeR(rows, cols), shapeC(rows, cols);
        std::fill(shapeR.begin(), shapeR.end(), NA_REAL);
        std::fill(shapeC.begin(), shapeC.end(), NA_REAL);
        NumericVector Rbdown(n, NA_REAL), Rddown(n, NA_REAL), Rdup(n, NA_REAL);
        NumericVector radLsw(n, NA_REAL), radLpar(n, NA_REAL), lwout(n, NA_REAL), zend(n, NA_REAL);
        NumericVector radCpar(n, NA_REAL), radGpar(n, NA_REAL);
        NumericVector uz(n, NA_REAL), uf(n, NA_REAL), rGz(n, NA_REAL), rGm(n, NA_REAL), rHa(n, NA_REAL), a2(n, NA_REAL), Lout(n, NA_REAL);
        NumericVector uzActual(n, NA_REAL);
        radGsw.attr("dim") = dim; radGlw.attr("dim") = dim; radCsw.attr("dim") = dim; radClw.attr("dim") = dim;
        emGround.attr("dim") = dim; emCanopy.attr("dim") = dim;
        rGh.attr("dim") = dim; rGreq.attr("dim") = dim;
        Rbdown.attr("dim") = dim; Rddown.attr("dim") = dim; Rdup.attr("dim") = dim;
        radLsw.attr("dim") = dim; radLpar.attr("dim") = dim; lwout.attr("dim") = dim; zend.attr("dim") = dim;
        radCpar.attr("dim") = dim; radGpar.attr("dim") = dim;
        uz.attr("dim") = dim; uzActual.attr("dim") = dim; uf.attr("dim") = dim; rGz.attr("dim") = dim; rGm.attr("dim") = dim;
        rHa.attr("dim") = dim; a2.attr("dim") = dim; Lout.attr("dim") = dim;

        // Aerodynamic geometry is fixed by each cell's canopy structure.
        NumericMatrix dOut(rows, cols), zmOut(rows, cols), zhOut(rows, cols);
        std::fill(dOut.begin(), dOut.end(), NA_REAL);
        std::fill(zmOut.begin(), zmOut.end(), NA_REAL);
        std::fill(zhOut.begin(), zhOut.end(), NA_REAL);

        // Solar geometry is common only when all cells share one location;
        // otherwise latitude/longitude give each cell its own sun path.
        bool latSpatial = (lats.nrow() != 1 || lats.ncol() != 1);
        SolarSeries sharedSolar;
        if (!latSpatial) {
            sharedSolar = computeSolarSeriesCpp(lats(0, 0), lons(0, 0), year, month, day, hour, tsteps);
        }

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (NumericMatrix::is_na(hgt(i, j))) continue;

                tirstruct tir = GridRadswabsSetupCpp(pai(i, j), paia(i, j), x(i, j), lref(i, j), ltra(i, j),
                    clump(i, j), gref(i, j));
                // The same canopy with PAR optical properties, derived as the
                // point model derives them, for the foliage PAR that drives
                // stomatal conductance.
                tirstruct tirPAR = GridRadswabsSetupCpp(pai(i, j), paia(i, j), x(i, j), 0.25 * lref(i, j),
                    0.25 * ltra(i, j), clump(i, j), 0.75 * gref(i, j));

                // This cell's own solar position -- a genuine per-cell
                // series where location varies spatially, otherwise the
                // shared series already computed above.
                SolarSeries cellSolar;
                const SolarSeries* solar = &sharedSolar;
                if (latSpatial) {
                    cellSolar = computeSolarSeriesCpp(lats(i, j), lons(i, j), year, month, day, hour, tsteps);
                    solar = &cellSolar;
                }

                // Canopy displacement and roughness set this cell's wind profile.
                double dCell = utils::zeroplanedisCpp(hgt(i, j), pai(i, j));
                double zmCell = utils::roughlengthCpp(hgt(i, j), pai(i, j), dCell);
                dOut(i, j) = dCell;
                zmOut(i, j) = zmCell;
                zhOut(i, j) = utils::scalarRoughlengthCpp(hgt(i, j), pai(i, j), dCell);

                // The transport column's per-canopy constants depend only on
                // this cell's vegetation, so they are built once here rather
                // than at every timestep below. The column exists at every
                // canopy height above zero; its resistances vanish continuously
                // as the canopy does.
                utils::ColumnStruct colCell{};
                bool hasColumn = hgt(i, j) > 0.0;
                if (hasColumn) colCell = utils::columnSetupCpp(hgt(i, j), pai(i, j), dCell);
                // The profile shape, from the column just built, for the one
                // requested height: the only quadrature in this path, once per cell.
                if (hasColumn && z > 0.0 && z < hgt(i, j)) {
                    utils::ProfileShape ps = utils::columnProfileShapeCpp(colCell, z);
                    shapeR(i, j) = ps.fR; shapeC(i, j) = ps.fC;
                }
                if (hasColumn) zsOut(i, j) = colCell.zs;

                // This cell's own reference wind height/geometry -- shared
                // by the whole grid or resolved per cell, resolved once
                // here since it too is time-invariant.
                double zrefCell = (zref.size() == 1) ? zref[0] : zref[i + rows * j];
                double dRefCell = (dRef.size() == 1) ? dRef[0] : dRef[i + rows * j];
                double zmRefCell = (zmRef.size() == 1) ? zmRef[0] : zmRef[i + rows * j];
                // The stable-branch peak of the stability inversion depends only
                // on this cell's geometry.
                double uPeakCell = utils::stableRecoveryPeakCpp(zrefCell, dCell, zmCell);

                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + rows * cols * k;

                    // Select the reference climate state covering this cell and timestep.
                    double tc_i = climSpatial ? tc[idx] : tc[k];
                    double Rsw_i = climSpatial ? Rsw[idx] : Rsw[k];
                    double Rdif_i = climSpatial ? Rdif[idx] : Rdif[k];
                    double Rlw_i = climSpatial ? Rlw[idx] : Rlw[k];
                    double uref_i = climSpatial ? uref[idx] : uref[k];
                    double ufRef_i = ufRefSpatial ? ufRef[idx] : ufRef[k];
                    double HRef_i = HRefSpatial ? HRef[idx] : HRef[k];
                    double pk_i = climSpatial ? pk[idx] : pk[k];
                    double wstar_i = utils::freeConvectiveVelocityCpp(HRef_i, tc_i, pk_i);

                    double solaralt = (pi / 2.0) - solar->zenr[k];
                    double ha = hor[solar->sindex[k] * rows * cols + j * rows + i];
                    bool shadow = (solaralt <= 0.0) || (ha > std::tan(solaralt));
                    double si = utils::solarindexCpp(slope(i, j), aspect(i, j), solar->zend[k], solar->azid[k], true, shadow);

                    radmodel2 radm = GridRadswabsStepCpp(tir, pai(i, j), clump(i, j), gref(i, j),
                        svfa(i, j), si, shadow, solar->zenr[k], x(i, j), Rsw_i, Rdif_i);
                    GridRadlwabsStepCpp(tir, svfa(i, j), tc_i, Rlw_i, radm);
                    radmodel2 radp = GridRadswabsStepCpp(tirPAR, pai(i, j), clump(i, j), 0.75 * gref(i, j),
                        svfa(i, j), si, shadow, solar->zenr[k], x(i, j), Rsw_i, Rdif_i);
                    radCpar[idx] = radp.radCsw; radGpar[idx] = radp.radGsw;
                    radm.zend = solar->zend[k];

                    radGsw[idx] = radm.radGsw; radGlw[idx] = radm.radGlw;
                    emGround[idx] = radm.emGround; emCanopy[idx] = radm.emCanopy;
                    radCsw[idx] = radm.radCsw; radClw[idx] = radm.radClw;
                    Rbdown[idx] = radm.Rbdown; Rddown[idx] = radm.Rddown; Rdup[idx] = radm.Rdup;
                    radLsw[idx] = radm.radLsw; radLpar[idx] = radm.radLpar;
                    lwout[idx] = radm.lwout; zend[idx] = radm.zend;

                    // Scale the reference aerodynamic state to this cell's canopy
                    // geometry and topographic shelter.
                    windresult2 windm = windCpp(z, zrefCell, hgt(i, j), pai(i, j), dCell, zmCell, uref_i,
                        shelterc[idx], ufRef_i, dRefCell, zmRefCell, colCell, hasColumn, uPeakCell, wstar_i);
                    uz[idx] = windm.uz; uzActual[idx] = windm.uzActual; uf[idx] = windm.uf; rGz[idx] = windm.rGz; rGm[idx] = windm.rGm;
                    rGh[idx] = windm.rGh; rGreq[idx] = windm.rGreq;
                    rHa[idx] = windm.rHa; a2[idx] = windm.a2; Lout[idx] = windm.L;
                }
            }
        }

        List out;
        out["radGsw"] = radGsw; out["radGlw"] = radGlw; out["radCsw"] = radCsw; out["radClw"] = radClw;
        out["emGround"] = emGround; out["emCanopy"] = emCanopy;
        out["Rbdown"] = Rbdown; out["Rddown"] = Rddown; out["Rdup"] = Rdup;
        out["radLsw"] = radLsw; out["radLpar"] = radLpar; out["lwout"] = lwout; out["zend"] = zend;
        out["radCpar"] = radCpar; out["radGpar"] = radGpar;
        out["uz"] = uz; out["uzActual"] = uzActual; out["uf"] = uf; out["rGz"] = rGz; out["rGm"] = rGm; out["rHa"] = rHa; out["a2"] = a2; out["L"] = Lout;
        out["rGh"] = rGh; out["shapeR"] = shapeR; out["shapeC"] = shapeC;
        out["rGreq"] = rGreq; out["zs"] = zsOut;
        out["d"] = dOut; out["zm"] = zmOut; out["zh"] = zhOut;
        return out;
    }

    // ========================================================================
    // Wind and turbulent exchange
    // ========================================================================

    windresult2 windCpp(double z, double zref, double hgt, double pai, double d, double zm,
        double uref, double shelterc, double ufRef, double dRef, double zmRef,
        const utils::ColumnStruct& col, bool hasColumn, double uPeak, double wstar)
    {
        // Topographic shelter reduces the atmospheric wind driving the whole
        // profile. A small positive floor avoids a physically unresolved
        // exactly-calm state.
        if (!std::isfinite(shelterc)) shelterc = 1.0;
        if (shelterc < 0.05) shelterc = 0.05;

        // The wind driving this cell: sheltered wind combined with the
        // free-convective velocity as the point model combines them, and
        // floored. The free-convective velocity is a boundary-layer scale set
        // by the reference run's sensible heat, so shelter does not reduce it.
        double Ueff = utils::drivingWindCpp(shelterc * uref, wstar);

        // Scale reference friction velocity for the target cell's own canopy
        // displacement and roughness geometry, and for the ratio of the wind
        // driving this cell to the wind that drove the reference, which is
        // unsheltered and formed the same way.
        double lnDenom = std::log((zref - d) / zm);
        if (lnDenom < 1e-6) lnDenom = 1e-6;
        double muU = std::log((zref - dRef) / zmRef) / lnDenom;
        double UeffRef = utils::drivingWindCpp(uref, wstar);
        double uf = ufRef * muU * (Ueff / UeffRef);
        if (uf < 1e-6) uf = 1e-6;

        // Infer a target-cell stability state from the scaled friction
        // velocity and driving wind. Because local sensible heat flux is not
        // independently solved here, constrain friction velocity around its
        // neutral value before this inversion so geometric differences are not
        // mistaken for unrealistically strong buoyancy effects.
        double ufNeutral = (mc::ka * Ueff) / lnDenom;
        if (uf < 0.5 * ufNeutral) uf = 0.5 * ufNeutral;
        if (uf > 1.5 * ufNeutral) uf = 1.5 * ufNeutral;
        double L = utils::recoverLClosedCpp(uf, Ueff, zref, d, zm, uPeak);

        // `rHa` links the canopy exchange surface to the atmosphere; `rGz` is
        // the whole column from the soil surface, which the point model
        // evaluates with the same relation from its own friction velocity and
        // Obukhov length. A canopy no taller than the soil's own heat
        // roughness height has no interior air column, and the two coincide.
        double zh = utils::scalarRoughlengthCpp(hgt, pai, d);
        double a2 = hasColumn ? utils::columnA2Cpp(col, L) : 0.0;
        double rGz = hasColumn ? utils::columnResistCpp(col, uf, L, zref)
                               : utils::rHaToHeightScalarCpp(zref, d, zh, L, uf);
        // Without a column there is no canopy source beneath the reference, and
        // the ground exchanges with reference air across the whole resistance.
        double rGm = hasColumn ? utils::columnMeanResistCpp(col, uf, L) : rGz;
        // Resistance from the soil surface to canopy top. The below-canopy
        // profile needs it every hour; the column is already built here, so it
        // is produced once with the rest rather than a second time there.
        double rGh = hasColumn ? utils::columnResistCpp(col, uf, L, hgt) : rGz;
        // Canopy exchange surface to the reference height, as the point model forms
        // it: the exchange surface's own empirical segment to canopy top, then the
        // column above it, so this resistance and the diffusivity field are one
        // description at every stability.
        double rHa = rGz;
        if (hasColumn) {
            double r1 = utils::rHaToHeightScalarCpp(hgt, d, zh, L, uf);
            if (r1 < 0.0) r1 = 0.0;
            rHa = r1 + (rGz - rGh);
        }

        // Above the momentum roughness sublayer, wind follows the
        // stability-corrected logarithmic profile; inside it the diffusivity is
        // constant, so wind rises linearly from canopy top to that profile at
        // the sublayer top (Raupach 1992). Inside vegetation, canopy drag
        // attenuates wind exponentially downward from the canopy top according
        // to plant area density and the canopy-top friction-velocity ratio.
        double uh = 0.0;
        if (hgt > 0.0) {
            double psiMh = utils::dpsimCpp(zm / L) - utils::dpsimCpp((hgt - d) / L);
            // Canopy top lies inside the roughness sublayer, whose influence
            // is folded into the roughness length, so it is taken back out
            // again here (Raupach 1992 Eq. 26b).
            uh = (uf / mc::ka) * (std::log((hgt - d) / zm) + utils::momentumInfluenceCpp(pai) + psiMh);
            if (uh < 1e-6) uh = 1e-6;
        }
        double uz;
        if (z >= hgt || hgt <= 0.0) {
            double psiMz = utils::dpsimCpp(zm / L) - utils::dpsimCpp((z - d) / L);
            uz = (uf / mc::ka) * (std::log((z - d) / zm) + psiMz);
            if (hgt > 0.0) {
                double zw = utils::momentumSublayerTopCpp(hgt, pai, d);
                if (z < zw) {
                    double psiMw = utils::dpsimCpp(zm / L) - utils::dpsimCpp((zw - d) / L);
                    double uw = (uf / mc::ka) * (std::log((zw - d) / zm) + psiMw);
                    uz = uh + (uw - uh) * (z - hgt) / (zw - hgt);
                }
            }
        }
        else {
            double Be = uf / uh;
            double Lc = 1.0 / (0.25 * pai / hgt);
            double Lm = 2.0 * Be * Be * Be * Lc;
            uz = uh * std::exp(Be * (z - hgt) / Lm);
        }
        // The profile is built on the effective wind, which carries the
        // free-convective velocity and the floor so that exchange stays well
        // behaved in light wind; it is not the actual wind. The reported wind is
        // the profile scaled by the cell's actual (sheltered) wind over the
        // profile's own wind at the reference height, so that it equals the
        // actual wind there, as in the point model.
        double uzActual = uz;
        {
            double psiMr = utils::dpsimCpp(zm / L) - utils::dpsimCpp((zref - d) / L);
            double uProfileRef = (uf / mc::ka) * (lnDenom + psiMr);
            if (uProfileRef > 0.0) uzActual = uz * (shelterc * uref) / uProfileRef;
        }
        // On the exchange scale, wind at any height cannot exceed the effective
        // driving wind that produces it.
        if (uz > Ueff) uz = Ueff;

        // Resistance from the soil surface to the requested height, wanted only
        // where that height is above canopy top: the above-canopy profile needs
        // the column inside the roughness sublayer, and the column is built here.
        // Closed form at every height, so no quadrature enters the hourly loop.
        double rGreq = (hasColumn && z >= hgt) ? utils::columnResistCpp(col, uf, L, z) : NA_REAL;

        windresult2 out;
        out.uz = uz;
        out.uzActual = uzActual;
        out.uf = uf;
        out.rGz = rGz;
        out.rGm = rGm;
        out.rGh = rGh;
        out.rGreq = rGreq;
        out.rHa = rHa;
        out.a2 = a2;
        out.L = L;
        return out;
    }

    static pointmodel::soilpstruct miniSoilpCpp(double Vq, double Vm, double Vo, double Mc,
        double thetaS, double psie, double b);

    NumericVector dampingDepthGridCpp(NumericMatrix Vq, NumericMatrix Vm, NumericMatrix Vo,
        NumericMatrix Mc, NumericMatrix thetaS, NumericMatrix psie, NumericMatrix b,
        NumericVector soilm, NumericVector Ts, NumericVector pk)
    {
        IntegerVector dim = soilm.attr("dim");
        int rows = dim[0], cols = dim[1], tsteps = dim[2];
        if (tsteps % 24 != 0) {
            Rcpp::stop("dampingDepthGridCpp: tsteps (%d) is not a multiple of 24", tsteps);
        }
        bool pkSpatial = (pk.size() != tsteps);
        NumericVector out(rows * cols * tsteps, NA_REAL);
        out.attr("dim") = dim;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (NumericMatrix::is_na(Vq(i, j)) || NumericMatrix::is_na(thetaS(i, j))) continue;
                pointmodel::soilpstruct sp = miniSoilpCpp(Vq(i, j), Vm(i, j), Vo(i, j), Mc(i, j),
                    thetaS(i, j), psie(i, j), b(i, j));
                for (int d = 0; d < tsteps / 24; ++d) {
                    double th = 0.0, tc = 0.0, pp = 0.0;
                    for (int h = 0; h < 24; ++h) {
                        int idx = i + rows * j + rows * cols * (d * 24 + h);
                        th += soilm[idx];
                        tc += Ts[idx];
                        pp += pkSpatial ? pk[idx] : pk[d * 24 + h];
                    }
                    double D = pointmodel::diurnalDampingDepthCpp(sp, th / 24.0, tc / 24.0, pp / 24.0);
                    for (int h = 0; h < 24; ++h) out[i + rows * j + rows * cols * (d * 24 + h)] = D;
                }
            }
        }
        return out;
    }

    NumericVector soilWaterPotentialGridCpp(NumericMatrix thetaS, NumericMatrix psie, NumericMatrix b,
        NumericVector soilm)
    {
        IntegerVector dim = soilm.attr("dim");
        int rows = dim[0], cols = dim[1], tsteps = dim[2];
        NumericVector out(rows * cols * tsteps, NA_REAL);
        out.attr("dim") = dim;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (NumericMatrix::is_na(thetaS(i, j))) continue;
                pointmodel::soilpstruct sp = miniSoilpCpp(0.0, 0.0, 0.0, 0.0, thetaS(i, j), psie(i, j), b(i, j));
                for (int t = 0; t < tsteps; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    if (NumericVector::is_na(soilm[idx])) continue;
                    out[idx] = pointmodel::waterPotentialCpp(sp, soilm[idx], 0) / 1000.0; // J/kg -> MPa
                }
            }
        }
        return out;
    }


    List surfaceResistGridCpp(List rad, List clim, NumericVector soilm, List soil, List veg,
        List ref)
    {
        NumericVector radCpar = rad["radCpar"], radGpar = rad["radGpar"], RabsCanopy = rad["RabsCanopy"];
        NumericVector zend = rad["zend"], rHa = rad["rHa"], rGz = rad["rGz"], rGm = rad["rGm"];
        NumericVector emCanopy = rad["emCanopy"];
        NumericVector Rsw = clim["Rsw"], Rdif = clim["Rdif"], Ta = clim["Ta"], rh = clim["rh"];
        NumericVector pk = clim["pk"], precip = clim["precip"], Ca = clim["Ca"];
        NumericMatrix thetaS = soil["thetaS"], psie = soil["psie"], b = soil["b"];
        NumericMatrix hgt = veg["hgt"], pai = veg["pai"], x = veg["x"], lref = veg["leafr"], ltra = veg["leaft"];
        NumericMatrix svfa = veg["svfa"], gsmax = veg["gsmax"], LfracM = veg["Lfrac"];
        NumericMatrix Vcmax25 = veg["Vcmax25"], Tup = veg["Tup"], Tlw = veg["Tlw"], Dcrit = veg["Dcrit"];
        NumericMatrix alpha = veg["alpha"], f0 = veg["f0"], fd = veg["fd"], psi50 = veg["psi50"];
        NumericMatrix apsi = veg["apsi"], rpmin = veg["rpmin"];
        LogicalMatrix isC3 = veg["isC3"];
        NumericVector wetRef = ref["wetShare"], availRef = ref["filmAvailable"], paiRef = ref["pai"];

        IntegerVector dim = soilm.attr("dim");
        int rows = dim[0], cols = dim[1], tsteps = dim[2];
        bool climSpatial = (Ta.size() != tsteps);
        bool refSpatial = (wetRef.size() != tsteps);
        int n = rows * cols * tsteps;
        NumericVector rSurf(n, NA_REAL), hSurf(n, NA_REAL), hFol(n, NA_REAL), hrOut(n, NA_REAL);
        rSurf.attr("dim") = dim;
        hSurf.attr("dim") = dim;
        hFol.attr("dim") = dim;
        hrOut.attr("dim") = dim;
        const double dT = 3600.0;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (NumericMatrix::is_na(hgt(i, j)) || NumericMatrix::is_na(thetaS(i, j))) continue;
                pointmodel::soilpstruct sp = miniSoilpCpp(0.0, 0.0, 0.0, 0.0, thetaS(i, j), psie(i, j), b(i, j));
                double paiC = pai(i, j);
                bool vegetated = (hgt(i, j) > 0.0 && paiC > 0.0);
                pointmodel::vegpstruct vp{};
                if (vegetated) {
                    vp.Vcmax25 = Vcmax25(i, j); vp.Tup = Tup(i, j); vp.Tlw = Tlw(i, j);
                    vp.Dcrit = Dcrit(i, j); vp.alpha = alpha(i, j); vp.f0 = f0(i, j); vp.fd = fd(i, j);
                    vp.psi50 = psi50(i, j); vp.apsi = apsi(i, j); vp.rpmin = rpmin(i, j);
                    vp.gsmaxCap = NumericMatrix::is_na(gsmax(i, j)) ? -1.0 : gsmax(i, j);
                }
                bool C3 = vegetated ? static_cast<bool>(isC3(i, j)) : true;
                double om = 0.25 * (lref(i, j) + ltra(i, j));
                double paiRefCell = (paiRef.size() == 1) ? paiRef[0] : paiRef[i + rows * j];
                double Lfrac = LfracM(i, j);

                for (int t = 0; t < tsteps; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    double sm = soilm[idx], rha = rHa[idx], rgz = rGz[idx];
                    if (NumericVector::is_na(sm) || NumericVector::is_na(rha) || NumericVector::is_na(rgz)) continue;
                    double Ta_t = climSpatial ? Ta[idx] : Ta[t];
                    double rh_t = climSpatial ? rh[idx] : rh[t];
                    double pk_t = climSpatial ? pk[idx] : pk[t];

                    // Soil surface humidity, as the grid's zero-flux pass forms it.
                    double hr = pointmodel::soilrelhumCpp(sp, Ta_t, sm);
                    if (hr < 1e-6) hr = 1e-6;
                    if (hr > 1.0) hr = 1.0;
                    hrOut[idx] = hr;
                    if (!vegetated) {
                        rSurf[idx] = 0.0;
                        hSurf[idx] = hr;
                        hFol[idx] = 0.0;
                        continue;
                    }

                    double Rsw_t = climSpatial ? Rsw[idx] : Rsw[t];
                    double Rdif_t = climSpatial ? Rdif[idx] : Rdif[t];
                    double P_t = climSpatial ? precip[idx] : precip[t];
                    double Ca_t = (Ca.size() == tsteps) ? Ca[t] : Ca[idx];

                    // Sunlit and shaded foliage and their absorbed PAR, as
                    // pointmodel::StomatalSetupCpp forms them, rescaled to this
                    // cell's foliage PAR.
                    double zd = zend[idx];
                    if (zd > 90.0) zd = 90.0;
                    double k = utils::canopyKCpp(zd * pi / 180.0, x(i, j));
                    double Lsun = (1.0 - std::exp(-k * paiC * Lfrac)) / k;
                    double Lshade = paiC * Lfrac - Lsun;
                    double shadeAbs = Rdif_t * ((1.0 - std::exp(-paiC)) / paiC) * (1.0 - om);
                    double sunAbs = (Rsw_t - Rdif_t) * k * (1.0 - om) + shadeAbs;
                    double unclumped = sunAbs * Lsun + shadeAbs * Lshade;
                    double target = radCpar[idx] - radGpar[idx];
                    double rescale = (unclumped > 0.0) ? target / unclumped : 1.0;
                    double PARsun = sunAbs * rescale, PARshade = shadeAbs * rescale;

                    // Wet share is copied from the reference state (and is 1
                    // whenever rain is falling); only the available film water
                    // is scaled by this cell's PAI relative to the reference PAI.
                    // With no reference plant area, both wet share and available
                    // film water come from this hour's rain interception.
                    double w0, avail;
                    if (paiRefCell > 0.0) {
                        w0 = (P_t > 0.0) ? 1.0 : (refSpatial ? wetRef[idx] : wetRef[t]);
                        avail = (refSpatial ? availRef[idx] : availRef[t]) * paiC / paiRefCell;
                    } else {
                        w0 = (P_t > 0.0) ? 1.0 : 0.0;
                        avail = (1.0 - std::exp(-0.5 * paiC)) * P_t;
                    }

                    double ph = utils::phairCpp(Ta_t, pk_t);
                    double Tk = Ta_t + 273.15;
                    double ea = 1000.0 * utils::satvapCpp(Ta_t) * rh_t / 100.0;

                    pointmodel::envstruct env{};
                    env.tair = Ta_t; env.rh = rh_t; env.pk = pk_t; env.Ca = Ca_t;
                    env.psi_r = pointmodel::waterPotentialCpp(sp, sm, 0) / 1000.0;

                    // The network's legs: the node at whichever of the mean source
                    // height and the exchange surface is nearer the reference, the
                    // foliage's leg through it, the soil's the rest of its path.
                    double rRise = rgz - rGm[idx];
                    double rsh = (rRise > 0.0) ? std::min(rRise, rha) : 0.0;
                    double rw = std::max(rha - rsh, 0.0);
                    double rbN = std::max(rgz - rsh, 1e-9);
                    double gWetLeg = 1.0 / std::max(rw, 0.1);   // numerical floor, as the point model's
                    double bN = 1.0 / rbN, cN = (rsh > 0.0) ? 1.0 / rsh : 0.0;
                    double eaK = ea / 1000.0;
                    double egK = hr * utils::satvapCpp(Ta_t);

                    // Two passes: the first at air temperature and with the film
                    // unlimited gives the canopy temperature at which the second
                    // evaluates stomata and the film's free evaporation.
                    double tc = Ta_t, rV = 0.0, hV = 1.0, hN = 0.0;
                    for (int pass = 0; pass < 2; ++pass) {
                        env.tcanopy = tc;
                        // Each leaf's conductance is at least its closed-stomata
                        // value, as in the point model.
                        const double gsMin = ph / mc::rStomClosed;
                        env.PARabs = PARsun;
                        double gsun = std::max(pointmodel::leafgsCpp(env, vp, 0.5 * hgt(i, j), C3), gsMin);
                        env.PARabs = PARshade;
                        double gshade = std::max(pointmodel::leafgsCpp(env, vp, 0.5 * hgt(i, j), C3), gsMin);
                        double Gs = gsun * Lsun + gshade * Lshade;
                        double gStom = (Gs > 0.0) ? Gs / ph : 0.0;
                        double gDryLeg = (gStom > 0.0) ? 1.0 / (1.0 / gStom + rw) : 0.0;
                        double ecK = utils::satvapCpp(tc);
                        // The wet share the film's water sustains over the step, in the
                        // point model's closed form, on the second pass.
                        double wtN = w0;
                        if (pass == 1 && w0 > 0.0) {
                            double cm = 1000.0 * mc::Mw * dT / (mc::RgasC * Tk);
                            auto nodeAt = [&](double aa) {
                                return (rsh > 0.0) ? (aa * ecK + bN * egK + cN * eaK) / (aa + bN + cN) : eaK; };
                            double aW = w0 * gWetLeg + (1.0 - w0) * gDryLeg;
                            double filmW = cm * w0 * gWetLeg * (ecK - nodeAt(aW));
                            if (filmW > 0.0 && filmW > avail) {
                                if (rsh > 0.0) {
                                    double N = bN * (ecK - egK) + cN * (ecK - eaK);
                                    double S0 = gDryLeg + bN + cN, dA = gWetLeg - gDryLeg;
                                    double A = avail / (cm * gWetLeg * N);
                                    double den = 1.0 - A * dA;
                                    wtN = (den > 0.0) ? A * S0 / den : w0;
                                } else {
                                    wtN = avail / (cm * gWetLeg * (ecK - eaK));
                                }
                                if (!(wtN >= 0.0) || wtN > w0) wtN = w0;
                            }
                        }
                        double aN = wtN * gWetLeg + (1.0 - wtN) * gDryLeg;
                        hN = aN / (aN + bN);
                        rV = 1.0 / (aN + bN) + rsh;
                        hV = hN + (1.0 - hN) * hr;
                        if (pass == 0) {
                            tc = utils::penmanMonteithCpp(RabsCanopy[idx], Ta_t, pk_t, rh_t,
                                emCanopy[idx], rha, rV, Ta_t, 0.0, hV);
                        }
                    }
                    rSurf[idx] = rV - rha;
                    hSurf[idx] = hV;
                    hFol[idx] = hN;
                }
            }
        }
        return List::create(Named("rSurf") = rSurf, Named("hSurf") = hSurf, Named("hFol") = hFol,
            Named("hr") = hrOut);
    }

    NumericVector referenceGroundResistCpp(NumericVector uref, NumericVector ufRef,
        NumericVector HRef, NumericVector Ta, NumericVector pk,
        double zref, double hRef, double paiRef)
    {
        int n = uref.size();
        NumericVector out(n, NA_REAL);
        if (!std::isfinite(hRef) || !std::isfinite(paiRef)) return out;
        double dRef = utils::zeroplanedisCpp(hRef, paiRef);
        double zmRef = utils::roughlengthCpp(hRef, paiRef, dRef);
        utils::ColumnStruct col{};
        bool hasColumn = hRef > 0.0;
        if (hasColumn) col = utils::columnSetupCpp(hRef, paiRef, dRef);
        double uPeak = utils::stableRecoveryPeakCpp(zref, dRef, zmRef);
        for (int t = 0; t < n; ++t) {
            double wstar = utils::freeConvectiveVelocityCpp(HRef[t], Ta[t], pk[t]);
            windresult2 w = windCpp(zref, zref, hRef, paiRef, dRef, zmRef, uref[t], 1.0,
                ufRef[t], dRef, zmRef, col, hasColumn, uPeak, wstar);
            out[t] = w.rGz;
        }
        return out;
    }

    // Distribute the reference surface-soil moisture across the landscape.
    // Each cell inherits the reference degree of wetness within its own soil-
    // specific residual-to-saturated range, then a topographic-wetness anomaly
    // shifts that fraction in logit space so the result remains physically
    // bounded between the cell's dry and saturated limits. Last, the cell is
    // raised toward saturation by its wet weight (0 to 1): the share of the
    // gap to saturation that is closed in every hour, 1 where the ground is
    // permanently wet, whatever the reference soil is doing.
    NumericVector soilmDistributeCpp(NumericMatrix Smin, NumericMatrix Smax, NumericMatrix tadd,
        NumericMatrix wet, NumericVector theta0, int tsteps)
    {
        int rows = Smin.nrow();
        int cols = Smin.ncol();
        bool theta0Spatial = (theta0.size() != tsteps);
        int n = rows * cols * tsteps;
        IntegerVector dim = { rows, cols, tsteps };
        NumericVector soilm(n, NA_REAL);
        soilm.attr("dim") = dim;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double smn = Smin(i, j);
                double smx = Smax(i, j);
                double ta = tadd(i, j);
                double w = wet(i, j);
                if (NumericMatrix::is_na(smn) || NumericMatrix::is_na(smx) || NumericMatrix::is_na(ta) ||
                    NumericMatrix::is_na(w)) continue;
                double rge = smx - smn;
                for (int k = 0; k < tsteps; ++k) {
                    int idx = i + rows * j + rows * cols * k;
                    // The reference soil moisture is expressed as a
                    // fraction of this cell's own physically plausible
                    // range, shifted in logit space by this cell's own
                    // topographic wetness anomaly (working in logit space
                    // keeps the shifted value within [0, 1] regardless of
                    // how large the anomaly is), then rescaled back into
                    // this cell's own physical units.
                    double theta0_i = theta0Spatial ? theta0[idx] : theta0[k];
                    // An evaporating surface dries below the residual water
                    // content, and the wetness anomaly has no meaning there:
                    // the reference's own value is carried across unchanged, so
                    // that a cell identical to the reference reproduces it.
                    double thetaCell;
                    if (theta0_i < smn) {
                        thetaCell = theta0_i;
                    } else {
                        double theta = (theta0_i - smn) / rge;
                        if (theta > 0.9999) theta = 0.9999;
                        if (theta < 0.0001) theta = 0.0001;
                        double lt = std::log(theta / (1.0 - theta));
                        double sm = lt + ta;
                        sm = 1.0 / (1.0 + std::exp(-sm));
                        thetaCell = sm * rge + smn;
                    }
                    soilm[idx] = thetaCell + w * (smx - thetaCell);
                }
            }
        }
        return soilm;
    }

    // ========================================================================
    // Ground heat storage and surface temperature
    // ========================================================================

    // Describe the strength of the diurnal surface-temperature cycle by its
    // first 24-hour Fourier harmonic. Its amplitude is used below as a measure
    // of how strongly a target cell's surface forcing differs from the
    // reference simulation.
    struct Harmonic { double mean, amp, phaseHr; };

    static std::vector<Harmonic> dailyHarmonicsCpp(const std::vector<double>& y, int nDays)
    {
        std::vector<Harmonic> out(nDays);
        const double w24 = 2.0 * pi / 24.0;
        for (int d = 0; d < nDays; ++d) {
            int base = d * 24;
            double c0 = 0.0, c1 = 0.0, c2 = 0.0;
            for (int h = 0; h < 24; ++h) {
                double a = w24 * h;
                c0 += y[base + h];
                c1 += y[base + h] * std::cos(a);
                c2 += y[base + h] * std::sin(a);
            }
            c0 /= 24.0; c1 *= 2.0 / 24.0; c2 *= 2.0 / 24.0;
            double amp = std::sqrt(c1 * c1 + c2 * c2);
            double phaseRad = std::atan2(c2, c1);
            double phaseHr = phaseRad * 24.0 / (2.0 * pi);
            out[d] = { c0, amp, phaseHr };
        }
        return out;
    }

    // ========================================================================
    // Diurnal timing adjustment for ground heat flux
    // ========================================================================
    // Slopes can shift the timing of direct solar heating relative to the flat
    // reference surface. Sunrise/sunset landmarks are therefore used to sample
    // the reference ground-heat-flux cycle at a corresponding phase of the
    // target cell's illuminated day rather than blindly at the same clock hour.

    struct Crossings { double rise, set; bool valid; };

    // Locate the first and last crossings of positive direct-beam exposure
    // within a day. Invalid crossings represent days without a usable
    // sunrise/sunset pair, such as polar night or effectively continuous day.
    static Crossings findCrossingsCpp(const std::vector<double>& si)
    {
        bool above[24];
        bool anyAbove = false;
        for (int h = 0; h < 24; ++h) {
            above[h] = si[h] > 1e-6;
            anyAbove = anyAbove || above[h];
        }
        if (!anyAbove) return { 0.0, 0.0, false };

        int riseIdx = -1;
        for (int h = 0; h < 24; ++h) {
            bool prevAbove = (h == 0) ? false : above[h - 1];
            if (above[h] && !prevAbove) { riseIdx = h; break; }
        }
        int setIdx = -1;
        for (int h = 0; h < 24; ++h) {
            bool nextAbove = (h == 23) ? false : above[h + 1];
            if (above[h] && !nextAbove) setIdx = h;
        }

        double rise;
        if (riseIdx > 0) {
            double h0 = riseIdx - 1, h1 = riseIdx;
            double s0 = si[riseIdx - 1], s1 = si[riseIdx];
            rise = h0 + (0.0 - s0) / (s1 - s0) * (h1 - h0);
        } else {
            rise = riseIdx;
        }
        double set;
        if (setIdx < 23) {
            double h0 = setIdx, h1 = setIdx + 1;
            double s0 = si[setIdx], s1 = si[setIdx + 1];
            set = h0 + (0.0 - s0) / (s1 - s0) * (h1 - h0);
        } else {
            set = setIdx;
        }

        double nightSpan = 24.0 - (set - rise);
        if (nightSpan <= 1.0) return { 0.0, 0.0, false };

        return { rise, set, true };
    }

    // Maps a target-cell clock hour onto the reference's own clock hour,
    // piecewise: night-before -> night-before, day -> day (rescaled to fit
    // between the two sides' own rise/set), night-after -> night-after.
    static double warpTimeCpp(double tTarget, double riseT, double setT, double riseR, double setR)
    {
        if (tTarget <= riseT) {
            double denom = std::max(riseT, 1e-6);
            return riseR * (tTarget / denom);
        } else if (tTarget >= setT) {
            double denom = std::max(24.0 - setT, 1e-6);
            return setR + (24.0 - setR) * ((tTarget - setT) / denom);
        } else {
            double denom = std::max(setT - riseT, 1e-6);
            return riseR + (setR - riseR) * ((tTarget - riseT) / denom);
        }
    }

    // Construct the surface-layer soil properties needed to evaluate local
    // surface humidity, thermal conductivity and thermal damping without
    // building a full multi-layer soil profile for every grid cell.
    static pointmodel::soilpstruct miniSoilpCpp(double Vq, double Vm, double Vo, double Mc,
        double thetaS, double psie, double b)
    {
        pointmodel::soilpstruct sp;
        sp.nLayers = 1;
        sp.FreeDrain = true;
        sp.gref = 0.0;
        sp.grefPAR = 0.0;
        sp.Vq = { Vq };
        sp.Vm = { Vm };
        sp.Vo = { Vo };
        sp.Mc = { Mc };
        // Convert the stored positive air-entry magnitude to the negative
        // matric potential used by the soil-water relations.
        sp.psie = { -std::abs(psie) };
        sp.b = { b };
        sp.thetaR = { 0.0 };
        sp.thetaS = { thetaS };
        sp.Ksat = { 0.0 };
        sp.psi_min = { 0.0 };
        return sp;
    }

    // Linearised exchange coefficient (W/m2/K) of a soil surface at temperature
    // Ts: emitted longwave, sensible heat and moisture-limited latent heat,
    // each per kelvin of surface temperature, across resistance r.
    static double surfaceExchangeCoefCpp(double Ts, double Ta, double pk, double em, double r, double hr)
    {
        double ph = utils::phairCpp(Ta, pk);
        double cp = utils::cpairCpp(Ta);
        double des = utils::satvapCpp(Ts + 0.5) - utils::satvapCpp(Ts - 0.5);
        return 4.0 * em * sb * std::pow(Ts + 273.15, 3.0) + ph * cp / r +
            utils::latentHeatCpp(Ts) * ph / (pk * r) * hr * des;
    }

    // Empirical envelope for plausible surface-air temperature differences as
    // a function of net radiation, with a physical floor below the dew point
    // against extreme nocturnal over-cooling. A backstop on the grid's one-shot
    // surface temperatures, applied alike to the ground and to the canopy
    // exchange surface; not an additional heat flux.
    //
    // The lines are the tightest straight bounds that contain every hour of
    // fully converged point-model runs over the regimes that set them: the
    // upper edge by bare and sparse ground, sheltered or in a hollow, at
    // midday; the lower by closed canopy on a slope or under a restricted sky
    // at night; both with warmer and drier variants of the same year. Nothing
    // the point model itself produces is clipped. `emView` is the coefficient
    // that surface's own balance puts on its emitted longwave (`emGround` or
    // `emCanopy`), so `Rnet` is the energy it actually has at air temperature.
    // Intercepts carry 2 K beyond the fitted bound, which is more than the
    // 1.3 K an 8 K warmer, 25% drier year moved those bounds by; the runs
    // behind them are one mid-latitude location, and conditions outside that
    // year's range are not evidence these lines cover them.
    //
    // The floor is physical, not fitted. A surface below the air's dew point
    // gains latent heat from condensation as well as sensible heat from the
    // air, so cooling beyond the dew point is divided by 1 + Delta/gamma, with
    // Delta the slope of saturation vapour pressure at the dew point and gamma
    // the psychrometric constant. The most a surface can cool below air at
    // night, about 10 K on observed surfaces worldwide, less the air's own
    // dew-point depression D, is what remains to take it below the dew point:
    // T_d - T_s <= max(10 - D, 0) / (1 + Delta/gamma). In saturated air that
    // is near 10 K in deep cold and 2-5 K in warm weather, tighter at altitude
    // where gamma is smaller; drier air leaves less room.
    static double boundSurfaceTempCpp(double ts, double Ta, double rh, double pk, double rabs,
        double emView)
    {
        const double SURFACE_AIR_MAX = 10.0;
        const double DEWPOINT_MARGIN = 1.5;
        double Rnet = rabs - emView * sb * utils::rademCpp(Ta);
        double tdifMax = 19.31 + 0.01086 * Rnet;
        double tdifMin = -10.98 + 0.0111425 * Rnet;
        if (ts > Ta + tdifMax) ts = Ta + tdifMax;
        if (ts < Ta + tdifMin) ts = Ta + tdifMin;
        double tdew = utils::dewpointCpp(Ta, rh);
        double slope = utils::satvapCpp(tdew + 0.5) - utils::satvapCpp(tdew - 0.5);
        double gamma = utils::cpairCpp(tdew) * pk / utils::latentHeatCpp(tdew);
        double room = SURFACE_AIR_MAX - (Ta - tdew);
        if (room < 0.0) room = 0.0;
        double tfloor = tdew - room / (1.0 + slope / gamma);
        // Observed surfaces during dew sit a median 1.3 K below the air's dew
        // point (interquartile range 0.7-2.0 K). The one-shot estimate is not
        // trusted beyond 1.5 K, which is well inside the physical bound.
        if (tfloor < tdew - DEWPOINT_MARGIN) tfloor = tdew - DEWPOINT_MARGIN;
        if (ts < tfloor) ts = tfloor;
        return ts;
    }

    // Daily swing of ground heat flux (W/m2) sustained by a zero-flux surface
    // temperature swing A0 (K), for exchange coefficient h and soil admittance Y.
    static double fluxSwingCpp(double A0, double h, double Y)
    {
        return A0 * h * Y / std::sqrt(h * h + std::sqrt(2.0) * h * Y + Y * Y);
    }

    // Estimate each cell's ground heat flux and surface temperature without a
    // full local soil-diffusion solve. First solve the local surface energy
    // balance with G = 0, then scale and phase-shift the reference model's
    // resolved diurnal ground-heat-flux cycle according to the target cell's
    // own surface-temperature amplitude and soil thermal properties. A second
    // surface-energy-balance calculation then includes that estimated G.
    List groundHeatFluxCpp(
        NumericVector RabsGround, NumericVector rGz, NumericVector rGm,
        NumericVector RabsCanopy, NumericVector rHa, NumericVector rSurf_ref,
        NumericVector hSurf_ref, NumericVector hFol_ref, NumericVector soilm,
        NumericVector Ta, NumericVector rh, NumericVector pk,
        NumericMatrix Vq, NumericMatrix Vm, NumericMatrix Vo,
        NumericMatrix Mc, NumericMatrix thetaS, NumericMatrix psie, NumericMatrix b,
        NumericMatrix svfa, NumericVector emGround, NumericVector emCanopy,
        NumericMatrix slope, NumericMatrix aspect,
        NumericVector Rsw, NumericVector Rdif,
        NumericMatrix lats, NumericMatrix lons,
        IntegerVector year, IntegerVector month, IntegerVector day, NumericVector hour,
        NumericVector G_ref, NumericVector RabsGround_ref, NumericVector rGz_ref,
        NumericVector theta0_ref, NumericVector emGround_ref,
        NumericVector VqRef, NumericVector VmRef, NumericVector VoRef, NumericVector McRef,
        NumericVector thetaSRef, NumericVector psieRef, NumericVector bRef)
    {
        IntegerVector dim = RabsGround.attr("dim");
        int rows = dim[0], cols = dim[1], tsteps = dim[2];
        int nHours = tsteps;
        if (nHours % 24 != 0) {
            Rcpp::stop("groundHeatFluxCpp: tsteps (%d) is not a multiple of 24 -- the daily "
                "harmonic fit needs complete real calendar-day blocks", nHours);
        }
        int nDays = nHours / 24;
        const double MU_G_MIN = 0.1;
        // Below this predicted daily swing of the reference ground heat flux
        // (W/m2) the reference carries no usable diurnal cycle to scale.
        const double G_SWING_MIN = 5.0;
        const double G_ABS_MAX = 500.0;
        const double TS_DEVIATION_MAX = 20.0;

        // Climate forcing and reference soil/flux diagnostics can each be
        // common to the domain or spatially varying by reference cell.
        bool climSpatial = (Ta.size() != nHours);
        bool refSpatial = (G_ref.size() != nHours);
        bool rSurfSpatial = (rSurf_ref.size() != nHours);

        // Establish the flat-surface solar day used as the timing reference
        // for each cell's slope/aspect-specific illumination cycle.
        bool latSpatial = (lats.nrow() != 1 || lats.ncol() != 1);
        SolarSeries sharedSolar;
        std::vector<double> siRefShared;
        std::vector<Crossings> crossRShared;
        if (!latSpatial) {
            sharedSolar = computeSolarSeriesCpp(lats(0, 0), lons(0, 0), year, month, day, hour, nHours);
            siRefShared.resize(nHours);
            for (int t = 0; t < nHours; ++t) {
                siRefShared[t] = utils::solarindexCpp(0.0, 0.0, sharedSolar.zend[t], sharedSolar.azid[t]);
            }
            crossRShared.resize(nDays);
            for (int d = 0; d < nDays; ++d) {
                std::vector<double> siRefDay(siRefShared.begin() + d * 24, siRefShared.begin() + d * 24 + 24);
                crossRShared[d] = findCrossingsCpp(siRefDay);
            }
        }

        int n = rows * cols * tsteps;
        NumericVector G_est(n, NA_REAL), Ts_est(n, NA_REAL), Tc_est(n, NA_REAL);
        IntegerVector outDim = { rows, cols, tsteps };
        G_est.attr("dim") = outDim;
        Ts_est.attr("dim") = outDim;
        Tc_est.attr("dim") = outDim;   // the bulk surface where it is solved with the ground
        NumericVector Tc_one(n, NA_REAL);
        Tc_one.attr("dim") = outDim;   // its one-surface estimate where the pass forms one instead

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double cVq = Vq(i, j), cVm = Vm(i, j), cVo = Vo(i, j), cMc = Mc(i, j);
                double cThetaS = thetaS(i, j), cPsie = psie(i, j), cB = b(i, j);
                double cSvfa = svfa(i, j);  // land mask only: sky exposure enters
                                            // through emGround / emCanopy
                if (NumericMatrix::is_na(cVq) || NumericMatrix::is_na(cVm) ||
                    NumericMatrix::is_na(cVo) || NumericMatrix::is_na(cMc) ||
                    NumericMatrix::is_na(cThetaS) || NumericMatrix::is_na(cPsie) ||
                    NumericMatrix::is_na(cB) || NumericMatrix::is_na(cSvfa)) continue;

                pointmodel::soilpstruct sp = miniSoilpCpp(cVq, cVm, cVo, cMc, cThetaS, cPsie, cB);

                // Gather the climate and reference-state series applying to this cell.
                std::vector<double> TaC(nHours), rhC(nHours), pkC(nHours), RswC(nHours), RdifC(nHours);
                std::vector<double> GrefC(nHours), RabsGrefC(nHours), rGzrefC(nHours), theta0refC(nHours);
                // The reference's own longwave emission coefficient, so its side of
                // the scaling is evaluated on the same terms as the cell's: a cell
                // identical to the reference must still scale by exactly one.
                std::vector<double> emGrefC(nHours);
                for (int t = 0; t < nHours; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    TaC[t] = climSpatial ? Ta[idx] : Ta[t];
                    rhC[t] = climSpatial ? rh[idx] : rh[t];
                    pkC[t] = climSpatial ? pk[idx] : pk[t];
                    RswC[t] = climSpatial ? Rsw[idx] : Rsw[t];
                    RdifC[t] = climSpatial ? Rdif[idx] : Rdif[t];
                    GrefC[t] = refSpatial ? G_ref[idx] : G_ref[t];
                    RabsGrefC[t] = refSpatial ? RabsGround_ref[idx] : RabsGround_ref[t];
                    rGzrefC[t] = refSpatial ? rGz_ref[idx] : rGz_ref[t];
                    theta0refC[t] = refSpatial ? theta0_ref[idx] : theta0_ref[t];
                    emGrefC[t] = refSpatial ? emGround_ref[idx] : emGround_ref[t];
                }
                double VqRefCell = (VqRef.size() == 1) ? VqRef[0] : VqRef[i + rows * j];
                double VmRefCell = (VmRef.size() == 1) ? VmRef[0] : VmRef[i + rows * j];
                double VoRefCell = (VoRef.size() == 1) ? VoRef[0] : VoRef[i + rows * j];
                double McRefCell = (McRef.size() == 1) ? McRef[0] : McRef[i + rows * j];
                double thetaSRefCell = (thetaSRef.size() == 1) ? thetaSRef[0] : thetaSRef[i + rows * j];
                double psieRefCell = (psieRef.size() == 1) ? psieRef[0] : psieRef[i + rows * j];
                double bRefCell = (bRef.size() == 1) ? bRef[0] : bRef[i + rows * j];
                pointmodel::soilpstruct spRef = miniSoilpCpp(VqRefCell, VmRefCell, VoRefCell, McRefCell,
                    thetaSRefCell, psieRefCell, bRefCell);

                // Only direct sunlight carries a slope/aspect-dependent timing
                // signal. On diffuse-dominated days progressively suppress the
                // solar time warp and retain the reference clock-hour cycle.
                std::vector<double> dayDirectFrac(nDays);
                for (int d = 0; d < nDays; ++d) {
                    double sum = 0.0; int cnt = 0;
                    for (int h = 0; h < 24; ++h) {
                        int t = d * 24 + h;
                        if (RswC[t] > 5.0) { sum += 1.0 - RdifC[t] / RswC[t]; ++cnt; }
                    }
                    dayDirectFrac[d] = (cnt > 0) ? (sum / cnt) : 0.0;
                }

                // The reference's side of the scaling, evaluated exactly as this
                // cell's side is below: its zero-flux surface temperature from the
                // same one-shot balance, on its own absorbed ground radiation,
                // grid-evaluated resistance and soil, then its exchange
                // coefficient, conductivity and damping depth. A cell identical to
                // the reference therefore scales by exactly one.
                std::vector<double> Tg0refC(nHours), hExchRef(nHours);
                for (int t = 0; t < nHours; ++t) {
                    double hrR = pointmodel::soilrelhumCpp(spRef, TaC[t], theta0refC[t]);
                    if (hrR < 0.0) hrR = 0.0;
                    if (hrR > 1.0) hrR = 1.0;
                    double rgR = rGzrefC[t];
                    Tg0refC[t] = utils::penmanMonteithCpp(RabsGrefC[t], TaC[t], pkC[t], rhC[t],
                        emGrefC[t], rgR, rgR, /*Ts=*/TaC[t], /*G=*/0.0, hrR);
                    hExchRef[t] = surfaceExchangeCoefCpp(Tg0refC[t], TaC[t], pkC[t], emGrefC[t], rgR, hrR);
                }
                std::vector<Harmonic> Href = dailyHarmonicsCpp(Tg0refC, nDays);
                std::vector<double> ksRefDay(nDays), DDRefDay(nDays);
                for (int d = 0; d < nDays; ++d) {
                    double thetaSum = 0.0, TcSum = 0.0;
                    for (int h = 0; h < 24; ++h) {
                        int t = d * 24 + h;
                        thetaSum += theta0refC[t];
                        TcSum += Tg0refC[t];
                    }
                    double thetaMean = thetaSum / 24.0, TcMean = TcSum / 24.0;
                    ksRefDay[d] = utils::thermalConductivityCpp(VqRefCell, VmRefCell, VoRefCell, thetaMean,
                        McRefCell, TcMean, pkC[d * 24]);
                    DDRefDay[d] = pointmodel::diurnalDampingDepthCpp(spRef, thetaMean, TcMean, pkC[d * 24]);
                }

                // This cell's own solar position -- a genuine per-cell
                // series where location varies spatially (in which case the
                // reference's own flat-surface curve is also recomputed
                // from that same per-cell location, so both sides of the
                // time-warp stay geographically consistent with each
                // other), otherwise the shared series already computed
                // above.
                SolarSeries cellSolar;
                const SolarSeries* solar = &sharedSolar;
                std::vector<Crossings> crossRCell;
                const std::vector<Crossings>* crossRP = &crossRShared;
                if (latSpatial) {
                    cellSolar = computeSolarSeriesCpp(lats(i, j), lons(i, j), year, month, day, hour, nHours);
                    solar = &cellSolar;
                    std::vector<double> siRefCell(nHours);
                    for (int t = 0; t < nHours; ++t) {
                        siRefCell[t] = utils::solarindexCpp(0.0, 0.0, cellSolar.zend[t], cellSolar.azid[t]);
                    }
                    crossRCell.resize(nDays);
                    for (int d = 0; d < nDays; ++d) {
                        std::vector<double> siRefDay(siRefCell.begin() + d * 24, siRefCell.begin() + d * 24 + 24);
                        crossRCell[d] = findCrossingsCpp(siRefDay);
                    }
                    crossRP = &crossRCell;
                }

                // This cell's own solar-index curve (real slope/aspect) and
                // per-day rise/set crossings, for the landmark time-warp
                // below.
                double cSlope = slope(i, j), cAspect = aspect(i, j);
                std::vector<double> siTgt(nHours);
                for (int t = 0; t < nHours; ++t) {
                    siTgt[t] = utils::solarindexCpp(cSlope, cAspect, solar->zend[t], solar->azid[t]);
                }
                std::vector<Crossings> crossT(nDays);
                for (int d = 0; d < nDays; ++d) {
                    std::vector<double> siTgtDay(siTgt.begin() + d * 24, siTgt.begin() + d * 24 + 24);
                    crossT[d] = findCrossingsCpp(siTgtDay);
                }

                // Each surface's emission coefficient comes from the radiation
                // step that supplied its absorbed longwave (`emGround`,
                // `emCanopy`; see GridRadlwabsStepCpp), per cell and hour.

                // First surface-energy-balance pass: the temperature this cell
                // would reach if no heat were conducted into or out of the soil.
                std::vector<double> Tground0(nHours), hExch(nHours);
                for (int t = 0; t < nHours; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    double sm = soilm[idx];
                    // Relative humidity is a physical fraction and must lie
                    // in [0,1] -- clamped defensively, since the underlying
                    // Kelvin-equation calculation can in principle overshoot
                    // slightly at the numerical edge of its valid range.
                    double hr = pointmodel::soilrelhumCpp(sp, TaC[t], sm);
                    if (hr < 0.0) hr = 0.0;
                    if (hr > 1.0) hr = 1.0;
                    double rabs = RabsGround[idx];
                    double rg = rGz[idx];
                    double emEffG = emGround[idx];
                    Tground0[t] = utils::penmanMonteithCpp(rabs, TaC[t], pkC[t], rhC[t],
                        emEffG, rg, rg, /*Ts=*/TaC[t], /*G=*/0.0, hr);
                    hExch[t] = surfaceExchangeCoefCpp(Tground0[t], TaC[t], pkC[t], emEffG, rg, hr);
                }

                // Scale the reference heat-flux cycle by the ratio of the daily
                // flux swings the two surfaces are expected to sustain. A zero-flux
                // surface temperature swing A0 drives a flux swing
                // A0 h Y / sqrt(h^2 + sqrt(2) h Y + Y^2), where h is the surface's
                // linearised exchange coefficient and Y = sqrt(2) k / D_D the soil's
                // thermal admittance: the flux itself damps the surface swing.
                std::vector<Harmonic> Htgt = dailyHarmonicsCpp(Tground0, nDays);
                std::vector<double> muGday(nDays);
                for (int d = 0; d < nDays; ++d) {
                    double thetaSum = 0.0, hSum = 0.0, hRefSum = 0.0;
                    for (int h = 0; h < 24; ++h) {
                        thetaSum += soilm[i + rows * j + rows * cols * (d * 24 + h)];
                        hSum += hExch[d * 24 + h];
                        hRefSum += hExchRef[d * 24 + h];
                    }
                    double thetaMean = thetaSum / 24.0;
                    double TcMean = 0.0;
                    for (int h = 0; h < 24; ++h) TcMean += Tground0[d * 24 + h];
                    TcMean /= 24.0;
                    double ksTgt = utils::thermalConductivityCpp(cVq, cVm, cVo, thetaMean, cMc,
                        TcMean, pkC[d * 24]);
                    double DDTgt = pointmodel::diurnalDampingDepthCpp(sp, thetaMean, TcMean, pkC[d * 24]);

                    double gTgt = fluxSwingCpp(Htgt[d].amp, hSum / 24.0, std::sqrt(2.0) * ksTgt / DDTgt);
                    double gRef = fluxSwingCpp(Href[d].amp, hRefSum / 24.0, std::sqrt(2.0) * ksRefDay[d] / DDRefDay[d]);
                    double muG = (gRef >= G_SWING_MIN && std::isfinite(gTgt)) ? gTgt / gRef : 1.0;
                    if (muG < MU_G_MIN) muG = MU_G_MIN;
                    muGday[d] = muG;
                }

                // Transfer the reference ground heat flux to this cell. Only its
                // 24-hour harmonic is the sun-driven cycle the factor and the
                // illumination timing describe, so only that part is scaled and
                // moved to the cell's own timing; the rest, the daily mean and
                // weather-driven departures, is kept at clock time unscaled.
                std::vector<Harmonic> HGref = dailyHarmonicsCpp(GrefC, nDays);
                std::vector<double> Gtgt(nHours);
                for (int d = 0; d < nDays; ++d) {
                    int base = d * 24;
                    const double w24 = 2.0 * pi / 24.0;
                    const double ampG = HGref[d].amp, phG = HGref[d].phaseHr;
                    for (int h = 0; h < 24; ++h) {
                        double tquery = h;
                        if (crossT[d].valid && (*crossRP)[d].valid) {
                            double tqueryFull = warpTimeCpp((double)h, crossT[d].rise, crossT[d].set,
                                (*crossRP)[d].rise, (*crossRP)[d].set);
                            tquery = dayDirectFrac[d] * tqueryFull + (1.0 - dayDirectFrac[d]) * h;
                        }
                        double cycleClock = ampG * std::cos(w24 * (h - phG));
                        double cycleCell = ampG * std::cos(w24 * (tquery - phG));
                        double g = GrefC[base + h] - cycleClock + muGday[d] * cycleCell;
                        if (g > G_ABS_MAX) g = G_ABS_MAX;
                        if (g < -G_ABS_MAX) g = -G_ABS_MAX;
                        Gtgt[base + h] = g;
                    }
                }

                // Second surface-energy-balance pass: include the estimated
                // ground heat flux to obtain the cell's final surface temperature.
                for (int t = 0; t < nHours; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    double sm = soilm[idx];
                    double hr = pointmodel::soilrelhumCpp(sp, TaC[t], sm);
                    if (hr < 0.0) hr = 0.0;
                    if (hr > 1.0) hr = 1.0;
                    double rabs = RabsGround[idx];
                    double rg = rGz[idx];
                    double rm = rGm[idx];
                    // The air the ground exchanges with. Canopy heat and vapour
                    // enter the column above the ground; with the sources spread
                    // uniformly with height the ground balance is exactly one
                    // with air raised by the whole surface's fluxes across the
                    // rise resistance rg - rm, exchanged across rm. Those fluxes
                    // come from the canopy balance canopyTempCpp solves, with the
                    // same inputs, so the pass stays closed form. Without a
                    // column rm equals rg and the air is the reference air.
                    double taG = TaC[t], rhG = rhC[t];
                    double rRise = rg - rm;
                    // With foliage, the bulk surface and the ground are solved
                    // together. Their two energy balances, each linearised about
                    // air temperature as the Penman-Monteith step does, are two
                    // linear equations in the two surface temperatures:
                    //   bulk:   Lc dTc = Fc - Kc (1 - h) hr sg dTg
                    //   ground: Lg dTg = Fg + beta Lg dTc
                    // The soil's vapour is at the ground's own temperature; the
                    // ground's air is raised by the bulk surface's sensible flux
                    // across the node's resistance, and its vapour boundary is the
                    // network seen from the soil. Substituting the second into the
                    // first gives dTc in closed form. Two passes, the second
                    // linearised about the first, as the bulk step elsewhere makes.
                    // The node is at whichever of the mean source height and the
                    // exchange surface is nearer the reference.
                    double rhaC = rHa[idx];
                    double rvC = rhaC + (rSurfSpatial ? rSurf_ref[idx] : rSurf_ref[t]);
                    double hF = rSurfSpatial ? hFol_ref[idx] : hFol_ref[t];
                    double rsh = (rRise > 0.0) ? std::min(rRise, rhaC) : 0.0;
                    bool coupled = rsh > 0.0 && std::isfinite(hF) && rvC - rsh > 0.0;
                    double ts;
                    if (coupled) {
                        const double Tair = TaC[t], pkt = pkC[t];
                        const double aN = hF / (rvC - rsh), cN = 1.0 / rsh;
                        const double rleg = rg - rsh, rTh = rleg + 1.0 / (aN + cN);
                        const double ea = utils::satvapCpp(Tair) * rhC[t] / 100.0;
                        const double esA = utils::satvapCpp(Tair);
                        const double phA = utils::phairCpp(Tair, pkt), cpA = utils::cpairCpp(Tair);
                        const double emC = emCanopy[idx], emG = emGround[idx], rabsC = RabsCanopy[idx];
                        const double TkA4 = utils::rademCpp(Tair), g = Gtgt[t];
                        auto slope = [](double te) { return utils::satvapCpp(te + 0.5) - utils::satvapCpp(te - 0.5); };
                        auto solvePair = [&](double gc, double gg, double& Tc, double& Tg) {
                            double Tec = 0.5 * (gc + Tair), Teg = 0.5 * (gg + Tair);
                            double sc = slope(Tec), sg = slope(Teg);
                            double Kc = utils::latentHeatCpp(gc) * phA / (pkt * rvC);
                            double Tsg = Tair + (gc - Tair) * rsh / rhaC;
                            double phG = utils::phairCpp(Tsg, pkt), cpG = utils::cpairCpp(Tsg);
                            double Kg = utils::latentHeatCpp(gg) * phG / (pkt * rTh);
                            double Lc = 4.0 * emC * sb * std::pow(Tec + 273.15, 3.0) + phA * cpA / rhaC + Kc * hF * sc;
                            double Fc = rabsC - emC * sb * TkA4 - g - Kc * (hF * esA + (1.0 - hF) * hr * esA - ea);
                            double Lg = 4.0 * emG * sb * std::pow(Teg + 273.15, 3.0) + phG * cpG / rleg + Kg * hr * sg;
                            double eTh0 = (aN * esA + cN * ea) / (aN + cN);
                            double Fg = rabs - emG * sb * TkA4 - g - Kg * (hr * esA - eTh0);
                            double beta = (phG * cpG * rsh / (rleg * rhaC) + Kg * aN * sc / (aN + cN)) / Lg;
                            double q = Kc * (1.0 - hF) * hr * sg;
                            double dTc = (Fc - q * Fg / Lg) / (Lc + q * beta);
                            Tc = Tair + dTc;
                            Tg = Tair + Fg / Lg + beta * dTc;
                        };
                        double Tc1, Tg1, Tc2, Tg2;
                        solvePair(Tair, Tair, Tc1, Tg1);
                        solvePair(Tc1, Tg1, Tc2, Tg2);
                        ts = Tg2;
                        Tc_est[idx] = boundSurfaceTempCpp(Tc2, Tair, rhC[t], pkt, rabsC, emC);
                    } else {
                    // Otherwise the ground alone, against air raised by the bulk
                    // surface's fluxes at the same node; without a column, reference air.
                    if (rRise > rhaC) { rRise = rhaC; rm = rg - rRise; }
                    if (rRise > 0.0) {
                        double rSurf_t = rSurfSpatial ? rSurf_ref[idx] : rSurf_ref[t];
                        double hSurf_t = rSurfSpatial ? hSurf_ref[idx] : hSurf_ref[t];
                        double rha = rHa[idx];
                        double rv = rha + rSurf_t;
                        double tc1 = utils::penmanMonteithCpp(RabsCanopy[idx], TaC[t], pkC[t], rhC[t],
                            emCanopy[idx], rha, rv, /*Ts=*/TaC[t], Gtgt[t], hSurf_t);
                        double tc2 = utils::penmanMonteithCpp(RabsCanopy[idx], TaC[t], pkC[t], rhC[t],
                            emCanopy[idx], rha, rv, /*Ts=*/tc1, Gtgt[t], hSurf_t);
                        tc2 = boundSurfaceTempCpp(tc2, TaC[t], rhC[t], pkC[t], RabsCanopy[idx], emCanopy[idx]);
                        // This is the one-surface bulk temperature, so the cell
                        // needs no separate estimate of it.
                        Tc_one[idx] = tc2;
                        double ph = utils::phairCpp(TaC[t], pkC[t]);
                        double cp = utils::cpairCpp(TaC[t]);
                        double ea = utils::satvapCpp(TaC[t]) * rhC[t] / 100.0;
                        // Sensible and latent heat leaving the surface, the latter
                        // from the linearisation the second canopy pass balanced.
                        double Hs = (ph * cp / rha) * (tc2 - TaC[t]);
                        double Te = (tc1 + TaC[t]) / 2.0;
                        double De = hSurf_t * (utils::satvapCpp(Te + 0.5) - utils::satvapCpp(Te - 0.5));
                        double Ls = (utils::latentHeatCpp(tc1) * ph / (pkC[t] * rv)) *
                            (hSurf_t * utils::satvapCpp(TaC[t]) - ea + De * (tc2 - TaC[t]));
                        taG = TaC[t] + Hs * rRise / (ph * cp);
                        double eG = ea + Ls * rRise * pkC[t] / (utils::latentHeatCpp(TaC[t]) * ph);
                        if (eG < 1e-6) eG = 1e-6;
                        rhG = 100.0 * eG / utils::satvapCpp(taG);
                        rg = rm;
                    }
                    ts = utils::penmanMonteithCpp(rabs, taG, pkC[t], rhG,
                        emGround[idx], rg, rg, /*Ts=*/Tground0[t], Gtgt[t], hr);
                    }
                    if (ts > Tground0[t] + TS_DEVIATION_MAX) ts = Tground0[t] + TS_DEVIATION_MAX;
                    if (ts < Tground0[t] - TS_DEVIATION_MAX) ts = Tground0[t] - TS_DEVIATION_MAX;
                    ts = boundSurfaceTempCpp(ts, TaC[t], rhC[t], pkC[t], rabs, emGround[idx]);

                    G_est[idx] = Gtgt[t];
                    Ts_est[idx] = ts;
                }
            }
        }

        return List::create(Named("G_est") = G_est, Named("Ts_est") = Ts_est, Named("Tcanopy_est") = Tc_est,
            Named("Tcanopy_one") = Tc_one);
    }

    // For an hourly series spanning whole calendar days, computes one
    // summary statistic (mean/max/min) per day and repeats it across that
    // day's own 24 hours.
    enum class DayStat { Mean, Max, Min };
    static std::vector<double> hourToDayCpp(const std::vector<double>& hourly, DayStat stat) {
        int n = (int)hourly.size();
        int nDays = n / 24;
        std::vector<double> daily(n);
        for (int d = 0; d < nDays; ++d) {
            double v = hourly[d * 24];
            for (int h = 1; h < 24; ++h) {
                double x = hourly[d * 24 + h];
                if (stat == DayStat::Mean) v += x;
                else if (stat == DayStat::Max) { if (x > v) v = x; }
                else { if (x < v) v = x; }
            }
            if (stat == DayStat::Mean) v /= 24.0;
            for (int h = 0; h < 24; ++h) daily[d * 24 + h] = v;
        }
        return daily;
    }

    // Estimate temperature below ground by transferring the reference point
    // model's resolved subsurface response to each cell. Very shallow depths
    // follow the cell's own surface temperature; with increasing depth the
    // diurnal anomaly is progressively taken from the reference soil profile,
    // scaled to the cell's surface diurnal range and mean. At depths beyond
    // the diurnal signal the estimate transitions toward mean annual
    // temperature using the annual thermal damping scale.
    NumericVector belowGroundShortcutCpp(
        NumericVector Ts_est, NumericVector Tgp_ref, NumericVector Tbp_ref,
        NumericVector DD, NumericVector mat, double reqhgt, double hiy)
    {
        if (reqhgt >= 0.0) {
            Rcpp::stop("belowGroundShortcutCpp: reqhgt must be negative (a below-ground depth)");
        }
        IntegerVector dim = Ts_est.attr("dim");
        int rows = dim[0], cols = dim[1], tsteps = dim[2];
        if (tsteps % 24 != 0) {
            Rcpp::stop("belowGroundShortcutCpp: tsteps (%d) is not a multiple of 24 -- the "
                "daily min/max/mean blend needs complete real calendar-day blocks", tsteps);
        }
        bool refSpatial = (Tgp_ref.size() != tsteps);

        int n = rows * cols * tsteps;
        NumericVector Tz_est(n, NA_REAL);
        IntegerVector outDim = { rows, cols, tsteps };
        Tz_est.attr("dim") = outDim;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                std::vector<double> Tg(tsteps);
                bool anyNA = false;
                for (int t = 0; t < tsteps; ++t) {
                    double v = Ts_est[i + rows * j + rows * cols * t];
                    if (NumericVector::is_na(v)) { anyNA = true; break; }
                    Tg[t] = v;
                }
                if (anyNA) continue; // non-land cell, matching Ts_est's own NA convention

                double matCell = (mat.size() == 1) ? mat[0] : mat[i + rows * j];

                std::vector<double> TgMax = hourToDayCpp(Tg, DayStat::Max);
                std::vector<double> TgMin = hourToDayCpp(Tg, DayStat::Min);
                std::vector<double> TgMean = hourToDayCpp(Tg, DayStat::Mean);

                // Separate each reference day into its mean and within-day
                // anomaly so the anomaly can be rescaled to this cell's own
                // surface range without losing the reference soil phase lag.
                std::vector<double> Tgp(tsteps), Tbp(tsteps);
                for (int t = 0; t < tsteps; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    Tgp[t] = refSpatial ? Tgp_ref[idx] : Tgp_ref[t];
                    Tbp[t] = refSpatial ? Tbp_ref[idx] : Tbp_ref[t];
                }
                std::vector<double> TbpDaily = hourToDayCpp(Tbp, DayStat::Mean);
                std::vector<double> Tbpa(tsteps);
                for (int t = 0; t < tsteps; ++t) Tbpa[t] = Tbp[t] - TbpDaily[t];
                std::vector<double> TgpMax = hourToDayCpp(Tgp, DayStat::Max);
                std::vector<double> TgpMin = hourToDayCpp(Tgp, DayStat::Min);
                std::vector<double> TgpMean = hourToDayCpp(Tgp, DayStat::Mean);

                for (int t = 0; t < tsteps; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    // Express depth relative to this cell's diurnal and annual
                    // thermal damping scales for the day. These determine whether
                    // local surface forcing, the reference subsurface anomaly, or
                    // mean annual temperature dominates at the requested depth.
                    double DDc = DD[idx];
                    double tz;
                    if (!(DDc > 0.0)) {
                        tz = NA_REAL;
                    } else {
                        double nb = -118.35 * reqhgt / DDc;
                        double nb_ann = -118.35 * reqhgt / (DDc * std::sqrt(365.0));
                        if (nb <= 1.0) {
                            tz = Tg[t];
                        } else {
                            // Scale the reference subsurface anomaly by the ratio of
                            // target-to-reference surface diurnal range; if the
                            // reference surface has no daily range, no anomaly is transferred.
                            double refRange = TgpMax[t] - TgpMin[t];
                            double rat = (refRange != 0.0) ? (TgMax[t] - TgMin[t]) / refRange : 0.0;
                            double dif = TgMean[t] - TgpMean[t];
                            double Tzd = rat * Tbpa[t] + TbpDaily[t] + dif;
                            if (nb <= 24.0) {
                                double w1 = 1.0 / nb, w2 = nb / 24.0;
                                double wgt = w1 / (w1 + w2);
                                tz = wgt * Tg[t] + (1.0 - wgt) * Tzd;
                            } else if (nb_ann < hiy) {
                                double w1 = 24.0 / nb_ann, w2 = nb_ann / hiy;
                                double wgt = w1 / (w1 + w2);
                                tz = wgt * Tzd + (1.0 - wgt) * matCell;
                            } else {
                                tz = matCell;
                            }
                        }
                    }
                    Tz_est[idx] = tz;
                }
            }
        }
        return Tz_est;
    }

    // Estimate the bulk canopy+ground exchange-surface temperature for each
    // cell from its absorbed radiation, aerodynamic resistance and ground
    // heat flux. The `_ref`-named resistance arguments are the cell's own
    // values from surfaceResistGridCpp's vapour network, supplied rather than
    // solved here; svfa is used only to mask non-land cells. The resulting
    // temperature is bounded by the empirical surface-air envelope and the
    // physical dew-point floor in boundSurfaceTempCpp.
    NumericVector canopyTempCpp(
        NumericVector RabsCanopy, NumericVector rHa, NumericVector G_est,
        NumericVector Ta, NumericVector rh, NumericVector pk,
        NumericMatrix svfa, NumericVector emCanopy,
        NumericVector rSurf_ref, NumericVector hSurf_ref)
    {
        IntegerVector dim = RabsCanopy.attr("dim");
        int rows = dim[0], cols = dim[1], tsteps = dim[2];
        bool climSpatial = (Ta.size() != tsteps);
        bool rSurfSpatial = (rSurf_ref.size() != tsteps);

        int n = rows * cols * tsteps;
        NumericVector Tc_est(n, NA_REAL);
        IntegerVector outDim = { rows, cols, tsteps };
        Tc_est.attr("dim") = outDim;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (NumericMatrix::is_na(svfa(i, j))) continue; // non-land cell

                for (int t = 0; t < tsteps; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    double rabs = RabsCanopy[idx];
                    double rha_cell = rHa[idx];
                    if (NumericVector::is_na(rabs) || NumericVector::is_na(rha_cell)) continue;

                    double Ta_t = climSpatial ? Ta[idx] : Ta[t];
                    double rh_t = climSpatial ? rh[idx] : rh[t];
                    double pk_t = climSpatial ? pk[idx] : pk[t];
                    double rSurf_t = rSurfSpatial ? rSurf_ref[idx] : rSurf_ref[t];
                    double hSurf_t = rSurfSpatial ? hSurf_ref[idx] : hSurf_ref[t];

                    double rv = rha_cell + rSurf_t;
                    double g = G_est[idx];

                    // Two successive Penman-Monteith evaluations reduce the
                    // error from linearising emitted longwave around an initial
                    // surface-temperature guess.
                    double tc1 = utils::penmanMonteithCpp(rabs, Ta_t, pk_t, rh_t,
                        emCanopy[idx], rha_cell, rv, /*Ts=*/Ta_t, g, hSurf_t);
                    double tc2 = utils::penmanMonteithCpp(rabs, Ta_t, pk_t, rh_t,
                        emCanopy[idx], rha_cell, rv, /*Ts=*/tc1, g, hSurf_t);

                    // The same backstop as the ground's: the canopy surface's
                    // one-shot temperature cannot know the cell's own heating.
                    Tc_est[idx] = boundSurfaceTempCpp(tc2, Ta_t, rh_t, pk_t, rabs, emCanopy[idx]);
                }
            }
        }
        return Tc_est;
    }

    // ========================================================================
    // Leaf temperature at the requested canopy height
    // ========================================================================
    // Solve the energy balance of an individual leaf, distinct from the bulk
    // canopy exchange surface above. The leaf receives the shortwave and
    // longwave radiation at its own height, exchanges sensible heat through a
    // leaf-scale boundary layer, and loses latent heat according to stomatal
    // conductance from the same photosynthesis-hydraulic optimisation used by
    // the point model. PAR drives stomatal physiology but is not added a second
    // time to the energy balance. Output is defined only where the requested
    // height lies within vegetation.

    NumericVector leafTempCpp(
        NumericVector radLpar, NumericVector radLsw, NumericVector radLlw, NumericVector uz,
        NumericVector Ta, NumericVector rh, NumericVector pk, NumericVector Ca, NumericVector psi_r,
        NumericMatrix hgt, NumericMatrix pai,
        NumericMatrix Vcmax25, NumericMatrix Tup, NumericMatrix Tlw, NumericMatrix Dcrit,
        NumericMatrix alpha, NumericMatrix f0, NumericMatrix fd, NumericMatrix psi50,
        NumericMatrix apsi, NumericMatrix rpmin, NumericMatrix leafd, LogicalMatrix isC3,
        double reqhgt)
    {
        IntegerVector dim = radLpar.attr("dim");
        int rows = dim[0], cols = dim[1], tsteps = dim[2];
        // Atmospheric and plant-water drivers may either be common to the
        // domain or vary from cell to cell; each leaf therefore uses the
        // local value where one is available and otherwise the shared forcing.
        bool taSpatial = (Ta.size() != tsteps);
        bool rhSpatial = (rh.size() != tsteps);
        bool pkSpatial = (pk.size() != tsteps);
        bool CaSpatial = (Ca.size() != tsteps);
        bool psirSpatial = (psi_r.size() != tsteps);

        int n = rows * cols * tsteps;
        NumericVector Tleaf(n, NA_REAL);
        IntegerVector outDim = { rows, cols, tsteps };
        Tleaf.attr("dim") = outDim;

        if (!(reqhgt > 0.0)) return Tleaf; // below-ground/surface: no leaf, NA throughout

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double hgtCell = hgt(i, j);
                if (NumericMatrix::is_na(hgtCell) || !(hgtCell > 0.0) || !(reqhgt < hgtCell)) continue; // bare ground, above canopy, or no land

                pointmodel::vegpstruct vp{};
                vp.Vcmax25 = Vcmax25(i, j);
                vp.Tup = Tup(i, j);
                vp.Tlw = Tlw(i, j);
                vp.Dcrit = Dcrit(i, j);
                vp.alpha = alpha(i, j);
                vp.f0 = f0(i, j);
                vp.fd = fd(i, j);
                vp.psi50 = psi50(i, j);
                vp.apsi = apsi(i, j);
                vp.rpmin = rpmin(i, j);
                vp.gsmaxCap = -1.0; // no extra empirical cap -- theoretical Vcmax25-derived ceiling only
                double leafdCell = leafd(i, j);
                bool C3 = isC3(i, j);

                for (int t = 0; t < tsteps; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    double parCell = radLpar[idx], swCell = radLsw[idx], lwCell = radLlw[idx], uzCell = uz[idx];
                    if (NumericVector::is_na(parCell) || NumericVector::is_na(swCell) ||
                        NumericVector::is_na(lwCell) || NumericVector::is_na(uzCell) || !(uzCell > 0.0)) continue;

                    double Ta_t = taSpatial ? Ta[idx] : Ta[t];
                    double rh_t = rhSpatial ? rh[idx] : rh[t];
                    double pk_t = pkSpatial ? pk[idx] : pk[t];
                    double Ca_t = CaSpatial ? Ca[idx] : Ca[t];
                    double psir_t = psirSpatial ? psi_r[idx] : psi_r[t];

                    double rHa = 190.9 * std::sqrt(leafdCell / uzCell);
                    if (rHa > 200.0) rHa = 200.0;

                    double Rabs = swCell + lwCell;

                    pointmodel::envstruct env{};
                    env.tair = Ta_t; env.rh = rh_t; env.pk = pk_t; env.Ca = Ca_t; env.psi_r = psir_t;
                    env.PARabs = parCell;

                    env.tcanopy = Ta_t;
                    double gs1 = pointmodel::leafgsCpp(env, vp, reqhgt, C3);
                    double ph = utils::phairCpp(Ta_t, pk_t);
                    // Closed stomata (no light, e.g. night-time) are given
                    // a large but finite resistance rather than infinite,
                    // a conventional order-of-magnitude value for
                    // cuticular/closed-stomata resistance. Conductance falls
                    // smoothly to zero as light fades, so that value also
                    // caps the resistance in dim light; without the cap a leaf
                    // in dim light would lose less water than one in the dark.
                    const double rStomMax = mc::rStomClosed;
                    double rStom1 = (gs1 > 0.0) ? std::min(ph / gs1, rStomMax) : rStomMax;
                    double rV1 = rHa + rStom1;
                    double tleaf1 = utils::penmanMonteithCpp(Rabs, Ta_t, pk_t, rh_t,
                        surfaceEmissivity, rHa, rV1, /*Ts=*/Ta_t, /*G=*/0.0);

                    env.tcanopy = tleaf1;
                    double gs2 = pointmodel::leafgsCpp(env, vp, reqhgt, C3);
                    double rStom2 = (gs2 > 0.0) ? std::min(ph / gs2, rStomMax) : rStomMax;
                    double rV2 = rHa + rStom2;
                    double tleaf2 = utils::penmanMonteithCpp(Rabs, Ta_t, pk_t, rh_t,
                        surfaceEmissivity, rHa, rV2, /*Ts=*/tleaf1, /*G=*/0.0);

                    Tleaf[idx] = tleaf2;
                }
            }
        }
        return Tleaf;
    }

    // Longwave radiation at the requested height: the flux arriving from above
    // and the flux arriving from below, as a radiometer at that height would
    // read them. Downward is open sky over `svfa` of the upper hemisphere and
    // surrounding terrain, at air temperature, over the rest, filtered by any
    // foliage above the height. Upward is the ground seen through the foliage
    // below the height, with that foliage's own emission filling the rest.
    // Emission is `em*sigma*T^4` whatever stands above it, so no sky-view
    // factor multiplies an emitted term here.
    //
    // Tcanopy_est is the bulk surface of canopy and ground together, so its
    // emission L_c is the whole surface's, reported as the upward flux at and
    // above canopy top. The foliage's own emission follows from the bulk
    // surface being foliage over 1 - exp(-P) and ground through the gaps,
    // L_c = exp(-P) L_g + (1 - exp(-P)) L_f. A layer of foliage transmitting
    // tau therefore contributes (1 - tau) L_f = w (L_c - exp(-P) L_g), with
    // w = (1 - tau) / (1 - exp(-P)) between 0 and 1, the form evaluated here so
    // that nothing is divided by a vanishing foliage fraction. The upward flux
    // then runs from L_g at the ground to exactly L_c at canopy top.
    List longwaveGridCpp(
        NumericVector Ts_est, NumericVector Tcanopy_est, NumericVector Ta,
        NumericMatrix hgt, NumericMatrix pai, NumericMatrix paia,
        NumericMatrix svfa, NumericVector lwdown, double reqhgt)
    {
        IntegerVector dim = Ts_est.attr("dim");
        int rows = dim[0], cols = dim[1], tsteps = dim[2];
        bool lwdownSpatial = (lwdown.size() != tsteps);
        bool TaSpatial = (Ta.size() != tsteps);

        int n = rows * cols * tsteps;
        NumericVector Rlwdown(n, NA_REAL), Rlwup(n, NA_REAL);
        IntegerVector outDim = { rows, cols, tsteps };
        Rlwdown.attr("dim") = outDim;
        Rlwup.attr("dim") = outDim;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double hgtCell = hgt(i, j);
                if (NumericMatrix::is_na(hgtCell)) continue; // non-land cell

                double svfaCell = svfa(i, j);
                if (NumericMatrix::is_na(svfaCell)) continue;

                double paiCell = pai(i, j);
                bool aboveCanopy = (reqhgt >= hgtCell) || !(paiCell > 0.0);
                double paiaCell = aboveCanopy ? 0.0 : std::min(std::max(paia(i, j), 0.0), paiCell);
                double tr = 0.0, trb = 0.0, trAll = 0.0, wa = 0.0, wb = 0.0;
                if (!aboveCanopy) {
                    // Plant area above and below the requested height attenuates
                    // the view of sky and ground respectively; the obscured
                    // fraction is filled by the foliage's own emission.
                    tr = std::exp(-paiaCell);
                    trb = std::exp(-(paiCell - paiaCell));
                    trAll = std::exp(-paiCell);
                    double folAll = -std::expm1(-paiCell);
                    wa = -std::expm1(-paiaCell) / folAll;
                    wb = -std::expm1(-(paiCell - paiaCell)) / folAll;
                }

                for (int t = 0; t < tsteps; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    double tsCell = Ts_est[idx], tcCell = Tcanopy_est[idx];
                    if (NumericVector::is_na(tsCell) || NumericVector::is_na(tcCell)) continue;

                    double lwdown_t = lwdownSpatial ? lwdown[idx] : lwdown[t];
                    double Ta_t = TaSpatial ? Ta[idx] : Ta[t];
                    // Above the canopy: sky over the open fraction, terrain
                    // over the blocked fraction.
                    double lwSky = svfaCell * lwdown_t +
                        (1.0 - svfaCell) * surfaceEmissivity * sb * utils::rademCpp(Ta_t);
                    double lwcan = surfaceEmissivity * sb * utils::rademCpp(tcCell);
                    if (aboveCanopy) {
                        Rlwdown[idx] = lwSky;
                        Rlwup[idx] = lwcan;
                    } else {
                        double lwgro = surfaceEmissivity * sb * utils::rademCpp(tsCell);
                        double lwfol = lwcan - trAll * lwgro; // (1 - exp(-P)) L_f
                        Rlwdown[idx] = tr * lwSky + wa * lwfol;
                        Rlwup[idx] = trb * lwgro + wb * lwfol;
                    }
                }
            }
        }
        return List::create(Named("Rlwdown") = Rlwdown, Named("Rlwup") = Rlwup);
    }

    // Extend the estimated canopy exchange-surface state to a requested
    // height above the canopy, for heat and moisture, through the transport
    // column's resistances where there is a canopy (roughness sublayer
    // included) and Monin-Obukhov similarity over bare ground. Temperature
    // and vapour pressure move continuously
    // from their surface values toward the reference atmosphere according to
    // the fraction of the full aerodynamic resistance accumulated by reqhgt.
    List aboveCanopyProfileGridCpp(
        NumericVector Tcanopy_est, NumericVector rHa, NumericVector L, NumericVector uf,
        NumericMatrix d, NumericMatrix zh,
        NumericVector Ta, NumericVector rh, NumericVector rSurf_ref, NumericVector hSurf_ref,
        NumericVector rGz, NumericVector rGreq, NumericMatrix zs,
        double reqhgt, double zref)
    {
        IntegerVector dim = Tcanopy_est.attr("dim");
        int rows = dim[0], cols = dim[1], tsteps = dim[2];
        bool climSpatial = (Ta.size() != tsteps);
        bool rSurfSpatial = (rSurf_ref.size() != tsteps);

        int n = rows * cols * tsteps;
        NumericVector Tz(n, NA_REAL), RHz(n, NA_REAL);
        IntegerVector outDim = { rows, cols, tsteps };
        Tz.attr("dim") = outDim;
        RHz.attr("dim") = outDim;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double dCell = d(i, j), zhCell = zh(i, j);
                if (NumericMatrix::is_na(dCell) || NumericMatrix::is_na(zhCell)) continue;
                // The top of this cell's roughness sublayer. NA where there is no
                // canopy and so no transport column, in which case the ordinary
                // logarithmic profile from the roughness length for heat is the
                // whole answer.
                double zsCell = zs(i, j);
                bool hasCol = !NumericMatrix::is_na(zsCell);

                for (int t = 0; t < tsteps; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    double Tsurface = Tcanopy_est[idx];
                    double rZref = rHa[idx];   // canopy-to-zref resistance, already computed by windCpp
                    double Lc = L[idx];
                    double ufc = uf[idx];
                    if (NumericVector::is_na(Tsurface) || NumericVector::is_na(rZref) ||
                        NumericVector::is_na(Lc) || NumericVector::is_na(ufc)) continue;

                    double Ta_t = climSpatial ? Ta[idx] : Ta[t];
                    double rh_t = climSpatial ? rh[idx] : rh[t];
                    double rSurf_t = rSurfSpatial ? rSurf_ref[idx] : rSurf_ref[t];
                    double hSurf_t = rSurfSpatial ? hSurf_ref[idx] : hSurf_ref[t];

                    // The resistance accumulated between the exchange surface
                    // and reqhgt determines how far the air has moved from the
                    // surface state toward the reference atmosphere. Above canopy
                    // top the exchange-to-reference resistance is the surface's
                    // own segment to canopy top plus the transport column above
                    // it, so the resistance still to cross between reqhgt and the
                    // reference height is the column's, at every height, inside
                    // the roughness sublayer and above it.
                    double ratio;
                    if (hasCol) {
                        double rc_zR = rGz[idx], rc_z = rGreq[idx];
                        if (NumericVector::is_na(rc_zR) || NumericVector::is_na(rc_z) || !(rZref > 0.0)) continue;
                        ratio = 1.0 - (rc_zR - rc_z) / rZref;
                    }
                    else {
                        ratio = utils::rHaToHeightScalarCpp(reqhgt, dCell, zhCell, Lc, ufc) / rZref;
                    }
                    double tz = Tsurface - (Tsurface - Ta_t) * ratio;

                    // Vapour pressure at the canopy exchange surface blends
                    // the surface value (saturation scaled by the cell's
                    // humidity factor) and the free-air value
                    // according to how much of the canopy's total vapour-
                    // transfer resistance is aerodynamic versus surface-
                    // controlled, then is extrapolated to the requested
                    // height by the same profile ratio used for temperature.
                    double rv = rZref + rSurf_t;
                    double es = 1000.0 * hSurf_t * utils::satvapCpp(Tsurface);
                    double ea = 1000.0 * utils::satvapCpp(Ta_t) * (rh_t / 100.0);
                    double eSurface = ea + (es - ea) * (rZref / rv);
                    double ez = eSurface - (eSurface - ea) * ratio;
                    double esatTz = 1000.0 * utils::satvapCpp(tz);
                    double rhz = 100.0 * ez / esatTz;
                    if (rhz < 0.0) rhz = 0.0;
                    if (rhz > 100.0) rhz = 100.0;

                    Tz[idx] = tz;
                    RHz[idx] = rhz;
                }
            }
        }
        return List::create(Named("Tz") = Tz, Named("RHz") = RHz);
    }

    // ========================================================================
    // Below-canopy air temperature and humidity
    // ========================================================================
    // A far-field canopy-scale gradient, constrained by both ground and
    // canopy-top conditions. No near-field term is added: it depends on the
    // vertical distribution of foliage, which canopy height and plant area
    // do not determine.

    // Air temperature and humidity at a height inside the canopy, for every
    // cell of the grid. The profile is the far-field reconstruction: the
    // canopy source is uniform with height, the column
    // supplies every resistance, and the two boundary values are the ground
    // surface and the air at the reference height. Its vertical shape at the
    // requested height is fixed by canopy geometry, so it is taken once per
    // cell from the neutral column and rescaled each hour by that hour's own
    // resistance totals.
    //
    // There is no blend towards a bare-ground profile and no bound on the
    // result: the column is already exact at zero plant area, and the profile
    // returns the ground surface at the column floor and the canopy-top air at
    // canopy top by construction.
    List belowCanopyProfileGridCpp(
        NumericVector Tcanopy_est, NumericVector Ts_est, NumericVector rBL,
        NumericVector rGh, NumericVector rGm, NumericVector rGz,
        NumericMatrix shapeR, NumericMatrix shapeC, NumericMatrix hgtCell, NumericMatrix paiCell,
        NumericVector Ta, NumericVector rh, NumericVector pk, NumericVector rSurf_ref,
        NumericVector hSurf_ref, double reqhgt, double zref)
    {
        IntegerVector dim = Tcanopy_est.attr("dim");
        int rows = dim[0], cols = dim[1], tsteps = dim[2];
        bool climSpatial = (Ta.size() != tsteps);
        bool rSurfSpatial = (rSurf_ref.size() != tsteps);

        int n = rows * cols * tsteps;
        NumericVector Tz(n, NA_REAL);
        NumericVector RHz(n, NA_REAL);
        IntegerVector outDim = { rows, cols, tsteps };
        Tz.attr("dim") = outDim;
        RHz.attr("dim") = outDim;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double hgt = hgtCell(i, j), pai = paiCell(i, j);
                double fR = shapeR(i, j), fC = shapeC(i, j);
                if (!(hgt > 0.0) || !(pai > 0.0) || !(reqhgt < hgt) ||
                    NumericMatrix::is_na(fR) || NumericMatrix::is_na(fC)) continue;

                for (int t = 0; t < tsteps; ++t) {
                    int idx = i + rows * j + rows * cols * t;
                    double tcEst = Tcanopy_est[idx], tsEst = Ts_est[idx], rbl_cell = rBL[idx];
                    // Every resistance this hour needs was produced by the wind
                    // step, which built the column once: A to canopy top, the
                    // column mean, and the resistance to the reference height.
                    double A = rGh[idx], Rbar = rGm[idx], RzR = (zref > hgt) ? rGz[idx] : rGh[idx];
                    if (NumericVector::is_na(tcEst) || NumericVector::is_na(tsEst) ||
                        NumericVector::is_na(rbl_cell) || NumericVector::is_na(A) ||
                        NumericVector::is_na(Rbar)) continue;

                    double Ta_t = climSpatial ? Ta[idx] : Ta[t];
                    double rh_t = climSpatial ? rh[idx] : rh[t];
                    double pk_t = climSpatial ? pk[idx] : pk[t];
                    double rSurf_t = rSurfSpatial ? rSurf_ref[idx] : rSurf_ref[t];
                    double hSurf_t = rSurfSpatial ? hSurf_ref[idx] : hSurf_ref[t];

                    double tz, rhz;
                    utils::farFieldPairCpp(Ta_t, rh_t, pk_t, tcEst, tsEst, 1.0, rbl_cell, rSurf_t, hSurf_t,
                        A, Rbar, RzR, A * fR, (A - Rbar) * fC, tz, rhz);
                    Tz[idx] = tz;
                    RHz[idx] = rhz;
                }
            }
        }
        return List::create(Named("Tz") = Tz, Named("RHz") = RHz);
    }

    // ========================================================================
    // Below-canopy air temperature and humidity for a point-model run
    // ========================================================================
    // The same profile as the grid's, evaluated for one location's time series.
    // The point model carries an explicit ground-surface humidity, so its
    // lower vapour boundary is generally unsaturated where the grid's is
    // saturated; that is the only difference between the two, and it is an
    // argument rather than a second set of equations. Meaningful only for a
    // height strictly inside a vegetated canopy; other regimes return NA and
    // belong to the above-canopy and surface routines.
    List belowCanopyProfilePointCpp(
        NumericVector Tcanopy, NumericVector Tground, NumericVector groundhr,
        NumericVector rBL, NumericVector L, NumericVector uf,
        double d, double zm, double hgt, double pai,
        NumericVector Ta, NumericVector rh, NumericVector pk, NumericVector rSurf,
        NumericVector hSurf, double reqhgt, double zref)
    {
        int n = Tcanopy.size();
        NumericVector Tz(n, NA_REAL), RHz(n, NA_REAL);
        if (!(hgt > 0.0) || !(pai > 0.0) || !(reqhgt < hgt) || !(reqhgt > 0.0))
            return List::create(Named("Tz") = Tz, Named("RHz") = RHz);

        utils::ColumnStruct col = utils::columnSetupCpp(hgt, pai, d);
        utils::ProfileShape shp = utils::columnProfileShapeCpp(col, reqhgt);

        for (int i = 0; i < n; ++i) {
            double tcEst = Tcanopy[i], tsEst = Tground[i], rbl_i = rBL[i];
            double Lc = L[i], ufc = uf[i];
            if (NumericVector::is_na(tcEst) || NumericVector::is_na(tsEst) ||
                NumericVector::is_na(rbl_i) || NumericVector::is_na(Lc) ||
                NumericVector::is_na(ufc)) continue;

            double A = utils::columnResistCpp(col, ufc, Lc, hgt);
            double Rbar = utils::columnMeanResistCpp(col, ufc, Lc);
            double RzR = (zref > hgt) ? utils::columnResistCpp(col, ufc, Lc, zref) : A;

            double hrGround = groundhr[i];
            if (!(hrGround > 0.0)) hrGround = 1e-6;
            if (hrGround > 1.0) hrGround = 1.0;

            double tz, rhz;
            utils::farFieldPairCpp(Ta[i], rh[i], pk[i], tcEst, tsEst, hrGround, rbl_i, rSurf[i], hSurf[i],
                A, Rbar, RzR, A * shp.fR, (A - Rbar) * shp.fC, tz, rhz);
            Tz[i] = tz;
            RHz[i] = rhz;
        }
        return List::create(Named("Tz") = Tz, Named("RHz") = RHz);
    }

    // ========================================================================
    // Smooth seasonal vegetation trajectories
    // ========================================================================
    // Natural cubic splines provide a continuous daily trajectory through
    // discrete vegetation snapshots (e.g. seasonal PAI or canopy height),
    // fitted independently for each grid cell with zero curvature at the
    // endpoints.

    List splineFitCpp(NumericVector knotX, NumericMatrix knotY)
    {
        int nCells = knotY.nrow();
        int n = knotX.size();
        NumericMatrix y2(nCells, n);
        std::vector<double> u(n);
        for (int c = 0; c < nCells; ++c) {
            // Natural boundary condition: zero second derivative at the
            // first knot.
            y2(c, 0) = 0.0;
            u[0] = 0.0;
            for (int i = 1; i < n - 1; ++i) {
                double sig = (knotX[i] - knotX[i - 1]) / (knotX[i + 1] - knotX[i - 1]);
                double p = sig * y2(c, i - 1) + 2.0;
                y2(c, i) = (sig - 1.0) / p;
                double uu = (knotY(c, i + 1) - knotY(c, i)) / (knotX[i + 1] - knotX[i])
                    - (knotY(c, i) - knotY(c, i - 1)) / (knotX[i] - knotX[i - 1]);
                u[i] = (6.0 * uu / (knotX[i + 1] - knotX[i - 1]) - sig * u[i - 1]) / p;
            }
            // Natural boundary condition at the last knot too.
            y2(c, n - 1) = 0.0;
            for (int k = n - 2; k >= 0; --k) {
                y2(c, k) = y2(c, k) * y2(c, k + 1) + u[k];
            }
        }
        return List::create(Named("knotX") = knotX, Named("knotY") = knotY, Named("y2") = y2);
    }

    NumericVector splineEvalCpp(List fit, double queryX)
    {
        NumericVector knotX = fit["knotX"];
        NumericMatrix knotY = fit["knotY"];
        NumericMatrix y2 = fit["y2"];
        int n = knotX.size();
        int nCells = knotY.nrow();

        // The knot positions are shared across every cell, so the segment
        // the query falls in is located once and reused for every row.
        double qx = queryX;
        if (qx < knotX[0]) qx = knotX[0];
        if (qx > knotX[n - 1]) qx = knotX[n - 1];
        int klo = 0, khi = n - 1;
        while (khi - klo > 1) {
            int k = (khi + klo) >> 1;
            if (knotX[k] > qx) khi = k; else klo = k;
        }
        double h = knotX[khi] - knotX[klo];
        double a = (knotX[khi] - qx) / h;
        double b = (qx - knotX[klo]) / h;

        NumericVector out(nCells);
        for (int c = 0; c < nCells; ++c) {
            out[c] = a * knotY(c, klo) + b * knotY(c, khi)
                + ((a * a * a - a) * y2(c, klo) + (b * b * b - b) * y2(c, khi)) * (h * h) / 6.0;
        }
        return out;
    }

} // namespace gridmodel

// [[Rcpp::export]]
Rcpp::List runmicro1Cpp(Rcpp::DataFrame obstime, Rcpp::DataFrame climdata, Rcpp::List vegp,
    Rcpp::List soilc, Rcpp::NumericMatrix lats, Rcpp::NumericMatrix lons, Rcpp::NumericVector zref, double z,
    Rcpp::NumericVector ufRef, Rcpp::NumericVector HRef, Rcpp::NumericVector dRef,
    Rcpp::NumericVector zmRef, Rcpp::NumericVector shelterc)
{
    return gridmodel::runmicro1Cpp(obstime, climdata, vegp, soilc, lats, lons,
        zref, z, ufRef, HRef, dRef, zmRef, shelterc);
}

// [[Rcpp::export]]
Rcpp::NumericVector dampingDepthGridCpp(Rcpp::NumericMatrix Vq, Rcpp::NumericMatrix Vm,
    Rcpp::NumericMatrix Vo, Rcpp::NumericMatrix Mc, Rcpp::NumericMatrix thetaS,
    Rcpp::NumericMatrix psie, Rcpp::NumericMatrix b,
    Rcpp::NumericVector soilm, Rcpp::NumericVector Ts, Rcpp::NumericVector pk)
{
    return gridmodel::dampingDepthGridCpp(Vq, Vm, Vo, Mc, thetaS, psie, b, soilm, Ts, pk);
}

// [[Rcpp::export]]
Rcpp::NumericVector soilWaterPotentialGridCpp(Rcpp::NumericMatrix thetaS, Rcpp::NumericMatrix psie,
    Rcpp::NumericMatrix b, Rcpp::NumericVector soilm)
{
    return gridmodel::soilWaterPotentialGridCpp(thetaS, psie, b, soilm);
}

// [[Rcpp::export]]
Rcpp::List surfaceResistGridCpp(Rcpp::List rad, Rcpp::List clim, Rcpp::NumericVector soilm,
    Rcpp::List soil, Rcpp::List veg, Rcpp::List ref)
{
    return gridmodel::surfaceResistGridCpp(rad, clim, soilm, soil, veg, ref);
}

// [[Rcpp::export]]
Rcpp::NumericVector referenceGroundResistCpp(Rcpp::NumericVector uref, Rcpp::NumericVector ufRef,
    Rcpp::NumericVector HRef, Rcpp::NumericVector Ta, Rcpp::NumericVector pk,
    double zref, double hRef, double paiRef)
{
    return gridmodel::referenceGroundResistCpp(uref, ufRef, HRef, Ta, pk, zref, hRef, paiRef);
}

// [[Rcpp::export]]
Rcpp::NumericVector soilmDistributeCpp(Rcpp::NumericMatrix Smin, Rcpp::NumericMatrix Smax,
    Rcpp::NumericMatrix tadd, Rcpp::NumericMatrix wet, Rcpp::NumericVector theta0, int tsteps)
{
    return gridmodel::soilmDistributeCpp(Smin, Smax, tadd, wet, theta0, tsteps);
}

// [[Rcpp::export]]
Rcpp::List groundHeatFluxCpp(
    Rcpp::NumericVector RabsGround, Rcpp::NumericVector rGz, Rcpp::NumericVector rGm,
    Rcpp::NumericVector RabsCanopy, Rcpp::NumericVector rHa, Rcpp::NumericVector rSurf_ref,
    Rcpp::NumericVector hSurf_ref, Rcpp::NumericVector hFol_ref,
    Rcpp::NumericVector soilm,
    Rcpp::NumericVector Ta, Rcpp::NumericVector rh, Rcpp::NumericVector pk,
    Rcpp::NumericMatrix Vq, Rcpp::NumericMatrix Vm, Rcpp::NumericMatrix Vo,
    Rcpp::NumericMatrix Mc, Rcpp::NumericMatrix thetaS, Rcpp::NumericMatrix psie,
    Rcpp::NumericMatrix b, Rcpp::NumericMatrix svfa,
    Rcpp::NumericVector emGround, Rcpp::NumericVector emCanopy,
    Rcpp::NumericMatrix slope, Rcpp::NumericMatrix aspect,
    Rcpp::NumericVector Rsw, Rcpp::NumericVector Rdif,
    Rcpp::NumericMatrix lats, Rcpp::NumericMatrix lons,
    Rcpp::IntegerVector year, Rcpp::IntegerVector month, Rcpp::IntegerVector day,
    Rcpp::NumericVector hour,
    Rcpp::NumericVector G_ref, Rcpp::NumericVector RabsGround_ref, Rcpp::NumericVector rGz_ref,
    Rcpp::NumericVector theta0_ref, Rcpp::NumericVector emGround_ref,
    Rcpp::NumericVector VqRef, Rcpp::NumericVector VmRef,
    Rcpp::NumericVector VoRef, Rcpp::NumericVector McRef,
    Rcpp::NumericVector thetaSRef, Rcpp::NumericVector psieRef, Rcpp::NumericVector bRef)
{
    return gridmodel::groundHeatFluxCpp(RabsGround, rGz, rGm, RabsCanopy, rHa, rSurf_ref, hSurf_ref, hFol_ref, soilm, Ta, rh, pk,
        Vq, Vm, Vo, Mc, thetaS, psie, b, svfa, emGround, emCanopy, slope, aspect, Rsw, Rdif, lats, lons,
        year, month, day, hour, G_ref, RabsGround_ref, rGz_ref, theta0_ref, emGround_ref,
        VqRef, VmRef, VoRef, McRef, thetaSRef, psieRef, bRef);
}

// [[Rcpp::export]]
Rcpp::NumericVector belowGroundShortcutCpp(
    Rcpp::NumericVector Ts_est, Rcpp::NumericVector Tgp_ref, Rcpp::NumericVector Tbp_ref,
    Rcpp::NumericVector DD, Rcpp::NumericVector mat, double reqhgt, double hiy = 8760.0)
{
    return gridmodel::belowGroundShortcutCpp(Ts_est, Tgp_ref, Tbp_ref, DD, mat, reqhgt, hiy);
}

// [[Rcpp::export]]
Rcpp::NumericVector canopyTempCpp(
    Rcpp::NumericVector RabsCanopy, Rcpp::NumericVector rHa, Rcpp::NumericVector G_est,
    Rcpp::NumericVector Ta, Rcpp::NumericVector rh, Rcpp::NumericVector pk,
    Rcpp::NumericMatrix svfa, Rcpp::NumericVector emCanopy,
    Rcpp::NumericVector rSurf_ref, Rcpp::NumericVector hSurf_ref)
{
    return gridmodel::canopyTempCpp(RabsCanopy, rHa, G_est, Ta, rh, pk, svfa, emCanopy, rSurf_ref, hSurf_ref);
}

// [[Rcpp::export]]
Rcpp::NumericVector leafTempCpp(
    Rcpp::NumericVector radLpar, Rcpp::NumericVector radLsw, Rcpp::NumericVector radLlw,
    Rcpp::NumericVector uz, Rcpp::NumericVector Ta, Rcpp::NumericVector rh,
    Rcpp::NumericVector pk, Rcpp::NumericVector Ca, Rcpp::NumericVector psi_r,
    Rcpp::NumericMatrix hgt, Rcpp::NumericMatrix pai,
    Rcpp::NumericMatrix Vcmax25, Rcpp::NumericMatrix Tup, Rcpp::NumericMatrix Tlw,
    Rcpp::NumericMatrix Dcrit, Rcpp::NumericMatrix alpha, Rcpp::NumericMatrix f0,
    Rcpp::NumericMatrix fd, Rcpp::NumericMatrix psi50, Rcpp::NumericMatrix apsi,
    Rcpp::NumericMatrix rpmin, Rcpp::NumericMatrix leafd, Rcpp::LogicalMatrix isC3,
    double reqhgt)
{
    return gridmodel::leafTempCpp(radLpar, radLsw, radLlw, uz, Ta, rh, pk, Ca, psi_r,
        hgt, pai, Vcmax25, Tup, Tlw, Dcrit, alpha, f0, fd, psi50, apsi, rpmin, leafd, isC3, reqhgt);
}

// [[Rcpp::export]]
Rcpp::List longwaveGridCpp(
    Rcpp::NumericVector Ts_est, Rcpp::NumericVector Tcanopy_est, Rcpp::NumericVector Ta,
    Rcpp::NumericMatrix hgt, Rcpp::NumericMatrix pai, Rcpp::NumericMatrix paia,
    Rcpp::NumericMatrix svfa, Rcpp::NumericVector lwdown, double reqhgt)
{
    return gridmodel::longwaveGridCpp(Ts_est, Tcanopy_est, Ta, hgt, pai, paia, svfa, lwdown, reqhgt);
}

// [[Rcpp::export]]
Rcpp::List aboveCanopyProfileGridCpp(
    Rcpp::NumericVector Tcanopy_est, Rcpp::NumericVector rHa, Rcpp::NumericVector L, Rcpp::NumericVector uf,
    Rcpp::NumericMatrix d, Rcpp::NumericMatrix zh,
    Rcpp::NumericVector Ta, Rcpp::NumericVector rh, Rcpp::NumericVector rSurf_ref,
    Rcpp::NumericVector hSurf_ref,
    Rcpp::NumericVector rGz, Rcpp::NumericVector rGreq, Rcpp::NumericMatrix zs,
    double reqhgt, double zref)
{
    return gridmodel::aboveCanopyProfileGridCpp(Tcanopy_est, rHa, L, uf, d, zh, Ta, rh, rSurf_ref, hSurf_ref,
        rGz, rGreq, zs, reqhgt, zref);
}

// [[Rcpp::export]]
Rcpp::List belowCanopyProfileGridCpp(
    Rcpp::NumericVector Tcanopy_est, Rcpp::NumericVector Ts_est, Rcpp::NumericVector rBL,
    Rcpp::NumericVector rGh, Rcpp::NumericVector rGm, Rcpp::NumericVector rGz,
    Rcpp::NumericMatrix shapeR, Rcpp::NumericMatrix shapeC,
    Rcpp::NumericMatrix hgtCell, Rcpp::NumericMatrix paiCell,
    Rcpp::NumericVector Ta, Rcpp::NumericVector rh, Rcpp::NumericVector pk, Rcpp::NumericVector rSurf_ref,
    Rcpp::NumericVector hSurf_ref,
    double reqhgt, double zref)
{
    return gridmodel::belowCanopyProfileGridCpp(Tcanopy_est, Ts_est, rBL, rGh, rGm, rGz,
        shapeR, shapeC, hgtCell, paiCell, Ta, rh, pk, rSurf_ref, hSurf_ref, reqhgt, zref);
}

// [[Rcpp::export]]
Rcpp::List belowCanopyProfilePointCpp(
    Rcpp::NumericVector Tcanopy, Rcpp::NumericVector Tground, Rcpp::NumericVector groundhr,
    Rcpp::NumericVector rBL, Rcpp::NumericVector L, Rcpp::NumericVector uf,
    double d, double zm, double hgt, double pai,
    Rcpp::NumericVector Ta, Rcpp::NumericVector rh, Rcpp::NumericVector pk, Rcpp::NumericVector rSurf,
    Rcpp::NumericVector hSurf,
    double reqhgt, double zref)
{
    return gridmodel::belowCanopyProfilePointCpp(Tcanopy, Tground, groundhr, rBL, L, uf,
        d, zm, hgt, pai, Ta, rh, pk, rSurf, hSurf, reqhgt, zref);
}

// [[Rcpp::export]]
Rcpp::List computeSolarSeriesRCpp(double lat, double lon, Rcpp::IntegerVector year,
    Rcpp::IntegerVector month, Rcpp::IntegerVector day, Rcpp::NumericVector hour)
{
    return gridmodel::computeSolarSeriesRCpp(lat, lon, year, month, day, hour);
}

// [[Rcpp::export]]
Rcpp::List radiationPointCpp(double pai, double paia, double x, double lref, double ltra,
    double clump, double gref, double svfa, double slope, double aspect,
    Rcpp::NumericVector zend, Rcpp::NumericVector azid, Rcpp::NumericVector shadow,
    Rcpp::NumericVector Rsw, Rcpp::NumericVector Rdif)
{
    return gridmodel::radiationPointCpp(pai, paia, x, lref, ltra, clump, gref, svfa,
        slope, aspect, zend, azid, shadow, Rsw, Rdif);
}

// [[Rcpp::export]]
Rcpp::List splineFitCpp(Rcpp::NumericVector knotX, Rcpp::NumericMatrix knotY)
{
    return gridmodel::splineFitCpp(knotX, knotY);
}

// [[Rcpp::export]]
Rcpp::NumericVector splineEvalCpp(Rcpp::List fit, double queryX)
{
    return gridmodel::splineEvalCpp(fit, queryX);
}
