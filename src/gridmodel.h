// gridmodel.h
// Landscape-scale microclimate model.
//
// A fully coupled point simulation supplies the temporally resolved reference
// surface state. The grid model then represents how fine-scale terrain, canopy
// structure and soil properties alter radiation, turbulent exchange, surface
// temperature, soil moisture and vertical microclimate from cell to cell. This
// keeps the expensive nonlinear energy/water balance at the reference scale
// while retaining the physical mechanisms that generate spatial heterogeneity.
//
// Two spatial effects require explicit treatment beyond the bulk point model:
// slope/aspect and terrain horizons alter the local radiation field, and a
// requested height can lie inside a real canopy, requiring vertical foliage and
// scalar-transport structure rather than a single bulk exchange surface.
#pragma once
#include <Rcpp.h>
#include "pointmodel.h"

namespace gridmodel {

    // Optical state of one cell's canopy, fixed while solar geometry changes.
    // Whole-canopy terms describe absorption/reflection to the ground and sky;
    // the `_z` terms describe the part of the canopy above a requested height.
    struct tirstruct {
        bool vegetated;            // whether the cell has any plant area (pai > 0)
        double pait;                // clumping-adjusted total plant area index
        utils::tsdifstruct tsd;     // diffuse two-stream parameters, whole canopy
        double amx;                  // maximum permissible albedo
        double trdn;                 // squared ground-level gap fraction (clump^2)
        double albd;                 // diffuse albedo, canopy top
        double Rddn_g;                // diffuse transmission fraction to the ground
        double gi;                   // downward gap fraction, canopy top -> requested height
        double trdu;                  // squared gap fraction, ground -> requested height
        double paiaa;                 // clumping-adjusted plant area index above requested height
        double Rddn_z;                // diffuse transmission fraction down to requested height
        double Rdup_z;                 // diffuse fraction reflected back up past requested height
    };

    // Radiation state for one cell and timestep: absorbed energy used by
    // surface energy balances plus the direct/diffuse/reflected streams and
    // leaf-scale absorption at the requested height.
    struct radmodel2 {
        double radGsw = 0.0;   // shortwave absorbed by the ground alone (W/m^2)
        double radGlw = 0.0;   // longwave absorbed by the ground alone (W/m^2)
        double radCsw = 0.0;   // shortwave absorbed, canopy + ground combined (W/m^2)
        double radClw = 0.0;   // longwave absorbed, canopy + ground combined (W/m^2)
        double Rbdown = 0.0;   // direct-beam radiation reaching the requested height (W/m^2)
        double Rddown = 0.0;   // diffuse (+ direct-beam-scattered) radiation reaching the requested height (W/m^2)
        double Rdup = 0.0;     // radiation reflected back upward past the requested height (W/m^2)
        double radLsw = 0.0;   // shortwave absorbed by a unit leaf at the requested height (W/m^2)
        double radLpar = 0.0;  // PAR absorbed by a unit leaf at the requested height (W/m^2)
        double lwout = 0.0;    // longwave emitted upward from the ground surface (W/m^2)
        double emGround = 0.0; // coefficient on the ground's own emitted longwave in its
                               // energy balance: not a bare emissivity, but net of the share
                               // canopy and terrain return (see GridRadlwabsStepCpp)
        double emCanopy = 0.0; // the same, for the canopy+ground surface

        double zend = 0.0;     // solar zenith angle for this timestep (degrees)
    };

    tirstruct GridRadswabsSetupCpp(double pai, double paia, double x, double lref, double ltra,
        double clump, double gref);

    radmodel2 GridRadswabsStepCpp(const tirstruct& tir, double pai, double clump, double gref,
        double svfa, double si, bool shadow, double zenr, double x, double Rsw, double Rdif);

    void GridRadlwabsStepCpp(const tirstruct& tir, double svfa, double tc, double Rlw, radmodel2& out);

    // Evaluate shortwave, PAR, longwave radiation and wind across the full
    // grid and time series. Forcing may be shared by the whole landscape or
    // supplied per coarse reference cell. Local latitude/longitude, terrain,
    // canopy structure and shelter then determine each fine cell's solar
    // exposure and aerodynamic state at the requested height `z`.
    Rcpp::List runmicro1Cpp(Rcpp::DataFrame obstime, Rcpp::DataFrame climdata, Rcpp::List vegp,
        Rcpp::List soilc, Rcpp::NumericMatrix lats, Rcpp::NumericMatrix lons, Rcpp::NumericVector zref, double z,
        Rcpp::NumericVector ufRef, Rcpp::NumericVector HRef, Rcpp::NumericVector dRef,
        Rcpp::NumericVector zmRef, Rcpp::NumericVector shelterc);

    // ========================================================================
    // Wind component.
    // ========================================================================

    // Aerodynamic state of a grid cell. `uf` and L describe the turbulent
    // surface layer; `uz` is wind at the requested height; `rHa` connects the
    // canopy exchange surface to the reference atmosphere; `rGz` additionally
    // includes transport through the canopy from the ground; `rGm` is the mean
    // resistance from the ground to canopy heights, across which the ground
    // exchanges with canopy-warmed air; `a2` controls within-canopy scalar
    // mixing.
    struct windresult2 {
        double uz = 0.0;        // wind at the requested height on the exchange scale (m/s)
        double uzActual = 0.0;  // the same wind as actual wind, for reporting (m/s)
        // Friction velocity (m/s).
        double uf = 0.0;
        // Soil surface to reference-height scalar resistance (s/m).
        double rGz = 0.0;
        // Canopy exchange surface to reference-height resistance (s/m).
        double rHa = 0.0;
        // Mean soil-to-canopy-source resistance (s/m).
        double rGm = 0.0;
        double rGh = 0.0;   // ground to canopy top, for the below-canopy profile (s/m)
        double rGreq = 0.0; // ground to the requested height, above canopy top only (s/m)
        // Within-canopy scalar-mixing parameter.
        double a2 = 0.0;
        // Monin-Obukhov length (m).
        double L = 0.0;
    };

    // Scale the converged reference aerodynamic state to one grid cell.
    // Canopy displacement/roughness and topographic shelter modify friction
    // velocity and the wind profile; the implied Monin-Obukhov length is then
    // recovered from that scaled state, in closed form: `uPeak` is the cell's
    // stable-branch peak from utils::stableRecoveryPeakCpp, fixed by its
    // geometry. `wstar` is the free-convective velocity of the reference
    // run's sensible heat; both driving winds include it as the point model
    // does, so the scaling leaves it unreduced by shelter. Above-canopy wind follows similarity theory, while
    // within-canopy wind decays exponentially with depth.
    windresult2 windCpp(double z, double zref, double hgt, double pai, double d, double zm,
        double uref, double shelterc, double ufRef, double dRef, double zmRef,
        const utils::ColumnStruct& col, bool hasColumn, double uPeak, double wstar);

    // The reference run's ground-to-reference resistance as this grid model
    // evaluates it for the reference's own geometry, unsheltered, from the
    // reference driving wind, friction velocity and sensible heat. Used so
    // the reference's side of the ground heat flux scaling is computed as a
    // cell's is.
    Rcpp::NumericVector referenceGroundResistCpp(Rcpp::NumericVector uref, Rcpp::NumericVector ufRef,
        Rcpp::NumericVector HRef, Rcpp::NumericVector Ta, Rcpp::NumericVector pk,
        double zref, double hRef, double paiRef);

    // ========================================================================
    // Soil-moisture spatial distribution.
    // ========================================================================

    // Distribute reference surface soil moisture according to local soil
    // capacity and topographic wetness. The reference moisture is first
    // expressed as a fraction between each cell's residual and saturated water
    // contents, shifted by the wetness anomaly in logit space, and converted
    // back to volumetric water content. Each cell is then raised toward
    // saturation by its wet weight, which is 1 on permanently wet ground.
    Rcpp::NumericVector soilmDistributeCpp(Rcpp::NumericMatrix Smin, Rcpp::NumericMatrix Smax,
        Rcpp::NumericMatrix tadd, Rcpp::NumericMatrix wet, Rcpp::NumericVector theta0, int tsteps);

    // ========================================================================
    // Ground heat flux / surface temperature shortcut
    // ========================================================================
    // Transfer the reference simulation's resolved diurnal ground-heat-flux
    // cycle to each cell rather than solving soil heat diffusion independently.
    // The amplitude is scaled by the target/reference ratio of the daily flux
    // swing each surface sustains: its zero-flux surface temperature swing,
    // its linearised surface exchange coefficient and its soil's thermal
    // admittance, for a periodic flux into a semi-infinite conducting medium
    // (de Vries 1963; Campbell 1985). The reference's side is evaluated as
    // the cell's is, from its absorbed ground radiation (`RabsGround_ref`),
    // its resistance as this grid evaluates it (`rGz_ref`, from
    // referenceGroundResistCpp), its surface moisture and its soil, so a
    // cell identical to the reference scales by one. Only the reference
    // flux's 24-hour harmonic is scaled and moved to the cell's timing; the
    // rest is kept at clock time. A reference day with too small a predicted
    // swing leaves the reference flux unscaled.
    // The phase is adjusted so slope/aspect-driven illumination changes occur
    // at the target cell's own sunrise/sunset, since a flat reference surface
    // peaks and crosses zero at a different clock time than a sloped,
    // differently-oriented cell's true peak would. A local Penman-Monteith
    // balance first estimates the zero-flux surface temperature and then
    // resolves the final surface temperature using the transferred heat flux.
    // Sky-view factor modifies radiative exposure and local soil moisture
    // modifies evaporative cooling. Where the cell has foliage, the final pass
    // solves the bulk canopy+ground surface and the ground together: their two
    // energy balances, linearised about air temperature, are two linear
    // equations in the two temperatures, with the soil's vapour at the ground's
    // own temperature (the vapour network, with the foliage's share `hFol_ref`
    // from surfaceResistGridCpp), the ground's air raised by the bulk surface's
    // sensible flux across the node's resistance, and its vapour boundary the
    // network seen from the soil. The node is at whichever of the mean source
    // height and the exchange surface is nearer the reference. The pair is
    // returned as `Ts_est` and `Tcanopy_est` (NA where not solved together).
    // Without foliage the ground alone is solved, against reference air. The
    // zero-flux pass keeps reference air across rGz, as the reference run's own
    // zero-flux temperature does.
    //
    // This is an approximation to the full soil heat solve. It is bounded by
    // empirically calibrated ground-air temperature limits and a dew-point
    // margin to prevent pathological values under weakly constrained calm/night
    // conditions. Complete 24-hour days are required because the transfer is
    // based on daily harmonics. A known, accepted limitation: the phase
    // correction relies on a direct-beam timing signal, so it is least
    // effective for a cell that receives little direct beam for much of the
    // day (e.g. a steep, persistently shaded slope near the summer solstice).
    Rcpp::List groundHeatFluxCpp(
        Rcpp::NumericVector RabsGround, Rcpp::NumericVector rGz, Rcpp::NumericVector rGm,
        Rcpp::NumericVector RabsCanopy, Rcpp::NumericVector rHa, Rcpp::NumericVector rSurf_ref, Rcpp::NumericVector hSurf_ref,
        Rcpp::NumericVector hFol_ref, Rcpp::NumericVector soilm,
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
        Rcpp::NumericVector thetaSRef, Rcpp::NumericVector psieRef, Rcpp::NumericVector bRef);

    // Estimate temperature at a requested soil depth by treating the soil as
    // a depth-dependent low-pass filter. Near the surface the estimate follows
    // the cell's own surface temperature. At intermediate depths it adopts the
    // reference model's resolved subsurface anomaly, rescaled to the target
    // cell's daily surface range and mean. At still greater depths the diurnal
    // signal is suppressed and the estimate tends toward mean annual
    // temperature according to the annual damping scale (diurnal damping depth
    // multiplied by sqrt(365)). `DD` is each cell's diurnal damping depth for
    // each hour (constant within a day), from dampingDepthGridCpp. Accuracy is necessarily weaker where target and
    // reference surface regimes differ strongly because no independent local
    // heat-diffusion equation is solved -- most noticeable in the 30cm-1m
    // range, of diminishing practical consequence by 2m, where the true
    // signal itself has largely flattened.
    Rcpp::NumericVector belowGroundShortcutCpp(
        Rcpp::NumericVector Ts_est, Rcpp::NumericVector Tgp_ref, Rcpp::NumericVector Tbp_ref,
        Rcpp::NumericVector DD, Rcpp::NumericVector mat, double reqhgt, double hiy = 8760.0);

    // Each cell's diurnal damping depth (m) for each day, from its own soil and
    // that day's mean surface soil moisture, surface temperature and pressure,
    // repeated for the day's hours.
    Rcpp::NumericVector dampingDepthGridCpp(Rcpp::NumericMatrix Vq, Rcpp::NumericMatrix Vm,
        Rcpp::NumericMatrix Vo, Rcpp::NumericMatrix Mc, Rcpp::NumericMatrix thetaS,
        Rcpp::NumericMatrix psie, Rcpp::NumericMatrix b,
        Rcpp::NumericVector soilm, Rcpp::NumericVector Ts, Rcpp::NumericVector pk);

    // Each cell's vapour network, for every hour, formed as the point model
    // forms it (pointmodel::canopyWaterBudgetCpp) but from the cell's own
    // stomata (sunlit and shaded foliage PAR from `rad`, its plant-type
    // physiology and soil water potential), its own column resistances and
    // soil humidity, and the reference run's wet share and film water per unit
    // plant area. Returns `rSurf` (the network's vapour resistance less rHa),
    // `hFol` (the foliage's share of the network's conductance) and `hSurf`
    // (a humidity factor with the soil at the surface temperature, h + (1 - h)
    // hr, used only for cells the coupled ground pass does not solve) and `hr`
    // (the soil surface's humidity, at air temperature). The
    // soil's vapour pressure, needed before any ground temperature exists for
    // the film's water limit, is taken at air temperature. Two fixed passes:
    // the first at air temperature with the film unlimited gives the canopy
    // temperature at which the second is evaluated.
    Rcpp::List surfaceResistGridCpp(Rcpp::List rad, Rcpp::List clim, Rcpp::NumericVector soilm,
        Rcpp::List soil, Rcpp::List veg, Rcpp::List ref);

    // Each cell's soil water potential (MPa) from its own surface soil moisture
    // and water-retention curve, used as the root-zone potential for leaves.
    Rcpp::NumericVector soilWaterPotentialGridCpp(Rcpp::NumericMatrix thetaS, Rcpp::NumericMatrix psie,
        Rcpp::NumericMatrix b, Rcpp::NumericVector soilm);

    // Bulk canopy+ground exchange-surface temperature from a local
    // Penman-Monteith balance. Absorbed radiation, aerodynamic resistance,
    // ground heat flux vary by cell, and sky exposure enters through the
    // cell's own longwave emission coefficient `emCanopy` (`svfa` itself is
    // used here only to mask non-land cells). The surface resistance
    // and humidity factor are the cell's own, from surfaceResistGridCpp. Two
    // temperature updates reduce the longwave linearisation error. The result
    // is bounded by the same empirical surface-air envelope and dew-point
    // floor as the ground temperature, since this one-shot balance cannot
    // know the cell's own heating either.
    Rcpp::NumericVector canopyTempCpp(
        Rcpp::NumericVector RabsCanopy, Rcpp::NumericVector rHa, Rcpp::NumericVector G_est,
        Rcpp::NumericVector Ta, Rcpp::NumericVector rh, Rcpp::NumericVector pk,
        Rcpp::NumericMatrix svfa, Rcpp::NumericVector emCanopy,
        Rcpp::NumericVector rSurf_ref, Rcpp::NumericVector hSurf_ref);

    // Temperature of an individual leaf at a requested height within the
    // canopy. This is distinct from the bulk canopy exchange temperature: the
    // leaf uses radiation at its own height, leaf-scale boundary-layer
    // resistance and stomatal conductance from the photosynthesis/hydraulic
    // model. PAR controls stomatal physiology while total absorbed radiation
    // drives the leaf energy balance. Values are defined only for
    // 0 < reqhgt < canopy height.
    Rcpp::NumericVector leafTempCpp(
        Rcpp::NumericVector radLpar, Rcpp::NumericVector radLsw, Rcpp::NumericVector radLlw,
        Rcpp::NumericVector uz, Rcpp::NumericVector Ta, Rcpp::NumericVector rh,
        Rcpp::NumericVector pk, Rcpp::NumericVector Ca, Rcpp::NumericVector psi_r,
        Rcpp::NumericMatrix hgt, Rcpp::NumericMatrix pai,
        Rcpp::NumericMatrix Vcmax25, Rcpp::NumericMatrix Tup, Rcpp::NumericMatrix Tlw,
        Rcpp::NumericMatrix Dcrit, Rcpp::NumericMatrix alpha, Rcpp::NumericMatrix f0,
        Rcpp::NumericMatrix fd, Rcpp::NumericMatrix psi50, Rcpp::NumericMatrix apsi,
        Rcpp::NumericMatrix rpmin, Rcpp::NumericMatrix leafd, Rcpp::LogicalMatrix isC3,
        double reqhgt);

    // Downward and upward longwave radiation at the requested height. Above
    // canopy this is atmospheric longwave from the visible sky versus surface
    // emission. Within canopy, plant area above and below the point determines
    // the mixture of sky, canopy and ground emission seen in each direction.
    // Sky-view factor reduces exchange with the unobstructed atmosphere.
    Rcpp::List longwaveGridCpp(
        Rcpp::NumericVector Ts_est, Rcpp::NumericVector Tcanopy_est, Rcpp::NumericVector Ta,
        Rcpp::NumericMatrix hgt, Rcpp::NumericMatrix pai, Rcpp::NumericMatrix paia,
        Rcpp::NumericMatrix svfa, Rcpp::NumericVector lwdown, double reqhgt);

    // Air temperature and humidity above the canopy, from the transport
    // column's resistances where there is a canopy (roughness sublayer
    // included) and from Monin-Obukhov similarity over bare ground. The
    // estimated exchange-surface temperature and
    // vapour pressure form the lower boundary; reqhgt samples the fraction of
    // the full surface-to-reference aerodynamic resistance accumulated by that
    // height, preserving continuity at canopy top and at the reference level.
    // This function itself has no notion of whether a given cell's requested
    // height genuinely sits above its own canopy -- it simply evaluates the
    // above-canopy profile everywhere; it is the caller's responsibility to
    // use the result only where reqhgt is at or above that cell's own canopy
    // height, and use the below-canopy profile elsewhere.
    Rcpp::List aboveCanopyProfileGridCpp(
        Rcpp::NumericVector Tcanopy_est, Rcpp::NumericVector rHa, Rcpp::NumericVector L, Rcpp::NumericVector uf,
        Rcpp::NumericMatrix d, Rcpp::NumericMatrix zh,
        Rcpp::NumericVector Ta, Rcpp::NumericVector rh, Rcpp::NumericVector rSurf_ref, Rcpp::NumericVector hSurf_ref,
        Rcpp::NumericVector rGz, Rcpp::NumericVector rGreq, Rcpp::NumericMatrix zs,
        double reqhgt, double zref);

    // Air temperature and relative humidity inside the canopy, from the
    // far-field closure: a smooth gradient constrained by the ground and
    // canopy-top boundary conditions and by the integrated canopy resistance
    // between them. The profile's vertical shape is taken from the column at
    // neutral; its resistances carry each hour's stability.
    Rcpp::List belowCanopyProfileGridCpp(
        Rcpp::NumericVector Tcanopy_est, Rcpp::NumericVector Ts_est, Rcpp::NumericVector rBL,
        Rcpp::NumericVector rGh, Rcpp::NumericVector rGm, Rcpp::NumericVector rGz,
        Rcpp::NumericMatrix shapeR, Rcpp::NumericMatrix shapeC,
        Rcpp::NumericMatrix hgtCell, Rcpp::NumericMatrix paiCell,
        Rcpp::NumericVector Ta, Rcpp::NumericVector rh, Rcpp::NumericVector pk, Rcpp::NumericVector rSurf_ref, Rcpp::NumericVector hSurf_ref,
        double reqhgt, double zref);

    // Single-location version of the same below-canopy far-field
    // closure, driven by the fully converged point-model state. Unlike the grid
    // shortcut it can use the point model's actual (generally unsaturated)
    // ground-surface humidity as the lower moisture boundary. As in the
    // gridded version, the vertical shape is taken from the column at
    // neutral, and the resistances from the point model's own friction
    // velocity and Obukhov length each hour. Returns NA when reqhgt is not
    // strictly inside a vegetated canopy.
    Rcpp::List belowCanopyProfilePointCpp(
        Rcpp::NumericVector Tcanopy, Rcpp::NumericVector Tground, Rcpp::NumericVector groundhr,
        Rcpp::NumericVector rBL, Rcpp::NumericVector L, Rcpp::NumericVector uf,
        double d, double zm, double hgt, double pai,
        Rcpp::NumericVector Ta, Rcpp::NumericVector rh, Rcpp::NumericVector pk, Rcpp::NumericVector rSurf, Rcpp::NumericVector hSurf,
        double reqhgt, double zref);

    // ========================================================================
    // Seasonal interpolation of vegetation properties
    // ========================================================================
    // Fit independent natural cubic splines through each cell's supplied
    // vegetation snapshots, allowing seasonal properties to vary smoothly
    // rather than jumping between layers.
    Rcpp::List splineFitCpp(Rcpp::NumericVector knotX, Rcpp::NumericMatrix knotY);

    // Evaluate every cell's fitted spline at one position; values outside
    // the supplied seasonal range are held at the nearest endpoint.
    Rcpp::NumericVector splineEvalCpp(Rcpp::List fit, double queryX);

} // namespace gridmodel
