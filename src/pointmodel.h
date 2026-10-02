// pointmodel.h
// Point-scale microclimate model.
//
// The model solves the coupled surface energy and water balances of a
// vegetation--soil patch through time. Radiation interception sets the energy
// available to canopy and ground; aerodynamic exchange links those surfaces to
// the atmosphere; photosynthesis, stomatal regulation and plant hydraulics set
// transpiration; canopy interception partitions rainfall; and multilayer soil
// heat and water transport update the subsurface state. These processes feed
// back on one another through surface temperature, humidity, soil moisture and
// atmospheric stability until each timestep converges.
//
// vegpstruct contains the patch's vegetation properties, whereas envstruct
// carries the meteorological forcing and the coupled surface state for the
// current timestep. Generic physical relationships shared with the grid model
// (e.g. two-stream radiation and Monin-Obukhov utilities) are declared in
// utils.h.
#pragma once
#include <vector>
#include <limits>
#include "utils.h"
#include "constants.h"

namespace pointmodel {

    // ========================================================================
    // Vegetation and environment: parameters and state shared across sections
    // ========================================================================

    // Vegetation parameters for a patch. Most fields are constant across an
    // entire model run; Lfrac is the exception, and may be updated between
    // timesteps (e.g. to track seasonal leaf-on/leaf-off state).
    struct vegpstruct {
        // Canopy structure
        double hgt;      // canopy height (m)
        double pai;      // total one-sided plant area index (m^2/m^2), whole canopy
        double x;        // Campbell foliage angle coefficient
        double clump;    // canopy clumping factor
        double Lfrac;    // fraction of pai that is living leaf material, whole canopy
        double len;      // leaf length (m)
        double wid;      // leaf width (m)

        // Radiative properties
        double lref;     // leaf reflectance, shortwave (0-1)
        double ltra;     // leaf transmittance, shortwave (0-1)
        double lrefp;    // leaf reflectance, PAR (0-1)
        double ltrap;    // leaf transmittance, PAR (0-1)
        // Longwave emissivity is represented by the shared surfaceEmissivity
        // constant rather than varying among vegetation types.
        double mwft;     // maximum water film thickness on the canopy (mm)

        // Photosynthesis (JULES-SOX)
        double Vcmax25;  // maximum carboxylation rate at 25 degC (micromol / m^2 / s)
        double Tup;      // upper temperature limit for photosynthesis (deg C)
        double Tlw;      // lower temperature limit for photosynthesis (deg C)
        double Dcrit;    // vapour pressure deficit at stomatal closure (kPa)
        double alpha;    // photosynthetic quantum efficiency (mol CO2 / mol PAR)
        double f0;       // ratio of ci / (ci - Gamma) under low vapour pressure deficit
        double fd;       // fraction of Vcmax lost to leaf dark respiration
        double gsmaxCap = -1.0; // optional empirical ceiling on stomatal conductance
                          // (mol m^-2 s^-1), applied in addition to the
                          // photosynthesis/hydraulic model's own limit; a
                          // negative value means no additional cap.

        // Xylem hydraulics
        double rpmin;    // minimum whole-plant hydraulic resistance (m^2 s MPa / mol),
                          // precomputed per plant-functional-type from stem/vessel
                          // geometry via the xylem-tapering model -- a constant for a
                          // given plant-functional-type, since it depends only on
                          // static structural parameters, not on run-time conditions
        double psi50;    // water potential at 50% loss of hydraulic conductance (MPa)
        double apsi = -1.0; // hydraulic vulnerability curve shape parameter; computed
                             // from psi50 on first use if left negative.

        // Root water uptake
        double pTAW;     // fraction of total available soil water depletable before stress
        double root50;   // depth above which 50% of roots lie (m)
        double root95;   // depth above which 95% of roots lie (m)
    };

    // Atmospheric forcing and coupled surface state for one timestep. In
    // addition to the external weather drivers, this structure carries the
    // radiation, turbulence, canopy temperature and soil-water demands that
    // link the model's component energy and water balances during convergence.
    struct envstruct {
        // Meteorological drivers
        double tair;    // air temperature (deg C)
        double rh;      // relative humidity (%)
        double pk;      // atmospheric pressure (kPa)
        double Ca;      // atmospheric CO2 concentration (ppm)
        double precip;  // precipitation (mm)
        double psi_r;   // mean water potential in the root zone (MPa)
        double Rsw;     // shortwave radiation flux density on the horizontal (W/m^2)
        double Rdif;    // diffuse radiation flux density on the horizontal (W/m^2)
        double Rlw;     // longwave radiation from the sky, on the horizontal (W/m^2)
        double uref;    // wind speed at the reference height (m/s)
        // Coefficients on each surface's own emitted longwave in its energy balance,
        // set by RadlwabsStepCpp with the absorbed longwave they belong to. Not bare
        // emissivities: they carry the share of that emission which terrain and canopy
        // return, so each pair of surfaces exchanges through one coefficient and two
        // surfaces at the same temperature exchange nothing. `emGround` is the soil
        // surface's, `emCanopy` the canopy+ground surface's; they coincide on bare ground.
        double emGround = mc::surfaceEmissivity;
        double emCanopy = mc::surfaceEmissivity;

        // Canopy surface state
        double tcanopy = 15.0; // combined canopy+ground "big leaf" surface temperature (deg C) --
                                // solved by penmanMonteithCpp against RabsCanopy (the combined
                                // canopy+ground absorption total). This is the whole surface's
                                // temperature (used as such throughout the energy balance and as
                                // the model's own Tcanopy output), not an individual leaf's -- it
                                // is also fed into the stomatal-conductance formulas, which do
                                // treat it as a single leaf's temperature, since the "big leaf"
                                // abstraction has no separate leaf-scale temperature of its own.
        double PARabs = 0.0;   // radiation absorbed at the leaf surface, used for photosynthesis (W/m^2)

        // Radiation outputs. Rground/Rcanopy/RgroundPAR/RcanopyPAR/albedo
        // are written by RadswabsStepCpp; RabsGround/RabsCanopy (shortwave +
        // longwave, the totals consumed by the surface energy balances)
        // are written by RadlwabsStepCpp, which should be called
        // immediately afterwards for the same timestep. For all of these,
        // "ground" is the ground surface only, while "canopy" is the
        // combined canopy + ground total (the convention a Big Leaf model
        // needs: a single combined-surface value driving the canopy-level
        // energy balance, and a ground-only value driving the soil model).
        double Rground = 0.0;     // shortwave absorbed by the ground alone (W/m^2)
        double Rcanopy = 0.0;     // shortwave absorbed, canopy + ground combined (W/m^2)
        double RgroundPAR = 0.0;  // PAR absorbed by the ground alone (W/m^2)
        double RcanopyPAR = 0.0;  // PAR absorbed, canopy + ground combined (W/m^2); canopy
                                   // (foliage) only PAR absorption, as needed for
                                   // photosynthesis, is RcanopyPAR - RgroundPAR
        double albedo = 0.0;      // effective albedo
        double RabsGround = 0.0;  // total (shortwave + longwave) radiation absorbed by the ground (W/m^2)
        double RabsCanopy = 0.0;  // total (shortwave + longwave) radiation absorbed, canopy + ground combined (W/m^2)

        // Wind / stability. H is an input: the current sensible heat flux
        // estimate from the outer energy-balance iteration. uf, LL, uh
        // and zm are outputs. psi_m and psi_h are both: read as the
        // previous iterate, then overwritten with the updated value.
        double H = 0.0;        // sensible heat flux (W/m^2)
        double uf = 0.0;       // friction velocity (m/s)
        double LL = 1e99;      // Monin-Obukhov length (m)
        double uh = 0.0;       // wind speed at the top of the canopy (m/s)
        double windScale = 1.0; // measured wind over the exchange profile's wind at zref: turns
                                // profile winds into actual winds (see windmodelCpp)
        double zm = 0.0;       // roughness length for momentum (m)
        double zh = 0.0;       // roughness length for heat (m); see utils::scalarRoughlengthCpp
        double psi_m = 0.0;    // diabatic correction for momentum
        double psi_h = 0.0;    // diabatic correction for heat
        // Low-wind stability damping (windmodelCpp), both reset at the start of
        // each timestep: whether an earlier pass of this timestep has set uf and
        // LL, and whether stability has since changed sign.
        bool stabPassDone = false;
        bool stabFlipped = false;

        // Soil energy/water balance drivers (read by SoilHeatCpp/SoilWaterCpp).
        // RabsGround itself is one of the radiation outputs above.
        double rHa = 0.0;   // aerodynamic resistance to heat transfer, ground to reference height (s/m); see groundaeroresistCpp
        double rVg = -1.0;  // the soil's vapour resistance where it differs from rHa (s/m); <= 0 means rHa
        // Conductances of the vapour network's foliage (a) and soil (b) legs
        // (m/s), set only on the copy that carries the whole surface's balance
        // in the soil's surface row; there the soil and foliage share one
        // temperature and the humidity factor is (a + b hr)/(a + b). < 0 = unused.
        double wholeA = -1.0;
        double wholeB = -1.0;
        double Et = 0.0;    // stomatal transpiration demand for this timestep (mm), driving root
                             // water uptake in SoilWaterCpp -- the stomatal pathway alone
                             // (aerodynamic + stomatal resistance in series), not the combined-
                             // surface evaporative flux, since only water actually drawn through
                             // the stomata depletes root-zone water.
        double precipGround = 0.0; // rainfall reaching the ground after canopy interception
                             // and throughfall (mm); equals precip over bare ground and supplies
                             // the surface input to the soil-water balance.
    };

    // ========================================================================
    // Radiation: two-stream shortwave absorption by canopy and ground
    // ========================================================================
    // Canopy structure and optical properties determine how incoming solar
    // radiation is reflected, transmitted and absorbed by foliage and ground.
    // Diffuse transfer depends only on canopy state, whereas the direct-beam
    // component changes with solar position and slope/aspect.

    struct solmodel {
        double zend;  // solar zenith angle (degrees)
        double zenr;  // solar zenith angle (radians)
        double azid;  // solar azimuth angle (degrees)
        double azir;  // solar azimuth angle (radians)
    };

    // Optical state of the canopy that is fixed while the sun moves. The
    // same canopy structure is used for total shortwave and PAR, but each
    // waveband has its own leaf/ground reflectance, transmittance and resulting
    // diffuse two-stream solution.
    struct radsetup {
        bool vegetated;    // whether the patch has any plant area (pai > 0)
        double pait;       // clumping-adjusted plant area index
        double trd;        // squared canopy gap fraction
        double amx;        // maximum permissible albedo for this patch
        double albd;       // diffuse albedo
        double groundRdd;  // fraction of diffuse radiation reaching the ground
        utils::tsdifstruct tsd;   // diffuse two-stream parameters
        double amxPAR;      // maximum permissible albedo for this patch, PAR
        double albdPAR;     // diffuse albedo, PAR
        double groundRddPAR; // fraction of diffuse PAR reaching the ground
        utils::tsdifstruct tsdPAR;  // diffuse two-stream parameters, PAR
    };

    // Solar zenith and azimuth for a location, date and local solar time.
    solmodel solpositionCpp(double lat, double lon, int year, int month, int day, double lt);

    // Solar zenith for flat-surface calculations, where azimuth does not
    // affect beam incidence. Azimuth fields are therefore left at zero.
    solmodel solzenithCpp(double lat, double lon, int year, int month, int day, double lt);

    // Establish the canopy's fixed diffuse optical state for total
    // shortwave and PAR, including clumping, albedo and transmission to the
    // ground.
    radsetup RadswabsSetupCpp(const vegpstruct& vegp, double gref, double grefPAR);

    // Shortwave balance for one timestep. Direct and diffuse radiation are
    // propagated through the clumped canopy to give absorption by the ground
    // and by the combined vegetation--ground surface, together with effective
    // albedo. The calculation is repeated with PAR optical properties for the
    // photosynthesis/stomatal model. Slope and aspect control direct-beam
    // incidence and sky-view factor controls exposure to diffuse sky
    // radiation. Terrain shadow is the caller's to apply, as forcing with no
    // direct component.
    void RadswabsStepCpp(const radsetup& s, const vegpstruct& vegp, double gref, double grefPAR,
        double slope, double aspect, double svfa, const solmodel& solp, envstruct& env);

    // Longwave balance of canopy and ground. Each surface exchanges with the
    // sky it can see, with the surrounding terrain filling the rest of the
    // upper hemisphere, and - for the ground - with the canopy over whatever
    // the canopy blocks. Vegetation and ground use the same thermal
    // emissivity. Writes the absorbed longwave of each surface and the
    // coefficient its own emission carries in its balance (env.emGround,
    // env.emCanopy); those are combined with the shortwave balance to give
    // the net radiative forcing of the surface energy balances.
    void RadlwabsStepCpp(const radsetup& s, double svfa, envstruct& env);

    // ========================================================================
    // Wind: within- and above-canopy wind speed and stability corrections
    // ========================================================================
    // Above-canopy aerodynamic state and canopy-top wind, with
    // Monin-Obukhov stability. Canopy height and plant area determine
    // displacement height and roughness, while sensible heat flux determines
    // atmospheric stability.
    //
    // Under unstable, weak-wind conditions a Beljaars free-convection velocity
    // scale is combined with the measured wind so turbulent exchange does not
    // collapse when buoyancy, rather than shear, is driving the flow. Stability
    // and sensible heat flux are coupled, so this state is updated as part of
    // the timestep's outer energy-balance iteration. `shelterc` optionally
    // reduces the driving wind for topographic shelter; bare ground has no
    // separate canopy-top or within-canopy mixing state.
    void windmodelCpp(const vegpstruct& vegp, double zref, double d, envstruct& env,
        double shelterc = 1.0, double zi = mc::freeConvZi, double beta = mc::freeConvBeta);

    // Aerodynamic resistance to heat transfer (s/m) from the surface
    // roughness length to an arbitrary height above the canopy: the
    // Monin-Obukhov-corrected logarithmic profile up to canopy top and the
    // transport column above it. This is the building block
    // for both canopy-to-atmosphere and ground-to-atmosphere exchange paths.
    double rHaToHeightCpp(const vegpstruct& vegp, double height, double d, const envstruct& env,
        const utils::ColumnStruct& col, bool hasColumn);

    // Bulk surface aerodynamic resistance to heat transfer (s/m), from the
    // reference height zref down to the roughness length for heat. Call
    // after windmodelCpp() has populated env.zm, env.LL and env.uf for this
    // timestep.
    double bulkaeroresistCpp(const vegpstruct& vegp, double zref, double d, const envstruct& env,
        const utils::ColumnStruct& col, bool hasColumn);

    // Total scalar resistance from the soil surface to the reference
    // atmosphere, through the transport column of utils.
    double groundaeroresistCpp(const vegpstruct& vegp, double zref, double d, const envstruct& env,
        const utils::ColumnStruct& col, bool hasColumn, double* rRise = nullptr);

    // Wind, air temperature and relative humidity at a requested height
    // above the canopy. Wind follows the Monin-Obukhov-corrected logarithmic
    // profile implied by the timestep's converged friction velocity and
    // stability state. Temperature is obtained by scaling the surface-to-air
    // temperature difference by the fraction of the aerodynamic resistance
    // accumulated below the requested height.
    //
    // Water vapour is treated with the same constant-flux resistance logic.
    // The exchange-surface vapour pressure is inferred from saturation at the
    // surface temperature and the combined aerodynamic + stomatal resistance,
    // then extrapolated to height with the aerodynamic resistance ratio. The
    // final conversion to relative humidity is bounded to 20--100% because
    // temperature and vapour pressure are extrapolated separately and can
    // otherwise combine to give unrealistic RH under strongly stratified
    // conditions. This routine is for z >= canopy height; within-canopy
    // profiles use the separate far-field calculation (utils::farFieldPairCpp).
    struct AboveCanopyPoint { double windz; double Tz; double RHz; };
    AboveCanopyPoint aboveCanopyProfileCpp(const vegpstruct& vegp, double z, double d,
        double zref, double Tsurface, const envstruct& env, double rHa, double rV, double hSurf,
        const utils::ColumnStruct& col, bool hasColumn, double LEtotal = std::numeric_limits<double>::quiet_NaN());

    // Temperature or volumetric water content at an arbitrary below-ground
    // depth, obtained by linear interpolation between the solved soil nodes,
    // which hold the state at depths z[0] (the surface) to z[n] (the fixed
    // lower boundary). Node depths are positive whereas reqhgt is negative
    // below ground, so the sign is converted internally. Requests at or above
    // the surface return the surface node, and requests below the deepest
    // node return that node, rather than extrapolating beyond the column.
    double interpSoilProfileCpp(const std::vector<double>& z,
        const std::vector<double>& values, double reqhgt);

    // ========================================================================
    // Stomatal conductance (JULES-SOX)
    // ========================================================================
    // Coupled photosynthesis, stomatal regulation and plant hydraulics.
    // Carbon gain favours opening stomata, while vapour-pressure deficit and
    // declining hydraulic conductance constrain water loss. Sunlit and shaded
    // foliage are solved separately because they experience different PAR,
    // then combined to give the canopy-scale stomatal resistance used by the
    // surface energy and water balances.

    // Stomatal conductance of a sunlit or shaded leaf fraction (mol m^-2
    // s^-1). Absorbed PAR and temperature determine photosynthetic demand,
    // while vapour-pressure deficit, plant height, root-zone water potential
    // and the xylem vulnerability curve determine the hydraulic cost of
    // supplying that transpiration.
    double leafgsCpp(const envstruct& env, vegpstruct& vegp, double z, bool C3 = true);

    // Radiation state needed to represent the canopy as sunlit and shaded
    // foliage. Solar geometry determines the exposed leaf-area fractions and
    // the canopy radiation balance determines absorbed PAR per unit leaf area;
    // these provide the contrasting light environments passed to leafgsCpp.
    struct stomatalsetup {
        double z;             // leaf height above ground for leafgsCpp's hydraulic head
                              // (mid-canopy, vegp.hgt / 2) (m)
        double L_sun;         // sunlit leaf area index (m^2/m^2)
        double L_shade;       // shaded leaf area index (m^2/m^2)
        double PARabs_sun;    // PAR absorbed per unit sunlit leaf area (W/m^2)
        double PARabs_shade;  // PAR absorbed per unit shaded leaf area (W/m^2)
    };

    // Partition the canopy into sunlit and shaded leaf area and assign PAR
    // absorption to each fraction for the current solar geometry. Defined only
    // for vegetated surfaces.
    stomatalsetup StomatalSetupCpp(const solmodel& solp, const envstruct& env, const vegpstruct& vegp);

    // Canopy-scale stomatal resistance (s/m), obtained by weighting the
    // conductance of sunlit and shaded foliage by their leaf areas. The light
    // partition is reconciled with the canopy's clumping-aware total PAR
    // absorption so photosynthetic water demand remains consistent with the
    // radiation balance.
    double bulkstomatalresistCpp(const stomatalsetup& setup, const envstruct& env, vegpstruct& vegp, bool C3 = true);

    // ========================================================================
    // Canopy surface energy balance and canopy water budget
    // ========================================================================
    // The bulk canopy exchange temperature is set by the surface energy
    // balance, while stomatal resistance and canopy wetness determine how the
    // available energy is partitioned between sensible and latent heat. Rain
    // interception changes that partition by allowing evaporation directly
    // from wet foliage as well as transpiration through stomata.

    // Canopy interception and wet-surface exchange for the big-leaf canopy.
    // Rain first fills the finite canopy water store; excess becomes
    // throughfall. The fraction of the surface that is wet lowers the bulk
    // surface resistance and therefore permits direct evaporation, while the
    // dry fraction retains the stomatal resistance used for transpiration.
    //
    // Vapour network. The foliage and the soil are two vapour sources, each at
    // its own temperature, exhaling into one node of canopy air that drains to
    // the reference height. The foliage's leg is the stomata (dry share) or
    // nothing (wet share) in series with the aerodynamic part rHa - rsh; the
    // soil's leg is RgR - rsh; the shared leg rsh runs from the node to the
    // reference, the node sitting at whichever of the mean source height and
    // the exchange surface is nearer the reference. Given the two surface
    // vapour pressures the node is a closed form. The whole surface's latent
    // flux is affine in es(Tc), (hSurf * es(Tc) - eV) / (rHa + rSurf), with
    // hSurf the foliage's share of the conductance and eV the air vapour
    // pressure less the soil's contribution; the soil solves see the network
    // as its Thevenin equivalent, air at eTh across rTh. With no foliage the
    // soil surface is the exchange surface and evaporates across rHa alone.
    struct canopywaterresult {
        double rSurf;        // resistance beyond rHa for the combined-surface
                              // Penman-Monteith solve (s/m) -- add to rHa for rV
        double hSurf;        // factor on the surface's saturation vapour pressure:
                              // the foliage's share of the network's conductance, or the
                              // soil's humidity over bare ground
        double eV = -999.0;  // reference vapour pressure the whole surface's balance sees (kPa);
                              // -999 over bare ground, where the air's own is used
        double eTh = -1.0;   // the network seen from the soil: air vapour pressure (kPa)
        double rTh = -1.0;   // and resistance (s/m); < 0 over bare ground
        double aNet = -1.0;  // foliage leg conductance a (m/s)
        double bNet = -1.0;  // soil leg conductance b (m/s)
        double swaterdepth;  // updated canopy water storage (mm per unit ground area),
                              // at most vegp.mwft per unit plant area
        double filmEvap = 0.0;      // evaporation of intercepted water over the step (mm)
        double filmAvailable = 0.0; // water the film could give up over the step (mm)
        double wetShare = 0.0;      // share of the surface that is wet at the start of the step
        double Et;           // stomatal transpiration (mm), kept separate from
                              // evaporation of intercepted canopy water and used as
                              // the demand for root water uptake
        double precipGround; // rainfall reaching the ground surface (mm), after
                              // interception/throughfall -- feeds env.precipGround
    };
    // `eg` is the soil surface vapour pressure (kPa); `rsh` is the
    // node-to-reference resistance (s/m; 0 = no column, so the node is the
    // reference air); `RgR` is the soil-to-reference resistance (s/m), so the
    // soil-to-node leg is `RgR - rsh`.
    canopywaterresult canopyWaterBudgetCpp(vegpstruct& vegp, const envstruct& env,
        const stomatalsetup& stomSetup, bool vegetated, double rHa, double hr,
        double swaterdepth, double dT, double eg, double rsh, double RgR);

    // ========================================================================
    // Soil: heat and water balance for a layered soil profile
    // ========================================================================
    // A layered soil profile is described by soilpstruct (time-invariant
    // physical properties per layer) together with two pieces of
    // time-varying state: soilheatmod (temperature profile, advanced by
    // SoilHeatCpp) and soilwatermod (water potential/content profile,
    // advanced by SoilWaterCpp). Both solvers read the meteorological
    // drivers and the ground-to-reference-height aerodynamic resistance
    // needed for the surface energy/water balance from envstruct
    // (env.RabsGround, env.tair, env.rh, env.pk, env.rHa, env.precipGround,
    // env.Et -- SoilWaterCpp's surface flux uses precipGround, the
    // post-interception throughfall, not the raw precip input; see
    // canopyWaterBudgetCpp above).

    // Physical properties of a layered soil profile, constant across an
    // entire model run.
    struct soilpstruct {
        // Terrain geometry does not enter the layered soil properties;
        // topographic effects on radiation and moisture are handled outside
        // this soil-parameter structure.
        double gref;     // ground reflectance, shortwave (0-1)
        double grefPAR;  // ground reflectance, PAR (0-1)
        // Ground longwave emissivity uses the shared surfaceEmissivity constant.
        int nLayers;     // number of soil layers
        bool FreeDrain;  // whether the bottom layer is free-draining

        std::vector<double> Vq;      // volumetric quartz fraction (m^3/m^3)
        std::vector<double> Vm;      // volumetric other-mineral fraction (m^3/m^3)
        std::vector<double> Vo;      // volumetric organic matter fraction (m^3/m^3)
        std::vector<double> Mc;      // mass fraction of clay (kg/kg)
        std::vector<double> psie;    // air-entry water potential (J/kg)
        std::vector<double> b;       // Campbell water retention curve shape parameter
        std::vector<double> thetaR;  // residual volumetric water fraction, used to set the initial water content
        std::vector<double> thetaS;  // volumetric water fraction at saturation
        std::vector<double> Ksat;    // saturated hydraulic conductivity (kg s / m^3)
        // In the Campbell formulation the hydraulic-conductivity exponent
        // is derived from b as 2*b + 3 wherever it is required.
        std::vector<double> psi_min; // oven-dry water potential, the solver's lower bound
    };

    // Time-varying state of the soil temperature profile.
    struct soilheatmod {
        int n;                       // number of layers
        std::vector<double> z;       // node depths (m); node i holds the state at z[i]
        std::vector<double> dz;      // spacing between node i and node i + 1 (m)
        std::vector<double> vol;     // control volume around each node (m^3, per unit ground area)
        std::vector<double> wc;      // volumetric water fraction at each node
        std::vector<double> Te;      // temperature at each node (deg C)
        std::vector<double> oldTe;   // temperature at each node at the previous timestep (deg C)
        int iters = 0;               // number of iterations taken to converge
        double lamS = 0.0;           // restoring slope of the surface row's balance, -dq/dTs (W/m^2/K)
        double condSlope = 0.0;      // surface row's conduction slope, C0/dT + f0 (W/m^2/K)
    };

    // Time-varying state of the soil water profile.
    struct soilwatermod {
        int n;                        // number of layers
        std::vector<double> z;        // node depths (m); node i holds the state at z[i]
        std::vector<double> dz;       // spacing between node i and node i + 1 (m)
        std::vector<double> vol;      // volume around each node (m^3, per unit ground area)
        std::vector<double> psiw;     // water potential of each layer (J/kg)
        std::vector<double> k;        // hydraulic conductivity of each layer
        std::vector<double> vapor;    // water vapour concentration of each layer
        std::vector<double> oldvapor; // water vapour concentration at the previous timestep
        std::vector<double> theta;    // volumetric water fraction of each layer
        std::vector<double> oldtheta; // volumetric water fraction at the previous timestep
        std::vector<double> Tc;       // temperature of each layer (deg C)
        std::vector<double> oldTc;    // temperature of each layer at the previous timestep (deg C)
        std::vector<double> rootfrac; // fraction of root water uptake taken from each layer
    };

    // Result of a single call to SoilWaterCpp.
    struct soilwaterresult {
        soilwatermod state;  // updated soil water state
        bool success;        // whether the Newton iteration converged
        int iterations;      // number of iterations taken
        double Evapmmhr;     // bare-soil evaporation for this timestep (mm)
        double surplus;      // water a saturated surface could not accept (mm)
    };

    // Relative humidity (0--1) of soil pore air at the surface. The top
    // layer's volumetric water content is converted to matric water potential
    // with the Campbell retention curve, and the Kelvin equation then gives
    // the equilibrium vapour pressure relative to saturation at the soil
    // temperature. This is the humidity boundary used for soil evaporation.
    double soilrelhumCpp(const soilpstruct& soilp, double Tsoil, double theta);

    // Net radiation balance (W/m^2, positive into the surface) of the soil
    // surface at temperature Tsurface (deg C) and water content theta
    // (m^3/m^3), combining absorbed radiation (env.RabsGround), emitted
    // longwave (through env.emGround, the coefficient that belongs with it),
    // sensible heat exchange with the reference-height air
    // (via env.rHa) and latent heat exchange with the surface pore-space
    // humidity (soilrelhumCpp) -- the ground-surface analogue of the
    // canopy's own surface energy balance, but returning the raw imbalance
    // Ba = Rnet - H - L rather than solving for Tsurface directly (that
    // solve is penmanMonteithCpp's job; this function is used instead
    // where the residual itself, not an updated temperature, is what's
    // needed, e.g. inside the soil heat solver).
    double soilsurfaceEBCpp(const soilpstruct& soilp, const envstruct& env, double Tsurface, double theta);

    // Ground-surface equilibrium temperature (deg C) if subsurface heat
    // storage were zero (G = 0). Holding the solved surface soil moisture and
    // atmospheric/aerodynamic drivers fixed, this finds the temperature at
    // which absorbed radiation is balanced by sensible and latent heat loss.
    // It provides an unbuffered surface-temperature forcing against which the
    // effect of soil heat storage can be scaled in the grid model.
    double groundTemp0Cpp(const soilpstruct& soilp, const envstruct& env, double theta,
        double Tguess, int maxIter = 50, double tolerance = 1e-4);

    // Ground-surface temperature (deg C) consistent with a prescribed ground
    // heat flux G (W/m^2, positive into the soil). The surface moisture and
    // atmospheric/aerodynamic state are held fixed while temperature is solved
    // so that the residual surface energy balance equals G. This allows a
    // spatially estimated heat-storage flux to be translated back into a
    // physically consistent surface temperature without re-solving the full
    // soil column.
    double groundTempGCpp(const soilpstruct& soilp, const envstruct& env, double theta,
        double G, double Tguess, int maxIter = 50, double tolerance = 1e-4);

    // Diurnal thermal damping depth (m), sqrt(2*kappa/omega), where kappa is
    // the top soil layer's thermal diffusivity and omega is the angular
    // frequency of the 24-hour cycle. It measures how deeply daily surface
    // temperature fluctuations penetrate and therefore provides the soil-
    // dependent thermal response scale used by the grid heat-storage
    // approximation.
    double diurnalDampingDepthCpp(const soilpstruct& soilp, double theta, double Tc, double pk);

    // Effective degree of saturation (0--1) from water potential using the
    // Campbell retention curve. At or wetter than the air-entry potential the
    // layer is treated as saturated (Se = 1); below air entry, saturation falls
    // as a power law controlled by the pore-size parameter b.
    double degreeOfSaturationCpp(const soilpstruct& soilp, double psiw, int i);

    // Soil matric water potential (J/kg) from volumetric water content using
    // the Campbell retention curve. Water content is capped at saturation
    // before applying the power law, so this relation represents the
    // unsaturated/saturated retention curve but not positive pressure heads.
    double waterPotentialCpp(const soilpstruct& soilp, double theta, int i);

    // Root-zone mean water potential (MPa), weighting each soil layer's
    // Campbell water potential by the fraction of roots assigned to that
    // layer. This is the plant-accessible water-status signal supplied to the
    // stomatal/hydraulic model rather than an unweighted profile average.
    double rootzonePsiCpp(const soilpstruct& soilp, const std::vector<double>& theta, const std::vector<double>& rootfrac);

    // Volumetric water content (m^3/m^3) from water potential, the inverse
    // Campbell retention relationship used by waterPotentialCpp.
    double thetaFromPsiCpp(const soilpstruct& soilp, double psiw, int i);

    // Unsaturated hydraulic conductivity of a soil layer from volumetric
    // water content, using the Campbell relationship
    // K = Ksat * (theta/thetaS)^(2*b+3). Deriving the exponent from the same
    // b used in the retention curve keeps water retention and conductivity
    // internally consistent.
    double hydraulicConductivityFromThetaCpp(const soilpstruct& soilp, double theta, int i);

    // Water-vapour mass per bulk soil volume from pore-space humidity and
    // air-filled porosity. Matric potential sets equilibrium relative humidity
    // through the Kelvin equation, temperature sets saturated vapour density,
    // and (thetaS - theta) supplies the volume available to the gas phase.
    double vaporFromPsiCpp(const soilpstruct& soilp, double psiw, double theta, double Tk, int i);

    // Differential soil-water capacity d(theta)/d(psi), obtained analytically
    // from the Campbell retention curve. This is the liquid-water storage
    // derivative used in the Newton solve; it is zero on the saturated side
    // where the model clamps effective saturation at one.
    double dThetaDPsiCpp(const soilpstruct& soilp, double psiw, int i);

    // Vapour-phase hydraulic conductivity (same units as
    // hydraulicConductivityFromThetaCpp's liquid-phase k) of soil layer i
    // from its water potential psiw (J/kg), water content theta (m^3/m^3)
    // and temperature Tk (Kelvin): k_vapour = 0.66 * (thetaS[i]-theta) *
    // Dv * rho_vs(Tk) * hr * Mw/(RgasC*Tk), where Dv (0.000024 m^2/s)
    // is the vapour diffusivity in soil air and 0.66 represents the reduced
    // effective diffusivity of vapour moving through tortuous pore space
    // (Campbell 1985; Bittelli et al. 2015). Added, in series/parallel
    // with the liquid-phase conductivity, to give SoilWaterCpp's total
    // layer conductivity for the mixed-form Richards equation.
    double vaporConductivityFromPsiThetaCpp(const soilpstruct& soilp, double psiw, double theta, double Tk, int i);

    // Sensitivity of soil-vapour storage to water potential. The derivative
    // combines the Kelvin-equation response of pore humidity with the change
    // in air-filled pore volume as liquid water content changes, and provides
    // the vapour-storage Jacobian term in the soil-water solve.
    double dVaporDPsiCpp(const soilpstruct& soilp, double psiw, double theta, double Tk, int i);

    // Bare-soil evaporation over the timestep (mm), driven by the vapour-
    // pressure difference between soil pore air and the reference atmosphere
    // across the ground aerodynamic resistance. Surface pore humidity follows
    // the Kelvin relation from the top-layer water potential, evaluated at the
    // surface temperature as the surface energy balance evaluates it.
    double evaporationFluxCpp(const soilpstruct& soilp, const envstruct& env, double theta, double Tsurface, double dT);

    // Distributes the canopy's total transpiration demand among soil layers.
    // Potential uptake is proportional to root fraction and is reduced when a
    // layer is either too wet or too dry for effective uptake; the remaining
    // weights are renormalised so available layers share the demand. If every
    // layer is fully stressed, uptake is zero rather than forcing extraction
    // from unavailable water.
    std::vector<double> transpirationDistributeCpp(const soilpstruct& soilp, const std::vector<double>& rootfrac,
        double totalTransp_mm, double dT, const std::vector<double>& psiw, double p = 0.5);

    // Advances the multilayer soil-temperature profile by one timestep.
    // Heat diffusion is solved implicitly with temperature- and moisture-
    // dependent thermal properties, while the upper boundary is coupled to
    // the radiative, sensible and latent surface energy balance. The nonlinear
    // problem is iterated to convergence; the deepest boundary temperature is
    // held fixed over the timestep. Conduction is weighted in time the same way
    // in every row, including the surface row, so that the heat one row sends
    // between two nodes is the heat the next row receives. Unlike the water solve
    // (which offers a free-drainage lower boundary via soilp.FreeDrain), there
    // is no free-flux option for heat -- the fixed deep temperature is the
    // only bottom boundary condition available, so a profile that is shallow
    // relative to the timestep can show an artificially damped deep-layer
    // response.
    //
    // The surface row is implicit by construction, so Fact must remain 1. Any
    // other interior weighting breaks row-to-row heat conservation and leaks
    // energy.
    //
    // For vegetation shorter than 10 mm the surface row solves a blend of two
    // balances (the matching constraint): the ground's own (env) and the whole
    // surface's (envS), weighted w and 1 - w. See matchResidualCpp.
    soilheatmod SoilHeatCpp(soilheatmod state, const soilpstruct& soilp, const envstruct& env,
        double dT = 3600.0, double Fact = 1.0, int maxNrIterations = 100, double tolerance = 1e-2,
        const envstruct* envS = nullptr, double wMatch = 1.0);

    // The surface row's flux under the matching constraint: the blend of two
    // energy-balance residuals, each zero at its own solution,
    //     w (q_g - C) + (1 - w) kappa (q_s - C) = 0 ,
    // with kappa = (lam_g + C')/(lam_s + C') equalising their restoring
    // strengths. Solved for C this is a flux whose weights sum to one, so a
    // flux level the two share passes through unchanged at every w; the slope
    // is the same average. q_g is the ground's balance (env), q_s the whole
    // surface's (envS), C' the row's conduction slope.
    struct matchresidual { double q = 0.0; double lam = 0.0; double kappa = 1.0; };
    matchresidual matchResidualCpp(const soilpstruct& soilp, const envstruct& env, const envstruct* envS,
        double w, double T, double theta, double Cp);

    // Advances the multilayer soil-water profile by one timestep using water
    // potential as the solved state variable. Liquid and vapour transport are
    // coupled to precipitation/throughfall, bare-soil evaporation and root
    // uptake, and the nonlinear storage/flux equations are solved iteratively.
    // With free drainage, water may leave the bottom under gravity; otherwise
    // the lower boundary is held at saturation, representing a fixed water-
    // table boundary. The Newton iteration uses the exact derivative of the
    // mass-balance residual and takes full steps, halving a step only when it
    // fails to reduce the residual. The tolerance is on the profile-wide
    // mass-balance residual (kg m-2 s-1); 1e-6 is 0.0036 mm per hour.
    soilwaterresult SoilWaterCpp(soilwatermod state, const soilpstruct& soilp, const envstruct& env,
        double dT = 3600.0, double pTAW = 0.5, int maxNrIterations = 100, double tolerance = 1e-6);

    // The timestep orchestration that couples these component models and
    // advances their shared state is implemented in pointmodel.cpp.

} // namespace pointmodel
