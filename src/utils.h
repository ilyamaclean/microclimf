// utils.h
// General-purpose numerical and physical helper functions shared across the
// point and grid model code.
#pragma once
#include <string>
#include <vector>

namespace utils {

    // Result of a tridiagonal solve via the Thomas algorithm: the
    // (forward-eliminated) diagonal and right-hand-side coefficients, and
    // the solution vector x.
    struct ThomasResult {
        std::vector<double> bb;
        std::vector<double> cc;
        std::vector<double> dd;
        std::vector<double> x;
    };

    // Zero-plane displacement height (m): the effective height at which a
    // canopy's drag on the wind above it appears to originate, so that the
    // above-canopy wind profile behaves as a standard logarithmic profile
    // measured from this height rather than from the ground. Depends on
    // canopy height (h, m) and plant area index (pai, m^2/m^2), following
    // Raupach (1994) Eq. 8. A denser or taller canopy pushes this
    // "effective ground" higher into the canopy. pai is floored at a small
    // positive value before use, so bare ground (pai = 0) still returns a
    // small, physically sensible displacement rather than triggering a
    // divide-by-zero.
    double zeroplanedisCpp(double h, double pai);

    // Ratio of friction velocity to wind speed at canopy top, from the
    // partition of surface drag between the substrate and the roughness
    // elements standing on it (Raupach 1994 Eq. 7). Rises with plant area
    // index until elements begin to shelter one another, beyond which it is
    // held at the limit his Sect. 4c takes from drag data.
    double dragPartitionCpp(double pai);

    // Aerodynamic closure of the canopy on the drag partition's own scale:
    // 0 over bare substrate, 1 once the partition saturates.
    double dragClosureCpp(double pai);

    // Momentum roughness length (m), including the momentum-sublayer influence
    // (Raupach 1992 Eq. 27); returns bareGroundZ0 for h <= 0.
    double roughlengthCpp(double h, double pai, double d);

    // Heat roughness length (m): 0.2 times the momentum construction without
    // any roughness-sublayer influence. The 0.2 factor is exp(-kB^-1), with
    // kB^-1 = ln 5. The sublayer is integrated above canopy top by the transport
    // column, so including it here would count it twice.
    double scalarRoughlengthCpp(double h, double pai, double d);

    // Raupach (1992) influence function Psi(c) = ln(c) - 1 + 1/c; zero for c <= 1.
    double sublayerInfluenceCpp(double c);

    // Influence function of the momentum roughness sublayer, whose depth
    // coefficient follows canopy closure as the scalar sublayer's does:
    // Psi(1 + (cW - 1) f). Zero over bare substrate, Psi(2) = 0.1931 once
    // closed. Folded into roughlengthCpp and taken back out at canopy top.
    double momentumInfluenceCpp(double pai);
    // Top of that sublayer, z_w = d + (1 + (cW - 1) f)(h - d).
    double momentumSublayerTopCpp(double h, double pai, double d);
    // Height at which the flow has cleared both roughness sublayers, the top
    // of the scalar one, d + (1 + (cH - 1) f)(h - d); zero without a canopy.
    double rslClearHeightCpp(double h, double pai);

    // Molar density of air (mol/m^3) at a given temperature (tc, deg C)
    // and atmospheric pressure (pk, kPa), from the ideal gas law
    // (Campbell & Norman 2012). Used wherever a flux needs to be expressed
    // per mole of air rather than per unit mass -- the two are not
    // interchangeable without this conversion, and should not be confused
    // with a mass density (kg/m^3) computed elsewhere in this model.
    double phairCpp(double tc, double pk);

    // Molar specific heat of air at constant pressure (J / mol / K), as a
    // function of air temperature (tc, deg C): how much energy is needed
    // to raise the temperature of one mole of air by one degree.
    // A quadratic fit to measured air properties (Campbell & Norman 2012).
    double cpairCpp(double tc);

    // Integrated stability correction for momentum transfer, as a function
    // of the dimensionless stability parameter z/L (height over the
    // Monin-Obukhov length, which measures how strongly buoyancy is
    // enhancing or suppressing turbulent mixing). This correction adjusts
    // the standard neutral logarithmic wind profile for the effect of
    // atmospheric stability: under unstable (convective) conditions
    // (z/L < 0) turbulence is enhanced and wind shear reduced, following
    // Businger et al. (1971); under stable conditions (z/L >= 0) turbulence
    // is suppressed and wind shear increased. The stable branch approaches
    // a fixed asymptotic correction rather than growing without bound as
    // stability increases, which keeps the model numerically well-behaved
    // under very calm, strongly stable conditions (e.g. clear, still
    // nights) without materially changing the correction under typical
    // conditions.
    double dpsimCpp(double ze);

    // Integrated stability correction for heat transfer, the heat-transfer
    // counterpart to dpsimCpp: same physical role (adjusting the neutral
    // profile for atmospheric stability), same Businger et al. (1971)
    // unstable branch and asymptotically-bounded stable branch, but with
    // its own turbulent Prandtl-number scaling (0.74) reflecting that heat
    // and momentum are not transported with identical efficiency by
    // turbulent eddies.
    double dpsihCpp(double ze);

    // Bounds the Monin-Obukhov length L (the height scale characterising
    // atmospheric stability) so that its associated momentum stability
    // correction cannot exceed a fixed multiple of the neutral logarithmic
    // profile's own magnitude, on either the unstable or the stable side.
    // This is a numerical safeguard, not a physical law: without it, an
    // outer wind/stability solve can enter a runaway feedback loop under
    // strongly unstable or strongly stable conditions (very calm wind,
    // strong heating or strong nocturnal cooling), where a small change in
    // one quantity drives an ever-larger change in the other. `beta`
    // controls how much stability correction is allowed, symmetrically,
    // before L is pulled back towards the neutral (unstressed) case.
    double clipMOlength(double L, double zref, double d, double zm, double beta = 0.9);

    // Exact, closed-form recovery of the bulk momentum stability
    // correction psi_m implied by an already-known friction velocity uf
    // and effective driving wind speed Ueff, by inverting the logarithmic
    // wind-profile equation itself. Used where a full iterative
    // stability solve has already been done elsewhere (so uf is already
    // known) and the corresponding correction just needs to be read back
    // out, rather than iterated for independently.
    double recoverPsiMCpp(double uf, double Ueff, double zref, double d, double zm);

    // Inverts the stable-branch (L > 0) relationship between the Monin-
    // Obukhov length L and the momentum stability correction psi_m it
    // implies, recovering L from a target psi_m. This inversion is not
    // one-to-one everywhere: because the stable correction approaches the
    // same asymptotic value from both the near-neutral and strongly-stable
    // limits, a given psi_m below a certain threshold can correspond to
    // two different L values, and above that threshold no exact L exists
    // at all. This function always returns the more conservative
    // (larger-L, closer to neutral) of the two possible roots, and
    // saturates to the boundary value if the requested psi_m cannot be
    // reached exactly. Returns +Infinity for a neutral-or-better (psi_m <= 0)
    // target, for which no finite stable-branch L applies.
    double lStableFinalCpp(double psi_m, double zref, double d, double zm);

    // Inverts the unstable-branch (L < 0) relationship between the
    // Monin-Obukhov length and the momentum stability correction, using a
    // closed-form starting estimate refined by a small number of
    // Newton-Raphson iterations. Unlike the stable branch, this
    // relationship is one-to-one, so no ambiguity in the recovered L
    // arises here. Returns -Infinity if the requested correction is too
    // small to distinguish from the neutral case.
    double lUnstableNewtonCpp(double psi_m, double zref, double d, double zm, int nNewton = 1);

    // Recover L from friction velocity and driving wind by inverting the bulk
    // momentum correction, not by an outer stability solve. The stable branch
    // is solved numerically; see recoverLClosedCpp for the closed form used by
    // the grid model.
    double recoverLCpp(double uf, double Ueff, double zref, double d, double zm);

    // Peak of the stable-branch bulk momentum correction, fixed by geometry
    // and returned as u = a/L, with a = 1/stableZetaMaxM. A closed-form start
    // is refined by two Newton steps.
    double stableRecoveryPeakCpp(double zref, double d, double zm);

    // Recover L on the stable branch from psi_m and the geometry's peak
    // (from stableRecoveryPeakCpp): a closed-form start that neglects the
    // roughness term, refined by two Newton steps kept between zero and the
    // peak. Returns the milder root, as lStableFinalCpp does, and the peak's
    // own L for a correction at or beyond it. Within 1% of the exact root in
    // the heat resistance for zr/zm >= 5, and within 0.14% up to 90% of
    // the peak.
    double lStableClosedCpp(double psi_m, double zref, double d, double zm, double uPeak);

    // As recoverLCpp, but with the closed-form stable branch above, for the
    // grid model.
    double recoverLClosedCpp(double uf, double Ueff, double zref, double d, double zm, double uPeak);

    // Free-convective velocity scale w* (Beljaars 1994) for a surface giving
    // sensible heat H (W/m2) to air at tc (deg C) and pk (kPa), over a
    // boundary layer of depth mc::freeConvZi; zero unless H > 0.
    double freeConvectiveVelocityCpp(double H, double tc, double pk);

    // The wind driving the surface layer, as windmodelCpp forms it: the
    // (sheltered) wind combined in quadrature with the free-convective
    // velocity, then floored at mc::minWindSpeed.
    double drivingWindCpp(double u, double wstar);

    // Canopy light-extinction coefficient for direct-beam radiation, as a
    // function of solar zenith angle (zenr, radians) and the leaf angle
    // distribution parameter x -- how strongly a canopy's own foliage
    // angle statistics concentrate its shading rather than passing light
    // straight through, following Campbell (1986)'s closed-form
    // approximation for an ellipsoidal leaf angle distribution. x = 1
    // describes a spherical (random) leaf orientation; x approaching
    // infinity describes horizontal (planophile) leaves; x = 0 describes
    // vertical (erectophile) leaves -- these three limiting cases are
    // handled explicitly, since the general formula is undefined there.
    double canopyKCpp(double zenr, double x);

    // ------------------------------------------------------------------
    // Two-stream canopy radiation transfer (Sellers 1985). This solves for
    // how much shortwave/PAR radiation a canopy of given optical properties
    // and structure absorbs, reflects and transmits, separately for the
    // diffuse (sky) and direct-beam (sun) components. Shared by the point
    // and grid models -- the physics here has no knowledge of either
    // model's own vegetation/soil representation.
    // ------------------------------------------------------------------

    // Canopy extinction coefficient (see canopyKCpp), together with two
    // path-length adjustments needed for tracking radiation travelling at
    // an angle through the canopy rather than straight down. kd rescales
    // the extinction coefficient for the true (possibly oblique) path
    // length of direct-beam radiation, for use in the direct-beam
    // two-stream solve. Kc is the equivalent path-length correction used
    // to convert plant area into an optical depth along that same path for
    // simpler gap-fraction (Beer's law) calculations. Both fall back to
    // fixed values when the sun is exactly on or below the horizon, rather
    // than dividing by zero.
    struct kstruct {
        double k;    // canopy extinction coefficient (vertical)
        double kd;   // extinction coefficient adjusted for path length (si)
        double Kc;   // path-length correction factor (1/si)
    };
    kstruct cankCpp(double zenr, double x, double si);

    // Two-stream parameters governing how diffuse (sky) radiation is
    // scattered, absorbed and transmitted through the canopy. These depend
    // only on the canopy's own optical properties and total leaf/plant
    // area, not on where the sun is, so they are computed once per canopy
    // state and reused for both the diffuse solve and, as inputs, for the
    // direct-beam solve below.
    struct tsdifstruct {
        double p1, p2, p3, p4;
        double om, a, gma, J, del, h, u1, S1, D1, D2;
    };
    tsdifstruct twostreamdifCpp(double pait, double x, double lref, double ltra, double gref);

    // Two-stream parameters governing how direct-beam solar radiation is
    // scattered, absorbed and transmitted through the canopy, building on
    // the diffuse optical solution and the solar-position-dependent beam
    // extinction coefficient. This implementation corrects a known error in
    // the denominator term as originally published by Sellers (1985): the
    // published form can allow beam-generated diffuse fluxes to become
    // unphysically negative, so the denominator here is amended to keep
    // them non-negative; downstream radiation fractions are additionally
    // bounded to their physically admissible range as a further safeguard.
    struct tsdirstruct {
        double sig, p5, p6, p7, p8, p9, p10;
    };
    tsdirstruct twostreamdirCpp(double pait, double om, double a, double gma, double J, double del, double h,
        double gref, double kd, double u1, double S1, double D1, double D2);

    // Combines a radiative or transmission quantity's two possible fates
    // in a canopy with structural gaps: a `gapfrac` share of it passes
    // straight through an open gap, retaining whatever value it already
    // had (`bypass` -- e.g. an unattenuated transmission, or radiation
    // reflected off the ground and straight back out through the gap);
    // the remaining share instead passes through the closed, "clean"
    // canopy and takes the value the ordinary two-stream calculation gives
    // it (`clean`). Every clumped-canopy radiation quantity (albedo,
    // ground transmission, angle of the beam reaching the ground, and
    // their longwave equivalents) follows this same gap/no-gap blend, with
    // a different gap fraction and pair of values supplied at each call
    // site.
    double canopyGapMixCpp(double gapfrac, double bypass, double clean);

    // Saturated vapour pressure of air (kPa) at a given temperature (tc,
    // deg C): the maximum water vapour pressure the air can hold before
    // condensation begins, following the Tetens (1930) formula, with
    // separate empirical constants above and below freezing to represent
    // saturation over liquid water versus over ice respectively. The
    // branch is chosen from the air temperature itself, so this always
    // represents saturation over liquid water above 0 deg C and over ice
    // at or below it.
    double satvapCpp(double tc);

    // Saturated water vapour density (kg/m^3) at a given temperature (Tk,
    // Kelvin): the mass of water vapour per unit volume of air at
    // saturation, obtained by applying the ideal gas law to the saturation
    // vapour pressure (satvapCpp). Because this varies strongly with
    // temperature (roughly ten-fold across a typical 0-40 deg C range), it
    // is computed at the actual surface/air temperature rather than
    // assumed constant, wherever it drives a vapour-flux calculation (e.g.
    // soil evaporation).
    double satVapDensityCpp(double Tk);

    // Astronomical Julian day number for a given calendar date -- a
    // continuous day count used as the common time base for solar-position
    // calculations in both the point and grid models.
    int juldayCpp(int year, int month, int day);

    // Apparent solar time (hours) from Julian day, clock time and longitude.
    // The clock time must be UTC, or use a timezone consistent with longitude;
    // otherwise the result can be wrong by whole hours.
    double soltimeCpp(int jd, double lt, double lon);

    // Fraction of direct-beam solar radiation intercepted by a surface of
    // given slope and aspect, for a given solar position: the standard
    // cosine-of-incidence-angle law for a tilted plane, which reduces to
    // the familiar cosine of the solar zenith angle for a flat surface.
    // Clamped to zero for a surface facing away from the sun, rather than
    // returning an unphysical negative flux.
    //
    // `shadow` supplies the result of an actual terrain/horizon shadow
    // test for this location and moment (the sun sitting below the local
    // skyline, as determined elsewhere from real topography) -- this
    // function has no terrain data of its own to test that against. When
    // `shadow` is true, the intercepted fraction is forced to exactly
    // zero, ensuring every downstream direct-beam quantity computed from
    // it is consistently zero for a genuinely shadowed surface, rather
    // than only the geometric incidence fraction being zeroed while other
    // direct-beam terms computed independently remain nonzero.
    double solarindexCpp(double slope, double aspect, double zend, double azid, bool shadowmask = false, bool shadow = false);

    // Combines two adjoining materials' thermal (or other) conductivities
    // into a single effective value at their shared interface, using
    // either the geometric mean or the logarithmic mean of the two.
    // Needed wherever two adjacent soil layers of different composition
    // meet, so that heat conduction across the boundary between them can
    // be represented by one effective value rather than an ill-defined
    // discontinuity (Campbell 1985; Bittelli et al. 2015).
    double kMeanCpp(const std::string& meanType, double k1, double k2);

    // Solves a tridiagonal system of linear equations (as arises from an
    // implicit finite-difference discretisation of a 1-D diffusion
    // process, e.g. heat or water flow through a layered soil column) via
    // the Thomas algorithm (Thomas 1949) -- an efficient, direct
    // (non-iterative) solution method exploiting the tridiagonal
    // structure. Used by both the soil heat and soil water solvers to
    // advance their respective layered profiles by one timestep.
    ThomasResult thomasSolveCpp(std::vector<double> aa, std::vector<double> bb,
        std::vector<double> cc, std::vector<double> dd, std::vector<double> x,
        int first, int last);

    // Node depths (m) for a layered soil profile of n layers spanning a
    // total depth, with layers becoming progressively thicker with depth
    // -- concentrating resolution near the surface, where temperature and
    // moisture vary fastest, while still representing the deeper profile
    // economically. Returns n + 2 depths: the surface, a thin near-surface
    // node, and the remaining layer boundaries.
    std::vector<double> geometricCpp(int n, double totalDepth);

    // Thermal conductivity (W/m/K) of a soil (or other porous medium),
    // from the volumetric fractions of its mineral, organic and water
    // content and the properties of any remaining air-filled pore space,
    // following the de Vries (1963)/Campbell (1985) mixing model for
    // granular media. Quartz is treated separately from other mineral
    // material because it conducts heat substantially better (roughly
    // 3.5-fold) than typical non-quartz minerals, so the two fractions are
    // not interchangeable -- a soil's actual quartz content materially
    // affects how efficiently heat moves through it.
    double thermalConductivityCpp(double Vq, double Vm, double Vo, double Vw,
        double Mc, double Tc, double pk);

    // Volumetric heat capacity (J/m^3/K) of a soil (or other porous
    // medium): how much energy is needed to raise a unit volume of the
    // soil's mineral, organic, water and air content by one degree,
    // combining each component's own specific heat weighted by its
    // volumetric share (Campbell & Norman 2012; de Vries 1963). Air's
    // contribution is computed from its actual temperature and pressure
    // rather than a fixed value, since air's heat capacity is comparatively
    // sensitive to both.
    double heatCapacityCpp(double Vq, double Vm, double Vo, double Vw, double Tc, double pk);

    // The (temperature in Kelvin)^4 term of the Stefan-Boltzmann emitted-
    // radiation law for a surface at temperature tc (deg C). Multiplying
    // this by the Stefan-Boltzmann constant and a surface's emissivity
    // gives the actual emitted longwave flux (W/m^2); this function alone
    // returns only the temperature term, not a flux.
    double rademCpp(double tc);

    // Latent heat of vaporisation above 0 deg C and sublimation below (J/mol),
    // at surface temperature Tsurface (deg C).
    double latentHeatCpp(double Tsurface);

    // One linearised Penman-Monteith surface-temperature update.
    // `rV` is vapour resistance; `hr` scales the surface saturation vapour
    // pressure. `eaRef` replaces the rh-derived vapour pressure only when
    // eaRef > -998. `gSlope` (W/m^2/K) is dG/dTs at the current Ts.
    double penmanMonteithCpp(double Rabs, double Ta, double pk, double rh, double em,
        double rHa, double rV, double Ts, double G = 0.0, double hr = 1.0, double eaRef = -999.0,
        double gSlope = 0.0);

    // Water stress factor (0-1) representing restricted gas exchange in a
    // waterlogged surface: unstressed (1) once soil water potential is at
    // or drier than the air-entry potential (the point at which the
    // largest soil pores begin to drain), falling to fully stressed (0) at
    // full saturation. A saturated, air-free soil restricts evaporation
    // and root gas exchange just as a very dry soil does, for a different
    // physical reason (no air-filled pore space rather than no available
    // water) -- this is the wet-side complement to alphaDryCpp below.
    double alphaWetCpp(double psiw, double psie);

    // Water stress factor (0-1) representing reduced plant/soil water
    // availability as soil dries: unstressed (1) at or above an
    // onset-of-stress water potential, declining linearly to fully
    // stressed (0) at the wilting point -- a simple, standard
    // computationally-cheap stand-in for a fuller plant hydraulic
    // response.
    double alphaDryCpp(double psiw, double psiDry, double psiWilt);

    // Root fractions from the Schenk & Jackson (2002) profile. `z` contains
    // node depths 0..n plus one depth below the boundary. Each node takes the
    // roots between halfway points to its neighbours; node n-1 takes all deeper
    // roots and boundary node n takes none. Fractions sum to one.
    std::vector<double> rootDistributeCpp(const std::vector<double>& z, double D50, double D95);

    // State carried between successive calls to aitken1d(), for damped
    // relaxation of a single scalar fixed-point iterate.
    struct Aitken1DState {
        double omega = 0.5;    // current relaxation factor, self-tuning
        double r_prev = 0.0;   // previous call's residual (newv - oldv)
        bool have_prev = false;
    };

    // Damped (Aitken-relaxed) update of a quantity being solved for by
    // fixed-point iteration: given the previously accepted value and a
    // fresh, undamped candidate, returns a blend of the two rather than
    // jumping straight to the new value. The blending strength is
    // re-learned on each call from how quickly the iteration is
    // converging, which stabilises quantities (such as ground heat flux
    // within the canopy energy-balance iteration) that would otherwise
    // oscillate under naive iterative updating.
    double aitken1d(double oldv, double newv, Aitken1DState& st);

    // Dew point temperature (deg C): the temperature to which air at a
    // given temperature and relative humidity would need to cool for its
    // water vapour content to reach saturation, found by inverting the
    // saturation vapour pressure relationship (satvapCpp). Used as a
    // physically-motivated lower bound on estimated ground-surface
    // temperature in the grid model's ground heat flux shortcut: a real
    // surface's temperature is not expected to fall much below the
    // dewpoint of the air above it.
    double dewpointCpp(double tc, double rh);

    // Fills missing (NaN) values in a gridded time series by iterative
    // inverse-distance-weighted averaging of each cell's valid neighbours
    // -- used to patch small coastal or otherwise masked data gaps in
    // climate-forcing grids before they are used to drive the model.
    // `data` is a flattened rows x cols x time array and `landMask` a
    // flattened rows x cols grid marking which cells should ever be
    // filled (e.g. land, as opposed to open sea for a land-only climate
    // product). Only cells missing in the very first time slice are
    // filled, on the assumption that the pattern of missing data (e.g. a
    // fixed land-sea boundary) does not itself change over time. Within
    // each timestep, the fill is applied repeatedly so that a cell can be
    // filled using a neighbour that was itself only just filled, letting
    // the fill spread inward from the nearest real data; a cell with no
    // valid data anywhere in its neighbourhood at a given timestep is left
    // unfilled.
    std::vector<double> fillNAIDWCpp(std::vector<double> data, int nrow, int ncol, int ntime,
        const std::vector<double>& landMask);

    // Aerodynamic resistance (s/m) to heat from the scalar roughness length
    // zh to a given height, by Monin-Obukhov similarity at Obukhov length LL
    // and friction velocity uf. The grid uses it over bare ground, where
    // there is no transport column, and for the canopy exchange surface's
    // own segment up to canopy top.
    double rHaToHeightScalarCpp(double height, double d, double zh, double LL, double uf);

    // ------------------------------------------------------------------
    // The transport column: scalar resistance from the soil surface to any
    // height within or above the canopy.
    //
    // Four layers, joined so that the diffusivity is continuous. At the
    // soil a wall layer, in which eddies can be no larger than their
    // distance from the ground, driven by the gust strength that survives
    // to the floor and weakens as the foliage above thickens. Above it the
    // canopy interior, Raupach's cosine gust profile with a canopy-scale
    // Lagrangian time scale, taking over at the height where it mixes more
    // effectively than the wall layer does. Above canopy top a sublayer, to
    // about twice canopy height, across which the canopy's diffusivity is
    // carried up and blended into the surface-layer value. Above that the
    // ordinary surface layer.
    //
    // With no foliage the wall layer is the bare-ground surface layer and
    // the other three terms cancel exactly, so bare ground is this same
    // formula at zero plant area rather than a separate branch.
    //
    // Stability enters through the wall layer, the sublayer and the surface
    // layer, and through the canopy time scale; the shape of the canopy
    // interior is neutral, which is the regime Raupach's constants were
    // measured in.
    // ------------------------------------------------------------------

    // Per-canopy constants of the column. Everything here depends only on
    // canopy height, plant area index and zero-plane displacement, so it is
    // built once and reused across the hours of a timestep loop.
    struct ColumnStruct {
        double h = 0.0;     // canopy height (m)
        double d = 0.0;     // zero-plane displacement (m)
        double z0h = 0.0;   // soil heat roughness length, the column's floor (m)
        double a0 = 0.0;    // gust strength at the floor, relative to friction velocity
        double a2n = 0.0;   // dimensionless Lagrangian time scale at neutral
        double Aw = 0.0;    // mean of the canopy gust profile
        double Bw = 0.0;    // half-amplitude of the canopy gust profile
        double Dw = 0.0;    // Aw^2 - Bw^2, equal to a1*a0
        double rw = 0.0;    // sqrt((Aw+Bw)/(Aw-Bw)), equal to sqrt(a1/a0)
        double Jg = 0.0;    // the canopy integral evaluated at the floor
        double zs = 0.0;    // top of the canopy-top sublayer (m)
        double cEff = 0.0;  // canopy-top diffusivity enhancement over the surface layer
        double IJ = 0.0;    // integral of (J - Jg) from the floor to canopy top, in units of h
    };

    // The gradient function implied by dpsihCpp, phi_h* = 1 - zeta
    // dpsi_h/dzeta. It is derived from the correction rather than taken
    // from a separate similarity expression so that the diffusivity the
    // column integrates and the resistance the rest of the model uses
    // describe the same atmosphere.
    double phihStarCpp(double ze);

    // Per-canopy setup. d is supplied rather than recomputed so that the
    // column and its caller cannot disagree about the displacement height.
    ColumnStruct columnSetupCpp(double h, double pai, double d);

    // Dimensionless Lagrangian time scale at the given Obukhov length.
    double columnA2Cpp(const ColumnStruct& c, double LL);

    // Height at which the canopy interior takes over from the wall layer.
    // The crossing z/phi_h*(z/L) = a0 a1 a2 h/kappa is independent of uf.
    // The neutral solution phi_h* = 1 is used at every stability, while a2
    // retains the current hour's stability dependence. The resistance is
    // stationary in this height, so the approximation is second order.
    double columnHandoverCpp(const ColumnStruct& c, double uf, double LL, double a2);

    // Scalar resistance from the soil surface up to z (s/m).
    double columnResistCpp(const ColumnStruct& c, double uf, double LL, double z);

    // Scalar diffusivity at z (m^2/s). The resistance above is its exact
    // integral, so this is the quantity to check continuity against.
    double columnDiffusivityCpp(const ColumnStruct& c, double uf, double LL, double z);

    // Resistance that turns a canopy source spread uniformly with height into
    // a rise in the air the ground exchanges with: the mean, over canopy height,
    // of the resistance from each height to z. A source of total S (W/m2) then
    // raises that air by S times this, over rho cp. The canopy-interior part of
    // the mean depends only on the canopy and is precomputed in setup; the wall
    // layer's part is closed form, its stability term integrated through
    // psihIntegralCpp. No iteration.
    double columnUniformRiseResistCpp(const ColumnStruct& c, double uf, double LL, double z);

    // Antiderivative of dpsihCpp, so a stability correction can be integrated
    // over height without quadrature:
    //   int_a^b psi_h(z/L) dz = L [psihIntegralCpp(b/L) - psihIntegralCpp(a/L)].
    // Exact except in the saturated tail below zeta = -1.05, which is integrated
    // on four points per unit of zeta; see the definition.
    double psihIntegralCpp(double zeta);

    // Mean, over canopy height, of the resistance from the soil surface to
    // each height: the resistance across which the ground exchanges with air
    // raised by canopy sources spread uniformly with height.
    double columnMeanResistCpp(const ColumnStruct& c, double uf, double LL);

    // Neutral profile-shape ratios at one height:
    //   fR = R(z)/A
    //   fC = r_C(z)/B, with B = A - Rbar = r_C(z0h)
    // Both are rescaled hourly. fC is normalised by the quadrature's own floor
    // value so it is exactly 1 at z0h.
    struct ProfileShape { double fR = 0.0; double fC = 1.0; };
    ProfileShape columnProfileShapeCpp(const ColumnStruct& c, double z);

    // One height of the far-field profile, for temperature or for vapour
    // pressure, from the boundary values, the whole surface's flux and the
    // column's resistances. See the definition for the arguments.
    double farFieldProfileCpp(double sRef, double sGround, double flux, double rc,
        double A, double Rbar, double RzR, double Rz, double rCz);

    // Temperature and relative humidity at one height inside the canopy.
    // `groundhr` is ground-surface relative humidity as a fraction; input `rh`
    // and returned `RHzOut` are percent. `rBL` is the resistance crossed by the
    // big-leaf fluxes.
    void farFieldPairCpp(double Ta, double rh, double pk, double Tc, double Tg,
        double groundhr, double rBL, double rSurf, double hSurf,
        double A, double Rbar, double RzR, double Rz, double rCz,
        double& TzOut, double& RHzOut);

    // Plant functional type whose reference canopy lies closest to a cell's
    // height, plant area index and leaf angle coefficient, as an index into
    // the reference vectors (one entry per type; heights as natural logs).
    // Returns -1 where the height is not positive and finite. See the
    // definition for the distance.
    int pftFromStructureCpp(double h, double pai, double x,
        const std::vector<double>& logRefh, const std::vector<double>& refpai,
        const std::vector<double>& refx);

} // namespace utils
