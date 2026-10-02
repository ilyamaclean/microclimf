// constants.h
// Physical, mathematical and model-calibration constants shared across the
// point and grid models.
#pragma once

namespace mc {
    constexpr double pi = 3.14159265358979323846;
    constexpr double torad = pi / 180.0;          // degrees -> radians
    constexpr double sb = 5.67e-8;                // Stefan-Boltzmann constant (W m^-2 K^-4): radiated
                                                   // flux is proportional to absolute temperature^4
    constexpr double surfaceEmissivity = 0.97;    // thermal emissivity, applied identically to canopy
                                                   // and ground (no separate vegetation/soil value)
    constexpr double ka = 0.41;                   // von Karman constant for turbulent surface-layer profiles
    constexpr double g = 9.80665;                 // gravitational acceleration (m / s^2)
    constexpr double omdy = (2.0 * pi) / (24.0 * 3600.0); // angular frequency of the 24-h soil-temperature cycle (s^-1)
    constexpr double Mw = 0.018015;               // molar mass of water (kg mol^-1), used to convert
                                                   // between vapour pressure and vapour density/mixing
                                                   // ratio in humidity calculations
    constexpr double RgasC = 8.314;               // universal gas constant (J / mol / K)
    constexpr double a1 = 1.25;                   // Massman & Weil (1999) constant governing how fast wind
                                                   // speed decays with depth inside the canopy; one
                                                   // literature value used for every canopy type, not a
                                                   // per-plant-functional-type parameter
    constexpr double cH = 3.0;                    // enhancement of scalar diffusivity over its surface-layer
                                                   // value at canopy top, at full canopy closure
                                                   // (Raupach, Finnigan and Brunet 1996, p. 359)
    constexpr double phiStarMin = 0.5;            // bounds on the gradient function that makes the canopy
    constexpr double phiStarMax = 1.5;            // time scale stability-dependent, keeping the canopy
                                                   // interior within the range the closure was calibrated in
    constexpr double psihKnee = 1.5;              // the heat correction saturates smoothly beyond this value
                                                   // instead of meeting a hard ceiling, so that the gradient
                                                   // function it implies stays continuous
    constexpr double cW = 2.0;                    // depth coefficient of the momentum roughness sublayer
                                                   // at full canopy closure, z_w - d = cW (h - d)
                                                   // (Raupach 1992; Verhoef et al. 1997 Eq. 7)
    constexpr double freeConvZi = 1000.0;         // boundary-layer depth scaling the free-convective velocity (m)
    constexpr double freeConvBeta = 1.0;          // weight of the free-convective velocity in the driving wind
    constexpr double minWindSpeed = 0.5;          // floor on the wind driving the surface layer (m/s).
                                                   // Below it the shear-based resistance formulations
                                                   // become singular; both models apply it
    constexpr double dampWindSpeed = 1.0;         // below this reference-height wind (m/s) the point
                                                   // model's hourly iteration damps stability and the
                                                   // ground heat flux with fixed half-way steps
    constexpr double bareGroundZ0 = 0.004;        // aerodynamic roughness length for bare ground (m) --
                                                   // used because the canopy-height-based roughness
                                                   // formula is undefined with no canopy present
    constexpr double psiOvenDry = -1.0e6;         // oven-dry soil water potential (J/kg), the lower
                                                   // bound of the soil water solve. Pore-air humidity
                                                   // there is 0.0006, far drier than any atmosphere
    constexpr double psiFieldCapacity = -33.0;    // field capacity (J/kg): water held against gravity,
                                                   // the wet end of the water available to roots
    constexpr double psiWiltingPoint = -1500.0;   // permanent wilting point (J/kg): the dry end of
                                                   // the water available to roots
    constexpr double psiUptakeLimit = -10000.0;   // water potential (J/kg) at which a layer stops
                                                   // supplying root uptake; low enough that plants
                                                   // have all but stopped drawing water before it
    constexpr double rStomClosed = 2000.0;        // leaf resistance with stomata closed (s/m), a
                                                   // conventional cuticular value; no leaf's stomatal
                                                   // resistance exceeds it, in the dark or in dim light
}
