#include <Rcpp.h>
#include "utils.h"
#include "constants.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

using mc::pi;
using mc::torad;
using mc::ka;
using mc::sb;
using mc::Mw;
using mc::RgasC;

namespace utils {

    double zeroplanedisCpp(double h, double pai)
    {
        if (pai < 0.001) pai = 0.001;
        return (1.0 - (1.0 - std::exp(-std::sqrt(7.5 * pai))) / std::sqrt(7.5 * pai)) * h;
    }

    // Substrate and roughness-element drag coefficients of Raupach (1994)
    // Sect. 4a. Shared with dragClosureCpp below, which measures how far the
    // canopy has closed on the same scale.
    namespace {
        constexpr double dragCS = 0.003;
        constexpr double dragCR = 0.3;
        // Once elements are close enough to shelter one another the derivation
        // ceases to hold, and drag data place the ratio's limit here
        // (Raupach 1994 Sect. 4c). The ceiling engages at a plant area index
        // of 0.58; below that the unbounded expression governs.
        constexpr double dragRatioMax = 0.3;
    }

    // Friction-velocity to canopy-top wind ratio u*/U_h (Raupach 1994).
    // Capped at 0.3 once sheltering dominates, at PAI >= 0.58.
    double dragPartitionCpp(double pai)
    {
        double ratio = std::sqrt(dragCS + (dragCR * pai) / 2.0);
        if (ratio > dragRatioMax) ratio = dragRatioMax;
        return ratio;
    }

    // How far the canopy has closed aerodynamically, on the drag partition's
    // own scale: zero over bare substrate, one once roughness elements shelter
    // one another and the partition saturates. The canopy-top mixing-layer
    // instability that the roughness sublayer rests on requires that closure,
    // so this is the variable the sublayer's strength follows.
    double dragClosureCpp(double pai)
    {
        // The bare-substrate value of the partition, formed from the same
        // expression so that zero plant area returns exactly zero.
        const double bare = std::sqrt(dragCS);
        double f = (dragPartitionCpp(pai) - bare) / (dragRatioMax - bare);
        if (f < 0.0) f = 0.0;
        if (f > 1.0) f = 1.0;
        return f;
    }

    // Resistance effect of the roughness sublayer, as Raupach's (1992) influence
    // function of the sublayer's own depth coefficient. Integrating 1/K_H across a
    // sublayer of depth coefficient c, from canopy top to any height above it,
    // takes Psi(c)/(kappa u*) less than surface-layer similarity would; Psi(1) = 0,
    // so a surface with no sublayer is untouched.
    double sublayerInfluenceCpp(double c)
    {
        if (!(c > 1.0)) return 0.0;
        return std::log(c) - 1.0 + 1.0 / c;
    }

    // Influence function of the momentum roughness sublayer,
    // Psi(1 + (cW - 1) f): zero over bare substrate and Psi(2) = 0.193 once closed.
    // Folded into roughlengthCpp and taken back out at canopy top.
    double momentumInfluenceCpp(double pai)
    {
        return sublayerInfluenceCpp(1.0 + (mc::cW - 1.0) * dragClosureCpp(pai));
    }

    // Top of the momentum roughness sublayer, z_w - d = (1 + (cW - 1) f)(h - d).
    // Canopy top itself over bare substrate, where the sublayer has no depth.
    double momentumSublayerTopCpp(double h, double pai, double d)
    {
        return d + (1.0 + (mc::cW - 1.0) * dragClosureCpp(pai)) * (h - d);
    }

    // Height at which the flow has cleared the canopy's roughness sublayers: the
    // top of the scalar sublayer, d + (1 + (cH - 1) f)(h - d). Both sublayers
    // deepen with closure, the scalar one to cH and the momentum one to cW < cH,
    // so the scalar sublayer's top clears both. Ground level without a canopy.
    double rslClearHeightCpp(double h, double pai)
    {
        if (h <= 0.0) return 0.0;
        double d = zeroplanedisCpp(h, pai);
        return d + (1.0 + (mc::cH - 1.0) * dragClosureCpp(pai)) * (h - d);
    }

    // Roughness length for heat: the momentum length's construction, carrying the
    // scalar-to-momentum offset kB^-1 = ln 5 as the empirical 0.2 factor, and **no**
    // roughness-sublayer term of either kind.
    //
    // It is the lower limit of the canopy exchange surface's own segment, from that
    // surface to canopy top. The sublayer lies above canopy top and is integrated
    // explicitly there, from the transport column (rHaToHeightCpp); folding an
    // influence function in here as well would count it twice. The momentum length
    // keeps its own influence function, which belongs to the momentum sublayer and
    // to the matched pair that takes it back out at canopy top.
    //
    // The bounds mirror roughlengthCpp's, in the same order, applied to this length's
    // own raw value rather than to an already-bounded momentum length.
    double scalarRoughlengthCpp(double h, double pai, double d)
    {
        if (h <= 0.0) return 0.2 * mc::bareGroundZ0;
        double Be = dragPartitionCpp(pai);
        double zh = (h - d) * std::exp(-ka / Be);
        if (zh > (0.9 * (h - d))) zh = 0.9 * (h - d);
        if (zh < mc::bareGroundZ0) zh = mc::bareGroundZ0;
        return 0.2 * zh;
    }

    double roughlengthCpp(double h, double pai, double d)
    {
        if (h <= 0.0) return mc::bareGroundZ0;
        double Be = dragPartitionCpp(pai);
        // The roughness sublayer mixes momentum more vigorously than
        // surface-layer theory allows, so the canopy absorbs momentum more
        // effectively than its geometry alone implies and presents a larger
        // roughness length to the flow above (Raupach 1992 Eq. 27). The
        // sublayer, and so this enhancement, grows with canopy closure.
        double zm = (h - d) * std::exp(-ka / Be + momentumInfluenceCpp(pai));
        if (zm > (0.9 * (h - d))) zm = 0.9 * (h - d);
        // The canopy expression tends to a fixed fraction of canopy height as
        // foliage vanishes, so it describes an ever-smoother surface as the
        // canopy shortens and thins. Vegetation cannot make a surface smoother
        // than the soil beneath it, so bare ground bounds it from below.
        if (zm < mc::bareGroundZ0) zm = mc::bareGroundZ0;
        return zm;
    }

    // Molar density of air (mol/m^3) at temperature tc (deg C) and pressure pk (kPa).
    double phairCpp(double tc, double pk)
    {
        double tk = tc + 273.15;
        double ph = 44.6 * (pk / 101.3) * (273.15 / tk);
        return ph;
    }

    // Molar heat capacity of air at constant pressure (J/mol/K) at tc (deg C).
    double cpairCpp(double tc)
    {
        double cp = 2e-05 * tc * tc + 0.0002 * tc + 29.119;
        return cp;
    }

    // Stable-branch constants of dpsimCpp, dpsim(ze) = -c tanh(ze / zetaMax),
    // shared with the closed-form stable inversion below.
    static const double stableZetaMaxM = 4.0 / 4.7;
    static const double stableScaleM = 4.7 * stableZetaMaxM;

    double dpsimCpp(double ze)
    {
        double psim;
        if (ze < 0.0) {
            double x = std::pow((1.0 - 15.0 * ze), 0.25);
            psim = std::log(std::pow((1.0 + x) / 2.0, 2.0) * (1.0 + x * x) / 2.0) - 2.0 * std::atan(x) + pi / 2.0;
            if (psim > 3.0) psim = 3.0;
        } else {
            psim = -4.7 * (stableZetaMaxM * std::tanh(ze / stableZetaMaxM));
        }
        return psim;
    }

    double dpsihCpp(double ze)
    {
        double psih;
        if (ze < 0.0) {
            double y = std::sqrt(1.0 - 9.0 * ze);
            psih = std::log(std::pow((1.0 + y) / 2.0, 2.0));
            // A hard ceiling would make the implied gradient function jump
            // where it bit, and that gradient sets the height at which the
            // transport column hands over from its wall layer to the canopy
            // interior. The correction therefore saturates smoothly instead,
            // unchanged below the knee and approaching the same asymptote.
            if (psih > mc::psihKnee) {
                psih = mc::psihKnee + mc::psihKnee * std::tanh((psih - mc::psihKnee) / mc::psihKnee);
            }
        } else {
            const double zetaMaxH = 4.0 / (4.7 / 0.74);
            psih = -(4.7 / 0.74) * (zetaMaxH * std::tanh(ze / zetaMaxH));
        }
        return psih;
    }

    double recoverPsiMCpp(double uf, double Ueff, double zref, double d, double zm)
    {
        return (ka * Ueff) / uf - std::log((zref - d) / zm);
    }

    // The bulk stability correction on the stable branch, psi_m(L) =
    // dpsimCpp(zm/L) - dpsimCpp((zref-d)/L), rises from 0 as L -> infinity
    // (neutral) to an interior peak at some finite L, then falls back
    // towards 0 again as L -> 0+ -- both terms share the same asymptotic
    // value in that limit, so their difference vanishes at both ends of
    // the range rather than growing monotonically. Inverting psi_m(L)
    // therefore requires first locating that peak (a coarse log-spaced
    // scan to bracket it, refined by golden-section search), then solving
    // by bisection on the branch beyond the peak, where psi_m(L) is
    // guaranteed to decrease monotonically towards 0 -- always the milder,
    // more conservative (larger-L) of the two roots when both exist. If
    // the requested correction exceeds what the equation can produce at
    // all (at or beyond the peak), the result saturates to the peak's own
    // L rather than failing.
    double lStableFinalCpp(double psi_m, double zref, double d, double zm)
    {
        if (psi_m <= 1e-12) return std::numeric_limits<double>::infinity();

        auto psiOfL = [&](double L) {
            return dpsimCpp(zm / L) - dpsimCpp((zref - d) / L);
        };

        // 1. Coarse log-spaced scan (1e-4 to 1e6 m, covering any
        // physically plausible L) to bracket the peak.
        const int NSCAN = 60;
        double bestL = 1.0, bestVal = psiOfL(1.0);
        double logLo = std::log(1e-4), logHi = std::log(1e6);
        for (int i = 0; i <= NSCAN; ++i) {
            double L = std::exp(logLo + (logHi - logLo) * i / NSCAN);
            double v = psiOfL(L);
            if (v > bestVal) { bestVal = v; bestL = L; }
        }
        int bestIdx = static_cast<int>(std::round(NSCAN * (std::log(bestL) - logLo) / (logHi - logLo)));
        double loL = std::exp(logLo + (logHi - logLo) * std::max(0, bestIdx - 1) / NSCAN);
        double hiL = std::exp(logLo + (logHi - logLo) * std::min(NSCAN, bestIdx + 1) / NSCAN);

        // 2. Golden-section search refines the peak within that bracket.
        const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
        double a = loL, b = hiL;
        double c = b - gr * (b - a), e = a + gr * (b - a);
        for (int i = 0; i < 60; ++i) {
            if (psiOfL(c) > psiOfL(e)) { b = e; } else { a = c; }
            c = b - gr * (b - a); e = a + gr * (b - a);
            if (std::abs(b - a) < 1e-9 * b) break;
        }
        double Lpeak = 0.5 * (a + b);
        double psiPeak = psiOfL(Lpeak);

        // 3. A target psi_m at or beyond what this equation can produce at
        // all (at or above the peak) saturates to Lpeak, the tightest
        // stability this taper can represent, rather than failing.
        if (psi_m >= psiPeak) return Lpeak;

        // 4. Bisection for the unique root on the branch L > Lpeak, where
        // psiOfL is guaranteed monotonically decreasing (from psiPeak down
        // to 0 as L -> infinity) -- the milder, more conservative root.
        double lo = Lpeak, hi = std::max(Lpeak * 100.0, 1e6);
        while (psiOfL(hi) > psi_m && hi < 1e12) hi *= 10.0;
        for (int i = 0; i < 80; ++i) {
            double mid = 0.5 * (lo + hi);
            if (psiOfL(mid) > psi_m) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    double lUnstableNewtonCpp(double psi_m, double zref, double d, double zm, int nNewton)
    {
        // psi_m and its derivative, re-expressed in terms of the
        // substitution x = (1-15*ze)^(1/4) used by dpsimCpp's own unstable
        // branch, so this root can be solved by Newton-Raphson on x
        // directly. Unlike the stable branch above, the unstable
        // correction grows without a shared asymptote, so no equivalent
        // fold-back arises here and a single Newton step is sufficient.
        auto psiOfX = [](double x) {
            return 2.0 * std::log((1.0 + x) / 2.0) + std::log((1.0 + x * x) / 2.0)
                 - 2.0 * std::atan(x) + pi / 2.0;
        };
        auto dpsiDx = [](double x) {
            return 2.0 / (1.0 + x) - 2.0 * (1.0 - x) / (1.0 + x * x);
        };

        double r = zm / (zref - d);
        double target = psi_m;

        auto fFun = [&](double xr) {
            double xm4 = 1.0 - r + r * xr * xr * xr * xr;
            double xm = (xm4 > 0.0) ? std::pow(xm4, 0.25) : 1e-6;
            return psiOfX(xm) - psiOfX(xr) - target;
        };
        auto fPrimeFun = [&](double xr) {
            double xm4 = 1.0 - r + r * xr * xr * xr * xr;
            double xm = (xm4 > 0.0) ? std::pow(xm4, 0.25) : 1e-6;
            double dxm_dxr = (xm > 0.0) ? (r * xr * xr * xr / (xm * xm * xm)) : 0.0;
            return dpsiDx(xm) * dxm_dxr - dpsiDx(xr);
        };

        double xr = 2.0 * std::sqrt(1.0 + std::max(-target, 0.0)) - 1.0;
        for (int i = 0; i < nNewton; ++i) {
            double fx = fFun(xr);
            double fp = fPrimeFun(xr);
            if (fp == 0.0 || !std::isfinite(fp)) break;
            double xrNew = xr - fx / fp;
            if (xrNew <= 1.0) xrNew = (xr + 1.0) / 2.0;
            xr = xrNew;
        }
        double zetaR = (1.0 - xr * xr * xr * xr) / 15.0;
        if (zetaR >= 0.0) return -std::numeric_limits<double>::infinity();
        return (zref - d) / zetaR;
    }

    double recoverLCpp(double uf, double Ueff, double zref, double d, double zm)
    {
        double psi_m = recoverPsiMCpp(uf, Ueff, zref, d, zm);
        if (psi_m >= 0.0) {
            return lStableFinalCpp(psi_m, zref, d, zm);
        }
        return lUnstableNewtonCpp(psi_m, zref, d, zm, 1);
    }

    double stableRecoveryPeakCpp(double zref, double d, double zm)
    {
        const double zr = zref - d;
        // Where the roughness term is negligible the peak satisfies
        // cosh^2(zr u) = zr / zm; Newton on the derivative then includes it.
        double u = std::acosh(std::sqrt(zr / zm)) / zr;
        for (int k = 0; k < 2; ++k) {
            double cr = std::cosh(zr * u), cm = std::cosh(zm * u);
            double d1 = zr / (cr * cr) - zm / (cm * cm);
            double d2 = -2.0 * (zr * zr * std::tanh(zr * u) / (cr * cr) - zm * zm * std::tanh(zm * u) / (cm * cm));
            u -= d1 / d2;
        }
        return u;
    }

    double lStableClosedCpp(double psi_m, double zref, double d, double zm, double uPeak)
    {
        if (psi_m <= 1e-12) return std::numeric_limits<double>::infinity();
        const double zr = zref - d;
        const double a = 1.0 / stableZetaMaxM;
        const double c = stableScaleM;
        const double psiPeak = c * (std::tanh(zr * uPeak) - std::tanh(zm * uPeak));
        if (psi_m >= psiPeak) return a / uPeak;
        double t = psi_m / c;
        if (t > 0.999999) t = 0.999999;
        double u = std::atanh(t) / zr;
        for (int k = 0; k < 2; ++k) {
            double cr = std::cosh(zr * u), cm = std::cosh(zm * u);
            double f = c * (std::tanh(zr * u) - std::tanh(zm * u)) - psi_m;
            double fp = c * (zr / (cr * cr) - zm / (cm * cm));
            u -= f / fp;
            if (u < 1e-12) u = 1e-12;
            if (u > uPeak) u = uPeak;
        }
        return a / u;
    }

    double freeConvectiveVelocityCpp(double H, double tc, double pk)
    {
        if (!(H > 0.0)) return 0.0;
        double Tk = tc + 273.15;
        double cp = cpairCpp(tc);
        double ph = phairCpp(tc, pk);
        return std::cbrt((mc::g / Tk) * mc::freeConvZi * (H / (ph * cp)));
    }

    double drivingWindCpp(double u, double wstar)
    {
        double U = u;
        if (wstar > 0.0) U = std::sqrt(u * u + std::pow(mc::freeConvBeta * wstar, 2.0));
        if (U < mc::minWindSpeed) U = mc::minWindSpeed;
        return U;
    }

    double recoverLClosedCpp(double uf, double Ueff, double zref, double d, double zm, double uPeak)
    {
        double psi_m = recoverPsiMCpp(uf, Ueff, zref, d, zm);
        if (psi_m >= 0.0) {
            return lStableClosedCpp(psi_m, zref, d, zm, uPeak);
        }
        return lUnstableNewtonCpp(psi_m, zref, d, zm, 1);
    }

    double clipMOlength(double L, double zref, double d, double zm, double beta)
    {
        const double ln_z = std::log((zref - d) / zm);
        const double psim_min = -beta * ln_z;
        const double tol = 1e-4;
        double psim = dpsimCpp(zm / L) - dpsimCpp((zref - d) / L);
        if (L < 0.0 && psim < psim_min) {
            double L_low = -500.0;
            double L_high = L;
            double L_mid = L;
            for (int i = 0; i < 30; ++i) {
                L_mid = 0.5 * (L_low + L_high);
                double psim_mid = dpsimCpp(zm / L_mid) - dpsimCpp((zref - d) / L_mid);
                if (psim_mid < psim_min)
                    L_low = L_mid;
                else
                    L_high = L_mid;
                if (std::abs(psim_mid - psim_min) < tol)
                    break;
            }
            return L_high;
        }
        if (L >= 0.0) {
            // The bound is applied by comparing L itself against the L
            // that would produce exactly the maximum allowed correction,
            // rather than comparing the correction psim(L) directly
            // against that maximum. This matters because psim(L) on the
            // stable branch is not monotonic (see lStableFinalCpp): it
            // peaks at some intermediate L and falls back towards zero
            // again as L becomes very small, so an extremely over-stable L
            // can otherwise look numerically indistinguishable from a
            // near-neutral one if the correction itself is compared
            // instead. Comparing L directly against the milder
            // (larger-L) bound catches this over-stable case correctly.
            const double psim_max = beta * ln_z;
            const double L_bound = lStableFinalCpp(psim_max, zref, d, zm);
            if (L < L_bound) {
                return L_bound;
            }
        }
        return L;
    }

    // Beam extinction coefficient of a homogeneous canopy from solar zenith
    // and the Campbell leaf-angle distribution parameter.
    double canopyKCpp(double zenr, double x)
    {
        double k;
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
        return k;
    }

    // Convert canopy beam extinction to the geometry of the actual ground
    // surface: `k` describes interception by foliage, while `kd` and `Kc`
    // rescale that for the longer, more oblique optical path a beam takes
    // over sloping ground compared to a flat reference. Canopy gaps are a
    // separate mechanism, handled downstream by canopyGapMixCpp.
    kstruct cankCpp(double zenr, double x, double si)
    {
        if (zenr > (pi / 2.0)) zenr = pi / 2.0;
        if (si < 0.0) si = 0.0;
        double k = canopyKCpp(zenr, x);
        if (k > 6000.0) k = 6000.0;
        double kd = k * std::cos(zenr) / si;
        if (si == 0) kd = 1.0;
        // A beam grazing the surface travels an unbounded path through the
        // canopy. Bound the optical depth per unit plant area as k itself is
        // bounded, so that the two-stream coefficients stay finite; the beam
        // flux crossing the canopy plane vanishes with si in any case.
        if (kd > 6000.0) kd = 6000.0;
        double Kc = 1.0 / si;
        if (si == 0.0) Kc = 600.0;
        kstruct kparams;
        kparams.k = k;
        kparams.kd = kd;
        kparams.Kc = Kc;
        return kparams;
    }

    // Diffuse two-stream solution for a canopy of plant area index `pait`.
    // Leaf absorption/scattering, leaf-angle distribution and ground
    // reflectance are reduced to coefficients giving upward and downward
    // diffuse flux at arbitrary canopy depth.
    tsdifstruct twostreamdifCpp(double pait, double x, double lref, double ltra, double gref)
    {
        tsdifstruct params;
        // Single-scattering albedo and true absorption by foliage.
        params.om = lref + ltra;
        params.a = 1.0 - params.om;
        params.del = lref - ltra;
        params.J = 1.0 / 3.0;
        if (x != 1.0) {
            double mla = 9.65 * std::pow((3.0 + x), -1.65);
            if (mla > pi / 2.0) mla = pi / 2.0;
            params.J = std::cos(mla) * std::cos(mla);
        }
        // Backscatter and attenuation parameters determine the coupled
        // upward/downward diffuse streams.
        params.gma = 0.5 * (params.om + params.J * params.del);
        params.h = std::sqrt(params.a * params.a + 2.0 * params.a * params.gma);
        params.S1 = std::exp(-params.h * pait);
        params.u1 = params.a + params.gma * (1.0 - 1.0 / gref);
        double u2 = params.a + params.gma * (1.0 - gref);
        params.D1 = (params.a + params.gma + params.h) * (params.u1 - params.h) *
            1.0 / params.S1 - (params.a + params.gma - params.h) * (params.u1 + params.h) * params.S1;
        params.D2 = (u2 + params.h) * 1.0 / params.S1 - (u2 - params.h) * params.S1;
        // Boundary-condition coefficients: p1/p2 describe upward diffuse
        // flux; p3/p4 describe downward diffuse flux.
        params.p1 = (params.gma / (params.D1 * params.S1)) * (params.u1 - params.h);
        params.p2 = (-params.gma * params.S1 / params.D1) * (params.u1 + params.h);
        params.p3 = (1.0 / (params.D2 * params.S1)) * (u2 + params.h);
        params.p4 = (-params.S1 / params.D2) * (u2 - params.h);
        return params;
    }

    // Direct-beam companion to the diffuse two-stream solution. It describes
    // how intercepted beam radiation is scattered upward and downward through
    // the canopy while the unscattered beam attenuates exponentially.
    tsdirstruct twostreamdirCpp(double pait, double om, double a, double gma, double J, double del, double h, double gref,
        double kd, double u1, double S1, double D1, double D2)
    {
        tsdirstruct params;
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
        // Source of the downward beam-scattered stream: the beam scatters
        // forward into it directly, and what it scatters upward is in turn
        // scattered back down by the foliage above, so the two terms add.
        params.p8 = sstr * (a + gma + kd) + gma * ss;
        double v3 = (sstr + gma * gref - (params.p8 / params.sig) * (u2 - kd)) * S2;
        params.p9 = (-1 / D2) * ((params.p8 / (params.sig * S1)) * (u2 + h) + v3);
        params.p10 = (1 / D2) * (((params.p8 * S1) / params.sig) * (u2 - h) + v3);
        return params;
    }

    // Mix radiation travelling through explicit canopy gaps with the
    // homogeneous-canopy two-stream solution.
    double canopyGapMixCpp(double gapfrac, double bypass, double clean)
    {
        return gapfrac * bypass + (1.0 - gapfrac) * clean;
    }

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

    // Apparent solar time (hours). `lt` must be UTC, or use a timezone
    // consistent with `lon`; otherwise the result can be wrong by whole hours.
    double soltimeCpp(int jd, double lt, double lon)
    {
        double m = 6.24004077 + 0.01720197 * (jd - 2451545.0);
        double eot = -7.659 * std::sin(m) + 9.863 * std::sin(2 * m + 3.5932);
        double st = lt + (4.0 * lon + eot) / 60.0;
        return st;
    }

    // Cosine of solar incidence on the local ground surface, including
    // slope/aspect and optional terrain shadowing. This is the direct-beam
    // projection factor used by the radiation model.
    double solarindexCpp(double slope, double aspect, double zend, double azid, bool shadowmask, bool shadow)
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
        if (shadow) si = 0.0;
        return si;
    }

    // Saturation vapour pressure (kPa) at tc (deg C).
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

    // Saturated water-vapour density (kg/m^3) at temperature Tk (Kelvin).
    double satVapDensityCpp(double Tk)
    {
        double tc = Tk - 273.15;
        double es_Pa = 1000.0 * satvapCpp(tc);
        return Mw * es_Pa / (RgasC * Tk);
    }

    double kMeanCpp(const std::string& meanType, double k1, double k2) {
        if (meanType == "GEOMETRIC") {
            return std::sqrt(k1 * k2);
        }
        else if (meanType == "LOGARITHMIC") {
            if (k1 == k2) {
                return k1;
            }
            else {
                return (k1 - k2) / std::log(k1 / k2);
            }
        }
        else {
            throw std::invalid_argument("Unknown mean type. Use 'GEOMETRIC' or 'LOGARITHMIC'.");
        }
    }

    ThomasResult thomasSolveCpp(std::vector<double> aa, std::vector<double> bb,
        std::vector<double> cc, std::vector<double> dd, std::vector<double> x,
        int first, int last)
    {
        for (int i = first; i < last; ++i) {
            cc[i] = cc[i] / bb[i];
            dd[i] = dd[i] / bb[i];
            bb[i + 1] -= aa[i + 1] * cc[i];
            dd[i + 1] -= aa[i + 1] * dd[i];
        }
        x[last] = dd[last] / bb[last];
        for (int i = (last - 1); i >= first; --i) {
            x[i] = dd[i] - cc[i] * x[i + 1];
        }
        ThomasResult out;
        out.bb = bb;
        out.cc = cc;
        out.dd = dd;
        out.x = x;
        return out;
    }

    // Return n + 2 soil-node depths (m), starting at 0 and thickening with depth.
    std::vector<double> geometricCpp(int n, double totalDepth) {
        std::vector<double> z(n + 2);
        double weightSum = 0.0;
        for (int i = 1; i <= n; ++i) {
            weightSum += static_cast<double>(i) * static_cast<double>(i);
        }
        double dz_unit = totalDepth / weightSum;
        z[0] = 0.0;
        z[1] = dz_unit;
        for (int i = 2; i <= n + 1; ++i) {
            z[i] = z[i - 1] + dz_unit * static_cast<double>(i) * static_cast<double>(i);
        }
        return z;
    }

    // Effective soil thermal conductivity from the volume fractions of
    // quartz, other minerals, organic matter, water and air. The pore-space
    // term includes water-vapour heat transport, so conductivity responds to
    // both soil wetness and temperature.
    double thermalConductivityCpp(double Vq, double Vm, double Vo, double Vw,
        double Mc, double Tc, double pk)
    {
        // Conductivity of the solid skeleton, weighted geometrically by its
        // mineral and organic constituents.
        double Vsolid = Vq + Vm + Vo;
        double kSolid = std::pow(8.8, Vq / Vsolid) * std::pow(2.5, Vm / Vsolid) * std::pow(0.25, Vo / Vsolid);
        double Ga = 0.088 * ((Vq + Vm) / Vsolid) + 0.5 * (Vo / Vsolid);
        double Q = 7.25 * Mc + 2.52;
        double xwo = 0.33 * Mc + 0.078;
        // Split pore space between liquid water and gas.
        double porosity = 1.0 - Vsolid;
        double gasPorosity = porosity - Vw;
        if (gasPorosity < 0.0)  gasPorosity = 0.0;
        double Tk = Tc + 273.15;
        double Lv = 45144.0 - 48.0 * Tc;
        double svp = 0.611 * std::exp(17.502 * Tc / (Tc + 240.97));
        double slope = 17.502 * 240.97 * svp / std::pow(240.97 + Tc, 2);
        double Dv = 0.0000212 * (101.3 / pk) * std::pow(Tk / 273.15, 1.75);
        double rhoAir = 44.65 * (pk / 101.3) * (273.15 / Tk);
        double stcor = 1.0 - svp / pk;
        if (stcor < 0.3) stcor = 0.3;
        double kWater = 0.56 + 0.0018 * Tc;
        double wf;
        if (Vw < 0.01 * xwo) {
            wf = 0.0;
        }
        else {
            wf = 1.0 / (1.0 + std::pow(Vw / xwo, -Q));
        }
        // Gas-phase conduction includes latent heat transported by vapour
        // diffusion once enough water is present to sustain that pathway.
        double kGas = 0.0242 + 0.00007 * Tc + wf * Lv * rhoAir * Dv * slope / (pk * stcor);
        double Gc = 1.0 - 2.0 * Ga;
        // Combine the water/gas pore phase and solid skeleton with geometry-
        // dependent weighting to obtain bulk conductivity.
        double kFluid = kGas + (kWater - kGas) * std::pow(Vw / porosity, 2);
        double wA = (2.0 / (1.0 + (kGas / kFluid - 1.0) * Ga) +
            1.0 / (1.0 + (kGas / kFluid - 1.0) * Gc)) / 3.0;
        double wW = (2.0 / (1.0 + (kWater / kFluid - 1.0) * Ga) +
            1.0 / (1.0 + (kWater / kFluid - 1.0) * Gc)) / 3.0;
        double wS = (2.0 / (1.0 + (kSolid / kFluid - 1.0) * Ga) +
            1.0 / (1.0 + (kSolid / kFluid - 1.0) * Gc)) / 3.0;
        double out = (wW * kWater * Vw + wA * kGas * gasPorosity + wS * kSolid * Vsolid)
            / (wW * Vw + wA * gasPorosity + wS * Vsolid);
        return out;
    }

    // Volumetric heat capacity of the soil mixture. Water contributes liquid
    // heat capacity above freezing and ice heat capacity below freezing; the
    // remaining pore space contributes the heat capacity of air.
    double heatCapacityCpp(double Vq, double Vm, double Vo, double Vw, double Tc, double pk)
    {
        double Vsum = Vq + Vm + Vo + Vw;
        double Va = 0.0;
        if (Vsum > 1) {
            Vq = Vq / Vsum;
            Vm = Vm / Vsum;
            Vo = Vo / Vsum;
            Vw = Vw / Vsum;
        }
        else {
            Va = 1.0 - Vsum;
        }
        double CH;
        double CHa = cpairCpp(Tc) * phairCpp(Tc, pk) / 1e6;
        if (Tc >= 0) {
            CH = Vq * 2.13 + Vm * 2.31 + Vo * 2.50 + Vw * 4.18 + Va * CHa;
        }
        else {
            double CHi = 1.93 + 0.0067 * Tc;
            CH = Vq * 2.13 + Vm * 2.31 + Vo * 2.50 + Vw * CHi + Va * CHa;
        }
        return CH * 1e6;
    }

    double rademCpp(double tc) {
        return std::pow(tc + 273.15, 4.0);
    }

    // Latent heat of vaporisation above 0 deg C and sublimation below (J/mol),
    // at surface temperature Tsurface (deg C).
    double latentHeatCpp(double Tsurface)
    {
        double la;
        if (Tsurface >= 0) {
            la = 45068.7 - 42.8428 * Tsurface;
        }
        else {
            la = 51078.69 - 4.338 * Tsurface - 0.06367 * Tsurface * Tsurface;
        }
        return la;
    }

    // One linearised Penman-Monteith surface-temperature update.
    // `rV` is vapour resistance; `hr` scales surface saturation vapour pressure.
    // `eaRef` replaces the rh-derived air vapour pressure when eaRef > -998.
    // `gSlope` (W/m^2/K) is dG/dTs at the current surface temperature.
    double penmanMonteithCpp(double Rabs, double Ta, double pk, double rh, double em,
        double rHa, double rV, double Ts, double G, double hr, double eaRef, double gSlope)
    {
        double Rema = em * sb * rademCpp(Ta);
        double cp = cpairCpp(Ta);
        double ph = phairCpp(Ta, pk);
        double Te = (Ts + Ta) / 2.0;
        double Rer = 4.0 * em * sb * std::pow(Te + 273.15, 3.0);
        double la = latentHeatCpp(Ts);
        // Vapour-pressure deficit is measured between ambient air and the
        // effective surface vapour pressure. `hr < 1` lowers the latter for
        // moisture-limited soil while leaving atmospheric vapour pressure unchanged.
        double ea = satvapCpp(Ta) * (rh / 100.0);
        if (eaRef > -998.0) ea = eaRef;
        double Da = hr * satvapCpp(Ta) - ea;
        double De = hr * (satvapCpp(Te + 0.5) - satvapCpp(Te - 0.5));
        // With the substrate flux linearised about Ts, its slope joins the
        // restoring coefficient and its offset from Ts the forcing.
        double TsNew = Ta + ((Rabs - Rema - ((la * ph) / (pk * rV)) * Da - G + gSlope * (Ts - Ta)) /
            (Rer + ph * (cp / rHa + ((la * De) / (pk * rV))) + gSlope));
        return TsNew;
    }

    double alphaWetCpp(double psiw, double psie)
    {
        if (psiw >= 0) return 0.0;
        if (psiw <= psie) return 1.0;
        return std::abs(psiw) / std::abs(psie);
    }

    double alphaDryCpp(double psiw, double psiDry, double psiWilt) {
        if (psiw <= psiWilt) return 0.0;
        if (psiw >= psiDry)  return 1.0;
        return (psiw - psiWilt) / (psiDry - psiWilt);
    }

    // Root fractions by soil node. `z` contains node depths 0..n plus one depth
    // below the boundary. Each solved node receives roots between halfway points
    // to its neighbours; node n-1 receives all remaining deeper roots and boundary
    // node n receives none.
    std::vector<double> rootDistributeCpp(const std::vector<double>& z, double D50, double D95)
    {
        int n = static_cast<int>(z.size()) - 2;
        std::vector<double> v(n + 1, 0.0);
        double p = std::log(19.0) / std::log(D95 / D50);
        auto cumulative = [&](double D) {
            if (D <= 0.0) return 0.0;
            double x = std::pow(D / D50, p);
            return x / (1.0 + x);
        };
        double above = 0.0;
        for (int i = 0; i < n; ++i) {
            double below = (i < n - 1) ? cumulative(0.5 * (z[i] + z[i + 1])) : 1.0;
            v[i] = below - above;
            above = below;
        }
        return v;
    }

    double aitken1d(double oldv, double newv, Aitken1DState& st)
    {
        double r = newv - oldv;
        if (!st.have_prev) {
            st.r_prev = r;
            st.have_prev = true;
            return newv;
        }
        double dr = r - st.r_prev;
        if (dr != 0.0) {
            st.omega = -st.omega * st.r_prev / dr;
        }
        if (st.omega < 0.05) st.omega = 0.05;
        if (st.omega > 0.9)  st.omega = 0.9;
        st.r_prev = r;
        return oldv + st.omega * r;
    }

    double dewpointCpp(double tc, double rh)
    {
        double ea = satvapCpp(tc) * (rh / 100.0);
        double gamma = std::log(ea / 0.61078);
        double td = (237.3 * gamma) / (17.27 - gamma);
        if (td < 0.0) {
            td = (265.5 * gamma) / (21.875 - gamma);
        }
        return td;
    }

    // Fill missing climate values over land from neighbouring valid cells.
    // Only cells missing in the first time slice are candidates for filling.
    // Orthogonal neighbours receive greater weight than diagonals; successive
    // sweeps propagate values across wider gaps. Non-land cells are never filled.
    std::vector<double> fillNAIDWCpp(std::vector<double> data, int nrow, int ncol, int ntime,
        const std::vector<double>& landMask)
    {
        const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        int plane = nrow * ncol;

        // Identify only land cells that need infilling.
        std::vector<int> fillCells;
        fillCells.reserve(plane / 4);
        for (int j = 0; j < ncol; ++j) {
            for (int i = 0; i < nrow; ++i) {
                int idx2d = i + nrow * j;
                if (!std::isnan(landMask[idx2d]) && std::isnan(data[idx2d])) {
                    fillCells.push_back(idx2d);
                }
            }
        }
        if (fillCells.empty()) return data;

        // Neighbour geometry and inverse-distance weights are time-invariant,
        // so calculate them once and reuse for every climate layer.
        int nFill = static_cast<int>(fillCells.size());
        std::vector<int> neighStart(nFill + 1);
        std::vector<int> neighIdx;
        std::vector<double> neighWeight;
        neighIdx.reserve(nFill * 8);
        neighWeight.reserve(nFill * 8);
        for (int k = 0; k < nFill; ++k) {
            int idx2d = fillCells[k];
            int i = idx2d % nrow;
            int j = idx2d / nrow;
            neighStart[k] = static_cast<int>(neighIdx.size());
            for (int dj = -1; dj <= 1; ++dj) {
                int jj = j + dj;
                if (jj < 0 || jj >= ncol) continue;
                for (int di = -1; di <= 1; ++di) {
                    if (di == 0 && dj == 0) continue;
                    int ii = i + di;
                    if (ii < 0 || ii >= nrow) continue;
                    double w = (di == 0 || dj == 0) ? 1.0 : inv_sqrt2;
                    neighIdx.push_back(ii + nrow * jj);
                    neighWeight.push_back(w);
                }
            }
        }
        neighStart[nFill] = static_cast<int>(neighIdx.size());

        // Successive sweeps allow values to propagate across gaps wider than
        // one cell while always using the nearest values currently available.
        for (int t = 0; t < ntime; ++t) {
            int base = plane * t;
            bool changed = true;
            while (changed) {
                changed = false;
                for (int k = 0; k < nFill; ++k) {
                    int idx3d = base + fillCells[k];
                    if (!std::isnan(data[idx3d])) continue;
                    double wsum = 0.0, vsum = 0.0;
                    for (int m = neighStart[k]; m < neighStart[k + 1]; ++m) {
                        double nv = data[base + neighIdx[m]];
                        if (!std::isnan(nv)) {
                            wsum += neighWeight[m];
                            vsum += neighWeight[m] * nv;
                        }
                    }
                    if (wsum > 0.0) {
                        data[idx3d] = vsum / wsum;
                        changed = true;
                    }
                }
            }
        }
        return data;
    }

    double rHaToHeightScalarCpp(double height, double d, double zh, double LL, double uf)
    {
        double psih = dpsihCpp(zh / LL) - dpsihCpp((height - d) / LL);
        return (std::log((height - d) / zh) + psih) / (mc::ka * uf);
    }


    // ------------------------------------------------------------------
    // Transport column
    // ------------------------------------------------------------------

    double phihStarCpp(double ze)
    {
        if (ze < 0.0) {
            double y = std::sqrt(1.0 - 9.0 * ze);
            // -zeta dpsi/dzeta for the uncorrected unstable branch reduces
            // to (1-y)/y, so the gradient function is 1/y.
            double excess = 1.0 / y - 1.0;
            double raw = std::log(std::pow((1.0 + y) / 2.0, 2.0));
            if (raw > mc::psihKnee) {
                // Beyond the knee the correction saturates, and the gradient
                // it implies is damped by the saturation's own slope.
                double s = 1.0 / std::cosh((raw - mc::psihKnee) / mc::psihKnee);
                excess *= s * s;
            }
            return 1.0 + excess;
        }
        const double S = 4.7 / 0.74;
        const double zetaMaxH = 4.0 / S;
        double s = 1.0 / std::cosh(ze / zetaMaxH);
        return 1.0 + S * ze * s * s;
    }

    // The canopy interior's dimensionless integral, J(x) = int_0^x dx'/g^2
    // with g the cosine gust profile. Closed form; at the canopy top the
    // arctangent reaches its half-turn and the first term vanishes.
    static double columnJCpp(const ColumnStruct& c, double x)
    {
        double root = c.Dw * std::sqrt(c.Dw);
        if (x >= 1.0) return c.Aw / root;
        if (x <= 0.0) return 0.0;
        double th = mc::pi * x;
        double first = c.Bw * std::sin(th) / (c.Dw * (c.Aw - c.Bw * std::cos(th)));
        // tan(theta/2) rather than sin/(1+cos): the same quantity, without
        // the cancellation as theta approaches its half-turn.
        double second = 2.0 * c.Aw * std::atan(c.rw * std::tan(0.5 * th)) / root;
        return (first + second) / mc::pi;
    }

    ColumnStruct columnSetupCpp(double h, double pai, double d)
    {
        ColumnStruct c;
        if (pai < 0.0) pai = 0.0;
        c.h = h;
        c.d = d;
        c.z0h = 0.2 * mc::bareGroundZ0;
        // Canopy-top gustiness reaching the floor, decaying with the drag
        // area that shelters it (Massman and Weil 1999 Eq. 10, with
        // Massman's sheltering). At zero plant area it is the canopy-top
        // value itself, which is what makes the bare-ground limit exact.
        double zeta = 0.25 * pai / (1.0 + 0.4 * pai);
        c.a0 = mc::a1 * std::exp(-4.18 * zeta);
        // Canopy-top mixing is enhanced over the surface-layer value only to the
        // extent the canopy has closed aerodynamically: the coherent canopy-scale
        // eddies responsible need an inflected wind profile at canopy top, which
        // sparse roughness does not produce. The enhancement therefore runs from
        // one over bare substrate to its full value once the drag partition
        // saturates, and carries through to the sublayer's depth below, which is
        // where the surface layer reaches this diffusivity.
        c.cEff = 1.0 + (mc::cH - 1.0) * dragClosureCpp(pai);
        c.a2n = c.cEff * ka * (1.0 - d / h) / (mc::a1 * mc::a1);
        c.Aw = 0.5 * (mc::a1 + c.a0);
        c.Bw = 0.5 * (mc::a1 - c.a0);
        c.Dw = mc::a1 * c.a0;           // = Aw^2 - Bw^2, formed this way to stay exact
        c.rw = std::sqrt(mc::a1 / c.a0);
        c.Jg = columnJCpp(c, c.z0h / h);
        // Mean of the canopy interior's dimensionless integral above the floor,
        // for the uniform-source rise. J is smooth, so composite Simpson is ample.
        {
            const int N = 64;
            double x0 = c.z0h / h;
            if (x0 < 1.0) {
                double hx = (1.0 - x0) / N, sum = 0.0;
                for (int k = 0; k <= N; ++k) {
                    double w = (k == 0 || k == N) ? 1.0 : ((k % 2 == 1) ? 4.0 : 2.0);
                    sum += w * (columnJCpp(c, x0 + k * hx) - c.Jg);
                }
                c.IJ = sum * hx / 3.0;
            }
        }
        // Depth of the canopy-top sublayer, from the neutral time scale so
        // that the geometry of the column does not move with stability.
        c.zs = d + (mc::a1 * mc::a1 * c.a2n * h) / ka;
        if (c.zs < h + 1e-9) c.zs = h + 1e-9;
        return c;
    }

    double columnA2Cpp(const ColumnStruct& c, double LL)
    {
        double ph = 1.0;
        if (LL != 0.0) ph = phihStarCpp((c.h - c.d) / LL);
        if (ph < mc::phiStarMin) ph = mc::phiStarMin;
        if (ph > mc::phiStarMax) ph = mc::phiStarMax;
        return c.a2n / ph;
    }

    double columnHandoverCpp(const ColumnStruct& c, double uf, double LL, double a2)
    {
        // Handover satisfies z / phi_h*(z/L) = a0 a1 a2 h / kappa and is
        // independent of uf. We use the neutral solution phi_h* = 1 at every
        // stability, but a2 retains the current hour's stability dependence.
        // The resistance is stationary in this height, so the approximation is
        // second order (<= 0.005 K ground and 0.026 K canopy in the bundled year).
        double A = c.a0 * mc::a1 * a2 * c.h / ka;
        if (A > c.h) A = c.h;
        if (A < c.z0h) A = c.z0h;
        return A;
    }

    // Diffusivity the column itself produces at canopy top, and the value the
    // sublayer above starts from. Where the canopy interior has taken over
    // below the top, this is the interior's own value and the wall term is
    // negative there, so nothing is added. Where the canopy is too sparse for
    // the hand-over to fall below its own top the wall layer still governs,
    // and the sublayer starts from that smaller value rather than from an
    // enhancement no canopy generated.
    static double columnKTopCpp(const ColumnStruct& c, double uf, double LL, double a2)
    {
        double Kc = mc::a1 * mc::a1 * a2 * c.h * uf;
        double ug = (c.a0 / mc::a1) * uf;
        double Kc0 = c.a0 * c.a0 * a2 * c.h * uf;
        double wall = phihStarCpp(c.h / LL) / (ka * ug * c.h) - 1.0 / Kc0;
        // Returned directly rather than through its own reciprocal, so that a
        // canopy whose interior governs at the top is left bit for bit alone.
        if (wall <= 0.0) return Kc;
        return 1.0 / (1.0 / Kc + wall);
    }

    double columnResistCpp(const ColumnStruct& c, double uf, double LL, double z)
    {
        double a2 = columnA2Cpp(c, LL);
        double ug = (c.a0 / mc::a1) * uf;
        double Kc0 = c.a0 * c.a0 * a2 * c.h * uf;
        double zw = columnHandoverCpp(c, uf, LL, a2);
        double zc = z;
        if (zc < c.z0h) zc = c.z0h;
        if (zc > c.h) zc = c.h;
        // The canopy integral already runs from the floor. The wall term therefore
        // adds only its excess resistance, hence subtraction of (zl - z0h) / Kc0.
        double R = (columnJCpp(c, zc / c.h) - c.Jg) / (a2 * uf);
        double zl = (zc < zw) ? zc : zw;
        if (zl > c.z0h) {
            double psi = dpsihCpp(c.z0h / LL) - dpsihCpp(zl / LL);
            R += (std::log(zl / c.z0h) + psi) / (ka * ug) - (zl - c.z0h) / Kc0;
        }
        if (z <= c.h) return R;
        double Ktop = columnKTopCpp(c, uf, LL, a2);
        double Ks = ka * uf * (c.zs - c.d) / phihStarCpp((c.zs - c.d) / LL);
        double m = (Ks - Ktop) / (c.zs - c.h);
        double zb = (z < c.zs) ? z : c.zs;
        double dz = zb - c.h;
        // Diffusivity rising linearly across the sublayer integrates to a
        // logarithm, except where the two ends coincide and it is uniform.
        if (std::abs(m) * (c.zs - c.h) < 1e-9 * Ktop) R += dz / Ktop;
        else R += std::log(1.0 + m * dz / Ktop) / m;
        if (z <= c.zs) return R;
        double psi = dpsihCpp((c.zs - c.d) / LL) - dpsihCpp((z - c.d) / LL);
        R += (std::log((z - c.d) / (c.zs - c.d)) + psi) / (ka * uf);
        return R;
    }

    double columnDiffusivityCpp(const ColumnStruct& c, double uf, double LL, double z)
    {
        double a2 = columnA2Cpp(c, LL);
        double ug = (c.a0 / mc::a1) * uf;
        double Kc0 = c.a0 * c.a0 * a2 * c.h * uf;
        if (z <= c.h) {
            double gx = c.Aw - c.Bw * std::cos(mc::pi * z / c.h);
            double inv = 1.0 / (a2 * c.h * uf * gx * gx);
            // The wall term is positive below the hand-over height and negative
            // above it, so its own sign says where the wall layer governs and
            // the height itself is not needed.
            double zz = (z > c.z0h) ? z : c.z0h;
            double wall = phihStarCpp(zz / LL) / (ka * ug * zz) - 1.0 / Kc0;
            if (wall > 0.0) inv += wall;
            return 1.0 / inv;
        }
        double Ktop = columnKTopCpp(c, uf, LL, a2);
        if (z <= c.zs) {
            double Ks = ka * uf * (c.zs - c.d) / phihStarCpp((c.zs - c.d) / LL);
            return Ktop + (Ks - Ktop) * (z - c.h) / (c.zs - c.h);
        }
        return ka * uf * (z - c.d) / phihStarCpp((z - c.d) / LL);
    }

    // Antiderivative of the heat stability correction, with Psi(0) = 0, so that
    //     int_a^b psi_h(z/L) dz = L [ Psi(b/L) - Psi(a/L) ].
    // Both branches of dpsihCpp integrate in closed form:
    //     stable    psi_h = -4 tanh(zeta/zm)   ->  -4 zm ln cosh(zeta/zm)
    //     unstable  psi_h = 2 ln((1+y)/2), y = sqrt(1-9 zeta), which under the
    //               substitution zeta = (1-y^2)/9 gives the expression below.
    // Above psihKnee the unstable branch saturates through a tanh, which has no
    // elementary integral. That knee sits at zeta = -1.05, reached only where
    // the Obukhov length is shorter than the hand-over height: a tall canopy in
    // strong instability. Only that tail is integrated numerically, on four
    // points per unit of zeta, so the common case is exact algebra.
    double psihIntegralCpp(double zeta)
    {
        if (zeta >= 0.0) {
            const double zm = 4.0 / (4.7 / 0.74);
            return -4.0 * zm * std::log(std::cosh(zeta / zm));
        }
        auto F = [](double y) {
            double l = std::log(0.5 * (1.0 + y));
            return -(4.0 / 9.0) * (0.5 * y * y * l - 0.5 * (0.5 * y * y - y + std::log(1.0 + y)));
        };
        const double yKnee = 2.0 * std::exp(0.5 * mc::psihKnee) - 1.0;
        const double zKnee = (1.0 - yKnee * yKnee) / 9.0;
        double zUnsat = (zeta >= zKnee) ? zeta : zKnee;
        double val = F(std::sqrt(1.0 - 9.0 * zUnsat)) - F(1.0);
        if (zeta < zKnee) {
            static const double gx[4] = { -0.8611363115940526, -0.3399810435848563,
                                           0.3399810435848563,  0.8611363115940526 };
            static const double gw[4] = { 0.3478548451374538, 0.6521451548625461,
                                          0.6521451548625461, 0.3478548451374538 };
            // In deep instability this tail is long -- it reaches zeta = -8 in a
            // tall canopy at a fifth of a metre per second -- so it is panelled
            // by its own length rather than taken whole. Four points per unit of
            // zeta holds it to the accuracy of the algebra above; the tail is
            // rare enough that the extra points cost nothing on average.
            double len = zKnee - zeta;
            int np = (int)std::ceil(len); if (np < 1) np = 1; if (np > 16) np = 16;
            double half = 0.5 * len / np, acc = 0.0;
            for (int q = 0; q < np; ++q) {
                double mid = zKnee - (2 * q + 1) * half;
                for (int k = 0; k < 4; ++k) acc += gw[k] * dpsihCpp(mid + half * gx[k]);
            }
            val -= acc * half;
        }
        return val;
    }

    double columnMeanResistCpp(const ColumnStruct& c, double uf, double LL)
    {
        double a2 = columnA2Cpp(c, LL);
        double ug = (c.a0 / mc::a1) * uf;
        double Kc0 = c.a0 * c.a0 * a2 * c.h * uf;
        double zw = columnHandoverCpp(c, uf, LL, a2);
        // Height integral of the resistance from the floor, over the canopy.
        // Below the floor the resistance is zero, so the interior integral
        // starts there.
        double total = c.h * c.IJ / (a2 * uf);
        if (zw > c.z0h) {
            double psi0 = dpsihCpp(c.z0h / LL);
            double ipsi = LL * (psihIntegralCpp(zw / LL) - psihIntegralCpp(c.z0h / LL));
            double span = zw - c.z0h;
            double intW = (zw * std::log(zw / c.z0h) - span + psi0 * span - ipsi) / (ka * ug)
                        - span * span / (2.0 * Kc0);
            double Ww = (std::log(zw / c.z0h) + psi0 - dpsihCpp(zw / LL)) / (ka * ug) - span / Kc0;
            total += intW + (c.h - zw) * Ww;
        }
        return total / c.h;
    }

    double columnUniformRiseResistCpp(const ColumnStruct& c, double uf, double LL, double z)
    {
        return columnResistCpp(c, uf, LL, z) - columnMeanResistCpp(c, uf, LL);
    }

    // Vertical shape of the within-canopy far-field profile at one height,
    // taken from the column at neutral. The profile needs two resistances at
    // that height: R(z), the resistance from the soil surface up to it, and
    // r_C(z), the resistance that turns a canopy source spread uniformly with
    // height into the rise at it,
    //
    //     r_C(z) = A - (z/h) R(z) - (1/h) int_z^h R(xi) dxi,   A = R(h).
    //
    // The integral has no closed form, so it is done once per canopy here, at
    // neutral, and kept as the two ratios R(z)/A and r_C(z)/B, where
    // B = A - Rbar = r_C(z0h). Each hour those ratios are rescaled by that
    // hour's own A and B, which are closed forms. Over the bundled year that
    // reproduces the same profile evaluated exactly to about 0.01 K, at a
    // sixteenth of the cost.
    //
    // Simpson in log height, which is where the resistance varies smoothly.
    // At z = z0h the quadrature must return B exactly, since the whole
    // integral is then Rbar; that identity is what the panel count is set by.
    ProfileShape columnProfileShapeCpp(const ColumnStruct& c, double z)
    {
        ProfileShape out;
        const double h = c.h;
        const double uf = 1.0, LL = 1e12;      // neutral; uf cancels from both ratios
        double zlo = z; if (zlo < c.z0h) zlo = c.z0h; if (zlo > h) zlo = h;
        double A = columnResistCpp(c, uf, LL, h);
        double B = A - columnMeanResistCpp(c, uf, LL);
        if (!(A > 0.0) || !(B > 0.0)) { out.fR = 0.0; out.fC = 1.0; return out; }
        double Rz = columnResistCpp(c, uf, LL, zlo);
        // (1/h) int_z^h R dxi, Simpson over an even number of panels in log height.
        // Evaluated twice: once from the requested height and once from the column
        // floor. The ratio below is taken against the floor's own value rather than
        // against the closed-form B, so that the shape is exactly 1 at the floor
        // whatever the quadrature's error, and the profile returns the ground
        // surface there identically.
        const int nPanel = 64;
        auto meanAbove = [&](double zFrom) {
            double u0 = std::log(zFrom), u1 = std::log(h), du = (u1 - u0) / nPanel;
            double acc = 0.0;
            for (int i = 0; i <= nPanel; ++i) {
                double zi = std::exp(u0 + i * du);
                double w = (i == 0 || i == nPanel) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
                acc += w * columnResistCpp(c, uf, LL, zi) * zi;   // R dxi = R xi dln xi
            }
            return acc * du / (3.0 * h);
        };
        double rCz = A - (zlo / h) * Rz - meanAbove(zlo);
        double rC0 = A - meanAbove(c.z0h);
        out.fR = Rz / A;
        out.fC = (rC0 > 0.0) ? rCz / rC0 : 1.0;
        return out;
    }

    // One height of the far-field profile, for any scalar carried by the same
    // transport: temperature with rc = rho cp, vapour pressure with
    // rc = rho lambda / p. `sRef` is the scalar at the reference height and
    // `flux` the whole surface's flux of it; `sGround` is its value at the
    // ground surface. The canopy-top value and the air the ground exchanges
    // with follow from the flux and the column, the ground's own flux is then
    // closed on `sGround`, and the profile carries the two sources across the
    // resistances they each cross.
    double farFieldProfileCpp(double sRef, double sGround, double flux, double rc,
        double A, double Rbar, double RzR, double Rz, double rCz)
    {
        if (!(rc > 0.0) || !(Rbar > 0.0)) return sRef;
        double sTop = sRef + flux * (RzR - A) / rc;
        double sStar = sRef + flux * (RzR - Rbar) / rc;
        double fluxGround = rc * (sGround - sStar) / Rbar;
        return sTop + (fluxGround * (A - Rz) + (flux - fluxGround) * rCz) / rc;
    }

    // One cell-hour of the within-canopy profile, shared by grid and point models.
    // `groundhr` is ground-surface relative humidity as a fraction; input `rh`
    // and returned `RHzOut` are percent. `rBL` is the resistance crossed by the
    // bulk canopy+ground fluxes to the reference height.
    void farFieldPairCpp(double Ta, double rh, double pk, double Tc, double Tg,
        double groundhr, double rBL, double rSurf, double hSurf,
        double A, double Rbar, double RzR, double Rz, double rCz,
        double& TzOut, double& RHzOut)
    {
        double ph = phairCpp(Ta, pk);
        double rhocp = ph * cpairCpp(Ta);
        double rcE = latentHeatCpp(Tc) * ph / pk;

        // The whole surface's two fluxes, on the resistances they cross.
        double HBL = rhocp * (Tc - Ta) / rBL;
        double ea = satvapCpp(Ta) * (rh / 100.0);
        double LBL = rcE * (hSurf * satvapCpp(Tc) - ea) / (rBL + rSurf);

        double Tz = farFieldProfileCpp(Ta, Tg, HBL, rhocp, A, Rbar, RzR, Rz, rCz);
        double eg = groundhr * satvapCpp(Tg);
        double ez = farFieldProfileCpp(ea, eg, LBL, rcE, A, Rbar, RzR, Rz, rCz);

        TzOut = Tz;
        double rhz = 100.0 * ez / satvapCpp(Tz);
        if (rhz < 0.0) rhz = 0.0;
        if (rhz > 100.0) rhz = 100.0;
        RHzOut = rhz;
    }

    // Plant functional type whose reference canopy lies closest to a cell's,
    // by a weighted sum of squared differences in height (0.7), plant area
    // index (0.2) and leaf angle coefficient (0.1). Height dominates, and is
    // compared on a log scale because it spans two orders of magnitude across
    // types, so a raw difference would let tall-canopy noise swamp everything
    // else. Leaf angle is a minor tie-breaker, useful mainly for separating
    // needleleaf from broadleaf trees. Each squared difference is divided by
    // its largest value over the types, which puts the three terms on a
    // common 0-1 scale despite their different units. A missing plant area or
    // leaf angle drops its term. The first type wins a tie.
    int pftFromStructureCpp(double h, double pai, double x,
        const std::vector<double>& logRefh, const std::vector<double>& refpai,
        const std::vector<double>& refx)
    {
        if (!(h > 0.0) || !std::isfinite(h)) return -1;
        const std::size_t n = logRefh.size();
        const bool usePai = !std::isnan(pai);
        const bool useX = !std::isnan(x);
        const double lh = std::log(h);
        double maxh = 0.0, maxp = 0.0, maxx = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            double dh = lh - logRefh[k];
            maxh = std::max(maxh, dh * dh);
            if (usePai) { double dp = pai - refpai[k]; maxp = std::max(maxp, dp * dp); }
            if (useX) { double dx = x - refx[k]; maxx = std::max(maxx, dx * dx); }
        }
        int best = -1;
        double bestScore = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0; k < n; ++k) {
            double dh = lh - logRefh[k];
            double score = 0.7 * (dh * dh / maxh);
            if (usePai) { double dp = pai - refpai[k]; score = score + 0.2 * (dp * dp / maxp); }
            if (useX) { double dx = x - refx[k]; score = score + 0.1 * (dx * dx / maxx); }
            if (score < bestScore) { best = static_cast<int>(k); bestScore = score; }
        }
        return best;
    }

} // namespace utils

// Thin global-namespace forwarders exposing selected utils:: functions to R --
// an Rcpp::export tag only works on a plain global-namespace function, not one
// nested inside namespace utils, so these wrappers exist purely to make the
// underlying implementation callable from R while the implementation itself
// stays organised inside the namespace.
// [[Rcpp::export]]
double satvapCpp(double tc) {
    return utils::satvapCpp(tc);
}

// [[Rcpp::export]]
std::vector<double> geometricCpp(int n, double totalDepth) {
    return utils::geometricCpp(n, totalDepth);
}

// [[Rcpp::export]]
std::vector<double> fillNAIDWCpp(std::vector<double> data, int nrow, int ncol, int ntime,
    std::vector<double> landMask) {
    return utils::fillNAIDWCpp(std::move(data), nrow, ncol, ntime, landMask);
}

// Canopy geometry for R, vectorised over cells: arguments are recycled to the
// longest, and a missing value in any argument gives a missing result.
namespace {
    template <typename F>
    Rcpp::NumericVector vectoriseCanopyCpp(F f, std::initializer_list<Rcpp::NumericVector> args)
    {
        R_xlen_t n = 0;
        for (const auto& a : args) {
            if (a.size() == 0) return Rcpp::NumericVector(0);
            if (a.size() > n) n = a.size();
        }
        Rcpp::NumericVector out(n);
        std::vector<double> x(args.size());
        for (R_xlen_t i = 0; i < n; ++i) {
            bool missing = false;
            std::size_t k = 0;
            for (const auto& a : args) {
                x[k] = a[i % a.size()];
                if (Rcpp::NumericVector::is_na(x[k])) missing = true;
                ++k;
            }
            out[i] = missing ? NA_REAL : f(x);
        }
        return out;
    }
}

// [[Rcpp::export]]
Rcpp::NumericVector zeroplanedisCpp(Rcpp::NumericVector h, Rcpp::NumericVector pai) {
    return vectoriseCanopyCpp([](const std::vector<double>& x) {
        return utils::zeroplanedisCpp(x[0], x[1]); }, { h, pai });
}

// [[Rcpp::export]]
Rcpp::NumericVector roughlengthCpp(Rcpp::NumericVector h, Rcpp::NumericVector pai, Rcpp::NumericVector d) {
    return vectoriseCanopyCpp([](const std::vector<double>& x) {
        return utils::roughlengthCpp(x[0], x[1], x[2]); }, { h, pai, d });
}

// [[Rcpp::export]]
Rcpp::NumericVector rslClearHeightCpp(Rcpp::NumericVector h, Rcpp::NumericVector pai) {
    return vectoriseCanopyCpp([](const std::vector<double>& x) {
        return utils::rslClearHeightCpp(x[0], x[1]); }, { h, pai });
}

// [[Rcpp::export]]
double scalarRoughlengthCpp(double h, double pai, double d) {
    return utils::scalarRoughlengthCpp(h, pai, d);
}

// [[Rcpp::export]]
double sublayerInfluenceCpp(double c) {
    return utils::sublayerInfluenceCpp(c);
}

// [[Rcpp::export]]
double momentumInfluenceCpp(double pai) {
    return utils::momentumInfluenceCpp(pai);
}


// The transport column, reachable from R for validation and diagnostics.
// [[Rcpp::export]]
double columnResistCpp(double h, double pai, double d, double uf, double LL, double z) {
    return utils::columnResistCpp(utils::columnSetupCpp(h, pai, d), uf, LL, z);
}

// [[Rcpp::export]]
double columnDiffusivityCpp(double h, double pai, double d, double uf, double LL, double z) {
    return utils::columnDiffusivityCpp(utils::columnSetupCpp(h, pai, d), uf, LL, z);
}

// [[Rcpp::export]]
double columnUniformRiseResistCpp(double h, double pai, double d, double uf, double LL, double z) {
    return utils::columnUniformRiseResistCpp(utils::columnSetupCpp(h, pai, d), uf, LL, z);
}

// [[Rcpp::export]]
std::vector<double> columnProfileShapeCpp(double h, double pai, double d, double z) {
    utils::ProfileShape sh = utils::columnProfileShapeCpp(utils::columnSetupCpp(h, pai, d), z);
    return std::vector<double>{ sh.fR, sh.fC };
}

// [[Rcpp::export]]
double columnHandoverCpp(double h, double pai, double d, double uf, double LL) {
    utils::ColumnStruct c = utils::columnSetupCpp(h, pai, d);
    return utils::columnHandoverCpp(c, uf, LL, utils::columnA2Cpp(c, LL));
}

// [[Rcpp::export]]
std::vector<double> columnGeometryCpp(double h, double pai, double d, double LL) {
    utils::ColumnStruct c = utils::columnSetupCpp(h, pai, d);
    return { c.z0h, c.a0, c.a2n, utils::columnA2Cpp(c, LL), c.zs, c.cEff };
}

// [[Rcpp::export]]
double phihStarCpp(double ze) { return utils::phihStarCpp(ze); }

// [[Rcpp::export]]
double dpsihTestCpp(double ze) { return utils::dpsihCpp(ze); }


// Vegetation fields for R. Each field arrives as a matrix whose rows are cells
// and whose columns are seasonal layers, and is walked once. A field with one
// layer is reused for every layer of a seasonal one. A field that does not
// change is not copied: NULL is returned in its place.

// A plant of zero height or zero plant area has no mass, so a cell where either
// is zero is bare ground and both are set to zero. A seasonal field is zeroed
// in the layers where the cell is bare; a single-layer field paired with a
// seasonal one only where the cell is bare in every layer, so a cell leafless
// for part of the year keeps its height. Fields with different numbers of
// seasonal layers are left unchanged. nBare is the number of cells altered;
// with countOnly, that count alone is produced.
// [[Rcpp::export]]
Rcpp::List zeroBareVegCpp(Rcpp::NumericMatrix hgt, Rcpp::NumericMatrix pai, bool countOnly = false) {
    const R_xlen_t n = hgt.nrow();
    const int nh = hgt.ncol(), np = pai.ncol();
    if (pai.nrow() != n) Rcpp::stop("hgt and pai must have the same number of cells");
    Rcpp::NumericMatrix h = hgt, p = pai;
    bool hChanged = false, pChanged = false;
    int nBare = 0;
    if (nh == np || nh == 1 || np == 1) {
        const int nl = std::max(nh, np);
        std::vector<char> bare(nl);
        for (R_xlen_t i = 0; i < n; ++i) {
            bool anyBare = false, allBare = true;
            for (int k = 0; k < nl; ++k) {
                bare[k] = hgt(i, k % nh) <= 0.0 || pai(i, k % np) <= 0.0;
                anyBare = anyBare || bare[k];
                allBare = allBare && bare[k];
            }
            if (!anyBare) continue;
            bool changed = false;
            for (int k = 0; k < nh; ++k) {
                double v = hgt(i, k);
                if (!(nh == nl ? bare[k] : allBare) || std::isnan(v) || v == 0.0) continue;
                changed = true;
                if (countOnly) continue;
                if (!hChanged) { h = Rcpp::clone(hgt); hChanged = true; }
                h(i, k) = 0.0;
            }
            for (int k = 0; k < np; ++k) {
                double v = pai(i, k);
                if (!(np == nl ? bare[k] : allBare) || std::isnan(v) || v == 0.0) continue;
                changed = true;
                if (countOnly) continue;
                if (!pChanged) { p = Rcpp::clone(pai); pChanged = true; }
                p(i, k) = 0.0;
            }
            if (changed) ++nBare;
        }
    }
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("hgt") = R_NilValue,
        Rcpp::Named("pai") = R_NilValue, Rcpp::Named("nBare") = nBare);
    if (hChanged) out["hgt"] = h;
    if (pChanged) out["pai"] = p;
    return out;
}

// The plant functional type each cell's canopy structure matches
// (utils::pftFromStructureCpp), for cells with a height and plant area in
// every layer and not marked in skip; x, the leaf angle coefficient, may have
// no columns. The match is made on the middle layer, or on the layer of
// largest plant area where the middle one is leafless. Returned per cell: pft,
// the 1-based index of the matched type, 0 for a cell bare in every layer and
// NA for a cell not matched; hgt, the height on the layer matched; and
// evergreen, whether a seasonal plant area never falls below half its maximum.
// [[Rcpp::export]]
Rcpp::List matchStructureCpp(Rcpp::NumericMatrix hgt, Rcpp::NumericMatrix pai, Rcpp::NumericMatrix x,
    Rcpp::LogicalVector skip, Rcpp::NumericVector refh, Rcpp::NumericVector refpai,
    Rcpp::NumericVector refx) {
    const R_xlen_t n = hgt.nrow();
    const int nh = hgt.ncol(), np = pai.ncol(), nx = x.ncol();
    const int nl = std::max(nh, np);
    if (pai.nrow() != n) Rcpp::stop("hgt and pai must have the same number of cells");
    const bool useSkip = skip.size() == n;
    std::vector<double> logRefh(refh.size());
    for (R_xlen_t k = 0; k < refh.size(); ++k) logRefh[k] = std::log(refh[k]);
    const std::vector<double> refp(refpai.begin(), refpai.end()), refxv(refx.begin(), refx.end());

    Rcpp::IntegerVector pft(n, NA_INTEGER);
    Rcpp::NumericVector hOut(n, NA_REAL);
    Rcpp::LogicalVector evergreen(n, NA_LOGICAL);
    for (R_xlen_t i = 0; i < n; ++i) {
        if (useSkip && skip[i] == TRUE) continue;
        bool known = true;
        for (int k = 0; k < nh && known; ++k) known = !std::isnan(hgt(i, k));
        for (int k = 0; k < np && known; ++k) known = !std::isnan(pai(i, k));
        if (!known) continue;
        // Middle layer, or the first layer of largest plant area if it is leafless.
        const int mid = (nl + 1) / 2 - 1;
        bool anyLeafy = false, midLeafy = false;
        int kmax = 0;
        for (int k = 0; k < nl; ++k) {
            bool leafy = hgt(i, k % nh) > 0.0 && pai(i, k % np) > 0.0;
            anyLeafy = anyLeafy || leafy;
            if (k == mid) midLeafy = leafy;
            if (pai(i, k % np) > pai(i, kmax % np)) kmax = k;
        }
        if (!anyLeafy) { pft[i] = 0; continue; }
        const int j = midLeafy ? mid : kmax;
        const double h = hgt(i, j % nh);
        int type = utils::pftFromStructureCpp(h, pai(i, j % np), nx > 0 ? x(i, j % nx) : NA_REAL,
            logRefh, refp, refxv);
        if (type < 0) Rcpp::stop("A cell's layer of largest plant area has no canopy height, so no plant functional type can be matched to it");
        bool ever = false;
        if (np > 1) {
            double pmin = pai(i, 0), pmax = pai(i, 0);
            for (int k = 1; k < np; ++k) {
                pmin = std::min(pmin, pai(i, k));
                pmax = std::max(pmax, pai(i, k));
            }
            ever = pmin >= 0.5 * pmax;
        }
        pft[i] = type + 1;
        hOut[i] = h;
        evergreen[i] = ever;
    }
    return Rcpp::List::create(Rcpp::Named("pft") = pft, Rcpp::Named("hgt") = hOut,
        Rcpp::Named("evergreen") = evergreen);
}

// Fills the missing entries of v, in cells whose type is given, with that
// type's value: type is a 1-based index into byType, and a cell of type NA or
// zero is left alone. Returns the filled matrix, NULL if nothing was missing,
// and n, the number of cells filled.
// [[Rcpp::export]]
Rcpp::List fillByTypeCpp(Rcpp::NumericMatrix v, Rcpp::IntegerVector type, Rcpp::NumericVector byType) {
    const R_xlen_t n = v.nrow();
    const int nc = v.ncol();
    if (type.size() != n) Rcpp::stop("type must have one entry per cell");
    Rcpp::NumericMatrix out = v;
    bool filled = false;
    int nFilled = 0;
    for (R_xlen_t i = 0; i < n; ++i) {
        const int t = type[i];
        if (t == NA_INTEGER || t < 1) continue;
        if (t > byType.size()) Rcpp::stop("type index outside byType");
        bool any = false;
        for (int k = 0; k < nc; ++k) {
            if (!std::isnan(v(i, k))) continue;
            if (!filled) {
                out = Rcpp::clone(v);
                for (double& e : out) if (std::isnan(e)) e = NA_REAL;
                filled = true;
            }
            out(i, k) = byType[t - 1];
            any = true;
        }
        if (any) ++nFilled;
    }
    Rcpp::List res = Rcpp::List::create(Rcpp::Named("values") = R_NilValue,
        Rcpp::Named("n") = nFilled);
    if (filled) res["values"] = out;
    return res;
}

// Whether each cell holds a positive value in any layer; missing values count
// as not positive.
// [[Rcpp::export]]
Rcpp::LogicalVector rowAnyPositiveCpp(Rcpp::NumericMatrix v) {
    const R_xlen_t n = v.nrow();
    Rcpp::LogicalVector out(n, FALSE);
    for (int k = 0; k < v.ncol(); ++k) {
        for (R_xlen_t i = 0; i < n; ++i) if (v(i, k) > 0.0) out[i] = TRUE;
    }
    return out;
}

// The highest roughness-sublayer clear height (utils::rslClearHeightCpp) over
// every cell and every seasonal layer, height and plant area paired cell by
// cell. Cells missing either are passed over; -Inf if there are no others.
// [[Rcpp::export]]
double maxClearHeightCpp(Rcpp::NumericMatrix hgt, Rcpp::NumericMatrix pai) {
    const R_xlen_t n = hgt.nrow();
    const int nh = hgt.ncol(), np = pai.ncol();
    if (pai.nrow() != n) Rcpp::stop("hgt and pai must have the same number of cells");
    double best = R_NegInf;
    for (int k = 0; k < std::max(nh, np); ++k) {
        const int kh = std::min(k, nh - 1), kp = std::min(k, np - 1);
        for (R_xlen_t i = 0; i < n; ++i) {
            const double h = hgt(i, kh), p = pai(i, kp);
            if (std::isnan(h) || std::isnan(p)) continue;
            const double z = utils::rslClearHeightCpp(h, p);
            if (std::isfinite(z) && z > best) best = z;
        }
    }
    return best;
}
