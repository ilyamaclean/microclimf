// pointmodel.cpp
// Point-scale microclimate model. For each timestep it couples canopy and
// ground radiation, turbulent exchange, plant water use, and layered soil
// heat/water dynamics to solve the surface temperatures and fluxes experienced
// by organisms at a representative location. See pointmodel.h for the model
// state variables and scientific definitions used by each section.
#include "pointmodel.h"
#include "constants.h"
#include "utils.h"
#include <Rcpp.h>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <numeric>

using mc::pi;
using mc::torad;
using mc::ka;
using mc::g;
using mc::sb;
using mc::surfaceEmissivity;
using mc::a1;
using mc::Mw;
using mc::RgasC;

namespace pointmodel {

    // ========================================================================
    // Radiation: two-stream shortwave absorption by canopy and ground
    // ========================================================================

    // Solar zenith and azimuth angle for a given location, date and local time.
    solmodel solpositionCpp(double lat, double lon, int year, int month, int day, double lt)
    {
        int jd = utils::juldayCpp(year, month, day);
        double st = utils::soltimeCpp(jd, lt, lon);
        // Solar zenith (degrees)
        double latr = lat * pi / 180.0;
        double tt = 0.261799 * (st - 12);
        double dec = (pi * 23.5 / 180) * std::cos(2 * pi * ((jd - 159.5) / 365.25));
        double coh = std::sin(dec) * std::sin(latr) + std::cos(dec) * std::cos(latr) * std::cos(tt);
        double z = std::acos(coh) * (180 / pi);
        // Solar azimuth (degrees)
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
        solmodel solpos;
        solpos.zend = z;
        solpos.zenr = z * torad;
        solpos.azid = azi;
        solpos.azir = azi * torad;
        return solpos;
    }

    // Solar zenith angle only -- see pointmodel.h for why this point model
    // skips the azimuth calculation solpositionCpp performs above.
    solmodel solzenithCpp(double lat, double lon, int year, int month, int day, double lt)
    {
        int jd = utils::juldayCpp(year, month, day);
        double st = utils::soltimeCpp(jd, lt, lon);
        double latr = lat * pi / 180.0;
        double tt = 0.261799 * (st - 12);
        double dec = (pi * 23.5 / 180) * std::cos(2 * pi * ((jd - 159.5) / 365.25));
        double coh = std::sin(dec) * std::sin(latr) + std::cos(dec) * std::cos(latr) * std::cos(tt);
        double z = std::acos(coh) * (180 / pi);
        solmodel solpos;
        solpos.zend = z;
        solpos.zenr = z * torad;
        solpos.azid = 0.0;
        solpos.azir = 0.0;
        return solpos;
    }

    // Generic Sellers (1985) two-stream radiation functions are shared with
    // the grid model through utils.cpp; the functions below supply the canopy
    // structure, optical properties and timestep-specific solar forcing.

    // Precompute the canopy optical state that is fixed for a patch: the
    // clumping-adjusted optical depth, diffuse two-stream solution, effective
    // albedo and diffuse transmission to the ground. Timestep-specific solar
    // angle and direct-beam terms are handled separately in RadswabsStepCpp.
    radsetup RadswabsSetupCpp(const vegpstruct& vegp, double gref, double grefPAR)
    {
        radsetup s;
        s.vegetated = (vegp.pai > 0.0);
        if (s.vegetated) {
            s.pait = vegp.pai;
            if (vegp.clump > 0.0) s.pait = vegp.pai / (1 - vegp.clump);
            s.tsd = utils::twostreamdifCpp(s.pait, vegp.x, vegp.lref, vegp.ltra, gref);
            s.tsdPAR = utils::twostreamdifCpp(s.pait, vegp.x, vegp.lrefp, vegp.ltrap, grefPAR);
            s.trd = vegp.clump * vegp.clump;
            s.amx = gref;
            if (s.amx < vegp.lref) s.amx = vegp.lref;
            s.amxPAR = grefPAR;
            if (s.amxPAR < vegp.lrefp) s.amxPAR = vegp.lrefp;
            // Diffuse albedo: a round-trip-through-gaps share (s.trd*s.trd)
            // reflects straight off the ground (bypass = gref); the rest
            // sees the clean two-stream diffuse albedo (tsd.p1+tsd.p2).
            s.albd = utils::canopyGapMixCpp(s.trd * s.trd, gref, s.tsd.p1 + s.tsd.p2);
            if (s.albd > s.amx) s.albd = s.amx;
            if (s.albd < 0.01) s.albd = 0.01;
            s.albdPAR = utils::canopyGapMixCpp(s.trd * s.trd, grefPAR, s.tsdPAR.p1 + s.tsdPAR.p2);
            if (s.albdPAR > s.amxPAR) s.albdPAR = s.amxPAR;
            if (s.albdPAR < 0.01) s.albdPAR = 0.01;
            // Diffuse radiation reaches the ground either directly through
            // canopy gaps or after attenuation/scattering by foliage. The
            // resulting transmission fraction is constrained to [0, 1].
            s.groundRdd = utils::canopyGapMixCpp(s.trd, 1.0,
                s.tsd.p3 * std::exp(-s.tsd.h * s.pait) + s.tsd.p4 * std::exp(s.tsd.h * s.pait));
            if (s.groundRdd > 1.0) s.groundRdd = 1.0;
            if (s.groundRdd < 0.0) s.groundRdd = 0.0;
            s.groundRddPAR = utils::canopyGapMixCpp(s.trd, 1.0,
                s.tsdPAR.p3 * std::exp(-s.tsdPAR.h * s.pait) + s.tsdPAR.p4 * std::exp(s.tsdPAR.h * s.pait));
            if (s.groundRddPAR > 1.0) s.groundRddPAR = 1.0;
            if (s.groundRddPAR < 0.0) s.groundRddPAR = 0.0;
        }
        else {
            s.pait = 0.0;
            s.trd = 0.0;
            s.amx = gref;
            s.albd = gref;
            s.groundRdd = 0.0;
            s.amxPAR = grefPAR;
            s.albdPAR = grefPAR;
            s.groundRddPAR = 0.0;
        }
        return s;
    }

    // Solve shortwave radiation for one timestep. Incoming radiation is split
    // into direct and diffuse components, adjusted for slope/aspect and sky
    // exposure, then propagated through the clumped canopy with the two-stream
    // solution. The function returns absorption by the ground alone and by the
    // combined canopy+ground surface, together with the same quantities for PAR
    // used by the photosynthesis/stomatal model.
    void RadswabsStepCpp(const radsetup& s, const vegpstruct& vegp, double gref, double grefPAR,
        double slope, double aspect, double svfa, const solmodel& solp, envstruct& env)
    {
        double Rsw = env.Rsw;
        double Rdif = env.Rdif;
        if (s.vegetated) {
            if (Rsw > 0.0) {
                double zenr = solp.zenr;
                if (zenr > pi / 2.0) zenr = pi / 2.0;
                double cosz = std::cos(zenr);
                // Direct-beam incidence on the local surface.
                double si = utils::solarindexCpp(slope, aspect, solp.zend, solp.azid);
                utils::kstruct kp = utils::cankCpp(zenr, vegp.x, si);
                utils::tsdirstruct tspdir = utils::twostreamdirCpp(s.pait, s.tsd.om, s.tsd.a, s.tsd.gma, s.tsd.J, s.tsd.del, s.tsd.h, gref,
                    kp.kd, s.tsd.u1, s.tsd.S1, s.tsd.D1, s.tsd.D2);
                utils::tsdirstruct tspdirPAR = utils::twostreamdirCpp(s.pait, s.tsdPAR.om, s.tsdPAR.a, s.tsdPAR.gma, s.tsdPAR.J, s.tsdPAR.del, s.tsdPAR.h, grefPAR,
                    kp.kd, s.tsdPAR.u1, s.tsdPAR.S1, s.tsdPAR.D1, s.tsdPAR.D2);
                // Convert horizontal direct radiation to beam-normal
                // irradiance.
                double Rbeam = (Rsw - Rdif) / cosz;
                if (Rbeam > 1352.0) Rbeam = 1352.0;
                double trb = std::pow(vegp.clump, kp.Kc);
                if (trb > 0.999) trb = 0.999;
                if (trb < 0.0) trb = 0.0;
                // Beam flux crossing the canopy plane. The two-stream solution
                // treats the canopy as slabs parallel to the ground surface --
                // that is what kd = k*cos(zenith)/si describes -- so every beam
                // fraction it returns is a fraction of this flux, not of the
                // beam-normal irradiance. On level ground si equals cos(zenith)
                // and the two coincide.
                double Rb = Rbeam * si;
                // Direct radiation reaches the ground through explicit canopy
                // gaps plus the attenuated beam passing through foliage.
                double trg = utils::canopyGapMixCpp(trb, 1.0, std::exp(-kp.kd * s.pait));
                // Direct-beam-viewed albedo: a combined diffuse-round-trip
                // and direct-beam gap share (s.trd*trb) reflects straight
                // off the ground (bypass = gref); the rest sees the clean
                // two-stream direct-beam albedo.
                double albb = utils::canopyGapMixCpp(s.trd * trb, gref, tspdir.p5 / -tspdir.sig + tspdir.p6 + tspdir.p7);
                if (albb > s.amx) albb = s.amx;
                if (albb < 0.01) albb = 0.01;
                double albbPAR = utils::canopyGapMixCpp(s.trd * trb, grefPAR, tspdirPAR.p5 / -tspdirPAR.sig + tspdirPAR.p6 + tspdirPAR.p7);
                if (albbPAR > s.amxPAR) albbPAR = s.amxPAR;
                if (albbPAR < 0.01) albbPAR = 0.01;
                // Direct-beam transmission to the ground via canopy
                // scattering: only the vegetated fraction of the ground can
                // produce this (a photon travelling straight through an
                // actual physical gap never hits a leaf to scatter off), so
                // bypass = 0.0 here -- the gap fraction's own direct-beam
                // contribution to the ground is credited separately below,
                // in full, via trg.
                double groundRbdd = utils::canopyGapMixCpp(trb, 0.0,
                    (tspdir.p8 / tspdir.sig) * std::exp(-kp.kd * s.pait) +
                    tspdir.p9 * std::exp(-s.tsd.h * s.pait) + tspdir.p10 * std::exp(s.tsd.h * s.pait));
                if (groundRbdd > s.amx) groundRbdd = s.amx;
                if (groundRbdd < 0.0) groundRbdd = 0.0;
                double groundRbddPAR = utils::canopyGapMixCpp(trb, 0.0,
                    (tspdirPAR.p8 / tspdirPAR.sig) * std::exp(-kp.kd * s.pait) +
                    tspdirPAR.p9 * std::exp(-s.tsdPAR.h * s.pait) + tspdirPAR.p10 * std::exp(s.tsdPAR.h * s.pait));
                if (groundRbddPAR > s.amxPAR) groundRbddPAR = s.amxPAR;
                if (groundRbddPAR < 0.0) groundRbddPAR = 0.0;
                env.Rcanopy = (1.0 - s.albd) * Rdif * svfa + (1.0 - albb) * Rb;
                env.RcanopyPAR = (1.0 - s.albdPAR) * Rdif * svfa + (1.0 - albbPAR) * Rb;
                double Rgdif = s.groundRdd * Rdif * svfa + groundRbdd * Rb;
                double RgdifPAR = s.groundRddPAR * Rdif * svfa + groundRbddPAR * Rb;
                // Direct-beam radiation reaching the ground still travelling
                // as a direct beam (as opposed to the canopy-scattered
                // contribution above) reuses trg (already computed above)
                // -- the gap fraction's own bypass = 1.0 share of the
                // direct beam. groundRbdd above already credits the gap
                // fraction's contribution to the canopy-scattered term with
                // bypass = 0.0, so between the two terms the gap fraction's
                // direct beam is counted exactly once, not twice.
                env.Rground = (1.0 - gref) * (Rgdif + trg * Rb);
                env.RgroundPAR = (1.0 - grefPAR) * (RgdifPAR + trg * Rb);
                // Effective albedo is diagnosed against the radiation actually
                // incident on this surface: diffuse irradiance is reduced by
                // sky-view factor, and the beam contributes the flux it
                // delivers across the ground surface.
                env.albedo = 1.0 - (env.Rcanopy / (Rdif * svfa + Rb));
                if (env.albedo > s.amx) env.albedo = s.amx;
                if (env.albedo < 0.01) env.albedo = 0.01;
            }
            else {
                env.Rground = 0.0;
                env.Rcanopy = 0.0;
                env.RgroundPAR = 0.0;
                env.RcanopyPAR = 0.0;
                env.albedo = vegp.lref;
            }
        }
        else {
            env.albedo = gref;
            if (Rsw > 0.0) {
                double zenr = solp.zenr;
                if (zenr > pi / 2.0) zenr = pi / 2.0;
                double cosz = std::cos(zenr);
                double si = utils::solarindexCpp(slope, aspect, solp.zend, solp.azid);
                // dirr is the beam-normal irradiance (the horizontal beam
                // divided by the zenith cosine); si*dirr is that beam's
                // incidence on the (possibly tilted) surface.
                double dirr = (Rsw - Rdif) / cosz;
                // As in the vegetated branch: near the horizon the division by
                // cos(zenith) amplifies the recorded beam without bound, which
                // a surface tilted toward the sun would otherwise absorb.
                if (dirr > 1352.0) dirr = 1352.0;
                env.Rground = (1 - gref) * (Rdif * svfa + si * dirr);
                env.Rcanopy = env.Rground;
                env.RgroundPAR = (1 - grefPAR) * (Rdif * svfa + si * dirr);
                env.RcanopyPAR = env.RgroundPAR;
            }
            else {
                env.Rground = 0.0;
                env.Rcanopy = 0.0;
                env.RgroundPAR = 0.0;
                env.RcanopyPAR = 0.0;
            }
        }
    }

    // Add longwave exchange to the shortwave solution.
    //
    // Each surface sees three things and exchanges with all of them: the sky,
    // the surrounding terrain that blocks the rest of the sky, and - for the
    // ground beneath vegetation - the canopy. The ground's hemisphere is
    // canopy over `1 - trd`, open sky over `trd*svfa` and terrain over
    // `trd*(1 - svfa)`; the whole canopy+ground surface sees sky over `svfa`
    // and terrain over the rest.
    //
    // Every surface emits `em*sigma*T^4` into that hemisphere and absorbs `em`
    // of what arrives. Terrain and canopy return `1 - em` of what reaches
    // them, so the share of a surface's own emission that comes straight back
    // never leaves it, and its net loss to those two carries `em^2` rather
    // than `em`. Collecting that reflected share into the emission term gives
    // one coefficient per pair, `env.emGround` and `env.emCanopy` below, which
    // the energy balances use in place of a bare emissivity. Two surfaces at
    // the same temperature then exchange nothing, whatever the sky view.
    //
    // Terrain radiates at air temperature. The sky's own emission is the
    // measured `env.Rlw` and needs no such closure.
    void RadlwabsStepCpp(const radsetup& s, double svfa, envstruct& env)
    {
        const double em2 = surfaceEmissivity * surfaceEmissivity;
        // The ground's view out past the foliage: the canopy's diffuse
        // transmission when vegetated, the whole hemisphere on bare ground.
        double trd = s.vegetated ? s.groundRdd : 1.0;
        double RemC = sb * utils::rademCpp(env.tcanopy);  // canopy emission, before emissivity
        double RemT = sb * utils::rademCpp(env.tair);     // terrain emission, before emissivity
        double skyG = trd * svfa, skyC = svfa;            // each surface's view of open sky

        // Absorbed longwave, and the matching coefficient on the surface's own
        // emission, for the ground alone and for the canopy+ground surface.
        double radGlw = surfaceEmissivity * skyG * env.Rlw
            + em2 * (trd * (1.0 - svfa) * RemT + (1.0 - trd) * RemC);
        double radClw = surfaceEmissivity * skyC * env.Rlw + em2 * (1.0 - svfa) * RemT;
        env.emGround = surfaceEmissivity * skyG + em2 * (1.0 - skyG);
        env.emCanopy = surfaceEmissivity * skyC + em2 * (1.0 - skyC);

        env.RabsGround = env.Rground + radGlw; // ground only
        // With no canopy `trd` is 1, so the ground's view and the whole
        // surface's are the same and both quantities above already coincide.
        env.RabsCanopy = s.vegetated ? (env.Rcanopy + radClw) : env.RabsGround; // combined canopy + ground
    }

    // ========================================================================
    // Wind: within- and above-canopy wind speed and stability corrections
    // ========================================================================

    // Solve the canopy aerodynamic state for the current sensible heat flux.
    // Canopy geometry sets displacement and roughness, buoyancy modifies the
    // surface-layer profile through Monin-Obukhov stability, and the resulting
    // friction velocity sets both canopy-top wind and within-canopy mixing.
    void windmodelCpp(const vegpstruct& vegp, double zref, double d, envstruct& env,
        double shelterc, double zi, double beta)
    {
        double hgt = vegp.hgt;
        double pai = vegp.pai;
        // Apply topographic shelter, with the shelter coefficient floored at
        // 0.05 so the driving wind never vanishes.
        if (!std::isfinite(shelterc)) shelterc = 1.0;
        if (shelterc < 0.05) shelterc = 0.05;
        double uref = shelterc * env.uref;
        double H = env.H;
        double tc = env.tair;
        double pk = env.pk;
        // Both corrections are carried between iterations, evaluated at the
        // reference height, and seed the first estimate of friction velocity
        // before the Obukhov length for this step is known.
        double psi_h = env.psi_h;
        double psi_m = env.psi_m;

        if (zref < hgt) {
            throw std::invalid_argument("Height of wind speed measurement must be above canopy");
        }
        double Tk = tc + 273.15;
        // Canopy geometry defines displacement and roughness, while air
        // temperature and pressure set the thermodynamic properties needed
        // to convert temperature differences into sensible heat flux.
        double cp = utils::cpairCpp(tc);
        double ph = utils::phairCpp(tc, pk);
        double zm = utils::roughlengthCpp(hgt, pai, d);
        double zh = utils::scalarRoughlengthCpp(hgt, pai, d);
        // Under free convection (H > 0, light wind), shear-generated
        // turbulence alone can no longer sustain the friction velocity
        // implied by similarity theory. Combine the reference wind speed
        // in quadrature with a free-convective velocity scale (Beljaars,
        // 1994) so uf stays physically bounded as uref -> 0, rather than
        // letting it collapse and aerodynamic resistance diverge.
        double Ueff = uref;
        if (H > 0.0) {
            double wstar = std::cbrt((g / Tk) * zi * (H / (ph * cp)));
            Ueff = std::sqrt(uref * uref + std::pow(beta * wstar, 2.0));
        }
        // The free-convective scale above only helps while the surface is
        // heating. On a cooling surface nothing bounds the driving wind, and
        // below roughly 0.2 m/s sustained the ground decouples far enough that
        // the surface balance becomes nearly singular in its exchange term and
        // stops converging. The same floor the grid model applies is therefore
        // applied here, which also makes the two models agree about what a calm
        // hour is.
        if (Ueff < mc::minWindSpeed) Ueff = mc::minWindSpeed;
        double uf = (ka * Ueff) / (std::log((zref - d) / zm) + psi_m);
        // Monin-Obukhov length expresses the balance between buoyant and
        // mechanically generated turbulence: negative under surface heating,
        // positive under stable cooling, and effectively infinite at neutrality.
        double LL = 1e99;
        if (H != 0.0) LL = (ph * cp * std::pow(uf, 3.0) * Tk) / (-ka * g * H);
        // The safeguard returns the length unchanged unless the stability
        // correction it implies exceeds the bound, and the milder length when
        // it does, so its result is taken as given. Comparing the two and
        // keeping the smaller would reject it: on the stable side the milder
        // length is the larger one.
        LL = utils::clipMOlength(LL, zref, d, zm);
        // Apply the stability corrections implied by the current heat flux;
        // sensible heat and turbulence converge together in the outer surface
        // energy-balance iteration. With no heat flux the surface layer is
        // neutral by definition and both corrections are zero, set here
        // rather than inferred from the magnitude of the neutral sentinel
        // held in LL.
        if (H == 0.0) {
            psi_m = 0.0;
            psi_h = 0.0;
        }
        else {
            psi_m = utils::dpsimCpp(zm / LL) - utils::dpsimCpp((zref - d) / LL);
            psi_h = utils::dpsihCpp(zh / LL) - utils::dpsihCpp((zref - d) / LL);
        }
        uf = (ka * Ueff) / (std::log((zref - d) / zm) + psi_m);
        // In light wind, stability can switch sign from one pass to the next
        // and throw the iteration between an unstable and a strongly stable
        // state. Once that has happened within a timestep, each later pass of
        // the timestep moves friction velocity and 1/L only half-way from the
        // previous pass's values, and the corrections follow the damped length.
        if (env.uref < mc::dampWindSpeed && env.stabPassDone) {
            const double prevL = env.LL, prevUf = env.uf;
            const bool prevFinite = std::abs(prevL) < 1e98 && prevUf > 0.0;
            if (prevFinite && H != 0.0 && ((prevL < 0.0) != (LL < 0.0))) env.stabFlipped = true;
            if (env.stabFlipped && prevFinite) {
                const double invL = 0.5 * (1.0 / prevL + (H != 0.0 ? 1.0 / LL : 0.0));
                LL = (invL != 0.0) ? utils::clipMOlength(1.0 / invL, zref, d, zm) : 1e99;
                uf = 0.5 * (prevUf + uf);
                if (std::abs(LL) < 1e98) {
                    psi_m = utils::dpsimCpp(zm / LL) - utils::dpsimCpp((zref - d) / LL);
                    psi_h = utils::dpsihCpp(zh / LL) - utils::dpsihCpp((zref - d) / LL);
                } else {
                    psi_m = 0.0;
                    psi_h = 0.0;
                }
            }
        }
        env.stabPassDone = true;
        // Canopy-top wind and the within-canopy mixing parameter exist only
        // for vegetated surfaces. Bare ground has no distinct canopy air
        // space, so both quantities are zero and scalar exchange is handled
        // entirely by the above-surface aerodynamic pathway.
        double uh = 0.0;
        if (hgt > 0.0) {
            // Canopy-top wind needs the momentum correction over the shorter
            // roughness-to-canopy-top path, which is a different quantity from
            // the reference-height correction that sets friction velocity. It
            // is held separately so the reference-height value survives to be
            // carried into the next iteration.
            double psi_m_hgt = 0.0;
            if (H != 0.0) {
                psi_m_hgt = utils::dpsimCpp(zm / LL) - utils::dpsimCpp((hgt - d) / LL);
            }
            // Canopy top lies inside the roughness sublayer, whose influence
            // is folded into the roughness length, so it is taken back out
            // again here (Raupach 1992 Eq. 26b). The two uses are a matched
            // pair, and together they realise a friction velocity to
            // canopy-top wind ratio equal to the drag partition itself.
            double lnprofile = std::log((hgt - d) / zm) + utils::momentumInfluenceCpp(pai) + psi_m_hgt;
            // Vegetation only a few centimetres tall can be shallower than the
            // bare-ground roughness length bounding zm from below, which would
            // otherwise drive the profile term to zero or below.
            if (lnprofile < 1e-6) lnprofile = 1e-6;
            uh = (uf / ka) * lnprofile;
        }

        // The profile winds above are built on the effective wind, which carries the
        // free-convective velocity and the floor so that exchange stays well behaved in
        // light wind; they are not the actual wind. The ratio of the measured wind to the
        // profile's own wind at the reference height turns them into actual winds. It is
        // applied only where a wind is reported or used to translate the forcing, never
        // to exchange, and is taken against the profile rather than the effective wind so
        // that it holds after the light-wind damper has moved the friction velocity.
        {
            double psiR = utils::dpsimCpp(zm / LL) - utils::dpsimCpp((zref - d) / LL);
            double uProfileRef = (uf / ka) * (std::log((zref - d) / zm) + psiR);
            env.windScale = (uProfileRef > 0.0) ? uref / uProfileRef : 1.0;
        }
        env.uf = uf;
        env.LL = LL;
        env.uh = uh;
        env.zm = zm;
        env.zh = zh;
        env.psi_m = psi_m;
        env.psi_h = psi_h;
    }

    // Similarity resistance from the canopy heat-exchange roughness height to an
    // arbitrary height, with no sublayer term: the scalar roughness length carries
    // none. Used for the exchange surface's own segment up to canopy top, and for
    // the whole path where there is no canopy column to describe.
    static double rHaSimilarityCpp(const vegpstruct& vegp, double height, double d,
        const envstruct& env)
    {
        // windmodelCpp produces the scalar roughness length for this timestep; it is
        // recomputed here only if this is reached before that has run.
        double zh = (env.zh > 0.0) ? env.zh : utils::scalarRoughlengthCpp(vegp.hgt, vegp.pai, d);
        double psih = utils::dpsihCpp(zh / env.LL) - utils::dpsihCpp((height - d) / env.LL);
        return (std::log((height - d) / zh) + psih) / (ka * env.uf);
    }

    // Aerodynamic resistance from the canopy heat-exchange surface to a height above
    // the canopy, in two parts: the exchange surface's own empirical segment up to
    // canopy top, and above canopy top the transport column's own resistance, which
    // carries the roughness sublayer and reduces to ordinary similarity above it.
    //
    // Taking the upper part from the column rather than continuing the similarity
    // profile is what makes this resistance and the scalar diffusivity field one
    // description. A similarity form cannot reproduce the column there at any
    // stability but neutral: the column carries the sublayer's stability dependence
    // as a factor on a constant-diffusivity layer, similarity as a difference of
    // correction functions, and no choice of roughness length reconciles the two.
    double rHaToHeightCpp(const vegpstruct& vegp, double height, double d, const envstruct& env,
        const utils::ColumnStruct& col, bool hasColumn)
    {
        if (!hasColumn || height <= vegp.hgt) return rHaSimilarityCpp(vegp, height, d, env);
        // The exchange surface sits below canopy top, so its own segment is never
        // negative in any canopy the column describes; the guard covers vegetation
        // barely taller than the soil's own heat roughness height.
        double r1 = rHaSimilarityCpp(vegp, vegp.hgt, d, env);
        if (r1 < 0.0) r1 = 0.0;
        double rTop = utils::columnResistCpp(col, env.uf, env.LL, vegp.hgt);
        double rZ = utils::columnResistCpp(col, env.uf, env.LL, height);
        return r1 + (rZ - rTop);
    }

    // Aerodynamic resistance between the bulk canopy exchange surface and
    // the reference atmosphere used by the canopy energy balance.
    double bulkaeroresistCpp(const vegpstruct& vegp, double zref, double d, const envstruct& env,
        const utils::ColumnStruct& col, bool hasColumn)
    {
        return rHaToHeightCpp(vegp, zref, d, env, col, hasColumn);
    }

    // Total scalar resistance from the soil surface to the reference
    // atmosphere, through the transport column (utils): a wall layer at the
    // soil, the canopy interior above it, the canopy-top sublayer, and the
    // surface layer. Bare ground is the same column at zero plant area, so
    // no separate blend is needed.
    //
    // If rRise is given it receives the column's uniform-source rise
    // resistance: canopy sources, spread uniformly with height, raise the air
    // the ground exchanges with by their total times this, over rho cp. It is
    // zero where there is no column.
    // The column's own constants depend only on canopy height, plant area and
    // displacement height, so they are built once for a run and passed in here
    // rather than rebuilt on every call. `hasColumn` is false only at zero
    // canopy height: the ground is then the whole exchange surface and the
    // single log-law profile is the answer.
    double groundaeroresistCpp(const vegpstruct& vegp, double zref, double d, const envstruct& env,
        const utils::ColumnStruct& col, bool hasColumn, double* rRise)
    {
        if (rRise) *rRise = 0.0;
        if (!hasColumn) {
            return rHaToHeightCpp(vegp, zref, d, env, col, hasColumn);
        }
        if (rRise) *rRise = utils::columnUniformRiseResistCpp(col, env.uf, env.LL, zref);
        return utils::columnResistCpp(col, env.uf, env.LL, zref);
    }

    // ========================================================================
    // Vertical profiles outside the bulk exchange surface
    // ========================================================================
    // Above the canopy, the converged turbulent state is extended vertically:
    // wind by Monin-Obukhov similarity, temperature and humidity by the
    // transport column's resistances, which include the roughness sublayer
    // and reduce to Monin-Obukhov similarity above it. Below ground,
    // temperature and moisture are read directly from the explicitly solved multilayer soil profile.
    // These diagnostics therefore use the converged timestep state rather than
    // introducing an additional energy- or water-balance calculation.

    AboveCanopyPoint aboveCanopyProfileCpp(const vegpstruct& vegp, double z, double d,
        double zref, double Tsurface, const envstruct& env, double rHa, double rV, double hSurf,
        const utils::ColumnStruct& col, bool hasColumn, double LEtotal)
    {
        // Extend the converged friction velocity and stability state to the
        // requested height with the same similarity profile used to relate
        // the canopy to the reference atmosphere.
        double psi_mz = utils::dpsimCpp(env.zm / env.LL) - utils::dpsimCpp((z - d) / env.LL);
        double windz = (env.uf / ka) * (std::log((z - d) / env.zm) + psi_mz);
        // Inside the momentum roughness sublayer the diffusivity is constant, so
        // wind rises linearly from its canopy-top value to the logarithmic
        // profile at the sublayer top (Raupach 1992), rather than following the
        // logarithmic profile, which the sublayer term in zm makes exact only
        // from the sublayer top upward.
        if (vegp.hgt > 0.0 && z >= vegp.hgt) {
            double zw = utils::momentumSublayerTopCpp(vegp.hgt, vegp.pai, d);
            if (z < zw) {
                double psi_mw = utils::dpsimCpp(env.zm / env.LL) - utils::dpsimCpp((zw - d) / env.LL);
                double uw = (env.uf / ka) * (std::log((zw - d) / env.zm) + psi_mw);
                windz = env.uh + (uw - env.uh) * (z - vegp.hgt) / (zw - vegp.hgt);
            }
        }
        double rZ = rHaToHeightCpp(vegp, z, d, env, col, hasColumn);
        double rZref = rHaToHeightCpp(vegp, zref, d, env, col, hasColumn);
        double ratio = rZ / rZref;
        // Above canopy top the exchange-to-reference resistance is the
        // surface's own segment to canopy top plus the transport column above
        // it, so the share of the surface-to-reference difference still to come
        // between z and the reference height is the column's own resistance
        // between the two, inside the roughness sublayer and above it.
        if (hasColumn) {
            double rcolZ = utils::columnResistCpp(col, env.uf, env.LL, z);
            double rcolR = utils::columnResistCpp(col, env.uf, env.LL, zref);
            if (rZref > 0.0) ratio = 1.0 - (rcolR - rcolZ) / rZref;
        }
        double Tz = Tsurface - (Tsurface - env.tair) * ratio;

        // Humidity at z: the vapour pressure at the canopy/ground exchange
        // surface, the same node Tsurface represents, is the reference value
        // raised by the surface's latent flux across rHa. Where the converged
        // total latent flux `LEtotal` (W/m2) is supplied it is used directly,
        // e = ea + LEtotal * rHa / k with k = lambda * rho / p: the vegetated
        // surface's flux has a soil part at the ground's own temperature,
        // which no factor on es(Tsurface) reproduces. Otherwise the flux is
        // (hSurf * es - ea) / rV through two resistances in series, and the
        // vapour pressure at their junction splits the surface value and ea
        // by rHa's share of rV. Extrapolated to z with the same resistance
        // ratio as Tz above, since it is the same constant-flux profile.
        double es = 1000.0 * hSurf * utils::satvapCpp(Tsurface);            // kPa -> Pa
        double ea = 1000.0 * utils::satvapCpp(env.tair) * (env.rh / 100.0); // kPa -> Pa
        double eSurface = ea + (es - ea) * (rHa / rV);
        if (std::isfinite(LEtotal)) {
            double kv = utils::latentHeatCpp(Tsurface) * utils::phairCpp(env.tair, env.pk) / env.pk; // W/m2 per kPa per s/m
            eSurface = ea + 1000.0 * LEtotal * rHa / kv;
        }
        double ez = eSurface - (eSurface - ea) * ratio;
        double esatTz = 1000.0 * utils::satvapCpp(Tz);
        double RHz = 100.0 * ez / esatTz;
        // Temperature and vapour pressure are extrapolated separately, so
        // their conversion back to relative humidity can exceed the intended
        // range under strong stratification. Bound it to the physically
        // admissible interval.
        if (RHz < 0.0) RHz = 0.0;
        if (RHz > 100.0) RHz = 100.0;

        AboveCanopyPoint out;
        // Returned as actual wind (see windmodelCpp); the profile itself is on the
        // exchange scale.
        out.windz = windz * env.windScale;
        out.Tz = Tz;
        out.RHz = RHz;
        return out;
    }

    // E-folding depth of a sinusoidal 24-hour soil-temperature signal,
    // determined by the top layer's thermal diffusivity. Larger thermal
    // diffusivity carries the diurnal temperature wave deeper into the soil.
    double diurnalDampingDepthCpp(const soilpstruct& soilp, double theta, double Tc, double pk)
    {
        double ks = utils::thermalConductivityCpp(soilp.Vq[0], soilp.Vm[0], soilp.Vo[0], theta,
            soilp.Mc[0], Tc, pk);
        double Cs = utils::heatCapacityCpp(soilp.Vq[0], soilp.Vm[0], soilp.Vo[0], theta, Tc, pk);
        double kappa = ks / Cs;
        return std::sqrt(2.0 * kappa / mc::omdy);
    }

    // Read a temperature or moisture value from the solved soil profile at
    // the requested depth, by linear interpolation between the two nodes that
    // bracket it. Node i holds the state at depth z[i]; the surface node
    // answers at or above the surface, and the deepest node at or below its
    // own depth.
    double interpSoilProfileCpp(const std::vector<double>& z,
        const std::vector<double>& values, double reqhgt)
    {
        double depth = -reqhgt; // z is positive-downward; reqhgt is negative below ground
        int n = static_cast<int>(values.size());
        if (depth <= z[0]) return values[0];
        if (depth >= z[n - 1]) return values[n - 1];
        for (int i = 0; i < n - 1; ++i) {
            if (depth >= z[i] && depth <= z[i + 1]) {
                double f = (depth - z[i]) / (z[i + 1] - z[i]);
                return values[i] + f * (values[i + 1] - values[i]);
            }
        }
        return values[n - 1];
    }

    // ========================================================================
    // Stomatal conductance (JULES-SOX)
    // ========================================================================

    // Leaf stomatal conductance from the JULES-SOX optimisation framework.
    // Photosynthetic carbon gain is balanced against the hydraulic cost of
    // supplying transpiration through a vulnerable xylem pathway. Light,
    // temperature, CO2, VPD, root-zone water potential and leaf height all
    // therefore influence the optimum stomatal opening.
    double leafgsCpp(const envstruct& env, vegpstruct& vegp, double z, bool C3)
    {
        double gs = 0.0;
        if (env.PARabs > 0.0) {
            double IPAR = env.PARabs * (0.48 / 0.219) * 1e-6; // conversion to mol photons / m^2 / s
            // Temperature-dependent photosynthetic capacity and dark respiration.
            double xx = 0.1 * (env.tcanopy - 25.0);
            double Vcmax = vegp.Vcmax25 * std::exp2(xx) / ((1.0 + std::exp(0.3 * (env.tcanopy - vegp.Tup))) *
                (1.0 + std::exp(0.3 * (vegp.Tlw - env.tcanopy))));
            double Rd = vegp.fd * Vcmax;
            // CO2 compensation point, including its temperature dependence.
            double Q10rs = 0.57;
            double Oa = 0.2095 * env.pk * 1000;
            double photocomp = Oa / (2.0 * 2600.0 * std::pow(Q10rs, 0.1 * (env.tcanopy - 25.0)));
            // Atmospheric/intercellular CO2 and leaf-to-air vapour pressure deficit.
            double ea = utils::satvapCpp(env.tair) * (env.rh / 100.0);
            double es = utils::satvapCpp(env.tcanopy);
            // Vapour pressure deficit. A leaf at or below the air's dew point
            // loses no water through its stomata, so a negative deficit is
            // taken as zero: the water cost of opening vanishes and conductance
            // is at its maximum, continuous with its limit as the deficit falls
            // to zero from above.
            double DD = std::max(es - ea, 0.0);
            double ca = env.Ca * env.pk / 1000.0;  // convert CO2 concentration to Pa
            double ci = vegp.f0 * (1.0 - DD / vegp.Dcrit) * (ca - photocomp); // in Pa
            double Aca;
            double Acol;
            double cicol;
            if (C3) {
                // C3 Rubisco kinetics and the three potential limitations on
                // assimilation: carboxylation, light and electron transport.
                double Q10Kc = 2.1;
                double Kc = 30.0 * std::pow(Q10Kc, 0.1 * (env.tair - 25.0));
                double Q10Ko = 1.2;
                double Ko = 30000.0 * std::pow(Q10Ko, 0.1 * (env.tair - 25.0));
                // Gross assimilation, limited by carbon, light or transport --
                // evaluated at ci (the Jacobs (1994) formula). This retains
                // the physically-meaningful ci-response of Wl (real leaves
                // draw ci down below ca as light/assimilation rises, which
                // is what gives the light-response curve its shape) and
                // feeds the co-limitation smoothing (Wcol/cicol) below.
                double Wc = Vcmax * ((ci - photocomp) / (ci + Kc * (1 + Oa / Ko)));
                double Wl = vegp.alpha * IPAR * ((ci - photocomp) / (ci + 2.0 * photocomp));
                double We = 0.5 * Vcmax;
                if (Wc < 0.0) Wc = 0.0;
                if (Wl < 0.0) Wl = 0.0;
                if (We < 0.0) We = 0.0;
                // Co-limiting assimilation and CO2 concentration -- still
                // derived from the ci-evaluated We/Wl (see below for why).
                double Wcol = ((We + Wl) - std::sqrt(std::pow(We + Wl, 2.0) - 4.0 * 0.93 * (We * Wl))) / (2.0 * 0.93);
                Acol = Wcol - Rd;
                cicol = (-Vcmax * photocomp - Kc * (1.0 + Oa / Ko) * Wcol) / (Wcol - Vcmax);
                // Assimilation at ambient CO2 represents the fully open-
                // stomata end member. Together with the co-limitation point
                // above it gives the marginal carbon gain from increasing
                // intercellular CO2, required by the SOX optimum.
                double Wc_ca = Vcmax * ((ca - photocomp) / (ca + Kc * (1 + Oa / Ko)));
                double Wl_ca = vegp.alpha * IPAR * ((ca - photocomp) / (ca + 2.0 * photocomp));
                double We_ca = We;
                if (Wc_ca < 0.0) Wc_ca = 0.0;
                if (Wl_ca < 0.0) Wl_ca = 0.0;
                double Wca = Wc_ca;
                if (Wl_ca < Wca) Wca = Wl_ca;
                if (We_ca < Wca) Wca = We_ca;
                Aca = Wca - Rd;
            }
            else {
                // C4 photosynthetic pathway: PEP-carboxylase-limited
                // assimilation, using a simpler light/Rubisco/CO2
                // co-limitation than the C3 branch above (C4 plants
                // concentrate CO2 internally, so their assimilation is not
                // CO2-limited in the same way C3 photosynthesis is).
                double k = 2.0e-4;
                double Wc = Vcmax;
                double Wl = vegp.alpha * IPAR;
                double We = k * Vcmax * (ci / (env.pk * 1000.0));
                double W = Wc;
                if (Wl < W) W = Wl;
                if (We < W) W = We;
                Aca = W - Rd;
                double Wcol = ((Wc + Wl) - std::sqrt(std::pow(Wc + Wl, 2.0) - 4.0 * 0.83 * (We * Wl))) / (2.0 * 0.83);
                Acol = Wcol - Rd;
                cicol = (Wcol * env.pk * 1000.0) / (k * Vcmax);
            }
            // Change in assimilation per unit change in intercellular CO2:
            // a finite difference between assimilation evaluated at ambient
            // CO2 (ca, the fully-open-stomata case) and at the co-limitation
            // point (cicol) -- the marginal carbon gain the stomatal-
            // optimisation scheme below weighs against the marginal water
            // cost of opening the stomata further.
            double dadc = (Aca - Acol) / (ca - cicol);
            // Hydraulic supply cost. Leaf water potential becomes more
            // negative with height, and the vulnerability curve translates
            // that water potential into the fraction of xylem conductivity
            // still available to support transpiration.
            if (vegp.apsi < 0.0) {
                double stem_slope = 65.15 * std::pow(-vegp.psi50, -1.25);
                vegp.apsi = -4.0 * stem_slope / 100.0 * vegp.psi50;
            }
            double rhow = 1000.0 * (1 - (env.tcanopy + 288.9414) * std::pow(env.tcanopy - 3.9863, 2.0) /
                (508929.2 * (env.tcanopy + 68.12963)));
            double psi_pd = env.psi_r - z * g * rhow * 1e-6;
            double K_psi_pd = 1.0 / (1.0 + std::pow(psi_pd / vegp.psi50, vegp.apsi));
            double K_50f = 0.5 / (1.0 + std::pow((psi_pd + vegp.psi50) / vegp.psi50, vegp.apsi));
            double psi_50f = (psi_pd + vegp.psi50) / 2.0;
            double dKdpKi = ((K_psi_pd - K_50f) / (psi_pd - psi_50f)) * (1.0 / K_psi_pd);
            // Whole-plant hydraulic resistance increases as xylem
            // conductivity is lost relative to its unstressed maximum.
            double rp = vegp.rpmin / K_psi_pd;
            // Solve the coupled carbon-gain / hydraulic-cost optimum for
            // stomatal conductance.
            double DDm = DD / env.pk; // divide by pk to convert to mol / mol
            double zeta = 2.0 / (dKdpKi * rp * 1.6 * DDm);
            if (dadc < 1e-99) dadc = 1e-99;
            double mu = 1.0 + (4.0 * zeta) / dadc;
            if (mu < 1.0) mu = 1.0;
            gs = 0.5 * dadc * (std::sqrt(mu) - 1.0);
            // Keep the optimum within the conductance supported by the
            // photosynthetic parameterisation.
            double gsmax = vegp.Vcmax25 * 1e4;
            gs = std::fmin(gs, gsmax);
            // An optional measured/empirical maximum can impose a stricter
            // ceiling; negative means no additional cap.
            if (vegp.gsmaxCap >= 0.0) gs = std::fmin(gs, vegp.gsmaxCap);
        }
        return gs;
    }

    // Partition living foliage into sunlit and shaded fractions and assign
    // absorbed PAR to each. The simple Beer's-law split is rescaled so that
    // its leaf-area-weighted total exactly matches the clumping-aware
    // two-stream PAR absorbed by foliage, preserving the correct canopy-scale
    // energy while retaining a cheap sunlit/shaded representation.
    stomatalsetup StomatalSetupCpp(const solmodel& solp, const envstruct& env, const vegpstruct& vegp)
    {
        stomatalsetup out;
        out.z = vegp.hgt / 2.0;
        double om = vegp.lrefp + vegp.ltrap;
        // Canopy extinction coefficient (zenith clamped to the horizon)
        double zend = solp.zend;
        if (zend > 90.0) zend = 90.0;
        double zenr = zend * torad;
        double k = utils::canopyKCpp(zenr, vegp.x);
        double PAI = vegp.pai;
        // Beer's-law sunlit/shaded leaf-area partition. Clumping is not
        // represented inside this local split; instead the resulting PAR
        // totals are reconciled with the full clumping-aware two-stream
        // solution below.
        double L_sun = (1.0 - std::exp(-k * PAI * vegp.Lfrac)) / k;
        double L_shade = (PAI * vegp.Lfrac) - L_sun;
        // First-pass PAR absorbed per unit sunlit and shaded leaf area.
        double Rshade_abs = env.Rdif * ((1.0 - std::exp(-PAI)) / PAI) * (1.0 - om);
        double Rsun_abs = (env.Rsw - env.Rdif) * k * (1.0 - om) + Rshade_abs;
        // Scale both leaf classes together so their weighted sum equals
        // the foliage-only PAR absorption from the full two-stream model.
        double unclumpedTotal = Rsun_abs * L_sun + Rshade_abs * L_shade;
        double targetTotal = env.RcanopyPAR - env.RgroundPAR;
        double rescale = (unclumpedTotal > 0.0) ? (targetTotal / unclumpedTotal) : 1.0;
        out.L_sun = L_sun;
        out.L_shade = L_shade;
        out.PARabs_sun = Rsun_abs * rescale;
        out.PARabs_shade = Rshade_abs * rescale;
        return out;
    }

    // Integrate leaf conductance across sunlit and shaded foliage to obtain
    // the bulk canopy stomatal resistance seen by the surface energy balance.
    double bulkstomatalresistCpp(const stomatalsetup& setup, const envstruct& env, vegpstruct& vegp, bool C3)
    {
        envstruct envSun = env;
        envSun.PARabs = setup.PARabs_sun;
        envstruct envShade = env;
        envShade.PARabs = setup.PARabs_shade;
        // Convert total canopy conductance (mol m^-2 s^-1) to an
        // aerodynamic-style resistance in s/m using molar air density. Each
        // leaf's conductance is at least its closed-stomata value, which it
        // reaches in the dark and approaches as light fades.
        double ph = utils::phairCpp(env.tair, env.pk);
        double gsMin = ph / mc::rStomClosed;
        double gs_sun = std::max(leafgsCpp(envSun, vegp, setup.z, C3), gsMin);
        double gs_shade = std::max(leafgsCpp(envShade, vegp, setup.z, C3), gsMin);
        double Gs = gs_sun * setup.L_sun + gs_shade * setup.L_shade;
        return ph / Gs;
    }

    // ========================================================================
    // Canopy surface energy balance and canopy water budget
    // ========================================================================

    // Couple canopy interception, wet-surface evaporation and stomatal
    // transpiration through the vapour network (see canopywaterresult). Rain is
    // split between direct throughfall and temporary canopy storage; stored
    // water evaporates without stomatal resistance, whereas dry-canopy
    // transpiration is controlled by stomata. Foliage and soil exhale into one
    // node of canopy air, so the foliage's evaporation, the soil's and the
    // whole surface's all come from the same closed form. precipGround is the
    // water reaching the soil.
    canopywaterresult canopyWaterBudgetCpp(vegpstruct& vegp, const envstruct& env,
        const stomatalsetup& stomSetup, bool vegetated, double rHa, double hr,
        double swaterdepth, double dT, double eg, double rsh, double RgR)
    {
        canopywaterresult out;
        if (!vegetated) {
            // Bare ground has no interception or stomatal pathway. The soil
            // surface is the exchange surface, so its vapour crosses rHa alone
            // at the pore air's humidity.
            out.rSurf = 0.0;
            out.hSurf = hr;
            out.swaterdepth = 0.0;
            out.Et = 0.0;
            out.precipGround = env.precip;
            return out;
        }

        double rStom = bulkstomatalresistCpp(stomSetup, env, vegp);
        double ec = utils::satvapCpp(env.tcanopy);                       // kPa
        double ea = utils::satvapCpp(env.tair) * (env.rh / 100.0);        // kPa
        double Tk = env.tair + 273.15;
        double cm = 1000.0 * Mw * dT / (RgasC * Tk);                      // kPa/(s/m) -> mm over the step

        // The three legs. The wet leg's aerodynamic part vanishes where the node
        // is the exchange surface; its floor is numerical and moves the film's
        // flux by about 0.1/rsh.
        const bool hasNode = rsh > 0.0;
        double rw = std::max(rHa - rsh, 0.0);
        double rb = std::max(RgR - rsh, 1e-9);
        double gDryLeg = std::isfinite(rStom) ? 1.0 / (rStom + rw) : 0.0;
        double gWetLeg = 1.0 / std::max(rw, 0.1);
        double b = 1.0 / rb;
        double c = hasNode ? 1.0 / rsh : 0.0;

        // Share of the surface that is wet: all of it while rain falls, and
        // afterwards a share that shrinks as the film held per unit plant area
        // evaporates. The store itself is held per unit ground area.
        double wetShare = (env.precip > 0.0) ? 1.0 : 1.0 - std::exp(-30.0 * swaterdepth / vegp.pai);

        // Interception of rainfall uses the same foliage-angle extinction
        // geometry as radiation, with the trajectory of falling drops set by
        // rain terminal velocity and canopy-top wind. The un-intercepted
        // fraction becomes direct throughfall.
        double tr = 1.0;
        if (env.precip > 0.0) {
            double vr = 3.78 * std::pow(env.precip, 0.067); // empirical rain terminal velocity (m/s)
            double rainZ = std::atan(env.uh / vr);           // rain angle from vertical
            double si = std::cos(rainZ);                     // the canopy is treated as horizontal for rain; slope does not enter
            utils::kstruct kp = utils::cankCpp(rainZ, vegp.x, si);
            tr = std::exp(-kp.kd * vegp.pai);
        }

        // Rain keeps arriving through the step, and whatever the film
        // evaporates makes room for more, so the water the film can give up
        // over the step is what it held plus all the rain it intercepted.
        // The film's capacity (mwft per unit plant area) limits only what is
        // left at the end of the step; the rest drips to the ground.
        double capacity = vegp.mwft * vegp.pai;
        double intercepted = (1.0 - tr) * env.precip;
        double available = swaterdepth + intercepted;

        // The film's evaporation at a wet share wt: with a(wt) the foliage
        // conductance, ec - e* = N / (a + b + c) where N does not depend on wt,
        // so the film's flux is wt*gWetLeg*N/(S0 + wt*dA) and the share the
        // water available over the step can sustain has a closed form.
        // Condensation keeps the wet share wet all step. Without a node the
        // film exhales to the reference air.
        auto nodeAt = [&](double a) {
            return hasNode ? (a * ec + b * eg + c * ea) / (a + b + c) : ea;
        };
        double wt = wetShare;
        if (wetShare > 0.0) {
            double aW = wetShare * gWetLeg + (1.0 - wetShare) * gDryLeg;
            double filmW = cm * wetShare * gWetLeg * (ec - nodeAt(aW));
            if (filmW > 0.0 && filmW > available) {
                if (hasNode) {
                    double N = b * (ec - eg) + c * (ec - ea);
                    double S0 = gDryLeg + b + c;
                    double dA = gWetLeg - gDryLeg;
                    double A = available / (cm * gWetLeg * N);
                    double den = 1.0 - A * dA;
                    wt = (den > 0.0) ? A * S0 / den : wetShare;
                } else {
                    wt = available / (cm * gWetLeg * (ec - ea));
                }
                if (!(wt >= 0.0) || wt > wetShare) wt = wetShare;
            }
        }

        double a = wt * gWetLeg + (1.0 - wt) * gDryLeg;
        double est = nodeAt(a);
        double filmEvap = cm * wt * gWetLeg * (ec - est);
        double Et = cm * (1.0 - wt) * gDryLeg * (ec - est);

        // The whole surface's latent flux in Penman-Monteith form,
        // (hSurf * es(Tc) - eV) / rV, which keeps its exact slope in Tc.
        double rV = 1.0 / (a + b) + rsh;
        out.rSurf = rV - rHa;
        out.hSurf = a / (a + b);
        out.eV = ea - (b / (a + b)) * eg;
        // The network seen from the soil surface: its Thevenin equivalent,
        // which gives the node's soil evaporation at any soil vapour pressure.
        out.eTh = hasNode ? (a * ec + c * ea) / (a + c) : ea;
        out.rTh = rb + (hasNode ? 1.0 / (a + c) : 0.0);
        out.aNet = a;
        out.bNet = b;

        // End-of-step store and drip. Rain = throughfall + drip + change in
        // store + film evaporation, exactly; condensed water is kept the same way.
        double remaining = available - filmEvap;
        if (remaining < 0.0) remaining = 0.0;
        double swaterdepthNew = std::min(remaining, capacity);
        double drip = remaining - swaterdepthNew;
        out.swaterdepth = swaterdepthNew;
        out.filmEvap = filmEvap;
        out.filmAvailable = available;
        out.wetShare = wetShare;
        out.Et = Et;
        out.precipGround = tr * env.precip + drip;
        return out;
    }

    // ========================================================================
    // Soil: heat and water balance for a layered soil profile
    // ========================================================================

    // Relative humidity of pore air at the soil surface from the Kelvin
    // relationship between matric water potential and vapour pressure.
    // Dry soil therefore presents a lower effective vapour pressure than a
    // saturated surface at the same temperature.
    double soilrelhumCpp(const soilpstruct& soilp, double Tsoil, double theta)
    {
        double psiw = waterPotentialCpp(soilp, theta, 0);
        double Tk = Tsoil + 273.15;
        double hr = std::exp(Mw * psiw / (RgasC * Tk));
        return hr;
    }

    // Net energy available to enter the soil at a trial surface temperature:
    // absorbed radiation minus emitted longwave, sensible heat and latent
    // heat. Positive values imply net downward heat flux into the soil.
    double soilsurfaceEBCpp(const soilpstruct& soilp, const envstruct& env, double Tsurface, double theta)
    {
        // Net radiation. `env.emGround` is the coefficient that belongs with
        // the absorbed longwave RadlwabsStepCpp supplied: the ground's emission
        // net of the share the canopy and terrain return to it.
        double Rnet = env.RabsGround - env.emGround * sb * utils::rademCpp(Tsurface);
        // Turbulent sensible heat exchange with the reference air.
        double cp = utils::cpairCpp(env.tair);
        double ph = utils::phairCpp(env.tair, env.pk);
        double H = ((ph * cp) / env.rHa) * (Tsurface - env.tair);
        // Moisture-limited evaporation from the soil surface.
        double hr = soilrelhumCpp(soilp, Tsurface, theta);
        // The whole surface's balance, set for the matching constraint (see
        // matchResidualCpp): foliage and soil as the network's two sources at
        // one temperature, in parallel, with factor (a + b hr)/(a + b).
        if (env.wholeB >= 0.0 && env.wholeA >= 0.0 && env.wholeA + env.wholeB > 0.0)
            hr = (env.wholeA + env.wholeB * hr) / (env.wholeA + env.wholeB);
        double es = utils::satvapCpp(Tsurface) * hr;
        double ea = utils::satvapCpp(env.tair) * env.rh / 100.0;
        double la = utils::latentHeatCpp(Tsurface);
        double rVs = (env.rVg > 0.0) ? env.rVg : env.rHa;
        double L = ((la * ph) / (rVs * env.pk)) * (es - ea);
        double Ba = Rnet - H - L;
        return Ba;
    }

    // Equilibrium ground-surface temperature if no heat is conducted into
    // or out of the deeper soil (G = 0). This isolates the instantaneous
    // surface energy balance from subsurface heat storage.
    double groundTemp0Cpp(const soilpstruct& soilp, const envstruct& env, double theta,
        double Tguess, int maxIter, double tolerance)
    {
        double Ts = Tguess;
        for (int iter = 0; iter < maxIter; ++iter) {
            double Ba = soilsurfaceEBCpp(soilp, env, Ts, theta);
            // Centred finite-difference derivative of Ba w.r.t. Ts (same
            // technique utils::penmanMonteithCpp uses for its own slope term).
            double BaPlus = soilsurfaceEBCpp(soilp, env, Ts + 0.5, theta);
            double BaMinus = soilsurfaceEBCpp(soilp, env, Ts - 0.5, theta);
            double dBa = BaPlus - BaMinus;
            if (std::abs(dBa) < 1e-9) break; // Ba(Ts) is smooth and monotonic over any
                                              // physically realistic range, so this shouldn't
                                              // trigger in practice -- guards only against a
                                              // degenerate input (e.g. a zero resistance).
            double step = Ba / dBa;
            // Clamp the step: Ba(Ts) is only close to linear over a modest range because of
            // the T^4 emitted-longwave term, so an unclamped Newton step can overshoot badly
            // starting from a poor Tguess.
            if (step > 10.0) step = 10.0;
            if (step < -10.0) step = -10.0;
            Ts -= step;
            if (std::abs(step) < tolerance) break;
        }
        return Ts;
    }

    // Ground-surface temperature consistent with a specified conductive
    // heat flux G into the soil, solved from the full surface energy balance.
    double groundTempGCpp(const soilpstruct& soilp, const envstruct& env, double theta,
        double G, double Tguess, int maxIter, double tolerance)
    {
        double Ts = Tguess;
        for (int iter = 0; iter < maxIter; ++iter) {
            double Ba = soilsurfaceEBCpp(soilp, env, Ts, theta) - G;
            double BaPlus = soilsurfaceEBCpp(soilp, env, Ts + 0.5, theta) - G;
            double BaMinus = soilsurfaceEBCpp(soilp, env, Ts - 0.5, theta) - G;
            double dBa = BaPlus - BaMinus;
            if (std::abs(dBa) < 1e-9) break;
            double step = Ba / dBa;
            if (step > 10.0) step = 10.0;
            if (step < -10.0) step = -10.0;
            Ts -= step;
            if (std::abs(step) < tolerance) break;
        }
        return Ts;
    }

    // Effective saturation of a soil layer from matric water potential,
    // using the Campbell retention curve and capping at full saturation.
    double degreeOfSaturationCpp(const soilpstruct& soilp, double psiw, int i)
    {
        if (psiw >= 0) return 1.0;
        double Se;
        if (psiw >= soilp.psie[i]) {
            Se = 1.0;
        }
        else {
            Se = std::pow(psiw / soilp.psie[i], -1.0 / soilp.b[i]);
        }
        return Se;
    }

    // Campbell retention curve in the forward direction: convert volumetric
    // water content to matric water potential for soil layer i.
    double waterPotentialCpp(const soilpstruct& soilp, double theta, int i)
    {
        double Se = 1.0;
        if (theta < soilp.thetaS[i]) Se = theta / soilp.thetaS[i];
        double psiw = soilp.psie[i] * std::pow(Se, -soilp.b[i]);
        return psiw;
    }

    // Root-zone water potential seen by the plant: a root-fraction-weighted
    // mean of layer water potentials, returned in MPa for the hydraulic model.
    double rootzonePsiCpp(const soilpstruct& soilp, const std::vector<double>& theta, const std::vector<double>& rootfrac)
    {
        double psir = 0.0;
        int n = static_cast<int>(rootfrac.size());
        for (int i = 0; i < n; ++i) {
            psir += rootfrac[i] * waterPotentialCpp(soilp, theta[i], i);
        }
        return psir / 1000.0; // J/kg -> MPa
    }

    // Inverse Campbell retention curve: water content corresponding to a
    // specified matric water potential in layer i.
    double thetaFromPsiCpp(const soilpstruct& soilp, double psiw, int i)
    {
        double Se = degreeOfSaturationCpp(soilp, psiw, i);
        double theta = Se * soilp.thetaS[i];
        return theta;
    }

    // Liquid hydraulic conductivity declines steeply as a Campbell soil
    // dries, scaling saturated conductivity by relative water content.
    double hydraulicConductivityFromThetaCpp(const soilpstruct& soilp, double theta, int i)
    {
        double n = 2.0 * soilp.b[i] + 3.0;
        double k = soilp.Ksat[i] * std::pow(theta / soilp.thetaS[i], n);
        return k;
    }

    // Water vapour stored in air-filled pore space. Matric potential lowers
    // pore-air humidity below saturation, while drier soil provides more
    // air-filled volume in which vapour can reside.
    double vaporFromPsiCpp(const soilpstruct& soilp, double psiw, double theta, double Tk, int i)
    {
        double humidity = std::exp(Mw * psiw / (RgasC * Tk));
        double vapor = (soilp.thetaS[i] - theta) * utils::satVapDensityCpp(Tk) * humidity;
        return vapor;
    }

    // Differential liquid-water storage capacity d(theta)/d(psi), used in
    // the Jacobian of the nonlinear soil-water solve.
    double dThetaDPsiCpp(const soilpstruct& soilp, double psiw, int i)
    {
        double psie = soilp.psie[i]; // assumed negative
        if (psiw >= psie) return 0.0;
        double Se = std::pow(psiw / psie, -1.0 / soilp.b[i]);
        double theta = soilp.thetaS[i] * Se;
        return -theta / (soilp.b[i] * psiw);
    }

    // Vapour-phase contribution to vertical water transport through the
    // air-filled pore network, added to the liquid hydraulic conductivity.
    double vaporConductivityFromPsiThetaCpp(const soilpstruct& soilp, double psiw, double theta, double Tk, int i)
    {
        const double dv = 0.000024;
        double humidity = std::exp(Mw * psiw / (RgasC * Tk));
        double vp = utils::satVapDensityCpp(Tk);
        double k = 0.66 * (soilp.thetaS[i] - theta) * dv * vp * humidity * Mw / (RgasC * Tk);
        return k;
    }

    // Differential vapour-storage capacity with respect to water potential,
    // including both humidity and changing air-filled pore volume.
    double dVaporDPsiCpp(const soilpstruct& soilp, double psiw, double theta, double Tk, int i)
    {
        double humidity = std::exp(Mw * psiw / (RgasC * Tk));
        double vp = utils::satVapDensityCpp(Tk);
        double capacity_vapor = (soilp.thetaS[i] - theta) * vp * humidity *
            (Mw / (RgasC * Tk)) - dThetaDPsiCpp(soilp, psiw, i) * vp * humidity;
        return capacity_vapor;
    }

    // Soil evaporation over one timestep, driven by the vapour-pressure
    // gradient from moisture-limited pore air at the surface to the air the
    // soil exchanges with, across the soil's vapour resistance. The pore air's
    // Kelvin humidity is taken at the surface temperature, as the surface
    // energy balance takes it.
    double evaporationFluxCpp(const soilpstruct& soilp, const envstruct& env, double theta, double Tsurface, double dT)
    {
        double psiw = waterPotentialCpp(soilp, theta, 0); // J/kg
        double Tk = env.tair + 273.15; // deg C -> K
        double hs = std::exp(Mw * psiw / (RgasC * (Tsurface + 273.15))); // soil effective relative humidity (0-1)
        double es = 1000.0 * utils::satvapCpp(Tsurface) * hs;              // kPa -> Pa
        double ea = 1000.0 * utils::satvapCpp(env.tair) * (env.rh / 100.0); // kPa -> Pa
        double rV = (env.rVg > 0.0) ? env.rVg : env.rHa;
        double Ev = (Mw / (rV * RgasC * Tk)) * (es - ea) * dT; // soil evaporation (kg/m^2/s -> mm)
        return Ev;
    }

    // Allocate the canopy's total transpiration demand among soil layers.
    // Potential uptake follows root abundance but is reduced where a layer
    // is either waterlogged or too dry; the remaining weights are normalised
    // so that available layers collectively meet the imposed demand.
    std::vector<double> transpirationDistributeCpp(const soilpstruct& soilp, const std::vector<double>& rootfrac,
        double totalTransp_mm, double dT, const std::vector<double>& psiw, double p)
    {
        int n = static_cast<int>(psiw.size());
        std::vector<double> S(n);
        if (totalTransp_mm > 0.0) {
            std::vector<double> w(n);
            double Trate = totalTransp_mm / dT;
            for (int i = 0; i < n; ++i) {
                // Water available to roots spans field capacity to the wilting
                // point, and `p` is the share of it a plant draws on before
                // uptake from a layer begins to fall away. It falls to zero only
                // at a much lower potential, so a layer drier than the wilting
                // point still supplies roots that reach it. All three limits are
                // potentials, the same in every soil; each soil's own retention
                // curve turns them into water contents.
                double theta_fc = thetaFromPsiCpp(soilp, mc::psiFieldCapacity, i);
                double theta_wilt = thetaFromPsiCpp(soilp, mc::psiWiltingPoint, i);
                double psi_dry = waterPotentialCpp(soilp, theta_fc - p * (theta_fc - theta_wilt), i);
                double aw = utils::alphaWetCpp(psiw[i], soilp.psie[i]);
                double ad = utils::alphaDryCpp(psiw[i], psi_dry, mc::psiUptakeLimit);
                double alpha = aw * ad;
                w[i] = rootfrac[i] * alpha;
            }
            double sumw = std::accumulate(w.begin(), w.end(), 0.0);
            for (int i = 0; i < n; ++i) {
                w[i] = (sumw > 0.0) ? (w[i] / sumw) : 0.0;
                S[i] = w[i] * Trate;
            }
        }
        else {
            for (int i = 0; i < n; ++i) S[i] = 0.0;
        }
        return S;
    }

    matchresidual matchResidualCpp(const soilpstruct& soilp, const envstruct& env, const envstruct* envS,
        double w, double T, double theta, double Cp)
    {
        matchresidual r;
        double qg = 0.0, lg = 0.0;
        if (w > 0.0) {
            qg = soilsurfaceEBCpp(soilp, env, T, theta);
            lg = -(soilsurfaceEBCpp(soilp, env, T + 0.5, theta) - soilsurfaceEBCpp(soilp, env, T - 0.5, theta));
        }
        if (envS == nullptr || w >= 1.0) { r.q = qg; r.lam = lg; return r; }
        double qs = soilsurfaceEBCpp(soilp, *envS, T, theta);
        double ls = -(soilsurfaceEBCpp(soilp, *envS, T + 0.5, theta) - soilsurfaceEBCpp(soilp, *envS, T - 0.5, theta));
        if (w <= 0.0) { r.q = qs; r.lam = ls; r.kappa = 0.0; return r; }
        double kappa = (lg + Cp) / (ls + Cp);
        double den = w + (1.0 - w) * kappa;
        r.q = (w * qg + (1.0 - w) * kappa * qs) / den;
        r.lam = (w * lg + (1.0 - w) * kappa * ls) / den;
        r.kappa = kappa;
        return r;
    }

    // Advance the layered soil temperature profile through one timestep.
    // Surface temperature is coupled to radiation, sensible and latent heat,
    // while heat diffuses vertically according to moisture- and temperature-
    // dependent conductivity and heat capacity. The nonlinear profile is
    // solved iteratively with a fixed deep-temperature boundary -- there is
    // no free-drainage-style flux option for heat as there is for water, so
    // a shallow profile relative to the timestep can show an artificially
    // damped deep-layer response.
    soilheatmod SoilHeatCpp(soilheatmod state, const soilpstruct& soilp, const envstruct& env,
        double dT, double Fact, int maxNrIterations, double tolerance, const envstruct* envS, double wMatch)
    {
        int n = state.n;
        double boundaryT = state.oldTe[n];
        std::vector<double> ff(n + 1), CT(n + 1), lambda(n + 1);
        std::vector<double> aa(n + 1), bb(n + 1), cc(n + 1), dd(n + 1);
        const double gg = 1.0 - Fact;
        // The old profile is the storage reference for the entire timestep;
        // inner iterations refine the new state rather than advancing time.
        const std::vector<double> oldTe_fixed = state.oldTe;
        std::vector<double> Te_new = state.Te;
        std::vector<double> Te_prev = Te_new;
        std::vector<double> wc = state.wc;
        std::vector<double> dz = state.dz;
        std::vector<double> Vq = soilp.Vq;
        std::vector<double> Vm = soilp.Vm;
        std::vector<double> Vo = soilp.Vo;
        std::vector<double> Mc = soilp.Mc;
        int nrIterations = 0;
        double maxdT = 1e99;
        double qsurface = 0.0;
        double Tsurf_iter = state.Te[0];
        double lamS = 0.0;
        while (maxdT > tolerance && nrIterations < maxNrIterations) {
            // Surface boundary heat flux at the current surface temperature
            // estimate, together with the rate at which the flux weakens as
            // that temperature rises. Carrying the slope into the surface row
            // below solves the boundary implicitly. Every term of the balance
            // -- emitted longwave, sensible heat, and latent heat through
            // saturation vapour pressure -- opposes a rise in surface
            // temperature, so holding the flux fixed while solving for that
            // temperature drives the iteration away from the solution wherever
            // the ground is tightly coupled to the air. Both added terms cancel
            // once the surface stops changing, so the converged state still
            // satisfies the original nonlinear balance exactly: they change how
            // the solution is reached, not what it converges to.
            for (int i = 0; i <= n; ++i) {
                lambda[i] = utils::thermalConductivityCpp(Vq[i], Vm[i], Vo[i], wc[i], Mc[i], Te_new[i], env.pk);
                CT[i] = utils::heatCapacityCpp(Vq[i], Vm[i], Vo[i], wc[i], Te_new[i], env.pk) * state.vol[i];
            }
            // Conductive coupling between adjacent soil layers.
            ff[0] = utils::kMeanCpp("LOGARITHMIC", lambda[0], lambda[1]) / dz[0];
            for (int i = 1; i < n; ++i) ff[i] = lambda[i] / dz[i];
            if (envS != nullptr && wMatch < 1.0) {
                // Short vegetation: the row solves the matching constraint, the
                // ground's balance blended toward the whole surface's.
                matchresidual mr = matchResidualCpp(soilp, env, envS, wMatch, Tsurf_iter, wc[0], CT[0] / dT + ff[0]);
                qsurface = mr.q;
                lamS = mr.lam;
            } else {
                qsurface = soilsurfaceEBCpp(soilp, env, Tsurf_iter, wc[0]);
                // Centred difference over 1 K, the convention groundTemp0Cpp and
                // utils::penmanMonteithCpp already use on this same balance.
                double qplus = soilsurfaceEBCpp(soilp, env, Tsurf_iter + 0.5, wc[0]);
                double qminus = soilsurfaceEBCpp(soilp, env, Tsurf_iter - 0.5, wc[0]);
                lamS = -(qplus - qminus);
            }
            // Emitted longwave alone guarantees this much restoring strength,
            // view-weighted exactly as the emitted term in the balance is.
            // The floor is a backstop against a non-finite or non-positive
            // difference, which would otherwise leave a zero slope and a lagged
            // boundary; the measured slope exceeds it in ordinary cases.
            double lamRad = 4.0 * env.emGround * sb *
                std::pow(Tsurf_iter + 273.15, 3.0);
            if (!(lamS > lamRad)) lamS = lamRad;
            // Assemble heat storage and conduction for the surface, interior
            // nodes 1..n-1 and the fixed-temperature lower boundary at node n.
            for (int i = 0; i <= n; ++i) {
                if (i == 0) {
                    aa[i] = 0.0;
                    bb[i] = CT[i] / dT + ff[i] + lamS;
                    cc[i] = -ff[i];
                    dd[i] = CT[i] / dT * oldTe_fixed[i] + qsurface + lamS * Tsurf_iter;
                }
                else if (i < n) {
                    aa[i] = -ff[i - 1] * Fact;
                    bb[i] = CT[i] / dT + (ff[i - 1] + ff[i]) * Fact;
                    cc[i] = -ff[i] * Fact;
                    dd[i] = CT[i] / dT * oldTe_fixed[i] + gg * (
                        ff[i - 1] * oldTe_fixed[i - 1] +
                        ff[i] * oldTe_fixed[i + 1] - (ff[i - 1] + ff[i]) * oldTe_fixed[i]);
                }
                else {
                    aa[i] = 0.0;
                    bb[i] = 1.0;
                    cc[i] = 0.0;
                    dd[i] = boundaryT;
                }
            }
            Te_prev = Te_new;
            utils::ThomasResult TBC = utils::thomasSolveCpp(aa, bb, cc, dd, Te_new, 0, n);
            Te_new = TBC.x;
            // The surface balance was linearised about Tsurf_iter, and the
            // penalty lamS*(Tsurf_iter - Te_new[0]) it carries into the surface
            // row vanishes only when the two agree. That mismatch is therefore
            // part of the convergence test, alongside the profile's own change.
            maxdT = std::abs(Te_new[0] - Tsurf_iter);
            Tsurf_iter = Te_new[0];
            // Converge when no soil node changes appreciably between iterations.
            for (int i = 0; i <= n; ++i) {
                double d = std::abs(Te_new[i] - Te_prev[i]);
                if (d > maxdT) maxdT = d;
            }
            ++nrIterations;
        }
        state.Te = Te_new;
        state.iters = nrIterations;
        state.lamS = lamS;
        state.condSlope = CT[0] / dT + ff[0];
        return state;
    }

    // Maximum rainfall that can enter the soil during this timestep without
    // forcing the surface layer beyond saturation. It evaluates the same
    // Darcy transport and root-uptake terms as the water solver at a saturated
    // surface; rainfall above this capacity must remain ponded or become runoff.
    static double maxSurfaceInfiltrationCpp(const soilpstruct& soilp, const soilwatermod& state,
        const envstruct& env, double dT, double pTAW)
    {
        double Evapmmhr = evaporationFluxCpp(soilp, env, state.theta[0], state.Tc[0], dT);
        double nPsi0 = 2.0 + 3.0 / soilp.b[0];
        double k0sat = soilp.Ksat[0];
        double kh1 = hydraulicConductivityFromThetaCpp(soilp, state.theta[1], 1);
        // A saturated surface holds no vapour pathway, so the face carries
        // half the vapour conductivity of the layer below.
        double kv1 = vaporConductivityFromPsiThetaCpp(soilp, state.psiw[1], state.theta[1], state.Tc[1] + 273.15, 1);
        double ff0max = (state.psiw[1] * kh1 - soilp.psie[0] * k0sat) / (state.dz[0] * (1.0 - nPsi0))
            + 0.5 * kv1 * (state.psiw[1] - soilp.psie[0]) / state.dz[0] - g * k0sat;
        std::vector<double> STr = transpirationDistributeCpp(soilp, state.rootfrac, env.Et, dT, state.psiw, pTAW);
        return Evapmmhr - dT * (ff0max - STr[0]);
    }

    // Advance the layered soil-water profile through one timestep by solving
    // the mixed liquid-plus-vapour water balance in water-potential space.
    // Vertical transport, storage, bare-soil evaporation and root uptake are
    // coupled nonlinearly, with either free drainage or a saturated lower
    // boundary controlling water loss at depth.
    soilwaterresult SoilWaterCpp(soilwatermod state, const soilpstruct& soilp, const envstruct& env,
        double dT, double pTAW, int maxNrIterations, double tolerance)
    {
        const double rho = 1000.0;
        const int n = soilp.nLayers;
        // Fixed start-of-timestep water and vapour stores provide the mass-
        // balance reference while the current state is iterated to convergence.
        const std::vector<double> oldtheta = state.oldtheta;
        const std::vector<double> oldvapor = state.oldvapor;
        std::vector<double> psiw = state.psiw;
        std::vector<double> theta = state.theta;
        std::vector<double> vapor = state.vapor;
        std::vector<double> k(n + 1, 0.0);
        std::vector<double> aa(n, 0.0), bb(n, 0.0), cc(n, 0.0), dd(n, 0.0);
        std::vector<double> ff(n, 0.0), Ca(n, 0.0), dpsi(n, 0.0);
        std::vector<double> kh(n + 1, 0.0), kv(n + 1, 0.0), dkh(n, 0.0), dkv(n, 0.0);
        std::vector<double> dFaceA(n, 0.0), dFaceB(n, 0.0);
        // Lower boundary: either a saturated fixed-head condition or, under
        // free drainage, continuation of the deepest resolved layer.
        if (!soilp.FreeDrain) {
            psiw[n] = soilp.psie[n - 1];
            theta[n] = soilp.thetaS[n - 1];
            k[n] = soilp.Ksat[n - 1];
        }
        // State and step of the last full Newton step, so that a step which
        // increases the residual can be retraced at half the length.
        std::vector<double> stepOrigin(n, 0.0), stepTaken(n, 0.0);
        double originResidual = 0.0, stepFraction = 1.0;
        bool haveStep = false;
        double surplusRate = 0.0;   // water a saturated surface cannot accept (kg m-2 s-1)
        int iter = 1;
        double massBalance = 1.0;
        double Evapmmhr = 0.0;
        while (iter < maxNrIterations) {
            // Net water flux at the surface combines rainfall reaching the
            // ground with evaporation from the current surface state.
            Evapmmhr = evaporationFluxCpp(soilp, env, theta[0], state.Tc[0], dT);
            // precipGround is throughfall after canopy interception, so the
            // soil surface receives only the water that actually reaches the ground.
            double surfaceFlux = (Evapmmhr - env.precipGround) / dT;
            // Rate at which that surface flux changes as the surface layer
            // wets. Adding it to the surface diagonal below solves the
            // evaporative boundary by the same Newton step the solver already
            // applies to its internal Darcy fluxes, rather than holding it at
            // the previous estimate. It cancels at the fixed point, so the
            // converged water profile is unchanged. The perturbation is a
            // fraction of the potential because psiw spans orders of magnitude
            // between air entry and the dry limit; the floor covers the
            // near-saturated end, where evaporation barely responds to wetness
            // and the derivative is near zero in any case.
            double dEvap_dpsi = 0.0;
            {
                constexpr double PSI_SLOPE_FRAC = 1e-4;
                constexpr double PSI_SLOPE_MIN = 1.0;
                double hp = PSI_SLOPE_FRAC * std::abs(psiw[0]);
                if (hp < PSI_SLOPE_MIN) hp = PSI_SLOPE_MIN;
                double theta_probe = thetaFromPsiCpp(soilp, psiw[0] + hp, 0);
                double Evprobe = evaporationFluxCpp(soilp, env, theta_probe, state.Tc[0], dT);
                dEvap_dpsi = ((Evprobe - Evapmmhr) / hp) / dT;
                if (!(dEvap_dpsi > 0.0)) dEvap_dpsi = 0.0;
            }
            // Distribute plant water uptake across roots using the current
            // layer water potentials.
            std::vector<double> STr = transpirationDistributeCpp(soilp, state.rootfrac, env.Et, dT, psiw, pTAW);

            // Liquid and vapour conductivities, their derivatives with respect
            // to water potential, and the derivative of stored water, for each
            // layer. Liquid conductivity, (theta/thetaS)^(2b+3) on the Campbell
            // retention curve, varies as psi^-(2+3/b) below air entry. Vapour
            // conductivity varies through the Kelvin humidity and through the
            // air-filled porosity.
            for (int i = 0; i < n; ++i) {
                double Tkelvin = state.Tc[i] + 273.15;
                kh[i] = hydraulicConductivityFromThetaCpp(soilp, theta[i], i);
                kv[i] = vaporConductivityFromPsiThetaCpp(soilp, psiw[i], theta[i], Tkelvin, i);
                k[i] = kh[i] + kv[i];
                double Cw = dThetaDPsiCpp(soilp, psiw[i], i);
                double Cv = dVaporDPsiCpp(soilp, psiw[i], theta[i], Tkelvin, i);
                Ca[i] = state.vol[i] * (rho * Cw + Cv) / dT;
                dkh[i] = -(2.0 + 3.0 / soilp.b[i]) * kh[i] / psiw[i];
                double airFilled = soilp.thetaS[i] - theta[i];
                dkv[i] = kv[i] * Mw / (RgasC * Tkelvin) - ((airFilled > 0.0) ? kv[i] / airFilled * Cw : 0.0);
            }
            // Lower boundary node: a saturated fixed head holds no vapour
            // pathway; free drainage continues the deepest layer.
            if (!soilp.FreeDrain) {
                psiw[n] = soilp.psie[n - 1];
                theta[n] = soilp.thetaS[n - 1];
                kh[n] = soilp.Ksat[n - 1];
                kv[n] = 0.0;
            }
            else {
                psiw[n] = psiw[n - 1];
                theta[n] = theta[n - 1];
                kh[n] = kh[n - 1];
                kv[n] = kv[n - 1];
            }
            k[n] = kh[n] + kv[n];
            // Upward flux across the face below each layer, and its derivatives
            // with respect to the potential of the node above (A) and below (B)
            // that face. Liquid water moves down the gradient of the matric flux
            // potential, psi*kh/(1 - n) with n = 2 + 3/b, which integrates the
            // Campbell conductivity exactly between the nodes, and falls under
            // gravity. Vapour diffuses down the gradient of water potential with
            // the mean of the two nodes' vapour conductivities.
            for (int i = 0; i < n; ++i) {
                double nPsi = 2.0 + 3.0 / soilp.b[i];
                double sh = state.dz[i] * (1.0 - nPsi);
                double dpsiFace = psiw[i + 1] - psiw[i];
                double kvFace = 0.5 * (kv[i] + kv[i + 1]);
                ff[i] = (psiw[i + 1] * kh[i + 1] - psiw[i] * kh[i]) / sh
                    + kvFace * dpsiFace / state.dz[i] - g * kh[i];
                double dPhi = kh[i] * (1.0 - nPsi);
                if (i == n - 1 && soilp.FreeDrain) {
                    // The node below moves with this one, so only gravity
                    // responds to this layer's potential.
                    dFaceA[i] = -g * dkh[i];
                    dFaceB[i] = 0.0;
                    continue;
                }
                double dPhiBelow = (i + 1 < n) ? kh[i + 1] * (1.0 - (2.0 + 3.0 / soilp.b[i + 1])) : 0.0;
                double dkvBelow = (i + 1 < n) ? dkv[i + 1] : 0.0;
                dFaceA[i] = -dPhi / sh - kvFace / state.dz[i] + 0.5 * dkv[i] * dpsiFace / state.dz[i] - g * dkh[i];
                dFaceB[i] = dPhiBelow / sh + kvFace / state.dz[i] + 0.5 * dkvBelow * dpsiFace / state.dz[i];
            }
            // Assemble each layer's mass-balance residual, gaining water from
            // the face below and losing it through the face above, and its
            // exact Jacobian.
            massBalance = 0.0;
            aa[0] = 0.0;
            cc[0] = -dFaceB[0];
            bb[0] = -dFaceA[0] + Ca[0] + dEvap_dpsi;
            dd[0] = surfaceFlux + STr[0] - ff[0]
                + state.vol[0] * (rho * (theta[0] - oldtheta[0]) + (vapor[0] - oldvapor[0])) / dT;
            massBalance += std::abs(dd[0]);

            for (int i = 1; i < n; ++i) {
                aa[i] = dFaceA[i - 1];
                cc[i] = -dFaceB[i];
                bb[i] = dFaceB[i - 1] - dFaceA[i] + Ca[i];
                dd[i] = ff[i - 1] + STr[i] - ff[i]
                    + state.vol[i] * (rho * (theta[i] - oldtheta[i]) + (vapor[i] - oldvapor[i])) / dT;
                massBalance += std::abs(dd[i]);
            }
            // A surface already at saturation can take no more water, however
            // much arrives: its potential is held there and the water it cannot
            // pass on to the layer below becomes surface excess, for the caller
            // to pond or shed. Without this the balance has no solution inside
            // the saturation bound and the excess would simply vanish.
            surplusRate = 0.0;
            if (psiw[0] >= soilp.psie[0] - 2e-8 && dd[0] < 0.0) {
                surplusRate = -dd[0];
                massBalance -= std::abs(dd[0]);
                dd[0] = 0.0;
                aa[0] = 0.0;
                bb[0] = 1.0;
                cc[0] = 0.0;
            }
            if (massBalance <= tolerance) break;

            // A full step that increased the residual is retraced from where
            // it started at half the length, down to 1/128 of it. Checking at
            // the next assembly costs nothing extra, since that assembly is
            // needed anyway.
            if (haveStep && massBalance > originResidual && stepFraction > 1.0 / 128.0) {
                stepFraction *= 0.5;
                for (int i = 0; i < n; ++i) {
                    psiw[i] = stepOrigin[i] - stepFraction * stepTaken[i];
                    psiw[i] = std::min(psiw[i], soilp.psie[i] - 1e-8);
                    psiw[i] = std::max(psiw[i], soilp.psi_min[i]);
                    theta[i] = thetaFromPsiCpp(soilp, psiw[i], i);
                    vapor[i] = vaporFromPsiCpp(soilp, psiw[i], theta[i], state.Tc[i] + 273.15, i);
                }
                ++iter;
                continue;
            }

            // Solve for the water-potential correction that removes the
            // profile-wide mass-balance residual, and take it in full, within
            // the physically allowed range between oven dryness and air
            // entry.
            utils::ThomasResult TBC = utils::thomasSolveCpp(aa, bb, cc, dd, dpsi, 0, n - 1);
            dpsi = TBC.x;
            for (int i = 0; i < n; ++i) {
                stepOrigin[i] = psiw[i];
                stepTaken[i] = dpsi[i];
                psiw[i] = psiw[i] - dpsi[i];
                psiw[i] = std::min(psiw[i], soilp.psie[i] - 1e-8);
                psiw[i] = std::max(psiw[i], soilp.psi_min[i]);

                theta[i] = thetaFromPsiCpp(soilp, psiw[i], i);
                vapor[i] = vaporFromPsiCpp(soilp, psiw[i], theta[i], state.Tc[i] + 273.15, i);
            }
            originResidual = massBalance;
            stepFraction = 1.0;
            haveStep = true;
            ++iter;
        }
        state.psiw = psiw;
        state.theta = theta;
        state.vapor = vapor;
        state.k = k;
        soilwaterresult out;
        out.state = state;
        out.success = (massBalance <= tolerance);
        out.iterations = iter;
        out.Evapmmhr = Evapmmhr;
        out.surplus = surplusRate * dT;
        return out;
    }

    // ========================================================================
    // Construct C++ model state from R-side parameter lists
    // ========================================================================
    // These adapters translate vegetation and soil parameters into the compact
    // structures used by the physical routines above. They also establish
    // derived quantities such as PAR optics, signed matric potentials, soil
    // layer geometry and root distribution before the timestep integration.

    // Vegetation state for the physical model. `Lfrac` is kept separate
    // because the fraction of plant area that is living foliage can change
    // seasonally while the remaining PFT parameters are structural traits.
    static vegpstruct toVegpstructCpp(Rcpp::List vegp, double Lfrac)
    {
        vegpstruct out;
        out.hgt = vegp["h"];
        out.pai = vegp["pai"];
        // A plant of zero height or zero plant area has no mass: either one
        // zero is bare ground, and both are set to zero so that bare ground is
        // the same surface whatever height or plant area was supplied with it.
        if (out.hgt <= 0.0 || out.pai <= 0.0) { out.hgt = 0.0; out.pai = 0.0; }
        out.x = vegp["x"];
        out.clump = vegp["clump"];
        out.Lfrac = Lfrac;
        out.len = vegp["len"];
        out.wid = vegp["wid"];
        out.lref = vegp["lref"];
        out.ltra = vegp["ltra"];
        // PAR reflectance/transmittance use the model's fixed quarter-of-
        // shortwave optical approximation.
        out.lrefp = 0.25 * out.lref;
        out.ltrap = 0.25 * out.ltra;
        // Longwave emissivity is shared by vegetation and ground and is
        // therefore not part of the vegetation parameter list.
        out.mwft = vegp["mwft"];
        out.Vcmax25 = vegp["Vcmx25"];
        out.Tup = vegp["Tup"];
        out.Tlw = vegp["Tlow"];
        out.Dcrit = vegp["Dcrit"];
        out.alpha = vegp["alpha"];
        out.f0 = vegp["f0"];
        out.fd = vegp["fd"];
        out.gsmaxCap = vegp["gsmaxCap"];
        // Minimum whole-plant hydraulic resistance is supplied directly
        // as the PFT-specific structural parameter rpmin.
        out.rpmin = vegp["rpmin"];
        out.psi50 = vegp["psi50"];
        out.apsi = vegp["apsi"];
        out.pTAW = vegp["pTAW"];
        out.root50 = vegp["root50"];
        out.root95 = vegp["root95"];
        if (!(out.root50 > 0.0) || !(out.root95 > out.root50)) {
            Rcpp::stop("vegp root50 must be positive and root95 greater than root50");
        }
        return out;
    }

    // Fixed physical properties of the layered soil profile.
    static soilpstruct toSoilpstructCpp(Rcpp::List soilc)
    {
        soilpstruct out;
        // Soil hydraulic and thermal properties are independent of terrain
        // slope/aspect; terrain affects the forcing rather than this profile.
        out.gref = Rcpp::as<double>(soilc["gref"]);
        out.grefPAR = Rcpp::as<double>(soilc["grefPAR"]);
        // Ground longwave emissivity is represented by the shared
        // surfaceEmissivity constant rather than a soil-specific field.
        out.nLayers = Rcpp::as<int>(soilc["nLayers"]);
        out.FreeDrain = Rcpp::as<bool>(soilc["FreeDrain"]);
        out.Vq = Rcpp::as<std::vector<double>>(soilc["Vq"]);
        out.Vm = Rcpp::as<std::vector<double>>(soilc["Vm"]);
        out.Vo = Rcpp::as<std::vector<double>>(soilc["Vo"]);
        out.Mc = Rcpp::as<std::vector<double>>(soilc["Mc"]);
        out.psie = Rcpp::as<std::vector<double>>(soilc["psi_e"]); // R element key is "psi_e"
        out.b = Rcpp::as<std::vector<double>>(soilc["b"]);
        out.thetaR = Rcpp::as<std::vector<double>>(soilc["Smin"]);
        out.thetaS = Rcpp::as<std::vector<double>>(soilc["Smax"]);
        out.Ksat = Rcpp::as<std::vector<double>>(soilc["Ksat"]);
        // Matric potentials are negative by convention; the R-side table
        // supplies magnitudes, so convert them to physical signed values.
        for (double& v : out.psie) v = -std::abs(v);
        // Lower bound on water potential for the water solve. Evaporation is
        // limited by the pore air's Kelvin humidity, which approaches the
        // atmosphere's own humidity well before the soil reaches oven dryness,
        // so the bound is a numerical backstop rather than a physical limit and
        // is placed where no soil in equilibrium with the air can go.
        out.psi_min.assign(out.psie.size(), mc::psiOvenDry);
        return out;
    }

    // Initialise the soil heat grid and state. Layer thickness increases with
    // depth so the near-surface thermal gradient is resolved most finely.
    static soilheatmod toSoilheatmodCpp(Rcpp::List soilc, std::vector<double> Te, std::vector<double> wc)
    {
        int n = Rcpp::as<int>(soilc["nLayers"]);
        double totalDepth = Rcpp::as<double>(soilc["totalDepth"]);
        std::vector<double> z = utils::geometricCpp(n, totalDepth);
        std::vector<double> dz(n + 1);
        for (int i = 0; i <= n; ++i) {
            dz[i] = z[i + 1] - z[i];
        }
        // Each node stores heat over the control volume around it: half of each
        // neighbouring layer, and half the top layer at the surface (Campbell
        // 1985, eqs 4.12 and 4.14). The deepest node is the fixed boundary.
        std::vector<double> vol(n + 1);
        vol[0] = 0.5 * (z[1] - z[0]);
        for (int i = 1; i <= n; ++i) vol[i] = 0.5 * (z[i + 1] - z[i - 1]);
        soilheatmod state;
        state.n = n;
        state.z = z;             // node depths; node i holds the state at z[i]
        state.dz = dz;           // spacing between node i and node i + 1
        state.vol = vol;         // control volume around each node
        state.wc = wc;           // volumetric water fraction at each node
        state.Te = Te;           // temperature at each node
        state.oldTe = Te;        // previous-timestep temperature at each node
        state.iters = 0;
        return state;
    }

    // Initialise soil water on the same vertical grid as soil heat, deriving
    // matric potential, hydraulic conductivity and pore-space vapour from the
    // starting water contents. `root50` and `root95` set the depth distribution
    // of root uptake.
    static soilwatermod toSoilwatermodCpp(const soilpstruct& soilp, const soilheatmod& heat, double root50, double root95)
    {
        soilwatermod out;
        out.n = heat.n;
        out.z = heat.z;
        out.dz = heat.dz;
        out.Tc = heat.Te;
        out.oldTc = heat.oldTe;
        out.theta = heat.wc;
        std::vector<double> vol(out.n + 1);
        std::vector<double> psiw(out.n + 1);
        std::vector<double> k(out.n + 1);
        std::vector<double> vapor(out.n + 1);
        // Node-centred control volumes extend halfway to neighbouring
        // nodes. The surface has only the downward half-cell, so it still
        // carries a finite water-storage capacity in the mass balance.
        vol[0] = (out.z[1] - out.z[0]) / 2.0;
        for (int i = 0; i <= out.n; ++i) {
            if (i > 0) vol[i] = (out.z[i + 1] - out.z[i - 1]) / 2.0;
            psiw[i] = waterPotentialCpp(soilp, out.theta[i], i);
            k[i] = hydraulicConductivityFromThetaCpp(soilp, out.theta[i], i);
            vapor[i] = vaporFromPsiCpp(soilp, psiw[i], out.theta[i], out.Tc[i] + 273.15, i);
        }
        out.vol = vol;
        out.psiw = psiw;
        out.vapor = vapor;
        out.k = k;
        out.oldvapor = vapor;
        out.oldtheta = out.theta;
        out.rootfrac = utils::rootDistributeCpp(out.z, root50, root95);
        return out;
    }

    // ========================================================================
    // Big Leaf timestep integration
    // ========================================================================
    // The point model represents the vegetation/ground system as a coupled
    // exchange surface above a layered soil profile. At each timestep the
    // solution links radiation, atmospheric stability and aerodynamic
    // resistance, stomatal/soil water loss, canopy interception, soil heat
    // and water transport, and the bulk surface temperature. Soil temperature,
    // soil water and intercepted canopy water persist between timesteps; the
    // atmospheric/surface energy-balance iteration is solved afresh each hour.

    // Time series returned by the integrated point model. The core outputs
    // describe radiation, turbulent exchange, surface/soil temperature and
    // water state; additional fields expose quantities needed by the grid
    // model to scale the reference solution across a landscape.
    struct BigLeafOutput {
        // The emission coefficients of the soil surface and of the whole
        // canopy+ground surface (see envstruct::emGround). Returned so either
        // energy balance can be reconstructed from the time series, and so the
        // grid model can evaluate this run as a reference on its own terms.
        std::vector<double> emGround;
        std::vector<double> emCanopy;
        // Surface energy balance and turbulent state
        std::vector<double> RabsGround;  // total radiation absorbed by the ground (W/m^2)
        std::vector<double> RabsCanopy;  // total radiation absorbed by canopy + ground (W/m^2)
        std::vector<double> uf;          // friction velocity (m/s)
        std::vector<double> LL;          // Monin-Obukhov length (m)
        std::vector<double> Hout;        // sensible heat flux from the bulk surface (W/m^2)
        std::vector<double> Tground;     // ground-surface temperature (deg C)
        std::vector<double> G;           // ground heat flux, positive into the soil (W/m^2)
        std::vector<double> Tground0;    // ground temperature implied by the same balance if G = 0 (deg C)
        std::vector<double> Tcanopy;     // bulk canopy + ground exchange-surface temperature (deg C)

        // Plant water state and solver diagnostics
        std::vector<double> psi_r;       // root-weighted soil water potential (MPa)
        std::vector<double> iters;       // outer surface-energy-balance iterations
        std::vector<double> witers;      // soil-water iterations on the final outer pass

        // Above-canopy profiles used either for a requested output height or
        // for moving meteorological forcing between measurement heights.
        std::vector<double> Tabove;      // air temperature at zAbove (deg C)
        std::vector<double> windAbove;   // wind speed at zAbove (m/s)
        std::vector<double> RHabove;     // relative humidity at zAbove (%)

        // Requested below-ground state, populated only when reqhgt <= 0. The
        // two fields differ in how their own deep boundary behaves: Tzbelow's
        // deep clamp is the fixed value seeded at the start of the run (see
        // SoilHeatCpp's own boundary-condition note), while thetazbelow's is
        // recomputed every timestep from soilp.FreeDrain.
        std::vector<double> Tzbelow;     // soil temperature at reqhgt (deg C)
        std::vector<double> thetazbelow; // volumetric soil water content at reqhgt

        // Requested above-canopy state, populated when reqhgt is above canopy.
        std::vector<double> Tzabove;     // air temperature at reqhgt (deg C)
        std::vector<double> windzabove;  // wind speed at reqhgt (m/s)
        std::vector<double> RHzabove;    // relative humidity at reqhgt (%)

        // Canopy-top and solar geometry retained for landscape-scale scaling.
        std::vector<double> uh;          // wind speed at canopy top (m/s), on the exchange scale
        std::vector<double> windScale;   // factor turning exchange-scale winds into actual winds
        std::vector<double> zenr;        // solar zenith angle (radians)
        std::vector<double> zend;        // solar zenith angle (degrees)
        std::vector<double> azid;        // solar azimuth (degrees)

        // Soil/surface state used to anchor landscape-scale spatial scaling.
        std::vector<double> D_D;         // diurnal thermal damping depth of the surface soil (m)
        std::vector<double> rHa;         // aerodynamic resistance from ground to reference atmosphere (s/m)
        std::vector<double> rBL;         // aerodynamic resistance from the canopy+ground exchange surface to the
                                         // reference atmosphere (s/m). This is the resistance the big leaf's own
                                         // sensible and latent fluxes cross, and the one `rSurf` is defined
                                         // against; `rHa` above starts at the ground instead and is the soil
                                         // balance's. The two are not interchangeable.
        std::vector<double> theta0;      // volumetric water content of the surface soil layer
        std::vector<std::vector<double>> TeProfile; // complete soil-temperature profile (deg C);
                                          // diagnostic only, kept for validating the grid model's
                                          // below-ground shortcut against this real multilayer solve
        std::vector<double> swaterdepth; // intercepted water stored on the canopy (mm)
        std::vector<double> rSurf;       // additional canopy surface/vapour resistance beyond aerodynamic resistance (s/m)
        std::vector<double> hSurf;       // factor on the bulk surface's saturation vapour pressure; can exceed 1
                                         // where the soil is warmer than the exchange surface
        std::vector<double> filmEvap;    // evaporation of intercepted water (mm)
        std::vector<double> filmAvailable; // intercepted water the film could give up in the hour (mm)
        std::vector<double> precipGround; // rain reaching the ground: throughfall plus drip (mm)
        std::vector<double> wetShare;    // share of the canopy surface wet at the start of the hour
        std::vector<double> groundhr;    // relative humidity of air in equilibrium with the soil surface (0-1)
        // Soil state at the end of the run, one value per soil node
        std::vector<double> TeEnd;       // temperature (deg C)
        std::vector<double> thetaEnd;    // volumetric water content
    };

    // Integrate the coupled point model through the supplied meteorological
    // time series. For each timestep:
    //   1. compute solar geometry and shortwave/PAR absorption;
    //   2. iterate longwave exchange, atmospheric stability/aerodynamic
    //      resistance, canopy water/stomatal fluxes, soil heat and soil water;
    //   3. close the bulk surface energy balance for canopy temperature and
    //      sensible heat flux;
    //   4. commit the converged soil/canopy water states and derive requested
    //      above- or below-ground diagnostics.
    // The outer fixed-point iteration is necessary because surface temperature,
    // sensible heat flux and atmospheric stability feed back on one another,
    // while soil moisture also feeds plant water status and evaporative cooling.
    // `vegp` is taken by non-const reference because leafgsCpp lazily caches a
    // derived hydraulic parameter (apsi) on it the first time it is needed,
    // rather than recomputing it every timestep. Atmospheric CO2 (Ca) is a
    // fixed placeholder for the whole run, not a real per-timestep driver the
    // way root-zone water potential (psi_r) genuinely is re-solved each hour.
    BigLeafOutput BigLeafCoreCpp(
        const std::vector<double>& year, const std::vector<double>& month,
        const std::vector<double>& day, const std::vector<double>& hour,
        const std::vector<double>& temp, const std::vector<double>& relhum,
        const std::vector<double>& pres, const std::vector<double>& swdown,
        const std::vector<double>& difrad, const std::vector<double>& lwdown,
        const std::vector<double>& windspeed, const std::vector<double>& precip,
        vegpstruct& vegp, const soilpstruct& soilp,
        soilheatmod soilheat, soilwatermod soilwater,
        double gref, double grefPAR,
        double lat, double lon,
        double zref, double reqhgt, double Ca, double zAbove,
        int maxIter = 100, double tolerance = 1e-2, bool pooling = false,
        double slope = 0.0, double aspect = 0.0, double svfa = 1.0,
        const std::vector<double>& shelterc = std::vector<double>{1.0})
    {
        int n = static_cast<int>(temp.size());
        // Per-timestep wind shelter coefficient -- real topographic wind
        // shelter depends on wind direction, which varies hourly, unlike
        // slope/aspect/sky-view factor, which are fixed terrain properties
        // of the location. Broadcasts a length-1 default (no shelter) to
        // every timestep when a real per-timestep series isn't supplied.
        bool shelterSpatial = (static_cast<int>(shelterc.size()) == n);
        BigLeafOutput out;
        out.rBL.resize(n);
        out.emGround.resize(n);
        out.emCanopy.resize(n);
        out.RabsGround.resize(n);
        out.RabsCanopy.resize(n);
        out.uf.resize(n);
        out.LL.resize(n);
        out.Hout.resize(n);
        out.Tground.resize(n);
        out.G.resize(n);
        out.Tground0.resize(n);
        out.Tcanopy.resize(n);
        out.psi_r.resize(n);
        out.iters.resize(n);
        double nanv = std::numeric_limits<double>::quiet_NaN();
        out.Tzbelow.assign(n, nanv);
        out.thetazbelow.assign(n, nanv);
        out.Tabove.assign(n, nanv);
        out.windAbove.assign(n, nanv);
        out.RHabove.assign(n, nanv);
        out.Tzabove.assign(n, nanv);
        out.windzabove.assign(n, nanv);
        out.RHzabove.assign(n, nanv);
        out.D_D.resize(n);
        out.rHa.resize(n);
        out.rSurf.resize(n);
        out.hSurf.resize(n);
        out.filmEvap.resize(n);
        out.filmAvailable.resize(n);
        out.precipGround.resize(n);
        out.wetShare.resize(n);
        out.theta0.resize(n);
        out.TeProfile.resize(n);
        out.witers.resize(n);
        out.groundhr.resize(n);
        out.swaterdepth.resize(n);
        out.uh.resize(n);
        out.windScale.resize(n);
        out.zenr.resize(n);
        out.zend.resize(n);
        out.azid.resize(n);

        // Canopy optical properties are fixed for the run; solar geometry and
        // incoming radiation are applied separately at each timestep.
        radsetup s = RadswabsSetupCpp(vegp, gref, grefPAR);
        // Aerodynamic displacement is fixed by canopy structure for this run.
        double d = utils::zeroplanedisCpp(vegp.hgt, vegp.pai);
        // The transport column's constants are fixed by canopy height, plant
        // area and that displacement height, none of which change during a
        // run, so it is built once here rather than on every convergence
        // iteration of every timestep. It exists at every canopy height above
        // zero; its resistances vanish continuously as the canopy does.
        const bool runHasColumn = vegp.hgt > 0.0;
        utils::ColumnStruct runCol{};
        if (runHasColumn) runCol = utils::columnSetupCpp(vegp.hgt, vegp.pai, d);
        // Weight of the ground's own balance in the soil's surface row (the
        // matching constraint, see matchResidualCpp): 0 at and below the soil's
        // heat roughness height, where ground and foliage are one surface, 1 at
        // and above 10 mm, where the ground's balance is used unaltered.
        const double z0hSoil = 0.2 * mc::bareGroundZ0;
        const double wMatch = std::min(1.0, std::max(0.0, (vegp.hgt - z0hSoil) / (0.010 - z0hSoil)));
        // Intercepted canopy water is a prognostic state carried between
        // timesteps. During the within-timestep energy-balance iteration the
        // start-of-hour storage is held fixed, so repeated evaluations solve
        // the same physical hour rather than repeatedly advancing the water
        // balance. The final candidate is committed once convergence is reached.
        double swaterdepth = 0.0;
        // Water that exceeds the soil's infiltration capacity is either stored
        // at the surface for later infiltration (pooling = true) or treated as
        // runoff (pooling = false).
        double surfacePool = 0.0;
        // Sensible and latent heat leaving the whole surface, carried between
        // iterations and hours to set the air the ground exchanges with.
        double Hsurf = 0.0, Lsurf = 0.0;

        for (int i = 0; i < n; ++i) {
            envstruct env{};
            env.tair = temp[i];
            env.rh = relhum[i];
            env.pk = pres[i];
            env.Ca = Ca;
            env.precip = precip[i];
            env.Rsw = swdown[i];
            env.Rdif = difrad[i];
            env.Rlw = lwdown[i];
            env.uref = windspeed[i];
            env.tcanopy = temp[i]; // initial guess for this timestep's iteration

            // Solar position is computed once per timestep and shared
            // between RadswabsStepCpp and the stomatal setup below. Full
            // solpositionCpp (zenith + azimuth), not the zenith-only
            // shortcut used elsewhere, since this point model can be run
            // for an inclined slope (slope/aspect default to flat), and
            // azimuth is needed for the direct-beam weighting on a tilted
            // surface.
            solmodel solp = solpositionCpp(lat, lon, static_cast<int>(year[i]), static_cast<int>(month[i]),
                static_cast<int>(day[i]), hour[i]);
            {
                RadswabsStepCpp(s, vegp, gref, grefPAR, slope, aspect, svfa, solp, env);
            }

            // Freeze the soil state at the start of the hour. The heat and
            // water solvers iterate away from this fixed initial state, so the
            // outer surface-energy-balance loop cannot accidentally advance
            // the soil through several fictitious timesteps.
            soilheat.oldTe = soilheat.Te;
            soilwater.oldtheta = soilwater.theta;
            soilwater.oldvapor = soilwater.vapor;
            env.psi_r = rootzonePsiCpp(soilp, soilwater.theta, soilwater.rootfrac);

            // Hold canopy and surface-pool storage at their start-of-hour
            // values while the coupled fluxes are iterated, only committing
            // once the outer iteration converges (a few lines below).
            // canopyWaterBudgetCpp itself does plain forward-Euler integration
            // over a full timestep on each call, unlike the soil heat/water
            // solvers, which manage their own old/new-state distinction
            // internally -- so this anchor-and-commit-once discipline has to
            // be imposed here, externally, rather than inside that function.
            double swaterdepthOld = swaterdepth;
            double swaterdepthCandidate = swaterdepth;
            double surfacePoolOld = surfacePool;
            double surfacePoolCandidate = surfacePool;

            // Sunlit/shaded leaf fractions and absorbed PAR are fixed by the
            // current radiation field, so they are established before solving
            // the temperature-water-stability feedbacks below.
            stomatalsetup stomSetup;
            if (s.vegetated) stomSetup = StomatalSetupCpp(solp, env, vegp);

            utils::Aitken1DState st_G;
            utils::Aitken1DState st_tcanopy;
            env.stabPassDone = false;
            env.stabFlipped = false;
            // The ground heat flux is relaxed between passes. Its self-tuning
            // Aitken step can stall after a reversal and then jump, which in light
            // wind feeds a limit cycle with the canopy balance; there a fixed
            // half-way step is taken after the first pass instead.
            const bool gFixedStep = env.uref < mc::dampWindSpeed;
            auto relaxG = [&](double oldv, double newv) {
                if (!gFixedStep) return utils::aitken1d(oldv, newv, st_G);
                if (!st_G.have_prev) { st_G.have_prev = true; return newv; }
                return 0.5 * (oldv + newv);
            };
            double G = 0.0;
            double err = 1e99;
            int iter = 0;
            int witers = 0;
            // Retain the converged canopy aerodynamic and vapour resistances
            // for any above-canopy profile requested after the coupled solve.
            double rHa_conv = 0.0, rV_conv = 0.0, hSurf_conv = 1.0;
            double filmEvap_conv = 0.0, filmAvail_conv = 0.0, precipGround_conv = 0.0, wetShare_conv = 0.0;
            envstruct envG = env;
            // Soil surface temperature from the previous pass, for the
            // convergence test: the ground lags the canopy by one pass.
            double tgroundPrev = soilheat.Te[0];

            // At least two passes, so that stability has been evaluated from a
            // non-zero sensible heat flux before the solution is accepted.
            while ((err > tolerance || iter < 2) && iter < maxIter) {
                // Longwave depends on the current env.tcanopy (the canopy's
                // own emission), so it's recomputed every iteration, unlike
                // the shortwave partition above which doesn't depend on
                // tcanopy or wind at all.
                RadlwabsStepCpp(s, svfa, env);

                // The sensible heat this pass's stability is formed from.
                const double Hused = env.H;
                windmodelCpp(vegp, zref, d, env, shelterSpatial ? shelterc[i] : shelterc[0]);
                double rRise = 0.0;
                double rGz = groundaeroresistCpp(vegp, zref, d, env, runCol, runHasColumn, &rRise); // ground to zref, for the soil model
                // Canopy/heat-exchange surface to zref; on bare ground the exchange
                // surface is the soil surface and the two resistances are one.
                double rHa = s.vegetated ? bulkaeroresistCpp(vegp, zref, d, env, runCol, runHasColumn) : rGz;

                // Soil drying lowers the vapour pressure at the ground surface,
                // through the Kelvin-equilibrium humidity of its pore air. Over
                // bare ground the soil surface is the exchange surface.
                double hr = std::max(soilrelhumCpp(soilp, env.tcanopy, soilwater.theta[0]), 1e-6);

                // Vegetated: the vapour network, with the soil's source at its
                // own temperature and moisture from the previous pass, and the
                // node at whichever of the mean source height and the exchange
                // surface is nearer the reference height.
                double eg = -1.0, rshV = 0.0;
                if (s.vegetated) {
                    double hg = std::max(soilrelhumCpp(soilp, soilheat.Te[0], soilwater.theta[0]), 1e-6);
                    eg = hg * utils::satvapCpp(soilheat.Te[0]);
                    rshV = (rRise > 0.0) ? std::min(rRise, rHa) : 0.0;
                }
                canopywaterresult cw = canopyWaterBudgetCpp(vegp, env, stomSetup, s.vegetated,
                    rHa, hr, swaterdepthOld, 3600.0, eg, rshV, rGz);
                double rV = rHa + cw.rSurf;
                double hSurf = cw.hSurf;
                swaterdepthCandidate = cw.swaterdepth;
                rHa_conv = rHa;
                rV_conv = rV;
                hSurf_conv = hSurf;
                filmEvap_conv = cw.filmEvap;
                filmAvail_conv = cw.filmAvailable;
                precipGround_conv = cw.precipGround;
                wetShare_conv = cw.wetShare;

                env.rHa = rGz;
                env.Et = cw.Et;
                env.precipGround = cw.precipGround;

                // The air the soil exchanges with. The canopy's heat and vapour
                // enter the column above the ground, so the ground does not see
                // reference air. With the canopy's sources spread uniformly with
                // height, the ground balance is one with air at T* across the
                // resistance from the ground to the level of that air, where T*
                // is the reference value raised by the whole surface's sensible
                // flux across the rest of the path. That level is the mean source
                // height, or the exchange surface where it is nearer the reference
                // (dense short canopies): a node below the exchange surface would
                // put the air the ground sees above the surface heating it. It is
                // the vapour network's node. The fluxes are the previous
                // iterate's, which converge with the canopy temperature. No
                // column, no rise: the ground then sees reference air across rGz.
                envG = env;
                envG.rVg = -1.0;
                double rShH = (rRise > 0.0) ? std::min(rRise, rHa) : rRise;
                if (rRise > 0.0) {
                    double phA = utils::phairCpp(env.tair, env.pk);
                    double cpA = utils::cpairCpp(env.tair);
                    double laA = utils::latentHeatCpp(env.tair);
                    double eaA = utils::satvapCpp(env.tair) * env.rh / 100.0;
                    double tStar = env.tair + Hsurf * rShH / (phA * cpA);
                    double eStar = eaA + Lsurf * rShH * env.pk / (laA * phA);
                    if (eStar < 1e-6) eStar = 1e-6;
                    envG.tair = tStar;
                    envG.rh = 100.0 * eStar / utils::satvapCpp(tStar);
                    envG.rHa = rGz - rShH;
                }
                if (s.vegetated) {
                    // The soil's vapour boundary is the network seen from the
                    // soil: its Thevenin air and resistance.
                    double eTh = std::max(cw.eTh, 1e-6);
                    envG.rh = 100.0 * eTh / utils::satvapCpp(envG.tair);
                    envG.rVg = cw.rTh;
                }
                // Limit rainfall entering the soil to what the surface layer can
                // accept without exceeding saturation. Any excess becomes
                // temporary ponded water or runoff according to `pooling`.
                {
                    double availablePrecip = env.precipGround + surfacePoolOld;
                    double infiltCapacity = maxSurfaceInfiltrationCpp(soilp, soilwater, envG, 3600.0, vegp.pTAW);
                    if (infiltCapacity < 0.0) infiltCapacity = 0.0;
                    if (availablePrecip > infiltCapacity) {
                        env.precipGround = infiltCapacity;
                        surfacePoolCandidate = pooling ? (availablePrecip - infiltCapacity) : 0.0;
                    } else {
                        env.precipGround = availablePrecip;
                        surfacePoolCandidate = 0.0;
                    }
                }
                envG.precipGround = env.precipGround;
                // Couple soil heat and water: moisture controls thermal
                // conductivity, heat capacity and surface humidity; the updated
                // temperature profile in turn controls water/vapour transport.
                soilheat.wc = soilwater.theta;
                // Below 10 mm the soil's surface row blends the ground's balance
                // toward the whole surface's, the equation the canopy step
                // solves, written in the row's form: the combined absorbed
                // radiation and emission, the big leaf's path, the network's
                // conductances with foliage and soil at one temperature.
                const bool matchHere = s.vegetated && wMatch < 1.0;
                envstruct envS = env;
                if (matchHere) {
                    envS.RabsGround = env.RabsCanopy;
                    envS.emGround = env.emCanopy;
                    envS.rHa = rHa;
                    envS.rVg = rV;
                    envS.wholeA = cw.aNet;
                    envS.wholeB = cw.bNet;
                }
                soilheat = SoilHeatCpp(soilheat, soilp, envG, 3600.0, 1.0, maxIter, 1e-2,
                    matchHere ? &envS : nullptr, wMatch);
                soilwater.Tc = soilheat.Te;
                soilwater.oldTc = soilheat.oldTe;
                soilwaterresult swres = SoilWaterCpp(soilwater, soilp, envG, 3600.0, vegp.pTAW, maxIter);
                soilwater = swres.state;
                witers = swres.iterations;
                // Water the saturated surface could not take after all, beyond
                // what the capacity estimate above withheld, joins the same
                // pool or runs off with it.
                if (pooling) surfacePoolCandidate += swres.surplus;
                env.psi_r = rootzonePsiCpp(soilp, soilwater.theta, soilwater.rootfrac);

                // Ground heat flux is the residual of the instantaneous soil-
                // surface energy balance (under the matching constraint, of the
                // blended balance the row solved). It therefore represents the
                // energy being conducted into or released from the soil at the
                // surface, including the phase lag created by diffusion through
                // the profile.
                double gWeightMatch = 1.0;
                if (matchHere) {
                    matchresidual mG = matchResidualCpp(soilp, envG, &envS, wMatch, soilheat.Te[0], soilwater.theta[0], soilheat.condSlope);
                    gWeightMatch = (wMatch > 0.0) ? wMatch / (wMatch + (1.0 - wMatch) * mG.kappa) : 0.0;
                    G = relaxG(G, mG.q);
                } else {
                    // On bare ground nothing within the pass reads G, so it needs
                    // no relaxation.
                    double Gnew = soilsurfaceEBCpp(soilp, envG, soilheat.Te[0], soilwater.theta[0]);
                    G = s.vegetated ? relaxG(G, Gnew) : Gnew;
                }

                double tcanopy_old = env.tcanopy;
                // The ground heat flux's response to the canopy temperature
                // through the air the soil sees: raising Tc raises T* by
                // rShH/rHa, the soil surface follows by the share
                // (rho cp / leg) / (C' + lamS), and conduction takes C' of that.
                // Carried into the canopy step as a slope, it turns the one-pass
                // lag between the canopy and the soil into a Newton step; the
                // fixed point is unchanged. Under the matching constraint the
                // ground's balance carries weight w/(w + (1 - w) kappa) of it.
                double gSlope = 0.0;
                if (s.vegetated && rRise > 0.0) {
                    double rhocp = utils::phairCpp(env.tair, env.pk) * utils::cpairCpp(env.tair);
                    double Cp = soilheat.condSlope, lS = soilheat.lamS, Rbar = rGz - rShH;
                    if (Cp > 0.0 && Rbar > 0.0 && Cp + lS > 0.0)
                        gSlope = Cp * (rhocp / Rbar) / (Cp + lS) * (rShH / rHa);
                }
                if (matchHere) gSlope *= gWeightMatch;
                // On bare ground the exchange surface is the soil surface, which
                // the soil's surface row has just solved; there is nothing more
                // to solve.
                double tcanopy_raw = soilheat.Te[0];
                if (s.vegetated) {
                    tcanopy_raw = utils::penmanMonteithCpp(env.RabsCanopy, env.tair, env.pk, env.rh, env.emCanopy,
                        rHa, rV, env.tcanopy, G, hSurf, cw.eV, gSlope);
                    // Canopy temperature feeds back on sensible heat, stability and
                    // aerodynamic resistance; Aitken relaxation damps oscillation of
                    // this nonlinear feedback while preserving the same fixed point.
                    env.tcanopy = utils::aitken1d(tcanopy_old, tcanopy_raw, st_tcanopy);
                } else {
                    env.tcanopy = tcanopy_raw;
                }

                // Latent heat leaving the whole surface, from the same
                // linearisation the Penman-Monteith solve just balanced. Read
                // only by the vegetated surface's network and humidity factor.
                if (s.vegetated) {
                    double phL = utils::phairCpp(env.tair, env.pk);
                    double Te = (tcanopy_old + env.tair) / 2.0;
                    double eaL = s.vegetated ? cw.eV : utils::satvapCpp(env.tair) * env.rh / 100.0;
                    double Da = hSurf * utils::satvapCpp(env.tair) - eaL;
                    double De = hSurf * (utils::satvapCpp(Te + 0.5) - utils::satvapCpp(Te - 0.5));
                    Lsurf =(utils::latentHeatCpp(tcanopy_old) * phL / (env.pk * rV)) * (Da + De * (env.tcanopy - env.tair));
                }

                double ph = utils::phairCpp(env.tair, env.pk);
                double cp = utils::cpairCpp(env.tair);
                env.H = (ph * cp / rHa) * (env.tcanopy - env.tair);
                Hsurf = env.H;

                // Convergence is judged on the unrelaxed residual of the canopy
                // balance, not on the relaxed step, which can be a small fraction
                // of it; and on the change in the soil surface over the pass. On
                // bare ground there is no canopy balance, and the remaining
                // coupling is stability: its residual is the difference between
                // the sensible heat stability was formed from and the sensible
                // heat that resulted, as a surface temperature.
                double groundStep = std::abs(soilheat.Te[0] - tgroundPrev);
                tgroundPrev = soilheat.Te[0];
                double stepC = s.vegetated ? std::abs(tcanopy_raw - tcanopy_old)
                                           : std::abs(env.H - Hused) * rHa / (ph * cp);
                err = std::max(stepC, groundStep);
                ++iter;
            }

            // The reported humidity factor for a vegetated surface is the one
            // that reproduces the solved whole-surface latent flux at the solved
            // surface temperature, k (hSurf es(Tc) - ea)/(rHa + rSurf) = Lsurf,
            // so that consumers of the pair (the within-canopy profile of
            // runpointmodelasgrid) see the flux the model balanced. It can exceed
            // one where the soil surface is warmer than the exchange surface.
            if (s.vegetated && rV_conv > 0.0) {
                double kL = utils::latentHeatCpp(env.tcanopy) * utils::phairCpp(env.tair, env.pk) / env.pk;
                double eaR = utils::satvapCpp(env.tair) * env.rh / 100.0;
                hSurf_conv = (eaR + Lsurf * rV_conv / kL) / utils::satvapCpp(env.tcanopy);
            }

            // Commit the converged water stores once for this real timestep.
            swaterdepth = swaterdepthCandidate;
            surfacePool = surfacePoolCandidate;

            // Zero-ground-heat-flux surface temperature: a diagnostic used by
            // the landscape shortcut to separate local radiative/aerodynamic
            // forcing from the soil's thermal-storage response. It takes reference
            // air across the full column, as the shortcut's own zero-flux pass
            // does, so that the two stay comparable.
            double Tground0 = groundTemp0Cpp(soilp, env, soilwater.theta[0], soilheat.Te[0]);
            out.rHa[i] = env.rHa;
            out.rBL[i] = rHa_conv;
            out.rSurf[i] = rV_conv - rHa_conv;
            out.hSurf[i] = hSurf_conv;
            out.filmEvap[i] = filmEvap_conv;
            out.filmAvailable[i] = filmAvail_conv;
            out.precipGround[i] = precipGround_conv;
            out.wetShare[i] = wetShare_conv;
            out.theta0[i] = soilwater.theta[0];
            out.TeProfile[i] = soilheat.Te;
            out.swaterdepth[i] = swaterdepth;
            // Surface-soil thermal damping depth for scaling the ground heat flux.
            out.D_D[i] = diurnalDampingDepthCpp(soilp, soilwater.theta[0], soilheat.Te[0], env.pk);

            out.emGround[i] = env.emGround;
            out.emCanopy[i] = env.emCanopy;
            out.RabsGround[i] = env.RabsGround;
            out.RabsCanopy[i] = env.RabsCanopy;
            out.uf[i] = env.uf;
            out.LL[i] = env.LL;
            out.Hout[i] = env.H;
            out.Tground[i] = soilheat.Te[0];
            out.G[i] = G;
            out.Tground0[i] = Tground0;
            out.Tcanopy[i] = env.tcanopy;
            out.psi_r[i] = env.psi_r;
            out.iters[i] = iter;
            out.witers[i] = witers;
            // Humidity in equilibrium with the actual ground temperature and
            // surface soil water content, used as the lower boundary for
            // below-canopy vapour transport.
            out.groundhr[i] = std::max(soilrelhumCpp(soilp, soilheat.Te[0], soilwater.theta[0]), 1e-6);
            out.uh[i] = env.uh;
            out.windScale[i] = env.windScale;
            out.zenr[i] = solp.zenr;
            out.zend[i] = solp.zend;
            out.azid[i] = solp.azid;

            // Resolve a requested height directly only where the point model
            // has an explicit profile: within the soil, or in the atmosphere
            // above the canopy. A height inside the canopy is a
            // per-cell quantity, not one this flat reference run generalises,
            // so it belongs to the grid model's own per-cell routines
            // (gridmodel::belowCanopyProfileGridCpp, and
            // belowCanopyProfilePointCpp for runpointmodelasgrid) rather than
            // here. Above a
            // vegetated surface the vapour profile is carried by the converged
            // total latent flux itself (see aboveCanopyProfileCpp).
            const double LEprof = s.vegetated ? Lsurf : std::numeric_limits<double>::quiet_NaN();
            if (!std::isnan(reqhgt) && reqhgt <= 0.0) {
                // Interpolate the resolved soil temperature and moisture
                // profiles to the requested depth.
                out.Tzbelow[i] = interpSoilProfileCpp(soilheat.z, soilheat.Te, reqhgt);
                out.thetazbelow[i] = interpSoilProfileCpp(soilwater.z, soilwater.theta, reqhgt);
            } else if (!std::isnan(reqhgt) && reqhgt >= vegp.hgt) {
                // Extend the converged surface state upward: wind by
                // Monin-Obukhov similarity, heat and moisture through the
                // transport column's resistances.
                AboveCanopyPoint ptR = aboveCanopyProfileCpp(vegp, reqhgt, d, zref, env.tcanopy, env,
                    rHa_conv, rV_conv, hSurf_conv, runCol, runHasColumn, LEprof);
                out.Tzabove[i] = ptR.Tz;
                out.windzabove[i] = ptR.windz;
                out.RHzabove[i] = ptR.RHz;
            }

            // A second above-canopy height can be evaluated from the same
            // converged profile when meteorological forcing must be translated
            // between measurement heights.
            if (!std::isnan(zAbove)) {
                AboveCanopyPoint pt = aboveCanopyProfileCpp(vegp, zAbove, d, zref, env.tcanopy, env,
                    rHa_conv, rV_conv, hSurf_conv, runCol, runHasColumn, LEprof);
                out.Tabove[i] = pt.Tz;
                out.windAbove[i] = pt.windz;
                out.RHabove[i] = pt.RHz;
            }
        }
        out.TeEnd = soilheat.Te;
        out.thetaEnd = soilwater.theta;
        return out;
    }

    // R-facing wrapper for the point model. It converts R vegetation/soil
    // parameter lists into the C++ state structures, seeds the soil profile,
    // runs the coupled time series, and returns the principal surface, soil
    // and profile diagnostics. The default soil initial condition is a
    // vertically uniform temperature (`matemp`, or the mean supplied air
    // temperature when absent) and water content midway between residual and
    // saturated values; callers should therefore allow an appropriate spin-up
    // before interpreting early-timestep soil state.
    //
    // `reqhgt` requests soil state at/below the surface or atmospheric state
    // above the canopy. `zAbove` is a separate above-canopy diagnostic used
    // when meteorological variables need to be translated between heights.
    // `Tinit` and `thetainit`, when given one value per soil node, replace the
    // default starting temperature and water profiles. The deepest node is the
    // fixed-temperature lower boundary for the whole run, so it starts at the
    // deep temperature (`matemp`) whatever `Tinit` holds. The soil state at the
    // end of the run is attached to the result as attributes TeEnd and thetaEnd.
    Rcpp::DataFrame BigLeafCpp(
        Rcpp::NumericVector year, Rcpp::NumericVector month,
        Rcpp::NumericVector day, Rcpp::NumericVector hour,
        Rcpp::NumericVector temp, Rcpp::NumericVector relhum, Rcpp::NumericVector pres,
        Rcpp::NumericVector swdown, Rcpp::NumericVector difrad, Rcpp::NumericVector lwdown,
        Rcpp::NumericVector windspeed, Rcpp::NumericVector winddir, Rcpp::NumericVector precip,
        Rcpp::List vegp, Rcpp::List soilc, double Lfrac, double lat, double lon,
        double zref, double reqhgt = NA_REAL, double Ca = 400.0, double matemp = NA_REAL,
        double zAbove = NA_REAL, int maxIter = 100, double tolerance = 1e-2, bool pooling = false,
        double slope = 0.0, double aspect = 0.0, double svfa = 1.0,
        Rcpp::NumericVector shelterc = Rcpp::NumericVector::create(1.0),
        Rcpp::NumericVector Tinit = Rcpp::NumericVector(0),
        Rcpp::NumericVector thetainit = Rcpp::NumericVector(0))
    {
        (void)winddir;

        vegpstruct vp = toVegpstructCpp(vegp, Lfrac);
        double gref = Rcpp::as<double>(soilc["gref"]);
        double grefPAR = Rcpp::as<double>(soilc["grefPAR"]);

        soilpstruct sp = toSoilpstructCpp(soilc);
        double Te0seed = matemp;
        if (std::isnan(Te0seed)) {
            std::vector<double> tempVec = Rcpp::as<std::vector<double>>(temp);
            Te0seed = std::accumulate(tempVec.begin(), tempVec.end(), 0.0) /
                static_cast<double>(tempVec.size());
        }
        std::vector<double> Te0(sp.nLayers + 1, Te0seed);
        std::vector<double> wc0(sp.nLayers + 1);
        for (int i = 0; i <= sp.nLayers; ++i) {
            wc0[i] = 0.5 * (sp.thetaR[i] + sp.thetaS[i]);
        }
        const R_xlen_t nNodes = sp.nLayers + 1;
        if (Tinit.size() == nNodes) {
            for (int i = 0; i < sp.nLayers; ++i) Te0[i] = Tinit[i];
            Te0[sp.nLayers] = Te0seed;
        }
        if (thetainit.size() == nNodes) {
            for (int i = 0; i <= sp.nLayers; ++i) wc0[i] = thetainit[i];
        }
        soilheatmod heat0 = toSoilheatmodCpp(soilc, Te0, wc0);
        soilwatermod water0 = toSoilwatermodCpp(sp, heat0, vp.root50, vp.root95);

        BigLeafOutput res = BigLeafCoreCpp(
            Rcpp::as<std::vector<double>>(year), Rcpp::as<std::vector<double>>(month),
            Rcpp::as<std::vector<double>>(day), Rcpp::as<std::vector<double>>(hour),
            Rcpp::as<std::vector<double>>(temp), Rcpp::as<std::vector<double>>(relhum),
            Rcpp::as<std::vector<double>>(pres), Rcpp::as<std::vector<double>>(swdown),
            Rcpp::as<std::vector<double>>(difrad), Rcpp::as<std::vector<double>>(lwdown),
            Rcpp::as<std::vector<double>>(windspeed), Rcpp::as<std::vector<double>>(precip),
            vp, sp, heat0, water0,
            gref, grefPAR, lat, lon,
            zref, reqhgt, Ca, zAbove, maxIter, tolerance, pooling, slope, aspect, svfa,
            Rcpp::as<std::vector<double>>(shelterc));

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("RabsGround") = res.RabsGround,
            Rcpp::Named("RabsCanopy") = res.RabsCanopy,
            Rcpp::Named("uf") = res.uf,
            Rcpp::Named("LL") = res.LL,
            Rcpp::Named("Hout") = res.Hout,
            Rcpp::Named("Tground") = res.Tground,
            Rcpp::Named("G") = res.G,
            Rcpp::Named("Tground0") = res.Tground0,
            Rcpp::Named("Tcanopy") = res.Tcanopy,
            Rcpp::Named("psi_r") = res.psi_r,
            Rcpp::Named("iters") = res.iters,
            Rcpp::Named("witers") = res.witers,
            Rcpp::Named("Tzbelow") = res.Tzbelow,
            Rcpp::Named("thetazbelow") = res.thetazbelow,
            Rcpp::Named("D_D") = res.D_D,
            Rcpp::Named("rHa") = res.rHa,
            Rcpp::Named("rBL") = res.rBL,
            Rcpp::Named("rSurf") = res.rSurf,
            Rcpp::Named("hSurf") = res.hSurf,
            Rcpp::Named("filmEvap") = res.filmEvap,
            Rcpp::Named("filmAvailable") = res.filmAvailable,
            Rcpp::Named("precipGround") = res.precipGround,
            Rcpp::Named("wetShare") = res.wetShare,
            Rcpp::Named("theta0") = res.theta0,
            Rcpp::Named("Tabove") = res.Tabove,
            Rcpp::Named("windAbove") = res.windAbove,
            Rcpp::Named("RHabove") = res.RHabove,
            Rcpp::Named("swaterdepth") = res.swaterdepth,
            Rcpp::Named("groundhr") = res.groundhr,
            Rcpp::Named("Tzabove") = res.Tzabove,
            Rcpp::Named("windzabove") = res.windzabove,
            Rcpp::Named("RHzabove") = res.RHzabove,
            Rcpp::Named("uh") = res.uh,
            Rcpp::Named("windScale") = res.windScale,
            Rcpp::Named("zenr") = res.zenr,
            Rcpp::Named("zend") = res.zend,
            Rcpp::Named("azid") = res.azid,
            Rcpp::Named("emGround") = res.emGround,
            Rcpp::Named("emCanopy") = res.emCanopy
        );
        df.attr("TeEnd") = res.TeEnd;
        df.attr("thetaEnd") = res.thetaEnd;
        return df;
    }

} // namespace pointmodel

// R-facing wrapper for the namespaced point-model implementation.
// [[Rcpp::export]]
Rcpp::DataFrame BigLeafCpp(
    Rcpp::NumericVector year, Rcpp::NumericVector month,
    Rcpp::NumericVector day, Rcpp::NumericVector hour,
    Rcpp::NumericVector temp, Rcpp::NumericVector relhum, Rcpp::NumericVector pres,
    Rcpp::NumericVector swdown, Rcpp::NumericVector difrad, Rcpp::NumericVector lwdown,
    Rcpp::NumericVector windspeed, Rcpp::NumericVector winddir, Rcpp::NumericVector precip,
    Rcpp::List vegp, Rcpp::List soilc, double Lfrac, double lat, double lon,
    double zref, double reqhgt = NA_REAL, double Ca = 400.0, double matemp = NA_REAL,
    double zAbove = NA_REAL, int maxIter = 100, double tolerance = 1e-2, bool pooling = false,
    double slope = 0.0, double aspect = 0.0, double svfa = 1.0,
    Rcpp::NumericVector shelterc = Rcpp::NumericVector::create(1.0),
    Rcpp::NumericVector Tinit = Rcpp::NumericVector::create(),
    Rcpp::NumericVector thetainit = Rcpp::NumericVector::create())
{
    return pointmodel::BigLeafCpp(year, month, day, hour, temp, relhum, pres,
        swdown, difrad, lwdown, windspeed, winddir, precip, vegp, soilc, Lfrac, lat, lon,
        zref, reqhgt, Ca, matemp, zAbove, maxIter, tolerance, pooling, slope, aspect, svfa, shelterc,
        Tinit, thetainit);
}
