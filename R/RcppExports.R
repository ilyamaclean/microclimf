# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

solpositionCpp <- function(lat, lon, year, month, day, lt) {
    .Call(`_microclimf_solpositionCpp`, lat, lon, year, month, day, lt)
}

solarindexCpp <- function(slope, aspect, zen, azi, shadowmask = FALSE) {
    .Call(`_microclimf_solarindexCpp`, slope, aspect, zen, azi, shadowmask)
}

clearskyradCpp <- function(year, month, day, lt, lat, lon, tc, rh, pk) {
    .Call(`_microclimf_clearskyradCpp`, year, month, day, lt, lat, lon, tc, rh, pk)
}

zeroplanedisCpp <- function(h, pai) {
    .Call(`_microclimf_zeroplanedisCpp`, h, pai)
}

roughlengthCpp <- function(h, pai, d, psi_h) {
    .Call(`_microclimf_roughlengthCpp`, h, pai, d, psi_h)
}

dpsimCpp <- function(ze) {
    .Call(`_microclimf_dpsimCpp`, ze)
}

dpsihCpp <- function(ze) {
    .Call(`_microclimf_dpsihCpp`, ze)
}

dphihCpp <- function(ze) {
    .Call(`_microclimf_dphihCpp`, ze)
}

satvapCpp <- function(tc) {
    .Call(`_microclimf_satvapCpp`, tc)
}

dewpointCpp <- function(tc, ea) {
    .Call(`_microclimf_dewpointCpp`, tc, ea)
}

hourtodayCpp <- function(hourly, stat) {
    .Call(`_microclimf_hourtodayCpp`, hourly, stat)
}

maCpp <- function(x, n) {
    .Call(`_microclimf_maCpp`, x, n)
}

mayCpp <- function(x) {
    .Call(`_microclimf_mayCpp`, x)
}

BigLeafCpp <- function(obstime, climdata, vegp, groundp, soilm, lat, lon, dTmx = 25.0, zref = 2.0, maxiter = 100L, bwgt = 0.5, tol = 0.5, gmn = 0.1, yearG = TRUE) {
    .Call(`_microclimf_BigLeafCpp`, obstime, climdata, vegp, groundp, soilm, lat, lon, dTmx, zref, maxiter, bwgt, tol, gmn, yearG)
}

weatherhgtCpp <- function(obstime, climdata, zin, uzin, zout, lat, lon) {
    .Call(`_microclimf_weatherhgtCpp`, obstime, climdata, zin, uzin, zout, lat, lon)
}

soilmCpp <- function(climdata, rmu, mult, pwr, Smax, Smin, Ksat, a) {
    .Call(`_microclimf_soilmCpp`, climdata, rmu, mult, pwr, Smax, Smin, Ksat, a)
}

pointmprocess <- function(pointvars, zref, h, pai, rho, Vm, Vq, Mc) {
    .Call(`_microclimf_pointmprocess`, pointvars, zref, h, pai, rho, Vm, Vq, Mc)
}

aperm3D2 <- function(Tz, rows, cols, tsteps) {
    .Call(`_microclimf_aperm3D2`, Tz, rows, cols, tsteps)
}

slice_2d <- function(a, k, dim) {
    .Call(`_microclimf_slice_2d`, a, k, dim)
}

applycpp3 <- function(a, fun_name) {
    .Call(`_microclimf_applycpp3`, a, fun_name)
}

manCpp <- function(x, n) {
    .Call(`_microclimf_manCpp`, x, n)
}

flowdirCpp <- function(md) {
    .Call(`_microclimf_flowdirCpp`, md)
}

flowaccCpp <- function(dm) {
    .Call(`_microclimf_flowaccCpp`, dm)
}

soildCppv <- function(soilm, Smin, Smax, tadd) {
    .Call(`_microclimf_soildCppv`, soilm, Smin, Smax, tadd)
}

soildCppm <- function(twi, Smin, Smax, tfact) {
    .Call(`_microclimf_soildCppm`, twi, Smin, Smax, tfact)
}

solargrid <- function(slope, aspect, obstime, micro) {
    .Call(`_microclimf_solargrid`, slope, aspect, obstime, micro)
}

twostreamgrid <- function(reqhgt, micro) {
    .Call(`_microclimf_twostreamgrid`, reqhgt, micro)
}

windgrid <- function(reqhgt, micro) {
    .Call(`_microclimf_windgrid`, reqhgt, micro)
}

PenmanMonteith2Cpp <- function(Rabs, gHa, gV, tc, mxtc, pk, ea, es, G, surfwet, tdew) {
    .Call(`_microclimf_PenmanMonteith2Cpp`, Rabs, gHa, gV, tc, mxtc, pk, ea, es, G, surfwet, tdew)
}

soiltempgrid <- function(micro) {
    .Call(`_microclimf_soiltempgrid`, micro)
}

rhcanopy <- function(uf, h, d, z) {
    .Call(`_microclimf_rhcanopy`, uf, h, d, z)
}

abovegrid <- function(reqhgt, micro) {
    .Call(`_microclimf_abovegrid`, reqhgt, micro)
}

belowgrid <- function(reqhgt, micro, hiy, complete) {
    .Call(`_microclimf_belowgrid`, reqhgt, micro, hiy, complete)
}

runmicro1Cpp <- function(obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lat, lon, Sminp, Smaxp, tfact, complete, mat, out) {
    .Call(`_microclimf_runmicro1Cpp`, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lat, lon, Sminp, Smaxp, tfact, complete, mat, out)
}

runmicro2Cpp <- function(obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lats, lons, Sminp, Smaxp, tfact, complete, mat, out) {
    .Call(`_microclimf_runmicro2Cpp`, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lats, lons, Sminp, Smaxp, tfact, complete, mat, out)
}

runmicro3Cpp <- function(dfsel, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lat, lon, Sminp, Smaxp, tfact, complete, mat, out) {
    .Call(`_microclimf_runmicro3Cpp`, dfsel, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lat, lon, Sminp, Smaxp, tfact, complete, mat, out)
}

runmicro4Cpp <- function(dfsel, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lats, lons, Sminp, Smaxp, tfact, complete, mat, out) {
    .Call(`_microclimf_runmicro4Cpp`, dfsel, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lats, lons, Sminp, Smaxp, tfact, complete, mat, out)
}

runbioclim1Cpp <- function(obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lat, lon, Sminp, Smaxp, tfact, mat, out, wetq, dryq, hotq, colq, air) {
    .Call(`_microclimf_runbioclim1Cpp`, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lat, lon, Sminp, Smaxp, tfact, mat, out, wetq, dryq, hotq, colq, air)
}

runbioclim2Cpp <- function(obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lats, lons, Sminp, Smaxp, tfact, mat, out, wetq, dryq, hotq, colq, air) {
    .Call(`_microclimf_runbioclim2Cpp`, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lats, lons, Sminp, Smaxp, tfact, mat, out, wetq, dryq, hotq, colq, air)
}

runbioclim3Cpp <- function(obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lat, lon, Sminp, Smaxp, tfact, mat, out, wetq, dryq, hotq, colq, air) {
    .Call(`_microclimf_runbioclim3Cpp`, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lat, lon, Sminp, Smaxp, tfact, mat, out, wetq, dryq, hotq, colq, air)
}

runbioclim4Cpp <- function(obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lats, lons, Sminp, Smaxp, tfact, mat, out, wetq, dryq, hotq, colq, air) {
    .Call(`_microclimf_runbioclim4Cpp`, obstime, climdata, pointm, vegp, soilc, reqhgt, zref, lats, lons, Sminp, Smaxp, tfact, mat, out, wetq, dryq, hotq, colq, air)
}

pointmodelsnow <- function(obstime, climdata, vegp, other, snowenv, tol = 0.5, maxiter = 100) {
    .Call(`_microclimf_pointmodelsnow`, obstime, climdata, vegp, other, snowenv, tol, maxiter)
}

gridmodelsnow1 <- function(obstime, climdata, pointm, vegp, other, snowenv) {
    .Call(`_microclimf_gridmodelsnow1`, obstime, climdata, pointm, vegp, other, snowenv)
}

snowdayan <- function(stempg) {
    .Call(`_microclimf_snowdayan`, stempg)
}

canintfrac <- function(hgt, pai, uf, prec, tc, Li) {
    .Call(`_microclimf_canintfrac`, hgt, pai, uf, prec, tc, Li)
}

meltmu <- function(mu, stemp, tc) {
    .Call(`_microclimf_meltmu`, mu, stemp, tc)
}

meltmu2 <- function(mu, stemp, tc) {
    .Call(`_microclimf_meltmu2`, mu, stemp, tc)
}

snowdaysfun <- function(maxsnowdepth, minsnowdepth) {
    .Call(`_microclimf_snowdaysfun`, maxsnowdepth, minsnowdepth)
}

gridmicrosnow1 <- function(reqhgt, obstime, climdata, snowm, micro, vegp, other, out) {
    .Call(`_microclimf_gridmicrosnow1`, reqhgt, obstime, climdata, snowm, micro, vegp, other, out)
}

gridmodelsnow2 <- function(obstime, climdata, pointm, vegp, other, snowenv) {
    .Call(`_microclimf_gridmodelsnow2`, obstime, climdata, pointm, vegp, other, snowenv)
}

gridmicrosnow2 <- function(reqhgt, obstime, climdata, snowm, micro, vegp, other, out) {
    .Call(`_microclimf_gridmicrosnow2`, reqhgt, obstime, climdata, snowm, micro, vegp, other, out)
}

leafrcpp <- function(om, pai, gref, albin, x, ltrr) {
    .Call(`_microclimf_leafrcpp`, om, pai, gref, albin, x, ltrr)
}

