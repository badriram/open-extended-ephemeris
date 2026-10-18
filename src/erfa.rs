// SPDX-License-Identifier: MIT
//! The few ERFA routines used, from the vendored C (vendor/erfa, ERFA 2.0.1, BSD-3-Clause):
//! IAU 2006/2000A CIO-based Earth orientation, the Earth rotation angle, polar motion, and the
//! Sun and Moon. Dates are two-part Julian dates as in ERFA.

pub type M3 = [[f64; 3]; 3];

extern "C" {
    fn eraC2i06a(date1: f64, date2: f64, rc2i: *mut [f64; 3]);
    fn eraEra00(dj1: f64, dj2: f64) -> f64;
    fn eraSp00(date1: f64, date2: f64) -> f64;
    fn eraPom00(xp: f64, yp: f64, sp: f64, rpom: *mut [f64; 3]);
    fn eraEpv00(date1: f64, date2: f64, pvh: *mut [f64; 3], pvb: *mut [f64; 3]) -> i32;
    fn eraMoon98(date1: f64, date2: f64, pv: *mut [f64; 3]);
}

/// Celestial-to-intermediate matrix (IAU 2006/2000A), TT.
pub fn c2i06a(d1: f64, d2: f64) -> M3 {
    let mut m = [[0.0; 3]; 3];
    unsafe { eraC2i06a(d1, d2, m.as_mut_ptr()) };
    m
}

/// Earth rotation angle, UT1.
pub fn era00(d1: f64, d2: f64) -> f64 {
    unsafe { eraEra00(d1, d2) }
}

/// The TIO locator s', TT.
pub fn sp00(d1: f64, d2: f64) -> f64 {
    unsafe { eraSp00(d1, d2) }
}

/// Polar motion matrix.
pub fn pom00(xp: f64, yp: f64, sp: f64) -> M3 {
    let mut m = [[0.0; 3]; 3];
    unsafe { eraPom00(xp, yp, sp, m.as_mut_ptr()) };
    m
}

/// Heliocentric position of the Earth, AU, BCRS axes (TDB, for which TT serves here).
pub fn earth_heliocentric(d1: f64, d2: f64) -> [f64; 3] {
    let mut pvh = [[0.0; 3]; 2];
    let mut pvb = [[0.0; 3]; 2];
    unsafe { eraEpv00(d1, d2, pvh.as_mut_ptr(), pvb.as_mut_ptr()) };
    pvh[0]
}

/// Geocentric position of the Moon, AU, GCRS axes (TT).
pub fn moon_geocentric(d1: f64, d2: f64) -> [f64; 3] {
    let mut pv = [[0.0; 3]; 2];
    unsafe { eraMoon98(d1, d2, pv.as_mut_ptr()) };
    pv[0]
}
