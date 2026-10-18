// SPDX-License-Identifier: MIT
//! Earth orientation (GCRS -> ITRS) and the Sun and Moon, from ERFA and IGS ERP rows. What
//! changes slowly is tabulated every 30 minutes and interpolated: calling ERFA at every
//! integration step would dominate the run time.
//!
//! The CIO-based celestial-to-intermediate matrix (c2i06a) pairs with the Earth rotation angle.
//! (The equinox-based matrix pairs with sidereal time instead; mixing it with the ERA turned the
//! Earth by ~0.5 degrees in an earlier version.)
use crate::erfa;
use crate::linalg::{mm, M3, V3};
use crate::time::{JD_GPS0, TT_GPS, UTC_GPS};

pub const AU: f64 = 1.495978707e11;
const STEP: f64 = 1800.0;

/// One ERP row: GPS seconds, x and y pole (rad), UT1 - UTC (s).
pub type Erp = (f64, f64, f64, f64);

/// IGS ERP text -> its rows, as given.
pub fn read_erp(text: &str) -> Result<Vec<Erp>, String> {
    let arc = std::f64::consts::PI / 180.0 / 3600.0 * 1e-6;
    let mut rows = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 5 || !f[0].starts_with(|c: char| c.is_ascii_digit()) { continue; }
        let num = |s: &str| s.parse::<f64>().map_err(|e| format!("ERP: {e}"));
        let t = (num(f[0])? - 44244.0) * 86400.0 - UTC_GPS;
        rows.push((t, num(f[1])? * arc, num(f[2])? * arc, num(f[3])? * 1e-7, num(f[4])? * 1e-7));
    }
    if rows.is_empty() { return Err("no ERP rows".into()); }
    Ok(rows.into_iter().map(|r| (r.0, r.1, r.2, r.3)).collect())
}

/// The last length of day (s) in an IGS ERP file.
pub fn last_lod(text: &str) -> f64 {
    text.lines().filter_map(|l| {
        let f: Vec<&str> = l.split_whitespace().collect();
        if f.len() >= 5 && f[0].starts_with(|c: char| c.is_ascii_digit()) { f[4].parse::<f64>().ok() } else { None }
    }).last().unwrap_or(0.0) * 1e-7
}

/// Extend ERP rows ten days ahead, UT1 drifting by the last length of day and the pole held.
pub fn extrapolate(rows: &[Erp], lod: f64) -> Vec<Erp> {
    let mut out = rows.to_vec();
    let (t, xp, yp, ut) = *rows.last().unwrap();
    for d in 1..=10 {
        out.push((t + d as f64 * 86400.0, xp, yp, ut - lod * d as f64));
    }
    out
}

/// IERS Bulletin A's prediction table (lines "year month day MJD x" y" UT1-UTC s") -> rows.
pub fn read_bulletin_a(text: &str) -> Vec<Erp> {
    let arc = std::f64::consts::PI / 180.0 / 3600.0;
    let mut rows = Vec::new();
    let mut in_pred = false;
    for line in text.lines() {
        if line.contains("PREDICTIONS:") { in_pred = true; continue; }
        if !in_pred { continue; }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() == 7 && f[0].len() == 4 && f[0].starts_with("20") {
            let v: Vec<f64> = f.iter().filter_map(|x| x.parse().ok()).collect();
            if v.len() == 7 {
                rows.push(((v[3] - 44244.0) * 86400.0 - UTC_GPS, v[4] * arc, v[5] * arc, v[6]));
            }
        }
    }
    rows
}

/// IGS ERP rows where they exist (they are measured, and more accurate), then Bulletin A's
/// predictions, shifted to join the IGS values at the last IGS row. Bulletin A's predictions
/// include tidal and atmospheric terms that extrapolating the length of day misses.
pub fn join_bulletin_a(igs: &[Erp], bul: &[Erp]) -> Vec<Erp> {
    let measured = igs.to_vec();
    let last = *measured.last().unwrap();
    if bul.len() < 2 || bul[0].0 > last.0 || bul.last().unwrap().0 < last.0 + 7.0 * 86400.0 {
        return vec![];
    }
    let off = (last.1 - interp(bul, last.0, 1), last.2 - interp(bul, last.0, 2), last.3 - interp(bul, last.0, 3));
    let mut rows = measured.clone();
    rows.extend(bul.iter().filter(|r| r.0 > last.0).map(|r| (r.0, r.1 + off.0, r.2 + off.1, r.3 + off.2)));
    rows
}

fn interp(rows: &[Erp], t: f64, k: usize) -> f64 {
    let get = |r: &Erp| match k { 1 => r.1, 2 => r.2, _ => r.3 };
    if t <= rows[0].0 { return get(&rows[0]); }
    for w in rows.windows(2) {
        if t <= w[1].0 {
            let f = (t - w[0].0) / (w[1].0 - w[0].0);
            return get(&w[0]) * (1.0 - f) + get(&w[1]) * f;
        }
    }
    get(rows.last().unwrap())
}

/// Natural cubic spline on a uniform grid, one per coordinate.
struct Spline { t0: f64, h: f64, y: Vec<V3>, m: Vec<V3> }

impl Spline {
    fn new(t0: f64, h: f64, y: Vec<V3>) -> Spline {
        let n = y.len();
        let mut m = vec![[0.0; 3]; n];
        for c in 0..3 {
            // tridiagonal system for second derivatives, natural ends
            let mut cp = vec![0.0; n];
            let mut dp = vec![0.0; n];
            for i in 1..n - 1 {
                let rhs = 6.0 * (y[i + 1][c] - 2.0 * y[i][c] + y[i - 1][c]) / (h * h);
                let den = 4.0 - cp[i - 1];
                cp[i] = 1.0 / den;
                dp[i] = (rhs - dp[i - 1]) / den;
            }
            for i in (1..n - 1).rev() {
                m[i][c] = dp[i] - cp[i] * m[i + 1][c];
            }
        }
        Spline { t0, h, y, m }
    }

    fn at(&self, t: f64) -> V3 {
        let k = (((t - self.t0) / self.h).floor() as isize).clamp(0, self.y.len() as isize - 2) as usize;
        let a = (self.t0 + (k + 1) as f64 * self.h - t) / self.h;
        let b = 1.0 - a;
        let h2 = self.h * self.h / 6.0;
        let mut out = [0.0; 3];
        for c in 0..3 {
            out[c] = a * self.y[k][c] + b * self.y[k + 1][c]
                + ((a * a * a - a) * self.m[k][c] + (b * b * b - b) * self.m[k + 1][c]) * h2;
        }
        out
    }
}

pub struct EarthOrientation {
    t0: f64,
    dut1: Vec<f64>,
    pom: Vec<M3>,
    c2i: Vec<M3>,
    sun: Spline,
    moon: Spline,
}

impl EarthOrientation {
    pub fn new(erp: &[Erp], t_start: f64, t_end: f64) -> EarthOrientation {
        let t0 = t_start - 3.0 * 3600.0;
        let n = ((t_end + 86400.0 - t0) / STEP).ceil() as usize + 1;
        let (mut dut1, mut pom, mut c2i, mut sun, mut moon) = (vec![], vec![], vec![], vec![], vec![]);
        for i in 0..n {
            let t = t0 + i as f64 * STEP;
            let d2 = (t + TT_GPS) / 86400.0;
            dut1.push(interp(erp, t, 3));
            pom.push(erfa::pom00(interp(erp, t, 1), interp(erp, t, 2), erfa::sp00(JD_GPS0, d2)));
            c2i.push(erfa::c2i06a(JD_GPS0, d2));
            let e = erfa::earth_heliocentric(JD_GPS0, d2);
            sun.push([-e[0] * AU, -e[1] * AU, -e[2] * AU]);
            let m = erfa::moon_geocentric(JD_GPS0, d2);
            moon.push([m[0] * AU, m[1] * AU, m[2] * AU]);
        }
        EarthOrientation { t0, dut1, pom, c2i, sun: Spline::new(t0, STEP, sun), moon: Spline::new(t0, STEP, moon) }
    }

    fn lin<T: Copy, F: Fn(T, T, f64) -> T>(&self, table: &[T], t: f64, mix: F) -> T {
        let k = (((t - self.t0) / STEP).floor() as isize).clamp(0, table.len() as isize - 2) as usize;
        let f = (t - (self.t0 + k as f64 * STEP)) / STEP;
        mix(table[k], table[k + 1], f)
    }

    /// GCRS -> ITRS rotation matrix at GPS time t.
    pub fn c2t(&self, t: f64) -> M3 {
        let mixm = |a: M3, b: M3, f: f64| {
            let mut c = [[0.0; 3]; 3];
            for i in 0..3 { for j in 0..3 { c[i][j] = a[i][j] * (1.0 - f) + b[i][j] * f; } }
            c
        };
        let dut1 = self.lin(&self.dut1, t, |a, b, f| a * (1.0 - f) + b * f);
        let era = erfa::era00(JD_GPS0, (t + UTC_GPS + dut1) / 86400.0);
        let (s, c) = era.sin_cos();
        let rz = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]];
        mm(&mm(&self.lin(&self.pom, t, mixm), &rz), &self.lin(&self.c2i, t, mixm))
    }

    /// Geocentric Sun and Moon (m, GCRS) at GPS time t.
    pub fn bodies(&self, t: f64) -> (V3, V3) {
        (self.sun.at(t), self.moon.at(t))
    }
}
