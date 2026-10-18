// SPDX-License-Identifier: MIT
//! A GPS orbit model for predicting days ahead, integrated in GCRS.
//!
//! Forces: Earth's gravity field (gravity.rs), the Sun and the Moon (ERFA, via eo.rs), and solar
//! radiation pressure as five ECOM-like terms (constant along the Sun direction, the solar panel
//! axis and the third axis, plus a once-per-revolution term on the third), faded through Earth's
//! penumbra. Integrated with 8th-order Dormand-Prince at a fixed 100 s step: adaptive steps made
//! the result jitter with its inputs by centimetres to metres, and the fit could not converge.
use crate::dop853::{A, B, C, STAGES};
use crate::eo::EarthOrientation;
use crate::gravity::{gravity_tide, solid_tide, R_E};
use crate::linalg::*;

pub const GM_SUN: f64 = 1.32712440018e20;
pub const GM_MOON: f64 = 4.9028e12;
const R_SUN: f64 = 6.957e8;
pub const STEP: f64 = 100.0;

/// A priori radiation-pressure terms, nm/s²: D0 (along the Sun direction) is left free; the
/// others are usually a few nm/s² for GPS and are held near zero unless the data insists.
/// Terms, nm/s²: D0 Y0 B0 Bc Bs (ECOM1) and D2c D2s D4c D4s (the extra ECOM2 terms along the Sun
/// direction, twice and four times per revolution).
pub const NSRP: usize = 9;
pub const SRP_PRIOR_MEAN: [f64; NSRP] = [-100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
/// A priori sigmas: ECOM2 (all nine terms; the default). In two backtests it cut the day-4-5 orbit
/// error by 6-11% against ECOM1; looser or tighter priors did not help.
pub const SRP_SIGMA_ECOM2: [f64; NSRP] = [1e6, 2.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0];
const NP: usize = 6 + NSRP;

pub type State = [f64; 6];
pub type Srp = [f64; NSRP];

/// Fraction of sunlight reaching a satellite: 1 in sunlight, 0 in Earth's umbra, and a smooth
/// ramp across the penumbra (~190 km wide at GPS distance).
pub fn sunlit(r: V3, sun: V3) -> f64 {
    let sn = norm(sun);
    let ss = scale(sun, 1.0 / sn);
    let proj = dot(r, ss);
    if proj >= 0.0 { return 1.0; }
    let rho = norm(sub(r, scale(ss, proj)));
    let w = 2.0 * (-proj) * R_SUN / sn;
    let x = (rho - (R_E - w / 2.0)) / w;
    if x <= 0.0 { 0.0 } else if x >= 1.0 { 1.0 } else { x * x * (3.0 - 2.0 * x) }
}

/// Which optional forces to include (all on by default; switches exist to measure each).
#[derive(Clone, Copy)]
pub struct Forces { pub tides: bool, pub relativity: bool, pub srp_sigma: [f64; NSRP], pub box_wing: bool, pub antenna_offset: bool }

impl Default for Forces {
    fn default() -> Self { Forces { tides: true, relativity: true, srp_sigma: SRP_SIGMA_ECOM2, box_wing: false, antenna_offset: true } }
}

pub struct Orbit<'a> {
    pub eo: &'a EarthOrientation,
    pub t0: f64,
    pub forces: Forces,
    /// An a priori radiation-pressure model (Galileo); the ECOM terms then fit what it misses.
    pub box_wing: Option<crate::boxwing::BoxWing>,
    /// Fit a constant radial offset between the orbit and the track: on for tracks built from
    /// broadcast orbits, which describe the antenna (~1 m Earthward of the centre of mass). An
    /// orbit cannot follow a radially shifted point at the true angular rate; without this term
    /// the fit bought the radius with a wrong period, worth ~70 m along-track by day 5 (Galileo).
    pub fit_antenna_offset: bool,
}

impl<'a> Orbit<'a> {
    pub fn with(eo: &'a EarthOrientation, t0: f64, forces: Forces) -> Orbit<'a> { Orbit { eo, t0, forces, box_wing: None, fit_antenna_offset: false } }

    pub fn accel(&self, t: f64, r: V3, v: V3, p: &Srp) -> V3 {
        let m = self.eo.c2t(t);
        let (sun, moon) = self.eo.bodies(t);
        let (dc, ds) = if self.forces.tides {
            solid_tide(&[(mv(&m, sun), GM_SUN), (mv(&m, moon), GM_MOON)])
        } else { ([0.0; 3], [0.0; 3]) };
        let mut a = mtv(&m, gravity_tide(mv(&m, r), &dc, &ds));
        if self.forces.relativity {
            // Schwarzschild term, IERS Conventions 2010 eq. 10.12 (the Lense-Thirring and de Sitter
            // terms are far smaller for GPS)
            const C2: f64 = 299792458.0 * 299792458.0;
            let gm = crate::gravity::GM;
            let rn = norm(r);
            let k = gm / (C2 * rn * rn * rn);
            a = add(a, scale(add(scale(r, 4.0 * gm / rn - dot(v, v)), scale(v, 4.0 * dot(r, v))), k));
        }
        for (s, gm) in [(sun, GM_SUN), (moon, GM_MOON)] {
            let d = sub(s, r);
            let (dn, sn) = (norm(d), norm(s));
            a = add(a, scale(sub(scale(d, 1.0 / (dn * dn * dn)), scale(s, 1.0 / (sn * sn * sn))), gm));
        }
        let lit = sunlit(r, sun);
        if let Some(bw) = &self.box_wing { a = add(a, bw.accel(r, sun, lit)); }
        if lit > 0.0 {
            let ss = unit(sun);
            let e_d = unit(sub(sun, r));
            let e_r = unit(r);
            let e_y = unit(cross(e_d, e_r));
            let e_b = cross(e_d, e_y);
            let h = unit(cross(r, v));
            let sp = unit(sub(ss, scale(h, dot(ss, h))));
            let u = dot(cross(sp, e_r), h).atan2(dot(sp, e_r));
            let bt = p[2] + p[3] * u.cos() + p[4] * u.sin();
            let dt = p[0] + p[5] * (2.0 * u).cos() + p[6] * (2.0 * u).sin() + p[7] * (4.0 * u).cos() + p[8] * (4.0 * u).sin();
            let f = 1e-9 * lit;
            a = add(a, scale(add(add(scale(e_d, dt), scale(e_y, p[1])), scale(e_b, bt)), f));
        }
        a
    }

    fn deriv(&self, t: f64, y: &State, p: &Srp) -> State {
        let a = self.accel(t, [y[0], y[1], y[2]], [y[3], y[4], y[5]], p);
        [y[3], y[4], y[5], a[0], a[1], a[2]]
    }

    fn step(&self, t: f64, y: &State, h: f64, p: &Srp) -> State {
        let mut k = [[0.0; 6]; STAGES];
        for s in 0..STAGES {
            let mut ys = *y;
            for j in 0..s {
                if A[s][j] != 0.0 {
                    for i in 0..6 { ys[i] += h * A[s][j] * k[j][i]; }
                }
            }
            k[s] = self.deriv(t + C[s] * h, &ys, p);
        }
        let mut out = *y;
        for s in 0..STAGES {
            if B[s] != 0.0 {
                for i in 0..6 { out[i] += h * B[s] * k[s][i]; }
            }
        }
        out
    }

    /// States at the given increasing times (all >= t0), from `state` at t0. The main track
    /// advances by fixed steps; a time between steps is reached by one short step off the track,
    /// so every output is a smooth function of the inputs.
    pub fn propagate_states(&self, state: &State, p: &Srp, times: &[f64]) -> Vec<State> {
        let mut out = Vec::with_capacity(times.len());
        let (mut k, mut y) = (0u64, *state);
        for &target in times {
            while self.t0 + (k + 1) as f64 * STEP <= target + 1e-6 {
                y = self.step(self.t0 + k as f64 * STEP, &y, STEP, p);
                k += 1;
            }
            let t = self.t0 + k as f64 * STEP;
            let rest = target - t;
            out.push(if rest.abs() < 1e-6 { y } else { self.step(t, &y, rest, p) });
        }
        out
    }

    pub fn propagate(&self, state: &State, p: &Srp, times: &[f64]) -> Vec<V3> {
        self.propagate_states(state, p, times).iter().map(|s| [s[0], s[1], s[2]]).collect()
    }

    /// The initial state and radiation-pressure terms that best reproduce positions (GCRS, m) at
    /// times starting at t0: Levenberg-Marquardt on scaled unknowns with a finite-difference
    /// Jacobian, priors on the secondary radiation terms, run until a step moves the orbit by
    /// under 0.1 mm. `start` changes how quickly the fit gets there, not where.
    /// Returns (state, srp, rms m, iterations, converged).
    pub fn fit(&self, ts: &[f64], xs: &[V3], start: Option<(State, Srp)>) -> (State, Srp, f64, usize, bool) {
        const NQ: usize = NP + 1;                       // + the radial antenna offset, m
        let mut steps = [10.0; NQ];
        steps[NP] = 0.1;
        for j in 3..6 { steps[j] = 1e-2; }
        let mut p = [0.0; NQ];
        match start {
            Some((st, q)) => { p[..6].copy_from_slice(&st); p[6..NP].copy_from_slice(&q); }
            None => {
                let v0 = scale(sub(lagrange(ts, xs, ts[0] + 1.0, 10), lagrange(ts, xs, ts[0] - 1.0, 10)), 0.5);
                p[..3].copy_from_slice(&xs[0]);
                p[3..6].copy_from_slice(&v0);
                p[6..NP].copy_from_slice(&SRP_PRIOR_MEAN);
            }
        }
        let sig = self.forces.srp_sigma;
        let lam: Vec<f64> = (0..NQ).map(|j| if j < 6 { 0.0 } else if j == NP { if self.fit_antenna_offset { 0.0 } else { 1e12 } } else { 1.0 / sig[j - 6].powi(2) }).collect();
        let prior: Vec<f64> = (0..NQ).map(|j| if j < 6 || j == NP { 0.0 } else { SRP_PRIOR_MEAN[j - 6] }).collect();
        let split = |q: &[f64; NQ]| -> (State, Srp) {
            let mut st = [0.0; 6]; let mut r = [0.0; NSRP];
            st.copy_from_slice(&q[..6]); r.copy_from_slice(&q[6..NP]); (st, r)
        };
        let resid = |q: &[f64; NQ]| -> Vec<f64> {
            let (st, r) = split(q);
            self.propagate(&st, &r, ts).iter().zip(xs).flat_map(|(a, b)| sub(add(*a, scale(unit(*a), q[NP])), *b)).collect()
        };
        let cost = |q: &[f64; NQ], res: &[f64]| -> f64 {
            res.iter().map(|x| x * x).sum::<f64>() + (0..NQ).map(|j| lam[j] * (q[j] - prior[j]).powi(2)).sum::<f64>()
        };
        let mut res = resid(&p);
        let mut c = cost(&p, &res);
        let (mut mu, mut converged, mut it) = (1e-3, false, 0);
        let n = res.len();
        while it < 25 {
            it += 1;
            // Jacobian, column by column
            let mut jac = vec![0.0; n * NQ];
            for j in 0..NQ {
                let mut q = p;
                q[j] += steps[j];
                let r2 = resid(&q);
                for i in 0..n { jac[i * NQ + j] = (r2[i] - res[i]) / steps[j]; }
            }
            // Scale every unknown to unit effect (velocity partials are ~1e5 times position ones).
            let d: Vec<f64> = (0..NQ).map(|j| (0..n).map(|i| jac[i * NQ + j].powi(2)).sum::<f64>().sqrt()).collect();
            let mut nm = vec![0.0; NQ * NQ];
            let mut g = vec![0.0; NQ];
            for a in 0..NQ {
                for b in 0..NQ {
                    nm[a * NQ + b] = (0..n).map(|i| jac[i * NQ + a] * jac[i * NQ + b]).sum::<f64>() / (d[a] * d[b]);
                }
                nm[a * NQ + a] += lam[a] / (d[a] * d[a]);
                g[a] = (0..n).map(|i| jac[i * NQ + a] * res[i]).sum::<f64>() / d[a] + lam[a] * (p[a] - prior[a]) / d[a];
            }
            let mut moved;
            loop {
                let mut damped = nm.clone();
                for a in 0..NQ { damped[a * NQ + a] += mu * nm[a * NQ + a]; }
                let step = solve(NQ, &damped, &g.iter().map(|x| -x).collect::<Vec<_>>());
                let dp: Vec<f64> = (0..NQ).map(|j| step[j] / d[j]).collect();
                moved = (0..n).map(|i| (0..NQ).map(|j| jac[i * NQ + j] * dp[j]).sum::<f64>().abs()).fold(0.0, f64::max);
                let mut q = p;
                for j in 0..NQ { q[j] += dp[j]; }
                let rq = resid(&q);
                let cq = cost(&q, &rq);
                if cq <= c {
                    p = q; res = rq; c = cq; mu = (mu / 10.0).max(1e-6);
                    break;
                }
                if moved < 1e-5 || mu > 1e6 { break; }
                mu *= 10.0;
            }
            if moved < 1e-4 { converged = true; break; }
        }
        let (s, r) = split(&p);
        let rms = (res.iter().map(|x| x * x).sum::<f64>() / n as f64).sqrt();
        (s, r, rms, it, converged)
    }
}
