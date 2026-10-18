// SPDX-License-Identifier: MIT
//! GPS LNAV broadcast ephemeris (IS-GPS-200): evaluate one, and fit one to a trajectory.
//! Positions are ECEF metres; times are continuous GPS seconds.
use crate::linalg::*;
use crate::time::WEEK;

pub const GM: f64 = 3.986005e14;          // IS-GPS-200's values, which receivers use
pub const OMEGA_E: f64 = 7.2921151467e-5;

/// Orbit parameters in this order: sqrtA, e, M0, omega, i0, Omega0, OmegaDot, deltaN, IDOT,
/// Cuc, Cus, Crc, Crs, Cic, Cis (angles in radians).
pub type Params = [f64; 15];

pub fn sat_pos(p: &Params, toe: f64, t: f64) -> V3 {
    let [sqrt_a, e, m0, w, i0, w0, wd, dn, idot, cuc, cus, crc, crs, cic, cis] = *p;
    let a = sqrt_a * sqrt_a;
    let tk = t - toe;
    let m = m0 + ((GM / (a * a * a)).sqrt() + dn) * tk;
    let mut ea = m;
    for _ in 0..12 { ea -= (ea - e * ea.sin() - m) / (1.0 - e * ea.cos()); }
    let phi = ((1.0 - e * e).sqrt() * ea.sin()).atan2(ea.cos() - e) + w;
    let (s2, c2) = (2.0 * phi).sin_cos();
    let u = phi + cus * s2 + cuc * c2;
    let r = a * (1.0 - e * ea.cos()) + crs * s2 + crc * c2;
    let i = i0 + idot * tk + cis * s2 + cic * c2;
    let (x, y) = (r * u.cos(), r * u.sin());
    let om = w0 + (wd - OMEGA_E) * tk - OMEGA_E * toe.rem_euclid(WEEK);
    let (so, co) = om.sin_cos();
    let ci = i.cos();
    [x * co - y * ci * so, x * so + y * ci * co, y * i.sin()]
}

/// Osculating elements, in broadcast form, from an ECEF state at toe: the fit's starting point.
fn elements(r: V3, v_ecef: V3, toe: f64) -> Params {
    let v = add(v_ecef, cross([0.0, 0.0, OMEGA_E], r));
    let h = cross(r, v);
    let hn = norm(h);
    let i = (h[2] / hn).acos();
    let node = cross([0.0, 0.0, 1.0], h);
    let nhat = unit(node);
    let ev = sub(scale(cross(v, h), 1.0 / GM), unit(r));
    let e = norm(ev);
    let a = 1.0 / (2.0 / norm(r) - dot(v, v) / GM);
    let w = (dot(cross(nhat, ev), h) / hn).atan2(dot(nhat, ev));
    let u = (dot(cross(nhat, r), h) / hn).atan2(dot(nhat, r));
    let nu = u - w;
    let ea = ((1.0 - e * e).sqrt() * nu.sin()).atan2(e + nu.cos());
    let pi = std::f64::consts::PI;
    let w0 = (node[1].atan2(node[0]) + OMEGA_E * toe.rem_euclid(WEEK) + pi).rem_euclid(2.0 * pi) - pi;
    [a.sqrt(), e, ea - e * ea.sin(), w, i, w0, -8.0e-9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
}

// The fit runs in non-singular variables (e cos w, e sin w, M0 + w), which stay well defined for
// near-circular orbits, scaled so that a unit step moves the satellite by roughly a metre.
const SCALE: Params = [1e-4, 1e-7, 1e-8, 1e-7, 1e-8, 1e-8, 1e-11, 1e-11, 1e-11, 1e-7, 1e-7, 1.0, 1.0, 1e-7, 1e-7];

fn to_ns(p: &Params) -> Params {
    let mut q = *p;
    q[1] = p[1] * p[3].cos();
    q[3] = p[1] * p[3].sin();
    q[2] = p[2] + p[3];
    q
}

fn from_ns(q: &Params) -> Params {
    let pi = std::f64::consts::PI;
    let mut p = *q;
    p[1] = q[1].hypot(q[3]);
    p[3] = q[3].atan2(q[1]);
    p[2] = (q[2] - p[3] + pi).rem_euclid(2.0 * pi) - pi;
    p[5] = (p[5] + pi).rem_euclid(2.0 * pi) - pi;
    p
}

/// Fit ephemeris parameters to positions (ECEF, m) at times spanning toe ± 2 h, by Levenberg-
/// Marquardt with a central-difference Jacobian (steps of ~1 cm: much smaller steps drown in
/// rounding at GPS distances). Returns (parameters, rms residual in metres).
pub fn fit(ts: &[f64], xs: &[V3], toe: f64) -> (Params, f64) {
    let n = ts.len();
    let k = (0..n).min_by(|&a, &b| (ts[a] - toe).abs().partial_cmp(&(ts[b] - toe).abs()).unwrap()).unwrap().clamp(1, n - 2);
    let vel = scale(sub(xs[k + 1], xs[k - 1]), 1.0 / (ts[k + 1] - ts[k - 1]));
    let mut p0 = elements(xs[k], vel, ts[k]);
    p0[2] += (GM / p0[0].powi(6)).sqrt() * (toe - ts[k]);
    let q0 = to_ns(&p0);
    let resid = |d: &Params| -> Vec<f64> {
        let mut q = q0;
        for j in 0..15 { q[j] += d[j] * SCALE[j]; }
        let p = from_ns(&q);
        ts.iter().zip(xs).flat_map(|(&t, &x)| sub(sat_pos(&p, toe, t), x)).collect()
    };
    let m = 3 * n;
    let mut d = [0.0; 15];
    let mut res = resid(&d);
    let mut c: f64 = res.iter().map(|x| x * x).sum();
    let mut mu = 1e-3;
    for _ in 0..60 {
        let mut jac = vec![0.0; m * 15];
        for j in 0..15 {
            let (mut a, mut b) = (d, d);
            a[j] += 1e-2;
            b[j] -= 1e-2;
            let (ra, rb) = (resid(&a), resid(&b));
            for i in 0..m { jac[i * 15 + j] = (ra[i] - rb[i]) / 2e-2; }
        }
        let mut nm = vec![0.0; 225];
        let mut g = vec![0.0; 15];
        for a in 0..15 {
            for b in a..15 {
                let v: f64 = (0..m).map(|i| jac[i * 15 + a] * jac[i * 15 + b]).sum();
                nm[a * 15 + b] = v;
                nm[b * 15 + a] = v;
            }
            g[a] = -(0..m).map(|i| jac[i * 15 + a] * res[i]).sum::<f64>();
        }
        let mut improved = false;
        let mut small = false;
        for _ in 0..12 {
            let mut damped = nm.clone();
            for a in 0..15 { damped[a * 15 + a] *= 1.0 + mu; }
            let step = solve(15, &damped, &g);
            let mut dn = d;
            for j in 0..15 { dn[j] += step[j]; }
            let rn = resid(&dn);
            let cn: f64 = rn.iter().map(|x| x * x).sum();
            small = step.iter().map(|x| x.abs()).fold(0.0, f64::max) < 1e-6;   // ~1 µm
            if cn <= c {
                improved = c - cn > 1e-14 * c.max(1e-30);
                d = dn; res = rn; c = cn; mu = (mu / 10.0).max(1e-9);
                break;
            }
            mu *= 10.0;
        }
        if !improved || small { break; }
    }
    let mut q = q0;
    for j in 0..15 { q[j] += d[j] * SCALE[j]; }
    (from_ns(&q), (c / m as f64).sqrt())
}
