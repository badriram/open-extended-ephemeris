// SPDX-License-Identifier: MIT
//! An a priori box-wing model of solar radiation pressure for Galileo satellites, from the
//! European GNSS Service Centre's satellite metadata (gsc-europa.eu, "Galileo Satellite Metadata",
//! section 6: geometry and optical properties; section 3: attitude; section 4: mass).
//!
//! Each surface of area A, outward normal n, lit at cos θ = n·e (e: unit vector to the Sun) feels
//!   F = -P A cos θ [ (α + δ) e + (2/3 δ + 2 ρ cos θ) n ]  (Montenbruck & Gill, eq. 3.75 with diffusion)
//! plus, for the multi-layer insulation, instantaneous re-emission of what it absorbs, -P A cos θ (2/3 α) n.
//! P is the solar flux over c at the satellite's distance from the Sun.
//!
//! Attitude, nominal yaw steering: +Z toward the Earth, Y perpendicular to the Sun, +X (the
//! radiators) toward deep space; the solar wings rotate about Y to face the Sun.
use crate::linalg::*;

const P0: f64 = 1367.0 / 299792458.0;          // N/m² at 1 AU
const AU: f64 = 1.495978707e11;

/// One surface: body-frame normal, area m², absorption α, specular ρ, diffuse δ, re-emits (MLI).
struct Face { n: V3, area: f64, a: f64, r: f64, d: f64, mli: bool }

pub struct BoxWing { mass: f64, body: Vec<Face>, wing: Vec<(f64, f64, f64, f64)> }

const MLI_IOV: (f64, f64, f64) = (0.94, 0.00, 0.06);   // carbon-filled Kapton
const RAD_IOV: (f64, f64, f64) = (0.10, 0.72, 0.18);   // optical surface radiator, beginning of life
const GE_KAPTON: (f64, f64, f64) = (0.57, 0.22, 0.21); // germanium-coated black Kapton
const MLI_FOC: (f64, f64, f64) = (0.93, 0.00, 0.07);   // material A
const RAD_FOC: (f64, f64, f64) = (0.08, 0.73, 0.19);   // material C

fn face(n: V3, area: f64, m: (f64, f64, f64), mli: bool) -> Face { Face { n, area, a: m.0, r: m.1, d: m.2, mli } }

impl BoxWing {
    /// Galileo In-Orbit Validation satellites (GSAT0101-0103: E11, E12, E19).
    pub fn galileo_iov(mass: f64) -> BoxWing {
        let (x, y, z) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        let neg = |v: V3| scale(v, -1.0);
        BoxWing { mass, body: vec![
            face(neg(x), 1.32, MLI_IOV, true),
            face(x, 0.54, MLI_IOV, true), face(x, 0.78, RAD_IOV, false),
            face(y, 1.00, MLI_IOV, true), face(y, 2.00, RAD_IOV, false),
            face(neg(y), 1.03, MLI_IOV, true), face(neg(y), 1.97, RAD_IOV, false),
            face(z, 1.72, MLI_IOV, true), face(z, 1.28, GE_KAPTON, false),
            face(neg(z), 3.00, MLI_IOV, true),
        ], wing: vec![(2.0 * 3.88, 0.92, 0.08, 0.00), (2.0 * 1.53, 0.90, 0.10, 0.00)] }
    }

    /// Galileo Full Operational Capability satellites (GSAT02xx).
    pub fn galileo_foc(mass: f64) -> BoxWing {
        let (x, y, z) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        let neg = |v: V3| scale(v, -1.0);
        BoxWing { mass, body: vec![
            face(x, 0.440, MLI_FOC, true), face(x, 0.880, RAD_FOC, false),
            face(neg(x), 1.320, MLI_FOC, true),
            face(y, 1.129, MLI_FOC, true), face(y, 1.654, RAD_FOC, false),
            face(neg(y), 1.244, MLI_FOC, true), face(neg(y), 1.539, RAD_FOC, false),
            face(z, 1.053, MLI_FOC, true), face(z, 1.969, GE_KAPTON, false),
            face(neg(z), 2.077, MLI_FOC, true), face(neg(z), 0.959, RAD_FOC, false),
        ], wing: vec![(2.0 * 3.880, 0.92, 0.08, 0.00), (2.0 * 1.530, 0.90, 0.10, 0.00)] }
    }

    /// The model for a Galileo PRN: IOV geometry for E11, E12, E19; the two eccentric-orbit FOC
    /// satellites (E14, E18) are lighter; GSC masses, rounded.
    pub fn for_galileo(prn: u8) -> BoxWing {
        match prn {
            11 | 12 | 19 => BoxWing::galileo_iov(697.0),
            14 | 18 => BoxWing::galileo_foc(661.5),
            _ => BoxWing::galileo_foc(708.0),
        }
    }

    /// Acceleration (m/s², inertial) at position r with the Sun at `sun` (both geocentric, m), in
    /// sunlight fraction `lit`.
    pub fn accel(&self, r: V3, sun: V3, lit: f64) -> V3 {
        if lit <= 0.0 { return [0.0; 3]; }
        let to_sun = sub(sun, r);
        let dist = norm(to_sun);
        let e = scale(to_sun, 1.0 / dist);
        let p = P0 * (AU / dist).powi(2) * lit;
        // body axes: +Z to the Earth, Y perpendicular to the Sun, +X away from the Sun
        let zb = scale(unit(r), -1.0);
        let c = cross(zb, e);
        if norm(c) < 1e-9 { return [0.0; 3]; }                   // Sun along Z: yaw undefined, skip
        let mut yb = unit(c);
        let mut xb = cross(yb, zb);
        if dot(xb, e) > 0.0 { xb = scale(xb, -1.0); yb = scale(yb, -1.0); }
        let mut f = [0.0; 3];
        for s in &self.body {
            let n = add(add(scale(xb, s.n[0]), scale(yb, s.n[1])), scale(zb, s.n[2]));
            let cos = dot(n, e);
            if cos <= 0.0 { continue; }
            let k = (2.0 / 3.0) * s.d + 2.0 * s.r * cos + if s.mli { (2.0 / 3.0) * s.a } else { 0.0 };
            f = add(f, scale(add(scale(e, s.a + s.d), scale(n, k)), -p * s.area * cos));
        }
        for &(area, a, rr, d) in &self.wing {                   // normal along e, cos θ = 1
            f = add(f, scale(e, -p * area * (a + d + (2.0 / 3.0) * d + 2.0 * rr)));
        }
        scale(f, 1.0 / self.mass)
    }
}
