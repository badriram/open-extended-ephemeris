// SPDX-License-Identifier: MIT
//! Earth's gravity field: EGM2008 to degree and order 4 (degree 6 changed the fit by under
//! 0.1 m), evaluated with the recursion of Montenbruck & Gill, *Satellite Orbits*, eq. 3.33.
use crate::linalg::V3;

pub const GM: f64 = 3.986004418e14;
pub const R_E: f64 = 6378136.3;
const NMAX: usize = 4;

// fully normalised (n, m, C, S)
const CS: [(usize, usize, f64, f64); 12] = [
    (2, 0, -0.484165143790815e-3, 0.0),
    (2, 1, -0.206615509074176e-9, 0.138441389137979e-8),
    (2, 2, 0.243938357328313e-5, -0.140027370385934e-5),
    (3, 0, 0.957161207093473e-6, 0.0),
    (3, 1, 0.203046201047864e-5, 0.248200415856872e-6),
    (3, 2, 0.904787894809528e-6, -0.619005475177618e-6),
    (3, 3, 0.721321757121568e-6, 0.141434926192941e-5),
    (4, 0, 0.539965866638991e-6, 0.0),
    (4, 1, -0.536157389388867e-6, -0.473567346518086e-6),
    (4, 2, 0.350501623962649e-6, 0.662480026275829e-6),
    (4, 3, 0.990856766672321e-6, -0.200956723567452e-6),
    (4, 4, -0.188519633023033e-6, 0.308803882149194e-6),
];

fn factorial(n: usize) -> f64 { (1..=n).map(|k| k as f64).product() }

/// Unnormalised coefficients, C[n][m] and S[n][m], C[0][0] = 1.
fn coefficients() -> ([[f64; NMAX + 1]; NMAX + 1], [[f64; NMAX + 1]; NMAX + 1]) {
    let mut c = [[0.0; NMAX + 1]; NMAX + 1];
    let mut s = [[0.0; NMAX + 1]; NMAX + 1];
    c[0][0] = 1.0;
    for &(n, m, cb, sb) in CS.iter() {
        let f = ((if m == 0 { 1.0 } else { 2.0 }) * (2 * n + 1) as f64 * factorial(n - m) / factorial(n + m)).sqrt();
        c[n][m] = cb * f;
        s[n][m] = sb * f;
    }
    (c, s)
}

thread_local! {
    static COEF: ([[f64; NMAX + 1]; NMAX + 1], [[f64; NMAX + 1]; NMAX + 1]) = coefficients();
}

/// Acceleration from the geopotential at an Earth-fixed position (m), m/s².
#[allow(dead_code)]
pub fn gravity(r: V3) -> V3 { gravity_tide(r, &[0.0; 3], &[0.0; 3]) }

/// Solid Earth tide, IERS Conventions 2010 section 6.2.1, step 1 (frequency-independent, elastic
/// Love numbers): the Sun and Moon raise a tide in the Earth that changes the fully normalised
/// degree-2 coefficients by dC2m - i dS2m = k2m/5 * sum (GMj/GM)(R/rj)^3 P2m(sin phi_j) e^(-i m lambda_j).
/// EGM2008 is a tide-free model, so the whole change is added. Bodies are Earth-fixed positions (m).
/// Returns (dC20, dC21, dC22), (0, dS21, dS22), normalised.
pub fn solid_tide(bodies: &[(V3, f64)]) -> ([f64; 3], [f64; 3]) {
    const K2: [f64; 3] = [0.29525, 0.29470, 0.29801];
    let (mut dc, mut ds) = ([0.0; 3], [0.0; 3]);
    for &(b, gm) in bodies {
        let rb = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
        let (sp, cp) = (b[2] / rb, (b[0] * b[0] + b[1] * b[1]).sqrt() / rb);
        let lam = b[1].atan2(b[0]);
        let f = gm / GM * (R_E / rb).powi(3);
        let p = [5f64.sqrt() * (3.0 * sp * sp - 1.0) / 2.0, 15f64.sqrt() * sp * cp, 15f64.sqrt() / 2.0 * cp * cp];
        for m in 0..3 {
            dc[m] += K2[m] / 5.0 * f * p[m] * (m as f64 * lam).cos();
            ds[m] += K2[m] / 5.0 * f * p[m] * (m as f64 * lam).sin();
        }
    }
    (dc, ds)
}

/// The geopotential with the degree-2 coefficients changed by (dc, ds), fully normalised.
pub fn gravity_tide(r: V3, dc: &[f64; 3], ds: &[f64; 3]) -> V3 {
    COEF.with(|(cc0, ss0)| {
        let (mut cc, mut ss) = (*cc0, *ss0);
        for m in 0..3 {
            let f = ((if m == 0 { 1.0 } else { 2.0 }) * 5.0 * factorial(2 - m) / factorial(2 + m)).sqrt();
            cc[2][m] += dc[m] * f;
            ss[2][m] += ds[m] * f;
        }
        const N: usize = NMAX + 2;
        let (x, y, z) = (r[0], r[1], r[2]);
        let r2 = x * x + y * y + z * z;
        let rho = R_E * R_E / r2;
        let (x0, y0, z0) = (R_E * x / r2, R_E * y / r2, R_E * z / r2);
        let mut v = [[0.0; N]; N];
        let mut w = [[0.0; N]; N];
        v[0][0] = R_E / r2.sqrt();
        v[1][0] = z0 * v[0][0];
        for n in 2..N {
            v[n][0] = ((2 * n - 1) as f64 * z0 * v[n - 1][0] - (n - 1) as f64 * rho * v[n - 2][0]) / n as f64;
        }
        for m in 1..N {
            let k = (2 * m - 1) as f64;
            v[m][m] = k * (x0 * v[m - 1][m - 1] - y0 * w[m - 1][m - 1]);
            w[m][m] = k * (x0 * w[m - 1][m - 1] + y0 * v[m - 1][m - 1]);
            if m + 1 < N {
                v[m + 1][m] = (2 * m + 1) as f64 * z0 * v[m][m];
                w[m + 1][m] = (2 * m + 1) as f64 * z0 * w[m][m];
            }
            for n in m + 2..N {
                let (a, b) = ((2 * n - 1) as f64, (n + m - 1) as f64);
                v[n][m] = (a * z0 * v[n - 1][m] - b * rho * v[n - 2][m]) / (n - m) as f64;
                w[n][m] = (a * z0 * w[n - 1][m] - b * rho * w[n - 2][m]) / (n - m) as f64;
            }
        }
        let (mut ax, mut ay, mut az) = (0.0, 0.0, 0.0);
        for n in 0..=NMAX {
            for m in 0..=n {
                let (c, s) = (cc[n][m], ss[n][m]);
                if m == 0 {
                    ax -= c * v[n + 1][1];
                    ay -= c * w[n + 1][1];
                    az -= (n + 1) as f64 * c * v[n + 1][0];
                } else {
                    let f = 0.5 * ((n - m + 1) * (n - m + 2)) as f64;
                    ax += 0.5 * (-c * v[n + 1][m + 1] - s * w[n + 1][m + 1]) + f * (c * v[n + 1][m - 1] + s * w[n + 1][m - 1]);
                    ay += 0.5 * (-c * w[n + 1][m + 1] + s * v[n + 1][m + 1]) + f * (-c * w[n + 1][m - 1] + s * v[n + 1][m - 1]);
                    az += (n - m + 1) as f64 * (-c * v[n + 1][m] - s * w[n + 1][m]);
                }
            }
        }
        let k = GM / (R_E * R_E);
        [k * ax, k * ay, k * az]
    })
}
