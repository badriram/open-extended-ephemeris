// SPDX-License-Identifier: MIT
//! Small fixed-size vectors and matrices, and a dense solver for the fits' normal equations.

pub type V3 = [f64; 3];
pub type M3 = [[f64; 3]; 3];

pub fn add(a: V3, b: V3) -> V3 { [a[0] + b[0], a[1] + b[1], a[2] + b[2]] }
pub fn sub(a: V3, b: V3) -> V3 { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
pub fn scale(a: V3, s: f64) -> V3 { [a[0] * s, a[1] * s, a[2] * s] }
pub fn dot(a: V3, b: V3) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
pub fn norm(a: V3) -> f64 { dot(a, a).sqrt() }
pub fn unit(a: V3) -> V3 { scale(a, 1.0 / norm(a)) }
pub fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
/// m · v
pub fn mv(m: &M3, v: V3) -> V3 { [dot(m[0], v), dot(m[1], v), dot(m[2], v)] }
/// mᵀ · v
pub fn mtv(m: &M3, v: V3) -> V3 {
    [m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
     m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
     m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2]]
}
pub fn mm(a: &M3, b: &M3) -> M3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 { for j in 0..3 { c[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum(); } }
    c
}

/// Solve the n×n system a·x = b (row-major a) by Gaussian elimination with partial pivoting.
pub fn solve(n: usize, a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut m = a.to_vec();
    let mut x = b.to_vec();
    for col in 0..n {
        let piv = (col..n).max_by(|&i, &j| m[i * n + col].abs().partial_cmp(&m[j * n + col].abs()).unwrap()).unwrap();
        if piv != col {
            for k in 0..n { m.swap(col * n + k, piv * n + k); }
            x.swap(col, piv);
        }
        let d = m[col * n + col];
        for row in col + 1..n {
            let f = m[row * n + col] / d;
            if f != 0.0 {
                for k in col..n { m[row * n + k] -= f * m[col * n + k]; }
                x[row] -= f * x[col];
            }
        }
    }
    for col in (0..n).rev() {
        let s: f64 = (col + 1..n).map(|k| m[col * n + k] * x[k]).sum();
        x[col] = (x[col] - s) / m[col * n + col];
    }
    x
}

pub fn median(v: &mut Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 { v[n / 2] } else { 0.5 * (v[n / 2 - 1] + v[n / 2]) }
}

pub fn percentile(v: &mut Vec<f64>, p: f64) -> f64 {
    // linear interpolation between closest ranks, as numpy's default
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let r = p / 100.0 * (v.len() - 1) as f64;
    let (lo, f) = (r.floor() as usize, r - r.floor());
    if lo + 1 < v.len() { v[lo] * (1.0 - f) + v[lo + 1] * f } else { v[lo] }
}

/// Lagrange interpolation of order `order` around t in the sorted samples (ts, ys).
pub fn lagrange(ts: &[f64], ys: &[V3], t: f64, order: usize) -> V3 {
    let i = ts.partition_point(|&x| x < t);
    let lo = (i as isize - (order / 2) as isize).clamp(0, (ts.len() - order) as isize) as usize;
    let mut out = [0.0; 3];
    for j in 0..order {
        let mut w = 1.0;
        for m in 0..order {
            if m != j { w *= (t - ts[lo + m]) / (ts[lo + j] - ts[lo + m]); }
        }
        out = add(out, scale(ys[lo + j], w));
    }
    out
}
