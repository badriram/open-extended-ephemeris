// SPDX-License-Identifier: MIT
//! SP3 orbit files (IGS ultra-rapid, CODE ultra-rapid and 5-day predictions): GPS and Galileo
//! positions and clocks, with the "predicted" flag set on the prediction part.
use crate::linalg::V3;
use crate::time::gps_seconds;
use std::collections::BTreeMap;

#[derive(Default, Clone)]
pub struct Track {
    pub ts: Vec<f64>,
    pub xs: Vec<V3>,      // ITRF, m
    pub cs: Vec<f64>,     // s, NaN where absent
    pub predicted: Vec<bool>,
    /// Built from broadcast orbits (antenna-referenced) rather than read from an IGS orbit.
    pub from_broadcast: bool,
}

pub fn read_sp3(text: &str) -> BTreeMap<u8, Track> {
    let mut out: BTreeMap<u8, Track> = BTreeMap::new();
    let mut t = f64::NAN;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('*') {
            let f: Vec<f64> = rest.split_whitespace().filter_map(|x| x.parse().ok()).collect();
            if f.len() >= 6 { t = gps_seconds(f[0] as i64, f[1] as i64, f[2] as i64, f[3] as i64, f[4] as i64, f[5]); }
        } else if (line.starts_with("PG") || line.starts_with("PE")) && !t.is_nan() && line.len() >= 60 {
            let v: Vec<f64> = (0..4).map(|k| line[4 + 14 * k..18 + 14 * k].trim().parse().unwrap_or(0.0)).collect();
            if v[0].abs() < 1e-6 { continue; }
            let prn: u8 = match line[2..4].trim().parse() { Ok(s) => s, Err(_) => continue };
            let sv = crate::rinex::sat_key(line.as_bytes()[1] as char, prn);
            let tr = out.entry(sv).or_default();
            tr.ts.push(t);
            tr.xs.push([v[0] * 1e3, v[1] * 1e3, v[2] * 1e3]);
            tr.cs.push(if v[3].abs() < 999999.0 { v[3] * 1e-6 } else { f64::NAN });
            tr.predicted.push(line.len() > 75 && &line[75..76] == "P");
        }
    }
    out
}
