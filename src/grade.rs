// SPDX-License-Identifier: MIT
//! Grade a predicted RINEX file against what the satellites actually broadcast: for every
//! predicted record, over its central two hours at 15-minute steps, compare with the broadcast
//! record a receiver would have used then. Reported by age (from the first predicted record):
//! signal-in-space range error (median, 95th percentile, worst), median orbit error, and the
//! worst satellites.
use crate::fetch::read_text;
use crate::linalg::{median, percentile};
use crate::predict::sisre;
use crate::rinex::{read_nav, Navs};

pub fn run(pred: &str, truth: &[String]) -> Result<(), String> {
    // "--system E" anywhere among the arguments grades one constellation
    let (mut sys, mut files) = (None, vec![]);
    let mut it = truth.iter();
    while let Some(a) = it.next() {
        if a == "--system" { sys = it.next().and_then(|s| s.chars().next()); } else { files.push(a.clone()); }
    }
    let truth = &files[..];
    let mut p = Navs::new();
    read_nav(&read_text(pred)?, &mut p)?;
    if let Some(c) = sys { p.retain(|sv, _| crate::rinex::sat_sys(*sv) == c); }
    let mut t = Navs::new();
    for f in truth { read_nav(&read_text(f)?, &mut t)?; }
    let t0 = p.values().flatten().map(|r| r.toe_abs).fold(f64::INFINITY, f64::min) - 7200.0;
    let mut rows: Vec<(f64, f64, f64, f64, bool)> = vec![];
    let mut worst: Vec<(u8, f64)> = vec![];
    for (sv, recs) in &p {
        let Some(tb) = t.get(sv) else { continue };
        let mut w = 0.0f64;
        for r in recs {
            for k in 0..8 {
                let tt = r.toe_abs - 3600.0 + 900.0 * k as f64;
                let Some(b) = tb.iter().filter(|b| (b.toe_abs - tt).abs() <= 7200.0)
                    .min_by(|a, b| (a.toe_abs - tt).abs().partial_cmp(&(b.toe_abs - tt).abs()).unwrap()) else { continue };
                let (s, d3) = sisre(r, b, tt);
                let dclk = 299792458.0 * (r.clock(tt) - b.clock(tt));
                rows.push(((tt - t0) / 86400.0, s, d3, dclk.abs(), s > r.svacc));
                w = w.max(s);
            }
        }
        worst.push((*sv, w));
    }
    if rows.is_empty() { return Err("no overlap between the prediction and the broadcasts".into()); }
    println!("days      | range error median / 95% / worst | orbit 3D median | clock |median| | over stated accuracy");
    let max_d = rows.iter().map(|r| r.0).fold(0.0, f64::max);
    let mut d = 0.0;
    while d < max_d {
        let sel: Vec<&(f64, f64, f64, f64, bool)> = rows.iter().filter(|r| r.0 >= d && r.0 < d + 0.5).collect();
        if !sel.is_empty() {
            let mut s: Vec<f64> = sel.iter().map(|r| r.1).collect();
            let mut o: Vec<f64> = sel.iter().map(|r| r.2).collect();
            let mut ck: Vec<f64> = sel.iter().map(|r| r.3).collect();
            let mx = s.iter().cloned().fold(0.0, f64::max);
            let over = 100.0 * sel.iter().filter(|r| r.4).count() as f64 / sel.len() as f64;
            println!("{d:3.1}-{:3.1}   | {:6.2} / {:6.2} / {:6.1} m        | {:6.2} m        | {:6.2} m       | {over:5.1}%", d + 0.5,
                     median(&mut s.clone()), percentile(&mut s, 95.0), mx, median(&mut o), median(&mut ck));
        }
        d += 0.5;
    }
    worst.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    println!("worst satellites: {}", worst.iter().take(5).map(|(s, v)| format!("{} {v:.1} m", crate::rinex::sat_name(*s))).collect::<Vec<_>>().join(", "));
    Ok(())
}
