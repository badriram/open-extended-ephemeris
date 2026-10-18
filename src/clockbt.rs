// SPDX-License-Identifier: MIT
//! clock-backtest: how well each clock model predicts the satellites' own broadcast clocks, over
//! as much broadcast history as given (a year of daily files, say). For every epoch (every
//! `every` hours) and satellite, fit the model to the broadcasts of the preceding days, then
//! compare hourly over the next seven days with the clock the satellite actually broadcast.
//! Only healthy records are used. Errors are in metres of range (c × clock error).
use crate::fetch::read_text;
use crate::linalg::{median, percentile};
use crate::predict::fit_clock_poly;
use crate::rinex::{read_nav, Navs};
use crate::time::{calendar, UTC_GPS};

const C: f64 = 299792458.0;

#[derive(Clone, Copy)]
struct Model { quad: bool, days: f64 }

fn name(m: &Model) -> String { format!("{}:{}", if m.quad { "quad" } else { "line" }, m.days) }

pub fn run(args: &[String]) -> Result<(), String> {
    let mut files = vec![];
    let mut models = vec![];
    let mut every = 6.0;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--models" => {
                for spec in args.get(i + 1).ok_or("--models needs a value")?.split(',') {
                    let (kind, d) = spec.split_once(':').ok_or("model as quad:DAYS or line:DAYS")?;
                    models.push(Model { quad: kind == "quad", days: d.parse().map_err(|_| "model days")? });
                }
                i += 2;
            }
            "--every" => { every = args.get(i + 1).ok_or("--every needs hours")?.parse().map_err(|_| "--every: hours")?; i += 2; }
            f => { files.push(f.to_string()); i += 1; }
        }
    }
    if models.is_empty() {
        models = [(false, 1.0), (false, 2.0), (true, 3.0), (true, 4.0), (true, 5.0), (true, 6.0), (true, 7.0), (true, 10.0)]
            .iter().map(|&(quad, days)| Model { quad, days }).collect();
    }
    let mut navs = Navs::new();
    for (k, f) in files.iter().enumerate() {
        read_nav(&read_text(f)?, &mut navs)?;
        if k % 50 == 49 { eprintln!("read {} of {} files", k + 1, files.len()); }
    }
    // per satellite: healthy records sorted by toe, and (toc, af0) samples
    let sats: Vec<(u8, Vec<crate::rinex::Nav>)> = navs.into_iter()
        .map(|(sv, v)| (sv, v.into_iter().filter(|r| r.health == 0.0).collect::<Vec<_>>()))
        .filter(|(_, v)| v.len() > 50).collect();
    let t_first = sats.iter().map(|(_, v)| v[0].toc_abs).fold(f64::INFINITY, f64::min);
    let t_last = sats.iter().map(|(_, v)| v.last().unwrap().toc_abs).fold(f64::NEG_INFINITY, f64::max);
    let max_days = models.iter().map(|m| m.days).fold(0.0, f64::max);
    let epochs: Vec<f64> = (0..).map(|k| (t_first / 3600.0).ceil() * 3600.0 + max_days * 86400.0 + k as f64 * every * 3600.0)
        .take_while(|&t| t + 7.0 * 86400.0 <= t_last).collect();
    eprintln!("{} satellites, {} epochs every {every} h, {} models", sats.len(), epochs.len(), models.len());

    // errors[model][day] and, for the first model listed with quad:5 if present, per month and PRN
    let mut errs: Vec<Vec<Vec<f32>>> = vec![vec![vec![]; 7]; models.len()];
    let focus = models.iter().position(|m| m.quad && m.days == 5.0).unwrap_or(0);
    let mut by_month: Vec<Vec<f32>> = vec![vec![]; 13];
    let mut by_sat: std::collections::BTreeMap<u8, Vec<f32>> = Default::default();
    let mut big = vec![0usize; models.len()];
    for (sv, recs) in &sats {
        let tocs: Vec<f64> = recs.iter().map(|r| r.toc_abs).collect();
        let toes: Vec<f64> = recs.iter().map(|r| r.toe_abs).collect();
        for &t0 in &epochs {
            // the broadcast in force at each future hour: nearest toe within 2 h
            let truth: Vec<Option<(f64, f64)>> = (1..=168).map(|h| {
                let t = t0 + h as f64 * 3600.0;
                let k = toes.partition_point(|&x| x < t);
                let cand = [k.wrapping_sub(1), k].into_iter().filter(|&j| j < toes.len())
                    .min_by(|&a, &b| (toes[a] - t).abs().partial_cmp(&(toes[b] - t).abs()).unwrap());
                cand.filter(|&j| (toes[j] - t).abs() <= 7200.0).map(|j| (t, recs[j].clock(t)))
            }).collect();
            if truth.iter().all(|x| x.is_none()) { continue; }
            for (mi, m) in models.iter().enumerate() {
                let lo = tocs.partition_point(|&x| x < t0 - m.days * 86400.0);
                let hi = tocs.partition_point(|&x| x <= t0);
                let pts: Vec<(f64, f64)> = (lo..hi).map(|j| (tocs[j] - t0, recs[j].af0)).collect();
                if pts.len() < 12 || pts.last().map_or(true, |p| p.0 < -6.0 * 3600.0) { continue; }
                let Some(c) = fit_clock_poly(&pts, m.quad) else { continue };
                for (h, tr) in truth.iter().enumerate() {
                    let Some((t, truth_clock)) = tr else { continue };
                    let dt = t - t0;
                    let e = (C * (c[0] + c[1] * dt + c[2] * dt * dt - truth_clock)).abs() as f32;
                    errs[mi][h / 24].push(e);
                    if e > 100.0 { big[mi] += 1; }
                    if mi == focus && (96..120).contains(&h) {
                        by_month[calendar(t0 + UTC_GPS).1 as usize].push(e);
                        by_sat.entry(*sv).or_default().push(e);
                    }
                }
            }
        }
    }
    let stat = |v: &[f32]| -> (f64, f64, f64, f64) {
        let mut x: Vec<f64> = v.iter().map(|&e| e as f64).collect();
        if x.is_empty() { return (f64::NAN, f64::NAN, f64::NAN, f64::NAN); }
        let mx = x.iter().cloned().fold(0.0, f64::max);
        (median(&mut x.clone()), percentile(&mut x.clone(), 95.0), percentile(&mut x, 99.0), mx)
    };
    println!("clock prediction error against broadcast clocks, metres (median / 95% / 99%), by forecast day");
    print!("{:<9}", "model");
    for d in 0..7 { print!(" | day {}-{}            ", d, d + 1); }
    println!(" | >100 m");
    for (mi, m) in models.iter().enumerate() {
        print!("{:<9}", name(m));
        for d in 0..7 {
            let (md, p95, p99, _) = stat(&errs[mi][d]);
            print!(" | {md:5.2} /{p95:6.2} /{p99:6.1}");
        }
        let n: usize = errs[mi].iter().map(|v| v.len()).sum();
        println!(" | {:.3}%", 100.0 * big[mi] as f64 / n.max(1) as f64);
    }
    println!("\n{} at day 4-5, by month of the prediction: median / 95% / 99% / worst (m)", name(&models[focus]));
    for (mo, v) in by_month.iter().enumerate().skip(1) {
        if v.is_empty() { continue; }
        let (md, p95, p99, mx) = stat(v);
        println!("  {mo:2}: {md:5.2} / {p95:6.2} / {p99:6.1} / {mx:7.1}");
    }
    println!("\n{} at day 4-5, by satellite (worst 95th percentiles first): median / 95% / worst (m)", name(&models[focus]));
    let mut rows: Vec<(f64, u8, f64, f64)> = by_sat.iter().map(|(sv, v)| { let (md, p95, _, mx) = stat(v); (p95, *sv, md, mx) }).collect();
    rows.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    for (p95, sv, md, mx) in rows.iter().take(12) { println!("  {}: {md:5.2} / {p95:6.2} / {mx:7.1}", crate::rinex::sat_name(*sv)); }
    println!("  ... best: {}", rows.iter().rev().take(5).map(|r| format!("{} {:.2}", crate::rinex::sat_name(r.1), r.0)).collect::<Vec<_>>().join(", "));
    Ok(())
}
