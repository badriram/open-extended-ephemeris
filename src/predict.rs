// SPDX-License-Identifier: MIT
//! The prediction, per GPS satellite:
//!
//! 1. Take the newest IGS ultra-rapid orbit (24 h measured, 24 h predicted by the IGS) and its
//!    Earth rotation parameters.
//! 2. Fit the force model (orbit.rs) to that 48 h arc and integrate it a week ahead.
//! 3. Every 2 h, fit a broadcast ephemeris (lnav.rs) to 4 h of trajectory: to the IGS's own
//!    prediction while it lasts, then to ours.
//! 4. Calibrate against the satellite's real broadcasts over the measured day: the antenna offset
//!    (IGS orbits are centre of mass) and the clock reference. The clock is a straight line
//!    through the IGS clock estimates.
//! 5. Give each record an accuracy (URA) that grows with the prediction's age.
use crate::eo::{read_erp, EarthOrientation};
use crate::fetch::{brdc_urls, latest_igu, read_text};
use crate::linalg::*;
use crate::lnav;
use crate::orbit::{Orbit, Srp, State};
use crate::rinex::{read_nav, Nav, Navs};
use crate::sp3::{read_sp3, Track};
use crate::time::*;
use crate::{casic, rinex};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

const STEP: f64 = 7200.0;             // one ephemeris every 2 h, as the satellites do
const C_LIGHT: f64 = 299792458.0;
/// 95th-percentile range error by day of age in the backtest, as metres for the URA index.
const URA_BY_DAY: [f64; 7] = [1.9, 5.1, 8.6, 12.1, 17.2, 25.0, 35.0];
/// A warm-started force-model fit worse than this is redone from scratch.
const WARM_MAX_RMS: f64 = 6.0;

pub struct Options {
    pub out: String,
    pub sp3: Option<String>,
    pub erp: Option<String>,
    pub rinex: Vec<String>,
    pub days: f64,
    pub warm: Option<String>,
    pub max_check_m: f64,
    pub jobs: usize,
}

#[derive(Clone, Copy)]
pub struct FitState { pub t0: f64, pub state: State, pub srp: Srp, pub rms: f64 }

struct SatResult {
    records: Vec<(f64, Nav)>,
    fit: FitState,
    started: &'static str,
    iterations: usize,
    converged: bool,
    seconds: f64,
}

fn nearest(recs: &[Nav], t: f64) -> Option<&Nav> {
    recs.iter().min_by(|a, b| (a.toe_abs - t).abs().partial_cmp(&(b.toe_abs - t).abs()).unwrap())
        .filter(|r| (r.toe_abs - t).abs() <= 7200.0)
}

/// Radial, along-track and cross-track unit vectors (rows).
fn rac(x: V3, v: V3) -> M3 {
    let r = unit(x);
    let c = unit(cross(x, v));
    [r, cross(c, r), c]
}

fn predict_sat(sv: u8, tr: &Track, recs: &[Nav], eo: &EarthOrientation, t_end: f64, warm: Option<FitState>)
               -> Option<SatResult> {
    let start = std::time::Instant::now();
    let ts = &tr.ts;
    let n = ts.len();
    let pred_start = tr.predicted.iter().position(|&p| p).map(|k| ts[k]).unwrap_or(ts[n / 2]);
    let model = Orbit::new(eo, ts[0]);
    let xg: Vec<V3> = ts.iter().zip(&tr.xs).map(|(&t, &x)| mtv(&eo.c2t(t), x)).collect();

    let mut started = "cold";
    let mut fit = None;
    if let Some(w) = warm {
        // The previous run's state, carried forward to this arc's start with its own forces.
        let y0 = if w.t0 < ts[0] { Orbit::new(eo, w.t0).propagate_states(&w.state, &w.srp, &[ts[0]])[0] } else { w.state };
        let f = model.fit(ts, &xg, Some((y0, w.srp)));
        started = if f.2 > WARM_MAX_RMS || !f.4 { "warm-rejected" } else { "warm" };
        fit = Some(f);
    }
    if started != "warm" { fit = Some(model.fit(ts, &xg, None)); }
    let (state, srp, fit_rms, iterations, converged) = fit.unwrap();

    // Calibration over the measured day: where the broadcasts put the satellite and its clock,
    // relative to the IGS centre of mass and clock.
    let (mut d_rac, mut d_clk): (Vec<V3>, Vec<f64>) = (vec![], vec![]);
    for k in (1..n - 1).step_by(4) {
        if tr.predicted[k] { continue; }
        let Some(b) = nearest(recs, ts[k]) else { continue };
        let v = scale(sub(tr.xs[k + 1], tr.xs[k - 1]), 1.0 / (ts[k + 1] - ts[k - 1]));
        d_rac.push(mv(&rac(tr.xs[k], v), sub(lnav::sat_pos(&b.p, b.toe_abs, ts[k]), tr.xs[k])));
        if !tr.cs[k].is_nan() { d_clk.push(b.clock(ts[k]) - tr.cs[k]); }
    }
    if d_rac.len() < 10 { return None; }
    let off_rac: V3 = [0, 1, 2].map(|c| median(&mut d_rac.iter().map(|d| d[c]).collect()));
    let latest = recs.last()?;
    let ok: Vec<usize> = (0..n).filter(|&k| !tr.cs[k].is_nan()).collect();
    let (c0, c1) = if ok.len() >= 20 && d_clk.len() >= 10 {
        // least-squares line through the IGS clock estimates
        let m = ok.len() as f64;
        let (sx, sy) = ok.iter().fold((0.0, 0.0), |(a, b), &k| (a + ts[k] - pred_start, b + tr.cs[k]));
        let (mx, my) = (sx / m, sy / m);
        let (sxy, sxx) = ok.iter().fold((0.0, 0.0), |(a, b), &k| {
            let dx = ts[k] - pred_start - mx;
            (a + dx * (tr.cs[k] - my), b + dx * dx)
        });
        let slope = sxy / sxx;
        (my - slope * mx + median(&mut d_clk), slope)
    } else {
        (latest.clock(pred_start), latest.af1)                    // no IGS clock: extend the latest broadcast's
    };

    let grid: Vec<f64> = (0..).map(|i| ts[0] + 300.0 * i as f64).take_while(|&t| t <= t_end + STEP + 1.0).collect();
    let xe: Vec<V3> = model.propagate(&state, &srp, &grid).iter().zip(&grid).map(|(x, &t)| mv(&eo.c2t(t), *x)).collect();

    let mut records = Vec::new();
    let mut toe = (pred_start / STEP).ceil() * STEP;
    while toe <= t_end {
        let win: Vec<f64> = (0..49).map(|i| toe - 7200.0 + 300.0 * i as f64).collect();
        let pos: Vec<V3> = if *win.last().unwrap() <= ts[n - 1] {
            win.iter().map(|&t| lagrange(ts, &tr.xs, t, 10)).collect()       // the IGS's own orbit
        } else {
            win.iter().map(|&t| xe[((t - grid[0]) / 300.0).round() as usize]).collect()
        };
        // move the centre of mass to where the broadcasts put the satellite
        let pos: Vec<V3> = (0..pos.len()).map(|i| {
            let v = sub(pos[(i + 1).min(pos.len() - 1)], pos[i.saturating_sub(1)]);
            add(pos[i], mtv(&rac(pos[i], v), off_rac))
        }).collect();
        let (p, _rms) = lnav::fit(&win, &pos, toe);
        let age_days = ((toe - pred_start) / 86400.0).max(0.0) as usize;
        let iodc = 768.0 + ((toe / STEP) as i64 % 240) as f64;
        records.push((toe, Nav {
            p, toe_abs: toe, toc_abs: toe, af0: c0 + c1 * (toe - pred_start), af1: c1, af2: 0.0,
            iode: (iodc as i64 & 0xFF) as f64, iodc, svacc: URA_BY_DAY[age_days.min(URA_BY_DAY.len() - 1)],
            health: 0.0, tgd: latest.tgd, ttx: (toe - 7200.0).rem_euclid(WEEK),
        }));
        toe += STEP;
    }
    let _ = sv;
    Some(SatResult {
        records, fit: FitState { t0: ts[0], state, srp, rms: fit_rms }, started, iterations, converged,
        seconds: start.elapsed().as_secs_f64(),
    })
}

/// Compare the first predicted records with the satellites' real broadcasts where both exist.
fn self_check(results: &BTreeMap<u8, SatResult>, navs: &Navs) -> Vec<f64> {
    let mut errs = vec![];
    for (sv, res) in results {
        for (toe, rec) in res.records.iter().take(3) {
            for k in 0..8 {
                let t = toe - 3600.0 + 900.0 * k as f64;
                let Some(b) = nearest(&navs[sv], t) else { continue };
                errs.push(sisre(rec, b, t).0);
            }
        }
    }
    errs
}

/// Signal-in-space range error and 3D orbit difference of a record against a reference record.
pub fn sisre(r: &Nav, b: &Nav, t: f64) -> (f64, f64) {
    let xb = lnav::sat_pos(&b.p, b.toe_abs, t);
    let d = sub(lnav::sat_pos(&r.p, r.toe_abs, t), xb);
    let dr = dot(d, unit(xb));
    let dc = C_LIGHT * (r.clock(t) - b.clock(t));
    (((0.98 * dr - dc).powi(2) + (dot(d, d) - dr * dr) / 49.0).sqrt(), norm(d))
}

fn read_warm(src: &str) -> BTreeMap<u8, FitState> {
    let mut out = BTreeMap::new();
    let Ok(text) = read_text(src) else { eprintln!("warm start unavailable: {src}"); return out };
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let v: Vec<f64> = line.split_whitespace().filter_map(|x| x.parse().ok()).collect();
        if v.len() == 14 {
            out.insert(v[0] as u8, FitState { t0: v[1], state: [v[2], v[3], v[4], v[5], v[6], v[7]],
                                              srp: [v[8], v[9], v[10], v[11], v[12]], rms: v[13] });
        }
    }
    out
}

fn doy(y: i64, m: i64, d: i64) -> i64 { days_from_civil(y, m, d) - days_from_civil(y, 1, 1) + 1 }

pub fn run(o: &Options) -> Result<(), String> {
    let t_run = std::time::Instant::now();
    let now = now_gps();
    let (sp3_src, erp_src) = match (&o.sp3, &o.erp) {
        (Some(s), Some(e)) => (s.clone(), e.clone()),
        (None, None) => latest_igu(now)?,
        _ => return Err("--sp3 and --erp go together".into()),
    };
    eprintln!("orbit: {sp3_src}\nERP:   {erp_src}");
    let sp3 = read_sp3(&read_text(&sp3_src)?);
    let erp = read_erp(&read_text(&erp_src)?)?;
    if sp3.is_empty() { return Err("no GPS satellites in the orbit file".into()); }
    let t0 = sp3.values().map(|t| t.ts[0]).fold(f64::INFINITY, f64::min);
    let t_arc_end = sp3.values().map(|t| *t.ts.last().unwrap()).fold(f64::NEG_INFINITY, f64::max);

    let mut navs = Navs::new();
    let mut used = vec![];
    if !o.rinex.is_empty() {
        for f in &o.rinex { read_nav(&read_text(f)?, &mut navs)?; used.push(f.clone()); }
    } else {
        let (mut day, last) = (((t0 + UTC_GPS + UNIX_GPS0) / 86400.0).floor() as i64,
                               ((now + UTC_GPS + UNIX_GPS0) / 86400.0).floor() as i64);
        while day <= last {
            let (y, m, d, ..) = calendar(day as f64 * 86400.0 - UNIX_GPS0);
            for url in brdc_urls(y, doy(y, m, d)) {
                if let Ok(text) = read_text(&url) {
                    if read_nav(&text, &mut navs).is_ok() { used.push(url); break; }
                }
            }
            day += 1;
        }
    }
    eprintln!("broadcasts:\n  {}", used.join("\n  "));
    let warm = o.warm.as_deref().map(read_warm).unwrap_or_default();

    let pred_start = sp3.values().filter_map(|t| t.predicted.iter().position(|&p| p).map(|k| t.ts[k]))
        .fold(f64::INFINITY, f64::min);
    let t_end = pred_start + o.days * 86400.0;
    let warm_t0 = warm.values().map(|w| w.t0).fold(t0, f64::min);
    let eo = EarthOrientation::new(&erp, warm_t0, t_end + STEP);
    let jobs: Vec<(u8, &Track)> = sp3.iter()
        .filter(|(sv, _)| navs.get(sv).and_then(|r| r.last()).map_or(false, |r| r.health == 0.0))
        .map(|(sv, t)| (*sv, t)).collect();

    let next = AtomicUsize::new(0);
    let results: Mutex<BTreeMap<u8, SatResult>> = Mutex::new(BTreeMap::new());
    std::thread::scope(|s| {
        for _ in 0..o.jobs.max(1) {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= jobs.len() { break; }
                let (sv, tr) = jobs[i];
                if let Some(r) = predict_sat(sv, tr, &navs[&sv], &eo, t_end, warm.get(&sv).copied()) {
                    results.lock().unwrap().insert(sv, r);
                }
            });
        }
    });
    let results = results.into_inner().unwrap();
    if results.is_empty() { return Err("no satellite could be predicted".into()); }

    let mut entries: Vec<(u8, Nav)> = results.iter().flat_map(|(sv, r)| r.records.iter().map(|(_, n)| (*sv, n.clone()))).collect();
    entries.sort_by(|a, b| a.1.toc_abs.partial_cmp(&b.1.toc_abs).unwrap().then(a.0.cmp(&b.0)));
    let mut blob = Vec::new();
    for (sv, n) in &entries { blob.extend(casic::gpseph_frame(n, *sv)?); }
    casic::check_frames(&blob)?;
    let mut errs = self_check(&results, &navs);
    let check = if errs.is_empty() { (0, f64::NAN, f64::NAN) } else {
        (errs.len(), median(&mut errs.clone()), percentile(&mut errs, 95.0))
    };
    let mut fits: Vec<f64> = results.values().map(|r| r.fit.rms).collect();
    let worst_fit = fits.iter().cloned().fold(0.0, f64::max);
    let count = |s: &str| results.values().filter(|r| r.started == s).count();
    let mut iters: Vec<f64> = results.values().map(|r| r.iterations as f64).collect();
    let mut secs: Vec<f64> = results.values().map(|r| r.seconds).collect();
    let base = |s: &str| s.rsplit('/').next().unwrap_or(s).to_string();

    std::fs::create_dir_all(&o.out).map_err(|e| e.to_string())?;
    let path = |f: &str| format!("{}/{f}", o.out.trim_end_matches('/'));
    std::fs::write(path("gps-7d-casic.bin"), &blob).map_err(|e| e.to_string())?;
    std::fs::write(path("gps-7d.rnx"), rinex::write_nav(&entries, &[
        "PREDICTED broadcast ephemeris (extended ephemeris), not".into(),
        "broadcast by the satellites. From IGS ultra-rapid orbit".into(),
        base(&sp3_src),
        "github.com/badriram/open-extended-ephemeris  MIT".into(),
    ])).map_err(|e| e.to_string())?;
    let mut state = String::from("# sv t0 x y z vx vy vz (GCRS m, m/s) D0 Y0 B0 Bc Bs (nm/s^2) rms (m)\n");
    for (sv, r) in &results {
        let f = &r.fit;
        state += &format!("{sv} {:.3} {} {} {:.6}\n", f.t0,
            f.state.iter().map(|v| format!("{v:.9}")).collect::<Vec<_>>().join(" "),
            f.srp.iter().map(|v| format!("{v:.9}")).collect::<Vec<_>>().join(" "), f.rms);
    }
    std::fs::write(path("fit-state.txt"), state).map_err(|e| e.to_string())?;
    let sats: Vec<String> = results.keys().map(|s| s.to_string()).collect();
    let json = format!(r#"{{
  "made": "{made}",
  "igs_orbit": "{sp3_src}",
  "igs_erp": "{erp_src}",
  "broadcasts": [{bc}],
  "prediction_starts": "{ps}",
  "valid_until": "{vu}",
  "igs_prediction_until": "{ipu}",
  "satellites": [{sats}],
  "records": {nrec},
  "casic_bytes": {nb},
  "force_model_fit_m": {{"median": {fm:.4}, "worst": {fw:.4}}},
  "fits": {{"cold": {cold}, "warm": {warm_}, "warm_rejected": {rej}, "not_converged": {nc}, "iterations_median": {im}}},
  "satellite_seconds": {{"median": {sm:.2}, "total": {st:.1}}},
  "run_seconds": {rs:.1},
  "self_check_vs_broadcast": {{"samples": {cn}, "median_m": {cm:.4}, "p95_m": {cp:.4}}},
  "backtest": "range error median/95%: 0.3/1.1 m first half day, 1.1/5.1 m at 2 days, 4.7/17 m at 4.5-5 days, worst 25 m (IGS orbit of 2026-10-13 06:00)"
}}
"#, made = iso_utc(now), bc = used.iter().map(|u| format!("\"{u}\"")).collect::<Vec<_>>().join(", "),
        ps = iso_utc(pred_start), vu = iso_utc(t_end), ipu = iso_utc(t_arc_end), sats = sats.join(", "),
        nrec = entries.len(), nb = blob.len(), fm = median(&mut fits), fw = worst_fit,
        cold = count("cold"), warm_ = count("warm"), rej = count("warm-rejected"),
        nc = results.values().filter(|r| !r.converged).count(), im = median(&mut iters),
        sm = median(&mut secs.clone()), st = secs.iter().sum::<f64>(), rs = t_run.elapsed().as_secs_f64(),
        cn = check.0, cm = check.1, cp = check.2);
    std::fs::write(path("gps-7d.json"), &json).map_err(|e| e.to_string())?;
    print!("{json}");
    if check.0 == 0 || !(check.1 <= o.max_check_m) {
        return Err(format!("self-check failed: median {:.3} m over {} samples", check.1, check.0));
    }
    let _ = &mut secs;
    Ok(())
}
