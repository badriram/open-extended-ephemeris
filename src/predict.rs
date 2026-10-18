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
use crate::eo::{extrapolate, join_bulletin_a, last_lod, read_bulletin_a, read_erp, EarthOrientation};
use crate::fetch::{brdc_urls, latest_igu, read_text};
use crate::linalg::*;
use crate::lnav;
use crate::orbit::{Forces, Orbit, Srp, State};
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
const URA_BY_DAY: [f64; 7] = [1.6, 4.3, 7.0, 9.8, 12.3, 15.0, 18.0];
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
    pub forces: Forces,
    pub clock: ClockModel,
    /// IERS Bulletin A (file or URL) for Earth rotation beyond the IGS values; "none" to
    /// extrapolate the IGS values with their last length of day instead.
    pub eop: String,
}

pub const BULLETIN_A: &str = "https://datacenter.iers.org/data/latestVersion/bulletinA.txt";

/// How the satellite clock is projected forward from the IGS clock estimates.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ClockModel {
    /// A straight line through 48 h.
    Line,
    /// A line plus the four harmonics of the orbital period (12, 6, 4, 3 h) that GPS clocks show,
    /// all fitted over 48 h and all carried forward.
    Harmonic,
    /// As Harmonic, but the harmonics are fitted and then left out of the prediction: they keep
    /// the drift estimate clean without extrapolating terms whose amplitude changes.
    HarmonicDrift,
    /// A straight line through the satellites' own broadcast clocks over the last N days.
    BroadcastLine(f64),
    /// A quadratic through the broadcast clocks over the last N days (for ageing clocks).
    BroadcastQuad(f64),
}

/// The IGS issues an ultra-rapid orbit about three hours after its prediction starts; a run
/// sees broadcasts up to about then and no later.
const LATENCY: f64 = 4.0 * 3600.0;

#[derive(Clone, Copy)]
pub struct FitState { pub t0: f64, pub state: State, pub srp: Srp, pub rms: f64 }

struct SatResult {
    records: Vec<(f64, Nav)>,
    fit: FitState,
    clocks: Vec<(f64, f64)>,
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

fn predict_sat(sv: u8, tr: &Track, recs: &[Nav], eo: &EarthOrientation, t_end: f64, warm: Option<FitState>,
               forces: Forces, clock_model: ClockModel, history: &[(f64, f64)]) -> Option<SatResult> {
    let start = std::time::Instant::now();
    let ts = &tr.ts;
    let n = ts.len();
    let pred_start = tr.predicted.iter().position(|&p| p).map(|k| ts[k]).unwrap_or(ts[n / 2]);
    let model = Orbit::with(eo, ts[0], forces);
    let xg: Vec<V3> = ts.iter().zip(&tr.xs).map(|(&t, &x)| mtv(&eo.c2t(t), x)).collect();
    let (ts_fit, xg_fit) = (&ts[..], &xg[..]);

    let mut started = "cold";
    let mut fit = None;
    if let Some(w) = warm {
        // The previous run's state, carried forward to this arc's start with its own forces.
        let y0 = if w.t0 < ts[0] { Orbit::with(eo, w.t0, forces).propagate_states(&w.state, &w.srp, &[ts[0]])[0] } else { w.state };
        let f = model.fit(ts_fit, xg_fit, Some((y0, w.srp)));
        started = if f.2 > WARM_MAX_RMS || !f.4 { "warm-rejected" } else { "warm" };
        fit = Some(f);
    }
    if started != "warm" { fit = Some(model.fit(ts_fit, xg_fit, None)); }
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
    // what a run made when this orbit appeared could have seen
    let seen: Vec<Nav> = recs.iter().filter(|r| r.toe_abs <= pred_start + LATENCY).cloned().collect();
    let recs = &seen[..];
    let off_rac: V3 = [0, 1, 2].map(|c| median(&mut d_rac.iter().map(|d| d[c]).collect()));
    let latest = recs.last()?;
    let ok: Vec<usize> = (0..n).filter(|&k| !tr.cs[k].is_nan()).collect();
    let broadcast_fit = |days: f64, quad: bool| -> Option<Box<dyn Fn(f64) -> (f64, f64)>> {
        let from = pred_start + LATENCY - days * 86400.0;
        let mut all: BTreeMap<i64, f64> = history.iter().map(|&(t, c)| (t.round() as i64, c)).collect();
        for r in recs { all.insert(r.toc_abs.round() as i64, r.af0); }
        let pts: Vec<(f64, f64)> = all.iter().map(|(&t, &c)| (t as f64, c))
            .filter(|&(t, _)| t >= from && t <= pred_start + LATENCY)
            .map(|(t, c)| (t - pred_start, c)).collect();
        let c = fit_clock_poly(&pts, quad)?;
        let c2 = if quad { c[2] } else { 0.0 };
        Some(Box::new(move |dt: f64| (c[0] + c[1] * dt + c2 * dt * dt, c[1] + 2.0 * c2 * dt)))
    };
    let from_broadcast = match clock_model {
        ClockModel::BroadcastLine(d) => broadcast_fit(d, false),
        ClockModel::BroadcastQuad(d) => broadcast_fit(d, true),
        _ => None,
    };
    let clock: Box<dyn Fn(f64) -> (f64, f64)> = if let Some(f) = from_broadcast { f }
    else if ok.len() >= 20 && d_clk.len() >= 10 {
        let off = median(&mut d_clk);
        let coef = fit_clock(&ok.iter().map(|&k| ts[k] - pred_start).collect::<Vec<_>>(),
                             &ok.iter().map(|&k| tr.cs[k]).collect::<Vec<_>>(), clock_model);
        Box::new(move |dt: f64| clock_at(&coef, dt, clock_model, off))
    } else {
        let (c0, c1) = (latest.clock(pred_start), latest.af1);   // no IGS clock: extend the latest broadcast's
        Box::new(move |dt: f64| (c0 + c1 * dt, c1))
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
        let (af0, af1) = clock(toe - pred_start);
        records.push((toe, Nav {
            p, toe_abs: toe, toc_abs: toe, af0, af1, af2: 0.0,
            iode: (iodc as i64 & 0xFF) as f64, iodc, svacc: URA_BY_DAY[age_days.min(URA_BY_DAY.len() - 1)],
            health: 0.0, tgd: latest.tgd, ttx: (toe - 7200.0).rem_euclid(WEEK),
        }));
        toe += STEP;
    }
    let _ = sv;
    // broadcast clock history to carry to the next run: the last seven days
    let mut all: BTreeMap<i64, f64> = history.iter().map(|&(t, c)| (t.round() as i64, c)).collect();
    for r in recs { all.insert(r.toc_abs.round() as i64, r.af0); }
    let clocks = all.iter().map(|(&t, &c)| (t as f64, c)).filter(|&(t, _)| t >= pred_start + LATENCY - 7.0 * 86400.0).collect();
    Some(SatResult {
        records, fit: FitState { t0: ts[0], state, srp, rms: fit_rms }, clocks, started, iterations, converged,
        seconds: start.elapsed().as_secs_f64(),
    })
}

/// A line (or, with `quad`, a quadratic) through clock samples (t, clock): least squares, then up
/// to three rounds dropping points more than four robust sigmas out (a clock jump or a bad upload)
/// and refitting. Coefficients c0, c1[, c2] in the samples' time origin; None with under 6 points.
pub fn fit_clock_poly(pts: &[(f64, f64)], quad: bool) -> Option<Vec<f64>> {
    if pts.len() < 6 { return None; }
    let m = if quad { 3 } else { 2 };
    let fit = |use_: &dyn Fn(usize) -> bool| -> Vec<f64> {
        let (mut nm, mut g) = (vec![0.0; m * m], vec![0.0; m]);
        for (k, &(t, y)) in pts.iter().enumerate() {
            if !use_(k) { continue; }
            let b = [1.0, t, t * t];
            for i in 0..m { g[i] += b[i] * y; for j in 0..m { nm[i * m + j] += b[i] * b[j]; } }
        }
        solve(m, &nm, &g)
    };
    let eval = |c: &[f64], t: f64| c[0] + c[1] * t + if m == 3 { c[2] * t * t } else { 0.0 };
    let mut c = fit(&|_| true);
    for _ in 0..3 {
        let mut r: Vec<f64> = pts.iter().map(|&(t, y)| (y - eval(&c, t)).abs()).collect();
        let lim = 4.0 * 1.4826 * median(&mut r.clone()).max(1e-10);
        if r.iter().all(|&x| x <= lim) { break; }
        let keep: Vec<bool> = r.iter_mut().map(|x| *x <= lim).collect();
        if keep.iter().filter(|&&k| k).count() < 6 { break; }
        c = fit(&|k| keep[k]);
    }
    if !quad { c.push(0.0); }
    Some(c)
}

/// Orbital frequency of GPS satellites: two revolutions per sidereal day.
const ORBIT_CPD: f64 = 2.0029;

fn clock_basis(dt: f64, model: ClockModel) -> Vec<f64> {
    let mut b = vec![1.0, dt];
    if matches!(model, ClockModel::Harmonic | ClockModel::HarmonicDrift) {
        for n in 1..=4 {
            let w = 2.0 * std::f64::consts::PI * ORBIT_CPD * n as f64 / 86400.0;
            b.push((w * dt).cos());
            b.push((w * dt).sin());
        }
    }
    b
}

/// Least-squares clock coefficients over times dt (s from prediction start) and clocks (s).
fn fit_clock(dt: &[f64], c: &[f64], model: ClockModel) -> Vec<f64> {
    let m = clock_basis(0.0, model).len();
    let (mut nm, mut g) = (vec![0.0; m * m], vec![0.0; m]);
    for (&t, &y) in dt.iter().zip(c) {
        let b = clock_basis(t, model);
        for i in 0..m {
            g[i] += b[i] * y;
            for j in 0..m { nm[i * m + j] += b[i] * b[j]; }
        }
    }
    solve(m, &nm, &g)
}

/// (af0, af1) at dt seconds from prediction start: the clock and its rate.
fn clock_at(coef: &[f64], dt: f64, model: ClockModel, offset: f64) -> (f64, f64) {
    let keep = if model == ClockModel::HarmonicDrift { 2 } else { coef.len() };
    let b = clock_basis(dt, model);
    let value: f64 = (0..keep).map(|i| coef[i] * b[i]).sum();
    let mut rate = coef[1];
    for n in 1..=((keep.saturating_sub(2)) / 2) {
        let w = 2.0 * std::f64::consts::PI * ORBIT_CPD * n as f64 / 86400.0;
        rate += w * (-coef[2 * n] * (w * dt).sin() + coef[2 * n + 1] * (w * dt).cos());
    }
    (value + offset, rate)
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

type ClockHistory = BTreeMap<u8, Vec<(f64, f64)>>;

/// A previous run's fit-state.txt: fitted states, and recent broadcast clock values ("clk" lines),
/// so a run downloads only the newest broadcasts.
fn read_warm(src: &str) -> (BTreeMap<u8, FitState>, ClockHistory) {
    let (mut fits, mut clocks) = (BTreeMap::new(), ClockHistory::new());
    let Ok(text) = read_text(src) else { eprintln!("warm start unavailable: {src}"); return (fits, clocks) };
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        if let Some(rest) = line.strip_prefix("clk ") {
            let v: Vec<f64> = rest.split_whitespace().filter_map(|x| x.parse().ok()).collect();
            if v.len() == 3 { clocks.entry(v[0] as u8).or_default().push((v[1], v[2])); }
            continue;
        }
        let v: Vec<f64> = line.split_whitespace().filter_map(|x| x.parse().ok()).collect();
        // sv t0 state(6) srp(5 or 9) rms
        if v.len() == 8 + 5 + 1 || v.len() == 8 + crate::orbit::NSRP + 1 {
            let mut srp = [0.0; crate::orbit::NSRP];
            for (k, x) in v[8..v.len() - 1].iter().enumerate() { srp[k] = *x; }
            fits.insert(v[0] as u8, FitState { t0: v[1], state: [v[2], v[3], v[4], v[5], v[6], v[7]],
                                               srp, rms: v[v.len() - 1] });
        }
    }
    (fits, clocks)
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
    let erp_text = read_text(&erp_src)?;
    let igs_erp = read_erp(&erp_text)?;
    let mut erp = extrapolate(&igs_erp, last_lod(&erp_text));
    let mut eop_used = String::from("IGS ERP, extrapolated with its last length of day");
    if o.eop != "none" {
        let joined = read_text(&o.eop).map(|t| join_bulletin_a(&igs_erp, &read_bulletin_a(&t))).unwrap_or_default();
        if joined.is_empty() {
            eprintln!("Bulletin A unavailable or not covering the week ({}); extrapolating the IGS ERP", o.eop);
        } else {
            eop_used = format!("IGS ERP, then IERS Bulletin A ({})", o.eop.rsplit('/').next().unwrap_or(&o.eop));
            erp = joined;
        }
    }
    eprintln!("Earth rotation: {eop_used}");
    if sp3.is_empty() { return Err("no GPS satellites in the orbit file".into()); }
    let t0 = sp3.values().map(|t| t.ts[0]).fold(f64::INFINITY, f64::min);
    let t_arc_end = sp3.values().map(|t| *t.ts.last().unwrap()).fold(f64::NEG_INFINITY, f64::max);

    let (warm, history) = o.warm.as_deref().map(read_warm).unwrap_or_default();
    // Broadcast history for the clock fit comes from the previous run when it has it; otherwise
    // download six days of it once.
    let history_from = history.values().flatten().map(|c| c.0).fold(f64::INFINITY, f64::min);
    let first_day = if history_from <= t0 - 4.0 * 86400.0 { t0 } else { t0 - 5.0 * 86400.0 };
    let mut navs = Navs::new();
    let mut used = vec![];
    if !o.rinex.is_empty() {
        for f in &o.rinex { read_nav(&read_text(f)?, &mut navs)?; used.push(f.clone()); }
    } else {
        let (mut day, last) = (((first_day + UTC_GPS + UNIX_GPS0) / 86400.0).floor() as i64,
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
                if let Some(r) = predict_sat(sv, tr, &navs[&sv], &eo, t_end, warm.get(&sv).copied(), o.forces, o.clock,
                                             history.get(&sv).map_or(&[][..], |v| &v[..])) {
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
    let mut state = String::from("# sv t0 x y z vx vy vz (GCRS m, m/s) D0 Y0 B0 Bc Bs D2c D2s D4c D4s (nm/s^2) rms (m)\n");
    for (sv, r) in &results {
        let f = &r.fit;
        state += &format!("{sv} {:.3} {} {} {:.6}\n", f.t0,
            f.state.iter().map(|v| format!("{v:.9}")).collect::<Vec<_>>().join(" "),
            f.srp.iter().map(|v| format!("{v:.9}")).collect::<Vec<_>>().join(" "), f.rms);
    }
    state += "# clk sv toc (GPS s) af0 (s): recent broadcast clocks, for the next run's clock fit\n";
    for (sv, r) in &results {
        for (t, c) in &r.clocks { state += &format!("clk {sv} {t:.0} {c:.15e}\n"); }
    }
    std::fs::write(path("fit-state.txt"), state).map_err(|e| e.to_string())?;
    let sats: Vec<String> = results.keys().map(|s| s.to_string()).collect();
    let json = format!(r#"{{
  "made": "{made}",
  "igs_orbit": "{sp3_src}",
  "igs_erp": "{erp_src}",
  "earth_rotation": "{eop_used}",
  "clock_model": "{clock_model:?}",
  "solid_tides": {tides},
  "relativity": {rel},
  "srp_prior_sigma_nm_s2": "{srp:?}",
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
  "backtest": "range error median/95%: 0.2/0.7 m first half day, 0.6/4.3 m at 2 days, 2.2/12.3 m at 4.5-5 days, worst 14.6 m (IGS orbit of 2026-10-13 06:00); second backtest from 2026-10-14 00:00: 1.3/7.5 m at 4-4.5 days"
}}
"#, made = iso_utc(now), bc = used.iter().map(|u| format!("\"{u}\"")).collect::<Vec<_>>().join(", "),
        ps = iso_utc(pred_start), vu = iso_utc(t_end), ipu = iso_utc(t_arc_end), sats = sats.join(", "),
        nrec = entries.len(), nb = blob.len(), fm = median(&mut fits), fw = worst_fit,
        cold = count("cold"), warm_ = count("warm"), rej = count("warm-rejected"),
        nc = results.values().filter(|r| !r.converged).count(), im = median(&mut iters),
        sm = median(&mut secs.clone()), st = secs.iter().sum::<f64>(), rs = t_run.elapsed().as_secs_f64(),
        cn = check.0, cm = check.1, cp = check.2, clock_model = o.clock, tides = o.forces.tides, rel = o.forces.relativity, srp = o.forces.srp_sigma);
    std::fs::write(path("gps-7d.json"), &json).map_err(|e| e.to_string())?;
    print!("{json}");
    if check.0 == 0 || !(check.1 <= o.max_check_m) {
        return Err(format!("self-check failed: median {:.3} m over {} samples", check.1, check.0));
    }
    let _ = &mut secs;
    Ok(())
}
