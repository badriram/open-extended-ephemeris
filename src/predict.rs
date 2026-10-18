// SPDX-License-Identifier: MIT
//! The prediction, per GPS or Galileo satellite:
//!
//! 1. Take the newest orbit prediction and its Earth rotation parameters: CODE's 5-day GPS and
//!    Galileo prediction by default (joined to the one before it while the newest starts in the
//!    future), or the IGS ultra-rapid orbit (24 h measured, 24 h predicted), GPS only.
//! 2. Fit the force model (orbit.rs) to that arc and integrate it a week ahead.
//! 3. Every 2 h, fit a broadcast ephemeris (lnav.rs) to 4 h of trajectory: to the source's own
//!    prediction while it lasts, then to ours.
//! 4. Calibrate against the satellite's real broadcasts over the known part of the arc: the
//!    antenna offset (precise orbits are centre of mass).
//! 5. Predict the clock from five days of the satellite's own broadcast clocks, as a quadratic.
//! 6. Give each record an accuracy (URA/SISA): the population figure for its age, or more for a
//!    satellite whose clock predicted its own last 1.5 days badly.
use crate::eo::{extrapolate, join_bulletin_a, last_lod, read_bulletin_a, read_erp, EarthOrientation};
use crate::fetch::{brdc_urls, latest_code_prd, latest_igu, read_text};
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
/// Galileo satellites left out: the In-Orbit Validation satellites, whose clocks were the least
/// predictable over 2025 (day 4-5 medians 11-35 m, 95th percentiles 100-570 m, against ~1 m for the
/// rest; results/2026-10-18-code-prediction.md).
const EXCLUDED: [u8; 3] = [111, 112, 119];   // E11, E12, E19
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
    /// Build the 48 h tracks from the satellites' broadcast orbits known at this GPS time, which
    /// is also where the prediction starts, instead of from an IGS ultra-rapid orbit.
    pub from_broadcast: Option<f64>,
    /// Length of the broadcast-built arc, hours.
    pub arc_hours: f64,
    /// Treat orbit-file samples after this GPS time as predicted, whatever their flags (for CODE's
    /// 5-day prediction, whose flags start minutes after its first epoch).
    pub prediction_start: Option<f64>,
    /// Where orbits come from when no --sp3 is given: "auto" (CODE's 5-day prediction, falling
    /// back to the IGS ultra-rapid orbit for GPS alone), "code" or "igs".
    pub source: String,
    /// Constellations to predict: 'G' (from the IGS orbit, or broadcasts with --from-broadcast) and
    /// 'E' (Galileo, always from broadcasts: there is no anonymous multi-GNSS ultra-rapid orbit).
    pub systems: Vec<char>,
    /// --systems was given: build broadcast arcs for systems the orbit file lacks (otherwise only
    /// with --from-broadcast; Galileo from broadcasts is too weak to publish by default).
    pub systems_explicit: bool,
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
    let pred_start = tr.predicted.iter().position(|&p| p).map(|k| ts[k]).unwrap_or(ts[n - 1]);
    let galileo_bw = || if forces.box_wing && crate::rinex::sat_sys(sv) == 'E' {
        Some(crate::boxwing::BoxWing::for_galileo(crate::rinex::sat_prn(sv)))
    } else { None };
    let mut model = Orbit::with(eo, ts[0], forces);
    model.box_wing = galileo_bw();
    model.fit_antenna_offset = tr.from_broadcast && forces.antenna_offset;
    let xg: Vec<V3> = ts.iter().zip(&tr.xs).map(|(&t, &x)| mtv(&eo.c2t(t), x)).collect();
    let (ts_fit, xg_fit) = (&ts[..], &xg[..]);

    let mut started = "cold";
    let mut fit = None;
    if let Some(w) = warm {
        // The previous run's state, carried forward to this arc's start with its own forces.
        let y0 = if w.t0 < ts[0] {
            let mut prev = Orbit::with(eo, w.t0, forces);
            prev.box_wing = galileo_bw();
            prev.propagate_states(&w.state, &w.srp, &[ts[0]])[0]
        } else { w.state };
        let f = model.fit(ts_fit, xg_fit, Some((y0, w.srp)));
        started = if f.2 > WARM_MAX_RMS || !f.4 { "warm-rejected" } else { "warm" };
        fit = Some(f);
    }
    if started != "warm" { fit = Some(model.fit(ts_fit, xg_fit, None)); }
    let (state, srp, fit_rms, iterations, converged) = fit.unwrap();

    // How predictable this satellite's clock is: fit the clock model to the 3.5 days of broadcast
    // history before the last 1.5 days, predict those, and take the median error. Each record's
    // accuracy (URA/SISA) is then the larger of the population figure for its age and this
    // satellite's own error grown with age (as the year of 2025 clocks grew, ~age^1.3, doubled
    // for a 95th percentile). Receivers can then weight a poor clock down.
    let hold_out_m: Option<f64> = {
        let end = pred_start + LATENCY;
        let mut all: BTreeMap<i64, f64> = history.iter().map(|&(t, c)| (t.round() as i64, c)).collect();
        for r in recs { all.insert(r.toc_abs.round() as i64, r.af0); }
        let pts: Vec<(f64, f64)> = all.iter().map(|(&t, &c)| (t as f64 - end, c)).filter(|&(t, _)| t >= -5.0 * 86400.0 && t <= 0.0).collect();
        let (fit_pts, test): (Vec<_>, Vec<_>) = pts.iter().partition(|&&(t, _)| t < -1.5 * 86400.0);
        fit_clock_poly(&fit_pts, true).filter(|_| test.len() >= 6).map(|c| {
            let mut e: Vec<f64> = test.iter().map(|&(t, y)| (C_LIGHT * (c[0] + c[1] * t + c[2] * t * t - y)).abs()).collect();
            median(&mut e)
        })
    };
    let accuracy = |age_days: f64| -> f64 {
        let table = URA_BY_DAY[(age_days.max(0.0) as usize).min(URA_BY_DAY.len() - 1)];
        let own = hold_out_m.map_or(0.0, |e| 2.0 * e * (age_days.max(0.25) / 0.75).powf(1.3));
        table.max(own)
    };
    let grid: Vec<f64> = (0..).map(|i| ts[0] + 300.0 * i as f64).take_while(|&t| t <= t_end + STEP + 1.0).collect();
    let xe: Vec<V3> = model.propagate(&state, &srp, &grid).iter().zip(&grid).map(|(x, &t)| mv(&eo.c2t(t), *x)).collect();
    // Calibration over the measured day: where the broadcasts put the satellite and its clock,
    // relative to the IGS centre of mass and clock.
    let (mut d_rac, mut d_clk): (Vec<V3>, Vec<f64>) = (vec![], vec![]);
    // Normally the orbit's measured part; an orbit that is all prediction (CODE's 5-day file from
    // its first epoch) falls back to whatever the broadcasts seen by the run cover.
    let measured = tr.predicted.iter().filter(|&&p| !p).count() >= 40;
    for k in (1..n - 1).step_by(4) {
        if if measured { tr.predicted[k] } else { ts[k] > pred_start + LATENCY } { continue; }
        let Some(b) = nearest(recs, ts[k]) else { continue };
        let v = scale(sub(tr.xs[k + 1], tr.xs[k - 1]), 1.0 / (ts[k + 1] - ts[k - 1]));
        // Against the IGS orbit (centre of mass) when there is one. Built from broadcasts, the
        // track is itself at the antenna, which no orbit can follow; the fitted orbit settles on
        // the centre of mass, so measure the offset against that instead (it was +1.4 m radial).
        let reference = if tr.from_broadcast { xe[((ts[k] - grid[0]) / 300.0).round() as usize] } else { tr.xs[k] };
        d_rac.push(mv(&rac(reference, v), sub(b.pos(ts[k]), reference)));
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


    let mut records = Vec::new();
    // The first record whose whole fit window (toe ± 2 h) the orbit covers. With CODE's 5-day
    // file, which starts at the prediction start, an earlier toe extrapolated the interpolation
    // off the start of the data (errors of 40-150 km in the self-check).
    let mut toe = ((pred_start / STEP).ceil() * STEP).max((ts[0] / STEP).ceil() * STEP + STEP);
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
        let gm = if crate::rinex::sat_sys(sv) == 'E' { lnav::GM_GAL } else { lnav::GM };
        let (p, _rms) = lnav::fit(&win, &pos, toe, gm);
        let age_days = ((toe - pred_start) / 86400.0).max(0.0);
        let iodc = 768.0 + ((toe / STEP) as i64 % 240) as f64;
        let (af0, af1) = clock(toe - pred_start);
        let sys = crate::rinex::sat_sys(sv);
        records.push((toe, Nav {
            sys, ds: if sys == 'E' { 517.0 } else { 0.0 }, p, toe_abs: toe, toc_abs: toe, af0, af1, af2: 0.0,
            iode: (iodc as i64 & 0xFF) as f64, iodc, svacc: accuracy(age_days),
            health: 0.0, tgd: latest.tgd, ttx: (toe - 7200.0).rem_euclid(WEEK),
        }));
        toe += STEP;
    }
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
                let e = sisre(rec, b, t).0;
                if e > 100.0 && std::env::var("OEE_DEBUG").is_ok() {
                    eprintln!("self-check {} at {}: {e:.0} m against the record of toe {} (health {})",
                              crate::rinex::sat_name(*sv), iso_utc(t), iso_utc(b.toe_abs), b.health);
                }
                errs.push(e);
            }
        }
    }
    errs
}

/// Signal-in-space range error and 3D orbit difference of a record against a reference record.
pub fn sisre(r: &Nav, b: &Nav, t: f64) -> (f64, f64) {
    let xb = b.pos(t);
    let d = sub(r.pos(t), xb);
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
    let mut systems = o.systems.clone();
    let mut prediction_start = o.prediction_start;
    let (sp3_src, erp_src) = match (&o.sp3, &o.erp) {
        (Some(s), Some(e)) => (s.clone(), e.clone()),
        (None, None) if o.from_broadcast.is_some() => return Err("--from-broadcast needs --erp (and --rinex)".into()),
        (None, None) => {
            let code = if o.source == "igs" { Err("not asked for".to_string()) } else { latest_code_prd(now) };
            match code {
                Ok(pair) => {
                    // CODE's file starts at midnight of its issue day: predict from now, and use the
                    // hours before as the orbit's known part
                    prediction_start.get_or_insert(now);
                    pair
                }
                Err(e) if o.source == "code" => return Err(format!("CODE: {e}")),
                Err(e) => {
                    if o.source != "igs" { eprintln!("CODE unavailable ({e}); falling back to the IGS ultra-rapid orbit, GPS only"); }
                    systems.retain(|&c| c == 'G');
                    latest_igu(now)?
                }
            }
        }
        (None, Some(e)) if o.from_broadcast.is_some() => (String::new(), e.clone()),
        _ => return Err("--sp3 and --erp go together".into()),
    };
    eprintln!("orbit: {}\nERP:   {erp_src}", if o.from_broadcast.is_some() { "broadcast orbits" } else { &sp3_src });
    // Several orbit files may be given, comma-separated and oldest first: each later one takes over
    // from its own first epoch.
    let mut sp3: BTreeMap<u8, Track> = BTreeMap::new();
    if o.from_broadcast.is_none() {
        for src in sp3_src.split(',') {
            for (sv, tr) in read_sp3(&read_text(src)?) {
                match sp3.get_mut(&sv) {
                    None => { sp3.insert(sv, tr); }
                    Some(old) => {
                        let keep = old.ts.iter().take_while(|&&t| t < tr.ts[0] - 1.0).count();
                        old.ts.truncate(keep); old.xs.truncate(keep); old.cs.truncate(keep); old.predicted.truncate(keep);
                        old.ts.extend(tr.ts); old.xs.extend(tr.xs); old.cs.extend(tr.cs); old.predicted.extend(tr.predicted);
                    }
                }
            }
        }
    }
    if let Some(ps) = prediction_start {
        for tr in sp3.values_mut() {
            tr.predicted = tr.ts.iter().map(|&t| t > ps + 1.0).collect();
        }
    }
    let igs_start = sp3.values().filter_map(|t| t.predicted.iter().position(|&p| p).map(|k| t.ts[k]))
        .fold(f64::INFINITY, f64::min);
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
    let t0 = if let Some(ps) = o.from_broadcast { ps - o.arc_hours * 3600.0 } else {
        if sp3.is_empty() { return Err("no GPS satellites in the orbit file".into()); }
        sp3.values().map(|t| t.ts[0]).fold(f64::INFINITY, f64::min)
    };

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
    sp3.retain(|sv, _| systems.contains(&crate::rinex::sat_sys(*sv)) && !EXCLUDED.contains(sv));
    // Tracks from broadcast orbits, for every system not read from an IGS orbit: all of them with
    // --from-broadcast, Galileo always. The arc ends where the prediction starts.
    let ps_bc = o.from_broadcast.unwrap_or(igs_start);
    for (sv, recs) in &navs {
        let sys = crate::rinex::sat_sys(*sv);
        if !systems.contains(&sys) || (sys == 'G' && o.from_broadcast.is_none()) || sp3.contains_key(sv)
            || EXCLUDED.contains(sv) || (o.from_broadcast.is_none() && !o.systems_explicit) { continue; }
        // a record is uploaded and broadcast from about two hours before its reference time
        let known: Vec<&Nav> = recs.iter().filter(|r| r.toe_abs <= ps_bc + 7200.0 && r.health == 0.0).collect();
        let mut tr = Track { from_broadcast: true, ..Track::default() };
        let nk = (o.arc_hours * 4.0).round() as usize;
        for k in 0..=nk {
            let t = ps_bc - o.arc_hours * 3600.0 + 900.0 * k as f64;
            let Some(r) = known.iter().filter(|r| (r.toe_abs - t).abs() <= 7200.0)
                .min_by(|a, b| (a.toe_abs - t).abs().partial_cmp(&(b.toe_abs - t).abs()).unwrap()) else { continue };
            tr.ts.push(t);
            tr.xs.push(r.pos(t));
            tr.cs.push(r.clock(t));
            tr.predicted.push(false);
        }
        if tr.ts.len() as f64 >= 0.75 * nk as f64 { sp3.insert(*sv, tr); }
    }
    if sp3.is_empty() { return Err("no satellites to predict".into()); }
    let t_arc_end = sp3.values().map(|t| *t.ts.last().unwrap()).fold(f64::NEG_INFINITY, f64::max);

    let pred_start = ps_bc;
    let t_end = pred_start + o.days * 86400.0;
    // The Earth-orientation table must cover every track: broadcast-built arcs start 48 h before
    // the prediction, IGS ones 24 h. (Sized from the IGS arc alone, it left the first 21 h of each
    // Galileo arc on a clamped, wrong Earth orientation: fits of ~27 m RMS.)
    let t0 = sp3.values().map(|t| t.ts[0]).fold(t0, f64::min);
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
    // A record that will not encode (a clock rate beyond the broadcast range, say, after a bad
    // fit) is dropped and reported; one odd satellite must never stop the whole publication.
    let before = entries.len();
    let mut dropped: Vec<String> = vec![];
    entries.retain(|(sv, n)| {
        let ok = if n.sys == 'E' { crate::ubx::galeph_payload(n, crate::rinex::sat_prn(*sv)).map(|_| ()) }
                 else { casic::gpseph_payload(n, *sv).and(crate::ubx::gpseph_payload(n, *sv)).map(|_| ()) };
        if let Err(e) = ok { dropped.push(format!("{} {}: {e}", crate::rinex::sat_name(*sv), iso_utc(n.toe_abs))); false } else { true }
    });
    if !dropped.is_empty() {
        eprintln!("dropped {} of {before} records that would not encode:\n  {}", dropped.len(), dropped.join("\n  "));
    }
    let mut blob = Vec::new();
    for (sv, n) in entries.iter().filter(|e| e.1.sys == 'G') { blob.extend(casic::gpseph_frame(n, *sv)?); }
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
    let mut ublob = Vec::new();
    for (sv, n) in &entries {
        ublob.extend(if n.sys == 'E' { crate::ubx::galeph_frame(n, crate::rinex::sat_prn(*sv))? } else { crate::ubx::gpseph_frame(n, *sv)? });
    }
    crate::ubx::check_frames(&ublob)?;
    std::fs::write(path("gps-7d-ubx.bin"), &ublob).map_err(|e| e.to_string())?;
    std::fs::write(path("gps-7d.rnx"), rinex::write_nav(&entries, &[
        "PREDICTED broadcast ephemeris (extended ephemeris), not".into(),
        "broadcast by the satellites. Orbits from CODE or IGS:".into(),
        format!("orbit {}", base(&sp3_src)),
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
    let sats: Vec<String> = results.keys().map(|s| format!("\"{}\"", crate::rinex::sat_name(*s))).collect();
    let backtest = if sp3_src.contains("COD0OPSPRD") {
        "range error median/95% at 4.5-5 days: GPS 1.8/12.2 m, Galileo 1.5/6.9 m (CODE 5-day orbit, predicting from 2026-10-13 06:00); from 2026-10-14 00:00, at 4-4.5 days: GPS 1.3/8.3 m, Galileo 0.7/5.8 m"
    } else {
        "range error median/95%: 0.2/0.7 m first half day, 0.6/4.3 m at 2 days, 2.2/12.3 m at 4.5-5 days, worst 14.6 m (IGS orbit of 2026-10-13 06:00); second backtest from 2026-10-14 00:00: 1.3/7.5 m at 4-4.5 days"
    };
    let json = format!(r#"{{
  "made": "{made}",
  "orbit": "{sp3_src}",
  "erp": "{erp_src}",
  "earth_rotation": "{eop_used}",
  "clock_model": "{clock_model:?}",
  "solid_tides": {tides},
  "relativity": {rel},
  "srp_prior_sigma_nm_s2": "{srp:?}",
  "broadcasts": [{bc}],
  "prediction_starts": "{ps}",
  "valid_until": "{vu}",
  "orbit_prediction_until": "{ipu}",
  "satellites": [{sats}],
  "records": {nrec},
  "dropped_records": {ndrop},
  "casic_bytes": {nb},
  "ubx_bytes": {nub},
  "force_model_fit_m": {{"median": {fm:.4}, "worst": {fw:.4}}},
  "fits": {{"cold": {cold}, "warm": {warm_}, "warm_rejected": {rej}, "not_converged": {nc}, "iterations_median": {im}}},
  "satellite_seconds": {{"median": {sm:.2}, "total": {st:.1}}},
  "run_seconds": {rs:.1},
  "self_check_vs_broadcast": {{"samples": {cn}, "median_m": {cm:.4}, "p95_m": {cp:.4}}},
  "backtest": "{backtest}"
}}
"#, made = iso_utc(now), bc = used.iter().map(|u| format!("\"{u}\"")).collect::<Vec<_>>().join(", "),
        ps = iso_utc(pred_start), vu = iso_utc(t_end), ipu = iso_utc(t_arc_end), sats = sats.join(", "),
        nrec = entries.len(), ndrop = dropped.len(), nb = blob.len(), nub = ublob.len(), fm = median(&mut fits), fw = worst_fit,
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
