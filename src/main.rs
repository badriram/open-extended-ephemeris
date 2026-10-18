// SPDX-License-Identifier: MIT
//! extended-ephemeris: offline rapid cold start for GPS. Predicts the satellites' orbit data a
//! week ahead from IGS ultra-rapid orbits and writes it as RINEX 3 and CASIC frames.
mod boxwing;
mod casic;
mod clockbt;
mod dop853;
mod eo;
mod erfa;
mod fetch;
mod grade;
mod gravity;
mod linalg;
mod lnav;
mod orbit;
mod predict;
mod rinex;
mod select;
mod sp3;
mod time;
mod ubx;
#[cfg(test)]
mod tests;

const USAGE: &str = "\
extended-ephemeris: a week of predicted GPS orbit data (extended ephemeris)

usage:
  extended-ephemeris [predict] [--out DIR] [--days N] [--warm fit-state.txt]
                     [--sp3 FILE|URL --erp FILE|URL] [--rinex FILE|URL ...]
                     [--max-check-m M] [--jobs N] [--clock line|harmonic|harmonic-drift]
                     [--no-tides] [--no-relativity] [--srp-sigma S1,...,S9]
                     [--eop BULLETIN-A-FILE|URL|none] [--systems G,E]
  extended-ephemeris grade PREDICTED.rnx BROADCAST.rnx[.gz] ...
  extended-ephemeris select --format ubx|casic --out FILE [--at YYYY-MM-DDTHH:MMZ]
                     [--position LAT,LON[,HEIGHT]] [PREDICTED.rnx|URL]
  extended-ephemeris clock-backtest [--models quad:5,line:1,...] [--every HOURS] BROADCAST.rnx[.gz] ...

predict (the default) fetches the newest IGS ultra-rapid orbit and the broadcasts it needs from
BKG, and writes gps-7d.rnx, gps-7d-casic.bin, gps-7d-ubx.bin, gps-7d.json and fit-state.txt to
--out (default .).
select writes a ready-to-send burst (time, optional position, the ephemerides for now) for a
u-blox or CASIC receiver, from a predicted RINEX file (default: the published one).
grade compares a predicted RINEX file with what the satellites broadcast.
clock-backtest tests clock models against a span of broadcast history (a year of daily files).
Needs curl and gzip on the path.";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("-h") | Some("--help") => { println!("{USAGE}"); return; }
        Some("grade") if args.len() >= 3 => grade::run(&args[1], &args[2..]),
        Some("grade") => Err(USAGE.into()),
        Some("clock-backtest") => clockbt::run(&args[1..]),
        Some("select") => select::run(&args[1..]),
        _ => parse(if args.first().map(String::as_str) == Some("predict") { &args[1..] } else { &args })
            .and_then(|o| predict::run(&o)),
    };
    if let Err(e) = result {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn parse(args: &[String]) -> Result<predict::Options, String> {
    let mut o = predict::Options {
        out: ".".into(), sp3: None, erp: None, rinex: vec![], days: 7.0, warm: None, max_check_m: 5.0,
        jobs: std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4),
        forces: orbit::Forces::default(), clock: predict::ClockModel::BroadcastQuad(5.0),
        eop: predict::BULLETIN_A.into(), from_broadcast: None, arc_hours: 48.0, systems: vec!['G'], prediction_start: None,
    };
    let mut i = 0;
    let val = |i: usize| args.get(i + 1).cloned().ok_or_else(|| format!("{} needs a value\n\n{USAGE}", args[i]));
    while i < args.len() {
        match args[i].as_str() {
            "--out" => { o.out = val(i)?; i += 2; }
            "--sp3" => { o.sp3 = Some(val(i)?); i += 2; }
            "--erp" => { o.erp = Some(val(i)?); i += 2; }
            "--warm" => { o.warm = Some(val(i)?); i += 2; }
            "--days" => { o.days = val(i)?.parse().map_err(|_| "--days: a number")?; i += 2; }
            "--max-check-m" => { o.max_check_m = val(i)?.parse().map_err(|_| "--max-check-m: a number")?; i += 2; }
            "--jobs" => { o.jobs = val(i)?.parse().map_err(|_| "--jobs: a number")?; i += 2; }
            "--no-tides" => { o.forces.tides = false; i += 1; }
            "--eop" => { o.eop = val(i)?; i += 2; }
            "--no-relativity" => { o.forces.relativity = false; i += 1; }
            "--box-wing" => { o.forces.box_wing = true; i += 1; }
            "--no-antenna-offset" => { o.forces.antenna_offset = false; i += 1; }
            "--systems" => {
                o.systems = val(i)?.split(',').filter_map(|x| x.trim().chars().next()).collect();
                if o.systems.iter().any(|c| *c != 'G' && *c != 'E') { return Err("--systems: G and/or E".into()); }
                i += 2;
            }
            "--arc-hours" => { o.arc_hours = val(i)?.parse().map_err(|_| "--arc-hours: a number")?; i += 2; }
            "--prediction-start" => { o.prediction_start = Some(select::parse_iso(&val(i)?)?); i += 2; }
            "--from-broadcast" => { o.from_broadcast = Some(select::parse_iso(&val(i)?)?); i += 2; }
            "--srp-sigma" => {
                // nine comma-separated a priori sigmas, nm/s²: D0 Y0 B0 Bc Bs D2c D2s D4c D4s
                let v: Vec<f64> = val(i)?.split(',').filter_map(|x| x.parse().ok()).collect();
                if v.len() != orbit::NSRP { return Err("--srp-sigma: nine comma-separated numbers".into()); }
                o.forces.srp_sigma.copy_from_slice(&v);
                i += 2;
            }
            "--clock" => {
                o.clock = match val(i)?.as_str() {
                    "line" => predict::ClockModel::Line,
                    "harmonic" => predict::ClockModel::Harmonic,
                    "harmonic-drift" => predict::ClockModel::HarmonicDrift,
                    v if v.starts_with("broadcast-line:") => predict::ClockModel::BroadcastLine(v[15..].parse().map_err(|_| "days")?),
                    v if v.starts_with("broadcast-quad:") => predict::ClockModel::BroadcastQuad(v[15..].parse().map_err(|_| "days")?),
                    v => return Err(format!("--clock {v}: line, harmonic, harmonic-drift, broadcast-line:DAYS or broadcast-quad:DAYS")),
                };
                i += 2;
            }
            "--rinex" => {
                i += 1;
                while i < args.len() && !args[i].starts_with("--") { o.rinex.push(args[i].clone()); i += 1; }
            }
            other => return Err(format!("unknown argument {other}\n\n{USAGE}")),
        }
    }
    Ok(o)
}
