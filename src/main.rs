// SPDX-License-Identifier: MIT
//! extended-ephemeris: offline rapid cold start for GPS. Predicts the satellites' orbit data a
//! week ahead from IGS ultra-rapid orbits and writes it as RINEX 3 and CASIC frames.
mod casic;
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
mod sp3;
mod time;
#[cfg(test)]
mod tests;

const USAGE: &str = "\
extended-ephemeris: a week of predicted GPS orbit data (extended ephemeris)

usage:
  extended-ephemeris [predict] [--out DIR] [--days N] [--warm fit-state.txt]
                     [--sp3 FILE|URL --erp FILE|URL] [--rinex FILE|URL ...]
                     [--max-check-m M] [--jobs N]
  extended-ephemeris grade PREDICTED.rnx BROADCAST.rnx[.gz] ...

predict (the default) fetches the newest IGS ultra-rapid orbit and the broadcasts it needs from
BKG, and writes gps-7d.rnx, gps-7d-casic.bin, gps-7d.json and fit-state.txt to --out (default .).
grade compares a predicted RINEX file with what the satellites broadcast.
Needs curl and gzip on the path.";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("-h") | Some("--help") => { println!("{USAGE}"); return; }
        Some("grade") if args.len() >= 3 => grade::run(&args[1], &args[2..]),
        Some("grade") => Err(USAGE.into()),
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
            "--rinex" => {
                i += 1;
                while i < args.len() && !args[i].starts_with("--") { o.rinex.push(args[i].clone()); i += 1; }
            }
            other => return Err(format!("unknown argument {other}\n\n{USAGE}")),
        }
    }
    Ok(o)
}
