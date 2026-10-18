// SPDX-License-Identifier: MIT
//! select: from a predicted week (RINEX, a local file or the published URL), pick each satellite's
//! record nearest a moment (now by default, within 2 h) and write a ready-to-send assistance burst
//! for a receiver: for u-blox, MGA-INI-TIME_UTC [+ MGA-INI-POS_LLH] then MGA-GPS-EPH; for CASIC
//! (Quectel L76K and other AT6558 receivers), AID-INI then MSG-GPSEPH. Galileo records, when the
//! file has them, go out as UBX-MGA-GAL-EPH (u-blox only).
use crate::fetch::read_text;
use crate::rinex::{read_nav, Navs};
use crate::time::{gps_seconds, iso_utc, now_gps, UTC_GPS};
use crate::{casic, ubx};

pub const PUBLISHED: &str = "https://raw.githubusercontent.com/badriram/open-extended-ephemeris/ephemeris-data/gps-7d.rnx";

pub fn run(args: &[String]) -> Result<(), String> {
    let (mut format, mut src, mut out, mut at) = (String::new(), PUBLISHED.to_string(), String::new(), None::<f64>);
    let mut pos: Option<(f64, f64, f64)> = None;
    let mut pos_acc = 10000.0f64;
    let mut time_acc = 2.0f64;
    let mut i = 0;
    let val = |i: usize| args.get(i + 1).cloned().ok_or_else(|| format!("{} needs a value", args[i]));
    while i < args.len() {
        match args[i].as_str() {
            "--format" => { format = val(i)?; i += 2; }
            "--out" | "-o" => { out = val(i)?; i += 2; }
            "--at" => { at = Some(parse_iso(&val(i)?)?); i += 2; }
            "--position" => {
                let v: Vec<f64> = val(i)?.split(',').filter_map(|x| x.trim().parse().ok()).collect();
                if v.len() < 2 { return Err("--position LAT,LON[,HEIGHT_M] in degrees".into()); }
                pos = Some((v[0], v[1], *v.get(2).unwrap_or(&0.0)));
                i += 2;
            }
            "--position-accuracy" => { pos_acc = val(i)?.parse().map_err(|_| "metres")?; i += 2; }
            "--time-accuracy" => { time_acc = val(i)?.parse().map_err(|_| "seconds")?; i += 2; }
            s if !s.starts_with("--") => { src = s.to_string(); i += 1; }
            other => return Err(format!("unknown argument {other}")),
        }
    }
    if format != "ubx" && format != "casic" { return Err("--format ubx or --format casic".into()); }
    if out.is_empty() { return Err("--out FILE".into()); }
    let t = at.unwrap_or_else(now_gps);
    let mut navs = Navs::new();
    read_nav(&read_text(&src)?, &mut navs)?;
    let mut burst = match (format.as_str(), pos) {
        ("ubx", p) => {
            let mut b = ubx::ini_time_utc(t, time_acc);
            if let Some((la, lo, h)) = p { b.extend(ubx::ini_pos_llh(la, lo, h, pos_acc)); }
            b
        }
        (_, p) => casic::aid_ini(t, time_acc as f32, p.map(|(la, lo, h)| (la, lo, h, pos_acc as f32))),
    };
    let mut n = 0;
    for (sv, recs) in &navs {
        let Some(r) = recs.iter().filter(|r| (r.toe_abs - t).abs() <= 7200.0)
            .min_by(|a, b| (a.toe_abs - t).abs().partial_cmp(&(b.toe_abs - t).abs()).unwrap()) else { continue };
        let prn = crate::rinex::sat_prn(*sv);
        match (format.as_str(), r.sys) {
            ("ubx", 'E') => burst.extend(ubx::galeph_frame(r, prn)?),
            ("ubx", _) => burst.extend(ubx::gpseph_frame(r, prn)?),
            (_, 'G') => burst.extend(casic::gpseph_frame(r, prn)?),
            _ => continue,                       // AT6558 receivers have no Galileo
        }
        n += 1;
    }
    if n == 0 {
        return Err(format!("no record within 2 h of {} in {src} (is it more than 7 days old?)", iso_utc(t)));
    }
    std::fs::write(&out, &burst).map_err(|e| format!("{out}: {e}"))?;
    eprintln!("{out}: {} bytes for {}: {} time{} + {n} ephemerides", burst.len(), iso_utc(t),
              if format == "ubx" { "MGA-INI" } else { "AID-INI" }, if pos.is_some() { "/position" } else { "" });
    Ok(())
}

/// "2026-10-18T12:00:00Z" (UTC) -> GPS seconds.
pub fn parse_iso(s: &str) -> Result<f64, String> {
    let f: Vec<f64> = s.trim_end_matches('Z').split(|c| c == '-' || c == 'T' || c == ':' || c == ' ')
        .filter_map(|x| x.parse().ok()).collect();
    if f.len() < 5 { return Err(format!("--at {s}: use YYYY-MM-DDTHH:MM[:SS]Z (UTC)")); }
    Ok(gps_seconds(f[0] as i64, f[1] as i64, f[2] as i64, f[3] as i64, f[4] as i64, *f.get(5).unwrap_or(&0.0)) - UTC_GPS)
}
