// SPDX-License-Identifier: MIT
//! RINEX 3 navigation files, GPS records only: read the satellites' broadcasts, write predictions.
//!
//! A GPS record is an epoch line (PRN, the clock reference time toc, af0, af1, af2) and seven lines
//! of four values: IODE Crs deltaN M0 / Cuc e Cus sqrtA / toe Cic Omega0 Cis / i0 Crc omega
//! OmegaDot / IDOT L2codes week L2P / accuracy health TGD IODC / transmission-time fit-interval.
use crate::lnav::Params;
use crate::time::{calendar, gps_seconds, WEEK};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct Nav {
    pub p: Params,       // sqrtA e M0 omega i0 Omega0 OmegaDot deltaN IDOT Cuc Cus Crc Crs Cic Cis
    pub toe_abs: f64,    // continuous GPS seconds
    pub toc_abs: f64,
    pub af0: f64,
    pub af1: f64,
    pub af2: f64,
    pub iode: f64,
    pub iodc: f64,
    pub svacc: f64,      // metres
    pub health: f64,
    pub tgd: f64,
    pub ttx: f64,        // seconds of week
}

impl Nav {
    pub fn week(&self) -> f64 { (self.toe_abs / WEEK).floor() }
    pub fn toe(&self) -> f64 { self.toe_abs.rem_euclid(WEEK) }
    pub fn toc(&self) -> f64 { self.toc_abs.rem_euclid(WEEK) }
    pub fn clock(&self, t: f64) -> f64 {
        let d = t - self.toc_abs;
        self.af0 + self.af1 * d + self.af2 * d * d
    }
}

pub type Navs = BTreeMap<u8, Vec<Nav>>;

fn num(s: &str) -> f64 {
    let s = s.trim().replace(['D', 'd'], "E");
    if s.is_empty() { 0.0 } else { s.parse().unwrap_or(0.0) }
}

fn field(line: &str, a: usize, b: usize) -> &str {
    let n = line.len();
    if a >= n { "" } else { &line[a..b.min(n)] }
}

/// Add the GPS records of RINEX 3 navigation text to `out`.
pub fn read_nav(text: &str, out: &mut Navs) -> Result<(), String> {
    let lines: Vec<&str> = text.lines().collect();
    let mut k = lines.iter().position(|l| l.contains("END OF HEADER")).ok_or("no RINEX header")? + 1;
    while k < lines.len() {
        let l = lines[k];
        if !l.starts_with('G') || k + 7 >= lines.len() || l.len() < 23 { k += 1; continue; }
        let mut v: Vec<f64> = (0..3).map(|j| num(field(l, 23 + 19 * j, 42 + 19 * j))).collect();
        for m in 1..8 {
            v.extend((0..4).map(|j| num(field(lines[k + m], 4 + 19 * j, 23 + 19 * j))));
        }
        let e: Vec<f64> = l[4..23].split_whitespace().map(|x| x.parse().unwrap_or(0.0)).collect();
        if e.len() < 6 { k += 1; continue; }
        let toc_abs = gps_seconds(e[0] as i64, e[1] as i64, e[2] as i64, e[3] as i64, e[4] as i64, e[5]);
        let sv: u8 = l[1..3].trim().parse().map_err(|_| "bad PRN")?;
        // af0 af1 af2 | IODE Crs dn M0 | Cuc e Cus sqrtA | toe Cic W0 Cis | i0 Crc w Wd |
        // IDOT L2 week L2P | acc health TGD IODC | ttx fit
        let p = [v[10], v[8], v[6], v[17], v[15], v[13], v[18], v[5], v[19], v[7], v[9], v[16], v[4], v[12], v[14]];
        let nav = Nav {
            p, toe_abs: v[21] * WEEK + v[11], toc_abs, af0: v[0], af1: v[1], af2: v[2], iode: v[3],
            iodc: v[26], svacc: v[23], health: v[24], tgd: v[25], ttx: v[27],
        };
        let list = out.entry(sv).or_default();
        if !list.iter().any(|x| x.toe_abs == nav.toe_abs && x.iode == nav.iode) { list.push(nav); }
        k += 8;
    }
    for list in out.values_mut() {
        list.sort_by(|a, b| a.toe_abs.partial_cmp(&b.toe_abs).unwrap());
    }
    Ok(())
}

fn f(v: f64) -> String {
    // RINEX D19.12, written with E: sign, d.dddddddddddd, E, sign, two digits
    let s = format!("{:.12E}", v);
    let (mant, exp) = s.split_once('E').unwrap();
    let e: i32 = exp.parse().unwrap();
    format!("{:>19}", format!("{mant}E{}{:02}", if e < 0 { '-' } else { '+' }, e.abs()))
}

/// RINEX 3.04 GPS navigation text from (PRN, record) pairs.
pub fn write_nav(entries: &[(u8, Nav)], comments: &[String]) -> String {
    let (y, mo, d, h, mi, s) = calendar(crate::time::now_gps() + crate::time::UTC_GPS);
    let mut out = String::new();
    out += "     3.04           N: GNSS NAV DATA    G: GPS              RINEX VERSION / TYPE\n";
    out += &format!("{:<20}{:<20}{:<20}PGM / RUN BY / DATE\n", "open-extended-ephem", "",
                    format!("{y:04}{mo:02}{d:02} {h:02}{mi:02}{:02} UTC", s as i64));
    for c in comments { out += &format!("{:<60}COMMENT\n", &c[..c.len().min(60)]); }
    out += &format!("{:<60}END OF HEADER\n", "");
    for (sv, r) in entries {
        let (y, mo, d, h, mi, s) = calendar(r.toc_abs);
        out += &format!("G{sv:02} {y:04} {mo:02} {d:02} {h:02} {mi:02} {:02}{}{}{}\n",
                        s.round() as i64, f(r.af0), f(r.af1), f(r.af2));
        let p = &r.p;
        let rows = [
            [r.iode, p[12], p[7], p[2]],
            [p[9], p[1], p[10], p[0]],
            [r.toe(), p[13], p[5], p[14]],
            [p[4], p[11], p[3], p[6]],
            [p[8], 0.0, r.week(), 0.0],
            [r.svacc, r.health, r.tgd, r.iodc],
        ];
        for row in rows {
            out += "    ";
            for v in row { out += &f(v); }
            out += "\n";
        }
        out += &format!("    {}{}\n", f(r.ttx), f(4.0));
    }
    out
}
