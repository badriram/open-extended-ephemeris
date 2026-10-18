// SPDX-License-Identifier: MIT
//! Time scales. Times are continuous GPS seconds since 1980-01-06 00:00 GPS.

pub const WEEK: f64 = 604800.0;
pub const JD_GPS0: f64 = 2444244.5;     // Julian date of the GPS epoch
pub const TT_GPS: f64 = 51.184;         // TT - GPS, s
pub const UTC_GPS: f64 = -18.0;         // UTC - GPS since 2017, s
pub const UNIX_GPS0: f64 = 315964800.0; // the GPS epoch in Unix seconds

/// Days from 1970-01-01 to a civil date (Howard Hinnant's algorithm).
pub fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = if y >= 0 { y } else { y - 399 } / 400;
    let yoe = y - era * 400;
    let doy = (153 * (m + if m > 2 { -3 } else { 9 }) + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146097 + doe - 719468
}

/// GPS seconds of a calendar time given in GPS time (as SP3 and RINEX GPS records are).
pub fn gps_seconds(y: i64, mo: i64, d: i64, h: i64, mi: i64, s: f64) -> f64 {
    (days_from_civil(y, mo, d) * 86400 + h * 3600 + mi * 60) as f64 + s - UNIX_GPS0
}

/// Calendar fields of a GPS-seconds time, in the same time scale.
pub fn calendar(t: f64) -> (i64, i64, i64, i64, i64, f64) {
    let u = t + UNIX_GPS0;
    let days = (u / 86400.0).floor() as i64;
    let sod = u - days as f64 * 86400.0;
    // civil_from_days
    let z = days + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + if m <= 2 { 1 } else { 0 };
    let h = (sod / 3600.0).floor() as i64;
    let mi = ((sod - h as f64 * 3600.0) / 60.0).floor() as i64;
    (y, m, d, h, mi, sod - h as f64 * 3600.0 - mi as f64 * 60.0)
}

/// ISO 8601 UTC string for a GPS-seconds time.
pub fn iso_utc(t: f64) -> String {
    let (y, m, d, h, mi, s) = calendar(t + UTC_GPS);
    format!("{y:04}-{m:02}-{d:02}T{h:02}:{mi:02}:{:02}Z", s.round() as i64 % 60)
}

pub fn now_gps() -> f64 {
    let u = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs_f64();
    u - UNIX_GPS0 - UTC_GPS
}
