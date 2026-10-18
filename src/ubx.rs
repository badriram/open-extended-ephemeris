// SPDX-License-Identifier: MIT
//! u-blox UBX messages for assistance (M8, M9, M10 receivers):
//!
//! - UBX-MGA-GPS-EPH (0x13 0x00, 68 bytes): one GPS satellite's ephemeris, at its broadcast
//!   resolution. It carries no week number: the receiver must know the date (MGA-INI-TIME_UTC).
//! - UBX-MGA-INI-TIME_UTC (0x13 0x40, 24 bytes) and UBX-MGA-INI-POS_LLH (0x13 0x40, 20 bytes): the
//!   time and a rough position, sent first.
//!
//! Layouts as in u-blox's interface descriptions, checked field by field against pyubx2
//! (semuconsulting/pyubx2, ubxtypes_set.py).
use crate::casic::ura_index;
use crate::rinex::Nav;
use crate::time::{calendar, UTC_GPS};

const PI: f64 = 3.1415926535898;   // IS-GPS-200's pi, which semicircles are defined with

/// A UBX frame: sync, class, id, length, payload, 8-bit Fletcher checksum over class..payload.
pub fn frame(cls: u8, id: u8, payload: &[u8]) -> Vec<u8> {
    let mut f = vec![0xB5, 0x62, cls, id];
    f.extend_from_slice(&(payload.len() as u16).to_le_bytes());
    f.extend_from_slice(payload);
    let (mut a, mut b) = (0u8, 0u8);
    for &x in &f[2..] {
        a = a.wrapping_add(x);
        b = b.wrapping_add(a);
    }
    f.push(a);
    f.push(b);
    f
}

/// Check every UBX frame's checksum; returns the number of frames.
pub fn check_frames(buf: &[u8]) -> Result<usize, String> {
    let (mut j, mut n) = (0, 0);
    while j + 8 <= buf.len() {
        if buf[j] != 0xB5 || buf[j + 1] != 0x62 { return Err(format!("bad sync at offset {j}")); }
        let ln = u16::from_le_bytes([buf[j + 4], buf[j + 5]]) as usize;
        let end = j + 8 + ln;
        if end > buf.len() { return Err(format!("short frame at offset {j}")); }
        let f = frame(buf[j + 2], buf[j + 3], &buf[j + 6..j + 6 + ln]);
        if f[..] != buf[j..end] { return Err(format!("bad checksum at offset {j}")); }
        j = end;
        n += 1;
    }
    Ok(n)
}

fn q(v: f64, lsb: f64, lo: f64, hi: f64) -> Result<i64, String> {
    let r = (v / lsb).round();
    if r < lo || r > hi { return Err(format!("{v} out of range at {lsb}")); }
    Ok(r as i64)
}

/// The 68-byte UBX-MGA-GPS-EPH payload for a record.
pub fn gpseph_payload(r: &Nav, sv: u8) -> Result<[u8; 68], String> {
    let s = 1.0 / PI;
    let p = &r.p;
    let i32r = (i32::MIN as f64, i32::MAX as f64);
    let u32r = (0.0, u32::MAX as f64);
    let i16r = (i16::MIN as f64, i16::MAX as f64);
    let u16r = (0.0, u16::MAX as f64);
    let i8r = (i8::MIN as f64, i8::MAX as f64);
    let mut b = [0u8; 68];
    let mut put = |off: usize, bytes: &[u8]| b[off..off + bytes.len()].copy_from_slice(bytes);
    put(0, &[0x01, 0x00, sv, 0, 0, ura_index(r.svacc), r.health as u8]);   // type, version, svId, -, fitInterval 4 h
    put(7, &[(q(r.tgd, 2f64.powi(-31), i8r.0, i8r.1)? as i8) as u8]);
    put(8, &(r.iodc as u16).to_le_bytes());
    put(10, &(q(r.toc(), 16.0, u16r.0, u16r.1)? as u16).to_le_bytes());
    put(13, &[(q(r.af2, 2f64.powi(-55), i8r.0, i8r.1)? as i8) as u8]);
    put(14, &(q(r.af1, 2f64.powi(-43), i16r.0, i16r.1)? as i16).to_le_bytes());
    put(16, &(q(r.af0, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());
    put(20, &(q(p[12], 2f64.powi(-5), i16r.0, i16r.1)? as i16).to_le_bytes());        // crs
    put(22, &(q(p[7] * s, 2f64.powi(-43), i16r.0, i16r.1)? as i16).to_le_bytes());    // deltaN
    put(24, &(q(p[2] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());    // M0
    put(28, &(q(p[9], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());        // cuc
    put(30, &(q(p[10], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());       // cus
    put(32, &(q(p[1], 2f64.powi(-33), u32r.0, u32r.1)? as u32).to_le_bytes());        // e
    put(36, &(q(p[0], 2f64.powi(-19), u32r.0, u32r.1)? as u32).to_le_bytes());        // sqrtA
    put(40, &(q(r.toe(), 16.0, u16r.0, u16r.1)? as u16).to_le_bytes());
    put(42, &(q(p[13], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());       // cic
    put(44, &(q(p[5] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());    // Omega0
    put(48, &(q(p[14], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());       // cis
    put(50, &(q(p[11], 2f64.powi(-5), i16r.0, i16r.1)? as i16).to_le_bytes());        // crc
    put(52, &(q(p[4] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());    // i0
    put(56, &(q(p[3] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());    // omega
    put(60, &(q(p[6] * s, 2f64.powi(-43), i32r.0, i32r.1)? as i32).to_le_bytes());    // OmegaDot
    put(64, &(q(p[8] * s, 2f64.powi(-43), i16r.0, i16r.1)? as i16).to_le_bytes());    // IDOT
    Ok(b)
}

pub fn gpseph_frame(r: &Nav, sv: u8) -> Result<Vec<u8>, String> {
    Ok(frame(0x13, 0x00, &gpseph_payload(r, sv)?))
}

/// UBX-MGA-INI-TIME_UTC for a GPS-seconds time, accurate to acc_s seconds.
pub fn ini_time_utc(t_gps: f64, acc_s: f64) -> Vec<u8> {
    let (y, mo, d, h, mi, s) = calendar(t_gps + UTC_GPS);
    let whole = s.floor();
    let ns = ((s - whole) * 1e9).round().min(999_999_999.0) as u32;
    let mut b = vec![0x10, 0x00, 0x00, (-UTC_GPS) as i8 as u8];   // type, version, ref: on receipt, leap seconds
    b.extend_from_slice(&(y as u16).to_le_bytes());
    b.extend_from_slice(&[mo as u8, d as u8, h as u8, mi as u8, whole as u8, 0]);
    b.extend_from_slice(&ns.to_le_bytes());
    let acc_whole = acc_s.floor();
    b.extend_from_slice(&(acc_whole.min(65535.0) as u16).to_le_bytes());
    b.extend_from_slice(&[0, 0]);
    b.extend_from_slice(&(((acc_s - acc_whole) * 1e9).round() as u32).to_le_bytes());
    frame(0x13, 0x40, &b)
}

/// UBX-MGA-INI-POS_LLH: latitude and longitude in degrees, altitude in metres, accuracy in metres.
pub fn ini_pos_llh(lat: f64, lon: f64, alt_m: f64, acc_m: f64) -> Vec<u8> {
    let mut b = vec![0x01, 0x00, 0x00, 0x00];
    b.extend_from_slice(&((lat * 1e7).round() as i32).to_le_bytes());
    b.extend_from_slice(&((lon * 1e7).round() as i32).to_le_bytes());
    b.extend_from_slice(&((alt_m * 100.0).round() as i32).to_le_bytes());
    b.extend_from_slice(&((acc_m * 100.0).round().max(0.0) as u32).to_le_bytes());
    frame(0x13, 0x40, &b)
}

/// Galileo SISA index for an accuracy in metres (Galileo OS SIS ICD): 1 cm steps to 0.5 m, 2 cm
/// to 1 m, 4 cm to 2 m, 16 cm to 6 m. Beyond 6 m the ICD has only "no accuracy prediction
/// available" (255), which receivers may treat as unusable, so predictions are capped at 6 m.
pub fn sisa_index(metres: f64) -> u8 {
    let cm = (metres * 100.0).max(0.0);
    (if cm < 50.0 { cm } else if cm < 100.0 { 50.0 + (cm - 50.0) / 2.0 }
     else if cm < 200.0 { 75.0 + (cm - 100.0) / 4.0 } else { (100.0 + (cm - 200.0) / 16.0).min(125.0) }).round() as u8
}

/// The 76-byte UBX-MGA-GAL-EPH payload (0x13 0x02) for a Galileo I/NAV record. Galileo's own
/// resolutions: toe and toc in 60 s, af0 2^-34 s, af1 2^-46, af2 2^-59, BGD 2^-32 s.
pub fn galeph_payload(r: &Nav, prn: u8) -> Result<[u8; 76], String> {
    let s = 1.0 / PI;
    let p = &r.p;
    let i32r = (i32::MIN as f64, i32::MAX as f64);
    let u32r = (0.0, u32::MAX as f64);
    let i16r = (i16::MIN as f64, i16::MAX as f64);
    let u16r = (0.0, u16::MAX as f64);
    let i8r = (i8::MIN as f64, i8::MAX as f64);
    let mut b = [0u8; 76];
    let mut put = |off: usize, bytes: &[u8]| b[off..off + bytes.len()].copy_from_slice(bytes);
    put(0, &[0x01, 0x00, prn, 0]);                                                   // type, version, svId
    put(4, &((r.iodc as u16) & 0x3FF).to_le_bytes());                                // IODnav
    put(6, &(q(p[7] * s, 2f64.powi(-43), i16r.0, i16r.1)? as i16).to_le_bytes());    // deltaN
    put(8, &(q(p[2] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());    // M0
    put(12, &(q(p[1], 2f64.powi(-33), u32r.0, u32r.1)? as u32).to_le_bytes());       // e
    put(16, &(q(p[0], 2f64.powi(-19), u32r.0, u32r.1)? as u32).to_le_bytes());       // sqrtA
    put(20, &(q(p[5] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());   // Omega0
    put(24, &(q(p[4] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());   // i0
    put(28, &(q(p[3] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());   // omega
    put(32, &(q(p[6] * s, 2f64.powi(-43), i32r.0, i32r.1)? as i32).to_le_bytes());   // OmegaDot
    put(36, &(q(p[8] * s, 2f64.powi(-43), i16r.0, i16r.1)? as i16).to_le_bytes());   // IDOT
    put(38, &(q(p[9], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());       // cuc
    put(40, &(q(p[10], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());      // cus
    put(42, &(q(p[11], 2f64.powi(-5), i16r.0, i16r.1)? as i16).to_le_bytes());       // crc
    put(44, &(q(p[12], 2f64.powi(-5), i16r.0, i16r.1)? as i16).to_le_bytes());       // crs
    put(46, &(q(p[13], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());      // cic
    put(48, &(q(p[14], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());      // cis
    put(50, &(q(r.toe(), 60.0, u16r.0, u16r.1)? as u16).to_le_bytes());
    put(52, &(q(r.af0, 2f64.powi(-34), i32r.0, i32r.1)? as i32).to_le_bytes());
    put(56, &(q(r.af1, 2f64.powi(-46), i32r.0, i32r.1)? as i32).to_le_bytes());
    put(60, &[(q(r.af2, 2f64.powi(-59), i8r.0, i8r.1)? as i8) as u8, sisa_index(r.svacc)]);
    put(62, &(q(r.toc(), 60.0, u16r.0, u16r.1)? as u16).to_le_bytes());
    put(64, &(q(r.tgd, 2f64.powi(-32), i16r.0, i16r.1)? as i16).to_le_bytes());      // BGD E1/E5b
    put(68, &[0, 0, 0, 0]);                                                          // health and data validity: OK
    Ok(b)
}

pub fn galeph_frame(r: &Nav, prn: u8) -> Result<Vec<u8>, String> {
    Ok(frame(0x13, 0x02, &galeph_payload(r, prn)?))
}
