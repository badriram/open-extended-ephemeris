// SPDX-License-Identifier: MIT
//! CASIC MSG-GPSEPH frames, for AT6558-based receivers (Quectel L76K, ATGM336H and others).
//!
//! The 72-byte payload carries the LNAV ephemeris at its broadcast resolution; its first word,
//! "reserved" in CASIC's specification, is the sum of the payload's other 32-bit words. Verified
//! byte for byte against CASIC's own AGNSS file, and on an L76K, by
//! https://github.com/badriram/l76k-at6558-agnss (which has casic_eph_select() in C to pick each
//! satellite's frame for now from a file of many).
use crate::rinex::Nav;

const PI: f64 = 3.1415926535898;   // IS-GPS-200's pi, which semicircles are defined with

pub fn checksum(cls: u8, id: u8, payload: &[u8]) -> u32 {
    let mut ck = ((id as u32) << 24).wrapping_add((cls as u32) << 16).wrapping_add(payload.len() as u32);
    for w in payload.chunks_exact(4) {
        ck = ck.wrapping_add(u32::from_le_bytes([w[0], w[1], w[2], w[3]]));
    }
    ck
}

pub fn frame(cls: u8, id: u8, payload: &[u8]) -> Vec<u8> {
    let mut f = vec![0xBA, 0xCE];
    f.extend_from_slice(&(payload.len() as u16).to_le_bytes());
    f.extend_from_slice(&[cls, id]);
    f.extend_from_slice(payload);
    f.extend_from_slice(&checksum(cls, id, payload).to_le_bytes());
    f
}

pub fn ura_index(metres: f64) -> u8 {
    let tops = [2.4, 3.4, 4.85, 6.85, 9.65, 13.65, 24.0, 48.0, 96.0, 192.0, 384.0, 768.0, 1536.0, 3072.0, 6144.0];
    tops.iter().position(|&t| metres <= t).unwrap_or(15) as u8
}

fn q(v: f64, lsb: f64, lo: f64, hi: f64) -> Result<i64, String> {
    let r = (v / lsb).round();
    if r < lo || r > hi { return Err(format!("{v} out of range at {lsb}")); }
    Ok(r as i64)
}

/// The 72-byte MSG-GPSEPH payload for a record.
pub fn gpseph_payload(r: &Nav, sv: u8) -> Result<[u8; 72], String> {
    let s = 1.0 / PI;
    let p = &r.p;
    let (i32r, u32r, i16r, u16r, i8r) = ((i32::MIN as f64, i32::MAX as f64), (0.0, u32::MAX as f64),
        (i16::MIN as f64, i16::MAX as f64), (0.0, u16::MAX as f64), (i8::MIN as f64, i8::MAX as f64));
    let mut b = [0u8; 72];
    let mut put = |off: usize, bytes: &[u8]| b[off..off + bytes.len()].copy_from_slice(bytes);
    put(4, &(q(p[0], 2f64.powi(-19), u32r.0, u32r.1)? as u32).to_le_bytes());
    put(8, &(q(p[1], 2f64.powi(-33), u32r.0, u32r.1)? as u32).to_le_bytes());
    put(12, &(q(p[3] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());
    put(16, &(q(p[2] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());
    put(20, &(q(p[4] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());
    put(24, &(q(p[5] * s, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());
    put(28, &(q(p[6] * s, 2f64.powi(-43), i32r.0, i32r.1)? as i32).to_le_bytes());
    put(32, &(q(p[7] * s, 2f64.powi(-43), i16r.0, i16r.1)? as i16).to_le_bytes());
    put(34, &(q(p[8] * s, 2f64.powi(-43), i16r.0, i16r.1)? as i16).to_le_bytes());
    put(36, &(q(p[9], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());
    put(38, &(q(p[10], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());
    put(40, &(q(p[11], 2f64.powi(-5), i16r.0, i16r.1)? as i16).to_le_bytes());
    put(42, &(q(p[12], 2f64.powi(-5), i16r.0, i16r.1)? as i16).to_le_bytes());
    put(44, &(q(p[13], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());
    put(46, &(q(p[14], 2f64.powi(-29), i16r.0, i16r.1)? as i16).to_le_bytes());
    put(48, &(q(r.toe(), 16.0, u16r.0, u16r.1)? as u16).to_le_bytes());
    put(50, &((r.week() as i64 % 1024) as u16).to_le_bytes());
    put(52, &(q(r.toc(), 16.0, u32r.0, u32r.1)? as u32).to_le_bytes());
    put(56, &(q(r.af0, 2f64.powi(-31), i32r.0, i32r.1)? as i32).to_le_bytes());
    put(60, &(q(r.af1, 2f64.powi(-43), i16r.0, i16r.1)? as i16).to_le_bytes());
    put(62, &[(q(r.af2, 2f64.powi(-55), i8r.0, i8r.1)? as i8) as u8]);
    put(63, &[(q(r.tgd, 2f64.powi(-31), i8r.0, i8r.1)? as i8) as u8]);
    put(64, &(r.iodc as u16).to_le_bytes());
    put(66, &[ura_index(r.svacc), r.health as u8, sv, if r.health == 0.0 { 3 } else { 1 }]);
    let sum = b[4..].chunks_exact(4).fold(0u32, |a, w| a.wrapping_add(u32::from_le_bytes([w[0], w[1], w[2], w[3]])));
    b[..4].copy_from_slice(&sum.to_le_bytes());
    Ok(b)
}

pub fn gpseph_frame(r: &Nav, sv: u8) -> Result<Vec<u8>, String> {
    Ok(frame(0x08, 0x07, &gpseph_payload(r, sv)?))
}

/// Check every frame's checksum; returns the number of frames.
pub fn check_frames(buf: &[u8]) -> Result<usize, String> {
    let (mut j, mut n) = (0, 0);
    while j + 6 <= buf.len() {
        let ln = u16::from_le_bytes([buf[j + 2], buf[j + 3]]) as usize;
        let f = &buf[j..j + 10 + ln];
        let ck = u32::from_le_bytes([f[6 + ln], f[7 + ln], f[8 + ln], f[9 + ln]]);
        if f[0] != 0xBA || f[1] != 0xCE || ck != checksum(f[4], f[5], &f[6..6 + ln]) {
            return Err(format!("bad frame at offset {j}"));
        }
        j += 10 + ln;
        n += 1;
    }
    Ok(n)
}
