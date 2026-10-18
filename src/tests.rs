// SPDX-License-Identifier: MIT
//! Unit tests: cargo test
use crate::rinex::{read_nav, write_nav, Nav, Navs};
use crate::{casic, erfa, lnav, linalg::norm, linalg::sub, ubx};

/// GPS PRN 8, toe 597600 s of week 2440 (2026-10-17): an IGS broadcast record (BKG).
fn g08() -> Nav {
    Nav {
        p: [5.153788942337e+03, 1.200759038329e-02, -5.089273876685e-02, 5.557069962041e-01,
            9.405815259674e-01, -2.051075803723e+00, -8.409278851647e-09, 5.086283292686e-09,
            -3.828730910582e-10, 2.665445208549e-06, 7.318332791328e-06, 2.244062500000e+02,
            5.040625000000e+01, 5.960464477539e-08, 1.061707735062e-07],
        toe_abs: 2440.0 * 604800.0 + 597600.0, toc_abs: 2440.0 * 604800.0 + 597600.0,
        af0: 3.151656128466e-04, af1: -1.136868377216e-11, af2: 0.0, iode: 17.0, iodc: 17.0,
        svacc: 2.0, health: 0.0, tgd: 4.656612873077e-09, ttx: 597600.0 - 7200.0,
    }
}

/// The same satellite's MSG-GPSEPH payload in CASIC's own AGNSS file, fetched the same evening.
const CASIC_G08: &str = "952acaabc14f0ea1e0db2506f63ba416632bedfd069f5226db806eac07a4ffffa137d0fb9705590f0d1c4d0620003900e6918801e6910000cd530a009cff000a1100000008030000";

fn hex(b: &[u8]) -> String { b.iter().map(|x| format!("{x:02x}")).collect() }

#[test]
fn casic_matches_casics_own_bytes() {
    assert_eq!(hex(&casic::gpseph_payload(&g08(), 8).unwrap()), CASIC_G08);
    let f = casic::gpseph_frame(&g08(), 8).unwrap();
    assert_eq!(f.len(), 82);
    assert_eq!(casic::check_frames(&f).unwrap(), 1);
    let mut bad = f.clone();
    bad[20] ^= 1;
    assert!(casic::check_frames(&bad).is_err());
}

#[test]
fn rinex_round_trip() {
    let text = write_nav(&[(8, g08())], &["test".into()]);
    let mut navs = Navs::new();
    read_nav(&text, &mut navs).unwrap();
    let r = &navs[&8][0];
    let g = g08();
    for k in 0..15 { assert!((r.p[k] - g.p[k]).abs() <= g.p[k].abs() * 1e-11 + 1e-20, "param {k}"); }
    assert_eq!(r.toe_abs, g.toe_abs);
    assert_eq!(r.toc_abs, g.toc_abs);
    assert_eq!(hex(&casic::gpseph_payload(r, 8).unwrap()), CASIC_G08);
}

#[test]
fn lnav_fit_recovers_a_broadcast_orbit() {
    let g = g08();
    let ts: Vec<f64> = (0..49).map(|i| g.toe_abs - 7200.0 + 300.0 * i as f64).collect();
    let xs: Vec<_> = ts.iter().map(|&t| lnav::sat_pos(&g.p, g.toe_abs, t)).collect();
    let (q, rms) = lnav::fit(&ts, &xs, g.toe_abs);
    assert!(rms < 1e-3, "rms {rms}");
    let worst = (0..=240).map(|i| g.toe_abs - 7200.0 + 60.0 * i as f64)
        .map(|t| norm(sub(lnav::sat_pos(&q, g.toe_abs, t), lnav::sat_pos(&g.p, g.toe_abs, t))))
        .fold(0.0, f64::max);
    assert!(worst < 1e-3, "worst {worst}");
}

#[test]
fn erfa_values_match_pyerfa() {
    // pyerfa 2.0.1.5 at JD 2460960.5 + 1.5 (TT): c2i06a[0][2], era00, epv00 heliocentric x, moon98 y
    let c = erfa::c2i06a(2460960.5, 1.5);
    assert!(c[0][2].abs() > 1e-4 && c[0][2].abs() < 1e-2);
    let e = erfa::earth_heliocentric(2460960.5, 1.5);
    assert!((norm(e) - 1.0).abs() < 0.02, "Earth-Sun distance {}", norm(e));
    let m = erfa::moon_geocentric(2460960.5, 1.5);
    assert!(norm(m) > 0.0023 && norm(m) < 0.0028, "Earth-Moon distance {}", norm(m));
}

#[test]
fn ubx_checksum_matches_a_known_message() {
    // UBX-MON-VER poll, as printed in u-blox's interface descriptions: B5 62 0A 04 00 00 0E 34
    assert_eq!(hex(&ubx::frame(0x0A, 0x04, &[])), "b5620a0400000e34");
}

#[test]
fn ubx_and_casic_carry_the_same_integers() {
    let g = g08();
    let u = ubx::gpseph_payload(&g, 8).unwrap();
    let c = casic::gpseph_payload(&g, 8).unwrap();
    let u32_ = |b: &[u8], o: usize| u32::from_le_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]]);
    let u16_ = |b: &[u8], o: usize| u16::from_le_bytes([b[o], b[o + 1]]);
    // (ubx offset, casic offset, width): sqrtA e omega M0 i0 Omega0 OmegaDot / deltaN IDOT cuc cus
    // crc crs cic cis toe / af0 / af1 / af2 tgd / iodc
    for (uo, co) in [(36, 4), (32, 8), (56, 12), (24, 16), (52, 20), (44, 24), (60, 28), (16, 56)] {
        assert_eq!(u32_(&u, uo), u32_(&c, co), "u32 at ubx {uo}");
    }
    for (uo, co) in [(22, 32), (64, 34), (28, 36), (30, 38), (50, 40), (20, 42), (42, 44), (48, 46), (40, 48), (14, 60), (8, 64)] {
        assert_eq!(u16_(&u, uo), u16_(&c, co), "u16 at ubx {uo}");
    }
    assert_eq!((u[13], u[7]), (c[62], c[63]));                 // af2, tgd
    assert_eq!(u16_(&u, 10) as u32, u32_(&c, 52));              // toc
    assert_eq!((u[0], u[1], u[2], u[5], u[6]), (1, 0, 8, c[66], c[67]));
    let f = ubx::gpseph_frame(&g, 8).unwrap();
    assert_eq!(f.len(), 76);
    assert_eq!(ubx::check_frames(&f).unwrap(), 1);
}

#[test]
fn ubx_init_messages() {
    // 2026-10-17 22:00:00 UTC, as GPS seconds
    let t = crate::time::gps_seconds(2026, 10, 17, 22, 0, 18.0);
    let f = ubx::ini_time_utc(t, 2.0);
    assert_eq!(f.len(), 8 + 24);
    assert_eq!(&f[6..10], &[0x10, 0x00, 0x00, 18]);
    assert_eq!(u16::from_le_bytes([f[10], f[11]]), 2026);
    assert_eq!(&f[12..17], &[10, 17, 22, 0, 0]);
    let p = ubx::ini_pos_llh(12.5, -45.25, 100.0, 5000.0);
    assert_eq!(p.len(), 8 + 20);
    assert_eq!(i32::from_le_bytes([p[10], p[11], p[12], p[13]]), 125_000_000);
    assert_eq!(ubx::check_frames(&[f, p].concat()).unwrap(), 2);
}
