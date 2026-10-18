// SPDX-License-Identifier: MIT
//! Unit tests: cargo test
use crate::rinex::{read_nav, write_nav, Nav, Navs};
use crate::{casic, erfa, lnav, linalg::norm, linalg::sub};

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
