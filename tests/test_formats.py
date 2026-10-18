# SPDX-License-Identifier: MIT
"""Format and fitting tests: python3 -m unittest discover tests"""
import struct
import unittest

import numpy as np

from extended_ephemeris import casic, lnav, rinex

# GPS PRN 8, toe 597600 s of week 2440 (2026-10-17): an IGS broadcast record (BKG), and the same
# satellite's MSG-GPSEPH payload in CASIC's own AGNSS file, fetched the same evening.
G08 = dict(af0=3.151656128466e-04, af1=-1.136868377216e-11, af2=0.0, IODE=17.0, Crs=5.040625000000e+01,
           deltaN=5.086283292686e-09, M0=-5.089273876685e-02, Cuc=2.665445208549e-06, e=1.200759038329e-02,
           Cus=7.318332791328e-06, sqrtA=5.153788942337e+03, toe=597600.0, Cic=5.960464477539e-08,
           Omega0=-2.051075803723e+00, Cis=1.061707735062e-07, i0=9.405815259674e-01, Crc=2.244062500000e+02,
           omega=5.557069962041e-01, OmegaDot=-8.409278851647e-09, IDOT=-3.828730910582e-10, week=2440.0,
           svacc=2.0, health=0.0, TGD=4.656612873077e-09, IODC=17.0, toc=597600.0,
           toc_abs=2440 * 604800.0 + 597600.0)
CASIC_G08 = bytes.fromhex(
    "952acaabc14f0ea1e0db2506f63ba416632bedfd069f5226db806eac07a4ffff"
    "a137d0fb9705590f0d1c4d0620003900e6918801e6910000cd530a009cff000a"
    "1100000008030000")


class Casic(unittest.TestCase):
    def test_matches_casics_own_bytes(self):
        self.assertEqual(casic.gpseph_payload(G08, 8), CASIC_G08)

    def test_frame_checksum(self):
        f = casic.gpseph_frame(G08, 8)
        self.assertEqual(len(f), 82)
        self.assertEqual(casic.split_frames(f), {(8, 7): 1})
        bad = bytearray(f)
        bad[20] ^= 1
        with self.assertRaises(ValueError):
            casic.split_frames(bytes(bad))


class Rinex(unittest.TestCase):
    def test_round_trip(self):
        text = rinex.write_nav([(8, G08)], comments=["test"])
        recs = rinex.read_nav(text)[8]
        self.assertEqual(len(recs), 1)
        r = recs[0]
        for k in lnav.KEYS + ["af0", "af1", "TGD", "IODC", "IODE", "week", "toe"]:
            self.assertAlmostEqual(r[k], G08[k], delta=abs(G08[k]) * 1e-11 + 1e-20, msg=k)
        self.assertEqual(r["toe_abs"], G08["toc_abs"])
        # the re-read record encodes to CASIC's bytes as well
        self.assertEqual(casic.gpseph_payload(r, 8), CASIC_G08)


class Fit(unittest.TestCase):
    def test_recovers_a_broadcast_orbit(self):
        p = np.array([G08[k] for k in lnav.KEYS])
        toe = G08["toc_abs"]
        ts = np.arange(toe - 7200, toe + 7201, 300.0)
        q, rms = lnav.fit(ts, [lnav.sat_pos(p, toe, t) for t in ts], toe)
        worst = max(np.linalg.norm(lnav.sat_pos(q, toe, t) - lnav.sat_pos(p, toe, t))
                    for t in np.arange(toe - 7200, toe + 7201, 60.0))
        self.assertLess(rms, 1e-3)
        self.assertLess(worst, 1e-3)


if __name__ == "__main__":
    unittest.main()
