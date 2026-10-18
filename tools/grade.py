#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Grade a predicted RINEX file against what the satellites actually broadcast.

    python3 tools/grade.py out/gps-7d.rnx BRDC-day1.rnx.gz BRDC-day2.rnx.gz ...

For every predicted record, over its central two hours at 15-minute steps, compare with the
broadcast record a receiver would have used at that moment, and report by age of prediction
(from the first predicted record): signal-in-space range error (median, 95th percentile, worst),
the median orbit error, and the worst satellites. Broadcast files: e.g. BKG's merged
https://igs.bkg.bund.de/root_ftp/IGS/BRDC/<year>/<doy>/BRDC00WRD_S_<year><doy>0000_01D_MN.rnx.gz
"""
import gzip
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from extended_ephemeris import lnav, rinex   # noqa: E402

C = 299792458.0


def text(path):
    data = open(path, "rb").read()
    return (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode("ascii", "replace")


def clock(r, t):
    d = t - r["toc_abs"]
    return r["af0"] + r["af1"] * d + r["af2"] * d * d


def main():
    pred = rinex.read_nav(text(sys.argv[1]))
    truth = {}
    for f in sys.argv[2:]:
        rinex.read_nav(text(f), truth)
    t0 = min(r["toe_abs"] for v in pred.values() for r in v) - 7200
    rows, worst = [], {}
    for sv, recs in pred.items():
        for r in recs:
            for t in np.arange(r["toe_abs"] - 3600, r["toe_abs"] + 3600, 900):
                cand = [x for x in truth.get(sv, []) if abs(x["toe_abs"] - t) <= 7200]
                if not cand:
                    continue
                b = min(cand, key=lambda x: abs(x["toe_abs"] - t))
                xb = lnav.sat_pos(b["p"], b["toe_abs"], t)
                d = lnav.sat_pos(r["p"], r["toe_abs"], t) - xb
                dr = float(d @ (xb / np.linalg.norm(xb)))
                dc = C * (clock(r, t) - clock(b, t))
                sisre = math.sqrt((0.98 * dr - dc) ** 2 + (float(d @ d) - dr * dr) / 49)
                rows.append(((t - t0) / 86400, sisre, float(np.linalg.norm(d))))
                worst[sv] = max(worst.get(sv, 0), sisre)
    if not rows:
        sys.exit("no overlap between the prediction and the broadcasts")
    a = np.array(rows)
    print("days      | range error median / 95% / worst | orbit 3D median")
    for d in np.arange(0, a[:, 0].max(), 0.5):
        s = a[(a[:, 0] >= d) & (a[:, 0] < d + 0.5)]
        if len(s):
            print(f"{d:3.1f}-{d + 0.5:3.1f}   | {np.median(s[:, 1]):6.2f} / {np.percentile(s[:, 1], 95):6.2f} / "
                  f"{s[:, 1].max():6.1f} m        | {np.median(s[:, 2]):6.2f} m")
    print("worst satellites:", ", ".join(f"G{sv:02d} {v:.1f} m" for sv, v in
                                         sorted(worst.items(), key=lambda x: -x[1])[:5]))


if __name__ == "__main__":
    main()
