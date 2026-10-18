# SPDX-License-Identifier: MIT
"""RINEX 3 navigation files, GPS records only: read the satellites' broadcasts, write predictions.

A GPS record is an epoch line (PRN, the clock reference time toc, af0, af1, af2) and seven lines of
four values: IODE Crs deltaN M0 / Cuc e Cus sqrtA / toe Cic Omega0 Cis / i0 Crc omega OmegaDot /
IDOT L2codes week L2P / accuracy health TGD IODC / transmission-time fit-interval."""
import datetime as dt

import numpy as np

from . import lnav

GPS_EPOCH = dt.datetime(1980, 1, 6)
NAMES = ["af0", "af1", "af2", "IODE", "Crs", "deltaN", "M0", "Cuc", "e", "Cus", "sqrtA", "toe",
         "Cic", "Omega0", "Cis", "i0", "Crc", "omega", "OmegaDot", "IDOT", "L2codes", "week",
         "L2P", "svacc", "health", "TGD", "IODC", "ttx", "fit"]


def _num(s):
    s = s.strip().replace("D", "E").replace("d", "e")
    return float(s) if s else 0.0


def read_nav(text, out=None):
    """GPS records from RINEX 3 navigation text: {prn: [record]}, each with continuous-time
    toe_abs and toc_abs (GPS seconds since 1980-01-06) and 'p', the orbit in lnav.KEYS order."""
    lines = text.splitlines()
    k = next(n for n, l in enumerate(lines) if "END OF HEADER" in l) + 1
    out = {} if out is None else out
    while k < len(lines):
        l = lines[k]
        if not l.startswith("G") or k + 7 >= len(lines):
            k += 1
            continue
        vals = [_num(l[23 + 19 * j:42 + 19 * j]) for j in range(3)]
        for m in range(1, 8):
            b = lines[k + m]
            vals += [_num(b[4 + 19 * j:23 + 19 * j]) for j in range(4)]
        r = dict(zip(NAMES, vals))
        epoch = dt.datetime(*map(int, l[4:23].split()))
        r["toc_abs"] = (epoch - GPS_EPOCH).total_seconds()
        r["toc"] = r["toc_abs"] % lnav.WEEK
        r["toe_abs"] = r["week"] * lnav.WEEK + r["toe"]
        r["p"] = np.array([r[key] for key in lnav.KEYS])
        sv = int(l[1:3])
        lst = out.setdefault(sv, [])
        if not any(x["toe_abs"] == r["toe_abs"] and x["IODE"] == r["IODE"] for x in lst):
            lst.append(r)
        k += 8
    for lst in out.values():
        lst.sort(key=lambda r: r["toe_abs"])
    return out


def _f(v):
    return f"{v:19.12E}"


def write_nav(entries, program="open-extended-ephemeris", comments=()):
    """RINEX 3.04 GPS navigation text from (prn, record) pairs, records with lnav.KEYS plus
    toc_abs, week, toe (seconds of week), af0, af1, af2, IODE, IODC, svacc, health, TGD."""
    now = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%d %H%M%S UTC")
    out = ["     3.04           N: GNSS NAV DATA    G: GPS              RINEX VERSION / TYPE",
           f"{program[:20]:<20}{'':<20}{now:<20}PGM / RUN BY / DATE"]
    out += [f"{c[:60]:<60}COMMENT" for c in comments]
    out.append(f"{'':<60}END OF HEADER")
    for sv, r in sorted(entries, key=lambda e: (e[1]["toc_abs"], e[0])):
        t = GPS_EPOCH + dt.timedelta(seconds=r["toc_abs"])
        out.append(f"G{sv:02d} {t:%Y %m %d %H %M %S}" + _f(r["af0"]) + _f(r["af1"]) + _f(r["af2"]))
        rows = [(r["IODE"], r["Crs"], r["deltaN"], r["M0"]),
                (r["Cuc"], r["e"], r["Cus"], r["sqrtA"]),
                (r["toe"], r["Cic"], r["Omega0"], r["Cis"]),
                (r["i0"], r["Crc"], r["omega"], r["OmegaDot"]),
                (r["IDOT"], 0.0, r["week"], 0.0),
                (r["svacc"], r["health"], r["TGD"], r["IODC"]),
                (r.get("ttx", r["toe"] - 7200.0) % lnav.WEEK, 4.0)]
        out += ["    " + "".join(_f(v) for v in row) for row in rows]
    return "\n".join(out) + "\n"
