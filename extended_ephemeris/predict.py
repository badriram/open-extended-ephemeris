# SPDX-License-Identifier: MIT
"""Extended ephemeris for GPS: predict the satellites' orbit data a week ahead, so a receiver
cold-starts in seconds for days without a network.

    python3 -m extended_ephemeris --out out/                      # newest IGS data, 7 days
    python3 -m extended_ephemeris --sp3 F --erp F --rinex F [F...] --out out/

Per GPS satellite:

  1. Take the newest IGS ultra-rapid orbit (24 h measured, 24 h predicted by the IGS) and its
     Earth rotation parameters, from BKG's data centre (anonymous HTTPS).
  2. Fit a force model (orbit.py) to that 48 h arc and integrate it a week ahead.
  3. Every 2 h, fit a broadcast ephemeris (lnav.py) to 4 h of trajectory: to the IGS's own
     prediction while it lasts, then to ours.
  4. Calibrate against the satellite's real broadcasts over the measured day: the antenna offset
     (IGS orbits are centre of mass) and the clock reference. The clock is a straight line
     through the IGS clock estimates.
  5. Write the records as RINEX 3 navigation data and as CASIC MSG-GPSEPH frames, with an
     accuracy (URA) that grows with the prediction's age.

Output, in --out: gps-7d.rnx (RINEX 3.04), gps-7d-casic.bin (CASIC frames, 223 kB), gps-7d.json
(when and from what it was made, and a self-check against the satellites' broadcasts) and
fit-state.json (for --warm next time). A receiver wants, for each satellite, the record whose
toe is nearest now.

Backtest (IGS file of 2026-10-13 06:00, against the satellites' broadcasts over the next 4.7
days): range error, median / 95th percentile, 0.3 / 1.1 m in the first half day, 1.1 / 5.1 m at
2 days, 4.7 / 17 m at 4.5-5 days; worst 25 m.
"""
import argparse, datetime as dt, gzip, io, json, math, os, re, sys, time, urllib.request
from multiprocessing import Pool

import numpy as np

from . import casic, lnav, orbit, rinex

BKG = "https://igs.bkg.bund.de/root_ftp/IGS"
GPS_EPOCH = dt.datetime(1980, 1, 6, tzinfo=dt.timezone.utc)
STEP = 7200                                   # one ephemeris every 2 h, as the satellites do
DAYS = 7
C = 299792458.0
# 95th-percentile range error by age in the backtest, as metres for the URA index.
URA_BY_DAY = [1.9, 5.1, 8.6, 12.1, 17.2, 25.0, 35.0]
WARM_MAX_RMS = 6.0       # a warm-started force-model fit worse than this is redone from scratch


def gps_s(t):
    return (t - GPS_EPOCH).total_seconds()


def utc_of(gps_seconds):
    return GPS_EPOCH + dt.timedelta(seconds=gps_seconds + orbit.UTC_GPS)


# ---------------------------------------------------------------- inputs

def fetch(url):
    data = urllib.request.urlopen(url, timeout=120).read()
    return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data


def read_text(where):
    data = fetch(where) if "://" in where else open(where, "rb").read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("ascii", "replace")


def latest_igu(now):
    """URLs of the newest IGS ultra-rapid orbit and ERP files on BKG."""
    week = int(gps_s(now) // lnav.WEEK)
    for w in (week, week - 1):
        try:
            listing = fetch(f"{BKG}/products/{w}/").decode("ascii", "replace")
        except Exception:                     # noqa: BLE001
            continue
        names = sorted(set(re.findall(r"IGS0OPSULT_\d{11}_02D_15M_ORB\.SP3\.gz", listing)))
        if names:
            orb = names[-1]
            return f"{BKG}/products/{w}/{orb}", f"{BKG}/products/{w}/{orb.replace('15M_ORB.SP3', '01D_ERP.ERP')}"
    raise RuntimeError("no IGS ultra-rapid orbit found on BKG")


def brdc_urls(day):
    y, doy = day.year, day.timetuple().tm_yday
    base = f"{BKG}/BRDC/{y}/{doy:03d}/"
    return [base + f"BRDC00WRD_{k}_{y}{doy:03d}0000_01D_MN.rnx.gz" for k in ("S", "R")] + \
           [base + f"BRDC00IGS_R_{y}{doy:03d}0000_01D_MN.rnx.gz"]


def read_sp3(text):
    """{prn: (times, positions m, clocks s or nan, predicted flags)} for GPS."""
    data, t = {}, None
    for line in text.splitlines():
        if line.startswith("*"):
            f = line[1:].split()
            t = gps_s(dt.datetime(*map(int, f[:5]), int(float(f[5])), tzinfo=dt.timezone.utc))
        elif line.startswith("PG") and t is not None:
            x, y, z, c = (float(line[4 + 14 * k:18 + 14 * k]) for k in range(4))
            if abs(x) < 1e-6:
                continue
            d = data.setdefault(int(line[2:4]), ([], [], [], []))
            d[0].append(t)
            d[1].append(np.array([x, y, z]) * 1e3)
            d[2].append(c * 1e-6 if abs(c) < 999999 else float("nan"))
            d[3].append(len(line) > 75 and line[75] == "P")
    return {sv: tuple(np.array(v) for v in d) for sv, d in data.items()}


def broadcasts(texts):
    out = {}
    for text in texts:
        rinex.read_nav(text, out)
    return out


def nearest(recs, t):
    r = min(recs, key=lambda r: abs(r["toe_abs"] - t))
    return r if abs(r["toe_abs"] - t) <= 7200 else None


def clock(r, t):
    d = t - r["toc_abs"]
    return r["af0"] + r["af1"] * d + r["af2"] * d * d


def rac(x, v):
    r = x / np.linalg.norm(x)
    c = np.cross(x, v)
    c /= np.linalg.norm(c)
    return np.array([r, np.cross(c, r), c])


def lagrange(ts, ys, t, order=10):
    i = int(np.searchsorted(ts, t))
    lo = max(0, min(len(ts) - order, i - order // 2))
    xs = ts[lo:lo + order]
    out = 0
    for j in range(order):
        w = 1.0
        for m in range(order):
            if m != j:
                w *= (t - xs[m]) / (xs[j] - xs[m])
        out = out + w * ys[lo + j]
    return out


# ---------------------------------------------------------------- one satellite

def predict_sat(job):
    sv, (ts, xs, cs, pf), recs, erp, t_end, warm = job
    t_begin = time.time()
    pred_start = ts[pf][0] if pf.any() else ts[len(ts) // 2]
    eo = orbit.EarthOrientation(erp, min(ts[0], warm["t0"]) if warm else ts[0], t_end)
    model = orbit.Orbit(eo, ts[0])
    xg = np.array([eo.c2t(t).T @ x for t, x in zip(ts, xs)])
    started = "cold"
    if warm:
        # The previous run's fitted state, carried forward to this arc's start with its own forces.
        prev = orbit.Orbit(eo, warm["t0"])
        y0 = np.array(warm["state"])
        if warm["t0"] < ts[0]:
            y0 = prev.state_at(y0, np.array(warm["srp"]), ts[0])
        state, srp, fit_rms, nfev, conv = model.fit(ts, xg, start=(y0, np.array(warm["srp"])))
        started = "warm"
        if fit_rms > WARM_MAX_RMS or not conv:    # a manoeuvre, or a bad start: fit from scratch
            started = "warm-rejected"
    if started != "warm":
        state, srp, fit_rms, nfev, conv = model.fit(ts, xg)

    # calibration over the measured day: where the broadcast puts the satellite and its clock,
    # relative to the IGS centre of mass and clock
    d_rac, d_clk = [], []
    for k in range(1, len(ts) - 1, 4):
        r = nearest(recs, ts[k]) if not pf[k] else None
        if r is None:
            continue
        v = (xs[k + 1] - xs[k - 1]) / (ts[k + 1] - ts[k - 1])
        d_rac.append(rac(xs[k], v) @ (lnav.sat_pos(r["p"], r["toe_abs"], ts[k]) - xs[k]))
        if not np.isnan(cs[k]):
            d_clk.append(clock(r, ts[k]) - cs[k])
    if len(d_rac) < 10:
        return sv, None
    off_rac = np.median(d_rac, axis=0)
    ok = ~np.isnan(cs)
    latest = recs[-1]
    if ok.sum() >= 20 and len(d_clk) >= 10:
        c1, c0 = np.polyfit(ts[ok] - pred_start, cs[ok], 1)
        c0 += float(np.median(d_clk))
    else:                                       # no IGS clock: extend the latest broadcast's
        c1, c0 = latest["af1"], clock(latest, pred_start)

    grid = np.arange(ts[0], t_end + STEP + 1, 300.0)
    xe = np.array([eo.c2t(t) @ x for t, x in zip(grid, model.propagate(state, srp, grid))])

    out, fits = [], []
    toe = math.ceil(pred_start / STEP) * STEP
    while toe <= t_end:
        win = np.arange(toe - 7200, toe + 7201, 300.0)
        if win[-1] <= ts[-1]:                   # inside the IGS arc: fit to the IGS's own orbit
            pos = [lagrange(ts, xs, t) for t in win]
        else:
            pos = [xe[int(round((t - grid[0]) / 300))] for t in win]
        # move centre of mass to where the broadcasts put the satellite
        pos = [p + rac(p, (pos[min(i + 1, len(pos) - 1)] - pos[max(i - 1, 0)])).T @ off_rac
               for i, p in enumerate(pos)]
        p, rms = lnav.fit(win, pos, toe)
        fits.append(rms)
        age_days = (toe - pred_start) / 86400.0
        rec = dict(zip(lnav.KEYS, p))
        iodc = 768 + int((toe // STEP) % 240)
        rec.update(toe=toe % lnav.WEEK, week=int(toe // lnav.WEEK), toc=toe % lnav.WEEK, toc_abs=toe,
                   af0=c0 + c1 * (toe - pred_start), af1=c1, af2=0.0, TGD=latest["TGD"],
                   IODC=iodc, IODE=iodc & 0xFF, health=0,
                   svacc=URA_BY_DAY[min(int(max(age_days, 0)), len(URA_BY_DAY) - 1)])
        out.append((toe, rec, p))
        toe += STEP
    return sv, dict(records=out, fit_rms=fit_rms, eph_rms=float(np.median(fits)), started=started,
                    iterations=nfev, converged=conv, seconds=time.time() - t_begin,
                    fit_state=dict(t0=float(ts[0]), state=[float(v) for v in state],
                                   srp=[float(v) for v in srp], rms=fit_rms))


# ---------------------------------------------------------------- self-check

def self_check(results, recs_by_sv):
    """Compare the first predicted records with the satellites' real broadcasts where both exist."""
    errs = []
    for sv, res in results.items():
        for toe, rec, p in res["records"][:3]:
            for t in np.arange(toe - 3600, toe + 3600, 900):
                r = nearest(recs_by_sv[sv], t)
                if r is None:
                    continue
                x_t = lnav.sat_pos(r["p"], r["toe_abs"], t)
                d = lnav.sat_pos(p, toe, t) - x_t
                dr = float(d @ (x_t / np.linalg.norm(x_t)))
                dc = C * ((rec["af0"] + rec["af1"] * (t - toe)) - clock(r, t))
                errs.append(math.sqrt((0.98 * dr - dc) ** 2 + (float(d @ d) - dr * dr) / 49))
    return errs


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=".", help="directory for the outputs")
    ap.add_argument("--sp3", help="IGS ultra-rapid orbit (path or URL; default: newest on BKG)")
    ap.add_argument("--erp", help="its Earth rotation parameters (path or URL)")
    ap.add_argument("--rinex", nargs="+", help="broadcast RINEX covering the orbit's measured day")
    ap.add_argument("--days", type=float, default=DAYS)
    ap.add_argument("--max-check-m", type=float, default=5.0,
                    help="fail if the self-check's median range error exceeds this")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--warm", help="a previous run's fit-state.json (path or URL) to start the fits "
                                   "from; anything missing or unreadable falls back to a full fit")
    a = ap.parse_args()

    now = dt.datetime.now(dt.timezone.utc)
    sp3_url, erp_url = (a.sp3, a.erp) if a.sp3 else latest_igu(now)
    if a.sp3 and not a.erp:
        sys.exit("--erp is needed with --sp3")
    print("orbit:", sp3_url, "\nERP:  ", erp_url, file=sys.stderr)
    sp3 = read_sp3(read_text(sp3_url))
    erp = orbit.read_erp(read_text(erp_url))
    t0 = min(v[0][0] for v in sp3.values())
    t_arc_end = max(v[0][-1] for v in sp3.values())

    texts, used = [], []
    if a.rinex:
        texts, used = [read_text(f) for f in a.rinex], a.rinex
    else:
        day = utc_of(t0).date()
        while day <= now.date():
            for url in brdc_urls(day):
                try:
                    texts.append(fetch(url).decode("ascii", "replace"))
                    used.append(url)
                    break
                except Exception:                 # noqa: BLE001 - try the next name
                    continue
            day += dt.timedelta(days=1)
    recs = broadcasts(texts)
    print("broadcasts:", *used, sep="\n  ", file=sys.stderr)

    warm = {}
    if a.warm:
        try:
            warm = {int(k): v for k, v in json.loads(read_text(a.warm))["satellites"].items()}
        except Exception as e:                    # noqa: BLE001 - a cold start is always possible
            print("warm start unavailable:", e, file=sys.stderr)
    pred_start = min(v[0][v[3]][0] for v in sp3.values() if v[3].any())
    t_end = pred_start + a.days * 86400
    jobs = [(sv, d, recs[sv], erp, t_end, warm.get(sv)) for sv, d in sorted(sp3.items())
            if sv in recs and recs[sv][-1]["health"] == 0]
    with Pool(a.jobs) as pool:
        results = {sv: r for sv, r in pool.map(predict_sat, jobs) if r is not None}

    entries = sorted(((sv, rec) for sv, r in results.items() for _, rec, _ in r["records"]),
                     key=lambda e: (e[1]["toc_abs"], e[0]))
    blob = b"".join(casic.gpseph_frame(rec, sv) for sv, rec in entries)
    kinds = casic.split_frames(blob)
    errs = self_check(results, recs)
    check = dict(samples=len(errs), median_m=float(np.median(errs)) if errs else None,
                 p95_m=float(np.percentile(errs, 95)) if errs else None)
    meta = dict(
        made=now.isoformat(timespec="seconds"),
        igs_orbit=sp3_url, igs_erp=erp_url, broadcasts=used,
        prediction_starts=utc_of(pred_start).isoformat(), valid_until=utc_of(t_end).isoformat(),
        igs_prediction_until=utc_of(t_arc_end).isoformat(),
        satellites=sorted(results), records=len(entries), casic_bytes=len(blob),
        force_model_fit_m=dict(median=float(np.median([r["fit_rms"] for r in results.values()])),
                               worst=float(max(r["fit_rms"] for r in results.values()))),
        self_check_vs_broadcast=check,
        backtest="range error median/95%: 0.3/1.1 m first half day, 1.1/5.1 m at 2 days, "
                 "4.7/17 m at 4.5-5 days, worst 25 m (IGS file 2026-10-13 06:00)")
    meta["fits"] = {s: sum(1 for r in results.values() if r["started"] == s)
                    for s in ("cold", "warm", "warm-rejected")}
    meta["fits"]["not_converged"] = sum(1 for r in results.values() if not r["converged"])
    meta["fits"]["iterations_median"] = float(np.median([r["iterations"] for r in results.values()]))
    meta["satellite_seconds"] = dict(median=float(np.median([r["seconds"] for r in results.values()])),
                                     total=float(sum(r["seconds"] for r in results.values())))
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "fit-state.json"), "w") as f:
        json.dump(dict(made=meta["made"], igs_orbit=sp3_url,
                       satellites={str(sv): r["fit_state"] for sv, r in results.items()}), f)
    with open(os.path.join(a.out, "gps-7d-casic.bin"), "wb") as f:
        f.write(blob)
    with open(os.path.join(a.out, "gps-7d.rnx"), "w") as f:
        f.write(rinex.write_nav(entries, comments=[
            "PREDICTED broadcast ephemeris (extended ephemeris), not",
            "broadcast by the satellites. From IGS ultra-rapid orbit",
            os.path.basename(sp3_url)[:60],
            "github.com/badriram/open-extended-ephemeris  MIT"]))
    with open(os.path.join(a.out, "gps-7d.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta, indent=2))
    if not errs or check["median_m"] > a.max_check_m:
        sys.exit(f"self-check failed: {check}")


if __name__ == "__main__":
    main()
