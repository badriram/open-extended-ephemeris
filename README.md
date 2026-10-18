# open-extended-ephemeris

**Offline rapid cold start for GPS.** A GPS receiver with no network spends 30-40 s on a cold start
downloading each satellite's orbit data at 50 bit/s. Give it that data and it fixes in a few
seconds, but broadcast orbit data expires after about four hours. This project **predicts it a
week ahead**: download one file while you have a network, and cold starts stay fast for days
without one.

This is what the industry calls **extended ephemeris**, sold as u-blox AssistNow Offline,
MediaTek EPO, Broadcom LTO, Qualcomm XTRA and SiRF InstantFix. All of those are closed services.
As far as we could find, this is the first open one: the method, the code and the measurements.

## How good it is

Backtest: predicted from the IGS ultra-rapid orbit of 2026-10-13 06:00 UTC, then graded against
what the satellites actually broadcast over the next 4.7 days ([results](results/2026-10-18-backtest.md),
[`tools/grade.py`](tools/grade.py)). Range error is the satellite's contribution to the error of a
distance measurement; a position is off by roughly twice it.

| days after the IGS orbit | range error, median | 95th percentile | worst | satellite position, median |
| --- | --- | --- | --- | --- |
| 0-1 (the IGS's own prediction) | 0.3-0.5 m | 1.1-1.9 m | 4.1 m | 0.7 m |
| 1-2 | 0.8-1.1 m | 3.2-5.1 m | 7.0 m | 0.9 m |
| 2-3 | 1.5-1.9 m | 6.8-8.6 m | 11 m | 1.3-1.6 m |
| 3-4 | 2.5-3.2 m | 9.9-12 m | 18 m | 2.2-3.1 m |
| 4-5 | 4.0-4.7 m | 15-17 m | 25 m | 4.1-5.0 m |

So after five days offline, a first fix is typically about 10 m out, and the receiver replaces the
prediction with the satellites' own data once it has tracked them for about 30 s. Satellite
clocks, not orbits, are now most of the error. Days 5-7 have not been graded yet, and no receiver
has yet been given predicted records (they are ordinary records with a current reference time,
which receivers accept from other sources).

## How it works

Per GPS satellite, with plain physics and least squares (no machine learning):

1. **Start from the IGS.** The International GNSS Service's ultra-rapid orbit gives 24 h of
   measured positions and clocks and 24 h of its own prediction; BKG serves it anonymously.
2. **Fit a force model** to those 48 h: Earth's gravity field (EGM2008 to degree and order 4), the
   Sun and Moon, and solar radiation pressure as five ECOM-style terms with Earth's shadow and a
   smooth penumbra; Earth orientation from the IAU 2006/2000A model (ERFA) with the IGS's polar
   motion and UT1. Fitted by Levenberg-Marquardt, with a priori limits on the minor radiation
   terms, to 0.1-0.26 m RMS.
3. **Integrate a week ahead**, with fixed steps so the result is a smooth function of the inputs.
4. **Fit an ordinary broadcast ephemeris** (IS-GPS-200 LNAV) to every 2 h of the trajectory: to the
   IGS's prediction while it lasts, then to ours.
5. **Calibrate against the satellites' own broadcasts** over the measured day: the antenna offset
   (IGS orbits refer to the centre of mass) and the clock reference. The clock is a straight line
   through the IGS's clock estimates.

Three faults found on the way are recorded in the [results](results/2026-10-18-backtest.md): an
Earth-orientation mix of two IAU conventions (0.49° of error), noisy adaptive integration, and an
unconverged fit. Fixing them cut the five-day error from 10.6 m to 4.7 m (median) and from 66 m to
25 m (worst).

## Output

`python3 -m extended_ephemeris --out out/` fetches the newest IGS data and writes:

| file | what |
| --- | --- |
| `gps-7d.rnx` | the predicted records as **RINEX 3.04** navigation data: 32 satellites, one every 2 h, 7 days |
| `gps-7d-casic.bin` | the same as **CASIC MSG-GPSEPH** frames (223 kB) for AT6558 receivers: Quectel L76K, ATGM336H and others |
| `gps-7d.json` | when and from what it was made, fit quality, and a self-check against the satellites' current broadcasts (the run fails if that check is off) |
| `fit-state.json` | each satellite's fitted state; `--warm fit-state.json` starts the next run from it |

A receiver wants, for each satellite, the record whose reference time (toe) is nearest now. For
the L76K, [l76k-at6558-agnss](https://github.com/badriram/l76k-at6558-agnss) has a C function,
`casic_eph_select()`, that picks them from `gps-7d-casic.bin`, and the measurements showing that
the receiver accepts injected orbit data and cold-starts in ~6 s with it. Other receivers' formats
(u-blox UBX-MGA-GPS-EPH and so on) can be written from the same records.

Other inputs: `--sp3 FILE --erp FILE --rinex FILE...` to predict from given files, `--days N`.
Needs Python 3.9+ with numpy, scipy and pyerfa (`pip install -r requirements.txt`). A run takes
about 6 minutes on 14 cores.

## Publishing

`.github/workflows/publish.yml` publishes the four files to the `ephemeris-data` branch. It is
manual for now: a run takes far longer on GitHub's 4-core runners than here, and the warm start
does not yet save time. Once it is faster it can run after each IGS issue, four times a day.

## Limits

- GPS only: the IGS ultra-rapid orbit has no BeiDou or Galileo. (MGEX products do; not yet used.)
- One backtest, 4.7 days. Days 5-7 are extrapolated.
- Satellites that manoeuvre during the week will be wrong until the next run; the receiver still
  recovers by downloading their real orbit data.
- No receiver test with predicted records yet.

## Credits and sources

- [IGS](https://igs.org/) ultra-rapid orbits and Earth rotation parameters, and merged broadcast
  ephemeris, via [BKG](https://igs.bkg.bund.de/). NASA CDDIS archives the same products.
- [ERFA](https://github.com/liberfa/erfa) (pyerfa), the open copy of the IAU's SOFA routines.
- Montenbruck & Gill, *Satellite Orbits*, for the gravity-field recursion; EGM2008 coefficients;
  the ECOM radiation-pressure model of Beutler, Springer and others; IS-GPS-200 for the broadcast
  ephemeris.

## Licence

MIT. See [LICENSE](LICENSE). The published data carries the same licence; the IGS asks that its
products be credited.
