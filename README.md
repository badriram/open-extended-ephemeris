# open-extended-ephemeris

**Offline rapid cold start for GPS and Galileo.** A GPS receiver with no network spends 30-40 s on a cold start
downloading each satellite's orbit data at 50 bit/s. Give it that data and it fixes in a few
seconds, but broadcast orbit data expires after about four hours. This project **predicts it a
week ahead**: download one file while you have a network, and cold starts stay fast for days
without one.

This is what the industry calls **extended ephemeris**, sold as u-blox AssistNow Offline,
MediaTek EPO, Broadcom LTO, Qualcomm XTRA and SiRF InstantFix. All of those are closed services.
As far as we could find, this is the first open one: the method, the code and the measurements.

## How good it is

Two backtests: predicted from 2026-10-13 06:00 and 2026-10-14 00:00 UTC with the data that was
current then, then graded against what the satellites actually broadcast over the following days
([results](results/2026-10-18-code-pipeline.md), `extended-ephemeris grade`). Range error is the
satellite's contribution to the error of a distance measurement; a position is off by roughly
twice it. Median / 95th percentile, both backtests:

| days ahead | GPS | Galileo |
| --- | --- | --- |
| 0-0.5 | 0.2-0.3 / 0.7-1.4 m | 0.1 / 0.5 m |
| 1-1.5 | 0.3-0.5 / 2.9-3.3 m | 0.2 / 1.6-2.4 m |
| 2-2.5 | 0.5-0.8 / 5.0-5.4 m | 0.3-0.4 / 2.3-5.9 m |
| 3-3.5 | 0.8-1.3 / 6.2-8.4 m | 0.5-0.7 / 4.4-7.8 m |
| 4-4.5 | 1.3-1.5 / 8.3-10.8 m | 0.7-1.1 / 5.8-6.0 m |
| 4.5-5 | 1.8 / 12.2 m | 1.5 / 6.9 m |

So after five days offline, a first fix is typically a few metres out, and the receiver replaces
the prediction with the satellites' own data once it has tracked them for about 30 s. Every
record also states its own accuracy (URA for GPS, SISA for Galileo), which the range error
exceeded in at most about 4% of samples, so a receiver can weight an aged or poorly behaved satellite
down. These figures are at least as good as every published consumer-grade extended-ephemeris
result we found. The clock model alone was also tested over the whole of 2025 (every 6 hours,
[results](results/2025-clock-backtest.md)): GPS day 4-5 clock error median 2.0 m, 95th
percentile 13 m, with no seasonal variation. Days 5-7 have not been graded yet, and no receiver
has yet been given predicted records (they are ordinary records with a current reference time,
which receivers accept from other sources).

## How it works

Per satellite, with plain physics and least squares (no machine learning):

1. **Start from a precise orbit prediction.** By default that is the 5-day GPS and Galileo
   prediction of the Center for Orbit Determination in Europe (CODE, University of Bern), served
   openly from its archive. The fallback is the International GNSS Service's ultra-rapid orbit
   (24 h measured, 24 h predicted, GPS only), which BKG serves anonymously.
2. **Fit a force model** to that arc: Earth's gravity field (EGM2008 to degree and order 4), the
   Sun and Moon, relativity, and solar radiation pressure as the nine ECOM2 terms with Earth's
   shadow and a smooth penumbra; Earth orientation from the IAU 2006/2000A model (ERFA) with the
   source's polar motion and UT1, continued by IERS Bulletin A's predictions. Fitted by
   Levenberg-Marquardt, with a priori limits on the minor radiation terms, to 0.1-0.3 m RMS.
3. **Integrate a week ahead** (8th-order Dormand-Prince, fixed 100 s steps, so the result is a
   smooth function of the inputs).
4. **Fit an ordinary broadcast ephemeris** (IS-GPS-200 LNAV; Galileo I/NAV uses the same
   Keplerian form) to every 2 h of the trajectory: to the source's prediction while it lasts,
   then to ours.
5. **Calibrate against the satellites' own broadcasts** over the known part of the arc: the
   antenna offset (precise orbits refer to the centre of mass).
6. **Predict each clock from the satellite's own broadcast clocks over the last five days**, as a
   quadratic: GPS clocks age, so their rate drifts. This was the biggest single gain: clocks were
   98% of the remaining error, and the quadratic cut the five-day median from 4.7 m to 2.2 m. The
   clock history travels in `fit-state.txt`, so a run downloads only the newest broadcasts.
7. **State each record's accuracy**: the population figure for its age, or more for a satellite
   whose clock failed to predict its own last 1.5 days.

Galileo's three remaining IOV satellites (E11, E12, E19) are left out: their clocks are ten times
less predictable than the rest ([results](results/2025-clock-backtest.md)).

What changed the numbers, and what did not, is recorded in the results, starting with
[the IGS-based backtest](results/2026-10-18-backtest.md):
fixing an Earth-orientation mix of two IAU conventions (0.49° of error), noisy adaptive
integration and an unconverged fit took the five-day median from 10.6 m to 4.7 m; the clock model
then took it to 2.2 m. The solid Earth tide and Bulletin A are physically right but changed the
results by less than the backtests' noise: Earth rotation is not what limits the orbits (even the
UT1 measured afterwards did no better), the force model is. Of the force model, ECOM2's extra
radiation-pressure terms cut the five-day orbit error by 6-11%; nothing else tried (relativity,
gravity to degree 6, the IGS's GM, a finer step, different priors, a Galileo box-wing model)
moved it. CODE's 5-day prediction then cut the five-day GPS orbit error from 5.0 m to 1.7 m and
made Galileo possible at all ([results](results/2026-10-18-code-prediction.md)); for range error
the clocks still decide how far out a first fix is.

## Running it

A single Rust program with no crate dependencies; ERFA's C routines are compiled in from
`vendor/erfa`. Needs Rust, a C compiler, and `curl` and `gzip` on the path.

```bash
cargo build --release
./target/release/extended-ephemeris --out out/
```

That fetches CODE's newest 5-day prediction (falling back to the IGS ultra-rapid orbit) and the
broadcasts it needs from BKG, and writes:

| file | what |
| --- | --- |
| `gps-7d.rnx` | the predicted records as **RINEX 3.04** navigation data: ~31 GPS and ~25 Galileo satellites, one every 2 h, 7 days |
| `gps-7d-ubx.bin` | the same as **u-blox UBX-MGA-GPS-EPH and UBX-MGA-GAL-EPH** frames (~370 kB) for u-blox M8, M9 and M10 receivers |
| `gps-7d-casic.bin` | the GPS records as **CASIC MSG-GPSEPH** frames (~215 kB) for AT6558 receivers: Quectel L76K, ATGM336H and others |
| `gps-7d.json` | when and from what it was made, fit quality, and a self-check against the satellites' current broadcasts (the run fails if that check is off) |
| `fit-state.txt` | each satellite's fitted state and recent broadcast clocks (~100 kB); `--warm fit-state.txt` starts the next run from it, and the run then needs only the newest broadcasts |

**A whole run takes 17-35 seconds** on 14 cores, downloads aside: 56 satellites, each fitted to
up to six days of 5-minute positions. From the IGS orbit (32 satellites, 15-minute positions) it
takes about 3 seconds. The first version, in Python, took 5.6 minutes for that.

The file names keep `gps-` so that existing links keep working. Other options: `--source
auto|code|igs`, `--systems G` (GPS only), `--sp3 FILE|URL[,FILE|URL...] --erp FILE|URL` with
`--prediction-start TIME` for an orbit that is all prediction, and `--rinex FILE|URL ...` to predict from given
files, `--eop FILE|URL|none` (Bulletin A; default the latest from the IERS), `--days N`, `--jobs N`,
and for experiments `--clock line|harmonic|harmonic-drift|broadcast-line:DAYS|broadcast-quad:DAYS`,
`--srp-sigma S1,...,S9`, `--no-tides`, `--no-relativity`, `--box-wing`
(a priori Galileo radiation pressure from GSC metadata), `--no-antenna-offset`, and
`--from-broadcast TIME [--arc-hours H]`
(predict from broadcast orbits, no precise orbit; needs `--erp` and `--rinex`). `extended-ephemeris grade PREDICTED.rnx BROADCAST.rnx.gz ...`
grades a prediction against what the satellites broadcast (`--system G|E` for one
constellation; it also reports how often the stated accuracy was exceeded), and `extended-ephemeris clock-backtest
BROADCAST.rnx.gz ...` tests clock models over any span of broadcast history (BKG keeps a week of daily files at
`https://igs.bkg.bund.de/root_ftp/IGS/BRDC/<year>/<day>/BRDC00WRD_S_<year><day>0000_01D_MN.rnx.gz`).

A receiver wants the time, ideally a rough position, and for each satellite the record whose
reference time (toe) is nearest now. `select` writes exactly that, ready to send down a serial
port:

```bash
# from the published week, for now (or --at 2026-10-23T12:00Z), with an optional rough position
./target/release/extended-ephemeris select --format ubx --out assist.ubx --position LAT,LON
./target/release/extended-ephemeris select --format casic --out assist.bin --position LAT,LON
```

For u-blox that is UBX-MGA-INI-TIME_UTC, UBX-MGA-INI-POS_LLH if a position was given, then one
UBX-MGA-GPS-EPH or UBX-MGA-GAL-EPH per satellite (about 4.5 kB); the receiver needs the time
first, since those messages carry no week number. For CASIC it is AID-INI then MSG-GPSEPH. For the L76K, [l76k-at6558-agnss](https://github.com/badriram/l76k-at6558-agnss) has a C function,
`casic_eph_select()`, that picks them from `gps-7d-casic.bin`, and the measurements showing that
the receiver accepts injected orbit data and cold-starts in ~6 s with it. The u-blox messages
follow u-blox's interface descriptions, checked field by field against
[pyubx2](https://github.com/semuconsulting/pyubx2); no u-blox receiver has been given them yet.

## Publishing

`.github/workflows/publish.yml` publishes the four files to the `ephemeris-data` branch four times
a day (04:30, 10:30, 16:30 and 22:30 UTC), each with the newest orbit and broadcasts. Each run
starts from the previous one's `fit-state.txt`, so it downloads only the newest broadcasts.
The branch holds only the latest files:

- `https://raw.githubusercontent.com/badriram/open-extended-ephemeris/ephemeris-data/gps-7d.rnx`
- `https://raw.githubusercontent.com/badriram/open-extended-ephemeris/ephemeris-data/gps-7d-ubx.bin`
- `https://raw.githubusercontent.com/badriram/open-extended-ephemeris/ephemeris-data/gps-7d-casic.bin`
- `https://raw.githubusercontent.com/badriram/open-extended-ephemeris/ephemeris-data/gps-7d.json`

## Limits

- GPS and Galileo. BeiDou and GLONASS are not predicted; CODE's file has GLONASS, but GLONASS
  broadcasts a different kind of ephemeris (positions and velocities, not Keplerian elements).
  CASIC output is GPS only.
- Two backtests, of 4.5 and 5 days. Days 5-7 are extrapolated and not yet graded.
- Satellites that manoeuvre or whose clocks jump during the week will be wrong until the next
  run; the receiver still recovers by downloading their real orbit data.
- No receiver test with predicted records yet.

## Credits and sources

- [IERS](https://www.iers.org/) Bulletin A Earth orientation predictions; IERS Conventions 2010
  for the solid Earth tide.
- CODE, the Center for Orbit
  Determination in Europe at the Astronomical Institute of the University of Bern: the 5-day GPS
  and Galileo orbit prediction and its Earth rotation parameters (product DOI
  [10.48620/97774](https://doi.org/10.48620/97774)).
- [IGS](https://igs.org/) ultra-rapid orbits and Earth rotation parameters, and merged broadcast
  ephemeris, via [BKG](https://igs.bkg.bund.de/). NASA CDDIS archives the same products.
- The European GNSS Service Centre's Galileo satellite metadata (geometry, optical properties,
  attitude, mass) for the box-wing model.
- [ERFA](https://github.com/liberfa/erfa) 2.0.1, the open copy of the IAU's SOFA routines; the 37
  files used are vendored unmodified in `vendor/erfa` under its BSD-3-Clause licence.
- The DOP853 coefficients of Hairer, Nørsett and Wanner, as tabulated in SciPy.
- Montenbruck & Gill, *Satellite Orbits*, for the gravity-field recursion; EGM2008 coefficients;
  the ECOM radiation-pressure model of Beutler, Springer and others; IS-GPS-200 for the broadcast
  ephemeris.

## Licence

MIT. See [LICENSE](LICENSE), except `vendor/erfa` (BSD-3-Clause, see its LICENSE). The published
data carries the MIT licence; CODE and the IGS ask that their products be credited.
