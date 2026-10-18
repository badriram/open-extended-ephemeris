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

Two backtests: predicted from the IGS ultra-rapid orbits of 2026-10-13 06:00 and 2026-10-14 00:00
UTC, then graded against what the satellites actually broadcast over the following days
([results](results/2026-10-18-backtest.md), `extended-ephemeris grade`). Range error is the
satellite's contribution to the error of a distance measurement; a position is off by roughly
twice it.

| days after the IGS orbit | range error, median | 95th percentile | worst |
| --- | --- | --- | --- |
| 0-1 (the IGS's own prediction) | 0.2-0.3 m | 0.7-1.9 m | 3.1 m |
| 1-2 | 0.3-0.6 m | 2.3-4.3 m | 5.9 m |
| 2-3 | 0.6-1.2 m | 4.6-7.0 m | 8.8 m |
| 3-4 | 1.0-1.6 m | 6.1-9.8 m | 11 m |
| 4-5 | 1.3-2.2 m | 7.5-12.3 m | 14.6 m |

So after five days offline, a first fix is typically a few metres out, and the receiver replaces
the prediction with the satellites' own data once it has tracked them for about 30 s. These
figures are at least as good as every published consumer-grade extended-ephemeris result we
found. The clock model alone was also tested over the whole of 2025 (every 6 hours, 32 satellites,
[results](results/2025-clock-backtest.md)): day 4-5 clock error median 2.0 m, 95th percentile 13 m,
with no seasonal variation. Days 5-7 have not been graded yet, and no receiver has yet been given predicted records
(they are ordinary records with a current reference time, which receivers accept from other
sources).

## How it works

Per GPS satellite, with plain physics and least squares (no machine learning):

1. **Start from the IGS.** The International GNSS Service's ultra-rapid orbit gives 24 h of
   measured positions and clocks and 24 h of its own prediction; BKG serves it anonymously.
2. **Fit a force model** to those 48 h: Earth's gravity field (EGM2008 to degree and order 4), the
   Sun and Moon, relativity, and solar radiation pressure as the nine ECOM2 terms with Earth's
   shadow and a smooth penumbra; Earth orientation from the IAU 2006/2000A model (ERFA) with the IGS's polar
   motion and UT1, continued by IERS Bulletin A's predictions. Fitted by Levenberg-Marquardt, with a priori limits on the minor radiation
   terms, to 0.1-0.26 m RMS.
3. **Integrate a week ahead** (8th-order Dormand-Prince, fixed 100 s steps, so the result is a
   smooth function of the inputs).
4. **Fit an ordinary broadcast ephemeris** (IS-GPS-200 LNAV) to every 2 h of the trajectory: to the
   IGS's prediction while it lasts, then to ours.
5. **Calibrate against the satellites' own broadcasts** over the measured day: the antenna offset
   (IGS orbits refer to the centre of mass).
6. **Predict each clock from the satellite's own broadcast clocks over the last five days**, as a
   quadratic: GPS clocks age, so their rate drifts. This was the biggest single gain: clocks were
   98% of the remaining error, and the quadratic cut the five-day median from 4.7 m to 2.2 m. The
   clock history travels in `fit-state.txt`, so a run downloads only the newest broadcasts.

What changed the numbers, and what did not, is recorded in the [results](results/2026-10-18-backtest.md):
fixing an Earth-orientation mix of two IAU conventions (0.49° of error), noisy adaptive
integration and an unconverged fit took the five-day median from 10.6 m to 4.7 m; the clock model
then took it to 2.2 m. The solid Earth tide and Bulletin A are physically right but changed the
results by less than the backtests' noise: Earth rotation is not what limits the orbits (even the
UT1 measured afterwards did no better), the force model is. Of the force model, ECOM2's extra
radiation-pressure terms cut the five-day orbit error by 6-11%; nothing else tried (relativity,
gravity to degree 6, the IGS's GM, a finer step, different priors) moved it. The orbit error left
is almost all along-track, which range error weights about seven times down, so the clocks still
decide how far out a first fix is.

## Running it

A single Rust program with no crate dependencies; ERFA's C routines are compiled in from
`vendor/erfa`. Needs Rust, a C compiler, and `curl` and `gzip` on the path.

```bash
cargo build --release
./target/release/extended-ephemeris --out out/
```

That fetches the newest IGS ultra-rapid orbit and the broadcasts it needs from BKG, and writes:

| file | what |
| --- | --- |
| `gps-7d.rnx` | the predicted records as **RINEX 3.04** navigation data: 32 satellites, one every 2 h, 7 days |
| `gps-7d-ubx.bin` | the same as **u-blox UBX-MGA-GPS-EPH** frames (207 kB) for u-blox M8, M9 and M10 receivers |
| `gps-7d-casic.bin` | the same as **CASIC MSG-GPSEPH** frames (223 kB) for AT6558 receivers: Quectel L76K, ATGM336H and others |
| `gps-7d.json` | when and from what it was made, fit quality, and a self-check against the satellites' current broadcasts (the run fails if that check is off) |
| `fit-state.txt` | each satellite's fitted state and recent broadcast clocks (~100 kB); `--warm fit-state.txt` starts the next run from it, and the run then needs only the newest broadcasts |

**A whole run takes about 2 seconds** (32 satellites on 14 cores, 22 CPU-seconds; downloads
aside). The first version, in Python, took 5.6 minutes for identical results.

Other options: `--sp3 FILE|URL --erp FILE|URL` and `--rinex FILE|URL ...` to predict from given
files, `--eop FILE|URL|none` (Bulletin A; default the latest from the IERS), `--days N`, `--jobs N`,
and for experiments `--clock line|harmonic|harmonic-drift|broadcast-line:DAYS|broadcast-quad:DAYS`,
`--srp-sigma S1,...,S9`, `--no-tides` and `--no-relativity`. `extended-ephemeris grade PREDICTED.rnx BROADCAST.rnx.gz ...`
grades a prediction against what the satellites broadcast, and `extended-ephemeris clock-backtest
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
UBX-MGA-GPS-EPH per satellite (about 2.5 kB); the receiver needs the time first, since MGA-GPS-EPH
carries no week number. For CASIC it is AID-INI then MSG-GPSEPH. For the L76K, [l76k-at6558-agnss](https://github.com/badriram/l76k-at6558-agnss) has a C function,
`casic_eph_select()`, that picks them from `gps-7d-casic.bin`, and the measurements showing that
the receiver accepts injected orbit data and cold-starts in ~6 s with it. The u-blox messages
follow u-blox's interface descriptions, checked field by field against
[pyubx2](https://github.com/semuconsulting/pyubx2); no u-blox receiver has been given them yet.

## Publishing

`.github/workflows/publish.yml` publishes the four files to the `ephemeris-data` branch four times
a day, about 4.5 hours after each IGS ultra-rapid issue (04:30, 10:30, 16:30 and 22:30 UTC). Each
run starts from the previous one's `fit-state.txt`, so it downloads only the newest broadcasts.
The branch holds only the latest files:

- `https://raw.githubusercontent.com/badriram/open-extended-ephemeris/ephemeris-data/gps-7d.rnx`
- `https://raw.githubusercontent.com/badriram/open-extended-ephemeris/ephemeris-data/gps-7d-ubx.bin`
- `https://raw.githubusercontent.com/badriram/open-extended-ephemeris/ephemeris-data/gps-7d-casic.bin`
- `https://raw.githubusercontent.com/badriram/open-extended-ephemeris/ephemeris-data/gps-7d.json`

## Limits

- GPS only: the IGS ultra-rapid orbit has no BeiDou or Galileo. (MGEX products do; not yet used.)
- One backtest, 4.7 days. Days 5-7 are extrapolated.
- Satellites that manoeuvre during the week will be wrong until the next run; the receiver still
  recovers by downloading their real orbit data.
- No receiver test with predicted records yet.

## Credits and sources

- [IERS](https://www.iers.org/) Bulletin A Earth orientation predictions; IERS Conventions 2010
  for the solid Earth tide.
- [IGS](https://igs.org/) ultra-rapid orbits and Earth rotation parameters, and merged broadcast
  ephemeris, via [BKG](https://igs.bkg.bund.de/). NASA CDDIS archives the same products.
- [ERFA](https://github.com/liberfa/erfa) 2.0.1, the open copy of the IAU's SOFA routines; the 37
  files used are vendored unmodified in `vendor/erfa` under its BSD-3-Clause licence.
- The DOP853 coefficients of Hairer, Nørsett and Wanner, as tabulated in SciPy.
- Montenbruck & Gill, *Satellite Orbits*, for the gravity-field recursion; EGM2008 coefficients;
  the ECOM radiation-pressure model of Beutler, Springer and others; IS-GPS-200 for the broadcast
  ephemeris.

## Licence

MIT. See [LICENSE](LICENSE), except `vendor/erfa` (BSD-3-Clause, see its LICENSE). The published
data carries the MIT licence; the IGS asks that its products be credited.
