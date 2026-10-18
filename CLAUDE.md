# Working in this repo

**open-extended-ephemeris** predicts GPS broadcast ephemeris a week ahead (extended ephemeris, the
open counterpart to u-blox AssistNow Offline / MediaTek EPO), so a receiver cold-starts in seconds
for days without a network. One Rust program, `extended-ephemeris`; MIT. README.md has the method
and the headline numbers; `results/` has every experiment, including the ones that did nothing.
Keep both current when a result changes.

## Structure

- `src/predict.rs` — the pipeline: IGS ultra-rapid orbit + ERP + broadcasts in; per satellite, fit
  the force model, propagate, fit an LNAV record every 2 h, calibrate, predict the clock; write
  `gps-7d.rnx`, `gps-7d-ubx.bin`, `gps-7d-casic.bin`, `gps-7d.json`, `fit-state.txt`.
- `src/orbit.rs` — force model (gravity + solid tide, Sun, Moon, relativity, ECOM2 radiation
  pressure with a smooth penumbra), fixed-step DOP853, Levenberg-Marquardt fit with priors.
- `src/eo.rs` — Earth orientation (ERFA, IGS ERP, IERS Bulletin A), Sun/Moon tables.
- `src/lnav.rs` — IS-GPS-200 broadcast ephemeris: evaluate and fit.
- `src/rinex.rs`, `src/sp3.rs`, `src/casic.rs`, `src/ubx.rs` — file formats. `src/fetch.rs` — curl/gzip.
- `src/select.rs` — `select`: the ready-to-send burst (time, position, current frames) per receiver.
- `src/grade.rs`, `src/clockbt.rs` — `grade` and `clock-backtest` subcommands.
- `vendor/erfa/` — 37 ERFA 2.0.1 files, **unmodified**, BSD-3-Clause, compiled by `build.rs`.
- `.github/workflows/publish.yml` — runs four times a day and force-pushes the outputs to the
  `ephemeris-data` branch, which holds only the latest files. Consumers fetch from there.

## Conventions

- **No crates.** Downloads go through the system `curl` and `gzip`; linear algebra, solvers and
  JSON are hand-written. Keep it that way unless there is a strong reason.
- **ERFA is vendored, not ported.** Never edit `vendor/erfa`; to need another routine, copy the
  file (and its dependencies) from the same release. Coefficient tables are never retyped by hand.
- **Measure every model change against both backtests** before keeping it, and record the numbers
  in `results/`, including changes that did nothing. Choosing settings on one backtest and
  confirming on the other is how the 5-day clock window was picked.
- **No look-ahead.** A run may use broadcasts only up to `pred_start + LATENCY` (4 h), as it would
  in practice. A backtest that lets later broadcasts in is not a backtest.
- **The metric is signal-in-space range error** (radial and clock, along/cross-track down-weighted
  by 1/7), not 3D orbit error. A change that only improves along-track barely helps a fix.

## Lessons that cost time

- Earth orientation once mixed the equinox-based `pnm06a` with the CIO-based Earth rotation angle:
  0.49° of error, a 2.6 m fit floor and a fake 80 nm/s² radiation term. Use `c2i06a` with the ERA.
- Adaptive integration made outputs jitter with their inputs (a 0.01 nm/s² change "moved" the
  orbit 0.5 m), so no fit could converge. Integration is fixed-step; keep it so.
- The fit needs scaled unknowns (velocity partials are ~1e5 × position ones), damping, and priors
  on the secondary radiation terms, or it stops wherever it happens to be.
- Clocks, not orbits, dominate range error. A quadratic through 5 days of the satellite's own
  broadcast clocks beat every line and harmonic model. The orbit error left is along-track.
- Bulletin A, the solid tide, relativity, degree-6 gravity and a finer step changed the results by
  less than noise. Even the UT1 measured afterwards did no better: Earth rotation is not the limit.

## Running

```bash
cargo test --release
cargo build --release
./target/release/extended-ephemeris --out out/          # newest IGS data from BKG
./target/release/extended-ephemeris grade out/gps-7d.rnx BRDC...rnx.gz ...
./target/release/extended-ephemeris clock-backtest --models quad:5,line:1 BRDC...rnx.gz ...
```

A run takes ~3 s of computing on 14 cores; downloads dominate. BKG keeps about a week of IGS
ultra-rapid products under `root_ftp/IGS/products/<gps week>/`; daily broadcast files are under
`root_ftp/IGS/BRDC/<year>/<doy>/` for the recent past and for past years (older years have also
lived under `root_ftp/IGS/obs/<year>/<doy>/`). Backtests need the IGS orbit of the start day plus
broadcasts from five days before it to the end of the graded span.

## Related

`github.com/badriram/l76k-at6558-agnss` is the receiver side for the Quectel L76K (C library with
`casic_eph_select()`, measurements). It still carries its own CASIC/IGS encoders; one day it may
consume this repo's output instead. Not yet.
