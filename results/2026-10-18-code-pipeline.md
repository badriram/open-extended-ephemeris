# CODE's 5-day orbit as the default source, Galileo published

Before this change the published files were GPS only, built on the IGS ultra-rapid orbit (24 h
measured plus 24 h predicted), with our force model doing all the work from day 1 to day 7.
[The CODE comparison](2026-10-18-code-prediction.md) showed that CODE's 5-day prediction
(COD0OPSPRD) beats our extension for GPS and is the only good orbit source for Galileo. So:

- **Orbits:** CODE's 5-day GPS and Galileo prediction for days 0-5, our force model fitted to it
  for days 5-7. The IGS ultra-rapid orbit is the fallback (`--source igs`, or automatically if
  CODE cannot be reached), and the fallback is GPS only.
- **Clocks:** unchanged: a quadratic through five days of each satellite's own broadcast clocks.
- **Galileo:** the three IOV satellites still transmitting (E11, E12, E19) are left out. Their
  clocks were 11-35 m median at day 4-5 over 2025, against 1.2 m for the rest, and E11 and E12
  were the tail of the CODE comparison.
- **Accuracy per satellite:** each record's URA (GPS) or SISA (Galileo) is the larger of the
  population figure for its age (1.6 m on day 0 up to 18 m by day 6) and the satellite's own clock
  error. That error comes from a hold-out test in every run: fit the clock model to the 3.5 days
  of broadcasts before the last 1.5, predict the last 1.5, take the median error, then double it
  and grow it with age as `(age / 0.75 d)^1.3`.
- **u-blox:** `gps-7d-ubx.bin` now carries UBX-MGA-GAL-EPH for Galileo next to UBX-MGA-GPS-EPH.
  CASIC stays GPS only: the base AT6558 does not track Galileo, and CASIC has no documented
  Galileo ephemeris message for variants that do.
  The file names keep `gps-` so that existing URLs keep working.

## Timing of CODE's files

The file named for day N starts at 00:00 on day N+1. On 2026-10-18 the one for day 291 was
uploaded at 13:08 UTC, and the one for day 290 at 06:31. So for part of each day the newest file
starts in the future and covers none of the present. The run then joins it to the one before:
the older file up to the newer file's first epoch, then the newer one. Calibration against
broadcasts uses the joined arc up to the prediction start plus four hours.

Both backtests use the file that would have been newest at their start time. A file uploaded
the afternoon before is in use at 06:00 the next morning, so in backtest 1 CODE's prediction was
already about 17 h old at "day 0". That is how a live run sees it too.

## Backtests

The same two starts as before. Range error is median / 95th percentile / worst, in metres.
"Over" is the share of samples whose range error exceeded the accuracy the record states.

Backtest 1, predicting from 2026-10-13 06:00 (CODE file of day 285):

| days | GPS | over | Galileo | over |
| --- | --- | --- | --- | --- |
| 0-0.5 | 0.26 / 0.68 / 1.4 | 0.0% | 0.12 / 0.48 / 1.3 | 0.0% |
| 1.0-1.5 | 0.46 / 3.29 / 4.3 | 0.6% | 0.22 / 1.55 / 2.7 | 0.0% |
| 2.0-2.5 | 0.82 / 5.39 / 7.3 | 0.7% | 0.40 / 2.34 / 7.1 | 0.0% |
| 3.0-3.5 | 1.32 / 8.42 / 9.7 | 0.3% | 0.71 / 4.35 / 8.7 | 0.0% |
| 4.0-4.5 | 1.50 / 10.75 / 12.1 | 0.3% | 1.13 / 6.01 / 8.3 | 0.0% |
| 4.5-5.0 | 1.81 / 12.16 / 15.2 | 1.1% | 1.47 / 6.94 / 7.7 | 0.0% |

Backtest 2, predicting from 2026-10-14 00:00 (CODE file of day 286):

| days | GPS | over | Galileo | over |
| --- | --- | --- | --- | --- |
| 0-0.5 | 0.23 / 1.38 / 2.2 | 3.3% | 0.12 / 0.47 / 1.8 | 0.5% |
| 1.0-1.5 | 0.34 / 2.90 / 4.1 | 0.0% | 0.20 / 2.37 / 6.5 | 4.0% |
| 2.0-2.5 | 0.54 / 5.00 / 5.7 | 0.0% | 0.34 / 5.94 / 8.5 | 1.7% |
| 3.0-3.5 | 0.83 / 6.24 / 6.8 | 0.0% | 0.51 / 7.82 / 13.2 | 0.0% |
| 4.0-4.5 | 1.28 / 8.29 / 10.2 | 0.0% | 0.68 / 5.83 / 17.3 | 0.0% |

Against the IGS-based pipeline it replaces, at day 4.5-5 of backtest 1:

- **GPS:** 1.81 / 12.2 m, against 2.19 / 12.3 m. The median orbit error falls from 4.98 m to
  1.65 m. Range error barely moves, because clocks decide it.
- **Galileo:** 1.47 / 6.9 / 7.7 m. Built from broadcast orbits it was 9.3 / 23.4 m; with CODE but
  the IOV satellites still in, 1.68 / 12.7 / 19.3 m.

Backtest 2's Galileo tail is E09, at 17 m. Its clock has jumped in the past (see the 2025 clock
backtest), and this is the same kind of event, which no smooth model predicts.

## Does the per-satellite accuracy help?

Both backtests were rerun with only the age table. The share of GPS samples over the stated
accuracy, per half-day bin:

- **Table only:** up to 5.6% in backtest 1 and 5.9% in backtest 2.
- **With the hold-out:** at most 3.2% and 3.3%.

G27, the worst GPS satellite in both backtests, is flagged: in the live run of 2026-10-18 it
states 23 m at day 4, against 12 m from the table. G07 and G12 are raised too. G29 and G03 are almost as bad and are not flagged. Their
error is not a steady clock drift the last 1.5 days would show. A modest gain, kept because it
is free and conservative.

The stated accuracy acts as roughly a 95-99% bound, not the one-sigma value IS-GPS-200 defines
for URA. That is deliberate: a receiver that over-trusts a prediction is worse off than one that
under-trusts it.

## Cost

A run now takes 17-35 s instead of 2-3 s. CODE's file has 5-minute samples over up to six days,
and there are 56 satellites instead of 32. `fit-state.txt` grows to ~750 kB, mostly Galileo's
10-minute broadcast clocks.
