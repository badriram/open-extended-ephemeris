# Predicting from broadcast orbits alone

The IGS ultra-rapid orbit the predictor starts from is GPS only, and its multi-GNSS equivalents
are not anonymously available (BKG's MGEX folder stopped in 2022; NASA CDDIS needs a login). The
daily broadcast files carry every constellation, so: how much is lost if the 48 h track is built
from the satellites' own broadcast orbits instead? `--from-broadcast TIME` does that, using only
records a receiver could have had by TIME (reference time at most 2 h after it), and starts the
prediction there. Everything else is unchanged; graded as before against later broadcasts.

**An antenna-offset trap.** Broadcast orbits describe the antenna, 1-2.5 m from the centre of
mass, and no orbit can follow a point offset radially, so the fit settles on the centre of mass
(fit RMS 0.87 m against 0.15 m from IGS orbits). The first version then measured the offset
against the track, which is the broadcast itself, found zero, and every prediction sat +1.37 m
radially out from the first record (range error 1.43 m at day 0). Measuring the offset against
the fitted orbit instead fixed it (0.38 m).

Range error median / 95th percentile / worst, and orbit 3D median:

| days | backtest 1, IGS orbit | backtest 1, broadcast only | backtest 2, IGS orbit | backtest 2, broadcast only |
| --- | --- | --- | --- | --- |
| 0.5-1.0 | 0.34 / 1.56 / 3.1 m (0.7 m) | 0.55 / 1.80 / 2.6 m (1.9 m) | 0.24 / 1.94 / 3.1 m (0.8 m) | 0.44 / 2.16 / 3.4 m (1.8 m) |
| 2.5-3.0 | 1.22 / 6.94 / 8.8 m (1.6 m) | 2.06 / 6.63 / 9.9 m (8.0 m) | 0.74 / 5.20 / 6.6 m (1.5 m) | 1.40 / 5.99 / 8.0 m (5.8 m) |
| 4.0-4.5 | 1.84 / 10.9 / 12.0 m (3.8 m) | 3.64 / 12.2 / 14.7 m (16.6 m) | 1.26 / 7.49 / 9.9 m (3.5 m) | 2.06 / 11.8 / 13.1 m (10.3 m) |
| 4.5-5.0 | 2.19 / 12.3 / 14.6 m (5.0 m) | 4.28 / 14.1 / 16.2 m (19.2 m) | | |

Broadcast-only roughly doubles the median range error and adds 15-55% to the 95th percentile;
orbits drift three to four times faster, since the arc has no IGS prediction half and broadcast
orbits are noisier. Still well within what a first fix needs, so it is a workable source for the
constellations IGS ultra-rapid orbits lack.

Longer broadcast arcs trade one error for another (backtest 1, 4.5-5 days):

| arc | fit RMS | range error | orbit 3D |
| --- | --- | --- | --- |
| 48 h | 0.87 m | 4.28 / 14.1 / 16.2 m | 19.2 m |
| 72 h | 0.97 m | 4.10 / 16.9 / 20.8 m | 15.3 m |
| 96 h | 1.11 m | 4.29 / 21.6 / 26.5 m | 13.9 m |

A 120 h arc produced one satellite whose clock rate would not encode, and the run aborted: now a
record that will not encode is dropped and counted (`dropped_records` in gps-7d.json) instead.
