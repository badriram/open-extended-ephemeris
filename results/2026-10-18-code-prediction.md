# Starting from CODE's 5-day prediction

CODE (Center for Orbit Determination in Europe, University of Bern) publishes, anonymously at
`https://zhw-b.s3.cloud.switch.ch/aiub/CODE/` (the old ftp.aiub.unibe.ch no longer answers):

- `COD0OPSULT.SP3`: ultra-rapid orbits, 24 h measured + 24 h predicted, GPS, GLONASS and Galileo,
  refreshed through the day (DOI 10.48620/97774);
- `COD0OPSPRD_<yyyy><ddd>0000_05D_05M_ORB.SP3`: a 0-5 day orbit prediction every 5 minutes from
  those ultra-rapid orbits, issued daily about 06:40 UTC for the following day, with a 21-day ERP.

Both carry Galileo, which the IGS ultra-rapid orbit does not. The predictor reads Galileo from SP3
now, and `--prediction-start TIME` treats samples after TIME as predicted (CODE's 5-day file flags
its samples predicted from minutes after its first epoch). Two fixes it needed: calibration falls
back to the samples the run's broadcasts cover when an orbit is all prediction, and no record is
fitted whose window (toe ± 2 h) reaches before the orbit's first epoch (the extrapolation gave
errors of 40-150 km in the self-check, which would have failed a publication).

Graded the usual way (CODE's files for days 285 and 286, the same start times as backtests 1 and
2, GPS clocks from our 5-day broadcast quadratic for both), range error median / 95% / worst and
orbit 3D median:

| | backtest 1, day 4.5-5 | backtest 2, day 4.0-4.5 |
| --- | --- | --- |
| GPS, ours (IGS ultra-rapid + our force model) | 2.19 / 12.3 / 14.6 m (4.98 m) | 1.25 / 7.49 / 9.9 m (3.07 m) |
| GPS, from CODE's prediction | 1.81 / 12.2 / 15.2 m (1.65 m) | 1.28 / 8.29 / 10.2 m (1.35 m) |
| Galileo, ours (broadcast orbits) | 9.33 / 23.4 / 26.8 m (53.1 m) | 9.04 / 17.5 / 42.9 m (51.9 m) |
| Galileo, from CODE's prediction | 1.68 / 12.7 / 19.3 m (3.72 m) | 0.72 / 15.4 / 35.2 m (1.97 m) |

CODE's orbit prediction is about three times better than ours for GPS and fourteen times better
for Galileo; with it, the error left is nearly all clock. Galileo's tail is clock, on a few
satellites (backtest 2, days 3-4.5: E11 25 m, E09 14 m, E12 13 m of clock against 0.6-2.6 m of
orbit).

**Galileo clocks over 2025** (`clock-backtest --system E` over the 365 daily broadcast files, 27
satellites, 1,392 epochs): the 5-day quadratic gives a day 4-5 median of 1.24 m (GPS: 2.04 m), but
the tail is wide: the IOV satellites E11, E12 and E19 have medians of 11-35 m and 95th
percentiles of 100-570 m, and clock switches make broadcast clocks jump by up to thousands of km
(1.4% of samples above 100 m). Straight lines are far worse for some Galileo clocks, and no window
length wins outright. Galileo wants per-satellite accuracy (or the IOV satellites left out).
