# Clock prediction over a year: 2025

`extended-ephemeris clock-backtest` over all 365 daily broadcast files of 2025 (BKG,
`root_ftp/IGS/BRDC/2025/<doy>/BRDC00IGS_R_2025<doy>0000_01D_MN.rnx.gz`, 0.74 GB): every 6 hours,
for every healthy GPS satellite, each clock model is fitted to the broadcast clocks of the
preceding days and compared hourly over the next seven days with the clock the satellite actually
broadcast. 32 satellites, 1,393 epochs. Errors in metres of range (c × clock error). Runs in about
9 seconds.

```
clock prediction error against broadcast clocks, metres (median / 95% / 99%), by forecast day
model     | day 0-1             | day 1-2             | day 2-3             | day 3-4             | day 4-5             | day 5-6             | day 6-7             | >100 m
line:1    |  0.27 /  1.41 /   2.9 |  1.05 /  4.46 /   9.7 |  2.19 /  9.45 /  19.8 |  3.66 / 16.45 /  32.5 |  5.41 / 25.37 /  49.2 |  7.43 / 36.10 /  69.8 |  9.77 / 48.25 /  94.4 | 0.194%
line:2    |  0.47 /  2.32 /   5.1 |  1.40 /  6.35 /  13.5 |  2.65 / 12.32 /  24.8 |  4.21 / 20.26 /  39.6 |  6.05 / 29.97 /  58.6 |  8.19 / 41.34 /  81.2 | 10.67 / 53.91 / 107.9 | 0.294%
quad:3    |  0.22 /  1.09 /   2.5 |  0.62 /  2.79 /   7.9 |  1.17 /  5.43 /  14.4 |  1.90 /  9.02 /  21.7 |  2.80 / 13.43 /  30.9 |  3.87 / 18.58 /  42.2 |  5.12 / 24.23 /  57.6 | 0.125%
quad:4    |  0.22 /  1.13 /   3.0 |  0.56 /  2.82 /   8.1 |  1.00 /  5.38 /  14.8 |  1.58 /  8.76 /  21.8 |  2.29 / 12.87 /  30.6 |  3.11 / 17.57 /  41.5 |  4.05 / 22.78 /  55.7 | 0.092%
quad:5    |  0.23 /  1.26 /   3.4 |  0.54 /  3.09 /   8.4 |  0.93 /  5.72 /  15.0 |  1.44 /  9.14 /  21.9 |  2.04 / 13.17 /  30.9 |  2.72 / 17.72 /  41.6 |  3.49 / 22.58 /  56.0 | 0.066%
quad:6    |  0.25 /  1.48 /   4.0 |  0.54 /  3.43 /   9.4 |  0.91 /  6.04 /  16.3 |  1.37 /  9.30 /  23.8 |  1.90 / 13.03 /  34.0 |  2.51 / 17.26 /  46.5 |  3.17 / 21.77 /  61.2 | 0.051%
quad:7    |  0.27 /  1.71 /   4.7 |  0.56 /  3.66 /  10.8 |  0.91 /  6.19 /  17.7 |  1.34 /  9.18 /  25.7 |  1.82 / 12.60 /  36.1 |  2.36 / 16.52 /  49.0 |  2.96 / 21.08 /  63.4 | 0.052%
quad:10   |  0.35 /  2.29 /   7.8 |  0.62 /  4.28 /  14.5 |  0.93 /  6.69 /  22.1 |  1.29 /  9.56 /  30.9 |  1.68 / 12.86 /  40.7 |  2.10 / 16.71 /  51.5 |  2.56 / 21.07 /  62.9 | 0.149%

quad:5 at day 4-5, by month of the prediction: median / 95% / 99% / worst (m)
   1:  2.38 /  11.76 /   29.1 /   122.8
   2:  2.23 /  12.28 /   25.6 /   126.9
   3:  2.49 /  11.65 /   27.8 /    71.6
   4:  2.23 /  15.10 /   29.2 /    77.5
   5:  1.80 /  13.10 /   26.6 /    70.2
   6:  2.16 /  13.54 /   42.7 /    83.6
   7:  2.28 /  13.44 /   28.8 /    70.7
   8:  1.92 /  13.25 /   35.8 /    78.6
   9:  1.59 /  12.66 /   33.5 /    84.4
  10:  1.92 /  13.56 /   50.4 /   267.9
  11:  1.57 /  13.01 /   39.4 /    73.1
  12:  2.25 /  13.77 /   23.9 /   113.9

quad:5 at day 4-5, by satellite (worst 95th percentiles first): median / 95% / worst (m)
  G03: 20.70 /  55.35 /    77.5
  G26:  2.36 /  20.13 /    83.6
  G13:  6.34 /  18.82 /    71.1
  G07:  6.46 /  17.29 /    35.4
  G29:  6.07 /  16.41 /    27.3
  G31:  6.09 /  13.84 /    23.0
  G27:  1.65 /  12.87 /   267.9
  G21:  2.62 /   9.88 /    25.2
  G02:  2.56 /   7.71 /    61.2
  G16:  1.84 /   7.08 /    84.4
  G20:  2.01 /   6.62 /    27.1
  G01:  1.96 /   6.57 /   122.8
  ... best: G23 4.05, G09 4.24, G25 4.36, G32 4.42, G06 4.42
```

Findings:

- **The 5-day quadratic (the predictor's default) holds over the whole year**, as the two October
  backtests suggested: day 4-5 median 2.0 m, 95th percentile 13.2 m; a 1-day line is 5.4 / 25.4 m.
- **The best window depends on the horizon.** 3-4 days of history are best for days 0-3 (quad:4
  at day 2-3: 1.00 / 5.4 m); 7 days for days 4-7 (quad:7 at day 6-7: 2.96 / 21.1 m against
  quad:5's 3.49 / 22.6 m). Choosing the window by the record's age would gain roughly 5-15%.
- **A few satellites dominate the tail.** G03 (day 4-5 median 20.7 m), then G13, G07, G29 and G31
  (about 6 m) against about 2 m for the rest. The URA sent today is the same for every satellite;
  setting it per satellite from a year of predictability would let receivers weight them down.
- **No seasonality**: the day 4-5 95th percentile stays between 11.7 and 15.1 m in every month.
  The rare large errors (up to 268 m) are clock events no model anticipates.
