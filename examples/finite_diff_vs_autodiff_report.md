# Finite-Difference vs ForwardDiff Jacobian — Benchmark Report

> Auto-generated on 2026-02-25 17:34:14
> Julia 1.11.9 — x86_64-w64-mingw32

## Default Tolerances

### Performance

| Scenario | FD Median (ms) | AD Median (ms) | Speedup | FD Allocs | AD Allocs | FD Mem (KiB) | AD Mem (KiB) |
|:---------|---------------:|---------------:|--------:|----------:|----------:|-------------:|-------------:|
| LEO (AMAZONIA 1) — no initial guess | 250.7 | 639.2 | 2.55x FD | 2 | 332334 | 1 | 1204704 |
| HEO (MOLNIYA 1-83) — no initial guess | 13290.4 | 130.2 | 102.12x AD | 2 | 31792 | 1 | 115240 |
| LEO (AMAZONIA 1) — TLE initial guess | 8.3 | 7.4 | 1.11x AD | 2 | 4006 | 1 | 14515 |
| HEO (MOLNIYA 1-83) — TLE initial guess | 5.3 | 4.7 | 1.13x AD | 2 | 1158 | 1 | 4191 |

### Accuracy (Absolute Error vs Reference TLE)

#### LEO (AMAZONIA 1) — no initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | 1.986000000000e-03 | 2.9098e-07 | 2.9669e-07 | FD |
| `eccentricity` | 1.597000000000e-04 | 2.1677e-11 | 2.0739e-11 | AD |
| `inclination` | 9.848890000000e+01 | 1.7788e-10 | 3.4106e-13 | AD |
| `raan` | 3.446059000000e+02 | 2.4824e-10 | 7.2760e-12 | AD |
| `argument_of_perigee` | 7.442440000000e+01 | 3.4058e-06 | 3.5994e-06 | FD |
| `mean_anomaly` | 2.857135000000e+02 | 3.4194e-06 | 3.6132e-06 | FD |
| `mean_motion` | 1.440801240000e+01 | 1.2466e-09 | 1.2702e-09 | FD |

#### HEO (MOLNIYA 1-83) — no initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | -1.352500000000e-04 | 6.8312e-08 | 1.2042e-08 | AD |
| `eccentricity` | 7.421690000000e-01 | 1.3213e-11 | 5.4731e-11 | FD |
| `inclination` | 6.217490000000e+01 | 1.1395e-08 | 1.3401e-11 | AD |
| `raan` | 1.980096000000e+02 | 5.9339e-08 | 4.0330e-11 | AD |
| `argument_of_perigee` | 2.530462000000e+02 | 2.5049e-08 | 2.8433e-09 | AD |
| `mean_anomaly` | 2.015610000000e+01 | 6.0885e-08 | 1.3307e-08 | AD |
| `mean_motion` | 2.012699940000e+00 | 3.2725e-10 | 4.5816e-11 | AD |

#### LEO (AMAZONIA 1) — TLE initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | 1.986000000000e-03 | 2.9090e-07 | 2.9669e-07 | FD |
| `eccentricity` | 1.597000000000e-04 | 2.1680e-11 | 2.0739e-11 | AD |
| `inclination` | 9.848890000000e+01 | 1.7843e-10 | 3.4106e-13 | AD |
| `raan` | 3.446059000000e+02 | 2.4886e-10 | 7.2760e-12 | AD |
| `argument_of_perigee` | 7.442440000000e+01 | 3.4039e-06 | 3.5993e-06 | FD |
| `mean_anomaly` | 2.857135000000e+02 | 3.4174e-06 | 3.6131e-06 | FD |
| `mean_motion` | 1.440801240000e+01 | 1.2463e-09 | 1.2702e-09 | FD |

#### HEO (MOLNIYA 1-83) — TLE initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | -1.352500000000e-04 | 4.0804e-09 | 1.2042e-08 | FD |
| `eccentricity` | 7.421690000000e-01 | 2.9982e-11 | 5.4731e-11 | FD |
| `inclination` | 6.217490000000e+01 | 5.2099e-09 | 1.3408e-11 | AD |
| `raan` | 1.980096000000e+02 | 2.6194e-08 | 4.0330e-11 | AD |
| `argument_of_perigee` | 2.530462000000e+02 | 1.3657e-08 | 2.8433e-09 | AD |
| `mean_anomaly` | 2.015610000000e+01 | 1.4843e-08 | 1.3307e-08 | AD |
| `mean_motion` | 2.012699940000e+00 | 3.8728e-11 | 4.5814e-11 | FD |

---

## Tight Tolerances (atol = rtol = 1e-14)

### Performance

| Scenario | FD Median (ms) | AD Median (ms) | Speedup | FD Allocs | AD Allocs | FD Mem (KiB) | AD Mem (KiB) |
|:---------|---------------:|---------------:|--------:|----------:|----------:|-------------:|-------------:|
| LEO (AMAZONIA 1) — no initial guess | 20753.8 | 19026.3 | 1.09x AD | 2 | 10010002 | 1 | 36286251 |
| HEO (MOLNIYA 1-83) — no initial guess | 28219.6 | 26958.2 | 1.05x AD | 2 | 5780002 | 1 | 20952501 |
| LEO (AMAZONIA 1) — TLE initial guess | 20853.7 | 19226.6 | 1.08x AD | 2 | 10010002 | 1 | 36286251 |
| HEO (MOLNIYA 1-83) — TLE initial guess | 26891.4 | 23947.5 | 1.12x AD | 2 | 5780002 | 1 | 20952501 |

### Accuracy (Absolute Error vs Reference TLE)

#### LEO (AMAZONIA 1) — no initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | 1.986000000000e-03 | 2.9099e-07 | 2.9669e-07 | FD |
| `eccentricity` | 1.597000000000e-04 | 2.1676e-11 | 2.0739e-11 | AD |
| `inclination` | 9.848890000000e+01 | 1.7789e-10 | 3.4106e-13 | AD |
| `raan` | 3.446059000000e+02 | 2.4824e-10 | 7.2760e-12 | AD |
| `argument_of_perigee` | 7.442440000000e+01 | 3.4057e-06 | 3.5993e-06 | FD |
| `mean_anomaly` | 2.857135000000e+02 | 3.4193e-06 | 3.6131e-06 | FD |
| `mean_motion` | 1.440801240000e+01 | 1.2466e-09 | 1.2702e-09 | FD |

#### HEO (MOLNIYA 1-83) — no initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | -1.352500000000e-04 | 1.1174e-06 | 1.2042e-08 | AD |
| `eccentricity` | 7.421690000000e-01 | 3.2745e-10 | 5.4731e-11 | AD |
| `inclination` | 6.217490000000e+01 | 3.9307e-08 | 1.3401e-11 | AD |
| `raan` | 1.980096000000e+02 | 2.1117e-07 | 4.0330e-11 | AD |
| `argument_of_perigee` | 2.530462000000e+02 | 3.0063e-08 | 2.8433e-09 | AD |
| `mean_anomaly` | 2.015610000000e+01 | 1.1598e-06 | 1.3307e-08 | AD |
| `mean_motion` | 2.012699940000e+00 | 5.5268e-09 | 4.5814e-11 | AD |

#### LEO (AMAZONIA 1) — TLE initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | 1.986000000000e-03 | 2.9098e-07 | 2.9669e-07 | FD |
| `eccentricity` | 1.597000000000e-04 | 2.1677e-11 | 2.0739e-11 | AD |
| `inclination` | 9.848890000000e+01 | 1.7789e-10 | 3.2685e-13 | AD |
| `raan` | 3.446059000000e+02 | 2.4824e-10 | 7.3328e-12 | AD |
| `argument_of_perigee` | 7.442440000000e+01 | 3.4058e-06 | 3.5994e-06 | FD |
| `mean_anomaly` | 2.857135000000e+02 | 3.4193e-06 | 3.6131e-06 | FD |
| `mean_motion` | 1.440801240000e+01 | 1.2466e-09 | 1.2702e-09 | FD |

#### HEO (MOLNIYA 1-83) — TLE initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | -1.352500000000e-04 | 1.1176e-06 | 1.2042e-08 | AD |
| `eccentricity` | 7.421690000000e-01 | 3.2748e-10 | 5.4730e-11 | AD |
| `inclination` | 6.217490000000e+01 | 3.9312e-08 | 1.3408e-11 | AD |
| `raan` | 1.980096000000e+02 | 2.1119e-07 | 4.0302e-11 | AD |
| `argument_of_perigee` | 2.530462000000e+02 | 3.0068e-08 | 2.8433e-09 | AD |
| `mean_anomaly` | 2.015610000000e+01 | 1.1599e-06 | 1.3307e-08 | AD |
| `mean_motion` | 2.012699940000e+00 | 5.5274e-09 | 4.5816e-11 | AD |

---

