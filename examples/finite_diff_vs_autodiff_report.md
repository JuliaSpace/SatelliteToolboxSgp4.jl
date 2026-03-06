# Finite-Difference vs ForwardDiff Jacobian — Benchmark Report

> Auto-generated on 2026-02-25 19:52:20
> Julia 1.11.9 — x86_64-w64-mingw32

## Default Tolerances

### Performance

| Scenario | FD Median (ms) | AD Median (ms) | Speedup | FD Allocs | AD Allocs | FD Mem (KiB) | AD Mem (KiB) |
|:---------|---------------:|---------------:|--------:|----------:|----------:|-------------:|-------------:|
| LEO (AMAZONIA 1) — no initial guess | 262.2 | 378.8 | 1.44x FD | 2 | 4 | 1 | 8 |
| HEO (MOLNIYA 1-83) — no initial guess | 13942.4 | 120.6 | 115.57x AD | 2 | 4 | 1 | 8 |
| LEO (AMAZONIA 1) — TLE initial guess | 8.8 | 5.9 | 1.49x AD | 2 | 4 | 1 | 8 |
| HEO (MOLNIYA 1-83) — TLE initial guess | 5.7 | 4.5 | 1.29x AD | 2 | 4 | 1 | 8 |

### Accuracy (Absolute Error vs Reference TLE)

#### LEO (AMAZONIA 1) — no initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | 1.986000000000e-03 | 2.9098e-07 | 2.9669e-07 | FD |
| `eccentricity` | 1.597000000000e-04 | 2.1677e-11 | 2.0739e-11 | AD |
| `inclination` | 9.848890000000e+01 | 1.7788e-10 | 3.2685e-13 | AD |
| `raan` | 3.446059000000e+02 | 2.4824e-10 | 7.3328e-12 | AD |
| `argument_of_perigee` | 7.442440000000e+01 | 3.4058e-06 | 3.5994e-06 | FD |
| `mean_anomaly` | 2.857135000000e+02 | 3.4194e-06 | 3.6132e-06 | FD |
| `mean_motion` | 1.440801240000e+01 | 1.2466e-09 | 1.2702e-09 | FD |

#### HEO (MOLNIYA 1-83) — no initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | -1.352500000000e-04 | 6.8312e-08 | 1.2042e-08 | AD |
| `eccentricity` | 7.421690000000e-01 | 1.3213e-11 | 5.4730e-11 | FD |
| `inclination` | 6.217490000000e+01 | 1.1395e-08 | 1.3408e-11 | AD |
| `raan` | 1.980096000000e+02 | 5.9339e-08 | 4.0330e-11 | AD |
| `argument_of_perigee` | 2.530462000000e+02 | 2.5049e-08 | 2.8433e-09 | AD |
| `mean_anomaly` | 2.015610000000e+01 | 6.0885e-08 | 1.3307e-08 | AD |
| `mean_motion` | 2.012699940000e+00 | 3.2725e-10 | 4.5818e-11 | AD |

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
| LEO (AMAZONIA 1) — no initial guess | 22737.5 | 15936.1 | 1.43x AD | 2 | 4 | 1 | 8 |
| HEO (MOLNIYA 1-83) — no initial guess | 29576.2 | 22522.1 | 1.31x AD | 2 | 4 | 1 | 8 |
| LEO (AMAZONIA 1) — TLE initial guess | 22739.7 | 15932.0 | 1.43x AD | 2 | 4 | 1 | 8 |
| HEO (MOLNIYA 1-83) — TLE initial guess | 29301.5 | 22458.9 | 1.30x AD | 2 | 4 | 1 | 8 |

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
| `mean_motion` | 2.012699940000e+00 | 5.5268e-09 | 4.5816e-11 | AD |

#### LEO (AMAZONIA 1) — TLE initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | 1.986000000000e-03 | 2.9098e-07 | 2.9669e-07 | FD |
| `eccentricity` | 1.597000000000e-04 | 2.1677e-11 | 2.0739e-11 | AD |
| `inclination` | 9.848890000000e+01 | 1.7789e-10 | 3.2685e-13 | AD |
| `raan` | 3.446059000000e+02 | 2.4824e-10 | 7.2760e-12 | AD |
| `argument_of_perigee` | 7.442440000000e+01 | 3.4058e-06 | 3.5993e-06 | FD |
| `mean_anomaly` | 2.857135000000e+02 | 3.4193e-06 | 3.6131e-06 | FD |
| `mean_motion` | 1.440801240000e+01 | 1.2466e-09 | 1.2702e-09 | FD |

#### HEO (MOLNIYA 1-83) — TLE initial guess

| Field | Reference | FD Error | AD Error | Winner |
|:------|----------:|---------:|---------:|:------:|
| `bstar` | -1.352500000000e-04 | 1.1176e-06 | 1.2042e-08 | AD |
| `eccentricity` | 7.421690000000e-01 | 3.2748e-10 | 5.4730e-11 | AD |
| `inclination` | 6.217490000000e+01 | 3.9312e-08 | 1.3401e-11 | AD |
| `raan` | 1.980096000000e+02 | 2.1119e-07 | 4.0330e-11 | AD |
| `argument_of_perigee` | 2.530462000000e+02 | 3.0068e-08 | 2.8434e-09 | AD |
| `mean_anomaly` | 2.015610000000e+01 | 1.1599e-06 | 1.3307e-08 | AD |
| `mean_motion` | 2.012699940000e+00 | 5.5274e-09 | 4.5811e-11 | AD |

---

