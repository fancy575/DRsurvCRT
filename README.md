# DRsurvCRT

Doubly-robust survival estimation for **cluster-randomized trials** (CRTs) with right censoring.

Supports:

- **SPCE** (survival probability contrast estimate)
- **RMST** (restricted mean survival time) at user-specified target times `tau`
- **Marginal** Cox or **frailty** (gamma) outcome/censoring models
- Optional **leave-one-cluster-out jackknife** variances on a shared grid

> Core entry point: `DRsurvfit()` with S3 methods `summary()` and `plot()`.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("fancy575/DRsurvCRT")
```

### Dependencies

- `survival` (Cox models, `Surv()`, `cluster()`)
- `frailtyEM` (gamma frailty fits; only when `method = "frailty"`)
- `Rcpp` (C++ backends; the package registers native routines)

---

## Example dataset: `dat`

A toy CRT survival dataset bundled with the package to demonstrate `DRsurvfit`.

**Columns**

- `M`: cluster ID (integer or factor)
- `A`: cluster-level treatment (0/1), **constant within** `M`
- `W1`, `W2`: cluster-level covariates, **constant within** `M`
- `Z1`, `Z2`: individual-level covariates
- `time`: follow-up time
- `event`: event indicator (1 = event, 0 = censored)

**Load & sanity checks**

```r
data(dat)  # from DRsurvCRT

# Check cluster-constancy of cluster-level variables (should be TRUE)
all(tapply(dat$A,  dat$M, function(x) length(unique(x)) == 1))
all(tapply(dat$W1, dat$M, function(x) length(unique(x)) == 1))
all(tapply(dat$W2, dat$M, function(x) length(unique(x)) == 1))
```

---

## Quick start

```r
library(DRsurvCRT)
library(survival)

data(dat)

# --- SPCE on the full grid (marginal Cox) ---
fit_spce <- DRsurvfit(
  data = dat,
  formula = Surv(time, event) ~ W1 + W2 + Z1 + Z2 + cluster(M),
  intv = "A",
  method = "marginal",
  estimand = "SPCE"
)

summary(fit_spce)                  # S1, S0, S1 - S0 at selected times
plot(fit_spce, level = "cluster")  # survival curves (S1 solid, S0 dashed)

# --- RMST at landmarks with jackknife variance (frailty models) ---
fit_rmst <- DRsurvfit(
  data = dat,
  formula = Surv(time, event) ~ W1 + W2 + Z1 + Z2 + cluster(M),
  intv = "A",
  method = "frailty",
  estimand = "RMST",
  tau = c(0.5, 1, 2),
  variance = "jackknife"
)

summary(fit_rmst)
plot(fit_rmst, level = "individual")
```

---

## API

### `DRsurvfit()`

```r
DRsurvfit(
  data,
  formula,                 # Surv(time, event) ~ covariates + cluster(M)
  cens_formula = NULL,     # optional; default uses RHS (no cluster()) with event==0
  intv,                    # cluster-level treatment (0/1), constant within cluster
  method   = c("marginal", "frailty"),
  estimand = c("SPCE", "RMST"),
  tau = NULL,              # numeric vector for RMST landmarks
  trt_prob = NULL,         # optional c(p0, p1); otherwise inferred from first row per cluster
  variance = c("none","jackknife"),
  fit_controls = NULL,     # frailtyEM::emfrail_control() when method = "frailty"
  strata = NULL,           # cluster id if not present via cluster() in formula
  spce_summary_times = NULL,
  verbose = FALSE
)
```

**Returns** (class `"DRsurvfit"`):

- Always computes on the **shared event-time grid**.
- For **SPCE**: `S_full_cluster`, `S_full_ind`, and displayed slices in `summary()`.
- For **RMST**: `RMST_cluster`, `RMST_ind` at each `tau` (rows named by `tau`).
- Optional jackknife: `var_cluster`, `var_ind` (variances for S1, S0, diff; plus `cov`).
- Common fields: `event_time`, `trt_prob = c(p0, p1)`, and `tau` (for RMST).

### `summary.DRsurvfit()`

- **SPCE**: prints `(S1, S0, S1 − S0)` at `times` (or `spce_summary_times`; else default grid quantiles).
- **RMST**: prints `(R1, R0, R1 − R0)` for each `tau`.

### `plot.DRsurvfit()`

- `plot(x, level = c("cluster","individual"))` draws S1 (solid) and S0 (dashed) across time; vertical lines at `tau` if present.

---

## Modeling details

- **Outcome/Censoring models**
  - `method = "marginal"`: `survival::coxph()`
  - `method = "frailty"`: `frailtyEM::emfrail()` with gamma frailty
- **Censoring formula**: if `cens_formula` is `NULL`, uses `Surv(time, event == 0) ~ <RHS without cluster()>`; appends `+ cluster(<id>)` if a cluster id is known.
- **Treatment probabilities**: set via `trt_prob = c(p0, p1)`; otherwise computed from the first observation per cluster.
- **Jackknife**: leave-one-cluster-out on the same shared event-time grid; RMST reuses the core’s RMST-at-τ (no re-patching in the variance step).

---

## Practical notes

- `intv` must be **binary** and **constant within clusters**.
- For `method = "frailty"`, include `cluster(<id>)` in `formula` or pass `strata = "<id>"`.
- Choose `tau` within the observed time range for stable RMST.
- You can pass `spce_summary_times` to control which times appear in `summary()` for SPCE.

---

## Citation

If you use **DRsurvCRT** in published work, please cite:


---

## License

MIT. See `LICENSE`.

---

## Contributing

Issues and pull requests are welcome. Please attach a minimal reproducible example and your `sessionInfo()` when reporting bugs.
