# DRsurvCRT 
Doubly-robust estimation of survival outcomes in **cluster randomized trials (CRTs)** with
support for **multi-state events** and **Restricted Mean Time in Favor (RMT-IF)** and
**state-specific prioritized composite estimands (SPCE)**.

The package is designed for CRTs where:

- Treatment is assigned at the **cluster level** (e.g., hospitals, practices).
- Outcomes are right-censored event times at the **individual level**.
- The event indicator can be
  - **binary** (single state: 0 = censored, 1 = event), or
  - **multi-state** (0 = censored / not yet failed, 1,2,3,… = ordered states).

Key features:

- Doubly-robust (DR) estimation via augmentation + inverse probability weighting.
- Support for **marginal Cox** and **shared frailty (gamma) Cox** nuisance models.
- Joint estimation of:
  - **SPCE(t)**: state-specific survival probabilities by arm over time.
  - **RMT-IF(t)**: cumulative win probabilities based on DR survival curves.
- Both **cluster-level** and **individual-level** estimands.
- Leave-one-cluster-out **jackknife variance** with CRT degrees of freedom.
- Convenient **summary** and **ggplot-based plot** methods.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("fancy575/DRsurvCRT")
```

Load the package:

```r
library(DRsurvCRT)
```

---

## Data Included in the Package

Two example datasets are bundled in the package.

### 1. `dats`: Single-state CRT survival data

A synthetic CRT dataset with **one event type** (`event` ∈ {0,1}).

- `cluster` : Cluster id (e.g., practice/hospital).
- `trt`     : Cluster-level treatment (0 = control, 1 = treated).
- `W1, W2`  : Cluster-level or baseline covariates.
- `Z1, Z2`  : Individual-level covariates.
- `time`    : Observed follow-up time.
- `event`   : Event indicator (0 = censored, 1 = event).

You can inspect it with:

```r
data("dats")
str(dats)
```

### 2. `datm`: Multi-state CRT survival data

A synthetic CRT dataset with **multi-state events** (`event` ∈ {0,1,2,3}). Each subject can
experience multiple transitions, stored in **long format**.

- `id`      : Individual identifier.
- `time`    : Time of transition or censoring.
- `event`   : State indicator:
  - 0 = still at risk / censored at that row,
  - 1,2,3 = different ordered event types (e.g., prioritized composite states).
- `cluster` : Cluster id.
- `trt`     : Cluster-level treatment (0/1).
- `W1, W2`  : Cluster-level covariates.
- `Z1, Z2`  : Individual-level covariates (may vary over time in general).

Example:

```r
data("datm")
str(datm)
```

Internally, the DRsurvCRT core loops over `s = 1, …, max(event)` and, for each `s`,
constructs a binary event indicator and fits DR survival curves in the same way as in the
single-state case. These per-state survival curves are then assembled into **state-specific
SPCE(t)** and used to construct **RMT-IF(t)**.

---

## Main Function: `DRsurvfit()`

The central user-facing function is:

```r
fit <- DRsurvfit(
  formula,
  data,
  intv,
  method   = c("marginal", "frailty"),
  tau      = NULL,
  variance = c("none", "jackknife"),
  surv_cov   = NULL,
  censor_cov = NULL,
  strata   = NULL,
  verbose  = FALSE
)
```

### Arguments

- `formula`  
  A survival formula of the form  
  `Surv(time, event) ~ trt + W1 + W2 + Z1 + Z2 + cluster(cluster)`  
  where:
  - `time`  is the follow-up time variable.
  - `event` is the event indicator:
    - For single state: 0 = censored, 1 = event.
    - For multi-state: 0 = censored/still at risk, 1,2,3,… = ordered states.
  - `cluster()` must appear on the RHS (or you can provide `strata=` separately).

- `data`  
  A data.frame containing all variables in the formula. For multi-state data (`datm`), each
  row is a time interval or transition for a subject (`id`).

- `intv`  
  Character name of the **treatment column** (e.g., `"trt"`), constant within each cluster.

- `method`  
  - `"marginal"`: Cox proportional hazards model at the individual level with robust/cluster
    variance.
  - `"frailty"`: Shared frailty Cox model (`frailtyEM::emfrail`) with gamma frailty.

- `tau`  
  Numeric vector of **truncation times** at which you want summaries of:
  - **SPCE(t)** (survival difference at time `t`), and
  - **RMT-IF(t)** (integrated win probabilities up to `t`).  
  If `tau = NULL`, summary and plot methods use internal defaults (25%, 50%, 75% quantiles of
  the observed event-time distribution).

- `variance`  
  - `"none"`: Point estimates only.
  - `"jackknife"`: Leave-one-cluster-out jackknife variances, with degrees of freedom
    `M - 1` where `M` is the number of clusters.

- `surv_cov`, `censor_cov`  
  Character vectors of covariate names used in the DR augmentation models for the survival
  and censoring processes, respectively. Defaults reuse all non-treatment covariates on the
  RHS of the formula.

- `strata`  
  Optional character giving the cluster-id column if `cluster()` is not in the formula.

- `verbose`  
  Logical; print progress messages during fitting.

### Returned Object

`DRsurvfit()` returns an object of class `"DRsurvfit"` with (main slots):

- `event_time`       : Shared event-time grid.
- `S_state_cluster`  : Array `[time, 2, S]` with cluster-level survival by state (S1,S0).
- `S_state_ind`      : Array `[time, 2, S]` with individual-level survival by state (S1,S0).
- `SPCE_cluster`     : Matrix `[time, S]` of state-specific survival differences S1–S0.
- `SPCE_ind`         : Same, individual-level.
- `RMTIF_cluster`    : Matrix `[time, 3]` with columns `(R1, R0, R1–R0)` aggregated across all states.
- `RMTIF_ind`        : Same, individual-level.
- `var_cluster`      : Jackknife variances (if requested) for all time points.
- `var_ind`          : Same, individual-level.
- `trt_prob`         : Treatment probabilities `(p0, p1)` used in DR weighting.
- `n_clusters`, `df_jackknife`, `n_obs`, `cluster_col`, etc.

These slots are primarily accessed via `summary()` and `plot()`.

---

## Basic Workflow

### 1. Single-State CRT Example (`dats`)

```r
library(DRsurvCRT)
data("dats")

# Inspect data
head(dats)

# Fit DRsurvfit using a marginal Cox model
fit_single <- DRsurvfit(
  formula  = Surv(time, event) ~ trt + W1 + W2 + Z1 + Z2 + cluster(cluster),
  data     = dats,
  intv     = "trt",
  method   = "marginal",
  tau      = c(1, 2, 3),        # truncation times for summaries
  variance = "jackknife",
  verbose  = TRUE
)
```

#### Summary for Single-State

The `summary()` method reports **SPCE(t)** and **RMT-IF(t)** at the truncation times `tau`,
with t-based confidence intervals using the CRT jackknife degrees of freedom.

```r
# Cluster-level SPCE and RMT-IF at t = 1, 2, 3
summary(
  fit_single,
  tau    = c(1, 2, 3),
  level  = "cluster",    # "cluster", "individual", or "both"
  type   = "both",       # "spce", "rmtif", or "both"
  alpha  = 0.05,
  digits = 3
)
```

If you call `summary(fit_single)` without `tau`, it will automatically choose the 25%, 50%, and
75% quantiles of the observed event-time grid.

#### Plot for Single-State

The `plot()` method uses **ggplot2** to display **only the treatment–control difference** with
confidence bands.

- For SPCE(t): difference of survival probabilities by arm.
- For RMT-IF(t): difference of cumulative “time in favor” (sum across states).

```r
# SPCE(t) difference (cluster-level) with 95% CI over the full time range
plot(
  fit_single,
  type  = "spce",
  level = "cluster",
  tau   = NULL,      # NULL = full event-time range
  alpha = 0.05
)

# RMT-IF(t) difference (cluster-level) truncated at t = 3 with 95% CI
plot(
  fit_single,
  type  = "rmtif",
  level = "cluster",
  tau   = 3,         # truncate the plot at t = 3
  alpha = 0.05
)
```

You can switch to individual-level curves via `level = "individual"`.

---

### 2. Multi-State CRT Example (`datm`)

```r
library(DRsurvCRT)
data("datm")

# Inspect structure
head(datm)

# Fit DRsurvfit using a shared gamma frailty model
fit_multi <- DRsurvfit(
  formula  = Surv(time, event) ~ trt + W1 + W2 + Z1 + Z2 + cluster(cluster),
  data     = datm,
  intv     = "trt",
  method   = "frailty",
  tau      = c(2, 5, 10),   # truncation times for multi-state RMT-IF and SPCE
  variance = "jackknife",
  verbose  = TRUE
)
```

Internally, `event` is treated as a multi-state indicator (`0,1,2,3` in `datm`).
The core constructs, for each state `s`:

- a state-specific binary event indicator,
- DR survival curves for treatment and control, and
- state-specific SPCE(t) and contributions to RMT-IF(t).

All states are combined for the **overall RMT-IF**; SPCE is reported **per state**.

#### Summary for Multi-State

```r
# Cluster-level, both SPCE and RMT-IF at t = 2, 5, 10
summary(
  fit_multi,
  tau    = c(2, 5, 10),
  level  = "cluster",
  type   = "both",
  alpha  = 0.05,
  digits = 3
)
```

For SPCE, the table shows (for each state and time `tau`):

- S1(t,s) : treatment survival probability at time t for state s,
- S0(t,s) : control survival probability at time t for state s,
- S1(t,s) − S0(t,s) and its CI.

For RMT-IF, the table reports:

- R1(t) : total time in favor of treatment up to t,
- R0(t) : total time in favor of control up to t,
- R1(t) − R0(t) and its CI (aggregated across all states).

#### Plot for Multi-State

```r
# SPCE(t) differences by state (cluster-level)
plot(
  fit_multi,
  type  = "spce",
  level = "cluster",
  tau   = NULL,      # full time range for each state
  alpha = 0.05
)

# RMT-IF(t) overall difference (sum across states), truncated at t = 10
plot(
  fit_multi,
  type  = "rmtif",
  level = "cluster",
  tau   = 10,
  alpha = 0.05
)
```

For `type = "spce"`, the “bouquet” plot overlays the SPCE(t) differences across states with
confidence ribbons (state as color).  
For `type = "rmtif"`, a single curve for the **overall** RMT-IF(t) difference is drawn with a
CI band.

---

## Reproducible Example

Below is a fully reproducible mini-example assuming the package is installed:

```r
library(DRsurvCRT)

## Single-state example
data("dats")
fit_single <- DRsurvfit(
  Surv(time, event) ~ trt + W1 + W2 + Z1 + Z2 + cluster(cluster),
  data     = dats,
  intv     = "trt",
  method   = "marginal",
  tau      = c(1, 2, 3),
  variance = "jackknife"
)

summary(fit_single, tau = c(1, 2, 3), level = "cluster", type = "both")
plot(fit_single, type = "spce",  level = "cluster", tau = NULL)
plot(fit_single, type = "rmtif", level = "cluster", tau = 3)

## Multi-state example
data("datm")
fit_multi <- DRsurvfit(
  Surv(time, event) ~ trt + W1 + W2 + Z1 + Z2 + cluster(cluster),
  data     = datm,
  intv     = "trt",
  method   = "frailty",
  tau      = c(2, 5, 10),
  variance = "jackknife"
)

summary(fit_multi, tau = c(2, 5, 10), level = "cluster", type = "both")
plot(fit_multi, type = "spce",  level = "cluster", tau = NULL)
plot(fit_multi, type = "rmtif", level = "cluster", tau = 10)
```

---

## Citation

If you use **DRsurvCRT** in your work, please cite it as:

> Fang et al. (2025). *DRsurvCRT: Doubly Robust Estimation for Survival Outcomes in Cluster Randomized Trials*. R package version 0.1.0.

---

## Issues and Contributions

Bug reports, feature requests, and pull requests are welcome at:

- GitHub: <https://github.com/fancy575/DRsurvCRT>
