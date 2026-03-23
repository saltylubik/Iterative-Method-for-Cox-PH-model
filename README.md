# Iterative Method for Modelling the Distribution of Regression Effects in the Cox Proportional Hazards Model

**Supplementary material for the SAfJR 2026 poster**

*Lubomír Seif & Ivana Malá — Prague University of Economics and Business*

---

## Overview

This repository contains the full R implementation accompanying the poster
*"Iterative Method for Modelling the Distribution of Regression Effects in the
Cox Proportional Hazards Model"*. The method constructs a smooth plausibility
distribution π(β) over log-hazard-ratio values by bootstrap-aggregating Wald
p-values and transforming them via exponential weighting — without requiring
Bayesian priors.

The core idea: for each point β_k on a fine grid, the mean bootstrap p-value
p̄(β_k) measures how consistent that value is with the data. The probability
measure is then

```
π(β_k) ∝ exp(λ · p̄(β_k))
```

where the concentration coefficient λ controls how sharply π(β) peaks around
the MLE. Two data-driven calibration strategies for λ are provided (see
below).

The method is a semi-parametric extension of iterative frequentist work
originally developed for linear regression.

---

## Repository structure

```
.
├── cox_iterative_supplement.R   # Complete implementation (run this)
├── README.md
```

---

## Requirements

R ≥ 4.2 is recommended. The script installs any missing packages automatically on first run. The packages used are:

| Package      | Purpose                                      |
|--------------|----------------------------------------------|
| `survival`   | Cox model fitting (`coxph`)                  |
| `ggplot2`    | All figures                                  |
| `dplyr`      | Data manipulation                            |
| `tidyr`      | Table reshaping                              |
| `kableExtra` | LaTeX table output                           |
| `scales`     | Axis formatting                              |
| `patchwork`  | Combining ggplot panels                      |

---

## Usage

Clone the repository and run the script from your R console or terminal. All
outputs are written to `./outputs/` (created if it does not exist).

```r
source("cox_iterative_supplement.R")
```

Or from the terminal:

```bash
Rscript cox_iterative_supplement.R
```

### Adjusting the simulation design

All tunable parameters are collected in the `CONFIG` list at the top of the
script (Section 1). You do not need to edit any other part of the code.

```r
CONFIG <- list(
  sim_n_vec    = c(200, 500, 1000),   # sample sizes
  sim_beta_vec = c(log(1.1), log(1.5), log(2.0)),  # true log-HRs
  sim_R_reps   = 100,    # replicates per condition (use 500 for final paper)
  sim_B        = 500,    # bootstrap replicates per replicate
  ...
)
```

A full run with the default settings (100 replicates, B = 500, m = 500) takes
roughly 30–90 minutes depending on hardware. Set `sim_R_reps = 10` and
`sim_B = 100` for a fast test run.

---

## Method summary

### Bootstrap p-value grid

For each bootstrap resample b = 1, …, B:
1. Resample the dataset with replacement (case resampling).
2. Refit the Cox model; record β̂⁽ᵇ⁾ and SE⁽ᵇ⁾.
3. For each grid point β_k compute the Wald p-value
   p_k⁽ᵇ⁾ = 2(1 − Φ(|β̂⁽ᵇ⁾ − β_k| / SE⁽ᵇ⁾)).

Aggregate: p̄_k = mean over b, with a Monte Carlo standard error for
uncertainty quantification.

### Probability measure

```
π(β_k) = exp(λ · p̄_k) / Σ_ℓ exp(λ · p̄_ℓ)
```

The formulation exp(λ p̄) is equivalent to exp(−λ(1 − p̄)) after
normalisation. The quantity (1 − p̄) plays the role of "energy" in a
Boltzmann analogy; 1/λ is the corresponding "temperature".

### Concentration coefficient λ

Two calibration strategies are implemented:

| Strategy | Function | Description |
|---|---|---|
| **CI matching** | `lambda_ci_match()` | Numerical: finds λ such that the mass of π(β) within the Wald 95% CI equals a target (default 0.95). Aligns the new measure with standard frequentist coverage. |
| **Contrast ratio** | `lambda_contrast()` | Closed form: fixes the plausibility ratio π(β̂) / π(β_ref) = R (default R = 50). Makes λ directly interpretable as a signal-to-noise scale. |

### Interval summaries

- **HDI** (`hdi_discrete`): highest density interval — the narrowest interval
  covering a given probability mass.
- **QI** (`quantile_interval`): equal-tailed quantile interval.

---

## Outputs

| File | Contents |
|------|----------|
| `fig_illustration.pdf` | p̄(β) curve with confidence band and π(β) overlay for a single dataset (n = 300, HR = 1.5) |
| `fig_simulation_panel.pdf` | Four-panel simulation summary: HDI coverage, HDI width, calibration accuracy, and λ vs. n |
| `sim_results.rds` | Raw simulation results (one row per replicate × method); load with `readRDS("outputs/sim_results.rds")` |

---

## Validity checks

Section 11 of the script runs 20 automated checks and prints a tally:

```
[ OK ]  A1: p̄ ∈ [0,1]
[ OK ]  A4: pi sums to 1 (ci_match)
...
──────────────────────────────────────────────────────────────────────
  Validity check tally:  20 PASS  |  0 WARN  |  0 FAIL  (of 20)
──────────────────────────────────────────────────────────────────────
```

Checks labelled `WARN` are expected to fail occasionally on a single sample
(e.g. the true β not being covered by one HDI) and do not indicate a bug.
`FAIL` indicates a mathematical or numerical problem that should be
investigated.

---

## Citation

If you use this code, please cite the accompanying poster:

```
Seif, L. & Malá, I. (2026). Iterative Method for Modelling the Distribution
of Regression Effects in the Cox Proportional Hazards Model.
Poster presented at SAfJR 2026, Prague University of Economics and Business.
```

---

## License

MIT — see `LICENSE` for details.
