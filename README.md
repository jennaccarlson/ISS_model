# ISS Model: Estimating Insulin Sensitivity and Beta-Cell Function from the OGTT

R implementation of the Insulin Sensitivity and Secretion (ISS) model for fitting Oral Glucose Tolerance Test (OGTT) data. This code estimates two key metabolic parameters for each subject:

- **σ (sigma)** — a unitless measure of beta-cell function (BCF)
- **SI** — peripheral insulin sensitivity (ml/µU·day)

The ISS model and its validation against IVGTT-derived measures are described in:

> Ha J, et al. *Estimating Insulin Sensitivity and Beta-Cell Function from the Oral Glucose Tolerance Test: Validation of a New Insulin Sensitivity and Secretion (ISS) Model.* Am J Physiol Endocrinol Metab. (2024) doi: 10.1152/ajpendo.00189.2023.

This R code was translated from the original MATLAB implementation by the authors. The original MATLAB implementation is available in [Figshare](https://figshare.com/articles/journal_contribution/Fitting_the_Insulin_Sensitivity_and_Secretion_ISS_Model_to_OGTTs/23535612).


## Model Overview

The model consists of four ordinary differential equations tracking:

| State Variable | Description |
|---|---|
| G (mg/dl) | Plasma glucose concentration |
| I (µU/ml) | Plasma insulin concentration |
| N5 | Docked insulin granule pool |
| N6 | Reserve insulin granule pool |

Key physiological sub-models include:

- **Hepatic glucose production (HGP):** suppressed by insulin, modulated by SI
- **Glucose rate of appearance (Ra):** piecewise linear absorption profile for the oral glucose load
- **Cytosolic and microdomain calcium:** Hill functions of a metabolic signal driving exocytosis
- **Insulin granule exocytosis:** based on [Chen et al. (2008)](https://pubmed.ncbi.nlm.nih.gov/18515381/), with fast granule states (N1–N4) reduced to quasi-steady state

## Repository Contents

```
├── GI_ode_universal.R               # ODE system definition (4 state variables)
├── OGTT_cost_function_universal.R    # Cost function: solves ODEs, computes WRSS
├── OGTT_fit.R                        # Main fitting script (loop over subjects)
├── Supplemental_Equations.pdf        # Full equation and parameter documentation
└── README.md
```

## Requirements

- **R** ≥ 4.0
- R packages:

```r
install.packages(c("deSolve", "readxl", "writexl"))
```

| Package | Purpose |
|---|---|
| `deSolve` | ODE integration (stiff solver) |
| `readxl` | Read input data from Excel |
| `writexl` | Write fitted results to Excel |

## Input Data Format

The script expects an Excel file (default: `Synthetic_Data.xlsx`) with the following structure:

| subject | ... | G0 | G30 | G60 | G90 | G120 | I0 | I30 | I60 | I90 | I120 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | ... | 85.2 | 142.1 | 130.5 | 110.3 | 95.0 | 6.1 | 45.2 | 38.7 | 25.1 | 12.3 |

- Glucose columns (G0–G120) are in **mg/dl**
- Insulin columns (I0–I120) are in **µU/ml**
- Time points are in **minutes** and are configured at the top of `OGTT_fit.R`
- Column indices (`gstart`, `istart`) should be adjusted to match your file layout
- Missing values (`NA`) are handled automatically; subjects with no valid data are skipped

## Usage

1. Place your data file in the working directory (or update the path in `OGTT_fit.R`).

2. Adjust configuration in `OGTT_fit.R` as needed:
   - `tdata` — your OGTT measurement time points (default: `c(0, 30, 60, 90, 120)`)
   - `gstart`, `istart` — column indices for glucose and insulin in your Excel file
   - `minmod_units` — set to `TRUE` to convert SI to standard MINMOD units (×10⁴/1440)

3. Run the fitting:

```r
source("OGTT_fit.R")
```

4. Results are written to `Synthetic_Fits.xlsx` (or your specified `outfile`) with columns:

| Column | Description |
|---|---|
| `subject` | Subject identifier |
| `sigma` | Fitted beta-cell function (unitless) |
| `si` | Fitted insulin sensitivity (ml/µU·day, or MINMOD units if converted) |
| `WRSS` | Weighted residual sum of squares |
| `true_error` | Unweighted sum of squared errors |
| `RSQ` | R-squared (1 − SSE/SST) |
| `AIC` | Akaike Information Criterion |
| `G0`–`G120` | Simulated glucose at protocol time points |
| `I0`–`I120` | Simulated insulin at protocol time points |

## Optimization Details

The fitting uses R's `optim()` with the **L-BFGS-B** method, which provides:

- Native box constraints on parameters (sigma ∈ [0.01, 10], SI ∈ [0.01, 100])
- Quasi-Newton gradient-based convergence

The ODE system is solved with `deSolve::ode()` using the `lsoda` method (automatic stiff/non-stiff switching) with tight tolerances (`rtol = 1e-8`, `atol = 1e-10`). A 1440-minute burn-in integration is run before each OGTT simulation to ensure the system reaches fasting steady state.

## Customization

**Different OGTT protocols:** Change `tdata` in `OGTT_fit.R` to match your sampling schedule (e.g., `c(0, 10, 20, 30, 60, 90, 120)` for a 7-point OGTT).

**Rate of appearance parameters:** The default `odeparams$a1`, `a2`, `a3`, `t1`, `t2`, `t3` values are population medians. These can be replaced with individualized estimates if available.

**Additional fitted parameters:** The current model fits only σ and SI. To fit additional parameters (e.g., Ra shape parameters), extend `theta0`, `LB`, `UB`, and the parameter assignment block in `OGTT_cost_function_universal.R`.

## Notes on the MATLAB-to-R Translation

This code was translated from MATLAB with the following substitutions:

| MATLAB | R |
|---|---|
| `ode15s` | `deSolve::ode(..., method = "lsoda")` |
| `fminsearch` (Nelder-Mead) | `optim(..., method = "L-BFGS-B")` |
| `interp1` | `approx()` |
| `readtable` / `writetable` | `readxl::read_excel` / `writexl::write_xlsx` |

The optimizer was changed from Nelder-Mead to L-BFGS-B and ODE tolerances were tightened to improve reproducibility across platforms. Results should agree with the MATLAB implementation to 2–4 significant figures on well-identified subjects; small differences are expected due to solver internals.

## Generative AI Disclosure

The R code in this repository was translated from the original MATLAB implementation with the assistance of Claude (Anthropic), a large language model. Specifically, Claude was used to:

- Convert MATLAB syntax and idioms to equivalent R code
- Select appropriate R packages and solver methods (e.g., `deSolve::ode` with `lsoda` as a substitute for MATLAB's `ode15s`)
- Add explanatory comments throughout, drawing on the supplemental equations document
- Draft this README

The model equations, parameter values, algorithmic structure, and scientific content are unchanged from the original MATLAB code and the published reference. All AI-generated code was reviewed and tested by the authors. The optimizer was changed from Nelder-Mead to L-BFGS-B and ODE tolerances were tightened on the recommendation of the AI assistant to improve numerical reproducibility across platforms.

