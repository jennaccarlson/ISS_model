#' OGTT_fit.R: Main script to fit the ISS model to OGTT data
#'
#' This script reads glucose and insulin time-series data from an Excel file,
#' fits the ISS (Insulin Sensitivity and Secretion) model to each subject
#' by optimizing two parameters:
#'   sigma = beta-cell function (unitless scaling factor)
#'   si    = peripheral insulin sensitivity (ml/µU·day)
#'
#' The fitting minimizes a weighted residual sum of squares (WRSS) using
#' Nelder-Mead optimization. Results (fitted parameters, goodness-of-fit
#' statistics, and simulated glucose/insulin curves) are saved to an Excel file.

# =========================================================================
# Required packages
# =========================================================================
library(deSolve)      # ODE solver (ode function with stiff methods)
library(readxl)       # Read Excel files
library(writexl)      # Write Excel files

# Source the model and cost function files
source("GI_ode_universal.R")
source("OGTT_cost_function_universal.R")

# =========================================================================
# Configuration
# =========================================================================

# Set to TRUE to convert SI to standard MINMOD units at the end
minmod_units <- FALSE

# OGTT measurement time points (minutes)
tdata   <- c(0, 30, 60, 90, 120)
outtime <- tdata
numout  <- length(tdata) * 2  # Total output columns: G and I at each time

# Build column names for the output: G0, G30, ..., I0, I30, ...
outnames <- character(numout)
for (i in seq_along(tdata)) {
  outnames[i]                    <- paste0("G", tdata[i])
  outnames[i + length(tdata)]    <- paste0("I", tdata[i])
}
GInames <- outnames

# =========================================================================
# Data input: column indices for glucose and insulin in the Excel file
# =========================================================================
# Adjust these indices to match your Excel file layout
gstart <- 4
gend   <- gstart + length(tdata) - 1

istart <- gend + 1
iend   <- istart + length(tdata) - 1

# Read data and specify output file
DATA    <- read_excel("Synthetic_Data.xlsx")
outfile <- "Synthetic_Fits.xlsx"

# Extract glucose and insulin matrices (subjects x time points)
GDATA <- as.matrix(DATA[, gstart:gend])
IDATA <- as.matrix(DATA[, istart:iend])

# =========================================================================
# Population-level parameters for the rate of appearance (Ra)
# =========================================================================
# These are median values from UPMC data, used for all subjects.
# They define the piecewise linear glucose absorption profile (Eq. A1).
odeparams <- list(
  a1 = 5.99,       # Glucose flux at t1 (mg/min, scaled)
  a2 = 2.14,       # Glucose flux at t2
  a3 = 0.0013,     # Glucose flux at t3 (~0)
  t1 = 15.60,      # First knot time (min)
  t2 = 137.2,      # Second knot time (min)
  t3 = 258.3,      # Third knot time (min)
  r20 = 0.006,     # Basal refilling rate for exocytosis (sec^-1)
  hepasi = 1,      # Hepatic insulin sensitivity multiplier
  gamma_bar = 1,   # GLP-1 signaling scaling factor
  gamma = 2.5e-6 * 17500  # GLP-1 signaling parameter (unitless)
)

# =========================================================================
# Initial conditions for the ODE system
# =========================================================================
# y = c(G, I, N5, N6) at fasting steady state
init <- c(80.1842, 5.8462, 60.9341, 443.7764)

# =========================================================================
# Initial guesses for the fitted parameters
# =========================================================================
sigma <- 1       # Beta-cell function (unitless)
si    <- 0.8     # Insulin sensitivity (ml/µU·day)
theta0 <- c(sigma, si)

# =========================================================================
# Parameter bounds for optimization
# =========================================================================
sigmalb <- 0.01;  sigmaub <- 10
silb    <- 0.01;  siub    <- 100

LB <- c(sigmalb, silb)
UB <- c(sigmaub, siub)

# =========================================================================
# Cost function configuration
# =========================================================================
costparams <- list(
  init  = init,
  LB    = LB,
  UB    = UB,
  doPlot = FALSE
)

# Compute population-level standard deviations for weighting
# (inverse-variance weighting normalizes glucose and insulin contributions)
var_G <- apply(GDATA, 2, sd, na.rm = TRUE)
var_I <- apply(IDATA, 2, sd, na.rm = TRUE)

# =========================================================================
# Fitting configuration
# =========================================================================
numpar   <- 2                # Number of fitted parameters (sigma, si)
numpts   <- length(tdata)    # Number of OGTT time points
numextra <- 4                # Extra output columns: WRSS, true_error, RSQ, AIC
numvals  <- length(theta0) + numpts * 2 + numextra

# Pre-allocate results matrices
RESULTS <- matrix(0, nrow = nrow(DATA), ncol = numvals)
PARAMS  <- matrix(0, nrow = nrow(DATA), ncol = numpar)

# =========================================================================
# Main fitting loop: iterate over each subject
# =========================================================================
for (patientID in 1:nrow(DATA)) {
  
  # Set up data for this subject
  data_subj <- list(
    t       = tdata,
    G       = as.numeric(GDATA[patientID, ]),
    I       = as.numeric(IDATA[patientID, ]),
    outtime = outtime   # Protocol time points (even if data is missing)
  )
  
  # Weights: inverse of population SD at each time point
  costparams$weights_G <- 1 / var_G
  costparams$weights_I <- 1 / var_I
  
  # ------------------------------------------------------------------
  # Handle missing data: remove time points with NA in either G or I
  # ------------------------------------------------------------------
  missing <- is.na(data_subj$G) | is.na(data_subj$I)
  
  if (all(missing)) {
    message("No data for #", patientID, " - skipping")
    RESULTS[patientID, ] <- NA
    next
  }
  
  data_subj$t <- data_subj$t[!missing]
  data_subj$G <- data_subj$G[!missing]
  data_subj$I <- data_subj$I[!missing]
  costparams$weights_G <- costparams$weights_G[!missing]
  costparams$weights_I <- costparams$weights_I[!missing]
  
  message("Patient ID: ", patientID)
  
  # ------------------------------------------------------------------
  # Optimization using L-BFGS-B with native box constraints
  # ------------------------------------------------------------------
  # L-BFGS-B enforces parameter bounds directly (no penalty needed),
  # and its quasi-Newton steps are more deterministic than Nelder-Mead,
  # improving cross-platform reproducibility with MATLAB's fminsearch.
  costfun <- function(theta) {
    result <- OGTT_cost_function_universal(theta, data_subj, costparams, odeparams)
    return(result$S)
  }
  
  opt_result <- optim(
    par     = theta0,
    fn      = costfun,
    method  = "L-BFGS-B",
    lower   = LB,                  # Box constraints enforced by the algorithm
    upper   = UB,
    control = list(
      factr   = 1e7,               # Convergence tolerance on function value
      #   (relative, in units of machine epsilon)
      pgtol   = 1e-6,              # Convergence tolerance on projected gradient
      maxit   = 1e4                # Maximum iterations
    )
  )
  
  paramin <- opt_result$par
  ERRmin  <- opt_result$value
  
  message("Estimated parameters sigma si: ", paste(round(paramin, 4), collapse = "  "))
  
  # ------------------------------------------------------------------
  # Evaluate the fitted solution to get all outputs
  # ------------------------------------------------------------------
  fit_result <- OGTT_cost_function_universal(paramin, data_subj, costparams, odeparams)
  WRSS       <- fit_result$S
  true_error <- fit_result$true_err
  Gsim       <- fit_result$Gsim
  Isim       <- fit_result$Isim
  
  # ------------------------------------------------------------------
  # Optionally convert SI to standard MINMOD units
  # ------------------------------------------------------------------
  # MINMOD units: multiply by 10^4 / 1440 to convert from ml/µU·day
  # to ml/µU·min (with the 10^4 convenience factor)
  if (minmod_units) {
    paramin[2] <- paramin[2] * 10000 / 1440
  }
  
  PARAMS[patientID, ] <- paramin
  
  # ------------------------------------------------------------------
  # Goodness-of-fit statistics
  # ------------------------------------------------------------------
  # R-squared: 1 - SSE/SST, where SST is total sum of squares of the data
  SS  <- sum(data_subj$G^2, na.rm = TRUE) + sum(data_subj$I^2, na.rm = TRUE)
  RSQ <- 1 - true_error / SS
  
  # Akaike Information Criterion (AIC)
  AIC_val <- 2 * numpar + numpts * log(numpts * WRSS)
  
  # ------------------------------------------------------------------
  # Store results for this subject
  # ------------------------------------------------------------------
  RESULTS[patientID, ] <- c(paramin, WRSS, true_error, RSQ, AIC_val, Gsim, Isim)
  
  message("Final results (WRSS, true_error, RSQ, AIC): ",
          paste(round(c(WRSS, true_error, RSQ, AIC_val), 4), collapse = "  "))
  message("")
}

# =========================================================================
# Format and save results
# =========================================================================

# Convert results matrix to data frame with named columns
colnames(RESULTS) <- c("sigma", "si", "WRSS", "true_error", "RSQ", "AIC", GInames)
RESULTS <- as.data.frame(RESULTS)

# Prepend subject identifiers from the input data
RESULTS <- cbind(subject = DATA$subject, RESULTS)

# Write to Excel
write_xlsx(RESULTS, outfile)
message("Results written to: ", outfile)