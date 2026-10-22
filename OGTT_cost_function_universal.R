#' OGTT_cost_function_universal: Cost function for fitting the ISS model to OGTT data
#'
#' Solves the ODE system and computes the weighted residual sum of squares (WRSS)
#' between model predictions and observed glucose/insulin data. Used as the
#' objective function for parameter optimization.
#'
#' The solver first runs a 1440-minute (1-day) burn-in to reach steady state,
#' then integrates over the OGTT time span.
#'
#' @param theta       Parameter vector to optimize: c(sigma, si)
#'                      sigma = beta-cell function scaling factor (unitless)
#'                      si    = peripheral insulin sensitivity (ml/µU·day)
#' @param data        List with fields:
#'                      t       = time points of observations (min)
#'                      G       = observed glucose values (mg/dl)
#'                      I       = observed insulin values (µU/ml)
#'                      outtime = protocol time points for output (may differ
#'                                from t if data points are missing)
#' @param costparams  List with fields:
#'                      init      = initial conditions c(G0, I0, N5_0, N6_0)
#'                      LB        = lower bounds for theta
#'                      UB        = upper bounds for theta
#'                      weights_G = weights for glucose residuals (typically 1/SD)
#'                      weights_I = weights for insulin residuals (typically 1/SD)
#'                      doPlot    = logical, whether to plot during optimization
#' @param odeparams   List of ODE parameters (passed to GI_ode_universal)
#'
#' @return A list with:
#'   S        = weighted residual sum of squares (WRSS)
#'   true_err = unweighted sum of squared errors
#'   Gsim     = simulated glucose at protocol time points
#'   Isim     = simulated insulin at protocol time points

library(deSolve)

OGTT_cost_function_universal <- function(theta, data, costparams, odeparams) {
  
  # ------------------------------------------------------------------
  # Unpack the initial conditions and define the integration time span
  # ------------------------------------------------------------------
  init  <- costparams$init
  tspan <- c(data$t[1], data$t[length(data$t)])
  
  # ------------------------------------------------------------------
  # Enforce parameter bounds (soft check for non-box-constrained solvers)
  # ------------------------------------------------------------------
  # When using L-BFGS-B, bounds are enforced by the optimizer itself,
  # but we keep this check for safety if called outside optimization.
  if (any(theta < costparams$LB) || any(theta > costparams$UB)) {
    return(list(S = 1e6, true_err = 1e6, Gsim = NULL, Isim = NULL))
  }
  
  # ------------------------------------------------------------------
  # Set the two fitted parameters in the ODE parameter list
  # ------------------------------------------------------------------
  odeparams$sigma <- theta[1]   # Beta-cell function (unitless)
  odeparams$si    <- theta[2]   # Insulin sensitivity (ml/µU·day)
  
  # ------------------------------------------------------------------
  # Burn-in integration: run 1440 min (1 day) to reach steady state
  # ------------------------------------------------------------------
  # This ensures the initial conditions for the OGTT simulation are
  
  # at physiological equilibrium before the glucose challenge begins.
  burnin_times <- seq(0, 1440, length.out = 500)
  burnin_sol <- ode(
    y      = init,
    times  = burnin_times,
    func   = GI_ode_universal,
    parms  = odeparams,
    method = "lsoda",               # Stiff/non-stiff auto-switching solver
    rtol   = 1e-8,                  # Tight tolerance for cross-platform agreement
    atol   = 1e-10
  )
  # Use the final state of the burn-in as the initial condition for the OGTT
  init <- as.numeric(burnin_sol[nrow(burnin_sol), -1])  # Drop the time column
  
  # ------------------------------------------------------------------
  # Main OGTT integration over the data time span
  # ------------------------------------------------------------------
  # Create a fine time grid that includes all data time points
  fine_times <- sort(unique(c(seq(tspan[1], tspan[2], by = 0.5), data$t, data$outtime)))
  sol <- ode(
    y      = init,
    times  = fine_times,
    func   = GI_ode_universal,
    parms  = odeparams,
    method = "lsoda",               # Stiff/non-stiff auto-switching solver
    rtol   = 1e-8,                  # Tight tolerance for cross-platform agreement
    atol   = 1e-10
  )
  
  t_sol <- sol[, 1]   # Time column
  G_sol <- sol[, 2]   # Glucose (state 1)
  I_sol <- sol[, 3]   # Insulin (state 2)
  
  # ------------------------------------------------------------------
  # Optional plotting of model fit vs. data
  # ------------------------------------------------------------------
  if (costparams$doPlot) {
    plot(data$t, data$G, pch = 4, cex = 1.5, col = "blue",
         xlab = "Time (min)", ylab = "Concentration",
         ylim = range(c(data$G, data$I, G_sol, I_sol), na.rm = TRUE),
         main = "OGTT Model Fit")
    points(data$t, data$I, pch = 4, cex = 1.5, col = "red")
    lines(t_sol, G_sol, col = "blue")
    lines(t_sol, I_sol, col = "red")
    legend("topright", legend = c("G data", "I data", "G model", "I model"),
           col = c("blue", "red", "blue", "red"),
           pch = c(4, 4, NA, NA), lty = c(NA, NA, 1, 1))
  }
  
  # ------------------------------------------------------------------
  # Interpolate model predictions at observed data time points
  # ------------------------------------------------------------------
  Gpred <- approx(t_sol, G_sol, xout = data$t)$y
  Ipred <- approx(t_sol, I_sol, xout = data$t)$y
  
  # ------------------------------------------------------------------
  # Interpolate model predictions at protocol output time points
  # (these may differ from data$t if some measurements are missing)
  # ------------------------------------------------------------------
  OGTTt <- data$outtime
  Gsim  <- approx(t_sol, G_sol, xout = OGTTt)$y
  Isim  <- approx(t_sol, I_sol, xout = OGTTt)$y
  
  # ------------------------------------------------------------------
  # Compute cost: weighted residual sum of squares (WRSS)
  # ------------------------------------------------------------------
  # WRSS uses per-time-point weights (typically inverse of population SD)
  # to normalize glucose and insulin contributions to the cost.
  S <- sum(costparams$weights_G * (Gpred - data$G)^2 +
             costparams$weights_I * (Ipred - data$I)^2)
  
  # Unweighted sum of squared errors (for diagnostics / R-squared)
  true_err <- sum((Gpred - data$G)^2 + (Ipred - data$I)^2)
  
  return(list(S = S, true_err = true_err, Gsim = Gsim, Isim = Isim))
}