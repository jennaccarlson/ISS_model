#' GI_ode_universal: ODE system for the Insulin Sensitivity and Secretion (ISS) model
#'
#' This function defines the system of ordinary differential equations for modeling
#' glucose-insulin dynamics during an Oral Glucose Tolerance Test (OGTT).
#' Based on the longitudinal diabetes progression model from Ha & Sherman (2020),
#' restricted for OGTT fitting.
#'
#' The system has 4 state variables:
#'   y[1] = G  : Plasma glucose concentration (mg/dl)
#'   y[2] = I  : Plasma insulin concentration (µU/ml)
#'   y[3] = N5 : Insulin granule pool (docked vesicles)
#'   y[4] = N6 : Insulin granule pool (reserve vesicles)
#'
#' @param t     Current time (minutes)
#' @param y     State vector: c(G, I, N5, N6)
#' @param p     List of parameters (must include: si, sigma, r20, hepasi,
#'              gamma_bar, gamma, a1, a2, a3, t1, t2, t3)
#' @return A list containing the vector of derivatives dydt

GI_ode_universal <- function(t, y, p) {
  
  # -----------------------------------------------------------------------
  # Fixed physiological parameters
  # -----------------------------------------------------------------------
  # Eg0: Glucose effectiveness (also called SG), rate of glucose disposal
  #      independent of insulin (min^-1)
  # k:   Insulin clearance rate (min^-1)
  # BV:  Volume of distribution for insulin (mL)
  Eg0 <- 0.0118
  k   <- 0.4861
  BV  <- 7200
  
  # b: Beta-cell mass (mg). Held fixed during OGTT fitting.
  b <- 1553.6
  
  # -----------------------------------------------------------------------
  # Metabolic rate M: Hill function of glucose
  # -----------------------------------------------------------------------
  # M represents the metabolic signal from glucose that drives calcium entry.
  # Modeled as a Hill function with max Mmax, half-activation alpha_M,
  # and Hill coefficient kM.
  Mmax    <- 1
  alpha_M <- 150
  kM      <- 2
  M <- Mmax * y[1]^kM / (alpha_M^kM + y[1]^kM)
  
  # -----------------------------------------------------------------------
  # Rate of appearance of glucose (Ra) during the OGTT
  # -----------------------------------------------------------------------
  # The OGTT glucose flux is modeled as a piecewise linear function (Eq. A1
  # in the supplement). Glucose enters the plasma in three linear segments:
  #   Segment 1: t in (0, t1]     — rising from 0 to a1
  #   Segment 2: t in (t1, t2]    — declining from a1 to a2
  #   Segment 3: t in (t2, t3]    — declining from a2 to a3 (~0)
  # Each segment is a linear interpolation between its endpoint flux values.
  OGTT_rate <- ((t > 0) - (t > p$t1)) * t * p$a1 / p$t1 +
    ((t > p$t1) - (t > p$t2)) * ((t - p$t2) * (p$a2 - p$a1) / (p$t2 - p$t1) + p$a2) +
    ((t > p$t2) - (t > p$t3)) * (t - p$t3) * (p$a3 - p$a2) / (p$t3 - p$t2)
  
  # -----------------------------------------------------------------------
  # Hepatic Glucose Production (HGP) — Eqs. A2, A3, A4
  # -----------------------------------------------------------------------
  # HGP is suppressed by insulin (I) and modulated by peripheral insulin
  # sensitivity (si). As insulin rises post-meal, HGP is suppressed.
  #
  # hepa_max (Eq. A3): maximum HGP suppression capacity, decreasing function of si
  # alpha_HGP (Eq. A4): half-activation for HGP suppression, also decreasing in si
  # HGP (Eq. A2): sigmoidal suppression of HGP by insulin, with basal rate HGP_b
  
  # Parameters for hepa_max (Eq. A3)
  hepa_bar <- 15.443
  hepa_k   <- 0.27
  hepa_b   <- -3.54277
  hepa_max <- hepa_bar / (hepa_k + p$si) + hepa_b
  
  # Parameters for alpha_HGP (Eq. A4)
  alpha_max <- 6
  alpha_k   <- 0.4
  alpha_b   <- -0.5
  alpha_HGP <- alpha_max / (alpha_k + p$si) + alpha_b
  
  # HGP equation (Eq. A2): suppressed by insulin, with basal production rate
  HGP_b <- 0.104166
  HGP   <- hepa_max / (alpha_HGP + y[2] * p$hepasi) + HGP_b
  
  # -----------------------------------------------------------------------
  # Glucose Amplifying Factor GF — Eq. A7
  # -----------------------------------------------------------------------
  # GF captures the glucose-dependent amplification of insulin secretion,
  # implicitly including the incretin effect. It is a Hill function of
  # glucose shifted by GF_sh.
  GF_bar  <- 4.45
  kGF     <- 16
  alpha_GF <- 260
  shGF    <- -89
  GF_b    <- 1.78
  GF <- GF_bar * (y[1] - shGF)^kGF / (alpha_GF^kGF + (y[1] - shGF)^kGF) + GF_b
  
  # -----------------------------------------------------------------------
  # Bulk Cytosolic Calcium Concentration (Ci) — Eq. A5
  # -----------------------------------------------------------------------
  # Ci is a Hill function of the metabolic signal M plus a contribution
  # from GLP-1 signaling (gamma_bar * gamma). Calcium drives exocytosis.
  ca_bar   <- 2
  kca      <- 4
  alpha_ca <- 0.62
  ca_b     <- 0.07
  ci <- ca_bar * (M + p$gamma_bar * p$gamma)^kca /
    (alpha_ca^kca + (M + p$gamma_bar * p$gamma)^kca) + ca_b
  
  # -----------------------------------------------------------------------
  # Microdomain Ca2+ Concentration (Cmd) — Eq. A6
  # -----------------------------------------------------------------------
  # Cmd represents the local calcium concentration near secretory granules,
  # a steep Hill function of bulk cytosolic calcium ci.
  cmd_factor <- 150
  cmd_b      <- 0.0635
  cik        <- 4
  cialpha    <- 1
  cmd <- cmd_factor * ci^cik / (cialpha^cik + ci^cik) + cmd_b
  
  # -----------------------------------------------------------------------
  # Insulin Granule Exocytosis Model — Eqs. A8–A11
  # -----------------------------------------------------------------------
  # Based on Chen, Wang & Sherman (2008, Biophys J 95:2226-41).
  # Models the docking, priming, and fusion of insulin granules.
  # Fast steps (N1–N4) are set to quasi-steady state; N5 (docked) and
  # N6 (reserve) are the slow dynamic variables.
  
  # Kinetic rate constants for granule trafficking (Table S6)
  k1  <- 20       # Binding rate (µM·sec)^-1
  km1 <- 100      # Unbinding rate (sec^-1)
  r1  <- 0.6      # Mobilization rate (sec^-1)
  rm1 <- 1        # Reverse mobilization (sec^-1)
  rm2 <- 0.001    # Reverse priming (sec^-1)
  
  r30 <- 1.205    # Basal delivery rate to membrane (sec^-1)
  rm3 <- 0.0001   # Reverse delivery rate (sec^-1)
  u1  <- 2000     # Fusion rate from docked to fused (sec^-1)
  u2  <- 3        # Rate from fused to released (sec^-1)
  u3  <- 0.02     # Final release rate (sec^-1)
  Kp2 <- 2.3      # Half-activation for Ca2+-dependent rates (µM)
  
  # r2: Ca2+-dependent refilling rate (Eq. A10)
  r2 <- p$r20 * ci / (ci + Kp2)
  
  # ts: time-scale conversion factor (sec/min)
  ts <- 60
  # unit_con: unit conversion factor for SI*I*G term
  unit_con <- 0.00069444
  
  # r3: glucose- and sigma-amplified delivery rate (Eq. A9)
  # sigma (p$sigma) is fit to data as a measure of beta-cell function
  r3 <- p$sigma * GF * r30 * ci / (ci + Kp2)
  
  # -----------------------------------------------------------------------
  # Quasi-steady-state reduction of fast granule states (N1–N4)
  # -----------------------------------------------------------------------
  # Express N1–N4 as algebraic functions of N5 (y[3]) using the
  # chain of quasi-steady-state relations (Eqs. A8).
  
  # Coefficients for N1 as a function of N5
  N1_C <- km1 / (3 * k1 * cmd + rm1)
  N1_D <- r1  / (3 * k1 * cmd + rm1)
  
  # Coefficients for N2 as a function of N1
  N2_E <- 3 * k1 * cmd / (2 * k1 * cmd + km1)
  N2_F <- 2 * km1      / (2 * k1 * cmd + km1)
  
  # Coefficients for N3 as a function of N2
  N3_L <- 2 * k1 * cmd / (2 * km1 + k1 * cmd)
  N3_N <- 3 * km1      / (2 * km1 + k1 * cmd)
  
  # Chain substitution to get compound coefficients CN1–CN4
  CN4 <- k1 * cmd / (3 * km1 + u1)
  CN3 <- N3_L / (1 - N3_N * CN4)
  CN2 <- N2_E / (1 - N2_F * CN3)
  CN1 <- N1_D / (1 - N1_C * CN2)
  
  # Compute granule pool sizes from N5 (y[3])
  N1 <- CN1 * y[3]
  N2 <- CN2 * N1
  N3 <- CN3 * N2
  N4 <- CN4 * N3
  NF <- u1 * N4 / u2      # Fused pool
  NR <- (u2 / u3) * NF    # Released pool
  
  # -----------------------------------------------------------------------
  # Insulin Secretion Rate (ISR) — Eq. A11
  # -----------------------------------------------------------------------
  # ISR = ts * r * u3 * NR, where r ≈ 9 µU/mg (see supplement p.6-7)
  # The factor 9 converts vesicles/cell/min to µU/mg/min
  ISR <- ts * 9 * (u3 * NR)
  
  # -----------------------------------------------------------------------
  # Differential equations
  # -----------------------------------------------------------------------
  
  # dG/dt: Glucose dynamics (Eq. 1 in supplement)
  # = HGP + glucose appearance - (glucose effectiveness + insulin-mediated uptake) * G
  dGdt <- HGP + OGTT_rate - (Eg0 + unit_con * p$si * y[2]) * y[1]
  
  # dI/dt: Insulin dynamics (Eq. 2 in supplement)
  # = beta-cell mass * ISR / distribution volume - insulin clearance
  dIdt <- b * ISR / BV - k * y[2]
  
  # dN5/dt: Docked granule pool dynamics (part of Eq. A8)
  dN5dt <- ts * (rm1 * CN1 * y[3] - (r1 + rm2) * y[3] + r2 * y[4])
  
  # dN6/dt: Reserve granule pool dynamics (part of Eq. A8)
  dN6dt <- ts * (r3 + rm2 * y[3] - (rm3 + r2) * y[4])
  
  # Return derivatives as a list (required by deSolve)
  return(list(c(dGdt, dIdt, dN5dt, dN6dt)))
}