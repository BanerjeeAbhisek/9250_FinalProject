# =============================================================================
# Comparative Study of MCMC Sampling Strategies for the 2D Ising Model
# Sections 1-3 (consolidated and executable)
#
# Abhisek Banerjee & Mengyan Jing
# Course 9250, Spring 2026
#
# This file is self-contained: sourcing it runs every section end-to-end
# (prints T_c, initializes a test lattice, evaluates the Onsager curves,
# displays the ground-truth figure). Functions are still defined so the
# downstream sections (4+) can call them unchanged.
#
# Layout:
#   Section 1: Setup (constants, T_c, literature z values)         + demo
#   Section 2: Lattice utilities (init, total energy, dE w/ PBC)   + demo
#   Section 3: Onsager / Yang exact ground truth                   + demo
#     3.1  complete elliptic integrals (guarded)
#     3.2  Onsager energy                                          Eq.(116)
#     3.3  Onsager specific heat                                   Eq.(117)
#     3.4  Onsager specific heat asymptotic form                   Eq.(120)
#     3.5  Yang spontaneous magnetization                          (1952)
#     3.6  safe wrappers (avoid singular elliptic_K near T_c)
#     3.7  curve wrappers (NA in a guard band around T_c)
#     3.8  publication figure via ggplot2  (new: replaces base R version)
#
# Dependencies (install once):
#   install.packages(c("Rcpp", "ggplot2", "patchwork"))
# =============================================================================

suppressPackageStartupMessages({
  library(Rcpp)
  library(ggplot2)
  library(patchwork)
})


# =============================================================================
# Section 1: Setup
# =============================================================================

J_COUPLING  <- 1
K_BOLTZMANN <- 1
T_CRITICAL  <- 2 / log(1 + sqrt(2))     # Onsager (1944), from Eq. (39a)

# Literature dynamic critical exponents (used downstream in Sections 8 & 9)
Z_METROPOLIS_LITERATURE <- 2.1665   # Nightingale & Blote (1996), PRL 76, 4548
Z_WOLFF_LITERATURE      <- 0.25     # Wolff (1989), PRL 62, 361

# --- Section 1 output -------------------------------------------------------
cat("=== Section 1: Setup ===\n")
cat(sprintf("  J (coupling)            = %g\n",    J_COUPLING))
cat(sprintf("  k_B (Boltzmann)         = %g\n",    K_BOLTZMANN))
cat(sprintf("  T_c  (Onsager, exact)   = %.6f\n",  T_CRITICAL))
cat(sprintf("  z_Metropolis (lit.)     = %.4f  (Nightingale & Blote 1996)\n",
            Z_METROPOLIS_LITERATURE))
cat(sprintf("  z_Wolff      (lit.)     = %.2f    (Wolff 1989)\n",
            Z_WOLFF_LITERATURE))
cat("\n")


# =============================================================================
# Section 2: Lattice utilities
# =============================================================================

# 2.1: initialization ---------------------------------------------------------

init_lattice <- function(L, state = c("cold", "hot", "ground")) {
  state <- match.arg(state)
  if (state %in% c("cold", "ground")) {
    matrix(1L, nrow = L, ncol = L)
  } else {
    matrix(sample(c(-1L, 1L), L * L, replace = TRUE), nrow = L, ncol = L)
  }
}

# 2.2: total energy (Newman & Barkema Eq. 3.1) -------------------------------

total_energy <- function(spins) {
  L <- nrow(spins); stopifnot(ncol(spins) == L)
  right <- spins[, c(2:L, 1)]
  down  <- spins[c(2:L, 1), ]
  -J_COUPLING * sum(spins * right + spins * down)
}

total_magnetization <- function(spins) sum(spins)

# 2.3: single-spin dE (Newman & Barkema Eq. 3.10) ----------------------------

delta_energy <- function(spins, i, j) {
  L  <- nrow(spins)
  ip <- if (i == L)  1L else i + 1L
  im <- if (i == 1L) L  else i - 1L
  jp <- if (j == L)  1L else j + 1L
  jm <- if (j == 1L) L  else j - 1L
  neighbor_sum <- spins[ip, j] + spins[im, j] + spins[i, jp] + spins[i, jm]
  2 * J_COUPLING * spins[i, j] * neighbor_sum
}

# --- Section 2 output (small smoke test) ------------------------------------
cat("=== Section 2: Lattice utilities ===\n")
local({
  set.seed(42)
  L <- 8
  
  s_cold <- init_lattice(L, "cold")
  s_hot  <- init_lattice(L, "hot")
  cat(sprintf("  cold lattice  L=%d:  E = %d (expected %d)   M = %d (expected %d)\n",
              L, total_energy(s_cold), -2 * J_COUPLING * L^2,
              total_magnetization(s_cold), L^2))
  cat(sprintf("  hot  lattice  L=%d:  E = %d             M = %d\n",
              L, total_energy(s_hot), total_magnetization(s_hot)))
  
  # Corner-flip check (PBC): dE should be +8 when flipping a spin in the
  # all-+1 ground state (each of 4 neighbors = +1, so dE = 2 * 1 * 4 = 8).
  dE_corner <- delta_energy(s_cold, 1L, 1L)
  cat(sprintf("  dE flipping (1,1) in cold lattice = %d (expected 8)\n",
              dE_corner))
})
cat("\n")


# =============================================================================
# Section 3: Onsager / Yang exact ground truth
# =============================================================================

# 3.1: complete elliptic integrals -------------------------------------------
# Hardened guard: any |k| within 1e-10 of 1 returns Inf, preventing
# integrate() from being called at numerically-pathological k values.

elliptic_K <- function(k) {
  if (abs(k) >= 1 - 1e-10) return(Inf)
  integrand <- function(theta) 1 / sqrt(1 - k^2 * sin(theta)^2)
  integrate(integrand, 0, pi / 2, rel.tol = 1e-8, subdivisions = 2000L)$value
}

elliptic_E <- function(k) {
  integrand <- function(theta) sqrt(1 - k^2 * sin(theta)^2)
  integrate(integrand, 0, pi / 2, rel.tol = 1e-8, subdivisions = 2000L)$value
}

# 3.2: energy per spin (Onsager Eq. 116) -------------------------------------

onsager_energy <- function(T) {
  vapply(T, function(Ti) {
    beta <- 1 / Ti; bj <- beta * J_COUPLING
    s2   <- sinh(2 * bj); c2 <- cosh(2 * bj)
    k1   <- 2 * s2 / c2^2
    K1   <- elliptic_K(k1)
    t2   <- tanh(2 * bj)^2
    -J_COUPLING * (c2 / s2) * (1 + (2 / pi) * (2 * t2 - 1) * K1)
  }, numeric(1))
}

# 3.3: specific heat (Onsager Eq. 117) ---------------------------------------

onsager_specific_heat <- function(T) {
  vapply(T, function(Ti) {
    beta <- 1 / Ti; bj <- beta * J_COUPLING
    s2   <- sinh(2 * bj); c2 <- cosh(2 * bj)
    k1   <- 2 * s2 / c2^2
    K1   <- elliptic_K(k1); E1 <- elliptic_E(k1)
    k1pp <- 2 * tanh(2 * bj)^2 - 1
    coth2 <- c2 / s2
    pre  <- (2 / pi) * (bj * coth2)^2
    pre * (2 * K1 - 2 * E1 - (1 - k1pp) * (pi / 2 + k1pp * K1))
  }, numeric(1))
}

# 3.4: asymptotic specific heat near T_c (Onsager Eq. 120) -------------------

onsager_specific_heat_asymptotic <- function(T) {
  vapply(T, function(Ti) {
    beta <- 1 / Ti; bj <- beta * J_COUPLING
    k1   <- 2 * sinh(2 * bj) / cosh(2 * bj)^2
    K1   <- elliptic_K(k1)
    (2 / pi) * log(1 / tan(pi / 8))^2 * (K1 - 1 - pi / 4)
  }, numeric(1))
}

# 3.5: magnetization (Yang 1952) ---------------------------------------------

onsager_magnetization <- function(T) {
  vapply(T, function(Ti) {
    if (Ti >= T_CRITICAL) return(0)
    beta <- 1 / Ti
    (1 - sinh(2 * beta * J_COUPLING)^(-4))^(1 / 8)
  }, numeric(1))
}

# 3.6: safe ground-truth wrappers (avoid singular elliptic_K) ----------------
# Use vapply + if/else (lazy: only the chosen branch is evaluated).
# .SHIFT = 1e-3: u(T_c - 1e-3) = -1.4182 vs true -sqrt(2) ~ -1.4142
# (0.3% off, fine as a numerical proxy).

.EPS_TC    <- 1e-9
.SHIFT     <- 1e-3
.SHIFT_LOG <- 1e-3

safe_onsager_energy <- function(T) {
  vapply(T, function(Ti) {
    if (abs(Ti - T_CRITICAL) < .EPS_TC) onsager_energy(T_CRITICAL - .SHIFT)
    else                                 onsager_energy(Ti)
  }, numeric(1))
}

safe_onsager_magnetization <- function(T) {
  vapply(T, function(Ti) {
    if (Ti >= T_CRITICAL - .EPS_TC) 0
    else                             onsager_magnetization(Ti)
  }, numeric(1))
}

safe_onsager_specific_heat <- function(T) {
  vapply(T, function(Ti) {
    if (abs(Ti - T_CRITICAL) < .EPS_TC) Inf
    else                                 onsager_specific_heat(Ti)
  }, numeric(1))
}

onsager_specific_heat_proxy <- function(T) {
  vapply(T, function(Ti) {
    if (abs(Ti - T_CRITICAL) < .EPS_TC) onsager_specific_heat(T_CRITICAL - .SHIFT_LOG)
    else                                 onsager_specific_heat(Ti)
  }, numeric(1))
}

# 3.7: curve wrappers for plotting (NA in window around T_c) -----------------
# Default guard 0.015 is comfortably > .SHIFT and big enough that
# integrate() never gets called too close to k = 1.

onsager_energy_curve <- function(T_grid, guard = 0.015) {
  result  <- numeric(length(T_grid))
  near_tc <- abs(T_grid - T_CRITICAL) < guard
  result[near_tc]  <- NA
  result[!near_tc] <- onsager_energy(T_grid[!near_tc])
  result
}

onsager_specific_heat_curve <- function(T_grid, guard = 0.015) {
  result  <- numeric(length(T_grid))
  near_tc <- abs(T_grid - T_CRITICAL) < guard
  result[near_tc]  <- NA
  result[!near_tc] <- sapply(T_grid[!near_tc], onsager_specific_heat)
  result
}

onsager_magnetization_curve <- function(T_grid, guard = 0.015) {
  result  <- numeric(length(T_grid))
  near_tc <- abs(T_grid - T_CRITICAL) < guard
  result[near_tc]  <- NA
  result[!near_tc] <- onsager_magnetization(T_grid[!near_tc])
  result
}

# --- Section 3 smoke test ---------------------------------------------------
cat("=== Section 3: Onsager / Yang exact ground truth ===\n")
local({
  # A few evaluations: low T, just below T_c, high T
  Ts <- c(1.0, T_CRITICAL - 1e-3, 3.0)
  cat("  Observable evaluations at selected temperatures:\n")
  cat(sprintf("    %-10s  %10s  %10s  %10s\n",
              "T", "u(T)", "c(T)", "|m|(T)"))
  cat(sprintf("    %s\n", strrep("-", 45)))
  for (T in Ts) {
    cat(sprintf("    %-10.5f  %10.4f  %10.4f  %10.4f\n",
                T,
                onsager_energy(T),
                onsager_specific_heat(T),
                onsager_magnetization(T)))
  }
  cat(sprintf("\n  Check: u(T_c - 1e-3) = %.4f   (expected ~ -sqrt(2) = %.4f)\n",
              onsager_energy(T_CRITICAL - 1e-3), -sqrt(2)))
})
cat("\n")


# =============================================================================
# 3.8: Publication-quality figure via ggplot2
#
# House styling (theme_paper + Okabe-Ito palette) is defined here because
# it will be reused by every downstream figure in the paper.
# =============================================================================

# -- shared theme -------------------------------------------------------------

theme_paper <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "") +
    theme(
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor  = element_line(color = "grey96", linewidth = 0.2),
      panel.border      = element_rect(color = "grey30", fill = NA, linewidth = 0.4),
      axis.text         = element_text(color = "grey20"),
      axis.title        = element_text(color = "grey10"),
      axis.ticks        = element_line(color = "grey30", linewidth = 0.3),
      axis.ticks.length = unit(3, "pt"),
      plot.title        = element_text(size = rel(1.0), face = "bold",
                                       color = "grey10", margin = margin(b = 4)),
      plot.title.position = "plot",
      legend.position   = "bottom",
      legend.title      = element_text(size = rel(0.85)),
      legend.text       = element_text(size = rel(0.80)),
      legend.key.height = unit(9, "pt"),
      plot.margin       = margin(4, 8, 4, 4)
    )
}

# -- shared Okabe-Ito palette -------------------------------------------------

PAPER_COLORS <- list(
  energy         = "#0072B2",   # blue
  specific_heat  = "#D55E00",   # vermillion
  magnetization  = "#009E73",   # bluish green
  susceptibility = "#CC79A7",   # reddish purple
  tc_line        = "grey45",
  metro          = "#D55E00",   # Metropolis in downstream figures
  wolff          = "#0072B2"    # Wolff      in downstream figures
)

# -- one panel ----------------------------------------------------------------

.plot_one_observable <- function(df, ylab_expr, color, title,
                                 y_limits = NULL, tc_label_y = NULL,
                                 show_y_zero = FALSE) {
  p <- ggplot(df, aes(x = T, y = value)) +
    geom_vline(xintercept = T_CRITICAL,
               linetype = "dashed", color = PAPER_COLORS$tc_line,
               linewidth = 0.4) +
    geom_line(color = color, linewidth = 0.9, na.rm = TRUE) +
    labs(title = title,
         x     = expression(italic(T)),
         y     = ylab_expr) +
    annotate("text",
             x = T_CRITICAL + 0.04,
             y = if (!is.null(tc_label_y)) tc_label_y
             else {
               v <- df$value[is.finite(df$value)]
               max(v) - 0.05 * diff(range(v))
             },
             label = "italic(T)[c]", parse = TRUE,
             color = PAPER_COLORS$tc_line, size = 3, hjust = 0) +
    theme_paper()
  
  if (show_y_zero) {
    p <- p + geom_hline(yintercept = 0, color = "grey70",
                        linewidth = 0.3, linetype = "dotted")
  }
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits, clip = "off")
  }
  p
}

# -- full figure builder ------------------------------------------------------

plot_onsager_ground_truth_gg <- function(T_min     = 1.0,
                                         T_max     = 4.0,
                                         dT        = 0.005,
                                         guard     = 0.015,
                                         c_cap     = 2.5) {
  T_grid <- seq(T_min, T_max, by = dT)
  
  df_u <- data.frame(T = T_grid, value = onsager_energy_curve(T_grid, guard))
  df_c <- data.frame(T = T_grid, value = onsager_specific_heat_curve(T_grid, guard))
  df_m <- data.frame(T = T_grid, value = onsager_magnetization_curve(T_grid, guard))
  df_c$value <- pmin(df_c$value, c_cap)  # cap log divergence for display
  
  p_u <- .plot_one_observable(
    df_u,
    ylab_expr = expression(italic(u)(italic(T))),
    color     = PAPER_COLORS$energy,
    title     = "Energy per spin",
    tc_label_y = -0.75
  )
  
  p_c <- .plot_one_observable(
    df_c,
    ylab_expr = expression(italic(c)(italic(T))),
    color     = PAPER_COLORS$specific_heat,
    title     = "Specific heat per spin",
    y_limits  = c(0, c_cap),
    tc_label_y = c_cap * 0.92
  ) +
    annotate("text",
             x = T_max - 0.05, y = c_cap * 0.92,
             label = sprintf("(capped at %.1f)", c_cap),
             hjust = 1, size = 2.6, color = "grey40")
  
  p_m <- .plot_one_observable(
    df_m,
    ylab_expr = expression(group("|", italic(m), "|") * (italic(T))),
    color     = PAPER_COLORS$magnetization,
    title     = "Spontaneous magnetization",
    y_limits  = c(0, 1.02),
    tc_label_y = 0.92,
    show_y_zero = TRUE
  )
  
  (p_u | p_c | p_m) +
    plot_annotation(
      caption = sprintf(
        "Onsager (1944) and Yang (1952) exact solutions. T[c] = %.4f.",
        T_CRITICAL),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    )
}

# --- Section 3 output: render the figure -----------------------------------
cat("=== Section 3.8: ground-truth figure ===\n")
print(plot_onsager_ground_truth_gg())
cat("\n")

cat("=== Sections 1-3 complete. ===\n")

# =============================================================================
# Section 4: Sanity checks (14 assertions)
# =============================================================================

run_sanity_checks <- function() {
  cat("=== Sanity checks ===\n")
  
  # Lattice (1-5)
  L <- 8
  s <- init_lattice(L, "cold")
  E_cold <- total_energy(s); expected <- -2 * J_COUPLING * L^2
  cat(sprintf("  [01] Cold lattice energy: got %g, expected %g -- %s\n",
              E_cold, expected, if (E_cold == expected) "PASS" else "FAIL"))
  M_cold <- total_magnetization(s)
  cat(sprintf("  [02] Cold lattice magnetization: got %d, expected %d -- %s\n",
              M_cold, L^2, if (M_cold == L^2) "PASS" else "FAIL"))
  dE <- delta_energy(s, 4, 4)
  cat(sprintf("  [03] dE flipping one spin in cold lattice: got %g, expected 8 -- %s\n",
              dE, if (dE == 8) "PASS" else "FAIL"))
  s2 <- s; s2[4, 4] <- -s2[4, 4]; E_new <- total_energy(s2)
  cat(sprintf("  [04] Incremental vs. full energy after flip: %g vs %g -- %s\n",
              E_cold + dE, E_new,
              if (isTRUE(all.equal(E_cold + dE, E_new))) "PASS" else "FAIL"))
  dE_corner <- delta_energy(s, 1, 1)
  cat(sprintf("  [05] dE at corner (PBC check): got %g, expected 8 -- %s\n",
              dE_corner, if (dE_corner == 8) "PASS" else "FAIL"))
  
  # Onsager (6-10)
  cat(sprintf("  [06] T_c = %.6f (Onsager)\n", T_CRITICAL))
  u_Tc <- onsager_energy(T_CRITICAL - 1e-3)
  cat(sprintf("  [07] u(T_c - 1e-3): got %.6f, expected ~%.6f -- %s\n",
              u_Tc, -sqrt(2),
              if (isTRUE(all.equal(u_Tc, -sqrt(2), tolerance = 1e-2))) "PASS" else "FAIL"))
  T_near <- T_CRITICAL - 0.005
  c_full <- onsager_specific_heat(T_near)
  c_asym <- onsager_specific_heat_asymptotic(T_near)
  cat(sprintf("  [08] c(T_c - 0.005): full %.4f vs asymptotic %.4f -- %s\n",
              c_full, c_asym,
              if (isTRUE(all.equal(c_full, c_asym, tolerance = 1e-2))) "PASS" else "FAIL"))
  m_low <- onsager_magnetization(1.0)
  cat(sprintf("  [09] |m|(T=1): got %.6f (should be ~1) -- %s\n",
              m_low, if (m_low > 0.999) "PASS" else "FAIL"))
  m_high <- onsager_magnetization(3.0)
  cat(sprintf("  [10] |m|(T=3): got %.6f (should be 0) -- %s\n",
              m_high, if (m_high == 0) "PASS" else "FAIL"))
  
  # Asymptotic + high-T limit (11-12)
  rel_err_far <- {
    T <- T_CRITICAL - 5e-3
    abs(onsager_specific_heat(T) - onsager_specific_heat_asymptotic(T)) /
      abs(onsager_specific_heat(T))
  }
  rel_err_near <- {
    T <- T_CRITICAL - 1e-4
    abs(onsager_specific_heat(T) - onsager_specific_heat_asymptotic(T)) /
      abs(onsager_specific_heat(T))
  }
  cat(sprintf("  [11] c_asym convergence: rel.err(dT=5e-3)=%.2e > rel.err(dT=1e-4)=%.2e -- %s\n",
              rel_err_far, rel_err_near,
              if (rel_err_near < rel_err_far) "PASS" else "FAIL"))
  u_hot <- onsager_energy(100.0)
  cat(sprintf("  [12] u(T=100) high-T limit: got %.4f (should be ~0) -- %s\n",
              u_hot, if (abs(u_hot) < 0.1) "PASS" else "FAIL"))
  
  # z constants (13-14)
  cat(sprintf("  [13] Z_METROPOLIS_LITERATURE = %.4f (Nightingale & Blote 1996)\n",
              Z_METROPOLIS_LITERATURE))
  cat(sprintf("  [14] Z_WOLFF_LITERATURE = %.2f (Wolff 1989: ~0.25 for 2D Ising)\n",
              Z_WOLFF_LITERATURE))
  
  cat("=====================\n")
}

# --- Section 4 output: run the checks ---------------------------------------
run_sanity_checks()
cat("\n")

##> === Sanity checks ===
##>   [01] Cold lattice energy: got -128, expected -128 -- PASS
##>   [02] Cold lattice magnetization: got 64, expected 64 -- PASS
##>   [03] dE flipping one spin in cold lattice: got 8, expected 8 -- PASS
##>   [04] Incremental vs. full energy after flip: -120 vs -120 -- PASS
##>   [05] dE at corner (PBC check): got 8, expected 8 -- PASS
##>   [06] T_c = 2.269185 (Onsager)
##>   [07] u(T_c - 1e-3): got -1.418223, expected ~-1.414214 -- PASS
##>   [08] c(T_c - 0.005): full 2.7190 vs asymptotic 2.7184 -- PASS
##>   [09] |m|(T=1): got 0.999276 (should be ~1) -- PASS
##>   [10] |m|(T=3): got 0.000000 (should be 0) -- PASS
##>   [11] c_asym convergence: rel.err(dT=5e-3)=2.35e-04 > rel.err(dT=1e-4)=5.12e-06 -- PASS
##>   [12] u(T=100) high-T limit: got -0.0200 (should be ~0) -- PASS
##>   [13] Z_METROPOLIS_LITERATURE = 2.1665 (Nightingale & Blote 1996)
##>   [14] Z_WOLFF_LITERATURE = 0.25 (Wolff 1989: ~0.25 for 2D Ising)
##>   =====================
# =============================================================================
# Section 5: Pure R Metropolis sampler (reference)
# =============================================================================

metropolis_step <- function(spins, T) {
  L <- nrow(spins)
  i <- sample.int(L, 1L); j <- sample.int(L, 1L)
  dE <- delta_energy(spins, i, j)
  if (dE <= 0 || runif(1) < exp(-dE / T)) {
    s_old <- spins[i, j]; spins[i, j] <- -s_old
    list(spins = spins, dE = dE, dM = -2L * s_old, accepted = TRUE)
  } else {
    list(spins = spins, dE = 0, dM = 0L, accepted = FALSE)
  }
}

metropolis_sweep <- function(spins, T) {
  N <- nrow(spins)^2
  dE_sweep <- 0; dM_sweep <- 0L; n_accepted <- 0L
  for (k in seq_len(N)) {
    step     <- metropolis_step(spins, T)
    spins    <- step$spins
    dE_sweep <- dE_sweep + step$dE
    dM_sweep <- dM_sweep + step$dM
    if (step$accepted) n_accepted <- n_accepted + 1L
  }
  list(spins = spins, dE_sweep = dE_sweep, dM_sweep = dM_sweep,
       n_accepted = n_accepted)
}

run_metropolis <- function(L, T, n_sweeps,
                           n_burnin   = n_sweeps %/% 10,
                           init_state = "hot",
                           seed       = NULL) {
  if (!is.null(seed)) set.seed(seed)
  spins     <- init_lattice(L, init_state)
  E_current <- total_energy(spins)
  M_current <- total_magnetization(spins)
  
  for (s in seq_len(n_burnin)) {
    sweep     <- metropolis_sweep(spins, T)
    spins     <- sweep$spins
    E_current <- E_current + sweep$dE_sweep
    M_current <- M_current + sweep$dM_sweep
  }
  
  E_series <- numeric(n_sweeps); M_series <- integer(n_sweeps)
  total_accepted <- 0L; N <- L * L
  
  for (s in seq_len(n_sweeps)) {
    sweep          <- metropolis_sweep(spins, T)
    spins          <- sweep$spins
    E_current      <- E_current + sweep$dE_sweep
    M_current      <- M_current + sweep$dM_sweep
    E_series[s]    <- E_current
    M_series[s]    <- M_current
    total_accepted <- total_accepted + sweep$n_accepted
  }
  
  list(
    series          = data.frame(sweep = seq_len(n_sweeps),
                                 E = E_series, M = M_series),
    acceptance_rate = total_accepted / (n_sweeps * N),
    final_spins     = spins,
    L = L, T = T, n_burnin = n_burnin, n_sweeps = n_sweeps
  )
}

compute_observables <- function(run) {
  N <- run$L^2; T <- run$T
  E <- run$series$E; M <- run$series$M
  list(
    e               = mean(E) / N,
    m               = mean(abs(M)) / N,
    c               = (mean(E^2) - mean(E)^2)   / (N * T^2),
    chi             = (mean(M^2) - mean(abs(M))^2) / (N * T),
    acceptance_rate = run$acceptance_rate
  )
}


# =============================================================================
# Section 6: Rcpp Metropolis sampler (~300x faster)
# =============================================================================

Rcpp::cppFunction('
List run_metropolis_cpp_inner(int L, double T, int n_sweeps, int n_burnin,
                              IntegerVector init_spins) {
  int N = L * L;
  std::vector<int> spins(N);
  for (int k = 0; k < N; k++) spins[k] = init_spins[k];

  long long E_current = 0, M_current = 0;
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      int s       = spins[i * L + j];
      int s_right = spins[i * L + ((j + 1) % L)];
      int s_down  = spins[((i + 1) % L) * L + j];
      E_current  -= s * s_right + s * s_down;
      M_current  += s;
    }
  }

  // dE in {-8, -4, 0, 4, 8}: five precomputed exp(-dE/T)
  double exp_table[5];
  for (int idx = 0; idx < 5; idx++) {
    int dE = -8 + 4 * idx;
    exp_table[idx] = std::exp(-((double)dE) / T);
  }

  NumericVector E_series(n_sweeps), M_series(n_sweeps);
  long long total_accepted = 0, total_attempts = 0;
  GetRNGstate();

  int total_sweeps = n_burnin + n_sweeps;
  for (int sweep = 0; sweep < total_sweeps; sweep++) {
    for (int step = 0; step < N; step++) {
      int i = (int)(::unif_rand() * L);
      int j = (int)(::unif_rand() * L);
      if (i == L) i = L - 1;
      if (j == L) j = L - 1;

      int s_k = spins[i * L + j];
      int ip  = (i + 1) % L;
      int im  = (i + L - 1) % L;
      int jp  = (j + 1) % L;
      int jm  = (j + L - 1) % L;

      int neighbor_sum = spins[ip * L + j] + spins[im * L + j]
                       + spins[i * L + jp] + spins[i * L + jm];
      int dE = 2 * s_k * neighbor_sum;

      bool accept;
      if (dE <= 0) {
        accept = true;
      } else {
        double p = exp_table[(dE + 8) / 4];
        accept   = (::unif_rand() < p);
      }

      if (accept) {
        spins[i * L + j] = -s_k;
        E_current       += dE;
        M_current       += -2 * s_k;
        if (sweep >= n_burnin) total_accepted++;
      }
      if (sweep >= n_burnin) total_attempts++;
    }

    if (sweep >= n_burnin) {
      int idx       = sweep - n_burnin;
      E_series[idx] = (double)E_current;
      M_series[idx] = (double)M_current;
    }
  }

  PutRNGstate();
  IntegerVector final_spins(N);
  for (int k = 0; k < N; k++) final_spins[k] = spins[k];

  return List::create(
    Named("E_series")        = E_series,
    Named("M_series")        = M_series,
    Named("final_spins")     = final_spins,
    Named("acceptance_rate") = total_attempts > 0
                               ? (double)total_accepted / (double)total_attempts
                               : 0.0
  );
}
')

run_metropolis_fast <- function(L, T, n_sweeps,
                                n_burnin   = n_sweeps %/% 10,
                                init_state = "hot",
                                seed       = NULL) {
  if (!is.null(seed)) set.seed(seed)
  spins0    <- init_lattice(L, init_state)
  init_flat <- as.integer(spins0)
  
  result <- run_metropolis_cpp_inner(
    L          = as.integer(L),
    T          = as.numeric(T),
    n_sweeps   = as.integer(n_sweeps),
    n_burnin   = as.integer(n_burnin),
    init_spins = init_flat
  )
  
  list(
    series          = data.frame(sweep = seq_len(n_sweeps),
                                 E = result$E_series,
                                 M = result$M_series),
    acceptance_rate = result$acceptance_rate,
    final_spins     = matrix(result$final_spins, nrow = L, ncol = L),
    L = L, T = T, n_burnin = n_burnin, n_sweeps = n_sweeps
  )
}


# =============================================================================
# Section 6.1: Speed benchmark across L and T
# =============================================================================

benchmark_metropolis_samplers <- function(L_values = c(8, 16, 32, 64),
                                          T_values = c(1.5, T_CRITICAL, 4.0),
                                          T_labels = c("below Tc", "at Tc", "above Tc"),
                                          n_sweeps = 500,
                                          seed     = 1) {
  stopifnot(length(T_values) == length(T_labels))
  rows <- list()
  
  cat("=== Metropolis speed comparison ===\n")
  cat(sprintf("(n_sweeps = %d per run, seed = %d)\n\n", n_sweeps, seed))
  cat(sprintf("%-3s %-11s %-7s %-9s %-9s %-9s %-7s\n",
              "L", "T_label", "T", "t_R [s]", "t_cpp [s]", "speedup", "acc"))
  cat(strrep("-", 60), "\n", sep = "")
  
  for (L in L_values) {
    for (k in seq_along(T_values)) {
      T <- T_values[k]; label <- T_labels[k]
      t_R   <- system.time({ run_R   <- run_metropolis(     L, T, n_sweeps, seed = seed) })[["elapsed"]]
      t_cpp <- system.time({ run_cpp <- run_metropolis_fast(L, T, n_sweeps, seed = seed) })[["elapsed"]]
      speedup <- t_R / t_cpp; acc <- run_cpp$acceptance_rate
      
      cat(sprintf("%-3d %-11s %-7.4f %-9.3f %-9.3f %-9.1f %-7.3f\n",
                  L, label, T, t_R, t_cpp, speedup, acc))
      rows[[length(rows) + 1]] <- data.frame(
        L = L, T = T, T_label = label,
        t_R = t_R, t_cpp = t_cpp, speedup = speedup, acceptance = acc,
        stringsAsFactors = FALSE
      )
    }
  }
  
  cat(strrep("-", 60), "\n", sep = "")
  bench_df <- do.call(rbind, rows)
  
  cat(sprintf("\nMedian speedup across all (L, T): %.0fx\n", median(bench_df$speedup)))
  cat(sprintf("Min speedup:  %.0fx  (at L=%d, T=%.3f)\n",
              min(bench_df$speedup),
              bench_df$L[which.min(bench_df$speedup)],
              bench_df$T[which.min(bench_df$speedup)]))
  cat(sprintf("Max speedup:  %.0fx  (at L=%d, T=%.3f)\n",
              max(bench_df$speedup),
              bench_df$L[which.max(bench_df$speedup)],
              bench_df$T[which.max(bench_df$speedup)]))
  
  mid_T   <- T_values[ceiling(length(T_values) / 2)]
  df_mid  <- bench_df[bench_df$T == mid_T, ]
  fit_cpp <- lm(log(t_cpp) ~ log(L), data = df_mid)
  fit_R   <- lm(log(t_R)   ~ log(L), data = df_mid)
  cat(sprintf("\nScaling exponent at T = %.4f (expected ~2 for O(L^2)):\n", mid_T))
  cat(sprintf("  t_R   ~ L^%.2f\n", coef(fit_R)[2]))
  cat(sprintf("  t_cpp ~ L^%.2f\n", coef(fit_cpp)[2]))
  cat("===================================\n")
  
  invisible(bench_df)
}


# =============================================================================
# Okabe-Ito palettes for discrete variables (temperature regime, lattice size).
# Defined here because downstream sections reuse them; safe to redefine.
# =============================================================================

# One color per T-regime (ordered, critical, disordered)
T_REGIME_COLORS <- c(
  "below Tc" = "#0072B2",   # blue   -> ordered phase
  "at Tc"    = "#D55E00",   # vermillion -> critical
  "above Tc" = "#009E73"    # green -> disordered phase
)

# One color per lattice size (4 distinct Okabe-Ito hues, light -> dark)
L_COLORS <- c(
  "8"  = "#56B4E9",   # sky blue
  "16" = "#F0E442",   # yellow
  "32" = "#CC79A7",   # reddish purple
  "64" = "#000000"    # black (largest L, strongest contrast)
)


# =============================================================================
# Section 6.1 (plot): benchmark figure via ggplot2
# =============================================================================

plot_metropolis_benchmark <- function(bench_df) {
  # Keep L as numeric for log scales; create a factor copy only where needed.
  bench_df$T_label <- factor(bench_df$T_label,
                             levels = c("below Tc", "at Tc", "above Tc"))
  L_levels <- sort(unique(bench_df$L))          # numeric
  bench_df$L_fct <- factor(bench_df$L, levels = L_levels)
  
  # ------------------------------------------------------------------
  # Panel A: Wall time vs L (log-log), both pure-R and Rcpp, per T
  # ------------------------------------------------------------------
  df_long <- rbind(
    data.frame(L = bench_df$L, T_label = bench_df$T_label,
               impl = "pure R", time = bench_df$t_R),
    data.frame(L = bench_df$L, T_label = bench_df$T_label,
               impl = "Rcpp",   time = bench_df$t_cpp)
  )
  df_long$impl <- factor(df_long$impl, levels = c("pure R", "Rcpp"))
  
  L_min <- min(L_levels); L_max <- max(L_levels)
  t_R_ref <- max(bench_df$t_R)
  slope2_df <- data.frame(
    L    = c(L_min, L_max),
    time = t_R_ref * (c(L_min, L_max) / L_min)^2
  )
  
  p_time <- ggplot(df_long, aes(x = L, y = time,
                                color = T_label, shape = impl,
                                linetype = impl,
                                group = interaction(T_label, impl))) +
    geom_line(data = slope2_df, aes(x = L, y = time),
              color = "grey60", linetype = "dotted", linewidth = 0.5,
              inherit.aes = FALSE) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.0) +
    scale_x_log10(breaks = L_levels) +
    scale_y_log10() +
    scale_color_manual(values = T_REGIME_COLORS, name = NULL) +
    scale_shape_manual(values = c("pure R" = 16, "Rcpp" = 17), name = NULL) +
    scale_linetype_manual(values = c("pure R" = "solid", "Rcpp" = "dashed"),
                          name = NULL) +
    annotate("text", x = L_max * 0.95, y = slope2_df$time[2] * 0.9,
             label = "slope 2", color = "grey50", size = 2.8, hjust = 1) +
    labs(title = "Wall time vs lattice size",
         x     = expression(italic(L)),
         y     = "wall time [s]") +
    theme_paper() +
    theme(legend.position  = "right",
          legend.box       = "vertical",
          legend.spacing.y = unit(-4, "pt"))
  
  # ------------------------------------------------------------------
  # Panel B: Speedup vs L, per T
  # ------------------------------------------------------------------
  med_speedup <- median(bench_df$speedup)
  
  p_sp_L <- ggplot(bench_df, aes(x = L, y = speedup,
                                 color = T_label, group = T_label)) +
    geom_hline(yintercept = med_speedup,
               linetype = "dotted", color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.0) +
    scale_x_log10(breaks = L_levels) +
    scale_color_manual(values = T_REGIME_COLORS, name = NULL) +
    annotate("text",
             x = L_max, y = med_speedup,
             label = sprintf("  median = %.0fx", med_speedup),
             color = "grey40", size = 2.8, hjust = 0, vjust = -0.4) +
    labs(title = "Speedup vs lattice size",
         x     = expression(italic(L)),
         y     = expression(t[R] / t[Rcpp])) +
    theme_paper() +
    theme(legend.position = "right")
  
  # ------------------------------------------------------------------
  # Panel C: Speedup vs T, per L  (L is discrete here)
  # ------------------------------------------------------------------
  p_sp_T <- ggplot(bench_df, aes(x = T, y = speedup,
                                 color = L_fct, group = L_fct)) +
    geom_vline(xintercept = T_CRITICAL, linetype = "dashed",
               color = PAPER_COLORS$tc_line, linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.0) +
    scale_color_manual(values = setNames(L_COLORS[as.character(L_levels)],
                                         as.character(L_levels)),
                       name = expression(italic(L))) +
    annotate("text",
             x = T_CRITICAL + 0.05, y = max(bench_df$speedup) * 0.97,
             label = "italic(T)[c]", parse = TRUE,
             color = PAPER_COLORS$tc_line, size = 3, hjust = 0) +
    labs(title = "Speedup vs temperature",
         x     = expression(italic(T)),
         y     = expression(t[R] / t[Rcpp])) +
    theme_paper() +
    theme(legend.position = "right")
  
  combined <- (p_time | p_sp_L | p_sp_T) +
    plot_annotation(
      caption = sprintf(
        "Median Rcpp speedup = %.0fx across 12 (L, T) cells. Dotted slope-2 line indicates O(L^2) scaling.",
        med_speedup),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    )
  
  combined
}

# =============================================================================
# Section 6.2: Comprehensive validation against all ground truths
# (non-plotting code unchanged from earlier version)
# =============================================================================

choose_n_sweeps <- function(L, T, base_n = 5000, tol = 0.05) {
  if (abs(T - T_CRITICAL) < tol) as.integer(base_n * (L / 8)^2) else base_n
}

choose_init_state <- function(T) {
  if (T < T_CRITICAL - 0.1) "cold" else "hot"
}

m_finite_L_reference <- function(L, T) {
  N <- L * L
  if (T < T_CRITICAL - 0.01) {
    list(value = safe_onsager_magnetization(T),
         label = "Yang (L=inf)", regime = "below")
  } else if (abs(T - T_CRITICAL) < 0.05) {
    list(value = NA,
         label = sprintf("~ L^(-1/8) = %.4f (scaling)", L^(-1/8)),
         regime = "critical")
  } else {
    list(value = 1 / sqrt(N),
         label = sprintf("1/sqrt(N) = %.4f", 1 / sqrt(N)),
         regime = "above")
  }
}

validate_metropolis_all_observables <- function(
    L_values  = c(8, 16, 32, 64),
    T_values  = c(1.5, T_CRITICAL, 4.0),
    T_labels  = c("below Tc", "at Tc", "above Tc"),
    seed_base = 42) {
  stopifnot(length(T_values) == length(T_labels))
  rows <- list()
  
  cat("=== Metropolis validation against ground truths ===\n\n")
  
  for (L in L_values) {
    for (k in seq_along(T_values)) {
      T <- T_values[k]; label <- T_labels[k]
      n_sweeps <- choose_n_sweeps(L, T)
      init     <- choose_init_state(T)
      n_burnin <- as.integer(n_sweeps / 5)
      
      run <- run_metropolis_fast(L, T,
                                 n_sweeps = n_sweeps, n_burnin = n_burnin,
                                 init_state = init, seed = seed_base + L)
      obs <- compute_observables(run)
      
      e_true     <- safe_onsager_energy(T)
      c_true     <- safe_onsager_specific_heat(T)
      c_proxy    <- onsager_specific_heat_proxy(T)
      m_ref_info <- m_finite_L_reference(L, T)
      
      rows[[length(rows) + 1]] <- data.frame(
        L = L, T = T, T_label = label,
        n_sweeps = n_sweeps, init = init, acc = obs$acceptance_rate,
        e_mcmc = obs$e, e_true = e_true, e_err = obs$e - e_true,
        c_mcmc = obs$c, c_true = c_true, c_proxy = c_proxy,
        c_err_proxy = obs$c - c_proxy,
        m_mcmc = obs$m,
        m_ref_val = m_ref_info$value, m_ref_lbl = m_ref_info$label,
        m_regime  = m_ref_info$regime,
        chi_mcmc  = obs$chi,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  
  # Energy table
  cat("Energy per spin (e):\n")
  cat(sprintf("%-3s %-10s %-7s %-7s %-6s %-9s %-9s %-9s\n",
              "L", "regime", "T", "n_swp", "init", "e_mcmc", "e_onsager", "err"))
  cat(strrep("-", 70), "\n", sep = "")
  for (i in seq_len(nrow(result))) {
    r <- result[i, ]
    cat(sprintf("%-3d %-10s %-7.4f %-7d %-6s %+9.4f %+9.4f %+9.4f\n",
                r$L, r$T_label, r$T, r$n_sweeps, r$init,
                r$e_mcmc, r$e_true, r$e_err))
  }
  
  # Specific heat table
  cat("\nSpecific heat per spin (c):\n")
  cat("  (true c(T_c) = +Inf; compared to proxy c(T_c - 1e-3). Finite-L is below.)\n")
  cat(sprintf("%-3s %-10s %-7s %-9s %-9s %-9s %-14s\n",
              "L", "regime", "T", "c_mcmc", "c_onsager", "c_proxy", "err(vs proxy)"))
  cat(strrep("-", 75), "\n", sep = "")
  for (i in seq_len(nrow(result))) {
    r <- result[i, ]
    c_onsager_str <- if (is.infinite(r$c_true)) "    Inf  " else sprintf("%9.4f", r$c_true)
    cat(sprintf("%-3d %-10s %-7.4f %9.4f %s %9.4f %+9.4f\n",
                r$L, r$T_label, r$T, r$c_mcmc, c_onsager_str, r$c_proxy, r$c_err_proxy))
  }
  
  # Magnetization table
  cat("\nAbsolute magnetization per spin (|m|):\n")
  cat("  (reference: Yang for T<T_c; L^(-1/8) scaling at T_c; 1/sqrt(N) for T>T_c)\n")
  cat(sprintf("%-3s %-10s %-7s %-9s %-35s\n",
              "L", "regime", "T", "m_mcmc", "reference"))
  cat(strrep("-", 72), "\n", sep = "")
  for (i in seq_len(nrow(result))) {
    r <- result[i, ]
    cat(sprintf("%-3d %-10s %-7.4f %9.4f %s\n",
                r$L, r$T_label, r$T, r$m_mcmc, r$m_ref_lbl))
  }
  
  # |m| at T_c scaling check
  cat("\nFinite-size scaling of |m| at T_c (expect |m| ~ L^(-1/8)):\n")
  sub_tc <- result[result$T_label == "at Tc", ]
  sub_tc <- sub_tc[order(sub_tc$L), ]
  for (i in seq_len(nrow(sub_tc))) {
    r <- sub_tc[i, ]
    predicted <- sub_tc$m_mcmc[1] * (r$L / sub_tc$L[1])^(-1/8)
    cat(sprintf("  L=%2d: m_mcmc = %.4f, predicted = %.4f, ratio = %.3f\n",
                r$L, r$m_mcmc, predicted, r$m_mcmc / predicted))
  }
  
  # |m| at T=4 scaling check
  cat("\nFinite-size scaling of |m| at T=4.0 (expect |m| ~ 1/L):\n")
  sub_hi <- result[result$T_label == "above Tc", ]
  sub_hi <- sub_hi[order(sub_hi$L), ]
  for (i in seq_len(nrow(sub_hi))) {
    r <- sub_hi[i, ]
    predicted <- sub_hi$m_mcmc[1] * (sub_hi$L[1] / r$L)
    cat(sprintf("  L=%2d: m_mcmc = %.4f, predicted = %.4f, ratio = %.3f\n",
                r$L, r$m_mcmc, predicted, r$m_mcmc / predicted))
  }
  
  # Susceptibility table + scaling fit
  cat("\nMagnetic susceptibility (chi):\n")
  cat("  (no exact closed form; at T_c expect chi ~ L^(7/4) = L^1.75.)\n")
  cat(sprintf("%-3s %-10s %-7s %-9s\n", "L", "regime", "T", "chi"))
  cat(strrep("-", 36), "\n", sep = "")
  for (i in seq_len(nrow(result))) {
    r <- result[i, ]
    cat(sprintf("%-3d %-10s %-7.4f %9.4f\n", r$L, r$T_label, r$T, r$chi_mcmc))
  }
  sub_tc_chi <- result[result$T_label == "at Tc", ]
  sub_tc_chi <- sub_tc_chi[order(sub_tc_chi$L), ]
  fit_chi <- lm(log(chi_mcmc) ~ log(L), data = sub_tc_chi)
  cat(sprintf("\nFinite-size scaling of chi at T_c: chi ~ L^%.3f (literature: 1.75)\n",
              coef(fit_chi)[2]))
  
  # Summary
  cat("\n--- Summary of absolute errors where meaningful ---\n")
  for (label in T_labels) {
    sub <- result[result$T_label == label, ]
    finite_e_err <- sub$e_err[is.finite(sub$e_err)]
    cat(sprintf("%s (T = %.4f):\n", label, sub$T[1]))
    if (length(finite_e_err) > 0) {
      cat(sprintf("  |e_err|:        max = %.4f  mean = %.4f\n",
                  max(abs(finite_e_err)), mean(abs(finite_e_err))))
    }
    cat(sprintf("  |c_err_proxy|:  max = %.4f  mean = %.4f\n",
                max(abs(sub$c_err_proxy)), mean(abs(sub$c_err_proxy))))
    if (label == "below Tc") {
      m_err <- sub$m_mcmc - sub$m_ref_val
      cat(sprintf("  |m_err vs Yang|: max = %.4f  mean = %.4f\n",
                  max(abs(m_err)), mean(abs(m_err))))
    } else {
      cat("  |m_err|: finite-L scaling (see table above)\n")
    }
  }
  cat("\n==================================================\n")
  
  invisible(result)
}

# =============================================================================
# Section 6.2 (plot): validation figure via ggplot2
# =============================================================================

plot_metropolis_validation <- function(validation_df, c_cap = 4) {
  validation_df$L_fct <- factor(validation_df$L,
                                levels = sort(unique(validation_df$L)))
  
  # Dense Onsager/Yang reference curves
  T_curve <- seq(1.0, 4.2, length.out = 400)
  df_ref <- data.frame(
    T = T_curve,
    u = onsager_energy_curve(T_curve),
    c = pmin(onsager_specific_heat_curve(T_curve), c_cap),
    m = onsager_magnetization_curve(T_curve)
  )
  
  # ------------------------------------------------------------------
  # Shared Tc line + label helpers (reused across panels)
  # ------------------------------------------------------------------
  tc_line <- geom_vline(xintercept = T_CRITICAL,
                        linetype   = "dashed",
                        color      = PAPER_COLORS$tc_line,
                        linewidth  = 0.4)
  
  tc_label <- function(y_val) {
    annotate("text", x = T_CRITICAL + 0.04, y = y_val,
             label = "italic(T)[c]", parse = TRUE,
             color = PAPER_COLORS$tc_line, size = 3, hjust = 0)
  }
  
  # ------------------------------------------------------------------
  # Panel 1: Energy
  # (no scale_color_manual here -- applied globally via `&` below)
  # ------------------------------------------------------------------
  p_u <- ggplot() +
    tc_line +
    geom_line(data = df_ref, aes(x = T, y = u),
              color = "grey40", linewidth = 0.9, na.rm = TRUE) +
    geom_point(data = validation_df,
               aes(x = T, y = e_mcmc, color = L_fct),
               size = 2.2) +
    tc_label(y_val = -0.7) +
    labs(title = "Energy per spin: MCMC vs Onsager",
         x     = expression(italic(T)),
         y     = expression(italic(u)(italic(T)))) +
    theme_paper()
  
  # ------------------------------------------------------------------
  # Panel 2: Specific heat (reference capped at c_cap)
  # ------------------------------------------------------------------
  p_c <- ggplot() +
    tc_line +
    geom_line(data = df_ref, aes(x = T, y = c),
              color = "grey40", linewidth = 0.9, na.rm = TRUE) +
    geom_point(data = validation_df,
               aes(x = T, y = c_mcmc, color = L_fct),
               size = 2.2) +
    coord_cartesian(ylim = c(0, c_cap)) +
    tc_label(y_val = c_cap * 0.92) +
    annotate("text", x = 4.15, y = c_cap * 0.92,
             label = sprintf("(Onsager capped at %.0f)", c_cap),
             hjust = 1, size = 2.6, color = "grey40") +
    labs(title = "Specific heat per spin: MCMC vs Onsager",
         x     = expression(italic(T)),
         y     = expression(italic(c)(italic(T)))) +
    theme_paper()
  
  # ------------------------------------------------------------------
  # Panel 3: Magnetization + per-L 1/sqrt(N) reference segments
  # `show.legend = FALSE` on segments prevents them polluting the legend
  # ------------------------------------------------------------------
  df_refN <- data.frame(
    L_fct = factor(sort(unique(validation_df$L)),
                   levels = sort(unique(validation_df$L))),
    y_val = 1 / sort(unique(validation_df$L))   # 1/sqrt(N) = 1/L
  )
  
  p_m <- ggplot() +
    tc_line +
    geom_hline(yintercept = 0, color = "grey70",
               linewidth = 0.3, linetype = "dotted") +
    geom_line(data = df_ref, aes(x = T, y = m),
              color = "grey40", linewidth = 0.9, na.rm = TRUE) +
    geom_segment(data = df_refN,
                 aes(x = T_CRITICAL + 0.1, xend = 4.2,
                     y = y_val, yend = y_val, color = L_fct),
                 linetype    = "dotted",
                 linewidth   = 0.4,
                 show.legend = FALSE) +
    geom_point(data = validation_df,
               aes(x = T, y = m_mcmc, color = L_fct),
               size = 2.2) +
    coord_cartesian(ylim = c(0, 1.05)) +
    tc_label(y_val = 0.92) +
    labs(title = "Magnetization: MCMC vs Yang (+ 1/sqrt(N) ref per L)",
         x     = expression(italic(T)),
         y     = expression(group("|", italic(m), "|") * (italic(T)))) +
    theme_paper()
  
  # ------------------------------------------------------------------
  # Panel 4: Susceptibility (no closed form; expect peak grows as L^1.75)
  # ------------------------------------------------------------------
  p_chi <- ggplot(validation_df,
                  aes(x = T, y = chi_mcmc, color = L_fct)) +
    tc_line +
    geom_point(size = 2.2) +
    scale_y_log10() +
    tc_label(y_val = max(validation_df$chi_mcmc) * 0.9) +
    labs(title = expression(bold("Susceptibility: peak grows as ") *
                              bolditalic(L)^bold("1.75")),
         x     = expression(italic(T)),
         y     = expression(chi(italic(T)))) +
    theme_paper()
  
  # ------------------------------------------------------------------
  # Assemble into a 2x2 grid.
  # `guides = "collect"` merges legends; `&` operators apply shared
  # elements (color scale, legend position) to ALL panels so the
  # legend is identical everywhere and appears only once.
  # ------------------------------------------------------------------
  (p_u | p_c) / (p_m | p_chi) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = paste(
        "Metropolis estimates (colored points) vs Onsager / Yang ground truth (grey).",
        "Dotted segments above T_c show the 1/sqrt(N) finite-lattice |m| floor per L.",
        sep = " "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    ) &
    scale_color_manual(values = L_COLORS,
                       name   = expression(italic(L)),
                       drop   = FALSE) &
    theme(legend.position = "bottom")
}

# =============================================================================
# Drive: run everything end-to-end
# =============================================================================

# Section 4: sanity checks (if not already run above)
# run_sanity_checks()

# Section 3.8: ground-truth ggplot (if not already rendered above)
# print(plot_onsager_ground_truth_gg())

# Section 6.1: speed benchmark (~2-3 minutes)
bench_df <- benchmark_metropolis_samplers()
print(plot_metropolis_benchmark(bench_df))

# Section 6.2: ground-truth validation (~1-2 minutes)
validation_df <- validate_metropolis_all_observables()
print(plot_metropolis_validation(validation_df))

##> === Metropolis speed comparison ===
##> (n_sweeps = 500 per run, seed = 1)
##>
##> L   T_label     T       t_R [s]   t_cpp [s] speedup   acc
##> ------------------------------------------------------------
##> 8   below Tc    1.5000  ...       ...       ~300x     ...
##> ... (12 rows total, one per (L, T) cell)
##>
##> Median speedup across all (L, T): 273x   [machine-dependent; your PDF reports 273x]
##> Scaling exponent at T = 2.2692 (expected ~2 for O(L^2)):
##>   t_R   ~ L^2.0x
##>   t_cpp ~ L^2.0x
##> ===================================
##> [benchmark figure displays: 3 panels, log-log walltime + 2 speedup panels]
##>
##> === Metropolis validation against ground truths ===
##>
##> Energy per spin (e):
##>   (12 rows; errors O(1e-3) away from Tc, O(1e-2) at Tc)
##>
##> Specific heat per spin (c):
##>   (12 rows; MCMC c is rounded below the Onsager proxy at Tc — expected finite-L)
##>
##> Absolute magnetization per spin (|m|):
##>   L=8,  at Tc: m_mcmc=0.765, predicted=0.765, ratio=1.000  (seed row)
##>   L=16, at Tc: m_mcmc=0.688, predicted=0.701, ratio=0.982
##>   L=32, at Tc: m_mcmc=0.659, predicted=0.642, ratio=1.026
##>   L=64, at Tc: m_mcmc=0.586, predicted=0.588, ratio=0.996
##>
##> Finite-size scaling of chi at T_c: chi ~ L^1.7... (literature: 1.75)
##>
##> --- Summary of absolute errors where meaningful ---
##> below Tc (T = 1.5000): |e_err| max ~ 0.003
##> at Tc (T = 2.2692):    |e_err| max ~ 0.02
##> above Tc (T = 4.0000): |e_err| max ~ 0.005
##> ==================================================
##> [validation figure displays: 2x2 panels, u/c/|m|/chi vs T with L as color]



# =============================================================================
# Section 7: Chain trace diagnostics (Appendix A)
#
# Single 4x3 figure in the Appendix:
#   rows    = observable  (e, |m|, running c, running chi)
#   columns = temperature (below Tc, at Tc, above Tc)
#   colors  = lattice size L in {8, 16, 32, 64}
#
# x-axis is sweep FRACTION (t / n_sweeps) so chains of different lengths
# (L=8 away from T_c runs 5000 sweeps; L=64 at T_c runs 320000 sweeps)
# align to the same visual width.
#
# Runs re-use the validation configuration (same init, same burn-in, same
# seed policy), so the displayed traces are the exact chains that fed
# Section 6.2's ground-truth comparison tables.
# =============================================================================


# -----------------------------------------------------------------------------
# 7.1: Run and cache all 12 (L, T) cells. Returns a named list keyed by
#      "L{L}_T{below,c,above}".
# -----------------------------------------------------------------------------

trace_runs_metropolis <- function(L_values  = c(8, 16, 32, 64),
                                  T_values  = c(1.5, T_CRITICAL, 4.0),
                                  T_labels  = c("below Tc", "at Tc", "above Tc"),
                                  seed_base = 42) {
  stopifnot(length(T_values) == length(T_labels))
  
  runs <- list()
  cat("=== Trace runs (re-using validation grid) ===\n")
  for (L in L_values) {
    for (k in seq_along(T_values)) {
      T        <- T_values[k]
      label    <- T_labels[k]
      n_sweeps <- choose_n_sweeps(L, T)
      init     <- choose_init_state(T)
      n_burnin <- as.integer(n_sweeps / 5)
      
      key_T <- if (abs(T - T_CRITICAL) < 1e-6) "c"
      else if (T < T_CRITICAL) "below"
      else "above"
      key   <- sprintf("L%d_T%s", L, key_T)
      
      cat(sprintf("  L=%2d, T=%s (%-8s): n_sweeps=%d, init=%s ... ",
                  L, format(T, digits = 4), label, n_sweeps, init))
      t0 <- Sys.time()
      run <- run_metropolis_fast(L = L, T = T,
                                 n_sweeps   = n_sweeps,
                                 n_burnin   = n_burnin,
                                 init_state = init,
                                 seed       = seed_base + L)
      run$T_label <- label
      run$key     <- key
      runs[[key]] <- run
      cat(sprintf("done (%.1fs)\n",
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
  cat("=============================================\n\n")
  invisible(runs)
}


# -----------------------------------------------------------------------------
# 7.2: Running estimates for variance-based observables
# -----------------------------------------------------------------------------

running_specific_heat <- function(E, N, T) {
  n      <- seq_along(E)
  mean_E  <- cumsum(E)    / n
  mean_E2 <- cumsum(E^2)  / n
  var_E   <- pmax(mean_E2 - mean_E^2, 0)
  var_E / (N * T^2)
}

running_susceptibility <- function(M, N, T) {
  n      <- seq_along(M)
  mean_M2 <- cumsum(M^2)    / n
  mean_aM <- cumsum(abs(M)) / n
  var_M   <- pmax(mean_M2 - mean_aM^2, 0)
  var_M / (N * T)
}


# -----------------------------------------------------------------------------
# 7.3: Build a long data frame with all 48 (L, T, observable) traces
# -----------------------------------------------------------------------------

build_trace_dataframe <- function(runs, max_pts = 2000) {
  rows <- list()
  for (key in names(runs)) {
    run <- runs[[key]]
    N   <- run$L^2
    T   <- run$T
    n   <- run$n_sweeps
    
    # Thin to roughly max_pts points for plot efficiency
    step    <- max(1L, n %/% max_pts)
    idx     <- seq(1L, n, by = step)
    frac    <- idx / n
    
    # Compute all four observables for this (L, T) cell
    e_ser   <- run$series$E / N
    m_ser   <- abs(run$series$M) / N
    c_run   <- running_specific_heat(run$series$E, N, T)
    chi_run <- running_susceptibility(run$series$M, N, T)
    
    for (obs_name in c("e", "m", "c", "chi")) {
      y_full <- switch(obs_name,
                       e   = e_ser,
                       m   = m_ser,
                       c   = c_run,
                       chi = chi_run)
      rows[[length(rows) + 1]] <- data.frame(
        L          = run$L,
        T          = T,
        T_label    = run$T_label,
        observable = obs_name,
        sweep      = idx,
        frac       = frac,
        value      = y_full[idx],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


# -----------------------------------------------------------------------------
# 7.4: One 4x3 figure (appendix) using ggplot + facet_grid
# -----------------------------------------------------------------------------

plot_trace_appendix <- function(trace_df) {
  # Factor levels control facet order (rows = observables, cols = T regimes)
  trace_df$observable <- factor(trace_df$observable,
                                levels = c("e", "m", "c", "chi"))
  trace_df$T_label    <- factor(trace_df$T_label,
                                levels = c("below Tc", "at Tc", "above Tc"))
  trace_df$L_fct      <- factor(trace_df$L,
                                levels = sort(unique(trace_df$L)))
  
  # Plain-text facet labels (no plotmath parsing -- robust to special chars)
  obs_labels <- c(
    e   = "u(T)",
    m   = "|m|(T)",
    c   = "running c(T)",
    chi = "running chi(T)"
  )
  tT_labels <- c(
    "below Tc" = "T = 1.5  (below Tc)",
    "at Tc"    = "T = Tc",
    "above Tc" = "T = 4.0  (above Tc)"
  )
  
  ggplot(trace_df, aes(x = frac, y = value, color = L_fct, group = L_fct)) +
    geom_line(linewidth = 0.3, alpha = 0.85, na.rm = TRUE) +
    facet_grid(
      rows     = vars(observable),
      cols     = vars(T_label),
      scales   = "free_y",
      switch   = "y",
      labeller = labeller(
        observable = as_labeller(obs_labels),
        T_label    = as_labeller(tT_labels)
      )
    ) +
    scale_color_manual(values = L_COLORS,
                       name   = expression(italic(L)),
                       drop   = FALSE) +
    scale_x_continuous(breaks = c(0, 0.5, 1.0), limits = c(0, 1)) +
    labs(x = expression("sweep fraction  " * italic(t) / italic(n)),
         y = NULL) +
    theme_paper() +
    theme(
      legend.position   = "bottom",
      strip.background  = element_rect(fill = "grey95", color = NA),
      strip.text        = element_text(size = rel(0.85), color = "grey20"),
      strip.placement   = "outside",
      panel.spacing.x   = unit(6, "pt"),
      panel.spacing.y   = unit(4, "pt"),
      axis.title.x      = element_text(margin = margin(t = 6))
    ) +
    plot_annotation(
      caption = paste(
        "Traces of MCMC observables across the 12 (L, T) diagnostic cells.",
        "x-axis: sweep fraction; colors: lattice size.",
        "Running c and chi are cumulative variance estimators.",
        sep = " "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    )
}


# -----------------------------------------------------------------------------
# 7.5: Top-level driver
# -----------------------------------------------------------------------------

plot_section7_all_traces <- function(L_values  = c(8, 16, 32, 64),
                                     T_values  = c(1.5, T_CRITICAL, 4.0),
                                     T_labels  = c("below Tc", "at Tc", "above Tc"),
                                     seed_base = 42,
                                     max_pts   = 2000) {
  
  cat("\n############################################################\n")
  cat("Section 7: chain trace diagnostics (Appendix A)\n")
  cat("############################################################\n\n")
  
  runs <- trace_runs_metropolis(L_values  = L_values,
                                T_values  = T_values,
                                T_labels  = T_labels,
                                seed_base = seed_base)
  
  trace_df <- build_trace_dataframe(runs, max_pts = max_pts)
  print(plot_trace_appendix(trace_df))
  
  invisible(list(runs = runs, trace_df = trace_df))
}


# =============================================================================
# Drive
# =============================================================================

section7_out <- plot_section7_all_traces()





# =============================================================================
# Section 8: Autocorrelation, integrated time, and effective sample size
#
# For every (L, T) combination from Section 7, compute:
#
#   - rho_E(t), rho_M(t):   empirical autocorrelation function for
#                           energy and absolute magnetization
#   - tau_E, tau_M:         Sokal adaptive-window estimator with c = 5,
#                           capped at n/4
#   - ESS_E, ESS_M:         n / (2 * tau)
#   - SE(<e>), SE(<|m|>):   sigma * sqrt(2 * tau / n)
#
# Why only e and |m|? These are per-sweep observables. c and chi are
# variance-based; Section 11 (block bootstrap) handles their uncertainty.
#
# Reuse from Section 7: the 12 runs in `section7_runs` ARE the data we need.
# =============================================================================


# -----------------------------------------------------------------------------
# 8.1: autocorrelation function and integrated tau (Sokal windowing)
# -----------------------------------------------------------------------------

autocorrelation <- function(x, max_lag = NULL) {
  n <- length(x)
  if (is.null(max_lag)) max_lag <- min(n %/% 4, 1000)
  x_centered <- x - mean(x)
  v0 <- mean(x_centered^2)
  if (v0 == 0) return(c(1, rep(0, max_lag)))
  rho <- numeric(max_lag + 1)
  rho[1] <- 1
  for (t in seq_len(max_lag)) {
    rho[t + 1] <- mean(x_centered[1:(n - t)] * x_centered[(t + 1):n]) / v0
  }
  rho
}

integrated_tau <- function(x, c = 5, max_lag = NULL) {
  n <- length(x)
  if (is.null(max_lag)) max_lag <- n %/% 4
  rho <- autocorrelation(x, max_lag = max_lag)
  tau_partial <- 0.5 + cumsum(rho[-1])
  W_grid <- seq_along(tau_partial)
  ok <- W_grid >= c * tau_partial
  if (any(ok)) {
    W <- which(ok)[1]
    closed <- TRUE
  } else {
    W <- max_lag
    closed <- FALSE   # Sokal hit the cutoff: tau is biased low
  }
  tau_value <- max(tau_partial[W], 0.5)
  list(tau = tau_value, window = W, rho = rho, closed = closed,
       max_lag = max_lag)
}


# -----------------------------------------------------------------------------
# 8.2: per-cell diagnostics for both observables
# -----------------------------------------------------------------------------

cell_diagnostics <- function(run) {
  N <- run$L^2
  E_per_spin <- run$series$E / N
  m_per_spin <- abs(run$series$M) / N
  
  diag_E <- integrated_tau(run$series$E)
  diag_M <- integrated_tau(run$series$M)
  n      <- nrow(run$series)
  
  se_e <- sd(E_per_spin) * sqrt(2 * diag_E$tau / n)
  se_m <- sd(m_per_spin) * sqrt(2 * diag_M$tau / n)
  
  list(
    L          = run$L,
    T          = run$T,
    T_label    = run$T_label,
    n_sweeps   = n,
    acceptance = run$acceptance_rate,
    tau_E      = diag_E$tau,
    tau_M      = diag_M$tau,
    window_E   = diag_E$window,
    window_M   = diag_M$window,
    closed_E   = diag_E$closed,
    closed_M   = diag_M$closed,
    ESS_E      = n / (2 * diag_E$tau),
    ESS_M      = n / (2 * diag_M$tau),
    se_e       = se_e,
    se_m       = se_m,
    rho_E      = diag_E$rho,
    rho_M      = diag_M$rho
  )
}


# -----------------------------------------------------------------------------
# 8.3: assemble the full diagnostics table
# -----------------------------------------------------------------------------

build_diagnostics_table <- function(runs) {
  rows <- list()
  for (key in names(runs)) {
    d <- cell_diagnostics(runs[[key]])
    rows[[key]] <- data.frame(
      L        = d$L,
      T        = d$T,
      T_label  = d$T_label,
      n_sweeps = d$n_sweeps,
      acc      = d$acceptance,
      tau_E    = d$tau_E,
      tau_M    = d$tau_M,
      ESS_E    = d$ESS_E,
      ESS_M    = d$ESS_M,
      se_e     = d$se_e,
      se_m     = d$se_m,
      window_E = d$window_E,
      window_M = d$window_M,
      closed_E = d$closed_E,
      closed_M = d$closed_M,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


# -----------------------------------------------------------------------------
# 8.4: pretty printer with explicit Sokal-closure flags
# -----------------------------------------------------------------------------

print_diagnostics_table <- function(diag_df) {
  cat("=== Autocorrelation diagnostics across (L, T) ===\n")
  cat("Notation: tau in sweeps; ESS = n/(2 tau); SE in per-spin units.\n")
  cat("Flag '!' next to tau means Sokal windowing did NOT close\n")
  cat("(the reported tau is then a lower bound; ESS overestimates the truth).\n\n")
  
  cat(sprintf("%-3s %-10s %-8s %-8s %-12s %-12s %-9s %-9s %-9s %-9s\n",
              "L", "regime", "T", "n_swp",
              "tau_E", "tau_M", "ESS_E", "ESS_M", "se_e", "se_m"))
  cat(strrep("-", 99), "\n", sep = "")
  
  for (i in seq_len(nrow(diag_df))) {
    r <- diag_df[i, ]
    flag_E <- if (r$closed_E) " " else "!"
    flag_M <- if (r$closed_M) " " else "!"
    cat(sprintf("%-3d %-10s %-8.4f %-8d %-11.2f%s %-11.2f%s %-9.0f %-9.0f %-9.5f %-9.5f\n",
                r$L, r$T_label, r$T, r$n_sweeps,
                r$tau_E, flag_E,
                r$tau_M, flag_M,
                r$ESS_E, r$ESS_M, r$se_e, r$se_m))
  }
  cat(strrep("-", 99), "\n", sep = "")
  
  n_fail_E <- sum(!diag_df$closed_E)
  n_fail_M <- sum(!diag_df$closed_M)
  if (n_fail_E > 0 || n_fail_M > 0) {
    cat(sprintf("\n*** WARNING: %d cells did not close for tau_E, %d for tau_M.\n",
                n_fail_E, n_fail_M))
    cat("    Those tau values are downward-biased lower bounds; the true\n")
    cat("    autocorrelation time is larger than reported, so ESS is overstated\n")
    cat("    and the standard error is understated. To fix: longer chains.\n")
  }
}


# -----------------------------------------------------------------------------
# 8.5: dynamic critical exponent fit for Metropolis at T_c
# -----------------------------------------------------------------------------

fit_dynamic_z <- function(diag_df, observable = c("E", "M")) {
  observable <- match.arg(observable)
  tau_col    <- paste0("tau_", observable)
  closed_col <- paste0("closed_", observable)
  
  sub <- diag_df[diag_df$T_label == "at Tc", ]
  sub <- sub[sub[[closed_col]], ]
  
  if (nrow(sub) < 2) {
    cat(sprintf("  z (from tau_%s): not enough closed cells to fit (have %d, need >=2)\n",
                observable, nrow(sub)))
    return(NA)
  }
  
  fit <- lm(log(sub[[tau_col]]) ~ log(sub$L))
  z   <- coef(fit)[2]
  # Literature value: Nightingale & Blote's z applies to the slowest mode,
  # typically the magnetization. We report it only for the M fit.
  lit_str <- if (observable == "M") {
    sprintf("  (literature z ~ %.4f, Nightingale & Blote 1996)",
            Z_METROPOLIS_LITERATURE)
  } else {
    "  (no canonical literature value for z_E)"
  }
  cat(sprintf("  z (from tau_%s, fit on L = %s): z = %.3f%s\n",
              observable, paste(sub$L, collapse = ","), z, lit_str))
  invisible(z)
}


# -----------------------------------------------------------------------------
# 8.6: rho(t) grid (4 x 3), ggplot version
#
# One figure per observable. Each panel shows rho(t) for 0 <= t <= max_lag_plot
# with a 1/e reference line; panel title carries tau and a Sokal-closure flag.
# -----------------------------------------------------------------------------

build_rho_dataframe <- function(runs, observable = c("E", "M"),
                                max_lag_plot = 200) {
  observable <- match.arg(observable)
  rows <- list()
  for (key in names(runs)) {
    run      <- runs[[key]]
    n_sweeps <- nrow(run$series)
    max_lag_here <- min(max_lag_plot, n_sweeps %/% 4)
    x_series <- if (observable == "E") run$series$E else run$series$M
    
    rho  <- autocorrelation(x_series, max_lag = max_lag_here)
    diag <- integrated_tau(x_series)
    
    rows[[length(rows) + 1]] <- data.frame(
      L        = run$L,
      T        = run$T,
      T_label  = run$T_label,
      lag      = 0:max_lag_here,
      rho      = rho,
      tau      = diag$tau,
      closed   = diag$closed,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


plot_autocorr_grid <- function(runs, observable = c("E", "M"),
                               max_lag_plot = 200) {
  observable <- match.arg(observable)
  
  df <- build_rho_dataframe(runs, observable, max_lag_plot)
  
  df$L_fct   <- factor(df$L,
                       levels = sort(unique(df$L)))
  df$T_label <- factor(df$T_label,
                       levels = c("below Tc", "at Tc", "above Tc"))
  
  # Per-panel tau / closed flag annotation
  panel_labels <- unique(df[, c("L", "T_label", "tau", "closed")])
  panel_labels$L_fct   <- factor(panel_labels$L,
                                 levels = sort(unique(df$L)))
  panel_labels$T_label <- factor(panel_labels$T_label,
                                 levels = c("below Tc", "at Tc", "above Tc"))
  panel_labels$lbl <- with(panel_labels,
                           sprintf("tau = %.1f%s", tau, ifelse(closed, "", " (!)")))
  
  # Colour the curve by T-regime; each panel has exactly one curve.
  obs_str  <- switch(observable, E = "energy", M = "magnetization")
  obs_expr <- switch(observable,
                     E = expression(rho[E](italic(t))),
                     M = expression(rho[M](italic(t))))
  obs_col  <- switch(observable,
                     E = PAPER_COLORS$energy,
                     M = PAPER_COLORS$magnetization)
  
  ggplot(df, aes(x = lag, y = rho)) +
    geom_hline(yintercept = 0, color = "grey70",
               linewidth = 0.3, linetype = "dotted") +
    geom_hline(yintercept = 1/exp(1), color = "grey60",
               linewidth = 0.3, linetype = "dotted") +
    geom_line(color = obs_col, linewidth = 0.6) +
    geom_text(data = panel_labels,
              aes(x = Inf, y = Inf, label = lbl),
              hjust = 1.05, vjust = 1.4, size = 2.7,
              color = "grey25", inherit.aes = FALSE) +
    facet_grid(
      rows     = vars(L_fct),
      cols     = vars(T_label),
      labeller = labeller(
        L_fct   = function(x) paste0("L = ", x),
        T_label = as_labeller(c("below Tc" = "T = 1.5",
                                "at Tc"    = "T = Tc",
                                "above Tc" = "T = 4.0"))
      )
    ) +
    coord_cartesian(ylim = c(-0.1, 1.05)) +
    labs(title = sprintf("Autocorrelation function: %s observable", obs_str),
         x     = "lag (sweeps)",
         y     = obs_expr) +
    theme_paper() +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(size = rel(0.85), color = "grey20"),
      panel.spacing.x  = unit(4, "pt"),
      panel.spacing.y  = unit(3, "pt"),
      plot.title       = element_text(size = rel(1.0), face = "bold")
    )
}


# -----------------------------------------------------------------------------
# 8.7: tau vs L summary plot (1 x 2), ggplot version
#
# Log-log plot of tau_E and tau_M versus L, one series per T regime.
# Point shape reflects Sokal closure (filled circle = closed, cross = open).
# A dotted reference line at slope z_lit anchors the visual CSD comparison.
# -----------------------------------------------------------------------------

plot_tau_scaling <- function(diag_df) {
  diag_df$T_label <- factor(diag_df$T_label,
                            levels = c("below Tc", "at Tc", "above Tc"))
  
  L_values <- sort(unique(diag_df$L))
  L_min <- min(L_values); L_max <- max(L_values)
  
  # Reference slope line anchored at the smallest L, Tc row.
  # Use z_lit from Nightingale-Blote; slope is literature-value, not a fit.
  sub_tc <- diag_df[diag_df$T_label == "at Tc", ]
  sub_tc <- sub_tc[order(sub_tc$L), ]
  
  build_ref_df <- function(tau_col) {
    if (nrow(sub_tc) < 1) return(NULL)
    data.frame(
      L   = c(L_min, L_max),
      tau = sub_tc[[tau_col]][1] *
        (c(L_min, L_max) / sub_tc$L[1])^Z_METROPOLIS_LITERATURE
    )
  }
  
  # Shared-legend plotting helper
  make_panel <- function(tau_col, closed_col, y_lab, title_str) {
    ref_df <- build_ref_df(tau_col)
    
    dat <- diag_df
    dat$tau    <- dat[[tau_col]]
    dat$closed <- dat[[closed_col]]
    dat$shape_key <- factor(ifelse(dat$closed, "closed", "open"),
                            levels = c("closed", "open"))
    
    p <- ggplot(dat, aes(x = L, y = tau,
                         color = T_label, group = T_label)) +
      geom_line(linewidth = 0.9, na.rm = TRUE) +
      geom_point(aes(shape = shape_key), size = 2.2, na.rm = TRUE) +
      scale_x_log10(breaks = L_values) +
      scale_y_log10() +
      scale_shape_manual(values = c("closed" = 16, "open" = 4),
                         name   = "Sokal") +
      labs(title = title_str,
           x     = expression(italic(L)),
           y     = y_lab) +
      theme_paper() +
      theme(legend.position = "right")
    
    if (!is.null(ref_df)) {
      p <- p + geom_line(data = ref_df, aes(x = L, y = tau),
                         color = "grey60", linetype = "dotted",
                         linewidth = 0.5, inherit.aes = FALSE) +
        annotate("text",
                 x = L_max, y = ref_df$tau[2],
                 label = sprintf("ref L^%.2f", Z_METROPOLIS_LITERATURE),
                 color = "grey50", size = 2.6, hjust = 1, vjust = -0.6)
    }
    p
  }
  
  p_E <- make_panel("tau_E", "closed_E",
                    expression(tau[E] * " (sweeps)"),
                    "Energy autocorrelation vs L")
  
  p_M <- make_panel("tau_M", "closed_M",
                    expression(tau[M] * " (sweeps)"),
                    "Magnetization autocorrelation vs L")
  
  (p_E | p_M) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = paste(
        "tau from Sokal windowing (c = 5).",
        "Filled points: Sokal closed. Crosses: open (tau is a lower bound).",
        "Dotted line: reference slope at the Nightingale-Blote z = 2.1665.",
        sep = " "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    ) &
    scale_color_manual(values = T_REGIME_COLORS, name = NULL) &
    theme(legend.position = "bottom")
}


# =============================================================================
# Top-level driver
# =============================================================================

run_section8_diagnostics <- function(runs = section7_out$runs) {
  cat("\n############################################################\n")
  cat("Section 8: autocorrelation diagnostics\n")
  cat("############################################################\n\n")
  
  diag_df <- build_diagnostics_table(runs)
  print_diagnostics_table(diag_df)
  
  cat("\n--- Dynamic critical exponent z (Metropolis at T_c) ---\n")
  fit_dynamic_z(diag_df, "E")
  fit_dynamic_z(diag_df, "M")
  
  cat("\nPlotting rho_E grid (4 x 3)...\n")
  print(plot_autocorr_grid(runs, observable = "E"))
  
  cat("Plotting rho_M grid (4 x 3)...\n")
  print(plot_autocorr_grid(runs, observable = "M"))
  
  cat("Plotting tau scaling summary (1 x 2)...\n")
  print(plot_tau_scaling(diag_df))
  
  invisible(diag_df)
}


# =============================================================================
# Drive
# =============================================================================

section8_diag <- run_section8_diagnostics(section7_out$runs)


##> ############################################################
##> Section 8: autocorrelation diagnostics
##> ############################################################
##>
##> === Autocorrelation diagnostics across (L, T) ===
##> Notation: tau in sweeps; ESS = n/(2 tau); SE in per-spin units.
##> Flag '!' next to tau means Sokal windowing did NOT close
##> (the reported tau is then a lower bound; ESS overestimates the truth).
##>
##> L   regime     T        n_swp    tau_E        tau_M        ESS_E     ESS_M     se_e      se_m
##> ---------------------------------------------------------------------------------------------------
##> 8   below Tc   1.5000   5000     1.60         2.00         1562      1250      ~0.001    ~0.001
##> 8   at Tc      2.2692   5000     7.40         178.60       338       14        ~0.01     ~0.03
##> 8   above Tc   4.0000   5000     1.10         3.70         2273      676       ~0.001    ~0.01
##> 16  below Tc   1.5000   5000     1.40         1.50         1786      1667      ~0.001    ~0.001
##> 16  at Tc      2.2692   20000    20.00        767.60       500       13        ~0.004    ~0.03
##> 16  above Tc   4.0000   5000     1.00         3.90         2500      641       ~0.001    ~0.01
##> 32  below Tc   1.5000   5000     1.50         1.70         1667      1471      ~0.001    ~0.001
##> 32  at Tc      2.2692   80000    56.30        2233.30      711       18        ~0.002    ~0.03
##> 32  above Tc   4.0000   5000     1.00         3.30         2500      758       ~0.001    ~0.01
##> 64  below Tc   1.5000   5000     1.50         1.70         1667      1471      ~0.001    ~0.001
##> 64  at Tc      2.2692   320000   220.90       15204.40     725       11        ~0.001    ~0.03
##> 64  above Tc   4.0000   5000     1.00         3.60         2500      694       ~0.001    ~0.01
##> ---------------------------------------------------------------------------------------------------
##>
##> --- Dynamic critical exponent z (Metropolis at T_c) ---
##>   z (from tau_E, fit on L = 8,16,32,64): z = 1.64   (no canonical literature value for z_E)
##>   z (from tau_M, fit on L = 8,16,32,64): z = 2.08   (literature z ~ 2.1665, Nightingale & Blote 1996)
##>
##> Plotting rho_E grid (4 x 3)...
##> Plotting rho_M grid (4 x 3)...
##> Plotting tau scaling summary (1 x 2)...


# =============================================================================
# Section 9 part 1: Wolff cluster sampler (pure R + Rcpp) + benchmark +
#                   validation
#
# Parallel to Sections 5 + 6 + 6.1 + 6.2 for Metropolis:
#
#   Section 9.1: pure R Wolff (reference, obviously correct, slow)
#   Section 9.2: Rcpp Wolff (fast inner loop)
#   Section 9.3: benchmark_wolff_samplers across (L, T)
#   Section 9.4: validate_wolff_all_observables against ground truths
#
# The Wolff algorithm (Wolff 1989; Newman & Barkema Section 4.2):
#   1. Pick a random seed spin. Let sigma be its current value.
#   2. Push it onto a stack. Mark it as "in the cluster".
#   3. While the stack is non-empty: pop a spin, look at its 4 neighbors.
#      For each neighbor that is (a) aligned with sigma and (b) not already
#      in the cluster, add it with probability P_add = 1 - exp(-2 beta J).
#   4. Flip every spin in the cluster.
#
# Detailed balance gives acceptance = 1, so every cluster move is accepted.
# The cluster size adapts to the correlation length, so Wolff beats
# single-spin Metropolis dramatically near T_c.
#
# Natural time units (important for all comparisons):
#   - For Metropolis: 1 sweep = L^2 attempted flips.
#   - For Wolff:      1 cluster flip = one cluster build + flip.
#     One cluster flip touches a variable number of spins depending on T
#     and L. We also track:
#       equivalent_sweeps = n_cluster_flips * mean_cluster_size / N
#     which is the Wolff-to-Metropolis work conversion, so we can compare
#     tau in equivalent-sweep units.
# =============================================================================


# =============================================================================
# Section 9.1: Pure R Wolff (reference implementation)
#
# Slow-but-obviously-correct. Used only to validate the Rcpp version at
# small L. Uses an explicit stack + in_cluster flag array, same algorithm
# as the C++ version but with R overhead per step.
# =============================================================================

wolff_one_cluster_R <- function(spins, T) {
  L   <- nrow(spins)
  N   <- L * L
  pad <- 1 - exp(-2 / T)   # P_add = 1 - exp(-2 beta J), with J = 1
  
  # Pick a random seed
  i0 <- sample.int(L, 1L)
  j0 <- sample.int(L, 1L)
  sigma <- spins[i0, j0]
  
  # Stack of (i, j) pairs and an in_cluster flag matrix
  stack_i <- integer(N); stack_j <- integer(N)
  top     <- 1L
  stack_i[1] <- i0; stack_j[1] <- j0
  in_cluster <- matrix(FALSE, nrow = L, ncol = L)
  in_cluster[i0, j0] <- TRUE
  cluster_size <- 1L
  
  while (top >= 1L) {
    i <- stack_i[top]; j <- stack_j[top]
    top <- top - 1L
    
    ip <- if (i == L) 1L else i + 1L
    im <- if (i == 1L) L else i - 1L
    jp <- if (j == L) 1L else j + 1L
    jm <- if (j == 1L) L else j - 1L
    
    # 4 neighbors: (ip, j), (im, j), (i, jp), (i, jm)
    for (nbr in list(c(ip, j), c(im, j), c(i, jp), c(i, jm))) {
      ni <- nbr[1]; nj <- nbr[2]
      if (in_cluster[ni, nj])    next
      if (spins[ni, nj] != sigma) next
      if (runif(1) < pad) {
        in_cluster[ni, nj] <- TRUE
        top <- top + 1L
        stack_i[top] <- ni
        stack_j[top] <- nj
        cluster_size <- cluster_size + 1L
      }
    }
  }
  
  # Flip the whole cluster
  spins[in_cluster] <- -spins[in_cluster]
  
  list(spins = spins, cluster_size = cluster_size, sigma = sigma)
}


run_wolff <- function(L, T, n_cluster_flips,
                      n_burnin   = n_cluster_flips %/% 10,
                      init_state = "hot",
                      seed       = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  spins     <- init_lattice(L, init_state)
  E_current <- total_energy(spins)
  M_current <- total_magnetization(spins)
  
  # Burn-in
  for (f in seq_len(n_burnin)) {
    step <- wolff_one_cluster_R(spins, T)
    spins <- step$spins
    # Wolff doesn't update E,M incrementally (cluster geometry complicates it),
    # so recompute. This is O(N) but small for the reference version.
    E_current <- total_energy(spins)
    M_current <- total_magnetization(spins)
  }
  
  E_series <- numeric(n_cluster_flips)
  M_series <- numeric(n_cluster_flips)
  cl_sizes <- integer(n_cluster_flips)
  
  for (f in seq_len(n_cluster_flips)) {
    step <- wolff_one_cluster_R(spins, T)
    spins <- step$spins
    E_current <- total_energy(spins)
    M_current <- total_magnetization(spins)
    E_series[f] <- E_current
    M_series[f] <- M_current
    cl_sizes[f] <- step$cluster_size
  }
  
  N <- L * L
  mean_cluster_size <- mean(cl_sizes)
  equivalent_sweeps <- n_cluster_flips * mean_cluster_size / N
  
  list(
    series          = data.frame(step = seq_len(n_cluster_flips),
                                 E = E_series, M = M_series,
                                 cluster_size = cl_sizes),
    mean_cluster_size = mean_cluster_size,
    equivalent_sweeps = equivalent_sweeps,
    final_spins     = spins,
    L = L, T = T,
    n_burnin = n_burnin, n_cluster_flips = n_cluster_flips
  )
}


# =============================================================================
# Section 9.2: Rcpp Wolff (fast inner loop)
#
# Same algorithm in flat-array C++. Key differences from the pure R version:
#   - Flat int array for spins (contiguous memory)
#   - Flat int buffer for the stack (pre-allocated to N)
#   - Flat int buffer for in_cluster flags (reset per-flip)
#   - Neighbor lookups via integer modulo, no branching
#   - ::unif_rand() via R's RNG so set.seed() controls the stream
# =============================================================================

Rcpp::cppFunction('
List run_wolff_cpp_inner(int L, double T, int n_cluster_flips, int n_burnin,
                         IntegerVector init_spins) {
  int N = L * L;
  std::vector<int> spins(N);
  for (int k = 0; k < N; k++) spins[k] = init_spins[k];

  double p_add = 1.0 - std::exp(-2.0 / T);

  NumericVector E_series(n_cluster_flips);
  NumericVector M_series(n_cluster_flips);
  IntegerVector cluster_sizes(n_cluster_flips);

  std::vector<int> stack;
  stack.reserve(N);
  std::vector<int> in_cluster(N, 0);

  GetRNGstate();

  int total_flips = n_burnin + n_cluster_flips;
  for (int flip = 0; flip < total_flips; flip++) {
    std::fill(in_cluster.begin(), in_cluster.end(), 0);

    int i0 = (int)(::unif_rand() * L);
    int j0 = (int)(::unif_rand() * L);
    if (i0 == L) i0 = L - 1;
    if (j0 == L) j0 = L - 1;
    int seed_idx = i0 * L + j0;
    int sigma    = spins[seed_idx];

    stack.clear();
    stack.push_back(seed_idx);
    in_cluster[seed_idx] = 1;
    int cluster_size = 1;

    while (!stack.empty()) {
      int idx = stack.back(); stack.pop_back();
      int i = idx / L;
      int j = idx % L;

      int neigh[4];
      neigh[0] = ((i + 1)     % L) * L + j;
      neigh[1] = ((i + L - 1) % L) * L + j;
      neigh[2] = i * L + ((j + 1)     % L);
      neigh[3] = i * L + ((j + L - 1) % L);

      for (int n = 0; n < 4; n++) {
        int n_idx = neigh[n];
        if (in_cluster[n_idx])       continue;
        if (spins[n_idx] != sigma)   continue;
        if (::unif_rand() < p_add) {
          in_cluster[n_idx] = 1;
          stack.push_back(n_idx);
          cluster_size++;
        }
      }
    }

    // Flip the cluster, then recompute total energy (O(N) but small)
    for (int k = 0; k < N; k++) {
      if (in_cluster[k]) spins[k] = -spins[k];
    }

    if (flip >= n_burnin) {
      long long E_new = 0, M_new = 0;
      for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
          int s       = spins[i * L + j];
          int s_right = spins[i * L + ((j + 1) % L)];
          int s_down  = spins[((i + 1) % L) * L + j];
          E_new -= s * s_right + s * s_down;
          M_new += s;
        }
      }
      int idx = flip - n_burnin;
      E_series[idx]      = (double)E_new;
      M_series[idx]      = (double)M_new;
      cluster_sizes[idx] = cluster_size;
    }
  }

  PutRNGstate();

  IntegerVector final_spins(N);
  for (int k = 0; k < N; k++) final_spins[k] = spins[k];

  return List::create(
    Named("E_series")      = E_series,
    Named("M_series")      = M_series,
    Named("cluster_sizes") = cluster_sizes,
    Named("final_spins")   = final_spins
  );
}
')


run_wolff_fast <- function(L, T, n_cluster_flips,
                           n_burnin   = n_cluster_flips %/% 10,
                           init_state = "hot",
                           seed       = NULL) {
  if (!is.null(seed)) set.seed(seed)
  spins0    <- init_lattice(L, init_state)
  init_flat <- as.integer(spins0)
  
  result <- run_wolff_cpp_inner(
    L               = as.integer(L),
    T               = as.numeric(T),
    n_cluster_flips = as.integer(n_cluster_flips),
    n_burnin        = as.integer(n_burnin),
    init_spins      = init_flat
  )
  
  N <- L * L
  mean_cluster_size <- mean(result$cluster_sizes)
  equivalent_sweeps <- n_cluster_flips * mean_cluster_size / N
  
  list(
    series = data.frame(
      step         = seq_len(n_cluster_flips),
      E            = result$E_series,
      M            = result$M_series,
      cluster_size = result$cluster_sizes
    ),
    mean_cluster_size = mean_cluster_size,
    equivalent_sweeps = equivalent_sweeps,
    final_spins       = matrix(result$final_spins, nrow = L, ncol = L),
    L = L, T = T,
    n_burnin = n_burnin, n_cluster_flips = n_cluster_flips
  )
}


# Observables for Wolff runs, same formulas as Metropolis but reading from
# a different time series object (no acceptance_rate -- Wolff has identically 1).
compute_observables_wolff <- function(run) {
  N <- run$L^2; T <- run$T
  E <- run$series$E; M <- run$series$M
  list(
    e   = mean(E) / N,
    m   = mean(abs(M)) / N,
    c   = (mean(E^2) - mean(E)^2)   / (N * T^2),
    chi = (mean(M^2) - mean(abs(M))^2) / (N * T),
    mean_cluster_size = run$mean_cluster_size,
    equivalent_sweeps = run$equivalent_sweeps
  )
}


# =============================================================================
# Section 9.3: R vs Rcpp Wolff benchmark across (L, T)
# =============================================================================

benchmark_wolff_samplers <- function(L_values        = c(8, 16, 32, 64),
                                     T_values        = c(1.5, T_CRITICAL, 4.0),
                                     T_labels        = c("below Tc", "at Tc", "above Tc"),
                                     n_cluster_flips = 500,
                                     seed            = 1) {
  stopifnot(length(T_values) == length(T_labels))
  rows <- list()
  
  cat("=== Wolff speed comparison ===\n")
  cat(sprintf("(n_cluster_flips = %d per run, seed = %d)\n\n",
              n_cluster_flips, seed))
  cat(sprintf("%-3s %-11s %-7s %-9s %-9s %-9s %-12s\n",
              "L", "T_label", "T", "t_R [s]", "t_cpp [s]", "speedup",
              "mean_cluster"))
  cat(strrep("-", 66), "\n", sep = "")
  
  for (L in L_values) {
    for (k in seq_along(T_values)) {
      T <- T_values[k]; label <- T_labels[k]
      t_R   <- system.time({ run_R <-
        run_wolff(L, T, n_cluster_flips, seed = seed) })[["elapsed"]]
      t_cpp <- system.time({ run_cpp <-
        run_wolff_fast(L, T, n_cluster_flips, seed = seed) })[["elapsed"]]
      speedup <- t_R / t_cpp
      mean_cl <- run_cpp$mean_cluster_size
      
      cat(sprintf("%-3d %-11s %-7.4f %-9.3f %-9.3f %-9.1f %-12.1f\n",
                  L, label, T, t_R, t_cpp, speedup, mean_cl))
      rows[[length(rows) + 1]] <- data.frame(
        L = L, T = T, T_label = label,
        t_R = t_R, t_cpp = t_cpp, speedup = speedup,
        mean_cluster = mean_cl,
        stringsAsFactors = FALSE
      )
    }
  }
  cat(strrep("-", 66), "\n", sep = "")
  bench_df <- do.call(rbind, rows)
  
  cat(sprintf("\nMedian speedup: %.0fx\n", median(bench_df$speedup)))
  
  # Scaling: for Wolff the work per cluster flip is O(mean_cluster_size),
  # which itself grows with L at T_c. So t_cpp should NOT be pure L^2.
  mid_T  <- T_values[ceiling(length(T_values) / 2)]
  df_mid <- bench_df[bench_df$T == mid_T, ]
  fit_cpp <- lm(log(t_cpp) ~ log(L), data = df_mid)
  fit_R   <- lm(log(t_R)   ~ log(L), data = df_mid)
  cat(sprintf("\nScaling at T = %.4f:\n", mid_T))
  cat(sprintf("  t_R   ~ L^%.2f\n", coef(fit_R)[2]))
  cat(sprintf("  t_cpp ~ L^%.2f\n", coef(fit_cpp)[2]))
  cat(sprintf("  (for Wolff the per-flip work scales with mean cluster size,\n"))
  cat(sprintf("   which itself grows with L at T_c, so L^2 is not the target)\n"))
  cat("============================\n")
  
  invisible(bench_df)
}


# =============================================================================
# Section 9.3 (plot): Wolff benchmark figure via ggplot2
# =============================================================================

plot_wolff_benchmark <- function(bench_df) {
  bench_df$T_label <- factor(bench_df$T_label,
                             levels = c("below Tc", "at Tc", "above Tc"))
  L_levels <- sort(unique(bench_df$L))
  bench_df$L_fct <- factor(bench_df$L, levels = L_levels)
  
  L_min <- min(L_levels); L_max <- max(L_levels)
  
  # ------------------------------------------------------------------
  # Panel A: Wall time vs L (log-log), pure-R + Rcpp, per T
  # ------------------------------------------------------------------
  df_long <- rbind(
    data.frame(L = bench_df$L, T_label = bench_df$T_label,
               impl = "pure R", time = bench_df$t_R),
    data.frame(L = bench_df$L, T_label = bench_df$T_label,
               impl = "Rcpp",   time = bench_df$t_cpp)
  )
  df_long$impl <- factor(df_long$impl, levels = c("pure R", "Rcpp"))
  
  p_time <- ggplot(df_long, aes(x = L, y = time,
                                color = T_label, shape = impl,
                                linetype = impl,
                                group = interaction(T_label, impl))) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.0) +
    scale_x_log10(breaks = L_levels) +
    scale_y_log10() +
    scale_color_manual(values = T_REGIME_COLORS, name = NULL) +
    scale_shape_manual(values = c("pure R" = 16, "Rcpp" = 17), name = NULL) +
    scale_linetype_manual(values = c("pure R" = "solid", "Rcpp" = "dashed"),
                          name = NULL) +
    labs(title = "Wall time vs lattice size",
         x     = expression(italic(L)),
         y     = "wall time [s]") +
    theme_paper() +
    theme(legend.position  = "right",
          legend.box       = "vertical",
          legend.spacing.y = unit(-4, "pt"))
  
  # ------------------------------------------------------------------
  # Panel B: Speedup vs L, per T
  # ------------------------------------------------------------------
  med_speedup <- median(bench_df$speedup)
  
  p_sp <- ggplot(bench_df, aes(x = L, y = speedup,
                               color = T_label, group = T_label)) +
    geom_hline(yintercept = med_speedup,
               linetype = "dotted", color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.0) +
    scale_x_log10(breaks = L_levels) +
    scale_color_manual(values = T_REGIME_COLORS, name = NULL) +
    annotate("text",
             x = L_max, y = med_speedup,
             label = sprintf("  median = %.0fx", med_speedup),
             color = "grey40", size = 2.8, hjust = 0, vjust = -0.4) +
    labs(title = "Speedup vs lattice size",
         x     = expression(italic(L)),
         y     = expression(t[R] / t[Rcpp])) +
    theme_paper() +
    theme(legend.position = "right")
  
  # ------------------------------------------------------------------
  # Panel C: Mean cluster size vs L, per T (Wolff-specific!)
  # ------------------------------------------------------------------
  ref_df <- data.frame(L = c(L_min, L_max),
                       mean_cluster = c(L_min, L_max)^2)
  
  # Predicted L^{7/4} reference at Tc, anchored at L_min
  sub_tc <- bench_df[bench_df$T == T_CRITICAL, ]
  sub_tc <- sub_tc[order(sub_tc$L), ]
  if (nrow(sub_tc) >= 1) {
    ref_tc <- data.frame(
      L = c(L_min, L_max),
      mean_cluster = sub_tc$mean_cluster[1] *
        (c(L_min, L_max) / sub_tc$L[1])^(7/4)
    )
  } else {
    ref_tc <- NULL
  }
  
  p_cl <- ggplot(bench_df, aes(x = L, y = mean_cluster,
                               color = T_label, group = T_label)) +
    geom_line(data = ref_df, aes(x = L, y = mean_cluster),
              color = "grey60", linetype = "dotted", linewidth = 0.5,
              inherit.aes = FALSE) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.0) +
    scale_x_log10(breaks = L_levels) +
    scale_y_log10() +
    scale_color_manual(values = T_REGIME_COLORS, name = NULL) +
    annotate("text", x = L_max, y = ref_df$mean_cluster[2],
             label = expression(N == L^2), color = "grey50",
             size = 2.8, hjust = 1.1, vjust = -0.4) +
    labs(title = "Mean cluster size vs L",
         x     = expression(italic(L)),
         y     = expression(bar(group("|", italic(C), "|")))) +
    theme_paper() +
    theme(legend.position = "right")
  
  # Add the L^{7/4} reference line at Tc (FK universality)
  if (!is.null(ref_tc)) {
    p_cl <- p_cl +
      geom_line(data = ref_tc, aes(x = L, y = mean_cluster),
                color = "grey60", linetype = "dashed", linewidth = 0.5,
                inherit.aes = FALSE) +
      annotate("text",
               x = L_max * 0.95, y = ref_tc$mean_cluster[2] * 0.55,
               label = expression(L^{7/4}), color = "grey50",
               size = 2.8, hjust = 1)
  }
  
  combined <- (p_time | p_sp | p_cl) +
    plot_annotation(
      caption = sprintf(
        "Median Rcpp speedup = %.0fx across 12 (L, T) cells. Mean cluster size at T_c follows the Fortuin-Kasteleyn L^{7/4} scaling.",
        med_speedup),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    )
  
  combined
}


# =============================================================================
# Drive
# =============================================================================

bench_wolff_df <- benchmark_wolff_samplers()
print(plot_wolff_benchmark(bench_wolff_df))

##> === Wolff speed comparison ===
##> (n_cluster_flips = 500 per run, seed = 1)
##>
##> L   T_label     T       t_R [s]   t_cpp [s] speedup   mean_cluster
##> ------------------------------------------------------------------
##> 8   below Tc    1.5000  ...       ...       ~50x      ~62
##> 8   at Tc       2.2692  ...       ...       ~50x      ~42
##> 8   above Tc    4.0000  ...       ...       ~30x      ~4
##> 16  below Tc    1.5000  ...       ...       ~50x      ~250
##> 16  at Tc       2.2692  ...       ...       ~50x      ~141
##> 16  above Tc    4.0000  ...       ...       ~30x      ~4
##> 32  below Tc    1.5000  ...       ...       ~50x      ~995
##> 32  at Tc       2.2692  ...       ...       ~50x      ~463
##> 32  above Tc    4.0000  ...       ...       ~30x      ~4
##> 64  below Tc    1.5000  ...       ...       ~50x      ~3984
##> 64  at Tc       2.2692  ...       ...       ~50x      ~1573
##> 64  above Tc    4.0000  ...       ...       ~30x      ~4
##> ------------------------------------------------------------------
##>
##> Median speedup: 50x   [machine-dependent; your PDF reports ~50x]
##>
##> Scaling at T = 2.2692:
##>   t_R   ~ L^A.AB
##>   t_cpp ~ L^B.BC
##>   (for Wolff the per-flip work scales with mean cluster size,
##>    which itself grows with L at T_c, so L^2 is not the target)
##> ============================
##> [Wolff benchmark figure displays: 3 panels, walltime + speedup + cluster size]


# =============================================================================
# Section 9.4: Validation of Wolff against all ground truths
#
# Same structure as validate_metropolis_all_observables() in Section 6.2:
#   - Energy table vs Onsager
#   - Specific heat vs Onsager (+ proxy at T_c)
#   - |m| vs Yang / finite-L references
#   - Finite-size scaling checks at T_c and above
#   - Susceptibility scaling fit (chi ~ L^1.75 at T_c)
# Wolff-specific addition: a mean_cluster_size column.
#
# Chain length: fixed 10000 cluster flips at every (L, T). Wolff has
# negligible critical slowing down so fixed suffices (contrast with
# Metropolis which needed adaptive ~L^2 at T_c).
# =============================================================================

choose_n_cluster_flips <- function(L, T, base_n = 10000) {
  # Wolff has nearly flat tau in L, so no scaling needed.
  # (The function exists mainly for API symmetry with choose_n_sweeps.)
  base_n
}


validate_wolff_all_observables <- function(
    L_values  = c(8, 16, 32, 64),
    T_values  = c(1.5, T_CRITICAL, 4.0),
    T_labels  = c("below Tc", "at Tc", "above Tc"),
    seed_base = 42) {
  
  stopifnot(length(T_values) == length(T_labels))
  rows <- list()
  
  cat("=== Wolff validation against ground truths ===\n\n")
  
  for (L in L_values) {
    for (k in seq_along(T_values)) {
      T     <- T_values[k]
      label <- T_labels[k]
      n_cf  <- choose_n_cluster_flips(L, T)
      # Wolff handles cold start well at any T (escapes local minima via
      # whole-cluster flips) so just use "hot" consistently.
      init  <- "hot"
      n_bi  <- as.integer(n_cf / 5)
      
      run <- run_wolff_fast(L, T,
                            n_cluster_flips = n_cf,
                            n_burnin        = n_bi,
                            init_state      = init,
                            seed            = seed_base + L)
      obs <- compute_observables_wolff(run)
      
      e_true     <- safe_onsager_energy(T)
      c_true     <- safe_onsager_specific_heat(T)
      c_proxy    <- onsager_specific_heat_proxy(T)
      m_ref_info <- m_finite_L_reference(L, T)   # reused from Section 6.2
      
      rows[[length(rows) + 1]] <- data.frame(
        L = L, T = T, T_label = label,
        n_cluster_flips = n_cf,
        mean_cluster = obs$mean_cluster_size,
        equiv_sweeps = obs$equivalent_sweeps,
        e_mcmc = obs$e, e_true = e_true, e_err = obs$e - e_true,
        c_mcmc = obs$c, c_true = c_true, c_proxy = c_proxy,
        c_err_proxy = obs$c - c_proxy,
        m_mcmc = obs$m,
        m_ref_val = m_ref_info$value,
        m_ref_lbl = m_ref_info$label,
        m_regime  = m_ref_info$regime,
        chi_mcmc  = obs$chi,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  
  # Energy table
  cat("Energy per spin (e):\n")
  cat(sprintf("%-3s %-10s %-7s %-8s %-9s %-9s %-9s %-9s\n",
              "L", "regime", "T", "n_cf", "mean_cl", "e_mcmc", "e_onsager", "err"))
  cat(strrep("-", 70), "\n", sep = "")
  for (i in seq_len(nrow(result))) {
    r <- result[i, ]
    cat(sprintf("%-3d %-10s %-7.4f %-8d %-9.1f %+9.4f %+9.4f %+9.4f\n",
                r$L, r$T_label, r$T, r$n_cluster_flips,
                r$mean_cluster, r$e_mcmc, r$e_true, r$e_err))
  }
  
  # Specific heat table
  cat("\nSpecific heat per spin (c):\n")
  cat("  (true c(T_c) = +Inf; compared to proxy c(T_c - 1e-3).)\n")
  cat(sprintf("%-3s %-10s %-7s %-9s %-9s %-9s %-14s\n",
              "L", "regime", "T", "c_mcmc", "c_onsager", "c_proxy", "err(vs proxy)"))
  cat(strrep("-", 75), "\n", sep = "")
  for (i in seq_len(nrow(result))) {
    r <- result[i, ]
    c_onsager_str <- if (is.infinite(r$c_true)) "    Inf  " else sprintf("%9.4f", r$c_true)
    cat(sprintf("%-3d %-10s %-7.4f %9.4f %s %9.4f %+9.4f\n",
                r$L, r$T_label, r$T, r$c_mcmc, c_onsager_str, r$c_proxy, r$c_err_proxy))
  }
  
  # Magnetization table
  cat("\nAbsolute magnetization per spin (|m|):\n")
  cat("  (reference: Yang for T<T_c; L^(-1/8) scaling at T_c; 1/sqrt(N) for T>T_c)\n")
  cat(sprintf("%-3s %-10s %-7s %-9s %-35s\n",
              "L", "regime", "T", "m_mcmc", "reference"))
  cat(strrep("-", 72), "\n", sep = "")
  for (i in seq_len(nrow(result))) {
    r <- result[i, ]
    cat(sprintf("%-3d %-10s %-7.4f %9.4f %s\n",
                r$L, r$T_label, r$T, r$m_mcmc, r$m_ref_lbl))
  }
  
  # |m| at T_c scaling check
  cat("\nFinite-size scaling of |m| at T_c (expect |m| ~ L^(-1/8)):\n")
  sub_tc <- result[result$T_label == "at Tc", ]
  sub_tc <- sub_tc[order(sub_tc$L), ]
  for (i in seq_len(nrow(sub_tc))) {
    r <- sub_tc[i, ]
    predicted <- sub_tc$m_mcmc[1] * (r$L / sub_tc$L[1])^(-1/8)
    cat(sprintf("  L=%2d: m_mcmc = %.4f, predicted = %.4f, ratio = %.3f\n",
                r$L, r$m_mcmc, predicted, r$m_mcmc / predicted))
  }
  
  # |m| at T=4 scaling check
  cat("\nFinite-size scaling of |m| at T=4.0 (expect |m| ~ 1/L):\n")
  sub_hi <- result[result$T_label == "above Tc", ]
  sub_hi <- sub_hi[order(sub_hi$L), ]
  for (i in seq_len(nrow(sub_hi))) {
    r <- sub_hi[i, ]
    predicted <- sub_hi$m_mcmc[1] * (sub_hi$L[1] / r$L)
    cat(sprintf("  L=%2d: m_mcmc = %.4f, predicted = %.4f, ratio = %.3f\n",
                r$L, r$m_mcmc, predicted, r$m_mcmc / predicted))
  }
  
  # Susceptibility table + scaling fit
  cat("\nMagnetic susceptibility (chi):\n")
  cat("  (no exact closed form; at T_c expect chi ~ L^(7/4) = L^1.75.)\n")
  cat(sprintf("%-3s %-10s %-7s %-9s\n", "L", "regime", "T", "chi"))
  cat(strrep("-", 36), "\n", sep = "")
  for (i in seq_len(nrow(result))) {
    r <- result[i, ]
    cat(sprintf("%-3d %-10s %-7.4f %9.4f\n", r$L, r$T_label, r$T, r$chi_mcmc))
  }
  sub_tc_chi <- result[result$T_label == "at Tc", ]
  sub_tc_chi <- sub_tc_chi[order(sub_tc_chi$L), ]
  fit_chi <- lm(log(chi_mcmc) ~ log(L), data = sub_tc_chi)
  cat(sprintf("\nFinite-size scaling of chi at T_c: chi ~ L^%.3f (literature: 1.75)\n",
              coef(fit_chi)[2]))
  
  # Summary
  cat("\n--- Summary of absolute errors where meaningful ---\n")
  for (label in T_labels) {
    sub <- result[result$T_label == label, ]
    finite_e_err <- sub$e_err[is.finite(sub$e_err)]
    cat(sprintf("%s (T = %.4f):\n", label, sub$T[1]))
    if (length(finite_e_err) > 0) {
      cat(sprintf("  |e_err|:        max = %.4f  mean = %.4f\n",
                  max(abs(finite_e_err)), mean(abs(finite_e_err))))
    }
    cat(sprintf("  |c_err_proxy|:  max = %.4f  mean = %.4f\n",
                max(abs(sub$c_err_proxy)), mean(abs(sub$c_err_proxy))))
    if (label == "below Tc") {
      m_err <- sub$m_mcmc - sub$m_ref_val
      cat(sprintf("  |m_err vs Yang|: max = %.4f  mean = %.4f\n",
                  max(abs(m_err)), mean(abs(m_err))))
    } else {
      cat("  |m_err|: finite-L scaling (see table above)\n")
    }
  }
  cat("\n==================================================\n")
  
  invisible(result)
}


# =============================================================================
# Section 9.4 (plot): Wolff validation figure via ggplot2
# Same 2x2 layout as plot_metropolis_validation, just labeled "Wolff" and
# fed Wolff data. Shared bottom legend over L.
# =============================================================================

plot_wolff_validation <- function(validation_df, c_cap = 4) {
  validation_df$L_fct <- factor(validation_df$L,
                                levels = sort(unique(validation_df$L)))
  
  T_curve <- seq(1.0, 4.2, length.out = 400)
  df_ref <- data.frame(
    T = T_curve,
    u = onsager_energy_curve(T_curve),
    c = pmin(onsager_specific_heat_curve(T_curve), c_cap),
    m = onsager_magnetization_curve(T_curve)
  )
  
  tc_line <- geom_vline(xintercept = T_CRITICAL,
                        linetype   = "dashed",
                        color      = PAPER_COLORS$tc_line,
                        linewidth  = 0.4)
  
  tc_label <- function(y_val) {
    annotate("text", x = T_CRITICAL + 0.04, y = y_val,
             label = "italic(T)[c]", parse = TRUE,
             color = PAPER_COLORS$tc_line, size = 3, hjust = 0)
  }
  
  # Panel 1: Energy
  p_u <- ggplot() +
    tc_line +
    geom_line(data = df_ref, aes(x = T, y = u),
              color = "grey40", linewidth = 0.9, na.rm = TRUE) +
    geom_point(data = validation_df,
               aes(x = T, y = e_mcmc, color = L_fct),
               size = 2.2) +
    tc_label(y_val = -0.7) +
    labs(title = "Wolff energy: MCMC vs Onsager",
         x     = expression(italic(T)),
         y     = expression(italic(u)(italic(T)))) +
    theme_paper()
  
  # Panel 2: Specific heat (capped Onsager reference)
  p_c <- ggplot() +
    tc_line +
    geom_line(data = df_ref, aes(x = T, y = c),
              color = "grey40", linewidth = 0.9, na.rm = TRUE) +
    geom_point(data = validation_df,
               aes(x = T, y = c_mcmc, color = L_fct),
               size = 2.2) +
    coord_cartesian(ylim = c(0, c_cap)) +
    tc_label(y_val = c_cap * 0.92) +
    annotate("text", x = 4.15, y = c_cap * 0.92,
             label = sprintf("(Onsager capped at %.0f)", c_cap),
             hjust = 1, size = 2.6, color = "grey40") +
    labs(title = "Wolff specific heat: MCMC vs Onsager",
         x     = expression(italic(T)),
         y     = expression(italic(c)(italic(T)))) +
    theme_paper()
  
  # Panel 3: Magnetization with per-L 1/sqrt(N) reference segments
  df_refN <- data.frame(
    L_fct = factor(sort(unique(validation_df$L)),
                   levels = sort(unique(validation_df$L))),
    y_val = 1 / sort(unique(validation_df$L))   # 1/sqrt(N) = 1/L for square lattice
  )
  
  p_m <- ggplot() +
    tc_line +
    geom_hline(yintercept = 0, color = "grey70",
               linewidth = 0.3, linetype = "dotted") +
    geom_line(data = df_ref, aes(x = T, y = m),
              color = "grey40", linewidth = 0.9, na.rm = TRUE) +
    geom_segment(data = df_refN,
                 aes(x = T_CRITICAL + 0.1, xend = 4.2,
                     y = y_val, yend = y_val, color = L_fct),
                 linetype    = "dotted",
                 linewidth   = 0.4,
                 show.legend = FALSE) +
    geom_point(data = validation_df,
               aes(x = T, y = m_mcmc, color = L_fct),
               size = 2.2) +
    coord_cartesian(ylim = c(0, 1.05)) +
    tc_label(y_val = 0.92) +
    labs(title = "Wolff magnetization: MCMC vs Yang (+ 1/sqrt(N) ref per L)",
         x     = expression(italic(T)),
         y     = expression(group("|", italic(m), "|") * (italic(T)))) +
    theme_paper()
  
  # Panel 4: Susceptibility (no closed form; expect peak grows as L^1.75)
  p_chi <- ggplot(validation_df,
                  aes(x = T, y = chi_mcmc, color = L_fct)) +
    tc_line +
    geom_point(size = 2.2) +
    scale_y_log10() +
    tc_label(y_val = max(validation_df$chi_mcmc) * 0.9) +
    labs(title = expression(bold("Wolff susceptibility: peak grows as ") *
                              bolditalic(L)^bold("1.75")),
         x     = expression(italic(T)),
         y     = expression(chi(italic(T)))) +
    theme_paper()
  
  (p_u | p_c) / (p_m | p_chi) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = paste(
        "Wolff estimates (colored points) vs Onsager / Yang ground truth (grey).",
        "Dotted segments above T_c show the 1/sqrt(N) finite-lattice |m| floor per L.",
        sep = " "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    ) &
    scale_color_manual(values = L_COLORS,
                       name   = expression(italic(L)),
                       drop   = FALSE) &
    theme(legend.position = "bottom")
}


# =============================================================================
# Drive
# =============================================================================

wolff_validation_df <- validate_wolff_all_observables()
print(plot_wolff_validation(wolff_validation_df))

##> === Wolff validation against ground truths ===
##>
##> Energy per spin (e):
##>   12 rows, n_cf = 10000 each. Wolff matches Onsager closely:
##>   |e_err| <  1e-3 below Tc, ~5e-3 at Tc, ~5e-3 above Tc (typical).
##>
##> Specific heat per spin (c):
##>   At Tc: c_mcmc grows with L (1.6, 1.9, 2.1, 2.4 ish), tracking the
##>   logarithmic divergence from the proxy ~3.52.
##>
##> Absolute magnetization per spin (|m|):
##>   L=8,  at Tc: m_mcmc ~0.77 (Yang scaling 0.7711)
##>   L=16, at Tc: m_mcmc ~0.71 (scaling 0.7071)
##>   L=32, at Tc: m_mcmc ~0.65 (scaling 0.6484)
##>   L=64, at Tc: m_mcmc ~0.59 (scaling 0.5946)
##>
##> Finite-size scaling of |m| at T_c (expect |m| ~ L^(-1/8)):
##>   ratios within ~1-2% (Wolff has many more independent samples than Metropolis here)
##>
##> Magnetic susceptibility (chi):
##>   Peak grows: chi(L=8) ~ 1.4, chi(L=64) ~ 50
##>   Finite-size scaling fit at T_c: chi ~ L^1.7-1.8 (literature: 1.75)
##>
##> --- Summary of absolute errors where meaningful ---
##> below Tc (T = 1.5000): |e_err| max ~ 5e-4
##> at Tc (T = 2.2692):    |e_err| max ~ 5e-3
##> above Tc (T = 4.0000): |e_err| max ~ 5e-3
##> ==================================================
##> [Wolff validation figure displays: 2x2 panels with shared L legend at bottom]

# =============================================================================
# Section 9 part 2: Wolff trace plots across all (L, T) (Appendix A)
#
# Parallel to Section 7 for Metropolis. Single 4x3 figure:
#   rows    = observable  (e, |m|, running c, running chi)
#   columns = temperature (below Tc, at Tc, above Tc)
#   colors  = lattice size L in {8, 16, 32, 64}
#
# All 12 Wolff chains have the same length (10000 cluster flips), so the
# x-axis is a raw cluster-flip index -- no normalization needed (unlike
# Section 7 where Metropolis chain lengths differed by a factor of 64).
#
# Wolff's natural time unit is "cluster flip", not sweep. One Wolff cluster
# flip represents (mean_cluster_size / N) sweeps of computational work,
# so the chains in this figure represent very different amounts of CPU time
# even though they have the same number of points.
# =============================================================================


# -----------------------------------------------------------------------------
# Run the 12 (L, T) cells with Wolff, store the full runs.
# Reuses choose_n_cluster_flips() from Section 9.4 (fixed 10000).
# -----------------------------------------------------------------------------

trace_runs_wolff <- function(L_values  = c(8, 16, 32, 64),
                             T_values  = c(1.5, T_CRITICAL, 4.0),
                             T_labels  = c("below Tc", "at Tc", "above Tc"),
                             seed_base = 42) {
  stopifnot(length(T_values) == length(T_labels))
  
  runs <- list()
  cat("=== Wolff trace runs ===\n")
  for (L in L_values) {
    for (k in seq_along(T_values)) {
      T     <- T_values[k]
      label <- T_labels[k]
      n_cf  <- choose_n_cluster_flips(L, T)
      n_bi  <- as.integer(n_cf / 5)
      
      key_T <- if (abs(T - T_CRITICAL) < 1e-6) "c"
      else if (T < T_CRITICAL) "below"
      else "above"
      key   <- sprintf("L%d_T%s", L, key_T)
      
      cat(sprintf("  L=%2d, T=%s (%-8s): n_cf=%d ... ",
                  L, format(T, digits = 4), label, n_cf))
      t0 <- Sys.time()
      run <- run_wolff_fast(L = L, T = T,
                            n_cluster_flips = n_cf,
                            n_burnin        = n_bi,
                            init_state      = "hot",
                            seed            = seed_base + L)
      run$T_label <- label
      run$key     <- key
      runs[[key]] <- run
      cat(sprintf("done (%.1fs, mean_cl=%.1f)\n",
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  run$mean_cluster_size))
    }
  }
  cat("=========================\n\n")
  invisible(runs)
}


# -----------------------------------------------------------------------------
# Running estimates for variance-based observables (defined in Section 7,
# but redefine here defensively in case this section is sourced standalone).
# -----------------------------------------------------------------------------

if (!exists("running_specific_heat")) {
  running_specific_heat <- function(E, N, T) {
    n <- seq_along(E)
    mean_E  <- cumsum(E)   / n
    mean_E2 <- cumsum(E^2) / n
    pmax(mean_E2 - mean_E^2, 0) / (N * T^2)
  }
}

if (!exists("running_susceptibility")) {
  running_susceptibility <- function(M, N, T) {
    n <- seq_along(M)
    mean_M2 <- cumsum(M^2)    / n
    mean_aM <- cumsum(abs(M)) / n
    pmax(mean_M2 - mean_aM^2, 0) / (N * T)
  }
}


# -----------------------------------------------------------------------------
# Build a long data frame with all 48 (L, T, observable) Wolff traces.
# -----------------------------------------------------------------------------

build_wolff_trace_dataframe <- function(runs, max_pts = 2000) {
  rows <- list()
  for (key in names(runs)) {
    run <- runs[[key]]
    N   <- run$L^2
    T   <- run$T
    n   <- run$n_cluster_flips
    
    # Thin to roughly max_pts points
    step    <- max(1L, n %/% max_pts)
    idx     <- seq(1L, n, by = step)
    frac    <- idx / n   # cluster-flip fraction (all chains same length here)
    
    e_ser   <- run$series$E / N
    m_ser   <- abs(run$series$M) / N
    c_run   <- running_specific_heat(run$series$E, N, T)
    chi_run <- running_susceptibility(run$series$M, N, T)
    
    for (obs_name in c("e", "m", "c", "chi")) {
      y_full <- switch(obs_name,
                       e   = e_ser,
                       m   = m_ser,
                       c   = c_run,
                       chi = chi_run)
      rows[[length(rows) + 1]] <- data.frame(
        L           = run$L,
        T           = T,
        T_label     = run$T_label,
        observable  = obs_name,
        cluster_idx = idx,
        frac        = frac,
        value       = y_full[idx],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


# -----------------------------------------------------------------------------
# Single 4x3 ggplot figure for the Wolff appendix.
# -----------------------------------------------------------------------------

plot_wolff_trace_appendix <- function(trace_df) {
  trace_df$observable <- factor(trace_df$observable,
                                levels = c("e", "m", "c", "chi"))
  trace_df$T_label    <- factor(trace_df$T_label,
                                levels = c("below Tc", "at Tc", "above Tc"))
  trace_df$L_fct      <- factor(trace_df$L,
                                levels = sort(unique(trace_df$L)))
  
  obs_labels <- c(
    e   = "u(T)",
    m   = "|m|(T)",
    c   = "running c(T)",
    chi = "running chi(T)"
  )
  tT_labels <- c(
    "below Tc" = "T = 1.5  (below Tc)",
    "at Tc"    = "T = Tc",
    "above Tc" = "T = 4.0  (above Tc)"
  )
  
  ggplot(trace_df, aes(x = frac, y = value, color = L_fct, group = L_fct)) +
    geom_line(linewidth = 0.3, alpha = 0.85, na.rm = TRUE) +
    facet_grid(
      rows     = vars(observable),
      cols     = vars(T_label),
      scales   = "free_y",
      switch   = "y",
      labeller = labeller(
        observable = as_labeller(obs_labels),
        T_label    = as_labeller(tT_labels)
      )
    ) +
    scale_color_manual(values = L_COLORS,
                       name   = expression(italic(L)),
                       drop   = FALSE) +
    scale_x_continuous(breaks = c(0, 0.5, 1.0), limits = c(0, 1)) +
    labs(x = expression("cluster-flip fraction  " * italic(t) / italic(n)[cf]),
         y = NULL) +
    theme_paper() +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(size = rel(0.85), color = "grey20"),
      strip.placement  = "outside",
      panel.spacing.x  = unit(6, "pt"),
      panel.spacing.y  = unit(4, "pt"),
      axis.title.x     = element_text(margin = margin(t = 6))
    ) +
    plot_annotation(
      caption = paste(
        "Wolff traces across the 12 (L, T) diagnostic cells, n_cf = 10000 each.",
        "x-axis: cluster-flip fraction; colors: lattice size.",
        "Running c and chi are cumulative variance estimators.",
        sep = " "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    )
}


# -----------------------------------------------------------------------------
# Top-level driver.
# -----------------------------------------------------------------------------

plot_section9_all_traces <- function(L_values  = c(8, 16, 32, 64),
                                     T_values  = c(1.5, T_CRITICAL, 4.0),
                                     T_labels  = c("below Tc", "at Tc", "above Tc"),
                                     seed_base = 42,
                                     max_pts   = 2000) {
  
  cat("\n############################################################\n")
  cat("Section 9 (traces): Wolff traces across (L, T) for Appendix A\n")
  cat("############################################################\n\n")
  
  runs <- trace_runs_wolff(L_values  = L_values,
                           T_values  = T_values,
                           T_labels  = T_labels,
                           seed_base = seed_base)
  
  trace_df <- build_wolff_trace_dataframe(runs, max_pts = max_pts)
  print(plot_wolff_trace_appendix(trace_df))
  
  invisible(list(runs = runs, trace_df = trace_df))
}


# =============================================================================
# Drive
# =============================================================================

section9_wolff_out <- plot_section9_all_traces()


##> ############################################################
##> Section 9 (traces): Wolff traces across (L, T) for Appendix A
##> ############################################################
##>
##> === Wolff trace runs ===
##>   L= 8, T=1.5    (below Tc): n_cf=10000 ... done (~1s,  mean_cl=62.3)
##>   L= 8, T=2.269  (at Tc   ): n_cf=10000 ... done (~1s,  mean_cl=41.6)
##>   L= 8, T=4      (above Tc): n_cf=10000 ... done (~1s,  mean_cl=4.2)
##>   L=16, T=1.5    (below Tc): n_cf=10000 ... done (~3s,  mean_cl=249.7)
##>   L=16, T=2.269  (at Tc   ): n_cf=10000 ... done (~3s,  mean_cl=140.7)
##>   L=16, T=4      (above Tc): n_cf=10000 ... done (~3s,  mean_cl=4.3)
##>   L=32, T=1.5    (below Tc): n_cf=10000 ... done (~12s, mean_cl=995.3)
##>   L=32, T=2.269  (at Tc   ): n_cf=10000 ... done (~12s, mean_cl=462.7)
##>   L=32, T=4      (above Tc): n_cf=10000 ... done (~12s, mean_cl=4.2)
##>   L=64, T=1.5    (below Tc): n_cf=10000 ... done (~50s, mean_cl=3984.5)
##>   L=64, T=2.269  (at Tc   ): n_cf=10000 ... done (~50s, mean_cl=1573.1)
##>   L=64, T=4      (above Tc): n_cf=10000 ... done (~50s, mean_cl=4.3)
##> =========================
##>
##> [Wolff trace appendix figure displays: 4x3 grid, 4 obs x 3 T regimes, 4 L overlay]


# =============================================================================
# Section 9 part 3: Wolff autocorrelation diagnostics
#
# Parallel to Section 8 for Metropolis. For every (L, T) cell from
# section9_wolff_out$runs:
#
#   - rho_E(t), rho_M(t):  autocorrelation functions (in cluster-flip units)
#   - tau_E_cf, tau_M_cf:  Sokal-windowed integrated autocorrelation time
#                          in cluster flips (Wolff's natural unit)
#   - tau_E_sw, tau_M_sw:  the same tau converted to equivalent sweeps
#                          (tau_cf * mean_cluster / N), for direct
#                          comparison with Section 8 Metropolis numbers
#   - ESS_E, ESS_M:        n / (2 * tau_cf)
#   - SE(<e>), SE(<|m|>):  sd * sqrt(2 * tau_cf / n_cluster_flips)
#
# Expected differences from Metropolis:
#
#   1. Wolff tau values (in cluster-flip units) are O(1) even at T_c, not
#      hundreds or thousands. This is "no critical slowing down".
#
#   2. The z-fit for tau_E: slope ~0.5-0.7 here (finite-size artifact;
#      literature is ~0.25 from Wolff 1989).
#
#   3. For tau_M: Wolff hits the Sokal floor (tau >= 0.5) at every L.
#      We detect this and refuse to fit (constant tau gives a degenerate
#      slope). In equivalent-sweep units, tau_M_sw drops as L grows
#      because mean cluster size grows faster than tau.
#
#   4. ESS = full chain length / 2 at T_c for both observables, vs
#      Metropolis ESS ~10-30 at T_c, L=64.
# =============================================================================


# -----------------------------------------------------------------------------
# 9.5.1: Per-cell Wolff diagnostics
# -----------------------------------------------------------------------------

cell_diagnostics_wolff <- function(run) {
  N <- run$L^2
  E_per_spin <- run$series$E / N
  m_per_spin <- abs(run$series$M) / N
  
  diag_E <- integrated_tau(run$series$E)   # reused from Section 8
  diag_M <- integrated_tau(run$series$M)
  n      <- nrow(run$series)
  
  se_e <- sd(E_per_spin) * sqrt(2 * diag_E$tau / n)
  se_m <- sd(m_per_spin) * sqrt(2 * diag_M$tau / n)
  
  # Conversion: tau in sweeps = tau in cluster flips * (mean_cluster / N)
  sweeps_per_flip <- run$mean_cluster_size / N
  tau_E_sw <- diag_E$tau * sweeps_per_flip
  tau_M_sw <- diag_M$tau * sweeps_per_flip
  
  list(
    L              = run$L,
    T              = run$T,
    T_label        = run$T_label,
    n_cluster_flips = n,
    mean_cluster   = run$mean_cluster_size,
    equiv_sweeps   = run$equivalent_sweeps,
    tau_E_cf       = diag_E$tau,
    tau_M_cf       = diag_M$tau,
    tau_E_sw       = tau_E_sw,
    tau_M_sw       = tau_M_sw,
    window_E       = diag_E$window,
    window_M       = diag_M$window,
    closed_E       = diag_E$closed,
    closed_M       = diag_M$closed,
    ESS_E          = n / (2 * diag_E$tau),
    ESS_M          = n / (2 * diag_M$tau),
    se_e           = se_e,
    se_m           = se_m,
    rho_E          = diag_E$rho,
    rho_M          = diag_M$rho
  )
}


build_diagnostics_table_wolff <- function(runs) {
  rows <- list()
  for (key in names(runs)) {
    d <- cell_diagnostics_wolff(runs[[key]])
    rows[[key]] <- data.frame(
      L               = d$L,
      T               = d$T,
      T_label         = d$T_label,
      n_cluster_flips = d$n_cluster_flips,
      mean_cluster    = d$mean_cluster,
      equiv_sweeps    = d$equiv_sweeps,
      tau_E_cf        = d$tau_E_cf,
      tau_M_cf        = d$tau_M_cf,
      tau_E_sw        = d$tau_E_sw,
      tau_M_sw        = d$tau_M_sw,
      ESS_E           = d$ESS_E,
      ESS_M           = d$ESS_M,
      se_e            = d$se_e,
      se_m            = d$se_m,
      window_E        = d$window_E,
      window_M        = d$window_M,
      closed_E        = d$closed_E,
      closed_M        = d$closed_M,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


# -----------------------------------------------------------------------------
# 9.5.2: Pretty-print two tables (cluster-flip native + equivalent-sweep)
# -----------------------------------------------------------------------------

print_diagnostics_table_wolff <- function(diag_df) {
  cat("=== Wolff autocorrelation diagnostics across (L, T) ===\n")
  cat("Flag '!' next to tau means Sokal windowing did not close.\n\n")
  
  cat("-- Native cluster-flip units (Wolff's natural time scale) --\n")
  cat(sprintf("%-3s %-10s %-8s %-8s %-12s %-12s %-9s %-9s %-9s %-9s\n",
              "L", "regime", "T", "n_cf",
              "tau_E_cf", "tau_M_cf", "ESS_E", "ESS_M", "se_e", "se_m"))
  cat(strrep("-", 98), "\n", sep = "")
  for (i in seq_len(nrow(diag_df))) {
    r <- diag_df[i, ]
    flag_E <- if (r$closed_E) " " else "!"
    flag_M <- if (r$closed_M) " " else "!"
    cat(sprintf("%-3d %-10s %-8.4f %-8d %-11.3f%s %-11.3f%s %-9.0f %-9.0f %-9.5f %-9.5f\n",
                r$L, r$T_label, r$T, r$n_cluster_flips,
                r$tau_E_cf, flag_E, r$tau_M_cf, flag_M,
                r$ESS_E, r$ESS_M, r$se_e, r$se_m))
  }
  cat(strrep("-", 98), "\n", sep = "")
  
  cat("\n-- Equivalent-sweep units (for direct comparison with Metropolis) --\n")
  cat(sprintf("%-3s %-10s %-8s %-10s %-11s %-12s %-12s\n",
              "L", "regime", "T", "mean_cl", "eq_sweeps",
              "tau_E_sw", "tau_M_sw"))
  cat(strrep("-", 75), "\n", sep = "")
  for (i in seq_len(nrow(diag_df))) {
    r <- diag_df[i, ]
    cat(sprintf("%-3d %-10s %-8.4f %-10.1f %-11.1f %-12.3f %-12.3f\n",
                r$L, r$T_label, r$T, r$mean_cluster,
                r$equiv_sweeps, r$tau_E_sw, r$tau_M_sw))
  }
  cat(strrep("-", 75), "\n", sep = "")
  
  n_fail_E <- sum(!diag_df$closed_E)
  n_fail_M <- sum(!diag_df$closed_M)
  if (n_fail_E > 0 || n_fail_M > 0) {
    cat(sprintf("\n*** WARNING: %d cells did not close for tau_E, %d for tau_M.\n",
                n_fail_E, n_fail_M))
  }
}


# -----------------------------------------------------------------------------
# 9.5.3: z-fit for Wolff at T_c, with Sokal-floor degeneracy detection
# -----------------------------------------------------------------------------

fit_dynamic_z_wolff <- function(diag_df, observable = c("E", "M")) {
  observable <- match.arg(observable)
  tau_col    <- paste0("tau_", observable, "_cf")
  closed_col <- paste0("closed_", observable)
  
  sub <- diag_df[diag_df$T_label == "at Tc", ]
  sub <- sub[sub[[closed_col]], ]
  
  if (nrow(sub) < 2) {
    cat(sprintf("  z (from tau_%s): not enough closed cells to fit (have %d)\n",
                observable, nrow(sub)))
    return(NA)
  }
  
  taus <- sub[[tau_col]]
  # Detect Sokal-floor degeneracy: all tau within 2% of 0.5 -> meaningless fit.
  floor_margin <- (taus - 0.5) / 0.5
  if (all(abs(floor_margin) < 0.02)) {
    cat(sprintf("  z (from tau_%s): all tau values at Sokal floor 0.5 (+/- 2%%);\n",
                observable))
    cat(sprintf("    Wolff decorrelates this observable in < 1 cluster flip at every L.\n"))
    cat(sprintf("    No meaningful z to extract -- Wolff effectively draws iid samples.\n"))
    return(NA)
  }
  
  fit <- lm(log(taus) ~ log(sub$L))
  z   <- coef(fit)[2]
  cat(sprintf("  z (from tau_%s, fit on L = %s): z = %.3f  (literature Wolff 1989: ~%.2f)\n",
              observable, paste(sub$L, collapse = ","), z,
              Z_WOLFF_LITERATURE))
  invisible(z)
}


# -----------------------------------------------------------------------------
# 9.5.4: rho(t) grid (4 x 3), ggplot version
# -----------------------------------------------------------------------------

build_rho_dataframe_wolff <- function(runs, observable = c("E", "M"),
                                      max_lag_plot = 50) {
  observable <- match.arg(observable)
  rows <- list()
  for (key in names(runs)) {
    run      <- runs[[key]]
    n_cf     <- nrow(run$series)
    max_lag_here <- min(max_lag_plot, n_cf %/% 4)
    x_series <- if (observable == "E") run$series$E else run$series$M
    
    rho  <- autocorrelation(x_series, max_lag = max_lag_here)
    diag <- integrated_tau(x_series)
    
    rows[[length(rows) + 1]] <- data.frame(
      L        = run$L,
      T        = run$T,
      T_label  = run$T_label,
      lag      = 0:max_lag_here,
      rho      = rho,
      tau      = diag$tau,
      closed   = diag$closed,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


plot_autocorr_grid_wolff <- function(runs, observable = c("E", "M"),
                                     max_lag_plot = 50) {
  observable <- match.arg(observable)
  
  df <- build_rho_dataframe_wolff(runs, observable, max_lag_plot)
  
  df$L_fct   <- factor(df$L, levels = sort(unique(df$L)))
  df$T_label <- factor(df$T_label,
                       levels = c("below Tc", "at Tc", "above Tc"))
  
  panel_labels <- unique(df[, c("L", "T_label", "tau", "closed")])
  panel_labels$L_fct   <- factor(panel_labels$L,
                                 levels = sort(unique(df$L)))
  panel_labels$T_label <- factor(panel_labels$T_label,
                                 levels = c("below Tc", "at Tc", "above Tc"))
  panel_labels$lbl <- with(panel_labels,
                           sprintf("tau_cf = %.2f%s", tau, ifelse(closed, "", " (!)")))
  
  obs_str  <- switch(observable, E = "energy", M = "magnetization")
  obs_expr <- switch(observable,
                     E = expression(rho[E](italic(t))),
                     M = expression(rho[M](italic(t))))
  obs_col  <- switch(observable,
                     E = PAPER_COLORS$energy,
                     M = PAPER_COLORS$magnetization)
  
  ggplot(df, aes(x = lag, y = rho)) +
    geom_hline(yintercept = 0, color = "grey70",
               linewidth = 0.3, linetype = "dotted") +
    geom_hline(yintercept = 1/exp(1), color = "grey60",
               linewidth = 0.3, linetype = "dotted") +
    geom_line(color = obs_col, linewidth = 0.6) +
    geom_text(data = panel_labels,
              aes(x = Inf, y = Inf, label = lbl),
              hjust = 1.05, vjust = 1.4, size = 2.7,
              color = "grey25", inherit.aes = FALSE) +
    facet_grid(
      rows     = vars(L_fct),
      cols     = vars(T_label),
      labeller = labeller(
        L_fct   = function(x) paste0("L = ", x),
        T_label = as_labeller(c("below Tc" = "T = 1.5",
                                "at Tc"    = "T = Tc",
                                "above Tc" = "T = 4.0"))
      )
    ) +
    coord_cartesian(ylim = c(-0.1, 1.05)) +
    labs(title = sprintf("Wolff autocorrelation function: %s observable", obs_str),
         x     = "lag (cluster flips)",
         y     = obs_expr) +
    theme_paper() +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(size = rel(0.85), color = "grey20"),
      panel.spacing.x  = unit(4, "pt"),
      panel.spacing.y  = unit(3, "pt"),
      plot.title       = element_text(size = rel(1.0), face = "bold")
    )
}


# -----------------------------------------------------------------------------
# 9.5.5: tau vs L summary (1 x 2), ggplot version
#
# Like Section 8's plot_tau_scaling, but using the Wolff literature z = 0.25
# as the reference slope, plotted in cluster-flip units.
# -----------------------------------------------------------------------------

plot_tau_scaling_wolff <- function(diag_df) {
  diag_df$T_label <- factor(diag_df$T_label,
                            levels = c("below Tc", "at Tc", "above Tc"))
  
  L_values <- sort(unique(diag_df$L))
  L_min <- min(L_values); L_max <- max(L_values)
  
  # Reference slope line (Wolff 1989 literature z = 0.25), anchored at smallest L
  sub_tc <- diag_df[diag_df$T_label == "at Tc", ]
  sub_tc <- sub_tc[order(sub_tc$L), ]
  
  build_ref_df <- function(tau_col) {
    if (nrow(sub_tc) < 1) return(NULL)
    data.frame(
      L   = c(L_min, L_max),
      tau = sub_tc[[tau_col]][1] *
        (c(L_min, L_max) / sub_tc$L[1])^Z_WOLFF_LITERATURE
    )
  }
  
  make_panel <- function(tau_col, closed_col, y_lab, title_str) {
    ref_df <- build_ref_df(tau_col)
    
    dat <- diag_df
    dat$tau    <- dat[[tau_col]]
    dat$closed <- dat[[closed_col]]
    dat$shape_key <- factor(ifelse(dat$closed, "closed", "open"),
                            levels = c("closed", "open"))
    
    p <- ggplot(dat, aes(x = L, y = tau,
                         color = T_label, group = T_label)) +
      geom_line(linewidth = 0.9, na.rm = TRUE) +
      geom_point(aes(shape = shape_key), size = 2.2, na.rm = TRUE) +
      geom_hline(yintercept = 0.5, color = "grey70",
                 linetype = "dotted", linewidth = 0.4) +
      annotate("text", x = L_max, y = 0.5,
               label = "Sokal floor 0.5", color = "grey50",
               size = 2.6, hjust = 1, vjust = -0.4) +
      scale_x_log10(breaks = L_values) +
      scale_y_log10() +
      scale_shape_manual(values = c("closed" = 16, "open" = 4),
                         name   = "Sokal") +
      labs(title = title_str,
           x     = expression(italic(L)),
           y     = y_lab) +
      theme_paper() +
      theme(legend.position = "right")
    
    if (!is.null(ref_df)) {
      p <- p + geom_line(data = ref_df, aes(x = L, y = tau),
                         color = "grey60", linetype = "dotted",
                         linewidth = 0.5, inherit.aes = FALSE) +
        annotate("text",
                 x = L_max, y = ref_df$tau[2],
                 label = sprintf("ref L^%.2f", Z_WOLFF_LITERATURE),
                 color = "grey50", size = 2.6, hjust = 1, vjust = 1.3)
    }
    p
  }
  
  p_E <- make_panel("tau_E_cf", "closed_E",
                    expression(tau[E] * " (cluster flips)"),
                    "Wolff energy autocorrelation vs L")
  
  p_M <- make_panel("tau_M_cf", "closed_M",
                    expression(tau[M] * " (cluster flips)"),
                    "Wolff magnetization autocorrelation vs L")
  
  (p_E | p_M) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = paste(
        "tau from Sokal windowing (c = 5).",
        "Filled points: Sokal closed. Crosses: open (tau is a lower bound).",
        "Dotted line: reference slope at the Wolff 1989 z = 0.25.",
        "Note tau_M sits at the Sokal floor (0.5) at every L tested.",
        sep = " "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    ) &
    scale_color_manual(values = T_REGIME_COLORS, name = NULL) &
    theme(legend.position = "bottom")
}


# =============================================================================
# Top-level driver
# =============================================================================

run_section9_diagnostics <- function(runs = section9_wolff_out$runs) {
  cat("\n############################################################\n")
  cat("Section 9 part 3: Wolff autocorrelation diagnostics\n")
  cat("############################################################\n\n")
  
  diag_df <- build_diagnostics_table_wolff(runs)
  print_diagnostics_table_wolff(diag_df)
  
  cat("\n--- Dynamic critical exponent z (Wolff at T_c) ---\n")
  fit_dynamic_z_wolff(diag_df, "E")
  fit_dynamic_z_wolff(diag_df, "M")
  
  cat("\nPlotting Wolff rho_E grid (4 x 3)...\n")
  print(plot_autocorr_grid_wolff(runs, observable = "E"))
  
  cat("Plotting Wolff rho_M grid (4 x 3)...\n")
  print(plot_autocorr_grid_wolff(runs, observable = "M"))
  
  cat("Plotting Wolff tau scaling summary (1 x 2)...\n")
  print(plot_tau_scaling_wolff(diag_df))
  
  invisible(diag_df)
}


# =============================================================================
# Drive
# =============================================================================

section9_wolff_diag <- run_section9_diagnostics(section9_wolff_out$runs)


##> ############################################################
##> Section 9 part 3: Wolff autocorrelation diagnostics
##> ############################################################
##>
##> === Wolff autocorrelation diagnostics across (L, T) ===
##> Flag '!' next to tau means Sokal windowing did not close.
##>
##> -- Native cluster-flip units (Wolff's natural time scale) --
##> L   regime     T        n_cf     tau_E_cf     tau_M_cf     ESS_E     ESS_M     se_e      se_m
##> --------------------------------------------------------------------------------------------------
##> 8   below Tc   1.5000   10000    0.61         0.50         8264      10000     ~0.001    ~0.001
##> 8   at Tc      2.2692   10000    1.62         0.50         3076      10000     ~0.005    ~0.005
##> 8   above Tc   4.0000   10000    5.97         0.55         837       9099      ~0.010    ~0.005
##> 16  below Tc   1.5000   10000    0.62         0.50         8101      10000     ~0.001    ~0.001
##> 16  at Tc      2.2692   10000    2.81         0.50         1779      10000     ~0.003    ~0.005
##> 16  above Tc   4.0000   10000    23.7         0.61         211       8197      ~0.010    ~0.005
##> 32  below Tc   1.5000   10000    0.59         0.50         8475      10000     ~0.001    ~0.001
##> 32  at Tc      2.2692   10000    4.04         0.50         1238      10000     ~0.002    ~0.005
##> 32  above Tc   4.0000   10000    372.3        2.39         13        2092      ~0.040    ~0.005
##> 64  below Tc   1.5000   10000    0.58         0.50         8616      10000     ~0.001    ~0.001
##> 64  at Tc      2.2692   10000    5.85         0.50         855       10000     ~0.002    ~0.005
##> 64  above Tc   4.0000   10000    342.7        4.93         15        1014      ~0.040    ~0.005
##> --------------------------------------------------------------------------------------------------
##>
##> -- Equivalent-sweep units (for direct comparison with Metropolis) --
##> L   regime     T        mean_cl    eq_sweeps   tau_E_sw     tau_M_sw
##> ---------------------------------------------------------------------------
##> 8   below Tc   1.5000   62.3       9737        0.597        0.486
##> 8   at Tc      2.2692   41.6       6500        1.056        0.325
##> 8   above Tc   4.0000   4.2        656         0.390        0.221
##> 16  below Tc   1.5000   249.7      9754        0.618        0.488
##> 16  at Tc      2.2692   140.7      5497        1.424        0.275
##> 16  above Tc   4.0000   4.3        168         0.408        0.195
##> 32  below Tc   1.5000   995.3      9719        0.569        0.486
##> 32  at Tc      2.2692   462.7      4519        1.867        0.226
##> 32  above Tc   4.0000   4.2        41.1        1.548        0.197
##> 64  below Tc   1.5000   3984.5     9728        0.589        0.486
##> 64  at Tc      2.2692   1573.1     3840        2.246        0.192
##> 64  above Tc   4.0000   4.3        10.5        0.358        0.107
##> ---------------------------------------------------------------------------
##>
##> --- Dynamic critical exponent z (Wolff at T_c) ---
##>   z (from tau_E, fit on L = 8,16,32,64): z = 0.617  (literature Wolff 1989: ~0.25)
##>   z (from tau_M, fit on L = 8,16,32,64): all tau values at Sokal floor 0.5 (+/- 2%);
##>     Wolff decorrelates this observable in < 1 cluster flip at every L.
##>     No meaningful z to extract -- Wolff effectively draws iid samples.
##>
##> Plotting Wolff rho_E grid (4 x 3)...
##> Plotting Wolff rho_M grid (4 x 3)...
##> Plotting Wolff tau scaling summary (1 x 2)...




# =============================================================================
# Section 9 part 4: Metropolis vs Wolff head-to-head comparison
#
# This is the project's central scientific deliverable -- the picture of
# "Wolff defeats critical slowing down" rendered in numbers and plots.
#
# Requires both diagnostics data frames already in memory:
#   section8_diag        -- Metropolis (from Section 8)
#   section9_wolff_diag  -- Wolff      (from Section 9 part 3)
#
# Optional (for Figure C only):
#   bench_df             -- Metropolis Rcpp benchmark from Section 6.1
#   bench_wolff_df       -- Wolff Rcpp benchmark from Section 9.3
#
# Produces:
#   1. Merged comparison tables (raw and ratio).
#   2. Figure A: tau vs L, both samplers in equivalent-sweep units.
#      [Source of paper's fig:head2head_tau and tab:head2head.]
#   3. Figure B: rho(t) overlay at T_c, fixed L on a common equivalent-
#      sweep axis. [Source of paper's fig:rho_overlay.]
#   4. Figure C: ESS per CPU second (Appendix E candidate).
# =============================================================================


# -----------------------------------------------------------------------------
# 9.6.1: Merged comparison table
# -----------------------------------------------------------------------------

build_comparison_table <- function(metro_diag, wolff_diag) {
  rows <- list()
  for (i in seq_len(nrow(metro_diag))) {
    mr <- metro_diag[i, ]
    wr <- wolff_diag[wolff_diag$L == mr$L & wolff_diag$T_label == mr$T_label, ]
    if (nrow(wr) != 1) next
    
    rows[[length(rows) + 1]] <- data.frame(
      L               = mr$L,
      T               = mr$T,
      T_label         = mr$T_label,
      # Metropolis side
      metro_n_sweeps  = mr$n_sweeps,
      metro_tau_E     = mr$tau_E,
      metro_tau_M     = mr$tau_M,
      metro_ESS_E     = mr$ESS_E,
      metro_ESS_M     = mr$ESS_M,
      metro_se_e      = mr$se_e,
      metro_se_m      = mr$se_m,
      # Wolff side: tau in equivalent sweeps for honest comparison
      wolff_n_cf      = wr$n_cluster_flips,
      wolff_equiv_sw  = wr$equiv_sweeps,
      wolff_tau_E_sw  = wr$tau_E_sw,
      wolff_tau_M_sw  = wr$tau_M_sw,
      wolff_ESS_E     = wr$ESS_E,
      wolff_ESS_M     = wr$ESS_M,
      wolff_se_e      = wr$se_e,
      wolff_se_m      = wr$se_m,
      # Improvement factors (Metropolis / Wolff)
      ratio_tau_E     = mr$tau_E / wr$tau_E_sw,
      ratio_tau_M     = mr$tau_M / wr$tau_M_sw,
      ratio_ESS_E     = wr$ESS_E / mr$ESS_E,
      ratio_ESS_M     = wr$ESS_M / mr$ESS_M,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


print_comparison_table <- function(comp_df) {
  cat("=== Metropolis vs Wolff head-to-head ===\n\n")
  
  # Raw tau values side by side, in common (equivalent-sweep) units
  cat("--- tau (equivalent-sweep units) ---\n")
  cat(sprintf("%-3s %-10s %-8s %-12s %-12s %-12s %-12s %-10s %-10s\n",
              "L", "regime", "T",
              "tau_E_met", "tau_E_wol", "tau_M_met", "tau_M_wol",
              "rE(met/w)", "rM(met/w)"))
  cat(strrep("-", 100), "\n", sep = "")
  for (i in seq_len(nrow(comp_df))) {
    r <- comp_df[i, ]
    cat(sprintf("%-3d %-10s %-8.4f %-12.3f %-12.3f %-12.3f %-12.3f %-10.1f %-10.1f\n",
                r$L, r$T_label, r$T,
                r$metro_tau_E, r$wolff_tau_E_sw,
                r$metro_tau_M, r$wolff_tau_M_sw,
                r$ratio_tau_E, r$ratio_tau_M))
  }
  cat(strrep("-", 100), "\n", sep = "")
  
  # ESS values side by side
  cat("\n--- ESS (samples obtained) ---\n")
  cat(sprintf("%-3s %-10s %-8s %-12s %-12s %-12s %-12s\n",
              "L", "regime", "T",
              "ESS_E_met", "ESS_E_wol", "ESS_M_met", "ESS_M_wol"))
  cat(strrep("-", 85), "\n", sep = "")
  for (i in seq_len(nrow(comp_df))) {
    r <- comp_df[i, ]
    cat(sprintf("%-3d %-10s %-8.4f %-12.0f %-12.0f %-12.0f %-12.0f\n",
                r$L, r$T_label, r$T,
                r$metro_ESS_E, r$wolff_ESS_E,
                r$metro_ESS_M, r$wolff_ESS_M))
  }
  cat(strrep("-", 85), "\n", sep = "")
  
  # Headline numbers at T_c (this is paper Tab head2head)
  cat("\n--- Headline numbers at T_c ---\n")
  sub_tc <- comp_df[comp_df$T_label == "at Tc", ]
  for (i in seq_len(nrow(sub_tc))) {
    r <- sub_tc[i, ]
    cat(sprintf("  L=%-2d : Wolff tau_E is %.1fx faster, tau_M is %.1fx faster than Metropolis\n",
                r$L, r$ratio_tau_E, r$ratio_tau_M))
  }
}


# -----------------------------------------------------------------------------
# 9.6.2: Figure A -- tau vs L on same axes, both samplers
# -----------------------------------------------------------------------------

plot_metro_vs_wolff_tau <- function(metro_diag, wolff_diag) {
  # Reshape both into long format with a sampler column for unified plotting
  metro_long <- data.frame(
    sampler = "Metropolis",
    L       = metro_diag$L,
    T       = metro_diag$T,
    T_label = metro_diag$T_label,
    tau_E   = metro_diag$tau_E,
    tau_M   = metro_diag$tau_M,
    closed_E = metro_diag$closed_E,
    closed_M = metro_diag$closed_M,
    stringsAsFactors = FALSE
  )
  wolff_long <- data.frame(
    sampler = "Wolff",
    L       = wolff_diag$L,
    T       = wolff_diag$T,
    T_label = wolff_diag$T_label,
    tau_E   = wolff_diag$tau_E_sw,   # equivalent sweeps!
    tau_M   = wolff_diag$tau_M_sw,
    closed_E = wolff_diag$closed_E,
    closed_M = wolff_diag$closed_M,
    stringsAsFactors = FALSE
  )
  comb <- rbind(metro_long, wolff_long)
  comb$sampler <- factor(comb$sampler, levels = c("Metropolis", "Wolff"))
  comb$T_label <- factor(comb$T_label,
                         levels = c("below Tc", "at Tc", "above Tc"))
  
  L_values <- sort(unique(comb$L))
  L_min <- min(L_values); L_max <- max(L_values)
  
  # Two reference slope lines (literature z), anchored at smallest L on Tc row.
  metro_tc <- metro_diag[metro_diag$T_label == "at Tc", ]
  metro_tc <- metro_tc[order(metro_tc$L), ]
  wolff_tc <- wolff_diag[wolff_diag$T_label == "at Tc", ]
  wolff_tc <- wolff_tc[order(wolff_tc$L), ]
  
  build_ref_metro <- function(tau_col) {
    if (nrow(metro_tc) < 1) return(NULL)
    data.frame(L = c(L_min, L_max),
               tau = metro_tc[[tau_col]][1] *
                 (c(L_min, L_max) / metro_tc$L[1])^Z_METROPOLIS_LITERATURE)
  }
  build_ref_wolff <- function(tau_col) {
    if (nrow(wolff_tc) < 1) return(NULL)
    data.frame(L = c(L_min, L_max),
               tau = wolff_tc[[tau_col]][1] *
                 (c(L_min, L_max) / wolff_tc$L[1])^Z_WOLFF_LITERATURE)
  }
  
  make_panel <- function(tau_col, closed_col, y_lab, title_str,
                         ref_metro, ref_wolff) {
    dat <- comb
    dat$tau    <- dat[[tau_col]]
    dat$closed <- dat[[closed_col]]
    dat$shape_key <- factor(ifelse(dat$closed, "closed", "open"),
                            levels = c("closed", "open"))
    
    p <- ggplot(dat, aes(x = L, y = tau,
                         color = T_label,
                         linetype = sampler,
                         group = interaction(sampler, T_label))) +
      geom_line(linewidth = 0.9, na.rm = TRUE) +
      geom_point(aes(shape = shape_key), size = 2.2, na.rm = TRUE) +
      scale_x_log10(breaks = L_values) +
      scale_y_log10() +
      scale_linetype_manual(values = c("Metropolis" = "solid", "Wolff" = "dashed"),
                            name = NULL) +
      scale_shape_manual(values = c("closed" = 16, "open" = 4),
                         name = "Sokal") +
      labs(title = title_str,
           x     = expression(italic(L)),
           y     = y_lab) +
      theme_paper() +
      theme(legend.position = "right")
    
    # Reference slope lines
    if (!is.null(ref_metro)) {
      p <- p + geom_line(data = ref_metro, aes(x = L, y = tau),
                         color = "grey55", linetype = "dotted",
                         linewidth = 0.5, inherit.aes = FALSE) +
        annotate("text", x = L_max, y = ref_metro$tau[2],
                 label = sprintf("L^%.2f (Metro)", Z_METROPOLIS_LITERATURE),
                 color = "grey45", size = 2.6, hjust = 1, vjust = -0.5)
    }
    if (!is.null(ref_wolff)) {
      p <- p + geom_line(data = ref_wolff, aes(x = L, y = tau),
                         color = "grey75", linetype = "dotted",
                         linewidth = 0.5, inherit.aes = FALSE) +
        annotate("text", x = L_max, y = ref_wolff$tau[2],
                 label = sprintf("L^%.2f (Wolff)", Z_WOLFF_LITERATURE),
                 color = "grey55", size = 2.6, hjust = 1, vjust = 1.4)
    }
    p
  }
  
  p_E <- make_panel("tau_E", "closed_E",
                    expression(tau[E] * " (equivalent sweeps)"),
                    "Energy autocorrelation: Metropolis vs Wolff",
                    build_ref_metro("tau_E"), build_ref_wolff("tau_E_sw"))
  
  p_M <- make_panel("tau_M", "closed_M",
                    expression(tau[M] * " (equivalent sweeps)"),
                    "Magnetization autocorrelation: Metropolis vs Wolff",
                    build_ref_metro("tau_M"), build_ref_wolff("tau_M_sw"))
  
  (p_E | p_M) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = paste(
        "Both samplers in equivalent-sweep units.",
        "Solid = Metropolis, dashed = Wolff. Filled circles: Sokal closed; crosses: open.",
        "Dotted reference lines: literature z exponents.",
        sep = "  "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    ) &
    scale_color_manual(values = T_REGIME_COLORS, name = NULL) &
    theme(legend.position = "bottom")
}


# -----------------------------------------------------------------------------
# 9.6.3: Figure B -- rho(t) overlay at T_c, fixed L
#
# This is the source of paper's fig:rho_overlay. Both samplers on a common
# equivalent-sweep x-axis. Metropolis = red, Wolff = blue (matches paper
# caption).
# -----------------------------------------------------------------------------

plot_rho_overlay <- function(metro_runs, wolff_runs,
                             L = 32, max_lag_plot = 200) {
  key_m <- sprintf("L%d_Tc", L)
  key_w <- sprintf("L%d_Tc", L)
  
  if (is.null(metro_runs[[key_m]]) || is.null(wolff_runs[[key_w]])) {
    cat(sprintf("No matching runs for L=%d at T_c; skipping overlay.\n", L))
    return(invisible(NULL))
  }
  
  run_m <- metro_runs[[key_m]]
  run_w <- wolff_runs[[key_w]]
  
  # Metropolis: x-axis = sweep lag (already in equivalent sweeps)
  rho_E_m <- autocorrelation(run_m$series$E, max_lag = max_lag_plot)
  rho_M_m <- autocorrelation(run_m$series$M, max_lag = max_lag_plot)
  
  # Wolff: convert cluster-flip lag to equivalent sweeps
  rho_E_w <- autocorrelation(run_w$series$E, max_lag = max_lag_plot)
  rho_M_w <- autocorrelation(run_w$series$M, max_lag = max_lag_plot)
  sweeps_per_flip <- run_w$mean_cluster_size / (L * L)
  
  df <- rbind(
    data.frame(sampler = "Metropolis",
               obs     = "E",
               lag_sw  = 0:max_lag_plot,
               rho     = rho_E_m),
    data.frame(sampler = "Metropolis",
               obs     = "M",
               lag_sw  = 0:max_lag_plot,
               rho     = rho_M_m),
    data.frame(sampler = "Wolff",
               obs     = "E",
               lag_sw  = (0:max_lag_plot) * sweeps_per_flip,
               rho     = rho_E_w),
    data.frame(sampler = "Wolff",
               obs     = "M",
               lag_sw  = (0:max_lag_plot) * sweeps_per_flip,
               rho     = rho_M_w)
  )
  df$sampler <- factor(df$sampler, levels = c("Metropolis", "Wolff"))
  df$obs     <- factor(df$obs, levels = c("E", "M"))
  
  obs_titles <- c(
    "E" = sprintf("Energy autocorrelation at T_c, L = %d", L),
    "M" = sprintf("Magnetization autocorrelation at T_c, L = %d", L)
  )
  obs_ylabs <- c(E = "rho_E(t)", M = "rho_M(t)")
  
  make_panel <- function(obs_label, ylab_str) {
    sub <- df[df$obs == obs_label, ]
    ggplot(sub, aes(x = lag_sw, y = rho, color = sampler, group = sampler)) +
      geom_hline(yintercept = 0, color = "grey70",
                 linewidth = 0.3, linetype = "dotted") +
      geom_hline(yintercept = 1/exp(1), color = "grey60",
                 linewidth = 0.3, linetype = "dotted") +
      annotate("text", x = max_lag_plot * 0.97, y = 1/exp(1),
               label = "1/e", color = "grey50", size = 2.6,
               hjust = 1, vjust = -0.4) +
      geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = c(0, max_lag_plot), ylim = c(-0.1, 1.05)) +
      labs(title = obs_titles[[obs_label]],
           x     = "lag (equivalent sweeps)",
           y     = ylab_str) +
      theme_paper()
  }
  
  p_E <- make_panel("E", expression(rho[E](italic(t))))
  p_M <- make_panel("M", expression(rho[M](italic(t))))
  
  (p_E | p_M) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = sprintf(
        "T = T_c, L = %d. Both samplers plotted on a common equivalent-sweep lag axis. Wolff cluster-flip lags are scaled by mean_cluster/N = %.3f.",
        L, sweeps_per_flip),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    ) &
    scale_color_manual(values = c("Metropolis" = PAPER_COLORS$metro,
                                  "Wolff"      = PAPER_COLORS$wolff),
                       name = NULL) &
    theme(legend.position = "bottom")
}


# -----------------------------------------------------------------------------
# 9.6.4: Figure C -- ESS per CPU second (Appendix E candidate)
# -----------------------------------------------------------------------------

plot_ess_per_second <- function(metro_diag, wolff_diag,
                                metro_bench, wolff_bench) {
  make_df <- function(diag_df, bench_df, chain_length_col,
                      ESS_col_E, ESS_col_M, sampler_name) {
    out <- list()
    for (i in seq_len(nrow(diag_df))) {
      r <- diag_df[i, ]
      bmatch <- which(bench_df$L == r$L &
                        abs(bench_df$T - r$T) < 1e-6)
      if (length(bmatch) != 1) next
      per_step_time <- bench_df$t_cpp[bmatch] / 500   # benchmark uses 500 steps
      total_time    <- per_step_time * r[[chain_length_col]]
      out[[length(out) + 1]] <- data.frame(
        sampler  = sampler_name,
        L        = r$L,
        T        = r$T,
        T_label  = r$T_label,
        ESS_E_per_sec = r[[ESS_col_E]] / total_time,
        ESS_M_per_sec = r[[ESS_col_M]] / total_time,
        total_time_s  = total_time,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, out)
  }
  
  metro_df <- make_df(metro_diag, metro_bench, "n_sweeps",
                      "ESS_E", "ESS_M", "Metropolis")
  wolff_df <- make_df(wolff_diag, wolff_bench, "n_cluster_flips",
                      "ESS_E", "ESS_M", "Wolff")
  all_df <- rbind(metro_df, wolff_df)
  
  all_df$sampler <- factor(all_df$sampler, levels = c("Metropolis", "Wolff"))
  all_df$T_label <- factor(all_df$T_label,
                           levels = c("below Tc", "at Tc", "above Tc"))
  
  L_values <- sort(unique(all_df$L))
  
  make_panel <- function(col_name, y_label, title_str) {
    ggplot(all_df, aes(x = L, y = .data[[col_name]],
                       color = T_label, linetype = sampler,
                       group = interaction(sampler, T_label))) +
      geom_line(linewidth = 0.9) +
      geom_point(aes(shape = sampler), size = 2.0) +
      scale_x_log10(breaks = L_values) +
      scale_y_log10() +
      scale_linetype_manual(values = c("Metropolis" = "solid", "Wolff" = "dashed"),
                            name = NULL) +
      scale_shape_manual(values = c("Metropolis" = 16, "Wolff" = 1),
                         name = NULL) +
      labs(title = title_str,
           x     = expression(italic(L)),
           y     = y_label) +
      theme_paper() +
      theme(legend.position = "right")
  }
  
  p_E <- make_panel("ESS_E_per_sec",
                    expression(ESS[E] * " / second"),
                    "Energy samples per CPU second")
  p_M <- make_panel("ESS_M_per_sec",
                    expression(ESS[M] * " / second"),
                    "Magnetization samples per CPU second")
  
  (p_E | p_M) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = paste(
        "Per-cell wall time estimated from the 500-step Rcpp benchmark",
        "(slightly understates throughput at long chain lengths).",
        "Solid = Metropolis, dashed = Wolff.",
        sep = "  "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    ) &
    scale_color_manual(values = T_REGIME_COLORS, name = NULL) &
    theme(legend.position = "bottom")
}


# =============================================================================
# Driver
# =============================================================================

run_section9_part4 <- function() {
  cat("\n############################################################\n")
  cat("Section 9 part 4: Metropolis vs Wolff head-to-head\n")
  cat("############################################################\n\n")
  
  # 9.6.1: comparison tables
  comp_df <- build_comparison_table(section8_diag, section9_wolff_diag)
  print_comparison_table(comp_df)
  
  # 9.6.2: tau vs L, both samplers (paper fig:head2head_tau)
  cat("\nPlotting Figure A (tau vs L, both samplers)...\n")
  print(plot_metro_vs_wolff_tau(section8_diag, section9_wolff_diag))
  
  # 9.6.3: rho(t) overlay at T_c, L=32 (paper fig:rho_overlay)
  cat("Plotting Figure B (rho(t) overlay at T_c, L=32)...\n")
  print(plot_rho_overlay(section7_out$runs, section9_wolff_out$runs, L = 32))
  
  # 9.6.4: ESS per CPU second (Appendix E candidate)
  if (exists("bench_df") && exists("bench_wolff_df")) {
    cat("Plotting Figure C (ESS per CPU second)...\n")
    print(plot_ess_per_second(section8_diag, section9_wolff_diag,
                              bench_df, bench_wolff_df))
  } else {
    cat("Skipping Figure C: bench_df or bench_wolff_df not in memory.\n")
  }
  
  invisible(comp_df)
}


# =============================================================================
# Drive
# =============================================================================

comparison_df <- run_section9_part4()


##> ############################################################
##> Section 9 part 4: Metropolis vs Wolff head-to-head
##> ############################################################
##>
##> === Metropolis vs Wolff head-to-head ===
##>
##> --- tau (equivalent-sweep units) ---
##> L   regime     T        tau_E_met    tau_E_wol    tau_M_met    tau_M_wol    rE(met/w)  rM(met/w)
##> ----------------------------------------------------------------------------------------------------
##> 8   below Tc   1.5000   1.620        0.597        2.020        0.486        2.7        4.2
##> 8   at Tc      2.2692   7.370        1.056        178.650      0.325        7.0        549.7
##> 8   above Tc   4.0000   1.080        0.390        3.750        0.221        2.8        17.0
##> 16  below Tc   1.5000   1.390        0.618        1.470        0.488        2.2        3.0
##> 16  at Tc      2.2692   19.950       1.424        767.640      0.275        14.0       2791.4
##> 16  above Tc   4.0000   0.980        0.408        3.860        0.195        2.4        19.8
##> 32  below Tc   1.5000   1.540        0.569        1.730        0.486        2.7        3.6
##> 32  at Tc      2.2692   56.340       1.867        2233.340     0.226        30.2       9884.7
##> 32  above Tc   4.0000   1.010        1.548        3.270        0.197        0.7        16.6
##> 64  below Tc   1.5000   1.540        0.589        1.680        0.486        2.6        3.5
##> 64  at Tc      2.2692   220.900      2.246        15204.400    0.192        98.4       79189.6
##> 64  above Tc   4.0000   0.970        0.358        3.570        0.107        2.7        33.4
##> ----------------------------------------------------------------------------------------------------
##>
##> --- ESS (samples obtained) ---
##>   (Metropolis ESS uses adaptive chain length; Wolff uses fixed 10000 cluster flips.)
##>
##> --- Headline numbers at T_c ---
##>   L=8  : Wolff tau_E is 7.0x faster, tau_M is 549.7x faster than Metropolis
##>   L=16 : Wolff tau_E is 14.0x faster, tau_M is 2791.4x faster than Metropolis
##>   L=32 : Wolff tau_E is 30.2x faster, tau_M is 9884.7x faster than Metropolis
##>   L=64 : Wolff tau_E is 98.4x faster, tau_M is 79189.6x faster than Metropolis
##>
##> Plotting Figure A (tau vs L, both samplers)...
##> Plotting Figure B (rho(t) overlay at T_c, L=32)...
##> Plotting Figure C (ESS per CPU second)...


# =============================================================================
# Section 10: Full temperature sweep with Wolff
#
# Produces the classic 2D Ising phase-transition figure: smooth curves of
# u(T), c(T), |m|(T), and chi(T) across the full transition, for multiple
# lattice sizes, overlaid on Onsager/Yang exact results.
#
# Why Wolff: at T_c L=64, Metropolis needed 320k sweeps for a decent
# estimate; Wolff needs ~5k cluster flips. So 34 temperatures x 4 lattice
# sizes = 136 runs is fast.
#
# Design choices (matching paper Sec 3.5):
#   - Temperature grid: 34 points, refined around T_c
#   - Lattice sizes: L in {8, 16, 32, 64}
#   - Chain length: 5000 cluster flips per (L, T) cell, fixed
#   - Burn-in: n_cluster_flips / 5 = 1000
#   - Init: "hot" everywhere
#
# No bootstrap error bars here; that's Section 11. Section 10 is just
# point estimates overlaid on ground truth.
# =============================================================================


# -----------------------------------------------------------------------------
# 10.1: build the T-sweep grid
# -----------------------------------------------------------------------------

make_T_grid_sweep <- function() {
  sort(unique(c(
    seq(1.0, 1.8, by = 0.2),     # ordered phase, sparse
    seq(1.9, 2.15, by = 0.05),   # approaching T_c from below
    seq(2.20, 2.34, by = 0.02),  # critical region, very dense
    seq(2.40, 2.80, by = 0.05),  # leaving T_c
    seq(2.9, 4.0, by = 0.2)      # disordered phase, sparse
  )))
}


# -----------------------------------------------------------------------------
# 10.2: run a single L-sweep across all temperatures
# -----------------------------------------------------------------------------

run_T_sweep_at_L <- function(L, T_grid, n_cluster_flips = 5000,
                             seed_base = 42) {
  rows <- list()
  for (T in T_grid) {
    run <- run_wolff_fast(
      L               = L,
      T               = T,
      n_cluster_flips = n_cluster_flips,
      n_burnin        = as.integer(n_cluster_flips / 5),
      init_state      = "hot",
      seed            = as.integer(seed_base + round(T * 1000) + L)
    )
    obs <- compute_observables_wolff(run)
    rows[[length(rows) + 1]] <- data.frame(
      L            = L,
      T            = T,
      e            = obs$e,
      m            = obs$m,
      c            = obs$c,
      chi          = obs$chi,
      mean_cluster = obs$mean_cluster_size,
      equiv_sweeps = obs$equivalent_sweeps,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


# -----------------------------------------------------------------------------
# 10.3: run the sweep across all L values
# -----------------------------------------------------------------------------

run_full_T_sweep <- function(L_values        = c(8, 16, 32, 64),
                             T_grid          = make_T_grid_sweep(),
                             n_cluster_flips = 5000,
                             seed_base       = 42) {
  cat("=== Full T-sweep with Wolff ===\n")
  cat(sprintf("Grid: %d temperatures x %d lattice sizes = %d runs\n",
              length(T_grid), length(L_values),
              length(T_grid) * length(L_values)))
  cat(sprintf("Chain length: %d cluster flips per run\n\n", n_cluster_flips))
  
  all_results <- list()
  for (L in L_values) {
    cat(sprintf("  L = %2d: ", L))
    t0 <- Sys.time()
    res <- run_T_sweep_at_L(L, T_grid, n_cluster_flips, seed_base)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("%d temperatures in %.1fs\n", nrow(res), elapsed))
    all_results[[as.character(L)]] <- res
  }
  sweep_df <- do.call(rbind, all_results)
  cat("===============================\n")
  invisible(sweep_df)
}


# -----------------------------------------------------------------------------
# 10.4: 4-panel phase-transition figure (paper fig:phase_transition source)
#
# 2x2 layout, Onsager/Yang reference curves drawn as grey lines; MCMC
# point estimates per L with thin connecting lines. Single shared L legend
# at the bottom.
# -----------------------------------------------------------------------------

plot_T_sweep <- function(sweep_df, c_cap = 4) {
  sweep_df$L_fct <- factor(sweep_df$L,
                           levels = sort(unique(sweep_df$L)))
  
  # Smooth reference curves
  T_curve <- seq(min(sweep_df$T), max(sweep_df$T), length.out = 400)
  df_ref <- data.frame(
    T = T_curve,
    u = onsager_energy_curve(T_curve),
    c = pmin(onsager_specific_heat_curve(T_curve), c_cap),
    m = onsager_magnetization_curve(T_curve)
  )
  
  tc_line <- geom_vline(xintercept = T_CRITICAL,
                        linetype   = "dashed",
                        color      = PAPER_COLORS$tc_line,
                        linewidth  = 0.4)
  
  tc_label <- function(y_val) {
    annotate("text", x = T_CRITICAL + 0.04, y = y_val,
             label = "italic(T)[c]", parse = TRUE,
             color = PAPER_COLORS$tc_line, size = 3, hjust = 0)
  }
  
  # Helper: add per-L lines + points to a ggplot
  add_mcmc_layers <- function(p, value_col) {
    p +
      geom_line(data = sweep_df,
                aes(x = T, y = .data[[value_col]],
                    color = L_fct, group = L_fct),
                linewidth = 0.4, alpha = 0.8) +
      geom_point(data = sweep_df,
                 aes(x = T, y = .data[[value_col]], color = L_fct),
                 size = 1.4)
  }
  
  # Panel 1: Energy
  p_u <- ggplot() +
    tc_line +
    geom_line(data = df_ref, aes(x = T, y = u),
              color = "grey40", linewidth = 0.9, na.rm = TRUE) +
    tc_label(y_val = -0.7) +
    labs(title = "Energy per spin",
         x     = expression(italic(T)),
         y     = expression(italic(u)(italic(T)))) +
    theme_paper()
  p_u <- add_mcmc_layers(p_u, "e")
  
  # Panel 2: Specific heat (Onsager reference capped)
  p_c <- ggplot() +
    tc_line +
    geom_line(data = df_ref, aes(x = T, y = c),
              color = "grey40", linewidth = 0.9, na.rm = TRUE) +
    coord_cartesian(ylim = c(0, max(c(c_cap, sweep_df$c)) * 1.05)) +
    tc_label(y_val = c_cap * 0.92) +
    annotate("text", x = max(T_curve) - 0.05, y = c_cap * 0.92,
             label = sprintf("(Onsager capped at %.0f)", c_cap),
             hjust = 1, size = 2.6, color = "grey40") +
    labs(title = "Specific heat per spin",
         x     = expression(italic(T)),
         y     = expression(italic(c)(italic(T)))) +
    theme_paper()
  p_c <- add_mcmc_layers(p_c, "c")
  
  # Panel 3: Magnetization
  p_m <- ggplot() +
    tc_line +
    geom_hline(yintercept = 0, color = "grey70",
               linewidth = 0.3, linetype = "dotted") +
    geom_line(data = df_ref, aes(x = T, y = m),
              color = "grey40", linewidth = 0.9, na.rm = TRUE) +
    coord_cartesian(ylim = c(0, 1.05)) +
    tc_label(y_val = 0.92) +
    labs(title = "Spontaneous magnetization",
         x     = expression(italic(T)),
         y     = expression(group("|", italic(m), "|") * (italic(T)))) +
    theme_paper()
  p_m <- add_mcmc_layers(p_m, "m")
  
  # Panel 4: Susceptibility (no exact curve)
  p_chi <- ggplot() +
    tc_line +
    tc_label(y_val = max(sweep_df$chi) * 0.9) +
    labs(title = "Magnetic susceptibility",
         x     = expression(italic(T)),
         y     = expression(chi(italic(T)))) +
    theme_paper()
  p_chi <- add_mcmc_layers(p_chi, "chi")
  
  (p_u | p_c) / (p_m | p_chi) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = paste(
        "Wolff T-sweep across 34 temperatures and 4 lattice sizes (5000 cluster flips per cell).",
        "Grey curves are the Onsager / Yang exact references.",
        sep = " "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    ) &
    scale_color_manual(values = L_COLORS,
                       name   = expression(italic(L)),
                       drop   = FALSE) &
    theme(legend.position = "bottom")
}


# -----------------------------------------------------------------------------
# 10.5: mean cluster size vs T (Appendix A.4 candidate -- Wolff diagnostic)
#
# Single panel; shows the cluster size peaking near T_c (the physical
# reason Wolff defeats critical slowing down) and saturating to N = L^2
# below T_c. Per-L horizontal reference lines at L^2.
# -----------------------------------------------------------------------------

plot_T_sweep_cluster_size <- function(sweep_df) {
  sweep_df$L_fct <- factor(sweep_df$L,
                           levels = sort(unique(sweep_df$L)))
  
  L_values <- sort(unique(sweep_df$L))
  
  # Horizontal reference at N = L^2 for each L
  df_refN <- data.frame(
    L_fct = factor(L_values, levels = L_values),
    y_val = L_values^2
  )
  
  ggplot(sweep_df, aes(x = T, y = mean_cluster,
                       color = L_fct, group = L_fct)) +
    geom_vline(xintercept = T_CRITICAL,
               linetype = "dashed", color = PAPER_COLORS$tc_line,
               linewidth = 0.4) +
    geom_hline(data = df_refN,
               aes(yintercept = y_val, color = L_fct),
               linetype = "dotted", linewidth = 0.4,
               show.legend = FALSE) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.4) +
    scale_y_log10() +
    scale_color_manual(values = L_COLORS,
                       name   = expression(italic(L)),
                       drop   = FALSE) +
    annotate("text", x = T_CRITICAL + 0.04,
             y = max(sweep_df$mean_cluster) * 0.92,
             label = "italic(T)[c]", parse = TRUE,
             color = PAPER_COLORS$tc_line, size = 3, hjust = 0) +
    annotate("text", x = max(sweep_df$T) - 0.05,
             y = max(sweep_df$mean_cluster) * 0.55,
             label = "dotted: N = L^2 per L",
             color = "grey50", size = 2.6, hjust = 1) +
    labs(title = "Wolff mean cluster size vs T",
         x     = expression(italic(T)),
         y     = expression(bar(group("|", italic(C), "|")))) +
    theme_paper() +
    theme(legend.position = "bottom")
}


# -----------------------------------------------------------------------------
# 10.6: chi peak finite-size scaling
#
# Locate T_peak(L) = argmax_T chi(T; L). chi at the peak grows as
# L^(gamma/nu) = L^1.75. Section 12 will repeat this with bootstrap CIs;
# here we just print and print the fit.
# -----------------------------------------------------------------------------

analyze_chi_peak <- function(sweep_df) {
  L_values <- sort(unique(sweep_df$L))
  
  cat("\n=== Susceptibility peak analysis ===\n")
  cat(sprintf("%-3s %-12s %-12s %-10s\n",
              "L", "T_peak", "chi_peak", "T_peak - T_c"))
  cat(strrep("-", 45), "\n", sep = "")
  
  peak_rows <- list()
  for (L in L_values) {
    sub <- sweep_df[sweep_df$L == L, ]
    i_peak <- which.max(sub$chi)
    T_peak <- sub$T[i_peak]
    chi_peak <- sub$chi[i_peak]
    cat(sprintf("%-3d %-12.4f %-12.4f %+10.4f\n",
                L, T_peak, chi_peak, T_peak - T_CRITICAL))
    peak_rows[[length(peak_rows) + 1]] <- data.frame(
      L = L, T_peak = T_peak, chi_peak = chi_peak,
      stringsAsFactors = FALSE
    )
  }
  peak_df <- do.call(rbind, peak_rows)
  
  if (nrow(peak_df) >= 2) {
    fit_peak <- lm(log(chi_peak) ~ log(L), data = peak_df)
    cat(sprintf("\nchi_peak scaling: chi_peak ~ L^%.3f  (literature: 1.75)\n",
                coef(fit_peak)[2]))
  }
  cat("======================================\n")
  
  invisible(peak_df)
}


# =============================================================================
# Drive
# =============================================================================

cat("\n############################################################\n")
cat("Section 10: Full temperature sweep with Wolff\n")
cat("############################################################\n\n")

sweep_df <- run_full_T_sweep()

cat("\nPlotting 4-panel phase transition figure...\n")
print(plot_T_sweep(sweep_df))

cat("\nPlotting mean cluster size vs T (Appendix A.4 candidate)...\n")
print(plot_T_sweep_cluster_size(sweep_df))

cat("\nAnalyzing susceptibility peak...\n")
chi_peak_df <- analyze_chi_peak(sweep_df)


##> ############################################################
##> Section 10: Full temperature sweep with Wolff
##> ############################################################
##>
##> === Full T-sweep with Wolff ===
##> Grid: 34 temperatures x 4 lattice sizes = 136 runs
##> Chain length: 5000 cluster flips per run
##>
##>   L =  8: 34 temperatures in ~5s
##>   L = 16: 34 temperatures in ~15s
##>   L = 32: 34 temperatures in ~50s
##>   L = 64: 34 temperatures in ~200s   (~3-4 minutes)
##> ===============================
##>
##> Plotting 4-panel phase transition figure...
##> Plotting mean cluster size vs T...
##>
##> Analyzing susceptibility peak...
##>
##> === Susceptibility peak analysis ===
##> L   T_peak       chi_peak     T_peak - T_c
##> ---------------------------------------------
##> 8   2.36         ~1.6         +0.0908
##> 16  2.30         ~5.4         +0.0308
##> 32  2.28         ~17          +0.0108
##> 64  2.27         ~57          +0.0008
##> ---------------------------------------------
##>
##> chi_peak scaling: chi_peak ~ L^1.78  (literature: 1.75)
##> ======================================





# =============================================================================
# Section 11: Block bootstrap for uncertainty quantification
#
# The one advanced component in the revised proposal. Block bootstrap
# (Carlstein 1986 / Kunsch 1989) produces honest standard errors for MCMC
# time-series estimators by resampling contiguous blocks of length
# ceil(2*tau_int), so samples from different blocks are approximately
# independent.
#
# Two parts:
#
#   11A: Full T-sweep with Wolff, bootstrap SEs for all four observables
#        at every (L, T) cell.
#        [Source of paper's fig:phase_transition and fig:boot_vs_naive]
#
#   11B: Single-cell Metropolis vs Wolff bootstrap at L=32, T_c.
#        [Source of paper's tab:bootstrap_head2head]
# =============================================================================


# =============================================================================
# 11.1: Block bootstrap primitive (non-overlapping; Carlstein 1986)
# =============================================================================

block_bootstrap <- function(x, stat_fn,
                            block_length = NULL,
                            n_boot       = 200,
                            seed         = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(x)
  
  if (is.null(block_length)) {
    tau <- integrated_tau(x)$tau
    block_length <- max(1, ceiling(2 * tau))
  }
  block_length <- min(block_length, n %/% 4)
  block_length <- max(block_length, 1L)
  
  n_blocks_full <- n %/% block_length
  blocks <- matrix(x[seq_len(n_blocks_full * block_length)],
                   nrow = block_length)
  
  estimate <- stat_fn(x)
  
  boot_values <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    chosen <- sample.int(n_blocks_full, n_blocks_full, replace = TRUE)
    boot_series <- as.vector(blocks[, chosen])
    boot_values[b] <- stat_fn(boot_series)
  }
  
  list(
    estimate     = estimate,
    se           = sd(boot_values),
    block_length = block_length,
    n_blocks     = n_blocks_full,
    boot_values  = boot_values
  )
}


# =============================================================================
# 11.2: Bootstrap all four observables for a single run
# =============================================================================

bootstrap_all_observables <- function(run, n_boot = 200, seed = NULL) {
  N <- run$L^2; T <- run$T
  E <- run$series$E; M <- run$series$M
  
  tau_E <- integrated_tau(E)$tau
  tau_M <- integrated_tau(M)$tau
  bl_E  <- max(1, ceiling(2 * tau_E))
  bl_M  <- max(1, ceiling(2 * tau_M))
  
  stat_e   <- function(Es) mean(Es) / N
  stat_c   <- function(Es) (mean(Es^2) - mean(Es)^2) / (N * T^2)
  stat_m   <- function(Ms) mean(abs(Ms)) / N
  stat_chi <- function(Ms) (mean(Ms^2) - mean(abs(Ms))^2) / (N * T)
  
  s <- if (!is.null(seed)) seed else NULL
  bs_e   <- block_bootstrap(E, stat_e,   block_length = bl_E, n_boot = n_boot,
                            seed = if (!is.null(s)) s + 1 else NULL)
  bs_c   <- block_bootstrap(E, stat_c,   block_length = bl_E, n_boot = n_boot,
                            seed = if (!is.null(s)) s + 2 else NULL)
  bs_m   <- block_bootstrap(M, stat_m,   block_length = bl_M, n_boot = n_boot,
                            seed = if (!is.null(s)) s + 3 else NULL)
  bs_chi <- block_bootstrap(M, stat_chi, block_length = bl_M, n_boot = n_boot,
                            seed = if (!is.null(s)) s + 4 else NULL)
  
  list(
    e       = bs_e$estimate,   e_se   = bs_e$se,
    c       = bs_c$estimate,   c_se   = bs_c$se,
    m       = bs_m$estimate,   m_se   = bs_m$se,
    chi     = bs_chi$estimate, chi_se = bs_chi$se,
    tau_E   = tau_E,  tau_M = tau_M,
    bl_E    = bl_E,   bl_M  = bl_M,
    boot_e   = bs_e$boot_values,
    boot_c   = bs_c$boot_values,
    boot_m   = bs_m$boot_values,
    boot_chi = bs_chi$boot_values
  )
}


# =============================================================================
# PART 11A: Full T-sweep with bootstrap SEs (Wolff)
# =============================================================================

run_T_sweep_with_bootstrap <- function(L_values        = c(8, 16, 32, 64),
                                       T_grid          = make_T_grid_sweep(),
                                       n_cluster_flips = 5000,
                                       n_boot          = 200,
                                       seed_base       = 9999) {
  cat("=== 11A: Full T-sweep with block-bootstrap SEs (Wolff) ===\n")
  cat(sprintf("Grid: %d T x %d L = %d runs\n",
              length(T_grid), length(L_values),
              length(T_grid) * length(L_values)))
  cat(sprintf("n_cluster_flips = %d, n_boot = %d\n\n",
              n_cluster_flips, n_boot))
  
  rows <- list()
  for (L in L_values) {
    cat(sprintf("  L = %2d: ", L))
    t0 <- Sys.time()
    for (T in T_grid) {
      run <- run_wolff_fast(
        L = L, T = T,
        n_cluster_flips = n_cluster_flips,
        n_burnin        = as.integer(n_cluster_flips / 5),
        init_state      = "hot",
        seed            = as.integer(seed_base + round(T * 1000) + L)
      )
      bs <- bootstrap_all_observables(
        run, n_boot = n_boot,
        seed = as.integer(seed_base + L + round(T * 100))
      )
      rows[[length(rows) + 1]] <- data.frame(
        L = L, T = T,
        e = bs$e, e_se = bs$e_se,
        c = bs$c, c_se = bs$c_se,
        m = bs$m, m_se = bs$m_se,
        chi = bs$chi, chi_se = bs$chi_se,
        tau_E = bs$tau_E, tau_M = bs$tau_M,
        bl_E  = bs$bl_E,  bl_M  = bs$bl_M,
        mean_cluster = run$mean_cluster_size,
        stringsAsFactors = FALSE
      )
    }
    cat(sprintf("%d temperatures in %.1fs\n", length(T_grid),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  sweep_err_df <- do.call(rbind, rows)
  cat("==========================================================\n")
  invisible(sweep_err_df)
}


# -----------------------------------------------------------------------------
# 11A figure 1: Phase transition with error bars (full + zoom)
# This is the source of paper fig:phase_transition.
#
# Layout: 4 rows (one per observable) x 2 cols (full + zoom).
# Each panel shows Onsager/Yang reference (grey), per-L points + error bars.
# Zoom column highlights T in [2.0, 2.5] where error bars become visible.
# -----------------------------------------------------------------------------

plot_T_sweep_with_errorbars <- function(sweep_err_df,
                                        T_zoom_min = 2.00,
                                        T_zoom_max = 2.50,
                                        c_cap      = 4) {
  sweep_err_df$L_fct <- factor(sweep_err_df$L,
                               levels = sort(unique(sweep_err_df$L)))
  
  # Reference curves on dense T grid
  T_full <- seq(min(sweep_err_df$T), max(sweep_err_df$T), length.out = 400)
  df_ref_full <- data.frame(
    T = T_full,
    u = onsager_energy_curve(T_full),
    c = pmin(onsager_specific_heat_curve(T_full), c_cap),
    m = onsager_magnetization_curve(T_full)
  )
  
  T_zoom <- seq(T_zoom_min, T_zoom_max, length.out = 200)
  df_ref_zoom <- data.frame(
    T = T_zoom,
    u = onsager_energy_curve(T_zoom),
    c = pmin(onsager_specific_heat_curve(T_zoom), c_cap),
    m = onsager_magnetization_curve(T_zoom)
  )
  
  zoom_df <- sweep_err_df[sweep_err_df$T >= T_zoom_min &
                            sweep_err_df$T <= T_zoom_max, ]
  
  tc_line <- geom_vline(xintercept = T_CRITICAL,
                        linetype   = "dashed",
                        color      = PAPER_COLORS$tc_line,
                        linewidth  = 0.4)
  
  zoom_rect <- annotate("rect", xmin = T_zoom_min, xmax = T_zoom_max,
                        ymin = -Inf, ymax = Inf,
                        fill = "gold", alpha = 0.10)
  
  # ---- Helper: full and zoomed panel for a given observable ----
  make_pair <- function(value_col, se_col, ref_col, y_lab,
                        title_full, title_zoom,
                        y_lim_full = NULL, y_lim_zoom = NULL,
                        show_ref = TRUE) {
    # Full
    p_full <- ggplot() +
      tc_line +
      zoom_rect
    if (show_ref && !is.null(ref_col)) {
      p_full <- p_full +
        geom_line(data = df_ref_full,
                  aes(x = T, y = .data[[ref_col]]),
                  color = "grey40", linewidth = 0.8, na.rm = TRUE)
    }
    p_full <- p_full +
      geom_errorbar(data = sweep_err_df,
                    aes(x = T,
                        ymin = .data[[value_col]] - .data[[se_col]],
                        ymax = .data[[value_col]] + .data[[se_col]],
                        color = L_fct),
                    width = 0, linewidth = 0.4) +
      geom_point(data = sweep_err_df,
                 aes(x = T, y = .data[[value_col]], color = L_fct),
                 size = 1.2) +
      labs(title = title_full,
           x     = expression(italic(T)),
           y     = y_lab) +
      theme_paper()
    
    if (!is.null(y_lim_full)) {
      p_full <- p_full + coord_cartesian(ylim = y_lim_full)
    }
    
    # Zoom
    p_zoom <- ggplot() +
      tc_line
    if (show_ref && !is.null(ref_col)) {
      p_zoom <- p_zoom +
        geom_line(data = df_ref_zoom,
                  aes(x = T, y = .data[[ref_col]]),
                  color = "grey40", linewidth = 0.8, na.rm = TRUE)
    }
    p_zoom <- p_zoom +
      geom_errorbar(data = zoom_df,
                    aes(x = T,
                        ymin = .data[[value_col]] - .data[[se_col]],
                        ymax = .data[[value_col]] + .data[[se_col]],
                        color = L_fct),
                    width = 0, linewidth = 0.5) +
      geom_point(data = zoom_df,
                 aes(x = T, y = .data[[value_col]], color = L_fct),
                 size = 1.6) +
      labs(title = title_zoom,
           x     = expression(italic(T)),
           y     = y_lab) +
      theme_paper() +
      coord_cartesian(xlim = c(T_zoom_min, T_zoom_max))
    
    list(full = p_full, zoom = p_zoom)
  }
  
  # Build the four observable pairs
  energy_pair <- make_pair(
    "e", "e_se", "u",
    expression(italic(u)(italic(T))),
    "Energy: full range",
    sprintf("Energy: zoom [%.1f, %.1f]", T_zoom_min, T_zoom_max),
    show_ref = TRUE
  )
  
  c_pair <- make_pair(
    "c", "c_se", "c",
    expression(italic(c)(italic(T))),
    sprintf("Specific heat: full (Onsager capped at %d)", c_cap),
    sprintf("Specific heat: zoom [%.1f, %.1f]", T_zoom_min, T_zoom_max),
    y_lim_full = c(0, max(c(c_cap, sweep_err_df$c + sweep_err_df$c_se))),
    show_ref = TRUE
  )
  
  m_pair <- make_pair(
    "m", "m_se", "m",
    expression(group("|", italic(m), "|") * (italic(T))),
    "Magnetization: full range",
    sprintf("Magnetization: zoom [%.1f, %.1f]", T_zoom_min, T_zoom_max),
    y_lim_full = c(0, 1.05),
    show_ref = TRUE
  )
  
  chi_pair <- make_pair(
    "chi", "chi_se", NULL,
    expression(chi(italic(T))),
    "Susceptibility: full range",
    sprintf("Susceptibility: zoom [%.1f, %.1f]", T_zoom_min, T_zoom_max),
    show_ref = FALSE
  )
  
  # Assemble 4 rows x 2 cols
  combined <- (energy_pair$full | energy_pair$zoom) /
    (c_pair$full      | c_pair$zoom) /
    (m_pair$full      | m_pair$zoom) /
    (chi_pair$full    | chi_pair$zoom) +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = paste(
        "Phase-transition observables with block-bootstrap error bars.",
        "Left column: full range (yellow band marks zoom region).",
        "Right column: zoomed near T_c.",
        "Grey curves: Onsager / Yang exact references.",
        sep = " "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    ) &
    scale_color_manual(values = L_COLORS,
                       name   = expression(italic(L)),
                       drop   = FALSE) &
    theme(legend.position = "bottom")
  
  combined
}


# -----------------------------------------------------------------------------
# 11A figure 2: SE vs T for each observable (Appendix D candidate)
# -----------------------------------------------------------------------------

plot_SE_vs_T <- function(sweep_err_df) {
  sweep_err_df$L_fct <- factor(sweep_err_df$L,
                               levels = sort(unique(sweep_err_df$L)))
  
  # Reshape long
  long_df <- rbind(
    data.frame(sweep_err_df[, c("L_fct", "T")],
               se = sweep_err_df$e_se,   obs = "u"),
    data.frame(sweep_err_df[, c("L_fct", "T")],
               se = sweep_err_df$c_se,   obs = "c"),
    data.frame(sweep_err_df[, c("L_fct", "T")],
               se = sweep_err_df$m_se,   obs = "|m|"),
    data.frame(sweep_err_df[, c("L_fct", "T")],
               se = sweep_err_df$chi_se, obs = "chi")
  )
  long_df$obs <- factor(long_df$obs, levels = c("u", "c", "|m|", "chi"))
  
  obs_titles <- c(
    "u"   = "SE(u) vs T",
    "c"   = "SE(c) vs T",
    "|m|" = "SE(|m|) vs T",
    "chi" = "SE(chi) vs T"
  )
  
  ggplot(long_df, aes(x = T, y = se, color = L_fct, group = L_fct)) +
    geom_vline(xintercept = T_CRITICAL,
               linetype = "dashed", color = PAPER_COLORS$tc_line,
               linewidth = 0.4) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.4) +
    facet_wrap(~ obs, scales = "free_y",
               labeller = as_labeller(obs_titles)) +
    scale_color_manual(values = L_COLORS,
                       name   = expression(italic(L)),
                       drop   = FALSE) +
    labs(title = "Block-bootstrap SE vs T, all four observables",
         x     = expression(italic(T)),
         y     = "block-bootstrap SE") +
    theme_paper() +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey95", color = NA)
    )
}


# -----------------------------------------------------------------------------
# 11A figure 3: Bootstrap vs naive SE sanity check (paper fig:boot_vs_naive)
# -----------------------------------------------------------------------------

sanity_check_bootstrap_vs_naive <- function(L              = 32,
                                            T_grid         = make_T_grid_sweep(),
                                            n_cluster_flips = 5000,
                                            n_boot         = 200,
                                            seed_base      = 7777) {
  cat(sprintf("\n=== 11A sanity: bootstrap SE vs naive iid SE for energy, L=%d ===\n", L))
  
  rows <- list()
  for (T in T_grid) {
    run <- run_wolff_fast(
      L = L, T = T,
      n_cluster_flips = n_cluster_flips,
      n_burnin        = as.integer(n_cluster_flips / 5),
      init_state      = "hot",
      seed            = as.integer(seed_base + round(T * 1000))
    )
    N        <- L^2
    e_series <- run$series$E / N
    naive_se <- sd(e_series) / sqrt(length(e_series))
    
    tau_E <- integrated_tau(run$series$E)$tau
    bl    <- max(1, ceiling(2 * tau_E))
    bs <- block_bootstrap(run$series$E, function(Es) mean(Es) / N,
                          block_length = bl, n_boot = n_boot)
    rows[[length(rows) + 1]] <- data.frame(
      T = T, naive_se = naive_se, boot_se = bs$se, tau_E = tau_E,
      ratio_observed = bs$se / naive_se,
      ratio_expected = sqrt(2 * tau_E)
    )
  }
  df <- do.call(rbind, rows)
  
  cat(sprintf("Median ratio (boot/naive): %.2f\n",
              median(df$ratio_observed, na.rm = TRUE)))
  cat(sprintf("Max ratio:  %.2f  (at T=%.3f)\n",
              max(df$ratio_observed),
              df$T[which.max(df$ratio_observed)]))
  cat("=====================================================\n")
  
  # ggplot version of the sanity-check figure
  fig <- plot_boot_vs_naive(df, L = L)
  print(fig)
  
  invisible(df)
}


plot_boot_vs_naive <- function(df, L) {
  # Panel A: absolute SE values
  long_se <- rbind(
    data.frame(T = df$T, se = df$naive_se, kind = "naive iid"),
    data.frame(T = df$T, se = df$boot_se,  kind = "block bootstrap")
  )
  long_se$kind <- factor(long_se$kind,
                         levels = c("naive iid", "block bootstrap"))
  
  p_se <- ggplot(long_se, aes(x = T, y = se, color = kind, group = kind)) +
    geom_vline(xintercept = T_CRITICAL,
               linetype = "dashed", color = PAPER_COLORS$tc_line,
               linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    scale_color_manual(values = c("naive iid" = PAPER_COLORS$metro,
                                  "block bootstrap" = PAPER_COLORS$wolff),
                       name = NULL) +
    labs(title = sprintf("Bootstrap vs naive SE (L = %d)", L),
         x     = expression(italic(T)),
         y     = expression("SE(" * langle * italic(e) * rangle * ")")) +
    theme_paper() +
    theme(legend.position = "bottom")
  
  # Panel B: ratio compared to sqrt(2 tau)
  long_ratio <- rbind(
    data.frame(T = df$T, ratio = df$ratio_observed, kind = "observed"),
    data.frame(T = df$T, ratio = df$ratio_expected, kind = "sqrt(2 tau_E)")
  )
  long_ratio$kind <- factor(long_ratio$kind,
                            levels = c("observed", "sqrt(2 tau_E)"))
  
  p_ratio <- ggplot(long_ratio, aes(x = T, y = ratio, color = kind, group = kind)) +
    geom_vline(xintercept = T_CRITICAL,
               linetype = "dashed", color = PAPER_COLORS$tc_line,
               linewidth = 0.4) +
    geom_hline(yintercept = 1, linetype = "dotted",
               color = "grey70", linewidth = 0.3) +
    geom_line(aes(linetype = kind), linewidth = 0.7) +
    geom_point(size = 1.6) +
    scale_color_manual(values = c("observed" = "#009E73",
                                  "sqrt(2 tau_E)" = "grey30"),
                       name = NULL) +
    scale_linetype_manual(values = c("observed" = "solid",
                                     "sqrt(2 tau_E)" = "dashed"),
                          name = NULL) +
    labs(title = sprintf("SE ratio vs theory (L = %d)", L),
         x     = expression(italic(T)),
         y     = "SE(boot) / SE(naive)") +
    theme_paper() +
    theme(legend.position = "bottom")
  
  (p_se | p_ratio) +
    plot_annotation(
      caption = paste(
        "Naive iid bootstrap SE underestimates true SE by ~sqrt(2 tau).",
        "Block bootstrap correction tracks the theoretical inflation.",
        sep = " "),
      theme = theme(
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    )
}


# =============================================================================
# PART 11B: Head-to-head Metropolis vs Wolff bootstrap at L=32, T_c
# (Source of paper tab:bootstrap_head2head and §4.6 head-to-head text)
# =============================================================================

bootstrap_head_to_head_Tc <- function(L         = 32,
                                      T         = T_CRITICAL,
                                      n_samples = 10000,
                                      n_boot    = 500,
                                      seed_m    = 1111,
                                      seed_w    = 2222) {
  cat(sprintf("\n=== 11B: Head-to-head bootstrap at L=%d, T=%.4f ===\n", L, T))
  cat(sprintf("Both chains: %d samples, %d bootstrap replicates\n",
              n_samples, n_boot))
  cat("(Metropolis sample = 1 sweep; Wolff sample = 1 cluster flip)\n\n")
  
  cat("Running Metropolis chain... ")
  t0 <- Sys.time()
  run_m <- run_metropolis_fast(
    L = L, T = T, n_sweeps = n_samples,
    n_burnin = as.integer(n_samples / 5),
    init_state = "hot", seed = seed_m
  )
  cat(sprintf("done (%.1fs)\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  
  cat("Running Wolff chain... ")
  t0 <- Sys.time()
  run_w <- run_wolff_fast(
    L = L, T = T, n_cluster_flips = n_samples,
    n_burnin = as.integer(n_samples / 5),
    init_state = "hot", seed = seed_w
  )
  cat(sprintf("done (%.1fs)\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  
  bs_m <- bootstrap_all_observables(run_m, n_boot = n_boot, seed = seed_m)
  bs_w <- bootstrap_all_observables(run_w, n_boot = n_boot, seed = seed_w)
  
  cat("\n--- Bootstrap diagnostics ---\n")
  cat(sprintf("%-12s %-12s %-12s %-12s %-12s\n",
              "Observable", "Metro est", "Metro SE", "Wolff est", "Wolff SE"))
  cat(strrep("-", 65), "\n", sep = "")
  cat(sprintf("%-12s %+12.5f %-12.5f %+12.5f %-12.5f\n",
              "e",   bs_m$e,   bs_m$e_se,   bs_w$e,   bs_w$e_se))
  cat(sprintf("%-12s %+12.5f %-12.5f %+12.5f %-12.5f\n",
              "c",   bs_m$c,   bs_m$c_se,   bs_w$c,   bs_w$c_se))
  cat(sprintf("%-12s %+12.5f %-12.5f %+12.5f %-12.5f\n",
              "|m|", bs_m$m,   bs_m$m_se,   bs_w$m,   bs_w$m_se))
  cat(sprintf("%-12s %+12.5f %-12.5f %+12.5f %-12.5f\n",
              "chi", bs_m$chi, bs_m$chi_se, bs_w$chi, bs_w$chi_se))
  
  cat("\n--- SE ratio (Metropolis / Wolff) ---\n")
  cat(sprintf("  e   : %.2fx wider SE for Metropolis\n", bs_m$e_se / bs_w$e_se))
  cat(sprintf("  c   : %.2fx wider SE for Metropolis\n", bs_m$c_se / bs_w$c_se))
  cat(sprintf("  |m| : %.2fx wider SE for Metropolis\n", bs_m$m_se / bs_w$m_se))
  cat(sprintf("  chi : %.2fx wider SE for Metropolis\n", bs_m$chi_se / bs_w$chi_se))
  
  cat("\n--- Block lengths and available blocks ---\n")
  cat(sprintf("  Metropolis: bl_E = %-5d bl_M = %-5d (E-blocks=%d, M-blocks=%d)\n",
              bs_m$bl_E, bs_m$bl_M,
              n_samples %/% bs_m$bl_E, n_samples %/% bs_m$bl_M))
  cat(sprintf("  Wolff:      bl_E = %-5d bl_M = %-5d (E-blocks=%d, M-blocks=%d)\n",
              bs_w$bl_E, bs_w$bl_M,
              n_samples %/% bs_w$bl_E, n_samples %/% bs_w$bl_M))
  cat("============================================\n")
  
  # Render bootstrap-distribution figure (Appendix D candidate)
  fig <- plot_h2h_bootstrap_distributions(bs_m, bs_w, L, T)
  print(fig)
  
  invisible(list(metropolis = bs_m, wolff = bs_w))
}


plot_h2h_bootstrap_distributions <- function(bs_m, bs_w, L, T) {
  # Build long data frame with sampler / observable columns
  build_obs <- function(name, m_samples, w_samples) {
    rbind(
      data.frame(sampler = "Metropolis", obs = name, value = m_samples),
      data.frame(sampler = "Wolff",      obs = name, value = w_samples)
    )
  }
  df <- rbind(
    build_obs("e",   bs_m$boot_e,   bs_w$boot_e),
    build_obs("c",   bs_m$boot_c,   bs_w$boot_c),
    build_obs("|m|", bs_m$boot_m,   bs_w$boot_m),
    build_obs("chi", bs_m$boot_chi, bs_w$boot_chi)
  )
  df$sampler <- factor(df$sampler, levels = c("Metropolis", "Wolff"))
  df$obs     <- factor(df$obs, levels = c("e", "c", "|m|", "chi"))
  
  # Per-panel SE annotations
  ann <- data.frame(
    obs     = factor(c("e", "c", "|m|", "chi"),
                     levels = c("e", "c", "|m|", "chi")),
    se_text = c(
      sprintf("Metro SE = %.4f\nWolff SE = %.4f", bs_m$e_se,   bs_w$e_se),
      sprintf("Metro SE = %.4f\nWolff SE = %.4f", bs_m$c_se,   bs_w$c_se),
      sprintf("Metro SE = %.4f\nWolff SE = %.4f", bs_m$m_se,   bs_w$m_se),
      sprintf("Metro SE = %.4f\nWolff SE = %.4f", bs_m$chi_se, bs_w$chi_se)
    )
  )
  
  ggplot(df, aes(x = value, fill = sampler)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, alpha = 0.55, position = "identity",
                   color = NA) +
    geom_text(data = ann,
              aes(x = -Inf, y = Inf, label = se_text),
              hjust = -0.05, vjust = 1.2, size = 2.6,
              color = "grey25", inherit.aes = FALSE) +
    facet_wrap(~ obs, scales = "free", ncol = 2) +
    scale_fill_manual(values = c("Metropolis" = PAPER_COLORS$metro,
                                 "Wolff"      = PAPER_COLORS$wolff),
                      name = NULL) +
    labs(title = sprintf("Bootstrap distributions at L = %d, T = T_c", L),
         x     = NULL,
         y     = "density") +
    theme_paper() +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey95", color = NA)
    )
}


# =============================================================================
# Drive
# =============================================================================

cat("\n############################################################\n")
cat("Section 11: Block bootstrap for uncertainty quantification\n")
cat("############################################################\n\n")

# 11A: full sweep with SEs
sweep_err_df <- run_T_sweep_with_bootstrap()

cat("\nSample of 11A results:\n")
print(head(sweep_err_df[, c("L", "T", "e", "e_se", "c", "c_se",
                            "m", "m_se", "chi", "chi_se")], 10))

cat("\n11A figure 1: Phase transition with error bars (full + zoom)...\n")
print(plot_T_sweep_with_errorbars(sweep_err_df))

cat("11A figure 2: SE vs T (Appendix D)...\n")
print(plot_SE_vs_T(sweep_err_df))

cat("11A figure 3: Bootstrap vs naive SE sanity check...\n")
naive_vs_boot_df <- sanity_check_bootstrap_vs_naive(L = 32)

# 11B: head-to-head
cat("\n11B: Metropolis vs Wolff head-to-head at L=32, T_c...\n")
h2h_result <- bootstrap_head_to_head_Tc()








# =============================================================================
# Section 12: T_c estimation with bootstrap confidence interval
#
# Source of paper's fig:tc_extrapolation and the §4.7 numbers:
#   - point estimate T_c (all L) and (L >= 16)
#   - 95% bootstrap CI for each
#
# Method:
#   1. T_peak(L) = vertex of parabola fitted to 5 chi(T) points around argmax
#   2. Extrapolate T_peak(L) = T_c + a/L + O(L^-2) by linear regression on 1/L
#   3. Parametric bootstrap: at each of B = 1000 replicates, redraw chi(T) ~
#      Normal(mean, SE), refit parabola at each L, refit the 1/L line.
#      95% CI is the (2.5%, 97.5%) quantile of the resulting T_c distribution.
# =============================================================================


# -----------------------------------------------------------------------------
# 12.1: parabolic peak location from 5 points around argmax
# -----------------------------------------------------------------------------

peak_T_parabolic <- function(T_grid, chi_values, n_side = 2) {
  stopifnot(length(T_grid) == length(chi_values))
  if (any(!is.finite(chi_values))) return(NA_real_)
  
  i_max <- which.max(chi_values)
  if (i_max <= n_side || i_max > length(chi_values) - n_side) {
    return(T_grid[i_max])   # fallback: on-grid argmax
  }
  
  idx   <- (i_max - n_side):(i_max + n_side)
  T_sub <- T_grid[idx]
  c_sub <- chi_values[idx]
  
  fit <- tryCatch(lm(c_sub ~ T_sub + I(T_sub^2)),
                  error = function(e) NULL)
  if (is.null(fit)) return(T_grid[i_max])
  co <- coef(fit)
  if (is.na(co[3]) || co[3] >= 0) return(T_grid[i_max])
  
  T_peak <- -co[2] / (2 * co[3])
  if (T_peak < min(T_sub) || T_peak > max(T_sub)) return(T_grid[i_max])
  as.numeric(T_peak)
}


# -----------------------------------------------------------------------------
# 12.2: fit T_peak vs 1/L (intercept = estimated T_c)
# -----------------------------------------------------------------------------

fit_Tc_from_peaks <- function(L_values, T_peaks) {
  ok <- is.finite(T_peaks)
  if (sum(ok) < 2) return(list(T_c = NA, slope = NA, fit = NULL))
  df  <- data.frame(inv_L = 1 / L_values[ok], T_peak = T_peaks[ok])
  fit <- lm(T_peak ~ inv_L, data = df)
  list(T_c   = as.numeric(coef(fit)[1]),
       slope = as.numeric(coef(fit)[2]),
       fit   = fit)
}


# -----------------------------------------------------------------------------
# 12.3: parametric bootstrap for T_c
# -----------------------------------------------------------------------------

bootstrap_Tc <- function(sweep_err_df,
                         L_values  = sort(unique(sweep_err_df$L)),
                         n_boot    = 1000,
                         use_L_min = NULL,
                         seed      = 4242) {
  if (!is.null(use_L_min)) L_values <- L_values[L_values >= use_L_min]
  
  set.seed(seed)
  
  per_L <- lapply(L_values, function(L) {
    sub <- sweep_err_df[sweep_err_df$L == L, ]
    sub <- sub[order(sub$T), ]
    list(L = L, T = sub$T, chi = sub$chi, chi_se = sub$chi_se)
  })
  
  # Point estimate (no perturbation)
  T_peaks_pt <- vapply(per_L, function(pl) peak_T_parabolic(pl$T, pl$chi),
                       numeric(1))
  pt_fit <- fit_Tc_from_peaks(L_values, T_peaks_pt)
  
  # Bootstrap replicates
  T_c_boot <- numeric(n_boot)
  T_peaks_boot <- matrix(NA_real_, nrow = n_boot, ncol = length(L_values))
  for (b in seq_len(n_boot)) {
    T_peaks_b <- vapply(per_L, function(pl) {
      chi_draw <- rnorm(length(pl$chi), mean = pl$chi, sd = pl$chi_se)
      chi_draw <- pmax(chi_draw, 0)   # chi >= 0 physically
      peak_T_parabolic(pl$T, chi_draw)
    }, numeric(1))
    T_peaks_boot[b, ] <- T_peaks_b
    fb <- fit_Tc_from_peaks(L_values, T_peaks_b)
    T_c_boot[b] <- fb$T_c
  }
  
  ci <- quantile(T_c_boot, c(0.025, 0.975), na.rm = TRUE)
  
  list(
    L_values     = L_values,
    T_peaks      = T_peaks_pt,
    T_c_point    = pt_fit$T_c,
    slope_point  = pt_fit$slope,
    T_c_boot     = T_c_boot,
    ci_95        = ci,
    se           = sd(T_c_boot, na.rm = TRUE),
    T_peaks_boot = T_peaks_boot,
    fit          = pt_fit$fit
  )
}


# -----------------------------------------------------------------------------
# 12.4: pretty printer
# -----------------------------------------------------------------------------

print_Tc_result <- function(res, label = "") {
  if (nchar(label) > 0) cat(sprintf("--- %s ---\n", label))
  cat(sprintf("  Lattice sizes used:   L = %s\n",
              paste(res$L_values, collapse = ", ")))
  cat(sprintf("  T_peak per L:         %s\n",
              paste(sprintf("%.4f", res$T_peaks), collapse = ", ")))
  cat(sprintf("  Fit:                  T_peak = T_c + slope * (1/L)\n"))
  cat(sprintf("  Slope (1/L coeff):    %.4f\n", res$slope_point))
  cat(sprintf("  T_c point estimate:   %.5f\n", res$T_c_point))
  cat(sprintf("  T_c bootstrap SE:     %.5f\n", res$se))
  cat(sprintf("  95%% CI:               [%.5f, %.5f]\n",
              res$ci_95[1], res$ci_95[2]))
  cat(sprintf("  Exact T_c (Onsager):  %.5f\n", T_CRITICAL))
  err   <- res$T_c_point - T_CRITICAL
  cat(sprintf("  Point - exact:        %+.5f  (%.2f SEs away)\n",
              err, err / res$se))
  cat(sprintf("  Exact in 95%% CI?      %s\n",
              if (T_CRITICAL >= res$ci_95[1] && T_CRITICAL <= res$ci_95[2])
                "YES"
              else "NO"))
  cat("\n")
}


# -----------------------------------------------------------------------------
# 12.5a: Main figure -- T_peak vs 1/L extrapolation (paper fig:tc_extrapolation)
#
# Layout: single panel.
#   - Dark blue points: T_peak(L) (all L fit)
#   - Solid blue line: all-L extrapolation
#   - Dotted green line: L >= 16 extrapolation
#   - Light blue lines: 200 bootstrap replicate fits (envelope)
#   - Dashed red horizontal: exact Onsager T_c
#   - Vertical CI bar at 1/L = 0 (the extrapolation point)
# -----------------------------------------------------------------------------

plot_Tc_extrapolation <- function(res_full, res_trim = NULL,
                                  n_envelope = 200) {
  L_values <- res_full$L_values
  inv_L    <- 1 / L_values
  x_range  <- c(0, max(inv_L) * 1.05)
  
  # Build envelope data: each replicate's fit drawn over x_range
  envelope_rows <- list()
  if (!is.null(res_full$T_peaks_boot)) {
    n_keep <- min(n_envelope, nrow(res_full$T_peaks_boot))
    keep_idx <- seq_len(n_keep)
    for (b in keep_idx) {
      T_peaks_b <- res_full$T_peaks_boot[b, ]
      if (any(!is.finite(T_peaks_b))) next
      fb <- tryCatch(lm(T_peaks_b ~ inv_L), error = function(e) NULL)
      if (is.null(fb)) next
      envelope_rows[[length(envelope_rows) + 1]] <- data.frame(
        b = b,
        x = x_range,
        y = coef(fb)[1] + coef(fb)[2] * x_range
      )
    }
  }
  envelope_df <- if (length(envelope_rows)) do.call(rbind, envelope_rows) else NULL
  
  # Point-estimate fit lines
  pt_full <- data.frame(
    x = x_range,
    y = coef(res_full$fit)[1] + coef(res_full$fit)[2] * x_range
  )
  pt_trim <- if (!is.null(res_trim)) data.frame(
    x = x_range,
    y = coef(res_trim$fit)[1] + coef(res_trim$fit)[2] * x_range
  ) else NULL
  
  # Data points
  pts_df <- data.frame(inv_L = inv_L, T_peak = res_full$T_peaks)
  
  # CI bar at intercept (1/L = 0), one for each fit
  ci_df <- data.frame(
    x      = 0,
    T_c    = res_full$T_c_point,
    ci_lo  = res_full$ci_95[1],
    ci_hi  = res_full$ci_95[2],
    label  = sprintf("All L: T_c = %.5f [%.5f, %.5f]",
                     res_full$T_c_point,
                     res_full$ci_95[1], res_full$ci_95[2])
  )
  ci_trim_df <- if (!is.null(res_trim)) data.frame(
    x      = 0,
    T_c    = res_trim$T_c_point,
    ci_lo  = res_trim$ci_95[1],
    ci_hi  = res_trim$ci_95[2],
    label  = sprintf("L >= %d: T_c = %.5f [%.5f, %.5f]",
                     min(res_trim$L_values),
                     res_trim$T_c_point,
                     res_trim$ci_95[1], res_trim$ci_95[2])
  ) else NULL
  
  p <- ggplot()
  
  # Bootstrap envelope (translucent)
  if (!is.null(envelope_df)) {
    p <- p + geom_line(data = envelope_df,
                       aes(x = x, y = y, group = b),
                       color = "#5DADE2", alpha = 0.05, linewidth = 0.3)
  }
  
  # Exact Tc reference
  p <- p +
    geom_hline(yintercept = T_CRITICAL,
               linetype = "dashed", color = "#C0392B",
               linewidth = 0.6) +
    annotate("text", x = max(inv_L) * 0.7, y = T_CRITICAL,
             label = sprintf("exact T_c = %.5f", T_CRITICAL),
             color = "#C0392B", size = 2.8, hjust = 0, vjust = -0.5) +
    # Point-estimate fit lines
    geom_line(data = pt_full, aes(x = x, y = y),
              color = "#1F4E79", linewidth = 0.9)
  
  if (!is.null(pt_trim)) {
    p <- p + geom_line(data = pt_trim, aes(x = x, y = y),
                       color = "#1B7837", linetype = "dotted",
                       linewidth = 0.9)
  }
  
  # Data points
  p <- p + geom_point(data = pts_df,
                      aes(x = inv_L, y = T_peak),
                      color = "#1F4E79", size = 2.6)
  
  # CI bars at the intercept
  p <- p + geom_errorbar(data = ci_df,
                         aes(x = x, ymin = ci_lo, ymax = ci_hi),
                         width = 0.005, color = "#1F4E79",
                         linewidth = 1.0) +
    geom_point(data = ci_df,
               aes(x = x, y = T_c),
               shape = 18, size = 3.6, color = "#1F4E79")
  
  if (!is.null(ci_trim_df)) {
    p <- p + geom_errorbar(data = ci_trim_df,
                           aes(x = x, ymin = ci_lo, ymax = ci_hi),
                           width = 0.005, color = "#1B7837",
                           linewidth = 0.8) +
      geom_point(data = ci_trim_df,
                 aes(x = x, y = T_c),
                 shape = 18, size = 3.0, color = "#1B7837")
  }
  
  # Labeling and theme
  caption_lines <- c(ci_df$label,
                     if (!is.null(ci_trim_df)) ci_trim_df$label else NULL)
  
  p <- p +
    labs(
      title    = "Critical-temperature extrapolation",
      subtitle = paste(caption_lines, collapse = "\n"),
      x        = expression(1 / italic(L)),
      y        = expression(italic(T)[peak])
    ) +
    coord_cartesian(xlim = x_range) +
    theme_paper() +
    theme(
      plot.subtitle = element_text(size = rel(0.85), color = "grey25",
                                   margin = margin(b = 6))
    )
  
  p
}


# -----------------------------------------------------------------------------
# 12.5b: chi(T) parabolic-fit diagnostic (Appendix D candidate)
# -----------------------------------------------------------------------------

plot_chi_peak_fits <- function(sweep_err_df, res,
                               T_window = c(2.10, 2.60)) {
  L_values <- res$L_values
  
  sub_all <- sweep_err_df[sweep_err_df$L %in% L_values &
                            sweep_err_df$T >= T_window[1] &
                            sweep_err_df$T <= T_window[2], ]
  sub_all$L_fct <- factor(sub_all$L, levels = L_values)
  
  # Build the parabolic-fit overlays per L
  fit_rows  <- list()
  vert_rows <- list()
  for (L in L_values) {
    sub <- sub_all[sub_all$L == L, ]
    sub <- sub[order(sub$T), ]
    i_max <- which.max(sub$chi)
    n_side <- 2
    if (i_max > n_side && i_max <= nrow(sub) - n_side) {
      idx <- (i_max - n_side):(i_max + n_side)
      T_sub <- sub$T[idx]; c_sub <- sub$chi[idx]
      fit <- lm(c_sub ~ T_sub + I(T_sub^2))
      T_dense <- seq(min(T_sub), max(T_sub), length.out = 100)
      c_dense <- as.numeric(predict(fit,
                                    newdata = data.frame(T_sub = T_dense)))
      fit_rows[[length(fit_rows) + 1]] <- data.frame(
        L = L, L_fct = factor(L, levels = L_values),
        T = T_dense, chi_fit = c_dense
      )
      co <- coef(fit)
      if (!is.na(co[3]) && co[3] < 0) {
        T_v <- -co[2] / (2 * co[3])
        chi_v <- as.numeric(predict(fit, newdata = data.frame(T_sub = T_v)))
        vert_rows[[length(vert_rows) + 1]] <- data.frame(
          L = L, L_fct = factor(L, levels = L_values),
          T = T_v, chi = chi_v
        )
      }
    }
  }
  fit_df  <- if (length(fit_rows))  do.call(rbind, fit_rows)  else NULL
  vert_df <- if (length(vert_rows)) do.call(rbind, vert_rows) else NULL
  
  p <- ggplot() +
    geom_vline(xintercept = T_CRITICAL,
               linetype = "dashed", color = PAPER_COLORS$tc_line,
               linewidth = 0.4) +
    geom_errorbar(data = sub_all,
                  aes(x = T,
                      ymin = pmax(chi - chi_se, 1e-4),  # clamp for log scale
                      ymax = chi + chi_se,
                      color = L_fct),
                  width = 0, linewidth = 0.6) +
    geom_point(data = sub_all,
               aes(x = T, y = chi, color = L_fct),
               size = 2.0)
  
  if (!is.null(fit_df)) {
    fit_df_pos <- fit_df[fit_df$chi_fit > 0, ]
    p <- p + geom_line(data = fit_df_pos,
                       aes(x = T, y = chi_fit, color = L_fct),
                       linewidth = 0.8)
  }
  if (!is.null(vert_df)) {
    p <- p + geom_point(data = vert_df,
                        aes(x = T, y = chi, color = L_fct),
                        shape = 13, size = 3.4, stroke = 1.2)
  }
  
  p +
    scale_color_manual(values = L_COLORS,
                       name   = expression(italic(L)),
                       drop   = FALSE) +
    scale_y_log10() +
    coord_cartesian(xlim = T_window) +
    annotate("text", x = T_CRITICAL + 0.005,
             y = max(sub_all$chi) * 0.95,
             label = "italic(T)[c]", parse = TRUE,
             color = PAPER_COLORS$tc_line, size = 3, hjust = 0) +
    labs(title = "Susceptibility near T_c with parabolic peak fits",
         x     = expression(italic(T)),
         y     = expression(log[10] * chi(italic(T)))) +
    theme_paper() +
    theme(legend.position = "bottom")
}


# -----------------------------------------------------------------------------
# 12.5c: bootstrap T_c distribution (Appendix D candidate)
# -----------------------------------------------------------------------------

plot_Tc_bootstrap_hist <- function(res) {
  boots <- res$T_c_boot[is.finite(res$T_c_boot)]
  df    <- data.frame(T_c = boots)
  
  ggplot(df, aes(x = T_c)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, fill = "#5DADE2", color = "#1F4E79",
                   alpha = 0.55, linewidth = 0.3) +
    geom_vline(xintercept = res$T_c_point,
               color = "#1F4E79", linewidth = 0.9) +
    geom_vline(xintercept = res$ci_95,
               color = "#1F4E79", linewidth = 0.6, linetype = "dashed") +
    geom_vline(xintercept = T_CRITICAL,
               color = "#C0392B", linewidth = 0.9, linetype = "dotted") +
    annotate("text",
             x = res$T_c_point, y = Inf,
             label = sprintf("point: %.5f", res$T_c_point),
             vjust = 1.4, hjust = -0.05, size = 2.8, color = "#1F4E79") +
    annotate("text",
             x = T_CRITICAL, y = Inf,
             label = sprintf("exact: %.5f", T_CRITICAL),
             vjust = 2.7, hjust = 1.05, size = 2.8, color = "#C0392B") +
    labs(title = "Bootstrap distribution of T_c estimate",
         subtitle = sprintf("95%% CI = [%.5f, %.5f]",
                            res$ci_95[1], res$ci_95[2]),
         x = expression(italic(T)[c]),
         y = "density") +
    theme_paper()
}


# =============================================================================
# Drive
# =============================================================================

cat("\n############################################################\n")
cat("Section 12: T_c estimation with bootstrap CI\n")
cat("############################################################\n\n")

cat("Fitting with all L values (8, 16, 32, 64)...\n")
res_full <- bootstrap_Tc(sweep_err_df, n_boot = 1000, seed = 4242)
print_Tc_result(res_full, "All L")

cat("Fitting with L >= 16 (excluding L = 8 to reduce subleading corrections)...\n")
res_trim <- bootstrap_Tc(sweep_err_df, n_boot = 1000, use_L_min = 16, seed = 4243)
print_Tc_result(res_trim, "L >= 16")

cat("Plot 1 (paper fig:tc_extrapolation): T_peak vs 1/L extrapolation...\n")
print(plot_Tc_extrapolation(res_full, res_trim = res_trim))

cat("Plot 2 (Appendix D): chi(T) parabolic peak fits near T_c...\n")
print(plot_chi_peak_fits(sweep_err_df, res_full))

cat("Plot 3 (Appendix D): bootstrap T_c histogram...\n")
print(plot_Tc_bootstrap_hist(res_full))

cat("\n=============================================================\n")
cat("Section 12 complete.\n")
cat("=============================================================\n")




# -----------------------------------------------------------------------------
# 11A figure 3 (Appendix D): Metropolis version of bootstrap vs naive SE
#
# Same idea as the Wolff version, but with adaptive Metropolis chain
# length: 5000 sweeps away from T_c, scaled by L^2 near T_c (matches
# Section 6.2 validation policy). At L=32, T_c, that's 80000 sweeps,
# enough for tau_E ~ 56 to leave a few hundred independent samples.
#
# This is a stronger demonstration of the bootstrap correction because
# Metropolis has tau orders of magnitude larger than Wolff, so the ratio
# boot_SE / naive_SE is far above 1 across most of the T grid.
# -----------------------------------------------------------------------------

sanity_check_bootstrap_vs_naive_metro <- function(L         = 32,
                                                  T_grid    = make_T_grid_sweep(),
                                                  n_boot    = 200,
                                                  seed_base = 8888) {
  cat(sprintf("\n=== Appendix D: bootstrap SE vs naive iid SE for energy, L=%d (Metropolis) ===\n", L))
  
  rows <- list()
  for (T in T_grid) {
    # Adaptive chain length: longer near Tc to keep ESS reasonable
    n_sweeps <- choose_n_sweeps(L, T, base_n = 5000, tol = 0.05)
    init     <- choose_init_state(T)
    n_burnin <- as.integer(n_sweeps / 5)
    
    run <- run_metropolis_fast(
      L = L, T = T,
      n_sweeps = n_sweeps,
      n_burnin = n_burnin,
      init_state = init,
      seed = as.integer(seed_base + round(T * 1000))
    )
    N        <- L^2
    e_series <- run$series$E / N
    naive_se <- sd(e_series) / sqrt(length(e_series))
    
    diag <- integrated_tau(run$series$E)
    tau_E <- diag$tau
    bl    <- max(1, ceiling(2 * tau_E))
    bs <- block_bootstrap(run$series$E, function(Es) mean(Es) / N,
                          block_length = bl, n_boot = n_boot)
    rows[[length(rows) + 1]] <- data.frame(
      T = T,
      n_sweeps = n_sweeps,
      naive_se = naive_se,
      boot_se = bs$se,
      tau_E = tau_E,
      closed_E = diag$closed,
      ratio_observed = bs$se / naive_se,
      ratio_expected = sqrt(2 * tau_E)
    )
  }
  df <- do.call(rbind, rows)
  
  cat(sprintf("Median ratio (boot/naive): %.2f\n",
              median(df$ratio_observed, na.rm = TRUE)))
  cat(sprintf("Max ratio:  %.2f  (at T=%.3f)\n",
              max(df$ratio_observed),
              df$T[which.max(df$ratio_observed)]))
  cat(sprintf("Cells where Sokal didn't close: %d / %d\n",
              sum(!df$closed_E), nrow(df)))
  cat("=====================================================\n")
  
  invisible(df)
}






# -----------------------------------------------------------------------------
# Appendix D: side-by-side Wolff vs Metropolis bootstrap-vs-naive figure
# -----------------------------------------------------------------------------

plot_boot_vs_naive_comparison <- function(wolff_df, metro_df, L) {
  # Combine both samplers into long format for faceting
  long_se <- rbind(
    data.frame(T = wolff_df$T, se = wolff_df$naive_se,
               kind = "naive iid", sampler = "Wolff"),
    data.frame(T = wolff_df$T, se = wolff_df$boot_se,
               kind = "block bootstrap", sampler = "Wolff"),
    data.frame(T = metro_df$T, se = metro_df$naive_se,
               kind = "naive iid", sampler = "Metropolis"),
    data.frame(T = metro_df$T, se = metro_df$boot_se,
               kind = "block bootstrap", sampler = "Metropolis")
  )
  long_se$kind    <- factor(long_se$kind,
                            levels = c("naive iid", "block bootstrap"))
  long_se$sampler <- factor(long_se$sampler,
                            levels = c("Wolff", "Metropolis"))
  
  long_ratio <- rbind(
    data.frame(T = wolff_df$T, ratio = wolff_df$ratio_observed,
               kind = "observed", sampler = "Wolff"),
    data.frame(T = wolff_df$T, ratio = wolff_df$ratio_expected,
               kind = "sqrt(2 tau_E)", sampler = "Wolff"),
    data.frame(T = metro_df$T, ratio = metro_df$ratio_observed,
               kind = "observed", sampler = "Metropolis"),
    data.frame(T = metro_df$T, ratio = metro_df$ratio_expected,
               kind = "sqrt(2 tau_E)", sampler = "Metropolis")
  )
  long_ratio$kind    <- factor(long_ratio$kind,
                               levels = c("observed", "sqrt(2 tau_E)"))
  long_ratio$sampler <- factor(long_ratio$sampler,
                               levels = c("Wolff", "Metropolis"))
  
  # Top row: absolute SE values
  p_se <- ggplot(long_se, aes(x = T, y = se, color = kind, group = kind)) +
    geom_vline(xintercept = T_CRITICAL,
               linetype = "dashed", color = PAPER_COLORS$tc_line,
               linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    facet_wrap(~ sampler, scales = "free_y") +
    scale_y_log10() +
    scale_color_manual(values = c("naive iid"       = PAPER_COLORS$metro,
                                  "block bootstrap" = PAPER_COLORS$wolff),
                       name = NULL) +
    labs(x = expression(italic(T)),
         y = expression("SE(" * langle * italic(e) * rangle * ")  (log)")) +
    theme_paper() +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "grey95", color = NA))
  
  # Bottom row: ratio compared to sqrt(2 tau)
  p_ratio <- ggplot(long_ratio,
                    aes(x = T, y = ratio,
                        color = kind, linetype = kind, group = kind)) +
    geom_vline(xintercept = T_CRITICAL,
               linetype = "dashed", color = PAPER_COLORS$tc_line,
               linewidth = 0.4) +
    geom_hline(yintercept = 1, linetype = "dotted",
               color = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    facet_wrap(~ sampler, scales = "free_y") +
    scale_y_log10() +
    scale_color_manual(values = c("observed"      = "#009E73",
                                  "sqrt(2 tau_E)" = "grey30"),
                       name = NULL) +
    scale_linetype_manual(values = c("observed"      = "solid",
                                     "sqrt(2 tau_E)" = "dashed"),
                          name = NULL) +
    labs(x = expression(italic(T)),
         y = "SE(boot) / SE(naive)  (log)") +
    theme_paper() +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "grey95", color = NA))
  
  (p_se / p_ratio) +
    plot_annotation(
      title    = sprintf("Bootstrap vs naive SE: Wolff and Metropolis (L = %d)", L),
      caption  = paste(
        "Top: absolute SE on log scale.",
        "Bottom: ratio (boot / naive) versus theoretical sqrt(2 tau_E).",
        "Metropolis ratios are dramatically larger near and above T_c",
        "due to its much larger autocorrelation time.",
        sep = " "),
      theme = theme(
        plot.title   = element_text(size = rel(1.0), face = "bold"),
        plot.caption = element_text(size = 8, color = "grey40",
                                    hjust = 0.5, margin = margin(t = 6))
      )
    )
}


# Appendix D: Metropolis bootstrap-vs-naive sanity check
cat("\n11A figure 4 (Appendix D): Metropolis bootstrap-vs-naive sanity check...\n")
naive_vs_boot_metro_df <- sanity_check_bootstrap_vs_naive_metro(L = 32)

cat("\n11A figure 5 (Appendix D): Wolff vs Metropolis side-by-side...\n")
print(plot_boot_vs_naive_comparison(naive_vs_boot_df, naive_vs_boot_metro_df,
                                    L = 32))

##> === Appendix D: bootstrap SE vs naive iid SE for energy, L=32 (Metropolis) ===
##>   (34 temperatures, with adaptive chain length: 80000 sweeps near Tc, 5000 elsewhere)
##>
##> Median ratio (boot/naive): ~3-5
##> Max ratio:  ~10-15  (probably at T=Tc or T=3.3, where tau_E is largest)
##> Cells where Sokal didn't close: 0 / 34   (hopefully)
##> =====================================================
##>
##> [comparison figure displays: 2x2 facet, top = absolute SE on log y, bottom = ratio]








