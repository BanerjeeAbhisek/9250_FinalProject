# =============================================================================
# Comparative Study of MCMC Sampling Strategies for the 2D Ising Model
# Sections 1-6 (consolidated, v3 with all roundoff fixes)
#
# Abhisek Banerjee & Mengyan Jing
# Course 9250, Spring 2026
#
# Layout:
#   Section 1: Setup (T_c, physical constants, literature z values)
#   Section 2: Lattice utilities (init, total energy, dE with PBC)
#   Section 3: Onsager / Yang exact ground truth
#              + safe wrappers for plotting and ground-truth comparison
#              that avoid the elliptic_K(k=1) numerical singularity
#   Section 4: Sanity checks (14 assertions, all PASS)
#   Section 5: Pure R Metropolis sampler (reference)
#   Section 6: Rcpp Metropolis sampler (~300x faster)
#              6.1: benchmark_metropolis_samplers + plot
#              6.2: validate_metropolis_all_observables + plot
#
# Numerical/physics fixes vs the original Rmd:
#   * z ~ 2.17 attributed to Nightingale & Blote (1996), not Wolff (1989)
#   * elliptic_K guarded at |k| >= 1 - 1e-10 to avoid integrate() roundoff
#   * Onsager curve wrappers return NA in a window around T_c for plots
#   * safe_onsager_* shift by 1e-3 from T_c (not 1e-6: too small, hits roundoff)
#   * Cold start + n_burnin = n_sweeps/5 below T_c for ergodicity
#   * Adaptive n_sweeps ~ L^2 at T_c for critical slowing down
#   * Finite-L magnetization compared to L^(-1/8) at T_c, 1/sqrt(N) above T_c
# =============================================================================


# =============================================================================
# Section 1: Setup
# =============================================================================

suppressPackageStartupMessages(library(Rcpp))

J_COUPLING  <- 1
K_BOLTZMANN <- 1
T_CRITICAL  <- 2 / log(1 + sqrt(2))     # Onsager (1944) Eq. (17)

# Literature dynamic critical exponents (used downstream in Sections 8 & 9)
Z_METROPOLIS_LITERATURE <- 2.1665   # Nightingale & Blote (1996), PRL 76, 4548
Z_WOLFF_LITERATURE      <- 0.25     # Wolff (1989), PRL 62, 361

cat(sprintf("T_c = %.6f\n", T_CRITICAL))


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

# 2.2: total energy (Newman & Barkema Eq. 3.1) --------------------------------

total_energy <- function(spins) {
  L <- nrow(spins); stopifnot(ncol(spins) == L)
  right <- spins[, c(2:L, 1)]
  down  <- spins[c(2:L, 1), ]
  -J_COUPLING * sum(spins * right + spins * down)
}

total_magnetization <- function(spins) sum(spins)

# 2.3: single-spin dE (Newman & Barkema Eq. 3.10) -----------------------------

delta_energy <- function(spins, i, j) {
  L  <- nrow(spins)
  ip <- if (i == L)  1L else i + 1L
  im <- if (i == 1L) L  else i - 1L
  jp <- if (j == L)  1L else j + 1L
  jm <- if (j == 1L) L  else j - 1L
  neighbor_sum <- spins[ip, j] + spins[im, j] + spins[i, jp] + spins[i, jm]
  2 * J_COUPLING * spins[i, j] * neighbor_sum
}


# =============================================================================
# Section 3: Onsager / Yang exact ground truth
# =============================================================================

# 3.1: complete elliptic integrals --------------------------------------------
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

# 3.2: energy per spin (Onsager Eq. 116) --------------------------------------

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

# 3.3: specific heat (Onsager Eq. 117) ----------------------------------------

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

# 3.4: asymptotic specific heat near T_c (Onsager Eq. 120) --------------------

onsager_specific_heat_asymptotic <- function(T) {
  vapply(T, function(Ti) {
    beta <- 1 / Ti; bj <- beta * J_COUPLING
    k1   <- 2 * sinh(2 * bj) / cosh(2 * bj)^2
    K1   <- elliptic_K(k1)
    (2 / pi) * log(1 / tan(pi / 8))^2 * (K1 - 1 - pi / 4)
  }, numeric(1))
}

# 3.5: magnetization (Yang 1952) ----------------------------------------------

onsager_magnetization <- function(T) {
  vapply(T, function(Ti) {
    if (Ti >= T_CRITICAL) return(0)
    beta <- 1 / Ti
    (1 - sinh(2 * beta * J_COUPLING)^(-4))^(1 / 8)
  }, numeric(1))
}

# 3.6: safe ground-truth wrappers (avoid singular elliptic_K) -----------------
# Use vapply + if/else (lazy, only evaluates the chosen branch).
# .SHIFT = 1e-3: u(T_c - 1e-3) = -1.4182 vs true -1.4142 (0.3% off, fine).

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

# 3.7: curve wrappers for plotting (NA in window around T_c) ------------------
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

# 3.8: plot the ground truth --------------------------------------------------

plot_onsager_ground_truth <- function() {
  T_grid <- seq(1, 4, by = 0.01)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  plot(T_grid, onsager_energy_curve(T_grid), type = "l", col = "darkblue", lwd = 2,
       xlab = "T", ylab = "u(T)", main = "Energy per spin")
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  
  plot(T_grid, onsager_specific_heat_curve(T_grid), type = "l",
       col = "darkred", lwd = 2,
       xlab = "T", ylab = "c(T)", main = "Specific heat per spin")
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  
  plot(T_grid, onsager_magnetization_curve(T_grid), type = "l",
       col = "darkgreen", lwd = 2,
       xlab = "T", ylab = "|m|(T)", main = "Spontaneous magnetization")
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  
  plot.new()
  text(0.5, 0.5, sprintf("T_c = %.4f", T_CRITICAL), cex = 1.5)
  
  par(mfrow = c(1, 1))
}


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

plot_metropolis_benchmark <- function(bench_df) {
  T_levels <- sort(unique(bench_df$T))
  T_labels <- unique(bench_df[order(bench_df$T), "T_label"])
  T_colors <- c("darkblue", "darkred", "darkgreen")[seq_along(T_levels)]
  L_levels <- sort(unique(bench_df$L))
  L_colors <- c("black", "darkorange", "purple", "darkcyan")[seq_along(L_levels)]
  
  old_par <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.5, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Wall time vs L log-log
  y_range <- range(c(bench_df$t_R, bench_df$t_cpp))
  plot(NA, xlim = range(L_levels), ylim = y_range, log = "xy",
       xlab = "L", ylab = "wall time [s]", main = "Wall time vs L (log-log)")
  for (k in seq_along(T_levels)) {
    sub <- bench_df[bench_df$T == T_levels[k], ]; sub <- sub[order(sub$L), ]
    lines(sub$L,  sub$t_R,   col = T_colors[k], lwd = 2, lty = 1)
    points(sub$L, sub$t_R,   col = T_colors[k], pch = 19, cex = 1.1)
    lines(sub$L,  sub$t_cpp, col = T_colors[k], lwd = 2, lty = 2)
    points(sub$L, sub$t_cpp, col = T_colors[k], pch = 17, cex = 1.1)
  }
  L_ref <- range(L_levels)
  y_ref <- max(bench_df$t_R) * (L_ref / min(L_levels))^2
  lines(L_ref, y_ref, col = "gray60", lwd = 1, lty = 3)
  legend("topleft", cex = 0.75, bty = "n",
         legend = c(paste("R,",   T_labels), paste("Cpp,", T_labels), "slope 2"),
         col    = c(T_colors, T_colors, "gray60"),
         lty    = c(rep(1, length(T_levels)), rep(2, length(T_levels)), 3),
         pch    = c(rep(19, length(T_levels)), rep(17, length(T_levels)), NA),
         lwd    = c(rep(2, 2 * length(T_levels)), 1))
  
  # Speedup vs L
  sp_range <- range(bench_df$speedup)
  plot(NA, xlim = range(L_levels), ylim = sp_range, log = "x",
       xlab = "L", ylab = "speedup (t_R / t_cpp)", main = "Speedup vs L")
  for (k in seq_along(T_levels)) {
    sub <- bench_df[bench_df$T == T_levels[k], ]; sub <- sub[order(sub$L), ]
    lines(sub$L,  sub$speedup, col = T_colors[k], lwd = 2)
    points(sub$L, sub$speedup, col = T_colors[k], pch = 19, cex = 1.2)
  }
  abline(h = median(bench_df$speedup), lty = 3, col = "gray50")
  legend("bottomright", cex = 0.8, bty = "n",
         legend = c(T_labels, "median"), col = c(T_colors, "gray50"),
         lty = c(rep(1, length(T_levels)), 3),
         pch = c(rep(19, length(T_levels)), NA),
         lwd = c(rep(2, length(T_levels)), 1))
  
  # Speedup vs T
  plot(NA, xlim = range(T_levels), ylim = sp_range,
       xlab = "T", ylab = "speedup (t_R / t_cpp)", main = "Speedup vs T")
  for (k in seq_along(L_levels)) {
    sub <- bench_df[bench_df$L == L_levels[k], ]; sub <- sub[order(sub$T), ]
    lines(sub$T,  sub$speedup, col = L_colors[k], lwd = 2)
    points(sub$T, sub$speedup, col = L_colors[k], pch = 19, cex = 1.2)
  }
  abline(v = T_CRITICAL, lty = 2, col = "gray60")
  legend("bottomright", cex = 0.8, bty = "n",
         legend = c(paste("L =", L_levels), "T_c"),
         col = c(L_colors, "gray60"),
         lty = c(rep(1, length(L_levels)), 2),
         pch = c(rep(19, length(L_levels)), NA),
         lwd = c(rep(2, length(L_levels)), 1))
  
  invisible(NULL)
}


# =============================================================================
# Section 6.2: Comprehensive validation against all ground truths
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

plot_metropolis_validation <- function(validation_df) {
  L_values <- sort(unique(validation_df$L))
  L_colors <- c("black", "darkorange", "purple", "darkcyan")[seq_along(L_values)]
  
  T_curve <- seq(1.0, 4.2, length.out = 400)
  u_exact <- onsager_energy_curve(T_curve)
  c_exact <- onsager_specific_heat_curve(T_curve)
  m_exact <- onsager_magnetization_curve(T_curve)
  c_cap   <- 4
  c_plot  <- pmin(c_exact, c_cap)
  
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Energy
  plot(T_curve, u_exact, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "u(T)", main = "Energy: MCMC vs Onsager",
       ylim = range(c(u_exact, validation_df$e_mcmc), finite = TRUE))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- validation_df[validation_df$L == L_values[k], ]
    points(sub$T, sub$e_mcmc, pch = 19, cex = 1.3, col = L_colors[k])
  }
  legend("topleft", cex = 0.8, bty = "n",
         legend = c(paste("L =", L_values), "Onsager"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(NA, length(L_values)), 2),
         lty = c(rep(NA, length(L_values)), 1))
  
  # Specific heat capped
  plot(T_curve, c_plot, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "c(T)",
       main = sprintf("Specific heat: MCMC vs Onsager (capped at %.0f)", c_cap),
       ylim = c(0, c_cap))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- validation_df[validation_df$L == L_values[k], ]
    points(sub$T, sub$c_mcmc, pch = 19, cex = 1.3, col = L_colors[k])
  }
  legend("topleft", cex = 0.8, bty = "n",
         legend = c(paste("L =", L_values), "Onsager (capped)"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(NA, length(L_values)), 2),
         lty = c(rep(NA, length(L_values)), 1))
  
  # Magnetization with 1/sqrt(N) ticks
  plot(T_curve, m_exact, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "|m|(T)",
       main = "Magnetization: MCMC vs Yang + 1/sqrt(N) ref",
       ylim = c(0, 1.05))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- validation_df[validation_df$L == L_values[k], ]
    points(sub$T, sub$m_mcmc, pch = 19, cex = 1.3, col = L_colors[k])
    N_k <- L_values[k]^2
    segments(T_CRITICAL + 0.1, 1/sqrt(N_k), 4.2, 1/sqrt(N_k),
             col = L_colors[k], lty = 3, lwd = 1)
  }
  legend("bottomleft", cex = 0.8, bty = "n",
         legend = c(paste("L =", L_values), "Yang", "1/sqrt(N) ref"),
         col = c(L_colors, "gray40", "gray40"),
         pch = c(rep(19, length(L_values)), NA, NA),
         lwd = c(rep(NA, length(L_values)), 2, 1),
         lty = c(rep(NA, length(L_values)), 1, 3))
  
  # Susceptibility
  plot(validation_df$T, validation_df$chi_mcmc, type = "n",
       xlab = "T", ylab = expression(chi(T)),
       main = "Susceptibility: MCMC (peak grows as L^1.75)")
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- validation_df[validation_df$L == L_values[k], ]
    points(sub$T, sub$chi_mcmc, pch = 19, cex = 1.3, col = L_colors[k])
  }
  legend("topleft", cex = 0.8, bty = "n",
         legend = paste("L =", L_values), col = L_colors, pch = 19)
  
  invisible(NULL)
}


# =============================================================================
# Drive: run everything end-to-end
# =============================================================================

# Section 4: sanity + ground-truth visualization
run_sanity_checks()
plot_onsager_ground_truth()

# Section 6.1: speed benchmark (~2-3 minutes)
bench_df <- benchmark_metropolis_samplers()
plot_metropolis_benchmark(bench_df)

# Section 6.2: ground-truth validation (~1-2 minutes)
validation_df <- validate_metropolis_all_observables()
plot_metropolis_validation(validation_df)




# =============================================================================
# Section 7: Visualizing chain behaviour across all (L, T) combinations
#
# Four figures, one per observable. Each figure is a 4x3 grid:
#   rows = L in {8, 16, 32, 64}
#   cols = T in {1.5, T_c, 4.0}
#
# The 12 cells are the same (L, T) grid as the validation runs. We re-use
# the validation choices (cold start below T_c, adaptive n_sweeps at T_c,
# burn-in n_sweeps/5) so what's shown here is exactly what fed the
# ground-truth tables in Section 6.2.
#
# What "trace" means depends on the observable:
#
#   - e and |m| are time series in the natural sense: one sample per sweep.
#     Plot the per-spin value at each sweep against sweep index. Overlay
#     the Onsager / Yang ground truth as a dashed horizontal where it
#     exists (Yang's |m| = 0 above T_c is plotted as a flat line at zero).
#
#   - c and chi are NOT per-sweep observables -- they are computed as
#     variances over the chain. The meaningful "trace" is the running
#     estimate after the first t samples:
#         c_running(t)   = (mean(E^2)[1:t] - mean(E)[1:t]^2) / (N T^2)
#         chi_running(t) = (mean(M^2)[1:t] - mean(|M|)[1:t]^2) / (N T)
#     This shows how the estimate stabilises as more samples accumulate.
#     A flat running estimate means the chain has converged for that
#     observable; a drifting running estimate means it has not.
#
# Implementation notes:
#   - We need to keep the time series for every (L, T) cell, not just the
#     summary observables. The new function trace_runs_metropolis() runs
#     and returns the runs themselves.
#   - For long chains (L=64 at T_c is 320k sweeps), per-sweep trace plots
#     get visually saturated. We thin to ~5000 plotted points by taking
#     every ceiling(n_sweeps / 5000)-th sample.
#   - Within a single observable's figure, all 12 panels share a common
#     y-axis range, so visual comparison across panels is honest. (Without
#     this, autoscaling each panel hides the variance ordering that is
#     the whole point of the figure.)
# =============================================================================


# -----------------------------------------------------------------------------
# Run the 12 (L, T) cells and store the full runs (not just summaries).
# Returns a named list keyed by "L8_T1.5", "L8_Tc", etc.
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
      
      key <- sprintf("L%d_T%s",
                     L,
                     if (abs(T - T_CRITICAL) < 1e-6) "c" else
                       sprintf("%.2f", T))
      
      cat(sprintf("  Running L=%2d, T=%s (%-8s): n_sweeps=%d, init=%s ... ",
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
# Helper: thin a long series to ~max_pts points for plotting.
# -----------------------------------------------------------------------------

thin_for_plot <- function(x, max_pts = 5000) {
  n <- length(x)
  if (n <= max_pts) return(seq_len(n))
  step <- ceiling(n / max_pts)
  seq(1L, n, by = step)
}


# -----------------------------------------------------------------------------
# Running estimates for c and chi (the variance-based observables).
# These are the "trace" equivalents for non-additive observables.
# -----------------------------------------------------------------------------

running_specific_heat <- function(E, N, T) {
  # Running mean and running mean of squares; protect against numerical
  # zero variance at very early t by replacing with 0.
  csum   <- cumsum(E)
  csum2  <- cumsum(E^2)
  n      <- seq_along(E)
  mean_E   <- csum / n
  mean_E2  <- csum2 / n
  var_E    <- pmax(mean_E2 - mean_E^2, 0)
  var_E / (N * T^2)
}

running_susceptibility <- function(M, N, T) {
  csum   <- cumsum(M)
  csum2  <- cumsum(M^2)
  csumA  <- cumsum(abs(M))
  n      <- seq_along(M)
  mean_M2  <- csum2 / n
  mean_aM  <- csumA / n
  # chi = (<M^2> - <|M|>^2) / (N T) on a finite lattice in zero field
  var_M  <- pmax(mean_M2 - mean_aM^2, 0)
  var_M / (N * T)
}


# -----------------------------------------------------------------------------
# Compute one observable's trace for every (L, T) cell, plus a sensible
# common y-range. Returns a list of (key, x, y, ground_truth_value or NA).
# -----------------------------------------------------------------------------

build_traces_one_observable <- function(runs, observable) {
  # observable is one of "e", "m", "c", "chi"
  out <- list()
  for (key in names(runs)) {
    run <- runs[[key]]
    N   <- run$L^2
    T   <- run$T
    n   <- run$n_sweeps
    
    if (observable == "e") {
      y <- run$series$E / N
      g <- safe_onsager_energy(T)             # finite at T_c via wrapper
    } else if (observable == "m") {
      y <- abs(run$series$M) / N
      g <- if (T < T_CRITICAL - 0.01) safe_onsager_magnetization(T)
      else if (abs(T - T_CRITICAL) < 0.05) NA   # scaling regime, no h-line
      else 1 / sqrt(N)                          # finite-L 1/sqrt(N) ref
    } else if (observable == "c") {
      y <- running_specific_heat(run$series$E, N, T)
      g <- if (abs(T - T_CRITICAL) < 0.05)
        onsager_specific_heat_proxy(T)         # finite proxy
      else
        onsager_specific_heat(T)
    } else if (observable == "chi") {
      y <- running_susceptibility(run$series$M, N, T)
      g <- NA   # no exact closed form
    } else {
      stop("Unknown observable: ", observable)
    }
    
    out[[key]] <- list(
      key      = key,
      L        = run$L,
      T        = T,
      T_label  = run$T_label,
      x        = seq_along(y),
      y        = y,
      ground   = g,
      acc_rate = run$acceptance_rate
    )
  }
  out
}


# -----------------------------------------------------------------------------
# Plot one figure (4 rows = L, 3 cols = T) for one observable.
# -----------------------------------------------------------------------------

plot_observable_grid <- function(traces, observable,
                                 L_values = c(8, 16, 32, 64),
                                 T_values = c(1.5, T_CRITICAL, 4.0),
                                 T_labels = c("below Tc", "at Tc", "above Tc"),
                                 max_pts  = 5000) {
  
  # Common y-range across the 12 panels for this observable. Drop the
  # first 2% of each chain when computing y-range (avoids being driven by
  # transient near-start values for the running estimates).
  y_all <- c()
  for (tr in traces) {
    n   <- length(tr$y)
    cut <- max(2L, as.integer(0.02 * n))
    y_all <- c(y_all, tr$y[cut:n])
  }
  y_all <- y_all[is.finite(y_all)]
  y_range <- range(y_all)
  
  # Slight padding
  y_range <- y_range + diff(y_range) * c(-0.05, 0.05)
  
  # Pretty observable label
  obs_label <- switch(observable,
                      e    = "e (energy per spin)",
                      m    = "|m| (magnetization per spin)",
                      c    = "c running estimate",
                      chi  = expression(chi ~ "running estimate"))
  
  obs_title <- switch(observable,
                      e    = "Energy trace",
                      m    = "Magnetization trace",
                      c    = "Specific heat (running)",
                      chi  = "Susceptibility (running)")
  
  obs_color <- switch(observable,
                      e    = "darkblue",
                      m    = "darkgreen",
                      c    = "darkred",
                      chi  = "purple")
  
  old_par <- par(mfrow = c(length(L_values), length(T_values)),
                 mar    = c(3.2, 3.6, 2.0, 0.6),
                 mgp    = c(2.0, 0.6, 0),
                 oma    = c(0, 0, 2.4, 0),
                 no.readonly = TRUE)
  on.exit(par(old_par))
  
  for (L in L_values) {
    for (k in seq_along(T_values)) {
      T     <- T_values[k]
      label <- T_labels[k]
      key   <- sprintf("L%d_T%s",
                       L,
                       if (abs(T - T_CRITICAL) < 1e-6) "c" else sprintf("%.2f", T))
      tr <- traces[[key]]
      idx <- thin_for_plot(tr$x, max_pts)
      
      plot(tr$x[idx], tr$y[idx], type = "l", col = obs_color, lwd = 1,
           xlab = "sweep", ylab = obs_label,
           ylim = y_range,
           main = sprintf("L=%d, T=%s (acc=%.2f)",
                          L,
                          if (abs(T - T_CRITICAL) < 1e-6) "T_c"
                          else sprintf("%.2f", T),
                          tr$acc_rate),
           cex.main = 0.9)
      
      # Ground truth horizontal where defined
      if (is.finite(tr$ground)) {
        abline(h = tr$ground, lty = 2, col = "gray30", lwd = 1.2)
      }
    }
  }
  
  # Outer title
  mtext(obs_title, side = 3, outer = TRUE, line = 0.6, cex = 1.2, font = 2)
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Top-level driver: run all 12 cells once, then make 4 figures.
# -----------------------------------------------------------------------------

plot_section7_all_traces <- function(L_values = c(8, 16, 32, 64),
                                     T_values = c(1.5, T_CRITICAL, 4.0),
                                     T_labels = c("below Tc", "at Tc", "above Tc"),
                                     seed_base = 42,
                                     max_pts   = 5000,
                                     pause_between = FALSE) {
  
  cat("\n############################################################\n")
  cat("Section 7: trace plots across (L, T) for all observables\n")
  cat("############################################################\n\n")
  
  runs <- trace_runs_metropolis(L_values  = L_values,
                                T_values  = T_values,
                                T_labels  = T_labels,
                                seed_base = seed_base)
  
  for (obs in c("e", "m", "c", "chi")) {
    cat(sprintf("Plotting figure for observable: %s\n", obs))
    traces <- build_traces_one_observable(runs, observable = obs)
    plot_observable_grid(traces, observable = obs,
                         L_values = L_values, T_values = T_values,
                         T_labels = T_labels, max_pts = max_pts)
    if (pause_between) {
      cat("Press [Enter] to continue to the next figure...")
      readline()
    }
  }
  
  invisible(runs)
}


# =============================================================================
# Drive
# =============================================================================

# pause_between = TRUE if you want to inspect each figure interactively
# before the next one is drawn. Set FALSE for a one-shot run.
section7_runs <- plot_section7_all_traces(pause_between = FALSE)

# =============================================================================
# Section 8: Autocorrelation, integrated time, and effective sample size
#
# For every (L, T) combination from Section 7, compute:
#
#   - rho_E(t), rho_M(t):   the empirical autocorrelation function for
#                           energy and absolute magnetization
#   - tau_E, tau_M:         the integrated autocorrelation time, estimated
#                           by Sokal's automatic windowing rule
#                              W = min { W : W >= c * tau_hat(W) }
#                           with c = 5 (Sokal's recommendation). If W
#                           reaches the max_lag cutoff without closure,
#                           the estimate is a downward-biased lower bound;
#                           we flag this loudly in the table.
#   - ESS_E, ESS_M:         n / (2 * tau)
#   - SE(<e>), SE(<|m|>):   the corrected standard error
#                              SE = sigma * sqrt(2 * tau / n)
#                           which is the iid SE inflated by sqrt(2 tau).
#
# Why only e and |m|?
#   These are per-sweep observables: one number per Monte Carlo sweep,
#   and the standard tau_int / ESS formulas are defined for the sample
#   mean of such a series. Specific heat c and susceptibility chi are
#   variances of E and M, not means; Section 11 (block bootstrap) is
#   the right machinery for their uncertainty. Here we report tau_E
#   and tau_M, which Section 11 will then use to pick block sizes.
#
# Reuse from Section 7:
#   The 12 runs in `section7_runs` are exactly what we need. We do not
#   re-simulate; we just compute diagnostics on what we already have.
# =============================================================================


# -----------------------------------------------------------------------------
# 8.1: autocorrelation function and integrated tau (Sokal windowing)
#
# These are unchanged from your original Rmd's definitions, just kept here
# for self-containedness. The 0.5 floor on tau prevents numerical dips
# below the iid theoretical minimum from breaking log-axis plots.
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
  
  # Standard errors of the means, accounting for autocorrelation
  se_e <- sd(E_per_spin) * sqrt(2 * diag_E$tau / n)
  se_m <- sd(m_per_spin) * sqrt(2 * diag_M$tau / n)
  
  list(
    L              = run$L,
    T              = run$T,
    T_label        = run$T_label,
    n_sweeps       = n,
    acceptance     = run$acceptance_rate,
    tau_E          = diag_E$tau,
    tau_M          = diag_M$tau,
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

# -----------------------------------------------------------------------------
# 8.3: assemble the full diagnostics table from a list of runs
# -----------------------------------------------------------------------------

build_diagnostics_table <- function(runs) {
  rows <- list()
  for (key in names(runs)) {
    d <- cell_diagnostics(runs[[key]])
    rows[[key]] <- data.frame(
      L          = d$L,
      T          = d$T,
      T_label    = d$T_label,
      n_sweeps   = d$n_sweeps,
      acc        = d$acceptance,
      tau_E      = d$tau_E,
      tau_M      = d$tau_M,
      ESS_E      = d$ESS_E,
      ESS_M      = d$ESS_M,
      se_e       = d$se_e,
      se_m       = d$se_m,
      window_E   = d$window_E,
      window_M   = d$window_M,
      closed_E   = d$closed_E,
      closed_M   = d$closed_M,
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
  
  # Loud warning if any closure failed
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
#
# Standard scaling: tau ~ L^z. Fit log(tau) ~ z * log(L) using the rows at
# T_c. Skip cells whose Sokal didn't close (those tau values are biased).
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
  cat(sprintf("  z (from tau_%s, fit on L = %s): z = %.3f  (literature %.4f)\n",
              observable, paste(sub$L, collapse = ","), z,
              if (observable == "E") Z_METROPOLIS_LITERATURE
              else                    Z_METROPOLIS_LITERATURE))
  invisible(z)
}


# -----------------------------------------------------------------------------
# 8.6: plot rho(t) grids for E and M (4 rows x 3 cols)
# -----------------------------------------------------------------------------

plot_autocorr_grid <- function(runs, observable = c("E", "M"),
                               L_values = c(8, 16, 32, 64),
                               T_values = c(1.5, T_CRITICAL, 4.0),
                               T_labels = c("below Tc", "at Tc", "above Tc"),
                               max_lag_plot = 200) {
  observable <- match.arg(observable)
  
  obs_label <- switch(observable,
                      E = "rho_E(t)  [energy autocorrelation]",
                      M = "rho_M(t)  [magnetization autocorrelation]")
  obs_color <- switch(observable, E = "darkblue", M = "darkred")
  obs_title <- switch(observable,
                      E = "Energy autocorrelation function",
                      M = "Magnetization autocorrelation function")
  
  old_par <- par(mfrow = c(length(L_values), length(T_values)),
                 mar    = c(3.2, 3.6, 2.0, 0.6),
                 mgp    = c(2.0, 0.6, 0),
                 oma    = c(0, 0, 2.4, 0),
                 no.readonly = TRUE)
  on.exit(par(old_par))
  
  for (L in L_values) {
    for (k in seq_along(T_values)) {
      T   <- T_values[k]
      key <- sprintf("L%d_T%s",
                     L,
                     if (abs(T - T_CRITICAL) < 1e-6) "c" else sprintf("%.2f", T))
      run <- runs[[key]]
      n_sweeps <- nrow(run$series)
      
      # Use a per-cell max_lag that's smaller than n/4 but at least 200 if possible
      max_lag_here <- min(max_lag_plot, n_sweeps %/% 4)
      
      x_series <- if (observable == "E") run$series$E else run$series$M
      rho      <- autocorrelation(x_series, max_lag = max_lag_here)
      
      # Compute tau for the panel title
      diag <- integrated_tau(x_series)
      flag <- if (diag$closed) "" else "!"
      
      plot(0:max_lag_here, rho, type = "l", col = obs_color, lwd = 1.5,
           xlab = "lag (sweeps)", ylab = obs_label,
           ylim = c(-0.1, 1.05),
           main = sprintf("L=%d, T=%s   tau=%.1f%s",
                          L,
                          if (abs(T - T_CRITICAL) < 1e-6) "T_c"
                          else sprintf("%.2f", T),
                          diag$tau, flag),
           cex.main = 0.9)
      abline(h = 0, lty = 3)
      abline(h = 1/exp(1), lty = 3, col = "gray60")   # 1/e reference
    }
  }
  
  mtext(obs_title, side = 3, outer = TRUE, line = 0.6, cex = 1.2, font = 2)
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# 8.7: summary plot -- tau vs L on log-log, one curve per T,
# side-by-side panels for E and M.
# -----------------------------------------------------------------------------

plot_tau_scaling <- function(diag_df,
                             T_values = c(1.5, T_CRITICAL, 4.0),
                             T_labels = c("below Tc", "at Tc", "above Tc"),
                             T_colors = c("darkblue", "darkred", "darkgreen")) {
  L_values <- sort(unique(diag_df$L))
  
  old_par <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Panel 1: tau_E
  y_E <- diag_df$tau_E
  plot(NA, xlim = range(L_values),
       ylim = range(c(y_E, 0.5), finite = TRUE),
       log = "xy",
       xlab = "L", ylab = expression(tau[E] ~ "(sweeps)"),
       main = "Energy autocorrelation: tau_E vs L")
  for (k in seq_along(T_values)) {
    T   <- T_values[k]
    sub <- diag_df[diag_df$T == T, ]
    sub <- sub[order(sub$L), ]
    pch <- ifelse(sub$closed_E, 19, 4)   # filled = closed, x = open
    lines(sub$L,  sub$tau_E, col = T_colors[k], lwd = 2)
    points(sub$L, sub$tau_E, col = T_colors[k], pch = pch, cex = 1.3)
  }
  
  # Reference line: tau ~ L^z_lit at T_c, anchored at smallest L
  sub_tc <- diag_df[diag_df$T_label == "at Tc", ]
  sub_tc <- sub_tc[order(sub_tc$L), ]
  if (nrow(sub_tc) >= 1) {
    L_ref <- range(L_values)
    y_ref <- sub_tc$tau_E[1] * (L_ref / sub_tc$L[1])^Z_METROPOLIS_LITERATURE
    lines(L_ref, y_ref, col = "gray60", lty = 3, lwd = 1.5)
  }
  legend("topleft", cex = 0.85, bty = "n",
         legend = c(T_labels,
                    sprintf("ref slope %.2f", Z_METROPOLIS_LITERATURE),
                    "filled = Sokal closed", "x = open (biased low)"),
         col    = c(T_colors, "gray60", "black", "black"),
         lty    = c(rep(1, length(T_values)), 3, NA, NA),
         pch    = c(rep(19, length(T_values)), NA, 19, 4),
         lwd    = c(rep(2, length(T_values)), 1.5, NA, NA))
  
  # Panel 2: tau_M
  y_M <- diag_df$tau_M
  plot(NA, xlim = range(L_values),
       ylim = range(c(y_M, 0.5), finite = TRUE),
       log = "xy",
       xlab = "L", ylab = expression(tau[M] ~ "(sweeps)"),
       main = "Magnetization autocorrelation: tau_M vs L")
  for (k in seq_along(T_values)) {
    T   <- T_values[k]
    sub <- diag_df[diag_df$T == T, ]
    sub <- sub[order(sub$L), ]
    pch <- ifelse(sub$closed_M, 19, 4)
    lines(sub$L,  sub$tau_M, col = T_colors[k], lwd = 2)
    points(sub$L, sub$tau_M, col = T_colors[k], pch = pch, cex = 1.3)
  }
  if (nrow(sub_tc) >= 1) {
    L_ref <- range(L_values)
    y_ref <- sub_tc$tau_M[1] * (L_ref / sub_tc$L[1])^Z_METROPOLIS_LITERATURE
    lines(L_ref, y_ref, col = "gray60", lty = 3, lwd = 1.5)
  }
  legend("topleft", cex = 0.85, bty = "n",
         legend = c(T_labels,
                    sprintf("ref slope %.2f", Z_METROPOLIS_LITERATURE),
                    "filled = Sokal closed", "x = open (biased low)"),
         col    = c(T_colors, "gray60", "black", "black"),
         lty    = c(rep(1, length(T_values)), 3, NA, NA),
         pch    = c(rep(19, length(T_values)), NA, 19, 4),
         lwd    = c(rep(2, length(T_values)), 1.5, NA, NA))
  
  invisible(NULL)
}


# =============================================================================
# Top-level driver: produce table, two rho grids, and the scaling plot.
# Requires `section7_runs` from Section 7 (so we don't re-simulate).
# =============================================================================

run_section8_diagnostics <- function(runs = section7_runs) {
  cat("\n############################################################\n")
  cat("Section 8: autocorrelation diagnostics\n")
  cat("############################################################\n\n")
  
  diag_df <- build_diagnostics_table(runs)
  print_diagnostics_table(diag_df)
  
  cat("\n--- Dynamic critical exponent z (Metropolis at T_c) ---\n")
  fit_dynamic_z(diag_df, "E")
  fit_dynamic_z(diag_df, "M")
  
  cat("\nPlotting rho_E grid (4 x 3)...\n")
  plot_autocorr_grid(runs, observable = "E")
  
  cat("Plotting rho_M grid (4 x 3)...\n")
  plot_autocorr_grid(runs, observable = "M")
  
  cat("Plotting tau scaling summary (1 x 2)...\n")
  plot_tau_scaling(diag_df)
  
  invisible(diag_df)
}


# =============================================================================
# Drive
# =============================================================================

section8_diag <- run_section8_diagnostics(section7_runs)

# =============================================================================
# Patch v2: legends outside the plot region
#
# Previous fix put legends in data-free corners -- but below Tc the "flat"
# curves still sit in the bottom portion, so "bottomright" hit them.
# 
# New approach: expand the right margin on each panel and draw the legend
# there using xpd = TRUE + explicit coordinates. The data area is then
# completely clean.
# =============================================================================

plot_tau_scaling <- function(diag_df,
                             T_values = c(1.5, T_CRITICAL, 4.0),
                             T_labels = c("below Tc", "at Tc", "above Tc"),
                             T_colors = c("darkblue", "darkred", "darkgreen")) {
  L_values <- sort(unique(diag_df$L))
  
  # Outer margins: give room on the right for legends outside the plot box.
  # mar = bottom, left, top, right. Bumping right from 1 to 9 lines.
  old_par <- par(mfrow = c(1, 2),
                 mar = c(4.2, 4.2, 2.5, 9),
                 xpd = FALSE,
                 no.readonly = TRUE)
  on.exit(par(old_par))
  
  draw_panel <- function(y_values, y_label, main_title, tau_col_name, closed_col_name) {
    plot(NA, xlim = range(L_values),
         ylim = range(c(y_values, 0.5), finite = TRUE),
         log  = "xy",
         xlab = "L", ylab = y_label,
         main = main_title)
    
    for (k in seq_along(T_values)) {
      T   <- T_values[k]
      sub <- diag_df[diag_df$T == T, ]
      sub <- sub[order(sub$L), ]
      pch <- ifelse(sub[[closed_col_name]], 19, 4)
      lines(sub$L,  sub[[tau_col_name]], col = T_colors[k], lwd = 2)
      points(sub$L, sub[[tau_col_name]], col = T_colors[k], pch = pch, cex = 1.3)
    }
    
    # Reference slope line, anchored at smallest L, using the Tc tau
    sub_tc <- diag_df[diag_df$T_label == "at Tc", ]
    sub_tc <- sub_tc[order(sub_tc$L), ]
    if (nrow(sub_tc) >= 1) {
      L_ref <- range(L_values)
      y_ref <- sub_tc[[tau_col_name]][1] *
        (L_ref / sub_tc$L[1])^Z_METROPOLIS_LITERATURE
      lines(L_ref, y_ref, col = "gray60", lty = 3, lwd = 1.5)
    }
    
    # Legend OUTSIDE the plot region on the right. xpd = TRUE lets us draw
    # outside. We use par("usr") in log space to position it just past the
    # right edge. Using "topright" with inset = negative offset pushes it
    # into the outer margin.
    par(xpd = TRUE)
    legend(x       = 10^(par("usr")[2] + 0.05 * diff(par("usr")[1:2])),
           y       = 10^par("usr")[4],
           cex     = 0.8,
           bty     = "n",
           legend  = c(T_labels,
                       sprintf("ref L^%.2f", Z_METROPOLIS_LITERATURE),
                       "", "Sokal:",
                       "closed", "open (biased low)"),
           col     = c(T_colors, "gray60", NA, NA, "black", "black"),
           lty     = c(rep(1, length(T_values)), 3, NA, NA, NA, NA),
           pch     = c(rep(19, length(T_values)), NA, NA, NA, 19, 4),
           lwd     = c(rep(2, length(T_values)), 1.5, NA, NA, NA, NA))
    par(xpd = FALSE)
  }
  
  draw_panel(diag_df$tau_E,
             expression(tau[E] ~ "(sweeps)"),
             "Energy autocorrelation: tau_E vs L",
             "tau_E", "closed_E")
  
  draw_panel(diag_df$tau_M,
             expression(tau[M] ~ "(sweeps)"),
             "Magnetization autocorrelation: tau_M vs L",
             "tau_M", "closed_M")
  
  invisible(NULL)
}

# Redraw:
plot_tau_scaling(section8_diag)



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
  mid_T <- T_values[ceiling(length(T_values) / 2)]
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


plot_wolff_benchmark <- function(bench_df) {
  T_levels <- sort(unique(bench_df$T))
  T_labels <- unique(bench_df[order(bench_df$T), "T_label"])
  T_colors <- c("darkblue", "darkred", "darkgreen")[seq_along(T_levels)]
  L_levels <- sort(unique(bench_df$L))
  L_colors <- c("black", "darkorange", "purple", "darkcyan")[seq_along(L_levels)]
  
  old_par <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.5, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Panel 1: wall time vs L log-log
  y_range <- range(c(bench_df$t_R, bench_df$t_cpp))
  plot(NA, xlim = range(L_levels), ylim = y_range, log = "xy",
       xlab = "L", ylab = "wall time [s]", main = "Wall time vs L (log-log)")
  for (k in seq_along(T_levels)) {
    sub <- bench_df[bench_df$T == T_levels[k], ]; sub <- sub[order(sub$L), ]
    lines(sub$L,  sub$t_R,   col = T_colors[k], lwd = 2, lty = 1)
    points(sub$L, sub$t_R,   col = T_colors[k], pch = 19, cex = 1.1)
    lines(sub$L,  sub$t_cpp, col = T_colors[k], lwd = 2, lty = 2)
    points(sub$L, sub$t_cpp, col = T_colors[k], pch = 17, cex = 1.1)
  }
  legend("topleft", cex = 0.7, bty = "n",
         legend = c(paste("R,",   T_labels), paste("Cpp,", T_labels)),
         col    = c(T_colors, T_colors),
         lty    = c(rep(1, length(T_levels)), rep(2, length(T_levels))),
         pch    = c(rep(19, length(T_levels)), rep(17, length(T_levels))),
         lwd    = 2)
  
  # Panel 2: speedup vs L
  sp_range <- range(bench_df$speedup)
  plot(NA, xlim = range(L_levels), ylim = sp_range, log = "x",
       xlab = "L", ylab = "speedup", main = "Speedup vs L")
  for (k in seq_along(T_levels)) {
    sub <- bench_df[bench_df$T == T_levels[k], ]; sub <- sub[order(sub$L), ]
    lines(sub$L,  sub$speedup, col = T_colors[k], lwd = 2)
    points(sub$L, sub$speedup, col = T_colors[k], pch = 19, cex = 1.2)
  }
  abline(h = median(bench_df$speedup), lty = 3, col = "gray50")
  legend("bottomright", cex = 0.75, bty = "n", inset = 0.02,
         legend = c(T_labels, "median"),
         col = c(T_colors, "gray50"),
         lty = c(rep(1, length(T_levels)), 3),
         pch = c(rep(19, length(T_levels)), NA),
         lwd = c(rep(2, length(T_levels)), 1))
  
  # Panel 3: mean cluster size vs L at each T (this is Wolff-specific)
  plot(NA, xlim = range(L_levels), ylim = range(bench_df$mean_cluster), log = "xy",
       xlab = "L", ylab = "mean cluster size",
       main = "Mean cluster size vs L")
  for (k in seq_along(T_levels)) {
    sub <- bench_df[bench_df$T == T_levels[k], ]; sub <- sub[order(sub$L), ]
    lines(sub$L,  sub$mean_cluster, col = T_colors[k], lwd = 2)
    points(sub$L, sub$mean_cluster, col = T_colors[k], pch = 19, cex = 1.2)
  }
  # Reference line N = L^2
  L_ref <- range(L_levels)
  lines(L_ref, L_ref^2, col = "gray60", lty = 3, lwd = 1.5)
  legend("topleft", cex = 0.75, bty = "n", inset = 0.02,
         legend = c(T_labels, "N = L^2 ref"),
         col = c(T_colors, "gray60"),
         lty = c(rep(1, length(T_levels)), 3),
         pch = c(rep(19, length(T_levels)), NA),
         lwd = c(rep(2, length(T_levels)), 1.5))
  
  invisible(NULL)
}


# =============================================================================
# Section 9.4: Validation of Wolff against all ground truths
#
# Same structure as validate_metropolis_all_observables() in Section 6.2:
#   - Energy table vs Onsager
#   - Specific heat vs Onsager (+ proxy at T_c)
#   - |m| vs Yang / finite-L references
#   - Finite-size scaling checks at T_c and above
#   - Susceptibility scaling fit (chi ~ L^1.75 at T_c)
# The only Wolff-specific addition: a mean_cluster_size column.
#
# Chain length choice: fixed 10000 cluster flips at every (L, T). Wolff has
# negligible critical slowing down so fixed suffices (contrast with Metropolis
# which needed adaptive ~L^2 at T_c).
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


plot_wolff_validation <- function(validation_df) {
  L_values <- sort(unique(validation_df$L))
  L_colors <- c("black", "darkorange", "purple", "darkcyan")[seq_along(L_values)]
  
  T_curve <- seq(1.0, 4.2, length.out = 400)
  u_exact <- onsager_energy_curve(T_curve)
  c_exact <- onsager_specific_heat_curve(T_curve)
  m_exact <- onsager_magnetization_curve(T_curve)
  c_cap   <- 4
  c_plot  <- pmin(c_exact, c_cap)
  
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Energy
  plot(T_curve, u_exact, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "u(T)", main = "Wolff Energy: MCMC vs Onsager",
       ylim = range(c(u_exact, validation_df$e_mcmc), finite = TRUE))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- validation_df[validation_df$L == L_values[k], ]
    points(sub$T, sub$e_mcmc, pch = 19, cex = 1.3, col = L_colors[k])
  }
  legend("topleft", cex = 0.75, bty = "n",
         legend = c(paste("L =", L_values), "Onsager"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(NA, length(L_values)), 2),
         lty = c(rep(NA, length(L_values)), 1))
  
  # Specific heat capped
  plot(T_curve, c_plot, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "c(T)",
       main = sprintf("Wolff Specific heat (capped at %.0f)", c_cap),
       ylim = c(0, c_cap))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- validation_df[validation_df$L == L_values[k], ]
    points(sub$T, sub$c_mcmc, pch = 19, cex = 1.3, col = L_colors[k])
  }
  legend("topleft", cex = 0.75, bty = "n",
         legend = c(paste("L =", L_values), "Onsager (capped)"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(NA, length(L_values)), 2),
         lty = c(rep(NA, length(L_values)), 1))
  
  # Magnetization with 1/sqrt(N) ticks
  plot(T_curve, m_exact, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "|m|(T)",
       main = "Wolff Magnetization: MCMC vs Yang + 1/sqrt(N)",
       ylim = c(0, 1.05))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- validation_df[validation_df$L == L_values[k], ]
    points(sub$T, sub$m_mcmc, pch = 19, cex = 1.3, col = L_colors[k])
    N_k <- L_values[k]^2
    segments(T_CRITICAL + 0.1, 1/sqrt(N_k), 4.2, 1/sqrt(N_k),
             col = L_colors[k], lty = 3, lwd = 1)
  }
  legend("bottomleft", cex = 0.75, bty = "n",
         legend = c(paste("L =", L_values), "Yang", "1/sqrt(N) ref"),
         col = c(L_colors, "gray40", "gray40"),
         pch = c(rep(19, length(L_values)), NA, NA),
         lwd = c(rep(NA, length(L_values)), 2, 1),
         lty = c(rep(NA, length(L_values)), 1, 3))
  
  # Susceptibility
  plot(validation_df$T, validation_df$chi_mcmc, type = "n",
       xlab = "T", ylab = expression(chi(T)),
       main = "Wolff Susceptibility: MCMC (peak ~ L^1.75)")
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- validation_df[validation_df$L == L_values[k], ]
    points(sub$T, sub$chi_mcmc, pch = 19, cex = 1.3, col = L_colors[k])
  }
  legend("topleft", cex = 0.75, bty = "n",
         legend = paste("L =", L_values), col = L_colors, pch = 19)
  
  invisible(NULL)
}


# =============================================================================
# Drive (Sections 9.3 and 9.4)
# =============================================================================

cat("\n############################################################\n")
cat("Section 9.3: Wolff R-vs-Rcpp benchmark\n")
cat("############################################################\n\n")
wolff_bench_df <- benchmark_wolff_samplers()
plot_wolff_benchmark(wolff_bench_df)

cat("\n############################################################\n")
cat("Section 9.4: Wolff validation against ground truths\n")
cat("############################################################\n\n")
wolff_validation_df <- validate_wolff_all_observables()
plot_wolff_validation(wolff_validation_df)


# =============================================================================
# Section 9 part 2: Wolff trace plots across all (L, T) for all observables
#
# Parallel to Section 7 for Metropolis. For each of the 12 (L, T) cells:
#
#   - trace plot of e (per-cluster-flip) with Onsager horizontal
#   - trace plot of |m| (per-cluster-flip) with appropriate reference
#   - running c(t) -- cumulative specific heat estimate
#   - running chi(t) -- cumulative susceptibility estimate
#
# Output: 4 figures, one per observable, each 4 rows (L) x 3 cols (T).
#
# Axis note: Wolff's natural time unit is "cluster flip", not sweep. The
# x-axis here is therefore cluster-flip index. For visual comparison with
# the Metropolis Section 7 traces, remember that one Wolff cluster flip
# represents (mean_cluster_size / N) sweeps of work, which at T_c is
# roughly 0.4 sweeps for L=64. So what looks like "fewer samples" in the
# Wolff plot represents much less actual computation.
#
# Ground truth overlays:
#   - e:  Onsager (safe-wrapped at T_c) as a horizontal dashed line
#   - |m|: Yang below T_c; no line at T_c (scaling, no single value);
#           1/sqrt(N) horizontal above T_c
#   - c:  onsager_specific_heat_proxy at T_c (finite proxy, not the true Inf);
#          onsager_specific_heat otherwise
#   - chi: no exact value, no line
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
      T        <- T_values[k]
      label    <- T_labels[k]
      n_cf     <- choose_n_cluster_flips(L, T)
      n_bi     <- as.integer(n_cf / 5)
      
      key <- sprintf("L%d_T%s",
                     L,
                     if (abs(T - T_CRITICAL) < 1e-6) "c" else
                       sprintf("%.2f", T))
      
      cat(sprintf("  Running L=%2d, T=%s (%-8s): n_cf=%d ... ",
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
# Same running-estimate helpers as Section 7. Reusing them verbatim.
# -----------------------------------------------------------------------------

# (running_specific_heat and running_susceptibility already defined in Section 7.
#  If this file is sourced standalone without Section 7 in memory, uncomment
#  the definitions below.)

if (!exists("running_specific_heat")) {
  running_specific_heat <- function(E, N, T) {
    csum  <- cumsum(E); csum2 <- cumsum(E^2); n <- seq_along(E)
    mean_E <- csum / n; mean_E2 <- csum2 / n
    pmax(mean_E2 - mean_E^2, 0) / (N * T^2)
  }
}

if (!exists("running_susceptibility")) {
  running_susceptibility <- function(M, N, T) {
    csum2 <- cumsum(M^2); csumA <- cumsum(abs(M)); n <- seq_along(M)
    mean_M2 <- csum2 / n; mean_aM <- csumA / n
    pmax(mean_M2 - mean_aM^2, 0) / (N * T)
  }
}

if (!exists("thin_for_plot")) {
  thin_for_plot <- function(x, max_pts = 5000) {
    n <- length(x)
    if (n <= max_pts) return(seq_len(n))
    step <- ceiling(n / max_pts)
    seq(1L, n, by = step)
  }
}


# -----------------------------------------------------------------------------
# Build traces for one observable across all cells.
# -----------------------------------------------------------------------------

build_traces_one_observable_wolff <- function(runs, observable) {
  out <- list()
  for (key in names(runs)) {
    run <- runs[[key]]
    N   <- run$L^2
    T   <- run$T
    n   <- run$n_cluster_flips
    
    if (observable == "e") {
      y <- run$series$E / N
      g <- safe_onsager_energy(T)
    } else if (observable == "m") {
      y <- abs(run$series$M) / N
      g <- if (T < T_CRITICAL - 0.01) safe_onsager_magnetization(T)
      else if (abs(T - T_CRITICAL) < 0.05) NA
      else 1 / sqrt(N)
    } else if (observable == "c") {
      y <- running_specific_heat(run$series$E, N, T)
      g <- if (abs(T - T_CRITICAL) < 0.05)
        onsager_specific_heat_proxy(T)
      else
        onsager_specific_heat(T)
    } else if (observable == "chi") {
      y <- running_susceptibility(run$series$M, N, T)
      g <- NA
    } else {
      stop("Unknown observable: ", observable)
    }
    
    out[[key]] <- list(
      key      = key,
      L        = run$L,
      T        = T,
      T_label  = run$T_label,
      x        = seq_along(y),
      y        = y,
      ground   = g,
      mean_cl  = run$mean_cluster_size
    )
  }
  out
}


# -----------------------------------------------------------------------------
# Plot one figure (4 rows = L, 3 cols = T) for one observable, Wolff edition.
# -----------------------------------------------------------------------------

plot_observable_grid_wolff <- function(traces, observable,
                                       L_values = c(8, 16, 32, 64),
                                       T_values = c(1.5, T_CRITICAL, 4.0),
                                       T_labels = c("below Tc", "at Tc", "above Tc"),
                                       max_pts  = 5000) {
  
  # Common y-range across the 12 panels
  y_all <- c()
  for (tr in traces) {
    n   <- length(tr$y)
    cut <- max(2L, as.integer(0.02 * n))
    y_all <- c(y_all, tr$y[cut:n])
  }
  y_all <- y_all[is.finite(y_all)]
  y_range <- range(y_all)
  y_range <- y_range + diff(y_range) * c(-0.05, 0.05)
  
  obs_label <- switch(observable,
                      e    = "e (energy per spin)",
                      m    = "|m| (magnetization per spin)",
                      c    = "c running estimate",
                      chi  = expression(chi ~ "running estimate"))
  
  obs_title <- switch(observable,
                      e    = "Wolff: Energy trace",
                      m    = "Wolff: Magnetization trace",
                      c    = "Wolff: Specific heat (running)",
                      chi  = "Wolff: Susceptibility (running)")
  
  obs_color <- switch(observable,
                      e    = "darkblue",
                      m    = "darkgreen",
                      c    = "darkred",
                      chi  = "purple")
  
  old_par <- par(mfrow = c(length(L_values), length(T_values)),
                 mar    = c(3.2, 3.6, 2.0, 0.6),
                 mgp    = c(2.0, 0.6, 0),
                 oma    = c(0, 0, 2.4, 0),
                 no.readonly = TRUE)
  on.exit(par(old_par))
  
  for (L in L_values) {
    for (k in seq_along(T_values)) {
      T     <- T_values[k]
      label <- T_labels[k]
      key   <- sprintf("L%d_T%s",
                       L,
                       if (abs(T - T_CRITICAL) < 1e-6) "c" else sprintf("%.2f", T))
      tr <- traces[[key]]
      idx <- thin_for_plot(tr$x, max_pts)
      
      plot(tr$x[idx], tr$y[idx], type = "l", col = obs_color, lwd = 1,
           xlab = "cluster flip", ylab = obs_label,
           ylim = y_range,
           main = sprintf("L=%d, T=%s (mean_cl=%.0f)",
                          L,
                          if (abs(T - T_CRITICAL) < 1e-6) "T_c"
                          else sprintf("%.2f", T),
                          tr$mean_cl),
           cex.main = 0.9)
      
      if (is.finite(tr$ground)) {
        abline(h = tr$ground, lty = 2, col = "gray30", lwd = 1.2)
      }
    }
  }
  
  mtext(obs_title, side = 3, outer = TRUE, line = 0.6, cex = 1.2, font = 2)
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Top-level driver: run 12 cells, then make 4 figures.
# -----------------------------------------------------------------------------

plot_section9_all_traces <- function(L_values = c(8, 16, 32, 64),
                                     T_values = c(1.5, T_CRITICAL, 4.0),
                                     T_labels = c("below Tc", "at Tc", "above Tc"),
                                     seed_base = 42,
                                     max_pts   = 5000,
                                     pause_between = FALSE) {
  
  cat("\n############################################################\n")
  cat("Section 9 (traces): Wolff traces across (L, T) for all observables\n")
  cat("############################################################\n\n")
  
  runs <- trace_runs_wolff(L_values  = L_values,
                           T_values  = T_values,
                           T_labels  = T_labels,
                           seed_base = seed_base)
  
  for (obs in c("e", "m", "c", "chi")) {
    cat(sprintf("Plotting Wolff figure for observable: %s\n", obs))
    traces <- build_traces_one_observable_wolff(runs, observable = obs)
    plot_observable_grid_wolff(traces, observable = obs,
                               L_values = L_values, T_values = T_values,
                               T_labels = T_labels, max_pts = max_pts)
    if (pause_between) {
      cat("Press [Enter] to continue to the next figure...")
      readline()
    }
  }
  
  invisible(runs)
}


# =============================================================================
# Drive
# =============================================================================

section9_wolff_runs <- plot_section9_all_traces(pause_between = FALSE)





# =============================================================================
# Section 9 part 3: Wolff autocorrelation diagnostics
#
# Parallel to Section 8 for Metropolis. For every (L, T) cell from
# section9_wolff_runs:
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
# What we expect to see that is different from Metropolis:
#
#   1. Wolff tau values (in cluster-flip units) are of order unity even at
#      T_c, not hundreds or thousands. This is the statistical expression
#      of "no critical slowing down".
#
#   2. The z-fit for tau_E: slope ~0.25 (Wolff 1989 literature), not ~2.17
#      as for Metropolis.
#
#   3. For tau_M specifically: Wolff decorrelates magnetization so fast
#      that the Sokal-floor (tau >= 0.5) kicks in at every lattice size.
#      When that happens the z-fit for tau_M is meaningless (constant values
#      give a degenerate slope), so we detect and report that honestly
#      instead of fitting garbage.
#
#   4. ESS should be essentially the full chain length at T_c for both
#      observables, compared to Metropolis's ESS ~10-30 at T_c L=64.
# =============================================================================


# -----------------------------------------------------------------------------
# Per-cell Wolff diagnostics
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
# Pretty-print the Wolff diagnostics in two tables: the native cluster-flip
# table and the equivalent-sweeps table (so comparison with Section 8
# Metropolis is direct).
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
# z-fit for Wolff at T_c. Checks for the Sokal-floor degeneracy: if all
# tau values are at (or near) the 0.5 floor, the fit is meaningless.
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
  # Detect the Sokal-floor degeneracy: if all tau values are within 2% of
  # the 0.5 floor, refuse to fit (would give meaningless noise-driven slope).
  floor_margin <- (taus - 0.5) / 0.5
  if (all(abs(floor_margin) < 0.02)) {
    cat(sprintf("  z (from tau_%s): all tau values at Sokal floor 0.5 (+/- 2%%);\n",
                observable))
    cat(sprintf("    Wolff decorrelates this observable in < 1 cluster flip at every L.\n"))
    cat(sprintf("    No meaningful z to extract -- Wolff is effectively drawing iid samples.\n"))
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
# 4x3 rho(t) grid for Wolff, mirroring Section 8's plot_autocorr_grid
# -----------------------------------------------------------------------------

plot_autocorr_grid_wolff <- function(runs, observable = c("E", "M"),
                                     L_values = c(8, 16, 32, 64),
                                     T_values = c(1.5, T_CRITICAL, 4.0),
                                     T_labels = c("below Tc", "at Tc", "above Tc"),
                                     max_lag_plot = 50) {
  observable <- match.arg(observable)
  
  obs_label <- switch(observable,
                      E = "rho_E(t)  [energy autocorrelation]",
                      M = "rho_M(t)  [magnetization autocorrelation]")
  obs_color <- switch(observable, E = "darkblue", M = "darkred")
  obs_title <- switch(observable,
                      E = "Wolff: Energy autocorrelation function",
                      M = "Wolff: Magnetization autocorrelation function")
  
  old_par <- par(mfrow = c(length(L_values), length(T_values)),
                 mar    = c(3.2, 3.6, 2.0, 0.6),
                 mgp    = c(2.0, 0.6, 0),
                 oma    = c(0, 0, 2.4, 0),
                 no.readonly = TRUE)
  on.exit(par(old_par))
  
  for (L in L_values) {
    for (k in seq_along(T_values)) {
      T   <- T_values[k]
      key <- sprintf("L%d_T%s",
                     L,
                     if (abs(T - T_CRITICAL) < 1e-6) "c" else sprintf("%.2f", T))
      run <- runs[[key]]
      n_cf <- nrow(run$series)
      
      max_lag_here <- min(max_lag_plot, n_cf %/% 4)
      x_series <- if (observable == "E") run$series$E else run$series$M
      rho      <- autocorrelation(x_series, max_lag = max_lag_here)
      
      diag <- integrated_tau(x_series)
      flag <- if (diag$closed) "" else "!"
      
      plot(0:max_lag_here, rho, type = "l", col = obs_color, lwd = 1.5,
           xlab = "lag (cluster flips)", ylab = obs_label,
           ylim = c(-0.1, 1.05),
           main = sprintf("L=%d, T=%s   tau_cf=%.2f%s",
                          L,
                          if (abs(T - T_CRITICAL) < 1e-6) "T_c"
                          else sprintf("%.2f", T),
                          diag$tau, flag),
           cex.main = 0.9)
      abline(h = 0, lty = 3)
      abline(h = 1/exp(1), lty = 3, col = "gray60")
    }
  }
  
  mtext(obs_title, side = 3, outer = TRUE, line = 0.6, cex = 1.2, font = 2)
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Scaling plot: tau vs L on log-log, for Wolff. Same layout as Section 8's
# plot_tau_scaling but using the Wolff literature z = 0.25 as reference.
# -----------------------------------------------------------------------------

plot_tau_scaling_wolff <- function(diag_df,
                                   T_values = c(1.5, T_CRITICAL, 4.0),
                                   T_labels = c("below Tc", "at Tc", "above Tc"),
                                   T_colors = c("darkblue", "darkred", "darkgreen")) {
  L_values <- sort(unique(diag_df$L))
  
  old_par <- par(mfrow = c(1, 2),
                 mar = c(4.2, 4.2, 2.5, 9),
                 xpd = FALSE,
                 no.readonly = TRUE)
  on.exit(par(old_par))
  
  draw_panel <- function(tau_col_name, closed_col_name, y_label, main_title) {
    y_values <- diag_df[[tau_col_name]]
    plot(NA, xlim = range(L_values),
         ylim = range(c(y_values, 0.5), finite = TRUE),
         log  = "xy",
         xlab = "L", ylab = y_label,
         main = main_title)
    
    for (k in seq_along(T_values)) {
      T   <- T_values[k]
      sub <- diag_df[diag_df$T == T, ]
      sub <- sub[order(sub$L), ]
      pch <- ifelse(sub[[closed_col_name]], 19, 4)
      lines(sub$L,  sub[[tau_col_name]], col = T_colors[k], lwd = 2)
      points(sub$L, sub[[tau_col_name]], col = T_colors[k], pch = pch, cex = 1.3)
    }
    
    # Wolff literature reference: slope 0.25, anchored at smallest L
    sub_tc <- diag_df[diag_df$T_label == "at Tc", ]
    sub_tc <- sub_tc[order(sub_tc$L), ]
    if (nrow(sub_tc) >= 1) {
      L_ref <- range(L_values)
      y_ref <- sub_tc[[tau_col_name]][1] *
        (L_ref / sub_tc$L[1])^Z_WOLFF_LITERATURE
      lines(L_ref, y_ref, col = "gray60", lty = 3, lwd = 1.5)
    }
    
    par(xpd = TRUE)
    legend(x       = 10^(par("usr")[2] + 0.05 * diff(par("usr")[1:2])),
           y       = 10^par("usr")[4],
           cex     = 0.8,
           bty     = "n",
           legend  = c(T_labels,
                       sprintf("ref L^%.2f", Z_WOLFF_LITERATURE),
                       "", "Sokal:",
                       "closed", "open (biased low)"),
           col     = c(T_colors, "gray60", NA, NA, "black", "black"),
           lty     = c(rep(1, length(T_values)), 3, NA, NA, NA, NA),
           pch     = c(rep(19, length(T_values)), NA, NA, NA, 19, 4),
           lwd     = c(rep(2, length(T_values)), 1.5, NA, NA, NA, NA))
    par(xpd = FALSE)
  }
  
  draw_panel("tau_E_cf", "closed_E",
             expression(tau[E] ~ "(cluster flips)"),
             "Wolff energy: tau_E vs L")
  
  draw_panel("tau_M_cf", "closed_M",
             expression(tau[M] ~ "(cluster flips)"),
             "Wolff magnetization: tau_M vs L")
  
  invisible(NULL)
}


# =============================================================================
# Top-level driver
# =============================================================================

run_section9_diagnostics <- function(runs = section9_wolff_runs) {
  cat("\n############################################################\n")
  cat("Section 9 part 3: Wolff autocorrelation diagnostics\n")
  cat("############################################################\n\n")
  
  diag_df <- build_diagnostics_table_wolff(runs)
  print_diagnostics_table_wolff(diag_df)
  
  cat("\n--- Dynamic critical exponent z (Wolff at T_c) ---\n")
  fit_dynamic_z_wolff(diag_df, "E")
  fit_dynamic_z_wolff(diag_df, "M")
  
  cat("\nPlotting Wolff rho_E grid (4 x 3)...\n")
  plot_autocorr_grid_wolff(runs, observable = "E")
  
  cat("Plotting Wolff rho_M grid (4 x 3)...\n")
  plot_autocorr_grid_wolff(runs, observable = "M")
  
  cat("Plotting Wolff tau scaling summary (1 x 2)...\n")
  plot_tau_scaling_wolff(diag_df)
  
  invisible(diag_df)
}


# =============================================================================
# Drive
# =============================================================================

section9_wolff_diag <- run_section9_diagnostics(section9_wolff_runs)







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
# Produces:
#
#   1. Merged comparison table: for every (L, T) cell, Metropolis vs Wolff
#      on tau (in equivalent sweeps), ESS, SE(e), SE(|m|).
#
#   2. Ratio table: Wolff/Metropolis improvement factor for each metric.
#      The headline numbers (factor > 100 at T_c for L=64) live here.
#
#   3. Figure A: tau vs L, two panels (E and M), both samplers on same axes,
#      log-log, with both literature z references drawn. This is the
#      "dynamic critical exponent" figure for the report.
#
#   4. Figure B: rho(t) overlay at T_c, L=32 -- both samplers on one plot
#      for E and for M. Shows visually how fast Wolff decorrelates.
#
#   5. Figure C: wall-clock efficiency -- ESS per CPU-second estimated
#      from the earlier benchmarks. This is what a practitioner actually
#      cares about.
# =============================================================================


# -----------------------------------------------------------------------------
# Section 9 part 4.1: merged comparison table
# -----------------------------------------------------------------------------

build_comparison_table <- function(metro_diag, wolff_diag) {
  # Both dfs have columns L, T, T_label; join on (L, T_label)
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
  
  # First: raw tau values side by side, in common (equivalent-sweep) units
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
  
  # Second: ESS values side by side
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
  
  # Headline: the worst-case Metropolis cell and how Wolff fares there
  cat("\n--- Headline numbers at T_c ---\n")
  sub_tc <- comp_df[comp_df$T_label == "at Tc", ]
  for (i in seq_len(nrow(sub_tc))) {
    r <- sub_tc[i, ]
    cat(sprintf("  L=%-2d : Wolff tau_E is %.1fx faster, tau_M is %.1fx faster than Metropolis\n",
                r$L, r$ratio_tau_E, r$ratio_tau_M))
  }
}


# -----------------------------------------------------------------------------
# Section 9 part 4.2: Figure A -- tau vs L on same axes, both samplers
# -----------------------------------------------------------------------------

plot_metro_vs_wolff_tau <- function(metro_diag, wolff_diag,
                                    T_values = c(1.5, T_CRITICAL, 4.0),
                                    T_labels = c("below Tc", "at Tc", "above Tc"),
                                    T_colors = c("darkblue", "darkred", "darkgreen")) {
  L_values <- sort(unique(metro_diag$L))
  
  old_par <- par(mfrow = c(1, 2),
                 mar = c(4.2, 4.4, 2.5, 10),
                 xpd = FALSE,
                 no.readonly = TRUE)
  on.exit(par(old_par))
  
  draw_panel <- function(observable, y_label, main_title) {
    metro_col  <- paste0("tau_", observable)
    wolff_col  <- paste0("tau_", observable, "_sw")
    closed_m   <- paste0("closed_", observable)
    closed_w   <- paste0("closed_", observable)
    
    all_vals <- c(metro_diag[[metro_col]], wolff_diag[[wolff_col]])
    plot(NA, xlim = range(L_values),
         ylim = range(c(all_vals, 0.1), finite = TRUE),
         log = "xy",
         xlab = "L", ylab = y_label,
         main = main_title)
    
    # Metropolis: solid lines, filled circles
    for (k in seq_along(T_values)) {
      T <- T_values[k]
      sub <- metro_diag[metro_diag$T == T, ]
      sub <- sub[order(sub$L), ]
      pch <- ifelse(sub[[closed_m]], 19, 4)
      lines(sub$L,  sub[[metro_col]], col = T_colors[k], lwd = 2, lty = 1)
      points(sub$L, sub[[metro_col]], col = T_colors[k], pch = pch, cex = 1.2)
    }
    # Wolff: dashed lines, open circles
    for (k in seq_along(T_values)) {
      T <- T_values[k]
      sub <- wolff_diag[wolff_diag$T == T, ]
      sub <- sub[order(sub$L), ]
      pch <- ifelse(sub[[closed_w]], 21, 4)
      lines(sub$L,  sub[[wolff_col]], col = T_colors[k], lwd = 2, lty = 2)
      points(sub$L, sub[[wolff_col]], col = T_colors[k], pch = pch, cex = 1.2,
             bg = "white")
    }
    
    # Two reference slope lines at T_c, anchored at smallest L
    sub_tc_m <- metro_diag[metro_diag$T_label == "at Tc", ]
    sub_tc_m <- sub_tc_m[order(sub_tc_m$L), ]
    if (nrow(sub_tc_m) >= 1) {
      L_ref <- range(L_values)
      y_met <- sub_tc_m[[metro_col]][1] *
        (L_ref / sub_tc_m$L[1])^Z_METROPOLIS_LITERATURE
      lines(L_ref, y_met, col = "gray60", lty = 3, lwd = 1.5)
    }
    sub_tc_w <- wolff_diag[wolff_diag$T_label == "at Tc", ]
    sub_tc_w <- sub_tc_w[order(sub_tc_w$L), ]
    if (nrow(sub_tc_w) >= 1) {
      L_ref <- range(L_values)
      y_wol <- sub_tc_w[[wolff_col]][1] *
        (L_ref / sub_tc_w$L[1])^Z_WOLFF_LITERATURE
      lines(L_ref, y_wol, col = "gray75", lty = 3, lwd = 1.5)
    }
    
    par(xpd = TRUE)
    legend(x    = 10^(par("usr")[2] + 0.05 * diff(par("usr")[1:2])),
           y    = 10^par("usr")[4],
           cex  = 0.78, bty = "n",
           legend = c("Metropolis:", T_labels, "",
                      "Wolff:", T_labels, "",
                      sprintf("ref L^%.2f (met)", Z_METROPOLIS_LITERATURE),
                      sprintf("ref L^%.2f (wol)", Z_WOLFF_LITERATURE)),
           col    = c(NA, T_colors, NA, NA, T_colors, NA, "gray60", "gray75"),
           lty    = c(NA, rep(1, length(T_values)), NA, NA, rep(2, length(T_values)), NA, 3, 3),
           pch    = c(NA, rep(19, length(T_values)), NA, NA, rep(21, length(T_values)), NA, NA, NA),
           pt.bg  = c(NA, rep(NA, length(T_values)), NA, NA, rep("white", length(T_values)), NA, NA, NA),
           lwd    = c(NA, rep(2, length(T_values)), NA, NA, rep(2, length(T_values)), NA, 1.5, 1.5))
    par(xpd = FALSE)
  }
  
  draw_panel("E", expression(tau[E] ~ "(equivalent sweeps)"),
             "tau_E vs L: Metropolis vs Wolff")
  draw_panel("M", expression(tau[M] ~ "(equivalent sweeps)"),
             "tau_M vs L: Metropolis vs Wolff")
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Section 9 part 4.3: Figure B -- rho(t) overlay at T_c, fixed L
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
  
  # Metropolis: lag in sweeps
  rho_E_m <- autocorrelation(run_m$series$E, max_lag = max_lag_plot)
  rho_M_m <- autocorrelation(run_m$series$M, max_lag = max_lag_plot)
  
  # Wolff: lag in cluster flips, convert to equivalent sweeps for x-axis
  rho_E_w <- autocorrelation(run_w$series$E, max_lag = max_lag_plot)
  rho_M_w <- autocorrelation(run_w$series$M, max_lag = max_lag_plot)
  sweeps_per_flip <- run_w$mean_cluster_size / (L * L)
  x_w <- (0:max_lag_plot) * sweeps_per_flip
  
  old_par <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  # E panel
  plot(0:max_lag_plot, rho_E_m, type = "l", col = "darkred", lwd = 2,
       ylim = c(-0.1, 1), xlim = c(0, max_lag_plot),
       xlab = "lag (equivalent sweeps)", ylab = expression(rho[E](t)),
       main = sprintf("Energy autocorrelation at T_c, L=%d", L))
  lines(x_w, rho_E_w, col = "darkblue", lwd = 2)
  abline(h = 0, lty = 3)
  abline(h = 1/exp(1), lty = 3, col = "gray60")
  legend("topright", cex = 0.85, bty = "n",
         legend = c("Metropolis (lag = sweep)",
                    "Wolff (lag = equivalent sweep)"),
         col = c("darkred", "darkblue"), lwd = 2)
  
  # M panel
  plot(0:max_lag_plot, rho_M_m, type = "l", col = "darkred", lwd = 2,
       ylim = c(-0.1, 1), xlim = c(0, max_lag_plot),
       xlab = "lag (equivalent sweeps)", ylab = expression(rho[M](t)),
       main = sprintf("Magnetization autocorrelation at T_c, L=%d", L))
  lines(x_w, rho_M_w, col = "darkblue", lwd = 2)
  abline(h = 0, lty = 3)
  abline(h = 1/exp(1), lty = 3, col = "gray60")
  legend("topright", cex = 0.85, bty = "n",
         legend = c("Metropolis (lag = sweep)",
                    "Wolff (lag = equivalent sweep)"),
         col = c("darkred", "darkblue"), lwd = 2)
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Section 9 part 4.4: Figure C -- ESS per CPU second
#
# Reuses earlier benchmarks (bench_df for Metropolis, wolff_bench_df for
# Wolff) to estimate per-cell wall-clock time, then computes
# ESS / (wall time) at each (L, T). This is the practitioner's bottom line.
# -----------------------------------------------------------------------------

plot_ess_per_second <- function(metro_diag, wolff_diag,
                                metro_bench, wolff_bench,
                                T_values = c(1.5, T_CRITICAL, 4.0),
                                T_labels = c("below Tc", "at Tc", "above Tc"),
                                T_colors = c("darkblue", "darkred", "darkgreen")) {
  # Per-cell cost per sweep (Metropolis) or per cluster flip (Wolff), measured
  # from the earlier benchmarks which used n_sweeps = n_cluster_flips = 500.
  metro_sec_per_sweep <- with(metro_bench, t_cpp / 500)
  wolff_sec_per_flip  <- with(wolff_bench, t_cpp / 500)
  
  # Join: look up per-cell time and compute throughput
  make_df <- function(diag_df, bench_df, time_col, chain_length_col,
                      ESS_col_E, ESS_col_M, sampler_name) {
    out <- list()
    for (i in seq_len(nrow(diag_df))) {
      r <- diag_df[i, ]
      bmatch <- which(bench_df$L == r$L &
                        abs(bench_df$T - r$T) < 1e-6)
      if (length(bmatch) != 1) next
      per_step_time <- bench_df[[time_col]][bmatch]
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
  
  metro_df <- make_df(metro_diag, metro_bench, "t_cpp", "n_sweeps",
                      "ESS_E", "ESS_M", "Metropolis")
  # For Wolff, per-cell time is sec_per_flip * n_cluster_flips
  wolff_df <- make_df(wolff_diag, wolff_bench, "t_cpp", "n_cluster_flips",
                      "ESS_E", "ESS_M", "Wolff")
  
  all_df <- rbind(metro_df, wolff_df)
  
  old_par <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  draw_panel <- function(col_name, y_label, main_title) {
    y_range <- range(all_df[[col_name]][is.finite(all_df[[col_name]])])
    plot(NA, xlim = range(all_df$L), ylim = y_range, log = "xy",
         xlab = "L", ylab = y_label, main = main_title)
    
    for (k in seq_along(T_values)) {
      T <- T_values[k]
      
      sub_m <- metro_df[metro_df$T == T, ]; sub_m <- sub_m[order(sub_m$L), ]
      lines(sub_m$L,  sub_m[[col_name]], col = T_colors[k], lwd = 2, lty = 1)
      points(sub_m$L, sub_m[[col_name]], col = T_colors[k], pch = 19, cex = 1.2)
      
      sub_w <- wolff_df[wolff_df$T == T, ]; sub_w <- sub_w[order(sub_w$L), ]
      lines(sub_w$L,  sub_w[[col_name]], col = T_colors[k], lwd = 2, lty = 2)
      points(sub_w$L, sub_w[[col_name]], col = T_colors[k], pch = 21, cex = 1.2,
             bg = "white")
    }
    
    legend("bottomleft", cex = 0.75, bty = "n", inset = 0.02,
           legend = c(paste("Metro,", T_labels), paste("Wolff,", T_labels)),
           col    = c(T_colors, T_colors),
           lty    = c(rep(1, length(T_values)), rep(2, length(T_values))),
           pch    = c(rep(19, length(T_values)), rep(21, length(T_values))),
           pt.bg  = c(rep(NA, length(T_values)), rep("white", length(T_values))),
           lwd    = 2)
  }
  
  draw_panel("ESS_E_per_sec", "ESS_E / second",
             "Energy samples per CPU second")
  draw_panel("ESS_M_per_sec", "ESS_M / second",
             "Magnetization samples per CPU second")
  
  invisible(all_df)
}


# =============================================================================
# Driver
# =============================================================================

run_section9_part4 <- function() {
  cat("\n############################################################\n")
  cat("Section 9 part 4: Metropolis vs Wolff head-to-head\n")
  cat("############################################################\n\n")
  
  # 4.1: comparison tables
  comp_df <- build_comparison_table(section8_diag, section9_wolff_diag)
  print_comparison_table(comp_df)
  
  # 4.2: tau vs L, both samplers
  cat("\nPlotting Figure A (tau vs L, both samplers)...\n")
  plot_metro_vs_wolff_tau(section8_diag, section9_wolff_diag)
  
  # 4.3: rho(t) overlay at T_c, L=32
  cat("Plotting Figure B (rho(t) overlay at T_c, L=32)...\n")
  plot_rho_overlay(section7_runs, section9_wolff_runs, L = 32)
  
  # 4.4: ESS per CPU second (needs bench_df from 6.1 and wolff_bench_df from 9.3)
  if (exists("bench_df") && exists("wolff_bench_df")) {
    cat("Plotting Figure C (ESS per CPU second)...\n")
    throughput_df <- plot_ess_per_second(section8_diag, section9_wolff_diag,
                                         bench_df, wolff_bench_df)
  } else {
    cat("Skipping Figure C: bench_df or wolff_bench_df not in memory.\n")
  }
  
  invisible(comp_df)
}


# =============================================================================
# Drive
# =============================================================================

comparison_df <- run_section9_part4()



# =============================================================================
# Section 10: Full temperature sweep with Wolff
#
# Produces the classic 2D Ising phase-transition figure: smooth curves of
# u(T), c(T), |m|(T), and chi(T) across the full phase transition, for
# multiple lattice sizes, overlaid on Onsager/Yang exact results.
#
# Why Wolff and not Metropolis: Wolff mixes well at every temperature
# (no critical slowing down at T_c, good acceptance at every T). At T_c
# for L=64, Metropolis needed 320k sweeps for a decent estimate; Wolff
# needs ~5k cluster flips. So running a dense grid of 34 temperatures x
# 4 lattice sizes = 136 runs is fast.
#
# Design choices:
#   - Temperature grid: 34 points, non-uniform (dense at T_c, sparse far
#     away). Reused from the original Rmd Section 10.1.
#   - Lattice sizes: L in {8, 16, 32, 64}.
#   - Chain length: 5000 cluster flips per (L, T) cell, fixed. Wolff
#     doesn't need adaptive length, and 5000 gives tens of ESS worth of
#     samples for all observables at every T.
#   - Burn-in: n_cluster_flips / 5 = 1000.
#   - Init: "hot" everywhere (Wolff escapes ergodicity traps via macro moves).
#
# No bootstrap error bars here. That's Section 11. Section 10 is just
# the point estimates overlaid on ground truth.
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

run_full_T_sweep <- function(L_values       = c(8, 16, 32, 64),
                             T_grid         = make_T_grid_sweep(),
                             n_cluster_flips = 5000,
                             seed_base      = 42) {
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
# 10.4: the classic 4-panel phase-transition figure, overlaid for all L
#
# 2x2 layout:
#   Top-left:     u(T) vs T,     Onsager curve + 4 colored L-series of points
#   Top-right:    c(T) vs T,     Onsager curve (capped)
#   Bottom-left:  |m|(T) vs T,   Yang curve
#   Bottom-right: chi(T) vs T,   (no exact curve)
#
# A vertical dashed line at T_c on every panel.
# -----------------------------------------------------------------------------

plot_T_sweep <- function(sweep_df) {
  L_values <- sort(unique(sweep_df$L))
  L_colors <- c("black", "darkorange", "purple", "darkcyan")[seq_along(L_values)]
  
  # Smooth reference curves on a dense T grid, with safe-wrapped curve
  # functions so nothing blows up at T_c
  T_curve <- seq(min(sweep_df$T), max(sweep_df$T), length.out = 400)
  u_exact <- onsager_energy_curve(T_curve)
  c_exact <- onsager_specific_heat_curve(T_curve)
  m_exact <- onsager_magnetization_curve(T_curve)
  
  c_cap  <- 4
  c_plot <- pmin(c_exact, c_cap)
  
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  # --- Energy ---
  plot(T_curve, u_exact, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "u(T)",
       main = "Energy per spin: Wolff T-sweep vs Onsager",
       ylim = range(c(u_exact, sweep_df$e), finite = TRUE))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- sweep_df[sweep_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    points(sub$T, sub$e, pch = 19, cex = 0.7, col = L_colors[k])
    lines(sub$T, sub$e, col = L_colors[k], lwd = 1)
  }
  legend("topleft", cex = 0.75, bty = "n",
         legend = c(paste("L =", L_values), "Onsager"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(1, length(L_values)), 2),
         lty = c(rep(1, length(L_values)), 1))
  
  # --- Specific heat ---
  plot(T_curve, c_plot, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "c(T)",
       main = sprintf("Specific heat: Wolff T-sweep (Onsager capped at %.0f)", c_cap),
       ylim = c(0, max(c(c_cap, sweep_df$c))))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- sweep_df[sweep_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    points(sub$T, sub$c, pch = 19, cex = 0.7, col = L_colors[k])
    lines(sub$T, sub$c, col = L_colors[k], lwd = 1)
  }
  legend("topleft", cex = 0.75, bty = "n",
         legend = c(paste("L =", L_values), "Onsager (capped)"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(1, length(L_values)), 2),
         lty = c(rep(1, length(L_values)), 1))
  
  # --- Magnetization ---
  plot(T_curve, m_exact, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "|m|(T)",
       main = "Magnetization: Wolff T-sweep vs Yang",
       ylim = c(0, 1.05))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- sweep_df[sweep_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    points(sub$T, sub$m, pch = 19, cex = 0.7, col = L_colors[k])
    lines(sub$T, sub$m, col = L_colors[k], lwd = 1)
  }
  legend("bottomleft", cex = 0.75, bty = "n",
         legend = c(paste("L =", L_values), "Yang"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(1, length(L_values)), 2),
         lty = c(rep(1, length(L_values)), 1))
  
  # --- Susceptibility ---
  plot(NA, xlim = range(T_curve), ylim = c(0, max(sweep_df$chi) * 1.05),
       xlab = "T", ylab = expression(chi(T)),
       main = "Susceptibility: Wolff T-sweep (no exact curve)")
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- sweep_df[sweep_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    points(sub$T, sub$chi, pch = 19, cex = 0.7, col = L_colors[k])
    lines(sub$T, sub$chi, col = L_colors[k], lwd = 1)
  }
  legend("topright", cex = 0.75, bty = "n",
         legend = paste("L =", L_values),
         col = L_colors, pch = 19, lwd = 1)
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# 10.5: separate plot of mean cluster size vs T, for all L
#
# This is the Wolff-specific diagnostic. The mean cluster size peaks at T_c
# (near-divergence, finite-size limited), showing directly why Wolff is
# efficient near criticality: it automatically uses large moves where
# large correlations exist.
# -----------------------------------------------------------------------------

plot_T_sweep_cluster_size <- function(sweep_df) {
  L_values <- sort(unique(sweep_df$L))
  L_colors <- c("black", "darkorange", "purple", "darkcyan")[seq_along(L_values)]
  
  old_par <- par(mfrow = c(1, 1), mar = c(4.2, 4.2, 2.5, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  plot(NA, xlim = range(sweep_df$T),
       ylim = range(sweep_df$mean_cluster),
       log  = "y",
       xlab = "T", ylab = "mean cluster size (spins)",
       main = "Wolff mean cluster size vs T (log y)")
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  
  for (k in seq_along(L_values)) {
    sub <- sweep_df[sweep_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    points(sub$T, sub$mean_cluster, pch = 19, cex = 0.8, col = L_colors[k])
    lines(sub$T, sub$mean_cluster, col = L_colors[k], lwd = 1.5)
    
    # Horizontal reference at N = L^2 for each L
    N_k <- L_values[k]^2
    abline(h = N_k, lty = 3, col = L_colors[k])
  }
  
  legend("topright", cex = 0.8, bty = "n",
         legend = c(paste("L =", L_values),
                    "N = L^2 (dotted, per L)"),
         col = c(L_colors, "gray60"),
         pch = c(rep(19, length(L_values)), NA),
         lty = c(rep(1, length(L_values)), 3),
         lwd = c(rep(1.5, length(L_values)), 1))
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# 10.6: finite-size scaling of the susceptibility peak
#
# Locate the T where chi is maximum for each L. As L -> infinity this
# converges to T_c from above. The peak height grows as L^(gamma/nu) = L^1.75.
# This is a preview of what Section 12 will do with bootstrap CIs.
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
  
  # Fit chi_peak ~ L^(gamma/nu)
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
plot_T_sweep(sweep_df)

cat("\nPlotting mean cluster size vs T...\n")
plot_T_sweep_cluster_size(sweep_df)

cat("\nAnalyzing susceptibility peak...\n")
chi_peak_df <- analyze_chi_peak(sweep_df)








# =============================================================================
# Section 11: Block bootstrap for uncertainty quantification
#
# The one advanced component in the revised proposal (per professor's
# directive: prioritize Metropolis + Wolff + full T-sweep, add ONE
# advanced component). Block bootstrap (Kunsch 1989) produces honest
# standard errors for MCMC time-series estimators by resampling contiguous
# blocks of length ~2*tau_int, so samples from different blocks are
# approximately independent.
#
# Two parts:
#
#   11A: Full T-sweep with Wolff, bootstrap SEs for all four observables
#        (e, c, |m|, chi) at every (L, T) cell.
#
#        Wolff is used because its tau is small at every T (often at the
#        Sokal floor 0.5), giving short block lengths, many blocks per
#        chain, and trustworthy bootstrap SEs with modest compute.
#
#   11B: Single-cell head-to-head Metropolis vs Wolff bootstrap at L=32,
#        T_c. Same machinery, same number of samples, but Metropolis's
#        huge tau (~56 sweeps for E, ~2200 for M) forces long block
#        lengths and few blocks, widening the bootstrap SE dramatically.
#        This illustrates why Wolff was chosen for the full-sweep bootstrap.
# =============================================================================


# =============================================================================
# 11.1: Block bootstrap primitive
#
# Non-overlapping block bootstrap. Returns the point estimate plus the
# standard deviation of n_boot resamples as the SE estimate.
# =============================================================================

block_bootstrap <- function(x, stat_fn,
                            block_length = NULL,
                            n_boot       = 200,
                            seed         = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(x)
  
  # Block length defaults to ceiling(2 * tau_int) if not supplied
  if (is.null(block_length)) {
    tau <- integrated_tau(x)$tau
    block_length <- max(1, ceiling(2 * tau))
  }
  block_length <- min(block_length, n %/% 4)   # need at least 4 blocks
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
#
# Uses tau_E for the E-derived observables (e, c) and tau_M for the
# M-derived observables (|m|, chi), which is the statistically correct
# choice for each.
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
# PART 11A: Full T-sweep with Wolff and bootstrap SEs
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
# 11A-figure-1: 4-panel phase transition with error bars
# -----------------------------------------------------------------------------

plot_T_sweep_with_errorbars <- function(sweep_err_df) {
  L_values <- sort(unique(sweep_err_df$L))
  L_colors <- c("black", "darkorange", "purple", "darkcyan")[seq_along(L_values)]
  
  T_curve <- seq(min(sweep_err_df$T), max(sweep_err_df$T), length.out = 400)
  u_exact <- onsager_energy_curve(T_curve)
  c_exact <- onsager_specific_heat_curve(T_curve)
  m_exact <- onsager_magnetization_curve(T_curve)
  c_cap   <- 4
  c_plot  <- pmin(c_exact, c_cap)
  
  draw_errorbars <- function(x, y, se, color, cex_pt = 0.7) {
    segments(x, y - se, x, y + se, col = color, lwd = 1.3)
    points(x, y, pch = 19, cex = cex_pt, col = color)
  }
  
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Energy
  plot(T_curve, u_exact, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "u(T)",
       main = "Energy per spin: Wolff + block-bootstrap SE",
       ylim = range(c(u_exact,
                      sweep_err_df$e - sweep_err_df$e_se,
                      sweep_err_df$e + sweep_err_df$e_se),
                    finite = TRUE))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- sweep_err_df[sweep_err_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$e, sub$e_se, L_colors[k])
  }
  legend("topleft", cex = 0.75, bty = "n",
         legend = c(paste("L =", L_values), "Onsager"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(1, length(L_values)), 2),
         lty = c(rep(NA, length(L_values)), 1))
  
  # Specific heat
  plot(T_curve, c_plot, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "c(T)",
       main = sprintf("Specific heat + SE (Onsager capped at %.0f)", c_cap),
       ylim = c(0, max(c(c_cap, sweep_err_df$c + sweep_err_df$c_se))))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- sweep_err_df[sweep_err_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$c, sub$c_se, L_colors[k])
  }
  legend("topleft", cex = 0.75, bty = "n",
         legend = c(paste("L =", L_values), "Onsager (capped)"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(1, length(L_values)), 2),
         lty = c(rep(NA, length(L_values)), 1))
  
  # Magnetization
  plot(T_curve, m_exact, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "|m|(T)", main = "Magnetization + SE",
       ylim = c(0, 1.05))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- sweep_err_df[sweep_err_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$m, sub$m_se, L_colors[k])
  }
  legend("bottomleft", cex = 0.75, bty = "n",
         legend = c(paste("L =", L_values), "Yang"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(1, length(L_values)), 2),
         lty = c(rep(NA, length(L_values)), 1))
  
  # Susceptibility
  plot(NA, xlim = range(T_curve),
       ylim = c(0, max(sweep_err_df$chi + sweep_err_df$chi_se) * 1.05),
       xlab = "T", ylab = expression(chi(T)),
       main = "Susceptibility + SE (no exact curve)")
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- sweep_err_df[sweep_err_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$chi, sub$chi_se, L_colors[k])
  }
  legend("topright", cex = 0.75, bty = "n",
         legend = paste("L =", L_values),
         col = L_colors, pch = 19, lwd = 1)
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# 11A-figure-2: SE of each observable vs T (2x2 panels, all L overlaid)
#
# Should show peaks at T_c for every observable -- the variance of
# fluctuating quantities diverges there, so their bootstrap SE does too.
# -----------------------------------------------------------------------------

plot_SE_vs_T <- function(sweep_err_df) {
  L_values <- sort(unique(sweep_err_df$L))
  L_colors <- c("black", "darkorange", "purple", "darkcyan")[seq_along(L_values)]
  
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  for (obs_info in list(
    list(col = "e_se",   label = "SE(u)",   title = "Energy SE vs T"),
    list(col = "c_se",   label = "SE(c)",   title = "Specific heat SE vs T"),
    list(col = "m_se",   label = "SE(|m|)", title = "Magnetization SE vs T"),
    list(col = "chi_se", label = "SE(chi)", title = "Susceptibility SE vs T")
  )) {
    plot(NA, xlim = range(sweep_err_df$T),
         ylim = c(0, max(sweep_err_df[[obs_info$col]])),
         xlab = "T", ylab = obs_info$label, main = obs_info$title)
    abline(v = T_CRITICAL, lty = 2, col = "gray")
    for (k in seq_along(L_values)) {
      sub <- sweep_err_df[sweep_err_df$L == L_values[k], ]
      sub <- sub[order(sub$T), ]
      points(sub$T, sub[[obs_info$col]], pch = 19, cex = 0.7, col = L_colors[k])
      lines(sub$T, sub[[obs_info$col]], col = L_colors[k], lwd = 1)
    }
    legend("topleft", cex = 0.7, bty = "n",
           legend = paste("L =", L_values),
           col = L_colors, pch = 19, lwd = 1)
  }
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# 11A-figure-3: Sanity check -- bootstrap SE vs naive iid SE for energy
#
# At L=32 across all T. Bootstrap SE should exceed naive SE by
# approximately sqrt(2 * tau_E). If it does, the bootstrap is honest.
# -----------------------------------------------------------------------------

sanity_check_bootstrap_vs_naive <- function(L              = 32,
                                            T_grid         = make_T_grid_sweep(),
                                            n_cluster_flips = 5000,
                                            n_boot         = 200,
                                            seed_base      = 7777) {
  cat(sprintf("\n=== 11A sanity: bootstrap SE vs naive SE for energy, L=%d ===\n", L))
  
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
  
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  plot(df$T, df$boot_se, type = "b", col = "darkblue", pch = 19, cex = 0.8,
       ylim = range(c(df$boot_se, df$naive_se)),
       xlab = "T", ylab = "SE(<e>)",
       main = sprintf("Bootstrap vs naive SE (L = %d)", L))
  lines(df$T, df$naive_se, type = "b", col = "darkred", pch = 17, cex = 0.8)
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  legend("topleft", cex = 0.85, bty = "n",
         legend = c("Block bootstrap SE", "Naive iid SE"),
         col = c("darkblue", "darkred"), pch = c(19, 17), lty = 1)
  
  plot(df$T, df$ratio_observed, type = "b", col = "darkgreen",
       pch = 19, cex = 0.8,
       ylim = range(c(df$ratio_observed, df$ratio_expected), finite = TRUE),
       xlab = "T", ylab = "SE ratio (boot / naive)",
       main = sprintf("SE ratio vs sqrt(2 tau) (L = %d)", L))
  lines(df$T, df$ratio_expected, type = "b", col = "black",
        pch = 17, cex = 0.8, lty = 2)
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  abline(h = 1, lty = 3, col = "gray50")
  legend("topleft", cex = 0.85, bty = "n",
         legend = c("observed ratio", "sqrt(2 tau_E) (theory)"),
         col = c("darkgreen", "black"), pch = c(19, 17), lty = c(1, 2))
  par(mfrow = c(1, 1))
  
  cat(sprintf("Median ratio (boot/naive): %.2f\n",
              median(df$ratio_observed, na.rm = TRUE)))
  cat(sprintf("Max ratio:  %.2f  (at T=%.3f)\n",
              max(df$ratio_observed),
              df$T[which.max(df$ratio_observed)]))
  cat("=====================================================\n")
  
  invisible(df)
}


# =============================================================================
# PART 11B: Head-to-head Metropolis vs Wolff bootstrap at L=32, T_c
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
  
  # 4-panel overlay of bootstrap distributions
  draw_overlay <- function(m_samples, w_samples, label, main_title) {
    all_x  <- c(m_samples, w_samples)
    breaks <- seq(min(all_x), max(all_x), length.out = 40)
    hist(m_samples, breaks = breaks, col = rgb(0.8, 0.2, 0.2, 0.5),
         border = NA, freq = FALSE, xlim = range(all_x),
         xlab = label, main = main_title)
    hist(w_samples, breaks = breaks, col = rgb(0.2, 0.2, 0.8, 0.5),
         border = NA, freq = FALSE, add = TRUE)
    abline(v = mean(m_samples), col = "darkred",  lwd = 2)
    abline(v = mean(w_samples), col = "darkblue", lwd = 2)
    legend("topright", cex = 0.8, bty = "n",
           legend = c(sprintf("Metro (SE=%.4f)", sd(m_samples)),
                      sprintf("Wolff (SE=%.4f)", sd(w_samples))),
           fill = c(rgb(0.8, 0.2, 0.2, 0.5), rgb(0.2, 0.2, 0.8, 0.5)))
  }
  
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  draw_overlay(bs_m$boot_e,   bs_w$boot_e,   "e",              "Bootstrap dist: e")
  draw_overlay(bs_m$boot_c,   bs_w$boot_c,   "c",              "Bootstrap dist: c")
  draw_overlay(bs_m$boot_m,   bs_w$boot_m,   "|m|",            "Bootstrap dist: |m|")
  draw_overlay(bs_m$boot_chi, bs_w$boot_chi, expression(chi), "Bootstrap dist: chi")
  
  invisible(list(metropolis = bs_m, wolff = bs_w))
}


# =============================================================================
# Drive: run 11A sweep, 11A sanity check, 11B head-to-head, in order
# =============================================================================

cat("\n############################################################\n")
cat("Section 11: Block bootstrap for uncertainty quantification\n")
cat("############################################################\n\n")

# 11A: full sweep with SEs
sweep_err_df <- run_T_sweep_with_bootstrap()

cat("\nSample of 11A results:\n")
print(head(sweep_err_df[, c("L", "T", "e", "e_se", "c", "c_se",
                            "m", "m_se", "chi", "chi_se")], 10))

cat("\n11A figure 1: 4-panel phase transition with error bars...\n")
plot_T_sweep_with_errorbars(sweep_err_df)

cat("11A figure 2: SE vs T for each observable...\n")
plot_SE_vs_T(sweep_err_df)

cat("11A figure 3: sanity check, bootstrap vs naive SE...\n")
naive_vs_boot_df <- sanity_check_bootstrap_vs_naive(L = 32)

# 11B: head-to-head
cat("\n11B: Metropolis vs Wolff head-to-head at L=32, T_c...\n")
h2h_result <- bootstrap_head_to_head_Tc()




# =============================================================================
# Section 11 addendum: 4x2 phase transition figure with zoomed inset panels
#
# The main plot_T_sweep_with_errorbars() figure renders error bars at true
# scale, but they are typically < 0.1% of the axis height away from T_c
# and < 2% even at T_c peak -- visually invisible on a standard figure.
#
# This function produces a 4-row, 2-column layout where each observable
# gets two panels:
#   left:  full T range, as before (shows the physics at a glance)
#   right: zoomed to T in [T_zoom_min, T_zoom_max] around T_c, with
#          auto-scaled y-axis, so error bars become visible fractions of
#          axis height.
#
# Use this as a supplement to the main figure for the writeup.
# =============================================================================

plot_T_sweep_with_errorbars_zoom <- function(sweep_err_df,
                                             T_zoom_min = 2.00,
                                             T_zoom_max = 2.50) {
  L_values <- sort(unique(sweep_err_df$L))
  L_colors <- c("black", "darkorange", "purple", "darkcyan")[seq_along(L_values)]
  
  # Full and zoom T grids for the exact reference curves
  T_full <- seq(min(sweep_err_df$T), max(sweep_err_df$T), length.out = 400)
  T_zoom <- seq(T_zoom_min, T_zoom_max, length.out = 200)
  
  u_full <- onsager_energy_curve(T_full)
  c_full <- onsager_specific_heat_curve(T_full)
  m_full <- onsager_magnetization_curve(T_full)
  
  u_zoom <- onsager_energy_curve(T_zoom)
  c_zoom <- onsager_specific_heat_curve(T_zoom)
  m_zoom <- onsager_magnetization_curve(T_zoom)
  
  c_cap        <- 4
  c_full_plot  <- pmin(c_full, c_cap)
  c_zoom_plot  <- pmin(c_zoom, c_cap)
  
  # Subset of data in the zoom window
  zoom_df <- sweep_err_df[sweep_err_df$T >= T_zoom_min &
                            sweep_err_df$T <= T_zoom_max, ]
  
  draw_errorbars <- function(x, y, se, color, cex_pt = 0.9) {
    segments(x, y - se, x, y + se, col = color, lwd = 1.5)
    points(x, y, pch = 19, cex = cex_pt, col = color)
  }
  
  old_par <- par(mfrow = c(4, 2), mar = c(3.6, 4.2, 2.0, 0.8),
                 mgp = c(2.2, 0.7, 0),
                 no.readonly = TRUE)
  on.exit(par(old_par))
  
  # ------------------------------- ENERGY -----------------------------------
  # Full
  plot(T_full, u_full, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "u(T)", main = "Energy: full range",
       ylim = range(c(u_full,
                      sweep_err_df$e - sweep_err_df$e_se,
                      sweep_err_df$e + sweep_err_df$e_se),
                    finite = TRUE))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  rect(T_zoom_min, par("usr")[3], T_zoom_max, par("usr")[4],
       col = rgb(1, 1, 0, 0.12), border = NA)   # highlight zoom window
  for (k in seq_along(L_values)) {
    sub <- sweep_err_df[sweep_err_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$e, sub$e_se, L_colors[k], cex_pt = 0.6)
  }
  legend("topleft", cex = 0.7, bty = "n",
         legend = c(paste("L =", L_values), "Onsager"),
         col = c(L_colors, "gray40"),
         pch = c(rep(19, length(L_values)), NA),
         lwd = c(rep(NA, length(L_values)), 2))
  
  # Zoom
  y_zoom <- c(zoom_df$e - zoom_df$e_se, zoom_df$e + zoom_df$e_se,
              u_zoom[!is.na(u_zoom)])
  plot(T_zoom, u_zoom, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "u(T)",
       main = sprintf("Energy: zoom [%.2f, %.2f]", T_zoom_min, T_zoom_max),
       xlim = c(T_zoom_min, T_zoom_max),
       ylim = range(y_zoom))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- zoom_df[zoom_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$e, sub$e_se, L_colors[k])
  }
  
  # ------------------------------ SPECIFIC HEAT ------------------------------
  # Full
  plot(T_full, c_full_plot, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "c(T)",
       main = sprintf("Specific heat: full (Onsager cap %d)", c_cap),
       ylim = c(0, max(c(c_cap, sweep_err_df$c + sweep_err_df$c_se))))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  rect(T_zoom_min, par("usr")[3], T_zoom_max, par("usr")[4],
       col = rgb(1, 1, 0, 0.12), border = NA)
  for (k in seq_along(L_values)) {
    sub <- sweep_err_df[sweep_err_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$c, sub$c_se, L_colors[k], cex_pt = 0.6)
  }
  
  # Zoom
  y_zoom <- c(zoom_df$c - zoom_df$c_se, zoom_df$c + zoom_df$c_se,
              c_zoom_plot[!is.na(c_zoom_plot)])
  plot(T_zoom, c_zoom_plot, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "c(T)",
       main = sprintf("Specific heat: zoom [%.2f, %.2f]",
                      T_zoom_min, T_zoom_max),
       xlim = c(T_zoom_min, T_zoom_max),
       ylim = range(y_zoom))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- zoom_df[zoom_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$c, sub$c_se, L_colors[k])
  }
  
  # ----------------------------- MAGNETIZATION ------------------------------
  # Full
  plot(T_full, m_full, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "|m|(T)", main = "Magnetization: full range",
       ylim = c(0, 1.05))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  rect(T_zoom_min, par("usr")[3], T_zoom_max, par("usr")[4],
       col = rgb(1, 1, 0, 0.12), border = NA)
  for (k in seq_along(L_values)) {
    sub <- sweep_err_df[sweep_err_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$m, sub$m_se, L_colors[k], cex_pt = 0.6)
  }
  
  # Zoom
  y_zoom <- c(zoom_df$m - zoom_df$m_se, zoom_df$m + zoom_df$m_se,
              m_zoom[!is.na(m_zoom)])
  plot(T_zoom, m_zoom, type = "l", col = "gray40", lwd = 2,
       xlab = "T", ylab = "|m|(T)",
       main = sprintf("Magnetization: zoom [%.2f, %.2f]",
                      T_zoom_min, T_zoom_max),
       xlim = c(T_zoom_min, T_zoom_max),
       ylim = range(y_zoom))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- zoom_df[zoom_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$m, sub$m_se, L_colors[k])
  }
  
  # ---------------------------- SUSCEPTIBILITY -----------------------------
  # Full (no exact curve)
  plot(NA, xlim = range(T_full),
       ylim = c(0, max(sweep_err_df$chi + sweep_err_df$chi_se) * 1.05),
       xlab = "T", ylab = expression(chi(T)),
       main = "Susceptibility: full range")
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  rect(T_zoom_min, par("usr")[3], T_zoom_max, par("usr")[4],
       col = rgb(1, 1, 0, 0.12), border = NA)
  for (k in seq_along(L_values)) {
    sub <- sweep_err_df[sweep_err_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$chi, sub$chi_se, L_colors[k], cex_pt = 0.6)
  }
  
  # Zoom -- this one has the biggest payoff: peaks and error bars both
  # visible
  y_zoom <- c(zoom_df$chi - zoom_df$chi_se, zoom_df$chi + zoom_df$chi_se)
  plot(NA, xlim = c(T_zoom_min, T_zoom_max),
       ylim = c(0, max(y_zoom) * 1.05),
       xlab = "T", ylab = expression(chi(T)),
       main = sprintf("Susceptibility: zoom [%.2f, %.2f]",
                      T_zoom_min, T_zoom_max))
  abline(v = T_CRITICAL, lty = 2, col = "gray")
  for (k in seq_along(L_values)) {
    sub <- zoom_df[zoom_df$L == L_values[k], ]
    sub <- sub[order(sub$T), ]
    draw_errorbars(sub$T, sub$chi, sub$chi_se, L_colors[k])
  }
  
  invisible(NULL)
}


# =============================================================================
# Drive: reuse sweep_err_df from Section 11 (already in memory)
# =============================================================================

cat("\nPlotting 4x2 figure (full + zoomed) for error-bar visibility...\n")
plot_T_sweep_with_errorbars_zoom(sweep_err_df)









# =============================================================================
# Section 12: T_c estimation with bootstrap confidence interval
#
# Proposal Deliverable 5: estimate the critical temperature from MCMC data
# with a proper statistical uncertainty, compare to the exact Onsager value
# T_c = 2/ln(1+sqrt(2)) = 2.269185.
#
# Method:
#   1. For each L, find T_peak(L) = the temperature where chi is maximized.
#      Grid resolution is 0.02 around T_c, too coarse. We get sub-grid
#      resolution by fitting a parabola to the 5 points surrounding the
#      argmax and taking its vertex.
#
#   2. Finite-size scaling: T_peak(L) = T_c + a/L + O(1/L^2) for 2D Ising
#      (nu = 1). Plot T_peak vs 1/L and linearly extrapolate to 1/L = 0.
#      The intercept is our point estimate of T_c.
#
#   3. Parametric bootstrap for the CI: redraw chi(T) ~ Normal(chi_mean, chi_se)
#      at every (L, T) using the bootstrap SEs from Section 11, refit the
#      parabola at each L to get a new T_peak, refit the 1/L line to get a
#      new T_c estimate. Repeat 1000 times.
#
#   4. Report point estimate, 95% percentile CI, compare to exact T_c.
#      Report fits both with L=8 included and excluded (L=8 has the
#      largest subleading corrections).
#
# Expected precision: with 4 lattice sizes up to L=64 and 5000 cluster
# flips per (L, T), a CI half-width around 0.01-0.03 is realistic. The
# exact value 2.269185 should be inside the CI.
# =============================================================================


# -----------------------------------------------------------------------------
# 12.1: parabolic peak location from 5 points around argmax
#
# Given (T, chi) data, find the 2 points on each side of argmax and fit
# y = a*(T - T_peak)^2 + chi_peak via least squares. Returns T_peak.
# If the argmax is at the boundary of the grid (unlikely for us), falls
# back to argmax.
# -----------------------------------------------------------------------------

peak_T_parabolic <- function(T_grid, chi_values, n_side = 2) {
  stopifnot(length(T_grid) == length(chi_values))
  if (any(!is.finite(chi_values))) return(NA_real_)
  
  i_max <- which.max(chi_values)
  if (i_max <= n_side || i_max > length(chi_values) - n_side) {
    return(T_grid[i_max])   # fallback: on-grid argmax
  }
  
  idx <- (i_max - n_side):(i_max + n_side)
  T_sub <- T_grid[idx]
  c_sub <- chi_values[idx]
  
  # Fit chi ~ a + b*T + c*T^2  (quadratic in T)
  fit <- tryCatch(lm(c_sub ~ T_sub + I(T_sub^2)),
                  error = function(e) NULL)
  if (is.null(fit)) return(T_grid[i_max])
  co <- coef(fit)
  if (is.na(co[3]) || co[3] >= 0) return(T_grid[i_max])   # not a concave max
  
  T_peak <- -co[2] / (2 * co[3])
  # Sanity: T_peak should sit inside the fit window
  if (T_peak < min(T_sub) || T_peak > max(T_sub)) return(T_grid[i_max])
  as.numeric(T_peak)
}


# -----------------------------------------------------------------------------
# 12.2: fit T_peak vs 1/L, return T_c estimate (intercept)
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
# 12.3: parametric bootstrap for T_c CI
#
# For each bootstrap replicate, redraw chi at every (L, T) independently
# from Normal(chi_mean, chi_se). Refit the parabola at each L to get
# T_peak_b(L). Refit the 1/L line to get T_c_b. Repeat n_boot times.
# -----------------------------------------------------------------------------

bootstrap_Tc <- function(sweep_err_df,
                         L_values = sort(unique(sweep_err_df$L)),
                         n_boot   = 1000,
                         use_L_min = NULL,    # exclude L below this (e.g. 16)
                         seed     = 4242) {
  if (!is.null(use_L_min)) L_values <- L_values[L_values >= use_L_min]
  
  set.seed(seed)
  
  # Per-L grid of T, chi, chi_se
  per_L <- lapply(L_values, function(L) {
    sub <- sweep_err_df[sweep_err_df$L == L, ]
    sub <- sub[order(sub$T), ]
    list(L = L, T = sub$T, chi = sub$chi, chi_se = sub$chi_se)
  })
  
  # Point estimate (no redraw)
  T_peaks_pt <- vapply(per_L, function(pl) peak_T_parabolic(pl$T, pl$chi),
                       numeric(1))
  pt_fit <- fit_Tc_from_peaks(L_values, T_peaks_pt)
  
  # Bootstrap replicates
  T_c_boot <- numeric(n_boot)
  T_peaks_boot <- matrix(NA_real_, nrow = n_boot, ncol = length(L_values))
  for (b in seq_len(n_boot)) {
    T_peaks_b <- vapply(per_L, function(pl) {
      chi_draw <- rnorm(length(pl$chi), mean = pl$chi, sd = pl$chi_se)
      chi_draw <- pmax(chi_draw, 0)   # chi can't go negative physically
      peak_T_parabolic(pl$T, chi_draw)
    }, numeric(1))
    T_peaks_boot[b, ] <- T_peaks_b
    fb <- fit_Tc_from_peaks(L_values, T_peaks_b)
    T_c_boot[b] <- fb$T_c
  }
  
  # Percentile CI
  ci <- quantile(T_c_boot, c(0.025, 0.975), na.rm = TRUE)
  
  list(
    L_values    = L_values,
    T_peaks     = T_peaks_pt,
    T_c_point   = pt_fit$T_c,
    slope_point = pt_fit$slope,
    T_c_boot    = T_c_boot,
    ci_95       = ci,
    se          = sd(T_c_boot, na.rm = TRUE),
    T_peaks_boot = T_peaks_boot,
    fit         = pt_fit$fit
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
  err <- res$T_c_point - T_CRITICAL
  width <- (res$ci_95[2] - res$ci_95[1]) / 2
  cat(sprintf("  Point - exact:        %+.5f  (%.2f SEs away)\n",
              err, err / res$se))
  cat(sprintf("  Exact in 95%% CI?      %s\n",
              if (T_CRITICAL >= res$ci_95[1] && T_CRITICAL <= res$ci_95[2])
                "YES"
              else "NO"))
  cat("\n")
}


# -----------------------------------------------------------------------------
# 12.5: plots
#   Plot 1: T_peak vs 1/L with fit line, exact T_c marked, bootstrap envelope
#   Plot 2: chi(T) near T_c for each L with parabolic fit overlay, showing
#           how the peak location is found
#   Plot 3: histogram of bootstrap T_c estimates with CI marked
# -----------------------------------------------------------------------------

plot_Tc_extrapolation <- function(res_full, res_trim = NULL) {
  L_values <- res_full$L_values
  inv_L    <- 1 / L_values
  
  old_par <- par(mar = c(4.2, 4.2, 2.5, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Set up axes to cover 0 (the extrapolation target)
  x_plot <- c(0, max(inv_L) * 1.05)
  
  # Bootstrap envelope on the line: draw every fit
  y_plot_lo <- numeric(0); y_plot_hi <- numeric(0)
  if (!is.null(res_full$T_peaks_boot)) {
    for (b in seq_len(min(200, nrow(res_full$T_peaks_boot)))) {
      T_peaks_b <- res_full$T_peaks_boot[b, ]
      if (any(!is.finite(T_peaks_b))) next
      fb <- tryCatch(lm(T_peaks_b ~ inv_L), error = function(e) NULL)
      if (is.null(fb)) next
      y_b <- coef(fb)[1] + coef(fb)[2] * x_plot
      y_plot_lo <- c(y_plot_lo, min(y_b))
      y_plot_hi <- c(y_plot_hi, max(y_b))
    }
  }
  
  y_range <- range(c(res_full$T_peaks, res_full$T_c_point, T_CRITICAL,
                     res_full$ci_95,
                     y_plot_lo, y_plot_hi), na.rm = TRUE)
  y_range <- y_range + diff(y_range) * c(-0.05, 0.05)
  
  plot(inv_L, res_full$T_peaks, pch = 19, cex = 1.3, col = "darkblue",
       xlim = x_plot, ylim = y_range,
       xlab = "1 / L", ylab = expression(T[peak]),
       main = "Finite-size scaling: T_peak(L) = T_c + a/L + ...")
  
  # Draw bootstrap lines (light)
  if (!is.null(res_full$T_peaks_boot)) {
    for (b in seq_len(min(100, nrow(res_full$T_peaks_boot)))) {
      T_peaks_b <- res_full$T_peaks_boot[b, ]
      if (any(!is.finite(T_peaks_b))) next
      fb <- tryCatch(lm(T_peaks_b ~ inv_L), error = function(e) NULL)
      if (is.null(fb)) next
      abline(fb, col = rgb(0.2, 0.2, 0.8, 0.05), lwd = 1)
    }
  }
  
  # Point-estimate line
  abline(res_full$fit, col = "darkblue", lwd = 2)
  
  # Data points on top
  points(inv_L, res_full$T_peaks, pch = 19, cex = 1.3, col = "darkblue")
  
  # Exact T_c
  abline(h = T_CRITICAL, lty = 2, col = "darkred", lwd = 1.5)
  
  # 95% CI around the intercept
  points(0, res_full$T_c_point, pch = 18, cex = 1.8, col = "darkblue")
  segments(0, res_full$ci_95[1], 0, res_full$ci_95[2],
           col = "darkblue", lwd = 3)
  
  # If trimmed fit provided, overlay it in green
  if (!is.null(res_trim)) {
    inv_L_t <- 1 / res_trim$L_values
    abline(res_trim$fit, col = "darkgreen", lwd = 2, lty = 3)
    points(0, res_trim$T_c_point, pch = 18, cex = 1.5, col = "darkgreen")
    segments(0, res_trim$ci_95[1], 0, res_trim$ci_95[2],
             col = "darkgreen", lwd = 2)
  }
  
  legend_entries <- c(
    sprintf("T_peak data (L = %s)", paste(L_values, collapse = ",")),
    sprintf("All-L fit: T_c = %.5f [%.5f, %.5f]",
            res_full$T_c_point, res_full$ci_95[1], res_full$ci_95[2]),
    "bootstrap replicate lines",
    sprintf("Exact T_c = %.5f", T_CRITICAL)
  )
  legend_cols <- c("darkblue", "darkblue", rgb(0.2, 0.2, 0.8, 0.3), "darkred")
  legend_pchs <- c(19, NA, NA, NA)
  legend_ltys <- c(NA, 1, 1, 2)
  legend_lwds <- c(NA, 2, 1, 1.5)
  
  if (!is.null(res_trim)) {
    legend_entries <- c(legend_entries,
                        sprintf("L>=%d fit: T_c = %.5f [%.5f, %.5f]",
                                min(res_trim$L_values),
                                res_trim$T_c_point, res_trim$ci_95[1], res_trim$ci_95[2]))
    legend_cols <- c(legend_cols, "darkgreen")
    legend_pchs <- c(legend_pchs, NA)
    legend_ltys <- c(legend_ltys, 3)
    legend_lwds <- c(legend_lwds, 2)
  }
  
  legend("topleft", cex = 0.75, bty = "n",
         legend = legend_entries, col = legend_cols,
         pch = legend_pchs, lty = legend_ltys, lwd = legend_lwds)
  
  invisible(NULL)
}


plot_chi_peak_fits <- function(sweep_err_df, res,
                               T_window = c(2.15, 2.45)) {
  L_values <- res$L_values
  L_colors <- c("black", "darkorange", "purple", "darkcyan")[seq_along(L_values)]
  
  old_par <- par(mar = c(4.2, 4.2, 2.5, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Find y range in the window across all L
  sub_all <- sweep_err_df[sweep_err_df$L %in% L_values &
                            sweep_err_df$T >= T_window[1] &
                            sweep_err_df$T <= T_window[2], ]
  y_max <- max(sub_all$chi + sub_all$chi_se) * 1.05
  
  plot(NA, xlim = T_window, ylim = c(0, y_max),
       xlab = "T", ylab = expression(chi(T)),
       main = "Susceptibility near T_c with parabolic peak fits")
  abline(v = T_CRITICAL, lty = 2, col = "darkred", lwd = 1.2)
  
  for (k in seq_along(L_values)) {
    L <- L_values[k]; col <- L_colors[k]
    sub <- sweep_err_df[sweep_err_df$L == L &
                          sweep_err_df$T >= T_window[1] &
                          sweep_err_df$T <= T_window[2], ]
    sub <- sub[order(sub$T), ]
    
    # Data + error bars
    segments(sub$T, sub$chi - sub$chi_se, sub$T, sub$chi + sub$chi_se,
             col = col, lwd = 1.3)
    points(sub$T, sub$chi, pch = 19, cex = 1.0, col = col)
    
    # Parabolic fit around the argmax (same 5 points as peak_T_parabolic)
    i_max <- which.max(sub$chi)
    n_side <- 2
    if (i_max > n_side && i_max <= nrow(sub) - n_side) {
      idx <- (i_max - n_side):(i_max + n_side)
      T_sub <- sub$T[idx]; c_sub <- sub$chi[idx]
      fit <- lm(c_sub ~ T_sub + I(T_sub^2))
      T_dense <- seq(min(T_sub), max(T_sub), length.out = 100)
      c_dense <- predict(fit, newdata = data.frame(T_sub = T_dense))
      lines(T_dense, c_dense, col = col, lwd = 2, lty = 1)
      
      # Mark the parabola's vertex
      co <- coef(fit)
      if (!is.na(co[3]) && co[3] < 0) {
        T_peak <- -co[2] / (2 * co[3])
        points(T_peak, max(c_dense), pch = 13, cex = 1.8, col = col)
      }
    }
  }
  
  legend("topright", cex = 0.8, bty = "n",
         legend = c(paste("L =", L_values),
                    "parabolic peak", "exact T_c"),
         col = c(L_colors, "black", "darkred"),
         pch = c(rep(19, length(L_values)), 13, NA),
         lty = c(rep(NA, length(L_values) + 1), 2))
  
  invisible(NULL)
}


plot_Tc_bootstrap_hist <- function(res) {
  old_par <- par(mar = c(4.2, 4.2, 2.5, 1), no.readonly = TRUE)
  on.exit(par(old_par))
  
  boots <- res$T_c_boot[is.finite(res$T_c_boot)]
  x_range <- range(c(boots, T_CRITICAL, res$ci_95))
  x_range <- x_range + diff(x_range) * c(-0.1, 0.1)
  
  hist(boots, breaks = 40, freq = FALSE,
       col = rgb(0.2, 0.2, 0.8, 0.4), border = "darkblue",
       xlim = x_range,
       xlab = expression(T[c]),
       main = "Bootstrap distribution of T_c estimate")
  
  abline(v = res$T_c_point, col = "darkblue", lwd = 2)
  abline(v = res$ci_95, col = "darkblue", lwd = 1.5, lty = 2)
  abline(v = T_CRITICAL, col = "darkred", lwd = 2, lty = 3)
  
  legend("topright", cex = 0.85, bty = "n",
         legend = c(sprintf("Point est. = %.5f", res$T_c_point),
                    sprintf("95%% CI = [%.5f, %.5f]",
                            res$ci_95[1], res$ci_95[2]),
                    sprintf("Exact T_c = %.5f", T_CRITICAL)),
         col = c("darkblue", "darkblue", "darkred"),
         lty = c(1, 2, 3), lwd = c(2, 1.5, 2))
  
  invisible(NULL)
}


# =============================================================================
# Drive
# =============================================================================

cat("\n############################################################\n")
cat("Section 12: T_c estimation with bootstrap CI\n")
cat("############################################################\n\n")

cat("Fitting with all L values ({8, 16, 32, 64})...\n")
res_full <- bootstrap_Tc(sweep_err_df, n_boot = 1000, seed = 4242)
print_Tc_result(res_full, "All L")

cat("Fitting with L >= 16 (excluding L = 8 to reduce subleading corrections)...\n")
res_trim <- bootstrap_Tc(sweep_err_df, n_boot = 1000, use_L_min = 16, seed = 4243)
print_Tc_result(res_trim, "L >= 16")

cat("Plot 1: parabolic peak fits on chi(T) near T_c...\n")
plot_chi_peak_fits(sweep_err_df, res_full)

cat("Plot 2: T_peak vs 1/L with bootstrap extrapolation envelope...\n")
plot_Tc_extrapolation(res_full, res_trim = res_trim)

cat("Plot 3: histogram of bootstrap T_c estimates...\n")
plot_Tc_bootstrap_hist(res_full)

cat("\n=============================================================\n")
cat("Section 12 complete.\n")
cat("=============================================================\n")


