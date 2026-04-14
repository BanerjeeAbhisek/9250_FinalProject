# =============================================================================
# ising2d.R
# Comparative study of MCMC sampling strategies for the 2D Ising model.
# Authors: Abhisek Banerjee & Mengyan Jing
# Course 9250, Spring 2026
# =============================================================================
#
# This script is organized in sections. Run top to bottom, or step through
# section by section. Each section is independent of later sections but may
# depend on earlier ones.
#
# Section 1: Setup (packages, constants)
# Section 2: Lattice utilities (R)
# Section 3: Onsager exact solution (R)
# Section 4: Sanity checks
# Section 5: Metropolis sampler (pure R reference implementation)
# Section 6: Metropolis sampler (Rcpp, ~100x faster)
# Section 7: Autocorrelation and effective sample size
# (Section 8+ will add Wolff, bootstrap as we go)
# =============================================================================


# -----------------------------------------------------------------------------
# Section 1: Setup
# -----------------------------------------------------------------------------

# We use Rcpp for the performance-critical inner loops of Metropolis and Wolff
# (added in later sections). For now, only base R is needed.
suppressPackageStartupMessages({
  library(Rcpp)
})

# Throughout the project we set J = 1 and k_B = 1, following Newman & Barkema
# (1999, p. 51). Temperature is therefore measured in units of J / k_B.
J_COUPLING <- 1
K_BOLTZMANN <- 1

# The exact 2D Ising critical temperature on the square lattice. Onsager
# (1944) gives the implicit condition sinh(2J/kT_c) sinh(2J'/kT_c) = 1
# in the abstract on p. 117 and as Eq. (17) on p. 119, and tabulates the
# numerical value H_c = J/kT_c = (1/2) log cot(pi/8) = 0.4406868 in Eq. (121),
# p. 141. For the isotropic case J = J' = 1 the condition reduces to
# sinh(2/T_c) = 1, i.e. 2/T_c = log(1 + sqrt(2)), giving:
#   T_c = 2 / log(1 + sqrt(2)) ~ 2.269185
T_CRITICAL <- 2 / log(1 + sqrt(2))


# -----------------------------------------------------------------------------
# Section 2: Lattice utilities
# -----------------------------------------------------------------------------
#
# A configuration is stored as an L x L integer matrix with entries in {-1, +1}.
# Periodic boundary conditions are handled by modular indexing helpers below.
# These functions are written in pure R for clarity and for use in tests; the
# samplers themselves will use Rcpp inner loops for speed.
# -----------------------------------------------------------------------------


#' Initialize an L x L Ising lattice.
#'
#' @param L Lattice side length.
#' @param state One of "cold" (all +1), "hot" (random uniform on {-1, +1}),
#'   or "ground" (alias for "cold").
#' @return An integer matrix of dimensions L x L with entries in {-1, +1}.
init_lattice <- function(L, state = c("cold", "hot", "ground")) {
  state <- match.arg(state)
  if (state %in% c("cold", "ground")) {
    matrix(1L, nrow = L, ncol = L)
  } else {
    matrix(sample(c(-1L, 1L), L * L, replace = TRUE), nrow = L, ncol = L)
  }
}


#' Total energy of a configuration under the zero-field Ising Hamiltonian
#' H = -J * sum_{<ij>} s_i s_j with periodic boundary conditions.
#'
#' This is the slow but obviously-correct reference implementation: it sums
#' each nearest-neighbor bond exactly once by counting only the right and down
#' neighbors of each site. We will use it as ground truth in unit tests for
#' the Rcpp version.
#'
#' @param spins An L x L integer matrix in {-1, +1}.
#' @return Total energy E (a scalar).
total_energy <- function(spins) {
  L <- nrow(spins)
  stopifnot(ncol(spins) == L)
  # Right neighbor of column j is column (j mod L) + 1
  right <- spins[, c(2:L, 1)]
  # Down neighbor of row i is row (i mod L) + 1
  down  <- spins[c(2:L, 1), ]
  -J_COUPLING * sum(spins * right + spins * down)
}


#' Total magnetization (signed sum of spins).
total_magnetization <- function(spins) {
  sum(spins)
}


#' Energy change ΔE that would result from flipping spin (i, j).
#' Newman & Barkema, Eq. 3.10:
#'   ΔE = 2 * J * s_k * sum_{n.n. of k} s_n
#'
#' @param spins An L x L integer matrix in {-1, +1}.
#' @param i,j 1-based row and column indices of the spin to (hypothetically) flip.
#' @return The scalar energy change.
delta_energy <- function(spins, i, j) {
  L <- nrow(spins)
  # 1-based modular arithmetic for periodic boundaries
  ip <- if (i == L) 1L else i + 1L
  im <- if (i == 1L) L else i - 1L
  jp <- if (j == L) 1L else j + 1L
  jm <- if (j == 1L) L else j - 1L
  neighbor_sum <- spins[ip, j] + spins[im, j] + spins[i, jp] + spins[i, jm]
  2 * J_COUPLING * spins[i, j] * neighbor_sum
}


# -----------------------------------------------------------------------------
# Section 3: Onsager exact solution
# -----------------------------------------------------------------------------
#
# Onsager (1944) derived the exact internal energy per spin and specific heat
# per spin for the 2D square Ising model in the thermodynamic limit (L -> inf).
# The closed-form expressions involve complete elliptic integrals of the first
# kind. We use the following standard form (see Newman & Barkema 1999, eqs.
# 3.20--3.21, or any statistical mechanics textbook):
#
#   beta * J = 1 / T  (with k = J = 1)
#   k1 = 2 * sinh(2 beta J) / cosh(2 beta J)^2          (modulus)
#   K(k1) = complete elliptic integral of the first kind, modulus k1
#
#   u(T)   = -J * coth(2 beta J) * [1 + (2/pi) * (2 tanh(2 beta J)^2 - 1) * K(k1)]
#
#   c(T)   = (2/pi) * (beta J coth(2 beta J))^2 *
#            { 2 K(k1) - 2 E(k1)
#              - (1 - tanh(2 beta J)^2) *
#                [pi/2 + (2 tanh(2 beta J)^2 - 1) K(k1)] }
#
# where E(k) is the complete elliptic integral of the second kind.
#
# Spontaneous magnetization (for T < T_c, in the thermodynamic limit):
#   |m|(T) = (1 - sinh(2 beta J)^(-4))^(1/8)        for T < T_c
#   |m|(T) = 0                                       for T >= T_c
#
# This is the famous Yang (1952) result, also reproduced in Onsager.
# -----------------------------------------------------------------------------


#' Complete elliptic integral of the first kind K(k) using base R's `integrate`.
#' We use the modulus convention (not the parameter m = k^2). This is fine for
#' a single ground-truth call per temperature; we don't need it to be fast.
#'
#' Near k = 1 (which corresponds to T = T_c in the Onsager formulas), the
#' integrand develops a near-singularity at theta = pi/2 and the default
#' subdivision limit of 100 is not enough. We bump it to 2000 and use a
#' modest relative tolerance, which is plenty for ground-truth comparisons.
elliptic_K <- function(k) {
  if (abs(k) >= 1) return(Inf)  # K diverges at k = 1 (the critical point)
  integrand <- function(theta) 1 / sqrt(1 - k^2 * sin(theta)^2)
  integrate(integrand, 0, pi / 2,
            rel.tol = 1e-8, subdivisions = 2000L)$value
}

#' Complete elliptic integral of the second kind E(k).
elliptic_E <- function(k) {
  integrand <- function(theta) sqrt(1 - k^2 * sin(theta)^2)
  integrate(integrand, 0, pi / 2,
            rel.tol = 1e-8, subdivisions = 2000L)$value
}


#' Onsager exact internal energy per spin in the thermodynamic limit.
#'
#' @param T Temperature (scalar or vector).
#' @return Energy per spin u(T) (same length as T).
onsager_energy <- function(T) {
  vapply(T, function(Ti) {
    beta <- 1 / Ti
    bj <- beta * J_COUPLING
    s2 <- sinh(2 * bj)
    c2 <- cosh(2 * bj)
    k1 <- 2 * s2 / c2^2
    K1 <- elliptic_K(k1)
    t2 <- tanh(2 * bj)^2
    -J_COUPLING * (c2 / s2) * (1 + (2 / pi) * (2 * t2 - 1) * K1)
  }, numeric(1))
}


#' Onsager exact specific heat per spin in the thermodynamic limit.
#'
#' Source: Onsager (1944), Eq. (117), p. 140. The formula is
#'   C/Nk = (2/pi) * (H * coth(2H))^2 *
#'          [ 2 K(k1) - 2 E(k1) - (1 - k1'') * ( pi/2 + k1'' * K(k1) ) ]
#' where H = beta * J, k1 = 2 sinh(2H) / cosh(2H)^2, and the auxiliary
#' modulus is k1'' = 2 tanh^2(2H) - 1 (Onsager Eq. 114, p. 139).
#'
#' Note: 1 - k1'' = 2 - 2 tanh^2(2H) = 2 * sech^2(2H), NOT 1 - tanh^2(2H).
#' Getting this wrong gives values that are systematically too large.
onsager_specific_heat <- function(T) {
  vapply(T, function(Ti) {
    beta <- 1 / Ti
    bj <- beta * J_COUPLING
    s2 <- sinh(2 * bj)
    c2 <- cosh(2 * bj)
    k1 <- 2 * s2 / c2^2
    K1 <- elliptic_K(k1)
    E1 <- elliptic_E(k1)
    k1pp <- 2 * tanh(2 * bj)^2 - 1   # Onsager Eq. (114): k_1''
    coth2 <- c2 / s2
    prefactor <- (2 / pi) * (bj * coth2)^2
    prefactor * (
      2 * K1 - 2 * E1
      - (1 - k1pp) * (pi / 2 + k1pp * K1)
    )
  }, numeric(1))
}


#' Yang (1952) exact spontaneous magnetization per spin.
#' Returns 0 for T >= T_c.
onsager_magnetization <- function(T) {
  vapply(T, function(Ti) {
    if (Ti >= T_CRITICAL) return(0)
    beta <- 1 / Ti
    (1 - sinh(2 * beta * J_COUPLING)^(-4))^(1 / 8)
  }, numeric(1))
}


# -----------------------------------------------------------------------------
# Section 4: Sanity checks
# -----------------------------------------------------------------------------
#
# Quick correctness checks. Run interactively when you change anything in the
# sections above. These will become formal testthat tests later if we want;
# for now they just print PASS/FAIL.
# -----------------------------------------------------------------------------

run_sanity_checks <- function() {
  cat("=== Sanity checks ===\n")
  
  # ---- Lattice utilities ----
  
  # 1. Cold lattice has all spins +1, energy = -2 * J * N (each site has 2
  #    bonds counted: right and down, both contribute -J)
  L <- 8
  s <- init_lattice(L, "cold")
  E_cold <- total_energy(s)
  expected <- -2 * J_COUPLING * L^2
  cat(sprintf("  Cold lattice energy: got %g, expected %g  -- %s\n",
              E_cold, expected, if (E_cold == expected) "PASS" else "FAIL"))
  
  # 2. Cold lattice magnetization = +N
  M_cold <- total_magnetization(s)
  cat(sprintf("  Cold lattice magnetization: got %d, expected %d  -- %s\n",
              M_cold, L^2, if (M_cold == L^2) "PASS" else "FAIL"))
  
  # 3. Flipping a single spin in the cold lattice should change energy by
  #    +8J (each of 4 neighbors contributes 2J).
  dE <- delta_energy(s, 4, 4)
  cat(sprintf("  ΔE flipping one spin in cold lattice: got %g, expected 8  -- %s\n",
              dE, if (dE == 8) "PASS" else "FAIL"))
  
  # 4. Flipping the spin and recomputing the total energy should agree with
  #    E_old + ΔE. This is the key invariant the Metropolis sampler relies on.
  s2 <- s; s2[4, 4] <- -s2[4, 4]
  E_new <- total_energy(s2)
  cat(sprintf("  Incremental vs. full energy after flip: %g vs %g  -- %s\n",
              E_cold + dE, E_new,
              if (isTRUE(all.equal(E_cold + dE, E_new))) "PASS" else "FAIL"))
  
  # 5. Periodic boundaries: flipping a corner spin should give the same ΔE
  #    as flipping any other spin in the cold lattice (by symmetry).
  dE_corner <- delta_energy(s, 1, 1)
  cat(sprintf("  ΔE at corner (PBC check): got %g, expected 8  -- %s\n",
              dE_corner, if (dE_corner == 8) "PASS" else "FAIL"))
  
  # ---- Onsager exact solution ----
  
  # 6. T_c value
  cat(sprintf("  T_c = %.6f (Onsager)\n", T_CRITICAL))
  
  # 7. Energy at T_c should be u(T_c) = -J * sqrt(2) ~ -1.41421
  #    (a known closed-form value at criticality). We evaluate slightly off
  #    T_c to avoid the 0 * Inf indeterminate form at exactly T_c (the
  #    elliptic integral diverges while the prefactor vanishes; the limit
  #    is finite but the formula is hard to evaluate at the singular point).
  u_Tc <- onsager_energy(T_CRITICAL - 1e-3)
  cat(sprintf("  u(T_c - 1e-3): got %.6f, expected ~%.6f  -- %s\n",
              u_Tc, -sqrt(2),
              if (isTRUE(all.equal(u_Tc, -sqrt(2), tolerance = 1e-2))) "PASS" else "FAIL"))
  
  # 8. Specific heat: cross-check against Onsager's own asymptotic form near
  #    T_c, his Eq. (120) on p. 140:
  #      C/Nk ~ (2/pi) * (log cot(pi/8))^2 * (K(k_1) - 1 - pi/4)
  #    valid as T -> T_c. We evaluate at T_c - 0.005, where the asymptotic
  #    form should agree with the full Eq. (117) to ~3-4 decimal places.
  T_near <- T_CRITICAL - 0.005
  c_full <- onsager_specific_heat(T_near)
  beta_n <- 1 / T_near
  bj_n <- beta_n * J_COUPLING
  k1_n <- 2 * sinh(2 * bj_n) / cosh(2 * bj_n)^2
  K1_n <- elliptic_K(k1_n)
  c_asym <- (2 / pi) * log(1 / tan(pi / 8))^2 * (K1_n - 1 - pi / 4)
  cat(sprintf("  c(T_c - 0.005): full %.4f vs asymptotic %.4f  -- %s\n",
              c_full, c_asym,
              if (isTRUE(all.equal(c_full, c_asym, tolerance = 1e-2))) "PASS" else "FAIL"))
  
  # 9. Magnetization at T = 1 should be very close to 1 (deep in ordered phase)
  m_low <- onsager_magnetization(1.0)
  cat(sprintf("  |m|(T=1): got %.6f  (should be ~1)  -- %s\n",
              m_low, if (m_low > 0.999) "PASS" else "FAIL"))
  
  # 10. Magnetization above T_c is exactly 0
  m_high <- onsager_magnetization(3.0)
  cat(sprintf("  |m|(T=3): got %.6f  (should be 0)  -- %s\n",
              m_high, if (m_high == 0) "PASS" else "FAIL"))
  
  cat("=====================\n")
}

# Run the checks now if the script is executed top-to-bottom.
run_sanity_checks()


# -----------------------------------------------------------------------------
# Section 5: Metropolis sampler (pure R reference implementation)
# -----------------------------------------------------------------------------
#
# Random-site, single-spin-flip Metropolis algorithm for the 2D Ising model
# with periodic boundary conditions, following Newman & Barkema (1999),
# Chapter 3, especially Eqs. (3.7) and (3.10).
#
# At each step:
#   1. Pick a site (i, j) uniformly at random.
#   2. Compute dE if we were to flip s_{i,j} (Newman & Barkema Eq. 3.10).
#   3. Accept the flip with probability A = min(1, exp(-dE / T)) (Eq. 3.7).
#      Note: when dE <= 0 we have exp(-dE/T) >= 1, so A = 1: always accept.
#   4. Repeat.
#
# A "sweep" is N = L^2 such steps, the conventional unit of Monte Carlo time
# (Newman & Barkema p. 49). Observables are recorded once per sweep, never
# more often --- successive single-step states are too correlated to count
# as independent measurements anyway.
#
# This implementation is written for clarity, not speed. It will be
# painfully slow for L > 32. Its purpose is to serve as the reference
# implementation that the Rcpp version (Section 6, when we add it) will be
# tested against. Both versions, run from the same RNG seed, MUST produce
# bit-identical trajectories.
# -----------------------------------------------------------------------------


#' One Metropolis step: pick a random site, propose a flip, accept or reject.
#'
#' Modifies `spins` and returns the (possibly updated) matrix along with the
#' energy change `dE` (zero if rejected) and the magnetization change `dM`
#' (zero if rejected). The caller is responsible for keeping running totals
#' of E and M --- it would be wasteful to recompute them from scratch every
#' step. This is the pattern Newman & Barkema use throughout Chapter 3.
#'
#' @param spins An L x L integer matrix of +/-1 spins.
#' @param T Temperature (positive scalar).
#' @return A list with elements `spins`, `dE`, `dM`, and `accepted`.
metropolis_step <- function(spins, T) {
  L <- nrow(spins)
  # Pick a uniformly random site. sample.int(L, 1) is the standard R idiom
  # and is what Newman & Barkema's "choose a spin at random" amounts to.
  i <- sample.int(L, 1L)
  j <- sample.int(L, 1L)
  
  # Compute dE for the proposed flip using the incremental formula
  # (Newman & Barkema Eq. 3.10). delta_energy() handles periodic boundaries.
  dE <- delta_energy(spins, i, j)
  
  # Metropolis acceptance rule (Newman & Barkema Eq. 3.7):
  #   A = min(1, exp(-dE / T))
  # If dE <= 0 the flip lowers (or doesn't change) the energy, so we accept
  # unconditionally --- equivalent to the rule above since exp(-dE/T) >= 1
  # in that case. We special-case it here to avoid an unnecessary RNG draw
  # and exp() call, which is the standard optimization.
  if (dE <= 0 || runif(1) < exp(-dE / T)) {
    # Accept the flip.
    s_old <- spins[i, j]
    spins[i, j] <- -s_old
    # The magnetization change is the difference of new and old values:
    # if s_old was +1 it becomes -1, change = -2; if s_old was -1, change = +2.
    dM <- -2L * s_old
    list(spins = spins, dE = dE, dM = dM, accepted = TRUE)
  } else {
    # Reject. State unchanged.
    list(spins = spins, dE = 0, dM = 0L, accepted = FALSE)
  }
}


#' Perform one full Metropolis sweep: N = L^2 single-spin steps.
#'
#' Returns the updated lattice and the cumulative changes in E and M over
#' the sweep, plus the number of accepted flips (so we can monitor the
#' acceptance rate, which is a useful sanity check by itself).
#'
#' @param spins An L x L integer matrix.
#' @param T Temperature.
#' @return A list with `spins`, `dE_sweep`, `dM_sweep`, `n_accepted`.
metropolis_sweep <- function(spins, T) {
  L <- nrow(spins)
  N <- L * L
  dE_sweep <- 0
  dM_sweep <- 0L
  n_accepted <- 0L
  for (k in seq_len(N)) {
    step <- metropolis_step(spins, T)
    spins <- step$spins
    dE_sweep <- dE_sweep + step$dE
    dM_sweep <- dM_sweep + step$dM
    if (step$accepted) n_accepted <- n_accepted + 1L
  }
  list(
    spins      = spins,
    dE_sweep   = dE_sweep,
    dM_sweep   = dM_sweep,
    n_accepted = n_accepted
  )
}


#' Run a Metropolis chain at fixed temperature, recording observables.
#'
#' Initializes a lattice in the requested state, runs `n_burnin` sweeps to
#' let the chain equilibrate (no observables recorded), then runs `n_sweeps`
#' more sweeps recording the energy E and magnetization M after each sweep.
#'
#' Returns a data frame with one row per recorded sweep. Energy and
#' magnetization are kept as TOTALS (not per-spin), so that downstream
#' code can choose how to normalize. The acceptance rate is also returned
#' as a sanity-check diagnostic.
#'
#' @param L Lattice side length.
#' @param T Temperature.
#' @param n_sweeps Number of measurement sweeps.
#' @param n_burnin Number of burn-in sweeps (default: n_sweeps / 10).
#' @param init_state Initial configuration: "cold", "hot", or "ground".
#' @param seed Optional integer RNG seed for reproducibility.
#' @return A list with `series` (data frame of E, M per sweep),
#'         `acceptance_rate` (over the measurement phase),
#'         and `final_spins` (the configuration at the end of the run).
run_metropolis <- function(L, T, n_sweeps,
                           n_burnin = n_sweeps %/% 10,
                           init_state = "hot",
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  spins <- init_lattice(L, init_state)
  
  # Compute initial E and M from scratch. After this, we maintain them
  # incrementally using the dE / dM returned by each sweep, which is
  # essential for performance: a from-scratch energy computation is O(N),
  # and we'd be doing it n_burnin + n_sweeps times = thousands of times
  # if we recomputed every sweep.
  E_current <- total_energy(spins)
  M_current <- total_magnetization(spins)
  
  # ---- Burn-in phase: run sweeps but record nothing ----
  for (s in seq_len(n_burnin)) {
    sweep <- metropolis_sweep(spins, T)
    spins <- sweep$spins
    E_current <- E_current + sweep$dE_sweep
    M_current <- M_current + sweep$dM_sweep
  }
  
  # ---- Measurement phase: record E and M after every sweep ----
  E_series <- numeric(n_sweeps)
  M_series <- integer(n_sweeps)
  total_accepted <- 0L
  N <- L * L
  for (s in seq_len(n_sweeps)) {
    sweep <- metropolis_sweep(spins, T)
    spins <- sweep$spins
    E_current <- E_current + sweep$dE_sweep
    M_current <- M_current + sweep$dM_sweep
    E_series[s] <- E_current
    M_series[s] <- M_current
    total_accepted <- total_accepted + sweep$n_accepted
  }
  
  list(
    series = data.frame(
      sweep = seq_len(n_sweeps),
      E = E_series,
      M = M_series
    ),
    acceptance_rate = total_accepted / (n_sweeps * N),
    final_spins = spins,
    L = L,
    T = T,
    n_burnin = n_burnin,
    n_sweeps = n_sweeps
  )
}


#' Compute the four observables from a metropolis run, normalized per spin.
#'
#' Takes the output of run_metropolis() and computes:
#'   - mean energy per spin <e>
#'   - mean absolute magnetization per spin <|m|>
#'   - specific heat per spin c
#'   - magnetic susceptibility per spin chi
#' These are point estimates only --- error bars come later via block
#' bootstrap. Note that for chi we use |M| in the second moment as
#' discussed in the theory write-up, to handle the symmetry-sector issue
#' on finite lattices in the ordered phase.
#'
#' @param run Output of run_metropolis().
#' @return A named list of the four observables.
compute_observables <- function(run) {
  N <- run$L^2
  T <- run$T
  E <- run$series$E
  M <- run$series$M
  
  e_mean    <- mean(E) / N
  abs_m     <- mean(abs(M)) / N
  c_value   <- (mean(E^2) - mean(E)^2) / (N * T^2)
  chi_value <- (mean(M^2) - mean(abs(M))^2) / (N * T)
  
  list(
    e   = e_mean,
    m   = abs_m,
    c   = c_value,
    chi = chi_value,
    acceptance_rate = run$acceptance_rate
  )
}


# -----------------------------------------------------------------------------
# Section 6: Metropolis sampler (Rcpp, ~100x faster)
# -----------------------------------------------------------------------------
#
# Same algorithm as Section 5 (random-site, single-spin-flip Metropolis with
# periodic boundaries), but the entire run --- burn-in, measurement, all of
# it --- happens inside one C++ function. R's only job is to call it once
# with the parameters and process the returned time series.
#
# Two performance tricks beyond just "compile to C++":
#
#   1. The five possible exp(-dE/T) values are precomputed once at the top
#      of the function and stored in a small lookup table indexed by dE.
#      Newman & Barkema discuss this on p. 53. dE can only be in
#      {-8, -4, 0, +4, +8}, so there are only two positive values that need
#      exp() calls.
#
#   2. The lattice is stored as a flat int* array of length L*L (row-major)
#      rather than as a 2D structure, so neighbor lookups are integer
#      arithmetic. Periodic wraparound is handled with a single modulo per
#      neighbor (no branching).
#
# Random numbers come from R's RNG via R::unif_rand(), so set.seed() in R
# before calling this function controls the C++ stream. Don't expect
# bit-identical trajectories with the pure-R version (the order of RNG
# calls differs slightly), but DO expect statistical agreement to within
# Monte Carlo error.
# -----------------------------------------------------------------------------

Rcpp::cppFunction('
List run_metropolis_cpp_inner(int L, double T, int n_sweeps, int n_burnin,
                              IntegerVector init_spins) {
  int N = L * L;

  // Copy the initial lattice into a local int array (row-major).
  // We could operate on init_spins directly, but copying makes it explicit
  // that we are not mutating the R object.
  std::vector<int> spins(N);
  for (int k = 0; k < N; k++) spins[k] = init_spins[k];

  // Compute initial total energy and magnetization from scratch.
  // Use the right+down trick to count each bond once, with periodic
  // boundaries via modulo arithmetic.
  long long E_current = 0;
  long long M_current = 0;
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      int s = spins[i * L + j];
      int s_right = spins[i * L + ((j + 1) % L)];
      int s_down  = spins[((i + 1) % L) * L + j];
      E_current -= s * s_right + s * s_down;  // -J = -1
      M_current += s;
    }
  }

  // Precompute exp(-dE/T) lookup table.
  // dE can only be in {-8, -4, 0, +4, +8}. We index it by (dE + 8) / 4
  // which maps -8 -> 0, -4 -> 1, 0 -> 2, +4 -> 3, +8 -> 4.
  // For dE <= 0 we never actually consult the table (we accept
  // unconditionally), so only the +4 and +8 entries matter, but we fill
  // all five for safety/clarity.
  double exp_table[5];
  for (int idx = 0; idx < 5; idx++) {
    int dE = -8 + 4 * idx;
    exp_table[idx] = std::exp(-((double)dE) / T);
  }

  // Storage for the measurement-phase time series.
  NumericVector E_series(n_sweeps);
  NumericVector M_series(n_sweeps);

  long long total_accepted = 0;
  long long total_attempts = 0;

  // Pull the RNG state from R once at the start. This is the standard
  // Rcpp idiom; without it R::unif_rand() will not respect set.seed().
  GetRNGstate();

  int total_sweeps = n_burnin + n_sweeps;
  for (int sweep = 0; sweep < total_sweeps; sweep++) {
    // One sweep = N attempted single-spin flips at randomly chosen sites.
    for (int step = 0; step < N; step++) {
      // Pick a random site. floor(unif * L) gives 0..L-1 uniformly.
      int i = (int)(::unif_rand() * L);
      int j = (int)(::unif_rand() * L);
      // Guard against the (vanishingly rare) case unif_rand() == 1.0
      if (i == L) i = L - 1;
      if (j == L) j = L - 1;

      int s_k = spins[i * L + j];

      // Sum of four nearest neighbors with periodic boundaries.
      int ip = (i + 1) % L;
      int im = (i + L - 1) % L;
      int jp = (j + 1) % L;
      int jm = (j + L - 1) % L;
      int neighbor_sum = spins[ip * L + j]
                       + spins[im * L + j]
                       + spins[i  * L + jp]
                       + spins[i  * L + jm];

      int dE = 2 * s_k * neighbor_sum;  // J = 1

      // Metropolis accept/reject.
      bool accept;
      if (dE <= 0) {
        accept = true;
      } else {
        // dE is +4 or +8; look up exp(-dE/T).
        double p = exp_table[(dE + 8) / 4];
        accept = (::unif_rand() < p);
      }

      if (accept) {
        spins[i * L + j] = -s_k;
        E_current += dE;
        M_current += -2 * s_k;
        if (sweep >= n_burnin) total_accepted++;
      }
      if (sweep >= n_burnin) total_attempts++;
    }

    // Record observables only during the measurement phase.
    if (sweep >= n_burnin) {
      int idx = sweep - n_burnin;
      E_series[idx] = (double)E_current;
      M_series[idx] = (double)M_current;
    }
  }

  PutRNGstate();

  // Return the time series, the final lattice, and the acceptance rate.
  // The final lattice is returned as an IntegerVector of length N (the
  // R caller can reshape it back to L x L if it wants).
  IntegerVector final_spins(N);
  for (int k = 0; k < N; k++) final_spins[k] = spins[k];

  return List::create(
    Named("E_series")        = E_series,
    Named("M_series")        = M_series,
    Named("final_spins")     = final_spins,
    Named("acceptance_rate") = (double)total_accepted / (double)total_attempts
  );
}
')


#' Run a Metropolis chain at fixed temperature using the Rcpp inner loop.
#'
#' Drop-in replacement for `run_metropolis()` that runs ~100x faster.
#' Same arguments, same return value structure.
#'
#' @param L Lattice side length.
#' @param T Temperature.
#' @param n_sweeps Number of measurement sweeps.
#' @param n_burnin Number of burn-in sweeps (default: n_sweeps / 10).
#' @param init_state Initial configuration: "cold", "hot", or "ground".
#' @param seed Optional integer RNG seed for reproducibility.
#' @return Same structure as `run_metropolis()`: a list with `series`,
#'         `acceptance_rate`, `final_spins`, etc.
run_metropolis_fast <- function(L, T, n_sweeps,
                                n_burnin = n_sweeps %/% 10,
                                init_state = "hot",
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  spins0 <- init_lattice(L, init_state)
  
  # Pass the lattice to C++ as a flat integer vector (row-major order).
  # R stores matrices in column-major, but for Ising on a square lattice
  # with PBC the choice doesn't matter --- the lattice is symmetric under
  # row/column swap. We use as.integer(t(spins0)) to get row-major if you
  # ever care; here we just use as.integer(spins0) and treat C++'s "row i,
  # column j" as "R's column i, row j". Either is fine and gives the
  # same physics.
  init_flat <- as.integer(spins0)
  
  result <- run_metropolis_cpp_inner(
    L = as.integer(L),
    T = as.numeric(T),
    n_sweeps = as.integer(n_sweeps),
    n_burnin = as.integer(n_burnin),
    init_spins = init_flat
  )
  
  # Reshape the final lattice back to L x L.
  final_spins <- matrix(result$final_spins, nrow = L, ncol = L)
  
  list(
    series = data.frame(
      sweep = seq_len(n_sweeps),
      E = result$E_series,
      M = result$M_series
    ),
    acceptance_rate = result$acceptance_rate,
    final_spins = final_spins,
    L = L,
    T = T,
    n_burnin = n_burnin,
    n_sweeps = n_sweeps
  )
}


# -----------------------------------------------------------------------------
# Section 7: Autocorrelation and effective sample size
# -----------------------------------------------------------------------------
#
# An MCMC chain produces samples that are autocorrelated --- successive draws
# are not independent. To honestly report Monte Carlo error bars, we need to
# know HOW correlated they are, which is captured by the integrated
# autocorrelation time tau_int.
#
# Definitions (Newman & Barkema, Sec. 3.3, pp. 60-65):
#
#   The normalized autocorrelation function of a stationary time series x_t is
#       rho(t) = (E[x_s x_{s+t}] - E[x]^2) / (E[x^2] - E[x]^2)
#   so rho(0) = 1 and rho(t) -> 0 as t -> infinity.
#
#   The integrated autocorrelation time is
#       tau_int = 1/2 + sum_{t=1}^{infinity} rho(t).
#   In practice we cannot sum to infinity --- the tail of rho(t) is dominated
#   by noise --- so we need a windowing rule.
#
#   The effective sample size after M correlated samples is
#       ESS = M / (2 * tau_int).
#   Roughly: every 2*tau_int sweeps gives you one "effectively independent"
#   measurement. The standard error of the sample mean is then
#       SE(mean) ~ sigma_x / sqrt(ESS),
#   not sigma_x / sqrt(M).
#
# We use Sokal's automatic windowing rule: truncate the sum at the smallest W
# satisfying W >= c * tau_int(W), with c = 5 a standard choice. The intuition:
# for an exponentially decaying rho(t) with decay constant tau, the cumulative
# sum reaches its asymptotic value by t ~ 5*tau, and beyond that point we are
# just adding noise to a converged sum.
# -----------------------------------------------------------------------------


#' Normalized autocorrelation function of a time series.
#'
#' @param x Numeric vector (the time series).
#' @param max_lag Maximum lag to compute. Default min(length(x)/4, 1000).
#' @return Numeric vector of length max_lag + 1, with rho[1] = 1 (lag 0).
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


#' Integrated autocorrelation time with Sokal's automatic windowing.
#'
#' @param x Numeric vector.
#' @param c Sokal cutoff parameter (default 5).
#' @param max_lag Maximum lag to consider. Default length(x) / 4.
#' @return A list with `tau`, `window`, `rho`.
integrated_tau <- function(x, c = 5, max_lag = NULL) {
  n <- length(x)
  if (is.null(max_lag)) max_lag <- n %/% 4
  
  rho <- autocorrelation(x, max_lag = max_lag)
  
  # Cumulative sum: tau_partial[W] = 1/2 + sum_{t=1}^{W} rho(t)
  tau_partial <- 0.5 + cumsum(rho[-1])
  
  # Find smallest W such that W >= c * tau_partial[W]
  W_grid <- seq_along(tau_partial)
  ok <- W_grid >= c * tau_partial
  if (any(ok)) {
    W <- which(ok)[1]
  } else {
    W <- max_lag
    warning(sprintf(
      "Sokal windowing not satisfied within max_lag = %d. Chain may be too short. Returned tau = %.2f at W = %d.",
      max_lag, tau_partial[W], W))
  }
  
  list(tau = tau_partial[W], window = W, rho = rho)
}


#' Effective sample size of a correlated time series.
#' ESS = M / (2 * tau_int).
effective_sample_size <- function(x, ...) {
  tau <- integrated_tau(x, ...)$tau
  length(x) / (2 * tau)
}


#' Convenience: tau_int and ESS for both energy and magnetization series of
#' a metropolis run. Returns a one-row data frame with diagnostics.
chain_diagnostics <- function(run) {
  tau_E <- integrated_tau(run$series$E)
  tau_M <- integrated_tau(run$series$M)
  M <- length(run$series$E)
  data.frame(
    L         = run$L,
    T         = run$T,
    n_sweeps  = M,
    tau_E     = tau_E$tau,
    tau_M     = tau_M$tau,
    ESS_E     = M / (2 * tau_E$tau),
    ESS_M     = M / (2 * tau_M$tau),
    window_E  = tau_E$window,
    window_M  = tau_M$window,
    acc_rate  = run$acceptance_rate
  )
}