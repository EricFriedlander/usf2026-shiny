# mutation_simulation.R
# Steady-state distribution and fit score for the mutation-rate estimation
# exercise (Tab 3).
#
# The stationary distribution of the WF model with symmetric mutation rate μ
# is exactly Beta(4Nμ, 4Nμ).  For large N (or any N given enough time),
# sampling from this Beta IS sampling from steady state.


#' Sample steady-state allele frequencies
#'
#' Returns n_pops draws from the exact stationary distribution Beta(4Nμ, 4Nμ),
#' which is the limiting allele-frequency distribution of the Wright-Fisher
#' model with symmetric mutation rate mu.
#'
#' @param N      Numeric. Effective population size.
#' @param mu     Numeric. Per-allele per-generation mutation rate (symmetric).
#' @param n_pops Integer. Number of independent populations to sample.
#' @param seed   Integer or NULL. Random seed.
#'
#' @return Numeric vector of length n_pops with frequencies in [0, 1].
simulate_steady_state <- function(N, mu, n_pops, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  theta <- 4 * N * mu
  if (!is.finite(theta) || theta <= 0) return(rep(0.5, n_pops))
  rbeta(n_pops, theta, theta)
}


#' Cosine-similarity fit score between empirical and theoretical distributions
#'
#' @param emp_freqs  Numeric vector of empirical allele frequencies.
#' @param N          Numeric. Effective population size (used internally).
#' @param mu         Numeric. Mutation rate to evaluate.
#' @param n_bins     Integer. Number of histogram bins for comparison.
#'
#' @return Numeric score in [0, 1], higher = better fit.
mutation_fit_score <- function(emp_freqs, N, mu, n_bins = 30) {
  alpha <- 4 * N * mu
  if (!is.finite(alpha) || alpha <= 0) return(0)

  breaks   <- seq(0, 1, length.out = n_bins + 1)
  emp_dens <- hist(emp_freqs, breaks = breaks, plot = FALSE)$density

  mids      <- (breaks[-length(breaks)] + breaks[-1]) / 2
  theo_dens <- dbeta(mids, alpha, alpha)
  theo_dens[!is.finite(theo_dens)] <- 0

  denom <- sqrt(sum(emp_dens^2)) * sqrt(sum(theo_dens^2))
  if (!is.finite(denom) || denom == 0) return(0)

  round(max(0, min(1, sum(emp_dens * theo_dens) / denom)), 3)
}
