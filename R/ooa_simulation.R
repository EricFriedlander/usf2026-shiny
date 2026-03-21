# ooa_simulation.R
# Out of Africa demographic model and Allele Frequency Spectrum simulation.
#
# Demographic model (going forward in time):
#   Phase 1: Ancestral population Ne=10,000 (before t_bottle)
#   Phase 2: Bottleneck: Africa=10,000, non-Africa shrinks to bottleneck_size
#             Duration: 800 generations (t_split to t_bottle = t_split + 800)
#   Phase 3: After split (t_split to present): Africa=10,000,
#             non-Africa grows exponentially at growth_rate% per generation
#
# True hidden parameters: t_split=2400, bottleneck_size=500, growth_rate=2.0


#' Compute Ne history for Africa and non-Africa
#'
#' Returns a data frame with columns: gen_ago, ne_africa, ne_nonafrica, phase
#' gen_ago goes from 0 (present) to max_time.
#'
#' @param t_split       Generations ago of Africa/non-Africa split
#' @param bottleneck_size   Ne of non-Africa during bottleneck
#' @param growth_rate   Post-split exponential growth rate (% per generation)
#' @param max_time      How far back to compute (default: t_split + 1200)
#'
#' @return data.frame with ne history
compute_ne_history <- function(t_split, bottleneck_size, growth_rate,
                                max_time = NULL) {
  if (is.null(max_time)) max_time <- t_split + 1200
  t_bottle <- t_split + 800

  gen_ago    <- seq(0, max_time, by = 50)
  ne_africa  <- rep(10000, length(gen_ago))
  ne_nonafrica <- numeric(length(gen_ago))

  r <- growth_rate / 100  # per-generation rate

  for (i in seq_along(gen_ago)) {
    t <- gen_ago[i]

    if (t <= t_split) {
      # Post-split: exponential growth backward (shrinking going back)
      # Ne(t) = Ne_present / (1 + r)^t, but Ne_present determined by growth from t_split
      # Ne at split = bottleneck_size * (1 + r)^0 = bottleneck_size (just exited bottleneck)
      # Ne at present = bottleneck_size * (1 + r)^t_split
      # Ne at time t ago = bottleneck_size * (1 + r)^(t_split - t)
      ne_nonafrica[i] <- bottleneck_size * (1 + r)^(t_split - t)
    } else if (t <= t_bottle) {
      # In bottleneck
      ne_nonafrica[i] <- bottleneck_size
    } else {
      # Ancestral (merged back): single population
      ne_nonafrica[i] <- 10000
    }
  }

  data.frame(
    gen_ago      = gen_ago,
    ne_africa    = ne_africa,
    ne_nonafrica = ne_nonafrica,
    phase        = ifelse(gen_ago > t_bottle, "Ancestral",
                   ifelse(gen_ago > t_split,  "Bottleneck",
                                              "Post-split growth"))
  )
}


#' Compute display dot count (scaled from true Ne)
#'
#' Maps Ne to a visible number of dots (max ~50).
#'
#' @param ne        True effective population size
#' @param ne_max    Maximum Ne in the scenario (for scaling)
#' @param max_dots  Maximum dots to display
#'
#' @return Integer number of dots to display
ne_to_dots <- function(ne, ne_max = 10000, max_dots = 50) {
  pmax(1L, round(ne / ne_max * max_dots))
}


#' Simulate Allele Frequency Spectrum (AFS)
#'
#' Uses a fast Beta-distribution approximation to the AFS under different
#' demographic histories. For each of n_snps SNPs, draws an allele frequency
#' and then samples n_sample chromosomes.
#'
#' @param t_split        Generations ago of split
#' @param bottleneck_size  Ne during non-Africa bottleneck
#' @param growth_rate    Post-split growth rate (% per generation)
#' @param n_snps         Number of SNPs to simulate (default 150)
#' @param n_sample       Number of chromosomes sampled per population (default 50)
#' @param seed           Random seed (NULL for no seed)
#'
#' @return list with:
#'   africa_counts:    integer vector of derived allele counts (length n_snps)
#'   nonafrica_counts: integer vector of derived allele counts (length n_snps)
simulate_afs <- function(t_split, bottleneck_size, growth_rate,
                          n_snps = 150, n_sample = 50, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # --- Africa: large stable population → AFS shaped by drift/mutation balance
  # Approx: Beta(0.3, 0.3) gives U-shaped spectrum (drift dominates at large Ne)
  # but with many singletons → Beta(0.1, 1.5) approximates site-freq spectrum
  africa_alpha <- 0.15
  africa_beta  <- 1.5
  africa_freqs <- rbeta(n_snps, africa_alpha, africa_beta)
  africa_counts <- rbinom(n_snps, n_sample, africa_freqs)

  # --- Non-Africa: bottleneck + growth → more rare variants, some intermediate
  # Harmonic mean Ne approximates effective bottleneck strength
  r <- growth_rate / 100
  t_bottle <- t_split + 800

  # Ne_effective for non-Africa (harmonic mean over all phases)
  # Phase 1: bottleneck for 800 gens
  # Phase 2: growth phase for t_split gens
  gens_bottle <- 800
  gens_growth <- t_split
  ne_bottle   <- bottleneck_size
  # Average Ne during growth phase (geometric series)
  ne_growth_avg <- if (r > 0) {
    ne_bottle * ((1 + r)^gens_growth - 1) / (gens_growth * log(1 + r))
  } else {
    ne_bottle
  }

  total_gens  <- gens_bottle + gens_growth
  ne_harmonic <- total_gens / (gens_bottle / ne_bottle + gens_growth / ne_growth_avg)

  # Bottleneck severity: stronger bottleneck → more rare variants in derived pop
  # Map ne_harmonic to Beta parameters: smaller Ne → more extreme spectrum
  strength <- pmax(0.01, pmin(1, ne_harmonic / 10000))
  na_alpha <- 0.1 + 0.4 * strength    # ranges from 0.1 (severe) to 0.5 (mild)
  na_beta  <- 1.0 + 2.0 * strength    # ranges from 1.0 (severe) to 3.0 (mild)

  nonafrica_freqs  <- rbeta(n_snps, na_alpha, na_beta)
  nonafrica_counts <- rbinom(n_snps, n_sample, nonafrica_freqs)

  list(
    africa_counts    = africa_counts,
    nonafrica_counts = nonafrica_counts,
    n_sample         = n_sample,
    n_snps           = n_snps
  )
}


#' Compute goodness of fit between simulated and real AFS
#'
#' @param sim_counts   Integer vector of simulated derived allele counts
#' @param real_counts  Integer vector of real (hidden truth) counts
#' @param n_sample     Sample size (chromosome count)
#'
#' @return Numeric score in [0, 1], higher = better fit
afs_fit_score <- function(sim_counts, real_counts, n_sample) {
  bins      <- 0:n_sample
  sim_hist  <- tabulate(sim_counts  + 1L, nbins = n_sample + 1) / length(sim_counts)
  real_hist <- tabulate(real_counts + 1L, nbins = n_sample + 1) / length(real_counts)
  ssd       <- sum((sim_hist - real_hist)^2)
  score     <- 1 / (1 + ssd * length(sim_counts))
  round(score, 3)
}
