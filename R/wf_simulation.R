# wf_simulation.R
# Core simulation functions for the Wright-Fisher genetic drift model.
#
# Allele encoding:
#   0 = Mario (red)
#   1 = Luigi (green)


#' Initialize a generation
#'
#' Creates the starting population as a vector of 0s and 1s, where each
#' individual is assigned allele 1 (Luigi) with probability p0.
#'
#' @param N    Integer. Population size.
#' @param p0   Numeric in [0, 1]. Initial frequency of allele 1 (Luigi). Defaults to 0.5.
#'
#' @return Integer vector of length N with values 0 or 1.
init_generation <- function(N, p0 = 0.5) {
  stopifnot(is.numeric(N), N >= 1)
  stopifnot(is.numeric(p0), p0 >= 0, p0 <= 1)
  as.integer(rbinom(N, size = 1, prob = p0))
}


#' Sample one individual from the previous generation
#'
#' Performs a single Wright-Fisher draw: picks one individual at random
#' (with replacement) from the previous complete generation.
#'
#' @param prev_gen  Integer vector. The previous complete generation (0s and 1s).
#'
#' @return A single integer, 0 or 1.
sample_one_individual <- function(prev_gen) {
  stopifnot(length(prev_gen) >= 1)
  sample(prev_gen, 1, replace = TRUE)
}


#' Sample a full generation from the previous generation
#'
#' Fills an entire new generation of the same size as prev_gen by sampling
#' with replacement — the standard Wright-Fisher model.
#'
#' @param prev_gen  Integer vector. The previous complete generation (0s and 1s).
#'
#' @return Integer vector of the same length as prev_gen, with values 0 or 1.
sample_generation <- function(prev_gen) {
  stopifnot(length(prev_gen) >= 1)
  sample(prev_gen, length(prev_gen), replace = TRUE)
}
