#!/usr/bin/Rscript --no-save --no-site-file

library(coala)
library(dplyr)
library(zoo)

source("scripts/simulate_data.R")

# number of loci
n_loci <- 1000

# locus length in kb 
locus_length <- 150

# Repetetions for smoothing
reps <- 100


# ------------------------------------------------------------------------------
# Transformations for summary statistics
# ------------------------------------------------------------------------------

# Returns the most exterme iHS value that reached by a consequitve pair of SNPs 
# on the locus
get_max_ihs <- function(ihh) {
  ihs <- ihh$iHS
  t(vapply(1:n_loci, function(i) {
    x <- ihs$iHS[ihs$CHR == i]
    x <- x[is.finite(x)]
    if (length(x) < 2) return(c(min = 0, max = 0))
    c(min = min(rollmax(x, 2)), max = max(-rollmax(-x, 2)))
  }, numeric(2)))
}

# Returns the most exterme omega value per locus
get_max_omega <- function(omega) {
  (omega %>% 
     group_by(locus) %>% 
     summarise(omega = max(omega)))$omega
}


# ------------------------------------------------------------------------------
# The demographic model
# ------------------------------------------------------------------------------

model <- coal_model(20, n_loci, locus_length * 1000, ploidy = 2) +
  feat_recombination(locus_length) +
  feat_mutation(locus_length) +
  feat_size_change(0.1, 1, time = 0.2) +
  feat_size_change(1, 1, time = 0.3) +
  par_range("s", 0, 1) +
  feat_selection(strength_A = par_zero_inflation(1000, par_expr(1 - s), FALSE),
                 population = 1, time = 0.03) +
  sumstat_ihh(calc_ihs = TRUE, transformation = get_max_ihs) + 
  sumstat_mcmf() + 
  sumstat_tajimas_d() +
  sumstat_omega(min_win = 100, max_win = 50000, grid = 1000, 
                transformation = get_max_omega)


set.seed(238549)
simulate_data("normal", list(bn_fixed_rate_ = model), reps,  cores = detectCores())
simulate_data("unphased", list(bn_fixed_rate_ = model + feat_unphased(2)), reps, cores = detectCores())
