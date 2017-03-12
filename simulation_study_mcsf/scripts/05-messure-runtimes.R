#!/usr/bin/Rscript --no-save --no-site-file

library(coala)
library(parallel)

n_loci <- 1000
locus_length <- 150

model <- coal_model(20, n_loci, locus_length * 1000, ploidy = 2) +
  feat_recombination(locus_length) +
  feat_mutation(locus_length) +
  feat_size_change(0.1, 1, time = 0.2) +
  feat_size_change(1, 1, time = 0.3) +
  sumstat_seg_sites("segsites")

runtimes <- do.call(rbind, mclapply(1:100, function(i) {
  set.seed(i)
  segsites <- simulate(model)$segsites
  stats <- list(iHS = sumstat_ihh(calc_ihs = TRUE),
                MCMF = sumstat_mcmf(),
                Omega = sumstat_omega(min_win = 100, max_win = 50000, grid = 1000))
  
  sapply(stats, function(stat) {
    system.time(
      temp <- calc_sumstats_from_data(model + stat, segsites, NULL)
    )[3]
  })
}, mc.cores = detectCores()))

print(colMeans(runtimes))

save(runtimes, file = "data/runtimes.Rda")
