#!/usr/bin/Rscript --no-site-file

library(coala)
snp_data <- readRDS("cache/snp_data.Rds")
locus_length <- readRDS("cache/locus_length.Rds")

message("Creating model...")
locus_list <- c(snp_data[[1]]$segsites, 
                snp_data[[2]]$segsites) 
class(locus_list) <- "list"
saveRDS(locus_list, "cache/locus_list.Rds")

sample_sizes <- sapply(snp_data[[1]]$populations, length)
message("Number of loci: ", length(locus_list))
n_snps <- sapply(locus_list, ncol)
message("Total number of SNPs: ", sum(n_snps))
message("Average number of SNPS per locus: ", mean(n_snps))
message("Standard deviation in number of SNPS: ", sqrt(var(n_snps)))

model <- coal_model(sample_sizes, length(locus_list), locus_length) +
  #feat_mutation(par_variation(par_range("theta", 10, 2000), var(n_snps))) +
  feat_mutation(par_range("theta", 10, 2000)) +
  feat_recombination(par_range("rho", 10, 6000)) +
  feat_pop_merge(par_range("tau", 0.11, 10), 2, 1) +
  feat_size_change(par_range("q_2", 0.01, 10), population = 2, time = 0) +
  feat_migration(par_range("m12", 0, 20), pop_from = 1, pop_to = 2, time = 0) +
  feat_migration(par_range("m21", 0, 20), pop_from = 2, pop_to = 1, time = 0) +
  feat_size_change(par_range("q_tau", 0.00001, 1), population = 2, time = par_expr(tau - 0.1)) +
  feat_migration(0.0, symmetric = TRUE, time = par_expr(tau - 0.1)) +
  feat_ignore_singletons() +
  
  # Selection Features
  par_range("s", 0, .25) +
  feat_selection(strength_A = par_zero_inflation(1000, par_expr(1 - s)),
                 population = 2, time = 0.02) +
  sumstat_jsfs() +
  sumstat_four_gamete("fgc1", population = 1) +
  sumstat_mcmf(population = 2)

# Additionally calculate and print pi as sanity check
model2 <- model +
  sumstat_nucleotide_div("pi_1", population = 1) +
  sumstat_nucleotide_div("pi_2", population = 2)
message("")

message("Calculating summary statistics...")
data_sumstats <- calc_sumstats_from_data(model2, locus_list)
message("Average pi in african population: ", round(mean(data_sumstats$pi_1)))
message("Average pi in european population: ", round(mean(data_sumstats$pi_2)))
print(data_sumstats$jsfs)

saveRDS(model, "cache/model.Rds")
saveRDS(data_sumstats, "cache/data_sumstats.Rds")
