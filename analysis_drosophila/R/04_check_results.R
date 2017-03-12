library(jaatha)
library(coala)

model <- readRDS("cache/model.Rds")
data_sumstats <- readRDS("cache/data_sumstats.Rds")
results <- readRDS("cache/jaatha_results.Rda")

message("Real data:")
message("Average pi in african population: ", round(mean(data_sumstats$pi_1)))
message("Average pi in european population: ", round(mean(data_sumstats$pi_2)))

# Additionally calculate PI
model2 <- model +
  sumstat_nucleotide_div("pi_1", population = 1) +
  sumstat_nucleotide_div("pi_2", population = 2)

data_simulated <- simulate(model2, pars = results$estimate)
message("Simulated data:")
message("Average pi in african population: ", round(mean(data_simulated$pi_1)))
message("Average pi in european population: ", round(mean(data_simulated$pi_2)))

