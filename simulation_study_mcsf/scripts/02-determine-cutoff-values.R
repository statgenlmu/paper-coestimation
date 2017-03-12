#!/usr/bin/Rscript --no-save --no-site-file

library(dplyr)
library(parallel)
options(stringsAsFactors = FALSE)

probs <- 0:200 / 200
cores <- detectCores()


cutoff_values <- do.call(rbind, lapply(list.files("data"), function(study) {
  print(study)
  do.call(rbind, mclapply(list.files(file.path("data", study)), function(simulation) {
    # Load simulation data
    sim_data <- NULL
    sim_data_neutr <- NULL
    load(file.path("data", study, simulation))
    stopifnot(!is.null(sim_data))
    stopifnot(!is.null(sim_data_neutr))
    
    sim_data$study <- study
    sim_data_neutr$study <- study
    
    # For tests:
    # Determine cutoff values from neutral simulations
    sim_data_neutr %>% 
      filter(stat %in% c("mcmf", "tajimas_d", "pi", "omega")) %>%
      group_by(model, study, seed, stat) -> sum_results_neutr_split
    cutoff_tests <- do.call(rbind, lapply(probs, function(prob) {
      sum_results_neutr_split %>%
        summarize(prob = prob,
                  cutoff_min = quantile(value_min, prob, na.rm = TRUE),
                  cutoff_max = quantile(value_max, 1 - prob, na.rm = TRUE))
    }))
    
    # For scans:
    # Use extrem values of distribution
    sim_data %>% 
      filter(stat %in% c("ihs")) %>%
      group_by(model, study, seed, stat) -> sum_results_split
    cutoff_scan <- do.call(rbind, lapply(probs, function(prob) {
      sum_results_split %>%
        summarize(prob = prob,
                  cutoff_min = quantile(value_min, prob, na.rm = TRUE),
                  cutoff_max = quantile(value_max, 1 - prob, na.rm = TRUE))
    }))
    
    rbind(cutoff_tests, cutoff_scan)
  }, mc.cores = cores))
}))

dir.create("cache", showWarnings = FALSE)
save(cutoff_values, file = "cache/cutoff_values.Rda")
