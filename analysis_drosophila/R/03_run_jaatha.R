#!/usr/bin/Rscript --no-site-file

library(coala)
library(jaatha)

coala_model <- readRDS("cache/model.Rds")
data_sumstats <- readRDS("cache/data_sumstats.Rds")

set.seed(2575932)
options(error = dump.frames)

jaatha_model <- create_jaatha_model(coala_model, 
                                    jsfs_summary = "sums",
                                    mcmf_breaks = c(.5, .75, .9),
                                    test = TRUE)

jaatha_data <- create_jaatha_data(data_sumstats, jaatha_model)

jaatha_results <- jaatha(jaatha_model, jaatha_data, 3, 200,
                         init_method = "random",
                         zoom_in_steps = 5,
                         final_sim = 500,
                         cores = 24)

print(jaatha_results)

saveRDS(jaatha_results, file = "cache/jaatha_results.Rda")
