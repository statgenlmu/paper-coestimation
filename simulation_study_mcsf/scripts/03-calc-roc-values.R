#!/usr/bin/Rscript --no-save --no-site-file

library(dplyr)
library(parallel)

options(stringsAsFactors = FALSE)
cores <- detectCores()

cutoff_values <- NULL
load("cache/cutoff_values.Rda")
stopifnot(!is.null(cutoff_values))


na_to_false <- function(x) ifelse(is.na(x), FALSE, x)

roc_values <- do.call(rbind, lapply(list.files("data"), function(study) {
  print(study)
  do.call(rbind, mclapply(list.files(file.path("data", study)), function(simulation) {
    sim_data <- NULL
    load(file.path(file.path("data", study, simulation)))
    stopifnot(!is.null(sim_data))
    sim_data$study <- study
    
    sim_data %>% 
      left_join(cutoff_values, by = c("model", "stat", "seed", "study")) %>%
      group_by(model, stat, prob, sel, study) %>%
      summarise(classified = mean(na_to_false(value_min <= cutoff_min) | na_to_false(value_max > cutoff_max))) %>%
      group_by(model, stat, prob, study) %>%
      summarise(true_pos = classified[sel], false_pos = classified[!sel])
  }, mc.cores = cores))
}))

get_demography <- function(model) {
  ret <- rep(NA, length(model))
  ret <- ifelse(grepl("^const_", model), "const", ret)
  ret <- ifelse(grepl("^bn_", model), "bn", ret)
  ret <- ifelse(grepl("^gro_", model), "growth", ret)
  if (any(is.na(ret))) stop("Unknown demographic model")
  ret
}

get_mut_model <- function(model) {
  ret <- rep(NA, length(model))
  ret <- ifelse(grepl("fixed_num", model), "fixed_num", ret)
  ret <- ifelse(grepl("fixed_rate", model), "fixed_rate", ret)
  ret <- ifelse(grepl("variable_rate", model), "variable_rate", ret)
  if (any(is.na(ret))) stop("Unknown mutation model")
  ret
}

roc_values %>%
  group_by(model, stat, study, prob) %>%
  summarise(true_pos = mean(true_pos), false_pos = mean(false_pos)) %>%
  mutate(demography = get_demography(model), mutation = get_mut_model(model)) ->
  roc_values
  
roc_values[roc_values$prob == 0, c("true_pos", "false_pos")] <- c(0, 0)
    
save(roc_values, file = "cache/roc_values.Rda")
