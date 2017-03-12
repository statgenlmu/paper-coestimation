#!/usr/bin/Rscript --no-site-file

library(reshape2)
library(dplyr)
options(stringsAsFactors = FALSE)

# --------------------------------------------------------------
# Read data
# --------------------------------------------------------------
models <- list.dirs('./data', FALSE, FALSE)

get_model <- function(name) {
  if (grepl("^bottleneck_", name)) return("Bottleneck")
  if (grepl("^iwm_", name)) return("IWM")
  if (grepl("^growth_1pop_", name)) return("Growth (1POP)")
  if (grepl("^growth_", name)) return("Growth")
  stop("Failed to determine model")
}


get_sel <- function(name) {
  if (grepl("_neutr", name)) return("Neutral")
  if (grepl("_sel", name)) return("Selection")
  stop("Failed to determine selection status for ", name)
}


sim_data <- do.call(rbind, lapply(models, function(model) {
  message(model)
  model_name = get_model(model)
  sel = get_sel(model)

  test_data <- NULL
  load(paste0("./data/model_", tolower(model_name), ".Rda"))
  stopifnot(!is.null(test_data))
  true.values <- test_data$par_grid
  rm(test_data)

  estimates <- do.call(rbind, lapply(1:nrow(true.values), function(run) {
    result_file <- paste0("./data/", model,
                          "/logs/run_", run, "_result.Rda")
    
    if (!file.exists(result_file)) {
      message(" Missing Results: ", result_file)
      return(NA)
    }

    load(result_file)
    stopifnot(!is.null(result))
    result$estimate
  }))

  if (all(is.na(estimates))) return(NA)

  melt(true.values, as.is = TRUE) %>%
    right_join(melt(estimates, as.is = TRUE), by = c("Var1", "Var2")) %>%
    transmute(run = Var1, par = Var2,
              model = model_name, selection = sel, s = true.values[run, "s"],
              true_value = value.x, estimate = value.y) %>%
    filter(!is.na(estimate))
}))

sim_data$s <- as.factor(sim_data$s)

dir.create("cache", showWarnings = FALSE)
save(sim_data, file = 'cache//sim-data-df.Rda')
