#!/usr/bin/Rscript --no-site-file

library(jaatha)
library(coala)
library(doMC)
library(parallel)

# Number of cores used for procssing the results
cores_runs <- ceiling(detectCores() / 4)
# Number of cores used for the classification of one run
cores_classifiction <- 4

options(stringsAsFactors = FALSE)

# Get information about the model for which we do the classifications
models_overview <- read.csv("data/overview.csv")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(sum(args[1] == models_overview$folder) == 1)
model_info <- models_overview[args[1] == models_overview$folder, ]
message("Classifying loci for ", model_info$folder)

# Load the datasets
test_data <- NULL
load(model_info$data)
stopifnot(!is.null(test_data))

# Load the coala model
load("cache/models.Rda")
coala_model <- get(model_info$coala_model)

# Function to get the true selected loci
get_sel_loci <- function(dataset) {
  sim_cmds <- dataset$cmds
  sel_group <- which(grepl("-Sc 0 2 1000", sim_cmds))
  loci_per_group <- sapply(strsplit(sim_cmds, " "), function(x) {
    as.numeric(x[3])
  })
  if (sel_group == 1) {
    sel_true <- 1:loci_per_group[1]
  } else {
    sel_true <- 1:loci_per_group[2] + loci_per_group[1]
  }
}

# The functions used for the classification
classify_functions <- list()
classify_functions$mcmf_outlier <- function(dataset, result,
                                            jaatha_model, jaatha_data,
                                            coala_model, cores) {
  mcmf <- dataset$mcmf
  # Return the most positive loci
  cutoff <- quantile(mcmf, 1 - result$estimate["s"])
  which(mcmf >= cutoff)
}

classify_svm <- NULL
source("R/locus_classification/classify_svm.R")
stopifnot((!is.null(classify_svm)))
classify_functions$svm <- classify_svm

# Process the runs in parallel
classification <- mclapply(seq(along = test_data$data), function(run) {
  message("Run ", run, " of ", length(test_data$data))

  # Load the results of the run
  jaatha_data <- NULL
  jaatha_model <- NULL
  result <- NULL
  setup_file <- file.path(model_info$folder, "logs", paste0("run_", run, "_setup.Rda"))
  result_file <- file.path(model_info$folder, "logs", paste0("run_", run, "_result.Rda"))
  if (!(file.exists(result_file) && file.exists(setup_file))) {
    warning("No results for run ", run, " found")
    return(NA)
  }
  load(setup_file)
  load(result_file)
  stopifnot(!is.null(jaatha_data))
  stopifnot(!is.null(jaatha_model))
  stopifnot(!is.null(result))

  # Estimate selection loci using the classify functions
  sel_estimated <- lapply(classify_functions, function(x) {
    x(dataset = test_data$data[[run]],
      result = result,
      jaatha_model = jaatha_model,
      jaatha_data = jaatha_data,
      coala_model = coala_model,
      cores = cores_classifiction)
  })

  # Get true selected loci
  sel_true <- get_sel_loci(test_data$data[[run]])

  list(estimated = sel_estimated, true = sel_true)
}, mc.cores = cores_runs, mc.preschedule = FALSE)

saveRDS(classification, file = file.path(model_info$folder, "classified_loci.Rds"))
