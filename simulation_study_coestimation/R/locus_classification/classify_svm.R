#!/usr/bin/Rscript --no-site-file

library(coala)
library(jaatha)
library(e1071)
library(parallel)

svm_n_sim <- 100

shuffle_rows <- function(dataset) {
  dataset[sample.int(nrow(dataset)), , drop = FALSE]
}

classify_svm <- function(dataset, result,
                         jaatha_model, jaatha_data,
                         coala_model, cores) {

  # Simulate learning data
  coala_model$sum_stats <- list()
  coala_model <- coala_model +
    sumstat_jsfs(per_locus = TRUE) +
    sumstat_mcmf(population = 2)

  data_learn <- simulate(coala_model, nsim = svm_n_sim,
                         pars = result$estimate,
                         cores = cores)


  # Convert it into a data.frame
  data_learn_df <- do.call(rbind, lapply(data_learn, function(x) {
    sel_loci <- get_sel_loci(x)
    selected <- seq_along(x$mcmf) %in% sel_loci
    jsfs <- do.call(rbind, lapply(x$jsfs, as.vector))
    data.frame(selected = selected, jsfs = jsfs[, 2:(ncol(jsfs)-1)], mcmf = x$mcmf)
  }))
  rm(data_learn)

  # SVM's parameters optimized for one test dataset
  cost_par <- 3
  gamma_par <- 1e-3

  # Train the svm on simulated data
  svm_trained <- svm(x = data_learn_df[, -1], y = as.factor(data_learn_df[, 1]),
                     kernel = "radial", cost = cost_par, gamma = gamma_par)

  # Classify the observed loci
  jsfs <- do.call(rbind, lapply(dataset$jsfs_per_locus, as.vector))
  real_data_df <- data.frame(locus = 1:nrow(jsfs),
                             jsfs = jsfs[, 2:(ncol(jsfs)-1)],
                             mcmf = dataset$mcmf)
  real_data_df <- shuffle_rows(real_data_df)
  real_data_df$locus[predict(svm_trained, real_data_df[, -1]) == "TRUE"]
}
