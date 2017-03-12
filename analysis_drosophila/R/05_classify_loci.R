library(coala)
library(jaatha)
library(e1071)
library(parallel)

set.seed(1234)
svm_n_sim <- 200
cores <- parallel::detectCores()

# Load model & data
coala_model <- readRDS("cache/model.Rds") 
coala_model$sum_stats[["jsfs"]] <- sumstat_jsfs(per_locus = TRUE)

snp_data <- readRDS("cache/snp_data.Rds")
locus_list <- c(snp_data[[1]]$segsites, 
                snp_data[[2]]$segsites) 
sum_stats <- calc_sumstats_from_data(coala_model, locus_list)
rm(snp_data, locus_list)

result <- readRDS("cache/jaatha_results.Rda")


# Simulate learning data
data_learn <- simulate(coala_model, nsim = svm_n_sim,
                       pars = result$estimate,
                       cores = cores)
saveRDS(data_learn, "cache/svm_data_learn.Rds")

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
                   kernel = "radial", cost = cost_par, gamma = gamma_par, probability = TRUE)
saveRDS(svm_trained, "cache/svm_trained.Rds")


# Classify the observed loci
jsfs <- do.call(rbind, lapply(sum_stats$jsfs, as.vector))
real_data_df <- data.frame(locus = 1:nrow(jsfs),
                           jsfs = jsfs[, 2:(ncol(jsfs)-1)],
                           mcmf = sum_stats$mcmf)
select_loci <- real_data_df$locus[predict(svm_trained, real_data_df[, -1]) == "TRUE"]

message("Selected Loci: ")
print(select_loci)
saveRDS(select_loci, "cache/selected_loci.Rds")

class_probs <- predict(svm_trained, real_data_df[, -1], probability = TRUE)
saveRDS(class_probs, "cache/class_probs.Rds")
