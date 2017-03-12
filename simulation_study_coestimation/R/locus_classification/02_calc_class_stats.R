#!/usr/bin/Rscript --no-site-file

options(stringsAsFactors = FALSE)

models_overview <- read.csv("data/overview.csv")

class_stats <- do.call(rbind, lapply(which(models_overview$selection),
  function(i) {
    class_file <- file.path(models_overview$folder[i], "classified_loci.Rds")
    if (!file.exists(class_file)) {
      warning("No classification found for: ", models_overview$folder[i], call. = FALSE)
      return(data.frame())
    }
    class_results <- readRDS(class_file)
    model <- models_overview$model[i]

    do.call(rbind, lapply(class_results, function(class_result) {
      do.call(rbind, lapply(seq_along(class_result$estimated), function(i) {
        estimated <- class_result$estimated[[i]]
        data.frame(method = names(class_result$estimated)[i],
                   model = model,
                   selected = length(class_result$true),
                   neutral = 200 - length(class_result$true),
                   true_postive = sum(class_result$true %in% estimated),
                   false_postive = sum(!(estimated %in% class_result$true)))
      }))
    }))
  })
)


# Create the summary
library(dplyr)
class_stats %>%
  group_by(method, model) %>%
  summarise(percent_true_postive = sum(true_postive) / sum(selected),
            percent_false_positve = sum(false_postive) / sum(neutral)) ->
  class_stats_summary

print(class_stats_summary)
saveRDS(class_stats_summary, file = "cache/class_stats.Rds")
