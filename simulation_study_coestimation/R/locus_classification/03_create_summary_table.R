#!/usr/bin/Rscript --no-site-file

library(stargazer)

class_stats <- readRDS("cache/class_stats.Rds")

model_labels <- list(iwm = "IWM",
                     bottleneck = "Bottleneck",
                     growth = "Growth")

class_stats$model <- sapply(class_stats$model,
                            function(x) model_labels[[x]])

colnames(class_stats) <- c("Method", "Model", "True Positives", "False Positives")


for (i in 3:4) {
  class_stats[ , i] <- sapply(class_stats[ , i], function(x) {
    paste0(format(x * 100, digits = 3), "%")
  })
}


# Create LaTeX-Tables
for (method in unique(class_stats$Method)) {
  tex_table <- stargazer(class_stats[class_stats$Method == method, 2:4],
                         summary=FALSE, rownames = FALSE, float = FALSE)
  write(tex_table, file = paste0("plots/class_stats_", method, ".tex"))
}
