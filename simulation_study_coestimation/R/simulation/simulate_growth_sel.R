#!/usr/bin/Rscript --no-site-file --no-save

library(coala)
library(jaatha)
library(testJaatha)

# Load the model
model_growth_sel <- NULL
load("cache/models.Rda")
stopifnot(!is.null(model_growth_sel))


# Create the test datasets
test_data_file <- "data/model_growth.Rda"
if (!file.exists(test_data_file)) {
  set.seed(1234321)
  test_data <- createTestData(model_growth_sel + sumstat_jsfs("jsfs_per_locus", per_locus = TRUE),
                              reps = 16,
                              grid.values = list(m = 1,
                                                 tau = .75,
                                                 theta = 1,
                                                 rho = 1,
                                                 s = c(.05, .1, .25, .5),
                                                 q = 10))
  save(test_data, file = test_data_file)
} else {
  load(test_data_file)
}


# Execute the jaatha analyses
testJaatha(dm = model_growth_sel, test_data = test_data,
           folder = "data/growth_sel", seed = 3303,
           sim = 160, init_method = "zoom-in",
           cores = c(8, 4))
