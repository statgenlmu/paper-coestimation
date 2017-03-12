#!/usr/bin/Rscript --no-site-file --no-save

library(coala)
library(jaatha)
library(testJaatha)

# Load the model
model_bn_sel <- NULL
load("cache/models.Rda")
stopifnot(!is.null(model_bn_sel))

test_data_file <- "data/model_bottleneck.Rda"
if (!file.exists(test_data_file)) {
  set.seed(43562234)
  test_data <- createTestData(model_bn_sel + sumstat_jsfs("jsfs_per_locus", per_locus = TRUE),
                              reps = 16,
                              grid.values = list(m = 1,
                                                 tau = .75,
                                                 theta = 1,
                                                 rho = 1,
                                                 s = c(.05, .1, .25, .5),
                                                 q = 0.1))
  save(test_data, file = test_data_file)
} else {
  load(test_data_file)
}

testJaatha(dm = model_bn_sel, test_data = test_data,
           folder = "data/bottleneck_sel", seed = 5505,
           sim = 160, init_method = "zoom-in")
